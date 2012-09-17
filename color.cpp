//===- color.cpp - Color matching using boost::units ------------*- C++ -*-===//
//
// Color
//
// This file is distributed under the Simplified BSD License. See LICENSE.TXT
// for details.
//
//===----------------------------------------------------------------------===//
///
/// \file
/// \brief This file does color matching in CIE L*C*h space using CMC (1984).
///
//===----------------------------------------------------------------------===//

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <cmath>
#include <fstream>
#include <iostream>

using namespace boost::numeric::ublas;

matrix<double> sRGB_D65_XYZ() {
  matrix<double> ret(3, 3);
  ret(0, 0) = 0.4124564; ret(0, 1) = 0.3575761; ret(0, 2) = 0.1804375;
  ret(1, 0) = 0.2126729; ret(1, 1) = 0.7151522; ret(1, 2) = 0.0721750;
  ret(2, 0) = 0.0193339; ret(2, 1) = 0.1191920; ret(2, 2) = 0.9503041;
  return ret;
}

double sRGBInverseCompand(double V) {
  if (V > 0.04045)
    return pow((V + 0.055) / 1.055, 2.4);
  return V / 12.92;
}

struct sRGB {
  double R, G, B;
};

struct CIEXYZ10degD65 {
  double X, Y, Z;
};

struct CIELab10degD65 {
  double L, a, b;
};

/// \brief CIE L*C*h with CIE 10 degree standard observer and D65 standard
/// illuminant.
struct CIELCh10degD65 {
  double L, C, h;
};

CIEXYZ10degD65 D65() {
  CIEXYZ10degD65 ret = {.95047, 1.00, 1.08883};
  return ret;
}

CIEXYZ10degD65 sRGBToCIEXYZ10degD65(const sRGB &color) {
  // sRGB to linear RGB.
  vector<double> rgb(3);
  rgb(0) = sRGBInverseCompand(color.R);
  rgb(1) = sRGBInverseCompand(color.G);
  rgb(2) = sRGBInverseCompand(color.B);

  // linear RGB times the transformation matrix.
  vector<double> XYZ = prod(sRGB_D65_XYZ(), rgb);

  CIEXYZ10degD65 ret = {XYZ(0), XYZ(1), XYZ(2)};
  return ret;
}

CIELab10degD65 CIEXYZ10degD65ToCIELab10degD65(const CIEXYZ10degD65 &color) {
  double e = 216. / 24389.;
  double k = 24389. / 27.;
  double xr = color.X / D65().X;
  double yr = color.Y / D65().Y;
  double zr = color.Z / D65().Z;

  auto f = [k, e](double r) -> double {
    if (r > e)
      return std::pow(r, 1. / 3.);
    return (k * r + 16.) / 116.;
  };

  CIELab10degD65 ret;
  ret.L = 116. * f(yr) - 16;
  ret.a = 500 * (f(xr) - f(yr));
  ret.b = 200 * (f(yr) - f(zr));
  return ret;
}

CIELCh10degD65 CIELab10degD65ToCIELCh10degD65(const CIELab10degD65 &Lab) {
  CIELCh10degD65 ret;
  ret.L = Lab.L;
  ret.C = std::pow(std::pow(Lab.a, 2) + std::pow(Lab.b, 2), 1. / 2.);
  ret.h = std::atan2(Lab.b, Lab.a) * 57.2957795;
  if (ret.h < 0)
    ret.h += 360;
  else if (ret.h >= 360)
    ret.h -= 360;
  return ret;
}

CIELCh10degD65 convert(sRGB other) {
  CIEXYZ10degD65 XYZ = sRGBToCIEXYZ10degD65(other);
  CIELab10degD65 Lab = CIEXYZ10degD65ToCIELab10degD65(XYZ);
  return CIELab10degD65ToCIELCh10degD65(Lab);
}

double cmc_distance(CIELCh10degD65 a, CIELab10degD65 al, CIELCh10degD65 b, CIELab10degD65 bl, double l, double c) {
  double deltaL = a.L - b.L;
  double deltaC = a.C - b.C;
  double deltaH = std::pow(std::pow(al.a - bl.a, 2) +
                           std::pow(al.b - bl.b, 2) -
                           std::pow(deltaC, 2), 0.5);
  double Sl = a.L < 16.
            ? 0.511
            : (0.040975 * a.L) / (1. + 0.01765 * a.L);
  double Sc = ((0.0638 * a.C) / (1 + 0.0131 * a.C)) + 0.638;
  double T = 164. <= a.h && a.h <= 345.
           ? 0.56 + std::abs(0.2 * std::cos(a.h + 168.))
           : 0.36 + std::abs(0.4 * std::cos(a.h + 35.));
  double C4 = std::pow(a.C, 4);
  double F = std::pow(C4 / (C4 + 1900), 1. / 2.);
  double Sh = Sc * (F * T + 1 - F);

  double deltaLcmc = deltaL / (l * Sl);
  double deltaCcmc = deltaC / (c * Sc);
  double deltaHcmc = deltaH / Sh;

  return std::pow(std::pow(deltaLcmc, 2) +
                  std::pow(deltaCcmc, 2) +
                  std::pow(deltaHcmc, 2), 0.5);
}

int main(int argc, const char **argv) {
  sRGB in;
  std::cout << "sRGB r g b: ";
  std::cin >> in.R >> in.G >> in.B;
  in.R /= 255.;
  in.G /= 255.;
  in.B /= 255.;
  CIEXYZ10degD65 inXYZ = sRGBToCIEXYZ10degD65(in);
  CIELab10degD65 inLab = CIEXYZ10degD65ToCIELab10degD65(inXYZ);
  CIELCh10degD65 inLCh = CIELab10degD65ToCIELCh10degD65(inLab);
  std::cout << "L*=" << inLCh.L << " C*=" << inLCh.C << " h=" << inLCh.h << "\n";

  if (argc < 2)
    return 0;

  // For each thread color, try to get a match.
  std::string bestNum, bestName, bestFileName;
  double bestDeltaE = 1000000.;
  CIELCh10degD65 bestLCh;
  sRGB bestColor;

  auto testColor = [&](boost::filesystem::path path) {
    std::ifstream colorFile(path.string());
    // Skip the first line.
    std::string line;
    if (!std::getline(colorFile, line))
      return;

    while (colorFile) {
      std::string num, r, g, b, name;
      std::getline(colorFile, num, ',');
      std::getline(colorFile, r, ',');
      std::getline(colorFile, g, ',');
      std::getline(colorFile, b, ',');
      std::getline(colorFile, name);
      if (r.empty() || g.empty() || b.empty())
        break;
      sRGB color = { boost::lexical_cast<double>(r) / 255.
                    , boost::lexical_cast<double>(g) / 255.
                    , boost::lexical_cast<double>(b) / 255.};
      CIEXYZ10degD65 XYZ = sRGBToCIEXYZ10degD65(color);
      CIELab10degD65 Lab = CIEXYZ10degD65ToCIELab10degD65(XYZ);
      CIELCh10degD65 LCh = CIELab10degD65ToCIELCh10degD65(Lab);
      double deltaE = cmc_distance(inLCh, inLab, LCh, Lab, 2, 1);
      if (deltaE < bestDeltaE) {
        bestDeltaE = deltaE;
        bestNum = num;
        bestName = name;
        bestLCh = LCh;
        bestColor = color;
        bestFileName = path.filename().string();
      }
    }
  };

  try {
    if (boost::filesystem::is_directory(argv[1])) {
      boost::filesystem::directory_iterator di(argv[1]);
      boost::filesystem::directory_iterator de;
      for (; di != de; ++di) {
        if (di->path().extension() != ".thr")
          continue;
        testColor(di->path());
      }
    } else {
      testColor(boost::filesystem::path(argv[1]));
    }
  } catch (const std::exception &e) {
    std::cerr << "Blarg!: " << e.what() << "\n";
  } catch (...) {
    std::cerr << "Unknown exception\n";
  }

  std::cout << "Nearest match: " << bestDeltaE << " " << bestNum << " " << bestName << " " << bestFileName << "\n";
  std::cout << "L*=" << bestLCh.L << " C*=" << bestLCh.C << " h=" << bestLCh.h << "\n";
  std::cout << "sR=" << bestColor.R << " sG=" << bestColor.G << " sB=" << bestColor.B << "\n";
}
