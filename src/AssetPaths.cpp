#include <FinancialMathematics/AssetPaths.h>
#include <FinancialMathematics/brownian_motion.h>

namespace FinancialMathematics {
  AssetPaths::AssetPaths(std::vector<std::vector<double>> x) : x_(std::move(x)) {}

  std::ostream &operator<<(std::ostream &os, const AssetPaths &asset_paths) {
    if (asset_paths.paths() < 10) {
      for (int i = 0; i < asset_paths.paths(); ++i) {
        os << "Path " << i << ": ";
        if(asset_paths.steps() < 10){
          for (int j = 0; j < asset_paths.steps(); ++j) {
            os << asset_paths.x_[i][j] << " ";
          }
        } else {
          for (int j = 0; j < 3; ++j) {
            os << asset_paths.x_[i][j] << " ";
          }
          os << "... ";
          for (int j = asset_paths.steps() - 3; j < asset_paths.steps(); ++j) {
            os << asset_paths.x_[i][j] << " ";
          }
        }
        os << std::endl;
      }
    } else {
        for (int i = 0; i < 3; ++i) {
            os << "Path " << i << ": ";
            if(asset_paths.steps() < 10){
            for (int j = 0; j < asset_paths.steps(); ++j) {
                os << asset_paths.x_[i][j] << " ";
            }
            } else {
            for (int j = 0; j < 3; ++j) {
                os << asset_paths.x_[i][j] << " ";
            }
            os << "... ";
            for (int j = asset_paths.steps() - 3; j < asset_paths.steps(); ++j) {
                os << asset_paths.x_[i][j] << " ";
            }
            }
            os << std::endl;
        }
        os << "..." << std::endl;
        for (int i = asset_paths.paths() - 3; i < asset_paths.paths(); ++i) {
            os << "Path " << i << ": ";
            if(asset_paths.steps() < 10){
            for (int j = 0; j < asset_paths.steps(); ++j) {
                os << asset_paths.x_[i][j] << " ";
            }
            } else {
            for (int j = 0; j < 3; ++j) {
                os << asset_paths.x_[i][j] << " ";
            }
            os << "... ";
            for (int j = asset_paths.steps() - 3; j < asset_paths.steps(); ++j) {
                os << asset_paths.x_[i][j] << " ";
            }
            }
            os << std::endl;
        }
    }
    return os;
  }
}