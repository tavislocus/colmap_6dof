// Copyright (c), ETH Zurich and UNC Chapel Hill.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the name of ETH Zurich and UNC Chapel Hill nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include "colmap/estimators/similarity_transform.h"

namespace colmap {
namespace {

template <bool kEstimateScale>
inline bool EstimateRigidOrSim3d(const std::vector<Eigen::Vector3d>& src,
                                 const std::vector<Eigen::Vector3d>& tgt,
                                 Eigen::Matrix3x4d& tgt_from_src) {
  std::vector<Eigen::Matrix3x4d> models;
  SimilarityTransformEstimator<3, kEstimateScale>().Estimate(src, tgt, &models);
  if (models.empty()) {
    return false;
  }
  THROW_CHECK_EQ(models.size(), 1);
  tgt_from_src = models[0];
  return true;
}

template <bool kEstimateScale>
inline typename RANSAC<SimilarityTransformEstimator<3, kEstimateScale>>::Report
EstimateRigidOrSim3dRobust(const std::vector<Eigen::Vector3d>& src,
                           const std::vector<Eigen::Vector3d>& tgt,
                           const RANSACOptions& options,
                           Eigen::Matrix3x4d& tgt_from_src) {
  LORANSAC<SimilarityTransformEstimator<3, kEstimateScale>,
           SimilarityTransformEstimator<3, kEstimateScale>>
      ransac(options);
  auto report = ransac.Estimate(src, tgt);
  if (report.success) {
    tgt_from_src = report.model;
  }
  return report;
}

}  // namespace

// Robust and non-robust estimators for similarity with rotation
namespace {

template <bool kEstimateScale>
inline bool EstimateRigidOrSim3dWithRotation(
    const std::vector<Eigen::Vector3d>& src_pos,
    const std::vector<Eigen::Vector3d>& tgt_pos,
    const std::vector<Eigen::Quaterniond>& src_rot,
    const std::vector<Eigen::Quaterniond>& tgt_rot,
    typename SimilarityPoseEstimator<3, kEstimateScale>::M_t& tgt_from_src) {
  using Estimator = SimilarityPoseEstimator<3, kEstimateScale>;
  std::vector<typename Estimator::M_t> models;
  Estimator::Estimate(src_pos, tgt_pos, src_rot, tgt_rot, &models);
  if (models.empty()) {
    return false;
  }
  THROW_CHECK_EQ(models.size(), 1);
  tgt_from_src = models[0];
  return true;
}

template <bool kEstimateScale>
inline typename RANSAC<SimilarityPoseEstimator<3, kEstimateScale>>::Report
EstimateRigidOrSim3dWithRotationRobust(
    const std::vector<Eigen::Vector3d>& src_pos,
    const std::vector<Eigen::Vector3d>& tgt_pos,
    const std::vector<Eigen::Quaterniond>& src_rot,
    const std::vector<Eigen::Quaterniond>& tgt_rot,
    const RANSACOptions& options,
    typename SimilarityPoseEstimator<3, kEstimateScale>::M_t& tgt_from_src) {
  using Estimator = SimilarityPoseEstimator<3, kEstimateScale>;
  LORANSAC<Estimator, Estimator> ransac(options);
  auto report = ransac.Estimate(src_pos, tgt_pos, src_rot, tgt_rot);
  if (report.success) {
    tgt_from_src = report.model;
  }
  return report;
}

}  // namespace


// Position
bool EstimateRigid3d(const std::vector<Eigen::Vector3d>& src,
                     const std::vector<Eigen::Vector3d>& tgt,
                     Rigid3d& tgt_from_src) {
  Eigen::Matrix3x4d tgt_from_src_mat = Eigen::Matrix3x4d::Zero();
  if (!EstimateRigidOrSim3d<false>(src, tgt, tgt_from_src_mat)) {
    return false;
  }
  tgt_from_src = Rigid3d::FromMatrix(tgt_from_src_mat);
  return true;
}

typename RANSAC<SimilarityTransformEstimator<3, false>>::Report
EstimateRigid3dRobust(const std::vector<Eigen::Vector3d>& src,
                      const std::vector<Eigen::Vector3d>& tgt,
                      const RANSACOptions& options,
                      Rigid3d& tgt_from_src) {
  Eigen::Matrix3x4d tgt_from_src_mat = Eigen::Matrix3x4d::Zero();
  auto report =
      EstimateRigidOrSim3dRobust<false>(src, tgt, options, tgt_from_src_mat);
  tgt_from_src = Rigid3d::FromMatrix(tgt_from_src_mat);
  return report;
}

bool EstimateSim3d(const std::vector<Eigen::Vector3d>& src,
                   const std::vector<Eigen::Vector3d>& tgt,
                   Sim3d& tgt_from_src) {
  Eigen::Matrix3x4d tgt_from_src_mat = Eigen::Matrix3x4d::Zero();
  if (!EstimateRigidOrSim3d<true>(src, tgt, tgt_from_src_mat)) {
    return false;
  }
  tgt_from_src = Sim3d::FromMatrix(tgt_from_src_mat);
  return true;
}

typename RANSAC<SimilarityTransformEstimator<3, true>>::Report
EstimateSim3dRobust(const std::vector<Eigen::Vector3d>& src,
                    const std::vector<Eigen::Vector3d>& tgt,
                    const RANSACOptions& options,
                    Sim3d& tgt_from_src) {
  Eigen::Matrix3x4d tgt_from_src_mat = Eigen::Matrix3x4d::Zero();
  auto report =
      EstimateRigidOrSim3dRobust<true>(src, tgt, options, tgt_from_src_mat);
  tgt_from_src = Sim3d::FromMatrix(tgt_from_src_mat);
  return report;
}


// Rotation
// For rigid transform (no scale)
bool EstimateRigid3dWithRotation(const std::vector<Eigen::Vector3d>& src_pos,
                                 const std::vector<Eigen::Vector3d>& tgt_pos,
                                 const std::vector<Eigen::Quaterniond>& src_rot,
                                 const std::vector<Eigen::Quaterniond>& tgt_rot,
                                 Rigid3d& tgt_from_src) {
  SimilarityPoseEstimator<3, false>::M_t model;
  if (!EstimateRigidOrSim3dWithRotation<false>(src_pos, tgt_pos, src_rot, tgt_rot, model)) {
    return false;
  }
  tgt_from_src = Rigid3d(model.rotation, model.translation);
  return true;
}

typename RANSAC<SimilarityPoseEstimator<3, false>>::Report
EstimateRigid3dRobustWithRotation(const std::vector<Eigen::Vector3d>& src_pos,
                                  const std::vector<Eigen::Vector3d>& tgt_pos,
                                  const std::vector<Eigen::Quaterniond>& src_rot,
                                  const std::vector<Eigen::Quaterniond>& tgt_rot,
                                  const RANSACOptions& options,
                                  Rigid3d& tgt_from_src) {
  SimilarityPoseEstimator<3, false>::M_t model;
  auto report = EstimateRigidOrSim3dWithRotationRobust<false>(
      src_pos, tgt_pos, src_rot, tgt_rot, options, model);
  tgt_from_src = Rigid3d(model.rotation, model.translation);
  return report;
}

// For similarity transform (with scale)
bool EstimateSim3dWithRotation(const std::vector<Eigen::Vector3d>& src_pos,
                               const std::vector<Eigen::Vector3d>& tgt_pos,
                               const std::vector<Eigen::Quaterniond>& src_rot,
                               const std::vector<Eigen::Quaterniond>& tgt_rot,
                               Sim3d& tgt_from_src) {
  SimilarityPoseEstimator<3, true>::M_t model;
  if (!EstimateRigidOrSim3dWithRotation<true>(src_pos, tgt_pos, src_rot, tgt_rot, model)) {
    return false;
  }
  tgt_from_src = Sim3d(model.scale, model.rotation, model.translation);
  return true;
}

typename RANSAC<SimilarityPoseEstimator<3, true>>::Report
EstimateSim3dRobustWithRotation(const std::vector<Eigen::Vector3d>& src_pos,
                                const std::vector<Eigen::Vector3d>& tgt_pos,
                                const std::vector<Eigen::Quaterniond>& src_rot,
                                const std::vector<Eigen::Quaterniond>& tgt_rot,
                                const RANSACOptions& options,
                                Sim3d& tgt_from_src) {
  SimilarityPoseEstimator<3, true>::M_t model;
  auto report = EstimateRigidOrSim3dWithRotationRobust<true>(
      src_pos, tgt_pos, src_rot, tgt_rot, options, model);
  tgt_from_src = Sim3d(model.scale, model.rotation, model.translation);
  return report;
}




















}  // namespace colmap
