% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4PRRR2
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4PRRR2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR2_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR2_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR2_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S4PRRR2_invdynm_fixb_reg2_snew_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:27:27
% EndTime: 2019-07-18 13:27:28
% DurationCPUTime: 1.22s
% Computational Cost: add. (2518->144), mult. (3473->125), div. (0->0), fcn. (2310->6), ass. (0->84)
t241 = sin(qJ(2));
t244 = cos(qJ(2));
t238 = g(2) + qJDD(1);
t239 = sin(qJ(4));
t240 = sin(qJ(3));
t242 = cos(qJ(4));
t243 = cos(qJ(3));
t249 = (t239 * t243 + t240 * t242) * t238;
t276 = (t239 * t240 - t242 * t243) * t238;
t287 = -t241 * t276 + t244 * t249;
t237 = qJD(2) + qJD(3);
t232 = qJD(4) + t237;
t230 = t232 ^ 2;
t236 = qJDD(2) + qJDD(3);
t231 = qJDD(4) + t236;
t208 = t242 * t230 + t239 * t231;
t211 = t239 * t230 - t242 * t231;
t183 = t243 * t208 - t240 * t211;
t187 = t240 * t208 + t243 * t211;
t293 = t241 * t183 + t244 * t187;
t296 = qJ(1) * t293 + t287;
t169 = t244 * t183 - t241 * t187;
t256 = t241 * t249 + t244 * t276;
t295 = qJ(1) * t169 - t256;
t248 = (t244 * t240 + t241 * t243) * t238;
t235 = t237 ^ 2;
t215 = t243 * t235 + t240 * t236;
t218 = t240 * t235 - t243 * t236;
t286 = t241 * t215 + t244 * t218;
t294 = qJ(1) * t286 + t248;
t225 = t244 * g(1) + t241 * g(3);
t220 = qJDD(2) * pkin(1) + t225;
t226 = -t241 * g(1) + t244 * g(3);
t245 = qJD(2) ^ 2;
t221 = -t245 * pkin(1) - t226;
t197 = -t243 * t220 + t240 * t221;
t190 = t236 * pkin(2) - t197;
t198 = t240 * t220 + t243 * t221;
t191 = -t235 * pkin(2) + t198;
t174 = -t242 * t190 + t239 * t191;
t175 = t239 * t190 + t242 * t191;
t163 = t242 * t174 - t239 * t175;
t264 = t239 * t174 + t242 * t175;
t159 = t243 * t163 - t240 * t264;
t289 = t240 * t163 + t243 * t264;
t154 = t241 * t159 + t244 * t289;
t155 = t244 * t159 - t241 * t289;
t192 = t244 * t215 - t241 * t218;
t290 = qJ(1) * t192;
t179 = t243 * t197 - t240 * t198;
t263 = t240 * t197 + t243 * t198;
t165 = t241 * t179 + t244 * t263;
t166 = t244 * t179 - t241 * t263;
t277 = (t240 * t241 - t243 * t244) * t238;
t222 = t241 * qJDD(2) + t244 * t245;
t229 = t244 * t238;
t275 = qJ(1) * t222 + t229;
t274 = pkin(2) * t163;
t272 = t239 * t238;
t271 = t240 * t238;
t270 = t241 * t238;
t269 = t242 * t238;
t268 = t243 * t238;
t267 = pkin(2) * t271;
t266 = pkin(1) * t229;
t265 = t240 * t270;
t219 = (-pkin(2) * t243 - pkin(1)) * t238;
t259 = -t241 * t219 + t244 * t267;
t258 = t243 * t229 - t265;
t223 = -t244 * qJDD(2) + t241 * t245;
t257 = qJ(1) * t223 + t270;
t200 = -t244 * t225 + t241 * t226;
t199 = -t241 * t225 - t244 * t226;
t250 = -pkin(2) * t211 - t174;
t247 = pkin(2) * t265 + t244 * t219;
t246 = -pkin(2) * t208 - t175;
t228 = pkin(1) * t270;
t182 = -pkin(1) * t215 - t198;
t181 = -pkin(1) * t218 - t197;
t176 = pkin(1) * t179;
t168 = -pkin(1) * t183 + t246;
t167 = -pkin(1) * t187 + t250;
t156 = -pkin(1) * t159 - t274;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, t238, 0, g(3), qJ(1) * g(3), 0, 0, -t222, 0, t223, 0, t275, -t257, -t199, -qJ(1) * t199, 0, 0, -t192, 0, t286, 0, -t277 + t290, -t294, -t165, -qJ(1) * t165 + t266, 0, 0, -t169, 0, t293, 0, t295, -t296, -t154, -qJ(1) * t154 - t247; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, -g(1), -g(3), 0, 0, 0, 0, 0, 0, 0, -qJDD(2), -t225, -t226, 0, 0, 0, 0, 0, 0, 0, -t236, -t181, -t182, 0, t176, 0, 0, 0, 0, 0, -t231, -t167, -t168, 0, -t156; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, -g(1), -qJ(1) * g(1), 0, 0, -t223, 0, -t222, 0, t257, t275, t200, qJ(1) * t200, 0, 0, -t286, 0, -t192, 0, t294, t258 + t290, t166, qJ(1) * t166 + t228, 0, 0, -t293, 0, -t169, 0, t296, t295, t155, qJ(1) * t155 + t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, -g(1), 0, 0, 0, -t223, 0, -t222, 0, t270, t229, t200, 0, 0, 0, -t286, 0, -t192, 0, t248, t258, t166, t228, 0, 0, -t293, 0, -t169, 0, t287, -t256, t155, t259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, 0, -g(3), 0, 0, 0, t222, 0, -t223, 0, -t229, t270, t199, 0, 0, 0, t192, 0, -t286, 0, t277, t248, t165, -t266, 0, 0, t169, 0, -t293, 0, t256, t287, t154, t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1), g(3), 0, 0, 0, 0, 0, 0, 0, qJDD(2), t225, t226, 0, 0, 0, 0, 0, 0, 0, t236, t181, t182, 0, -t176, 0, 0, 0, 0, 0, t231, t167, t168, 0, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t245, 0, 0, t238, -t225, 0, 0, 0, -t218, 0, -t215, 0, t271, t268, t179, 0, 0, 0, -t187, 0, -t183, 0, t249, -t276, t159, t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, 0, qJDD(2), 0, -t238, 0, -t226, 0, 0, 0, t215, 0, -t218, 0, -t268, t271, t263, -pkin(1) * t238, 0, 0, t183, 0, -t187, 0, t276, t249, t289, t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t225, t226, 0, 0, 0, 0, 0, 0, 0, t236, t181, t182, 0, -t176, 0, 0, 0, 0, 0, t231, t167, t168, 0, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, 0, -t235, 0, 0, t238, t197, 0, 0, 0, -t211, 0, -t208, 0, t272, t269, t163, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, 0, t236, 0, -t238, 0, t198, 0, 0, 0, t208, 0, -t211, 0, -t269, t272, t264, -pkin(2) * t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236, -t197, -t198, 0, 0, 0, 0, 0, 0, 0, t231, t250, t246, 0, -t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, 0, -t230, 0, 0, t238, t174, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, 0, t231, 0, -t238, 0, t175, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t231, -t174, -t175, 0, 0;];
m_new_reg  = t1;
