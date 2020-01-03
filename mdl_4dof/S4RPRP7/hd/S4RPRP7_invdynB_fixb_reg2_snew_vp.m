% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RPRP7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynB_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:18
% EndTime: 2019-12-31 16:47:21
% DurationCPUTime: 1.37s
% Computational Cost: add. (1743->218), mult. (3607->252), div. (0->0), fcn. (1705->4), ass. (0->156)
t226 = qJD(3) ^ 2;
t224 = cos(qJ(3));
t220 = t224 ^ 2;
t227 = qJD(1) ^ 2;
t261 = t220 * t227;
t204 = t226 + t261;
t222 = sin(qJ(3));
t251 = t222 * t224 * t227;
t198 = qJDD(3) + t251;
t266 = t198 * t222;
t160 = t204 * t224 + t266;
t258 = qJD(1) * qJD(3);
t248 = t222 * t258;
t254 = qJDD(1) * t224;
t192 = -0.2e1 * t248 + t254;
t223 = sin(qJ(1));
t225 = cos(qJ(1));
t132 = t160 * t225 + t192 * t223;
t295 = pkin(4) * t132;
t137 = t160 * t223 - t192 * t225;
t294 = pkin(4) * t137;
t285 = -pkin(5) - pkin(1);
t293 = t285 * t160;
t265 = t198 * t224;
t165 = -t204 * t222 + t265;
t292 = t285 * t165;
t247 = t224 * t258;
t256 = qJDD(1) * t222;
t189 = 0.2e1 * t247 + t256;
t145 = t189 * t222 - t192 * t224;
t219 = t222 ^ 2;
t197 = (-t219 + t220) * t227;
t291 = t145 * t223 - t197 * t225;
t290 = t145 * t225 + t197 * t223;
t262 = t219 * t227;
t202 = -t226 + t262;
t159 = t202 * t222 + t265;
t253 = qJDD(1) * t225;
t289 = t159 * t223 - t222 * t253;
t255 = qJDD(1) * t223;
t288 = t159 * t225 + t222 * t255;
t287 = pkin(2) * t160 - qJ(2) * t165;
t203 = -t226 - t262;
t199 = qJDD(3) - t251;
t263 = t199 * t224;
t158 = t203 * t222 + t263;
t264 = t199 * t222;
t163 = t203 * t224 - t264;
t286 = pkin(2) * t158 - qJ(2) * t163;
t260 = t219 + t220;
t196 = t260 * t227;
t284 = pkin(2) * t196;
t190 = -t247 - t256;
t283 = pkin(3) * t190;
t282 = pkin(3) * t222;
t281 = pkin(3) * t224;
t131 = -t158 * t225 + t189 * t223;
t280 = pkin(4) * t131;
t193 = t260 * qJDD(1);
t269 = t193 * t225;
t149 = -t196 * t223 + t269;
t279 = pkin(4) * t149;
t277 = qJ(4) * t222;
t276 = qJDD(1) * pkin(1);
t201 = g(1) * t225 + g(2) * t223;
t217 = qJDD(1) * qJ(2);
t234 = t201 - t217;
t230 = -t285 * t227 + t234;
t257 = (qJD(2) * qJD(1));
t171 = t230 - (2 * t257);
t275 = t171 * t222;
t274 = t171 * t224;
t273 = t189 * t224;
t270 = t193 * t223;
t200 = g(1) * t223 - t225 * g(2);
t239 = qJDD(2) - t200;
t233 = -qJ(2) * t227 + t239;
t176 = t285 * qJDD(1) + t233;
t147 = t222 * g(3) + t224 * t176;
t237 = -qJ(4) * t224 + t282;
t259 = t227 * t237;
t252 = t285 * t163;
t246 = g(3) * t224 - t222 * t176;
t213 = 2 * t257;
t177 = -pkin(1) * t227 + t213 - t234;
t178 = -t233 + t276;
t134 = t225 * t177 - t178 * t223;
t154 = -t200 * t223 - t225 * t201;
t245 = t223 * t251;
t244 = t225 * t251;
t243 = qJ(2) * t189 + t285 * t158;
t242 = -qJ(2) * t196 - t285 * t193;
t194 = -t223 * t227 + t253;
t241 = pkin(4) * t194 + g(3) * t223;
t195 = t225 * t227 + t255;
t240 = -pkin(4) * t195 + g(3) * t225;
t238 = -t277 - t281;
t125 = t147 * t224 - t222 * t246;
t126 = -t147 * t222 - t224 * t246;
t130 = t177 * t223 + t178 * t225;
t236 = t192 * t222 + t273;
t235 = -t202 * t224 + t266;
t153 = t200 * t225 - t201 * t223;
t232 = qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) - t222 * t259 - t246;
t231 = -qJDD(3) * pkin(3) - t226 * qJ(4) + t224 * t259 + qJDD(4) - t147;
t229 = 0.2e1 * (qJD(4) * t224 - qJD(2)) * qJD(1) - qJ(4) * t248 - pkin(3) * t247 + t230;
t191 = -t248 + t254;
t228 = qJ(4) * t191 + t229;
t205 = t226 - t261;
t184 = t260 * t258;
t175 = qJDD(3) * t225 - t184 * t223;
t174 = qJDD(3) * t223 + t184 * t225;
t173 = -t191 * t222 - t220 * t258;
t172 = -t190 * t224 - t219 * t258;
t166 = -t205 * t222 + t263;
t164 = (t191 - t248) * t224;
t161 = -t205 * t224 - t264;
t156 = (-t190 + t247) * t222;
t155 = -pkin(2) * t193 + t238 * qJDD(1);
t150 = -t196 * t225 - t270;
t146 = pkin(4) * t150;
t143 = -t172 * t223 - t244;
t142 = -t173 * t223 + t244;
t141 = t172 * t225 - t245;
t140 = t173 * t225 + t245;
t139 = -t161 * t223 + t224 * t253;
t138 = t161 * t225 + t223 * t254;
t135 = t158 * t223 + t189 * t225;
t129 = pkin(4) * t135;
t127 = -pkin(3) * t226 + t232;
t124 = qJ(4) * t196 + t231;
t123 = (t196 - t226) * pkin(3) + t232;
t122 = t228 + t283;
t121 = -t126 - t284;
t120 = t246 - t287;
t119 = t147 + t286;
t118 = (-t189 + t190) * pkin(3) + t228;
t117 = t283 + (t191 + t192) * qJ(4) + t229;
t116 = pkin(2) * t189 + t252 - t274;
t115 = pkin(2) * t192 + t275 - t292;
t114 = t125 * t223 - t171 * t225;
t113 = -t125 * t225 - t171 * t223;
t112 = t127 * t224 + t222 * t231;
t111 = t127 * t222 - t224 * t231;
t110 = pkin(3) * t199 + qJ(4) * t203 - t231 + t286;
t109 = qJ(4) * t198 + (t204 - t226) * pkin(3) + t232 + t287;
t108 = pkin(2) * t125 - qJ(2) * t126;
t107 = -t123 * t224 - t124 * t222 - t284;
t106 = -pkin(2) * t171 + t285 * t126;
t105 = -t117 * t222 + (-pkin(2) - t281) * t192 + t292;
t104 = -t118 * t224 + (pkin(2) + t277) * t189 + t252;
t103 = t111 * t223 - t122 * t225;
t102 = -t111 * t225 - t122 * t223;
t101 = pkin(2) * t111 - pkin(3) * t231 - qJ(2) * t112 + qJ(4) * t127;
t100 = t285 * t112 + (-pkin(2) + t238) * t122;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t195, -t194, 0, t154, 0, 0, 0, 0, 0, 0, 0, t195, t194, t134, 0, 0, 0, 0, 0, 0, t135, -t137, t150, t114, 0, 0, 0, 0, 0, 0, t135, t150, t137, t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t194, -t195, 0, t153, 0, 0, 0, 0, 0, 0, 0, -t194, t195, t130, 0, 0, 0, 0, 0, 0, t131, t132, t149, t113, 0, 0, 0, 0, 0, 0, t131, t149, -t132, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t163, -t165, 0, t126, 0, 0, 0, 0, 0, 0, t163, 0, t165, t112; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t194, 0, -t195, 0, -t241, -t240, -t153, -pkin(4) * t153, 0, -t194, t195, 0, 0, 0, -t130, t241, t240, -pkin(4) * t130 + (-pkin(1) * t223 + qJ(2) * t225) * g(3), t142, -t291, t139, t143, t289, t175, -t116 * t223 + t119 * t225 - t280, -t115 * t223 + t120 * t225 - t295, -pkin(2) * t269 - t121 * t223 - t279, -pkin(4) * t113 - t106 * t223 + t108 * t225, t142, t139, t291, t175, -t289, t143, -t104 * t223 + t110 * t225 - t280, -t107 * t223 + t155 * t225 - t279, -t105 * t223 + t109 * t225 + t295, -pkin(4) * t102 - t100 * t223 + t101 * t225; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t195, 0, t194, 0, t240, -t241, t154, pkin(4) * t154, 0, -t195, -t194, 0, 0, 0, t134, -t240, t241, pkin(4) * t134 + (pkin(1) * t225 + qJ(2) * t223) * g(3), t140, t290, t138, t141, -t288, t174, t116 * t225 + t119 * t223 + t129, t115 * t225 + t120 * t223 - t294, -pkin(2) * t270 + t121 * t225 + t146, pkin(4) * t114 + t106 * t225 + t108 * t223, t140, t138, -t290, t174, t288, t141, t104 * t225 + t110 * t223 + t129, t107 * t225 + t155 * t223 + t146, t105 * t225 + t109 * t223 + t294, pkin(4) * t103 + t100 * t225 + t101 * t223; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t200, t201, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t239 - 0.2e1 * t276, -t201 + t213 + 0.2e1 * t217, pkin(1) * t178 + qJ(2) * t177, t164, -t236, t166, t156, -t235, 0, t243 - t275, qJ(2) * t192 - t274 - t293, -t125 + t242, -qJ(2) * t171 + t285 * t125, t164, t166, t236, 0, t235, t156, -qJ(4) * t273 - t118 * t222 + t243, -t123 * t222 + t124 * t224 + t242, t117 * t224 + (-qJ(2) - t282) * t192 + t293, t285 * t111 + (-qJ(2) - t237) * t122;];
tauB_reg = t1;
