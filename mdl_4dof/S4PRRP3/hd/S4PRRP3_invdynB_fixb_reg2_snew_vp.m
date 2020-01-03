% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRP3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:57
% EndTime: 2019-12-31 16:27:00
% DurationCPUTime: 1.96s
% Computational Cost: add. (3325->259), mult. (7253->352), div. (0->0), fcn. (4336->6), ass. (0->200)
t274 = sin(pkin(6));
t278 = sin(qJ(2));
t280 = cos(qJ(2));
t275 = cos(pkin(6));
t251 = t274 * g(1) - t275 * g(2);
t252 = t275 * g(1) + t274 * g(2);
t288 = t280 * t251 + t278 * t252;
t310 = -t278 * t251 + t280 * t252;
t290 = -t278 * t288 - t280 * t310;
t162 = t278 * t310 - t280 * t288;
t319 = t275 * t162;
t344 = -t274 * t290 + t319;
t323 = t274 * t162;
t121 = t275 * t290 + t323;
t282 = qJD(2) ^ 2;
t305 = t278 * qJDD(2);
t247 = t280 * t282 + t305;
t303 = t280 * qJDD(2);
t248 = -t278 * t282 + t303;
t197 = -t274 * t247 + t275 * t248;
t272 = g(3) - qJDD(1);
t226 = pkin(4) * t247 - t280 * t272;
t285 = -pkin(4) * t248 - t278 * t272;
t343 = -qJ(1) * t197 + t274 * t226 + t275 * t285;
t336 = t275 * t247 + t274 * t248;
t341 = qJ(1) * t336 + t275 * t226 - t274 * t285;
t277 = sin(qJ(3));
t279 = cos(qJ(3));
t260 = t277 * t282 * t279;
t302 = qJDD(3) + t260;
t340 = t302 * pkin(3);
t190 = -t282 * pkin(2) + qJDD(2) * pkin(5) - t310;
t173 = t277 * t190 + t279 * t272;
t307 = qJD(2) * qJD(3);
t294 = t279 * t307;
t306 = t277 * qJDD(2);
t242 = t294 + t306;
t229 = t242 * qJ(4);
t335 = -t229 - t173 + t340;
t281 = qJD(3) ^ 2;
t270 = t277 ^ 2;
t324 = t270 * t282;
t257 = -t281 - t324;
t255 = qJDD(3) - t260;
t314 = t277 * t255;
t214 = t279 * t257 - t314;
t333 = pkin(2) * t214;
t271 = t279 ^ 2;
t268 = t271 * t282;
t259 = -t268 - t281;
t315 = t277 * t302;
t216 = t279 * t259 - t315;
t296 = t277 * t307;
t304 = t279 * qJDD(2);
t244 = -0.2e1 * t296 + t304;
t174 = t278 * t216 + t280 * t244;
t332 = pkin(4) * t174;
t311 = t279 * t255;
t218 = -t277 * t257 - t311;
t241 = 0.2e1 * t294 + t306;
t175 = t278 * t218 - t280 * t241;
t331 = pkin(4) * t175;
t309 = t270 + t271;
t246 = t309 * qJDD(2);
t249 = t268 + t324;
t200 = t278 * t246 + t280 * t249;
t330 = pkin(4) * t200;
t239 = t279 * t302;
t212 = t277 * t259 + t239;
t329 = pkin(5) * t212;
t328 = pkin(5) * t214;
t177 = t280 * t216 - t278 * t244;
t131 = t275 * t174 + t274 * t177;
t327 = qJ(1) * t131;
t178 = t280 * t218 + t278 * t241;
t132 = t275 * t175 + t274 * t178;
t326 = qJ(1) * t132;
t201 = t280 * t246 - t278 * t249;
t156 = t275 * t200 + t274 * t201;
t325 = qJ(1) * t156;
t320 = t274 * t272;
t318 = t275 * t272;
t140 = (qJ(4) * qJD(3) * t279 - 0.2e1 * qJD(4) * t277) * qJD(2) + t335;
t317 = t277 * t140;
t189 = -qJDD(2) * pkin(2) - t282 * pkin(5) - t288;
t316 = t277 * t189;
t313 = t279 * t140;
t312 = t279 * t189;
t176 = t279 * t190 - t277 * t272;
t308 = qJD(2) * t277;
t301 = 0.2e1 * qJD(2) * qJD(4);
t300 = pkin(1) * t174 + pkin(2) * t244 + pkin(5) * t216;
t299 = pkin(1) * t175 - pkin(2) * t241 + pkin(5) * t218;
t297 = pkin(1) * t200 + pkin(2) * t249 + pkin(5) * t246;
t295 = t277 * t305;
t293 = t277 * t303;
t292 = -pkin(1) * t212 + pkin(4) * t177;
t291 = -pkin(1) * t214 + pkin(4) * t178;
t137 = t277 * t173 + t279 * t176;
t203 = -t274 * t251 - t275 * t252;
t287 = t278 * t260;
t286 = t280 * t260;
t151 = -pkin(2) * t212 + t173;
t136 = t279 * t173 - t277 * t176;
t202 = t275 * t251 - t274 * t252;
t243 = -t296 + t304;
t253 = qJD(3) * pkin(3) - qJ(4) * t308;
t283 = t243 * qJ(4) - qJD(3) * t253 + t279 * t301 + t176;
t153 = -t243 * pkin(3) - qJ(4) * t268 + t253 * t308 + qJDD(4) + t189;
t262 = t277 * t301;
t258 = t268 - t281;
t256 = t281 - t324;
t250 = t268 - t324;
t238 = t309 * t307;
t222 = t278 * qJDD(3) + t280 * t238;
t221 = t279 * t242 - t270 * t307;
t220 = -t280 * qJDD(3) + t278 * t238;
t219 = -t277 * t243 - t271 * t307;
t217 = -t277 * t256 + t239;
t215 = t279 * t258 - t314;
t213 = t279 * t256 + t315;
t211 = t277 * t258 + t311;
t210 = (t242 + t294) * t277;
t209 = (t243 - t296) * t279;
t208 = -pkin(3) * t241 - qJ(4) * t255;
t194 = pkin(4) * t201;
t192 = -t277 * t241 + t279 * t244;
t191 = t279 * t241 + t277 * t244;
t188 = t280 * t221 - t287;
t187 = t280 * t219 + t287;
t186 = t278 * t221 + t286;
t185 = t278 * t219 - t286;
t184 = t280 * t217 + t295;
t183 = t280 * t215 + t278 * t304;
t182 = t278 * t217 - t293;
t181 = t278 * t215 - t279 * t303;
t167 = t280 * t192 - t278 * t250;
t166 = t278 * t192 + t280 * t250;
t165 = -t274 * t220 + t275 * t222;
t164 = t275 * t220 + t274 * t222;
t159 = t312 - t328;
t158 = t316 - t329;
t157 = -t274 * t200 + t275 * t201;
t155 = qJ(1) * t157;
t154 = pkin(1) * t272 + pkin(4) * t290;
t152 = t176 - t333;
t150 = -qJ(4) * t257 + t153;
t149 = -t274 * t186 + t275 * t188;
t148 = -t274 * t185 + t275 * t187;
t147 = t275 * t186 + t274 * t188;
t146 = t275 * t185 + t274 * t187;
t145 = -t274 * t182 + t275 * t184;
t144 = -t274 * t181 + t275 * t183;
t143 = t275 * t182 + t274 * t184;
t142 = t275 * t181 + t274 * t183;
t141 = -pkin(3) * t268 + t283;
t139 = t262 + (-t294 + t306) * qJ(4) - t335;
t138 = pkin(3) * t244 + qJ(4) * t259 - t153;
t134 = -t274 * t175 + t275 * t178;
t133 = -t274 * t174 + t275 * t177;
t130 = qJ(1) * t134;
t129 = qJ(1) * t133;
t128 = qJ(4) * t304 + (t249 - t268) * pkin(3) + t283;
t127 = -t274 * t166 + t275 * t167;
t126 = t275 * t166 + t274 * t167;
t125 = -t333 + (-t257 - t268) * pkin(3) + t283;
t124 = -qJ(4) * t294 + t151 + t229 + t262 - 0.2e1 * t340;
t123 = t280 * t136 - t330;
t122 = t278 * t136 + t194;
t119 = t280 * t137 + t278 * t189;
t118 = t278 * t137 - t280 * t189;
t117 = -qJ(4) * t239 - t277 * t138 - t329;
t116 = t279 * t150 - t277 * t208 - t328;
t115 = -pkin(3) * t153 + qJ(4) * t141;
t114 = t279 * t141 - t317;
t113 = t277 * t141 + t313;
t112 = -t278 * t152 + t280 * t159 - t331;
t111 = -t278 * t151 + t280 * t158 - t332;
t110 = -t277 * t128 + t279 * t139;
t109 = t280 * t152 + t278 * t159 + t291;
t108 = t280 * t151 + t278 * t158 + t292;
t107 = -pkin(3) * t295 + t280 * t110 - t330;
t106 = pkin(3) * t293 + t278 * t110 + t194;
t105 = t280 * t114 + t278 * t153;
t104 = t278 * t114 - t280 * t153;
t103 = -pkin(2) * t113 - pkin(3) * t140;
t102 = -t274 * t118 + t275 * t119;
t101 = t275 * t118 + t274 * t119;
t100 = t280 * t116 - t278 * t125 - t331;
t99 = t280 * t117 - t278 * t124 - t332;
t98 = t278 * t116 + t280 * t125 + t291;
t97 = t278 * t117 + t280 * t124 + t292;
t96 = -pkin(4) * t118 - (pkin(2) * t278 - pkin(5) * t280) * t136;
t95 = -pkin(5) * t113 - qJ(4) * t313 - t277 * t115;
t94 = pkin(4) * t119 - (-pkin(2) * t280 - pkin(5) * t278 - pkin(1)) * t136;
t93 = -t274 * t104 + t275 * t105;
t92 = t275 * t104 + t274 * t105;
t91 = -pkin(4) * t104 - t278 * t103 + t280 * t95;
t90 = -pkin(1) * t113 + pkin(4) * t105 + t280 * t103 + t278 * t95;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, 0, 0, 0, 0, 0, 0, -t336, -t197, 0, t121, 0, 0, 0, 0, 0, 0, t133, t134, t157, t102, 0, 0, 0, 0, 0, 0, t133, t134, t157, t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t202, 0, 0, 0, 0, 0, 0, t197, -t336, 0, -t344, 0, 0, 0, 0, 0, 0, t131, t132, t156, t101, 0, 0, 0, 0, 0, 0, t131, t132, t156, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t272, 0, 0, 0, 0, 0, 0, t212, t214, 0, -t136, 0, 0, 0, 0, 0, 0, t212, t214, 0, t113; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t320, -t318, -t202, -qJ(1) * t202, 0, 0, t197, 0, -t336, 0, t343, t341, t344, pkin(4) * t319 + qJ(1) * t344 - t274 * t154, t149, t127, t145, t148, t144, t165, -t274 * t108 + t275 * t111 - t327, -t274 * t109 + t275 * t112 - t326, -t274 * t122 + t275 * t123 - t325, -qJ(1) * t101 - t274 * t94 + t275 * t96, t149, t127, t145, t148, t144, t165, -t274 * t97 + t275 * t99 - t327, t275 * t100 - t274 * t98 - t326, -t274 * t106 + t275 * t107 - t325, -qJ(1) * t92 - t274 * t90 + t275 * t91; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t318, -t320, t203, qJ(1) * t203, 0, 0, t336, 0, t197, 0, -t341, t343, t121, pkin(4) * t323 + qJ(1) * t121 + t275 * t154, t147, t126, t143, t146, t142, t164, t275 * t108 + t274 * t111 + t129, t275 * t109 + t274 * t112 + t130, t275 * t122 + t274 * t123 + t155, qJ(1) * t102 + t274 * t96 + t275 * t94, t147, t126, t143, t146, t142, t164, t274 * t99 + t275 * t97 + t129, t274 * t100 + t275 * t98 + t130, t275 * t106 + t274 * t107 + t155, qJ(1) * t93 + t274 * t91 + t275 * t90; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t251, t252, 0, 0, 0, 0, 0, 0, 0, qJDD(2), pkin(1) * t248 + t288, -pkin(1) * t247 + t310, 0, -pkin(1) * t162, t210, t191, t213, t209, t211, 0, t300 - t312, t299 + t316, t137 + t297, pkin(1) * t118 - pkin(2) * t189 + pkin(5) * t137, t210, t191, t213, t209, t211, 0, -qJ(4) * t315 + t279 * t138 + t300, t277 * t150 + t279 * t208 + t299, t279 * t128 + t277 * t139 + t297, pkin(1) * t104 - pkin(2) * t153 + pkin(5) * t114 - qJ(4) * t317 + t279 * t115;];
tauB_reg = t1;
