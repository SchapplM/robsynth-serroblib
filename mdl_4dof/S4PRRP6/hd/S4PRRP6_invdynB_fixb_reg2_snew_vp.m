% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4PRRP6
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
% Datum: 2019-12-31 16:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4PRRP6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP6_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP6_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP6_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP6_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:30:44
% EndTime: 2019-12-31 16:30:48
% DurationCPUTime: 2.44s
% Computational Cost: add. (2929->265), mult. (6131->352), div. (0->0), fcn. (3595->6), ass. (0->201)
t277 = qJD(3) ^ 2;
t273 = sin(qJ(3));
t267 = t273 ^ 2;
t278 = qJD(2) ^ 2;
t330 = t267 * t278;
t252 = t277 + t330;
t275 = cos(qJ(3));
t257 = t273 * t278 * t275;
t251 = qJDD(3) - t257;
t315 = t275 * t251;
t218 = -t273 * t252 + t315;
t306 = qJD(2) * qJD(3);
t298 = t275 * t306;
t304 = t273 * qJDD(2);
t241 = 0.2e1 * t298 + t304;
t274 = sin(qJ(2));
t276 = cos(qJ(2));
t173 = t276 * t218 - t274 * t241;
t321 = t273 * t251;
t212 = t275 * t252 + t321;
t270 = sin(pkin(6));
t271 = cos(pkin(6));
t141 = t270 * t173 - t271 * t212;
t363 = qJ(1) * t141;
t144 = t271 * t173 + t270 * t212;
t362 = qJ(1) * t144;
t169 = t274 * t218 + t276 * t241;
t361 = pkin(4) * t169;
t360 = -pkin(1) * t169 - pkin(5) * t218;
t359 = -pkin(1) * t212 + pkin(4) * t173;
t323 = t273 * t241;
t300 = t273 * t306;
t302 = t275 * qJDD(2);
t283 = 0.2e1 * t300 - t302;
t343 = t283 * t275;
t191 = t343 + t323;
t268 = t275 ^ 2;
t247 = (t267 - t268) * t278;
t162 = t276 * t191 - t274 * t247;
t189 = t275 * t241 - t273 * t283;
t358 = t270 * t162 + t271 * t189;
t357 = t271 * t162 - t270 * t189;
t329 = t268 * t278;
t254 = -t277 + t329;
t216 = -t275 * t254 + t321;
t177 = t276 * t216 - t274 * t302;
t210 = t273 * t254 + t315;
t356 = t270 * t177 + t271 * t210;
t355 = t271 * t177 - t270 * t210;
t354 = 2 * qJD(4);
t352 = pkin(2) * t212;
t351 = pkin(5) * t212;
t342 = t274 * t191 + t276 * t247;
t301 = t276 * qJDD(2);
t341 = t274 * t216 + t275 * t301;
t312 = g(3) - qJDD(1);
t340 = t270 * t312;
t339 = t271 * t312;
t248 = t270 * g(1) - t271 * g(2);
t232 = t271 * t248;
t249 = t271 * g(1) + t270 * g(2);
t197 = -t270 * t249 + t232;
t331 = qJ(4) * t273;
t337 = pkin(3) * t275;
t289 = -t331 - t337;
t240 = t289 * qJD(2);
t226 = -t276 * t249 - t274 * t312;
t200 = -t278 * pkin(2) + qJDD(2) * pkin(5) + t226;
t310 = -t275 * t200 + t273 * t248;
t287 = t275 * qJD(2) * t240 + qJDD(3) * qJ(4) + (qJD(3) * t354) - t310;
t307 = qJD(2) * t273;
t338 = t240 * t307 + qJDD(4);
t255 = -t277 - t329;
t250 = qJDD(3) + t257;
t322 = t273 * t250;
t215 = t275 * t255 - t322;
t168 = t274 * t215 - t276 * t283;
t336 = pkin(4) * t168;
t308 = t267 + t268;
t243 = t308 * qJDD(2);
t246 = t308 * t278;
t195 = t274 * t243 + t276 * t246;
t335 = pkin(4) * t195;
t316 = t275 * t250;
t209 = t273 * t255 + t316;
t334 = pkin(5) * t209;
t172 = t276 * t215 + t274 * t283;
t140 = t270 * t172 - t271 * t209;
t333 = qJ(1) * t140;
t196 = t276 * t243 - t274 * t246;
t332 = qJ(1) * t196;
t245 = -t274 * t278 + t301;
t328 = t270 * t245;
t327 = t270 * t248;
t187 = t271 * t196;
t325 = t271 * t245;
t294 = t274 * t249 - t276 * t312;
t286 = qJDD(2) * pkin(2) + t278 * pkin(5) + t294;
t324 = t273 * t286;
t318 = t275 * t286;
t311 = -pkin(1) * t209 + pkin(4) * t172;
t167 = t273 * t200 + t275 * t248;
t309 = t246 - t277;
t303 = t274 * qJDD(2);
t244 = t276 * t278 + t303;
t296 = pkin(1) * t245 + qJ(1) * t244 + t294;
t295 = -pkin(1) * t244 + qJ(1) * t245 - t226;
t165 = t276 * t226 - t274 * t294;
t198 = -t271 * t249 - t327;
t293 = t274 * t257;
t292 = t276 * t257;
t148 = -pkin(2) * t209 + t167;
t290 = -pkin(1) * t168 - pkin(5) * t215;
t288 = pkin(3) * t273 - qJ(4) * t275;
t203 = pkin(4) * t244 - t276 * t248;
t202 = -pkin(4) * t245 - t274 * t248;
t133 = t275 * t167 + t273 * t310;
t134 = t273 * t167 - t275 * t310;
t164 = -t274 * t226 - t276 * t294;
t285 = t298 + t304;
t284 = -t300 + t302;
t282 = -pkin(1) * t195 - pkin(2) * t246 - pkin(5) * t243;
t281 = -qJDD(3) * pkin(3) + t167 + t338;
t280 = -t284 * pkin(3) - t286 + (-t285 - t298) * qJ(4);
t279 = t307 * t354 - t280;
t253 = t277 - t330;
t239 = t288 * qJDD(2);
t238 = t308 * t306;
t230 = t271 * t244;
t229 = t270 * t244;
t224 = t274 * qJDD(3) + t276 * t238;
t223 = -t267 * t306 + t275 * t285;
t222 = -t276 * qJDD(3) + t274 * t238;
t221 = -t268 * t306 - t273 * t284;
t217 = -t273 * t253 + t316;
t211 = -t275 * t253 - t322;
t205 = t271 * t224;
t204 = t270 * t224;
t194 = pkin(4) * t196;
t186 = t270 * t196;
t183 = qJ(1) * t187;
t182 = t276 * t223 - t293;
t181 = t276 * t221 + t293;
t180 = t274 * t223 + t292;
t179 = t274 * t221 - t292;
t178 = t276 * t217 + t273 * t303;
t175 = t274 * t217 - t273 * t301;
t160 = -t318 + t351;
t159 = -t324 - t334;
t158 = t271 * t165 - t327;
t157 = t270 * t165 + t232;
t156 = t271 * t182 + t270 * t323;
t155 = t271 * t181 - t270 * t343;
t154 = t270 * t182 - t271 * t323;
t153 = t270 * t181 + t271 * t343;
t152 = t271 * t178 - t270 * t211;
t151 = t270 * t178 + t271 * t211;
t150 = t277 * qJ(4) - t281;
t149 = -t310 + t352;
t147 = -t277 * pkin(3) + t287;
t146 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t307 + t280;
t143 = t271 * t172 + t270 * t209;
t139 = qJ(1) * t143;
t138 = t309 * qJ(4) + t281;
t137 = t309 * pkin(3) + t287;
t136 = (-t283 - t300) * pkin(3) + t279;
t135 = -pkin(3) * t300 + qJ(4) * t241 + t279;
t131 = pkin(2) * t283 + t290 - t318;
t130 = pkin(2) * t241 + t324 - t360;
t129 = (-t255 - t277) * qJ(4) + (-qJDD(3) - t250) * pkin(3) + t148 + t338;
t128 = -t352 - qJ(4) * t251 + (-t252 + t277) * pkin(3) - t287;
t127 = t276 * t133 - t335;
t126 = t276 * t134 - t274 * t286;
t125 = t274 * t134 + t276 * t286;
t124 = -pkin(3) * t323 + t275 * t135 - t351;
t123 = -qJ(4) * t343 - t273 * t136 - t334;
t122 = t275 * t147 - t273 * t150;
t121 = t273 * t147 + t275 * t150;
t120 = -t134 + t282;
t119 = -t273 * t137 + t275 * t138;
t118 = -t274 * t149 + t276 * t160 + t361;
t117 = -t274 * t148 + t276 * t159 - t336;
t116 = -t273 * t135 + (-pkin(2) - t337) * t241 + t360;
t115 = -t275 * t136 - (-pkin(2) - t331) * t283 + t290;
t114 = t276 * t119 - t274 * t239 - t335;
t113 = -t275 * t137 - t273 * t138 + t282;
t112 = t276 * t122 + t274 * t146;
t111 = t274 * t122 - t276 * t146;
t110 = t271 * t126 - t133 * t270;
t109 = t270 * t126 + t133 * t271;
t108 = -pkin(1) * t125 - pkin(2) * t286 - pkin(5) * t134;
t107 = -pkin(2) * t121 - pkin(3) * t150 - qJ(4) * t147;
t106 = t276 * t123 - t274 * t129 - t336;
t105 = t276 * t124 - t274 * t128 - t361;
t104 = -pkin(5) * t121 + t146 * t288;
t103 = -pkin(4) * t125 - (pkin(2) * t274 - pkin(5) * t276) * t133;
t102 = t271 * t112 + t270 * t121;
t101 = t270 * t112 - t271 * t121;
t100 = -pkin(1) * t111 - pkin(5) * t122 + (pkin(2) - t289) * t146;
t99 = -pkin(4) * t111 + t276 * t104 - t274 * t107;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, 0, 0, 0, 0, 0, 0, -t230, -t325, 0, t158, 0, 0, 0, 0, 0, 0, t143, -t144, t187, t110, 0, 0, 0, 0, 0, 0, t143, t187, t144, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, 0, 0, 0, 0, 0, 0, -t229, -t328, 0, t157, 0, 0, 0, 0, 0, 0, t140, -t141, t186, t109, 0, 0, 0, 0, 0, 0, t140, t186, t141, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t312, 0, 0, 0, 0, 0, 0, t245, -t244, 0, -t164, 0, 0, 0, 0, 0, 0, t168, -t169, t195, t125, 0, 0, 0, 0, 0, 0, t168, t195, t169, t111; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t340, -t339, -t197, -qJ(1) * t197, 0, 0, t325, 0, -t230, t270 * qJDD(2), t271 * t202 + t296 * t270, t271 * t203 + t295 * t270, t271 * t164, -qJ(1) * t157 - (pkin(1) * t270 - pkin(4) * t271) * t164, t156, -t357, t152, t155, -t355, t205, t271 * t117 - t270 * t131 - t333, t271 * t118 - t270 * t130 + t363, t271 * t127 + (-t120 - t332) * t270, -qJ(1) * t109 + t271 * t103 - t270 * t108, t156, t152, t357, t205, t355, t155, t271 * t106 - t270 * t115 - t333, t271 * t114 + (-t113 - t332) * t270, t271 * t105 - t270 * t116 - t363, -qJ(1) * t101 - t270 * t100 + t271 * t99; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t339, -t340, t198, qJ(1) * t198, 0, 0, t328, 0, -t229, -t271 * qJDD(2), t270 * t202 - t296 * t271, t270 * t203 - t271 * t295, t270 * t164, qJ(1) * t158 - (-pkin(1) * t271 - pkin(4) * t270) * t164, t154, -t358, t151, t153, -t356, t204, t270 * t117 + t271 * t131 + t139, t270 * t118 + t271 * t130 - t362, t271 * t120 + t270 * t127 + t183, qJ(1) * t110 + t270 * t103 + t271 * t108, t154, t151, t358, t204, t356, t153, t270 * t106 + t271 * t115 + t139, t271 * t113 + t270 * t114 + t183, t270 * t105 + t271 * t116 + t362, qJ(1) * t102 + t271 * t100 + t270 * t99; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t248, t249, 0, 0, 0, 0, t244, 0, t245, 0, -t203, t202, t165, pkin(1) * t248 + pkin(4) * t165, t180, -t342, t175, t179, -t341, t222, t276 * t148 + t274 * t159 + t311, t276 * t149 + t274 * t160 - t359, t274 * t133 + t194, pkin(4) * t126 - (-pkin(2) * t276 - pkin(5) * t274 - pkin(1)) * t133, t180, t175, t342, t222, t341, t179, t274 * t123 + t276 * t129 + t311, t274 * t119 + t276 * t239 + t194, t274 * t124 + t276 * t128 + t359, -pkin(1) * t121 + pkin(4) * t112 + t274 * t104 + t276 * t107;];
tauB_reg = t1;
