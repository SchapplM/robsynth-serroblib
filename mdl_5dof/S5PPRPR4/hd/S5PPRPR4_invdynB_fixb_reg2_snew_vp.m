% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:26
% EndTime: 2019-12-31 17:32:31
% DurationCPUTime: 3.49s
% Computational Cost: add. (8748->349), mult. (18537->510), div. (0->0), fcn. (13033->8), ass. (0->226)
t332 = cos(qJ(3));
t330 = sin(qJ(3));
t350 = t330 * qJDD(3);
t372 = qJD(3) ^ 2;
t301 = t332 * t372 + t350;
t349 = t332 * qJDD(3);
t354 = t330 * t372;
t302 = -t349 + t354;
t326 = sin(pkin(7));
t328 = cos(pkin(7));
t256 = t328 * t301 + t326 * t302;
t323 = g(3) - qJDD(1);
t275 = pkin(5) * t302 + t330 * t323;
t341 = pkin(5) * t301 + t332 * t323;
t390 = -qJ(1) * t256 + t275 * t326 + t328 * t341;
t329 = sin(qJ(5));
t327 = cos(pkin(8));
t331 = cos(qJ(5));
t325 = sin(pkin(8));
t362 = t325 * t329;
t284 = (-t327 * t331 + t362) * qJD(3);
t338 = t325 * t331 + t327 * t329;
t286 = t338 * qJD(3);
t248 = t286 * t284;
t377 = qJDD(5) - t248;
t389 = t329 * t377;
t388 = t331 * t377;
t304 = t328 * g(1) + t326 * g(2);
t303 = t326 * g(1) - t328 * g(2);
t342 = qJDD(2) - t303;
t260 = -t332 * t304 + t330 * t342;
t242 = -t372 * pkin(3) + qJDD(3) * qJ(4) + t260;
t353 = qJD(3) * qJD(4);
t387 = t242 + 0.2e1 * t353;
t320 = t325 ^ 2;
t321 = t327 ^ 2;
t376 = t372 * (t320 + t321);
t292 = t327 * t376;
t267 = t332 * t292 + t327 * t350;
t345 = t327 * t349;
t337 = -t330 * t292 + t345;
t386 = t267 * t326 + t328 * t337;
t385 = t267 * t328 - t326 * t337;
t343 = -t326 * t301 + t328 * t302;
t383 = -qJ(1) * t343 + t275 * t328 - t326 * t341;
t259 = -t330 * t304 - t332 * t342;
t213 = t259 * t332 - t260 * t330;
t214 = t259 * t330 + t260 * t332;
t177 = t213 * t328 + t214 * t326;
t382 = t213 * t326 - t214 * t328;
t317 = t321 * t372;
t355 = t320 * t372;
t298 = t317 + t355;
t378 = t338 * qJDD(3);
t293 = t326 * t304;
t375 = t328 * t342 + t293;
t294 = t328 * t304;
t254 = t326 * t342 - t294;
t280 = t284 ^ 2;
t281 = t286 ^ 2;
t373 = -0.2e1 * t325;
t371 = pkin(1) + pkin(2);
t324 = qJDD(3) * pkin(3);
t311 = t327 * t323;
t207 = t353 * t373 + t311 + (pkin(4) * t372 * t327 - pkin(6) * qJDD(3) - t242) * t325;
t232 = t325 * t323 + t387 * t327;
t352 = qJDD(3) * t327;
t210 = -pkin(4) * t317 + pkin(6) * t352 + t232;
t174 = -t331 * t207 + t210 * t329;
t175 = t329 * t207 + t331 * t210;
t143 = -t174 * t331 + t175 * t329;
t370 = t143 * t325;
t369 = t143 * t327;
t237 = -t372 * qJ(4) + qJDD(4) + t259 - t324;
t229 = -pkin(4) * t352 - t298 * pkin(6) + t237;
t368 = t229 * t329;
t367 = t229 * t331;
t240 = qJDD(5) + t248;
t366 = t240 * t329;
t365 = t240 * t331;
t364 = t325 * t237;
t363 = t325 * t327;
t359 = t326 * t323;
t358 = t327 * t237;
t312 = t328 * t323;
t357 = t284 * qJD(5);
t356 = t286 * qJD(5);
t351 = t320 * qJDD(3);
t347 = t330 * t248;
t346 = t332 * t248;
t344 = -t237 + t324;
t144 = t174 * t329 + t331 * t175;
t261 = t328 * t303 - t293;
t262 = -t326 * t303 - t294;
t282 = qJDD(3) * t362 - t331 * t352;
t231 = t387 * t325 - t311;
t191 = t231 * t327 - t232 * t325;
t192 = t325 * t231 + t327 * t232;
t269 = t301 * t363;
t270 = t325 * t345 - t354 * t363;
t340 = t269 * t328 - t270 * t326;
t339 = t269 * t326 + t270 * t328;
t333 = qJD(5) ^ 2;
t316 = t321 * qJDD(3);
t299 = t317 - t355;
t297 = t316 - t351;
t296 = t316 + t351;
t291 = t325 * t376;
t273 = -t281 - t333;
t272 = -t281 + t333;
t271 = t280 - t333;
t266 = t332 * t291 + t325 * t350;
t263 = t330 * t291 - t325 * t349;
t252 = t332 * t297 - t330 * t299;
t251 = t332 * t296 - t330 * t298;
t250 = -t330 * t297 - t332 * t299;
t249 = t330 * t296 + t332 * t298;
t247 = -t281 + t280;
t246 = t378 - t357;
t245 = t378 - 0.2e1 * t357;
t244 = -t282 - t356;
t243 = t282 + 0.2e1 * t356;
t238 = -t333 - t280;
t235 = (-t284 * t331 + t286 * t329) * qJD(5);
t234 = (-t284 * t329 - t286 * t331) * qJD(5);
t233 = -t280 - t281;
t228 = t331 * t246 - t329 * t356;
t226 = t326 * t263 + t328 * t266;
t225 = t329 * t246 + t331 * t356;
t223 = -t328 * t263 + t326 * t266;
t222 = -t329 * t244 + t331 * t357;
t221 = t331 * t244 + t329 * t357;
t220 = -t273 * t329 - t365;
t219 = -t272 * t329 + t388;
t218 = t271 * t331 - t366;
t217 = t273 * t331 - t366;
t216 = t272 * t331 + t389;
t215 = t271 * t329 + t365;
t209 = t249 * t326 + t251 * t328;
t208 = -t249 * t328 + t251 * t326;
t203 = pkin(5) * t213 + qJ(2) * t323;
t202 = -pkin(5) * t214 + t371 * t323;
t201 = -t243 * t331 - t245 * t329;
t200 = -t282 * t331 + t329 * t378;
t199 = -t243 * t329 + t245 * t331;
t198 = -t282 * t329 - t331 * t378;
t197 = t238 * t331 - t389;
t196 = t238 * t329 + t388;
t195 = -t234 * t325 + t235 * t327;
t194 = qJDD(5) * t330 + t195 * t332;
t193 = qJDD(5) * t332 - t195 * t330;
t189 = -pkin(6) * t217 + t367;
t188 = -pkin(5) * t263 - t232 * t330 + t332 * t358;
t187 = -pkin(5) * t337 - t231 * t330 + t332 * t364;
t186 = -pkin(5) * t266 - t232 * t332 - t330 * t358;
t185 = pkin(5) * t267 - t231 * t332 - t330 * t364;
t184 = -t225 * t325 + t228 * t327;
t183 = -t221 * t325 + t222 * t327;
t182 = -t217 * t325 + t220 * t327;
t181 = -t216 * t325 + t219 * t327;
t180 = -t215 * t325 + t218 * t327;
t179 = t217 * t327 + t220 * t325;
t176 = -pkin(6) * t196 + t368;
t172 = -pkin(5) * t249 + t191 * t332;
t171 = -pkin(5) * t251 - t191 * t330;
t170 = t181 * t332 + t330 * t378;
t169 = t180 * t332 - t282 * t330;
t168 = -t181 * t330 + t332 * t378;
t167 = -t180 * t330 - t282 * t332;
t166 = t192 * t332 + t237 * t330;
t165 = t192 * t330 - t237 * t332;
t164 = -pkin(4) * t245 + pkin(6) * t220 + t368;
t163 = -t199 * t325 + t201 * t327;
t162 = -t198 * t325 + t200 * t327;
t161 = t198 * t327 + t200 * t325;
t160 = -t196 * t325 + t197 * t327;
t159 = t196 * t327 + t197 * t325;
t158 = t184 * t332 + t347;
t157 = t183 * t332 - t347;
t156 = -t184 * t330 + t346;
t155 = -t183 * t330 - t346;
t154 = t182 * t332 + t245 * t330;
t153 = t182 * t330 - t245 * t332;
t152 = -pkin(4) * t243 + pkin(6) * t197 - t367;
t151 = t163 * t332 - t247 * t330;
t150 = -t163 * t330 - t247 * t332;
t149 = t160 * t332 + t243 * t330;
t148 = t160 * t330 - t243 * t332;
t147 = t162 * t332 + t233 * t330;
t146 = t162 * t330 - t233 * t332;
t145 = -pkin(3) * t161 - pkin(4) * t198;
t142 = t165 * t326 + t166 * t328;
t141 = -t165 * t328 + t166 * t326;
t140 = -pkin(3) * t179 - pkin(4) * t217 + t175;
t139 = -pkin(4) * t229 + pkin(6) * t144;
t138 = t153 * t326 + t154 * t328;
t137 = -t153 * t328 + t154 * t326;
t136 = -pkin(6) * t198 - t143;
t135 = -pkin(3) * t159 - pkin(4) * t196 + t174;
t134 = -pkin(4) * t233 + pkin(6) * t200 + t144;
t133 = t148 * t326 + t149 * t328;
t132 = -t148 * t328 + t149 * t326;
t131 = -qJ(4) * t179 - t164 * t325 + t189 * t327;
t130 = t146 * t326 + t147 * t328;
t129 = -t146 * t328 + t147 * t326;
t128 = -qJ(4) * t159 - t152 * t325 + t176 * t327;
t127 = -pkin(5) * t165 - (pkin(3) * t330 - qJ(4) * t332 + qJ(2)) * t191;
t126 = t144 * t327 - t370;
t125 = t144 * t325 + t369;
t124 = -pkin(5) * t166 - (pkin(3) * t332 + qJ(4) * t330 + t371) * t191;
t123 = t126 * t332 + t229 * t330;
t122 = t126 * t330 - t229 * t332;
t121 = -pkin(3) * t125 - pkin(4) * t143;
t120 = -qJ(4) * t161 - t134 * t325 + t136 * t327;
t119 = -pkin(5) * t153 + qJ(2) * t179 + t131 * t332 - t140 * t330;
t118 = -pkin(5) * t154 - t330 * t131 - t332 * t140 + t371 * t179;
t117 = -pkin(6) * t369 - qJ(4) * t125 - t139 * t325;
t116 = -pkin(5) * t148 + qJ(2) * t159 + t128 * t332 - t135 * t330;
t115 = t122 * t326 + t123 * t328;
t114 = -t122 * t328 + t123 * t326;
t113 = -pkin(5) * t149 - t330 * t128 - t332 * t135 + t371 * t159;
t112 = -pkin(5) * t146 + qJ(2) * t161 + t120 * t332 - t145 * t330;
t111 = -pkin(5) * t147 - t330 * t120 - t332 * t145 + t371 * t161;
t110 = -pkin(5) * t122 + qJ(2) * t125 + t117 * t332 - t121 * t330;
t109 = -pkin(5) * t123 - t330 * t117 - t332 * t121 + t371 * t125;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254, 0, 0, 0, 0, 0, 0, -t256, t343, 0, -t382, 0, 0, 0, 0, 0, 0, -t385, t226, t209, t142, 0, 0, 0, 0, 0, 0, t133, t138, t130, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t375, 0, 0, 0, 0, 0, 0, t343, t256, 0, t177, 0, 0, 0, 0, 0, 0, -t386, t223, t208, t141, 0, 0, 0, 0, 0, 0, t132, t137, t129, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t323, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t323, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t323, 0, 0, 0, 0, 0, 0, 0, 0, 0, t191, 0, 0, 0, 0, 0, 0, -t159, -t179, -t161, -t125; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t359, -t312, -t261, -qJ(1) * t261, 0, 0, 0, 0, 0, 0, -t359, t375, t312, qJ(1) * t375 + (-pkin(1) * t326 + qJ(2) * t328) * t323, 0, 0, -t343, 0, -t256, 0, t383, t390, t177, -qJ(1) * t177 - t202 * t326 + t203 * t328, t339, -t250 * t326 + t252 * t328, t226, -t339, t385, 0, qJ(1) * t386 - t185 * t326 + t187 * t328, -qJ(1) * t223 - t186 * t326 + t188 * t328, -qJ(1) * t208 - t171 * t326 + t172 * t328, -qJ(1) * t141 - t124 * t326 + t127 * t328, -t156 * t326 + t158 * t328, -t150 * t326 + t151 * t328, -t168 * t326 + t170 * t328, -t155 * t326 + t157 * t328, -t167 * t326 + t169 * t328, -t193 * t326 + t194 * t328, -qJ(1) * t132 - t113 * t326 + t116 * t328, -qJ(1) * t137 - t118 * t326 + t119 * t328, -qJ(1) * t129 - t111 * t326 + t112 * t328, -qJ(1) * t114 - t109 * t326 + t110 * t328; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t312, -t359, t262, qJ(1) * t262, 0, 0, 0, 0, 0, 0, t312, t254, t359, qJ(1) * t254 + (pkin(1) * t328 + qJ(2) * t326) * t323, 0, 0, -t256, 0, t343, 0, t390, -t383, t382, -qJ(1) * t382 + t202 * t328 + t203 * t326, -t340, t250 * t328 + t252 * t326, t223, t340, t386, 0, -qJ(1) * t385 + t185 * t328 + t187 * t326, qJ(1) * t226 + t186 * t328 + t188 * t326, qJ(1) * t209 + t171 * t328 + t172 * t326, qJ(1) * t142 + t124 * t328 + t127 * t326, t156 * t328 + t158 * t326, t150 * t328 + t151 * t326, t168 * t328 + t170 * t326, t155 * t328 + t157 * t326, t167 * t328 + t169 * t326, t193 * t328 + t194 * t326, qJ(1) * t133 + t113 * t328 + t116 * t326, qJ(1) * t138 + t118 * t328 + t119 * t326, qJ(1) * t130 + t111 * t328 + t112 * t326, qJ(1) * t115 + t109 * t328 + t110 * t326; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t303, t304, 0, 0, 0, 0, 0, 0, 0, 0, -t342, 0, -t304, -pkin(1) * t342 - qJ(2) * t304, 0, 0, 0, 0, 0, -qJDD(3), -qJ(2) * t301 + t371 * t302 + t259, qJ(2) * t302 + t371 * t301 + t260, 0, qJ(2) * t214 + t213 * t371, -t351, t352 * t373, 0, -t316, 0, 0, -qJ(2) * t267 + qJ(4) * t292 - t344 * t327 - t337 * t371, qJ(2) * t266 - qJ(4) * t291 - t371 * t263 + t344 * t325, -pkin(3) * t298 + qJ(2) * t251 - qJ(4) * t296 - t371 * t249 - t192, pkin(3) * t237 + qJ(2) * t166 - qJ(4) * t192 - t371 * t165, -t225 * t327 - t228 * t325, -t199 * t327 - t201 * t325, -t216 * t327 - t219 * t325, -t221 * t327 - t222 * t325, -t215 * t327 - t218 * t325, -t234 * t327 - t235 * t325, pkin(3) * t243 + qJ(2) * t149 - qJ(4) * t160 - t371 * t148 - t327 * t152 - t325 * t176, pkin(3) * t245 + qJ(2) * t154 - qJ(4) * t182 - t371 * t153 - t327 * t164 - t325 * t189, pkin(3) * t233 + qJ(2) * t147 - qJ(4) * t162 - t327 * t134 - t325 * t136 - t371 * t146, pkin(3) * t229 + pkin(6) * t370 + qJ(2) * t123 - qJ(4) * t126 - t371 * t122 - t327 * t139;];
tauB_reg = t1;
