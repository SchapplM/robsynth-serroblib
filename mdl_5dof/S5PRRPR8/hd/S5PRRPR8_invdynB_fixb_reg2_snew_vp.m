% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PRRPR8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PRRPR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:44
% EndTime: 2019-12-31 17:42:52
% DurationCPUTime: 5.53s
% Computational Cost: add. (20301->379), mult. (27918->563), div. (0->0), fcn. (19096->10), ass. (0->252)
t377 = qJD(2) + qJD(3);
t375 = t377 ^ 2;
t376 = qJDD(2) + qJDD(3);
t381 = sin(pkin(9));
t383 = cos(pkin(9));
t334 = t383 * t375 + t381 * t376;
t337 = t381 * t375 - t383 * t376;
t386 = sin(qJ(3));
t389 = cos(qJ(3));
t281 = t389 * t334 - t386 * t337;
t285 = t386 * t334 + t389 * t337;
t387 = sin(qJ(2));
t390 = cos(qJ(2));
t245 = t387 * t281 + t390 * t285;
t382 = sin(pkin(8));
t461 = t382 * t245;
t384 = cos(pkin(8));
t460 = t384 * t245;
t358 = t382 * g(1) - t384 * g(2);
t354 = -qJDD(4) + t358;
t297 = qJ(4) * t334 - t383 * t354;
t446 = qJ(4) * t337 - t381 * t354;
t223 = pkin(6) * t281 + t389 * t297 - t386 * t446;
t453 = pkin(6) * t285 + t386 * t297 + t389 * t446;
t171 = pkin(5) * t245 + t387 * t223 + t390 * t453;
t444 = t390 * t281 - t387 * t285;
t172 = pkin(5) * t444 + t390 * t223 - t387 * t453;
t359 = t384 * g(1) + t382 * g(2);
t380 = g(3) - qJDD(1);
t325 = -t387 * t359 + t390 * t380;
t394 = qJDD(2) * pkin(2) - t325;
t326 = -t390 * t359 - t387 * t380;
t437 = qJD(2) ^ 2;
t396 = -t437 * pkin(2) + t326;
t261 = t386 * t394 + t389 * t396;
t258 = -t375 * pkin(3) + t261;
t393 = -t386 * t396 + t389 * t394;
t392 = t376 * pkin(3) + t393;
t214 = t381 * t258 - t383 * t392;
t215 = t383 * t258 + t381 * t392;
t402 = t381 * t214 + t383 * t215;
t183 = t383 * t214 - t381 * t215;
t416 = t389 * t183;
t156 = -t386 * t402 + t416;
t422 = t386 * t183;
t448 = t389 * t402 + t422;
t142 = t387 * t156 + t390 * t448;
t141 = t390 * t156 - t387 * t448;
t339 = t389 * t375 + t386 * t376;
t342 = t386 * t375 - t389 * t376;
t292 = t387 * t339 + t390 * t342;
t457 = t382 * t292;
t456 = t384 * t292;
t305 = pkin(6) * t339 - t389 * t358;
t447 = pkin(6) * t342 - t386 * t358;
t230 = pkin(5) * t292 + t387 * t305 + t390 * t447;
t395 = t390 * t339 - t387 * t342;
t231 = pkin(5) * t395 + t390 * t305 - t387 * t447;
t401 = t389 * t261 - t386 * t393;
t218 = -t386 * t261 - t389 * t393;
t414 = t390 * t218;
t188 = -t387 * t401 + t414;
t420 = t387 * t218;
t189 = t390 * t401 + t420;
t349 = t384 * t358;
t316 = -t382 * t359 + t349;
t385 = sin(qJ(5));
t378 = t385 ^ 2;
t436 = t378 * t375;
t357 = t390 * qJDD(2) - t387 * t437;
t434 = t382 * t357;
t433 = t382 * t358;
t431 = t382 * t376;
t430 = t382 * t380;
t428 = t384 * t357;
t427 = t384 * t380;
t207 = -t376 * pkin(4) - t375 * pkin(7) + t214;
t426 = t385 * t207;
t388 = cos(qJ(5));
t364 = t385 * t375 * t388;
t352 = qJDD(5) + t364;
t425 = t385 * t352;
t353 = qJDD(5) - t364;
t424 = t385 * t353;
t423 = t385 * t376;
t419 = t388 * t207;
t418 = t388 * t352;
t417 = t388 * t353;
t367 = t388 * t376;
t208 = -t375 * pkin(4) + t376 * pkin(7) + t215;
t205 = t388 * t208 - t385 * t354;
t379 = t388 ^ 2;
t413 = t378 + t379;
t412 = qJD(5) * t377;
t411 = t385 * t412;
t410 = t388 * t412;
t204 = t385 * t208 + t388 * t354;
t170 = t385 * t204 + t388 * t205;
t338 = t413 * t376;
t368 = t379 * t375;
t345 = t368 + t436;
t287 = t381 * t338 + t383 * t345;
t288 = t383 * t338 - t381 * t345;
t247 = t389 * t287 + t386 * t288;
t248 = -t386 * t287 + t389 * t288;
t200 = t390 * t247 + t387 * t248;
t201 = -t387 * t247 + t390 * t248;
t409 = -pkin(1) * t200 - pkin(2) * t247 - pkin(3) * t287 - pkin(4) * t345 - pkin(7) * t338 + qJ(1) * t201 - t170;
t408 = pkin(1) * t444 + pkin(2) * t281 + pkin(3) * t334 + qJ(1) * t245 + t215;
t407 = pkin(1) * t245 + pkin(2) * t285 + pkin(3) * t337 - qJ(1) * t444 + t214;
t406 = pkin(1) * t395 + pkin(2) * t339 + qJ(1) * t292 + t261;
t405 = pkin(1) * t292 + pkin(2) * t342 - qJ(1) * t395 - t393;
t356 = t387 * qJDD(2) + t390 * t437;
t404 = -pkin(1) * t356 + qJ(1) * t357 - t326;
t403 = pkin(1) * t357 + qJ(1) * t356 - t325;
t268 = t387 * t325 + t390 * t326;
t317 = -t384 * t359 - t433;
t398 = t381 * t364;
t397 = t383 * t364;
t323 = pkin(5) * t356 - t390 * t358;
t322 = -pkin(5) * t357 - t387 * t358;
t169 = t388 * t204 - t385 * t205;
t267 = t390 * t325 - t387 * t326;
t391 = qJD(5) ^ 2;
t366 = t384 * t376;
t363 = -t368 - t391;
t362 = t368 - t391;
t361 = -t391 - t436;
t360 = t391 - t436;
t355 = pkin(1) * t358;
t348 = t384 * t356;
t347 = t382 * t356;
t346 = t368 - t436;
t333 = t367 - 0.2e1 * t411;
t332 = t367 - t411;
t331 = t410 + t423;
t330 = 0.2e1 * t410 + t423;
t329 = t413 * t412;
t319 = t381 * qJDD(5) + t383 * t329;
t318 = -t383 * qJDD(5) + t381 * t329;
t315 = t388 * t331 - t378 * t412;
t314 = -t385 * t332 - t379 * t412;
t313 = -t385 * t361 - t417;
t312 = -t385 * t360 + t418;
t311 = t388 * t363 - t425;
t310 = t388 * t362 - t424;
t309 = t388 * t361 - t424;
t308 = -t388 * t360 - t425;
t307 = t385 * t363 + t418;
t306 = -t385 * t362 - t417;
t301 = (-t331 - t410) * t385;
t300 = (-t332 + t411) * t388;
t280 = -t385 * t330 + t388 * t333;
t279 = -t388 * t330 - t385 * t333;
t278 = t384 * t395;
t277 = t382 * t395;
t276 = t383 * t312 + t381 * t423;
t275 = t383 * t310 + t381 * t367;
t274 = t381 * t312 - t383 * t423;
t273 = t381 * t310 - t383 * t367;
t272 = t383 * t315 - t398;
t271 = t383 * t314 + t398;
t270 = t381 * t315 + t397;
t269 = t381 * t314 - t397;
t265 = t383 * t313 + t381 * t330;
t264 = t383 * t311 - t381 * t333;
t263 = t381 * t313 - t383 * t330;
t262 = t381 * t311 + t383 * t333;
t257 = -t386 * t318 + t389 * t319;
t256 = t389 * t318 + t386 * t319;
t255 = t383 * t280 - t381 * t346;
t254 = t381 * t280 + t383 * t346;
t250 = t384 * t268 - t433;
t249 = t382 * t268 + t349;
t241 = t384 * t444;
t240 = t382 * t444;
t239 = -t386 * t274 + t389 * t276;
t238 = -t386 * t273 + t389 * t275;
t237 = t389 * t274 + t386 * t276;
t236 = t389 * t273 + t386 * t275;
t235 = -t386 * t270 + t389 * t272;
t234 = -t386 * t269 + t389 * t271;
t233 = t389 * t270 + t386 * t272;
t232 = t389 * t269 + t386 * t271;
t227 = -t386 * t263 + t389 * t265;
t226 = -t386 * t262 + t389 * t264;
t225 = t389 * t263 + t386 * t265;
t224 = t389 * t262 + t386 * t264;
t213 = pkin(2) * t358 + pkin(6) * t401;
t212 = -t387 * t256 + t390 * t257;
t210 = -t386 * t254 + t389 * t255;
t209 = t389 * t254 + t386 * t255;
t203 = -pkin(7) * t309 + t419;
t202 = -pkin(7) * t307 + t426;
t199 = -pkin(4) * t309 + t205;
t198 = -pkin(4) * t307 + t204;
t197 = -t387 * t237 + t390 * t239;
t196 = -t387 * t236 + t390 * t238;
t195 = -t387 * t233 + t390 * t235;
t194 = -t387 * t232 + t390 * t234;
t193 = -t387 * t225 + t390 * t227;
t192 = -t387 * t224 + t390 * t226;
t191 = t390 * t225 + t387 * t227;
t190 = t390 * t224 + t387 * t226;
t186 = t384 * t189 - t433;
t185 = t382 * t189 + t349;
t178 = -t387 * t209 + t390 * t210;
t177 = t384 * t193 + t382 * t309;
t176 = t384 * t192 + t382 * t307;
t175 = t382 * t193 - t384 * t309;
t174 = t382 * t192 - t384 * t307;
t173 = pkin(3) * t354 + qJ(4) * t402;
t167 = pkin(1) * t188 + pkin(2) * t218;
t166 = -qJ(4) * t287 + t383 * t169;
t165 = qJ(4) * t288 + t381 * t169;
t164 = -qJ(4) * t263 - t381 * t199 + t383 * t203;
t163 = -qJ(4) * t262 - t381 * t198 + t383 * t202;
t162 = -pkin(3) * t309 + qJ(4) * t265 + t383 * t199 + t381 * t203;
t161 = -pkin(3) * t307 + qJ(4) * t264 + t383 * t198 + t381 * t202;
t160 = t383 * t170 + t381 * t207;
t159 = t381 * t170 - t383 * t207;
t158 = pkin(5) * t188 + pkin(6) * t414 - t387 * t213;
t153 = -pkin(1) * t191 - pkin(2) * t225 - pkin(3) * t263 + pkin(4) * t330 - pkin(7) * t313 - t426;
t152 = -pkin(1) * t190 - pkin(2) * t224 - pkin(3) * t262 - pkin(4) * t333 - pkin(7) * t311 + t419;
t150 = -pkin(6) * t247 - t386 * t165 + t389 * t166;
t149 = pkin(6) * t248 + t389 * t165 + t386 * t166;
t148 = -t386 * t159 + t389 * t160;
t147 = t389 * t159 + t386 * t160;
t146 = -pkin(6) * t225 - t386 * t162 + t389 * t164;
t145 = -pkin(6) * t224 - t386 * t161 + t389 * t163;
t144 = -pkin(2) * t309 + pkin(6) * t227 + t389 * t162 + t386 * t164;
t143 = -pkin(2) * t307 + pkin(6) * t226 + t389 * t161 + t386 * t163;
t139 = pkin(6) * t156 + qJ(4) * t416 - t386 * t173;
t138 = t384 * t142 - t382 * t354;
t137 = t382 * t142 + t384 * t354;
t136 = pkin(2) * t354 + pkin(6) * t448 + qJ(4) * t422 + t389 * t173;
t135 = -qJ(4) * t159 - (pkin(4) * t381 - pkin(7) * t383) * t169;
t134 = qJ(4) * t160 - (-pkin(4) * t383 - pkin(7) * t381 - pkin(3)) * t169;
t133 = -pkin(5) * t200 - t387 * t149 + t390 * t150;
t132 = pkin(1) * t141 + pkin(2) * t156 + pkin(3) * t183;
t131 = -t387 * t147 + t390 * t148;
t130 = t390 * t147 + t387 * t148;
t129 = t384 * t131 - t169 * t382;
t128 = t382 * t131 + t169 * t384;
t127 = -pkin(5) * t191 - t387 * t144 + t390 * t146;
t126 = -pkin(5) * t190 - t387 * t143 + t390 * t145;
t125 = pkin(5) * t141 - t387 * t136 + t390 * t139;
t124 = -pkin(6) * t147 - t386 * t134 + t389 * t135;
t123 = -pkin(1) * t130 - pkin(2) * t147 - pkin(3) * t159 + pkin(4) * t207 - pkin(7) * t170;
t122 = pkin(2) * t169 + pkin(6) * t148 + t389 * t134 + t386 * t135;
t121 = -pkin(5) * t130 - t387 * t122 + t390 * t124;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t317, 0, 0, 0, 0, 0, 0, -t348, -t428, 0, t250, 0, 0, 0, 0, 0, 0, -t278, t456, 0, t186, 0, 0, 0, 0, 0, 0, -t241, t460, 0, t138, 0, 0, 0, 0, 0, 0, t176, t177, t384 * t201, t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t316, 0, 0, 0, 0, 0, 0, -t347, -t434, 0, t249, 0, 0, 0, 0, 0, 0, -t277, t457, 0, t185, 0, 0, 0, 0, 0, 0, -t240, t461, 0, t137, 0, 0, 0, 0, 0, 0, t174, t175, t382 * t201, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t380, 0, 0, 0, 0, 0, 0, t357, -t356, 0, -t267, 0, 0, 0, 0, 0, 0, -t292, -t395, 0, -t188, 0, 0, 0, 0, 0, 0, -t245, -t444, 0, -t141, 0, 0, 0, 0, 0, 0, t190, t191, t200, t130; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t430, -t427, -t316, -qJ(1) * t316, 0, 0, t428, 0, -t348, t382 * qJDD(2), t384 * t322 + t403 * t382, t384 * t323 + t404 * t382, t384 * t267, -qJ(1) * t249 - (pkin(1) * t382 - pkin(5) * t384) * t267, 0, 0, -t456, 0, -t278, t431, t384 * t230 - t405 * t382, t384 * t231 - t406 * t382, t384 * t188, -qJ(1) * t185 + t384 * t158 - t382 * t167, 0, 0, -t460, 0, -t241, t431, t384 * t171 - t382 * t407, t384 * t172 - t382 * t408, t384 * t141, -qJ(1) * t137 + t384 * t125 - t382 * t132, t384 * t195 - t382 * t301, t384 * t178 - t382 * t279, t384 * t197 - t382 * t308, t384 * t194 - t382 * t300, t384 * t196 - t382 * t306, t384 * t212, -qJ(1) * t174 + t384 * t126 - t382 * t152, -qJ(1) * t175 + t384 * t127 - t382 * t153, t384 * t133 - t382 * t409, -qJ(1) * t128 + t384 * t121 - t382 * t123; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t427, -t430, t317, qJ(1) * t317, 0, 0, t434, 0, -t347, -t384 * qJDD(2), t382 * t322 - t403 * t384, t382 * t323 - t404 * t384, t382 * t267, qJ(1) * t250 - (-pkin(1) * t384 - pkin(5) * t382) * t267, 0, 0, -t457, 0, -t277, -t366, t382 * t230 + t405 * t384, t382 * t231 + t384 * t406, t382 * t188, qJ(1) * t186 + t382 * t158 + t384 * t167, 0, 0, -t461, 0, -t240, -t366, t382 * t171 + t384 * t407, t382 * t172 + t384 * t408, t382 * t141, qJ(1) * t138 + t382 * t125 + t384 * t132, t382 * t195 + t384 * t301, t382 * t178 + t384 * t279, t382 * t197 + t384 * t308, t382 * t194 + t384 * t300, t382 * t196 + t384 * t306, t382 * t212, qJ(1) * t176 + t382 * t126 + t384 * t152, qJ(1) * t177 + t382 * t127 + t384 * t153, t382 * t133 + t384 * t409, qJ(1) * t129 + t382 * t121 + t384 * t123; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t358, t359, 0, 0, 0, 0, t356, 0, t357, 0, -t323, t322, t268, pkin(5) * t268 + t355, 0, 0, t395, 0, -t292, 0, -t231, t230, t189, pkin(5) * t189 + pkin(6) * t420 + t390 * t213 + t355, 0, 0, t444, 0, -t245, 0, -t172, t171, t142, pkin(1) * t354 + pkin(5) * t142 + t390 * t136 + t387 * t139, t390 * t233 + t387 * t235, t390 * t209 + t387 * t210, t390 * t237 + t387 * t239, t390 * t232 + t387 * t234, t390 * t236 + t387 * t238, t390 * t256 + t387 * t257, -pkin(1) * t307 + pkin(5) * t192 + t390 * t143 + t387 * t145, -pkin(1) * t309 + pkin(5) * t193 + t390 * t144 + t387 * t146, pkin(5) * t201 + t390 * t149 + t387 * t150, pkin(1) * t169 + pkin(5) * t131 + t390 * t122 + t387 * t124;];
tauB_reg = t1;
