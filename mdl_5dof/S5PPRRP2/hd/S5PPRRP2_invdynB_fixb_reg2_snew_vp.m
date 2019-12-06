% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5PPRRP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:13
% EndTime: 2019-12-05 15:09:22
% DurationCPUTime: 5.75s
% Computational Cost: add. (8363->362), mult. (15751->503), div. (0->0), fcn. (10560->8), ass. (0->267)
t366 = qJD(4) ^ 2;
t362 = sin(qJ(4));
t353 = t362 ^ 2;
t367 = qJD(3) ^ 2;
t425 = t353 * t367;
t340 = t366 + t425;
t364 = cos(qJ(4));
t345 = t362 * t367 * t364;
t339 = qJDD(4) - t345;
t406 = t364 * t339;
t302 = -t362 * t340 + t406;
t398 = qJD(3) * qJD(4);
t390 = t364 * t398;
t396 = t362 * qJDD(3);
t328 = 0.2e1 * t390 + t396;
t363 = sin(qJ(3));
t365 = cos(qJ(3));
t253 = t363 * t302 + t365 * t328;
t256 = t365 * t302 - t363 * t328;
t357 = sin(pkin(8));
t359 = cos(pkin(8));
t208 = t357 * t253 - t359 * t256;
t412 = t362 * t339;
t296 = t364 * t340 + t412;
t358 = sin(pkin(7));
t360 = cos(pkin(7));
t176 = t358 * t208 + t360 * t296;
t470 = qJ(1) * t176;
t179 = t360 * t208 - t358 * t296;
t469 = qJ(1) * t179;
t204 = t359 * t253 + t357 * t256;
t468 = qJ(2) * t204;
t467 = -pkin(1) * t204 - pkin(2) * t253 - pkin(6) * t302;
t466 = -pkin(1) * t296 - qJ(2) * t208;
t414 = t362 * t328;
t392 = t362 * t398;
t394 = t364 * qJDD(3);
t373 = 0.2e1 * t392 - t394;
t442 = t373 * t364;
t273 = t442 + t414;
t354 = t364 ^ 2;
t335 = (t353 - t354) * t367;
t248 = t363 * t273 + t365 * t335;
t250 = t365 * t273 - t363 * t335;
t195 = t357 * t248 - t359 * t250;
t271 = t364 * t328 - t362 * t373;
t465 = t358 * t195 - t360 * t271;
t464 = t360 * t195 + t358 * t271;
t424 = t354 * t367;
t342 = -t366 + t424;
t300 = -t364 * t342 + t412;
t393 = t365 * qJDD(3);
t259 = t363 * t300 + t364 * t393;
t262 = t365 * t300 - t363 * t394;
t216 = t357 * t259 - t359 * t262;
t294 = t362 * t342 + t406;
t463 = t358 * t216 - t360 * t294;
t462 = t360 * t216 + t358 * t294;
t461 = pkin(5) * t253;
t459 = -pkin(2) * t296 + pkin(5) * t256;
t458 = t359 * t248 + t357 * t250;
t457 = t359 * t259 + t357 * t262;
t395 = t363 * qJDD(3);
t332 = t365 * t367 + t395;
t333 = -t363 * t367 + t393;
t384 = -t357 * t332 + t359 * t333;
t456 = t358 * t384;
t455 = t360 * t384;
t337 = t360 * g(1) + t358 * g(2);
t355 = g(3) - qJDD(1);
t311 = -t357 * t337 + t359 * t355;
t312 = -t359 * t337 - t357 * t355;
t245 = t365 * t311 + t363 * t312;
t246 = -t363 * t311 + t365 * t312;
t385 = t363 * t245 + t365 * t246;
t193 = t365 * t245 - t363 * t246;
t418 = t359 * t193;
t158 = -t357 * t385 + t418;
t423 = t357 * t193;
t159 = t359 * t385 + t423;
t336 = t358 * g(1) - t360 * g(2);
t330 = -qJDD(2) + t336;
t287 = pkin(5) * t332 - t365 * t330;
t377 = -pkin(5) * t333 - t363 * t330;
t209 = -qJ(2) * t384 + t357 * t287 + t359 * t377;
t454 = 2 * qJD(5);
t451 = pkin(3) * t296;
t450 = pkin(6) * t296;
t437 = t359 * t332 + t357 * t333;
t210 = qJ(2) * t437 + t359 * t287 - t357 * t377;
t426 = qJ(5) * t362;
t434 = pkin(4) * t364;
t379 = -t426 - t434;
t327 = t379 * qJD(3);
t236 = -t367 * pkin(3) + qJDD(3) * pkin(6) + t246;
t402 = -t364 * t236 + t362 * t330;
t376 = t364 * qJD(3) * t327 + qJDD(4) * qJ(5) + (qJD(4) * t454) - t402;
t399 = qJD(3) * t362;
t436 = t327 * t399 + qJDD(5);
t343 = -t366 - t424;
t338 = qJDD(4) + t345;
t413 = t362 * t338;
t299 = t364 * t343 - t413;
t252 = t363 * t299 - t365 * t373;
t433 = pkin(5) * t252;
t400 = t353 + t354;
t331 = t400 * qJDD(3);
t334 = t400 * t367;
t280 = t363 * t331 + t365 * t334;
t432 = pkin(5) * t280;
t407 = t364 * t338;
t293 = t362 * t343 + t407;
t431 = pkin(6) * t293;
t255 = t365 * t299 + t363 * t373;
t206 = -t357 * t252 + t359 * t255;
t174 = t358 * t206 - t360 * t293;
t430 = qJ(1) * t174;
t281 = t365 * t331 - t363 * t334;
t232 = -t357 * t280 + t359 * t281;
t429 = qJ(1) * t232;
t203 = t359 * t252 + t357 * t255;
t428 = qJ(2) * t203;
t231 = t359 * t280 + t357 * t281;
t427 = qJ(2) * t231;
t420 = t358 * t330;
t419 = t358 * t355;
t417 = t359 * t330;
t229 = t360 * t232;
t315 = t360 * t330;
t416 = t360 * t355;
t235 = -qJDD(3) * pkin(3) - t367 * pkin(6) + t245;
t415 = t362 * t235;
t409 = t364 * t235;
t403 = -pkin(1) * t293 + qJ(2) * t206;
t223 = t362 * t236 + t364 * t330;
t401 = t334 - t366;
t388 = -pkin(2) * t293 + pkin(5) * t255;
t387 = pkin(1) * t437 + pkin(2) * t332 - qJ(1) * t384 + t246;
t386 = -pkin(1) * t384 - pkin(2) * t333 - qJ(1) * t437 + t245;
t244 = t357 * t311 + t359 * t312;
t283 = -t358 * t336 - t360 * t337;
t382 = t363 * t345;
t381 = t365 * t345;
t199 = -pkin(3) * t293 + t223;
t378 = pkin(4) * t362 - qJ(5) * t364;
t169 = t364 * t223 + t362 * t402;
t170 = t362 * t223 - t364 * t402;
t243 = t359 * t311 - t357 * t312;
t282 = t360 * t336 - t358 * t337;
t375 = t390 + t396;
t374 = -t392 + t394;
t372 = -pkin(1) * t203 - pkin(2) * t252 - pkin(6) * t299;
t371 = -qJDD(4) * pkin(4) + t223 + t436;
t370 = -t374 * pkin(4) + t235 + (-t375 - t390) * qJ(5);
t369 = -pkin(1) * t231 - pkin(2) * t280 - pkin(3) * t334 - pkin(6) * t331;
t368 = t399 * t454 - t370;
t341 = t366 - t425;
t326 = t378 * qJDD(3);
t325 = pkin(1) * t330;
t324 = t400 * t398;
t310 = t363 * qJDD(4) + t365 * t324;
t309 = -t353 * t398 + t364 * t375;
t308 = -t365 * qJDD(4) + t363 * t324;
t307 = -t354 * t398 - t362 * t374;
t301 = -t362 * t341 + t407;
t295 = -t364 * t341 - t413;
t274 = pkin(5) * t281;
t269 = t360 * t437;
t268 = t358 * t437;
t267 = t365 * t309 - t382;
t266 = t365 * t307 + t382;
t265 = t363 * t309 + t381;
t264 = t363 * t307 - t381;
t263 = t365 * t301 + t362 * t395;
t260 = t363 * t301 - t362 * t393;
t241 = -t357 * t308 + t359 * t310;
t240 = t359 * t308 + t357 * t310;
t238 = t360 * t241;
t237 = t358 * t241;
t230 = qJ(2) * t232;
t228 = t358 * t232;
t227 = qJ(1) * t229;
t226 = t360 * t244 - t420;
t225 = t358 * t244 + t315;
t222 = -t357 * t265 + t359 * t267;
t221 = -t357 * t264 + t359 * t266;
t220 = t359 * t265 + t357 * t267;
t219 = t359 * t264 + t357 * t266;
t218 = -t357 * t260 + t359 * t263;
t215 = t359 * t260 + t357 * t263;
t214 = t409 + t450;
t213 = t415 - t431;
t201 = t366 * qJ(5) - t371;
t200 = -t402 + t451;
t198 = -t366 * pkin(4) + t376;
t197 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t399 + t370;
t190 = t401 * qJ(5) + t371;
t189 = t360 * t222 + t358 * t414;
t188 = t360 * t221 - t358 * t442;
t187 = t358 * t222 - t360 * t414;
t186 = t358 * t221 + t360 * t442;
t185 = t360 * t218 - t358 * t295;
t184 = t358 * t218 + t360 * t295;
t183 = t401 * pkin(4) + t376;
t182 = (-t373 - t392) * pkin(4) + t368;
t181 = -pkin(4) * t392 + qJ(5) * t328 + t368;
t180 = pkin(2) * t330 + pkin(5) * t385;
t177 = t360 * t206 + t358 * t293;
t173 = qJ(1) * t177;
t172 = (-t343 - t366) * qJ(5) + (-qJDD(4) - t338) * pkin(4) + t199 + t436;
t171 = -t451 - qJ(5) * t339 + (-t340 + t366) * pkin(4) - t376;
t167 = -pkin(4) * t414 + t364 * t181 - t450;
t166 = -qJ(5) * t442 - t362 * t182 - t431;
t165 = t365 * t169 - t432;
t164 = t363 * t169 + t274;
t163 = t364 * t198 - t362 * t201;
t162 = t362 * t198 + t364 * t201;
t161 = t365 * t170 + t363 * t235;
t160 = t363 * t170 - t365 * t235;
t156 = -t362 * t183 + t364 * t190;
t155 = -t363 * t200 + t365 * t214 + t461;
t154 = -t363 * t199 + t365 * t213 - t433;
t153 = t360 * t159 - t420;
t152 = t358 * t159 + t315;
t151 = pkin(3) * t328 - t415 - t467;
t150 = pkin(3) * t373 + t372 + t409;
t149 = t365 * t200 + t363 * t214 - t459;
t148 = t365 * t199 + t363 * t213 + t388;
t147 = -t170 + t369;
t146 = t365 * t156 - t363 * t326 - t432;
t145 = t363 * t156 + t365 * t326 + t274;
t144 = -t362 * t181 + (-pkin(3) - t434) * t328 + t467;
t143 = -t364 * t182 - (-pkin(3) - t426) * t373 + t372;
t142 = t365 * t163 + t363 * t197;
t141 = t363 * t163 - t365 * t197;
t140 = pkin(1) * t158 + pkin(2) * t193;
t139 = -t364 * t183 - t362 * t190 + t369;
t138 = t365 * t166 - t363 * t172 - t433;
t137 = t365 * t167 - t363 * t171 - t461;
t136 = t363 * t166 + t365 * t172 + t388;
t135 = t363 * t167 + t365 * t171 + t459;
t134 = -pkin(3) * t162 - pkin(4) * t201 - qJ(5) * t198;
t133 = -pkin(6) * t162 + t378 * t197;
t132 = -t357 * t160 + t359 * t161;
t131 = t359 * t160 + t357 * t161;
t130 = -t357 * t164 + t359 * t165 - t427;
t129 = pkin(5) * t418 + qJ(2) * t158 - t357 * t180;
t128 = -pkin(5) * t160 - (pkin(3) * t363 - pkin(6) * t365) * t169;
t127 = -t357 * t149 + t359 * t155 + t468;
t126 = -t357 * t148 + t359 * t154 - t428;
t125 = t360 * t132 - t169 * t358;
t124 = t358 * t132 + t169 * t360;
t123 = -t357 * t145 + t359 * t146 - t427;
t122 = -t357 * t141 + t359 * t142;
t121 = t359 * t141 + t357 * t142;
t120 = pkin(5) * t161 - (-pkin(3) * t365 - pkin(6) * t363 - pkin(2)) * t169;
t119 = -t357 * t136 + t359 * t138 - t428;
t118 = -t357 * t135 + t359 * t137 - t468;
t117 = t360 * t122 + t358 * t162;
t116 = t358 * t122 - t360 * t162;
t115 = -pkin(1) * t131 - pkin(2) * t160 + pkin(3) * t235 - pkin(6) * t170;
t114 = -pkin(5) * t141 + t365 * t133 - t363 * t134;
t113 = -pkin(2) * t162 + pkin(5) * t142 + t363 * t133 + t365 * t134;
t112 = -pkin(1) * t121 - pkin(2) * t141 - pkin(6) * t163 + (pkin(3) - t379) * t197;
t111 = -qJ(2) * t131 - t357 * t120 + t359 * t128;
t110 = -qJ(2) * t121 - t357 * t113 + t359 * t114;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t283, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, 0, 0, 0, 0, 0, 0, -t269, -t455, 0, t153, 0, 0, 0, 0, 0, 0, t177, t179, t229, t125, 0, 0, 0, 0, 0, 0, t177, t229, -t179, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t282, 0, 0, 0, 0, 0, 0, 0, 0, 0, t225, 0, 0, 0, 0, 0, 0, -t268, -t456, 0, t152, 0, 0, 0, 0, 0, 0, t174, t176, t228, t124, 0, 0, 0, 0, 0, 0, t174, t228, -t176, t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243, 0, 0, 0, 0, 0, 0, t384, -t437, 0, -t158, 0, 0, 0, 0, 0, 0, t203, -t204, t231, t131, 0, 0, 0, 0, 0, 0, t203, t231, t204, t121; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t419, -t416, -t282, -qJ(1) * t282, 0, 0, 0, 0, 0, 0, -t358 * t311 - t357 * t315, -t358 * t312 - t359 * t315, t360 * t243, -qJ(1) * t225 - (pkin(1) * t358 - qJ(2) * t360) * t243, 0, 0, t455, 0, -t269, t358 * qJDD(3), t360 * t209 - t386 * t358, t360 * t210 - t387 * t358, t360 * t158, -qJ(1) * t152 + t360 * t129 - t358 * t140, t189, t464, t185, t188, t462, t238, t360 * t126 - t358 * t150 - t430, t360 * t127 - t358 * t151 - t470, t360 * t130 + (-t147 - t429) * t358, -qJ(1) * t124 + t360 * t111 - t358 * t115, t189, t185, -t464, t238, -t462, t188, t360 * t119 - t358 * t143 - t430, t360 * t123 + (-t139 - t429) * t358, t360 * t118 - t358 * t144 + t470, -qJ(1) * t116 + t360 * t110 - t358 * t112; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t416, -t419, t283, qJ(1) * t283, 0, 0, 0, 0, 0, 0, t360 * t311 - t357 * t420, t360 * t312 - t358 * t417, t358 * t243, qJ(1) * t226 - (-pkin(1) * t360 - qJ(2) * t358) * t243, 0, 0, t456, 0, -t268, -t360 * qJDD(3), t358 * t209 + t386 * t360, t358 * t210 + t387 * t360, t358 * t158, qJ(1) * t153 + t358 * t129 + t360 * t140, t187, t465, t184, t186, t463, t237, t358 * t126 + t360 * t150 + t173, t358 * t127 + t360 * t151 + t469, t358 * t130 + t360 * t147 + t227, qJ(1) * t125 + t358 * t111 + t360 * t115, t187, t184, -t465, t237, -t463, t186, t358 * t119 + t360 * t143 + t173, t358 * t123 + t360 * t139 + t227, t358 * t118 + t360 * t144 - t469, qJ(1) * t117 + t358 * t110 + t360 * t112; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t336, t337, 0, 0, 0, 0, 0, 0, 0, 0, t417, -t357 * t330, t244, qJ(2) * t244 + t325, 0, 0, t437, 0, t384, 0, -t210, t209, t159, pkin(5) * t423 + qJ(2) * t159 + t359 * t180 + t325, t220, -t458, t215, t219, -t457, t240, t359 * t148 + t357 * t154 + t403, t359 * t149 + t357 * t155 - t466, t359 * t164 + t357 * t165 + t230, pkin(1) * t169 + qJ(2) * t132 + t359 * t120 + t357 * t128, t220, t215, t458, t240, t457, t219, t359 * t136 + t357 * t138 + t403, t359 * t145 + t357 * t146 + t230, t359 * t135 + t357 * t137 + t466, -pkin(1) * t162 + qJ(2) * t122 + t359 * t113 + t357 * t114;];
tauB_reg = t1;
