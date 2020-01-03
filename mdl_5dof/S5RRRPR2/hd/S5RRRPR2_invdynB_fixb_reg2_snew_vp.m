% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRRPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:40
% EndTime: 2020-01-03 12:07:49
% DurationCPUTime: 7.32s
% Computational Cost: add. (37767->383), mult. (48058->558), div. (0->0), fcn. (28562->10), ass. (0->254)
t427 = qJD(1) + qJD(2);
t420 = qJD(3) + t427;
t418 = t420 ^ 2;
t426 = qJDD(1) + qJDD(2);
t419 = qJDD(3) + t426;
t431 = sin(pkin(9));
t432 = cos(pkin(9));
t376 = t418 * t432 + t419 * t431;
t379 = t418 * t431 - t419 * t432;
t434 = sin(qJ(3));
t438 = cos(qJ(3));
t320 = t376 * t438 - t379 * t434;
t430 = g(1) - qJDD(4);
t358 = qJ(4) * t376 - t430 * t432;
t490 = qJ(4) * t379 - t430 * t431;
t270 = pkin(7) * t320 + t358 * t438 - t434 * t490;
t324 = t376 * t434 + t379 * t438;
t435 = sin(qJ(2));
t439 = cos(qJ(2));
t279 = t320 * t439 - t324 * t435;
t504 = pkin(7) * t324 + t358 * t434 + t438 * t490;
t214 = pkin(6) * t279 + t270 * t439 - t435 * t504;
t283 = t320 * t435 + t324 * t439;
t436 = sin(qJ(1));
t440 = cos(qJ(1));
t240 = t279 * t440 - t283 * t436;
t519 = pkin(6) * t283 + t270 * t435 + t439 * t504;
t528 = -pkin(5) * t240 - t214 * t440 + t436 * t519;
t508 = t279 * t436 + t283 * t440;
t527 = pkin(5) * t508 + t214 * t436 + t440 * t519;
t411 = g(2) * t440 + g(3) * t436;
t400 = qJDD(1) * pkin(1) - t411;
t410 = g(2) * t436 - g(3) * t440;
t477 = qJD(1) ^ 2;
t401 = -pkin(1) * t477 - t410;
t345 = -t400 * t439 + t401 * t435;
t336 = pkin(2) * t426 - t345;
t346 = t400 * t435 + t401 * t439;
t425 = t427 ^ 2;
t337 = -pkin(2) * t425 + t346;
t297 = t336 * t434 + t337 * t438;
t294 = -pkin(3) * t418 + t297;
t445 = t336 * t438 - t337 * t434;
t443 = pkin(3) * t419 + t445;
t246 = t294 * t431 - t432 * t443;
t247 = t294 * t432 + t431 * t443;
t452 = t246 * t431 + t247 * t432;
t205 = t246 * t432 - t247 * t431;
t460 = t438 * t205;
t179 = -t434 * t452 + t460;
t465 = t434 * t205;
t494 = t438 * t452 + t465;
t166 = t179 * t439 - t435 * t494;
t514 = t179 * t435 + t439 * t494;
t524 = t166 * t436 + t440 * t514;
t151 = -t166 * t440 + t436 * t514;
t382 = t418 * t438 + t419 * t434;
t385 = t418 * t434 - t419 * t438;
t328 = t382 * t439 - t385 * t435;
t364 = pkin(7) * t382 - g(1) * t438;
t491 = pkin(7) * t385 - g(1) * t434;
t278 = pkin(6) * t328 + t364 * t439 - t435 * t491;
t332 = t382 * t435 + t385 * t439;
t290 = t328 * t440 - t332 * t436;
t505 = pkin(6) * t332 + t364 * t435 + t439 * t491;
t521 = -pkin(5) * t290 - t278 * t440 + t436 * t505;
t488 = t328 * t436 + t332 * t440;
t520 = pkin(5) * t488 + t278 * t436 + t440 * t505;
t451 = t297 * t438 - t434 * t445;
t251 = -t297 * t434 - t438 * t445;
t459 = t439 * t251;
t209 = -t435 * t451 + t459;
t464 = t435 * t251;
t495 = t439 * t451 + t464;
t513 = t209 * t436 + t440 * t495;
t182 = -t209 * t440 + t436 * t495;
t393 = t425 * t439 + t426 * t435;
t396 = t425 * t435 - t426 * t439;
t341 = t393 * t440 - t396 * t436;
t369 = pkin(6) * t393 - g(1) * t439;
t492 = pkin(6) * t396 - g(1) * t435;
t507 = -pkin(5) * t341 - t369 * t440 + t436 * t492;
t444 = t393 * t436 + t396 * t440;
t506 = pkin(5) * t444 + t369 * t436 + t440 * t492;
t450 = t345 * t435 + t346 * t439;
t303 = t345 * t439 - t346 * t435;
t458 = t303 * t440;
t254 = t436 * t450 - t458;
t463 = t303 * t436;
t493 = t440 * t450 + t463;
t433 = sin(qJ(5));
t428 = t433 ^ 2;
t472 = t428 * t418;
t241 = -pkin(4) * t419 - pkin(8) * t418 + t246;
t469 = t433 * t241;
t437 = cos(qJ(5));
t406 = t437 * t418 * t433;
t397 = qJDD(5) + t406;
t468 = t433 * t397;
t398 = qJDD(5) - t406;
t467 = t433 * t398;
t466 = t433 * t419;
t462 = t437 * t241;
t461 = t437 * t398;
t412 = t437 * t419;
t242 = -pkin(4) * t418 + pkin(8) * t419 + t247;
t235 = t242 * t437 - t430 * t433;
t429 = t437 ^ 2;
t457 = t428 + t429;
t456 = qJD(5) * t420;
t455 = t433 * t456;
t454 = t437 * t456;
t407 = -qJDD(1) * t436 - t440 * t477;
t453 = pkin(5) * t407 + g(1) * t440;
t234 = t242 * t433 + t430 * t437;
t199 = t234 * t433 + t235 * t437;
t359 = -t410 * t436 - t411 * t440;
t447 = t431 * t406;
t446 = t432 * t406;
t198 = t234 * t437 - t235 * t433;
t360 = t410 * t440 - t411 * t436;
t442 = qJD(5) ^ 2;
t441 = pkin(1) * g(1);
t413 = t429 * t418;
t408 = qJDD(1) * t440 - t436 * t477;
t405 = -t413 - t442;
t404 = t413 - t442;
t403 = -t442 - t472;
t402 = t442 - t472;
t389 = t437 * t397;
t388 = pkin(5) * t408 + g(1) * t436;
t387 = t413 - t472;
t386 = t413 + t472;
t380 = t457 * t419;
t375 = t412 - 0.2e1 * t455;
t374 = t412 - t455;
t373 = t454 + t466;
t372 = 0.2e1 * t454 + t466;
t371 = t457 * t456;
t354 = qJDD(5) * t431 + t371 * t432;
t353 = -qJDD(5) * t432 + t371 * t431;
t352 = -t403 * t433 - t461;
t351 = -t402 * t433 + t389;
t350 = t405 * t437 - t468;
t349 = t404 * t437 - t467;
t348 = t403 * t437 - t467;
t347 = t405 * t433 + t389;
t344 = t373 * t437 - t428 * t456;
t343 = -t374 * t433 - t429 * t456;
t327 = t380 * t432 - t386 * t431;
t326 = t380 * t431 + t386 * t432;
t319 = -t372 * t433 + t375 * t437;
t318 = t351 * t432 + t431 * t466;
t317 = t349 * t432 + t412 * t431;
t316 = t351 * t431 - t432 * t466;
t315 = t349 * t431 - t412 * t432;
t314 = t344 * t432 - t447;
t313 = t343 * t432 + t447;
t312 = t344 * t431 + t446;
t311 = t343 * t431 - t446;
t310 = t352 * t432 + t372 * t431;
t309 = t350 * t432 - t375 * t431;
t308 = t352 * t431 - t372 * t432;
t307 = t350 * t431 + t375 * t432;
t306 = -t353 * t434 + t354 * t438;
t305 = t353 * t438 + t354 * t434;
t302 = t319 * t432 - t387 * t431;
t299 = t319 * t431 + t387 * t432;
t298 = pkin(6) * t450 + t441;
t286 = -t326 * t434 + t327 * t438;
t285 = t326 * t438 + t327 * t434;
t274 = -t316 * t434 + t318 * t438;
t273 = -t315 * t434 + t317 * t438;
t272 = t316 * t438 + t318 * t434;
t271 = t315 * t438 + t317 * t434;
t266 = -t312 * t434 + t314 * t438;
t265 = -t311 * t434 + t313 * t438;
t264 = t312 * t438 + t314 * t434;
t263 = t311 * t438 + t313 * t434;
t262 = -t308 * t434 + t310 * t438;
t261 = -t307 * t434 + t309 * t438;
t260 = t308 * t438 + t310 * t434;
t259 = t307 * t438 + t309 * t434;
t258 = -t305 * t435 + t306 * t439;
t257 = t305 * t439 + t306 * t435;
t255 = -t299 * t434 + t302 * t438;
t253 = t299 * t438 + t302 * t434;
t248 = pkin(2) * g(1) + pkin(7) * t451;
t244 = -t285 * t435 + t286 * t439;
t243 = t285 * t439 + t286 * t435;
t232 = -t272 * t435 + t274 * t439;
t231 = -t271 * t435 + t273 * t439;
t230 = t272 * t439 + t274 * t435;
t229 = t271 * t439 + t273 * t435;
t228 = -t264 * t435 + t266 * t439;
t227 = -t263 * t435 + t265 * t439;
t226 = t264 * t439 + t266 * t435;
t225 = t263 * t439 + t265 * t435;
t224 = -pkin(8) * t348 + t462;
t223 = -pkin(8) * t347 + t469;
t222 = -pkin(4) * t348 + t235;
t221 = -pkin(4) * t347 + t234;
t220 = -t260 * t435 + t262 * t439;
t219 = -t259 * t435 + t261 * t439;
t218 = t260 * t439 + t262 * t435;
t217 = t259 * t439 + t261 * t435;
t216 = -t253 * t435 + t255 * t439;
t215 = t253 * t439 + t255 * t435;
t202 = pkin(3) * t430 + qJ(4) * t452;
t201 = t243 * t436 - t244 * t440;
t200 = t243 * t440 + t244 * t436;
t196 = -qJ(4) * t326 + t198 * t432;
t195 = qJ(4) * t327 + t198 * t431;
t194 = t218 * t436 - t220 * t440;
t193 = t217 * t436 - t219 * t440;
t192 = t218 * t440 + t220 * t436;
t191 = t217 * t440 + t219 * t436;
t190 = -qJ(4) * t308 - t222 * t431 + t224 * t432;
t189 = -qJ(4) * t307 - t221 * t431 + t223 * t432;
t188 = -pkin(3) * t348 + qJ(4) * t310 + t222 * t432 + t224 * t431;
t187 = -pkin(3) * t347 + qJ(4) * t309 + t221 * t432 + t223 * t431;
t186 = t199 * t432 + t241 * t431;
t185 = t199 * t431 - t241 * t432;
t184 = pkin(6) * t209 + pkin(7) * t459 - t248 * t435;
t181 = pkin(6) * t495 + pkin(7) * t464 + t248 * t439 + t441;
t176 = -pkin(7) * t285 - t195 * t434 + t196 * t438;
t175 = pkin(7) * t286 + t195 * t438 + t196 * t434;
t174 = -pkin(7) * t260 - t188 * t434 + t190 * t438;
t173 = -pkin(7) * t259 - t187 * t434 + t189 * t438;
t172 = -t185 * t434 + t186 * t438;
t171 = t185 * t438 + t186 * t434;
t170 = -pkin(2) * t348 + pkin(7) * t262 + t188 * t438 + t190 * t434;
t169 = -pkin(2) * t347 + pkin(7) * t261 + t187 * t438 + t189 * t434;
t168 = -qJ(4) * t185 - (pkin(4) * t431 - pkin(8) * t432) * t198;
t163 = pkin(7) * t179 + qJ(4) * t460 - t202 * t434;
t162 = pkin(2) * t430 + pkin(7) * t494 + qJ(4) * t465 + t202 * t438;
t161 = qJ(4) * t186 - (-pkin(4) * t432 - pkin(8) * t431 - pkin(3)) * t198;
t160 = -pkin(6) * t243 - t175 * t435 + t176 * t439;
t159 = pkin(6) * t244 + t175 * t439 + t176 * t435;
t158 = -t171 * t435 + t172 * t439;
t157 = t171 * t439 + t172 * t435;
t156 = -pkin(6) * t218 - t170 * t435 + t174 * t439;
t155 = -pkin(6) * t217 - t169 * t435 + t173 * t439;
t154 = -pkin(1) * t348 + pkin(6) * t220 + t170 * t439 + t174 * t435;
t153 = -pkin(1) * t347 + pkin(6) * t219 + t169 * t439 + t173 * t435;
t150 = pkin(6) * t166 - t162 * t435 + t163 * t439;
t149 = pkin(1) * t430 + pkin(6) * t514 + t162 * t439 + t163 * t435;
t148 = -pkin(7) * t171 - t161 * t434 + t168 * t438;
t147 = t157 * t436 - t158 * t440;
t146 = t157 * t440 + t158 * t436;
t145 = pkin(2) * t198 + pkin(7) * t172 + t161 * t438 + t168 * t434;
t144 = -pkin(6) * t157 - t145 * t435 + t148 * t439;
t143 = pkin(1) * t198 + pkin(6) * t158 + t145 * t439 + t148 * t435;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t430, 0, 0, 0, 0, 0, 0, t347, t348, 0, -t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t408, t407, 0, t359, 0, 0, 0, 0, 0, 0, -t444, -t341, 0, t254, 0, 0, 0, 0, 0, 0, -t488, -t290, 0, t182, 0, 0, 0, 0, 0, 0, -t508, -t240, 0, t151, 0, 0, 0, 0, 0, 0, t191, t192, t200, t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t407, t408, 0, t360, 0, 0, 0, 0, 0, 0, t341, -t444, 0, -t493, 0, 0, 0, 0, 0, 0, t290, -t488, 0, -t513, 0, 0, 0, 0, 0, 0, t240, -t508, 0, -t524, 0, 0, 0, 0, 0, 0, t193, t194, t201, t147; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t411, t410, 0, 0, 0, 0, 0, 0, 0, t426, -pkin(1) * t396 - t345, -pkin(1) * t393 - t346, 0, -pkin(1) * t303, 0, 0, 0, 0, 0, t419, -pkin(1) * t332 - pkin(2) * t385 + t445, -pkin(1) * t328 - pkin(2) * t382 - t297, 0, -pkin(1) * t209 - pkin(2) * t251, 0, 0, 0, 0, 0, t419, -pkin(1) * t283 - pkin(2) * t324 - pkin(3) * t379 - t246, -pkin(1) * t279 - pkin(2) * t320 - pkin(3) * t376 - t247, 0, -pkin(1) * t166 - pkin(2) * t179 - pkin(3) * t205, (t373 + t454) * t433, t372 * t437 + t375 * t433, t402 * t437 + t468, (t374 - t455) * t437, t404 * t433 + t461, 0, pkin(1) * t217 + pkin(2) * t259 + pkin(3) * t307 + pkin(4) * t375 + pkin(8) * t350 - t462, pkin(1) * t218 + pkin(2) * t260 + pkin(3) * t308 - pkin(4) * t372 + pkin(8) * t352 + t469, pkin(1) * t243 + pkin(2) * t285 + pkin(3) * t326 + pkin(4) * t386 + pkin(8) * t380 + t199, pkin(1) * t157 + pkin(2) * t171 + pkin(3) * t185 - pkin(4) * t241 + pkin(8) * t199; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t407, 0, t408, 0, t453, -t388, -t360, -pkin(5) * t360, 0, 0, t341, 0, -t444, 0, t507, t506, t493, pkin(5) * t493 + pkin(6) * t463 + t298 * t440, 0, 0, t290, 0, -t488, 0, t521, t520, t513, pkin(5) * t513 + t181 * t440 + t184 * t436, 0, 0, t240, 0, -t508, 0, t528, t527, t524, pkin(5) * t524 + t149 * t440 + t150 * t436, t226 * t440 + t228 * t436, t215 * t440 + t216 * t436, t230 * t440 + t232 * t436, t225 * t440 + t227 * t436, t229 * t440 + t231 * t436, t257 * t440 + t258 * t436, -pkin(5) * t193 + t153 * t440 + t155 * t436, -pkin(5) * t194 + t154 * t440 + t156 * t436, -pkin(5) * t201 + t159 * t440 + t160 * t436, -pkin(5) * t147 + t143 * t440 + t144 * t436; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t408, 0, -t407, 0, t388, t453, t359, pkin(5) * t359, 0, 0, t444, 0, t341, 0, -t506, t507, t254, pkin(5) * t254 - pkin(6) * t458 + t298 * t436, 0, 0, t488, 0, t290, 0, -t520, t521, t182, pkin(5) * t182 + t181 * t436 - t184 * t440, 0, 0, t508, 0, t240, 0, -t527, t528, t151, pkin(5) * t151 + t149 * t436 - t150 * t440, t226 * t436 - t228 * t440, t215 * t436 - t216 * t440, t230 * t436 - t232 * t440, t225 * t436 - t227 * t440, t229 * t436 - t231 * t440, t257 * t436 - t258 * t440, pkin(5) * t191 + t153 * t436 - t155 * t440, pkin(5) * t192 + t154 * t436 - t156 * t440, pkin(5) * t200 + t159 * t436 - t160 * t440, pkin(5) * t146 + t143 * t436 - t144 * t440;];
tauB_reg = t1;
