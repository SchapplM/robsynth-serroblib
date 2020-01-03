% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:46
% EndTime: 2019-12-31 19:49:55
% DurationCPUTime: 6.78s
% Computational Cost: add. (17756->392), mult. (25393->520), div. (0->0), fcn. (14266->8), ass. (0->278)
t460 = qJD(4) ^ 2;
t447 = qJD(1) + qJD(2);
t445 = t447 ^ 2;
t454 = sin(qJ(4));
t449 = t454 ^ 2;
t516 = t449 * t445;
t427 = t460 + t516;
t457 = cos(qJ(4));
t432 = t457 * t445 * t454;
t424 = qJDD(4) - t432;
t498 = t457 * t424;
t375 = -t454 * t427 + t498;
t491 = qJD(4) * t447;
t486 = t457 * t491;
t446 = qJDD(1) + qJDD(2);
t502 = t454 * t446;
t403 = 0.2e1 * t486 + t502;
t452 = sin(pkin(8));
t453 = cos(pkin(8));
t326 = t452 * t375 + t453 * t403;
t329 = t453 * t375 - t452 * t403;
t455 = sin(qJ(2));
t458 = cos(qJ(2));
t275 = t458 * t326 + t455 * t329;
t279 = t455 * t326 - t458 * t329;
t456 = sin(qJ(1));
t459 = cos(qJ(1));
t227 = t459 * t275 - t456 * t279;
t573 = pkin(5) * t227;
t231 = t456 * t275 + t459 * t279;
t572 = pkin(5) * t231;
t571 = pkin(6) * t275;
t570 = pkin(1) * t275 + pkin(2) * t326 + pkin(7) * t375;
t505 = t454 * t424;
t369 = t457 * t427 + t505;
t569 = -pkin(1) * t369 - pkin(6) * t279;
t487 = t454 * t491;
t497 = t457 * t446;
t404 = -0.2e1 * t487 + t497;
t364 = t457 * t404;
t365 = t454 * t403;
t345 = -t364 + t365;
t450 = t457 ^ 2;
t421 = (t449 - t450) * t445;
t318 = t452 * t345 + t453 * t421;
t320 = t453 * t345 - t452 * t421;
t265 = t458 * t318 + t455 * t320;
t266 = t455 * t318 - t458 * t320;
t568 = t459 * t265 - t456 * t266;
t567 = t456 * t265 + t459 * t266;
t515 = t450 * t445;
t429 = -t460 + t515;
t373 = -t457 * t429 + t505;
t336 = t452 * t373 + t453 * t497;
t339 = t453 * t373 - t452 * t497;
t294 = t458 * t336 + t455 * t339;
t296 = t455 * t336 - t458 * t339;
t566 = t459 * t294 - t456 * t296;
t565 = t456 * t294 + t459 * t296;
t407 = t453 * t445 + t452 * t446;
t410 = t452 * t445 - t453 * t446;
t348 = t458 * t407 - t455 * t410;
t451 = g(3) - qJDD(3);
t386 = qJ(3) * t407 - t453 * t451;
t538 = qJ(3) * t410 - t452 * t451;
t290 = pkin(6) * t348 + t458 * t386 - t455 * t538;
t352 = t455 * t407 + t458 * t410;
t306 = t456 * t348 + t459 * t352;
t552 = pkin(6) * t352 + t455 * t386 + t458 * t538;
t564 = pkin(5) * t306 + t456 * t290 + t459 * t552;
t537 = t459 * t348 - t456 * t352;
t563 = pkin(5) * t537 + t459 * t290 - t456 * t552;
t562 = qJ(3) * t326;
t558 = -pkin(2) * t369 + qJ(3) * t329;
t434 = t459 * g(1) + t456 * g(2);
t528 = qJD(1) ^ 2;
t465 = -t528 * pkin(1) - t434;
t433 = t456 * g(1) - t459 * g(2);
t469 = qJDD(1) * pkin(1) + t433;
t363 = t455 * t469 + t458 * t465;
t356 = -t445 * pkin(2) + t363;
t462 = -t455 * t465 + t458 * t469;
t461 = t446 * pkin(2) + t462;
t310 = t452 * t356 - t453 * t461;
t311 = t453 * t356 + t452 * t461;
t483 = t452 * t310 + t453 * t311;
t252 = t453 * t310 - t452 * t311;
t496 = t458 * t252;
t214 = -t455 * t483 + t496;
t501 = t455 * t252;
t546 = t458 * t483 + t501;
t195 = t456 * t214 + t459 * t546;
t557 = t459 * t214 - t456 * t546;
t416 = t458 * t445 + t455 * t446;
t419 = t455 * t445 - t458 * t446;
t359 = t456 * t416 + t459 * t419;
t392 = pkin(6) * t416 - t458 * g(3);
t539 = pkin(6) * t419 - t455 * g(3);
t554 = pkin(5) * t359 + t456 * t392 + t459 * t539;
t470 = t459 * t416 - t456 * t419;
t553 = pkin(5) * t470 + t459 * t392 - t456 * t539;
t482 = t458 * t363 - t455 * t462;
t315 = -t455 * t363 - t458 * t462;
t500 = t456 * t315;
t259 = t459 * t482 + t500;
t495 = t459 * t315;
t545 = -t456 * t482 + t495;
t544 = 2 * qJD(5);
t541 = pkin(3) * t369;
t540 = pkin(7) * t369;
t527 = pkin(4) * t457;
t474 = -qJ(5) * t454 - t527;
t402 = t474 * t447;
t300 = -t445 * pkin(3) + t446 * pkin(7) + t311;
t493 = -t457 * t300 + t454 * t451;
t472 = t457 * t447 * t402 + qJDD(4) * qJ(5) + (qJD(4) * t544) - t493;
t530 = t454 * t429 + t498;
t517 = t447 * t454;
t529 = t402 * t517 + qJDD(5);
t430 = -t460 - t515;
t423 = qJDD(4) + t432;
t506 = t454 * t423;
t372 = t457 * t430 - t506;
t325 = t452 * t372 + t453 * t404;
t328 = t453 * t372 - t452 * t404;
t274 = t458 * t325 + t455 * t328;
t277 = -t455 * t325 + t458 * t328;
t226 = t459 * t274 + t456 * t277;
t526 = pkin(5) * t226;
t492 = t449 + t450;
t412 = t492 * t446;
t420 = t492 * t445;
t354 = t452 * t412 + t453 * t420;
t355 = t453 * t412 - t452 * t420;
t308 = t458 * t354 + t455 * t355;
t309 = -t455 * t354 + t458 * t355;
t251 = t459 * t308 + t456 * t309;
t525 = pkin(5) * t251;
t524 = pkin(6) * t274;
t523 = pkin(6) * t308;
t413 = t457 * t423;
t367 = t454 * t430 + t413;
t522 = pkin(7) * t367;
t519 = qJ(3) * t325;
t518 = qJ(3) * t354;
t299 = -t446 * pkin(3) - t445 * pkin(7) + t310;
t508 = t454 * t299;
t507 = t454 * t404;
t499 = t457 * t299;
t281 = t454 * t300 + t457 * t451;
t494 = t420 - t460;
t485 = -pkin(1) * t367 + pkin(6) * t277;
t484 = -pkin(2) * t367 + qJ(3) * t328;
t236 = t454 * t281 - t457 * t493;
t382 = -t456 * t433 - t459 * t434;
t480 = t452 * t432;
t479 = t453 * t432;
t478 = pkin(1) * t274 + pkin(2) * t325 + pkin(3) * t404 + pkin(7) * t372;
t477 = pkin(1) * t308 + pkin(2) * t354 + pkin(3) * t420 + pkin(7) * t412;
t262 = -pkin(3) * t367 + t281;
t426 = t459 * qJDD(1) - t456 * t528;
t475 = -pkin(5) * t426 - t456 * g(3);
t473 = pkin(4) * t454 - qJ(5) * t457;
t235 = t457 * t281 + t454 * t493;
t471 = t457 * t403 + t507;
t381 = t459 * t433 - t456 * t434;
t468 = t486 + t502;
t467 = -t487 + t497;
t466 = -qJDD(4) * pkin(4) + t281 + t529;
t464 = -t467 * pkin(4) + t299 + (-t468 - t486) * qJ(5);
t463 = t517 * t544 - t464;
t428 = t460 - t516;
t425 = t456 * qJDD(1) + t459 * t528;
t400 = -pkin(5) * t425 + t459 * g(3);
t399 = t473 * t446;
t398 = t492 * t491;
t380 = t452 * qJDD(4) + t453 * t398;
t379 = -t453 * qJDD(4) + t452 * t398;
t378 = -t449 * t491 + t457 * t468;
t377 = -t450 * t491 - t454 * t467;
t374 = -t454 * t428 + t413;
t368 = t457 * t428 + t506;
t346 = qJ(3) * t355;
t340 = t453 * t374 + t452 * t502;
t337 = t452 * t374 - t453 * t502;
t334 = t453 * t378 - t480;
t333 = t453 * t377 + t480;
t332 = t452 * t378 + t479;
t331 = t452 * t377 - t479;
t322 = -t455 * t379 + t458 * t380;
t321 = t458 * t379 + t455 * t380;
t312 = pkin(1) * g(3) + pkin(6) * t482;
t302 = pkin(6) * t309;
t298 = -t455 * t337 + t458 * t340;
t295 = t458 * t337 + t455 * t340;
t286 = -t455 * t332 + t458 * t334;
t285 = -t455 * t331 + t458 * t333;
t284 = t458 * t332 + t455 * t334;
t283 = t458 * t331 + t455 * t333;
t271 = t499 + t540;
t270 = t508 - t522;
t269 = -t456 * t321 + t459 * t322;
t268 = t459 * t321 + t456 * t322;
t263 = -t493 + t541;
t261 = t460 * qJ(5) - t466;
t260 = -t460 * pkin(4) + t472;
t257 = t494 * qJ(5) + t466;
t256 = t494 * pkin(4) + t472;
t255 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t517 + t464;
t254 = -t456 * t308 + t459 * t309;
t248 = (t404 - t487) * pkin(4) + t463;
t247 = -pkin(4) * t487 + qJ(5) * t403 + t463;
t246 = pkin(5) * t254;
t245 = pkin(2) * t451 + qJ(3) * t483;
t244 = (-t430 - t460) * qJ(5) + (-qJDD(4) - t423) * pkin(4) + t262 + t529;
t243 = -t541 - qJ(5) * t424 + (-t427 + t460) * pkin(4) - t472;
t242 = -t456 * t295 + t459 * t298;
t241 = t459 * t295 + t456 * t298;
t240 = -t456 * t284 + t459 * t286;
t239 = -t456 * t283 + t459 * t285;
t238 = t459 * t284 + t456 * t286;
t237 = t459 * t283 + t456 * t285;
t233 = -pkin(4) * t365 + t457 * t247 - t540;
t232 = qJ(5) * t364 - t454 * t248 - t522;
t229 = -t456 * t274 + t459 * t277;
t225 = pkin(5) * t229;
t224 = t453 * t235 - t518;
t223 = t452 * t235 + t346;
t222 = t457 * t260 - t454 * t261;
t221 = t454 * t260 + t457 * t261;
t220 = -t454 * t256 + t457 * t257;
t219 = t453 * t236 + t452 * t299;
t218 = t452 * t236 - t453 * t299;
t217 = -t452 * t263 + t453 * t271 + t562;
t216 = -t452 * t262 + t453 * t270 - t519;
t211 = t453 * t263 + t452 * t271 - t558;
t210 = t453 * t262 + t452 * t270 + t484;
t209 = t453 * t220 - t452 * t399 - t518;
t208 = t452 * t220 + t453 * t399 + t346;
t207 = t453 * t222 + t452 * t255;
t206 = t452 * t222 - t453 * t255;
t205 = t453 * t232 - t452 * t244 - t519;
t204 = t453 * t233 - t452 * t243 - t562;
t203 = t452 * t232 + t453 * t244 + t484;
t202 = t452 * t233 + t453 * t243 + t558;
t201 = -pkin(3) * t221 - pkin(4) * t261 - qJ(5) * t260;
t200 = -pkin(7) * t221 + t255 * t473;
t199 = -t455 * t223 + t458 * t224 - t523;
t198 = t458 * t223 + t455 * t224 + t302;
t197 = -t455 * t218 + t458 * t219;
t196 = t458 * t218 + t455 * t219;
t193 = pkin(6) * t214 + qJ(3) * t496 - t455 * t245;
t192 = pkin(1) * t451 + pkin(6) * t546 + qJ(3) * t501 + t458 * t245;
t191 = -qJ(3) * t218 - (pkin(3) * t452 - pkin(7) * t453) * t235;
t190 = -t455 * t211 + t458 * t217 + t571;
t189 = -t455 * t210 + t458 * t216 - t524;
t188 = t458 * t211 + t455 * t217 - t569;
t187 = t458 * t210 + t455 * t216 + t485;
t186 = -t455 * t208 + t458 * t209 - t523;
t185 = t458 * t208 + t455 * t209 + t302;
t184 = -t455 * t206 + t458 * t207;
t183 = t458 * t206 + t455 * t207;
t182 = qJ(3) * t219 - (-pkin(3) * t453 - pkin(7) * t452 - pkin(2)) * t235;
t181 = -t455 * t203 + t458 * t205 - t524;
t180 = -t455 * t202 + t458 * t204 - t571;
t179 = t458 * t203 + t455 * t205 + t485;
t178 = t458 * t202 + t455 * t204 + t569;
t177 = -t456 * t196 + t459 * t197;
t176 = t459 * t196 + t456 * t197;
t175 = -qJ(3) * t206 + t453 * t200 - t452 * t201;
t174 = -t456 * t183 + t459 * t184;
t173 = t459 * t183 + t456 * t184;
t172 = -pkin(2) * t221 + qJ(3) * t207 + t452 * t200 + t453 * t201;
t171 = -pkin(6) * t196 - t455 * t182 + t458 * t191;
t170 = pkin(1) * t235 + pkin(6) * t197 + t458 * t182 + t455 * t191;
t169 = -pkin(6) * t183 - t455 * t172 + t458 * t175;
t168 = -pkin(1) * t221 + pkin(6) * t184 + t458 * t172 + t455 * t175;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t425, -t426, 0, t382, 0, 0, 0, 0, 0, 0, -t470, t359, 0, t259, 0, 0, 0, 0, 0, 0, -t537, t306, 0, t195, 0, 0, 0, 0, 0, 0, t229, t231, t254, t177, 0, 0, 0, 0, 0, 0, t229, t254, -t231, t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t426, -t425, 0, t381, 0, 0, 0, 0, 0, 0, -t359, -t470, 0, -t545, 0, 0, 0, 0, 0, 0, -t306, -t537, 0, -t557, 0, 0, 0, 0, 0, 0, t226, -t227, t251, t176, 0, 0, 0, 0, 0, 0, t226, t251, t227, t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451, 0, 0, 0, 0, 0, 0, t367, -t369, 0, -t235, 0, 0, 0, 0, 0, 0, t367, 0, t369, t221; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t426, 0, -t425, 0, t475, -t400, -t381, -pkin(5) * t381, 0, 0, -t359, 0, -t470, 0, t554, t553, t545, pkin(5) * t545 + pkin(6) * t495 - t456 * t312, 0, 0, -t306, 0, -t537, 0, t564, t563, t557, pkin(5) * t557 - t456 * t192 + t459 * t193, t240, t567, t242, t239, t565, t269, -t456 * t187 + t459 * t189 - t526, -t456 * t188 + t459 * t190 + t573, -t456 * t198 + t459 * t199 - t525, -pkin(5) * t176 - t456 * t170 + t459 * t171, t240, t242, -t567, t269, -t565, t239, -t456 * t179 + t459 * t181 - t526, -t456 * t185 + t459 * t186 - t525, -t456 * t178 + t459 * t180 - t573, -pkin(5) * t173 - t456 * t168 + t459 * t169; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t425, 0, t426, 0, t400, t475, t382, pkin(5) * t382, 0, 0, t470, 0, -t359, 0, -t553, t554, t259, pkin(5) * t259 + pkin(6) * t500 + t459 * t312, 0, 0, t537, 0, -t306, 0, -t563, t564, t195, pkin(5) * t195 + t459 * t192 + t456 * t193, t238, -t568, t241, t237, -t566, t268, t459 * t187 + t456 * t189 + t225, t459 * t188 + t456 * t190 + t572, t459 * t198 + t456 * t199 + t246, pkin(5) * t177 + t459 * t170 + t456 * t171, t238, t241, t568, t268, t566, t237, t459 * t179 + t456 * t181 + t225, t459 * t185 + t456 * t186 + t246, t459 * t178 + t456 * t180 - t572, pkin(5) * t174 + t459 * t168 + t456 * t169; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t433, t434, 0, 0, 0, 0, 0, 0, 0, t446, -pkin(1) * t419 + t462, -pkin(1) * t416 - t363, 0, -pkin(1) * t315, 0, 0, 0, 0, 0, t446, -pkin(1) * t352 - pkin(2) * t410 - t310, -pkin(1) * t348 - pkin(2) * t407 - t311, 0, -pkin(1) * t214 - pkin(2) * t252, t365, t471, t368, t364, t530, 0, t478 - t499, -pkin(3) * t403 + t508 - t570, t236 + t477, pkin(1) * t196 + pkin(2) * t218 - pkin(3) * t299 + pkin(7) * t236, t365, t368, -t471, 0, -t530, t364, qJ(5) * t507 + t457 * t248 + t478, t457 * t256 + t454 * t257 + t477, t454 * t247 + (pkin(3) + t527) * t403 + t570, pkin(1) * t183 + pkin(2) * t206 + pkin(7) * t222 + (-pkin(3) + t474) * t255;];
tauB_reg = t1;
