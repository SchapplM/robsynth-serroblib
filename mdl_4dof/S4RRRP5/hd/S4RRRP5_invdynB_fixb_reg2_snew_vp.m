% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRP5
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% tauB_reg [6x(5*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRP5_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP5_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP5_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:17:10
% EndTime: 2019-12-31 17:17:16
% DurationCPUTime: 5.98s
% Computational Cost: add. (8860->356), mult. (19053->474), div. (0->0), fcn. (12378->6), ass. (0->269)
t463 = sin(qJ(3));
t464 = sin(qJ(2));
t466 = cos(qJ(3));
t467 = cos(qJ(2));
t422 = (t463 * t467 + t464 * t466) * qJD(1);
t419 = t422 ^ 2;
t460 = qJD(2) + qJD(3);
t536 = t460 ^ 2;
t375 = t536 + t419;
t503 = qJD(1) * t464;
t420 = -t466 * t467 * qJD(1) + t463 * t503;
t385 = t422 * t420;
t459 = qJDD(2) + qJDD(3);
t545 = t385 + t459;
t523 = t463 * t545;
t310 = t466 * t375 + t523;
t510 = t466 * t545;
t334 = t463 * t375 - t510;
t297 = t464 * t310 + t467 * t334;
t465 = sin(qJ(1));
t468 = cos(qJ(1));
t500 = qJD(1) * qJD(2);
t492 = t467 * t500;
t499 = t464 * qJDD(1);
t431 = t492 + t499;
t454 = t467 * qJDD(1);
t493 = t464 * t500;
t432 = t454 - t493;
t473 = t420 * qJD(3) - t466 * t431 - t463 * t432;
t532 = t460 * t420;
t553 = -t473 - t532;
t268 = t465 * t297 - t468 * t553;
t611 = pkin(4) * t268;
t270 = t468 * t297 + t465 * t553;
t610 = pkin(4) * t270;
t277 = t467 * t310 - t464 * t334;
t609 = pkin(5) * t277;
t608 = -pkin(1) * t277 - pkin(2) * t310;
t607 = pkin(1) * t553 - pkin(5) * t297;
t489 = t463 * t431 - t466 * t432;
t324 = (qJD(3) + t460) * t422 + t489;
t289 = -t463 * t324 + t466 * t553;
t526 = t463 * t553;
t291 = t466 * t324 + t526;
t253 = t464 * t289 + t467 * t291;
t537 = t420 ^ 2;
t382 = t419 - t537;
t606 = t465 * t253 + t468 * t382;
t605 = t468 * t253 - t465 * t382;
t406 = t537 - t536;
t338 = t463 * t406 + t510;
t342 = t466 * t406 - t523;
t302 = t464 * t338 - t467 * t342;
t325 = (qJD(3) - t460) * t422 + t489;
t604 = t465 * t302 - t468 * t325;
t603 = t468 * t302 + t465 * t325;
t360 = -t537 - t419;
t543 = -t532 + t473;
t562 = -t466 * t325 - t463 * t543;
t563 = -t463 * t325 + t466 * t543;
t579 = -t464 * t563 + t467 * t562;
t591 = t465 * t360 + t468 * t579;
t602 = pkin(4) * t591;
t593 = -t468 * t360 + t465 * t579;
t601 = pkin(4) * t593;
t578 = t464 * t562 + t467 * t563;
t599 = pkin(5) * t578;
t598 = pkin(6) * t310;
t597 = pkin(6) * t334;
t238 = -pkin(1) * t578 - pkin(2) * t563;
t596 = -pkin(1) * t360 + pkin(5) * t579;
t595 = -t467 * t289 + t464 * t291;
t407 = -t419 + t536;
t546 = -t385 + t459;
t522 = t463 * t546;
t564 = t466 * t407 + t522;
t373 = t466 * t546;
t565 = -t463 * t407 + t373;
t576 = -t464 * t564 + t467 * t565;
t594 = t465 * t576 + t468 * t543;
t592 = -t465 * t543 + t468 * t576;
t590 = t467 * t338 + t464 * t342;
t542 = -t536 - t537;
t549 = t466 * t542 - t522;
t552 = t463 * t542 + t373;
t560 = t464 * t549 + t467 * t552;
t588 = pkin(5) * t560;
t561 = -t464 * t552 + t467 * t549;
t587 = pkin(5) * t561;
t586 = pkin(6) * t563;
t584 = t465 * t561;
t582 = t468 * t561;
t581 = -pkin(1) * t560 - pkin(2) * t552;
t580 = -pkin(2) * t360 + pkin(6) * t562;
t577 = t464 * t565 + t467 * t564;
t572 = pkin(6) * t549;
t571 = pkin(6) * t552;
t566 = t553 * qJ(4);
t559 = 2 * qJD(4);
t474 = (-t420 * t463 - t422 * t466) * t460;
t530 = t460 * t463;
t403 = t422 * t530;
t529 = t460 * t466;
t495 = t420 * t529;
t480 = t403 - t495;
t540 = -t464 * t474 + t467 * t480;
t551 = -t468 * t459 + t465 * t540;
t494 = t468 * t385;
t358 = -t422 * qJD(3) - t489;
t476 = -t463 * t358 + t495;
t481 = t466 * t358 + t420 * t530;
t539 = -t464 * t481 + t467 * t476;
t550 = t465 * t539 + t494;
t496 = t465 * t385;
t548 = t468 * t539 - t496;
t547 = t465 * t459 + t468 * t540;
t470 = qJD(1) ^ 2;
t515 = t464 * t470;
t442 = t468 * g(1) + t465 * g(2);
t424 = -t470 * pkin(1) + qJDD(1) * pkin(5) - t442;
t518 = t464 * t424;
t348 = qJDD(2) * pkin(2) - t431 * pkin(6) - t518 + (pkin(2) * t515 + pkin(6) * t500 - g(3)) * t467;
t402 = -t464 * g(3) + t467 * t424;
t462 = t467 ^ 2;
t456 = t462 * t470;
t478 = qJD(2) * pkin(2) - pkin(6) * t503;
t352 = -pkin(2) * t456 + t432 * pkin(6) - qJD(2) * t478 + t402;
t307 = t463 * t348 + t466 * t352;
t381 = t420 * pkin(3) - t422 * qJ(4);
t479 = t459 * qJ(4) - t420 * t381 + t460 * t559 + t307;
t541 = t464 * t480 + t467 * t474;
t538 = t464 * t476 + t467 * t481;
t535 = pkin(3) * t466;
t534 = t358 * pkin(3);
t533 = qJ(4) * t466;
t531 = t460 * t422;
t461 = t464 ^ 2;
t528 = t461 * t470;
t441 = t465 * g(1) - t468 * g(2);
t477 = qJDD(1) * pkin(1) + t441;
t362 = t432 * pkin(2) - t478 * t503 + (pkin(6) * t462 + pkin(5)) * t470 + t477;
t524 = t463 * t362;
t306 = -t466 * t348 + t463 * t352;
t264 = -t466 * t306 + t463 * t307;
t520 = t464 * t264;
t423 = t470 * pkin(5) + t477;
t519 = t464 * t423;
t449 = t467 * t515;
t439 = qJDD(2) + t449;
t517 = t464 * t439;
t440 = qJDD(2) - t449;
t516 = t464 * t440;
t511 = t466 * t362;
t509 = t467 * t264;
t508 = t467 * t423;
t507 = t467 * t440;
t505 = -t360 - t536;
t504 = t461 + t462;
t498 = t465 * qJDD(1);
t497 = t468 * qJDD(1);
t491 = -qJ(4) * t463 - pkin(2);
t265 = t463 * t306 + t466 * t307;
t401 = t467 * g(3) + t518;
t351 = t464 * t401 + t467 * t402;
t391 = -t465 * t441 - t468 * t442;
t488 = t465 * t449;
t487 = t468 * t449;
t485 = t422 * t381 + qJDD(4) + t306;
t436 = -t465 * t470 + t497;
t484 = -pkin(4) * t436 - t465 * g(3);
t320 = t422 * t529 - t463 * t473;
t321 = -t466 * t473 - t403;
t286 = -t464 * t320 + t467 * t321;
t483 = t465 * t286 - t494;
t482 = t468 * t286 + t496;
t350 = t467 * t401 - t464 * t402;
t390 = t468 * t441 - t465 * t442;
t475 = -t459 * pkin(3) + t485;
t472 = -pkin(3) * t531 + t422 * t559 + t362;
t471 = t472 + t566;
t469 = qJD(2) ^ 2;
t447 = -t456 - t469;
t446 = t456 - t469;
t445 = -t469 - t528;
t444 = t469 - t528;
t438 = t456 - t528;
t437 = t456 + t528;
t435 = t468 * t470 + t498;
t434 = t504 * qJDD(1);
t433 = t454 - 0.2e1 * t493;
t430 = 0.2e1 * t492 + t499;
t428 = t467 * t439;
t427 = t504 * t500;
t416 = -pkin(4) * t435 + t468 * g(3);
t405 = t467 * t431 - t461 * t500;
t404 = -t464 * t432 - t462 * t500;
t397 = -t464 * t445 - t507;
t396 = -t464 * t444 + t428;
t395 = t467 * t447 - t517;
t394 = t467 * t446 - t516;
t393 = t467 * t445 - t516;
t392 = t464 * t447 + t428;
t388 = t468 * t434 - t465 * t437;
t387 = t465 * t434 + t468 * t437;
t386 = -t464 * t430 + t467 * t433;
t372 = t468 * t397 + t465 * t430;
t371 = t468 * t395 - t465 * t433;
t370 = t465 * t397 - t468 * t430;
t369 = t465 * t395 + t468 * t433;
t364 = -pkin(5) * t393 - t508;
t363 = -pkin(5) * t392 - t519;
t357 = -pkin(1) * t393 + t402;
t356 = -pkin(1) * t392 + t401;
t323 = -t358 + t531;
t315 = t468 * t351 - t465 * t423;
t314 = t465 * t351 + t468 * t423;
t304 = -t511 + t598;
t303 = -t524 - t571;
t283 = t467 * t320 + t464 * t321;
t276 = qJ(4) * t536 - t475;
t275 = -pkin(3) * t536 + t479;
t274 = -pkin(2) * t553 - t524 + t597;
t273 = -pkin(2) * t324 + t511 + t572;
t272 = t471 + t534;
t271 = t465 * t323 + t582;
t269 = -t468 * t323 + t584;
t267 = t505 * qJ(4) + t475;
t266 = t505 * pkin(3) + t479;
t263 = t465 * t324 + t582;
t261 = -t468 * t324 + t584;
t259 = (-t323 + t358) * pkin(3) + t471;
t258 = t472 + t534 + 0.2e1 * t566;
t257 = pkin(2) * t362 + pkin(6) * t265;
t250 = t307 - t608;
t245 = t306 + t581;
t244 = t466 * t275 - t463 * t276;
t243 = t463 * t275 + t466 * t276;
t242 = -t264 - t586;
t241 = -t463 * t259 - t323 * t533 - t571;
t240 = -pkin(3) * t526 + t466 * t258 - t598;
t239 = t265 + t580;
t237 = (-t542 - t536) * qJ(4) + (-t546 - t459) * pkin(3) + t485 + t581;
t236 = -t464 * t274 + t467 * t304 + t609;
t235 = t467 * t265 - t520;
t234 = t464 * t265 + t509;
t233 = t466 * t259 + t491 * t323 + t572;
t232 = -t597 + t463 * t258 + (pkin(2) + t535) * t553;
t231 = -qJ(4) * t545 + (-t375 + t536) * pkin(3) - t479 + t608;
t230 = -t464 * t273 + t467 * t303 - t588;
t229 = t468 * t235 - t465 * t362;
t228 = t465 * t235 + t468 * t362;
t227 = -pkin(3) * t543 + qJ(4) * t325 + t238;
t226 = -t463 * t266 + t466 * t267 - t586;
t225 = t466 * t266 + t463 * t267 + t580;
t224 = -pkin(1) * t234 - pkin(2) * t264;
t223 = -t464 * t243 + t467 * t244;
t222 = t467 * t243 + t464 * t244;
t221 = -pkin(6) * t243 + (-pkin(3) * t463 + t533) * t272;
t220 = t468 * t223 - t465 * t272;
t219 = t465 * t223 + t468 * t272;
t218 = pkin(6) * t244 + (-t491 + t535) * t272;
t217 = -pkin(5) * t234 - pkin(6) * t509 - t464 * t257;
t216 = -t464 * t233 + t467 * t241 - t588;
t215 = -t464 * t232 + t467 * t240 - t609;
t214 = -t464 * t239 + t467 * t242 - t599;
t213 = -pkin(1) * t222 - pkin(2) * t243 - pkin(3) * t276 - qJ(4) * t275;
t212 = -t464 * t225 + t467 * t226 - t599;
t211 = -pkin(5) * t222 - t464 * t218 + t467 * t221;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t435, -t436, 0, t391, 0, 0, 0, 0, 0, 0, t371, t372, t388, t315, 0, 0, 0, 0, 0, 0, t263, t270, t591, t229, 0, 0, 0, 0, 0, 0, t271, t591, -t270, t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t436, -t435, 0, t390, 0, 0, 0, 0, 0, 0, t369, t370, t387, t314, 0, 0, 0, 0, 0, 0, t261, t268, t593, t228, 0, 0, 0, 0, 0, 0, t269, t593, -t268, t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t392, t393, 0, -t350, 0, 0, 0, 0, 0, 0, t560, -t277, t578, t234, 0, 0, 0, 0, 0, 0, t560, t578, t277, t222; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t436, 0, -t435, 0, t484, -t416, -t390, -pkin(4) * t390, t468 * t405 - t488, t468 * t386 - t465 * t438, t468 * t396 + t464 * t498, t468 * t404 + t488, t468 * t394 + t465 * t454, t465 * qJDD(2) + t468 * t427, -pkin(4) * t369 - t465 * t356 + t468 * t363, -pkin(4) * t370 - t465 * t357 + t468 * t364, -pkin(4) * t387 + t468 * t350, -pkin(4) * t314 - (pkin(1) * t465 - pkin(5) * t468) * t350, t482, -t605, t592, t548, -t603, t547, -pkin(4) * t261 + t468 * t230 - t465 * t245, t468 * t236 - t465 * t250 - t611, t468 * t214 - t465 * t238 - t601, -pkin(4) * t228 + t468 * t217 - t465 * t224, t482, t592, t605, t547, t603, t548, -pkin(4) * t269 + t468 * t216 - t465 * t237, t468 * t212 - t465 * t227 - t601, t468 * t215 - t465 * t231 + t611, -pkin(4) * t219 + t468 * t211 - t465 * t213; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t435, 0, t436, 0, t416, t484, t391, pkin(4) * t391, t465 * t405 + t487, t465 * t386 + t468 * t438, t465 * t396 - t464 * t497, t465 * t404 - t487, t465 * t394 - t467 * t497, -t468 * qJDD(2) + t465 * t427, pkin(4) * t371 + t468 * t356 + t465 * t363, pkin(4) * t372 + t468 * t357 + t465 * t364, pkin(4) * t388 + t465 * t350, pkin(4) * t315 - (-pkin(1) * t468 - pkin(5) * t465) * t350, t483, -t606, t594, t550, -t604, t551, pkin(4) * t263 + t465 * t230 + t468 * t245, t465 * t236 + t468 * t250 + t610, t465 * t214 + t468 * t238 + t602, pkin(4) * t229 + t465 * t217 + t468 * t224, t483, t594, t606, t551, t604, t550, pkin(4) * t271 + t465 * t216 + t468 * t237, t465 * t212 + t468 * t227 + t602, t465 * t215 + t468 * t231 - t610, pkin(4) * t220 + t465 * t211 + t468 * t213; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t441, t442, 0, 0, (t431 + t492) * t464, t467 * t430 + t464 * t433, t467 * t444 + t517, (t432 - t493) * t467, t464 * t446 + t507, 0, pkin(1) * t433 + pkin(5) * t395 + t508, -pkin(1) * t430 + pkin(5) * t397 - t519, pkin(1) * t437 + pkin(5) * t434 + t351, pkin(1) * t423 + pkin(5) * t351, t283, -t595, t577, t538, t590, t541, -pkin(1) * t324 + t467 * t273 + t464 * t303 + t587, t467 * t274 + t464 * t304 - t607, t467 * t239 + t464 * t242 + t596, pkin(1) * t362 + pkin(5) * t235 - pkin(6) * t520 + t467 * t257, t283, t577, t595, t541, -t590, t538, -pkin(1) * t323 + t467 * t233 + t464 * t241 + t587, t467 * t225 + t464 * t226 + t596, t467 * t232 + t464 * t240 + t607, pkin(1) * t272 + pkin(5) * t223 + t467 * t218 + t464 * t221;];
tauB_reg = t1;
