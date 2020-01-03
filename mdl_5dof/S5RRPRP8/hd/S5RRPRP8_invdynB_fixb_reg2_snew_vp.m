% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRP8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynB_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:32
% EndTime: 2019-12-31 20:04:40
% DurationCPUTime: 5.81s
% Computational Cost: add. (13415->431), mult. (29091->556), div. (0->0), fcn. (17213->6), ass. (0->327)
t526 = qJD(2) ^ 2;
t521 = sin(qJ(2));
t518 = t521 ^ 2;
t527 = qJD(1) ^ 2;
t599 = t518 * t527;
t496 = t526 + t599;
t524 = cos(qJ(2));
t501 = t524 * t527 * t521;
t492 = qJDD(2) - t501;
t578 = t524 * t492;
t440 = -t521 * t496 + t578;
t564 = qJD(1) * qJD(2);
t555 = t524 * t564;
t562 = t521 * qJDD(1);
t481 = 0.2e1 * t555 + t562;
t522 = sin(qJ(1));
t525 = cos(qJ(1));
t404 = t522 * t440 + t525 * t481;
t639 = pkin(5) * t404;
t407 = t525 * t440 - t522 * t481;
t638 = pkin(5) * t407;
t514 = -qJDD(2) + qJDD(4);
t520 = sin(qJ(4));
t523 = cos(qJ(4));
t465 = (-t520 * t521 - t523 * t524) * qJD(1);
t567 = qJD(1) * t524;
t568 = qJD(1) * t521;
t467 = -t520 * t567 + t523 * t568;
t603 = t467 * t465;
t541 = t514 + t603;
t637 = pkin(4) * t541;
t636 = pkin(6) * t440;
t620 = 2 * qJD(3);
t619 = pkin(2) + pkin(3);
t589 = t521 * t492;
t434 = t524 * t496 + t589;
t635 = pkin(1) * t434;
t634 = pkin(6) * t434;
t482 = t555 + t562;
t508 = t524 * qJDD(1);
t557 = t521 * t564;
t483 = t508 - t557;
t396 = t465 * qJD(4) + t523 * t482 - t520 * t483;
t515 = qJD(2) - qJD(4);
t605 = t465 * t515;
t633 = t396 - t605;
t484 = t508 - 0.2e1 * t557;
t579 = t524 * t484;
t592 = t521 * t481;
t422 = -t579 + t592;
t519 = t524 ^ 2;
t490 = (t518 - t519) * t527;
t632 = t522 * t422 + t525 * t490;
t631 = t525 * t422 - t522 * t490;
t598 = t519 * t527;
t498 = -t526 + t598;
t438 = -t524 * t498 + t589;
t560 = t525 * qJDD(1);
t630 = t522 * t438 + t524 * t560;
t629 = t525 * t438 - t522 * t508;
t594 = t520 * t541;
t581 = t523 * t541;
t543 = t482 + t555;
t628 = t543 * qJ(3);
t565 = qJD(5) * t467;
t627 = -t396 * qJ(5) - 0.2e1 * t565;
t626 = t483 * pkin(3) - pkin(7) * t598;
t618 = pkin(2) * t524;
t548 = -qJ(3) * t521 - t618;
t479 = t548 * qJD(1);
t495 = t525 * g(1) + t522 * g(2);
t469 = -t527 * pkin(1) + qJDD(1) * pkin(6) - t495;
t553 = t521 * g(3) - t524 * t469;
t536 = qJDD(2) * qJ(3) + qJD(2) * t620 + t479 * t567 - t553;
t625 = t521 * t498 + t578;
t624 = t479 * t568 + qJDD(3);
t552 = t520 * t482 + t523 * t483;
t395 = -t467 * qJD(4) - t552;
t445 = -t515 * pkin(4) - t467 * qJ(5);
t623 = t395 * qJ(5) + 0.2e1 * qJD(5) * t465 + t515 * t445;
t622 = qJ(5) * t605 - t637;
t463 = t465 ^ 2;
t621 = -t395 * pkin(4) - t463 * qJ(5) + t467 * t445 + qJDD(5);
t464 = t467 ^ 2;
t513 = t515 ^ 2;
t375 = t396 + t605;
t532 = (-qJD(4) - t515) * t467 - t552;
t329 = -t523 * t375 + t520 * t532;
t331 = t520 * t375 + t523 * t532;
t293 = t521 * t329 + t524 * t331;
t397 = -t463 - t464;
t279 = t522 * t293 + t525 * t397;
t617 = pkin(5) * t279;
t409 = -t513 - t463;
t357 = t520 * t409 + t581;
t358 = t523 * t409 - t594;
t321 = t521 * t357 + t524 * t358;
t370 = (qJD(4) - t515) * t467 + t552;
t298 = t522 * t321 + t525 * t370;
t616 = pkin(5) * t298;
t442 = -t464 - t513;
t416 = -t514 + t603;
t595 = t520 * t416;
t376 = t523 * t442 + t595;
t582 = t523 * t416;
t377 = -t520 * t442 + t582;
t333 = t521 * t376 + t524 * t377;
t303 = t522 * t333 + t525 * t633;
t615 = pkin(5) * t303;
t499 = -t526 - t598;
t491 = qJDD(2) + t501;
t590 = t521 * t491;
t437 = t524 * t499 - t590;
t403 = t522 * t437 + t525 * t484;
t614 = pkin(5) * t403;
t569 = t518 + t519;
t486 = t569 * qJDD(1);
t489 = t569 * t527;
t424 = t522 * t486 + t525 * t489;
t613 = pkin(5) * t424;
t291 = -t524 * t329 + t521 * t331;
t612 = pkin(6) * t291;
t320 = -t524 * t357 + t521 * t358;
t611 = pkin(6) * t320;
t332 = -t524 * t376 + t521 * t377;
t610 = pkin(6) * t332;
t474 = t524 * t491;
t432 = t521 * t499 + t474;
t609 = pkin(6) * t432;
t608 = t483 * pkin(2);
t607 = qJ(3) * t524;
t601 = t515 * t520;
t600 = t515 * t523;
t388 = -t526 * pkin(2) + t536;
t493 = -qJD(2) * pkin(3) - pkin(7) * t568;
t354 = -pkin(3) * t598 - t483 * pkin(7) + qJD(2) * t493 + t388;
t443 = t524 * g(3) + t521 * t469;
t533 = -qJDD(2) * pkin(2) + t443 + t624;
t394 = t526 * qJ(3) - t533;
t356 = (-t482 + t555) * pkin(7) - t491 * pkin(3) - t394;
t317 = t520 * t354 - t523 * t356;
t534 = t317 + t622;
t289 = -t534 + t627;
t597 = t520 * t289;
t494 = t522 * g(1) - t525 * g(2);
t468 = qJDD(1) * pkin(1) + t527 * pkin(6) + t494;
t531 = -pkin(2) * t557 + t468;
t349 = t608 + t482 * qJ(3) + (qJD(2) * t607 + (t620 + t493) * t521) * qJD(1) + t531 + t626;
t596 = t520 * t349;
t593 = t521 * t468;
t591 = t521 * t484;
t584 = t523 * t289;
t583 = t523 * t349;
t580 = t524 * t468;
t575 = pkin(1) * t397 + pkin(6) * t293;
t574 = pkin(1) * t370 + pkin(6) * t321;
t573 = pkin(1) * t633 + pkin(6) * t333;
t318 = t523 * t354 + t520 * t356;
t572 = pkin(1) * t484 + pkin(6) * t437;
t571 = pkin(1) * t489 + pkin(6) * t486;
t570 = t489 - t526;
t561 = t522 * qJDD(1);
t559 = t522 * t603;
t558 = t525 * t603;
t387 = t521 * t443 - t524 * t553;
t429 = -t522 * t494 - t525 * t495;
t551 = t522 * t501;
t550 = t525 * t501;
t392 = -pkin(1) * t432 + t443;
t488 = -t522 * t527 + t560;
t549 = -pkin(5) * t488 - t522 * g(3);
t547 = pkin(2) * t521 - t607;
t546 = -pkin(7) * t329 + qJ(3) * t397;
t545 = -pkin(7) * t357 + qJ(3) * t370;
t544 = -pkin(7) * t376 + qJ(3) * t633;
t281 = -t523 * t317 + t520 * t318;
t282 = t520 * t317 + t523 * t318;
t386 = t524 * t443 + t521 * t553;
t542 = t524 * t481 + t591;
t428 = t525 * t494 - t522 * t495;
t539 = -pkin(7) * t358 + t619 * t370;
t538 = -pkin(7) * t377 + t619 * t633;
t537 = -pkin(7) * t331 + t619 * t397;
t535 = t318 + t623;
t258 = -pkin(1) * t291 - qJ(3) * t331 + t619 * t329;
t271 = -pkin(1) * t332 - qJ(3) * t377 + t619 * t376 - t318;
t266 = -pkin(1) * t320 - qJ(3) * t358 + t619 * t357 - t317;
t530 = t568 * t620 + t531;
t529 = t530 + t608;
t528 = t530 + t628;
t310 = t349 + t621;
t497 = t526 - t599;
t487 = t525 * t527 + t561;
t477 = t547 * qJDD(1);
t473 = t569 * t564;
t462 = -pkin(5) * t487 + t525 * g(3);
t451 = -t464 + t513;
t450 = t463 - t513;
t449 = t522 * qJDD(2) + t525 * t473;
t448 = t524 * t482 - t518 * t564;
t447 = -t525 * qJDD(2) + t522 * t473;
t446 = -t521 * t483 - t519 * t564;
t439 = -t521 * t497 + t474;
t433 = t524 * t497 + t590;
t431 = t543 * t521;
t430 = (t483 - t557) * t524;
t425 = t525 * t486 - t522 * t489;
t423 = pkin(5) * t425;
t419 = t464 - t463;
t415 = t525 * t448 - t551;
t414 = t525 * t446 + t551;
t413 = t522 * t448 + t550;
t412 = t522 * t446 - t550;
t411 = t525 * t439 + t521 * t561;
t410 = t522 * t439 - t521 * t560;
t406 = t525 * t437 - t522 * t484;
t402 = pkin(5) * t406;
t401 = (-t465 * t523 - t467 * t520) * t515;
t400 = (t465 * t520 - t467 * t523) * t515;
t399 = -t580 + t634;
t398 = -t593 - t609;
t393 = -t553 + t635;
t384 = t570 * qJ(3) + t533;
t383 = t570 * pkin(2) + t536;
t382 = t528 + t608;
t381 = t523 * t450 + t595;
t380 = -t520 * t451 + t581;
t379 = -t520 * t450 + t582;
t378 = -t523 * t451 - t594;
t366 = t523 * t396 + t467 * t601;
t365 = -t520 * t396 + t467 * t600;
t364 = -t520 * t395 + t465 * t600;
t363 = -t523 * t395 - t465 * t601;
t362 = (t483 + t484) * pkin(2) + t528;
t361 = (t481 + t543) * qJ(3) + t529;
t360 = t525 * t387 - t522 * t468;
t359 = t522 * t387 + t525 * t468;
t350 = (-t499 - t526) * qJ(3) + (-qJDD(2) - t491) * pkin(2) + t392 + t624;
t348 = -t635 - qJ(3) * t492 + (-t496 + t526) * pkin(2) - t536;
t347 = -t521 * t400 + t524 * t401;
t346 = t524 * t400 + t521 * t401;
t345 = t524 * t388 - t521 * t394;
t344 = t521 * t388 + t524 * t394;
t343 = -pkin(4) * t633 + qJ(5) * t416;
t342 = t525 * t347 - t522 * t514;
t341 = t522 * t347 + t525 * t514;
t340 = -pkin(2) * t592 + t524 * t361 - t634;
t339 = qJ(3) * t579 - t521 * t362 - t609;
t338 = -t521 * t383 + t524 * t384;
t337 = -t521 * t379 + t524 * t381;
t336 = -t521 * t378 + t524 * t380;
t335 = t524 * t379 + t521 * t381;
t334 = t524 * t378 + t521 * t380;
t330 = -t523 * t370 - t520 * t633;
t328 = t520 * t370 - t523 * t633;
t325 = -t521 * t365 + t524 * t366;
t324 = -t521 * t363 + t524 * t364;
t323 = t524 * t365 + t521 * t366;
t322 = t524 * t363 + t521 * t364;
t316 = t525 * t345 - t522 * t382;
t315 = t522 * t345 + t525 * t382;
t314 = t525 * t325 + t559;
t313 = t525 * t324 - t559;
t312 = t522 * t325 - t558;
t311 = t522 * t324 + t558;
t309 = t525 * t337 - t522 * t532;
t308 = t525 * t336 - t522 * t375;
t307 = t522 * t337 + t525 * t532;
t306 = t522 * t336 + t525 * t375;
t305 = -pkin(1) * t344 - pkin(2) * t394 - qJ(3) * t388;
t304 = t525 * t333 - t522 * t633;
t302 = pkin(5) * t304;
t301 = -qJ(5) * t442 + t310;
t300 = t544 + t583;
t299 = t525 * t321 - t522 * t370;
t297 = pkin(5) * t299;
t296 = -pkin(6) * t344 - t547 * t382;
t295 = t545 + t596;
t294 = -t463 * pkin(4) + t535;
t292 = -t521 * t328 + t524 * t330;
t290 = t524 * t328 + t521 * t330;
t287 = t538 - t596;
t286 = -pkin(4) * t370 + qJ(5) * t409 - t493 * t568 - t529 - t621 - t626 - t628;
t285 = t539 + t583;
t284 = t525 * t292 - t522 * t419;
t283 = t522 * t292 + t525 * t419;
t280 = t525 * t293 - t522 * t397;
t278 = pkin(5) * t280;
t277 = 0.2e1 * t565 + (t375 + t396) * qJ(5) + t534;
t276 = qJ(5) * t532 + (-t397 - t463) * pkin(4) + t535;
t275 = -pkin(7) * t281 + qJ(3) * t349;
t274 = t523 * t301 - t520 * t343 + t544;
t273 = -pkin(4) * t310 + qJ(5) * t294;
t272 = -qJ(5) * t581 - t520 * t286 + t545;
t270 = -t281 + t546;
t269 = -pkin(7) * t282 + t619 * t349;
t268 = -t520 * t301 - t523 * t343 + t538;
t267 = qJ(5) * t594 - t523 * t286 + t539;
t265 = t523 * t294 - t597;
t264 = t520 * t294 + t584;
t263 = -t282 + t537;
t262 = (t442 + t463) * pkin(4) + t271 - t623;
t261 = t521 * t281 + t524 * t282;
t260 = -t524 * t281 + t521 * t282;
t259 = -t521 * t287 + t524 * t300 - t610;
t257 = t266 - t622 + t627 + t637;
t256 = -t521 * t285 + t524 * t295 - t611;
t255 = t525 * t261 - t522 * t349;
t254 = t522 * t261 + t525 * t349;
t253 = -pkin(4) * t375 + t258;
t252 = -t520 * t276 + t523 * t277 + t546;
t251 = -t523 * t276 - t520 * t277 + t537;
t250 = t521 * t264 + t524 * t265;
t249 = -t524 * t264 + t521 * t265;
t248 = -t521 * t268 + t524 * t274 - t610;
t247 = -t521 * t267 + t524 * t272 - t611;
t246 = t525 * t250 - t522 * t310;
t245 = t522 * t250 + t525 * t310;
t244 = -t521 * t263 + t524 * t270 - t612;
t243 = -pkin(7) * t264 + qJ(3) * t310 - qJ(5) * t584 - t520 * t273;
t242 = -pkin(6) * t260 - t521 * t269 + t524 * t275;
t241 = -pkin(1) * t260 - qJ(3) * t282 + t619 * t281;
t240 = -pkin(7) * t265 + qJ(5) * t597 - t523 * t273 + t619 * t310;
t239 = -t521 * t251 + t524 * t252 - t612;
t238 = -pkin(1) * t249 + pkin(4) * t289 - qJ(3) * t265 + t619 * t264;
t237 = -pkin(6) * t249 - t521 * t240 + t524 * t243;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t487, -t488, 0, t429, 0, 0, 0, 0, 0, 0, t406, -t407, t425, t360, 0, 0, 0, 0, 0, 0, t406, t425, t407, t316, 0, 0, 0, 0, 0, 0, t299, t304, t280, t255, 0, 0, 0, 0, 0, 0, t299, t304, t280, t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t488, -t487, 0, t428, 0, 0, 0, 0, 0, 0, t403, -t404, t424, t359, 0, 0, 0, 0, 0, 0, t403, t424, t404, t315, 0, 0, 0, 0, 0, 0, t298, t303, t279, t254, 0, 0, 0, 0, 0, 0, t298, t303, t279, t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t432, -t434, 0, -t386, 0, 0, 0, 0, 0, 0, t432, 0, t434, t344, 0, 0, 0, 0, 0, 0, t320, t332, t291, t260, 0, 0, 0, 0, 0, 0, t320, t332, t291, t249; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t488, 0, -t487, 0, t549, -t462, -t428, -pkin(5) * t428, t415, -t631, t411, t414, -t629, t449, -t522 * t392 + t525 * t398 - t614, -t522 * t393 + t525 * t399 + t639, t525 * t386 - t613, -pkin(5) * t359 - (pkin(1) * t522 - pkin(6) * t525) * t386, t415, t411, t631, t449, t629, t414, t525 * t339 - t522 * t350 - t614, t525 * t338 - t522 * t477 - t613, t525 * t340 - t522 * t348 - t639, -pkin(5) * t315 + t525 * t296 - t522 * t305, t314, t284, t308, t313, t309, t342, t525 * t256 - t522 * t266 - t616, t525 * t259 - t522 * t271 - t615, t525 * t244 - t522 * t258 - t617, -pkin(5) * t254 - t522 * t241 + t525 * t242, t314, t284, t308, t313, t309, t342, t525 * t247 - t522 * t257 - t616, t525 * t248 - t522 * t262 - t615, t525 * t239 - t522 * t253 - t617, -pkin(5) * t245 + t525 * t237 - t522 * t238; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t487, 0, t488, 0, t462, t549, t429, pkin(5) * t429, t413, -t632, t410, t412, -t630, t447, t525 * t392 + t522 * t398 + t402, t525 * t393 + t522 * t399 - t638, t522 * t386 + t423, pkin(5) * t360 - (-pkin(1) * t525 - pkin(6) * t522) * t386, t413, t410, t632, t447, t630, t412, t522 * t339 + t525 * t350 + t402, t522 * t338 + t525 * t477 + t423, t522 * t340 + t525 * t348 + t638, pkin(5) * t316 + t522 * t296 + t525 * t305, t312, t283, t306, t311, t307, t341, t522 * t256 + t525 * t266 + t297, t522 * t259 + t525 * t271 + t302, t522 * t244 + t525 * t258 + t278, pkin(5) * t255 + t525 * t241 + t522 * t242, t312, t283, t306, t311, t307, t341, t522 * t247 + t525 * t257 + t297, t522 * t248 + t525 * t262 + t302, t522 * t239 + t525 * t253 + t278, pkin(5) * t246 + t522 * t237 + t525 * t238; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t494, t495, 0, 0, t431, t542, t433, t430, t625, 0, t572 + t580, -pkin(1) * t481 - t593 - t636, t387 + t571, pkin(1) * t468 + pkin(6) * t387, t431, t433, -t542, 0, -t625, t430, qJ(3) * t591 + t524 * t362 + t572, t524 * t383 + t521 * t384 + t571, t636 + t521 * t361 + (pkin(1) + t618) * t481, pkin(6) * t345 + (pkin(1) - t548) * t382, t323, t290, t334, t322, t335, t346, t524 * t285 + t521 * t295 + t574, t524 * t287 + t521 * t300 + t573, t524 * t263 + t521 * t270 + t575, pkin(1) * t349 + pkin(6) * t261 + t524 * t269 + t521 * t275, t323, t290, t334, t322, t335, t346, t524 * t267 + t521 * t272 + t574, t524 * t268 + t521 * t274 + t573, t524 * t251 + t521 * t252 + t575, pkin(1) * t310 + pkin(6) * t250 + t524 * t240 + t521 * t243;];
tauB_reg = t1;
