% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RPRPPR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:23:37
% EndTime: 2019-05-05 17:23:51
% DurationCPUTime: 7.44s
% Computational Cost: add. (11746->497), mult. (24310->596), div. (0->0), fcn. (11966->6), ass. (0->321)
t502 = qJD(3) ^ 2;
t497 = sin(qJ(3));
t492 = t497 ^ 2;
t503 = qJD(1) ^ 2;
t573 = t492 * t503;
t468 = t502 + t573;
t500 = cos(qJ(3));
t547 = t497 * t500 * t503;
t463 = qJDD(3) - t547;
t576 = t463 * t500;
t406 = -t468 * t497 + t576;
t561 = qJD(1) * qJD(3);
t542 = t500 * t561;
t557 = qJDD(1) * t497;
t451 = 0.2e1 * t542 + t557;
t498 = sin(qJ(1));
t501 = cos(qJ(1));
t354 = t406 * t501 - t451 * t498;
t604 = pkin(6) * t354;
t493 = t500 ^ 2;
t572 = t493 * t503;
t470 = t502 + t572;
t462 = qJDD(3) + t547;
t579 = t462 * t497;
t401 = t470 * t500 + t579;
t481 = t497 * t561;
t555 = qJDD(1) * t500;
t454 = -0.2e1 * t481 + t555;
t355 = t401 * t501 + t454 * t498;
t603 = pkin(6) * t355;
t358 = t406 * t498 + t451 * t501;
t350 = pkin(6) * t358;
t361 = t401 * t498 - t454 * t501;
t351 = pkin(6) * t361;
t578 = t462 * t500;
t408 = -t470 * t497 + t578;
t610 = -pkin(7) - pkin(1);
t551 = t610 * t408;
t633 = pkin(2) * t454 - t551;
t577 = t463 * t497;
t410 = t468 * t500 + t577;
t550 = t410 * t610;
t632 = pkin(2) * t451 - t550;
t536 = -qJ(2) * t454 + t401 * t610;
t535 = qJ(2) * t451 + t406 * t610;
t567 = t492 + t493;
t457 = t567 * qJDD(1);
t460 = t567 * t503;
t534 = -qJ(2) * t460 - t457 * t610;
t582 = t457 * t501;
t387 = -t460 * t498 + t582;
t602 = pkin(6) * t387;
t583 = t457 * t498;
t388 = t460 * t501 + t583;
t383 = pkin(6) * t388;
t496 = sin(qJ(6));
t499 = cos(qJ(6));
t566 = qJD(1) * t497;
t444 = qJD(3) * t499 + t496 * t566;
t446 = -t496 * qJD(3) + t499 * t566;
t391 = t446 * t444;
t453 = -t481 + t555;
t441 = qJDD(6) + t453;
t622 = -t391 + t441;
t631 = t496 * t622;
t630 = t499 * t622;
t464 = g(1) * t498 - t501 * g(2);
t529 = qJDD(2) - t464;
t518 = -qJ(2) * t503 + t529;
t425 = qJDD(1) * t610 + t518;
t384 = t497 * g(3) + t500 * t425;
t527 = pkin(3) * t497 - qJ(4) * t500;
t448 = t527 * qJD(1);
t565 = qJD(1) * t500;
t520 = -qJDD(3) * pkin(3) - t502 * qJ(4) + t448 * t565 + qJDD(4) - t384;
t621 = pkin(2) * t406 + qJ(2) * t410;
t294 = pkin(3) * t463 - qJ(4) * t468 - t520 + t621;
t380 = t451 * t497 - t454 * t500;
t461 = (-t492 + t493) * t503;
t629 = t380 * t498 - t461 * t501;
t346 = t380 * t501 + t461 * t498;
t469 = -t502 + t573;
t405 = t469 * t497 + t578;
t554 = qJDD(1) * t501;
t369 = t405 * t498 - t497 * t554;
t556 = qJDD(1) * t498;
t628 = t405 * t501 + t497 * t556;
t627 = -pkin(2) * t401 + qJ(2) * t408;
t452 = t542 + t557;
t626 = (t451 + t452) * pkin(3);
t625 = (t460 - t502) * pkin(3);
t624 = t503 * t610;
t526 = t453 + t481;
t623 = t526 * qJ(5);
t546 = qJ(5) * t565;
t521 = -qJD(3) * pkin(4) - t546;
t620 = -t452 * qJ(5) - qJD(3) * t521;
t619 = -pkin(3) * t542 - qJ(4) * t481;
t541 = g(3) * t500 - t497 * t425;
t618 = -qJDD(3) * qJ(4) - 0.2e1 * qJD(4) * qJD(3) + t541;
t560 = qJD(1) * qJD(5);
t617 = 0.2e1 * t497 * t560 - t620;
t616 = -t452 * pkin(4) - qJ(5) * t573 + qJDD(5);
t439 = t444 ^ 2;
t440 = t446 ^ 2;
t475 = qJD(6) + t565;
t472 = t475 ^ 2;
t615 = -2 * qJD(2);
t614 = 0.2e1 * qJD(4);
t613 = pkin(3) + pkin(4);
t612 = pkin(3) + pkin(8);
t611 = -pkin(4) - pkin(8);
t607 = pkin(2) * t457;
t606 = pkin(2) * t460;
t605 = pkin(3) * t500;
t601 = t452 * pkin(3);
t600 = -qJ(4) - pkin(5);
t598 = qJ(4) * t460;
t597 = qJ(4) * t497;
t596 = qJDD(1) * pkin(1);
t512 = t618 + t620;
t532 = pkin(5) * t500 - pkin(8) * t497;
t539 = (-0.2e1 * qJD(5) + t448) * qJD(1);
t552 = pkin(4) * t573;
t302 = -qJDD(3) * pkin(5) + t612 * t502 + t512 + t552 + (-t503 * t532 + t539) * t497;
t595 = t302 * t496;
t594 = t302 * t499;
t365 = t391 + t441;
t593 = t365 * t499;
t465 = t501 * g(1) + t498 * g(2);
t490 = qJDD(1) * qJ(2);
t522 = t465 - t490;
t513 = t522 - t624;
t559 = qJD(2) * qJD(1);
t418 = t513 - 0.2e1 * t559;
t592 = t418 * t497;
t591 = t444 * t475;
t589 = t451 * t500;
t587 = t453 * qJ(4);
t586 = t454 * t497;
t575 = t475 * t496;
t574 = t475 * t499;
t571 = t496 * t365;
t570 = t500 * t418;
t564 = qJD(3) * t497;
t563 = qJD(4) * t500;
t562 = pkin(3) - t611;
t553 = -t440 - t472;
t549 = t497 * t391;
t548 = t500 * t391;
t545 = t500 * t560;
t540 = (t454 + t453) * qJ(4);
t486 = 0.2e1 * t559;
t519 = t486 - t522;
t426 = -t503 * pkin(1) + t519;
t427 = -t518 + t596;
t357 = t501 * t426 - t427 * t498;
t398 = -t464 * t498 - t501 * t465;
t538 = t498 * t547;
t537 = t501 * t547;
t458 = -t498 * t503 + t554;
t531 = pkin(6) * t458 + g(3) * t498;
t459 = t501 * t503 + t556;
t530 = -pkin(6) * t459 + g(3) * t501;
t528 = -t597 - t605;
t376 = -t444 * qJD(6) - qJDD(3) * t496 + t499 * t452;
t510 = t513 + t619;
t507 = t510 + t616;
t283 = -t600 * t453 - t612 * t452 + (-pkin(5) * t564 + t615 + (qJD(3) * t611 - t546 + t614) * t500) * qJD(1) + t507;
t479 = -0.2e1 * t545;
t467 = pkin(4) * t547;
t517 = -t467 - t520;
t303 = -t502 * pkin(5) - t453 * qJ(5) + t479 + t611 * qJDD(3) + (-qJ(5) * t564 - t532 * t565) * qJD(1) - t517;
t260 = -t283 * t499 + t303 * t496;
t261 = t283 * t496 + t303 * t499;
t245 = -t260 * t499 + t261 * t496;
t246 = t496 * t260 + t261 * t499;
t331 = t384 * t500 - t497 * t541;
t332 = -t384 * t497 - t500 * t541;
t352 = t426 * t498 + t427 * t501;
t525 = t586 + t589;
t524 = -t469 * t500 + t579;
t397 = t464 * t501 - t465 * t498;
t523 = qJDD(3) * t499 + t452 * t496;
t516 = t376 - t591;
t514 = -t448 * t566 - t618;
t511 = (-qJD(6) + t475) * t446 - t523;
t344 = -pkin(3) * t502 + t514;
t508 = qJDD(3) * pkin(4) + t517 + 0.2e1 * t545;
t292 = pkin(3) * t470 + qJ(4) * t462 + t344 - t627;
t506 = 0.2e1 * (-qJD(2) + t563) * qJD(1) + t510;
t505 = t506 + t587;
t504 = (t615 + (t614 + t521) * t500) * qJD(1) - t601 + t507;
t471 = t502 - t572;
t443 = t567 * t561;
t424 = -t440 + t472;
t423 = t439 - t472;
t422 = qJDD(3) * t501 - t443 * t498;
t421 = qJDD(3) * t498 + t443 * t501;
t420 = -t453 * t497 - t493 * t561;
t419 = t452 * t500 - t492 * t561;
t413 = (t453 - t481) * t500;
t409 = -t471 * t497 + t576;
t402 = -t471 * t500 - t577;
t400 = (t452 + t542) * t497;
t399 = qJDD(1) * t528 - t607;
t396 = qJ(4) * t451 - qJ(5) * t463;
t390 = t440 - t439;
t382 = t607 + (t500 * t613 + t597) * qJDD(1);
t377 = -t472 - t439;
t375 = -qJD(6) * t446 - t523;
t374 = -t419 * t498 - t537;
t373 = -t420 * t498 + t537;
t372 = t419 * t501 - t538;
t371 = t420 * t501 + t538;
t370 = -t402 * t498 + t500 * t554;
t368 = t402 * t501 + t498 * t555;
t366 = -qJ(5) * t462 + t454 * t613;
t363 = -t439 - t440;
t349 = (t444 * t499 - t446 * t496) * t475;
t348 = (t444 * t496 + t446 * t499) * t475;
t343 = t376 + t591;
t338 = (qJD(6) + t475) * t446 + t523;
t336 = -t376 * t499 + t446 * t575;
t335 = -t376 * t496 - t446 * t574;
t334 = t375 * t496 - t444 * t574;
t333 = -t375 * t499 - t444 * t575;
t330 = t520 + t598;
t329 = -t349 * t500 - t441 * t497;
t328 = t514 + t625;
t327 = -t423 * t499 + t571;
t326 = t424 * t496 - t630;
t325 = -t423 * t496 - t593;
t324 = -t424 * t499 - t631;
t323 = t505 - t601;
t322 = -t496 * t553 - t593;
t321 = t499 * t553 - t571;
t320 = -t332 - t606;
t319 = t541 + t627;
t318 = t384 + t621;
t317 = t377 * t499 - t631;
t316 = t496 * t377 + t630;
t315 = t505 - t626;
t314 = t506 + t540 - t601;
t313 = -t570 + t632;
t312 = t592 + t633;
t311 = t508 + t623;
t310 = t344 - t552 + t617;
t309 = t331 * t498 - t418 * t501;
t308 = -t331 * t501 - t418 * t498;
t307 = -t334 * t500 + t549;
t306 = -t336 * t500 - t549;
t305 = t504 + t587;
t304 = -t598 + (t526 + t555) * qJ(5) + t508;
t301 = t344 * t500 + t497 * t520;
t300 = t344 * t497 - t500 * t520;
t299 = t496 * t343 + t499 * t511;
t298 = t338 * t499 + t496 * t516;
t297 = -t343 * t499 + t496 * t511;
t296 = t338 * t496 - t499 * t516;
t295 = (-t460 + t573) * pkin(4) - t625 + (-qJ(5) * qJDD(1) + t539) * t497 + t512;
t293 = qJ(5) * t470 + t504 + t540;
t291 = -t326 * t500 - t343 * t497;
t290 = -t327 * t500 - t497 * t511;
t289 = pkin(2) * t331 - qJ(2) * t332;
t288 = t322 * t497 + t500 * t516;
t287 = -t322 * t500 + t497 * t516;
t286 = pkin(4) * t451 - qJ(5) * t468 - 0.2e1 * qJD(1) * t563 - t521 * t565 + t519 - t587 - t616 - t619 + t624 + t626;
t285 = t317 * t497 + t338 * t500;
t284 = -t317 * t500 + t338 * t497;
t282 = -t328 * t500 - t330 * t497 - t606;
t281 = -t298 * t500 - t390 * t497;
t280 = -pkin(2) * t418 + t332 * t610;
t279 = t299 * t497 + t363 * t500;
t278 = -t299 * t500 + t363 * t497;
t277 = -t314 * t497 + (-pkin(2) - t605) * t454 + t551;
t276 = -t315 * t500 + (pkin(2) + t597) * t451 - t550;
t275 = t292 + (t470 - t573) * pkin(4) + t617;
t274 = t467 + t479 - t623 + (-qJDD(3) - t463) * pkin(4) - t294;
t273 = t310 * t500 - t311 * t497;
t272 = t310 * t497 + t311 * t500;
t271 = t300 * t498 - t323 * t501;
t270 = -t300 * t501 - t323 * t498;
t269 = qJ(4) * t305 + qJ(5) * t311;
t268 = t287 * t498 - t321 * t501;
t267 = -t287 * t501 - t321 * t498;
t266 = t284 * t498 - t316 * t501;
t265 = -t284 * t501 - t316 * t498;
t264 = -t286 * t500 - t396 * t497 - t632;
t263 = -t293 * t497 - t366 * t500 - t633;
t262 = -t295 * t500 - t304 * t497 + t606;
t258 = t278 * t498 - t297 * t501;
t257 = -t278 * t501 - t297 * t498;
t256 = t272 * t498 - t305 * t501;
t255 = -t272 * t501 - t305 * t498;
t254 = -qJ(5) * t310 + t305 * t613;
t253 = pkin(2) * t300 - pkin(3) * t520 - qJ(2) * t301 + qJ(4) * t344;
t252 = -qJ(5) * t299 - t297 * t600;
t251 = -qJ(5) * t516 + t321 * t562 + t594;
t250 = -qJ(5) * t338 + t316 * t562 + t595;
t249 = t610 * t301 + (-pkin(2) + t528) * t323;
t248 = -qJ(5) * t322 - t321 * t600 - t261;
t247 = -qJ(5) * t317 - t316 * t600 - t260;
t244 = t246 * t497 - t302 * t500;
t243 = -t246 * t500 - t302 * t497;
t242 = pkin(2) * t272 - qJ(2) * t273 + qJ(4) * t310 + t311 * t613;
t241 = pkin(2) * t287 - qJ(2) * t288 - t322 * t562 - t516 * t600 + t595;
t240 = pkin(2) * t284 - qJ(2) * t285 - t317 * t562 - t338 * t600 - t594;
t239 = -qJ(5) * t363 + t297 * t562 + t245;
t238 = -pkin(2) * t305 - t254 * t500 - t269 * t497 + t273 * t610;
t237 = t243 * t498 - t245 * t501;
t236 = -t243 * t501 - t245 * t498;
t235 = -pkin(2) * t321 - t248 * t497 - t251 * t500 + t288 * t610;
t234 = pkin(2) * t278 - qJ(2) * t279 - t299 * t562 - t363 * t600 - t246;
t233 = -pkin(2) * t316 - t247 * t497 - t250 * t500 + t285 * t610;
t232 = -qJ(5) * t246 - t245 * t600;
t231 = qJ(5) * t302 + t245 * t562;
t230 = -pkin(2) * t297 - t239 * t500 - t252 * t497 + t279 * t610;
t229 = pkin(2) * t243 - qJ(2) * t244 - t246 * t562 + t302 * t600;
t228 = -pkin(2) * t245 - t231 * t500 - t232 * t497 + t244 * t610;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t459, -t458, 0, t398, 0, 0, 0, 0, 0, 0, 0, t459, t458, t357, 0, 0, 0, 0, 0, 0, t358, -t361, -t388, t309, 0, 0, 0, 0, 0, 0, t358, -t388, t361, t271, 0, 0, 0, 0, 0, 0, t361, -t358, t388, t256, 0, 0, 0, 0, 0, 0, t266, t268, t258, t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t458, -t459, 0, t397, 0, 0, 0, 0, 0, 0, 0, -t458, t459, t352, 0, 0, 0, 0, 0, 0, -t354, t355, t387, t308, 0, 0, 0, 0, 0, 0, -t354, t387, -t355, t270, 0, 0, 0, 0, 0, 0, -t355, t354, -t387, t255, 0, 0, 0, 0, 0, 0, t265, t267, t257, t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t410, -t408, 0, t332, 0, 0, 0, 0, 0, 0, -t410, 0, t408, t301, 0, 0, 0, 0, 0, 0, t408, t410, 0, t273, 0, 0, 0, 0, 0, 0, t285, t288, t279, t244; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t458, 0, -t459, 0, -t531, -t530, -t397, -pkin(6) * t397, 0, -t458, t459, 0, 0, 0, -t352, t531, t530, -pkin(6) * t352 + (-pkin(1) * t498 + qJ(2) * t501) * g(3), t373, -t629, t370, t374, t369, t422, -t313 * t498 + t318 * t501 + t604, -t312 * t498 + t319 * t501 - t603, -pkin(2) * t582 - t320 * t498 - t602, -pkin(6) * t308 - t280 * t498 + t289 * t501, t373, t370, t629, t422, -t369, t374, -t276 * t498 + t294 * t501 + t604, -t282 * t498 + t399 * t501 - t602, -t277 * t498 + t292 * t501 + t603, -pkin(6) * t270 - t249 * t498 + t253 * t501, t374, -t629, t369, t373, t370, t422, -t263 * t498 + t275 * t501 + t603, -t264 * t498 + t274 * t501 - t604, -t262 * t498 + t382 * t501 + t602, -pkin(6) * t255 - t238 * t498 + t242 * t501, -t306 * t498 + t335 * t501, -t281 * t498 + t296 * t501, -t291 * t498 + t324 * t501, -t307 * t498 + t333 * t501, -t290 * t498 + t325 * t501, -t329 * t498 + t348 * t501, -pkin(6) * t265 - t233 * t498 + t240 * t501, -pkin(6) * t267 - t235 * t498 + t241 * t501, -pkin(6) * t257 - t230 * t498 + t234 * t501, -pkin(6) * t236 - t228 * t498 + t229 * t501; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t459, 0, t458, 0, t530, -t531, t398, pkin(6) * t398, 0, -t459, -t458, 0, 0, 0, t357, -t530, t531, pkin(6) * t357 + (pkin(1) * t501 + qJ(2) * t498) * g(3), t371, t346, t368, t372, -t628, t421, t313 * t501 + t318 * t498 + t350, t312 * t501 + t319 * t498 - t351, -pkin(2) * t583 + t320 * t501 - t383, pkin(6) * t309 + t280 * t501 + t289 * t498, t371, t368, -t346, t421, t628, t372, t276 * t501 + t294 * t498 + t350, t282 * t501 + t399 * t498 - t383, t277 * t501 + t292 * t498 + t351, pkin(6) * t271 + t249 * t501 + t253 * t498, t372, t346, -t628, t371, t368, t421, t263 * t501 + t275 * t498 + t351, t264 * t501 + t274 * t498 - t350, t262 * t501 + t382 * t498 + t383, pkin(6) * t256 + t238 * t501 + t242 * t498, t306 * t501 + t335 * t498, t281 * t501 + t296 * t498, t291 * t501 + t324 * t498, t307 * t501 + t333 * t498, t290 * t501 + t325 * t498, t329 * t501 + t348 * t498, pkin(6) * t266 + t233 * t501 + t240 * t498, pkin(6) * t268 + t235 * t501 + t241 * t498, pkin(6) * t258 + t230 * t501 + t234 * t498, pkin(6) * t237 + t228 * t501 + t229 * t498; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t464, t465, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t529 - 0.2e1 * t596, -t465 + t486 + 0.2e1 * t490, pkin(1) * t427 + qJ(2) * t426, t413, -t525, t409, t400, -t524, 0, t535 - t592, -t536 - t570, -t331 + t534, -qJ(2) * t418 + t331 * t610, t413, t409, t525, 0, t524, t400, -qJ(4) * t589 - t315 * t497 + t535, -t328 * t497 + t330 * t500 + t534, -pkin(3) * t586 + t314 * t500 + t536, t610 * t300 + (-qJ(2) - t527) * t323, t400, -t525, -t524, t413, t409, 0, t293 * t500 - t366 * t497 + t536, -t286 * t497 + t396 * t500 - t535, -t295 * t497 + t304 * t500 - t534, -qJ(2) * t305 - t254 * t497 + t269 * t500 + t272 * t610, -t336 * t497 + t548, -t298 * t497 + t390 * t500, -t326 * t497 + t343 * t500, -t334 * t497 - t548, -t327 * t497 + t500 * t511, -t349 * t497 + t441 * t500, -qJ(2) * t316 + t247 * t500 - t250 * t497 + t284 * t610, -qJ(2) * t321 + t248 * t500 - t251 * t497 + t287 * t610, -qJ(2) * t297 - t239 * t497 + t252 * t500 + t278 * t610, -qJ(2) * t245 - t231 * t497 + t232 * t500 + t243 * t610;];
tauB_reg  = t1;
