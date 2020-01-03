% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S4RRRP7
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S4RRRP7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynB_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynB_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynB_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynB_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:21:09
% EndTime: 2019-12-31 17:21:15
% DurationCPUTime: 5.73s
% Computational Cost: add. (7950->362), mult. (16109->481), div. (0->0), fcn. (10272->6), ass. (0->281)
t471 = sin(qJ(3));
t474 = cos(qJ(3));
t472 = sin(qJ(2));
t517 = qJD(1) * t472;
t435 = t471 * qJD(2) + t474 * t517;
t475 = cos(qJ(2));
t512 = qJD(1) * qJD(2);
t503 = t475 * t512;
t511 = t472 * qJDD(1);
t439 = t503 + t511;
t501 = -t474 * qJDD(2) + t471 * t439;
t515 = t475 * qJD(1);
t458 = -qJD(3) + t515;
t514 = qJD(3) - t458;
t347 = t514 * t435 + t501;
t433 = -t474 * qJD(2) + t471 * t517;
t485 = -t471 * qJDD(2) - t474 * t439;
t479 = t433 * qJD(3) + t485;
t547 = t433 * t458;
t565 = t547 - t479;
t540 = t471 * t565;
t304 = t474 * t347 + t540;
t429 = t435 ^ 2;
t555 = t433 ^ 2;
t394 = t429 - t555;
t284 = t475 * t304 - t472 * t394;
t300 = -t471 * t347 + t474 * t565;
t473 = sin(qJ(1));
t476 = cos(qJ(1));
t635 = t473 * t284 + t476 * t300;
t634 = t476 * t284 - t473 * t300;
t554 = t458 ^ 2;
t412 = t555 - t554;
t461 = t472 * t512;
t509 = t475 * qJDD(1);
t440 = -t461 + t509;
t430 = -qJDD(3) + t440;
t546 = t435 * t433;
t480 = t430 - t546;
t537 = t471 * t480;
t330 = t474 * t412 + t537;
t348 = (qJD(3) + t458) * t435 + t501;
t297 = t475 * t330 - t472 * t348;
t525 = t474 * t480;
t327 = t471 * t412 - t525;
t633 = t473 * t297 - t476 * t327;
t632 = t476 * t297 + t473 * t327;
t564 = t547 + t479;
t589 = -t471 * t348 + t474 * t564;
t566 = t429 + t555;
t588 = -t474 * t348 - t471 * t564;
t604 = -t472 * t566 + t475 * t588;
t613 = t473 * t589 + t476 * t604;
t631 = pkin(4) * t613;
t614 = t473 * t604 - t476 * t589;
t630 = pkin(4) * t614;
t387 = t554 + t429;
t314 = t474 * t387 - t537;
t629 = pkin(1) * t314;
t628 = pkin(2) * t314;
t627 = pkin(5) * t604;
t606 = t472 * t588 + t475 * t566;
t626 = pkin(5) * t606;
t625 = pkin(6) * t314;
t322 = t471 * t387 + t525;
t624 = pkin(6) * t322;
t623 = t472 * t322;
t622 = t473 * t314;
t620 = t475 * t322;
t619 = t476 * t314;
t617 = -pkin(1) * t606 - pkin(2) * t566 - pkin(6) * t588;
t616 = t472 * t304 + t475 * t394;
t615 = t472 * t330 + t475 * t348;
t612 = pkin(6) * t589;
t413 = -t429 + t554;
t481 = -t430 - t546;
t536 = t471 * t481;
t578 = -t474 * t413 - t536;
t524 = t474 * t481;
t577 = -t471 * t413 + t524;
t586 = -t472 * t564 + t475 * t577;
t605 = t473 * t586 + t476 * t578;
t603 = -t473 * t578 + t476 * t586;
t563 = -t554 - t555;
t576 = t471 * t563 + t524;
t602 = pkin(1) * t576;
t601 = pkin(2) * t576;
t573 = t474 * t563 - t536;
t600 = pkin(6) * t573;
t599 = pkin(6) * t576;
t598 = qJ(4) * t565;
t596 = t472 * t573;
t595 = t473 * t576;
t592 = t475 * t573;
t591 = t476 * t576;
t587 = t472 * t577 + t475 * t564;
t585 = -2 * qJD(4);
t451 = t476 * g(1) + t473 * g(2);
t477 = qJD(1) ^ 2;
t424 = -t477 * pkin(1) + qJDD(1) * pkin(5) - t451;
t551 = pkin(2) * t475;
t494 = -pkin(6) * t472 - t551;
t437 = t494 * qJD(1);
t549 = t475 * g(3);
t553 = qJD(2) ^ 2;
t366 = -qJDD(2) * pkin(2) - t553 * pkin(6) + (qJD(1) * t437 + t424) * t472 + t549;
t385 = -t435 * qJD(3) - t501;
t579 = -t385 * pkin(3) + t366 - t598;
t543 = t458 * t474;
t409 = t435 * t543;
t544 = t458 * t471;
t507 = t433 * t544;
t489 = -t409 - t507;
t408 = t435 * t544;
t505 = t433 * t543;
t490 = -t408 + t505;
t557 = -t472 * t430 + t475 * t490;
t575 = t473 * t557 + t476 * t489;
t483 = -t474 * t385 + t507;
t484 = -t471 * t385 - t505;
t506 = t472 * t546;
t558 = t475 * t484 - t506;
t574 = t473 * t558 + t476 * t483;
t572 = -t473 * t489 + t476 * t557;
t571 = -t473 * t483 + t476 * t558;
t450 = t473 * g(1) - t476 * g(2);
t423 = qJDD(1) * pkin(1) + t477 * pkin(5) + t450;
t487 = -t440 + t461;
t488 = t439 + t503;
t345 = pkin(2) * t487 - pkin(6) * t488 - t423;
t407 = -t472 * g(3) + t475 * t424;
t367 = -t553 * pkin(2) + qJDD(2) * pkin(6) + t437 * t515 + t407;
t309 = t471 * t345 + t474 * t367;
t389 = t433 * pkin(3) - t435 * qJ(4);
t486 = -t430 * qJ(4) - t433 * t389 + t458 * t585 + t309;
t341 = t471 * t479 + t409;
t342 = -t474 * t479 + t408;
t491 = t475 * t342 + t506;
t561 = t476 * t341 + t473 * t491;
t560 = t475 * t430 + t472 * t490;
t504 = t475 * t546;
t559 = t472 * t484 + t504;
t556 = -t473 * t341 + t476 * t491;
t552 = pkin(2) * t472;
t550 = pkin(3) * t474;
t548 = qJ(4) * t474;
t545 = t435 * t458;
t467 = t472 ^ 2;
t542 = t467 * t477;
t538 = t471 * t366;
t533 = t472 * t423;
t457 = t475 * t477 * t472;
t447 = qJDD(2) + t457;
t531 = t472 * t447;
t448 = qJDD(2) - t457;
t530 = t472 * t448;
t526 = t474 * t366;
t521 = t475 * t423;
t520 = t475 * t448;
t308 = -t474 * t345 + t471 * t367;
t519 = t566 - t554;
t468 = t475 ^ 2;
t518 = t467 + t468;
t510 = t473 * qJDD(1);
t508 = t476 * qJDD(1);
t502 = qJ(4) * t471 + pkin(2);
t406 = t472 * t424 + t549;
t357 = t472 * t406 + t475 * t407;
t399 = -t473 * t450 - t476 * t451;
t499 = t473 * t457;
t498 = t476 * t457;
t496 = t435 * t389 + qJDD(4) + t308;
t444 = -t473 * t477 + t508;
t493 = -pkin(4) * t444 - t473 * g(3);
t492 = t472 * t342 - t504;
t268 = -t474 * t308 + t471 * t309;
t269 = t471 * t308 + t474 * t309;
t356 = t475 * t406 - t472 * t407;
t398 = t476 * t450 - t473 * t451;
t482 = t430 * pkin(3) + t496;
t478 = 0.2e1 * qJD(4) * t435 - t579;
t465 = t468 * t477;
t455 = -t465 - t553;
t454 = t465 - t553;
t453 = -t542 - t553;
t452 = -t542 + t553;
t446 = t465 - t542;
t445 = t465 + t542;
t443 = t476 * t477 + t510;
t442 = t518 * qJDD(1);
t441 = -0.2e1 * t461 + t509;
t438 = 0.2e1 * t503 + t511;
t432 = t475 * t447;
t431 = t518 * t512;
t420 = -pkin(4) * t443 + t476 * g(3);
t411 = t475 * t439 - t467 * t512;
t410 = -t472 * t440 - t468 * t512;
t405 = -t472 * t453 - t520;
t404 = -t472 * t452 + t432;
t403 = t475 * t455 - t531;
t402 = t475 * t454 - t530;
t401 = t475 * t453 - t530;
t400 = t472 * t455 + t432;
t393 = t476 * t442 - t473 * t445;
t392 = t473 * t442 + t476 * t445;
t388 = -t472 * t438 + t475 * t441;
t373 = t476 * t405 + t473 * t438;
t372 = t476 * t403 - t473 * t441;
t371 = t473 * t405 - t476 * t438;
t370 = t473 * t403 + t476 * t441;
t369 = -pkin(5) * t401 - t521;
t368 = -pkin(5) * t400 - t533;
t359 = -pkin(1) * t401 + t407;
t358 = -pkin(1) * t400 + t406;
t353 = t514 * t433 + t485;
t346 = -t385 - t545;
t319 = t476 * t357 - t473 * t423;
t318 = t473 * t357 + t476 * t423;
t299 = t526 + t625;
t294 = t538 - t599;
t293 = t472 * t346 + t592;
t292 = -t472 * t353 + t620;
t291 = -t475 * t346 + t596;
t290 = t475 * t353 + t623;
t289 = t472 * t347 + t592;
t288 = -t472 * t565 - t620;
t287 = -t475 * t347 + t596;
t286 = t475 * t565 - t623;
t279 = qJ(4) * t554 - t482;
t278 = (-pkin(3) * t458 + t585) * t435 + t579;
t277 = -pkin(3) * t554 + t486;
t276 = t309 + t628;
t275 = t308 - t601;
t274 = t519 * qJ(4) + t482;
t273 = t476 * t293 + t595;
t272 = t476 * t292 - t622;
t271 = t473 * t293 - t591;
t270 = t473 * t292 + t619;
t267 = t519 * pkin(3) + t486;
t266 = t476 * t289 + t595;
t265 = t476 * t288 + t622;
t264 = t473 * t289 - t591;
t263 = t473 * t288 - t619;
t262 = (-t346 + t545) * pkin(3) + t478;
t261 = pkin(3) * t545 + t478 + t598;
t260 = -pkin(2) * t589 - pkin(3) * t564 + qJ(4) * t348;
t259 = t475 * t269 + t472 * t366;
t258 = t472 * t269 - t475 * t366;
t253 = -pkin(1) * t290 - pkin(2) * t353 - t538 - t624;
t252 = -t601 + (-t563 - t554) * qJ(4) + (-t481 + t430) * pkin(3) + t496;
t251 = -pkin(1) * t287 + pkin(2) * t347 + t526 - t600;
t250 = -t628 + qJ(4) * t480 + (-t387 + t554) * pkin(3) - t486;
t249 = -t268 - t612;
t248 = t474 * t277 - t471 * t279;
t247 = t471 * t277 + t474 * t279;
t246 = -t471 * t262 - t346 * t548 - t599;
t245 = -pkin(3) * t540 + t474 * t261 - t625;
t244 = -pkin(5) * t290 - t472 * t276 + t475 * t299;
t243 = -pkin(5) * t287 - t472 * t275 + t475 * t294;
t242 = t476 * t259 + t473 * t268;
t241 = t473 * t259 - t476 * t268;
t240 = t475 * t248 + t472 * t278;
t239 = t472 * t248 - t475 * t278;
t238 = -t269 + t617;
t237 = -pkin(1) * t258 + pkin(2) * t366 - pkin(6) * t269;
t236 = -pkin(1) * t291 - t474 * t262 + t502 * t346 - t600;
t235 = -t471 * t267 + t474 * t274 - t612;
t234 = -pkin(1) * t286 + t624 - t471 * t261 + (-pkin(2) - t550) * t565;
t233 = t475 * t249 + t552 * t589 - t626;
t232 = -pkin(6) * t247 + (pkin(3) * t471 - t548) * t278;
t231 = -pkin(2) * t247 - pkin(3) * t279 - qJ(4) * t277;
t230 = -pkin(5) * t258 + (-pkin(6) * t475 + t552) * t268;
t229 = -t474 * t267 - t471 * t274 + t617;
t228 = -pkin(5) * t291 + t475 * t246 - t472 * t252;
t227 = -pkin(5) * t286 + t475 * t245 - t472 * t250;
t226 = t476 * t240 + t473 * t247;
t225 = t473 * t240 - t476 * t247;
t224 = t475 * t235 - t472 * t260 - t626;
t223 = -pkin(1) * t239 - pkin(6) * t248 + (t502 + t550) * t278;
t222 = -pkin(5) * t239 - t472 * t231 + t475 * t232;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t443, -t444, 0, t399, 0, 0, 0, 0, 0, 0, t372, t373, t393, t319, 0, 0, 0, 0, 0, 0, t266, t272, t613, t242, 0, 0, 0, 0, 0, 0, t273, t613, t265, t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t444, -t443, 0, t398, 0, 0, 0, 0, 0, 0, t370, t371, t392, t318, 0, 0, 0, 0, 0, 0, t264, t270, t614, t241, 0, 0, 0, 0, 0, 0, t271, t614, t263, t225; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t400, t401, 0, -t356, 0, 0, 0, 0, 0, 0, t287, t290, t606, t258, 0, 0, 0, 0, 0, 0, t291, t606, t286, t239; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t444, 0, -t443, 0, t493, -t420, -t398, -pkin(4) * t398, t476 * t411 - t499, t476 * t388 - t473 * t446, t476 * t404 + t472 * t510, t476 * t410 + t499, t476 * t402 + t473 * t509, t473 * qJDD(2) + t476 * t431, -pkin(4) * t370 - t473 * t358 + t476 * t368, -pkin(4) * t371 - t473 * t359 + t476 * t369, -pkin(4) * t392 + t476 * t356, -pkin(4) * t318 - (pkin(1) * t473 - pkin(5) * t476) * t356, t556, -t634, t603, t571, t632, t572, -pkin(4) * t264 + t476 * t243 - t473 * t251, -pkin(4) * t270 + t476 * t244 - t473 * t253, t476 * t233 - t473 * t238 - t630, -pkin(4) * t241 + t476 * t230 - t473 * t237, t556, t603, t634, t572, -t632, t571, -pkin(4) * t271 + t476 * t228 - t473 * t236, t476 * t224 - t473 * t229 - t630, -pkin(4) * t263 + t476 * t227 - t473 * t234, -pkin(4) * t225 + t476 * t222 - t473 * t223; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t443, 0, t444, 0, t420, t493, t399, pkin(4) * t399, t473 * t411 + t498, t473 * t388 + t476 * t446, t473 * t404 - t472 * t508, t473 * t410 - t498, t473 * t402 - t475 * t508, -t476 * qJDD(2) + t473 * t431, pkin(4) * t372 + t476 * t358 + t473 * t368, pkin(4) * t373 + t476 * t359 + t473 * t369, pkin(4) * t393 + t473 * t356, pkin(4) * t319 - (-pkin(1) * t476 - pkin(5) * t473) * t356, t561, -t635, t605, t574, t633, t575, pkin(4) * t266 + t473 * t243 + t476 * t251, pkin(4) * t272 + t473 * t244 + t476 * t253, t473 * t233 + t476 * t238 + t631, pkin(4) * t242 + t473 * t230 + t476 * t237, t561, t605, t635, t575, -t633, t574, pkin(4) * t273 + t473 * t228 + t476 * t236, t473 * t224 + t476 * t229 + t631, pkin(4) * t265 + t473 * t227 + t476 * t234, pkin(4) * t226 + t473 * t222 + t476 * t223; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t450, t451, 0, 0, t488 * t472, t475 * t438 + t472 * t441, t475 * t452 + t531, -t487 * t475, t472 * t454 + t520, 0, pkin(1) * t441 + pkin(5) * t403 + t521, -pkin(1) * t438 + pkin(5) * t405 - t533, pkin(1) * t445 + pkin(5) * t442 + t357, pkin(1) * t423 + pkin(5) * t357, t492, -t616, t587, t559, t615, t560, pkin(5) * t289 + t475 * t275 + t472 * t294 - t602, pkin(5) * t292 + t475 * t276 + t472 * t299 + t629, t627 + t472 * t249 + (-pkin(1) - t551) * t589, pkin(5) * t259 + (-pkin(1) + t494) * t268, t492, t587, t616, t560, -t615, t559, pkin(5) * t293 + t472 * t246 + t475 * t252 - t602, -pkin(1) * t589 + t472 * t235 + t475 * t260 + t627, pkin(5) * t288 + t472 * t245 + t475 * t250 - t629, -pkin(1) * t247 + pkin(5) * t240 + t475 * t231 + t472 * t232;];
tauB_reg = t1;
