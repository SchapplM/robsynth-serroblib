% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPRPR14_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:35:24
% EndTime: 2019-12-31 18:35:32
% DurationCPUTime: 7.81s
% Computational Cost: add. (21518->488), mult. (47267->712), div. (0->0), fcn. (30944->8), ass. (0->328)
t501 = sin(pkin(8));
t502 = cos(pkin(8));
t505 = sin(qJ(3));
t508 = cos(qJ(3));
t463 = (t501 * t508 + t502 * t505) * qJD(1);
t550 = qJD(1) * t508;
t465 = -qJD(1) * t501 * t505 + t502 * t550;
t562 = t465 * t463;
t586 = qJDD(3) - t562;
t588 = t501 * t586;
t587 = t502 * t586;
t543 = qJD(1) * qJD(3);
t529 = t508 * t543;
t541 = qJDD(1) * t505;
t471 = -t529 - t541;
t530 = t505 * t543;
t539 = qJDD(1) * t508;
t472 = -t530 + t539;
t421 = t471 * t501 + t472 * t502;
t549 = qJD(3) * t463;
t398 = t421 - t549;
t504 = sin(qJ(5));
t507 = cos(qJ(5));
t432 = -t507 * qJD(3) + t465 * t504;
t434 = qJD(3) * t504 + t465 * t507;
t394 = t434 * t432;
t526 = -t502 * t471 + t472 * t501;
t519 = qJDD(5) + t526;
t583 = -t394 + t519;
t585 = t504 * t583;
t584 = t507 * t583;
t548 = qJD(3) * t465;
t395 = t526 + t548;
t373 = -t432 * qJD(5) + t504 * qJDD(3) + t507 * t421;
t455 = qJD(5) + t463;
t410 = t455 * t432;
t348 = -t410 + t373;
t511 = qJD(1) ^ 2;
t506 = sin(qJ(1));
t509 = cos(qJ(1));
t481 = g(1) * t506 - t509 * g(2);
t520 = qJDD(2) - t481;
t514 = -qJ(2) * t511 + t520;
t580 = pkin(6) + pkin(1);
t446 = -t580 * qJDD(1) + t514;
t424 = -t508 * g(3) + t505 * t446;
t516 = qJD(3) * pkin(3) - qJ(4) * t550;
t498 = t505 ^ 2;
t555 = t498 * t511;
t386 = -pkin(3) * t555 + t471 * qJ(4) - qJD(3) * t516 + t424;
t553 = t508 * t511;
t565 = t446 * t508;
t512 = qJDD(3) * pkin(3) - qJ(4) * t472 + t565 + (-pkin(3) * t553 - qJ(4) * t543 + g(3)) * t505;
t325 = -0.2e1 * qJD(4) * t463 + t502 * t386 + t501 * t512;
t527 = -t507 * qJDD(3) + t421 * t504;
t345 = (qJD(5) - t455) * t434 + t527;
t482 = g(1) * t509 + g(2) * t506;
t497 = qJDD(1) * qJ(2);
t517 = t482 - t497;
t582 = -pkin(3) * t471 - (qJ(4) * t498 + t580) * t511 + t516 * t550 + qJDD(4) - t517;
t430 = t432 ^ 2;
t431 = t434 ^ 2;
t454 = t455 ^ 2;
t461 = t463 ^ 2;
t462 = t465 ^ 2;
t581 = 0.2e1 * qJD(4);
t579 = pkin(4) * t501;
t578 = qJDD(1) * pkin(1);
t528 = t386 * t501 - t502 * t512;
t324 = t465 * t581 + t528;
t278 = -t324 * t502 + t325 * t501;
t577 = t278 * t505;
t576 = t278 * t508;
t411 = pkin(4) * t463 - pkin(7) * t465;
t510 = qJD(3) ^ 2;
t303 = -qJDD(3) * pkin(4) - t510 * pkin(7) + (t581 + t411) * t465 + t528;
t575 = t303 * t504;
t574 = t303 * t507;
t363 = t394 + t519;
t573 = t363 * t504;
t572 = t363 * t507;
t542 = qJD(2) * qJD(1);
t537 = -0.2e1 * t542;
t387 = t537 - t582;
t571 = t387 * t501;
t570 = t387 * t502;
t414 = qJDD(3) + t562;
t569 = t414 * t501;
t568 = t414 * t502;
t443 = t580 * t511 + t517 + t537;
t567 = t443 * t505;
t566 = t443 * t508;
t564 = t455 * t504;
t563 = t455 * t507;
t499 = t508 ^ 2;
t551 = t498 + t499;
t474 = t551 * qJDD(1);
t561 = t474 * t506;
t560 = t474 * t509;
t532 = t505 * t553;
t479 = qJDD(3) + t532;
t559 = t479 * t505;
t558 = t479 * t508;
t480 = qJDD(3) - t532;
t557 = t480 * t505;
t556 = t480 * t508;
t554 = t499 * t511;
t304 = -pkin(4) * t510 + qJDD(3) * pkin(7) - t411 * t463 + t325;
t495 = 0.2e1 * t542;
t317 = t395 * pkin(4) - t398 * pkin(7) + t495 + t582;
t270 = t507 * t304 + t504 * t317;
t547 = qJD(3) * t501;
t546 = qJD(3) * t502;
t540 = qJDD(1) * t506;
t538 = qJDD(1) * t509;
t536 = t501 * t394;
t535 = t502 * t394;
t534 = t506 * t562;
t533 = t509 * t562;
t531 = -pkin(4) * t502 - pkin(3);
t269 = t304 * t504 - t507 * t317;
t241 = t269 * t504 + t507 * t270;
t279 = t324 * t501 + t502 * t325;
t450 = -pkin(1) * t511 + t495 - t517;
t456 = -t514 + t578;
t405 = t509 * t450 - t456 * t506;
t429 = -t481 * t506 - t509 * t482;
t525 = t506 * t532;
t524 = t509 * t532;
t475 = -t506 * t511 + t538;
t522 = pkin(5) * t475 + g(3) * t506;
t476 = t509 * t511 + t540;
t521 = -pkin(5) * t476 + g(3) * t509;
t423 = g(3) * t505 + t565;
t240 = -t269 * t507 + t270 * t504;
t376 = t423 * t508 + t424 * t505;
t377 = -t423 * t505 + t424 * t508;
t402 = t450 * t506 + t456 * t509;
t428 = t481 * t509 - t482 * t506;
t515 = -t526 + t548;
t492 = t509 * qJDD(3);
t490 = t506 * qJDD(3);
t487 = -t510 - t554;
t486 = t510 - t554;
t485 = -t510 - t555;
t484 = -t510 + t555;
t478 = (-t498 + t499) * t511;
t477 = t551 * t511;
t473 = -0.2e1 * t530 + t539;
t470 = 0.2e1 * t529 + t541;
t468 = t551 * t543;
t449 = -t462 - t510;
t448 = -t462 + t510;
t447 = t461 - t510;
t445 = -t472 * t505 - t499 * t543;
t444 = -t471 * t508 - t498 * t543;
t440 = -t487 * t505 - t558;
t439 = t485 * t508 - t557;
t438 = t487 * t508 - t559;
t437 = -t486 * t508 - t557;
t436 = t485 * t505 + t556;
t435 = -t484 * t505 - t558;
t426 = -t477 * t509 - t561;
t425 = -t477 * t506 + t560;
t422 = t470 * t505 - t473 * t508;
t417 = t462 - t461;
t412 = -t510 - t461;
t409 = (-t463 * t502 + t465 * t501) * qJD(3);
t408 = (-t463 * t501 - t465 * t502) * qJD(3);
t407 = t438 * t506 + t473 * t509;
t406 = t436 * t506 + t470 * t509;
t404 = -t438 * t509 + t473 * t506;
t403 = -t436 * t509 + t470 * t506;
t401 = -t431 + t454;
t400 = t430 - t454;
t399 = t421 + t549;
t393 = -t461 - t462;
t392 = -t431 + t430;
t391 = t421 * t502 - t465 * t547;
t390 = t421 * t501 + t465 * t546;
t389 = t463 * t546 + t501 * t526;
t388 = t463 * t547 - t502 * t526;
t385 = -t431 - t454;
t384 = -t449 * t501 - t568;
t383 = -t448 * t501 + t587;
t382 = t447 * t502 - t569;
t381 = t449 * t502 - t569;
t380 = t448 * t502 + t588;
t379 = t447 * t501 + t568;
t372 = -qJD(5) * t434 - t527;
t371 = -t454 - t430;
t370 = -pkin(2) * t477 - t377;
t369 = t430 + t431;
t368 = t412 * t502 - t588;
t367 = t412 * t501 + t587;
t366 = pkin(2) * t438 - qJ(2) * t440 - t424;
t365 = pkin(2) * t436 - qJ(2) * t439 + t423;
t361 = pkin(2) * t470 - t580 * t439 - t566;
t360 = pkin(2) * t473 - t580 * t440 + t567;
t359 = (-t432 * t507 + t434 * t504) * t455;
t358 = (-t432 * t504 - t434 * t507) * t455;
t357 = -t408 * t508 - t409 * t505;
t356 = t376 * t506 - t443 * t509;
t355 = -t376 * t509 - t443 * t506;
t354 = t399 * t501 + t502 * t515;
t353 = -t395 * t502 - t398 * t501;
t352 = -t399 * t502 + t501 * t515;
t351 = -t395 * t501 + t398 * t502;
t349 = -t410 - t373;
t346 = (-qJD(5) - t455) * t434 - t527;
t344 = t373 * t507 - t434 * t564;
t343 = t373 * t504 + t434 * t563;
t342 = -t372 * t504 + t432 * t563;
t341 = t372 * t507 + t432 * t564;
t340 = -t390 * t508 - t391 * t505;
t339 = -t388 * t508 - t389 * t505;
t338 = -qJ(4) * t381 - t570;
t337 = pkin(2) * t376 - qJ(2) * t377;
t336 = -t381 * t505 + t384 * t508;
t335 = t381 * t508 + t384 * t505;
t334 = -t380 * t508 - t383 * t505;
t333 = -t379 * t508 - t382 * t505;
t332 = t359 * t502 + t501 * t519;
t331 = t359 * t501 - t502 * t519;
t330 = t400 * t507 - t573;
t329 = -t401 * t504 + t584;
t328 = t400 * t504 + t572;
t327 = t401 * t507 + t585;
t326 = -qJ(4) * t367 - t571;
t322 = -t385 * t504 - t572;
t321 = t385 * t507 - t573;
t320 = -pkin(2) * t443 - t580 * t377;
t319 = t371 * t507 - t585;
t318 = t371 * t504 + t584;
t314 = -t367 * t505 + t368 * t508;
t313 = t367 * t508 + t368 * t505;
t312 = t344 * t502 + t536;
t311 = t342 * t502 - t536;
t310 = t344 * t501 - t535;
t309 = t342 * t501 + t535;
t308 = -pkin(3) * t398 + qJ(4) * t384 - t571;
t307 = t335 * t506 + t398 * t509;
t306 = -t335 * t509 + t398 * t506;
t305 = -pkin(3) * t395 + qJ(4) * t368 + t570;
t301 = t313 * t506 + t395 * t509;
t300 = -t352 * t505 + t354 * t508;
t299 = -t313 * t509 + t395 * t506;
t298 = t352 * t508 + t354 * t505;
t297 = -t351 * t508 - t353 * t505;
t296 = -t345 * t507 - t349 * t504;
t295 = t346 * t507 - t348 * t504;
t294 = -t345 * t504 + t349 * t507;
t293 = t346 * t504 + t348 * t507;
t292 = t330 * t502 - t345 * t501;
t291 = t329 * t502 - t349 * t501;
t290 = t330 * t501 + t345 * t502;
t289 = t329 * t501 + t349 * t502;
t288 = t322 * t502 + t348 * t501;
t287 = t322 * t501 - t348 * t502;
t286 = -t331 * t508 - t332 * t505;
t285 = t319 * t502 - t346 * t501;
t284 = t319 * t501 + t346 * t502;
t283 = t298 * t506 + t393 * t509;
t282 = -t298 * t509 + t393 * t506;
t281 = t295 * t502 - t392 * t501;
t280 = t295 * t501 + t392 * t502;
t277 = t296 * t502 - t369 * t501;
t276 = t296 * t501 + t369 * t502;
t275 = -t310 * t508 - t312 * t505;
t274 = -t309 * t508 - t311 * t505;
t273 = pkin(3) * t387 + qJ(4) * t279;
t272 = -pkin(7) * t321 + t574;
t271 = -pkin(7) * t318 + t575;
t267 = -qJ(4) * t352 - t278;
t266 = -pkin(3) * t393 + qJ(4) * t354 + t279;
t265 = pkin(2) * t335 + pkin(3) * t381 - qJ(2) * t336 - t325;
t264 = -t290 * t508 - t292 * t505;
t263 = -t289 * t508 - t291 * t505;
t262 = -pkin(4) * t321 + t270;
t261 = -pkin(4) * t318 + t269;
t260 = -t287 * t505 + t288 * t508;
t259 = t287 * t508 + t288 * t505;
t258 = pkin(2) * t298 + pkin(3) * t352 - qJ(2) * t300;
t257 = pkin(2) * t313 + pkin(3) * t367 - qJ(2) * t314 - t324;
t256 = -t284 * t505 + t285 * t508;
t255 = t284 * t508 + t285 * t505;
t254 = -t280 * t508 - t281 * t505;
t253 = t279 * t508 - t577;
t252 = t279 * t505 + t576;
t251 = pkin(2) * t398 - t308 * t508 - t580 * t336 - t338 * t505;
t250 = -t276 * t505 + t277 * t508;
t249 = t276 * t508 + t277 * t505;
t248 = t252 * t506 - t387 * t509;
t247 = -t252 * t509 - t387 * t506;
t246 = t259 * t506 + t321 * t509;
t245 = -t259 * t509 + t321 * t506;
t244 = pkin(2) * t395 - t305 * t508 - t580 * t314 - t326 * t505;
t243 = t255 * t506 + t318 * t509;
t242 = -t255 * t509 + t318 * t506;
t239 = t249 * t506 + t294 * t509;
t238 = -t249 * t509 + t294 * t506;
t237 = t241 * t502 + t303 * t501;
t236 = t241 * t501 - t303 * t502;
t235 = -pkin(7) * t294 - t240;
t234 = -qJ(4) * t287 - t262 * t501 + t272 * t502;
t233 = -qJ(4) * t284 - t261 * t501 + t271 * t502;
t232 = -pkin(3) * t321 + qJ(4) * t288 + t262 * t502 + t272 * t501;
t231 = -pkin(3) * t318 + qJ(4) * t285 + t261 * t502 + t271 * t501;
t230 = pkin(2) * t393 - t266 * t508 - t267 * t505 - t580 * t300;
t229 = pkin(2) * t252 + pkin(3) * t278 - qJ(2) * t253;
t228 = -qJ(4) * t276 + t235 * t502 + t294 * t579;
t227 = qJ(4) * t277 + t235 * t501 + t531 * t294;
t226 = pkin(2) * t259 + pkin(3) * t287 - pkin(4) * t348 + pkin(7) * t322 - qJ(2) * t260 + t575;
t225 = pkin(2) * t255 + pkin(3) * t284 + pkin(4) * t346 + pkin(7) * t319 - qJ(2) * t256 - t574;
t224 = -t236 * t505 + t237 * t508;
t223 = t236 * t508 + t237 * t505;
t222 = -pkin(2) * t387 + qJ(4) * t577 - t580 * t253 - t273 * t508;
t221 = -qJ(4) * t236 + (-pkin(7) * t502 + t579) * t240;
t220 = pkin(2) * t249 + pkin(3) * t276 + pkin(4) * t369 + pkin(7) * t296 - qJ(2) * t250 + t241;
t219 = t223 * t506 + t240 * t509;
t218 = -t223 * t509 + t240 * t506;
t217 = qJ(4) * t237 + (-pkin(7) * t501 + t531) * t240;
t216 = pkin(2) * t321 - t232 * t508 - t234 * t505 - t580 * t260;
t215 = pkin(2) * t318 - t231 * t508 - t233 * t505 - t580 * t256;
t214 = pkin(2) * t294 - t227 * t508 - t228 * t505 - t580 * t250;
t213 = pkin(2) * t223 + pkin(3) * t236 - pkin(4) * t303 + pkin(7) * t241 - qJ(2) * t224;
t212 = pkin(2) * t240 - t217 * t508 - t221 * t505 - t580 * t224;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t476, -t475, 0, t429, 0, 0, 0, 0, 0, 0, 0, t476, t475, t405, 0, 0, 0, 0, 0, 0, t406, t407, t426, t356, 0, 0, 0, 0, 0, 0, t301, t307, t283, t248, 0, 0, 0, 0, 0, 0, t243, t246, t239, t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t475, -t476, 0, t428, 0, 0, 0, 0, 0, 0, 0, -t475, t476, t402, 0, 0, 0, 0, 0, 0, t403, t404, t425, t355, 0, 0, 0, 0, 0, 0, t299, t306, t282, t247, 0, 0, 0, 0, 0, 0, t242, t245, t238, t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t439, t440, 0, t377, 0, 0, 0, 0, 0, 0, t314, t336, t300, t253, 0, 0, 0, 0, 0, 0, t256, t260, t250, t224; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t475, 0, -t476, 0, -t522, -t521, -t428, -pkin(5) * t428, 0, -t475, t476, 0, 0, 0, -t402, t522, t521, -pkin(5) * t402 + (-pkin(1) * t506 + qJ(2) * t509) * g(3), -t445 * t506 + t524, -t422 * t506 + t478 * t509, -t437 * t506 + t508 * t538, -t444 * t506 - t524, -t435 * t506 - t505 * t538, -t468 * t506 + t492, -pkin(5) * t403 - t361 * t506 + t365 * t509, -pkin(5) * t404 - t360 * t506 + t366 * t509, -pkin(2) * t560 - pkin(5) * t425 - t370 * t506, -pkin(5) * t355 - t320 * t506 + t337 * t509, -t340 * t506 + t533, -t297 * t506 + t417 * t509, -t334 * t506 + t399 * t509, -t339 * t506 - t533, -t333 * t506 + t509 * t515, -t357 * t506 + t492, -pkin(5) * t299 - t244 * t506 + t257 * t509, -pkin(5) * t306 - t251 * t506 + t265 * t509, -pkin(5) * t282 - t230 * t506 + t258 * t509, -pkin(5) * t247 - t222 * t506 + t229 * t509, -t275 * t506 + t343 * t509, -t254 * t506 + t293 * t509, -t263 * t506 + t327 * t509, -t274 * t506 + t341 * t509, -t264 * t506 + t328 * t509, -t286 * t506 + t358 * t509, -pkin(5) * t242 - t215 * t506 + t225 * t509, -pkin(5) * t245 - t216 * t506 + t226 * t509, -pkin(5) * t238 - t214 * t506 + t220 * t509, -pkin(5) * t218 - t212 * t506 + t213 * t509; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t476, 0, t475, 0, t521, -t522, t429, pkin(5) * t429, 0, -t476, -t475, 0, 0, 0, t405, -t521, t522, pkin(5) * t405 + (pkin(1) * t509 + qJ(2) * t506) * g(3), t445 * t509 + t525, t422 * t509 + t478 * t506, t437 * t509 + t506 * t539, t444 * t509 - t525, t435 * t509 - t505 * t540, t468 * t509 + t490, pkin(5) * t406 + t361 * t509 + t365 * t506, pkin(5) * t407 + t360 * t509 + t366 * t506, -pkin(2) * t561 + pkin(5) * t426 + t370 * t509, pkin(5) * t356 + t320 * t509 + t337 * t506, t340 * t509 + t534, t297 * t509 + t417 * t506, t334 * t509 + t399 * t506, t339 * t509 - t534, t333 * t509 + t506 * t515, t357 * t509 + t490, pkin(5) * t301 + t244 * t509 + t257 * t506, pkin(5) * t307 + t251 * t509 + t265 * t506, pkin(5) * t283 + t230 * t509 + t258 * t506, pkin(5) * t248 + t222 * t509 + t229 * t506, t275 * t509 + t343 * t506, t254 * t509 + t293 * t506, t263 * t509 + t327 * t506, t274 * t509 + t341 * t506, t264 * t509 + t328 * t506, t286 * t509 + t358 * t506, pkin(5) * t243 + t215 * t509 + t225 * t506, pkin(5) * t246 + t216 * t509 + t226 * t506, pkin(5) * t239 + t214 * t509 + t220 * t506, pkin(5) * t219 + t212 * t509 + t213 * t506; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t481, t482, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t520 - 0.2e1 * t578, -t482 + t495 + 0.2e1 * t497, pkin(1) * t456 + qJ(2) * t450, (t472 - t530) * t508, -t470 * t508 - t473 * t505, -t486 * t505 + t556, (-t471 + t529) * t505, t484 * t508 - t559, 0, qJ(2) * t470 - t580 * t436 - t567, qJ(2) * t473 - t580 * t438 - t566, -qJ(2) * t477 + t580 * t474 - t376, -qJ(2) * t443 - t580 * t376, -t390 * t505 + t391 * t508, -t351 * t505 + t353 * t508, -t380 * t505 + t383 * t508, -t388 * t505 + t389 * t508, -t379 * t505 + t382 * t508, -t408 * t505 + t409 * t508, qJ(2) * t395 - t305 * t505 - t580 * t313 + t326 * t508, qJ(2) * t398 - t308 * t505 - t580 * t335 + t338 * t508, qJ(2) * t393 - t266 * t505 + t267 * t508 - t580 * t298, -qJ(2) * t387 - qJ(4) * t576 - t580 * t252 - t273 * t505, -t310 * t505 + t312 * t508, -t280 * t505 + t281 * t508, -t289 * t505 + t291 * t508, -t309 * t505 + t311 * t508, -t290 * t505 + t292 * t508, -t331 * t505 + t332 * t508, qJ(2) * t318 - t231 * t505 + t233 * t508 - t580 * t255, qJ(2) * t321 - t232 * t505 + t234 * t508 - t580 * t259, qJ(2) * t294 - t227 * t505 + t228 * t508 - t580 * t249, qJ(2) * t240 - t217 * t505 + t221 * t508 - t580 * t223;];
tauB_reg = t1;
