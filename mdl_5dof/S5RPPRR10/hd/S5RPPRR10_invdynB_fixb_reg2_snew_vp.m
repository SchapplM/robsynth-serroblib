% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RPPRR10_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynB_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:24
% EndTime: 2019-12-31 18:04:32
% DurationCPUTime: 7.53s
% Computational Cost: add. (23650->485), mult. (58087->707), div. (0->0), fcn. (40609->8), ass. (0->330)
t505 = sin(qJ(4));
t502 = sin(pkin(8));
t503 = cos(pkin(8));
t508 = cos(qJ(4));
t523 = t502 * t505 + t503 * t508;
t463 = t523 * qJD(1);
t548 = qJD(1) * t503;
t549 = qJD(1) * t502;
t465 = -t505 * t548 + t508 * t549;
t572 = t465 * t463;
t588 = qJDD(4) - t572;
t590 = t505 * t588;
t589 = t508 * t588;
t500 = t503 ^ 2;
t511 = qJD(1) ^ 2;
t512 = t502 ^ 2;
t477 = (t500 + t512) * t511;
t468 = t502 * t477;
t506 = sin(qJ(1));
t509 = cos(qJ(1));
t541 = t509 * qJDD(1);
t437 = -t506 * t468 + t502 * t541;
t587 = pkin(5) * t437;
t542 = t506 * qJDD(1);
t439 = t509 * t468 + t502 * t542;
t586 = pkin(5) * t439;
t504 = sin(qJ(5));
t507 = cos(qJ(5));
t412 = t507 * t463 + t504 * t465;
t414 = -t504 * t463 + t507 * t465;
t367 = t414 * t412;
t499 = qJDD(4) + qJDD(5);
t583 = -t367 + t499;
t585 = t504 * t583;
t584 = t507 * t583;
t543 = qJDD(1) * t503;
t544 = qJDD(1) * t502;
t462 = -t505 * t543 + t508 * t544;
t547 = t463 * qJD(4);
t426 = t462 - t547;
t519 = t523 * qJDD(1);
t546 = t465 * qJD(4);
t516 = t519 + t546;
t354 = -t412 * qJD(5) + t507 * t426 - t504 * t516;
t501 = qJD(4) + qJD(5);
t407 = t501 * t412;
t582 = -t407 + t354;
t483 = t509 * g(1) + t506 * g(2);
t467 = -t511 * pkin(1) + qJDD(1) * qJ(2) - t483;
t545 = qJD(1) * qJD(2);
t581 = (-t467 - 0.2e1 * t545) * t503;
t410 = t412 ^ 2;
t411 = t414 ^ 2;
t580 = t463 ^ 2;
t460 = t465 ^ 2;
t498 = t501 ^ 2;
t579 = pkin(2) + pkin(3);
t578 = pkin(2) * t503;
t493 = t512 * qJDD(1);
t494 = t500 * qJDD(1);
t475 = t494 + t493;
t433 = t506 * t475 + t509 * t477;
t577 = pkin(5) * t433;
t469 = t503 * t477;
t532 = t503 * t541;
t438 = -t506 * t469 + t532;
t576 = pkin(5) * t438;
t575 = qJ(2) * t468;
t574 = qJ(3) * t502;
t573 = qJDD(1) * pkin(1);
t571 = t500 * t511;
t570 = t501 * t504;
t569 = t501 * t507;
t568 = t502 * t503;
t567 = t503 * t511;
t482 = t506 * g(1) - t509 * g(2);
t518 = t511 * qJ(2) - qJDD(2) + t482;
t515 = 0.2e1 * qJD(3) * t549 + t518;
t531 = pkin(1) + t574;
t553 = t511 * t512;
t396 = (t579 * t503 + t531) * qJDD(1) + t515 + (-t553 - t571) * pkin(6);
t446 = qJD(4) * pkin(4) - t465 * pkin(7);
t336 = t516 * pkin(4) - t580 * pkin(7) + t465 * t446 + t396;
t566 = t504 * t336;
t364 = t367 + t499;
t565 = t504 * t364;
t535 = t502 * t545;
t487 = 0.2e1 * t535;
t527 = -t574 - t578;
t472 = t527 * qJD(1);
t552 = t503 * g(3) + t502 * t467;
t529 = t472 * t549 + qJDD(3) + t552;
t522 = t487 + t529;
t380 = (-pkin(3) * t567 - pkin(6) * qJDD(1)) * t502 + t522;
t431 = -t502 * g(3) - t581;
t456 = t472 * t548;
t400 = t431 + t456;
t391 = -pkin(3) * t571 - pkin(6) * t543 + t400;
t343 = -t508 * t380 + t505 * t391;
t301 = (-t426 - t547) * pkin(7) + t588 * pkin(4) - t343;
t344 = t505 * t380 + t508 * t391;
t305 = -t580 * pkin(4) - t516 * pkin(7) - qJD(4) * t446 + t344;
t273 = -t507 * t301 + t504 * t305;
t274 = t504 * t301 + t507 * t305;
t243 = -t507 * t273 + t504 * t274;
t564 = t505 * t243;
t563 = t505 * t396;
t422 = qJDD(4) + t572;
t562 = t505 * t422;
t458 = t518 + t573;
t561 = t506 * t458;
t560 = t506 * t511;
t559 = t507 * t336;
t558 = t507 * t364;
t557 = t508 * t243;
t556 = t508 * t396;
t555 = t508 * t422;
t554 = t509 * t458;
t551 = pkin(1) * t477 + qJ(2) * t475;
t550 = pkin(1) * t543 - qJ(2) * t469;
t540 = t506 * t367;
t539 = t506 * t572;
t538 = t509 * t367;
t537 = t509 * t572;
t536 = pkin(1) + t578;
t533 = t502 * t543;
t244 = t504 * t273 + t507 * t274;
t530 = t504 * t426 + t507 * t516;
t430 = t487 + t552;
t376 = t502 * t430 + t503 * t431;
t443 = -t506 * t482 - t509 * t483;
t481 = t541 - t560;
t528 = -pkin(5) * t481 - t506 * g(3);
t526 = pkin(2) * t502 - qJ(3) * t503;
t298 = -t508 * t343 + t505 * t344;
t299 = t505 * t343 + t508 * t344;
t375 = t503 * t430 - t502 * t431;
t476 = t494 - t493;
t478 = (-t500 + t512) * t511;
t525 = t509 * t476 + t506 * t478;
t524 = t506 * t476 - t509 * t478;
t442 = t509 * t482 - t506 * t483;
t480 = t509 * t511 + t542;
t521 = pkin(1) - t527;
t520 = t509 * t469 + t503 * t542;
t517 = (-qJD(5) + t501) * t414 - t530;
t510 = qJD(4) ^ 2;
t470 = t526 * qJDD(1);
t459 = -pkin(5) * t480 + t509 * g(3);
t449 = -t460 - t510;
t448 = -t460 + t510;
t447 = -t510 + t580;
t445 = t502 * t532 - t560 * t568;
t444 = t480 * t568;
t435 = pkin(5) * t520;
t434 = t509 * t475 - t506 * t477;
t432 = pkin(5) * t434;
t428 = t460 - t580;
t425 = t462 - 0.2e1 * t547;
t424 = 0.2e1 * t546 + t519;
t420 = -t510 - t580;
t418 = t521 * qJDD(1) + t515;
t417 = (t531 + 0.2e1 * t578) * qJDD(1) + t515;
t416 = (t536 + 0.2e1 * t574) * qJDD(1) + t515;
t409 = (-t508 * t463 + t505 * t465) * qJD(4);
t408 = (t505 * t463 + t508 * t465) * qJD(4);
t405 = -t411 + t498;
t404 = t410 - t498;
t403 = -pkin(2) * t493 + t503 * t416;
t402 = qJ(3) * t494 - t502 * t417;
t401 = -t411 - t498;
t399 = -t460 - t580;
t398 = -t529 - 0.2e1 * t535;
t395 = t508 * t426 - t505 * t546;
t394 = -t505 * t426 - t508 * t546;
t393 = t505 * t516 + t508 * t547;
t392 = -t505 * t547 + t508 * t516;
t390 = pkin(2) * t477 + t400;
t389 = -t505 * t449 - t555;
t388 = -t505 * t448 + t589;
t387 = t508 * t447 - t562;
t386 = t508 * t449 - t562;
t385 = -t508 * t448 - t590;
t384 = -t505 * t447 - t555;
t383 = qJ(3) * t477 + t522;
t382 = -pkin(2) * t553 - t456 + (qJ(3) * t567 + g(3)) * t502 + t581;
t379 = (-pkin(2) * t568 + qJ(3) * t500) * t511 + t522;
t373 = -t508 * t424 - t505 * t425;
t372 = t505 * t462 - t508 * t519;
t371 = t505 * t424 - t508 * t425;
t370 = -t508 * t462 - t505 * t519;
t369 = t508 * t420 - t590;
t368 = t505 * t420 + t589;
t366 = t411 - t410;
t362 = -t498 - t410;
t361 = -t502 * t408 + t503 * t409;
t360 = (-t412 * t507 + t414 * t504) * t501;
t359 = (-t412 * t504 - t414 * t507) * t501;
t358 = t509 * t376 - t561;
t357 = t506 * t376 + t554;
t356 = -t502 * t398 + t503 * t400;
t355 = t503 * t398 + t502 * t400;
t353 = -t414 * qJD(5) - t530;
t352 = -t410 - t411;
t351 = -t502 * t394 + t503 * t395;
t350 = -t502 * t392 + t503 * t393;
t349 = t502 * t386 + t503 * t389;
t348 = -t502 * t385 + t503 * t388;
t347 = -t502 * t384 + t503 * t387;
t346 = -t503 * t386 + t502 * t389;
t345 = t503 * t383 - t502 * t390;
t342 = t507 * t404 - t565;
t341 = -t504 * t405 + t584;
t340 = t504 * t404 + t558;
t339 = t507 * t405 + t585;
t338 = -t504 * t401 - t558;
t337 = t507 * t401 - t565;
t335 = t509 * t356 - t506 * t418;
t334 = t506 * t356 + t509 * t418;
t333 = -pkin(6) * t386 + qJ(3) * t425 + t556;
t332 = -t502 * t371 + t503 * t373;
t331 = t502 * t370 + t503 * t372;
t330 = -t503 * t370 + t502 * t372;
t329 = t407 + t354;
t324 = (qJD(5) + t501) * t414 + t530;
t323 = t502 * t368 + t503 * t369;
t322 = -t503 * t368 + t502 * t369;
t321 = t507 * t354 - t414 * t570;
t320 = t504 * t354 + t414 * t569;
t319 = -t504 * t353 + t412 * t569;
t318 = t507 * t353 + t412 * t570;
t317 = t509 * t349 - t506 * t425;
t316 = t506 * t349 + t509 * t425;
t315 = -pkin(6) * t368 + qJ(3) * t424 + t563;
t314 = t507 * t362 - t585;
t313 = t504 * t362 + t584;
t312 = -qJ(2) * t355 - t526 * t418;
t311 = -t505 * t359 + t508 * t360;
t310 = -t508 * t359 - t505 * t360;
t309 = -pkin(6) * t389 + t579 * t425 - t563;
t308 = t509 * t323 - t506 * t424;
t307 = t506 * t323 + t509 * t424;
t306 = -pkin(6) * t369 + t579 * t424 + t556;
t304 = -pkin(1) * t355 - pkin(2) * t398 - qJ(3) * t400;
t303 = t509 * t331 - t506 * t399;
t302 = t506 * t331 + t509 * t399;
t297 = -t505 * t340 + t508 * t342;
t296 = -t505 * t339 + t508 * t341;
t295 = -t508 * t340 - t505 * t342;
t294 = -t508 * t339 - t505 * t341;
t293 = -t505 * t337 + t508 * t338;
t292 = t508 * t337 + t505 * t338;
t291 = -pkin(7) * t337 + t559;
t290 = -pkin(7) * t313 + t566;
t289 = t504 * t329 + t507 * t517;
t288 = -t507 * t324 - t504 * t582;
t287 = -t507 * t329 + t504 * t517;
t286 = -t504 * t324 + t507 * t582;
t285 = -pkin(6) * t298 + qJ(3) * t396;
t284 = -t505 * t320 + t508 * t321;
t283 = -t505 * t318 + t508 * t319;
t282 = -t508 * t320 - t505 * t321;
t281 = -t508 * t318 - t505 * t319;
t280 = -t505 * t313 + t508 * t314;
t279 = t508 * t313 + t505 * t314;
t278 = -pkin(6) * t299 + t579 * t396;
t277 = -t502 * t310 + t503 * t311;
t276 = -pkin(1) * t330 - qJ(3) * t372 + t579 * t370;
t275 = -pkin(6) * t370 + qJ(3) * t399 - t298;
t271 = -pkin(4) * t582 + pkin(7) * t338 + t566;
t270 = -pkin(6) * t372 + t579 * t399 - t299;
t269 = -pkin(1) * t346 - qJ(3) * t389 + t579 * t386 - t344;
t268 = -pkin(4) * t324 + pkin(7) * t314 - t559;
t267 = -qJ(2) * t346 - t502 * t309 + t503 * t333;
t266 = -pkin(1) * t322 - qJ(3) * t369 + t579 * t368 - t343;
t265 = t502 * t298 + t503 * t299;
t264 = -t503 * t298 + t502 * t299;
t263 = -qJ(2) * t322 - t502 * t306 + t503 * t315;
t262 = -t502 * t295 + t503 * t297;
t261 = -t502 * t294 + t503 * t296;
t260 = t502 * t292 + t503 * t293;
t259 = -t503 * t292 + t502 * t293;
t258 = t509 * t265 - t506 * t396;
t257 = t506 * t265 + t509 * t396;
t256 = -t505 * t287 + t508 * t289;
t255 = -t505 * t286 + t508 * t288;
t254 = t508 * t287 + t505 * t289;
t253 = -t508 * t286 - t505 * t288;
t252 = -t502 * t282 + t503 * t284;
t251 = -t502 * t281 + t503 * t283;
t250 = t502 * t279 + t503 * t280;
t249 = -t503 * t279 + t502 * t280;
t248 = t509 * t260 - t506 * t582;
t247 = t506 * t260 + t509 * t582;
t246 = t509 * t250 - t506 * t324;
t245 = t506 * t250 + t509 * t324;
t242 = -qJ(2) * t330 - t502 * t270 + t503 * t275;
t241 = -pkin(4) * t336 + pkin(7) * t244;
t240 = -pkin(6) * t292 + qJ(3) * t582 - t505 * t271 + t508 * t291;
t239 = -pkin(7) * t287 - t243;
t238 = -qJ(2) * t264 - t502 * t278 + t503 * t285;
t237 = t502 * t254 + t503 * t256;
t236 = -t502 * t253 + t503 * t255;
t235 = -t503 * t254 + t502 * t256;
t234 = -pkin(6) * t279 + qJ(3) * t324 - t505 * t268 + t508 * t290;
t233 = -pkin(4) * t352 + pkin(7) * t289 + t244;
t232 = -pkin(6) * t293 - t508 * t271 - t505 * t291 + t579 * t582;
t231 = -pkin(1) * t264 - qJ(3) * t299 + t579 * t298;
t230 = t509 * t237 - t506 * t352;
t229 = t506 * t237 + t509 * t352;
t228 = -pkin(6) * t280 - t508 * t268 - t505 * t290 + t579 * t324;
t227 = t508 * t244 - t564;
t226 = t505 * t244 + t557;
t225 = -pkin(1) * t259 + pkin(4) * t337 - qJ(3) * t293 + t579 * t292 - t274;
t224 = -pkin(1) * t249 + pkin(4) * t313 - qJ(3) * t280 + t579 * t279 - t273;
t223 = -qJ(2) * t259 - t502 * t232 + t503 * t240;
t222 = -pkin(6) * t254 + qJ(3) * t352 - t505 * t233 + t508 * t239;
t221 = -pkin(6) * t256 - t508 * t233 - t505 * t239 + t579 * t352;
t220 = -qJ(2) * t249 - t502 * t228 + t503 * t234;
t219 = t502 * t226 + t503 * t227;
t218 = -t503 * t226 + t502 * t227;
t217 = -pkin(1) * t235 + pkin(4) * t287 - qJ(3) * t256 + t579 * t254;
t216 = t509 * t219 - t506 * t336;
t215 = t506 * t219 + t509 * t336;
t214 = -pkin(6) * t226 - pkin(7) * t557 + qJ(3) * t336 - t505 * t241;
t213 = -pkin(6) * t227 + pkin(7) * t564 - t508 * t241 + t579 * t336;
t212 = -qJ(2) * t235 - t502 * t221 + t503 * t222;
t211 = -pkin(1) * t218 + pkin(4) * t243 - qJ(3) * t227 + t579 * t226;
t210 = -qJ(2) * t218 - t502 * t213 + t503 * t214;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t480, -t481, 0, t443, 0, 0, 0, 0, 0, 0, -t520, t439, t434, t358, 0, 0, 0, 0, 0, 0, -t520, t434, -t439, t335, 0, 0, 0, 0, 0, 0, t308, t317, t303, t258, 0, 0, 0, 0, 0, 0, t246, t248, t230, t216; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t481, -t480, 0, t442, 0, 0, 0, 0, 0, 0, t438, -t437, t433, t357, 0, 0, 0, 0, 0, 0, t438, t433, t437, t334, 0, 0, 0, 0, 0, 0, t307, t316, t302, t257, 0, 0, 0, 0, 0, 0, t245, t247, t229, t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t375, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355, 0, 0, 0, 0, 0, 0, t322, t346, t330, t264, 0, 0, 0, 0, 0, 0, t249, t259, t235, t218; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t481, 0, -t480, 0, t528, -t459, -t442, -pkin(5) * t442, t445, t525, t439, -t445, t520, 0, -t506 * t430 - t502 * t554 - t576, -t506 * t431 - t503 * t554 + t587, t509 * t375 - t577, -pkin(5) * t357 - (pkin(1) * t506 - qJ(2) * t509) * t375, t445, t439, -t525, 0, -t520, -t445, -t506 * t379 + t509 * t402 - t576, t509 * t345 - t506 * t470 - t577, -t506 * t382 + t509 * t403 - t587, -pkin(5) * t334 - t506 * t304 + t509 * t312, t509 * t351 - t539, t509 * t332 - t506 * t428, t509 * t348 - t506 * t462, t509 * t350 + t539, t509 * t347 + t506 * t519, -t506 * qJDD(4) + t509 * t361, -pkin(5) * t307 + t509 * t263 - t506 * t266, -pkin(5) * t316 + t509 * t267 - t506 * t269, -pkin(5) * t302 + t509 * t242 - t506 * t276, -pkin(5) * t257 - t506 * t231 + t509 * t238, t509 * t252 - t540, t509 * t236 - t506 * t366, t509 * t261 - t506 * t329, t509 * t251 + t540, t509 * t262 - t506 * t517, t509 * t277 - t506 * t499, -pkin(5) * t245 + t509 * t220 - t506 * t224, -pkin(5) * t247 + t509 * t223 - t506 * t225, -pkin(5) * t229 + t509 * t212 - t506 * t217, -pkin(5) * t215 + t509 * t210 - t506 * t211; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t480, 0, t481, 0, t459, t528, t443, pkin(5) * t443, t444, t524, -t437, -t444, -t438, 0, t509 * t430 - t502 * t561 - t435, t509 * t431 - t503 * t561 + t586, t506 * t375 + t432, pkin(5) * t358 - (-pkin(1) * t509 - qJ(2) * t506) * t375, t444, -t437, -t524, 0, t438, -t444, t509 * t379 + t506 * t402 - t435, t506 * t345 + t509 * t470 + t432, t509 * t382 + t506 * t403 - t586, pkin(5) * t335 + t509 * t304 + t506 * t312, t506 * t351 + t537, t506 * t332 + t509 * t428, t506 * t348 + t509 * t462, t506 * t350 - t537, t506 * t347 - t509 * t519, t509 * qJDD(4) + t506 * t361, pkin(5) * t308 + t506 * t263 + t509 * t266, pkin(5) * t317 + t506 * t267 + t509 * t269, pkin(5) * t303 + t506 * t242 + t509 * t276, pkin(5) * t258 + t509 * t231 + t506 * t238, t506 * t252 + t538, t506 * t236 + t509 * t366, t506 * t261 + t509 * t329, t506 * t251 - t538, t506 * t262 + t509 * t517, t506 * t277 + t509 * t499, pkin(5) * t246 + t506 * t220 + t509 * t224, pkin(5) * t248 + t506 * t223 + t509 * t225, pkin(5) * t230 + t506 * t212 + t509 * t217, pkin(5) * t216 + t506 * t210 + t509 * t211; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t482, t483, 0, 0, t493, 0.2e1 * t533, 0, t494, 0, 0, t503 * t458 + t550, t575 + (-t458 - t573) * t502, t376 + t551, pkin(1) * t458 + qJ(2) * t376, t493, 0, -0.2e1 * t533, 0, 0, t494, (qJ(3) * t544 + t417) * t503 + t550, t502 * t383 + t503 * t390 + t551, -t575 + (qJDD(1) * t536 + t416) * t502, qJ(2) * t356 + t418 * t521, t503 * t394 + t502 * t395, t503 * t371 + t502 * t373, t503 * t385 + t502 * t388, t503 * t392 + t502 * t393, t503 * t384 + t502 * t387, t503 * t408 + t502 * t409, pkin(1) * t424 + qJ(2) * t323 + t503 * t306 + t502 * t315, pkin(1) * t425 + qJ(2) * t349 + t503 * t309 + t502 * t333, pkin(1) * t399 + qJ(2) * t331 + t503 * t270 + t502 * t275, pkin(1) * t396 + qJ(2) * t265 + t503 * t278 + t502 * t285, t503 * t282 + t502 * t284, t503 * t253 + t502 * t255, t503 * t294 + t502 * t296, t503 * t281 + t502 * t283, t503 * t295 + t502 * t297, t503 * t310 + t502 * t311, pkin(1) * t324 + qJ(2) * t250 + t503 * t228 + t502 * t234, pkin(1) * t582 + qJ(2) * t260 + t503 * t232 + t502 * t240, pkin(1) * t352 + qJ(2) * t237 + t503 * t221 + t502 * t222, pkin(1) * t336 + qJ(2) * t219 + t503 * t213 + t502 * t214;];
tauB_reg = t1;
