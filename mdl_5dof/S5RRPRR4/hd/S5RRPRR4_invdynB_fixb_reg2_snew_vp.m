% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauB_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S5RRPRR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynB_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynB_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:39
% EndTime: 2022-01-20 10:48:52
% DurationCPUTime: 11.77s
% Computational Cost: add. (45949->511), mult. (61110->765), div. (0->0), fcn. (38139->10), ass. (0->343)
t549 = qJD(1) + qJD(2);
t545 = t549 ^ 2;
t547 = qJDD(1) + qJDD(2);
t553 = sin(pkin(9));
t554 = cos(pkin(9));
t507 = t554 * t545 + t553 * t547;
t510 = t553 * t545 - t554 * t547;
t557 = sin(qJ(2));
t561 = cos(qJ(2));
t450 = t561 * t507 - t557 * t510;
t552 = g(3) - qJDD(3);
t484 = qJ(3) * t507 - t554 * t552;
t619 = qJ(3) * t510 - t553 * t552;
t379 = pkin(6) * t450 + t561 * t484 - t557 * t619;
t454 = t557 * t507 + t561 * t510;
t558 = sin(qJ(1));
t562 = cos(qJ(1));
t402 = t558 * t450 + t562 * t454;
t629 = pkin(6) * t454 + t557 * t484 + t561 * t619;
t638 = pkin(5) * t402 + t558 * t379 + t562 * t629;
t618 = t562 * t450 - t558 * t454;
t637 = pkin(5) * t618 + t562 * t379 - t558 * t629;
t533 = t562 * g(1) + t558 * g(2);
t609 = qJD(1) ^ 2;
t567 = -t609 * pkin(1) - t533;
t532 = t558 * g(1) - t562 * g(2);
t568 = qJDD(1) * pkin(1) + t532;
t465 = t557 * t568 + t561 * t567;
t458 = -t545 * pkin(2) + t465;
t565 = -t557 * t567 + t561 * t568;
t564 = t547 * pkin(2) + t565;
t406 = t553 * t458 - t554 * t564;
t407 = t554 * t458 + t553 * t564;
t576 = t553 * t406 + t554 * t407;
t341 = t554 * t406 - t553 * t407;
t586 = t561 * t341;
t293 = -t557 * t576 + t586;
t593 = t557 * t341;
t624 = t561 * t576 + t593;
t259 = t558 * t293 + t562 * t624;
t634 = t562 * t293 - t558 * t624;
t514 = t561 * t545 + t557 * t547;
t517 = t557 * t545 - t561 * t547;
t461 = t558 * t514 + t562 * t517;
t490 = pkin(6) * t514 - t561 * g(3);
t620 = pkin(6) * t517 - t557 * g(3);
t631 = pkin(5) * t461 + t558 * t490 + t562 * t620;
t569 = t562 * t514 - t558 * t517;
t630 = pkin(5) * t569 + t562 * t490 - t558 * t620;
t575 = t561 * t465 - t557 * t565;
t417 = -t557 * t465 - t561 * t565;
t592 = t558 * t417;
t349 = t562 * t575 + t592;
t585 = t562 * t417;
t623 = -t558 * t575 + t585;
t555 = sin(qJ(5));
t559 = cos(qJ(5));
t560 = cos(qJ(4));
t556 = sin(qJ(4));
t604 = t549 * t556;
t493 = -t559 * t560 * t549 + t555 * t604;
t495 = (t555 * t560 + t556 * t559) * t549;
t448 = t495 * t493;
t581 = qJDD(4) + qJDD(5);
t611 = -t448 + t581;
t622 = t555 * t611;
t621 = t559 * t611;
t548 = qJD(4) + qJD(5);
t486 = t548 * t493;
t583 = qJD(4) * t549;
t577 = t560 * t583;
t594 = t556 * t547;
t503 = t577 + t594;
t536 = t560 * t547;
t578 = t556 * t583;
t504 = t536 - t578;
t566 = t493 * qJD(5) - t559 * t503 - t555 * t504;
t610 = -t486 - t566;
t574 = t555 * t503 - t559 * t504;
t395 = (qJD(5) - t548) * t495 + t574;
t491 = t493 ^ 2;
t492 = t495 ^ 2;
t544 = t548 ^ 2;
t606 = t548 * t555;
t605 = t548 * t559;
t550 = t556 ^ 2;
t603 = t550 * t545;
t551 = t560 ^ 2;
t537 = t551 * t545;
t387 = -t547 * pkin(3) - t545 * pkin(7) + t406;
t523 = qJD(4) * pkin(4) - pkin(8) * t604;
t355 = -t504 * pkin(4) - pkin(8) * t537 + t523 * t604 + t387;
t600 = t555 * t355;
t438 = t448 + t581;
t599 = t555 * t438;
t388 = -t545 * pkin(3) + t547 * pkin(7) + t407;
t370 = t556 * t388 + t560 * t552;
t531 = t560 * t545 * t556;
t521 = qJDD(4) + t531;
t353 = (-t503 + t577) * pkin(8) + t521 * pkin(4) - t370;
t371 = t560 * t388 - t556 * t552;
t354 = -pkin(4) * t537 + t504 * pkin(8) - qJD(4) * t523 + t371;
t303 = -t559 * t353 + t555 * t354;
t305 = t555 * t353 + t559 * t354;
t265 = -t559 * t303 + t555 * t305;
t598 = t556 * t265;
t597 = t556 * t387;
t596 = t556 * t521;
t522 = qJDD(4) - t531;
t595 = t556 * t522;
t591 = t559 * t355;
t590 = t559 * t438;
t589 = t560 * t265;
t588 = t560 * t387;
t587 = t560 * t522;
t584 = t550 + t551;
t580 = t553 * t448;
t579 = t554 * t448;
t267 = t555 * t303 + t559 * t305;
t328 = t556 * t370 + t560 * t371;
t478 = -t558 * t532 - t562 * t533;
t572 = t553 * t531;
t571 = t554 * t531;
t525 = t562 * qJDD(1) - t558 * t609;
t570 = -pkin(5) * t525 - t558 * g(3);
t327 = t560 * t370 - t556 * t371;
t477 = t562 * t532 - t558 * t533;
t563 = qJD(4) ^ 2;
t529 = -t537 - t563;
t528 = t537 - t563;
t527 = -t563 - t603;
t526 = t563 - t603;
t524 = t558 * qJDD(1) + t562 * t609;
t519 = t537 - t603;
t518 = t537 + t603;
t513 = t560 * t521;
t512 = t584 * t547;
t505 = t536 - 0.2e1 * t578;
t502 = 0.2e1 * t577 + t594;
t501 = -pkin(5) * t524 + t562 * g(3);
t500 = t584 * t583;
t480 = -t492 + t544;
t479 = t491 - t544;
t476 = t553 * qJDD(4) + t554 * t500;
t475 = -t554 * qJDD(4) + t553 * t500;
t474 = t560 * t503 - t550 * t583;
t473 = -t556 * t504 - t551 * t583;
t472 = -t492 - t544;
t471 = -t556 * t527 - t587;
t470 = -t556 * t526 + t513;
t469 = t560 * t529 - t596;
t468 = t560 * t528 - t595;
t467 = t560 * t527 - t595;
t466 = t556 * t529 + t513;
t457 = t554 * t512 - t553 * t518;
t456 = t553 * t512 + t554 * t518;
t449 = -t556 * t502 + t560 * t505;
t444 = -t492 + t491;
t443 = t554 * t470 + t553 * t594;
t442 = t554 * t468 + t553 * t536;
t441 = t553 * t470 - t554 * t594;
t440 = t553 * t468 - t554 * t536;
t436 = -t544 - t491;
t435 = t554 * t474 - t572;
t434 = t554 * t473 + t572;
t433 = t553 * t474 + t571;
t432 = t553 * t473 - t571;
t431 = t554 * t471 + t553 * t502;
t430 = t554 * t469 - t553 * t505;
t429 = t553 * t471 - t554 * t502;
t428 = t553 * t469 + t554 * t505;
t427 = (-t493 * t559 + t495 * t555) * t548;
t426 = (-t493 * t555 - t495 * t559) * t548;
t425 = -t491 - t492;
t424 = -t557 * t475 + t561 * t476;
t423 = t561 * t475 + t557 * t476;
t421 = -t495 * qJD(5) - t574;
t420 = t554 * t449 - t553 * t519;
t419 = t553 * t449 + t554 * t519;
t414 = pkin(1) * g(3) + pkin(6) * t575;
t413 = t559 * t479 - t599;
t412 = -t555 * t480 + t621;
t411 = t555 * t479 + t590;
t410 = t559 * t480 + t622;
t409 = -t555 * t472 - t590;
t408 = t559 * t472 - t599;
t405 = -t557 * t456 + t561 * t457;
t404 = t561 * t456 + t557 * t457;
t399 = -t486 + t566;
t394 = (qJD(5) + t548) * t495 + t574;
t392 = -t495 * t606 - t559 * t566;
t391 = t495 * t605 - t555 * t566;
t390 = -t555 * t421 + t493 * t605;
t389 = t559 * t421 + t493 * t606;
t386 = -t557 * t441 + t561 * t443;
t385 = -t557 * t440 + t561 * t442;
t384 = t561 * t441 + t557 * t443;
t383 = t561 * t440 + t557 * t442;
t381 = t559 * t436 - t622;
t380 = t555 * t436 + t621;
t375 = -t557 * t433 + t561 * t435;
t374 = -t557 * t432 + t561 * t434;
t373 = t561 * t433 + t557 * t435;
t372 = t561 * t432 + t557 * t434;
t368 = -t557 * t429 + t561 * t431;
t367 = -t557 * t428 + t561 * t430;
t366 = t561 * t429 + t557 * t431;
t365 = t561 * t428 + t557 * t430;
t364 = -t556 * t426 + t560 * t427;
t363 = t554 * t364 + t553 * t581;
t362 = t553 * t364 - t554 * t581;
t361 = -pkin(7) * t467 + t588;
t360 = -pkin(7) * t466 + t597;
t359 = -t557 * t419 + t561 * t420;
t358 = t561 * t419 + t557 * t420;
t357 = -pkin(3) * t467 + t371;
t356 = -pkin(3) * t466 + t370;
t347 = -t556 * t411 + t560 * t413;
t346 = -t556 * t410 + t560 * t412;
t345 = -t556 * t408 + t560 * t409;
t344 = t560 * t408 + t556 * t409;
t343 = -t558 * t404 + t562 * t405;
t340 = t562 * t404 + t558 * t405;
t337 = -t395 * t559 - t555 * t399;
t336 = -t559 * t394 - t555 * t610;
t335 = -t395 * t555 + t559 * t399;
t334 = -t555 * t394 + t559 * t610;
t333 = pkin(2) * t552 + qJ(3) * t576;
t332 = -t556 * t391 + t560 * t392;
t331 = -t556 * t389 + t560 * t390;
t330 = -t556 * t380 + t560 * t381;
t329 = t560 * t380 + t556 * t381;
t325 = -pkin(8) * t408 + t591;
t324 = -t558 * t366 + t562 * t368;
t323 = -t558 * t365 + t562 * t367;
t322 = t562 * t366 + t558 * t368;
t321 = t562 * t365 + t558 * t367;
t320 = -pkin(8) * t380 + t600;
t319 = t554 * t332 + t580;
t318 = t554 * t331 - t580;
t317 = t553 * t332 - t579;
t316 = t553 * t331 + t579;
t315 = t554 * t347 - t553 * t395;
t314 = t554 * t346 - t553 * t399;
t313 = t553 * t347 + t554 * t395;
t312 = t553 * t346 + t554 * t399;
t311 = -qJ(3) * t456 + t554 * t327;
t310 = qJ(3) * t457 + t553 * t327;
t309 = -t557 * t362 + t561 * t363;
t308 = t554 * t345 + t553 * t610;
t307 = t561 * t362 + t557 * t363;
t306 = t553 * t345 - t554 * t610;
t304 = t554 * t330 + t553 * t394;
t302 = t553 * t330 - t554 * t394;
t300 = -pkin(4) * t610 + pkin(8) * t409 + t600;
t299 = t554 * t328 + t553 * t387;
t298 = t553 * t328 - t554 * t387;
t297 = -qJ(3) * t429 - t553 * t357 + t554 * t361;
t296 = -qJ(3) * t428 - t553 * t356 + t554 * t360;
t295 = -pkin(4) * t394 + pkin(8) * t381 - t591;
t290 = -pkin(2) * t467 + qJ(3) * t431 + t554 * t357 + t553 * t361;
t289 = -pkin(2) * t466 + qJ(3) * t430 + t554 * t356 + t553 * t360;
t288 = -t556 * t335 + t560 * t337;
t287 = -t556 * t334 + t560 * t336;
t286 = t560 * t335 + t556 * t337;
t285 = t554 * t287 - t553 * t444;
t284 = t553 * t287 + t554 * t444;
t283 = t554 * t288 + t553 * t425;
t282 = t553 * t288 - t554 * t425;
t281 = -t557 * t317 + t561 * t319;
t280 = -t557 * t316 + t561 * t318;
t279 = t561 * t317 + t557 * t319;
t278 = t561 * t316 + t557 * t318;
t277 = -t557 * t313 + t561 * t315;
t276 = -t557 * t312 + t561 * t314;
t275 = t561 * t313 + t557 * t315;
t274 = t561 * t312 + t557 * t314;
t273 = -t557 * t306 + t561 * t308;
t272 = t561 * t306 + t557 * t308;
t271 = -pkin(3) * t344 - pkin(4) * t408 + t305;
t270 = -pkin(3) * t286 - pkin(4) * t335;
t269 = -pkin(3) * t329 - pkin(4) * t380 + t303;
t268 = -t557 * t302 + t561 * t304;
t266 = t561 * t302 + t557 * t304;
t264 = -pkin(6) * t404 - t557 * t310 + t561 * t311;
t263 = pkin(6) * t405 + t561 * t310 + t557 * t311;
t262 = -t557 * t298 + t561 * t299;
t261 = t561 * t298 + t557 * t299;
t260 = -pkin(7) * t344 - t556 * t300 + t560 * t325;
t257 = pkin(6) * t293 + qJ(3) * t586 - t557 * t333;
t256 = pkin(1) * t552 + pkin(6) * t624 + qJ(3) * t593 + t561 * t333;
t255 = -pkin(4) * t355 + pkin(8) * t267;
t254 = -qJ(3) * t298 - (pkin(3) * t553 - pkin(7) * t554) * t327;
t253 = -pkin(7) * t329 - t556 * t295 + t560 * t320;
t252 = -pkin(8) * t335 - t265;
t251 = -pkin(6) * t366 - t557 * t290 + t561 * t297;
t250 = -pkin(6) * t365 - t557 * t289 + t561 * t296;
t249 = -pkin(1) * t467 + pkin(6) * t368 + t561 * t290 + t557 * t297;
t248 = -pkin(1) * t466 + pkin(6) * t367 + t561 * t289 + t557 * t296;
t247 = -t557 * t284 + t561 * t285;
t246 = t561 * t284 + t557 * t285;
t245 = -pkin(4) * t425 + pkin(8) * t337 + t267;
t244 = -t557 * t282 + t561 * t283;
t243 = t561 * t282 + t557 * t283;
t242 = qJ(3) * t299 - (-pkin(3) * t554 - pkin(7) * t553 - pkin(2)) * t327;
t241 = -t558 * t272 + t562 * t273;
t240 = t562 * t272 + t558 * t273;
t239 = -t558 * t266 + t562 * t268;
t238 = t560 * t267 - t598;
t237 = t562 * t266 + t558 * t268;
t236 = t556 * t267 + t589;
t235 = -t558 * t261 + t562 * t262;
t234 = t562 * t261 + t558 * t262;
t233 = t554 * t238 + t553 * t355;
t232 = t553 * t238 - t554 * t355;
t231 = -qJ(3) * t306 + t554 * t260 - t553 * t271;
t230 = -qJ(3) * t302 + t554 * t253 - t553 * t269;
t229 = -pkin(2) * t344 + qJ(3) * t308 + t553 * t260 + t554 * t271;
t228 = -t558 * t243 + t562 * t244;
t227 = t562 * t243 + t558 * t244;
t226 = -pkin(3) * t236 - pkin(4) * t265;
t225 = -pkin(2) * t329 + qJ(3) * t304 + t553 * t253 + t554 * t269;
t224 = -pkin(7) * t286 - t556 * t245 + t560 * t252;
t223 = -pkin(6) * t261 - t557 * t242 + t561 * t254;
t222 = -pkin(7) * t236 - pkin(8) * t589 - t556 * t255;
t221 = pkin(1) * t327 + pkin(6) * t262 + t561 * t242 + t557 * t254;
t220 = -t557 * t232 + t561 * t233;
t219 = t561 * t232 + t557 * t233;
t218 = -qJ(3) * t282 + t554 * t224 - t553 * t270;
t217 = -pkin(2) * t286 + qJ(3) * t283 + t553 * t224 + t554 * t270;
t216 = -pkin(6) * t272 - t557 * t229 + t561 * t231;
t215 = -pkin(1) * t344 + pkin(6) * t273 + t561 * t229 + t557 * t231;
t214 = -pkin(6) * t266 - t557 * t225 + t561 * t230;
t213 = -pkin(1) * t329 + pkin(6) * t268 + t561 * t225 + t557 * t230;
t212 = -t558 * t219 + t562 * t220;
t211 = t562 * t219 + t558 * t220;
t210 = -qJ(3) * t232 + t554 * t222 - t553 * t226;
t209 = -pkin(2) * t236 + qJ(3) * t233 + t553 * t222 + t554 * t226;
t208 = -pkin(6) * t243 - t557 * t217 + t561 * t218;
t207 = -pkin(1) * t286 + pkin(6) * t244 + t561 * t217 + t557 * t218;
t206 = -pkin(6) * t219 - t557 * t209 + t561 * t210;
t205 = -pkin(1) * t236 + pkin(6) * t220 + t561 * t209 + t557 * t210;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t524, -t525, 0, t478, 0, 0, 0, 0, 0, 0, -t569, t461, 0, t349, 0, 0, 0, 0, 0, 0, -t618, t402, 0, t259, 0, 0, 0, 0, 0, 0, t323, t324, t343, t235, 0, 0, 0, 0, 0, 0, t239, t241, t228, t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t525, -t524, 0, t477, 0, 0, 0, 0, 0, 0, -t461, -t569, 0, -t623, 0, 0, 0, 0, 0, 0, -t402, -t618, 0, -t634, 0, 0, 0, 0, 0, 0, t321, t322, t340, t234, 0, 0, 0, 0, 0, 0, t237, t240, t227, t211; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t552, 0, 0, 0, 0, 0, 0, t466, t467, 0, -t327, 0, 0, 0, 0, 0, 0, t329, t344, t286, t236; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t525, 0, -t524, 0, t570, -t501, -t477, -pkin(5) * t477, 0, 0, -t461, 0, -t569, 0, t631, t630, t623, pkin(5) * t623 + pkin(6) * t585 - t558 * t414, 0, 0, -t402, 0, -t618, 0, t638, t637, t634, pkin(5) * t634 - t558 * t256 + t562 * t257, -t558 * t373 + t562 * t375, -t558 * t358 + t562 * t359, -t558 * t384 + t562 * t386, -t558 * t372 + t562 * t374, -t558 * t383 + t562 * t385, -t558 * t423 + t562 * t424, -pkin(5) * t321 - t558 * t248 + t562 * t250, -pkin(5) * t322 - t558 * t249 + t562 * t251, -pkin(5) * t340 - t558 * t263 + t562 * t264, -pkin(5) * t234 - t558 * t221 + t562 * t223, -t558 * t279 + t562 * t281, -t558 * t246 + t562 * t247, -t558 * t274 + t562 * t276, -t558 * t278 + t562 * t280, -t558 * t275 + t562 * t277, -t558 * t307 + t562 * t309, -pkin(5) * t237 - t558 * t213 + t562 * t214, -pkin(5) * t240 - t558 * t215 + t562 * t216, -pkin(5) * t227 - t558 * t207 + t562 * t208, -pkin(5) * t211 - t558 * t205 + t562 * t206; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t524, 0, t525, 0, t501, t570, t478, pkin(5) * t478, 0, 0, t569, 0, -t461, 0, -t630, t631, t349, pkin(5) * t349 + pkin(6) * t592 + t562 * t414, 0, 0, t618, 0, -t402, 0, -t637, t638, t259, pkin(5) * t259 + t562 * t256 + t558 * t257, t562 * t373 + t558 * t375, t562 * t358 + t558 * t359, t562 * t384 + t558 * t386, t562 * t372 + t558 * t374, t562 * t383 + t558 * t385, t562 * t423 + t558 * t424, pkin(5) * t323 + t562 * t248 + t558 * t250, pkin(5) * t324 + t562 * t249 + t558 * t251, pkin(5) * t343 + t562 * t263 + t558 * t264, pkin(5) * t235 + t562 * t221 + t558 * t223, t562 * t279 + t558 * t281, t562 * t246 + t558 * t247, t562 * t274 + t558 * t276, t562 * t278 + t558 * t280, t562 * t275 + t558 * t277, t562 * t307 + t558 * t309, pkin(5) * t239 + t562 * t213 + t558 * t214, pkin(5) * t241 + t562 * t215 + t558 * t216, pkin(5) * t228 + t562 * t207 + t558 * t208, pkin(5) * t212 + t562 * t205 + t558 * t206; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t532, t533, 0, 0, 0, 0, 0, 0, 0, t547, -pkin(1) * t517 + t565, -pkin(1) * t514 - t465, 0, -pkin(1) * t417, 0, 0, 0, 0, 0, t547, -pkin(1) * t454 - pkin(2) * t510 - t406, -pkin(1) * t450 - pkin(2) * t507 - t407, 0, -pkin(1) * t293 - pkin(2) * t341, (t503 + t577) * t556, t560 * t502 + t556 * t505, t560 * t526 + t596, (t504 - t578) * t560, t556 * t528 + t587, 0, pkin(1) * t365 + pkin(2) * t428 + pkin(3) * t505 + pkin(7) * t469 - t588, pkin(1) * t366 + pkin(2) * t429 - pkin(3) * t502 + pkin(7) * t471 + t597, pkin(1) * t404 + pkin(2) * t456 + pkin(3) * t518 + pkin(7) * t512 + t328, pkin(1) * t261 + pkin(2) * t298 - pkin(3) * t387 + pkin(7) * t328, t560 * t391 + t556 * t392, t560 * t334 + t556 * t336, t560 * t410 + t556 * t412, t560 * t389 + t556 * t390, t560 * t411 + t556 * t413, t560 * t426 + t556 * t427, pkin(1) * t266 + pkin(2) * t302 - pkin(3) * t394 + pkin(7) * t330 + t560 * t295 + t556 * t320, pkin(1) * t272 + pkin(2) * t306 - pkin(3) * t610 + pkin(7) * t345 + t560 * t300 + t556 * t325, pkin(1) * t243 + pkin(2) * t282 - pkin(3) * t425 + pkin(7) * t288 + t560 * t245 + t556 * t252, pkin(1) * t219 + pkin(2) * t232 - pkin(3) * t355 + pkin(7) * t238 - pkin(8) * t598 + t560 * t255;];
tauB_reg = t1;
