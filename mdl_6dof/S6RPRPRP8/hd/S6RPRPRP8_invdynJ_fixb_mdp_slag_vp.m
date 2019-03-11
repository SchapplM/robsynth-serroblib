% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:14
% EndTime: 2019-03-09 03:26:22
% DurationCPUTime: 5.29s
% Computational Cost: add. (5058->492), mult. (9858->610), div. (0->0), fcn. (6630->10), ass. (0->198)
t444 = sin(qJ(1));
t434 = g(1) * t444;
t447 = cos(qJ(1));
t574 = -g(2) * t447 + t434;
t558 = sin(pkin(9));
t495 = qJD(3) * t558;
t559 = cos(pkin(9));
t577 = qJD(1) * t495 - qJDD(1) * t559;
t496 = qJD(3) * t559;
t576 = qJD(1) * t496 + qJDD(1) * t558;
t450 = qJD(1) ^ 2;
t459 = -qJ(2) * t450 - t574;
t443 = sin(qJ(3));
t446 = cos(qJ(3));
t457 = t443 * t559 + t446 * t558;
t379 = t457 * qJD(1);
t575 = qJD(5) + t379;
t386 = -t558 * t443 + t559 * t446;
t431 = t443 * pkin(3);
t537 = qJ(2) + t431;
t352 = pkin(4) * t457 - pkin(8) * t386 + t537;
t448 = -pkin(1) - pkin(7);
t536 = qJ(4) - t448;
t392 = t536 * t443;
t393 = t536 * t446;
t356 = -t392 * t559 - t393 * t558;
t442 = sin(qJ(5));
t445 = cos(qJ(5));
t530 = t442 * t352 + t445 * t356;
t507 = MDP(21) + MDP(23);
t573 = MDP(22) - MDP(25);
t436 = qJDD(1) * qJ(2);
t483 = g(1) * t447 + g(2) * t444;
t437 = qJD(1) * qJD(2);
t504 = 0.2e1 * t437;
t572 = 0.2e1 * t436 + t504 - t483;
t400 = qJD(1) * t448 + qJD(2);
t523 = qJD(1) * t443;
t371 = -qJ(4) * t523 + t400 * t443;
t365 = t558 * t371;
t522 = qJD(1) * t446;
t372 = -qJ(4) * t522 + t446 * t400;
t368 = qJD(3) * pkin(3) + t372;
t332 = t368 * t559 - t365;
t325 = -qJD(3) * pkin(4) - t332;
t382 = t386 * qJD(1);
t516 = t445 * qJD(3);
t358 = t382 * t442 - t516;
t360 = qJD(3) * t442 + t382 * t445;
t306 = t358 * pkin(5) - t360 * qJ(6) + t325;
t353 = t443 * t577 - t446 * t576;
t351 = -qJDD(5) + t353;
t414 = pkin(3) * t558 + pkin(8);
t542 = t414 * t351;
t571 = t306 * t575 + t542;
t435 = qJ(3) + pkin(9);
t424 = cos(t435);
t519 = qJD(5) * t414;
t502 = t575 * t519;
t423 = sin(t435);
t563 = g(3) * t423;
t570 = -t424 * t574 - t502 + t563;
t499 = t559 * t371;
t333 = t558 * t368 + t499;
t326 = qJD(3) * pkin(8) + t333;
t391 = pkin(3) * t523 + qJD(1) * qJ(2) + qJD(4);
t334 = pkin(4) * t379 - pkin(8) * t382 + t391;
t304 = -t326 * t442 + t334 * t445;
t513 = qJD(6) - t304;
t295 = -pkin(5) * t575 + t513;
t305 = t326 * t445 + t334 * t442;
t296 = qJ(6) * t575 + t305;
t476 = t295 * t445 - t296 * t442;
t545 = t360 * t445;
t549 = t358 * t442;
t568 = -(t545 + t549) * MDP(24) - t476 * MDP(26);
t567 = t360 ^ 2;
t565 = pkin(5) * t351;
t564 = pkin(8) * t424;
t562 = g(3) * t424;
t561 = g(3) * t443;
t560 = g(3) * t445;
t557 = pkin(1) * qJDD(1);
t555 = qJ(6) * t351;
t399 = qJDD(1) * t448 + qJDD(2);
t388 = t446 * t399;
t509 = qJDD(1) * t446;
t511 = qJD(1) * qJD(4);
t512 = qJD(1) * qJD(3);
t521 = qJD(3) * t443;
t330 = -t446 * t511 - t400 * t521 + qJDD(3) * pkin(3) + t388 + (t443 * t512 - t509) * qJ(4);
t520 = qJD(3) * t446;
t338 = (-qJ(4) * qJD(1) + t400) * t520 + (-qJ(4) * qJDD(1) + t399 - t511) * t443;
t307 = t330 * t559 - t558 * t338;
t302 = -qJDD(3) * pkin(4) - t307;
t456 = t443 * t576 + t446 * t577;
t518 = qJD(5) * t442;
t317 = -qJD(5) * t516 - t442 * qJDD(3) + t382 * t518 + t445 * t456;
t454 = -t445 * qJDD(3) - t442 * t456;
t318 = qJD(5) * t360 + t454;
t291 = t318 * pkin(5) + t317 * qJ(6) - t360 * qJD(6) + t302;
t554 = t291 * t386;
t553 = t305 * t575;
t552 = t317 * t442;
t551 = t318 * t445;
t550 = t358 * t379;
t548 = t358 * t445;
t547 = t360 * t358;
t546 = t360 * t442;
t544 = t575 * t379;
t543 = t386 * t445;
t346 = t442 * t351;
t541 = t442 * t444;
t540 = t442 * t447;
t539 = t444 * t445;
t347 = t445 * t351;
t538 = t447 * t445;
t517 = qJD(5) * t445;
t534 = -t442 * t318 - t358 * t517;
t308 = t558 * t330 + t559 * t338;
t341 = t372 * t559 - t365;
t343 = pkin(3) * t522 + pkin(4) * t382 + pkin(8) * t379;
t533 = t445 * t341 + t442 * t343;
t532 = t517 * t575 - t346;
t531 = -t442 * t544 - t347;
t340 = t372 * t558 + t499;
t478 = pkin(5) * t442 - qJ(6) * t445;
t529 = -qJD(6) * t442 + t478 * t575 - t340;
t528 = g(2) * t424 * t538 + t423 * t560;
t527 = t447 * pkin(1) + t444 * qJ(2);
t440 = t446 ^ 2;
t525 = t443 ^ 2 - t440;
t449 = qJD(3) ^ 2;
t524 = -t449 - t450;
t515 = pkin(3) * t520 + qJD(2);
t510 = qJDD(1) * t443;
t508 = qJDD(3) * t443;
t505 = t424 * t434;
t503 = t559 * pkin(3);
t501 = t446 * t512;
t500 = -pkin(1) * t444 + t447 * qJ(2);
t303 = qJDD(3) * pkin(8) + t308;
t472 = qJDD(4) + t436 + t437 + (t501 + t510) * pkin(3);
t311 = -t353 * pkin(4) + pkin(8) * t456 + t472;
t462 = t445 * t303 + t442 * t311 - t326 * t518 + t334 * t517;
t289 = qJD(6) * t575 + t462 - t555;
t494 = t295 * t379 + t289;
t488 = t442 * t303 - t445 * t311 + t326 * t517 + t334 * t518;
t290 = qJDD(6) + t488 + t565;
t493 = -t296 * t379 + t290;
t490 = t360 * t575;
t489 = qJDD(2) - t557;
t486 = pkin(4) * t423 - t564;
t374 = t423 * t541 - t538;
t376 = t423 * t540 + t539;
t485 = g(1) * t376 + g(2) * t374;
t375 = t423 * t539 + t540;
t377 = t423 * t538 - t541;
t484 = -g(1) * t377 - g(2) * t375;
t479 = t445 * pkin(5) + t442 * qJ(6);
t369 = -qJD(4) * t446 + t521 * t536;
t370 = -qJD(3) * t393 - qJD(4) * t443;
t336 = -t559 * t369 + t370 * t558;
t355 = -t392 * t558 + t559 * t393;
t477 = -t289 * t445 - t290 * t442;
t475 = -t295 * t442 - t296 * t445;
t441 = -qJ(4) - pkin(7);
t469 = t447 * t431 + t444 * t441 + t500;
t468 = pkin(4) + t479;
t467 = t444 * t431 - t441 * t447 + t527;
t381 = -t443 * t496 - t446 * t495;
t465 = t381 * t442 + t386 * t517;
t464 = -t381 * t445 + t386 * t518;
t463 = 0.2e1 * qJ(2) * t512 + qJDD(3) * t448;
t337 = t369 * t558 + t370 * t559;
t380 = t443 * t495 - t446 * t496;
t342 = -pkin(4) * t380 - pkin(8) * t381 + t515;
t461 = t445 * t337 + t442 * t342 + t352 * t517 - t356 * t518;
t460 = t325 * t575 + t542;
t458 = t307 * t386 + t332 * t381 - t574;
t455 = g(1) * t374 - g(2) * t376 + t442 * t562 - t488;
t453 = -t448 * t449 + t572;
t452 = t306 * t360 + qJDD(6) - t455;
t451 = -g(1) * t375 + g(2) * t377 - t424 * t560 + t462;
t427 = qJDD(3) * t446;
t415 = -t503 - pkin(4);
t384 = -t503 - t468;
t321 = pkin(5) * t360 + qJ(6) * t358;
t319 = t386 * t478 + t355;
t313 = -pkin(5) * t457 - t352 * t445 + t356 * t442;
t312 = qJ(6) * t457 + t530;
t299 = -pkin(5) * t382 + t341 * t442 - t343 * t445;
t298 = qJ(6) * t382 + t533;
t297 = t358 * t575 - t317;
t294 = t478 * t381 + (qJD(5) * t479 - qJD(6) * t445) * t386 + t336;
t293 = pkin(5) * t380 + qJD(5) * t530 + t337 * t442 - t342 * t445;
t292 = -qJ(6) * t380 + qJD(6) * t457 + t461;
t1 = [t574 * MDP(2) + (qJDD(2) - t574 - 0.2e1 * t557) * MDP(4) + (-t308 * t457 + t333 * t380 + t336 * t382 - t337 * t379 + t356 * t353 - t355 * t456 - t458) * MDP(14) + (qJDD(1) * t440 - 0.2e1 * t443 * t501) * MDP(7) + (-t292 * t358 + t293 * t360 - t312 * t318 - t313 * t317 + t483 * t424 + t476 * t381 + (qJD(5) * t475 - t289 * t442 + t290 * t445) * t386) * MDP(24) + t483 * MDP(3) + (t289 * t312 + t296 * t292 + t291 * t319 + t306 * t294 + t290 * t313 + t295 * t293 - g(1) * (pkin(5) * t377 + qJ(6) * t376 + t447 * t486 + t469) - g(2) * (pkin(5) * t375 + qJ(6) * t374 + t444 * t486 + t467)) * MDP(26) + t572 * MDP(5) + (t443 * t453 + t446 * t463) * MDP(12) + (-t443 * t463 + t446 * t453) * MDP(13) + (-t443 * t449 + t427) * MDP(9) + (t302 * t543 + t305 * t380 - t355 * t317 - t325 * t464 + t336 * t360 + t351 * t530 - t457 * t462 - t461 * t575 + t485) * MDP(22) + (-t290 * t457 - t293 * t575 + t294 * t358 + t295 * t380 + t306 * t465 + t313 * t351 + t318 * t319 + t442 * t554 + t484) * MDP(23) + (-t488 * t457 - t304 * t380 + t336 * t358 + t355 * t318 + ((-qJD(5) * t356 + t342) * t575 - t352 * t351 + t325 * qJD(5) * t386) * t445 + ((-qJD(5) * t352 - t337) * t575 + t356 * t351 + t302 * t386 + t325 * t381) * t442 + t484) * MDP(21) + (-t351 * t457 - t380 * t575) * MDP(20) + (-t317 * t457 - t347 * t386 - t360 * t380 - t464 * t575) * MDP(18) + (-t318 * t457 + t346 * t386 + t358 * t380 - t465 * t575) * MDP(19) + (t289 * t457 - t291 * t543 + t292 * t575 - t294 * t360 - t296 * t380 + t306 * t464 - t312 * t351 + t317 * t319 - t485) * MDP(25) + qJDD(1) * MDP(1) + (-t446 * t449 - t508) * MDP(10) + 0.2e1 * (-t443 * t509 + t512 * t525) * MDP(8) + (-t489 * pkin(1) - g(1) * t500 - g(2) * t527 + (t504 + t436) * qJ(2)) * MDP(6) + (-g(1) * t469 - g(2) * t467 - t307 * t355 + t308 * t356 - t332 * t336 + t333 * t337 + t391 * t515 + t472 * t537) * MDP(15) + (-t317 * t543 - t360 * t464) * MDP(16) + ((-t546 - t548) * t381 + (t552 - t551 + (-t545 + t549) * qJD(5)) * t386) * MDP(17); qJDD(1) * MDP(4) - t450 * MDP(5) + (t489 + t459) * MDP(6) + (t443 * t524 + t427) * MDP(12) + (t446 * t524 - t508) * MDP(13) + (-t381 * t382 + t386 * t456) * MDP(14) + t458 * MDP(15) + (-t306 * t381 - t554 - t574) * MDP(26) + (MDP(14) * t379 - t333 * MDP(15) + (-t546 + t548) * MDP(24) + t475 * MDP(26)) * t380 + (-t391 * MDP(15) - t568) * qJD(1) + (t507 * (-qJD(1) * t445 + t380 * t442) + t573 * (qJD(1) * t442 + t380 * t445)) * t575 - (-t353 * MDP(14) - t308 * MDP(15) + (t551 + t552) * MDP(24) + t477 * MDP(26) + (-t442 * t507 - t445 * t573) * t351 + ((-t442 * t573 + t445 * t507) * t575 + t568) * qJD(5)) * t457 + t507 * (-t386 * t318 - t381 * t358) + t573 * (t317 * t386 - t360 * t381); MDP(9) * t509 - MDP(10) * t510 + qJDD(3) * MDP(11) + (t446 * t459 + t388 + t561) * MDP(12) + (g(3) * t446 + (-t399 - t459) * t443) * MDP(13) + ((t333 - t340) * t382 - (-t341 + t332) * t379 + (t353 * t558 + t456 * t559) * pkin(3)) * MDP(14) + (t332 * t340 - t333 * t341 + (t558 * t308 + t559 * t307 + t561 + (-qJD(1) * t391 - t574) * t446) * pkin(3)) * MDP(15) + (t445 * t490 - t552) * MDP(16) + ((-t317 - t550) * t445 - t575 * t546 + t534) * MDP(17) + (-t360 * t382 + t445 * t544 + t532) * MDP(18) + (t358 * t382 - t518 * t575 + t531) * MDP(19) - t575 * t382 * MDP(20) + (-t304 * t382 + t415 * t318 - t340 * t358 + (-t505 - t302 + (-t343 - t519) * t575) * t445 + (t341 * t575 + t460) * t442 + t528) * MDP(21) + (-t415 * t317 + t533 * t575 + t305 * t382 - t340 * t360 + t460 * t445 + (t302 - t570) * t442) * MDP(22) + (t295 * t382 + t299 * t575 + t318 * t384 + t529 * t358 + (-t291 - t502 - t505) * t445 + t571 * t442 + t528) * MDP(23) + (-t562 + t298 * t358 - t299 * t360 - t574 * t423 + (-t318 * t414 + (t360 * t414 + t295) * qJD(5) + t494) * t445 + (-t317 * t414 + (t358 * t414 - t296) * qJD(5) + t493) * t442) * MDP(24) + (-t296 * t382 - t298 * t575 + t317 * t384 - t529 * t360 - t571 * t445 + (-t291 + t570) * t442) * MDP(25) + (t291 * t384 - t296 * t298 - t295 * t299 - g(3) * (-t431 + t564) + t468 * t563 + t529 * t306 + (qJD(5) * t476 - t477) * t414 - t574 * (pkin(3) * t446 + pkin(8) * t423 + t424 * t468)) * MDP(26) + (MDP(7) * t443 * t446 - MDP(8) * t525) * t450; -t379 ^ 2 * MDP(14) + (t333 * t379 + t472 - t483) * MDP(15) + t531 * MDP(21) + t534 * MDP(24) + t532 * MDP(25) - t483 * MDP(26) + (-MDP(14) * t382 + t332 * MDP(15) - t306 * MDP(26) - t358 * t507 - t360 * t573) * t382 + (-t351 * MDP(23) + (t317 - t550) * MDP(24) + (qJD(5) * t296 - t493) * MDP(26) + (-MDP(22) * t575 + t379 * MDP(25)) * t575) * t445 + (t351 * MDP(22) + (qJD(5) * t295 + t494) * MDP(26) + MDP(24) * t490 + (-qJD(5) * MDP(21) - MDP(23) * t575) * t575) * t442; MDP(16) * t547 + (-t358 ^ 2 + t567) * MDP(17) + t297 * MDP(18) + (-t454 + (-qJD(5) + t575) * t360) * MDP(19) - t351 * MDP(20) + (-t325 * t360 + t455 + t553) * MDP(21) + (t304 * t575 + t325 * t358 - t451) * MDP(22) + (-t321 * t358 - t452 + t553 - 0.2e1 * t565) * MDP(23) + (pkin(5) * t317 - qJ(6) * t318 + (t296 - t305) * t360 + (t295 - t513) * t358) * MDP(24) + (-0.2e1 * t555 - t306 * t358 + t321 * t360 + (0.2e1 * qJD(6) - t304) * t575 + t451) * MDP(25) + (t289 * qJ(6) - t290 * pkin(5) - t306 * t321 - t295 * t305 - g(1) * (-pkin(5) * t374 + qJ(6) * t375) - g(2) * (pkin(5) * t376 - qJ(6) * t377) + t478 * t562 + t513 * t296) * MDP(26); (t351 + t547) * MDP(23) + t297 * MDP(24) + (-t575 ^ 2 - t567) * MDP(25) + (-t296 * t575 + t452 + t565) * MDP(26);];
tau  = t1;
