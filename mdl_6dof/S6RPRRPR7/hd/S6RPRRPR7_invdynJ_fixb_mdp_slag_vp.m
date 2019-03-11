% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRPR7_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:21:11
% EndTime: 2019-03-09 05:21:18
% DurationCPUTime: 5.22s
% Computational Cost: add. (6015->447), mult. (12242->573), div. (0->0), fcn. (8728->14), ass. (0->213)
t518 = qJD(1) ^ 2;
t515 = cos(qJ(1));
t495 = g(2) * t515;
t511 = sin(qJ(1));
t496 = g(1) * t511;
t588 = t496 - t495;
t529 = -qJ(2) * t518 - t588;
t505 = qJ(3) + qJ(4);
t489 = sin(t505);
t490 = cos(t505);
t634 = g(3) * t489 - t490 * t588;
t509 = sin(qJ(4));
t513 = cos(qJ(4));
t514 = cos(qJ(3));
t584 = qJD(1) * t514;
t510 = sin(qJ(3));
t585 = qJD(1) * t510;
t434 = t509 * t585 - t513 * t584;
t445 = t509 * t514 + t510 * t513;
t435 = t445 * qJD(1);
t506 = sin(pkin(10));
t507 = cos(pkin(10));
t541 = -t507 * t434 - t435 * t506;
t499 = qJD(3) + qJD(4);
t403 = t499 * t445;
t574 = qJDD(1) * t514;
t575 = qJDD(1) * t510;
t380 = -qJD(1) * t403 - t509 * t575 + t513 * t574;
t581 = qJD(4) * t509;
t568 = t510 * t581;
t576 = qJD(1) * qJD(3);
t628 = t510 * t576 - t574;
t535 = -qJD(1) * t568 - t509 * t628;
t553 = t499 * t514;
t381 = (qJD(1) * t553 + t575) * t513 + t535;
t357 = t380 * t507 - t381 * t506;
t497 = qJDD(3) + qJDD(4);
t508 = sin(qJ(6));
t512 = cos(qJ(6));
t578 = qJD(6) * t512;
t569 = t512 * t357 + t508 * t497 + t499 * t578;
t579 = qJD(6) * t508;
t339 = -t541 * t579 + t569;
t446 = -t509 * t510 + t513 * t514;
t399 = t445 * t506 - t507 * t446;
t633 = t339 * t399;
t383 = -t512 * t499 + t508 * t541;
t562 = t434 * t506 - t507 * t435;
t622 = qJD(6) - t562;
t632 = t383 * t622;
t385 = t499 * t508 + t512 * t541;
t631 = t385 * t622;
t516 = -pkin(1) - pkin(7);
t460 = qJDD(1) * t516 + qJDD(2);
t447 = t514 * t460;
t462 = qJD(1) * t516 + qJD(2);
t583 = qJD(3) * t510;
t389 = qJDD(3) * pkin(3) + pkin(8) * t628 - t462 * t583 + t447;
t582 = qJD(3) * t514;
t566 = t514 * t576;
t629 = t566 + t575;
t395 = -pkin(8) * t629 + t460 * t510 + t462 * t582;
t423 = -pkin(8) * t585 + t462 * t510;
t409 = t513 * t423;
t424 = -pkin(8) * t584 + t514 * t462;
t411 = qJD(3) * pkin(3) + t424;
t542 = -t411 * t509 - t409;
t525 = qJD(4) * t542 + t513 * t389 - t509 * t395;
t335 = pkin(4) * t497 - qJ(5) * t380 + qJD(5) * t434 + t525;
t620 = (qJD(4) * t411 + t395) * t513 + t509 * t389 - t423 * t581;
t337 = -qJ(5) * t381 - qJD(5) * t435 + t620;
t323 = t335 * t507 - t337 * t506;
t321 = -pkin(5) * t497 - t323;
t488 = pkin(10) + t505;
t474 = cos(t488);
t630 = t474 * t496 + t321;
t473 = sin(t488);
t627 = g(3) * t473 + t474 * t495;
t560 = t622 * t512;
t356 = -t380 * t506 - t507 * t381;
t355 = qJDD(6) - t356;
t605 = t355 * t508;
t626 = -t622 * t560 - t605;
t324 = t506 * t335 + t507 * t337;
t427 = t434 * qJ(5);
t408 = t509 * t423;
t564 = t513 * t411 - t408;
t376 = t427 + t564;
t372 = pkin(4) * t499 + t376;
t609 = qJ(5) * t435;
t377 = -t542 - t609;
t603 = t377 * t506;
t345 = t372 * t507 - t603;
t373 = t507 * t377;
t346 = t506 * t372 + t373;
t404 = -t509 * t583 + t513 * t553 - t568;
t369 = t507 * t403 + t404 * t506;
t540 = t507 * t445 + t446 * t506;
t544 = -t403 * t506 + t507 * t404;
t625 = -t323 * t399 + t324 * t540 - t345 * t369 + t346 * t544 - t588;
t616 = pkin(4) * t434;
t360 = pkin(5) * t541 - pkin(9) * t562 - t616;
t475 = pkin(4) * t506 + pkin(9);
t624 = t622 * (qJD(6) * t475 + t360);
t613 = pkin(8) - t516;
t451 = t613 * t510;
t452 = t613 * t514;
t590 = -t513 * t451 - t509 * t452;
t500 = qJDD(1) * qJ(2);
t548 = g(1) * t515 + g(2) * t511;
t501 = qJD(1) * qJD(2);
t570 = 0.2e1 * t501;
t621 = 0.2e1 * t500 + t570 - t548;
t441 = t613 * t583;
t442 = qJD(3) * t452;
t580 = qJD(4) * t513;
t530 = t509 * t441 - t513 * t442 + t451 * t581 - t452 * t580;
t354 = -qJ(5) * t404 - qJD(5) * t445 + t530;
t524 = -qJD(4) * t590 + t513 * t441 + t442 * t509;
t521 = qJ(5) * t403 - qJD(5) * t446 + t524;
t333 = t507 * t354 + t506 * t521;
t343 = -pkin(5) * t499 - t345;
t388 = -qJ(5) * t445 + t590;
t561 = t451 * t509 - t513 * t452;
t533 = -qJ(5) * t446 + t561;
t362 = t507 * t388 + t506 * t533;
t493 = t510 * pkin(3);
t477 = qJ(2) + t493;
t550 = pkin(4) * t445 + t477;
t363 = pkin(5) * t540 + pkin(9) * t399 + t550;
t322 = pkin(9) * t497 + t324;
t454 = pkin(3) * t585 + qJD(1) * qJ(2);
t405 = pkin(4) * t435 + qJD(5) + t454;
t358 = -pkin(5) * t562 - pkin(9) * t541 + t405;
t558 = qJD(6) * t358 + t322;
t618 = -t321 * t399 - t343 * t369 - t362 * t355 - t622 * (qJD(6) * t363 + t333) - t540 * t558;
t614 = g(3) * t474;
t612 = pkin(3) * qJD(4);
t611 = pkin(1) * qJDD(1);
t608 = t339 * t508;
t607 = t343 * t562;
t606 = t343 * t399;
t604 = t363 * t355;
t602 = t383 * t541;
t601 = t385 * t541;
t600 = t506 * t509;
t599 = t507 * t509;
t598 = t508 * t511;
t597 = t508 * t515;
t596 = t511 * t512;
t352 = t512 * t355;
t595 = t512 * t515;
t594 = -t403 * t499 + t446 * t497;
t593 = t513 * t424 - t408;
t378 = t427 + t593;
t563 = -t424 * t509 - t409;
t534 = t563 + t609;
t592 = -t378 * t506 + t507 * t534 + (t506 * t513 + t599) * t612;
t591 = -t507 * t378 - t506 * t534 + (t507 * t513 - t600) * t612;
t481 = pkin(3) * t513 + pkin(4);
t426 = pkin(3) * t599 + t506 * t481;
t589 = t515 * pkin(1) + t511 * qJ(2);
t504 = t514 ^ 2;
t587 = t510 ^ 2 - t504;
t517 = qJD(3) ^ 2;
t586 = -t517 - t518;
t463 = pkin(3) * t582 + qJD(2);
t573 = qJDD(3) * t510;
t412 = pkin(3) * t629 + t500 + t501;
t528 = pkin(4) * t381 + qJDD(5) + t412;
t329 = -pkin(5) * t356 - pkin(9) * t357 + t528;
t344 = pkin(9) * t499 + t346;
t557 = qJD(6) * t344 - t329;
t420 = pkin(9) + t426;
t485 = pkin(3) * t584;
t555 = qJD(6) * t420 + t360 + t485;
t552 = qJDD(2) - t611;
t549 = pkin(4) * t404 + t463;
t546 = -t355 * t420 - t607;
t545 = t345 * t562 + t346 * t541;
t543 = -t404 * t499 - t445 * t497;
t331 = t344 * t512 + t358 * t508;
t538 = t331 * t541 + t343 * t578 + t508 * t630;
t330 = -t344 * t508 + t358 * t512;
t537 = -t330 * t541 + t343 * t579 + t512 * t627;
t536 = t352 + (t508 * t562 - t579) * t622;
t425 = -pkin(3) * t600 + t481 * t507;
t532 = -t369 * t512 + t399 * t579;
t531 = 0.2e1 * qJ(2) * t576 + qJDD(3) * t516;
t348 = t376 * t507 - t603;
t526 = t348 * t622 - t355 * t475 - t607;
t523 = -t516 * t517 + t621;
t470 = t512 * t497;
t340 = qJD(6) * t385 + t357 * t508 - t470;
t522 = -t434 * t435 * MDP(14) - t622 * t541 * MDP(27) + ((t339 - t632) * t512 + (-t340 - t631) * t508) * MDP(24) + (t536 + t602) * MDP(26) + (-t601 - t626) * MDP(25) + (t385 * t560 + t608) * MDP(23) + (t435 * t499 + t380) * MDP(16) + (-t434 * t499 + (-t499 * t584 - t575) * t513 - t535) * MDP(17) + (t434 ^ 2 - t435 ^ 2) * MDP(15) + t497 * MDP(18);
t520 = g(3) * t490 + t454 * t435 + t489 * t588 - t620;
t519 = t454 * t434 + t525 + t634;
t498 = -qJ(5) - pkin(8) - pkin(7);
t492 = t515 * qJ(2);
t487 = qJDD(3) * t514;
t476 = -pkin(4) * t507 - pkin(5);
t449 = pkin(4) * t489 + t493;
t419 = -pkin(5) - t425;
t418 = t473 * t595 - t598;
t417 = t473 * t597 + t596;
t416 = t473 * t596 + t597;
t415 = -t473 * t598 + t595;
t361 = t388 * t506 - t507 * t533;
t347 = t376 * t506 + t373;
t338 = pkin(5) * t544 + pkin(9) * t369 + t549;
t332 = t354 * t506 - t507 * t521;
t328 = t512 * t329;
t1 = [(t332 * t541 + t333 * t562 + t356 * t362 + t357 * t361 - t625) * MDP(21) + t621 * MDP(5) + qJDD(1) * MDP(1) + (-t514 * t517 - t573) * MDP(10) + (qJDD(1) * t504 - 0.2e1 * t510 * t566) * MDP(7) + t588 * MDP(2) + (t324 * t362 + t346 * t333 - t323 * t361 - t345 * t332 + t528 * t550 + t405 * t549 - g(1) * (t449 * t515 + t492 + (-pkin(1) + t498) * t511) - g(2) * (t449 * t511 - t498 * t515 + t589)) * MDP(22) + (-t552 * pkin(1) - g(1) * (-pkin(1) * t511 + t492) - g(2) * t589 + (t570 + t500) * qJ(2)) * MDP(6) + t594 * MDP(16) + (t510 * t523 + t514 * t531) * MDP(12) + (-t510 * t531 + t514 * t523) * MDP(13) + (-t510 * t517 + t487) * MDP(9) + 0.2e1 * (-t510 * t574 + t576 * t587) * MDP(8) + (t385 * t532 - t512 * t633) * MDP(23) + t548 * MDP(3) + (qJDD(2) - t588 - 0.2e1 * t611) * MDP(4) + (-t380 * t445 - t381 * t446 + t403 * t435 + t404 * t434) * MDP(15) + t543 * MDP(17) + (g(1) * t417 - g(2) * t415 - t331 * t544 + t332 * t385 + t361 * t339 + (-(-qJD(6) * t362 + t338) * t622 - t604 + t557 * t540 + qJD(6) * t606) * t508 + t618 * t512) * MDP(29) + (-g(1) * t418 - g(2) * t416 + t328 * t540 + t330 * t544 + t332 * t383 + t361 * t340 + (t338 * t622 + t604 + (-t344 * t540 - t362 * t622 - t606) * qJD(6)) * t512 + t618 * t508) * MDP(28) + (t355 * t540 + t544 * t622) * MDP(27) + (t399 * t605 - t340 * t540 - t544 * t383 + (t369 * t508 + t399 * t578) * t622) * MDP(26) + (t339 * t540 - t352 * t399 + t385 * t544 + t532 * t622) * MDP(25) + (-(-t383 * t512 - t385 * t508) * t369 - (-t608 - t340 * t512 + (t383 * t508 - t385 * t512) * qJD(6)) * t399) * MDP(24) + (t380 * t446 + t403 * t434) * MDP(14) + (t477 * t380 - t454 * t403 + t412 * t446 - t463 * t434 - t490 * t548 - t497 * t590 - t499 * t530) * MDP(20) + (t477 * t381 + t454 * t404 + t412 * t445 + t463 * t435 - t489 * t548 + t497 * t561 + t499 * t524) * MDP(19); qJDD(1) * MDP(4) - t518 * MDP(5) + (t552 + t529) * MDP(6) + (t510 * t586 + t487) * MDP(12) + (t514 * t586 - t573) * MDP(13) + (-qJD(1) * t435 + t594) * MDP(19) + (qJD(1) * t434 + t543) * MDP(20) + (t356 * t540 + t357 * t399 + t369 * t541 + t544 * t562) * MDP(21) + (-qJD(1) * t405 + t625) * MDP(22) + (t340 * t399 + t369 * t383 - t540 * t605) * MDP(28) + (-t352 * t540 + t369 * t385 + t633) * MDP(29) + ((-qJD(1) * t512 - t508 * t544 - t540 * t578) * MDP(28) + (qJD(1) * t508 - t512 * t544 + t540 * t579) * MDP(29)) * t622; (t419 * t339 + t546 * t512 - t627 * t508 + t592 * t385 + (t508 * t555 - t512 * t591) * t622 + t538) * MDP(29) + (t324 * t426 + t323 * t425 - t405 * (t485 - t616) + g(3) * t449 - t588 * (pkin(3) * t514 + pkin(4) * t490) + t591 * t346 - t592 * t345) * MDP(22) + (t356 * t426 - t357 * t425 + t541 * t592 + t562 * t591 + t545) * MDP(21) + qJDD(3) * MDP(11) + (t419 * t340 - t630 * t512 + t546 * t508 + t592 * t383 + (-t508 * t591 - t512 * t555) * t622 + t537) * MDP(28) + (-t563 * t499 + (-t435 * t584 + t497 * t513 - t499 * t581) * pkin(3) + t519) * MDP(19) + (t593 * t499 + (t434 * t584 - t497 * t509 - t499 * t580) * pkin(3) + t520) * MDP(20) + t522 + MDP(9) * t574 - MDP(10) * t575 + (g(3) * t510 + t514 * t529 + t447) * MDP(12) + (g(3) * t514 + (-t460 - t529) * t510) * MDP(13) + (MDP(7) * t510 * t514 - MDP(8) * t587) * t518; (t345 * t347 - t346 * t348 + (t323 * t507 + t324 * t506 + t405 * t434 + t634) * pkin(4)) * MDP(22) + (t476 * t340 - t347 * t383 + t526 * t508 + (-t630 - t624) * t512 + t537) * MDP(28) + (t476 * t339 - t347 * t385 + t526 * t512 + (-t627 + t624) * t508 + t538) * MDP(29) + (-t499 * t542 + t519) * MDP(19) + (t499 * t564 + t520) * MDP(20) + t522 + (-t347 * t541 - t348 * t562 + (t356 * t506 - t357 * t507) * pkin(4) + t545) * MDP(21); (-t541 ^ 2 - t562 ^ 2) * MDP(21) + (t345 * t541 - t346 * t562 + t528 - t548) * MDP(22) + (t536 - t602) * MDP(28) + (-t601 + t626) * MDP(29); t385 * t383 * MDP(23) + (-t383 ^ 2 + t385 ^ 2) * MDP(24) + (t569 + t632) * MDP(25) + (t470 + t631) * MDP(26) + t355 * MDP(27) + (-g(1) * t415 - g(2) * t417 + t331 * t622 - t343 * t385 + t328) * MDP(28) + (g(1) * t416 - g(2) * t418 + t330 * t622 + t343 * t383) * MDP(29) + ((-t322 + t614) * MDP(29) + (-MDP(26) * t541 - MDP(28) * t344 - MDP(29) * t358) * qJD(6)) * t512 + (-qJD(6) * t541 * MDP(25) + (-qJD(6) * t499 - t357) * MDP(26) + (-t558 + t614) * MDP(28) + t557 * MDP(29)) * t508;];
tau  = t1;
