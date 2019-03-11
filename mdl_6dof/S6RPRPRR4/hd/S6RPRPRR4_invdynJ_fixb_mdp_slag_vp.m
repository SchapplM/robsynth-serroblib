% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:39
% EndTime: 2019-03-09 03:46:46
% DurationCPUTime: 5.97s
% Computational Cost: add. (3064->475), mult. (6010->615), div. (0->0), fcn. (3878->14), ass. (0->237)
t489 = sin(pkin(10));
t463 = pkin(1) * t489 + pkin(7);
t446 = t463 * qJDD(1);
t643 = qJD(2) * qJD(3) + t446;
t493 = sin(qJ(3));
t569 = qJD(1) * qJD(3);
t554 = t493 * t569;
t497 = cos(qJ(3));
t565 = qJDD(1) * t497;
t639 = t554 - t565;
t449 = t463 * qJD(1);
t589 = t497 * qJD(2) - t493 * t449;
t636 = qJD(4) - t589;
t642 = qJDD(1) * MDP(12);
t622 = pkin(4) + t463;
t492 = sin(qJ(5));
t496 = cos(qJ(5));
t584 = qJD(1) * t497;
t436 = qJD(3) * t492 + t496 * t584;
t495 = cos(qJ(6));
t555 = t492 * t584;
t580 = qJD(3) * t496;
t438 = -t555 + t580;
t491 = sin(qJ(6));
t610 = t438 * t491;
t367 = t495 * t436 + t610;
t585 = qJD(1) * t493;
t462 = qJD(5) + t585;
t456 = qJD(6) + t462;
t641 = t367 * t456;
t533 = t436 * t491 - t495 * t438;
t640 = t456 * t533;
t439 = t491 * t496 + t492 * t495;
t523 = t439 * t493;
t630 = qJD(5) + qJD(6);
t595 = -qJD(1) * t523 - t630 * t439;
t388 = -qJD(3) * pkin(3) + t636;
t574 = qJD(5) * t497;
t556 = t492 * t574;
t638 = t493 * t580 + t556;
t553 = t497 * t569;
t566 = qJDD(1) * t493;
t521 = t553 + t566;
t431 = qJDD(5) + t521;
t418 = t496 * t431;
t582 = qJD(3) * t436;
t637 = (t582 - t418) * MDP(21);
t499 = -pkin(3) - pkin(8);
t570 = pkin(4) * t585 + t636;
t372 = qJD(3) * t499 + t570;
t490 = cos(pkin(10));
t626 = pkin(1) * t490;
t464 = -pkin(2) - t626;
t620 = qJ(4) * t493;
t403 = t497 * t499 + t464 - t620;
t379 = t403 * qJD(1);
t337 = t372 * t492 + t379 * t496;
t334 = -pkin(9) * t436 + t337;
t572 = qJD(6) * t491;
t332 = t334 * t572;
t406 = t493 * qJD(2) + t497 * t449;
t383 = pkin(4) * t584 + t406;
t485 = qJD(3) * qJ(4);
t376 = t485 + t383;
t353 = pkin(5) * t436 + t376;
t482 = qJ(1) + pkin(10);
t472 = sin(t482);
t473 = cos(t482);
t488 = qJ(5) + qJ(6);
t479 = cos(t488);
t478 = sin(t488);
t609 = t478 * t493;
t385 = t472 * t479 + t473 * t609;
t387 = -t472 * t609 + t473 * t479;
t481 = g(3) * t497;
t635 = g(1) * t385 - g(2) * t387 + t353 * t367 - t478 * t481 + t332;
t608 = t479 * t493;
t384 = -t472 * t478 + t473 * t608;
t386 = t472 * t608 + t473 * t478;
t360 = -t436 * qJD(5) + t496 * qJDD(3) + t492 * t639;
t564 = qJDD(2) * t497;
t579 = qJD(3) * t497;
t629 = t449 * t579 + t493 * t643;
t516 = qJDD(4) - t564 + t629;
t343 = t521 * pkin(4) + qJDD(3) * t499 + t516;
t461 = pkin(3) * t554;
t619 = qJ(4) * t497;
t536 = pkin(8) * t493 - t619;
t578 = qJD(4) * t493;
t508 = qJD(3) * t536 - t578;
t347 = qJD(1) * t508 + qJDD(1) * t403 + t461;
t549 = t496 * t343 - t347 * t492;
t506 = -t337 * qJD(5) + t549;
t321 = pkin(5) * t431 - pkin(9) * t360 + t506;
t361 = -qJD(5) * t555 + qJDD(3) * t492 + (qJD(3) * qJD(5) - t639) * t496;
t575 = qJD(5) * t496;
t562 = -t492 * t343 - t496 * t347 - t372 * t575;
t576 = qJD(5) * t492;
t322 = -pkin(9) * t361 - t379 * t576 - t562;
t550 = t495 * t321 - t491 * t322;
t634 = -g(1) * t384 - g(2) * t386 + t353 * t533 + t479 * t481 + t550;
t421 = qJDD(6) + t431;
t633 = t421 * MDP(27) + (-t367 ^ 2 + t533 ^ 2) * MDP(24) - t367 * MDP(23) * t533;
t423 = t622 * t493;
t410 = t492 * t423;
t593 = t496 * t403 + t410;
t391 = -t485 - t406;
t631 = g(1) * t473 + g(2) * t472;
t628 = t376 * t462 + t499 * t431;
t548 = t360 * t491 + t495 * t361;
t327 = -qJD(6) * t533 + t548;
t627 = t497 * t630;
t623 = g(3) * t493;
t621 = pkin(9) - t499;
t618 = qJDD(3) * pkin(3);
t617 = t327 * t493;
t336 = t496 * t372 - t379 * t492;
t333 = -pkin(9) * t438 + t336;
t331 = pkin(5) * t462 + t333;
t616 = t331 * t495;
t615 = t334 * t495;
t614 = t360 * t496;
t613 = t361 * t493;
t612 = t436 * t462;
t611 = t438 * t462;
t607 = t491 * t492;
t606 = t492 * t431;
t605 = t492 * t493;
t604 = t492 * t497;
t603 = t493 * t496;
t602 = t495 * t496;
t601 = t496 * t497;
t501 = qJD(1) ^ 2;
t600 = t497 * t501;
t571 = qJD(6) * t495;
t561 = t495 * t360 - t491 * t361 - t436 * t571;
t326 = -t438 * t572 + t561;
t598 = t326 * t493 - t533 * t579;
t532 = -t602 + t607;
t581 = qJD(3) * t493;
t342 = t439 * t627 - t532 * t581;
t408 = t532 * t497;
t597 = t342 * t456 + t408 * t421;
t596 = t360 * t493 + t438 * t579;
t559 = t496 * t585;
t594 = -t491 * t576 - t492 * t572 + t495 * t559 - t585 * t607 + t630 * t602;
t471 = pkin(3) * t585;
t407 = qJD(1) * t536 + t471;
t592 = t492 * t383 + t496 * t407;
t591 = t638 * t462;
t424 = t622 * t497;
t560 = -pkin(5) * t496 - pkin(4);
t590 = pkin(5) * t575 - t560 * t585 + t636;
t486 = t493 ^ 2;
t487 = t497 ^ 2;
t588 = t486 - t487;
t538 = pkin(3) * t497 + t620;
t528 = pkin(2) + t538;
t422 = -t528 - t626;
t392 = qJD(1) * t422;
t441 = -qJ(4) * t584 + t471;
t586 = qJD(1) * t441;
t450 = qJD(1) * t464;
t583 = qJD(3) * t367;
t577 = qJD(5) * t379;
t573 = qJD(5) * t499;
t567 = qJDD(1) * t422;
t563 = qJDD(3) * t463;
t558 = t492 * t581;
t443 = t621 * t496;
t552 = pkin(9) * t497 - t403;
t470 = pkin(3) * t581;
t389 = t470 + t508;
t415 = t622 * t579;
t547 = -t389 * t492 + t496 * t415;
t351 = t516 - t618;
t546 = -qJD(3) * t391 - t351;
t545 = qJD(6) * t331 + t322;
t544 = t493 * qJDD(2) - t449 * t581 + t497 * t643;
t494 = sin(qJ(1));
t498 = cos(qJ(1));
t542 = g(1) * t494 - g(2) * t498;
t541 = -t532 * t421 + t456 * t595;
t378 = t496 * t383;
t442 = t621 * t492;
t529 = pkin(5) * t497 - pkin(9) * t605;
t540 = qJD(1) * t529 - qJD(6) * t442 - t407 * t492 - t621 * t576 + t378;
t539 = pkin(9) * t559 + t443 * t630 + t592;
t537 = pkin(3) * t493 - t619;
t324 = t331 * t491 + t615;
t341 = qJD(3) * t523 + t532 * t627;
t409 = t439 * t497;
t534 = -t341 * t456 + t409 * t421;
t531 = t462 * t492;
t526 = t542 * pkin(1);
t483 = qJDD(3) * qJ(4);
t484 = qJD(3) * qJD(4);
t350 = -t483 - t484 - t544;
t525 = t462 * t575 + t606;
t524 = -qJ(4) * t579 - t578;
t522 = -t421 * t439 - t456 * t594;
t520 = t496 * t389 - t403 * t576 + t492 * t415 + t423 * t575;
t519 = qJD(3) * t589 - t544;
t518 = qJD(3) * t406 - t481 - t629;
t500 = qJD(3) ^ 2;
t517 = g(1) * t472 - g(2) * t473 - t463 * t500;
t515 = -qJD(1) * t450 + t631;
t514 = -t496 * t574 + t558;
t513 = 0.2e1 * qJD(3) * t450 - t563;
t512 = -0.2e1 * qJD(3) * t392 + t563;
t510 = -t497 * t631 - t623;
t344 = -pkin(4) * t639 - t350;
t507 = t344 + t510;
t505 = -0.2e1 * qJDD(1) * t464 + t517;
t504 = -t350 * t497 + t351 * t493 + (t388 * t497 + t391 * t493) * qJD(3);
t352 = qJD(1) * t524 + t461 + t567;
t416 = t470 + t524;
t503 = -qJD(1) * t416 - t352 + t517 - t567;
t465 = pkin(5) * t492 + qJ(4);
t445 = qJDD(3) * t497 - t493 * t500;
t444 = qJDD(3) * t493 + t497 * t500;
t414 = t622 * t581;
t411 = t496 * t423;
t402 = -t472 * t605 + t473 * t496;
t401 = t472 * t603 + t473 * t492;
t400 = t472 * t496 + t473 * t605;
t399 = -t472 * t492 + t473 * t603;
t393 = pkin(5) * t601 + t424;
t380 = t392 * t585;
t364 = -pkin(5) * t556 + (-t463 + t560) * t581;
t349 = -pkin(9) * t601 + t593;
t348 = pkin(5) * t493 + t492 * t552 + t411;
t330 = pkin(5) * t361 + t344;
t329 = pkin(9) * t638 + t520;
t328 = t529 * qJD(3) + (t496 * t552 - t410) * qJD(5) + t547;
t323 = -t334 * t491 + t616;
t1 = [(g(1) * t498 + g(2) * t494) * MDP(3) + (-t324 * t579 + g(1) * t386 - g(2) * t384 + t393 * t326 - t330 * t409 + t332 * t493 + t353 * t341 - t364 * t533 + (-(-qJD(6) * t349 + t328) * t456 - t348 * t421 - t321 * t493) * t491 + (-(qJD(6) * t348 + t329) * t456 - t349 * t421 - t545 * t493) * t495) * MDP(29) + (t326 * t408 + t327 * t409 - t341 * t367 - t342 * t533) * MDP(24) + (-t326 * t409 - t341 * t533) * MDP(23) + t444 * MDP(7) + t445 * MDP(8) + ((t486 + t487) * t446 + t504 - t631) * MDP(12) + (t547 * t462 + (-t403 * t492 + t411) * t431 + t549 * t493 - t414 * t436 + t424 * t361 + t344 * t601 - g(1) * t402 - g(2) * t400 + (t336 * t497 - t376 * t603) * qJD(3) + (-t337 * t493 - t376 * t604 - t593 * t462) * qJD(5)) * MDP(21) + (-t367 * t579 + t597 - t617) * MDP(26) + (-t613 + (-t582 - t418) * t497 + t591) * MDP(19) + ((-t436 * t492 + t438 * t496) * t581 + (-t614 + t361 * t492 + (t436 * t496 + t438 * t492) * qJD(5)) * t497) * MDP(17) + (-t431 * t604 + t462 * t514 + t596) * MDP(18) + (-t360 * t604 + t438 * t514) * MDP(16) + (-t534 + t598) * MDP(25) + 0.2e1 * (t493 * t565 - t569 * t588) * MDP(6) + (t421 * t493 + t456 * t579) * MDP(27) + ((t328 * t495 - t329 * t491) * t456 + (t348 * t495 - t349 * t491) * t421 + t550 * t493 + t323 * t579 + t364 * t367 + t393 * t327 - t330 * t408 - t353 * t342 - g(1) * t387 - g(2) * t385 + ((-t348 * t491 - t349 * t495) * t456 - t324 * t493) * qJD(6)) * MDP(28) + (t431 * t493 + t462 * t579) * MDP(20) + qJDD(1) * MDP(1) + (-t520 * t462 - t593 * t431 - t414 * t438 + t424 * t360 + g(1) * t401 - g(2) * t399 + ((qJD(3) * t376 + t577) * t492 + t562) * t493 + (-qJD(3) * t337 - t344 * t492 - t376 * t575) * t497) * MDP(22) + (qJDD(1) * t486 + 0.2e1 * t493 * t553) * MDP(5) + t542 * MDP(2) + (t352 * t422 + t392 * t416 + t526 + (-g(1) * pkin(7) - g(2) * t528) * t473 + (-g(2) * pkin(7) + g(1) * t528) * t472 + t504 * t463) * MDP(15) + ((t489 ^ 2 + t490 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t526) * MDP(4) + (t493 * t503 + t497 * t512) * MDP(14) + (t493 * t512 - t497 * t503) * MDP(13) + (t493 * t513 + t497 * t505) * MDP(10) + (-t493 * t505 + t497 * t513) * MDP(11); (qJDD(2) - g(3)) * MDP(4) + (-t350 * t493 + t388 * t581 - g(3)) * MDP(15) + (t591 + t613) * MDP(21) + (-t462 * t558 + t596) * MDP(22) + (t597 + t617) * MDP(28) + (t534 + t598) * MDP(29) + (MDP(10) - MDP(13)) * t445 + (-MDP(11) + MDP(14)) * t444 + (t546 * MDP(15) + t525 * MDP(22) + MDP(28) * t583 + t637) * t497; (-t351 * pkin(3) - g(3) * t538 - t350 * qJ(4) - t388 * t406 - t391 * t636 - t392 * t441 + t631 * t537) * MDP(15) + (qJ(4) * t361 - t378 * t462 + t570 * t436 + t628 * t496 + ((t407 - t573) * t462 + t507) * t492) * MDP(21) + (qJ(4) * t360 + t592 * t462 + t570 * t438 - t628 * t492 + (-t462 * t573 + t507) * t496) * MDP(22) - t537 * t642 + qJDD(3) * MDP(9) + (t497 * t515 + t519 + t623) * MDP(11) + (-0.2e1 * t618 + qJDD(4) + t380 + (-qJDD(2) - t586) * t497 - t631 * t493 - t518) * MDP(13) + ((-t361 - t611) * t496 + (-t360 + t612) * t492) * MDP(17) + (-t462 * t576 + t418 + (-t438 * t497 - t462 * t605) * qJD(1)) * MDP(18) + ((t436 * t497 - t462 * t603) * qJD(1) - t525) * MDP(19) + t522 * MDP(26) + (0.2e1 * t483 + 0.2e1 * t484 + (-g(3) + t586) * t493 + (qJD(1) * t392 - t631) * t497 - t519) * MDP(14) + (-t326 * t532 - t533 * t595) * MDP(23) + (-t326 * t439 + t327 * t532 - t367 * t595 + t533 * t594) * MDP(24) + ((t442 * t491 - t443 * t495) * t421 + t465 * t327 + t330 * t439 + (t491 * t539 - t495 * t540) * t456 + t590 * t367 + t594 * t353 + t510 * t478) * MDP(28) + (-(-t442 * t495 - t443 * t491) * t421 + t465 * t326 - t330 * t532 + (t491 * t540 + t495 * t539) * t456 - t590 * t533 + t595 * t353 + t510 * t479) * MDP(29) + (-t438 * t531 + t614) * MDP(16) + t541 * MDP(25) + (t493 * t515 + t518 + t564) * MDP(10) + t588 * MDP(6) * t501 - t493 * MDP(5) * t600 + MDP(8) * t565 + MDP(7) * t566 + (-t462 * MDP(20) - t336 * MDP(21) + t337 * MDP(22) + MDP(25) * t533 + t367 * MDP(26) - t456 * MDP(27) - t323 * MDP(28) + t324 * MDP(29)) * t584; qJDD(3) * MDP(13) + (-t486 * t501 - t500) * MDP(14) + (t380 + t481 - t546) * MDP(15) - t637 + (-qJD(3) * t438 - t606) * MDP(22) + (t541 - t583) * MDP(28) + (qJD(3) * t533 + t522) * MDP(29) + (MDP(13) * t600 - MDP(15) * t631 + t642) * t493 + (-MDP(22) * t462 * t496 - MDP(21) * t531) * t462; t438 * t436 * MDP(16) + (-t436 ^ 2 + t438 ^ 2) * MDP(17) + (t360 + t612) * MDP(18) + (-t361 + t611) * MDP(19) + t431 * MDP(20) + (-g(1) * t399 - g(2) * t401 + g(3) * t601 + t337 * t462 - t376 * t438 + t506) * MDP(21) + (g(1) * t400 - g(2) * t402 + t336 * t462 + t376 * t436 + (t577 - t481) * t492 + t562) * MDP(22) + (t326 + t641) * MDP(25) + (-t327 - t640) * MDP(26) + (-(-t333 * t491 - t615) * t456 - t324 * qJD(6) + (-t367 * t438 + t421 * t495 - t456 * t572) * pkin(5) + t634) * MDP(28) + ((-t334 * t456 - t321) * t491 + (t333 * t456 - t545) * t495 + (-t421 * t491 + t438 * t533 - t456 * t571) * pkin(5) + t635) * MDP(29) + t633; (t561 + t641) * MDP(25) + (-t548 - t640) * MDP(26) + (t324 * t456 + t634) * MDP(28) + (-t491 * t321 - t495 * t322 + t323 * t456 + t635) * MDP(29) + (-MDP(25) * t610 + MDP(26) * t533 - MDP(28) * t324 - MDP(29) * t616) * qJD(6) + t633;];
tau  = t1;
