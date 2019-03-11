% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRPRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:54:20
% EndTime: 2019-03-08 21:54:29
% DurationCPUTime: 7.22s
% Computational Cost: add. (5573->480), mult. (13288->642), div. (0->0), fcn. (11023->16), ass. (0->221)
t537 = cos(qJ(6));
t605 = qJD(6) * t537;
t527 = sin(pkin(12));
t530 = cos(pkin(12));
t535 = sin(qJ(3));
t539 = cos(qJ(3));
t492 = -t527 * t535 + t530 * t539;
t485 = t492 * qJD(2);
t538 = cos(qJ(5));
t472 = t538 * t485;
t493 = t527 * t539 + t530 * t535;
t487 = t493 * qJD(2);
t534 = sin(qJ(5));
t428 = -t487 * t534 + t472;
t669 = t428 * t537;
t673 = t605 - t669;
t533 = sin(qJ(6));
t604 = -qJD(6) + t428;
t579 = t604 * t533;
t644 = qJ(4) + pkin(8);
t588 = qJD(3) * t644;
t478 = qJD(4) * t539 - t535 * t588;
t479 = -qJD(4) * t535 - t539 * t588;
t529 = sin(pkin(6));
t540 = cos(qJ(2));
t624 = t529 * t540;
t595 = qJD(1) * t624;
t617 = -t478 * t527 + t530 * t479 + t493 * t595;
t616 = t530 * t478 + t527 * t479 - t492 * t595;
t565 = t485 * t534 + t538 * t487;
t486 = t493 * qJD(3);
t441 = -qJD(2) * t486 + qJDD(2) * t492;
t602 = qJD(2) * qJD(3);
t590 = t539 * t602;
t591 = t535 * t602;
t442 = qJDD(2) * t493 - t527 * t591 + t530 * t590;
t607 = qJD(5) * t534;
t377 = qJD(5) * t472 + t534 * t441 + t538 * t442 - t487 * t607;
t523 = qJDD(3) + qJDD(5);
t524 = qJD(3) + qJD(5);
t597 = t537 * t377 + t533 * t523 + t524 * t605;
t606 = qJD(6) * t533;
t365 = -t565 * t606 + t597;
t362 = t365 * t537;
t414 = t524 * t533 + t537 * t565;
t581 = t377 * t533 - t537 * t523;
t366 = t414 * qJD(6) + t581;
t630 = t565 * t533;
t412 = -t537 * t524 + t630;
t672 = -t533 * t366 - t673 * t412 + t362;
t361 = t365 * t533;
t378 = qJD(5) * t565 - t538 * t441 + t442 * t534;
t375 = qJDD(6) + t378;
t372 = t533 * t375;
t618 = -t604 * t605 + t372;
t632 = t428 * t524;
t634 = t565 * t524;
t636 = t414 * t565;
t671 = t523 * MDP(18) + (-t378 + t634) * MDP(17) - t428 ^ 2 * MDP(15) + (-MDP(14) * t428 + MDP(15) * t565 + MDP(25) * t604) * t565 + (t377 - t632) * MDP(16) + (t673 * t414 + t361) * MDP(21) + (t604 * t669 + t618 - t636) * MDP(23);
t536 = sin(qJ(2));
t608 = qJD(1) * t536;
t596 = t529 * t608;
t578 = t644 * qJD(2) + t596;
t531 = cos(pkin(6));
t609 = qJD(1) * t531;
t453 = t535 * t609 + t539 * t578;
t445 = t527 * t453;
t452 = -t535 * t578 + t539 * t609;
t642 = qJD(3) * pkin(3);
t449 = t452 + t642;
t400 = t530 * t449 - t445;
t649 = pkin(9) * t487;
t384 = qJD(3) * pkin(4) + t400 - t649;
t623 = t530 * t453;
t401 = t527 * t449 + t623;
t650 = pkin(9) * t485;
t387 = t401 + t650;
t363 = t384 * t538 - t387 * t534;
t358 = -pkin(5) * t524 - t363;
t670 = t358 * t428;
t599 = t531 * qJDD(1);
t508 = t539 * t599;
t603 = qJD(1) * qJD(2);
t464 = qJDD(2) * pkin(8) + (qJDD(1) * t536 + t540 * t603) * t529;
t548 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t609 + t464;
t563 = t578 * qJD(3);
t395 = qJDD(3) * pkin(3) - t535 * t548 - t539 * t563 + t508;
t396 = (-t563 + t599) * t535 + t548 * t539;
t370 = t530 * t395 - t396 * t527;
t354 = qJDD(3) * pkin(4) - pkin(9) * t442 + t370;
t371 = t527 * t395 + t530 * t396;
t355 = pkin(9) * t441 + t371;
t364 = t384 * t534 + t387 * t538;
t653 = -qJD(5) * t364 + t538 * t354 - t534 * t355;
t344 = -pkin(5) * t523 - t653;
t641 = cos(pkin(11));
t586 = t641 * t536;
t528 = sin(pkin(11));
t627 = t528 * t540;
t482 = t531 * t586 + t627;
t585 = t641 * t540;
t628 = t528 * t536;
t484 = -t531 * t628 + t585;
t522 = qJ(3) + pkin(12) + qJ(5);
t515 = sin(t522);
t516 = cos(t522);
t587 = t529 * t641;
t626 = t529 * t536;
t629 = t528 * t529;
t556 = -g(3) * (-t515 * t626 + t516 * t531) - g(2) * (-t482 * t515 - t516 * t587) - g(1) * (-t484 * t515 + t516 * t629);
t553 = -t344 + t556;
t489 = t492 * qJD(3);
t668 = pkin(9) * t489 - t617;
t667 = -pkin(9) * t486 + t616;
t444 = t492 * t534 + t493 * t538;
t519 = pkin(3) * t539 + pkin(2);
t466 = -pkin(4) * t492 - t519;
t564 = t538 * t492 - t493 * t534;
t380 = -pkin(5) * t564 - pkin(10) * t444 + t466;
t481 = -t531 * t585 + t628;
t483 = t531 * t627 + t586;
t574 = g(1) * t483 + g(2) * t481;
t555 = g(3) * t624 - t574;
t550 = t555 * t516;
t665 = t380 * t375 - t550;
t393 = pkin(5) * t565 - pkin(10) * t428;
t434 = t482 * t516 - t515 * t587;
t436 = t484 * t516 + t515 * t629;
t475 = -qJD(2) * t519 + qJD(4) - t595;
t437 = -pkin(4) * t485 + t475;
t462 = t515 * t531 + t516 * t626;
t654 = (qJD(5) * t384 + t355) * t538 + t534 * t354 - t387 * t607;
t664 = g(1) * t436 + g(2) * t434 + g(3) * t462 - t437 * t428 - t654;
t637 = t412 * t565;
t521 = t535 * t642;
t571 = pkin(4) * t486 + t521 - t596;
t651 = pkin(3) * t530;
t517 = pkin(4) + t651;
t652 = pkin(3) * t527;
t613 = t534 * t517 + t538 * t652;
t373 = t537 * t375;
t558 = -t604 * t606 - t373;
t359 = pkin(10) * t524 + t364;
t376 = -pkin(5) * t428 - pkin(10) * t565 + t437;
t568 = t359 * t533 - t376 * t537;
t660 = t358 * t606 + t565 * t568;
t346 = t359 * t537 + t376 * t533;
t659 = t346 * t565 + t358 * t605 - t553 * t533;
t658 = -t437 * t565 + t556 + t653;
t343 = pkin(10) * t523 + t654;
t500 = t644 * t535;
t501 = t644 * t539;
t457 = -t530 * t500 - t501 * t527;
t418 = -pkin(9) * t493 + t457;
t458 = -t527 * t500 + t530 * t501;
t419 = pkin(9) * t492 + t458;
t383 = t418 * t534 + t419 * t538;
t397 = qJD(5) * t564 - t486 * t534 + t489 * t538;
t573 = g(1) * t484 + g(2) * t482;
t554 = -g(3) * t626 - t573;
t567 = t418 * t538 - t419 * t534;
t621 = -qJD(5) * t567 + t668 * t534 - t667 * t538;
t657 = (qJD(6) * t376 + t343) * t564 + t344 * t444 + t358 * t397 - (-qJD(6) * t380 + t621) * t604 - t383 * t375 + t554;
t592 = t536 * t603;
t506 = t529 * t592;
t541 = qJD(3) ^ 2;
t589 = qJDD(1) * t624;
t655 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t541 + t529 * (-g(3) * t540 + t592) - t506 + t574 + t589;
t490 = t531 * t539 - t535 * t626;
t645 = g(3) * t490;
t643 = qJD(2) * pkin(2);
t639 = t358 * t444;
t635 = t414 * t533;
t625 = t529 * t539;
t622 = qJDD(1) - g(3);
t620 = qJD(5) * t383 + t667 * t534 + t668 * t538;
t405 = t530 * t452 - t445;
t403 = -t452 * t527 - t623;
t390 = t403 - t650;
t391 = t405 - t649;
t561 = t517 * t538 - t534 * t652;
t615 = -t561 * qJD(5) + t390 * t534 + t391 * t538;
t614 = qJD(5) * t613 + t390 * t538 - t391 * t534;
t525 = t535 ^ 2;
t612 = -t539 ^ 2 + t525;
t601 = qJDD(2) * t535;
t600 = qJDD(2) * t539;
t520 = t535 * qJD(2) * pkin(3);
t594 = qJD(2) * t624;
t459 = pkin(4) * t487 + t520;
t584 = t529 * t622;
t477 = pkin(10) + t613;
t575 = qJD(6) * t477 + t393 + t459;
t398 = qJD(5) * t444 + t538 * t486 + t489 * t534;
t572 = pkin(5) * t398 - pkin(10) * t397 + t571;
t569 = -t477 * t375 - t670;
t491 = t531 * t535 + t536 * t625;
t423 = t490 * t530 - t491 * t527;
t424 = t490 * t527 + t491 * t530;
t566 = t423 * t538 - t424 * t534;
t386 = t423 * t534 + t424 * t538;
t562 = -t428 * t579 - t558;
t560 = -g(1) * t528 + g(2) * t641;
t557 = t397 * t537 - t444 * t606;
t499 = -t595 - t643;
t552 = -qJD(2) * t499 - t464 + t573;
t551 = pkin(3) * t591 - qJDD(2) * t519 + qJDD(4) + t506;
t546 = -pkin(8) * qJDD(3) + (t499 + t595 - t643) * qJD(3);
t440 = t551 - t589;
t399 = -pkin(4) * t441 + t440;
t542 = qJD(2) ^ 2;
t476 = -pkin(5) - t561;
t451 = -qJD(3) * t491 - t535 * t594;
t450 = qJD(3) * t490 + t539 * t594;
t404 = t450 * t530 + t451 * t527;
t402 = -t450 * t527 + t451 * t530;
t350 = qJD(5) * t386 - t402 * t538 + t404 * t534;
t349 = qJD(5) * t566 + t402 * t534 + t404 * t538;
t348 = pkin(5) * t378 - pkin(10) * t377 + t399;
t347 = t537 * t348;
t1 = [t622 * MDP(1) + (qJD(3) * t451 + qJDD(3) * t490) * MDP(10) + (-qJD(3) * t450 - qJDD(3) * t491) * MDP(11) + (-t402 * t487 + t404 * t485 - t423 * t442 + t424 * t441) * MDP(12) + (t370 * t423 + t371 * t424 + t400 * t402 + t401 * t404 - g(3)) * MDP(13) + (-t350 * t524 + t523 * t566) * MDP(19) + (-t349 * t524 - t386 * t523) * MDP(20) + (-(-t349 * t533 - t386 * t605) * t604 - t386 * t372 + t350 * t412 - t566 * t366) * MDP(26) + ((t349 * t537 - t386 * t606) * t604 - t386 * t373 + t350 * t414 - t566 * t365) * MDP(27) + ((-qJDD(2) * MDP(4) + (-t539 * MDP(10) + t535 * MDP(11) - MDP(3)) * t542 + (MDP(13) * t475 - MDP(19) * t428 + MDP(20) * t565 - (MDP(26) * t537 - MDP(27) * t533) * t604) * qJD(2)) * t536 + (qJDD(2) * MDP(3) - t542 * MDP(4) + (-t591 + t600) * MDP(10) + (-t590 - t601) * MDP(11) - t440 * MDP(13) - t378 * MDP(19) - t377 * MDP(20) + t558 * MDP(26) + t618 * MDP(27)) * t540) * t529; qJDD(2) * MDP(2) + (t622 * t624 + t574) * MDP(3) + (-t536 * t584 + t573) * MDP(4) + (qJDD(2) * t525 + 0.2e1 * t535 * t590) * MDP(5) + 0.2e1 * (t535 * t600 - t602 * t612) * MDP(6) + (qJDD(3) * t535 + t539 * t541) * MDP(7) + (qJDD(3) * t539 - t535 * t541) * MDP(8) + (t546 * t535 + t539 * t655) * MDP(10) + (-t535 * t655 + t546 * t539) * MDP(11) + (-t370 * t493 + t371 * t492 - t400 * t489 - t401 * t486 + t441 * t458 - t442 * t457 + t485 * t616 - t487 * t617 + t554) * MDP(12) + (t371 * t458 + t370 * t457 - t440 * t519 + t475 * t521 - g(1) * (-t483 * t519 + t484 * t644) - g(2) * (-t481 * t519 + t482 * t644) + t616 * t401 + t617 * t400 + (-t475 * t608 - g(3) * (t519 * t540 + t536 * t644)) * t529) * MDP(13) + (t377 * t444 + t397 * t565) * MDP(14) + (t377 * t564 - t378 * t444 + t397 * t428 - t398 * t565) * MDP(15) + (t397 * t524 + t444 * t523) * MDP(16) + (-t398 * t524 + t523 * t564) * MDP(17) + (t378 * t466 + t398 * t437 - t399 * t564 - t428 * t571 + t523 * t567 - t524 * t620 - t550) * MDP(19) + (t377 * t466 - t383 * t523 + t397 * t437 + t399 * t444 + t515 * t555 + t524 * t621 + t565 * t571) * MDP(20) + (t362 * t444 + t557 * t414) * MDP(21) + ((-t412 * t537 - t635) * t397 + (-t361 - t366 * t537 + (t412 * t533 - t414 * t537) * qJD(6)) * t444) * MDP(22) + (-t365 * t564 + t373 * t444 + t398 * t414 - t557 * t604) * MDP(23) + (-t444 * t372 + t366 * t564 - t398 * t412 - (-t397 * t533 - t444 * t605) * t604) * MDP(24) + (-t375 * t564 - t398 * t604) * MDP(25) + (-t568 * t398 - t347 * t564 - t567 * t366 + t620 * t412 + (-t572 * t604 + (t359 * t564 + t383 * t604 + t639) * qJD(6) + t665) * t537 + t657 * t533) * MDP(26) + (-t346 * t398 - t567 * t365 + t620 * t414 + ((-qJD(6) * t359 + t348) * t564 - qJD(6) * t639 - (qJD(6) * t383 - t572) * t604 - t665) * t533 + t657 * t537) * MDP(27); (-MDP(5) * t535 * t539 + MDP(6) * t612) * t542 + (t604 * t635 + t672) * MDP(22) + (t562 + t637) * MDP(24) + (-t459 * t565 - t523 * t613 + t524 * t615 + t664) * MDP(20) + MDP(8) * t600 + MDP(7) * t601 + (t535 * t552 + t560 * t625 + t508 - t645) * MDP(10) + (t370 * t651 + t371 * t652 - t400 * t403 - t401 * t405 - t475 * t520 + (-g(1) * (-t484 * t535 + t528 * t625) - g(2) * (-t482 * t535 - t539 * t587) - t645) * pkin(3)) * MDP(13) + ((t401 + t403) * t487 + (t400 - t405) * t485 + (t441 * t527 - t442 * t530) * pkin(3)) * MDP(12) + (g(3) * t491 + (-t529 * t560 - t599) * t535 + t552 * t539) * MDP(11) + (t476 * t366 + t614 * t412 + (-t604 * t615 + t569) * t533 + (t575 * t604 + t553) * t537 + t660) * MDP(26) + (t476 * t365 + t569 * t537 + t614 * t414 - (t533 * t575 + t537 * t615) * t604 + t659) * MDP(27) + (t428 * t459 + t523 * t561 - t524 * t614 + t658) * MDP(19) + qJDD(3) * MDP(9) + t671; (-t485 ^ 2 - t487 ^ 2) * MDP(12) + (t400 * t487 - t401 * t485 - t540 * t584 + t551 - t574) * MDP(13) + (t378 + t634) * MDP(19) + (t377 + t632) * MDP(20) + (t562 - t637) * MDP(26) + (-t537 * t604 ^ 2 - t372 - t636) * MDP(27); (t364 * t524 + t658) * MDP(19) + (t363 * t524 + t664) * MDP(20) + (t414 * t579 + t672) * MDP(22) + (-t579 * t604 + t373 + t637) * MDP(24) + (-pkin(5) * t366 - t364 * t412 + (-pkin(10) * t375 - t363 * t604 - t670) * t533 + (-(-pkin(10) * qJD(6) - t393) * t604 + t553) * t537 + t660) * MDP(26) + (-pkin(5) * t365 - (t363 * t537 + t393 * t533) * t604 - t364 * t414 - t358 * t669 + t558 * pkin(10) + t659) * MDP(27) + t671; t414 * t412 * MDP(21) + (-t412 ^ 2 + t414 ^ 2) * MDP(22) + (-t412 * t604 + t597) * MDP(23) + (-t414 * t604 - t581) * MDP(24) + t375 * MDP(25) + (-t533 * t343 + t347 - t346 * t604 - t358 * t414 - g(1) * (-t436 * t533 + t483 * t537) - g(2) * (-t434 * t533 + t481 * t537) - g(3) * (-t462 * t533 - t537 * t624)) * MDP(26) + (-t537 * t343 - t533 * t348 + t568 * t604 + t358 * t412 - g(1) * (-t436 * t537 - t483 * t533) - g(2) * (-t434 * t537 - t481 * t533) - g(3) * (-t462 * t537 + t533 * t624)) * MDP(27) + (-MDP(23) * t630 - MDP(24) * t414 - MDP(26) * t346 + MDP(27) * t568) * qJD(6);];
tau  = t1;
