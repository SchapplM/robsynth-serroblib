% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPPRR8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:22
% EndTime: 2019-03-09 09:26:35
% DurationCPUTime: 8.38s
% Computational Cost: add. (4224->539), mult. (10469->716), div. (0->0), fcn. (7535->8), ass. (0->225)
t538 = cos(qJ(2));
t608 = qJD(1) * t538;
t517 = qJD(5) + t608;
t506 = qJD(6) + t517;
t536 = cos(qJ(6));
t530 = sin(pkin(10));
t535 = sin(qJ(2));
t609 = qJD(1) * t535;
t586 = t530 * t609;
t531 = cos(pkin(10));
t594 = t531 * qJD(2);
t475 = t586 - t594;
t584 = t531 * t609;
t607 = qJD(2) * t530;
t477 = t584 + t607;
t534 = sin(qJ(5));
t537 = cos(qJ(5));
t555 = -t537 * t475 + t477 * t534;
t533 = sin(qJ(6));
t554 = t475 * t534 + t477 * t537;
t634 = t554 * t533;
t352 = t536 * t555 + t634;
t658 = t352 * t506;
t556 = t533 * t555 - t536 * t554;
t657 = t506 * t556;
t482 = t530 * t534 + t531 * t537;
t550 = t482 * t538;
t616 = -qJD(1) * t550 - t482 * qJD(5);
t583 = t531 * t608;
t585 = t530 * t608;
t597 = qJD(5) * t537;
t598 = qJD(5) * t534;
t615 = t530 * t597 - t531 * t598 - t534 * t583 + t537 * t585;
t639 = qJ(3) * t535;
t493 = -pkin(2) * t538 - pkin(1) - t639;
t467 = t493 * qJD(1);
t521 = pkin(7) * t608;
t500 = qJD(2) * qJ(3) + t521;
t424 = t467 * t531 - t530 * t500;
t407 = pkin(3) * t608 + qJD(4) - t424;
t370 = pkin(4) * t608 - pkin(8) * t477 + t407;
t425 = t530 * t467 + t531 * t500;
t410 = -qJ(4) * t608 + t425;
t376 = pkin(8) * t475 + t410;
t343 = t370 * t534 + t376 * t537;
t339 = -pkin(9) * t555 + t343;
t596 = qJD(6) * t533;
t337 = t339 * t596;
t520 = pkin(7) * t609;
t532 = qJD(2) * pkin(2);
t492 = qJD(3) + t520 - t532;
t400 = t475 * pkin(3) - t477 * qJ(4) + t492;
t373 = -pkin(4) * t475 - t400;
t349 = pkin(5) * t555 + t373;
t656 = t349 * t352 + t337;
t593 = qJD(1) * qJD(2);
t580 = t538 * t593;
t502 = t530 * t580;
t566 = t531 * t580;
t371 = t475 * t597 - t477 * t598 + t534 * t502 + t537 * t566;
t559 = pkin(2) * t535 - qJ(3) * t538;
t458 = qJD(2) * t559 - qJD(3) * t535;
t448 = t458 * qJD(1);
t490 = (qJD(3) - t520) * qJD(2);
t404 = t448 * t531 - t530 * t490;
t625 = t531 * t538;
t592 = pkin(8) * t625;
t641 = -pkin(3) - pkin(4);
t362 = (t535 * t641 - t592) * t593 - t404;
t405 = t530 * t448 + t531 * t490;
t581 = t535 * t593;
t589 = qJ(4) * t581 + t405;
t363 = (pkin(8) * t607 - qJD(4)) * t608 + t589;
t574 = t537 * t362 - t363 * t534;
t543 = -t343 * qJD(5) + t574;
t330 = -pkin(5) * t581 - pkin(9) * t371 + t543;
t372 = qJD(5) * t554 - t537 * t502 + t534 * t566;
t549 = t534 * t362 + t537 * t363 + t370 * t597 - t376 * t598;
t331 = -pkin(9) * t372 + t549;
t575 = t536 * t330 - t533 * t331;
t655 = t349 * t556 + t575;
t578 = MDP(30) * t609;
t654 = (-t352 ^ 2 + t556 ^ 2) * MDP(27) - qJD(2) * t578 - t352 * MDP(26) * t556;
t653 = -0.2e1 * t593;
t652 = MDP(4) * t535;
t529 = t538 ^ 2;
t651 = MDP(5) * (t535 ^ 2 - t529);
t650 = t517 * t554;
t649 = t517 * t555;
t627 = t530 * t538;
t514 = pkin(7) * t627;
t527 = t538 * pkin(3);
t413 = pkin(4) * t538 + t514 + t527 + (-pkin(8) * t535 - t493) * t531;
t515 = pkin(7) * t625;
t444 = t530 * t493 + t515;
t435 = -qJ(4) * t538 + t444;
t629 = t530 * t535;
t423 = pkin(8) * t629 + t435;
t618 = t534 * t413 + t537 * t423;
t600 = qJD(4) * t530;
t613 = qJ(4) * t583 - t521;
t569 = -t585 * t641 + t600 - t613;
t640 = -pkin(8) + qJ(3);
t497 = t640 * t530;
t498 = t640 * t531;
t614 = t534 * t497 + t537 * t498;
t588 = -pkin(7) * t530 - pkin(3);
t544 = -t592 + (-pkin(4) + t588) * t535;
t486 = t559 * qJD(1);
t632 = t486 * t531;
t390 = qJD(1) * t544 - t632;
t601 = qJD(3) * t537;
t646 = -t537 * t390 + t530 * t601;
t599 = qJD(4) * t538;
t606 = qJD(2) * t535;
t645 = qJ(4) * t606 - t599;
t463 = t530 * t486;
t518 = qJ(4) * t609;
t626 = t531 * t535;
t551 = -pkin(7) * t626 + pkin(8) * t627;
t408 = qJD(1) * t551 + t463 + t518;
t602 = qJD(3) * t534;
t644 = t534 * t390 + t537 * t408 - t497 * t597 + t498 * t598 - t530 * t602 - t531 * t601;
t643 = qJD(1) * t606;
t573 = t533 * t371 + t536 * t372;
t542 = qJD(6) * t556 - t573;
t642 = t506 ^ 2;
t473 = t477 ^ 2;
t342 = t537 * t370 - t376 * t534;
t338 = -pkin(9) * t554 + t342;
t336 = pkin(5) * t517 + t338;
t638 = t336 * t536;
t637 = t339 * t536;
t636 = t407 * t535;
t635 = t410 * t535;
t633 = t458 * t531;
t630 = t517 * t538;
t628 = t530 * t537;
t539 = qJD(2) ^ 2;
t624 = t535 * t539;
t623 = t538 * t539;
t540 = qJD(1) ^ 2;
t622 = t538 * t540;
t483 = -t531 * t534 + t628;
t421 = t536 * t482 + t483 * t533;
t621 = -qJD(6) * t421 - t533 * t615 + t536 * t616;
t422 = -t482 * t533 + t483 * t536;
t620 = qJD(6) * t422 + t533 * t616 + t536 * t615;
t617 = pkin(5) * t615 + t569;
t582 = t538 * t594;
t612 = -qJ(4) * t582 - qJD(4) * t626;
t605 = qJD(2) * t538;
t604 = qJD(3) * t477;
t603 = qJD(3) * t531;
t595 = qJD(6) * t536;
t491 = -t531 * pkin(3) - t530 * qJ(4) - pkin(2);
t591 = pkin(7) * t606;
t590 = t536 * t371 - t533 * t372 - t555 * t595;
t587 = pkin(3) * t530 + pkin(7);
t579 = MDP(23) * t609;
t577 = -t492 - t532;
t576 = pkin(1) * t653;
t383 = qJD(2) * t544 - t633;
t449 = t530 * t458;
t384 = qJD(2) * t551 + t449 + t645;
t572 = t537 * t383 - t384 * t534;
t571 = t537 * t413 - t423 * t534;
t570 = t537 * t497 - t498 * t534;
t443 = t493 * t531 - t514;
t464 = t531 * pkin(4) - t491;
t441 = pkin(3) * t585 - t613;
t568 = t441 + t600;
t567 = qJD(6) * t336 + t331;
t516 = pkin(7) * t580;
t393 = pkin(3) * t502 - qJ(4) * t566 - t477 * qJD(4) + t516;
t565 = t530 * t641 - pkin(7);
t563 = t535 * t588;
t396 = -pkin(9) * t482 + t614;
t561 = -pkin(5) * t609 + pkin(9) * t616 + qJD(5) * t614 + qJD(6) * t396 - t408 * t534 + t531 * t602 - t646;
t395 = -pkin(9) * t483 + t570;
t560 = pkin(9) * t615 - qJD(6) * t395 + t644;
t329 = t336 * t533 + t637;
t454 = t482 * t535;
t347 = pkin(5) * t538 - pkin(9) * t454 + t571;
t453 = t534 * t626 - t535 * t628;
t348 = -pkin(9) * t453 + t618;
t558 = t347 * t533 + t348 * t536;
t391 = t536 * t453 + t454 * t533;
t392 = -t453 * t533 + t454 * t536;
t553 = t533 * t537 + t534 * t536;
t552 = t533 * t534 - t536 * t537;
t438 = -pkin(7) * t584 + t463;
t432 = -t531 * t591 + t449;
t548 = t534 * t383 + t537 * t384 + t413 * t597 - t423 * t598;
t546 = t554 * t596 - t590;
t512 = qJ(4) * t626;
t434 = t535 * t565 + t512;
t406 = t565 * t605 - t612;
t377 = -pkin(4) * t502 - t393;
t501 = qJD(3) * t585;
t456 = t475 * t608;
t455 = t475 * t603;
t450 = t535 * t587 - t512;
t437 = pkin(7) * t586 + t632;
t436 = -t443 + t527;
t431 = t530 * t591 + t633;
t430 = pkin(5) * t482 + t464;
t429 = qJD(1) * t563 - t632;
t428 = t438 + t518;
t427 = t587 * t605 + t612;
t415 = qJD(2) * t563 - t633;
t403 = t432 + t645;
t399 = qJD(5) * t483 * t535 + qJD(2) * t550;
t398 = qJD(5) * t454 + t534 * t582 - t605 * t628;
t388 = pkin(5) * t453 + t434;
t387 = -pkin(3) * t581 - t404;
t375 = -qJD(1) * t599 + t589;
t350 = pkin(5) * t398 + t406;
t345 = pkin(5) * t372 + t377;
t341 = qJD(6) * t392 + t536 * t398 + t399 * t533;
t340 = -qJD(6) * t391 - t398 * t533 + t399 * t536;
t333 = -pkin(9) * t398 + t548;
t332 = -pkin(5) * t606 - pkin(9) * t399 - qJD(5) * t618 + t572;
t328 = -t339 * t533 + t638;
t1 = [(-t371 * t453 - t372 * t454 - t398 * t554 - t399 * t555) * MDP(20) + (t371 * t454 + t399 * t554) * MDP(19) + (t434 * t371 + t373 * t399 + t377 * t454 + t406 * t554 - t548 * t517 - t549 * t538) * MDP(25) + (t375 * t435 + t387 * t436 + t393 * t450 + t400 * t427 + t403 * t410 + t407 * t415) * MDP(18) + (t340 * t506 - t538 * t546) * MDP(28) + ((-t517 - t608) * MDP(23) + (-t506 - t608) * MDP(30) + (-(t347 * t536 - t348 * t533) * qJD(1) - t328) * MDP(31) + (-qJD(1) * t454 - t554) * MDP(21) + (qJD(1) * t453 + t555) * MDP(22) + (-qJD(1) * t392 + t556) * MDP(28) + (qJD(1) * t391 + t352) * MDP(29) + (-qJD(1) * t571 - t342) * MDP(24) + (qJD(1) * t558 + t329) * MDP(32) + (qJD(1) * t618 + t343) * MDP(25)) * t606 + (-t340 * t556 - t392 * t546) * MDP(26) + (-t388 * t546 + t337 * t538 + t349 * t340 + t345 * t392 - t350 * t556 + (-(-qJD(6) * t348 + t332) * t506 - t330 * t538) * t533 + (-(qJD(6) * t347 + t333) * t506 - t567 * t538) * t536) * MDP(32) + (-t340 * t352 + t341 * t556 + t391 * t546 + t392 * t542) * MDP(27) + ((t332 * t536 - t333 * t533) * t506 + t575 * t538 + t350 * t352 - t388 * t542 + t345 * t391 + t349 * t341 + (-t329 * t538 - t506 * t558) * qJD(6)) * MDP(31) + (-t341 * t506 + t538 * t542) * MDP(29) + (t572 * t517 + t574 * t538 + t406 * t555 + t434 * t372 + t377 * t453 + t373 * t398 + (-t343 * t538 - t517 * t618) * qJD(5)) * MDP(24) + 0.2e1 * t580 * t652 + (-t372 * t538 - t398 * t517) * MDP(22) + (t371 * t538 + t399 * t517) * MDP(21) + t651 * t653 + MDP(6) * t623 + ((qJD(1) * t432 + t405) * t538 + ((pkin(7) * t477 + t492 * t531) * t538 + (-t425 + (-t444 + 0.2e1 * t515) * qJD(1)) * t535) * qJD(2)) * MDP(12) + ((-qJD(1) * t431 - t404) * t538 + ((pkin(7) * t475 + t492 * t530) * t538 + (t424 + (t443 + 0.2e1 * t514) * qJD(1)) * t535) * qJD(2)) * MDP(11) + (-t431 * t477 - t432 * t475 + (-t404 * t531 - t405 * t530) * t535 + (-t424 * t531 - t425 * t530 + (-t443 * t531 - t444 * t530) * qJD(1)) * t605) * MDP(13) + (-t403 * t475 + t415 * t477 + (-t375 * t530 + t387 * t531) * t535 + (t407 * t531 - t410 * t530 + (-t435 * t530 + t436 * t531) * qJD(1)) * t605) * MDP(16) + (t404 * t443 + t405 * t444 + t424 * t431 + t425 * t432 + (t492 + t520) * pkin(7) * t605) * MDP(14) + (-pkin(7) * t623 + t535 * t576) * MDP(9) - MDP(7) * t624 + (pkin(7) * t624 + t538 * t576) * MDP(10) + (-t393 * t626 - t427 * t477 + (-qJD(1) * t403 - t375) * t538 + (-t400 * t625 + t635 + (t435 * t535 - t450 * t625) * qJD(1)) * qJD(2)) * MDP(17) + (t393 * t629 + t427 * t475 + (qJD(1) * t415 + t387) * t538 + (t400 * t627 - t636 + (-t436 * t535 + t450 * t627) * qJD(1)) * qJD(2)) * MDP(15); (-t615 * t517 + (qJD(2) * t482 - t555) * t609) * MDP(22) + (-t371 * t482 - t372 * t483 - t554 * t615 - t555 * t616) * MDP(20) + (t371 * t483 + t554 * t616) * MDP(19) + (t616 * t517 + (-qJD(2) * t483 + t554) * t609) * MDP(21) + (t621 * t506 + (-qJD(2) * t422 - t556) * t609) * MDP(28) + (-t422 * t546 - t556 * t621) * MDP(26) + (-t430 * t546 + t345 * t422 + (t533 * t561 + t536 * t560) * t506 - t617 * t556 + t621 * t349 + ((t395 * t533 + t396 * t536) * qJD(2) - t329) * t609) * MDP(32) + (-t430 * t542 + t345 * t421 + (t533 * t560 - t536 * t561) * t506 + t617 * t352 + t620 * t349 + (-(t395 * t536 - t396 * t533) * qJD(2) + t328) * t609) * MDP(31) + (-t352 * t621 + t421 * t546 + t422 * t542 + t556 * t620) * MDP(27) - t622 * t652 + (MDP(9) * t535 * t540 + MDP(10) * t622) * pkin(1) + ((-qJ(3) * t594 + t425) * t535 + (-t438 + (-t477 + t607) * pkin(7) + (qJD(3) + t577) * t531) * t538) * qJD(1) * MDP(12) + t540 * t651 + (t464 * t371 + t377 * t483 + t644 * t517 + t569 * t554 + t616 * t373 + (qJD(2) * t614 - t343) * t609) * MDP(25) + (t464 * t372 + t377 * t482 + (-t498 * t597 + (-qJD(5) * t497 + t408 - t603) * t534 + t646) * t517 + t569 * t555 + t615 * t373 + (-qJD(2) * t570 + t342) * t609) * MDP(24) + t506 * t578 + t517 * t579 + (-t424 * t437 - t425 * t438 + (-t424 * t530 + t425 * t531) * qJD(3) + (-t404 * t530 + t405 * t531) * qJ(3) + t577 * t521) * MDP(14) + (qJ(3) * t375 * t531 + t393 * t491 - t400 * t441 - t407 * t429 + (-t428 + t603) * t410 + (qJ(3) * t387 + qJD(3) * t407 - qJD(4) * t400) * t530) * MDP(18) + (t501 + ((-qJ(3) * t607 - t424) * t535 + (t437 + t577 * t530 + (-t475 - t594) * pkin(7)) * t538) * qJD(1)) * MDP(11) + (t437 * t477 + t438 * t475 - t455 + (t424 * t608 + t405) * t531 + (t425 * t608 - t404 + t604) * t530) * MDP(13) + (t428 * t475 - t429 * t477 - t455 + (-t407 * t608 + t375) * t531 + (t410 * t608 + t387 + t604) * t530) * MDP(16) + (-t620 * t506 + (qJD(2) * t421 - t352) * t609) * MDP(29) + (-t393 * t530 + t568 * t477 + (-t635 + t428 * t538 + (qJ(3) * t606 + (-qJD(2) * t491 - qJD(3) + t400) * t538) * t531) * qJD(1)) * MDP(17) + (-t393 * t531 + t501 - t568 * t475 + (t636 - t429 * t538 + (-t400 * t538 + (t491 * t538 - t639) * qJD(2)) * t530) * qJD(1)) * MDP(15); (t424 * t477 + t425 * t475 + t516) * MDP(14) + (-t407 * t477 + t410 * t475 + t393) * MDP(18) + (-t372 - t650) * MDP(24) + (-t371 + t649) * MDP(25) + (t542 + t657) * MDP(31) + (t546 + t658) * MDP(32) + (MDP(11) + MDP(15)) * (-t477 * t608 + t502) + (MDP(12) - MDP(17)) * (t456 + t566) + (MDP(13) + MDP(16)) * (-t475 ^ 2 - t473); (t475 * t477 - t581) * MDP(15) + (-t456 + t566) * MDP(16) + (-t529 * t540 - t473) * MDP(17) + (t400 * t477 + (-pkin(3) * t606 + t410 * t538) * qJD(1) - t404) * MDP(18) + (-t517 * t598 - t555 * t477 + (-t534 * t630 - t537 * t606) * qJD(1)) * MDP(24) + (-t517 * t597 - t554 * t477 + (t534 * t606 - t537 * t630) * qJD(1)) * MDP(25) + (-t477 * t352 + t552 * t643 - t553 * t642) * MDP(31) + (t477 * t556 + t552 * t642 + t553 * t643) * MDP(32); t554 * t555 * MDP(19) + (t554 ^ 2 - t555 ^ 2) * MDP(20) + (t371 + t649) * MDP(21) + (-t372 + t650) * MDP(22) - qJD(2) * t579 + (t343 * t517 - t373 * t554 + t543) * MDP(24) + (t342 * t517 + t373 * t555 - t549) * MDP(25) + (-t546 + t658) * MDP(28) + (t542 - t657) * MDP(29) + (-(-t338 * t533 - t637) * t506 - t329 * qJD(6) + (-t352 * t554 - t506 * t596 - t536 * t581) * pkin(5) + t655) * MDP(31) + ((-t339 * t506 - t330) * t533 + (t338 * t506 - t567) * t536 + (-t506 * t595 + t533 * t581 + t554 * t556) * pkin(5) + t656) * MDP(32) + t654; (t590 + t658) * MDP(28) + (-t573 - t657) * MDP(29) + (t329 * t506 + t655) * MDP(31) + (t328 * t506 - t533 * t330 - t536 * t331 + t656) * MDP(32) + (-MDP(28) * t634 + MDP(29) * t556 - MDP(31) * t329 - MDP(32) * t638) * qJD(6) + t654;];
tauc  = t1;
