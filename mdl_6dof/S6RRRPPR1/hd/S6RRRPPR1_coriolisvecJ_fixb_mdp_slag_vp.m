% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:23:15
% EndTime: 2019-03-09 15:23:29
% DurationCPUTime: 6.89s
% Computational Cost: add. (9637->470), mult. (25178->640), div. (0->0), fcn. (19073->10), ass. (0->232)
t552 = sin(qJ(3));
t553 = sin(qJ(2));
t620 = qJD(1) * t553;
t609 = t552 * t620;
t555 = cos(qJ(3));
t556 = cos(qJ(2));
t619 = qJD(1) * t556;
t610 = t555 * t619;
t499 = -t609 + t610;
t500 = -t552 * t619 - t555 * t620;
t549 = sin(pkin(10));
t650 = cos(pkin(10));
t471 = t650 * t499 + t500 * t549;
t465 = qJD(6) - t471;
t545 = qJD(2) + qJD(3);
t548 = sin(pkin(11));
t550 = cos(pkin(11));
t572 = t549 * t499 - t500 * t650;
t454 = -t550 * t545 + t548 * t572;
t554 = cos(qJ(6));
t664 = t554 * t454;
t551 = sin(qJ(6));
t514 = t548 * t554 + t550 * t551;
t663 = t465 * t514;
t630 = t554 * t550;
t632 = t548 * t551;
t513 = -t630 + t632;
t662 = t465 * t513;
t456 = t545 * t548 + t550 * t572;
t661 = t454 * t551 - t456 * t554;
t613 = qJD(1) * qJD(2);
t660 = -0.2e1 * t613;
t659 = MDP(5) * (t553 ^ 2 - t556 ^ 2);
t658 = t553 * MDP(4);
t654 = pkin(7) + pkin(8);
t524 = t654 * t553;
t517 = qJD(1) * t524;
t651 = qJD(2) * pkin(2);
t509 = -t517 + t651;
t611 = qJD(2) * t654;
t593 = qJD(1) * t611;
t510 = t553 * t593;
t657 = t555 * (qJD(3) * t509 - t510);
t516 = t552 * t556 + t553 * t555;
t656 = qJD(1) * t516;
t560 = t545 * t656;
t511 = t556 * t593;
t525 = t654 * t556;
t519 = qJD(1) * t525;
t618 = qJD(3) * t552;
t597 = -t552 * t511 - t519 * t618;
t397 = -qJ(4) * t560 + t499 * qJD(4) + t597 + t657;
t608 = t556 * t613;
t473 = qJD(3) * t610 - t545 * t609 + t555 * t608;
t507 = t555 * t519;
t582 = -t509 * t552 - t507;
t598 = t552 * t510 - t555 * t511;
t564 = qJD(3) * t582 + t598;
t561 = -qJ(4) * t473 + qJD(4) * t500 + t564;
t365 = t650 * t397 + t549 * t561;
t362 = qJD(5) * t545 + t365;
t437 = t473 * t549 + t560 * t650;
t438 = t473 * t650 - t549 * t560;
t542 = pkin(2) * t620;
t559 = pkin(3) * t560 + qJD(2) * t542;
t371 = t437 * pkin(4) - t438 * qJ(5) - qJD(5) * t572 + t559;
t342 = -t362 * t548 + t550 * t371;
t343 = t550 * t362 + t548 * t371;
t590 = -t342 * t548 + t343 * t550;
t596 = t517 * t552 - t507;
t649 = qJ(4) * t499;
t459 = t596 - t649;
t494 = t500 * qJ(4);
t503 = t552 * t519;
t622 = -t555 * t517 - t503;
t460 = t494 + t622;
t420 = t549 * t459 + t460 * t650;
t653 = pkin(3) * t500;
t424 = pkin(4) * t572 - qJ(5) * t471 - t653;
t421 = t424 + t542;
t382 = -t420 * t548 + t550 * t421;
t532 = t549 * t552 * pkin(2);
t617 = qJD(3) * t555;
t491 = t650 * pkin(2) * t617 - qJD(3) * t532;
t483 = qJD(5) + t491;
t601 = t483 * t548 + t382;
t383 = t550 * t420 + t548 * t421;
t600 = t483 * t550 - t383;
t604 = t650 * t552;
t625 = -t650 * t459 + t460 * t549 - (t549 * t555 + t604) * qJD(3) * pkin(2);
t599 = t555 * t509 - t503;
t452 = t494 + t599;
t655 = t437 * t514 - t465 * t662;
t652 = t550 * pkin(5);
t544 = t550 * pkin(9);
t364 = t397 * t549 - t650 * t561;
t466 = -qJ(4) * t516 - t524 * t555 - t525 * t552;
t515 = t552 * t553 - t555 * t556;
t581 = t524 * t552 - t525 * t555;
t467 = -qJ(4) * t515 - t581;
t432 = -t650 * t466 + t467 * t549;
t647 = t364 * t432;
t445 = pkin(3) * t545 + t452;
t453 = -t582 + t649;
t447 = t549 * t453;
t402 = t445 * t650 - t447;
t400 = -t545 * pkin(4) + qJD(5) - t402;
t646 = t400 * t471;
t412 = t456 * t551 + t664;
t645 = t412 * t572;
t644 = t661 * t572;
t642 = t438 * t548;
t641 = t438 * t550;
t477 = t545 * t515;
t570 = t516 * qJD(3);
t478 = qJD(2) * t516 + t570;
t440 = -t477 * t650 - t549 * t478;
t640 = t440 * t548;
t639 = t471 * t548;
t638 = t471 * t550;
t476 = -t549 * t515 + t516 * t650;
t637 = t476 * t548;
t636 = t476 * t550;
t541 = -pkin(2) * t556 - pkin(1);
t523 = t541 * qJD(1);
t633 = t523 * t500;
t557 = qJD(2) ^ 2;
t631 = t553 * t557;
t629 = t556 * t557;
t558 = qJD(1) ^ 2;
t628 = t556 * t558;
t518 = t553 * t611;
t520 = t556 * t611;
t571 = -t555 * t518 - t552 * t520 - t524 * t617 - t525 * t618;
t417 = -qJ(4) * t478 - qJD(4) * t515 + t571;
t563 = qJD(3) * t581 + t518 * t552 - t555 * t520;
t418 = qJ(4) * t477 - qJD(4) * t516 + t563;
t379 = t417 * t650 + t549 * t418;
t439 = -t477 * t549 + t478 * t650;
t543 = t553 * t651;
t606 = pkin(3) * t478 + t543;
t386 = pkin(4) * t439 - qJ(5) * t440 - qJD(5) * t476 + t606;
t346 = t550 * t379 + t548 * t386;
t605 = t650 * t453;
t403 = t549 * t445 + t605;
t401 = qJ(5) * t545 + t403;
t479 = -pkin(3) * t499 + qJD(4) + t523;
t410 = -pkin(4) * t471 - qJ(5) * t572 + t479;
t373 = t550 * t401 + t548 * t410;
t408 = t452 * t650 - t447;
t381 = t550 * t408 + t548 * t424;
t475 = t515 * t650 + t516 * t549;
t579 = pkin(3) * t515 + t541;
t431 = pkin(4) * t475 - qJ(5) * t476 + t579;
t433 = t549 * t466 + t467 * t650;
t390 = t548 * t431 + t550 * t433;
t614 = qJD(6) * t554;
t626 = t438 * t630 - t454 * t614;
t462 = pkin(5) * t639;
t624 = -t462 - t625;
t623 = t491 - t420;
t540 = pkin(2) * t555 + pkin(3);
t493 = pkin(2) * t604 + t549 * t540;
t356 = -pkin(9) * t454 + t373;
t616 = qJD(6) * t356;
t615 = qJD(6) * t476;
t612 = pkin(9) * t639;
t607 = -pkin(2) * t545 - t509;
t603 = pkin(1) * t660;
t602 = t364 * t548 + t373 * t572;
t345 = -t379 * t548 + t550 * t386;
t372 = -t401 * t548 + t550 * t410;
t380 = -t408 * t548 + t550 * t424;
t378 = t417 * t549 - t650 * t418;
t389 = t550 * t431 - t433 * t548;
t407 = t452 * t549 + t605;
t339 = -pkin(9) * t642 + t343;
t352 = -pkin(5) * t471 - pkin(9) * t456 + t372;
t594 = -qJD(6) * t352 - t339;
t592 = pkin(5) * t572 - pkin(9) * t638;
t538 = -pkin(3) * t650 - pkin(4);
t591 = -t513 * t437 - t663 * t465;
t492 = t540 * t650 - t532;
t336 = t352 * t554 - t356 * t551;
t337 = t352 * t551 + t356 * t554;
t589 = -t364 * t550 - t372 * t572;
t588 = t364 * t476 + t432 * t438;
t368 = pkin(5) * t475 - pkin(9) * t636 + t389;
t376 = -pkin(9) * t637 + t390;
t587 = t368 * t554 - t376 * t551;
t586 = t368 * t551 + t376 * t554;
t585 = t372 * t548 - t373 * t550;
t584 = t402 * t471 + t403 * t572;
t583 = t454 * t550 - t456 * t548;
t580 = -qJD(6) * t456 - t642;
t578 = t372 * t638 + t373 * t639 + t590;
t487 = -pkin(4) - t492;
t577 = -t523 * t499 - t597;
t486 = qJ(5) + t493;
t481 = t486 * t550 + t544;
t576 = qJD(6) * t481 + t592 + t601;
t480 = (-pkin(9) - t486) * t548;
t575 = -qJD(6) * t480 - t600 - t612;
t536 = pkin(3) * t549 + qJ(5);
t502 = t536 * t550 + t544;
t574 = qJD(5) * t548 + qJD(6) * t502 + t380 + t592;
t501 = (-pkin(9) - t536) * t548;
t573 = -qJD(5) * t550 - qJD(6) * t501 + t381 - t612;
t351 = pkin(5) * t642 + t364;
t391 = t454 * pkin(5) + t400;
t569 = -t336 * t572 + t351 * t513 + t663 * t391;
t568 = t337 * t572 + t351 * t514 - t662 * t391;
t567 = t400 * t440 + t588;
t566 = -t437 * t486 + t438 * t487 + t471 * t483 - t646;
t565 = qJD(5) * t471 - t437 * t536 + t438 * t538 - t646;
t357 = t580 * t551 + t626;
t358 = -qJD(6) * t661 + t438 * t514;
t562 = t500 * t499 * MDP(11) - t465 * t572 * MDP(28) + (-t357 * t513 - t358 * t514 + t412 * t662 + t661 * t663) * MDP(25) + (t357 * t514 + t661 * t662) * MDP(24) + (t644 + t655) * MDP(26) + (t591 + t645) * MDP(27) + (-t499 * t545 + t473) * MDP(13) + (-t500 * t545 - t560) * MDP(14) + (-t499 ^ 2 + t500 ^ 2) * MDP(12);
t521 = t538 - t652;
t482 = t487 - t652;
t436 = t513 * t476;
t435 = t514 * t476;
t398 = pkin(5) * t637 + t432;
t392 = t407 + t462;
t375 = t440 * t514 + t614 * t636 - t615 * t632;
t374 = -t440 * t513 - t514 * t615;
t355 = pkin(5) * t640 + t378;
t344 = -pkin(9) * t640 + t346;
t340 = pkin(5) * t439 - t440 * t544 + t345;
t335 = pkin(5) * t437 - pkin(9) * t641 + t342;
t334 = t554 * t335;
t1 = [(-t345 * t456 - t346 * t454 + (-t342 * t476 - t372 * t440 - t389 * t438) * t550 + (-t343 * t476 - t373 * t440 - t390 * t438) * t548) * MDP(22) + (t342 * t475 - t345 * t471 + t372 * t439 + t378 * t454 + t389 * t437 + t548 * t567) * MDP(20) + (-t343 * t475 + t346 * t471 - t373 * t439 + t378 * t456 - t390 * t437 + t550 * t567) * MDP(21) + (-t473 * t515 - t477 * t499 + t500 * t478 - t516 * t560) * MDP(12) + (t541 * t473 - t523 * t477 + (-t500 + t656) * t543) * MDP(17) + (-t499 * t543 + t523 * t478 + (t541 * t570 + (t553 * pkin(2) * t515 + t516 * t541) * qJD(2)) * qJD(1)) * MDP(16) + t659 * t660 + 0.2e1 * t608 * t658 + MDP(6) * t629 + (t365 * t433 - t402 * t378 + t403 * t379 + t479 * t606 + t559 * t579 + t647) * MDP(19) + (t342 * t389 + t343 * t390 + t345 * t372 + t346 * t373 + t378 * t400 + t647) * MDP(23) + (-t365 * t475 + t378 * t572 + t379 * t471 - t402 * t440 - t403 * t439 - t433 * t437 + t588) * MDP(18) + ((t340 * t554 - t344 * t551) * t465 + t587 * t437 + (-t339 * t551 + t334) * t475 + t336 * t439 + t355 * t412 + t398 * t358 + t351 * t435 + t391 * t375 + (-t337 * t475 - t465 * t586) * qJD(6)) * MDP(29) + (-(t340 * t551 + t344 * t554) * t465 - t586 * t437 - (t335 * t551 + t339 * t554) * t475 - t337 * t439 - t355 * t661 + t398 * t357 - t351 * t436 + t391 * t374 + (-t336 * t475 - t465 * t587) * qJD(6)) * MDP(30) - MDP(7) * t631 + (pkin(7) * t631 + t556 * t603) * MDP(10) + (-pkin(7) * t629 + t553 * t603) * MDP(9) + (t473 * t516 + t477 * t500) * MDP(11) + (-t358 * t475 - t375 * t465 - t412 * t439 - t435 * t437) * MDP(27) + (t357 * t475 + t374 * t465 - t436 * t437 - t439 * t661) * MDP(26) + (t437 * t475 + t439 * t465) * MDP(28) + (-t357 * t435 + t358 * t436 - t374 * t412 + t375 * t661) * MDP(25) + (-t357 * t436 - t374 * t661) * MDP(24) + (-t477 * MDP(13) - t478 * MDP(14) + MDP(16) * t563 - MDP(17) * t571) * t545; (t364 * t487 - t372 * t601 + t373 * t600 - t400 * t625 + t486 * t590) * MDP(23) + ((t480 * t554 - t481 * t551) * t437 + t482 * t358 + (t551 * t575 - t554 * t576) * t465 + t624 * t412 + t569) * MDP(29) + (-t437 * t493 - t438 * t492 + t471 * t623 - t572 * t625 + t584) * MDP(18) + (t500 * t542 + t622 * t545 + (qJD(3) * t607 + t510) * t555 + t577) * MDP(17) + (t499 * t542 + t633 - t596 * t545 + (t552 * t607 - t507) * qJD(3) + t598) * MDP(16) + (t382 * t471 - t454 * t625 + t548 * t566 + t589) * MDP(20) + (-t383 * t471 - t456 * t625 + t550 * t566 + t602) * MDP(21) + (t365 * t493 - t364 * t492 - t479 * (t542 - t653) + t623 * t403 + t625 * t402) * MDP(19) + t558 * t659 + t562 + (-t454 * t600 + t456 * t601 + t578) * MDP(22) + (-(t480 * t551 + t481 * t554) * t437 + t482 * t357 + (t551 * t576 + t554 * t575) * t465 - t624 * t661 + t568) * MDP(30) - t628 * t658 + (MDP(9) * t553 * t558 + MDP(10) * t628) * pkin(1); ((t501 * t554 - t502 * t551) * t437 + t521 * t358 - t392 * t412 + (t551 * t573 - t554 * t574) * t465 + t569) * MDP(29) + (t380 * t471 - t407 * t454 + t548 * t565 + t589) * MDP(20) + (-t381 * t471 - t407 * t456 + t550 * t565 + t602) * MDP(21) + (-qJD(5) * t585 + t364 * t538 - t372 * t380 - t373 * t381 - t400 * t407 + t536 * t590) * MDP(23) + (-t407 * t572 - t408 * t471 + (-t437 * t549 - t438 * t650) * pkin(3) + t584) * MDP(18) + t562 + (-t545 * t582 + t564 + t633) * MDP(16) + (-qJD(5) * t583 + t380 * t456 + t381 * t454 + t578) * MDP(22) + (t545 * t599 + t577 - t657) * MDP(17) + (t402 * t407 - t403 * t408 + (-t364 * t650 + t365 * t549 + t479 * t500) * pkin(3)) * MDP(19) + (-(t501 * t551 + t502 * t554) * t437 + t521 * t357 + t392 * t661 + (t551 * t574 + t554 * t573) * t465 + t568) * MDP(30); (-t471 ^ 2 - t572 ^ 2) * MDP(18) + (t402 * t572 - t403 * t471 + t559) * MDP(19) + (t437 * t550 - t454 * t572 - t471 * t639) * MDP(20) + (-t437 * t548 - t456 * t572 - t471 * t638) * MDP(21) + (t583 * t471 + (-t548 ^ 2 - t550 ^ 2) * t438) * MDP(22) + (t342 * t550 + t343 * t548 - t400 * t572 + t471 * t585) * MDP(23) + (t591 - t645) * MDP(29) + (t644 - t655) * MDP(30); (-t456 * t471 + t642) * MDP(20) + (t454 * t471 + t641) * MDP(21) + (-t454 ^ 2 - t456 ^ 2) * MDP(22) + (t372 * t456 + t373 * t454 + t364) * MDP(23) + (-t661 * t465 + t358) * MDP(29) + (-t465 * t664 + (-t456 * t465 + t580) * t551 + t626) * MDP(30); -t412 ^ 2 * MDP(25) + (t412 * t465 + t626) * MDP(26) + t437 * MDP(28) + (t337 * t465 + t334) * MDP(29) + (t336 * t465 + t391 * t412) * MDP(30) - (MDP(24) * t412 - MDP(25) * t661 + t465 * MDP(27) - t391 * MDP(29)) * t661 + (MDP(27) * t580 - MDP(29) * t616 + MDP(30) * t594) * t554 + (t580 * MDP(26) + (qJD(6) * t454 - t641) * MDP(27) + t594 * MDP(29) + (-t335 + t616) * MDP(30)) * t551;];
tauc  = t1;
