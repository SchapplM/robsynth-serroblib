% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRPRP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:41:43
% EndTime: 2019-03-09 16:41:55
% DurationCPUTime: 7.51s
% Computational Cost: add. (10483->504), mult. (25356->637), div. (0->0), fcn. (18512->8), ass. (0->218)
t592 = sin(qJ(2));
t696 = cos(qJ(3));
t641 = qJD(1) * t696;
t591 = sin(qJ(3));
t594 = cos(qJ(2));
t674 = t591 * t594;
t544 = -qJD(1) * t674 - t592 * t641;
t585 = qJD(2) + qJD(3);
t588 = sin(pkin(10));
t589 = cos(pkin(10));
t521 = t544 * t589 - t585 * t588;
t590 = sin(qJ(5));
t593 = cos(qJ(5));
t530 = t588 * t544;
t704 = t585 * t589 + t530;
t628 = t593 * t704;
t466 = t521 * t590 + t628;
t720 = t466 ^ 2;
t656 = qJD(1) * t592;
t717 = -t591 * t656 + t594 * t641;
t537 = qJD(5) - t717;
t719 = t466 * t537;
t504 = t717 * t585;
t714 = -t521 * t593 + t590 * t704;
t698 = t714 ^ 2;
t557 = t588 * t593 + t589 * t590;
t541 = t557 * qJD(5);
t715 = -t557 * t717 + t541;
t672 = t593 * t589;
t556 = t588 * t590 - t672;
t653 = qJD(5) * t593;
t654 = qJD(5) * t590;
t702 = -t588 * t654 + t589 * t653;
t659 = -t556 * t717 - t702;
t651 = qJD(1) * qJD(2);
t713 = -0.2e1 * t651;
t712 = MDP(5) * (t592 ^ 2 - t594 ^ 2);
t711 = t504 * t557;
t559 = t592 * t696 + t674;
t517 = t585 * t559;
t505 = t517 * qJD(1);
t639 = t592 * t651;
t433 = pkin(2) * t639 + pkin(3) * t505 - qJ(4) * t504 + qJD(4) * t544;
t697 = -pkin(8) - pkin(7);
t569 = t697 * t592;
t561 = qJD(1) * t569;
t692 = qJD(2) * pkin(2);
t551 = t561 + t692;
t646 = qJD(2) * t697;
t633 = qJD(1) * t646;
t552 = t592 * t633;
t553 = t594 * t633;
t570 = t697 * t594;
t563 = qJD(1) * t570;
t640 = qJD(3) * t696;
t655 = qJD(3) * t591;
t601 = t551 * t640 + t552 * t696 + t591 * t553 + t563 * t655;
t450 = t585 * qJD(4) + t601;
t400 = t589 * t433 - t450 * t588;
t401 = t588 * t433 + t589 * t450;
t626 = -t400 * t588 + t401 * t589;
t583 = -pkin(2) * t594 - pkin(1);
t615 = -t591 * t592 + t594 * t696;
t509 = -pkin(3) * t615 - qJ(4) * t559 + t583;
t525 = t591 * t569 - t570 * t696;
t458 = t589 * t509 - t525 * t588;
t584 = t589 * pkin(9);
t435 = -pkin(4) * t615 - t559 * t584 + t458;
t459 = t588 * t509 + t589 * t525;
t678 = t559 * t588;
t445 = -pkin(9) * t678 + t459;
t710 = t590 * t435 + t593 * t445;
t506 = -pkin(3) * t544 - qJ(4) * t717;
t489 = pkin(2) * t656 + t506;
t547 = t591 * t563;
t514 = t561 * t696 + t547;
t448 = t589 * t489 - t514 * t588;
t679 = t717 * t589;
t632 = -t544 * pkin(4) - pkin(9) * t679;
t420 = t448 + t632;
t449 = t588 * t489 + t589 * t514;
t680 = t717 * t588;
t650 = pkin(9) * t680;
t434 = -t650 + t449;
t573 = pkin(2) * t640 + qJD(4);
t579 = pkin(2) * t591 + qJ(4);
t549 = (-pkin(9) - t579) * t588;
t550 = t579 * t589 + t584;
t620 = t549 * t593 - t550 * t590;
t709 = -qJD(5) * t620 + t590 * t420 + t593 * t434 + t556 * t573;
t501 = t549 * t590 + t550 * t593;
t708 = -qJD(5) * t501 - t420 * t593 + t434 * t590 - t557 * t573;
t510 = t551 * t696 + t547;
t451 = t589 * t506 - t510 * t588;
t423 = t451 + t632;
t452 = t588 * t506 + t589 * t510;
t436 = -t650 + t452;
t566 = (-pkin(9) - qJ(4)) * t588;
t567 = qJ(4) * t589 + t584;
t619 = t566 * t593 - t567 * t590;
t707 = qJD(4) * t556 - qJD(5) * t619 + t590 * t423 + t593 * t436;
t520 = t566 * t590 + t567 * t593;
t706 = -qJD(4) * t557 - qJD(5) * t520 - t423 * t593 + t436 * t590;
t548 = t696 * t563;
t513 = t591 * t561 - t548;
t705 = -pkin(2) * t655 + t513;
t528 = pkin(4) * t680;
t703 = -pkin(5) * t715 - qJ(6) * t659 + qJD(6) * t557 + t528;
t701 = t696 * t569 + t591 * t570;
t700 = qJD(1) * t559;
t684 = t504 * t588;
t405 = t590 * (-qJD(5) * t521 + t684) - qJD(5) * t628 - t504 * t672;
t516 = t585 * t615;
t648 = t592 * t692;
t446 = pkin(3) * t517 - qJ(4) * t516 - qJD(4) * t559 + t648;
t562 = t592 * t646;
t564 = t594 * t646;
t467 = qJD(3) * t701 + t696 * t562 + t591 * t564;
t403 = t589 * t446 - t467 * t588;
t392 = pkin(4) * t517 - t516 * t584 + t403;
t404 = t588 * t446 + t589 * t467;
t683 = t516 * t588;
t402 = -pkin(9) * t683 + t404;
t699 = -qJD(5) * t710 + t392 * t593 - t402 * t590;
t695 = pkin(5) * t505;
t694 = pkin(5) * t544;
t693 = t589 * pkin(4);
t691 = qJ(6) * t505;
t568 = t583 * qJD(1);
t484 = -pkin(3) * t717 + qJ(4) * t544 + t568;
t511 = t591 * t551 - t548;
t492 = qJ(4) * t585 + t511;
t438 = t589 * t484 - t492 * t588;
t411 = -pkin(4) * t717 + pkin(9) * t521 + t438;
t439 = t588 * t484 + t589 * t492;
t418 = pkin(9) * t704 + t439;
t380 = t411 * t590 + t418 * t593;
t690 = t380 * t537;
t688 = t714 * t466;
t686 = t620 * t505;
t685 = t501 * t505;
t682 = t619 * t505;
t681 = t520 * t505;
t675 = t589 * t504;
t595 = qJD(2) ^ 2;
t673 = t592 * t595;
t671 = t594 * t595;
t596 = qJD(1) ^ 2;
t670 = t594 * t596;
t534 = t544 * qJ(6);
t669 = -t534 + t709;
t668 = t694 + t708;
t667 = -t534 + t707;
t666 = t694 + t706;
t665 = t511 + t703;
t664 = t703 + t705;
t379 = t411 * t593 - t418 * t590;
t652 = qJD(6) - t379;
t649 = t696 * pkin(2);
t580 = -pkin(3) - t693;
t638 = pkin(1) * t713;
t454 = t551 * t655 + t591 * t552 - t696 * t553 - t563 * t640;
t637 = -t439 * t544 + t454 * t588;
t636 = t573 * t588 + t448;
t378 = pkin(4) * t505 - pkin(9) * t675 + t400;
t387 = -pkin(9) * t684 + t401;
t634 = -t593 * t378 + t590 * t387 + t411 * t654 + t418 * t653;
t582 = -t649 - pkin(3);
t631 = -t528 - t705;
t424 = pkin(4) * t684 + t454;
t629 = t589 * t704;
t623 = t435 * t593 - t445 * t590;
t621 = t438 * t544 - t454 * t589;
t617 = t438 * t679 + t439 * t680 + t626;
t486 = pkin(4) * t678 - t701;
t488 = -t585 * pkin(3) + qJD(4) - t510;
t455 = -pkin(4) * t704 + t488;
t397 = -pkin(5) * t466 - qJ(6) * t714 + t455;
t614 = t397 * t714 + t634;
t613 = t544 * t568 - t454;
t468 = t591 * t562 - t564 * t696 + t569 * t655 - t570 * t640;
t612 = t590 * t378 + t593 * t387 + t411 * t653 - t418 * t654;
t611 = t590 * t392 + t593 * t402 + t435 * t653 - t445 * t654;
t406 = qJD(5) * t714 + t711;
t369 = pkin(5) * t406 + qJ(6) * t405 - qJD(6) * t714 + t424;
t374 = -pkin(5) * t537 + t652;
t610 = t369 * t556 - t374 * t544 + t397 * t715;
t375 = qJ(6) * t537 + t380;
t609 = -t369 * t557 + t375 * t544 + t397 * t659;
t608 = t379 * t544 + t424 * t556 + t455 * t715;
t607 = -t380 * t544 + t424 * t557 - t455 * t659;
t437 = pkin(4) * t683 + t468;
t606 = t454 * t559 + t488 * t516 - t504 * t701;
t503 = t556 * pkin(5) - t557 * qJ(6) + t580;
t605 = -pkin(3) * t504 - qJ(4) * t505 - (-qJD(4) + t488) * t717;
t361 = qJD(6) * t537 + t612 + t691;
t363 = t634 - t695;
t604 = -t361 * t556 + t363 * t557 - t374 * t659 - t375 * t715;
t603 = t582 * t504 - t579 * t505 - (t488 - t573) * t717;
t598 = t585 * t700;
t600 = (t405 * t556 - t406 * t557 - t466 * t659 - t714 * t715) * MDP(23) + (-t405 * t557 - t659 * t714) * MDP(22) + (t505 * t557 - t537 * t659 + t544 * t714) * MDP(24) + (t466 * t544 - t505 * t556 - t537 * t715) * MDP(25) + (-t544 * t585 - t598) * MDP(14) + (t544 ^ 2 - t717 ^ 2) * MDP(12) + (MDP(11) * t717 + MDP(26) * t537) * t544;
t599 = -t568 * t717 - t601;
t565 = t582 - t693;
t496 = t556 * t559;
t495 = t557 * t559;
t490 = -t649 + t503;
t472 = t528 + t511;
t422 = t516 * t557 + t559 * t702;
t421 = t516 * t556 + t541 * t559;
t416 = t495 * pkin(5) + t496 * qJ(6) + t486;
t414 = pkin(5) * t714 - qJ(6) * t466;
t396 = pkin(5) * t615 - t623;
t395 = -qJ(6) * t615 + t710;
t388 = -t405 - t719;
t370 = t422 * pkin(5) + t421 * qJ(6) + t496 * qJD(6) + t437;
t366 = -pkin(5) * t517 - t699;
t365 = qJ(6) * t517 - qJD(6) * t615 + t611;
t1 = [(t405 * t615 - t421 * t537 - t496 * t505 + t517 * t714) * MDP(24) + (-t361 * t615 + t365 * t537 + t369 * t496 - t370 * t714 + t375 * t517 + t395 * t505 + t397 * t421 + t405 * t416) * MDP(31) + (-t380 * t517 - t486 * t405 - t455 * t421 - t424 * t496 + t437 * t714 - t505 * t710 - t537 * t611 + t612 * t615) * MDP(28) + (t405 * t496 - t421 * t714) * MDP(22) + (t400 * t458 + t401 * t459 + t403 * t438 + t404 * t439 - t454 * t701 + t468 * t488) * MDP(21) + (t403 * t521 + t404 * t530 + (-t400 * t559 + t404 * t585 - t438 * t516 - t458 * t504) * t589 + (-t401 * t559 - t439 * t516 - t459 * t504) * t588) * MDP(20) + (-t505 * t615 + t517 * t537) * MDP(26) + (-t467 * t585 + t504 * t583 + t516 * t568 + (-t544 + t700) * t648) * MDP(17) + (-pkin(7) * t671 + t592 * t638) * MDP(9) - MDP(7) * t673 + (pkin(7) * t673 + t594 * t638) * MDP(10) + 0.2e1 * t594 * MDP(4) * t639 + t516 * t585 * MDP(13) - t517 * t585 * MDP(14) + (t504 * t559 - t516 * t544) * MDP(11) + (t405 * t495 + t406 * t496 - t421 * t466 - t422 * t714) * MDP(23) + (-t361 * t495 - t363 * t496 + t365 * t466 + t366 * t714 - t374 * t421 - t375 * t422 - t395 * t406 - t396 * t405) * MDP(30) + (t379 * t517 + t486 * t406 + t455 * t422 + t424 * t495 - t437 * t466 + t623 * t505 + t537 * t699 + t615 * t634) * MDP(27) + (t406 * t615 - t422 * t537 + t466 * t517 - t495 * t505) * MDP(25) + (t363 * t615 - t366 * t537 + t369 * t495 - t370 * t466 - t374 * t517 - t396 * t505 + t397 * t422 + t406 * t416) * MDP(29) + (t401 * t615 + t404 * t717 - t439 * t517 - t459 * t505 - t468 * t521 + t589 * t606) * MDP(19) + (-t468 * t585 + t505 * t583 + t517 * t568 + (-qJD(1) * t615 - t717) * t648) * MDP(16) + (t504 * t615 - t505 * t559 + t516 * t717 + t517 * t544) * MDP(12) + (-t400 * t615 - t403 * t717 + t438 * t517 + t458 * t505 - t468 * t704 + t588 * t606) * MDP(18) + t712 * t713 + MDP(6) * t671 + (t361 * t395 + t363 * t396 + t365 * t375 + t366 * t374 + t369 * t416 + t370 * t397) * MDP(32); (t513 * t585 + (-t585 * t655 + t656 * t717) * pkin(2) + t613) * MDP(16) + (t454 * t582 + t626 * t579 - t705 * t488 + (t573 * t589 - t449) * t439 - t636 * t438) * MDP(21) - t592 * MDP(4) * t670 + (-t449 * t717 + t521 * t705 + t589 * t603 + t637) * MDP(19) + (t406 * t490 + t466 * t664 + t537 * t668 + t610 + t686) * MDP(29) + (t405 * t490 - t537 * t669 + t664 * t714 + t609 + t685) * MDP(31) + (t405 * t620 - t406 * t501 - t466 * t669 - t668 * t714 + t604) * MDP(30) + (-t565 * t405 + t537 * t709 + t631 * t714 + t607 - t685) * MDP(28) + (-t449 * t704 - t521 * t636 + t573 * t629 + t617) * MDP(20) + (t448 * t717 + t588 * t603 + t704 * t705 + t621) * MDP(18) + t596 * t712 + (t565 * t406 - t466 * t631 + t537 * t708 + t608 + t686) * MDP(27) + (t361 * t501 - t363 * t620 + t369 * t490 - t374 * t668 - t375 * t669 - t397 * t664) * MDP(32) + t600 + (t514 * t585 + (t544 * t656 - t585 * t640) * pkin(2) + t599) * MDP(17) + (MDP(9) * t592 * t596 + MDP(10) * t670) * pkin(1); (t361 * t520 - t363 * t619 + t369 * t503 - t374 * t666 - t375 * t667 - t397 * t665) * MDP(32) + (-pkin(3) * t454 - t438 * t451 - t439 * t452 - t488 * t511 + (-t438 * t588 + t439 * t589) * qJD(4) + t626 * qJ(4)) * MDP(21) + (t405 * t503 - t537 * t667 + t665 * t714 + t609 + t681) * MDP(31) + (t405 * t619 - t406 * t520 - t466 * t667 - t666 * t714 + t604) * MDP(30) + (-t580 * t405 - t472 * t714 + t537 * t707 + t607 - t681) * MDP(28) + (t451 * t717 + t511 * t704 + t588 * t605 + t621) * MDP(18) + (-t452 * t704 - t451 * t521 + (-t521 * t588 + t629) * qJD(4) + t617) * MDP(20) + (-t452 * t717 + t511 * t521 + t589 * t605 + t637) * MDP(19) + (t511 * t585 + t613) * MDP(16) + (t406 * t503 + t466 * t665 + t537 * t666 + t610 + t682) * MDP(29) + t600 + (t580 * t406 + t466 * t472 + t537 * t706 + t608 + t682) * MDP(27) + (t510 * t585 + t599) * MDP(17); (t521 * t717 + t684) * MDP(18) + (-t704 * t717 + t675) * MDP(19) + (-t521 ^ 2 - t704 ^ 2) * MDP(20) + (-t438 * t521 - t439 * t704 + t454) * MDP(21) + (-t698 - t720) * MDP(30) + (-t374 * t714 - t375 * t466 + t369) * MDP(32) + (MDP(27) + MDP(29)) * (t537 * t714 + t406) + (-MDP(28) + MDP(31)) * (t405 - t719); -MDP(22) * t688 + (t698 - t720) * MDP(23) + t388 * MDP(24) + (-t711 + (-qJD(5) + t537) * t714) * MDP(25) + t505 * MDP(26) + (-t455 * t714 - t634 + t690) * MDP(27) + (t379 * t537 - t455 * t466 - t612) * MDP(28) + (t414 * t466 - t614 + t690 + 0.2e1 * t695) * MDP(29) + (pkin(5) * t405 - qJ(6) * t406 + (t375 - t380) * t714 - (t374 - t652) * t466) * MDP(30) + (0.2e1 * t691 + t397 * t466 + t414 * t714 + (0.2e1 * qJD(6) - t379) * t537 + t612) * MDP(31) + (-pkin(5) * t363 + qJ(6) * t361 - t374 * t380 + t375 * t652 - t397 * t414) * MDP(32); t388 * MDP(30) + (-t537 ^ 2 - t698) * MDP(31) + (-t375 * t537 + t614 - t695) * MDP(32) + (-t688 - t598) * MDP(29);];
tauc  = t1;
