% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR10_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:18
% EndTime: 2019-03-09 09:37:33
% DurationCPUTime: 11.80s
% Computational Cost: add. (5684->624), mult. (12089->786), div. (0->0), fcn. (8530->14), ass. (0->268)
t626 = sin(qJ(2));
t721 = qJD(1) * t626;
t587 = qJD(5) + t721;
t576 = qJD(6) + t587;
t622 = cos(pkin(10));
t630 = cos(qJ(2));
t720 = qJD(1) * t630;
t697 = t622 * t720;
t621 = sin(pkin(10));
t719 = qJD(2) * t621;
t552 = t697 + t719;
t695 = t621 * t720;
t718 = qJD(2) * t622;
t554 = -t695 + t718;
t625 = sin(qJ(5));
t629 = cos(qJ(5));
t481 = t629 * t552 + t554 * t625;
t628 = cos(qJ(6));
t480 = t552 * t625 - t554 * t629;
t624 = sin(qJ(6));
t745 = t480 * t624;
t773 = -t628 * t481 + t745;
t772 = t576 * t773;
t557 = t621 * t629 + t622 * t625;
t540 = t557 * qJD(5);
t652 = t626 * t557;
t727 = -qJD(1) * t652 - t540;
t759 = pkin(3) + pkin(7);
t771 = t480 * t587;
t770 = t481 * t587;
t663 = t480 * t628 + t481 * t624;
t769 = t576 * t663;
t698 = t622 * t721;
t713 = qJD(5) * t629;
t714 = qJD(5) * t625;
t740 = t621 * t625;
t726 = -t621 * t714 + t622 * t713 + t629 * t698 - t721 * t740;
t660 = -t622 * t629 + t740;
t598 = pkin(7) * t721;
t767 = qJD(3) + t598;
t607 = t626 * qJ(3);
t690 = -pkin(1) - t607;
t753 = pkin(2) + qJ(4);
t646 = -t630 * t753 + t690;
t513 = t646 * qJD(1);
t709 = pkin(3) * t721 + t767;
t523 = -qJD(2) * t753 + t709;
t450 = -t513 * t621 + t622 * t523;
t430 = pkin(4) * t721 - pkin(8) * t554 + t450;
t451 = t622 * t513 + t621 * t523;
t437 = -pkin(8) * t552 + t451;
t404 = t430 * t625 + t437 * t629;
t399 = -pkin(9) * t481 + t404;
t712 = qJD(6) * t624;
t397 = t399 * t712;
t599 = pkin(7) * t720;
t564 = pkin(3) * t720 + t599;
t618 = qJD(2) * qJ(3);
t535 = qJD(4) + t618 + t564;
t489 = pkin(4) * t552 + t535;
t431 = pkin(5) * t481 + t489;
t615 = pkin(10) + qJ(5);
t605 = qJ(6) + t615;
t593 = sin(t605);
t594 = cos(t605);
t627 = sin(qJ(1));
t631 = cos(qJ(1));
t736 = t626 * t631;
t505 = t593 * t736 + t594 * t627;
t737 = t626 * t627;
t507 = -t593 * t737 + t594 * t631;
t614 = g(3) * t630;
t766 = g(1) * t505 - g(2) * t507 - t431 * t773 - t593 * t614 + t397;
t504 = -t593 * t627 + t594 * t736;
t506 = t593 * t631 + t594 * t737;
t707 = qJD(1) * qJD(2);
t694 = t626 * t707;
t705 = qJDD(1) * t630;
t502 = qJDD(2) * t621 + (-t694 + t705) * t622;
t679 = qJDD(2) * t622 - t621 * t705;
t503 = t621 * t694 + t679;
t419 = -t625 * t502 + t629 * t503 - t552 * t713 - t554 * t714;
t693 = t630 * t707;
t706 = qJDD(1) * t626;
t649 = t693 + t706;
t556 = qJDD(5) + t649;
t586 = pkin(2) * t694;
t750 = qJ(3) * t630;
t666 = qJ(4) * t626 - t750;
t715 = qJD(3) * t626;
t634 = qJD(2) * t666 - qJD(4) * t630 - t715;
t443 = qJD(1) * t634 + qJDD(1) * t646 + t586;
t585 = pkin(7) * t693;
t595 = pkin(7) * t706;
t691 = qJDD(3) + t585 + t595;
t474 = pkin(3) * t649 - qJD(2) * qJD(4) - qJDD(2) * t753 + t691;
t414 = -t443 * t621 + t622 * t474;
t405 = pkin(4) * t649 - pkin(8) * t503 + t414;
t415 = t622 * t443 + t621 * t474;
t411 = -pkin(8) * t502 + t415;
t687 = t629 * t405 - t411 * t625;
t638 = -qJD(5) * t404 + t687;
t388 = pkin(5) * t556 - pkin(9) * t419 + t638;
t420 = -qJD(5) * t480 + t629 * t502 + t503 * t625;
t648 = t625 * t405 + t629 * t411 + t430 * t713 - t437 * t714;
t389 = -pkin(9) * t420 + t648;
t688 = t628 * t388 - t624 * t389;
t765 = -g(1) * t504 - g(2) * t506 + t431 * t663 + t594 * t614 + t688;
t550 = qJDD(6) + t556;
t764 = t550 * MDP(30) + (t663 ^ 2 - t773 ^ 2) * MDP(27) + t773 * MDP(26) * t663;
t611 = t630 * pkin(2);
t724 = t611 + t607;
t749 = qJ(4) * t630;
t672 = t724 + t749;
t551 = -pkin(1) - t672;
t573 = t759 * t626;
t560 = t622 * t573;
t464 = pkin(4) * t626 + t560 + (pkin(8) * t630 - t551) * t621;
t488 = t622 * t551 + t621 * t573;
t738 = t622 * t630;
t470 = -pkin(8) * t738 + t488;
t729 = t625 * t464 + t629 * t470;
t763 = t727 * t628;
t486 = -t557 * t624 - t628 * t660;
t752 = -pkin(8) - t753;
t566 = t752 * t621;
t567 = t752 * t622;
t725 = t629 * t566 + t625 * t567;
t602 = pkin(2) * t721;
t530 = qJD(1) * t666 + t602;
t472 = -t530 * t621 + t622 * t564;
t658 = -pkin(8) * t621 * t626 + pkin(4) * t630;
t446 = qJD(1) * t658 + t472;
t473 = t622 * t530 + t621 * t564;
t455 = pkin(8) * t698 + t473;
t762 = qJD(4) * t660 - qJD(5) * t725 - t629 * t446 + t455 * t625;
t674 = g(1) * t631 + g(2) * t627;
t700 = -pkin(4) * t622 - pkin(3);
t710 = -t700 * t721 + t767;
t716 = qJD(2) * t630;
t761 = t753 * t716 + t626 * (qJD(4) - t535);
t686 = t624 * t419 + t628 * t420;
t394 = -qJD(6) * t663 + t686;
t760 = qJD(4) * t557 + t625 * t446 + t629 * t455 + t566 * t714 - t567 * t713;
t758 = g(1) * t627;
t755 = g(2) * t631;
t754 = g(3) * t626;
t751 = pkin(7) * qJDD(2);
t748 = qJDD(2) * pkin(2);
t403 = t629 * t430 - t437 * t625;
t398 = pkin(9) * t480 + t403;
t396 = pkin(5) * t587 + t398;
t747 = t396 * t628;
t746 = t399 * t628;
t485 = t628 * t557 - t624 * t660;
t744 = t485 * t550;
t743 = t486 * t550;
t619 = t626 ^ 2;
t633 = qJD(1) ^ 2;
t741 = t619 * t633;
t735 = t626 * t633;
t734 = t627 * t630;
t733 = t630 * t631;
t592 = t621 * pkin(4) + qJ(3);
t732 = -qJD(6) * t485 - t624 * t726 + t763;
t731 = qJD(6) * t486 + t624 * t727 + t628 * t726;
t717 = qJD(2) * t626;
t601 = pkin(2) * t717;
t491 = t601 + t634;
t565 = t759 * t716;
t454 = t622 * t491 + t621 * t565;
t728 = pkin(5) * t726 + t710;
t574 = t759 * t630;
t620 = t630 ^ 2;
t723 = t619 - t620;
t711 = qJD(6) * t628;
t704 = t630 * t735;
t703 = t628 * t419 - t624 * t420 - t481 * t711;
t702 = g(1) * t736 + g(2) * t737 - t614;
t596 = pkin(7) * t705;
t616 = qJDD(2) * qJ(3);
t617 = qJD(2) * qJD(3);
t701 = t596 + t616 + t617;
t539 = pkin(4) * t738 + t574;
t699 = t759 * qJD(2);
t692 = t753 * t706;
t689 = -qJD(2) * pkin(2) + qJD(3);
t453 = -t491 * t621 + t622 * t565;
t438 = qJD(2) * t658 + t453;
t442 = pkin(8) * t622 * t717 + t454;
t685 = t629 * t438 - t442 * t625;
t683 = t629 * t464 - t470 * t625;
t681 = -t566 * t625 + t629 * t567;
t680 = qJD(6) * t396 + t389;
t678 = t631 * pkin(1) + pkin(2) * t733 + t627 * pkin(7) + qJ(3) * t736;
t677 = -t595 + t702;
t676 = qJD(1) * t699;
t632 = qJD(2) ^ 2;
t675 = pkin(7) * t632 + t755;
t673 = -t755 + t758;
t671 = -t660 * t556 + t587 * t727;
t459 = -pkin(9) * t557 + t725;
t669 = pkin(5) * t720 + pkin(9) * t727 + qJD(6) * t459 - t762;
t458 = pkin(9) * t660 + t681;
t668 = pkin(9) * t726 - qJD(6) * t458 + t760;
t667 = qJD(6) * t660 - t726;
t391 = t396 * t624 + t746;
t665 = t414 * t622 + t415 * t621;
t526 = t660 * t630;
t527 = t557 * t630;
t662 = t628 * t526 + t527 * t624;
t462 = t526 * t624 - t527 * t628;
t568 = t598 + t689;
t572 = -t599 - t618;
t661 = t568 * t630 + t572 * t626;
t659 = pkin(3) * t705 + qJDD(4) + t701;
t657 = t690 - t611;
t656 = t674 * t630;
t655 = -0.2e1 * pkin(1) * t707 - t751;
t546 = t657 * qJD(1);
t654 = t546 * t721 + qJDD(3) - t677;
t516 = (-pkin(7) + t700) * t717;
t653 = -qJ(3) * t716 - t715;
t651 = (-t450 * t621 + t451 * t622) * t626;
t650 = -t556 * t557 - t587 * t726;
t647 = t625 * t438 + t629 * t442 + t464 * t713 - t470 * t714;
t393 = t480 * t712 + t703;
t645 = 0.2e1 * qJDD(1) * pkin(1) - t675;
t479 = -t626 * t676 + t659;
t644 = -t479 * t630 + t674;
t569 = -pkin(1) - t724;
t643 = t751 + (-qJD(1) * t569 - t546) * qJD(2);
t641 = -t656 - t754;
t639 = t479 + t641;
t478 = qJD(1) * t653 + qJDD(1) * t657 + t586;
t532 = t601 + t653;
t637 = qJD(1) * t532 + qJDD(1) * t569 + t478 + t675;
t514 = pkin(7) * t694 - t701;
t529 = t691 - t748;
t635 = qJD(2) * t661 - t514 * t630 + t529 * t626;
t441 = pkin(4) * t502 + t479;
t612 = t631 * pkin(7);
t604 = cos(t615);
t603 = sin(t615);
t590 = g(1) * t734;
t584 = qJ(3) * t733;
t581 = qJ(3) * t734;
t563 = t626 * t699;
t561 = -qJ(3) * t720 + t602;
t521 = -t603 * t737 + t604 * t631;
t520 = t603 * t631 + t604 * t737;
t519 = t603 * t736 + t604 * t627;
t518 = -t603 * t627 + t604 * t736;
t510 = pkin(5) * t557 + t592;
t487 = -t551 * t621 + t560;
t477 = -pkin(5) * t526 + t539;
t467 = t540 * t630 - t660 * t717;
t466 = qJD(2) * t652 + qJD(5) * t526;
t436 = -pkin(5) * t467 + t516;
t413 = pkin(9) * t526 + t729;
t412 = pkin(5) * t626 + pkin(9) * t527 + t683;
t407 = qJD(6) * t462 + t466 * t624 - t628 * t467;
t406 = qJD(6) * t662 + t466 * t628 + t467 * t624;
t400 = pkin(5) * t420 + t441;
t395 = pkin(9) * t467 + t647;
t392 = pkin(5) * t716 - pkin(9) * t466 - qJD(5) * t729 + t685;
t390 = -t399 * t624 + t747;
t1 = [(t685 * t587 + t683 * t556 + t687 * t626 + t403 * t716 + t516 * t481 + t539 * t420 - t441 * t526 - t489 * t467 - g(1) * t521 - g(2) * t519 + (-t404 * t626 - t587 * t729) * qJD(5)) * MDP(24) + (qJDD(2) * t626 + t630 * t632) * MDP(6) + (qJDD(2) * t630 - t626 * t632) * MDP(7) + qJDD(1) * MDP(1) + (t626 * t643 + t630 * t637 - t590) * MDP(12) + (pkin(7) * t635 - g(1) * t612 - g(2) * t678 + t478 * t569 + t546 * t532 - t657 * t758) * MDP(14) + (t655 * t630 + (-t645 - t758) * t626) * MDP(10) + (t643 * t630 + (-t637 + t758) * t626) * MDP(13) + (-t453 * t554 - t454 * t552 - t487 * t503 - t488 * t502 + t590 + qJD(2) * t651 + (t414 * t621 - t415 * t622 - t755) * t630) * MDP(17) + 0.2e1 * (t626 * t705 - t707 * t723) * MDP(5) + (-t420 * t626 + t467 * t587 - t481 * t716 + t526 * t556) * MDP(22) + (t556 * t626 + t587 * t716) * MDP(23) + (t550 * t626 + t576 * t716) * MDP(30) + (qJDD(1) * t619 + 0.2e1 * t626 * t693) * MDP(4) + (t393 * t662 - t394 * t462 + t406 * t773 + t407 * t663) * MDP(27) + (-t394 * t626 - t407 * t576 + t550 * t662 + t716 * t773) * MDP(29) + ((t392 * t628 - t395 * t624) * t576 + (t412 * t628 - t413 * t624) * t550 + t688 * t626 + t390 * t716 - t436 * t773 + t477 * t394 - t400 * t662 + t431 * t407 - g(1) * t507 - g(2) * t505 + ((-t412 * t624 - t413 * t628) * t576 - t391 * t626) * qJD(6)) * MDP(31) + (g(1) * t520 - g(2) * t518 - t404 * t716 + t539 * t419 - t441 * t527 + t489 * t466 - t480 * t516 - t556 * t729 - t587 * t647 - t626 * t648) * MDP(25) + (t419 * t626 + t466 * t587 - t480 * t716 - t527 * t556) * MDP(21) + (t419 * t526 + t420 * t527 - t466 * t481 - t467 * t480) * MDP(20) + (-t419 * t527 - t466 * t480) * MDP(19) + (t415 * t488 + t451 * t454 + t414 * t487 + t450 * t453 + t479 * t574 - t535 * t563 - g(1) * (pkin(3) * t631 + t612) - g(2) * (qJ(4) * t733 + t678) + (-g(1) * (t657 - t749) - g(2) * pkin(3)) * t627) * MDP(18) + (t574 * t502 - t563 * t552 + (qJD(1) * t487 + t450) * t716 - t644 * t622 + (t453 * qJD(1) + t487 * qJDD(1) - t535 * t718 + t621 * t673 + t414) * t626) * MDP(15) + (t574 * t503 - t563 * t554 + (-qJD(1) * t488 - t451) * t716 + t644 * t621 + (-t454 * qJD(1) - t488 * qJDD(1) + t535 * t719 + t622 * t673 - t415) * t626) * MDP(16) + (t626 * t655 + t630 * t645 + t590) * MDP(9) + t673 * MDP(2) + ((t619 + t620) * qJDD(1) * pkin(7) + t635 - t674) * MDP(11) + t674 * MDP(3) + (t393 * t462 - t406 * t663) * MDP(26) + (-t391 * t716 + g(1) * t506 - g(2) * t504 + t477 * t393 + t397 * t626 + t400 * t462 + t431 * t406 - t436 * t663 + (-(-qJD(6) * t413 + t392) * t576 - t412 * t550 - t388 * t626) * t624 + (-(qJD(6) * t412 + t395) * t576 - t413 * t550 - t680 * t626) * t628) * MDP(32) + (t393 * t626 + t406 * t576 + t462 * t550 - t663 * t716) * MDP(28); (t592 * t420 + t441 * t557 + t710 * t481 + t726 * t489 + t681 * t556 + t587 * t762 + t641 * t603) * MDP(24) + t650 * MDP(22) + (-t576 * t731 - t744) * MDP(29) + t671 * MDP(21) + (t576 * t732 + t743) * MDP(28) + (t654 - 0.2e1 * t748) * MDP(12) + qJDD(2) * MDP(8) + (t479 * qJ(3) - t451 * t473 - t450 * t472 - g(1) * t584 - g(2) * t581 - g(3) * t672 + t709 * t535 + (-t450 * t622 - t451 * t621) * qJD(4) + (t626 * t674 - t665) * t753) * MDP(18) + (t754 - t596 + (pkin(1) * t633 + t674) * t630) * MDP(10) + ((-pkin(2) * t626 + t750) * qJDD(1) + ((-t572 - t618) * t626 + (-t568 + t689) * t630) * qJD(1)) * MDP(11) + MDP(6) * t706 + (pkin(1) * t735 + t677) * MDP(9) + (-t514 * qJ(3) - t572 * qJD(3) - t529 * pkin(2) - t546 * t561 - g(1) * (-pkin(2) * t736 + t584) - g(2) * (-pkin(2) * t737 + t581) - g(3) * t724 - t661 * qJD(1) * pkin(7)) * MDP(14) - MDP(4) * t704 + ((t458 * t628 - t459 * t624) * t550 + t510 * t394 + t400 * t485 + (t624 * t668 - t628 * t669) * t576 + t731 * t431 - t728 * t773 + t641 * t593) * MDP(31) + (-t561 * MDP(12) + MDP(21) * t480 + t481 * MDP(22) - t587 * MDP(23) - t403 * MDP(24) + t404 * MDP(25) + MDP(28) * t663 - MDP(29) * t773 - t576 * MDP(30) - t390 * MDP(31) + t391 * MDP(32)) * t720 + (-t393 * t485 - t394 * t486 + t663 * t731 + t732 * t773) * MDP(27) + t723 * MDP(5) * t633 + (t592 * t419 - t441 * t660 - t480 * t710 + t727 * t489 - t725 * t556 + t587 * t760 + t641 * t604) * MDP(25) + (-t419 * t660 - t480 * t727) * MDP(19) + (-t419 * t557 + t420 * t660 + t480 * t726 - t481 * t727) * MDP(20) + (t472 * t554 + t473 * t552 + (qJD(4) * t554 - t451 * t721 + t503 * t753 - t414) * t622 + (qJD(4) * t552 + t450 * t721 + t502 * t753 - t415) * t621 + t702) * MDP(17) + (t596 + 0.2e1 * t616 + 0.2e1 * t617 + (qJD(1) * t561 - g(3)) * t626 + (qJD(1) * t546 - t674) * t630) * MDP(13) + (t621 * t692 + qJ(3) * t503 + t709 * t554 + t639 * t622 + (t451 * t630 + t473 * t626 + t621 * t761) * qJD(1)) * MDP(16) + (-t622 * t692 + qJ(3) * t502 + t709 * t552 + t639 * t621 + (-t450 * t630 - t472 * t626 - t622 * t761) * qJD(1)) * MDP(15) + MDP(7) * t705 + (-(t458 * t624 + t459 * t628) * t550 + t510 * t393 + t400 * t486 + (t624 * t669 + t628 * t668) * t576 + t732 * t431 - t728 * t663 + t641 * t594) * MDP(32) + (t393 * t486 - t663 * t732) * MDP(26); MDP(11) * t706 + (qJDD(2) + t704) * MDP(12) + (-t632 - t741) * MDP(13) + (qJD(2) * t572 + t585 + t654 - t748) * MDP(14) + (t622 * t706 - t621 * t741 + (-t552 + t697) * qJD(2)) * MDP(15) + (-t621 * t706 - t622 * t741 + (-t554 - t695) * qJD(2)) * MDP(16) + (-t502 * t621 - t503 * t622 + (-t552 * t622 + t554 * t621) * t721) * MDP(17) + (qJD(1) * t651 - qJD(2) * t535 + t665 - t702) * MDP(18) + (-qJD(2) * t481 + t671) * MDP(24) + (qJD(2) * t480 + t650) * MDP(25) + (t743 + qJD(2) * t773 + (-t557 * t711 + t624 * t667 + t763) * t576) * MDP(31) + (-t744 + qJD(2) * t663 + (t667 * t628 + (qJD(6) * t557 - t727) * t624) * t576) * MDP(32); (t554 * t721 + t502) * MDP(15) + ((-t552 + t719) * t721 + t679) * MDP(16) + (-t552 ^ 2 - t554 ^ 2) * MDP(17) + (t450 * t554 + t451 * t552 - t656 + (-g(3) - t676) * t626 + t659) * MDP(18) + (t420 - t771) * MDP(24) + (t419 - t770) * MDP(25) + (t394 - t769) * MDP(31) + (t393 + t772) * MDP(32); -t480 * t481 * MDP(19) + (t480 ^ 2 - t481 ^ 2) * MDP(20) + (t419 + t770) * MDP(21) + (-t420 - t771) * MDP(22) + t556 * MDP(23) + (-g(1) * t518 - g(2) * t520 + t404 * t587 + t480 * t489 + t604 * t614 + t638) * MDP(24) + (g(1) * t519 - g(2) * t521 + t403 * t587 + t481 * t489 - t603 * t614 - t648) * MDP(25) + (t393 - t772) * MDP(28) + (-t394 - t769) * MDP(29) + (-(-t398 * t624 - t746) * t576 - t391 * qJD(6) + (-t480 * t773 + t550 * t628 - t576 * t712) * pkin(5) + t765) * MDP(31) + ((-t399 * t576 - t388) * t624 + (t398 * t576 - t680) * t628 + (-t480 * t663 - t550 * t624 - t576 * t711) * pkin(5) + t766) * MDP(32) + t764; (t703 - t772) * MDP(28) + (-t686 - t769) * MDP(29) + (t391 * t576 + t765) * MDP(31) + (-t624 * t388 - t628 * t389 + t390 * t576 + t766) * MDP(32) + (MDP(28) * t745 + MDP(29) * t663 - MDP(31) * t391 - MDP(32) * t747) * qJD(6) + t764;];
tau  = t1;
