% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% MDP [34x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:25:39
% EndTime: 2019-03-09 07:25:59
% DurationCPUTime: 13.24s
% Computational Cost: add. (6968->609), mult. (14138->809), div. (0->0), fcn. (10117->14), ass. (0->255)
t613 = sin(qJ(4));
t618 = cos(qJ(4));
t694 = t618 * qJD(3);
t619 = cos(qJ(3));
t709 = qJD(1) * t619;
t558 = t613 * t709 - t694;
t706 = qJD(3) * t613;
t560 = t618 * t709 + t706;
t612 = sin(qJ(5));
t617 = cos(qJ(5));
t490 = t617 * t558 + t560 * t612;
t616 = cos(qJ(6));
t611 = sin(qJ(6));
t649 = t558 * t612 - t617 * t560;
t747 = t649 * t611;
t436 = -t616 * t490 + t747;
t693 = qJD(1) * qJD(3);
t675 = t619 * t693;
t614 = sin(qJ(3));
t691 = qJDD(1) * t614;
t555 = qJDD(4) + t675 + t691;
t550 = qJDD(5) + t555;
t531 = qJDD(6) + t550;
t650 = t490 * t611 + t616 * t649;
t784 = t531 * MDP(32) + (-t436 ^ 2 + t650 ^ 2) * MDP(29) + t436 * MDP(28) * t650;
t699 = qJD(4) * t619;
t638 = -t613 * t699 - t614 * t694;
t690 = qJDD(1) * t619;
t479 = qJD(1) * t638 + qJD(4) * t694 + t613 * qJDD(3) + t618 * t690;
t705 = qJD(3) * t614;
t680 = t613 * t705;
t480 = -qJD(1) * t680 + qJD(4) * t560 - t618 * qJDD(3) + t613 * t690;
t697 = qJD(5) * t617;
t698 = qJD(5) * t612;
t418 = t617 * t479 - t612 * t480 - t558 * t697 - t560 * t698;
t627 = qJD(5) * t649 - t479 * t612 - t617 * t480;
t695 = qJD(6) * t616;
t685 = t616 * t418 - t490 * t695 + t611 * t627;
t696 = qJD(6) * t611;
t394 = t649 * t696 + t685;
t710 = qJD(1) * t614;
t592 = qJD(4) + t710;
t587 = qJD(5) + t592;
t670 = t418 * t611 - t616 * t627;
t628 = qJD(6) * t650 - t670;
t579 = qJD(6) + t587;
t774 = t579 * t650;
t775 = t436 * t579;
t783 = t550 * MDP(25) + (-t490 ^ 2 + t649 ^ 2) * MDP(22) + (t490 * t587 + t418) * MDP(23) + (-t587 * t649 + t627) * MDP(24) - t490 * MDP(21) * t649 + (t628 - t774) * MDP(31) + (t394 - t775) * MDP(30) + t784;
t571 = pkin(3) * t614 - pkin(8) * t619 + qJ(2);
t537 = t571 * qJD(1);
t621 = -pkin(1) - pkin(7);
t589 = qJD(1) * t621 + qJD(2);
t570 = t614 * t589;
t547 = qJD(3) * pkin(8) + t570;
t474 = t618 * t537 - t547 * t613;
t454 = -pkin(9) * t560 + t474;
t449 = pkin(4) * t592 + t454;
t475 = t537 * t613 + t547 * t618;
t455 = -pkin(9) * t558 + t475;
t453 = t617 * t455;
t413 = t449 * t612 + t453;
t777 = pkin(10) * t490;
t407 = t413 - t777;
t405 = t407 * t696;
t741 = t589 * t619;
t548 = -qJD(3) * pkin(3) - t741;
t502 = pkin(4) * t558 + t548;
t442 = pkin(5) * t490 + t502;
t610 = qJ(4) + qJ(5);
t606 = qJ(6) + t610;
t596 = sin(t606);
t597 = cos(t606);
t620 = cos(qJ(1));
t615 = sin(qJ(1));
t736 = t614 * t615;
t508 = t596 * t620 + t597 * t736;
t734 = t614 * t620;
t510 = -t596 * t615 + t597 * t734;
t752 = g(3) * t619;
t768 = g(1) * t508 - g(2) * t510 - t442 * t436 + t597 * t752 + t405;
t507 = -t596 * t736 + t597 * t620;
t509 = t596 * t734 + t597 * t615;
t658 = pkin(3) * t619 + pkin(8) * t614;
t556 = qJD(3) * t658 + qJD(2);
t488 = qJD(1) * t556 + qJDD(1) * t571;
t482 = t618 * t488;
t585 = qJDD(1) * t621 + qJDD(2);
t704 = qJD(3) * t619;
t504 = qJDD(3) * pkin(8) + t585 * t614 + t589 * t704;
t404 = pkin(4) * t555 - pkin(9) * t479 - qJD(4) * t475 - t504 * t613 + t482;
t700 = qJD(4) * t618;
t684 = -t613 * t488 - t618 * t504 - t537 * t700;
t702 = qJD(4) * t613;
t641 = -t547 * t702 - t684;
t409 = -pkin(9) * t480 + t641;
t671 = t617 * t404 - t612 * t409;
t629 = -qJD(5) * t413 + t671;
t390 = pkin(5) * t550 - pkin(10) * t418 + t629;
t660 = -t612 * t404 - t617 * t409 - t449 * t697 + t455 * t698;
t391 = pkin(10) * t627 - t660;
t672 = t616 * t390 - t611 * t391;
t767 = -g(1) * t507 - g(2) * t509 + t442 * t650 + t596 * t752 + t672;
t567 = t658 * qJD(1);
t546 = t618 * t567;
t756 = pkin(8) + pkin(9);
t682 = qJD(4) * t756;
t735 = t614 * t618;
t688 = pkin(9) * t735;
t738 = t613 * t619;
t780 = -t589 * t738 + t546 + (pkin(4) * t619 + t688) * qJD(1) + t618 * t682;
t681 = t613 * t710;
t730 = t618 * t619;
t717 = t613 * t567 + t589 * t730;
t779 = pkin(9) * t681 + t613 * t682 + t717;
t656 = g(1) * t615 - g(2) * t620;
t503 = -qJDD(3) * pkin(3) - t585 * t619 + t589 * t705;
t753 = g(3) * t614;
t634 = t619 * t656 - t753;
t778 = qJD(4) * pkin(8) * t592 + t503 + t634;
t776 = pkin(10) * t649;
t563 = t612 * t618 + t613 * t617;
t757 = qJD(4) + qJD(5);
t772 = t757 * t563;
t562 = t612 * t613 - t617 * t618;
t762 = t614 * t562;
t720 = -qJD(1) * t762 - t757 * t562;
t536 = t563 * qJD(1);
t719 = t614 * t536 + t772;
t770 = qJDD(2) - t656;
t769 = t618 * t699 - t680;
t604 = sin(t610);
t605 = cos(t610);
t518 = t604 * t620 + t605 * t736;
t520 = -t604 * t615 + t605 * t734;
t766 = g(1) * t518 - g(2) * t520 + t490 * t502 + t605 * t752 + t660;
t517 = -t604 * t736 + t605 * t620;
t519 = t604 * t734 + t605 * t615;
t765 = -g(1) * t517 - g(2) * t519 + t502 * t649 + t604 * t752 + t629;
t554 = t618 * t571;
t737 = t613 * t621;
t674 = pkin(4) - t737;
t486 = -pkin(9) * t730 + t614 * t674 + t554;
t586 = t621 * t735;
t715 = t613 * t571 + t586;
t501 = -pkin(9) * t738 + t715;
t721 = t612 * t486 + t617 * t501;
t659 = -t570 + (t681 + t702) * pkin(4);
t760 = t780 * t617;
t580 = t756 * t613;
t581 = t756 * t618;
t716 = -t612 * t580 + t617 * t581;
t759 = t580 * t697 + t581 * t698 + t612 * t780 + t779 * t617;
t751 = pkin(1) * qJDD(1);
t758 = t751 - t770;
t755 = pkin(4) * t612;
t623 = qJD(1) ^ 2;
t750 = qJ(2) * t623;
t451 = t612 * t455;
t412 = t617 * t449 - t451;
t406 = t412 + t776;
t403 = pkin(5) * t587 + t406;
t749 = t403 * t616;
t748 = t479 * t613;
t746 = t555 * t613;
t745 = t555 * t618;
t744 = t558 * t592;
t743 = t560 * t592;
t742 = t560 * t618;
t740 = t611 * t531;
t739 = t612 * t616;
t733 = t615 * t618;
t732 = t616 * t407;
t731 = t616 * t531;
t729 = t618 * t620;
t496 = t616 * t562 + t563 * t611;
t728 = -qJD(6) * t496 - t611 * t719 + t616 * t720;
t497 = -t562 * t611 + t563 * t616;
t727 = qJD(6) * t497 + t611 * t720 + t616 * t719;
t526 = t562 * t619;
t726 = qJD(3) * t526 + t614 * t772 + t536;
t725 = t562 * qJD(1) - t563 * t704 + t757 * t762;
t724 = t617 * t454 - t451;
t722 = pkin(5) * t719 + t659;
t609 = t619 ^ 2;
t714 = t614 ^ 2 - t609;
t622 = qJD(3) ^ 2;
t713 = -t622 - t623;
t708 = qJD(3) * t558;
t707 = qJD(3) * t560;
t703 = qJD(3) * t621;
t701 = qJD(4) * t614;
t692 = qJDD(1) * qJ(2);
t689 = qJDD(3) * t614;
t678 = t619 * t703;
t683 = t613 * t556 + t571 * t700 + t618 * t678;
t599 = -pkin(4) * t618 - pkin(3);
t530 = t618 * t556;
t430 = t530 + (-t586 + (pkin(9) * t619 - t571) * t613) * qJD(4) + (t619 * t674 + t688) * qJD(3);
t439 = -pkin(9) * t769 - t701 * t737 + t683;
t669 = t617 * t430 - t439 * t612;
t668 = -t454 * t612 - t453;
t666 = t617 * t486 - t501 * t612;
t665 = t592 * t621 + t547;
t664 = -t617 * t580 - t581 * t612;
t557 = pkin(4) * t738 - t619 * t621;
t663 = -qJD(4) * t537 - t504;
t662 = qJD(6) * t403 + t391;
t661 = qJD(1) + t701;
t657 = g(1) * t620 + g(2) * t615;
t470 = -pkin(10) * t562 + t716;
t655 = pkin(5) * t709 + pkin(10) * t720 + qJD(5) * t716 + qJD(6) * t470 - t612 * t779 + t760;
t469 = -pkin(10) * t563 + t664;
t654 = pkin(10) * t719 - qJD(6) * t469 + t759;
t523 = t563 * t614;
t653 = qJD(6) * t523 + t726;
t652 = -qJD(6) * t762 - t725;
t393 = t611 * t403 + t732;
t524 = t563 * t619;
t462 = t616 * t524 - t526 * t611;
t463 = -t524 * t611 - t526 * t616;
t646 = t592 * t700 + t746;
t645 = -t592 * t702 + t745;
t505 = pkin(4) * t769 + t614 * t703;
t642 = 0.2e1 * qJ(2) * t693 + qJDD(3) * t621;
t640 = t612 * t430 + t617 * t439 + t486 * t697 - t501 * t698;
t637 = 0.2e1 * qJD(1) * qJD(2) - t657;
t636 = -t585 + t656 + t750;
t635 = -pkin(8) * t555 + t548 * t592;
t440 = pkin(4) * t480 + t503;
t630 = t637 + 0.2e1 * t692;
t626 = -t621 * t622 + t630;
t602 = qJDD(3) * t619;
t598 = pkin(4) * t617 + pkin(5);
t541 = -t613 * t615 + t614 * t729;
t540 = t613 * t734 + t733;
t539 = t613 * t620 + t614 * t733;
t538 = -t613 * t736 + t729;
t515 = pkin(5) * t562 + t599;
t495 = pkin(5) * t524 + t557;
t456 = pkin(4) * t560 - pkin(5) * t649;
t448 = -t698 * t738 + (t730 * t757 - t680) * t617 + t638 * t612;
t446 = qJD(3) * t762 - t619 * t772;
t427 = pkin(5) * t448 + t505;
t426 = -pkin(10) * t524 + t721;
t423 = pkin(5) * t614 + pkin(10) * t526 + t666;
t411 = t724 + t776;
t410 = t668 + t777;
t400 = qJD(6) * t463 + t446 * t611 + t616 * t448;
t399 = -qJD(6) * t462 + t446 * t616 - t448 * t611;
t398 = -pkin(5) * t627 + t440;
t397 = -pkin(10) * t448 + t640;
t396 = pkin(5) * t704 - pkin(10) * t446 - qJD(5) * t721 + t669;
t392 = -t407 * t611 + t749;
t1 = [(qJDD(1) * t609 - 0.2e1 * t614 * t675) * MDP(7) + (-0.2e1 * t751 + t770) * MDP(4) + (t614 * t626 + t619 * t642) * MDP(12) + (-t614 * t642 + t619 * t626) * MDP(13) + (-t393 * t704 + g(1) * t509 - g(2) * t507 + t495 * t394 + t398 * t463 + t442 * t399 + t405 * t614 - t427 * t650 + (-(-qJD(6) * t426 + t396) * t579 - t423 * t531 - t390 * t614) * t611 + (-(qJD(6) * t423 + t397) * t579 - t426 * t531 - t662 * t614) * t616) * MDP(34) + (t394 * t614 + t399 * t579 + t463 * t531 - t650 * t704) * MDP(30) + (t394 * t463 - t399 * t650) * MDP(28) + (g(1) * t519 - g(2) * t517 - t413 * t704 + t557 * t418 - t440 * t526 + t502 * t446 - t505 * t649 - t550 * t721 - t587 * t640 + t614 * t660) * MDP(27) + (t418 * t614 + t446 * t587 - t526 * t550 - t649 * t704) * MDP(23) + (-t418 * t526 - t446 * t649) * MDP(21) + (t758 * pkin(1) + (t637 + t692) * qJ(2)) * MDP(6) + t656 * MDP(2) + t630 * MDP(5) + t657 * MDP(3) + 0.2e1 * (-t614 * t690 + t693 * t714) * MDP(8) + (-t683 * t592 - t715 * t555 + g(1) * t540 - g(2) * t538 + (t665 * t702 + (-t548 * t618 + t560 * t621) * qJD(3) + t684) * t614 + (-qJD(3) * t475 - t479 * t621 + t503 * t618 - t548 * t702) * t619) * MDP(20) + ((-t592 * t694 + t479) * t614 + (t645 + t707) * t619) * MDP(16) + ((t592 * t706 - t480) * t614 + (-t646 - t708) * t619) * MDP(17) + (t479 * t730 + t560 * t638) * MDP(14) + (-g(1) * t541 - g(2) * t539 + t530 * t592 + t554 * t555 + (t558 * t703 - t665 * t700 + t482) * t614 + (qJD(3) * t474 - t480 * t621 + t548 * t700) * t619 + ((-qJD(4) * t571 - t678) * t592 + t503 * t619 + (-qJD(3) * t548 - t555 * t621 + t663) * t614) * t613) * MDP(19) + (t555 * t614 + t592 * t704) * MDP(18) + (t550 * t614 + t587 * t704) * MDP(25) + (t531 * t614 + t579 * t704) * MDP(32) + (-t619 * t622 - t689) * MDP(10) + (-t394 * t462 + t399 * t436 + t400 * t650 + t463 * t628) * MDP(29) + (-t400 * t579 + t436 * t704 - t462 * t531 + t614 * t628) * MDP(31) + ((t396 * t616 - t397 * t611) * t579 + (t423 * t616 - t426 * t611) * t531 + t672 * t614 + t392 * t704 - t427 * t436 - t495 * t628 + t398 * t462 + t442 * t400 - g(1) * t510 - g(2) * t508 + ((-t423 * t611 - t426 * t616) * t579 - t393 * t614) * qJD(6)) * MDP(33) + ((t558 * t618 + t560 * t613) * t705 + (-t748 - t480 * t618 + (t558 * t613 - t742) * qJD(4)) * t619) * MDP(15) + qJDD(1) * MDP(1) + (-t418 * t524 - t446 * t490 + t448 * t649 - t526 * t627) * MDP(22) + (-t448 * t587 - t490 * t704 - t524 * t550 + t614 * t627) * MDP(24) + (t669 * t587 + t666 * t550 + t671 * t614 + t412 * t704 + t505 * t490 - t557 * t627 + t440 * t524 + t502 * t448 - g(1) * t520 - g(2) * t518 + (-t413 * t614 - t587 * t721) * qJD(5)) * MDP(26) + (-t614 * t622 + t602) * MDP(9); qJDD(1) * MDP(4) - t623 * MDP(5) + (-t750 - t758) * MDP(6) + (t614 * t713 + t602) * MDP(12) + (t619 * t713 - t689) * MDP(13) + (-t480 * t619 + (t708 - t746) * t614 + (-t613 * t704 - t618 * t661) * t592) * MDP(19) + (-t479 * t619 + (t707 - t745) * t614 + (t613 * t661 - t619 * t694) * t592) * MDP(20) + (t490 * t705 - t523 * t550 + t587 * t725 + t619 * t627) * MDP(26) + (-t418 * t619 + t550 * t762 + t587 * t726 - t649 * t705) * MDP(27) + ((-t523 * t616 + t611 * t762) * t531 - t436 * t705 + t619 * t628 + (t611 * t653 - t616 * t652) * t579) * MDP(33) + (-(-t523 * t611 - t616 * t762) * t531 - t650 * t705 - t619 * t394 + (t611 * t652 + t616 * t653) * t579) * MDP(34); ((t479 - t744) * t618 + (-t480 - t743) * t613) * MDP(15) + ((-t592 * t613 * t614 + t558 * t619) * qJD(1) + t645) * MDP(17) + (t394 * t497 - t650 * t728) * MDP(28) + (-(t469 * t611 + t470 * t616) * t531 + t515 * t394 + t398 * t497 + (t611 * t655 + t616 * t654) * t579 + t728 * t442 - t722 * t650 + t634 * t596) * MDP(34) + (t418 * t563 - t649 * t720) * MDP(21) + (t599 * t418 + t440 * t563 + t720 * t502 - t716 * t550 + t587 * t759 + t634 * t604 - t649 * t659) * MDP(27) + (MDP(7) * t614 * t619 - MDP(8) * t714) * t623 + ((-t560 * t619 + t592 * t735) * qJD(1) + t646) * MDP(16) - MDP(10) * t691 + (t614 * t636 + t752) * MDP(13) + (-t619 * t636 + t753) * MDP(12) + (t550 * t563 + t587 * t720) * MDP(23) + (-pkin(3) * t479 - t560 * t570 + t717 * t592 + t613 * t778 + t635 * t618) * MDP(20) + (-t558 * t570 - pkin(3) * t480 - t546 * t592 + (t592 * t741 + t635) * t613 - t778 * t618) * MDP(19) + (-t592 * MDP(18) - MDP(19) * t474 + t475 * MDP(20) + MDP(23) * t649 + t490 * MDP(24) - t587 * MDP(25) - t412 * MDP(26) + t413 * MDP(27) + MDP(30) * t650 - MDP(31) * t436 - t579 * MDP(32) - t392 * MDP(33) + t393 * MDP(34)) * t709 + (-t394 * t496 + t436 * t728 + t497 * t628 + t650 * t727) * MDP(29) + ((t469 * t616 - t470 * t611) * t531 - t515 * t628 + t398 * t496 + (t611 * t654 - t616 * t655) * t579 + t727 * t442 - t722 * t436 - t634 * t597) * MDP(33) + (t664 * t550 - t599 * t627 + t440 * t562 + (-t581 * t697 + (qJD(5) * t580 + t779) * t612 - t760) * t587 + t719 * t502 + t659 * t490 - t634 * t605) * MDP(26) + MDP(9) * t690 + (t592 * t742 + t748) * MDP(14) + (-t550 * t562 - t587 * t719) * MDP(24) + qJDD(3) * MDP(11) + (t497 * t531 + t579 * t728) * MDP(30) + (-t496 * t531 - t579 * t727) * MDP(31) + (-t418 * t562 - t490 * t720 + t563 * t627 + t649 * t719) * MDP(22); (t479 + t744) * MDP(16) + (t456 * t650 + (-t598 * t531 - t390 + (t410 - (-qJD(5) - qJD(6)) * t755) * t579) * t611 + (-t531 * t755 + (-pkin(4) * t697 - qJD(6) * t598 + t411) * t579 - t662) * t616 + t768) * MDP(34) + (-t668 * t587 + (-t490 * t560 + t550 * t617 - t587 * t698) * pkin(4) + t765) * MDP(26) + (t724 * t587 + (-t550 * t612 + t560 * t649 - t587 * t697) * pkin(4) + t766) * MDP(27) + (t598 * t731 - (t410 * t616 - t411 * t611) * t579 + t456 * t436 + (-t612 * t740 + (-t611 * t617 - t739) * t579 * qJD(5)) * pkin(4) + ((-pkin(4) * t739 - t598 * t611) * t579 - t393) * qJD(6) + t767) * MDP(33) + (-t480 + t743) * MDP(17) + (g(1) * t539 - g(2) * t541 + g(3) * t730 + t474 * t592 + t548 * t558 - t641) * MDP(20) + (-t547 * t700 - g(1) * t538 - g(2) * t540 + t475 * t592 - t548 * t560 + t482 + (t663 + t752) * t613) * MDP(19) + t560 * t558 * MDP(14) + t555 * MDP(18) + (-t558 ^ 2 + t560 ^ 2) * MDP(15) + t783; (t413 * t587 + t765) * MDP(26) + (t412 * t587 + t766) * MDP(27) + (-(-t406 * t611 - t732) * t579 - t393 * qJD(6) + (-t436 * t649 - t579 * t696 + t731) * pkin(5) + t767) * MDP(33) + ((-t407 * t579 - t390) * t611 + (t406 * t579 - t662) * t616 + (-t579 * t695 - t649 * t650 - t740) * pkin(5) + t768) * MDP(34) + t783; (t685 - t775) * MDP(30) + (-t670 - t774) * MDP(31) + (t393 * t579 + t767) * MDP(33) + (-t611 * t390 - t616 * t391 + t392 * t579 + t768) * MDP(34) + (MDP(30) * t747 + MDP(31) * t650 - MDP(33) * t393 - MDP(34) * t749) * qJD(6) + t784;];
tau  = t1;
