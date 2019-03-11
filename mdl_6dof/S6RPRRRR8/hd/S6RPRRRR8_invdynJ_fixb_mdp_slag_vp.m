% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR8
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
%   see S6RPRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(34,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [34 1]), ...
  'S6RPRRRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [34x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:21:29
% EndTime: 2019-03-09 07:21:40
% DurationCPUTime: 7.71s
% Computational Cost: add. (6680->559), mult. (13076->718), div. (0->0), fcn. (9548->14), ass. (0->258)
t618 = sin(qJ(4));
t623 = cos(qJ(3));
t753 = cos(qJ(4));
t688 = t753 * t623;
t619 = sin(qJ(3));
t711 = qJD(1) * t619;
t536 = -qJD(1) * t688 + t618 * t711;
t608 = qJD(3) + qJD(4);
t617 = sin(qJ(5));
t622 = cos(qJ(5));
t506 = -t536 * t617 - t622 * t608;
t647 = -t618 * t623 - t619 * t753;
t537 = t647 * qJD(1);
t769 = qJD(5) - t537;
t781 = t506 * t769;
t653 = t536 * t622 - t608 * t617;
t780 = t653 * t769;
t624 = cos(qJ(1));
t606 = g(2) * t624;
t620 = sin(qJ(1));
t752 = g(1) * t620;
t764 = -t606 + t752;
t615 = qJ(3) + qJ(4);
t603 = cos(t615);
t568 = t603 * t606;
t601 = sin(t615);
t652 = -g(3) * t601 - t568;
t694 = t603 * t752;
t773 = t652 + t694;
t625 = -pkin(1) - pkin(7);
t572 = qJD(1) * t625 + qJD(2);
t527 = -pkin(8) * t711 + t572 * t619;
t519 = t618 * t527;
t710 = qJD(1) * t623;
t528 = -pkin(8) * t710 + t623 * t572;
t481 = t528 * t753 - t519;
t683 = t753 * qJD(4);
t772 = -pkin(3) * t683 + t481;
t493 = -pkin(4) * t536 - pkin(9) * t537;
t484 = pkin(3) * t710 + t493;
t779 = -t622 * t484 + t617 * t772;
t621 = cos(qJ(6));
t616 = sin(qJ(6));
t735 = t653 * t616;
t452 = t621 * t506 - t735;
t526 = qJD(6) + t769;
t778 = t452 * t526;
t654 = t506 * t616 + t621 * t653;
t777 = t526 * t654;
t607 = qJDD(3) + qJDD(4);
t570 = qJDD(1) * t625 + qJDD(2);
t552 = t623 * t570;
t700 = qJD(1) * qJD(3);
t682 = t619 * t700;
t697 = qJDD(1) * t623;
t709 = qJD(3) * t619;
t483 = -t572 * t709 + qJDD(3) * pkin(3) + t552 + (t682 - t697) * pkin(8);
t708 = qJD(3) * t623;
t681 = t623 * t700;
t698 = qJDD(1) * t619;
t771 = t681 + t698;
t486 = -pkin(8) * t771 + t570 * t619 + t572 * t708;
t521 = qJD(3) * pkin(3) + t528;
t707 = qJD(4) * t618;
t669 = -t753 * t483 + t618 * t486 + t521 * t707 + t527 * t683;
t417 = -pkin(4) * t607 + t669;
t776 = t417 + t694;
t727 = t616 * t622;
t548 = t617 * t621 + t727;
t758 = qJD(5) + qJD(6);
t499 = t758 * t548;
t719 = t548 * t537 - t499;
t545 = t616 * t617 - t621 * t622;
t775 = (t537 - t758) * t545;
t774 = qJDD(2) - t764;
t609 = qJDD(1) * qJ(2);
t610 = qJD(1) * qJD(2);
t664 = g(1) * t624 + g(2) * t620;
t637 = -t664 + 0.2e1 * t610;
t770 = 0.2e1 * t609 + t637;
t520 = t753 * t527;
t472 = t618 * t521 + t520;
t465 = pkin(9) * t608 + t472;
t562 = pkin(3) * t711 + qJD(1) * qJ(2);
t475 = -pkin(4) * t537 + pkin(9) * t536 + t562;
t432 = t465 * t622 + t475 * t617;
t420 = -pkin(10) * t506 + t432;
t704 = qJD(6) * t616;
t418 = t420 * t704;
t471 = t521 * t753 - t519;
t464 = -t608 * pkin(4) - t471;
t440 = t506 * pkin(5) + t464;
t614 = qJ(5) + qJ(6);
t600 = sin(t614);
t602 = cos(t614);
t729 = t602 * t620;
t516 = t600 * t624 + t601 * t729;
t728 = t602 * t624;
t518 = -t600 * t620 + t601 * t728;
t588 = g(3) * t603;
t768 = g(1) * t516 - g(2) * t518 + t440 * t452 + t602 * t588 + t418;
t731 = t601 * t620;
t515 = -t600 * t731 + t728;
t730 = t601 * t624;
t517 = t600 * t730 + t729;
t633 = t618 * t483 + t486 * t753 + t521 * t683 - t527 * t707;
t416 = t607 * pkin(9) + t633;
t679 = qJDD(1) * t753;
t765 = t537 * t608;
t462 = -t618 * t698 + t623 * t679 + t765;
t636 = t753 * qJD(3) + t683;
t687 = t619 * t707;
t463 = -qJD(1) * t687 - t618 * t682 + t619 * t679 + t623 * (qJD(1) * t636 + qJDD(1) * t618);
t522 = pkin(3) * t771 + t609 + t610;
t424 = pkin(4) * t463 - pkin(9) * t462 + t522;
t422 = t622 * t424;
t705 = qJD(5) * t622;
t706 = qJD(5) * t617;
t437 = t622 * t462 + t536 * t706 + t617 * t607 + t608 * t705;
t461 = qJDD(5) + t463;
t393 = pkin(5) * t461 - pkin(10) * t437 - qJD(5) * t432 - t416 * t617 + t422;
t438 = -qJD(5) * t653 + t462 * t617 - t622 * t607;
t641 = t622 * t416 + t617 * t424 - t465 * t706 + t475 * t705;
t394 = -pkin(10) * t438 + t641;
t677 = t621 * t393 - t616 * t394;
t767 = -g(1) * t515 - g(2) * t517 + t440 * t654 + t600 * t588 + t677;
t460 = qJDD(6) + t461;
t766 = t460 * MDP(32) + (-t452 ^ 2 + t654 ^ 2) * MDP(29) - t452 * MDP(28) * t654;
t480 = t618 * t528 + t520;
t666 = pkin(3) * t707 - t480;
t747 = pkin(1) * qJDD(1);
t763 = t747 - t774;
t762 = t617 * t484 + t622 * t772;
t761 = t616 * t706 + t617 * t704;
t734 = t537 * t617;
t760 = (t706 - t734) * pkin(5);
t759 = -qJD(6) * t622 - t705;
t676 = t437 * t616 + t621 * t438;
t404 = -qJD(6) * t654 + t676;
t754 = -pkin(9) - pkin(10);
t751 = g(3) * t622;
t750 = t622 * pkin(5);
t605 = t622 * pkin(10);
t749 = pkin(8) - t625;
t592 = pkin(3) * t618 + pkin(9);
t748 = -pkin(10) - t592;
t627 = qJD(1) ^ 2;
t746 = qJ(2) * t627;
t431 = -t465 * t617 + t622 * t475;
t419 = pkin(10) * t653 + t431;
t413 = pkin(5) * t769 + t419;
t745 = t413 * t621;
t744 = t420 * t621;
t743 = t437 * t617;
t742 = t460 * t545;
t741 = t460 * t548;
t740 = t461 * t622;
t739 = t464 * t537;
t500 = -t618 * t708 - t619 * t636 - t623 * t707;
t738 = t500 * t617;
t737 = t500 * t622;
t501 = -t618 * t709 + t623 * t636 - t687;
t736 = t501 * t526;
t546 = t618 * t619 - t688;
t733 = t546 * t617;
t732 = t546 * t622;
t726 = t617 * t461;
t725 = t617 * t620;
t724 = t617 * t624;
t723 = t620 * t622;
t556 = t749 * t619;
t557 = t749 * t623;
t505 = -t556 * t753 - t618 * t557;
t496 = t622 * t505;
t722 = t622 * t624;
t586 = t619 * pkin(3) + qJ(2);
t721 = t622 * t471 + t617 * t493;
t717 = t500 * t608 - t546 * t607;
t494 = -pkin(4) * t647 + pkin(9) * t546 + t586;
t716 = t617 * t494 + t496;
t715 = t760 + t666;
t613 = t623 ^ 2;
t714 = t619 ^ 2 - t613;
t626 = qJD(3) ^ 2;
t713 = -t626 - t627;
t703 = qJD(6) * t621;
t573 = pkin(3) * t708 + qJD(2);
t696 = qJDD(3) * t619;
t695 = pkin(10) * t734;
t691 = pkin(9) * qJD(5) * t769;
t690 = t621 * t437 - t616 * t438 - t506 * t703;
t689 = qJD(5) * t754;
t458 = t464 * t705;
t680 = qJD(5) * t748;
t674 = t622 * t769;
t673 = -qJD(5) * t475 - t416;
t672 = qJD(6) * t413 + t394;
t671 = -qJD(5) * t647 + qJD(1);
t593 = -pkin(3) * t753 - pkin(4);
t668 = -t536 * pkin(5) - t537 * t605;
t665 = -t472 + t760;
t663 = -t465 * t705 + t422;
t540 = t592 * t622 + t605;
t662 = qJD(6) * t540 - t622 * t680 + t668 - t779;
t488 = t622 * t493;
t564 = pkin(9) * t622 + t605;
t661 = qJD(6) * t564 - t471 * t617 - t622 * t689 + t488 + t668;
t539 = t748 * t617;
t660 = -qJD(6) * t539 - t617 * t680 - t695 + t762;
t563 = t754 * t617;
t659 = -qJD(6) * t563 - t617 * t689 - t695 + t721;
t658 = -pkin(9) * t461 - t739;
t399 = t413 * t616 + t744;
t656 = -t461 * t592 - t739;
t655 = -t501 * t608 + t607 * t647;
t651 = -t432 * t536 + t776 * t617 + t458;
t650 = t431 * t536 + t464 * t706 + t622 * t568 + t601 * t751;
t648 = t618 * t556 - t557 * t753;
t646 = t546 * t705 - t738;
t645 = t546 * t706 + t737;
t643 = t526 * t545;
t642 = 0.2e1 * qJ(2) * t700 + qJDD(3) * t625;
t445 = pkin(4) * t501 - pkin(9) * t500 + t573;
t541 = t749 * t709;
t542 = qJD(3) * t557;
t448 = qJD(4) * t648 + t618 * t541 - t542 * t753;
t640 = t617 * t445 + t622 * t448 + t494 * t705 - t505 * t706;
t403 = t653 * t704 + t690;
t638 = -t746 - t764;
t634 = t562 * t536 - t669 - t773;
t632 = -t625 * t626 + t770;
t398 = -t420 * t616 + t745;
t407 = pkin(5) * t438 + t417;
t631 = t398 * t536 + t407 * t545 - t719 * t440 - t602 * t773;
t449 = qJD(4) * t505 - t541 * t753 - t618 * t542;
t630 = -t399 * t536 + t407 * t548 + t775 * t440 + t773 * t600;
t629 = g(1) * t731 - g(2) * t730 - t562 * t537 + t588 - t633;
t628 = (-t403 * t545 - t404 * t548 - t452 * t775 - t654 * t719) * MDP(29) + (t403 * t548 - t654 * t775) * MDP(28) + ((t437 - t781) * t622 + (-t438 + t780) * t617) * MDP(22) + (t526 * t775 - t536 * t654 + t741) * MDP(30) + (-t452 * t536 + t526 * t719 - t742) * MDP(31) + (-t653 * t674 + t743) * MDP(21) + (-t617 * t769 ^ 2 - t506 * t536 + t740) * MDP(24) + (-t536 * t653 + t674 * t769 + t726) * MDP(23) + (t462 - t765) * MDP(16) + (-t536 * t608 - t463) * MDP(17) + (t536 ^ 2 - t537 ^ 2) * MDP(15) + t607 * MDP(18) + (MDP(14) * t537 + MDP(25) * t769 + MDP(32) * t526) * t536;
t599 = qJDD(3) * t623;
t594 = -pkin(4) - t750;
t561 = t593 - t750;
t532 = t601 * t722 - t725;
t531 = t601 * t724 + t723;
t530 = t601 * t723 + t724;
t529 = -t601 * t725 + t722;
t492 = t622 * t494;
t490 = t545 * t546;
t489 = t548 * t546;
t473 = -pkin(5) * t733 - t648;
t442 = t622 * t445;
t439 = pkin(10) * t733 + t716;
t433 = -pkin(5) * t647 + pkin(10) * t732 - t505 * t617 + t492;
t426 = -pkin(5) * t646 + t449;
t412 = t500 * t727 + (-t732 * t758 + t738) * t621 + t761 * t546;
t411 = t499 * t546 - t500 * t545;
t400 = pkin(10) * t646 + t640;
t397 = -pkin(10) * t737 + pkin(5) * t501 - t448 * t617 + t442 + (-t496 + (-pkin(10) * t546 - t494) * t617) * qJD(5);
t1 = [(t438 * t647 - t501 * t506 + t546 * t726 + t646 * t769) * MDP(24) + (-t461 * t647 + t501 * t769) * MDP(25) + ((-t505 * t705 + t442) * t769 + t492 * t461 - t663 * t647 + t431 * t501 + t449 * t506 - t648 * t438 - t546 * t458 - g(1) * t532 - g(2) * t530 + ((-qJD(5) * t494 - t448) * t769 - t505 * t461 - t673 * t647 - t417 * t546 + t464 * t500) * t617) * MDP(26) + (-t437 * t647 - t461 * t732 - t501 * t653 + t645 * t769) * MDP(23) + (g(1) * t531 - g(2) * t529 - t417 * t732 - t432 * t501 - t437 * t648 - t449 * t653 - t461 * t716 + t464 * t645 - t640 * t769 + t641 * t647) * MDP(27) + (t763 * pkin(1) + (t637 + t609) * qJ(2)) * MDP(6) + t764 * MDP(2) + t770 * MDP(5) + (-0.2e1 * t747 + t774) * MDP(4) + (t403 * t490 - t411 * t654) * MDP(28) + (t403 * t489 - t404 * t490 - t411 * t452 + t412 * t654) * MDP(29) + (-t437 * t732 - t645 * t653) * MDP(21) + ((-t506 * t622 + t617 * t653) * t500 + (t743 + t438 * t622 + (-t506 * t617 - t622 * t653) * qJD(5)) * t546) * MDP(22) + ((t397 * t621 - t400 * t616) * t526 + (t433 * t621 - t439 * t616) * t460 - t677 * t647 + t398 * t501 + t426 * t452 + t473 * t404 - t407 * t489 + t440 * t412 - g(1) * t518 - g(2) * t516 + ((-t433 * t616 - t439 * t621) * t526 + t399 * t647) * qJD(6)) * MDP(33) + (t404 * t647 - t412 * t526 - t452 * t501 + t460 * t489) * MDP(31) + (t462 * t647 + t463 * t546 + t500 * t537 + t501 * t536) * MDP(15) + (-t460 * t647 + t736) * MDP(32) + (g(1) * t517 - g(2) * t515 - t399 * t501 + t473 * t403 + t407 * t490 + t440 * t411 - t418 * t647 - t426 * t654 + (-(-qJD(6) * t439 + t397) * t526 - t433 * t460 + t393 * t647) * t616 + (-(qJD(6) * t433 + t400) * t526 - t439 * t460 + t672 * t647) * t621) * MDP(34) + (-t403 * t647 + t411 * t526 + t460 * t490 - t501 * t654) * MDP(30) + (-t449 * t608 + t463 * t586 + t501 * t562 - t522 * t647 - t537 * t573 - t601 * t664 + t607 * t648) * MDP(19) + t655 * MDP(17) + (t619 * t632 + t623 * t642) * MDP(12) + (-t619 * t642 + t623 * t632) * MDP(13) + t664 * MDP(3) + (-t448 * t608 + t462 * t586 + t500 * t562 - t505 * t607 - t522 * t546 - t536 * t573 - t603 * t664) * MDP(20) + qJDD(1) * MDP(1) + (qJDD(1) * t613 - 0.2e1 * t619 * t681) * MDP(7) + (-t623 * t626 - t696) * MDP(10) + 0.2e1 * (-t619 * t697 + t700 * t714) * MDP(8) + t717 * MDP(16) + (-t462 * t546 - t500 * t536) * MDP(14) + (-t619 * t626 + t599) * MDP(9); qJDD(1) * MDP(4) - t627 * MDP(5) + (-t746 - t763) * MDP(6) + (t619 * t713 + t599) * MDP(12) + (t623 * t713 - t696) * MDP(13) + (qJD(1) * t537 + t717) * MDP(19) + (qJD(1) * t536 + t655) * MDP(20) + (t647 * t726 + t438 * t546 - t500 * t506 + (-t501 * t617 - t622 * t671) * t769) * MDP(26) + (t647 * t740 + t437 * t546 + t500 * t653 + (-t501 * t622 + t617 * t671) * t769) * MDP(27) + (t546 * t404 - t500 * t452 - t548 * t736 + qJD(1) * t643 - ((t621 * t759 + t761) * t526 - t741) * t647) * MDP(33) + (t546 * t403 + t500 * t654 + t501 * t643 + t548 * t526 * qJD(1) - (-(t616 * t759 - t617 * t703 - t621 * t706) * t526 + t742) * t647) * MDP(34); (g(3) * t619 + t623 * t638 + t552) * MDP(12) + (t480 * t608 + (t537 * t710 + t607 * t753 - t608 * t707) * pkin(3) + t634) * MDP(19) + (t593 * t438 - t776 * t622 + t656 * t617 + t666 * t506 + (-t592 * t705 + t779) * t769 + t650) * MDP(26) + (t481 * t608 + (t536 * t710 - t607 * t618 - t608 * t683) * pkin(3) + t629) * MDP(20) + (t593 * t437 + t656 * t622 + t652 * t617 - t666 * t653 + (t592 * t706 + t762) * t769 + t651) * MDP(27) + ((t539 * t621 - t540 * t616) * t460 + t561 * t404 + (t616 * t660 - t621 * t662) * t526 + t715 * t452 + t631) * MDP(33) + (-(t539 * t616 + t540 * t621) * t460 + t561 * t403 + (t616 * t662 + t621 * t660) * t526 - t715 * t654 + t630) * MDP(34) + t628 + MDP(9) * t697 - MDP(10) * t698 + qJDD(3) * MDP(11) + (g(3) * t623 + (-t570 - t638) * t619) * MDP(13) + (MDP(7) * t619 * t623 - MDP(8) * t714) * t627; (t471 * t608 + t629) * MDP(20) + (-pkin(4) * t438 - t472 * t506 - t488 * t769 + (t471 * t769 + t658) * t617 + (-t776 - t691) * t622 + t650) * MDP(26) + (-pkin(4) * t437 + t721 * t769 + t472 * t653 + t658 * t622 + (t652 + t691) * t617 + t651) * MDP(27) + ((t563 * t621 - t564 * t616) * t460 + t594 * t404 + (t616 * t659 - t621 * t661) * t526 + t665 * t452 + t631) * MDP(33) + (-(t563 * t616 + t564 * t621) * t460 + t594 * t403 + (t616 * t661 + t621 * t659) * t526 - t665 * t654 + t630) * MDP(34) + t628 + (t472 * t608 + t634) * MDP(19); -t653 * t506 * MDP(21) + (-t506 ^ 2 + t653 ^ 2) * MDP(22) + (t437 + t781) * MDP(23) + (-t438 - t780) * MDP(24) + t461 * MDP(25) + (-g(1) * t529 - g(2) * t531 + t432 * t769 + t464 * t653 + (t673 + t588) * t617 + t663) * MDP(26) + (g(1) * t530 - g(2) * t532 + t431 * t769 + t464 * t506 + t603 * t751 - t641) * MDP(27) + (t403 + t778) * MDP(30) + (-t404 - t777) * MDP(31) + (-(-t419 * t616 - t744) * t526 - t399 * qJD(6) + (t452 * t653 + t621 * t460 - t526 * t704) * pkin(5) + t767) * MDP(33) + ((-t420 * t526 - t393) * t616 + (t419 * t526 - t672) * t621 + (-t616 * t460 - t526 * t703 - t653 * t654) * pkin(5) + t768) * MDP(34) + t766; (t690 + t778) * MDP(30) + (-t676 - t777) * MDP(31) + (t399 * t526 + t767) * MDP(33) + (-t616 * t393 - t621 * t394 + t398 * t526 + t768) * MDP(34) + (MDP(30) * t735 + MDP(31) * t654 - MDP(33) * t399 - MDP(34) * t745) * qJD(6) + t766;];
tau  = t1;
