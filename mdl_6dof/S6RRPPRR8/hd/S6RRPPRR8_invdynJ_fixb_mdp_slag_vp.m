% Calculate vector of inverse dynamics joint torques for
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:28
% EndTime: 2019-03-09 09:26:45
% DurationCPUTime: 13.47s
% Computational Cost: add. (5479->677), mult. (12387->856), div. (0->0), fcn. (9133->12), ass. (0->273)
t664 = cos(qJ(2));
t758 = qJD(1) * t664;
t624 = qJD(5) + t758;
t610 = qJD(6) + t624;
t662 = cos(qJ(6));
t655 = sin(pkin(10));
t660 = sin(qJ(2));
t759 = qJD(1) * t660;
t736 = t655 * t759;
t656 = cos(pkin(10));
t745 = t656 * qJD(2);
t582 = t736 - t745;
t734 = t656 * t759;
t757 = qJD(2) * t655;
t584 = t734 + t757;
t659 = sin(qJ(5));
t663 = cos(qJ(5));
t699 = -t663 * t582 + t584 * t659;
t658 = sin(qJ(6));
t698 = t582 * t659 + t584 * t663;
t787 = t698 * t658;
t449 = t662 * t699 + t787;
t821 = t449 * t610;
t702 = t658 * t699 - t662 * t698;
t820 = t610 * t702;
t590 = t655 * t659 + t656 * t663;
t685 = t664 * t590;
t768 = -qJD(1) * t685 - t590 * qJD(5);
t733 = t656 * t758;
t735 = t655 * t758;
t748 = qJD(5) * t663;
t749 = qJD(5) * t659;
t767 = t655 * t748 - t656 * t749 - t659 * t733 + t663 * t735;
t650 = g(3) * t664;
t661 = sin(qJ(1));
t665 = cos(qJ(1));
t709 = g(1) * t665 + g(2) * t661;
t811 = t709 * t660;
t819 = t650 - t811;
t743 = qJDD(1) * t660;
t631 = pkin(7) * t743;
t744 = qJD(1) * qJD(2);
t729 = t664 * t744;
t558 = -qJDD(2) * pkin(2) + pkin(7) * t729 + qJDD(3) + t631;
t727 = -t558 - t650;
t668 = -t811 - t727;
t763 = t664 * pkin(2) + t660 * qJ(3);
t815 = -pkin(1) - t763;
t575 = t815 * qJD(1);
t634 = pkin(7) * t758;
t604 = qJD(2) * qJ(3) + t634;
t516 = t575 * t656 - t655 * t604;
t495 = pkin(3) * t758 + qJD(4) - t516;
t461 = pkin(4) * t758 - pkin(8) * t584 + t495;
t517 = t655 * t575 + t656 * t604;
t501 = -qJ(4) * t758 + t517;
t464 = pkin(8) * t582 + t501;
t427 = t461 * t659 + t464 * t663;
t423 = -pkin(9) * t699 + t427;
t747 = qJD(6) * t658;
t421 = t423 * t747;
t596 = -qJD(2) * pkin(2) + pkin(7) * t759 + qJD(3);
t489 = t582 * pkin(3) - t584 * qJ(4) + t596;
t462 = -pkin(4) * t582 - t489;
t441 = pkin(5) * t699 + t462;
t776 = t661 * t664;
t567 = t655 * t776 + t656 * t665;
t774 = t665 * t655;
t568 = t656 * t776 - t774;
t653 = qJ(5) + qJ(6);
t640 = sin(t653);
t641 = cos(t653);
t480 = -t567 * t640 - t568 * t641;
t777 = t661 * t656;
t569 = t664 * t774 - t777;
t775 = t664 * t665;
t570 = t655 * t661 + t656 * t775;
t482 = t569 * t640 + t570 * t641;
t696 = t640 * t655 + t641 * t656;
t794 = g(3) * t660;
t818 = g(1) * t482 - g(2) * t480 + t441 * t449 + t696 * t794 + t421;
t479 = t567 * t641 - t568 * t640;
t481 = t569 * t641 - t570 * t640;
t697 = t640 * t656 - t641 * t655;
t636 = t656 * qJDD(2);
t684 = t729 + t743;
t541 = t655 * t684 - t636;
t542 = qJDD(2) * t655 + t656 * t684;
t442 = t659 * t541 + t663 * t542 + t582 * t748 - t584 * t749;
t639 = t664 * qJDD(1);
t730 = t660 * t744;
t683 = -t730 + t639;
t589 = -qJDD(5) - t683;
t706 = pkin(2) * t660 - qJ(3) * t664;
t561 = qJD(2) * t706 - qJD(3) * t660;
t506 = qJD(1) * t561 + qJDD(1) * t815;
t551 = pkin(7) * t683 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t456 = t506 * t656 - t655 * t551;
t689 = pkin(3) * t639 + qJDD(4) - t456;
t447 = -pkin(3) * t730 + t689;
t431 = pkin(4) * t683 - pkin(8) * t542 + t447;
t457 = t655 * t506 + t656 * t551;
t444 = qJ(4) * t730 + (-qJ(4) * qJDD(1) - qJD(1) * qJD(4)) * t664 + t457;
t434 = pkin(8) * t541 + t444;
t725 = t663 * t431 - t659 * t434;
t672 = -qJD(5) * t427 + t725;
t411 = -pkin(5) * t589 - pkin(9) * t442 + t672;
t443 = qJD(5) * t698 - t663 * t541 + t542 * t659;
t682 = -t659 * t431 - t663 * t434 - t461 * t748 + t464 * t749;
t412 = -pkin(9) * t443 - t682;
t726 = t662 * t411 - t658 * t412;
t817 = -g(1) * t481 - g(2) * t479 + t441 * t702 + t697 * t794 + t726;
t577 = -qJDD(6) + t589;
t816 = (-t449 ^ 2 + t702 ^ 2) * MDP(27) - t577 * MDP(30) - t449 * MDP(26) * t702;
t814 = t624 * t698;
t813 = t624 * t699;
t782 = t655 * t664;
t618 = pkin(7) * t782;
t646 = t664 * pkin(3);
t504 = pkin(4) * t664 + t618 + t646 + (-pkin(8) * t660 - t815) * t656;
t780 = t656 * t664;
t546 = pkin(7) * t780 + t655 * t815;
t792 = qJ(4) * t664;
t527 = t546 - t792;
t784 = t655 * t660;
t515 = pkin(8) * t784 + t527;
t770 = t659 * t504 + t663 * t515;
t751 = qJD(4) * t655;
t765 = qJ(4) * t733 - t634;
t800 = -pkin(3) - pkin(4);
t719 = -t735 * t800 + t751 - t765;
t793 = -pkin(8) + qJ(3);
t600 = t793 * t655;
t601 = t793 * t656;
t766 = t659 * t600 + t663 * t601;
t740 = -pkin(7) * t655 - pkin(3);
t674 = -pkin(8) * t780 + (-pkin(4) + t740) * t660;
t594 = t706 * qJD(1);
t785 = t594 * t656;
t476 = qJD(1) * t674 - t785;
t752 = qJD(3) * t663;
t809 = -t663 * t476 + t655 * t752;
t756 = qJD(2) * t660;
t808 = qJ(4) * t756 - qJD(4) * t664;
t571 = t655 * t594;
t629 = qJ(4) * t759;
t781 = t656 * t660;
t686 = -pkin(7) * t781 + pkin(8) * t782;
t497 = qJD(1) * t686 + t571 + t629;
t753 = qJD(3) * t659;
t807 = t659 * t476 + t663 * t497 - t600 * t748 + t601 * t749 - t655 * t753 - t656 * t752;
t806 = t656 * pkin(3) + t655 * qJ(4) + pkin(2);
t573 = t584 * qJD(4);
t805 = t541 * pkin(3) - t542 * qJ(4) - t573;
t724 = t658 * t442 + t662 * t443;
t670 = qJD(6) * t702 - t724;
t802 = t610 ^ 2;
t801 = -0.2e1 * pkin(1);
t580 = t584 ^ 2;
t799 = pkin(7) * t584;
t798 = g(1) * t661;
t795 = g(2) * t665;
t426 = t663 * t461 - t464 * t659;
t422 = -pkin(9) * t698 + t426;
t420 = pkin(5) * t624 + t422;
t791 = t420 * t662;
t790 = t423 * t662;
t789 = t444 * t656;
t788 = t457 * t656;
t786 = t561 * t656;
t783 = t655 * t663;
t779 = t660 * t661;
t778 = t660 * t665;
t591 = -t656 * t659 + t783;
t513 = t662 * t590 + t591 * t658;
t773 = -qJD(6) * t513 - t658 * t767 + t662 * t768;
t514 = -t590 * t658 + t591 * t662;
t772 = qJD(6) * t514 + t658 * t768 + t662 * t767;
t769 = pkin(5) * t767 + t719;
t732 = t664 * t745;
t764 = -qJ(4) * t732 - qJD(4) * t781;
t651 = t660 ^ 2;
t652 = t664 ^ 2;
t762 = t651 - t652;
t760 = qJD(1) * t656;
t755 = qJD(2) * t664;
t754 = qJD(3) * t656;
t746 = qJD(6) * t662;
t742 = pkin(7) * t756;
t741 = t662 * t442 - t658 * t443 - t699 * t746;
t739 = pkin(3) * t655 + pkin(7);
t738 = qJ(3) * t756;
t731 = qJ(3) * t639;
t446 = t558 + t805;
t728 = -t446 - t650;
t470 = qJD(2) * t674 - t786;
t550 = t655 * t561;
t471 = qJD(2) * t686 + t550 + t808;
t723 = t663 * t470 - t471 * t659;
t722 = t663 * t504 - t515 * t659;
t720 = t663 * t600 - t601 * t659;
t545 = t656 * t815 - t618;
t572 = t656 * pkin(4) + t806;
t718 = qJD(6) * t420 + t412;
t717 = g(1) * t656 * t778 + g(2) * t660 * t777 + qJD(3) * t735 + t655 * t731;
t716 = t665 * pkin(1) + pkin(2) * t775 + t661 * pkin(7) + qJ(3) * t778;
t715 = t656 * t731;
t714 = t655 * t800 - pkin(7);
t713 = -t656 * qJ(3) * t541 - t582 * t754 - t794;
t712 = t740 * t660;
t711 = -g(1) * t567 + g(2) * t569;
t710 = g(1) * t568 - g(2) * t570;
t485 = -pkin(9) * t590 + t766;
t708 = -pkin(5) * t759 + pkin(9) * t768 + qJD(5) * t766 + qJD(6) * t485 - t497 * t659 + t656 * t753 - t809;
t484 = -pkin(9) * t591 + t720;
t707 = pkin(9) * t767 - qJD(6) * t484 + t807;
t705 = pkin(7) * t582 + t596 * t655;
t414 = t420 * t658 + t790;
t555 = t659 * t781 - t660 * t783;
t556 = t590 * t660;
t477 = t662 * t555 + t556 * t658;
t478 = -t555 * t658 + t556 * t662;
t701 = t567 * t663 - t568 * t659;
t700 = t567 * t659 + t568 * t663;
t695 = t658 * t663 + t659 * t662;
t694 = t658 * t659 - t662 * t663;
t693 = qJ(3) * t542 + qJD(3) * t584;
t692 = t624 ^ 2;
t533 = -pkin(7) * t734 + t571;
t524 = -t656 * t742 + t550;
t687 = -pkin(7) * qJDD(2) + t744 * t801;
t681 = t659 * t470 + t663 * t471 + t504 * t748 - t515 * t749;
t679 = t698 * t747 - t741;
t667 = qJD(1) ^ 2;
t678 = pkin(1) * t667 + t709;
t666 = qJD(2) ^ 2;
t677 = pkin(7) * t666 + qJDD(1) * t801 + t795;
t616 = qJ(4) * t781;
t526 = t660 * t714 + t616;
t676 = t815 * t798;
t494 = t714 * t755 - t764;
t435 = -pkin(4) * t541 - t446;
t648 = t665 * pkin(7);
t626 = g(1) * t779;
t622 = qJ(3) * t775;
t619 = qJ(3) * t776;
t559 = t582 * t758;
t552 = t660 * t739 - t616;
t539 = pkin(3) * t735 - t765;
t532 = pkin(7) * t736 + t785;
t531 = -t545 + t646;
t523 = t655 * t742 + t786;
t522 = pkin(5) * t590 + t572;
t521 = qJD(1) * t712 - t785;
t520 = t533 + t629;
t519 = t739 * t755 + t764;
t507 = qJD(2) * t712 - t786;
t499 = t569 * t659 + t570 * t663;
t498 = t569 * t663 - t570 * t659;
t493 = t524 + t808;
t488 = qJD(5) * t591 * t660 + qJD(2) * t685;
t487 = qJD(5) * t556 + t659 * t732 - t755 * t783;
t474 = pkin(5) * t555 + t526;
t445 = pkin(5) * t487 + t494;
t437 = -pkin(9) * t555 + t770;
t436 = pkin(5) * t664 - pkin(9) * t556 + t722;
t425 = qJD(6) * t478 + t662 * t487 + t488 * t658;
t424 = -qJD(6) * t477 - t487 * t658 + t488 * t662;
t419 = pkin(5) * t443 + t435;
t418 = -pkin(9) * t487 + t681;
t417 = -pkin(5) * t756 - pkin(9) * t488 - qJD(5) * t770 + t723;
t413 = -t423 * t658 + t791;
t1 = [t709 * MDP(3) + (-t443 * t664 - t487 * t624 + t555 * t589 + t699 * t756) * MDP(22) + (t442 * t556 + t488 * t698) * MDP(19) + (-t442 * t555 - t443 * t556 - t487 * t698 - t488 * t699) * MDP(20) + (t442 * t664 + t488 * t624 - t556 * t589 - t698 * t756) * MDP(21) + (g(1) * t701 - g(2) * t498 + t427 * t756 + t435 * t556 + t526 * t442 + t462 * t488 + t494 * t698 + t589 * t770 - t624 * t681 + t664 * t682) * MDP(25) + (t660 * t677 + t664 * t687 - t626) * MDP(10) + (-t424 * t702 - t478 * t679) * MDP(26) + (t424 * t610 - t478 * t577 - t664 * t679 + t702 * t756) * MDP(28) + (t421 * t664 + t414 * t756 - t445 * t702 - t474 * t679 + t419 * t478 + t441 * t424 + g(1) * t479 - g(2) * t481 + (-(-qJD(6) * t437 + t417) * t610 + t436 * t577 - t411 * t664) * t658 + (-(qJD(6) * t436 + t418) * t610 + t437 * t577 - t718 * t664) * t662) * MDP(32) + ((pkin(7) * t542 + t558 * t656 + (-qJD(1) * t546 - t517) * qJD(2)) * t660 + (qJD(1) * t524 + qJDD(1) * t546 + t457 + (t596 * t656 + t799) * qJD(2)) * t664 + t711) * MDP(12) + ((pkin(7) * t541 + t558 * t655 + (qJD(1) * t545 + t516) * qJD(2)) * t660 + (-qJD(1) * t523 + qJD(2) * t705 - qJDD(1) * t545 - t456) * t664 + t710) * MDP(11) + (-t523 * t584 - t524 * t582 - t541 * t546 - t542 * t545 + t626 + (-t516 * t656 - t517 * t655) * t755 + (-t456 * t656 - t457 * t655 - t795) * t660) * MDP(13) + (-t493 * t582 + t507 * t584 - t527 * t541 + t531 * t542 + t626 + (t495 * t656 - t501 * t655) * t755 + (-t444 * t655 + t447 * t656 - t795) * t660) * MDP(16) + (t687 * t660 + (-t677 + t798) * t664) * MDP(9) + (-t795 + t798) * MDP(2) + (qJDD(2) * t660 + t664 * t666) * MDP(6) + (qJDD(2) * t664 - t660 * t666) * MDP(7) + (t444 * t527 + t501 * t493 + t446 * t552 + t489 * t519 + t447 * t531 + t495 * t507 - g(1) * (-pkin(3) * t568 - qJ(4) * t567 + t648) - g(2) * (pkin(3) * t570 + qJ(4) * t569 + t716) - t676) * MDP(18) + (-t424 * t449 + t425 * t702 + t477 * t679 + t478 * t670) * MDP(27) + ((t417 * t662 - t418 * t658) * t610 - (t436 * t662 - t437 * t658) * t577 + t726 * t664 - t413 * t756 + t445 * t449 - t474 * t670 + t419 * t477 + t441 * t425 - g(1) * t480 - g(2) * t482 + ((-t436 * t658 - t437 * t662) * t610 - t414 * t664) * qJD(6)) * MDP(31) + (-t425 * t610 + t449 * t756 + t477 * t577 + t664 * t670) * MDP(29) + (t723 * t624 - t722 * t589 + t725 * t664 - t426 * t756 + t494 * t699 + t526 * t443 + t435 * t555 + t462 * t487 + g(1) * t700 - g(2) * t499 + (-t427 * t664 - t624 * t770) * qJD(5)) * MDP(24) + (qJDD(1) * t651 + 0.2e1 * t660 * t729) * MDP(4) + qJDD(1) * MDP(1) + (-t519 * t584 - t542 * t552 + (-t446 * t656 + (qJD(1) * t527 + t501) * qJD(2)) * t660 + (-qJD(1) * t493 - qJDD(1) * t527 - t489 * t745 - t444) * t664 - t711) * MDP(17) + (t457 * t546 + t517 * t524 + t456 * t545 + t516 * t523 - g(1) * t648 - g(2) * t716 - t676 + (t558 * t660 + t596 * t755) * pkin(7)) * MDP(14) + (-t589 * t664 - t624 * t756) * MDP(23) + (-t577 * t664 - t610 * t756) * MDP(30) + (t519 * t582 + t541 * t552 + (t446 * t655 + (-qJD(1) * t531 - t495) * qJD(2)) * t660 + (qJD(1) * t507 + qJDD(1) * t531 + t489 * t757 + t447) * t664 + t710) * MDP(15) + 0.2e1 * (t639 * t660 - t744 * t762) * MDP(5); (-t514 * t577 + t610 * t773) * MDP(28) + (-t442 * t590 - t443 * t591 - t698 * t767 - t699 * t768) * MDP(20) + (t442 * t591 + t698 * t768) * MDP(19) + (-t589 * t591 + t624 * t768) * MDP(21) + MDP(7) * t639 + MDP(6) * t743 + (-MDP(4) * t660 * t664 + MDP(5) * t762) * t667 + (t715 - pkin(2) * t542 + t668 * t655 + ((-qJ(3) * t745 + t517) * t660 + (-t799 - t533 + (qJD(3) - t596) * t656) * t664) * qJD(1)) * MDP(12) + ((t484 * t658 + t485 * t662) * t577 - t522 * t679 + t419 * t514 + (t658 * t708 + t662 * t707) * t610 - t769 * t702 + t773 * t441 + t819 * t697) * MDP(32) + (t766 * t589 + t572 * t442 + t807 * t624 + t719 * t698 + t768 * t462 + (t435 - t819) * t591) * MDP(25) + (-(t484 * t662 - t485 * t658) * t577 - t522 * t670 + t419 * t513 + (t658 * t707 - t662 * t708) * t610 + t769 * t449 + t772 * t441 - t819 * t696) * MDP(31) + (t513 * t577 - t610 * t772) * MDP(29) + (t660 * t678 - t631 - t650) * MDP(9) + (-t720 * t589 + t572 * t443 - g(3) * t685 + (-t601 * t748 + (-qJD(5) * t600 + t497 - t754) * t659 + t809) * t624 + t719 * t699 + t767 * t462 + (t435 + t811) * t590) * MDP(24) + (t698 * MDP(21) - t699 * MDP(22) + t624 * MDP(23) + t426 * MDP(24) - t427 * MDP(25) - MDP(28) * t702 - t449 * MDP(29) + t610 * MDP(30) + t413 * MDP(31) - t414 * MDP(32)) * t759 + (-t514 * t679 - t702 * t773) * MDP(26) + (t794 + (-pkin(7) * qJDD(1) + t678) * t664) * MDP(10) + (t789 + t520 * t582 - t521 * t584 + (-t495 * t760 - t709) * t664 + (t501 * t758 + t447 + t693) * t655 + t713) * MDP(16) + (-t558 * pkin(2) - t517 * t533 - t516 * t532 - t596 * t634 - g(1) * (-pkin(2) * t778 + t622) - g(2) * (-pkin(2) * t779 + t619) - g(3) * t763 + (-t516 * t655 + t517 * t656) * qJD(3) + (-t456 * t655 + t788) * qJ(3)) * MDP(14) + (t788 + t532 * t584 + t533 * t582 + (t516 * t760 - t709) * t664 + (t517 * t758 - t456 + t693) * t655 + t713) * MDP(13) + (t589 * t590 - t624 * t767) * MDP(22) + (-t449 * t773 + t513 * t679 + t514 * t670 + t702 * t772) * MDP(27) + (-t715 + t539 * t584 + t542 * t806 + (t573 + t728 + t811) * t655 + (-t501 * t660 + t520 * t664 + (t738 + (-qJD(3) + t489) * t664) * t656) * qJD(1)) * MDP(17) + (-t541 * t806 + t728 * t656 + (-t539 - t751) * t582 + (t495 * t660 - t521 * t664 + (-t489 * t664 - t738) * t655) * qJD(1) + t717) * MDP(15) + qJDD(2) * MDP(8) + (-pkin(2) * t541 + t727 * t656 + ((-qJ(3) * t757 - t516) * t660 + (t532 - t705) * t664) * qJD(1) + t717) * MDP(11) + (qJ(3) * t789 - t489 * t539 - t495 * t521 - g(1) * t622 - g(2) * t619 - g(3) * (pkin(3) * t780 + t763) + (-t520 + t754) * t501 + (-g(3) * t792 + qJ(3) * t447 + qJD(3) * t495 - qJD(4) * t489) * t655 + (-t446 + t811) * t806) * MDP(18); (t516 * t584 + t517 * t582 + t668) * MDP(14) + (-t495 * t584 + t501 * t582 + t668 + t805) * MDP(18) + (-t443 - t814) * MDP(24) + (-t442 + t813) * MDP(25) + (t670 + t820) * MDP(31) + (t679 + t821) * MDP(32) + (MDP(11) + MDP(15)) * (t655 * t743 - t636 + (-t584 + t757) * t758) + (MDP(12) - MDP(17)) * (t559 + t542) + (MDP(16) + MDP(13)) * (-t582 ^ 2 - t580); (t582 * t584 + t683) * MDP(15) + (-t559 + t542) * MDP(16) + (-t652 * t667 - t580) * MDP(17) + (-g(3) * t784 - g(1) * t569 - g(2) * t567 + t489 * t584 + (-pkin(3) * t756 + t501 * t664) * qJD(1) + t689) * MDP(18) + (-t584 * t699 - t589 * t663 - t659 * t692) * MDP(24) + (-t584 * t698 + t589 * t659 - t663 * t692) * MDP(25) + (-t584 * t449 + t694 * t577 - t695 * t802) * MDP(31) + (t695 * t577 + t584 * t702 + t694 * t802) * MDP(32); t698 * t699 * MDP(19) + (t698 ^ 2 - t699 ^ 2) * MDP(20) + (t442 + t813) * MDP(21) + (-t443 + t814) * MDP(22) - t589 * MDP(23) + (-g(1) * t498 - g(2) * t701 + g(3) * t555 + t427 * t624 - t462 * t698 + t672) * MDP(24) + (g(1) * t499 + g(2) * t700 + g(3) * t556 + t426 * t624 + t462 * t699 + t682) * MDP(25) + (-t679 + t821) * MDP(28) + (t670 - t820) * MDP(29) + (-(-t422 * t658 - t790) * t610 - t414 * qJD(6) + (-t449 * t698 - t577 * t662 - t610 * t747) * pkin(5) + t817) * MDP(31) + ((-t423 * t610 - t411) * t658 + (t422 * t610 - t718) * t662 + (t577 * t658 - t610 * t746 + t698 * t702) * pkin(5) + t818) * MDP(32) + t816; (t741 + t821) * MDP(28) + (-t724 - t820) * MDP(29) + (t414 * t610 + t817) * MDP(31) + (-t658 * t411 - t662 * t412 + t413 * t610 + t818) * MDP(32) + (-MDP(28) * t787 + MDP(29) * t702 - MDP(31) * t414 - MDP(32) * t791) * qJD(6) + t816;];
tau  = t1;
