% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6PRRRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:18:07
% EndTime: 2019-03-09 00:18:23
% DurationCPUTime: 11.75s
% Computational Cost: add. (7490->647), mult. (16650->841), div. (0->0), fcn. (12816->14), ass. (0->278)
t638 = sin(qJ(3));
t641 = cos(qJ(3));
t701 = pkin(3) * t638 - pkin(9) * t641;
t584 = t701 * qJD(3);
t637 = sin(qJ(4));
t639 = sin(qJ(2));
t640 = cos(qJ(4));
t753 = qJD(3) * t638;
t634 = sin(pkin(6));
t760 = qJD(1) * t634;
t642 = cos(qJ(2));
t778 = t641 * t642;
t806 = pkin(8) * t637;
t842 = (-t637 * t778 + t639 * t640) * t760 - t640 * t584 - t753 * t806;
t591 = -pkin(3) * t641 - pkin(9) * t638 - pkin(2);
t747 = qJD(4) * t640;
t841 = -(t637 * t639 + t640 * t778) * t760 + t637 * t584 + t591 * t747;
t752 = qJD(3) * t640;
t757 = qJD(2) * t638;
t577 = -t637 * t757 + t752;
t754 = qJD(3) * t637;
t578 = t640 * t757 + t754;
t636 = sin(qJ(5));
t808 = cos(qJ(5));
t506 = -t808 * t577 + t578 * t636;
t683 = t636 * t577 + t578 * t808;
t840 = t506 * t683;
t779 = t640 * t641;
t618 = pkin(8) * t779;
t695 = pkin(4) * t638 - pkin(10) * t779;
t839 = t695 * qJD(3) + (-t618 + (pkin(10) * t638 - t591) * t637) * qJD(4) - t842;
t748 = qJD(4) * t637;
t751 = qJD(3) * t641;
t724 = t637 * t751;
t832 = t638 * t747 + t724;
t838 = -t832 * pkin(10) + (-t638 * t752 - t641 * t748) * pkin(8) + t841;
t581 = t701 * qJD(2);
t568 = t640 * t581;
t643 = -pkin(10) - pkin(9);
t731 = qJD(4) * t643;
t586 = qJD(2) * pkin(8) + t639 * t760;
t635 = cos(pkin(6));
t784 = t635 * t641;
t821 = qJD(1) * t784 - t638 * t586;
t837 = qJD(2) * t695 - t637 * t821 - t640 * t731 + t568;
t756 = qJD(2) * t641;
t725 = t637 * t756;
t767 = t637 * t581 + t640 * t821;
t836 = pkin(10) * t725 + t637 * t731 - t767;
t709 = qJD(4) + t756;
t737 = qJDD(2) * t638;
t663 = qJD(3) * t709 + t737;
t740 = qJD(2) * qJD(4);
t716 = t638 * t740;
t697 = -qJDD(3) + t716;
t677 = t697 * t637;
t648 = t663 * t640 - t677;
t694 = qJD(3) * qJD(4) + t737;
t741 = qJD(2) * qJD(3);
t718 = t641 * t741;
t662 = t694 + t718;
t706 = t637 * t662 + t640 * t716;
t668 = qJDD(3) * t640 - t706;
t720 = t808 * qJD(5);
t746 = qJD(5) * t636;
t442 = -t577 * t720 + t578 * t746 - t636 * t668 - t808 * t648;
t616 = -qJD(4) + t756;
t604 = -qJD(5) + t616;
t434 = -t506 * t604 - t442;
t443 = qJD(5) * t683 + t636 * t648 - t808 * t668;
t626 = t641 * qJDD(2);
t574 = t638 * t741 + qJDD(4) - t626;
t570 = qJDD(5) + t574;
t809 = t683 ^ 2;
t835 = t570 * MDP(23) + (-t604 * t683 - t443) * MDP(22) + MDP(19) * t840 + (-t506 ^ 2 + t809) * MDP(20) + t434 * MDP(21);
t529 = -qJD(3) * pkin(3) - t821;
t497 = -pkin(4) * t577 + t529;
t444 = pkin(5) * t506 - qJ(6) * t683 + t497;
t834 = t444 * t506;
t833 = t497 * t506;
t782 = t636 * t637;
t681 = t640 * t808 - t782;
t818 = qJD(4) + qJD(5);
t819 = t808 * qJD(4) + t720;
t769 = -t640 * t819 + t681 * t756 + t782 * t818;
t580 = t636 * t640 + t637 * t808;
t517 = t818 * t580;
t768 = -t580 * t756 + t517;
t461 = pkin(5) * t683 + qJ(6) * t506;
t828 = pkin(4) * t637 + pkin(8);
t598 = t643 * t637;
t599 = t643 * t640;
t682 = t598 * t808 + t636 * t599;
t827 = qJD(5) * t682 - t837 * t636 + t836 * t808;
t533 = t636 * t598 - t599 * t808;
t826 = qJD(5) * t533 + t836 * t636 + t837 * t808;
t759 = qJD(1) * t638;
t611 = t635 * t759;
t540 = t641 * t586 + t611;
t530 = qJD(3) * pkin(9) + t540;
t729 = t642 * t760;
t542 = qJD(2) * t591 - t729;
t484 = t640 * t530 + t637 * t542;
t460 = pkin(10) * t577 + t484;
t576 = t640 * t591;
t780 = t638 * t640;
t513 = -pkin(10) * t780 + t576 + (-pkin(4) - t806) * t641;
t763 = t637 * t591 + t618;
t781 = t637 * t638;
t526 = -pkin(10) * t781 + t763;
t825 = t513 * t720 - t526 * t746 + t839 * t636 + t838 * t808;
t703 = -t540 + (-t725 + t748) * pkin(4);
t824 = t636 * t513 + t808 * t526;
t541 = pkin(4) * t832 + pkin(8) * t751;
t705 = t638 * t729;
t823 = t541 - t705;
t555 = t570 * qJ(6);
t590 = t604 * qJD(6);
t822 = t555 - t590;
t797 = cos(pkin(11));
t712 = t797 * t639;
t633 = sin(pkin(11));
t788 = t633 * t642;
t558 = t635 * t712 + t788;
t713 = t634 * t797;
t519 = t558 * t641 - t638 * t713;
t711 = t797 * t642;
t789 = t633 * t639;
t560 = -t635 * t789 + t711;
t521 = t633 * t634 * t638 + t560 * t641;
t786 = t634 * t641;
t565 = t635 * t638 + t639 * t786;
t785 = t634 * t642;
t524 = -t565 * t637 - t640 * t785;
t557 = -t635 * t711 + t789;
t559 = t635 * t788 + t712;
t820 = -g(1) * (-t521 * t637 + t559 * t640) - g(2) * (-t519 * t637 + t557 * t640) - g(3) * t524;
t556 = t570 * pkin(5);
t817 = t556 - qJDD(6);
t787 = t634 * t639;
t564 = t638 * t787 - t784;
t816 = g(3) * t564 - g(2) * (-t558 * t638 - t641 * t713) - g(1) * (-t560 * t638 + t633 * t786);
t632 = qJ(4) + qJ(5);
t627 = sin(t632);
t628 = cos(t632);
t478 = t519 * t627 - t557 * t628;
t480 = t521 * t627 - t559 * t628;
t511 = t565 * t627 + t628 * t785;
t742 = qJD(1) * qJD(2);
t548 = qJDD(2) * pkin(8) + (qJDD(1) * t639 + t642 * t742) * t634;
t738 = qJDD(1) * t635;
t715 = t638 * t738;
t469 = qJDD(3) * pkin(9) + qJD(3) * t821 + t548 * t641 + t715;
t719 = t639 * t742;
t698 = -qJDD(1) * t785 + t634 * t719;
t492 = qJD(2) * t584 + qJDD(2) * t591 + t698;
t491 = t640 * t492;
t710 = -t637 * t469 + t491;
t735 = t637 * qJDD(3);
t736 = qJDD(2) * t640;
t423 = -(t638 * t736 + t640 * t718 + t735) * pkin(10) + t574 * pkin(4) - t460 * qJD(4) + t710;
t732 = t640 * t469 + t637 * t492 + t542 * t747;
t675 = -t530 * t748 + t732;
t430 = pkin(10) * t668 + t675;
t483 = -t530 * t637 + t640 * t542;
t459 = -pkin(10) * t578 + t483;
t450 = -pkin(4) * t616 + t459;
t707 = -t808 * t423 + t636 * t430 + t450 * t746 + t460 * t720;
t657 = g(1) * t480 + g(2) * t478 + g(3) * t511 - t707;
t650 = t444 * t683 - t657 - t817;
t815 = -t497 * t683 + t657;
t814 = -t533 * t570 - t627 * t816;
t644 = qJD(3) ^ 2;
t801 = g(2) * t557;
t803 = g(1) * t559;
t700 = t801 + t803;
t813 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t644 + t634 * (-g(3) * t642 + t719) - t698 + t700;
t666 = g(3) * t785 - t700;
t812 = qJD(4) * (pkin(8) * t616 + t530) - t666;
t811 = qJD(5) * t824 + t838 * t636 - t839 * t808;
t798 = qJD(2) * pkin(2);
t795 = qJDD(3) * pkin(3);
t730 = t808 * t460;
t436 = t636 * t450 + t730;
t794 = t436 * t604;
t792 = t578 * t616;
t791 = t627 * t641;
t790 = t628 * t641;
t783 = t636 * t460;
t777 = qJDD(1) - g(3);
t776 = qJ(6) * t753 - qJD(6) * t641 + t825;
t775 = -pkin(5) * t753 + t811;
t774 = pkin(5) * t768 + qJ(6) * t769 - qJD(6) * t580 + t703;
t772 = -qJ(6) * t757 + t827;
t771 = pkin(5) * t757 + t826;
t438 = t459 * t808 - t783;
t764 = pkin(4) * t720 + qJD(6) - t438;
t585 = pkin(4) * t781 + t638 * pkin(8);
t630 = t638 ^ 2;
t762 = -t641 ^ 2 + t630;
t758 = qJD(2) * t634;
t755 = qJD(3) * t577;
t750 = qJD(4) * t577;
t749 = qJD(4) * t616;
t745 = t529 * qJD(4);
t744 = t578 * qJD(3);
t435 = t450 * t808 - t783;
t743 = qJD(6) - t435;
t733 = t627 * t785;
t624 = pkin(4) * t640 + pkin(3);
t727 = t639 * t758;
t726 = t642 * t758;
t723 = t638 * t748;
t717 = t642 * t741;
t708 = t636 * t423 + t808 * t430 + t450 * t720 - t460 * t746;
t704 = t808 * t751;
t437 = t636 * t459 + t730;
t702 = pkin(4) * t746 - t437;
t699 = g(1) * t560 + g(2) * t558;
t645 = qJD(2) ^ 2;
t696 = qJDD(2) * t642 - t639 * t645;
t471 = t517 * t638 + t636 * t724 - t640 * t704;
t472 = t637 * t704 - t636 * t723 - t746 * t781 + (t636 * t751 + t638 * t819) * t640;
t552 = t681 * t638;
t432 = pkin(5) * t472 + qJ(6) * t471 - qJD(6) * t552 + t541;
t693 = -t432 + t705;
t692 = t624 * t641 - t638 * t643 + pkin(2);
t689 = -t565 * t640 + t637 * t785;
t687 = t513 * t808 - t636 * t526;
t685 = t524 * t808 + t636 * t689;
t465 = t636 * t524 - t689 * t808;
t679 = t574 * t637 - t616 * t747;
t678 = t640 * t574 + t616 * t748;
t676 = -t637 * t740 + t736;
t493 = -t557 * t791 - t558 * t628;
t495 = -t559 * t791 - t560 * t628;
t534 = -t628 * t787 + t641 * t733;
t673 = g(1) * t495 + g(2) * t493 + g(3) * t534;
t494 = -t557 * t790 + t558 * t627;
t496 = -t559 * t790 + t560 * t627;
t535 = (t627 * t639 + t628 * t778) * t634;
t672 = -g(1) * t496 - g(2) * t494 - g(3) * t535;
t670 = g(1) * t521 + g(2) * t519 + g(3) * t565;
t669 = qJD(3) * t611 + t638 * t548 + t586 * t751 - t641 * t738;
t665 = -g(3) * t787 - t699;
t664 = -pkin(9) * t574 - t529 * t616;
t470 = t669 - t795;
t479 = t519 * t628 + t557 * t627;
t481 = t521 * t628 + t559 * t627;
t512 = t565 * t628 - t733;
t658 = g(1) * t481 + g(2) * t479 + g(3) * t512 - t708;
t656 = -pkin(9) * t749 + t470 - t816;
t655 = t570 * t682 + t628 * t816;
t587 = -t729 - t798;
t654 = -pkin(8) * qJDD(3) + (t587 + t729 - t798) * qJD(3);
t652 = -t435 * t604 + t658;
t649 = -g(1) * (-t480 * pkin(5) + qJ(6) * t481) - g(2) * (-t478 * pkin(5) + qJ(6) * t479) - g(3) * (-t511 * pkin(5) + qJ(6) * t512);
t445 = -pkin(4) * t668 + t470;
t623 = -pkin(4) * t808 - pkin(5);
t620 = pkin(4) * t636 + qJ(6);
t551 = t580 * t638;
t523 = qJD(3) * t565 + t638 * t726;
t522 = -qJD(3) * t564 + t641 * t726;
t503 = -pkin(5) * t681 - qJ(6) * t580 - t624;
t488 = pkin(5) * t551 - qJ(6) * t552 + t585;
t458 = t641 * pkin(5) - t687;
t457 = -qJ(6) * t641 + t824;
t455 = qJD(4) * t524 + t522 * t640 + t637 * t727;
t454 = qJD(4) * t689 - t522 * t637 + t640 * t727;
t447 = pkin(4) * t578 + t461;
t433 = -t604 * qJ(6) + t436;
t431 = t604 * pkin(5) + t743;
t425 = qJD(5) * t465 - t454 * t808 + t636 * t455;
t424 = qJD(5) * t685 + t636 * t454 + t455 * t808;
t420 = t443 * pkin(5) + t442 * qJ(6) - qJD(6) * t683 + t445;
t419 = t707 - t817;
t418 = t708 + t822;
t1 = [t777 * MDP(1) + (-qJD(3) * t523 - qJDD(3) * t564) * MDP(10) + (-qJD(3) * t522 - qJDD(3) * t565) * MDP(11) + (-t454 * t616 - t523 * t577 + t524 * t574 - t564 * t668) * MDP(17) + (t455 * t616 + t523 * t578 + t564 * t648 + t574 * t689) * MDP(18) + (-t424 * t506 + t425 * t683 + t442 * t685 - t443 * t465) * MDP(27) + (t418 * t465 - t419 * t685 + t420 * t564 + t424 * t433 + t425 * t431 + t444 * t523 - g(3)) * MDP(29) + (MDP(24) + MDP(26)) * (t425 * t604 + t564 * t443 + t523 * t506 + t570 * t685) + (MDP(25) - MDP(28)) * (t424 * t604 - t442 * t564 - t465 * t570 + t523 * t683) + (t696 * MDP(3) + (-qJDD(2) * t639 - t642 * t645) * MDP(4) + (-t638 * t717 + t641 * t696) * MDP(10) + (-t638 * t696 - t641 * t717) * MDP(11)) * t634; (-t436 * t753 - t585 * t442 + t445 * t552 - t497 * t471 - t570 * t824 + t604 * t825 + t708 * t641 + t683 * t823 + t673) * MDP(25) + (t418 * t457 + t420 * t488 + t444 * t432 + t419 * t458 - g(1) * (pkin(5) * t496 + qJ(6) * t495 + t560 * t828) - g(2) * (pkin(5) * t494 + qJ(6) * t493 + t558 * t828) - g(3) * (pkin(5) * t535 + qJ(6) * t534) + t692 * t803 + t692 * t801 + t776 * t433 + t775 * t431 + (-g(3) * t828 * t639 + (-g(3) * t692 - t444 * t759) * t642) * t634) * MDP(29) + (t654 * t638 + t641 * t813) * MDP(10) + (-t638 * t813 + t654 * t641) * MDP(11) + (t435 * t753 + t585 * t443 + t445 * t551 + t497 * t472 + t506 * t823 + t687 * t570 + t604 * t811 + t707 * t641 + t672) * MDP(24) + (-t763 * t574 + t841 * t616 + t665 * t640 + ((pkin(8) * t578 + t529 * t640) * qJD(3) - t812 * t637 + t732) * t641 + (-t578 * t729 - t637 * t745 - t484 * qJD(3) + t470 * t640 + (t735 + t676 * t638 + (-t616 + t709) * t752) * pkin(8)) * t638) * MDP(18) + (t576 * t574 + t842 * t616 + (t591 * t749 + t665) * t637 + (-pkin(8) * t755 - t491 + (-pkin(8) * t574 + qJD(3) * t529 + qJD(4) * t542 + t469) * t637 + t812 * t640) * t641 + (-pkin(8) * t668 + t483 * qJD(3) + t470 * t637 + t577 * t729 + t640 * t745) * t638) * MDP(17) + (-t777 * t787 + t699) * MDP(4) + ((t577 * t640 - t578 * t637) * t751 + ((-t578 * qJD(4) + t668) * t640 + (-t750 - t648) * t637) * t638) * MDP(13) + ((t616 * t754 - t668) * t641 + (-t679 + t755) * t638) * MDP(15) + ((-t735 + (-t616 - t709) * t752) * t641 + (-t641 * t676 + t678 + t744) * t638) * MDP(14) + (-t574 * t641 - t616 * t753) * MDP(16) + (-t570 * t641 - t604 * t753) * MDP(23) + (t443 * t641 + t472 * t604 - t506 * t753 - t551 * t570) * MDP(22) + (t777 * t785 + t700) * MDP(3) + (qJDD(2) * t630 + 0.2e1 * t638 * t718) * MDP(5) + (t442 * t641 + t471 * t604 + t552 * t570 + t683 * t753) * MDP(21) + (-t418 * t641 - t420 * t552 + t433 * t753 + t442 * t488 + t444 * t471 + t457 * t570 - t604 * t776 + t683 * t693 - t673) * MDP(28) + (-t418 * t551 + t419 * t552 - t431 * t471 - t433 * t472 - t442 * t458 - t443 * t457 - t506 * t776 - t638 * t666 + t683 * t775) * MDP(27) + (t442 * t551 - t443 * t552 + t471 * t506 - t472 * t683) * MDP(20) + (-t442 * t552 - t471 * t683) * MDP(19) + 0.2e1 * (t626 * t638 - t741 * t762) * MDP(6) + (-t578 * t723 + (-t638 * t677 + t641 * t744 + t662 * t780) * t640) * MDP(12) + (t419 * t641 + t420 * t551 - t431 * t753 + t443 * t488 + t444 * t472 - t458 * t570 - t506 * t693 + t604 * t775 + t672) * MDP(26) + qJDD(2) * MDP(2) + (qJDD(3) * t638 + t641 * t644) * MDP(7) + (qJDD(3) * t641 - t638 * t644) * MDP(8); (-pkin(3) * t706 + t568 * t616 + t540 * t577 + (-t616 * t821 + t664) * t637 + (-t656 + t795) * t640) * MDP(17) + (-t697 * t637 ^ 2 + (t637 * t663 - t792) * t640) * MDP(12) + ((-t706 + t792) * t637 + (t750 + 0.2e1 * t735 + t694 * t640 + (-t723 + (-t577 + t752) * t641) * qJD(2)) * t640) * MDP(13) + (t624 * t442 + t445 * t580 - t769 * t497 + t604 * t827 + t703 * t683 + t814) * MDP(25) + (-t420 * t580 + t442 * t503 + t444 * t769 - t604 * t772 - t683 * t774 - t814) * MDP(28) + (-t715 + (-qJD(2) * t587 - t548) * t641 + t670) * MDP(11) + (qJD(3) * t540 - t669 + t816) * MDP(10) + ((-t616 * t637 * t641 - t577 * t638) * qJD(2) + t678) * MDP(15) + (-t767 * t616 - t540 * t578 + (-pkin(3) * t662 + t664) * t640 + (pkin(3) * t697 + t656) * t637) * MDP(18) + (-t624 * t443 - t445 * t681 + t768 * t497 + t703 * t506 + t604 * t826 + t655) * MDP(24) + (t570 * t681 + t604 * t768) * MDP(22) + (t570 * t580 + t604 * t769) * MDP(21) + (-t442 * t681 - t443 * t580 + t506 * t769 - t683 * t768) * MDP(20) + (-t442 * t580 - t683 * t769) * MDP(19) + (t418 * t681 + t419 * t580 - t431 * t769 - t433 * t768 + t442 * t682 - t443 * t533 - t506 * t772 + t683 * t771 - t670) * MDP(27) + (-t420 * t681 + t443 * t503 + t444 * t768 + t506 * t774 + t604 * t771 + t655) * MDP(26) + ((-t578 * t638 + t616 * t779) * qJD(2) + t679) * MDP(14) + qJDD(3) * MDP(9) + MDP(8) * t626 + MDP(7) * t737 + (-MDP(10) * t587 + t616 * MDP(16) - t483 * MDP(17) + t484 * MDP(18) - MDP(21) * t683 + t506 * MDP(22) + t604 * MDP(23) - t435 * MDP(24) + t436 * MDP(25) + t431 * MDP(26) - t433 * MDP(28)) * t757 + (-MDP(5) * t638 * t641 + MDP(6) * t762) * t645 + (t418 * t533 - t419 * t682 + t420 * t503 + t431 * t771 + t433 * t772 + t444 * t774 + t643 * t670 + t816 * (pkin(5) * t628 + qJ(6) * t627 + t624)) * MDP(29); -t578 * t577 * MDP(12) + (-t577 ^ 2 + t578 ^ 2) * MDP(13) + (t577 * t616 + t648) * MDP(14) + (t668 - t792) * MDP(15) + t574 * MDP(16) + (-t529 * t578 + t710 + (-qJD(4) - t616) * t484 + t820) * MDP(17) + (-t483 * t616 - t529 * t577 - g(1) * (-t521 * t640 - t559 * t637) - g(2) * (-t519 * t640 - t557 * t637) - g(3) * t689 - t675) * MDP(18) + (-t437 * t604 + (-t506 * t578 + t570 * t808 + t604 * t746) * pkin(4) + t815) * MDP(24) + (-t438 * t604 + t833 + (-t570 * t636 - t578 * t683 + t604 * t720) * pkin(4) + t658) * MDP(25) + (-t447 * t506 - t570 * t623 + t604 * t702 - t650) * MDP(26) + (-t442 * t623 - t443 * t620 + (t433 + t702) * t683 + (t431 - t764) * t506) * MDP(27) + (t447 * t683 + t570 * t620 - t604 * t764 - t658 + t822 - t834) * MDP(28) + (t418 * t620 + t419 * t623 - t444 * t447 - t431 * t437 + t764 * t433 + (t431 * t746 + t820) * pkin(4) + t649) * MDP(29) + t835; (-t794 + t815) * MDP(24) + (t652 + t833) * MDP(25) + (-t461 * t506 + t556 - t650 - t794) * MDP(26) + (pkin(5) * t442 - qJ(6) * t443 + (t433 - t436) * t683 + (t431 - t743) * t506) * MDP(27) + (t461 * t683 + 0.2e1 * t555 - 0.2e1 * t590 - t652 - t834) * MDP(28) + (-t419 * pkin(5) + t418 * qJ(6) - t431 * t436 + t433 * t743 - t444 * t461 + t649) * MDP(29) + t835; (-t570 + t840) * MDP(26) + t434 * MDP(27) + (-t604 ^ 2 - t809) * MDP(28) + (t433 * t604 + t650) * MDP(29);];
tau  = t1;
