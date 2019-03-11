% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRPRRP8_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:25:39
% EndTime: 2019-03-09 12:25:58
% DurationCPUTime: 14.73s
% Computational Cost: add. (12544->691), mult. (28157->846), div. (0->0), fcn. (20951->14), ass. (0->283)
t698 = sin(pkin(10));
t701 = sin(qJ(4));
t699 = cos(pkin(10));
t704 = cos(qJ(4));
t826 = t699 * t704;
t641 = t698 * t701 - t826;
t705 = cos(qJ(2));
t740 = t641 * t705;
t807 = qJD(1) * t740 - t641 * qJD(4);
t702 = sin(qJ(2));
t799 = qJD(1) * t702;
t779 = t698 * t799;
t790 = t699 * qJD(2);
t633 = t779 - t790;
t778 = t699 * t799;
t797 = qJD(2) * t698;
t635 = t778 + t797;
t561 = t633 * t701 - t635 * t704;
t700 = sin(qJ(5));
t756 = -t633 * t704 - t701 * t635;
t846 = cos(qJ(5));
t516 = t700 * t561 + t756 * t846;
t885 = t516 ^ 2;
t798 = qJD(1) * t705;
t672 = -qJD(4) + t798;
t665 = -qJD(5) + t672;
t884 = t516 * t665;
t642 = t698 * t704 + t699 * t701;
t866 = t642 * t705;
t595 = qJD(1) * t866;
t732 = t642 * qJD(4);
t857 = t595 - t732;
t758 = pkin(2) * t702 - qJ(3) * t705;
t647 = t758 * qJD(1);
t580 = pkin(7) * t779 + t699 * t647;
t825 = t699 * t705;
t754 = pkin(3) * t702 - pkin(8) * t825;
t550 = qJD(1) * t754 + t580;
t624 = t698 * t647;
t827 = t699 * t702;
t828 = t698 * t705;
t743 = -pkin(7) * t827 - pkin(8) * t828;
t567 = qJD(1) * t743 + t624;
t835 = pkin(8) + qJ(3);
t656 = t835 * t698;
t657 = t835 * t699;
t805 = -t701 * t656 + t704 * t657;
t883 = t642 * qJD(3) + qJD(4) * t805 + t704 * t550 - t567 * t701;
t792 = qJD(4) * t704;
t882 = qJD(3) * t826 - t701 * t550 - t704 * t567 - t656 * t792;
t822 = t701 * t657;
t768 = -t704 * t656 - t822;
t843 = pkin(9) * t642;
t544 = t768 - t843;
t545 = -pkin(9) * t641 + t805;
t507 = t700 * t544 + t545 * t846;
t687 = t705 * qJDD(1);
t788 = qJD(1) * qJD(2);
t737 = t702 * t788 - t687;
t640 = qJDD(4) + t737;
t627 = qJDD(5) + t640;
t695 = pkin(10) + qJ(4);
t688 = qJ(5) + t695;
t674 = sin(t688);
t693 = g(3) * t705;
t703 = sin(qJ(1));
t706 = cos(qJ(1));
t761 = g(1) * t706 + g(2) * t703;
t751 = t761 * t702;
t725 = t751 - t693;
t881 = t507 * t627 + t674 * t725;
t870 = -t561 * t846 + t700 * t756;
t847 = t870 ^ 2;
t880 = t516 * t870;
t879 = t665 * t870;
t878 = pkin(4) * t799 + pkin(9) * t807 + t883;
t794 = qJD(3) * t698;
t877 = pkin(9) * t595 - t701 * t794 + (-t822 - t843) * qJD(4) + t882;
t683 = t699 * qJDD(2);
t775 = t705 * t788;
t787 = qJDD(1) * t702;
t738 = t775 + t787;
t864 = t698 * t738;
t587 = -t683 + t864;
t786 = qJDD(2) * t698;
t716 = t699 * t738 + t786;
t793 = qJD(4) * t701;
t509 = -t701 * t587 - t633 * t792 - t635 * t793 + t704 * t716;
t753 = t587 * t704 - t633 * t793 + t635 * t792 + t701 * t716;
t777 = qJD(5) * t846;
t791 = qJD(5) * t700;
t460 = -t846 * t509 - t561 * t791 + t700 * t753 - t756 * t777;
t453 = -t460 + t884;
t461 = qJD(5) * t870 + t700 * t509 + t846 * t753;
t876 = t627 * MDP(26) + (-t461 - t879) * MDP(25) - MDP(22) * t880 + (t847 - t885) * MDP(23) + t453 * MDP(24);
t867 = t641 * t702;
t875 = pkin(9) * t867;
t774 = -qJD(2) * pkin(2) + qJD(3);
t650 = pkin(7) * t799 + t774;
t579 = pkin(3) * t633 + t650;
t525 = -pkin(4) * t756 + t579;
t466 = -pkin(5) * t516 - qJ(6) * t870 + t525;
t874 = t466 * t516;
t873 = t516 * t525;
t872 = t561 * t672;
t811 = t641 * t777 + t642 * t791 - t700 * t857 - t807 * t846;
t566 = -t700 * t641 + t642 * t846;
t810 = qJD(5) * t566 + t700 * t807 - t846 * t857;
t477 = pkin(5) * t870 - qJ(6) * t516;
t865 = t672 * t756;
t745 = t544 * t846 - t700 * t545;
t861 = qJD(5) * t745 - t700 * t878 + t846 * t877;
t860 = qJD(5) * t507 + t700 * t877 + t846 * t878;
t759 = pkin(2) * t705 + qJ(3) * t702;
t652 = -pkin(1) - t759;
t632 = t699 * t652;
t568 = -pkin(8) * t827 + t632 + (-pkin(7) * t698 - pkin(3)) * t705;
t593 = pkin(7) * t825 + t698 * t652;
t829 = t698 * t702;
t575 = -pkin(8) * t829 + t593;
t823 = t701 * t575;
t769 = t704 * t568 - t823;
t499 = -pkin(4) * t705 + t769 + t875;
t610 = t642 * t702;
t808 = t701 * t568 + t704 * t575;
t508 = -pkin(9) * t610 + t808;
t859 = t700 * t499 + t846 * t508;
t681 = pkin(7) * t798;
t784 = pkin(3) * t798;
t628 = t698 * t784 + t681;
t858 = -pkin(4) * t857 - t628;
t620 = t627 * qJ(6);
t655 = t665 * qJD(6);
t856 = t620 - t655;
t855 = pkin(7) * t775 + qJDD(3);
t622 = t627 * pkin(5);
t854 = t622 - qJDD(6);
t685 = sin(t695);
t686 = cos(t695);
t819 = t703 * t705;
t601 = t685 * t819 + t686 * t706;
t818 = t705 * t706;
t603 = -t685 * t818 + t686 * t703;
t836 = g(3) * t702;
t853 = -g(1) * t603 + g(2) * t601 + t685 * t836;
t675 = cos(t688);
t588 = t674 * t819 + t675 * t706;
t817 = t706 * t674;
t820 = t703 * t675;
t590 = t705 * t817 - t820;
t626 = t652 * qJD(1);
t658 = qJD(2) * qJ(3) + t681;
t569 = t699 * t626 - t658 * t698;
t529 = -pkin(8) * t635 + t569 - t784;
t570 = t698 * t626 + t699 * t658;
t532 = -pkin(8) * t633 + t570;
t491 = t529 * t701 + t532 * t704;
t617 = qJD(2) * t758 - qJD(3) * t702;
t560 = qJD(1) * t617 + qJDD(1) * t652;
t600 = -pkin(7) * t737 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t521 = t699 * t560 - t698 * t600;
t496 = pkin(3) * t737 - pkin(8) * t716 + t521;
t522 = t698 * t560 + t699 * t600;
t510 = -pkin(8) * t587 + t522;
t773 = t704 * t496 - t510 * t701;
t717 = -qJD(4) * t491 + t773;
t448 = pkin(4) * t640 - pkin(9) * t509 + t717;
t735 = t701 * t496 + t704 * t510 + t529 * t792 - t532 * t793;
t451 = -pkin(9) * t753 + t735;
t490 = t704 * t529 - t532 * t701;
t481 = pkin(9) * t561 + t490;
t478 = -pkin(4) * t672 + t481;
t482 = pkin(9) * t756 + t491;
t766 = -t846 * t448 + t700 * t451 + t478 * t791 + t482 * t777;
t831 = t674 * t702;
t719 = g(1) * t590 + g(2) * t588 + g(3) * t831 - t766;
t712 = t466 * t870 - t719 - t854;
t852 = -t525 * t870 + t719;
t546 = -qJD(2) * t740 - t702 * t732;
t796 = qJD(2) * t702;
t783 = pkin(7) * t796;
t573 = t699 * t617 + t698 * t783;
t540 = qJD(2) * t754 + t573;
t598 = t698 * t617;
t552 = qJD(2) * t743 + t598;
t771 = t704 * t540 - t552 * t701;
t468 = pkin(4) * t796 - pkin(9) * t546 - qJD(4) * t808 + t771;
t728 = qJD(2) * t866;
t782 = t701 * t540 + t704 * t552 + t568 * t792;
t471 = -pkin(9) * t728 + (-t823 + t875) * qJD(4) + t782;
t850 = -qJD(5) * t859 + t468 * t846 - t700 * t471;
t848 = -0.2e1 * pkin(1);
t845 = pkin(3) * t698;
t844 = pkin(7) * t633;
t841 = g(1) * t703;
t837 = g(2) * t706;
t834 = qJDD(2) * pkin(2);
t780 = t846 * t482;
t457 = t700 * t478 + t780;
t833 = t457 * t665;
t830 = t675 * t702;
t824 = t700 * t482;
t821 = t702 * t706;
t816 = -qJ(6) * t799 + t861;
t815 = pkin(5) * t799 + t860;
t814 = pkin(5) * t810 + qJ(6) * t811 - t566 * qJD(6) + t858;
t459 = t481 * t846 - t824;
t804 = pkin(4) * t777 + qJD(6) - t459;
t679 = pkin(7) * t787;
t803 = -t679 - t693;
t795 = qJD(2) * t705;
t629 = (pkin(7) + t845) * t795;
t648 = pkin(3) * t829 + t702 * pkin(7);
t802 = t706 * pkin(1) + t703 * pkin(7);
t696 = t702 ^ 2;
t801 = -t705 ^ 2 + t696;
t456 = t478 * t846 - t824;
t789 = qJD(6) - t456;
t781 = pkin(3) * t699 + pkin(2);
t767 = t700 * t448 + t846 * t451 + t478 * t777 - t482 * t791;
t673 = t702 * t841;
t765 = -g(2) * t821 + t673;
t572 = pkin(4) * t610 + t648;
t458 = t700 * t481 + t780;
t764 = pkin(4) * t791 - t458;
t763 = -g(1) * t588 + g(2) * t590;
t589 = t675 * t819 - t817;
t591 = t674 * t703 + t675 * t818;
t762 = g(1) * t589 - g(2) * t591;
t760 = -t837 + t841;
t597 = pkin(4) * t641 - t781;
t646 = pkin(4) * t686 + t781;
t752 = pkin(5) * t675 + qJ(6) * t674 + t646;
t612 = t679 - t834 + t855;
t750 = -pkin(7) * qJDD(2) + t788 * t848;
t748 = t499 * t846 - t700 * t508;
t539 = -t700 * t610 - t846 * t867;
t736 = t699 * t787 + t786;
t734 = t700 * t468 + t846 * t471 + t499 * t777 - t508 * t791;
t708 = qJD(1) ^ 2;
t731 = pkin(1) * t708 + t761;
t707 = qJD(2) ^ 2;
t730 = pkin(7) * t707 + qJDD(1) * t848 + t837;
t727 = qJD(4) * t867;
t726 = g(2) * t702 * t820 + t745 * t627 + (g(1) * t821 - t693) * t675;
t541 = t587 * pkin(3) + t612;
t723 = -t705 * t761 - t836;
t722 = -t751 - t834;
t720 = g(1) * t591 + g(2) * t589 + g(3) * t830 - t767;
t718 = -t612 + t725;
t713 = -t456 * t665 + t720;
t486 = pkin(4) * t753 + t541;
t710 = t727 - t728;
t709 = -g(1) * (-t590 * pkin(5) + qJ(6) * t591) - g(2) * (-t588 * pkin(5) + qJ(6) * t589) - g(3) * (-pkin(5) * t831 + qJ(6) * t830);
t526 = -pkin(4) * t710 + t629;
t443 = t461 * pkin(5) + t460 * qJ(6) - qJD(6) * t870 + t486;
t694 = -pkin(9) - t835;
t691 = t706 * pkin(7);
t678 = -pkin(4) * t846 - pkin(5);
t676 = pkin(4) * t700 + qJ(6);
t649 = pkin(4) * t685 + t845;
t604 = t685 * t703 + t686 * t818;
t602 = t685 * t706 - t686 * t819;
t592 = -pkin(7) * t828 + t632;
t581 = -pkin(7) * t778 + t624;
t574 = -t699 * t783 + t598;
t565 = t641 * t846 + t642 * t700;
t538 = t610 * t846 - t700 * t867;
t505 = pkin(5) * t565 - qJ(6) * t566 + t597;
t487 = pkin(5) * t538 - qJ(6) * t539 + t572;
t485 = qJD(5) * t539 + t700 * t546 - t710 * t846;
t484 = -t546 * t846 + t610 * t777 - t700 * t710 - t791 * t867;
t474 = -pkin(4) * t561 + t477;
t473 = t705 * pkin(5) - t748;
t472 = -qJ(6) * t705 + t859;
t455 = -t665 * qJ(6) + t457;
t454 = t665 * pkin(5) + t789;
t452 = t485 * pkin(5) + t484 * qJ(6) - t539 * qJD(6) + t526;
t445 = -pkin(5) * t796 - t850;
t444 = qJ(6) * t796 - qJD(6) * t705 + t734;
t442 = t766 - t854;
t441 = t767 + t856;
t1 = [(t456 * t796 + t572 * t461 + t525 * t485 + t486 * t538 - t516 * t526 + t748 * t627 - t665 * t850 + t766 * t705 + t762) * MDP(27) + (-t441 * t538 + t442 * t539 + t444 * t516 + t445 * t870 - t454 * t484 - t455 * t485 - t460 * t473 - t461 * t472 + t765) * MDP(30) + t760 * MDP(2) + t761 * MDP(3) + (t441 * t472 + t455 * t444 + t443 * t487 + t466 * t452 + t442 * t473 + t454 * t445 - g(1) * (-pkin(5) * t589 - qJ(6) * t588 + t649 * t706 + t691) - g(2) * (pkin(5) * t591 + qJ(6) * t590 + t646 * t818 - t694 * t821 + t802) + (-g(1) * (-t646 * t705 + t694 * t702 - pkin(1)) - g(2) * t649) * t703) * MDP(32) + (-t761 * t699 + (t612 * t699 + (-qJD(1) * t593 - t570) * qJD(2) + t736 * pkin(7)) * t702 + (t574 * qJD(1) + t593 * qJDD(1) + t522 - t760 * t698 + (t650 * t699 + (t635 + t778) * pkin(7)) * qJD(2)) * t705) * MDP(12) + (qJDD(1) * t696 + 0.2e1 * t702 * t775) * MDP(4) + (-t761 * t698 + (pkin(7) * t587 + t612 * t698 + (qJD(1) * t592 + t569) * qJD(2)) * t702 + (-t573 * qJD(1) - t592 * qJDD(1) - t521 + t760 * t699 + (t650 * t698 + t844) * qJD(2)) * t705) * MDP(11) + (t460 * t538 - t461 * t539 - t484 * t516 - t485 * t870) * MDP(23) + (t461 * t705 + t485 * t665 + t516 * t796 - t538 * t627) * MDP(25) + (t442 * t705 + t443 * t538 + t445 * t665 - t452 * t516 - t454 * t796 + t461 * t487 + t466 * t485 - t473 * t627 + t762) * MDP(29) + qJDD(1) * MDP(1) + (t750 * t702 + (-t730 + t841) * t705) * MDP(9) + (t522 * t593 + t570 * t574 + t521 * t592 + t569 * t573 - g(1) * t691 - g(2) * (t706 * t759 + t802) - t652 * t841 + (t612 * t702 + t650 * t795) * pkin(7)) * MDP(14) + (t702 * t730 + t705 * t750 - t673) * MDP(10) + (-t460 * t539 - t484 * t870) * MDP(22) + (t460 * t705 + t484 * t665 + t539 * t627 + t796 * t870) * MDP(24) + (-t441 * t705 - t443 * t539 - t444 * t665 - t452 * t870 + t455 * t796 + t460 * t487 + t466 * t484 + t472 * t627 - t763) * MDP(31) + (-t457 * t796 - t572 * t460 - t525 * t484 + t486 * t539 + t526 * t870 - t627 * t859 + t665 * t734 + t705 * t767 + t763) * MDP(28) + (-t771 * t672 + t769 * t640 - t773 * t705 - t629 * t756 + t648 * t753 + t541 * t610 - g(1) * t602 - g(2) * t604 + (t490 * t702 + t579 * t866) * qJD(2) + (t491 * t705 - t579 * t867 + t672 * t808) * qJD(4)) * MDP(20) + (qJDD(2) * t702 + t705 * t707) * MDP(6) + (qJDD(2) * t705 - t702 * t707) * MDP(7) + (-t509 * t610 + t546 * t756 - t561 * t710 + t753 * t867) * MDP(16) + ((-t575 * t793 + t782) * t672 - t808 * t640 + t735 * t705 - t491 * t796 - t629 * t561 + t648 * t509 - t541 * t867 + t579 * t546 - g(1) * t601 - g(2) * t603) * MDP(21) + (-t509 * t867 - t546 * t561) * MDP(15) + (-t509 * t705 - t546 * t672 - t561 * t796 - t640 * t867) * MDP(17) + (-t610 * t640 + t753 * t705 - t672 * t727 + (t672 * t866 + t702 * t756) * qJD(2)) * MDP(18) + 0.2e1 * (t687 * t702 - t788 * t801) * MDP(5) + (-t573 * t635 - t574 * t633 - t593 * t587 + (-qJDD(2) * t592 - t522 * t702 - t570 * t795) * t698 + (-t521 * t702 - t569 * t795 - t592 * t738) * t699 + t765) * MDP(13) + (-t640 * t705 - t672 * t796) * MDP(19) + (-t627 * t705 - t665 * t796) * MDP(26); (-MDP(4) * t702 * t705 + MDP(5) * t801) * t708 + (-t641 * t640 - t672 * t857) * MDP(18) + (t597 * t461 + t486 * t565 - t516 * t858 + t810 * t525 + t665 * t860 + t726) * MDP(27) + (t443 * t565 + t461 * t505 + t466 * t810 - t516 * t814 + t665 * t815 + t726) * MDP(29) + (t460 * t565 - t461 * t566 - t516 * t811 - t810 * t870) * MDP(23) + (-t441 * t565 + t442 * t566 - t454 * t811 - t455 * t810 + t460 * t745 - t461 * t507 + t516 * t816 + t815 * t870 + t723) * MDP(30) + (MDP(17) * t561 - MDP(18) * t756 + t672 * MDP(19) - t490 * MDP(20) + t491 * MDP(21) - MDP(24) * t870 - MDP(25) * t516 + t665 * MDP(26) - t456 * MDP(27) + t457 * MDP(28) + t454 * MDP(29) - t455 * MDP(31)) * t799 + (-t758 * t699 * qJDD(1) + (t612 + t722 + t693) * t698 + ((-qJ(3) * t790 + t570) * t702 + (-pkin(7) * t635 - t581 + (-t650 + t774) * t699) * t705) * qJD(1)) * MDP(12) + (t836 + (-pkin(7) * qJDD(1) + t731) * t705) * MDP(10) + (qJ(3) * t698 * t687 - pkin(2) * t587 + t718 * t699 + ((-qJ(3) * t797 - t569) * t702 + (-t844 + t580 + (qJD(3) - t650) * t698) * t705) * qJD(1)) * MDP(11) + qJDD(2) * MDP(8) + (-t650 * t681 - t569 * t580 - t570 * t581 + (-t569 * t698 + t570 * t699) * qJD(3) + t718 * pkin(2) + (-t521 * t698 + t522 * t699 + t723) * qJ(3)) * MDP(14) + (-t460 * t566 - t811 * t870) * MDP(22) + (t640 * t642 - t672 * t807) * MDP(17) + (t441 * t507 - t442 * t745 + t443 * t505 + t814 * t466 + t816 * t455 + t815 * t454 + (-g(3) * t752 + t694 * t761) * t705 + (g(3) * t694 + t752 * t761) * t702) * MDP(32) + (t509 * t642 - t561 * t807) * MDP(15) + (-t509 * t641 - t561 * t857 - t642 * t753 + t756 * t807) * MDP(16) + (t702 * t731 + t803) * MDP(9) + (t566 * t627 + t665 * t811) * MDP(24) + (-t565 * t627 + t665 * t810) * MDP(25) + (t580 * t635 + t581 * t633 + (qJ(3) * t786 + qJD(3) * t635 + t570 * t798 - t521) * t698 + (t569 * t798 - qJD(3) * t633 + t522 + (-t587 + t864) * qJ(3)) * t699 + t723) * MDP(13) + MDP(7) * t687 + MDP(6) * t787 + (-t597 * t460 + t486 * t566 - t811 * t525 + t665 * t861 + t858 * t870 - t881) * MDP(28) + (-t443 * t566 + t460 * t505 + t466 * t811 - t665 * t816 - t814 * t870 + t881) * MDP(31) + (-t805 * t640 - t781 * t509 + t541 * t642 + t628 * t561 + ((-qJD(4) * t657 - t794) * t701 + t882) * t672 + t807 * t579 - t725 * t685) * MDP(21) + (t541 * t641 - t857 * t579 + t628 * t756 + t768 * t640 + t672 * t883 + t725 * t686 - t781 * t753) * MDP(20); (t698 * t787 - t683 + (-t635 + t797) * t798) * MDP(11) + ((t633 + t790) * t798 + t736) * MDP(12) + (-t633 ^ 2 - t635 ^ 2) * MDP(13) + (t569 * t635 + t570 * t633 + t722 - t803 + t855) * MDP(14) + (t753 + t872) * MDP(20) + (t509 - t865) * MDP(21) + (-t847 - t885) * MDP(30) + (-t454 * t870 - t455 * t516 + t443 - t725) * MDP(32) + (-MDP(28) + MDP(31)) * (t460 + t884) + (MDP(27) + MDP(29)) * (t461 - t879); t561 * t756 * MDP(15) + (t561 ^ 2 - t756 ^ 2) * MDP(16) + (t509 + t865) * MDP(17) + (-t753 + t872) * MDP(18) + t640 * MDP(19) + (-t491 * t672 + t561 * t579 + t717 + t853) * MDP(20) + (g(1) * t604 - g(2) * t602 - t490 * t672 - t579 * t756 + t686 * t836 - t735) * MDP(21) + (-t458 * t665 + (-t516 * t561 + t627 * t846 + t665 * t791) * pkin(4) + t852) * MDP(27) + (-t459 * t665 - t873 + (t561 * t870 - t627 * t700 + t665 * t777) * pkin(4) + t720) * MDP(28) + (t474 * t516 - t627 * t678 + t665 * t764 - t712) * MDP(29) + (-t460 * t678 - t461 * t676 + (t455 + t764) * t870 - (t454 - t804) * t516) * MDP(30) + (t474 * t870 + t627 * t676 - t665 * t804 - t720 + t856 + t874) * MDP(31) + (t441 * t676 + t442 * t678 - t466 * t474 - t454 * t458 + t804 * t455 + (t454 * t791 + t853) * pkin(4) + t709) * MDP(32) + t876; (-t833 + t852) * MDP(27) + (t713 - t873) * MDP(28) + (t477 * t516 + t622 - t712 - t833) * MDP(29) + (pkin(5) * t460 - qJ(6) * t461 + (t455 - t457) * t870 - (t454 - t789) * t516) * MDP(30) + (t477 * t870 + 0.2e1 * t620 - 0.2e1 * t655 - t713 + t874) * MDP(31) + (-t442 * pkin(5) + t441 * qJ(6) - t454 * t457 + t455 * t789 - t466 * t477 + t709) * MDP(32) + t876; (-t627 - t880) * MDP(29) + t453 * MDP(30) + (-t665 ^ 2 - t847) * MDP(31) + (t455 * t665 + t712) * MDP(32);];
tau  = t1;
