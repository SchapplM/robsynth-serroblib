% Calculate vector of inverse dynamics joint torques for
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRPPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRRPPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:36:23
% EndTime: 2019-03-09 15:36:41
% DurationCPUTime: 14.24s
% Computational Cost: add. (7582->675), mult. (16921->855), div. (0->0), fcn. (12047->12), ass. (0->283)
t721 = sin(qJ(3));
t722 = sin(qJ(2));
t824 = qJD(1) * t722;
t797 = t721 * t824;
t725 = cos(qJ(3));
t810 = t725 * qJD(2);
t649 = t797 - t810;
t820 = qJD(2) * t721;
t651 = t725 * t824 + t820;
t717 = sin(pkin(10));
t718 = cos(pkin(10));
t583 = t718 * t649 + t651 * t717;
t720 = sin(qJ(6));
t724 = cos(qJ(6));
t758 = -t649 * t717 + t718 * t651;
t531 = t583 * t720 + t724 * t758;
t726 = cos(qJ(2));
t706 = t726 * qJDD(1);
t806 = qJD(1) * qJD(2);
t885 = -t722 * t806 + t706;
t641 = qJDD(3) - t885;
t638 = -qJDD(6) + t641;
t787 = t726 * t806;
t805 = qJDD(1) * t722;
t816 = qJD(3) * t722;
t894 = -qJD(1) * t816 + qJDD(2);
t574 = qJD(3) * t810 + (t787 + t805) * t725 + t894 * t721;
t823 = qJD(1) * t726;
t575 = t721 * (qJD(2) * (qJD(3) + t823) + t805) - t894 * t725;
t515 = t574 * t717 + t718 * t575;
t516 = t574 * t718 - t575 * t717;
t812 = qJD(6) * t724;
t813 = qJD(6) * t720;
t778 = -t720 * t515 - t724 * t516 - t583 * t812 + t758 * t813;
t683 = -qJD(3) + t823;
t807 = -qJD(6) - t683;
t893 = -t724 * t583 + t720 * t758;
t864 = t893 * t807;
t907 = MDP(24) * t531 * t893 + (t531 ^ 2 - t893 ^ 2) * MDP(25) - (t778 + t864) * MDP(26) - t638 * MDP(28);
t865 = t531 * t807;
t663 = -qJD(2) * pkin(2) + pkin(7) * t824;
t741 = -pkin(3) * t649 - qJD(4) - t663;
t733 = qJ(5) * t758 + t741;
t876 = pkin(4) + pkin(5);
t493 = -t583 * t876 + t733;
t714 = qJ(3) + pkin(10);
t703 = cos(t714);
t723 = sin(qJ(1));
t702 = sin(t714);
t727 = cos(qJ(1));
t842 = t727 * t702;
t613 = -t723 * t703 + t726 * t842;
t843 = t726 * t727;
t614 = t702 * t723 + t703 * t843;
t554 = t613 * t724 - t614 * t720;
t754 = t702 * t724 - t703 * t720;
t776 = pkin(2) * t722 - pkin(8) * t726;
t653 = t776 * qJD(2);
t660 = -pkin(2) * t726 - pkin(8) * t722 - pkin(1);
t595 = qJD(1) * t653 + qJDD(1) * t660;
t589 = t725 * t595;
t639 = t660 * qJD(1);
t700 = pkin(7) * t823;
t664 = qJD(2) * pkin(8) + t700;
t594 = t639 * t721 + t664 * t725;
t622 = pkin(7) * t885 + qJDD(2) * pkin(8);
t484 = pkin(3) * t641 - qJ(4) * t574 - qJD(3) * t594 - qJD(4) * t651 - t622 * t721 + t589;
t815 = qJD(3) * t725;
t817 = qJD(3) * t721;
t743 = t721 * t595 + t725 * t622 + t639 * t815 - t664 * t817;
t488 = -qJ(4) * t575 - qJD(4) * t649 + t743;
t468 = t718 * t484 - t717 * t488;
t785 = -qJDD(5) + t468;
t460 = -pkin(9) * t516 - t641 * t876 - t785;
t667 = t683 * qJD(5);
t469 = t717 * t484 + t718 * t488;
t799 = t641 * qJ(5) + t469;
t464 = -t667 + t799;
t462 = pkin(9) * t515 + t464;
t783 = t724 * t460 - t720 * t462;
t867 = g(3) * t722;
t846 = t723 * t726;
t611 = t702 * t846 + t703 * t727;
t612 = t703 * t846 - t842;
t884 = t611 * t724 - t612 * t720;
t904 = g(1) * t554 + g(2) * t884 + t493 * t531 + t754 * t867 - t783;
t789 = t717 * t817;
t796 = t721 * t823;
t856 = t718 * t725;
t833 = -t717 * t796 - t718 * t815 + t823 * t856 + t789;
t555 = t613 * t720 + t614 * t724;
t753 = t702 * t720 + t703 * t724;
t761 = t611 * t720 + t612 * t724;
t593 = t725 * t639 - t664 * t721;
t556 = -qJ(4) * t651 + t593;
t544 = -pkin(3) * t683 + t556;
t557 = -qJ(4) * t649 + t594;
t858 = t717 * t557;
t500 = t544 * t718 - t858;
t770 = qJD(5) - t500;
t891 = pkin(9) * t758;
t479 = t683 * t876 + t770 - t891;
t800 = t720 * t460 + t724 * t462 + t479 * t812;
t901 = -g(1) * t555 - g(2) * t761 - t493 * t893 - t753 * t867 + t800;
t899 = pkin(9) * t583;
t898 = t583 * t683;
t643 = t717 * t725 + t718 * t721;
t627 = t643 * qJD(3);
t834 = t643 * t823 - t627;
t713 = g(3) * t726;
t775 = g(1) * t727 + g(2) * t723;
t889 = t722 * t775;
t897 = -t713 + t889;
t818 = qJD(2) * t726;
t794 = t721 * t818;
t896 = t722 * t815 + t794;
t895 = -t700 + (-t796 + t817) * pkin(3);
t892 = t758 ^ 2;
t845 = t725 * t726;
t750 = pkin(3) * t722 - qJ(4) * t845;
t652 = t776 * qJD(1);
t830 = pkin(7) * t797 + t725 * t652;
t569 = qJD(1) * t750 + t830;
t633 = t721 * t652;
t849 = t722 * t725;
t852 = t721 * t726;
t587 = t633 + (-pkin(7) * t849 - qJ(4) * t852) * qJD(1);
t719 = -qJ(4) - pkin(8);
t784 = qJD(3) * t719;
t814 = qJD(4) * t725;
t620 = t721 * t784 + t814;
t621 = -qJD(4) * t721 + t725 * t784;
t837 = (-t569 + t621) * t718 + (t587 - t620) * t717;
t519 = t717 * t569 + t718 * t587;
t509 = qJ(5) * t824 + t519;
t564 = t718 * t620 + t717 * t621;
t835 = t564 - t509;
t886 = qJ(5) * t833 - qJD(5) * t643 + t895;
t504 = t718 * t556 - t858;
t808 = qJD(5) - t504;
t698 = pkin(7) * t805;
t774 = qJDD(2) * pkin(2) - pkin(7) * t787 - t698;
t883 = qJD(3) * pkin(8) * t683 + t774 + t897;
t506 = pkin(4) * t583 - t733;
t882 = g(1) * t613 + g(2) * t611 - t506 * t758 + t702 * t867 + t785;
t881 = t726 * t775 + t867;
t782 = -t724 * t515 + t516 * t720;
t472 = qJD(6) * t531 + t782;
t877 = -0.2e1 * pkin(1);
t875 = pkin(4) * t641;
t874 = pkin(7) * t721;
t871 = g(1) * t723;
t868 = g(2) * t727;
t857 = t718 * t557;
t503 = t556 * t717 + t857;
t866 = t503 * t758;
t863 = t574 * t721;
t861 = t649 * t683;
t860 = t651 * t683;
t859 = t651 * t725;
t501 = t717 * t544 + t857;
t497 = -t683 * qJ(5) + t501;
t485 = t497 + t899;
t855 = t720 * t485;
t854 = t721 * t722;
t853 = t721 * t723;
t851 = t721 * t727;
t850 = t722 * t723;
t848 = t722 * t727;
t847 = t723 * t725;
t844 = t725 * t727;
t697 = pkin(3) * t725 + pkin(2);
t665 = t726 * t697;
t642 = t717 * t721 - t856;
t759 = t724 * t642 - t643 * t720;
t841 = qJD(6) * t759 - t720 * t834 - t724 * t833;
t581 = t642 * t720 + t643 * t724;
t840 = qJD(6) * t581 - t720 * t833 + t724 * t834;
t688 = pkin(7) * t845;
t819 = qJD(2) * t722;
t831 = t725 * t653 + t819 * t874;
t520 = -t722 * t814 + t750 * qJD(2) + (-t688 + (qJ(4) * t722 - t660) * t721) * qJD(3) + t831;
t832 = t721 * t653 + t660 * t815;
t535 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t849 + (-qJD(4) * t722 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t726) * t721 + t832;
t483 = t717 * t520 + t718 * t535;
t839 = t834 * t876 - t886;
t838 = -pkin(4) * t834 + t886;
t836 = pkin(4) * t824 - t837;
t645 = t725 * t660;
t590 = -qJ(4) * t849 + t645 + (-pkin(3) - t874) * t726;
t828 = t721 * t660 + t688;
t596 = -qJ(4) * t854 + t828;
t537 = t717 * t590 + t718 * t596;
t661 = t719 * t721;
t662 = t719 * t725;
t598 = t717 * t661 - t718 * t662;
t685 = pkin(3) * t854;
t827 = t722 * pkin(7) + t685;
t715 = t722 ^ 2;
t826 = -t726 ^ 2 + t715;
t822 = qJD(2) * t649;
t821 = qJD(2) * t651;
t809 = -t808 + t891;
t802 = t876 * t722;
t801 = t721 * t843;
t798 = pkin(3) * t896 + pkin(7) * t818;
t693 = -pkin(3) * t718 - pkin(4);
t795 = t683 * t810;
t793 = t726 * t810;
t792 = t683 * t817;
t791 = t683 * t815;
t482 = t520 * t718 - t717 * t535;
t536 = t590 * t718 - t717 * t596;
t597 = -t718 * t661 - t662 * t717;
t780 = -qJD(3) * t639 - t622;
t690 = g(1) * t850;
t777 = -g(2) * t848 + t690;
t532 = -qJ(5) * t726 + t537;
t618 = -t717 * t854 + t718 * t849;
t772 = qJ(5) * t618 - t827;
t534 = t726 * pkin(4) - t536;
t771 = t664 * t815 - t589;
t567 = pkin(9) * t642 + t598;
t769 = -pkin(9) * t833 - qJD(1) * t802 + qJD(6) * t567 + t837;
t566 = -pkin(9) * t643 + t597;
t768 = pkin(9) * t834 - qJD(6) * t566 - t835;
t767 = -pkin(3) * t651 - qJ(5) * t583;
t766 = pkin(4) * t703 + qJ(5) * t702;
t765 = -pkin(8) * t641 + qJD(3) * t663;
t466 = t720 * t479 + t724 * t485;
t505 = pkin(5) * t726 - pkin(9) * t618 + t534;
t617 = t643 * t722;
t507 = pkin(9) * t617 + t532;
t764 = t505 * t724 - t507 * t720;
t763 = t505 * t720 + t507 * t724;
t760 = t724 * t617 - t618 * t720;
t560 = t617 * t720 + t618 * t724;
t681 = -pkin(5) + t693;
t691 = pkin(3) * t717 + qJ(5);
t757 = t681 * t724 - t691 * t720;
t756 = t681 * t720 + t691 * t724;
t751 = qJ(5) * t643 + t697;
t476 = qJ(5) * t819 - qJD(5) * t726 + t483;
t746 = -pkin(7) * qJDD(2) + t806 * t877;
t629 = t721 * t846 + t844;
t745 = t641 * t721 - t791;
t744 = t641 * t725 + t792;
t729 = qJD(1) ^ 2;
t742 = pkin(1) * t729 + t775;
t728 = qJD(2) ^ 2;
t740 = pkin(7) * t728 + qJDD(1) * t877 + t868;
t739 = t727 * pkin(1) + pkin(3) * t853 + t723 * pkin(7) + t697 * t843 - t719 * t848;
t738 = -pkin(3) * t575 - qJDD(4) + t774;
t562 = t627 * t722 + t717 * t794 - t718 * t793;
t737 = -qJ(5) * t562 + qJD(5) * t618 - t798;
t736 = t719 * t850 + pkin(3) * t851 + t727 * pkin(7) + (-pkin(1) - t665) * t723;
t731 = qJ(5) * t516 + qJD(5) * t758 + t738;
t730 = -t598 * t515 + t516 * t597 - t564 * t583 - t881;
t470 = pkin(4) * t515 - t731;
t687 = pkin(3) * t847;
t632 = t725 * t843 + t853;
t631 = -t801 + t847;
t630 = -t723 * t845 + t851;
t578 = pkin(4) * t642 - t751;
t561 = -t717 * t793 - t718 * t896 + t722 * t789;
t546 = -t642 * t876 + t751;
t545 = pkin(4) * t617 - t772;
t533 = -t617 * t876 + t772;
t508 = pkin(4) * t758 - t767;
t496 = pkin(4) * t683 + t770;
t495 = -t758 * t876 + t767;
t494 = -pkin(4) * t561 - t737;
t491 = t503 + t899;
t490 = qJD(6) * t560 + t724 * t561 - t562 * t720;
t489 = qJD(6) * t760 - t561 * t720 - t562 * t724;
t478 = -pkin(4) * t819 - t482;
t475 = t561 * t876 + t737;
t474 = -pkin(9) * t561 + t476;
t473 = pkin(9) * t562 - qJD(2) * t802 - t482;
t467 = -t785 - t875;
t465 = t479 * t724 - t855;
t463 = -t515 * t876 + t731;
t1 = [t775 * MDP(3) + (t489 * t531 - t560 * t778) * MDP(24) + (-g(1) * t736 - g(2) * t739 + t468 * t536 + t469 * t537 + t500 * t482 + t501 * t483 - t738 * t827 - t741 * t798) * MDP(19) + (-(-t660 * t817 + t831) * t683 + t645 * t641 - g(1) * t630 - g(2) * t632 + ((t791 + t822) * pkin(7) + (-pkin(7) * t641 + qJD(2) * t663 - t780) * t721 + t771) * t726 + (pkin(7) * t575 + qJD(2) * t593 + t663 * t815 - t721 * t774) * t722) * MDP(16) + (t832 * t683 - t828 * t641 - g(1) * t629 - g(2) * t631 + (t663 * t810 + (-t792 + t821) * pkin(7) + t743) * t726 + (-t663 * t817 - t594 * qJD(2) - t774 * t725 + (t574 - t795) * pkin(7)) * t722) * MDP(17) + (-t489 * t807 - t531 * t819 - t560 * t638 - t726 * t778) * MDP(26) + (-t638 * t726 + t807 * t819) * MDP(28) + ((-t574 - t795) * t726 + (t744 + t821) * t722) * MDP(13) + ((t683 * t820 + t575) * t726 + (-t745 - t822) * t722) * MDP(14) + (g(1) * t611 - g(2) * t613 - t464 * t726 - t470 * t618 - t476 * t683 - t494 * t758 + t497 * t819 + t506 * t562 - t516 * t545 + t532 * t641) * MDP(22) + (-t468 * t618 - t469 * t617 - t482 * t758 - t483 * t583 + t500 * t562 + t501 * t561 - t515 * t537 - t516 * t536 + t777) * MDP(18) + (-t464 * t617 + t467 * t618 - t476 * t583 + t478 * t758 - t496 * t562 + t497 * t561 - t515 * t532 + t516 * t534 + t777) * MDP(21) + (-t472 * t560 - t489 * t893 - t490 * t531 - t760 * t778) * MDP(25) + (-(t473 * t724 - t474 * t720) * t807 - t764 * t638 + t783 * t726 - t465 * t819 + t475 * t893 + t533 * t472 - t463 * t760 + t493 * t490 + g(1) * t761 - g(2) * t555 + (-t466 * t726 + t763 * t807) * qJD(6)) * MDP(29) + (-t472 * t726 + t490 * t807 - t638 * t760 + t819 * t893) * MDP(27) + (qJDD(1) * t715 + 0.2e1 * t722 * t787) * MDP(4) + (-t868 + t871) * MDP(2) + (t746 * t722 + (-t740 + t871) * t726) * MDP(9) + (t574 * t849 + (-t721 * t816 + t793) * t651) * MDP(11) + ((qJD(6) * t764 + t473 * t720 + t474 * t724) * t807 + t763 * t638 - (-t485 * t813 + t800) * t726 + t466 * t819 + t475 * t531 - t533 * t778 + t463 * t560 + t493 * t489 + g(1) * t884 - g(2) * t554) * MDP(30) + (t722 * t740 + t726 * t746 - t690) * MDP(10) + (-t641 * t726 - t683 * t819) * MDP(15) + (g(1) * t612 - g(2) * t614 + t467 * t726 + t470 * t617 + t478 * t683 + t494 * t583 - t496 * t819 - t506 * t561 + t515 * t545 - t534 * t641) * MDP(20) + (t464 * t532 + t497 * t476 + t470 * t545 + t506 * t494 + t467 * t534 + t496 * t478 - g(1) * (-pkin(4) * t612 - qJ(5) * t611 + t736) - g(2) * (pkin(4) * t614 + qJ(5) * t613 + t739)) * MDP(23) + ((-t649 * t725 - t651 * t721) * t818 + (-t863 - t575 * t725 + (t649 * t721 - t859) * qJD(3)) * t722) * MDP(12) + qJDD(1) * MDP(1) + 0.2e1 * (t706 * t722 - t806 * t826) * MDP(5) + (qJDD(2) * t726 - t722 * t728) * MDP(7) + (qJDD(2) * t722 + t726 * t728) * MDP(6); ((t574 + t861) * t725 + (-t575 + t860) * t721) * MDP(12) + (t722 * t742 - t698 - t713) * MDP(9) + (t531 * t841 - t581 * t778) * MDP(24) + (-t638 * t759 + t807 * t840) * MDP(27) + (-t581 * t638 - t807 * t841) * MDP(26) + (-MDP(4) * t722 * t726 + MDP(5) * t826) * t729 + (-t464 * t642 + t467 * t643 - t496 * t833 + t497 * t834 + t509 * t583 + t758 * t836 + t730) * MDP(21) + (-t468 * t643 - t469 * t642 + t500 * t833 + t501 * t834 + t519 * t583 - t758 * t837 + t730) * MDP(18) + (-t472 * t581 - t531 * t840 - t759 * t778 - t841 * t893) * MDP(25) + (t683 * MDP(15) + t496 * MDP(20) - t497 * MDP(22) + t531 * MDP(26) - MDP(27) * t893 - MDP(28) * t807 + t465 * MDP(29) - t466 * MDP(30)) * t824 + (-t470 * t643 + t506 * t833 - t516 * t578 + t598 * t641 - t683 * t835 + t702 * t897 - t758 * t838) * MDP(22) + (t470 * t642 - t506 * t834 + t515 * t578 + t583 * t838 - t597 * t641 + t683 * t836 + t703 * t897) * MDP(20) + ((-t651 * t722 + t683 * t845) * qJD(1) + t745) * MDP(13) + (t867 + (-pkin(7) * qJDD(1) + t742) * t726) * MDP(10) + ((t649 * t722 - t683 * t852) * qJD(1) + t744) * MDP(14) + ((t566 * t720 + t567 * t724) * t638 - t546 * t778 + t463 * t581 - (t720 * t769 + t724 * t768) * t807 + t839 * t531 + t841 * t493 + t897 * t754) * MDP(30) + (-(t566 * t724 - t567 * t720) * t638 + t546 * t472 - t463 * t759 - (t720 * t768 - t724 * t769) * t807 + t839 * t893 + t840 * t493 + t897 * t753) * MDP(29) + (-g(3) * t665 + t464 * t598 + t467 * t597 + t470 * t578 + t838 * t506 + t835 * t497 + t836 * t496 + (-g(3) * t766 + t719 * t775) * t726 + (g(3) * t719 + t775 * (t697 + t766)) * t722) * MDP(23) + (-pkin(2) * t575 + t830 * t683 + t765 * t721 + (-t593 * t722 + (-pkin(7) * t649 - t663 * t721) * t726) * qJD(1) + t883 * t725) * MDP(16) + (-pkin(2) * t574 - t633 * t683 + t765 * t725 + (-t663 * t845 + t594 * t722 + (-t651 * t726 + t683 * t849) * pkin(7)) * qJD(1) - t883 * t721) * MDP(17) + qJDD(2) * MDP(8) + (-t683 * t859 + t863) * MDP(11) + (t469 * t598 - t468 * t597 + t738 * t697 - g(3) * (-t719 * t722 + t665) - t895 * t741 + (t564 - t519) * t501 + t837 * t500 + t775 * (t697 * t722 + t719 * t726)) * MDP(19) + MDP(7) * t706 + MDP(6) * t805; t651 * t649 * MDP(11) + (-t649 ^ 2 + t651 ^ 2) * MDP(12) + (t574 - t861) * MDP(13) + (-t575 - t860) * MDP(14) + t641 * MDP(15) + (-g(1) * t631 + g(2) * t629 - t594 * t683 - t651 * t663 + (t780 + t867) * t721 - t771) * MDP(16) + (g(1) * t632 - g(2) * t630 + g(3) * t849 - t593 * t683 + t649 * t663 - t743) * MDP(17) + (t501 * t758 - t866 + (-t515 * t717 - t516 * t718) * pkin(3) + (-t500 + t504) * t583) * MDP(18) + (-g(1) * t687 + t500 * t503 - t501 * t504 + (g(2) * t844 + t468 * t718 + t469 * t717 + t651 * t741 + t721 * t881) * pkin(3)) * MDP(19) + (-t503 * t683 - t508 * t583 + (pkin(4) - t693) * t641 + t882) * MDP(20) + (t497 * t758 - t515 * t691 + t516 * t693 - t866 + (t496 - t808) * t583) * MDP(21) + (-g(1) * t614 - g(2) * t612 + t504 * t683 - t506 * t583 + t508 * t758 + t641 * t691 - t703 * t867 - 0.2e1 * t667 + t799) * MDP(22) + (t464 * t691 + t467 * t693 - t506 * t508 - t496 * t503 - g(1) * (-pkin(3) * t801 - pkin(4) * t613 + qJ(5) * t614 + t687) - g(2) * (-pkin(3) * t629 - pkin(4) * t611 + qJ(5) * t612) - g(3) * (-t685 + (-pkin(4) * t702 + qJ(5) * t703) * t722) + t808 * t497) * MDP(23) + (t472 + t865) * MDP(27) + (-t757 * t638 - t495 * t893 - (-t724 * t491 + t720 * t809) * t807 + (t756 * t807 + t466) * qJD(6) + t904) * MDP(29) + (t756 * t638 - t495 * t531 - (t720 * t491 + t724 * t809) * t807 + (t757 * t807 - t855) * qJD(6) + t901) * MDP(30) - t907; (t500 * t758 + t501 * t583 + t713 - t738) * MDP(19) + (-t683 * t758 + t515) * MDP(20) + (-t516 - t898) * MDP(22) + (-t496 * t758 + t497 * t583 + t470 + t713) * MDP(23) + (-t472 + t865) * MDP(29) + (t778 - t864) * MDP(30) + (-MDP(19) - MDP(23)) * t889 + (MDP(18) + MDP(21)) * (-t583 ^ 2 - t892); (t583 * t758 - t641) * MDP(20) + (t516 - t898) * MDP(21) + (-t683 ^ 2 - t892) * MDP(22) + (t497 * t683 - t875 - t882) * MDP(23) + (-t638 * t724 - t758 * t893) * MDP(29) + (-t531 * t758 + t638 * t720) * MDP(30) - (MDP(29) * t720 + MDP(30) * t724) * t807 ^ 2; (-t782 - t865) * MDP(27) + (-t466 * t807 - t904) * MDP(29) + (-t465 * t807 - t901) * MDP(30) + (-MDP(27) * t531 - MDP(29) * t466 + MDP(30) * t855) * qJD(6) + t907;];
tau  = t1;
