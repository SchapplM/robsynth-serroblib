% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:30:10
% EndTime: 2019-03-10 03:30:29
% DurationCPUTime: 14.12s
% Computational Cost: add. (17332->574), mult. (41984->739), div. (0->0), fcn. (33222->18), ass. (0->271)
t752 = sin(qJ(3));
t753 = sin(qJ(2));
t906 = cos(qJ(3));
t838 = qJDD(1) * t906;
t758 = cos(qJ(2));
t850 = qJDD(1) * t758;
t745 = qJD(2) + qJD(3);
t843 = qJD(1) * t906;
t862 = qJD(1) * t753;
t941 = -t752 * t862 + t758 * t843;
t955 = t745 * t941;
t607 = t752 * t850 + t753 * t838 + t955;
t882 = t752 * t758;
t678 = t753 * t906 + t882;
t640 = t745 * t678;
t851 = qJDD(1) * t753;
t807 = t752 * t851 - t758 * t838;
t608 = qJD(1) * t640 + t807;
t666 = -qJD(1) * t882 - t753 * t843;
t751 = sin(qJ(4));
t757 = cos(qJ(4));
t859 = qJD(4) * t757;
t860 = qJD(4) * t751;
t553 = t757 * t607 - t751 * t608 + t666 * t860 + t859 * t941;
t750 = sin(qJ(5));
t756 = cos(qJ(5));
t798 = t757 * t666 - t751 * t941;
t778 = qJD(4) * t798 - t607 * t751 - t757 * t608;
t799 = t751 * t666 + t757 * t941;
t857 = qJD(5) * t756;
t858 = qJD(5) * t750;
t510 = t756 * t553 + t750 * t778 + t798 * t858 + t799 * t857;
t744 = qJDD(2) + qJDD(3);
t739 = qJDD(4) + t744;
t722 = qJDD(5) + t739;
t740 = qJD(4) + t745;
t725 = qJD(5) + t740;
t749 = sin(qJ(6));
t755 = cos(qJ(6));
t855 = qJD(6) * t755;
t847 = t755 * t510 + t749 * t722 + t725 * t855;
t856 = qJD(6) * t749;
t927 = t750 * t799 - t756 * t798;
t502 = -t856 * t927 + t847;
t500 = t502 * t749;
t501 = t502 * t755;
t569 = t725 * t749 + t755 * t927;
t693 = t755 * t722;
t503 = qJD(6) * t569 + t510 * t749 - t693;
t511 = qJD(5) * t927 + t553 * t750 - t756 * t778;
t509 = qJDD(6) + t511;
t506 = t749 * t509;
t507 = t755 * t509;
t567 = -t755 * t725 + t749 * t927;
t582 = t750 * t798 + t756 * t799;
t929 = qJD(6) - t582;
t948 = t929 * t749;
t957 = -t582 * t755 + t855;
t958 = t722 * MDP(29) + (t569 * t957 + t500) * MDP(32) + (-t569 * t927 + t929 * t957 + t506) * MDP(34) - t511 * MDP(28) - t582 ^ 2 * MDP(26) + (-t582 * t725 + t510) * MDP(27) + (-MDP(25) * t582 + MDP(26) * t927 + MDP(28) * t725 - MDP(36) * t929) * t927 + (-t749 * t503 - t567 * t957 - t569 * t948 + t501) * MDP(33) + (t567 * t927 - t929 * t948 + t507) * MDP(35);
t956 = t739 * MDP(22) + t799 * MDP(18) * t798 + (-t740 * t799 + t553) * MDP(20) + (t798 ^ 2 - t799 ^ 2) * MDP(19) + (-t740 * t798 + t778) * MDP(21) + t958;
t624 = t798 * pkin(10);
t659 = t666 * pkin(9);
t907 = pkin(7) + pkin(8);
t701 = t907 * t758;
t684 = qJD(1) * t701;
t667 = t752 * t684;
t700 = t907 * t753;
t682 = qJD(1) * t700;
t899 = qJD(2) * pkin(2);
t673 = -t682 + t899;
t827 = t906 * t673 - t667;
t605 = t659 + t827;
t593 = pkin(3) * t745 + t605;
t671 = t906 * t684;
t794 = -t752 * t673 - t671;
t900 = t941 * pkin(9);
t606 = -t794 + t900;
t595 = t751 * t606;
t832 = t757 * t593 - t595;
t545 = t832 + t624;
t540 = pkin(4) * t740 + t545;
t597 = t757 * t606;
t802 = -t593 * t751 - t597;
t904 = pkin(10) * t799;
t546 = -t802 + t904;
t894 = t546 * t756;
t520 = t540 * t750 + t894;
t518 = pkin(11) * t725 + t520;
t735 = -pkin(2) * t758 - pkin(1);
t699 = t735 * qJD(1);
t643 = -pkin(3) * t941 + t699;
t589 = -pkin(4) * t799 + t643;
t531 = -pkin(5) * t582 - pkin(11) * t927 + t589;
t497 = -t518 * t749 + t531 * t755;
t895 = t546 * t750;
t519 = t540 * t756 - t895;
t517 = -pkin(5) * t725 - t519;
t748 = qJ(2) + qJ(3);
t743 = qJ(4) + t748;
t736 = qJ(5) + t743;
t720 = sin(t736);
t754 = sin(qJ(1));
t759 = cos(qJ(1));
t809 = g(1) * t759 + g(2) * t754;
t949 = t809 * t720;
t954 = -t497 * t927 + t517 * t856 + t755 * t949;
t498 = t518 * t755 + t531 * t749;
t852 = qJD(1) * qJD(2);
t840 = t758 * t852;
t641 = qJDD(2) * pkin(2) + t907 * (-t840 - t851);
t841 = t753 * t852;
t642 = t907 * (-t841 + t850);
t776 = qJD(3) * t794 + t906 * t641 - t752 * t642;
t543 = t744 * pkin(3) - t607 * pkin(9) + t776;
t842 = qJD(3) * t906;
t861 = qJD(3) * t752;
t774 = t752 * t641 + t642 * t906 + t673 * t842 - t684 * t861;
t550 = -t608 * pkin(9) + t774;
t781 = qJD(4) * t802 + t757 * t543 - t751 * t550;
t495 = pkin(4) * t739 - pkin(10) * t553 + t781;
t912 = (qJD(4) * t593 + t550) * t757 + t751 * t543 - t606 * t860;
t496 = pkin(10) * t778 + t912;
t911 = qJD(5) * t520 - t756 * t495 + t750 * t496;
t485 = -pkin(5) * t722 + t911;
t721 = cos(t736);
t901 = g(3) * t721;
t933 = t485 + t901;
t953 = t498 * t927 + t517 * t855 + t749 * t933;
t951 = t517 * t582;
t771 = -t749 * t949 + t953;
t777 = -t755 * t933 + t954;
t922 = pkin(5) * t927;
t947 = -pkin(11) * t582 + t922;
t711 = g(3) * t720;
t913 = (qJD(5) * t540 + t496) * t756 + t750 * t495 - t546 * t858;
t946 = -t589 * t582 + t721 * t809 + t711 - t913;
t826 = t682 * t752 - t671;
t610 = t826 - t900;
t869 = -t906 * t682 - t667;
t611 = t659 + t869;
t734 = pkin(2) * t906 + pkin(3);
t885 = t751 * t752;
t943 = -t734 * t859 - (-t752 * t860 + (t757 * t906 - t885) * qJD(3)) * pkin(2) + t751 * t610 + t757 * t611;
t883 = t752 * t757;
t942 = -t734 * t860 + (-t752 * t859 + (-t751 * t906 - t883) * qJD(3)) * pkin(2) - t757 * t610 + t611 * t751;
t939 = -t589 * t927 - t901 - t911 + t949;
t936 = pkin(4) * t798;
t932 = -t904 - t942;
t931 = t624 + t943;
t729 = sin(t743);
t730 = cos(t743);
t764 = -g(3) * t730 + t643 * t798 + t729 * t809 + t781;
t919 = t929 * (pkin(11) * t929 + t922);
t825 = -t906 * t700 - t701 * t752;
t622 = -pkin(9) * t678 + t825;
t793 = -t752 * t753 + t758 * t906;
t868 = -t752 * t700 + t906 * t701;
t623 = pkin(9) * t793 + t868;
t873 = t751 * t622 + t757 * t623;
t812 = -pkin(2) * t885 + t757 * t734;
t657 = pkin(4) + t812;
t661 = pkin(2) * t883 + t734 * t751;
t870 = t750 * t657 + t756 * t661;
t769 = g(3) * t729 - t643 * t799 + t730 * t809 - t912;
t638 = t678 * t757 + t751 * t793;
t639 = t745 * t793;
t564 = qJD(4) * t638 + t639 * t751 + t757 * t640;
t845 = qJD(2) * t907;
t683 = t753 * t845;
t685 = t758 * t845;
t789 = -t906 * t683 - t752 * t685 - t700 * t842 - t701 * t861;
t573 = -pkin(9) * t640 + t789;
t775 = -qJD(3) * t868 + t752 * t683 - t906 * t685;
t574 = -t639 * pkin(9) + t775;
t788 = t757 * t573 + t751 * t574 + t622 * t859 - t623 * t860;
t512 = -pkin(10) * t564 + t788;
t637 = t678 * t751 - t757 * t793;
t563 = -qJD(4) * t637 + t639 * t757 - t640 * t751;
t780 = -qJD(4) * t873 - t573 * t751 + t757 * t574;
t513 = -pkin(10) * t563 + t780;
t828 = t757 * t622 - t623 * t751;
t560 = -pkin(10) * t638 + t828;
t561 = -pkin(10) * t637 + t873;
t803 = t560 * t756 - t561 * t750;
t486 = qJD(5) * t803 + t512 * t756 + t513 * t750;
t584 = t756 * t637 + t638 * t750;
t527 = -qJD(5) * t584 + t563 * t756 - t564 * t750;
t530 = t560 * t750 + t561 * t756;
t585 = -t637 * t750 + t638 * t756;
t647 = -pkin(3) * t793 + t735;
t598 = pkin(4) * t637 + t647;
t536 = pkin(5) * t584 - pkin(11) * t585 + t598;
t484 = pkin(11) * t722 + t913;
t819 = qJD(6) * t531 + t484;
t910 = t485 * t585 - t530 * t509 + t517 * t527 - t929 * (qJD(6) * t536 + t486) - t584 * t819;
t905 = pkin(3) * t666;
t897 = t517 * t585;
t896 = t536 * t509;
t742 = cos(t748);
t890 = t742 * t754;
t889 = t742 * t759;
t888 = t749 * t754;
t887 = t749 * t759;
t886 = t750 * t751;
t884 = t751 * t756;
t881 = t754 * t755;
t880 = t755 * t759;
t800 = t657 * t756 - t661 * t750;
t877 = -qJD(5) * t800 + t750 * t932 + t756 * t931;
t876 = t870 * qJD(5) - t750 * t931 + t756 * t932;
t875 = t757 * t605 - t595;
t831 = -t605 * t751 - t597;
t551 = t831 - t904;
t552 = t624 + t875;
t733 = pkin(3) * t757 + pkin(4);
t872 = t551 * t750 + t552 * t756 - t733 * t857 - (-t751 * t858 + (t756 * t757 - t886) * qJD(4)) * pkin(3);
t871 = t551 * t756 - t552 * t750 + t733 * t858 + (t751 * t857 + (t750 * t757 + t884) * qJD(4)) * pkin(3);
t866 = pkin(3) * t884 + t750 * t733;
t746 = t753 ^ 2;
t865 = -t758 ^ 2 + t746;
t738 = t753 * t899;
t629 = pkin(3) * t640 + t738;
t660 = pkin(2) * t841 + qJDD(1) * t735;
t588 = pkin(3) * t608 + t660;
t534 = -pkin(4) * t778 + t588;
t489 = pkin(5) * t511 - pkin(11) * t510 + t534;
t817 = qJD(6) * t518 - t489;
t592 = -t905 - t936;
t533 = t592 + t947;
t618 = pkin(11) + t870;
t737 = pkin(2) * t862;
t816 = qJD(6) * t618 + t533 + t737;
t658 = pkin(11) + t866;
t815 = qJD(6) * t658 + t533;
t731 = pkin(4) * t750 + pkin(11);
t814 = qJD(6) * t731 - t936 + t947;
t521 = t545 * t750 + t894;
t811 = pkin(4) * t858 - t521;
t522 = t545 * t756 - t895;
t810 = -pkin(4) * t857 + t522;
t808 = g(1) * t754 - g(2) * t759;
t806 = -t509 * t731 - t951;
t805 = -t509 * t618 - t951;
t804 = -t509 * t658 - t951;
t557 = pkin(4) * t564 + t629;
t797 = -pkin(3) * t886 + t733 * t756;
t795 = -0.2e1 * pkin(1) * t852 - pkin(7) * qJDD(2);
t792 = t527 * t755 - t585 * t856;
t787 = -pkin(11) * t509 + t519 * t929 - t951;
t760 = qJD(2) ^ 2;
t784 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t760 + t808;
t761 = qJD(1) ^ 2;
t783 = pkin(1) * t761 - pkin(7) * qJDD(1) + t809;
t741 = sin(t748);
t766 = g(1) * t889 + g(2) * t890 + g(3) * t741 - t699 * t941 - t774;
t763 = -g(3) * t742 + t699 * t666 + t741 * t809 + t776;
t762 = t666 * t941 * MDP(11) + (t607 - t955) * MDP(13) + (-t807 + (-qJD(1) * t678 - t666) * t745) * MDP(14) + (t666 ^ 2 - t941 ^ 2) * MDP(12) + t744 * MDP(15) + t956;
t732 = -pkin(4) * t756 - pkin(5);
t656 = -pkin(5) - t797;
t652 = t721 * t880 + t888;
t651 = -t721 * t887 + t881;
t650 = -t721 * t881 + t887;
t649 = t721 * t888 + t880;
t644 = t737 - t905;
t617 = -pkin(5) - t800;
t590 = t592 + t737;
t528 = qJD(5) * t585 + t563 * t750 + t756 * t564;
t493 = pkin(5) * t528 - pkin(11) * t527 + t557;
t488 = t755 * t489;
t487 = qJD(5) * t530 + t512 * t750 - t513 * t756;
t1 = [(-t553 * t637 + t563 * t799 + t564 * t798 + t638 * t778) * MDP(19) + (t643 * t564 + t588 * t637 - t629 * t799 - t647 * t778 + t730 * t808 + t739 * t828 + t740 * t780) * MDP(23) + (t510 * t585 + t527 * t927) * MDP(25) + (-t486 * t725 + t510 * t598 + t527 * t589 - t530 * t722 + t534 * t585 + t557 * t927 - t720 * t808) * MDP(31) + (-t510 * t584 - t511 * t585 + t527 * t582 - t528 * t927) * MDP(26) + (-t487 * t725 + t511 * t598 + t528 * t589 + t534 * t584 - t557 * t582 + t721 * t808 + t722 * t803) * MDP(30) + (-g(1) * t649 - g(2) * t651 + t487 * t569 - t498 * t528 - t803 * t502 + (-(-qJD(6) * t530 + t493) * t929 - t896 + t817 * t584 - qJD(6) * t897) * t749 + t910 * t755) * MDP(38) + (-g(1) * t650 - g(2) * t652 + t487 * t567 + t488 * t584 + t497 * t528 - t803 * t503 + (t493 * t929 + t896 + (-t518 * t584 - t530 * t929 + t897) * qJD(6)) * t755 + t910 * t749) * MDP(37) + (t509 * t584 + t528 * t929) * MDP(36) + (t502 * t584 + t507 * t585 + t528 * t569 + t792 * t929) * MDP(34) + (-t585 * t506 - t503 * t584 - t528 * t567 + (-t527 * t749 - t585 * t855) * t929) * MDP(35) + (t553 * t638 - t563 * t798) * MDP(18) + (t647 * t553 + t643 * t563 + t588 * t638 - t629 * t798 - t729 * t808 - t739 * t873 - t740 * t788) * MDP(24) + (-t640 * t745 + t744 * t793) * MDP(14) + (t607 * t793 - t608 * t678 + t639 * t941 + t640 * t666) * MDP(12) + (g(1) * t890 - g(2) * t889 + t735 * t608 + t699 * t640 - t660 * t793 - t738 * t941 + t744 * t825 + t745 * t775) * MDP(16) + (t639 * t745 + t678 * t744) * MDP(13) + (-t564 * t740 - t637 * t739) * MDP(21) + (t563 * t740 + t638 * t739) * MDP(20) + (t527 * t725 + t585 * t722) * MDP(27) + (-t528 * t725 - t584 * t722) * MDP(28) + (t607 * t678 - t639 * t666) * MDP(11) + (qJDD(2) * t753 + t758 * t760) * MDP(6) + (qJDD(2) * t758 - t753 * t760) * MDP(7) + (t501 * t585 + t569 * t792) * MDP(32) + ((-t567 * t755 - t569 * t749) * t527 + (-t500 - t503 * t755 + (t567 * t749 - t569 * t755) * qJD(6)) * t585) * MDP(33) + (t753 * t795 + t758 * t784) * MDP(9) + (-t753 * t784 + t758 * t795) * MDP(10) + t808 * MDP(2) + t809 * MDP(3) + (qJDD(1) * t746 + 0.2e1 * t753 * t840) * MDP(4) + qJDD(1) * MDP(1) + 0.2e1 * (t753 * t850 - t852 * t865) * MDP(5) + (t735 * t607 + t699 * t639 + t660 * t678 - t666 * t738 - t741 * t808 - t744 * t868 - t745 * t789) * MDP(17); (t617 * t502 + t805 * t755 + t876 * t569 + (t749 * t816 + t755 * t877) * t929 + t771) * MDP(38) + t762 + (-g(3) * t758 + t753 * t783) * MDP(9) + (t617 * t503 + t805 * t749 + t876 * t567 + (t749 * t877 - t755 * t816) * t929 + t777) * MDP(37) + (g(3) * t753 + t758 * t783) * MDP(10) + qJDD(2) * MDP(8) + MDP(7) * t850 + (t644 * t799 + t812 * t739 + t740 * t942 + t764) * MDP(23) + (t644 * t798 - t661 * t739 + t740 * t943 + t769) * MDP(24) + MDP(6) * t851 + (-t590 * t927 - t722 * t870 + t725 * t877 + t946) * MDP(31) + (t869 * t745 + (t666 * t862 - t752 * t744 - t745 * t842) * pkin(2) + t766) * MDP(17) + (-t826 * t745 + (t744 * t906 - t745 * t861 + t862 * t941) * pkin(2) + t763) * MDP(16) + (t582 * t590 + t722 * t800 - t725 * t876 + t939) * MDP(30) + (-MDP(4) * t753 * t758 + MDP(5) * t865) * t761; (t656 * t502 + t804 * t755 + t871 * t569 + (t749 * t815 + t755 * t872) * t929 + t771) * MDP(38) + t762 + (t656 * t503 + t804 * t749 + t871 * t567 + (t749 * t872 - t755 * t815) * t929 + t777) * MDP(37) + (-t745 * t794 + t763) * MDP(16) + (t745 * t827 + t766) * MDP(17) + (t875 * t740 + (-t666 * t798 - t739 * t751 - t740 * t859) * pkin(3) + t769) * MDP(24) + (-t831 * t740 + (-t666 * t799 + t739 * t757 - t740 * t860) * pkin(3) + t764) * MDP(23) + (-t592 * t927 - t722 * t866 + t725 * t872 + t946) * MDP(31) + (t582 * t592 + t722 * t797 - t725 * t871 + t939) * MDP(30); (-t740 * t802 + t764) * MDP(23) + (t740 * t832 + t769) * MDP(24) + (t521 * t725 + (-t582 * t798 + t722 * t756 - t725 * t858) * pkin(4) + t939) * MDP(30) + (t522 * t725 + (-t722 * t750 - t725 * t857 + t798 * t927) * pkin(4) + t946) * MDP(31) + (t732 * t503 + t806 * t749 + t811 * t567 + (t749 * t810 - t755 * t814) * t929 + t777) * MDP(37) + (t732 * t502 + t806 * t755 + t811 * t569 + (t749 * t814 + t755 * t810) * t929 + t771) * MDP(38) + t956; (t520 * t725 + t939) * MDP(30) + (t519 * t725 + t946) * MDP(31) + (-pkin(5) * t503 - t520 * t567 + t787 * t749 + (-t933 - t919) * t755 + t954) * MDP(37) + (-pkin(5) * t502 - t520 * t569 + t787 * t755 + (-t949 + t919) * t749 + t953) * MDP(38) + t958; t569 * t567 * MDP(32) + (-t567 ^ 2 + t569 ^ 2) * MDP(33) + (t567 * t929 + t847) * MDP(34) + (t569 * t929 + t693) * MDP(35) + t509 * MDP(36) + (-g(1) * t651 + g(2) * t649 + t498 * t929 - t517 * t569 + t488) * MDP(37) + (g(1) * t652 - g(2) * t650 + t497 * t929 + t517 * t567) * MDP(38) + ((-t484 + t711) * MDP(38) + (-MDP(35) * t927 - MDP(37) * t518 - MDP(38) * t531) * qJD(6)) * t755 + (-qJD(6) * t927 * MDP(34) + (-qJD(6) * t725 - t510) * MDP(35) + (-t819 + t711) * MDP(37) + t817 * MDP(38)) * t749;];
tau  = t1;
