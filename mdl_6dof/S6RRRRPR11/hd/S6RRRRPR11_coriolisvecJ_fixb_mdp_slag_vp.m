% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:24:05
% EndTime: 2019-03-09 23:24:33
% DurationCPUTime: 19.22s
% Computational Cost: add. (16997->693), mult. (43794->951), div. (0->0), fcn. (35202->12), ass. (0->288)
t803 = sin(pkin(6));
t813 = cos(qJ(2));
t893 = qJD(1) * t813;
t872 = t803 * t893;
t782 = -qJD(3) + t872;
t809 = sin(qJ(2));
t805 = cos(pkin(6));
t894 = qJD(1) * t805;
t880 = pkin(1) * t894;
t749 = pkin(8) * t872 + t809 * t880;
t808 = sin(qJ(3));
t812 = cos(qJ(3));
t973 = t749 + t782 * (pkin(3) * t808 - pkin(10) * t812);
t810 = cos(qJ(6));
t854 = qJD(2) + t894;
t895 = qJD(1) * t803;
t873 = t809 * t895;
t726 = t808 * t854 + t812 * t873;
t807 = sin(qJ(4));
t811 = cos(qJ(4));
t678 = t726 * t807 + t811 * t782;
t680 = t726 * t811 - t782 * t807;
t802 = sin(pkin(12));
t804 = cos(pkin(12));
t836 = -t678 * t804 - t680 * t802;
t919 = t810 * t836;
t613 = t678 * t802 - t680 * t804;
t806 = sin(qJ(6));
t936 = t613 * t806;
t560 = t919 + t936;
t955 = -t808 * t873 + t812 * t854;
t719 = qJD(4) - t955;
t707 = qJD(6) + t719;
t940 = t560 * t707;
t746 = -pkin(8) * t873 + t813 * t880;
t704 = -pkin(2) * t854 - t746;
t630 = -pkin(3) * t955 - t726 * pkin(10) + t704;
t705 = pkin(9) * t854 + t749;
t742 = (-pkin(2) * t813 - pkin(9) * t809 - pkin(1)) * t803;
t718 = qJD(1) * t742;
t649 = t812 * t705 + t808 * t718;
t633 = -pkin(10) * t782 + t649;
t576 = t630 * t807 + t633 * t811;
t828 = t803 * (pkin(2) * t809 - pkin(9) * t813);
t748 = qJD(2) * t828;
t734 = qJD(1) * t748;
t925 = t803 * t809;
t791 = pkin(8) * t925;
t946 = pkin(1) * t813;
t750 = (t805 * t946 - t791) * qJD(2);
t735 = qJD(1) * t750;
t888 = qJD(3) * t812;
t890 = qJD(3) * t808;
t822 = -t705 * t890 + t718 * t888 + t808 * t734 + t812 * t735;
t882 = qJD(1) * qJD(2);
t864 = t803 * t882;
t848 = t809 * t864;
t597 = pkin(10) * t848 + t822;
t847 = t813 * t864;
t682 = qJD(3) * t955 + t812 * t847;
t832 = t808 * t847;
t683 = qJD(3) * t726 + t832;
t924 = t803 * t813;
t881 = pkin(8) * t924;
t947 = pkin(1) * t809;
t751 = (t805 * t947 + t881) * qJD(2);
t736 = qJD(1) * t751;
t605 = pkin(3) * t683 - pkin(10) * t682 + t736;
t528 = -qJD(4) * t576 - t597 * t807 + t811 * t605;
t885 = qJD(4) * t811;
t887 = qJD(4) * t807;
t606 = t811 * t682 - t726 * t887 - t782 * t885 + t807 * t848;
t504 = pkin(4) * t683 - qJ(5) * t606 - qJD(5) * t680 + t528;
t527 = t811 * t597 + t807 * t605 + t630 * t885 - t633 * t887;
t607 = qJD(4) * t680 + t682 * t807 - t811 * t848;
t507 = -qJ(5) * t607 - qJD(5) * t678 + t527;
t493 = t802 * t504 + t804 * t507;
t553 = -t606 * t802 - t607 * t804;
t491 = pkin(11) * t553 + t493;
t648 = -t808 * t705 + t718 * t812;
t632 = pkin(3) * t782 - t648;
t599 = pkin(4) * t678 + qJD(5) + t632;
t545 = -pkin(5) * t836 + t599;
t492 = t804 * t504 - t507 * t802;
t554 = t606 * t804 - t607 * t802;
t490 = pkin(5) * t683 - pkin(11) * t554 + t492;
t575 = t811 * t630 - t633 * t807;
t555 = -qJ(5) * t680 + t575;
t541 = pkin(4) * t719 + t555;
t556 = -qJ(5) * t678 + t576;
t923 = t804 * t556;
t519 = t802 * t541 + t923;
t952 = pkin(11) * t836;
t509 = t519 + t952;
t883 = qJD(6) * t806;
t861 = t806 * t490 - t509 * t883;
t972 = -t810 * t491 - t545 * t560 - t861;
t954 = -t810 * t613 + t806 * t836;
t971 = t683 * MDP(31) + (-t560 ^ 2 + t954 ^ 2) * MDP(28) - t560 * MDP(27) * t954;
t943 = pkin(9) * t807;
t970 = -t811 * t973 + t890 * t943;
t711 = t807 * t812 * t872 - t811 * t873;
t962 = -t807 * t888 + t711;
t917 = t812 * t813;
t712 = (t807 * t809 + t811 * t917) * t895;
t969 = -t811 * t888 + t712;
t747 = qJD(1) * t828;
t902 = t812 * t746 + t808 * t747;
t661 = pkin(10) * t873 + t902;
t781 = -pkin(3) * t812 - pkin(10) * t808 - pkin(2);
t968 = t811 * t661 - t781 * t885 + t807 * t973;
t941 = t954 * t707;
t664 = pkin(3) * t726 - pkin(10) * t955;
t663 = t811 * t664;
t942 = -qJ(5) - pkin(10);
t863 = qJD(4) * t942;
t966 = -pkin(4) * t726 - t663 + (qJ(5) * t955 + t863) * t811 + (-qJD(5) + t648) * t807;
t884 = qJD(5) * t811;
t909 = t811 * t648 + t807 * t664;
t929 = t955 * t807;
t965 = -qJ(5) * t929 - t807 * t863 - t884 + t909;
t918 = t811 * t812;
t794 = pkin(9) * t918;
t849 = t808 * t872;
t944 = pkin(4) * t808;
t964 = -pkin(4) * t849 + qJ(5) * t712 + t661 * t807 - t808 * t884 + (-qJ(5) * t918 + t944) * qJD(3) + (-t794 + (qJ(5) * t808 - t781) * t807) * qJD(4) + t970;
t921 = t808 * t811;
t963 = -qJ(5) * t711 - (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t921 - (-qJD(5) * t808 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t812) * t807 + t968;
t862 = t810 * t490 - t806 * t491;
t961 = -t545 * t954 + t862;
t960 = pkin(11) * t613;
t764 = t802 * t811 + t804 * t807;
t959 = t719 * t764;
t833 = t802 * t807 - t804 * t811;
t958 = t719 * t833;
t886 = qJD(4) * t808;
t957 = -t711 * t804 - t712 * t802 + t764 * t888 - t833 * t886;
t956 = t764 * t886 - t802 * t962 + t804 * t969;
t799 = t803 ^ 2;
t953 = -0.2e1 * t799 * t882;
t951 = MDP(5) * (t809 ^ 2 - t813 ^ 2);
t916 = t802 * t963 + t804 * t964;
t915 = t802 * t964 - t804 * t963;
t740 = t791 + (-pkin(2) - t946) * t805;
t757 = -t805 * t812 + t808 * t925;
t758 = t805 * t808 + t812 * t925;
t656 = pkin(3) * t757 - pkin(10) * t758 + t740;
t741 = t881 + (pkin(9) + t947) * t805;
t903 = t812 * t741 + t808 * t742;
t658 = -pkin(10) * t924 + t903;
t908 = t807 * t656 + t811 * t658;
t906 = t802 * t965 + t804 * t966;
t905 = t802 * t966 - t804 * t965;
t856 = -t808 * t746 + t747 * t812;
t660 = -pkin(3) * t873 - t856;
t949 = -pkin(4) * t962 + pkin(9) * t888 + t885 * t944 - t660;
t948 = -t649 + (t887 - t929) * pkin(4);
t860 = -t810 * t553 + t554 * t806;
t501 = qJD(6) * t954 + t860;
t814 = qJD(1) ^ 2;
t945 = pkin(4) * t802;
t851 = t705 * t888 + t718 * t890 - t812 * t734 + t808 * t735;
t598 = -pkin(3) * t848 + t851;
t939 = t598 * t807;
t938 = t598 * t811;
t937 = t606 * t807;
t935 = t678 * t719;
t934 = t680 * t719;
t933 = t683 * t807;
t932 = t683 * t811;
t931 = t683 * t812;
t930 = t955 * t782;
t928 = t726 * t782;
t824 = t782 * t808;
t927 = t782 * t812;
t926 = t799 * t814;
t549 = t802 * t556;
t922 = t807 * t808;
t518 = t804 * t541 - t549;
t505 = pkin(5) * t719 + t518 + t960;
t920 = t810 * t505;
t870 = qJD(2) * t924;
t695 = -qJD(3) * t757 + t812 * t870;
t696 = t758 * t807 + t811 * t924;
t892 = qJD(2) * t809;
t871 = t803 * t892;
t629 = -qJD(4) * t696 + t695 * t811 + t807 * t871;
t694 = qJD(3) * t758 + t808 * t870;
t878 = t807 * t924;
t697 = t758 * t811 - t878;
t821 = -t741 * t890 + t742 * t888 + t808 * t748 + t812 * t750;
t601 = pkin(10) * t871 + t821;
t621 = pkin(3) * t694 - pkin(10) * t695 + t751;
t816 = -qJD(4) * t908 - t601 * t807 + t811 * t621;
t516 = pkin(4) * t694 - qJ(5) * t629 - qJD(5) * t697 + t816;
t628 = -qJD(4) * t878 + t695 * t807 + t758 * t885 - t811 * t871;
t820 = t811 * t601 + t807 * t621 + t656 * t885 - t658 * t887;
t520 = -qJ(5) * t628 - qJD(5) * t696 + t820;
t499 = t802 * t516 + t804 * t520;
t525 = t804 * t555 - t549;
t859 = t811 * t656 - t658 * t807;
t566 = pkin(4) * t757 - qJ(5) * t697 + t859;
t577 = -qJ(5) * t696 + t908;
t534 = t802 * t566 + t804 * t577;
t744 = t764 * t808;
t745 = t833 * t808;
t835 = -t810 * t744 + t745 * t806;
t914 = qJD(6) * t835 - t806 * t957 - t810 * t956;
t667 = -t744 * t806 - t745 * t810;
t913 = qJD(6) * t667 - t806 * t956 + t810 * t957;
t834 = -t764 * t806 - t810 * t833;
t912 = qJD(6) * t834 - t806 * t959 - t810 * t958;
t687 = t764 * t810 - t806 * t833;
t911 = qJD(6) * t687 - t806 * t958 + t810 * t959;
t910 = pkin(5) * t957 + t949;
t766 = t811 * t781;
t688 = -qJ(5) * t921 + t766 + (-pkin(4) - t943) * t812;
t899 = t807 * t781 + t794;
t700 = -qJ(5) * t922 + t899;
t636 = t802 * t688 + t804 * t700;
t904 = pkin(5) * t959 + t948;
t783 = t942 * t807;
t784 = t942 * t811;
t703 = t802 * t783 - t804 * t784;
t898 = pkin(4) * t922 + t808 * pkin(9);
t891 = qJD(2) * t812;
t889 = qJD(3) * t811;
t877 = qJD(6) * t919 + t806 * t553 + t810 * t554;
t874 = -pkin(4) * t811 - pkin(3);
t867 = t719 * t887;
t498 = t804 * t516 - t520 * t802;
t524 = -t555 * t802 - t923;
t533 = t804 * t566 - t577 * t802;
t635 = t804 * t688 - t700 * t802;
t857 = -t808 * t741 + t742 * t812;
t702 = t804 * t783 + t784 * t802;
t855 = t719 * t811;
t852 = t799 * t809 * t813 * MDP(4);
t846 = pkin(1) * t953;
t657 = pkin(3) * t924 - t857;
t675 = -pkin(11) * t833 + t703;
t843 = pkin(5) * t726 - pkin(11) * t958 + qJD(6) * t675 - t906;
t674 = -pkin(11) * t764 + t702;
t842 = pkin(11) * t959 - qJD(6) * t674 - t905;
t612 = -pkin(11) * t744 + t636;
t841 = qJD(6) * t612 - t916 - t956 * pkin(11) + (t849 - t890) * pkin(5);
t608 = -pkin(5) * t812 + pkin(11) * t745 + t635;
t840 = pkin(11) * t957 - qJD(6) * t608 - t915;
t496 = t806 * t505 + t810 * t509;
t639 = -t696 * t802 + t697 * t804;
t522 = pkin(5) * t757 - pkin(11) * t639 + t533;
t638 = -t696 * t804 - t697 * t802;
t526 = pkin(11) * t638 + t534;
t839 = t522 * t810 - t526 * t806;
t838 = t522 * t806 + t526 * t810;
t837 = t810 * t638 - t639 * t806;
t583 = t638 * t806 + t639 * t810;
t831 = -t741 * t888 - t742 * t890 + t748 * t812 - t808 * t750;
t796 = pkin(4) * t804 + pkin(5);
t830 = t796 * t806 + t810 * t945;
t829 = t796 * t810 - t806 * t945;
t826 = -t719 * t885 - t933;
t823 = -pkin(10) * t683 + t632 * t719;
t500 = t613 * t883 + t877;
t819 = pkin(4) * t696 + t657;
t818 = pkin(1) * (-t805 * t882 + t926);
t602 = -pkin(3) * t871 - t831;
t542 = pkin(4) * t607 + t598;
t815 = pkin(4) * t628 + t602;
t733 = pkin(5) * t833 + t874;
t689 = pkin(5) * t744 + t898;
t659 = t683 * t757;
t590 = pkin(4) * t680 - pkin(5) * t613;
t573 = -pkin(5) * t638 + t819;
t572 = -t628 * t802 + t629 * t804;
t570 = -t628 * t804 - t629 * t802;
t531 = -pkin(5) * t570 + t815;
t521 = -pkin(5) * t553 + t542;
t513 = qJD(6) * t583 - t810 * t570 + t572 * t806;
t512 = qJD(6) * t837 + t570 * t806 + t572 * t810;
t511 = t525 + t960;
t510 = t524 - t952;
t497 = pkin(11) * t570 + t499;
t495 = -t509 * t806 + t920;
t494 = pkin(5) * t694 - pkin(11) * t572 + t498;
t488 = -t496 * qJD(6) + t862;
t487 = (qJD(6) * t505 + t491) * t810 + t861;
t1 = [(t500 * t757 + t512 * t707 + t583 * t683 + t694 * t954) * MDP(29) + (-(qJD(6) * t839 + t494 * t806 + t497 * t810) * t707 - t838 * t683 - t487 * t757 - t496 * t694 + t531 * t954 + t573 * t500 + t521 * t583 + t545 * t512) * MDP(33) + (t500 * t583 + t512 * t954) * MDP(27) + (-t682 * t757 - t683 * t758 - t694 * t726 + t695 * t955) * MDP(12) + (t694 * t782 + (t683 * t813 + (-qJD(1) * t757 + t955) * t892) * t803) * MDP(14) + (-t831 * t782 - t751 * t955 + t740 * t683 + t736 * t757 + t704 * t694 + (t851 * t813 + (qJD(1) * t857 + t648) * t892) * t803) * MDP(16) + (t500 * t837 - t501 * t583 + t512 * t560 - t513 * t954) * MDP(28) + (-t501 * t757 - t513 * t707 + t560 * t694 + t683 * t837) * MDP(30) + ((-qJD(6) * t838 + t494 * t810 - t497 * t806) * t707 + t839 * t683 + t488 * t757 + t495 * t694 - t531 * t560 + t573 * t501 - t521 * t837 + t545 * t513) * MDP(32) + (MDP(6) * t870 - MDP(7) * t871) * (qJD(2) + 0.2e1 * t894) + 0.2e1 * t852 * t882 + (t821 * t782 + t751 * t726 + t740 * t682 + t736 * t758 + t704 * t695 + (t822 * t813 + (-qJD(1) * t903 - t649) * t892) * t803) * MDP(17) + (-t492 * t639 + t493 * t638 + t498 * t613 + t499 * t836 - t518 * t572 + t519 * t570 - t533 * t554 + t534 * t553) * MDP(25) + (t528 * t757 + t575 * t694 + t598 * t696 + t602 * t678 + t657 * t607 + t632 * t628 + t683 * t859 + t719 * t816) * MDP(23) + (-t736 * t805 - t751 * t854 + t809 * t846) * MDP(9) + (-t735 * t805 - t750 * t854 + t813 * t846) * MDP(10) + t951 * t953 + (t682 * t758 + t695 * t726) * MDP(11) + (-t607 * t757 - t628 * t719 - t678 * t694 - t683 * t696) * MDP(21) + (t606 * t757 + t629 * t719 + t680 * t694 + t683 * t697) * MDP(20) + (t694 * t719 + t659) * MDP(22) + (t694 * t707 + t659) * MDP(31) + (-t606 * t696 - t607 * t697 - t628 * t680 - t629 * t678) * MDP(19) + (t606 * t697 + t629 * t680) * MDP(18) + (-t527 * t757 - t576 * t694 + t598 * t697 + t602 * t680 + t657 * t606 + t632 * t629 - t683 * t908 - t719 * t820) * MDP(24) + (-t782 * t803 - t799 * t893) * MDP(15) * t892 + (t492 * t533 + t493 * t534 + t518 * t498 + t519 * t499 + t542 * t819 + t599 * t815) * MDP(26) + (-t695 * t782 + (-t682 * t813 + (qJD(1) * t758 + t726) * t892) * t803) * MDP(13); ((-MDP(6) * t813 + MDP(7) * t809) * t803 * t805 - t852) * t814 + (-t500 * t812 + t667 * t683 + t707 * t914 - t824 * t954) * MDP(29) + (t500 * t667 + t914 * t954) * MDP(27) + (-(t608 * t806 + t612 * t810) * t683 + t487 * t812 + t689 * t500 + t521 * t667 + (t806 * t841 + t810 * t840) * t707 + t910 * t954 + t914 * t545 + t496 * t824) * MDP(33) + ((t682 - t930) * t812 + (-t683 + t928) * t808) * MDP(12) + (-pkin(2) * t683 - t736 * t812 + t856 * t782 + t749 * t955 + (pkin(9) * t927 + t704 * t808) * qJD(3) + (-t648 * t809 + (-pkin(9) * t892 - t704 * t813) * t808) * t895) * MDP(16) + (t782 * t890 + (-t813 * t824 + (-t955 + t891) * t809) * t895) * MDP(14) + (t492 * t745 - t493 * t744 + t518 * t956 - t519 * t957 + t553 * t636 - t554 * t635 + t916 * t613 + t915 * t836) * MDP(25) + (t500 * t835 - t501 * t667 + t560 * t914 - t913 * t954) * MDP(28) + (t501 * t812 - t560 * t824 + t683 * t835 - t707 * t913) * MDP(30) + ((t608 * t810 - t612 * t806) * t683 - t488 * t812 + t689 * t501 - t521 * t835 + (t806 * t840 - t810 * t841) * t707 - t910 * t560 + t913 * t545 - t495 * t824) * MDP(32) + (-t632 * t711 - t660 * t678 + t766 * t683 + ((-qJD(4) * t781 + t661) * t807 + t970) * t719 + (t632 * t807 * qJD(3) - t528 + (qJD(3) * t678 + t826) * pkin(9)) * t812 + (pkin(9) * t607 - t575 * t782 + t632 * t885 + t939) * t808) * MDP(23) + (t606 * t921 + (-t807 * t886 - t969) * t680) * MDP(18) + (-t606 * t812 - t969 * t719 + (-t680 * t782 - t867 + t932) * t808) * MDP(20) + (t607 * t812 + t962 * t719 + (t678 * t782 + t826) * t808) * MDP(21) + (-t899 * t683 - t660 * t680 - t632 * t712 + t968 * t719 + (t632 * t889 + t527 + (qJD(3) * t680 + t867) * pkin(9)) * t812 + (-t632 * t887 + t938 + t782 * t576 + (t719 * t889 + t606) * pkin(9)) * t808) * MDP(24) + (t492 * t635 + t493 * t636 + t916 * t518 + t915 * t519 + t542 * t898 + t599 * t949) * MDP(26) + t782 * MDP(15) * t873 + (-t719 * t824 - t931) * MDP(22) + (-t707 * t824 - t931) * MDP(31) + (-pkin(8) * t847 + t749 * t854 + t809 * t818) * MDP(9) + (pkin(8) * t848 + t746 * t854 + t813 * t818) * MDP(10) + t926 * t951 + (-t782 * t888 + (t782 * t917 + (qJD(2) * t808 - t726) * t809) * t895) * MDP(13) + (t682 * t808 - t726 * t927) * MDP(11) + (-pkin(2) * t682 + t736 * t808 - t902 * t782 - t749 * t726 + (-pkin(9) * t824 + t704 * t812) * qJD(3) + (-t704 * t917 + (-pkin(9) * t891 + t649) * t809) * t895) * MDP(17) + (t678 * t712 + t680 * t711 + (-t678 * t811 - t680 * t807) * t888 + (-t937 - t607 * t811 + (t678 * t807 - t680 * t811) * qJD(4)) * t808) * MDP(19); -t955 ^ 2 * MDP(12) + (t682 + t930) * MDP(13) + (-t832 - t928) * MDP(14) + MDP(15) * t848 + (-t649 * t782 - t851) * MDP(16) + (-t648 * t782 - t704 * t955 - t822) * MDP(17) + (t680 * t855 + t937) * MDP(18) + ((t606 - t935) * t811 + (-t607 - t934) * t807) * MDP(19) + (t719 * t855 + t933) * MDP(20) + (-t719 ^ 2 * t807 + t932) * MDP(21) + (-pkin(3) * t607 - t938 - t649 * t678 + (-pkin(10) * t885 - t663) * t719 + (t648 * t719 + t823) * t807) * MDP(23) + (-pkin(3) * t606 + t939 - t649 * t680 + (pkin(10) * t887 + t909) * t719 + t823 * t811) * MDP(24) + (-t492 * t764 - t493 * t833 + t518 * t958 - t519 * t959 + t553 * t703 - t554 * t702 + t906 * t613 + t905 * t836) * MDP(25) + (t492 * t702 + t493 * t703 + t906 * t518 + t905 * t519 + t542 * t874 + t599 * t948) * MDP(26) + (t500 * t687 + t912 * t954) * MDP(27) + (t500 * t834 - t501 * t687 + t560 * t912 - t911 * t954) * MDP(28) + (t683 * t687 + t707 * t912) * MDP(29) + (t683 * t834 - t707 * t911) * MDP(30) + ((t674 * t810 - t675 * t806) * t683 + t733 * t501 - t521 * t834 + (t806 * t842 - t810 * t843) * t707 - t904 * t560 + t911 * t545) * MDP(32) + (-(t674 * t806 + t675 * t810) * t683 + t733 * t500 + t521 * t687 + (t806 * t843 + t810 * t842) * t707 + t904 * t954 + t912 * t545) * MDP(33) + (-MDP(11) * t955 + t726 * MDP(12) - MDP(14) * qJD(3) - t704 * MDP(16) - t680 * MDP(20) + t678 * MDP(21) - t719 * MDP(22) - t575 * MDP(23) + t576 * MDP(24) - MDP(29) * t954 - MDP(30) * t560 - t707 * MDP(31) - t495 * MDP(32) + t496 * MDP(33)) * t726; t680 * t678 * MDP(18) + (-t678 ^ 2 + t680 ^ 2) * MDP(19) + (t606 + t935) * MDP(20) + (-t607 + t934) * MDP(21) + t683 * MDP(22) + (t576 * t719 - t632 * t680 + t528) * MDP(23) + (t575 * t719 + t632 * t678 - t527) * MDP(24) + ((t553 * t802 - t554 * t804) * pkin(4) + (t518 - t525) * t836 + (-t519 - t524) * t613) * MDP(25) + (-t518 * t524 - t519 * t525 + (t492 * t804 + t493 * t802 - t599 * t680) * pkin(4)) * MDP(26) + (t500 - t940) * MDP(29) + (-t501 + t941) * MDP(30) + (t829 * t683 - (t510 * t810 - t511 * t806) * t707 + t590 * t560 + (-t707 * t830 - t496) * qJD(6) + t961) * MDP(32) + (-t830 * t683 + (t510 * t806 + t511 * t810) * t707 - t590 * t954 + (-t707 * t829 - t920) * qJD(6) + t972) * MDP(33) + t971; (-t613 ^ 2 - t836 ^ 2) * MDP(25) + (-t518 * t613 - t519 * t836 + t542) * MDP(26) + (t501 + t941) * MDP(32) + (t500 + t940) * MDP(33); (t877 - t940) * MDP(29) + (-t860 + t941) * MDP(30) + (t496 * t707 + t961) * MDP(32) + (t495 * t707 + t972) * MDP(33) + (MDP(29) * t936 - MDP(30) * t954 - MDP(32) * t496 - MDP(33) * t920) * qJD(6) + t971;];
tauc  = t1;
