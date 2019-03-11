% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% MDP [45x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S7RRRRRRR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [7x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 08:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S7RRRRRRR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1),zeros(45,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [4x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [45 1]), ...
  'S7RRRRRRR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [45x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 07:44:11
% EndTime: 2019-03-10 07:45:45
% DurationCPUTime: 49.54s
% Computational Cost: add. (20897->1070), mult. (49401->1520), div. (0->0), fcn. (40367->12), ass. (0->373)
t764 = sin(qJ(4));
t766 = sin(qJ(2));
t770 = cos(qJ(3));
t769 = cos(qJ(4));
t771 = cos(qJ(2));
t981 = t769 * t771;
t710 = (t764 * t766 + t770 * t981) * qJD(1);
t765 = sin(qJ(3));
t954 = qJD(4) * t764;
t911 = t765 * t954;
t942 = t769 * qJD(3);
t914 = t770 * t942;
t816 = -t911 + t914;
t1051 = -t710 - t816;
t960 = qJD(1) * t771;
t918 = t770 * t960;
t956 = qJD(3) * t770;
t935 = pkin(2) * t956;
t1050 = -pkin(2) * t918 + pkin(3) * t1051 - t935;
t932 = -pkin(2) * t769 - pkin(3);
t876 = t765 * t932;
t915 = t765 * t942;
t951 = qJD(4) * t770;
t1049 = pkin(3) * qJD(3) * t765 + (t764 * t951 + t915) * pkin(2) - t876 * t960;
t763 = sin(qJ(5));
t768 = cos(qJ(5));
t920 = t765 * t960;
t953 = qJD(4) * t768;
t986 = t768 * t770;
t967 = (-qJD(5) - t942) * t986 + (t764 * t953 + (qJD(5) * t769 + qJD(3)) * t763) * t765 - t710 * t768 + t763 * t920;
t962 = qJD(1) * t766;
t709 = t764 * t918 - t769 * t962;
t952 = qJD(4) * t769;
t807 = t764 * t956 + t765 * t952 + t709;
t727 = t932 * t770;
t729 = (pkin(3) * t769 + pkin(2)) * t765;
t672 = t727 * t763 - t729 * t768;
t1031 = qJD(5) * t672 - t1049 * t768 + t1050 * t763;
t756 = t765 * qJD(2);
t961 = qJD(1) * t770;
t829 = -t766 * t961 + t756;
t884 = qJD(3) + t960;
t849 = t769 * t884;
t687 = -t764 * t829 + t849;
t836 = t769 * t829;
t851 = t764 * t884;
t689 = -t836 - t851;
t762 = sin(qJ(6));
t1020 = cos(qJ(6));
t925 = t768 * t1020;
t619 = -t687 * t925 + t762 * t689;
t899 = t1020 * qJD(5);
t878 = t768 * t899;
t1048 = t878 - t619;
t1000 = t763 * t764;
t1041 = -t687 - qJD(5);
t1022 = qJD(5) * t1041;
t822 = qJD(4) * t836;
t824 = t849 + qJD(5);
t939 = qJD(1) * qJD(2);
t898 = t766 * t939;
t964 = qJD(4) * t851 + t769 * t898;
t897 = t771 * t939;
t744 = t770 * t897;
t746 = t765 * t962;
t757 = t770 * qJD(2);
t828 = t746 + t757;
t690 = -qJD(3) * t828 + t744;
t997 = t764 * t690;
t874 = -t964 + t997;
t776 = t764 * (-0.2e1 * t822 + t874) + t824 * t952;
t1047 = t1022 * t1000 + t768 * t776;
t941 = qJD(4) - t828;
t647 = t768 * t689 + t763 * t941;
t604 = t1020 * t647 - t1041 * t762;
t711 = t768 * t941;
t645 = t689 * t763 - t711;
t639 = qJD(6) + t645;
t761 = sin(qJ(7));
t767 = cos(qJ(7));
t562 = t604 * t761 + t639 * t767;
t669 = t1020 * t1041;
t602 = t647 * t762 + t669;
t940 = qJD(7) - t602;
t1046 = t562 * t940;
t564 = t604 * t767 - t761 * t639;
t1045 = t564 * t940;
t987 = t768 * t769;
t837 = t763 * t770 + t765 * t987;
t900 = qJD(6) * t1020;
t877 = t764 * t900;
t854 = t765 * t877;
t974 = -t854 + (qJD(6) * t837 - t807) * t762 + t967 * t1020;
t947 = qJD(6) * t762;
t906 = t763 * t947;
t1044 = -t906 + t1048;
t1001 = t762 * t764;
t659 = t763 * t829 - t828 * t987;
t879 = t763 * t899;
t901 = qJD(4) * t1020;
t946 = qJD(6) * t768;
t969 = t1001 * t828 - t1020 * t659 + (-t768 * t901 + t900) * t769 + (t879 + (-qJD(4) + t946) * t762) * t764;
t982 = t769 * t770;
t992 = t765 * t768;
t966 = t992 * qJD(3) + t837 * qJD(5) + t768 * t920 + (qJD(3) * t982 + t710 - t911) * t763;
t1005 = t828 * t769;
t658 = -t1005 * t763 - t768 * t829;
t948 = qJD(5) * t768;
t813 = t763 * t952 + t764 * t948;
t806 = t658 + t813;
t759 = t766 ^ 2;
t848 = t771 * t884;
t1043 = -t759 * qJD(1) + t848;
t804 = t762 * t952 + t877;
t1023 = qJD(3) * t941;
t1042 = t1023 * t769;
t1036 = t940 ^ 2;
t817 = t756 * t771 + t766 * t956;
t691 = t817 * qJD(1) - t756 * qJD(3);
t1039 = pkin(2) * t691;
t724 = t828 * pkin(2);
t983 = t769 * t724;
t660 = pkin(3) * t941 - t983;
t723 = t829 * pkin(2);
t794 = pkin(3) * t689 - t723;
t600 = t763 * t660 + t768 * t794;
t965 = t769 * t690 + t764 * t898;
t614 = -qJD(4) * t687 + t965;
t686 = qJD(2) * t935 + (qJD(3) * t746 - t744) * pkin(2);
t592 = -pkin(3) * t614 + t686;
t985 = t769 * t1039;
t606 = -pkin(3) * t691 + t724 * t954 - t985;
t857 = t763 * t592 + t768 * t606;
t544 = -qJD(5) * t600 + t857;
t601 = t768 * t660 - t763 * t794;
t823 = t1001 * t724 - t1020 * t601;
t880 = t769 * t901;
t927 = t764 * t1020;
t532 = qJD(6) * t823 - t1039 * t927 - t762 * t544 - t724 * t880;
t1040 = pkin(4) * t1036 + t532;
t1038 = MDP(4) * t766;
t1027 = -t771 ^ 2 + t759;
t1037 = MDP(5) * t1027;
t1035 = t645 * t1041;
t1034 = t687 * t884;
t1033 = t764 * t1041;
t798 = -t822 - t964;
t791 = t798 + t997;
t1032 = t768 * t791;
t673 = t727 * t768 + t729 * t763;
t1030 = qJD(5) * t673 + t1049 * t763 + t1050 * t768;
t656 = (t757 - qJD(4)) * t981 + (-t915 + (qJD(2) - t951) * t764) * t766;
t993 = t765 * t766;
t1029 = -qJD(5) * t993 + t656;
t1019 = pkin(2) * t770;
t860 = pkin(2) * t764 * t920;
t995 = t764 * t765;
t933 = t762 * t995;
t890 = pkin(2) * t933;
t1026 = -qJD(3) * t890 + t1019 * t804 + t1020 * t1031 + t673 * t947 - t762 * t860;
t830 = t829 * qJD(3);
t835 = t771 * t829;
t1025 = qJD(1) * t835 + t830;
t1024 = qJD(3) * t884;
t1021 = t724 * t829;
t949 = qJD(5) * t763;
t560 = qJD(5) * t711 + t768 * t614 - t689 * t949 - t763 * t691;
t540 = -qJD(6) * t669 + t1020 * t560 - t647 * t947 + t762 * t791;
t539 = t767 * t540;
t894 = t763 * t614 + t768 * t691;
t561 = qJD(5) * t647 + t894;
t527 = -qJD(7) * t562 - t561 * t761 + t539;
t1018 = t527 * t761;
t1017 = t540 * t762;
t895 = -t768 * t592 + t763 * t606;
t545 = qJD(5) * t601 + t895;
t1016 = t545 * t762;
t1015 = t560 * t763;
t1014 = t561 * t768;
t1013 = t564 * t767;
t1012 = t604 * t639;
t1011 = t614 * t764;
t1010 = t647 * t768;
t1009 = t687 * t763;
t1008 = t687 * t768;
t1007 = t691 * t770;
t1006 = t724 * t764;
t896 = -t1020 * t791 + t762 * t560;
t541 = qJD(6) * t604 + t896;
t1004 = t761 * t541;
t1003 = t762 * t561;
t1002 = t762 * t763;
t999 = t763 * t767;
t998 = t764 * t1039;
t996 = t764 * t691;
t994 = t764 * t768;
t991 = t765 * t769;
t990 = t766 * t770;
t989 = t767 * t541;
t988 = t767 * t768;
t984 = t769 * t691;
t679 = -t1020 * t837 - t933;
t717 = t763 * t991 - t986;
t635 = t679 * t761 - t717 * t767;
t980 = -qJD(7) * t635 + t761 * t966 + t767 * t974;
t945 = qJD(7) * t761;
t979 = t717 * t945 + (qJD(7) * t679 - t966) * t767 + t974 * t761;
t902 = qJD(3) * t1020;
t978 = -t673 * t900 + (-t770 * t880 + (t765 * t902 + t770 * t947) * t764) * pkin(2) + t1020 * t860 + t1031 * t762;
t926 = t767 * t1020;
t720 = t761 * t768 + t763 * t926;
t977 = -t720 * qJD(7) + t1041 * t999 + (-qJD(5) * t925 + t619 + t906) * t761;
t721 = t762 * t769 - t764 * t925;
t976 = (qJD(7) * t1000 + t969) * t767 + (-qJD(7) * t721 + t806) * t761;
t678 = t1000 * t761 + t721 * t767;
t975 = qJD(7) * t678 + t761 * t969 - t767 * t806;
t853 = t765 * t880;
t881 = t770 * t902;
t905 = t764 * t947;
t973 = t1020 * t709 + t762 * t967 + t764 * t881 - t765 * t905 - t837 * t900 + t853;
t663 = pkin(3) * t829 + t723 * t769;
t670 = pkin(3) * t1005 + t724;
t610 = t663 * t768 + t670 * t763;
t972 = t723 * t1001 + t1020 * t610;
t633 = pkin(3) * t1009 + t724 * t994;
t971 = t1020 * t633 - t762 * t983;
t970 = -t1009 * t761 - t619 * t767 + (t899 + qJD(7)) * t988 + (-t767 * t947 + (-qJD(7) * t1020 - qJD(5)) * t761) * t763;
t909 = t764 * t949;
t968 = t828 * t927 - t764 * t901 + t768 * t804 - t769 * t947 + (t659 - t909) * t762;
t959 = qJD(2) * t766;
t958 = qJD(2) * t771;
t957 = qJD(3) * t766;
t955 = qJD(3) * t771;
t950 = qJD(5) * t762;
t944 = qJD(7) * t767;
t938 = 0.2e1 * t983;
t937 = pkin(2) * t993;
t936 = t1020 * pkin(4);
t934 = t762 * t1000;
t931 = t541 * t1020;
t930 = t545 * t1020;
t929 = t761 * t1020;
t928 = t763 * t1020;
t924 = t1020 * t561;
t923 = t1020 * t639;
t922 = t1020 * t645;
t921 = t1020 * t769;
t917 = t771 * t757;
t916 = t765 * t957;
t910 = t766 * t952;
t889 = -qJD(3) * t723 - t1039;
t530 = t540 * pkin(4) + t545;
t557 = -t639 * pkin(4) - t823;
t888 = qJD(7) * t557 + t530;
t887 = qJD(7) * t604 + t561;
t886 = t940 * t936;
t885 = pkin(3) * t925;
t883 = t765 * t927;
t882 = pkin(3) * t1020 + pkin(4);
t632 = pkin(3) * t1008 - t1000 * t724;
t873 = t771 * t941;
t871 = t639 * t950 + t541;
t618 = -t1008 * t762 - t1020 * t689;
t870 = -t762 * t948 + t618;
t621 = pkin(4) * t679 + t672;
t869 = pkin(4) * t966 - qJD(7) * t621 - t1026;
t852 = -t1001 * t1019 + t1020 * t673;
t620 = pkin(4) * t717 + t852;
t868 = pkin(4) * t974 + qJD(7) * t620 + t1030;
t692 = -pkin(3) * t994 + pkin(4) * t721;
t867 = -qJD(7) * t692 + (t764 * t878 + (t880 - t905) * t763) * pkin(3) - t972 + t806 * pkin(4);
t609 = t663 * t763 - t670 * t768;
t712 = t882 * t1000;
t814 = -t768 * t952 + t909;
t866 = pkin(3) * t814 + pkin(4) * t969 + qJD(7) * t712 - t609;
t726 = t882 * t768;
t865 = pkin(3) * t948 + pkin(4) * t1044 + qJD(7) * t726 + t632;
t728 = (t936 + pkin(3)) * t763;
t864 = -qJD(7) * t728 + (-t762 * t946 - t879) * pkin(3) - t971 + (-t1009 - t949) * pkin(4);
t862 = qJD(4) * t941;
t556 = pkin(4) * t604 + t600;
t537 = -t556 * t767 - t557 * t761;
t719 = -t764 * t771 + t766 * t982;
t681 = t719 * t768 - t763 * t993;
t718 = t764 * t990 + t981;
t638 = t1020 * t681 + t718 * t762;
t693 = -pkin(2) * t990 - pkin(3) * t719;
t713 = t766 * t876;
t856 = t693 * t768 - t713 * t763;
t581 = pkin(4) * t638 - t856;
t680 = t719 * t763 + t766 * t992;
t643 = t693 * t763 + t713 * t768;
t831 = t1020 * t643 - t766 * t890;
t593 = -pkin(4) * t680 + t831;
t859 = -t581 * t767 - t593 * t761;
t858 = t581 * t761 - t593 * t767;
t598 = t638 * t761 + t680 * t767;
t855 = t723 * t765 + t724 * t770;
t850 = t765 * t884;
t847 = t884 * t689;
t846 = qJD(1) * t884;
t843 = t770 * t1023;
t842 = t941 * t954;
t841 = t769 * t862;
t840 = qJD(2) * t873;
t839 = t639 * t899 + t540;
t575 = -t761 * t647 - t767 * t922;
t838 = -t767 * t900 + t575;
t834 = t1020 * t718 - t681 * t762;
t833 = t940 * t944 - t1004;
t832 = t940 * t945 + t989;
t827 = t766 * t1024;
t826 = t766 * t843;
t825 = t724 * t884 - t686;
t821 = qJD(2) * t835;
t819 = -t639 * t900 - t1003;
t818 = -t639 * t947 + t924;
t811 = t765 * t827;
t810 = pkin(4) * t541;
t809 = t900 + t922;
t531 = t1020 * t544 - t601 * t947 - t724 * t804 - t762 * t998;
t808 = -t836 + t689;
t802 = t841 - t996;
t800 = t763 * t900 - t870;
t797 = qJD(5) * t604 - t818;
t796 = qJD(5) * t602 + t819;
t795 = qJD(5) * t719 + t817;
t793 = -t828 * t941 + t862;
t792 = qJD(1) * t873 + t1023;
t617 = -pkin(3) * t656 + (t916 - t917) * pkin(2);
t641 = -t817 * pkin(3) + (-t766 * t914 + (t766 * t954 - t769 * t958) * t765) * pkin(2);
t555 = qJD(5) * t643 - t617 * t768 + t641 * t763;
t789 = t771 * t846 + t1024;
t787 = t763 * t791;
t786 = -t764 * t817 - t765 * t910;
t784 = -t762 * t813 - t763 * t877;
t780 = t768 * t1022;
t554 = qJD(5) * t856 + t617 * t763 + t641 * t768;
t775 = -t643 * t947 + t1020 * t554 + (t762 * t786 - t766 * t854) * pkin(2);
t774 = -t1041 * t687 - t1022;
t758 = t764 ^ 2;
t716 = t762 * t994 + t921;
t715 = t761 * t928 - t988;
t677 = -t762 * t837 + t883;
t676 = t721 * t761 - t764 * t999;
t657 = -t1019 * t927 - t762 * t673;
t655 = -t764 * t916 - t771 * t954 - t769 * t959 + (t764 * t958 + t910) * t770;
t636 = t679 * t767 + t717 * t761;
t630 = t639 * t945;
t615 = -pkin(2) * t766 * t883 - t762 * t643;
t605 = -t762 * t633 - t724 * t921;
t599 = t638 * t767 - t680 * t761;
t591 = -t762 * t610 + t723 * t927;
t587 = t1029 * t768 - t795 * t763;
t586 = t1029 * t763 + t768 * t795;
t574 = t767 * t647 - t761 * t922;
t571 = -t762 * t601 - t724 * t927;
t566 = -pkin(4) * t922 + t601;
t565 = -t647 * pkin(4) - t1020 * t600;
t552 = t556 * t945;
t550 = qJD(6) * t834 + t1020 * t587 + t655 * t762;
t549 = qJD(6) * t638 - t1020 * t655 + t587 * t762;
t542 = -t643 * t900 - t762 * t554 + (-t766 * t853 + (-t766 * t881 + (-t1020 * t958 + t766 * t947) * t765) * t764) * pkin(2);
t538 = -t556 * t761 + t557 * t767;
t536 = pkin(4) * t550 + t555;
t535 = -t586 * pkin(4) + t775;
t534 = -qJD(7) * t598 + t550 * t767 - t586 * t761;
t533 = t550 * t761 - t680 * t945 + (qJD(7) * t638 + t586) * t767;
t529 = -pkin(4) * t561 + t531;
t528 = t540 * t761 + t767 * t887 - t630;
t526 = -t529 * t761 - t767 * t888 + t552;
t525 = qJD(7) * t537 + t529 * t767 - t530 * t761;
t1 = [(t1041 * t586 - t561 * t718 - t645 * t655 - t680 * t791) * MDP(28) + (-t1041 * t587 + t560 * t718 + t647 * t655 + t681 * t791) * MDP(27) + (-t1041 * t655 + t718 * t791) * MDP(29) + (-t770 * t827 - t691 * t771 + (-t1043 * t765 + t766 * t828) * qJD(2)) * MDP(14) + (-t723 * t959 + t686 * t771 + (-t1043 * t757 + t811) * pkin(2)) * MDP(16) + (-t614 * t993 + t656 * t941 - t689 * t817 - t719 * t691) * MDP(20) + (-t655 * t941 + t687 * t817 + t718 * t691 + t791 * t993) * MDP(21) + ((-t828 * t958 + (-t691 + t830) * t766) * t770 + (-t690 * t766 + t828 * t957 + t821) * t765) * MDP(12) + (t690 * t990 - t770 * t821 + t830 * t993) * MDP(11) + (t723 * t656 + t686 * t719 + (-t765 * t983 + (-t770 * t689 + t941 * t991) * pkin(2)) * t958 + (-t765 * t985 - t816 * t724 + ((-t614 + t1042) * t770 + (t689 * qJD(3) - t842 - t984) * t765) * pkin(2)) * t766) * MDP(24) + (-t811 + t690 * t771 + (t756 * t766 + (t955 + (-t759 - t1027) * qJD(1)) * t770) * qJD(2)) * MDP(13) + (MDP(6) * t771 - MDP(7) * t766) * qJD(2) ^ 2 + (-t993 * t998 + t723 * t655 + t686 * t718 + t786 * t724 + (t764 * t826 - t687 * t917 - t791 * t990 + (t764 * t840 + (t687 * qJD(3) + t802) * t766) * t765) * pkin(2)) * MDP(23) + (-t614 * t718 - t689 * t655 - t656 * t687 - t719 * t791) * MDP(19) + (-t826 + (t691 * t766 - t840) * t765) * MDP(22) - 0.2e1 * t939 * t1037 + 0.2e1 * t897 * t1038 + (-t724 * t959 + t1039 * t771 + (t1024 * t990 + t1043 * t756) * pkin(2)) * MDP(17) + (-t541 * t680 - t549 * t639 + t561 * t834 - t586 * t602) * MDP(35) + (t540 * t834 - t541 * t638 - t549 * t604 - t550 * t602) * MDP(33) + (-t541 * t834 - t549 * t940) * MDP(43) + (-t528 * t834 - t533 * t940 + t541 * t598 + t549 * t562) * MDP(42) + (t527 * t834 + t534 * t940 - t541 * t599 - t549 * t564) * MDP(41) + ((qJD(7) * t858 - t535 * t761 - t536 * t767) * t940 - t859 * t541 + t526 * t834 - t537 * t549 + t542 * t562 + t615 * t528 + t532 * t598 + t571 * t533) * MDP(44) + (-(qJD(7) * t859 + t535 * t767 - t536 * t761) * t940 - t858 * t541 - t525 * t834 + t538 * t549 + t542 * t564 + t615 * t527 + t532 * t599 + t571 * t534) * MDP(45) + (-t531 * t680 - t540 * t856 + t545 * t638 + t600 * t550 + t555 * t604 - t561 * t831 + t586 * t823 - t639 * t775) * MDP(38) + (t532 * t680 - t541 * t856 + t542 * t639 - t545 * t834 + t549 * t600 + t555 * t602 + t561 * t615 + t571 * t586) * MDP(37) + (t614 * t719 + t656 * t689) * MDP(18) + (t540 * t680 + t550 * t639 + t561 * t638 + t586 * t604) * MDP(34) + (t561 * t680 + t586 * t639) * MDP(36) + (-t560 * t680 - t561 * t681 - t586 * t647 - t587 * t645) * MDP(26) + (t560 * t681 + t587 * t647) * MDP(25) + (t540 * t638 + t550 * t604) * MDP(32) + (-t527 * t598 - t528 * t599 - t533 * t564 - t534 * t562) * MDP(40) + (t527 * t599 + t534 * t564) * MDP(39) + (-qJD(3) - 0.2e1 * t960) * MDP(15) * t959 + (-t555 * t824 - t856 * t964 - t545 * t718 - t600 * t655 + (-t645 * t937 - t724 * t680 - t829 * t856) * t952 + (-t1039 * t680 - t724 * t586 + t555 * t829 + t856 * t690 + (-t561 * t993 - t645 * t817) * pkin(2)) * t764) * MDP(30) + (-t554 * t824 + t643 * t964 - t544 * t718 - t601 * t655 + (t643 * t829 - t647 * t937 - t724 * t681) * t952 + (-t1039 * t681 - t724 * t587 + t554 * t829 - t643 * t690 + (-t560 * t993 - t647 * t817) * pkin(2)) * t764) * MDP(31); (-t710 * t941 + (-t614 - t1042) * t770 + (t764 * t862 + t847 + t984) * t765) * MDP(20) + (-t1041 * t967 - t560 * t995 - t647 * t807 - t791 * t837) * MDP(27) + (-t1041 * t966 + t561 * t995 + t645 * t807 + t717 * t791) * MDP(28) + (t1033 * t956 + t1041 * t709 - t765 * t776) * MDP(29) + (-t723 * t710 + t855 * t954 + (t765 * t614 + (-t842 + t847) * t770) * pkin(2) + ((t889 - t1039) * t770 + ((-t941 * t960 - t1023) * pkin(2) + t825) * t765) * t769) * MDP(24) + (t764 * t843 + t791 * t770 + t709 * t941 + (t802 - t1034) * t765) * MDP(21) + (-qJD(3) ^ 2 * t770 + (qJD(1) * t1027 - 0.2e1 * t955) * t961) * MDP(13) + (-t560 * t837 + t647 * t967) * MDP(25) + (t560 * t717 + t561 * t837 - t645 * t967 + t647 * t966) * MDP(26) + ((-t828 + t757) * t962 + t789 * t765) * MDP(14) + ((t828 * t884 - t690) * t770 + (t691 - t1025) * t765) * MDP(12) + (t724 * t962 + (-t770 * t898 + (-qJD(1) * t848 - t1024) * t765) * pkin(2)) * MDP(17) + (-t1038 * t771 + t1037) * qJD(1) ^ 2 + (t1030 * t602 - t532 * t717 + t541 * t672 + t545 * t677 + t561 * t657 - t966 * t571 + t973 * t600 + t978 * t639) * MDP(37) + (-t723 * t709 - t855 * t952 + (t765 * t825 + t770 * t889) * t764 + (t765 * t798 + (-t1007 + (t690 - t792) * t765) * t764 + (t841 + t1034) * t770) * pkin(2)) * MDP(23) + (t1025 * t770 - t690 * t765) * MDP(11) + (t1026 * t639 + t1030 * t604 + t531 * t717 + t672 * t540 + t545 * t679 - t561 * t852 + t600 * t974 - t823 * t966) * MDP(38) + ((t620 * t767 - t621 * t761) * t541 + t525 * t677 + t657 * t527 + t532 * t636 - (-t761 * t868 + t767 * t869) * t940 + t980 * t571 + t978 * t564 + t973 * t538) * MDP(45) + (-t527 * t677 - t541 * t636 - t564 * t973 + t940 * t980) * MDP(41) + (-(-t620 * t761 - t621 * t767) * t541 - t526 * t677 + t657 * t528 + t532 * t635 - (t761 * t869 + t767 * t868) * t940 + t979 * t571 + t978 * t562 - t973 * t537) * MDP(44) + (t528 * t677 + t541 * t635 + t562 * t973 - t940 * t979) * MDP(42) + (t541 * t677 - t940 * t973) * MDP(43) + (t710 * t687 + t689 * t709 + (t687 * t769 + t689 * t764) * t956 + (t769 * t874 + t1011 + (-t687 * t764 + t769 * t808) * qJD(4)) * t765) * MDP(19) + (t765 * t792 + t1007) * MDP(22) + (t527 * t636 + t564 * t980) * MDP(39) + (-t527 * t635 - t528 * t636 - t562 * t980 - t564 * t979) * MDP(40) + (t541 * t717 - t561 * t677 + t602 * t966 - t639 * t973) * MDP(35) + (t540 * t679 + t604 * t974) * MDP(32) + (-t540 * t717 + t561 * t679 - t604 * t966 + t639 * t974) * MDP(34) + (-t540 * t677 - t541 * t679 - t602 * t974 - t604 * t973) * MDP(33) + (t723 * t962 + (-t765 * t898 + t770 * t789) * pkin(2)) * MDP(16) + (-t561 * t717 - t639 * t966) * MDP(36) + (t673 * t964 + t601 * t709 + (-t1019 * t647 + t601 * t765 + t673 * t829 + t724 * t837) * t952 + (t1039 * t837 - t673 * t690 + t544 * t765 + t601 * t956 - t967 * t724 + (-t770 * t560 + t647 * t850) * pkin(2) - t1031 * t829) * t764 + t1031 * t824) * MDP(31) + (t672 * t964 + t600 * t709 + (-t1019 * t645 + t600 * t765 + t672 * t829 + t724 * t717) * t952 + (t1039 * t717 - t672 * t690 + t545 * t765 + t600 * t956 + t966 * t724 + (-t770 * t561 + t645 * t850) * pkin(2) + t1030 * t829) * t764 - t1030 * t824) * MDP(30) + (t1051 * t689 - t614 * t991) * MDP(18) + t766 * MDP(15) * t846; (-t647 * t764 * t941 + t1041 * t659 + t560 * t769 - t1047) * MDP(27) + (-t544 * t769 + t610 * t824 + (t1039 * t768 - t724 * t949) * t758 + (-t723 * t647 + t724 * t659 - t610 * t829 - t601 * t828 + (t768 * t938 + t601) * qJD(4)) * t764 + (t1022 * t994 + (t1041 * t952 - t764 * t791) * t763) * pkin(3)) * MDP(31) + (-t561 * t769 - t658 * t1041 + t763 * t776 + (t645 * t941 - t780) * t764) * MDP(28) + (-t1021 * t769 - t686 * t764 - t724 * t689) * MDP(24) + (-t1021 * t764 + t686 * t769 - t724 * t687) * MDP(23) + (t828 * t960 + t744) * MDP(13) + ((t687 * t941 - t614) * t769 + (qJD(4) * t808 - t689 * t828 + t874) * t764) * MDP(19) + (t723 * t884 + t1039) * MDP(17) - t941 * t829 * MDP(22) + (-t545 * t769 + t609 * t824 + (t1039 * t763 + t724 * t948) * t758 + (-t723 * t645 + t724 * t658 - t609 * t829 - t600 * t828 + (t763 * t938 + t600) * qJD(4)) * t764 + t1047 * pkin(3)) * MDP(30) + (t1033 * t941 + t769 * t791) * MDP(29) + (-t756 - t829) * MDP(14) * t960 - t825 * MDP(16) - t829 * t828 * MDP(11) + (t531 * t1000 + t545 * t721 + t972 * t639 - t609 * t604 + t969 * t600 - t806 * t823 + ((-t604 * t768 - t763 * t923) * t952 + (t763 * t797 - t768 * t839) * t764) * pkin(3)) * MDP(38) + ((-t692 * t761 + t712 * t767) * t541 - t525 * t716 + t532 * t678 - t591 * t564 - (-t761 * t866 + t767 * t867) * t940 + t976 * t571 - t968 * t538 + (-t527 * t934 + t564 * t784) * pkin(3)) * MDP(45) + (t527 * t716 - t541 * t678 + t564 * t968 + t940 * t976) * MDP(41) + (-(-t692 * t767 - t712 * t761) * t541 + t526 * t716 + t532 * t676 - t591 * t562 - (t761 * t867 + t767 * t866) * t940 + t975 * t571 + t968 * t537 + (-t528 * t934 + t562 * t784) * pkin(3)) * MDP(44) + (-t528 * t716 + t541 * t676 - t562 * t968 - t940 * t975) * MDP(42) + (-t541 * t716 + t940 * t968) * MDP(43) + (t645 * t659 + t647 * t658 + (t645 * t768 + t647 * t763) * t952 + (t1015 + t1014 + (-t645 * t763 + t1010) * qJD(5)) * t764) * MDP(26) + (-t689 * t769 * t941 - t1011) * MDP(18) + (-t1000 * t561 - t639 * t806) * MDP(36) + (-t1000 * t540 + t561 * t721 - t604 * t806 + t639 * t969) * MDP(34) + (t1000 * t541 + t561 * t716 + t602 * t806 + t639 * t968) * MDP(35) + (-t532 * t1000 - t545 * t716 - t591 * t639 - t609 * t602 - t968 * t600 - t806 * t571 + ((-t1002 * t639 - t602 * t768) * t952 + (t763 * t796 - t768 * t871) * t764) * pkin(3)) * MDP(37) + (-t560 * t994 + (-t659 + t814) * t647) * MDP(25) + (t527 * t678 + t564 * t976) * MDP(39) + (-t527 * t676 - t528 * t678 - t562 * t976 - t564 * t975) * MDP(40) + (t540 * t721 + t604 * t969) * MDP(32) + (t540 * t716 - t541 * t721 - t602 * t969 + t604 * t968) * MDP(33) + (-t689 * t829 - t769 * t793 + t996) * MDP(20) + (t687 * t829 + t764 * t793 + t984) * MDP(21) + (-t828 ^ 2 + t829 ^ 2) * MDP(12) - MDP(15) * t898; -t791 * MDP(21) + ((t560 + t1035) * t768 + (t1041 * t647 - t561) * t763) * MDP(26) + (t971 * t639 + t632 * t604 - t600 * t619 + (pkin(3) * t797 + t600 * t899 + t531) * t768 + (pkin(3) * t839 - t1041 * t823 - t600 * t947 + t930) * t763) * MDP(38) + (-t1041 * t639 * t763 - t1014) * MDP(36) + (-t1010 * t1041 + t1015) * MDP(25) + (-t600 * t618 + t632 * t602 - t605 * t639 + (pkin(3) * t796 + t600 * t950 - t532) * t768 + (pkin(3) * t871 - t1041 * t571 + t600 * t900 + t1016) * t763) * MDP(37) + (t541 * t768 + t870 * t639 + (t1041 * t602 + t819) * t763) * MDP(35) + (MDP(18) * t687 + MDP(19) * t689 + MDP(21) * t941 - t723 * MDP(23) - t647 * MDP(27) + t645 * MDP(28) + MDP(29) * t1041 + t600 * MDP(30) + t601 * MDP(31)) * t689 + (t1039 * t994 + t632 * t1041 + (t780 - t787) * pkin(3) + ((t645 + t953) * t769 + t763 * t1033) * t724) * MDP(30) + (t768 * t774 + t787) * MDP(27) + (-t1006 * t828 + t723 * t687 + t985) * MDP(24) + (t828 * t983 + t998) * MDP(23) + (-t540 * t768 + t1048 * t639 + (-t1041 * t604 + t818) * t763) * MDP(34) + (-t763 * t998 - t633 * t1041 + (-t1022 * t763 - t1032) * pkin(3) + ((-qJD(4) * t763 + t647) * t769 + t768 * t1033) * t724) * MDP(31) - t687 ^ 2 * MDP(19) + (-t1002 * t527 - t720 * t541 - t564 * t800 + t940 * t970) * MDP(41) + (t1002 * t541 - t800 * t940) * MDP(43) + (t1002 * t528 + t715 * t541 + t562 * t800 + t940 * t977) * MDP(42) + (-(-t726 * t761 - t728 * t767) * t541 + t532 * t715 + t537 * t618 - t605 * t562 - (t761 * t864 + t767 * t865) * t940 - t977 * t571 + (-t537 * t928 - t562 * t885) * qJD(6) + (-t537 * t948 - t526 * t763 + (-t528 * t768 + t562 * t949) * pkin(3)) * t762) * MDP(44) + ((t726 * t767 - t728 * t761) * t541 + t532 * t720 - t538 * t618 - t605 * t564 - (-t761 * t865 + t767 * t864) * t940 + t970 * t571 + (t538 * t928 - t564 * t885) * qJD(6) + (t538 * t948 + t525 * t763 + (-t527 * t768 + t564 * t949) * pkin(3)) * t762) * MDP(45) + (-t763 * t774 + t1032) * MDP(28) - t691 * MDP(22) + (t1044 * t604 + t540 * t928) * MDP(32) + (t619 * t602 + t604 * t618 + (-t1020 * t602 - t604 * t762) * t948 + (-t931 - t1017 + (-t1020 * t604 + t602 * t762) * qJD(6)) * t763) * MDP(33) + (-t527 * t715 - t528 * t720 - t562 * t970 + t564 * t977) * MDP(40) + (t527 * t720 + t564 * t970) * MDP(39) + (-t687 * t828 + t965) * MDP(20); -t645 ^ 2 * MDP(26) + (t560 - t1035) * MDP(27) - t894 * MDP(28) + t791 * MDP(29) + (t601 * t687 - t895) * MDP(30) + (-t645 * t1006 - t600 * t687 - t857) * MDP(31) + (t604 * t809 + t1017) * MDP(32) + (t540 * t1020 - t809 * t602 + (-t541 - t1012) * t762) * MDP(33) + (t639 * t809 + t1003) * MDP(34) + (-t639 ^ 2 * t762 + t924) * MDP(35) + (-t601 * t602 - t930) * MDP(37) + (t1016 - t601 * t604 + (-t923 + t809) * t600) * MDP(38) + (t527 * t767 * t762 + (-t762 * t945 - t838) * t564) * MDP(39) + (t575 * t562 + t564 * t574 + (-t562 * t926 - t564 * t929) * qJD(6) + (-t1018 - t528 * t767 + (t562 * t761 - t1013) * qJD(7)) * t762) * MDP(40) + (t527 * t1020 - t838 * t940 + (-t564 * t639 - t832) * t762) * MDP(41) + (-t528 * t1020 - (t761 * t900 - t574) * t940 + (t562 * t639 - t833) * t762) * MDP(42) + (-t639 * t762 * t940 - t931) * MDP(43) + (t526 * t1020 - (-t565 * t761 - t566 * t767) * t940 - t571 * t574 + (t571 * t929 - t767 * t886) * qJD(6) + (pkin(4) * t832 + t761 * t532 - t537 * t639 - t600 * t562 + t571 * t944) * t762) * MDP(44) + (-t525 * t1020 + (t565 * t767 - t566 * t761) * t940 - t571 * t575 + (t571 * t926 + t761 * t886) * qJD(6) + (pkin(4) * t833 + t532 * t767 + t538 * t639 - t600 * t564 - t571 * t945) * t762) * MDP(45) + (MDP(25) * t645 + MDP(26) * t647 + MDP(28) * t687 + MDP(30) * t1006 - MDP(34) * t604 + MDP(35) * t602 - MDP(36) * t639 - MDP(37) * t571 - MDP(38) * t823) * t647; -t602 ^ 2 * MDP(33) + (t602 * t639 + t540) * MDP(34) + (t1012 - t896) * MDP(35) + t561 * MDP(36) + (-t639 * t823 + t532) * MDP(37) + (t571 * t639 + t600 * t602 - t531) * MDP(38) + (-t1013 * t940 - t1018) * MDP(39) + ((-t527 + t1046) * t767 + (t528 + t1045) * t761) * MDP(40) + (-t1036 * t767 + t1004) * MDP(41) + (t1036 * t761 + t989) * MDP(42) + (t1040 * t767 - t823 * t562 - t810 * t761) * MDP(44) + (-t1040 * t761 - t823 * t564 - t810 * t767) * MDP(45) + (MDP(32) * t602 + MDP(33) * t604 - MDP(35) * qJD(6) - MDP(37) * t600 + MDP(41) * t564 - MDP(42) * t562 + MDP(43) * t940 + MDP(44) * t537 - MDP(45) * t538) * t604; t564 * t562 * MDP(39) + (-t562 ^ 2 + t564 ^ 2) * MDP(40) + (t539 + t1046) * MDP(41) + (t630 + t1045) * MDP(42) - t541 * MDP(43) + (t538 * t940 - t564 * t571 + t552) * MDP(44) + (t537 * t940 + t562 * t571) * MDP(45) + (-MDP(41) * t887 - t540 * MDP(42) - t529 * MDP(44) + MDP(45) * t888) * t761 + (-qJD(7) * t639 * MDP(41) - t887 * MDP(42) - t888 * MDP(44) + (qJD(7) * t556 - t529) * MDP(45)) * t767;];
tauc  = t1;
