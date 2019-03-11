% Calculate inertial parameters regressor of coriolis matrix for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRRR6_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:48:52
% EndTime: 2019-03-08 20:49:22
% DurationCPUTime: 23.25s
% Computational Cost: add. (14345->772), mult. (32761->1096), div. (0->0), fcn. (35416->10), ass. (0->559)
t689 = cos(pkin(6));
t692 = sin(qJ(4));
t695 = cos(qJ(4));
t688 = sin(pkin(6));
t696 = cos(qJ(2));
t950 = t688 * t696;
t591 = t689 * t692 + t695 * t950;
t1033 = t591 / 0.2e1;
t1017 = -t692 / 0.2e1;
t1052 = -pkin(10) - pkin(9);
t694 = cos(qJ(5));
t634 = t1052 * t694;
t690 = sin(qJ(6));
t1012 = cos(qJ(6));
t691 = sin(qJ(5));
t850 = t1012 * t691;
t807 = -t1052 * t850 - t690 * t634;
t1073 = t807 / 0.2e1;
t693 = sin(qJ(2));
t934 = t693 * t694;
t865 = t688 * t934;
t592 = t689 * t695 - t692 * t950;
t959 = t592 * t691;
t474 = -t865 + t959;
t951 = t688 * t693;
t866 = t691 * t951;
t958 = t592 * t694;
t475 = t866 + t958;
t808 = t1012 * t474 + t690 * t475;
t1049 = -t808 / 0.2e1;
t942 = t690 * t694;
t613 = t850 + t942;
t929 = t695 * t613;
t1038 = t929 / 0.2e1;
t873 = qJD(5) + qJD(6);
t1006 = t694 * pkin(4);
t793 = t1052 * t695 + qJ(3);
t697 = -pkin(2) - pkin(8);
t938 = t691 * t697;
t812 = pkin(5) - t938;
t701 = t793 * t694 + (t812 + t1006) * t692;
t394 = t1012 * t701;
t1009 = t692 * pkin(4);
t792 = -t695 * pkin(9) + t1009;
t759 = qJ(3) + t792;
t931 = t694 * t697;
t861 = t692 * t931;
t545 = t691 * t759 + t861;
t939 = t691 * t695;
t471 = -pkin(10) * t939 + t545;
t947 = t690 * t471;
t251 = -t394 + t947;
t935 = t692 * t697;
t862 = t691 * t935;
t470 = -t862 + (t793 + t1009) * t694;
t858 = t1012 * t470;
t272 = t858 - t947;
t1091 = t251 + t272;
t700 = t690 * t701;
t857 = t1012 * t471;
t252 = t857 + t700;
t948 = t690 * t470;
t270 = -t857 - t948;
t1090 = t252 + t270;
t283 = t1012 * t475 - t690 * t474;
t637 = t690 * t939;
t676 = t1012 * t694;
t579 = t676 * t695 - t637;
t1085 = t283 * t1017 + t1033 * t579;
t926 = t696 * t694;
t937 = t692 * t693;
t557 = (-t691 * t937 + t926) * t688;
t855 = t1012 * t557;
t927 = t696 * t691;
t558 = (t692 * t934 + t927) * t688;
t944 = t690 * t558;
t731 = -t944 / 0.2e1 + t855 / 0.2e1;
t1076 = t731 - t1085;
t1089 = qJD(1) * t1076;
t1086 = t1038 * t591 + t1049 * t692;
t854 = t1012 * t558;
t945 = t690 * t557;
t732 = -t945 / 0.2e1 - t854 / 0.2e1;
t1079 = t732 + t1086;
t1088 = qJD(1) * t1079;
t1029 = -t613 / 0.2e1;
t1035 = -t579 / 0.2e1;
t943 = t690 * t691;
t510 = -t1012 * t634 + t1052 * t943;
t1011 = pkin(5) * t691;
t811 = -t697 + t1011;
t603 = t811 * t695;
t679 = -pkin(5) * t694 - pkin(4);
t1087 = t1017 * t510 - t1029 * t603 - t1035 * t679;
t608 = -t676 + t943;
t1031 = -t608 / 0.2e1;
t1037 = -t929 / 0.2e1;
t1075 = t1031 * t603 + t679 * t1037 + t1073 * t692;
t730 = t942 / 0.2e1 + t850 / 0.2e1;
t716 = t1029 + t730;
t237 = t716 * t591;
t1084 = qJD(2) * t1076 + t237 * qJD(4);
t1030 = t608 / 0.2e1;
t1058 = t676 / 0.2e1 - t943 / 0.2e1;
t238 = (t1030 + t1058) * t591;
t1083 = qJD(2) * t1079 + t238 * qJD(4);
t1077 = t731 + t1085;
t1028 = t613 / 0.2e1;
t1057 = t1028 + t730;
t240 = t1057 * t591;
t1082 = qJD(2) * t1077 + t240 * qJD(4) - t283 * t873;
t1078 = t732 - t1086;
t239 = (t1031 + t1058) * t591;
t1081 = qJD(2) * t1078 + t239 * qJD(4) + t808 * t873;
t955 = t613 * t929;
t964 = t579 * t608;
t285 = -t955 / 0.2e1 - t964 / 0.2e1;
t1080 = t873 * t285;
t1005 = t695 * pkin(4);
t1007 = t692 * pkin(9);
t633 = t1005 + t1007;
t615 = t691 * t633;
t928 = t695 * t697;
t640 = t694 * t928;
t555 = t640 + t615;
t940 = t691 * t692;
t482 = pkin(10) * t940 + t555;
t856 = t1012 * t482;
t795 = -t856 / 0.2e1;
t617 = t694 * t633;
t936 = t692 * t694;
t441 = pkin(10) * t936 + t695 * t812 + t617;
t949 = t690 * t441;
t734 = -t949 / 0.2e1 + t795;
t93 = t734 - t1075;
t859 = t1012 * t441;
t946 = t690 * t482;
t733 = -t946 / 0.2e1 + t859 / 0.2e1;
t95 = t733 + t1087;
t92 = t733 - t1087;
t1046 = t283 / 0.2e1;
t1074 = -t510 / 0.2e1;
t903 = qJD(4) * t613;
t1068 = -qJD(2) * t285 + t608 * t903;
t910 = qJD(2) * t579;
t1067 = qJD(4) * t285 - t929 * t910;
t1026 = t637 / 0.2e1;
t794 = -t676 / 0.2e1;
t406 = t1026 + (t794 + t1031) * t695;
t920 = t238 * qJD(1) + t406 * qJD(3);
t1066 = t807 * t873 - t920;
t404 = t716 * t695;
t921 = t237 * qJD(1) - t404 * qJD(3);
t1063 = -t510 * t873 - t921;
t1062 = 0.2e1 * t691;
t1050 = t252 / 0.2e1;
t805 = t873 * t613;
t1061 = t608 * t805;
t685 = t692 ^ 2;
t687 = t695 ^ 2;
t1059 = t685 + t687;
t1056 = t794 + t1030;
t684 = t691 ^ 2;
t686 = t694 ^ 2;
t655 = t686 - t684;
t932 = t694 * t695;
t809 = t932 * t1062;
t727 = qJD(2) * t809 - qJD(4) * t655;
t1051 = t251 / 0.2e1;
t1048 = t808 / 0.2e1;
t1047 = -t283 / 0.2e1;
t1045 = t474 / 0.2e1;
t1044 = t475 / 0.2e1;
t575 = t613 * t692;
t1039 = -t575 / 0.2e1;
t636 = t690 * t940;
t578 = t676 * t692 - t636;
t1036 = t578 / 0.2e1;
t1034 = t579 / 0.2e1;
t1032 = t592 / 0.2e1;
t1027 = t636 / 0.2e1;
t677 = t687 * t694;
t1024 = -t677 / 0.2e1;
t1023 = t677 / 0.2e1;
t1022 = t685 / 0.2e1;
t1021 = -t686 / 0.2e1;
t1020 = t690 / 0.2e1;
t1019 = -t691 / 0.2e1;
t1018 = t691 / 0.2e1;
t1016 = t692 / 0.2e1;
t1015 = -t694 / 0.2e1;
t1014 = t694 / 0.2e1;
t1013 = -t695 / 0.2e1;
t1010 = t690 * pkin(5);
t1008 = t692 * pkin(5);
t1004 = t695 * pkin(5);
t818 = t1047 + t1046;
t819 = t1049 + t1048;
t37 = t579 * t818 - t819 * t929;
t1003 = t37 * qJD(5);
t989 = t808 * t608;
t215 = t989 / 0.2e1;
t834 = -t989 / 0.2e1;
t42 = t215 + t834;
t43 = -t608 * t819 + t613 * t818;
t1002 = t43 * qJD(5) + t42 * qJD(6);
t1000 = pkin(5) * qJD(5);
t566 = t592 * t1013;
t960 = t591 * t692;
t982 = t475 * t694;
t984 = t474 * t691;
t77 = t566 + (t982 / 0.2e1 + t984 / 0.2e1) * t695 + (t1021 - t684 / 0.2e1 + 0.1e1 / 0.2e1) * t960;
t999 = t77 * qJD(4);
t995 = qJD(4) * t42;
t994 = t251 * t578;
t993 = t251 * t608;
t992 = t252 * t575;
t991 = t252 * t613;
t986 = t283 * t613;
t983 = t475 * t692;
t980 = t807 * t929;
t979 = t807 * t578;
t976 = t510 * t575;
t744 = t694 * t759;
t544 = -t744 + t862;
t974 = t544 * t695;
t973 = t545 * t695;
t327 = t855 - t944;
t328 = t854 + t945;
t864 = t695 * t951;
t511 = t591 * t864;
t55 = t283 * t328 - t327 * t808 - t511;
t972 = t55 * qJD(1);
t971 = t557 * t692;
t970 = t575 * t690;
t969 = t575 * t692;
t968 = t929 * t695;
t967 = t578 * t692;
t966 = t579 * t510;
t965 = t579 * t603;
t963 = t579 * t690;
t364 = t591 * t691;
t366 = t591 * t694;
t954 = t613 * t690;
t941 = t691 * t929;
t933 = t694 * t685;
t930 = t695 * t579;
t678 = t695 * t692;
t860 = t691 * t928;
t554 = t617 - t860;
t133 = -t692 * t928 + (t555 * t1016 + t973 / 0.2e1) * t694 + (t554 * t1017 + t974 / 0.2e1) * t691;
t924 = t133 * qJD(4);
t167 = -t613 * t579 + t608 * t929;
t923 = t873 * t167;
t914 = t684 + t686;
t654 = t685 - t687;
t293 = -t955 + t964;
t913 = qJD(2) * t293;
t911 = qJD(2) * t929;
t612 = t654 * t691;
t909 = qJD(2) * t612;
t614 = -t677 + t933;
t908 = qJD(2) * t614;
t907 = qJD(2) * t688;
t905 = qJD(3) * t692;
t904 = qJD(4) * t608;
t902 = qJD(4) * t679;
t901 = qJD(4) * t691;
t900 = qJD(4) * t694;
t899 = qJD(5) * t691;
t898 = qJD(5) * t694;
t897 = qJD(6) * t679;
t822 = t932 / 0.2e1;
t824 = t935 / 0.2e1;
t101 = (-t1008 + t691 * t824 - t744 / 0.2e1 + pkin(10) * t822 + t470 / 0.2e1) * t690;
t896 = t101 * qJD(2);
t868 = -t1008 / 0.2e1;
t745 = -t394 / 0.2e1 + t1012 * t868;
t103 = t858 / 0.2e1 + t745;
t895 = t103 * qJD(2);
t131 = -t474 * t557 + t475 * t558 - t511;
t894 = t131 * qJD(1);
t248 = t575 * t579 + t578 * t929;
t891 = t248 * qJD(2);
t305 = (-t591 * t695 + t592 * t692 + t950) * t951;
t890 = t305 * qJD(1);
t749 = -t967 / 0.2e1 - t930 / 0.2e1;
t308 = -t749 + t1058;
t889 = t308 * qJD(2);
t750 = t969 / 0.2e1 + t968 / 0.2e1;
t309 = -t730 - t750;
t888 = t309 * qJD(2);
t357 = -t968 + t969;
t887 = t357 * qJD(2);
t358 = t930 - t967;
t886 = t358 * qJD(2);
t403 = t1057 * t692;
t375 = t403 * qJD(2);
t405 = t1056 * t692 + t1027;
t377 = t405 * qJD(2);
t815 = t1022 + t687 / 0.2e1;
t580 = (-0.1e1 / 0.2e1 - t815) * t691;
t883 = t580 * qJD(2);
t581 = t1023 + (0.1e1 / 0.2e1 + t1022) * t694;
t882 = t581 * qJD(2);
t593 = (t684 / 0.2e1 + t1021) * t695;
t881 = t593 * qJD(5);
t611 = t914 * t695;
t880 = t611 * qJD(2);
t879 = t654 * qJD(2);
t878 = t692 * qJD(2);
t877 = t692 * qJD(4);
t876 = t695 * qJD(2);
t875 = t695 * qJD(4);
t874 = t695 * qJD(5);
t38 = t575 * t818 + t578 * t819;
t351 = t591 * t613;
t352 = t591 * t608;
t39 = t352 * t1036 + t283 * t1034 + t351 * t1039 + t808 * t1038 + t566 + t960 / 0.2e1;
t872 = t39 * qJD(4) + t38 * qJD(5);
t871 = pkin(5) * t932;
t870 = t1012 / 0.2e1;
t869 = t1011 / 0.2e1;
t867 = -t1004 / 0.2e1;
t863 = t697 * t951;
t853 = t1012 * t929;
t852 = t1012 * t578;
t851 = t1012 * t608;
t849 = qJ(3) * t878;
t848 = qJ(3) * t876;
t845 = t691 * t900;
t844 = t691 * t875;
t843 = t694 * t875;
t842 = t692 * t899;
t841 = t691 * t874;
t840 = t692 * t898;
t839 = t694 * t874;
t838 = t608 * t878;
t837 = t613 * t878;
t652 = t693 * t907;
t836 = t696 * t907;
t835 = t691 * t898;
t670 = t692 * t875;
t669 = t692 * t876;
t833 = t474 * t1017;
t832 = t474 * t1015;
t831 = t558 * t1014;
t830 = -t364 / 0.2e1;
t829 = -t951 / 0.2e1;
t828 = t951 / 0.2e1;
t826 = t939 / 0.2e1;
t825 = -t937 / 0.2e1;
t823 = -t932 / 0.2e1;
t821 = t1050 + t270 / 0.2e1;
t820 = t272 / 0.2e1 + t1051;
t817 = -t807 / 0.2e1 + t1073;
t816 = t510 / 0.2e1 + t1074;
t814 = t1012 * qJD(5);
t813 = t1012 * qJD(6);
t810 = t914 * t591;
t491 = (-0.1e1 / 0.2e1 + t815) * t951;
t682 = qJD(2) * qJ(3);
t806 = qJD(1) * t491 - t682;
t804 = t873 * t692;
t803 = pkin(5) * t822;
t802 = -qJD(5) - t878;
t800 = t687 * t835;
t799 = t591 * t822;
t798 = t687 * t828;
t797 = t695 * t829;
t796 = t695 * t828;
t791 = qJD(4) * t809;
t266 = t859 - t946;
t267 = t856 + t949;
t602 = t811 * t692;
t698 = t266 * t1049 + t267 * t1046 - t351 * t251 / 0.2e1 + t352 * t1050 - t602 * t1033 + t603 * t1032;
t708 = t1073 * t327 + t1074 * t328 + t679 * t796;
t2 = t698 + t708;
t44 = -t251 * t266 + t252 * t267 - t602 * t603;
t788 = t2 * qJD(1) + t44 * qJD(2);
t706 = t272 * t1047 + t270 * t1048 + t1050 * t808 - t1051 * t283;
t738 = t1020 * t328 + t327 * t870;
t3 = (t591 * t823 + t738) * pkin(5) + t706;
t49 = -t251 * t270 + t252 * t272 + t603 * t871;
t787 = -t3 * qJD(1) + t49 * qJD(2);
t786 = -qJD(6) + t802;
t785 = qJD(1) * t42;
t704 = t1035 * t351 + t1037 * t352 + t1046 * t575 + t1049 * t578;
t755 = t1029 * t327 + t1031 * t328;
t18 = t704 - t755;
t28 = -t266 * t579 - t267 * t929 + t992 - t994;
t784 = t18 * qJD(1) + t28 * qJD(2);
t735 = t954 / 0.2e1 - t851 / 0.2e1;
t19 = -t820 * t578 + t821 * t575 + (t1023 + t735) * pkin(5);
t783 = t38 * qJD(1) - t19 * qJD(2);
t29 = -t1090 * t579 - t1091 * t929;
t782 = t37 * qJD(1) + t29 * qJD(2);
t711 = t1036 * t328 + t1039 * t327 + t798;
t45 = -t986 / 0.2e1 + t834 + t711;
t75 = t991 + t993;
t781 = -t45 * qJD(1) + t75 * qJD(2);
t715 = t1017 * t351 - t1032 * t929 + t1033 * t575;
t61 = (t608 * t829 + t1048) * t695 + t715;
t65 = -t251 * t695 + t266 * t692 - t575 * t603 - t602 * t929;
t780 = -t61 * qJD(1) + t65 * qJD(2);
t714 = t1016 * t352 + t1033 * t578 + t1035 * t592;
t62 = (t613 * t829 + t1046) * t695 + t714;
t66 = t252 * t695 + t267 * t692 + t603 * t578 + t602 * t579;
t779 = -t62 * qJD(1) - t66 * qJD(2);
t778 = t38 * qJD(3);
t56 = t283 * t352 - t351 * t808 + t591 * t592;
t777 = t56 * qJD(1) + t39 * qJD(3);
t775 = qJD(2) * t37 + qJD(4) * t43;
t774 = -t554 * t691 + t555 * t694;
t129 = t270 * t692 + t871 * t929 + t965;
t773 = -qJD(2) * t129 + t1089;
t454 = t603 * t929;
t130 = t272 * t692 - t579 * t871 + t454;
t772 = qJD(2) * t130 + t1088;
t145 = -t251 * t692 + t454;
t771 = qJD(2) * t145 + t1088;
t146 = -t252 * t692 + t965;
t770 = -qJD(2) * t146 + t1089;
t127 = (t592 - t982 - t984) * t591;
t769 = t127 * qJD(1) + t77 * qJD(3);
t395 = t1011 * t608 + t613 * t679;
t70 = (-t941 / 0.2e1 + (t1015 * t608 + t870) * t695) * pkin(5) + t92;
t768 = qJD(2) * t70 - qJD(4) * t395;
t396 = t1011 * t613 - t608 * t679;
t69 = (t579 * t1019 + (-t690 / 0.2e1 + t613 * t1015) * t695) * pkin(5) + t93;
t767 = qJD(2) * t69 - qJD(4) * t396;
t766 = t802 * t695;
t751 = t1019 * t557 + t831;
t753 = t1018 * t475 + t832;
t121 = t692 * t751 - t753 + t798;
t322 = -t544 * t694 + t545 * t691;
t765 = -qJD(1) * t121 + qJD(2) * t322;
t135 = (t833 - t558 / 0.2e1) * t694 + (t983 / 0.2e1 + t557 / 0.2e1) * t691;
t144 = t544 * t936 - t545 * t940 + (t554 * t694 + t555 * t691) * t695;
t764 = t135 * qJD(1) - t144 * qJD(2);
t140 = (t691 * t829 + t1044 - t958 / 0.2e1) * t695;
t295 = t973 + (t555 - 0.2e1 * t640) * t692;
t763 = -t140 * qJD(1) - t295 * qJD(2);
t141 = (t694 * t828 + t1045 - t959 / 0.2e1) * t695;
t294 = -t974 + (t554 + 0.2e1 * t860) * t692;
t762 = -t141 * qJD(1) + t294 * qJD(2);
t718 = (t691 * t825 + t926 / 0.2e1) * t688;
t740 = -t983 / 0.2e1 + t799;
t221 = t718 - t740;
t419 = -t545 * t692 - t687 * t931;
t761 = qJD(1) * t221 - qJD(2) * t419;
t717 = (t694 * t825 - t927 / 0.2e1) * t688;
t741 = t591 * t826 + t833;
t222 = t717 + t741;
t418 = -t544 * t692 - t687 * t938;
t760 = qJD(1) * t222 + qJD(2) * t418;
t314 = -t579 ^ 2 + t929 ^ 2;
t80 = qJD(2) * t314 + qJD(4) * t167;
t370 = t608 ^ 2 - t613 ^ 2;
t107 = qJD(2) * t167 + qJD(4) * t370;
t758 = t1007 / 0.2e1 + t1005 / 0.2e1;
t742 = t758 * t691;
t472 = t615 / 0.2e1 + t742;
t757 = pkin(4) * t900 - qJD(2) * t472;
t743 = t758 * t694;
t473 = -t617 / 0.2e1 - t743;
t756 = pkin(4) * t901 - qJD(2) * t473;
t754 = t554 * t1045 - t475 * t555 / 0.2e1;
t752 = t1028 * t510 + t1030 * t807;
t748 = qJD(2) * t93 + t608 * t902;
t747 = qJD(2) * t92 - t613 * t902;
t746 = t694 * t766;
t514 = -qJD(2) * t593 + t845;
t739 = t1020 * t267 + t266 * t870;
t737 = t1020 * t352 + t351 * t870;
t736 = t963 / 0.2e1 - t853 / 0.2e1;
t495 = qJD(2) * t677 * t691 + qJD(4) * t593;
t610 = t655 * t687;
t728 = qJD(2) * t610 + t791;
t726 = t751 * pkin(9);
t719 = (t970 / 0.2e1 + t852 / 0.2e1) * pkin(5);
t7 = t579 * t816 + t608 * t820 + t613 * t821 + t817 * t929 + t719;
t725 = t43 * qJD(1) - t7 * qJD(2);
t699 = -t1013 * t602 + t1016 * t603 + t1034 * t252 + t1036 * t267 + t1038 * t251 + t1039 * t266;
t15 = -t699 + t752;
t220 = t575 * t929 + t578 * t579 - t678;
t724 = t39 * qJD(1) - t15 * qJD(2) + t220 * qJD(3);
t154 = -t678 * t697 ^ 2 - t544 * t554 + t545 * t555;
t47 = (pkin(4) * t828 + t1032 * t697) * t695 + t726 + (t545 * t1014 + t544 * t1018 - t935 / 0.2e1) * t591 + t754;
t723 = -t47 * qJD(1) + t154 * qJD(2) + t133 * qJD(3);
t516 = t611 * t692 - t678;
t722 = t77 * qJD(1) + t133 * qJD(2) + t516 * qJD(3);
t13 = (t830 + t737) * pkin(5);
t139 = t1011 * t679;
t707 = t1090 * t1073 + t1074 * t1091;
t5 = (t1019 * t603 + t679 * t823 + t739) * pkin(5) + t707;
t57 = -t817 * t578 + t816 * t575 + (t826 + t736) * pkin(5);
t712 = -t13 * qJD(1) - t5 * qJD(2) - t57 * qJD(3) + t139 * qJD(4);
t683 = qJ(3) * qJD(3);
t675 = -t876 / 0.2e1;
t674 = t876 / 0.2e1;
t673 = t875 / 0.2e1;
t668 = t694 * t878;
t667 = t691 * t877;
t666 = t691 * t878;
t653 = qJ(3) * t950;
t619 = t687 * t863;
t601 = t669 + t874 / 0.2e1;
t583 = -t933 / 0.2e1 + t1024 + t1014;
t582 = t1018 * t1059 + t1019;
t568 = t669 + (qJD(5) / 0.2e1 + qJD(6) / 0.2e1) * t695;
t492 = t815 * t951 + t828;
t410 = t692 * t730 + t1039;
t409 = t1056 * t695 + t1026;
t408 = t1017 * t608 + t692 * t794 + t1027;
t407 = -t730 * t695 + t1037;
t391 = -t860 + t617 / 0.2e1 - t743;
t390 = -t640 - t615 / 0.2e1 + t742;
t311 = t749 + t1058;
t310 = -t730 + t750;
t303 = t311 * qJD(3);
t302 = t310 * qJD(3);
t301 = t309 * qJD(3);
t300 = t308 * qJD(3);
t299 = qJD(4) * t403 + t579 * t878;
t298 = qJD(4) * t405 + t878 * t929;
t269 = -t805 - t375;
t268 = -t608 * t873 - t377;
t224 = t718 + t740;
t223 = t717 - t741;
t158 = qJD(4) * t406 + t888;
t157 = -qJD(4) * t404 + t889;
t153 = t410 * qJD(4) + t579 * t786;
t152 = t408 * qJD(4) + t786 * t929;
t143 = t1013 * t475 + t592 * t822 + t691 * t797;
t142 = t1013 * t474 + t592 * t826 + t694 * t796;
t134 = t692 * t753 + t751;
t122 = t832 + t692 * t831 + t798 + (t1044 - t971 / 0.2e1) * t691;
t110 = t407 * qJD(4) - t578 * t873 - t889;
t109 = t409 * qJD(4) + t575 * t873 - t888;
t104 = t947 - t858 / 0.2e1 + t745;
t102 = t690 * t868 - t857 - t700 / 0.2e1 - t948 / 0.2e1;
t94 = t734 + t1075;
t72 = t579 * t869 + t613 * t803 + t795 + (t867 - t441 / 0.2e1) * t690 + t1075;
t71 = pkin(5) * t941 / 0.2e1 + t608 * t803 + t870 * t1004 + t95;
t64 = t1013 * t283 + t613 * t797 - t714;
t63 = t1013 * t808 + t608 * t797 - t715;
t58 = pkin(5) * t736 - t1036 * t807 - t1039 * t510 + t691 * t867 - t976 / 0.2e1 + t979 / 0.2e1;
t48 = -t545 * t366 / 0.2e1 + t544 * t830 + t697 * t566 + t591 * t824 + pkin(4) * t796 + t726 - t754;
t46 = t986 / 0.2e1 + t215 + t711;
t20 = t272 * t1036 + t270 * t1039 + t994 / 0.2e1 - t992 / 0.2e1 + (t1024 + t735) * pkin(5);
t17 = t704 + t755;
t16 = t699 + t752;
t14 = pkin(5) * t737 + t591 * t869;
t8 = t1029 * t270 + t1031 * t272 - t1035 * t510 - t1037 * t807 + t719 - t993 / 0.2e1 - t991 / 0.2e1 - t980 / 0.2e1 - t966 / 0.2e1;
t6 = pkin(5) * t739 + t603 * t869 + t679 * t803 - t707;
t4 = -t706 + (t738 + t799) * pkin(5);
t1 = t698 - t708;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t305, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t131 + qJD(4) * t127, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t55 + qJD(4) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t652, -t836, 0, 0, 0, 0, 0, 0, 0, 0, 0, t652, t836 (-pkin(2) * t951 + t653) * qJD(2) + qJD(3) * t951, 0, 0, 0, 0, 0, 0 (t693 * t875 + t696 * t878) * t688 (-t693 * t877 + t696 * t876) * t688, -t1059 * t652, t890 + (t685 * t863 + t619 + t653) * qJD(2) + t492 * qJD(3), 0, 0, 0, 0, 0, 0 (-t687 * t866 + t971) * qJD(2) + t142 * qJD(4) + t224 * qJD(5) (-t558 * t692 - t687 * t865) * qJD(2) + t143 * qJD(4) + t223 * qJD(5), t134 * qJD(4) + (-t557 * t694 - t558 * t691) * t876, t894 + (-t544 * t557 + t545 * t558 + t619) * qJD(2) + t122 * qJD(3) + t48 * qJD(4), 0, 0, 0, 0, 0, 0 (t327 * t692 - t864 * t929) * qJD(2) + t63 * qJD(4) + t873 * t1077 (-t328 * t692 - t579 * t864) * qJD(2) + t64 * qJD(4) + t873 * t1078 (-t327 * t579 - t328 * t929) * qJD(2) + t17 * qJD(4) + t1003, t972 + (-t251 * t327 + t252 * t328 - t603 * t864) * qJD(2) + t46 * qJD(3) + t1 * qJD(4) + t4 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t652, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t492, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t122 + t999, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t46 + t872; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t592 + t652 * t695, qJD(4) * t591 - t652 * t692, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t142 + qJD(5) * t364 - t592 * t900, qJD(2) * t143 + qJD(5) * t366 + t592 * t901, t134 * qJD(2) - qJD(4) * t810, t48 * qJD(2) + (-t592 * pkin(4) - pkin(9) * t810) * qJD(4) + t769, 0, 0, 0, 0, 0, 0, t63 * qJD(2) + t240 * t873 + t592 * t904, t64 * qJD(2) + t239 * t873 + t592 * t903, t17 * qJD(2) + (-t351 * t613 - t352 * t608) * qJD(4) + t1002, t1 * qJD(2) + (-t351 * t807 + t352 * t510 + t592 * t679) * qJD(4) + t14 * qJD(5) + t777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t224 + qJD(4) * t364 - qJD(5) * t475, qJD(2) * t223 + qJD(4) * t366 + qJD(5) * t474, 0, 0, 0, 0, 0, 0, 0, 0, t1082, t1081, t775, t4 * qJD(2) + t14 * qJD(4) + (-t1012 * t283 - t690 * t808) * t1000 + t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1082, t1081, t995, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t491 - t890, 0, 0, 0, 0, 0, 0, -qJD(4) * t141 - qJD(5) * t221, -qJD(4) * t140 - qJD(5) * t222, qJD(4) * t135, -qJD(3) * t121 - qJD(4) * t47 - t894, 0, 0, 0, 0, 0, 0, -qJD(4) * t61 - t1076 * t873, -qJD(4) * t62 - t1079 * t873, qJD(4) * t18 + t1003, -qJD(3) * t45 + qJD(4) * t2 - qJD(5) * t3 - t972; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t683, -t670, t654 * qJD(4), 0, t670, 0, 0, qJ(3) * t875 + t905, -qJ(3) * t877 + qJD(3) * t695, 0, t683, -t670 * t686 - t800, -qJD(5) * t610 + t692 * t791, -qJD(4) * t614 - t692 * t841, -t670 * t684 + t800, qJD(4) * t612 - t692 * t839, t670, qJD(4) * t294 + qJD(5) * t419 + t694 * t905, -qJD(4) * t295 - qJD(5) * t418 - t691 * t905, -qJD(3) * t611 - qJD(4) * t144, qJD(3) * t322 + qJD(4) * t154 (-qJD(4) * t578 - t873 * t929) * t579, t248 * qJD(4) + t314 * t873, t358 * qJD(4) - t804 * t929 (-qJD(4) * t575 + t579 * t873) * t929, t357 * qJD(4) - t579 * t804, t670, qJD(4) * t65 + qJD(5) * t129 + qJD(6) * t146 - t608 * t905, -qJD(4) * t66 - qJD(5) * t130 - qJD(6) * t145 - t613 * t905, qJD(3) * t293 + qJD(4) * t28 + qJD(5) * t29, qJD(3) * t75 + qJD(4) * t44 + qJD(5) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t682, 0, 0, 0, 0, 0, 0, t878, t876, 0, -t806, 0, 0, 0, 0, 0, 0, qJD(5) * t583 + t668, qJD(5) * t582 - t666, -t880, t765 + t924, 0, 0, 0, 0, 0, 0, t311 * t873 - t838, t310 * t873 - t837, t913 (t575 * t608 + t613 * t578) * qJD(3) + t16 * qJD(4) + t20 * qJD(5) + t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t669, t879, -t877, t669, -t875, 0, -t697 * t877 + t848, -t697 * t875 - t849, 0, 0, -t881 + (-t686 * t876 - t845) * t692, t692 * t727 - 0.2e1 * t695 * t835, t844 - t908, t881 + (-t684 * t876 + t845) * t692, t843 + t909, t601 (t691 * t792 - t861) * qJD(4) + t391 * qJD(5) + t762 (-pkin(9) * t932 + (t938 + t1006) * t692) * qJD(4) + t390 * qJD(5) + t763, qJD(4) * t774 + t764 (-pkin(4) * t935 + pkin(9) * t774) * qJD(4) + t723 (-t903 - t910) * t578 + t1080, t891 + (t575 * t613 + t578 * t608) * qJD(4) + t923, t408 * t873 + t613 * t875 + t886 (-t904 - t911) * t575 - t1080, t410 * t873 - t608 * t875 + t887, t568 (-t575 * t679 - t602 * t608 - t695 * t807) * qJD(4) + t71 * qJD(5) + t95 * qJD(6) + t780 (-t510 * t695 - t578 * t679 - t602 * t613) * qJD(4) + t72 * qJD(5) + t94 * qJD(6) + t779 (-t266 * t613 - t267 * t608 + t976 - t979) * qJD(4) + t8 * qJD(5) + t784, t16 * qJD(3) + (-t266 * t807 + t267 * t510 - t602 * t679) * qJD(4) + t6 * qJD(5) + t788; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t495, -t728, t691 * t766, t495, t746, t673, qJD(3) * t583 + qJD(4) * t391 - qJD(5) * t545 - t761, qJD(3) * t582 + qJD(4) * t390 + qJD(5) * t544 - t760, 0, 0, t1067, t80, t152, -t1067, t153, t673, qJD(4) * t71 + qJD(5) * t270 + qJD(6) * t102 + t303 - t773, qJD(4) * t72 - qJD(5) * t272 + qJD(6) * t104 + t302 - t772, t8 * qJD(4) + (t853 - t963) * t1000 + t782, t20 * qJD(3) + t6 * qJD(4) + (t1012 * t270 + t272 * t690) * t1000 + t787; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1067, t80, t152, -t1067, t153, t673, qJD(4) * t95 + qJD(5) * t102 - qJD(6) * t252 + t303 - t770, qJD(4) * t94 + qJD(5) * t104 + qJD(6) * t251 + t302 - t771, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t491, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t121 + t999, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t45 + t872; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t682, 0, 0, 0, 0, 0, 0, -t878, -t876, 0, t806, 0, 0, 0, 0, 0, 0, -qJD(5) * t581 - t668, -qJD(5) * t580 + t666, t880, -t765 + t924, 0, 0, 0, 0, 0, 0, -t308 * t873 + t838, -t309 * t873 + t837, -t913, -qJD(4) * t15 - qJD(5) * t19 - t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t516 * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t877, -t875, 0, 0, 0, 0, 0, 0, 0, 0, -t694 * t877 - t841, t667 - t839, t611 * qJD(4) (pkin(9) * t611 - t1009) * qJD(4) + t722, 0, 0, 0, 0, 0, 0, t407 * t873 + t608 * t877, t409 * t873 + t613 * t877, -qJD(4) * t293 (t679 * t692 + t966 + t980) * qJD(4) + t58 * qJD(5) + t724; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t840 - t844 - t882, t842 - t843 - t883, 0, 0, 0, 0, 0, 0, 0, 0, t110, t109, 0, t58 * qJD(4) + (-t852 - t970) * t1000 + t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t109, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t141, qJD(2) * t140, -qJD(2) * t135, qJD(2) * t47 - t769, 0, 0, 0, 0, 0, 0, t61 * qJD(2) - t237 * t873, t62 * qJD(2) - t238 * t873, -qJD(2) * t18 + t1002, -qJD(2) * t2 - qJD(5) * t13 - t777; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t669, -t879, 0, -t669, 0, 0, -t848, t849, 0, 0, t669 * t686 - t881, t746 * t1062, t840 + t908, t669 * t684 + t881, -t842 - t909, -t601, qJD(5) * t473 - t762, qJD(5) * t472 - t763, -t764, -t723, t578 * t910 + t1080, -t891 + t923, -t405 * t873 - t886, t575 * t911 - t1080, -t403 * t873 - t887, -t568, -qJD(5) * t70 - qJD(6) * t92 - t780, -qJD(5) * t69 - qJD(6) * t93 - t779, -qJD(5) * t7 - t784, qJD(3) * t15 - qJD(5) * t5 - t788; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t722, 0, 0, 0, 0, 0, 0, t873 * t404, -t873 * t406, 0, -qJD(5) * t57 - t724; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t835, t655 * qJD(5), 0, -t835, 0, 0, -pkin(4) * t899, -pkin(4) * t898, 0, 0, -t1061, t873 * t370, 0, t1061, 0, 0, qJD(5) * t395 + t613 * t897, qJD(5) * t396 - t608 * t897, 0, qJD(5) * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t514, -t727, t668 + t898, -t514, -t666 - t899, t675, -pkin(9) * t898 - t756, pkin(9) * t899 - t757, 0, 0, -t1068, t107, t268, t1068, t269, t675, t1063 - t768, t1066 - t767 (t851 - t954) * t1000 + t725 (-t1012 * t510 - t690 * t807) * t1000 + t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1068, t107, t268, t1068, t269, t675, t1063 - t747, t1066 - t748, t785, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t221, qJD(2) * t222, 0, 0, 0, 0, 0, 0, 0, 0, t1084, t1083, -t775, qJD(2) * t3 + qJD(4) * t13 - t778; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t495, t728 (t691 * t876 - t900) * t692, -t495, t669 * t694 + t667, t673, qJD(3) * t581 - qJD(4) * t473 + t761, qJD(3) * t580 - qJD(4) * t472 + t760, 0, 0, -t1067, -t80, t298, t1067, t299, t673, qJD(4) * t70 + qJD(6) * t101 + t300 + t773, qJD(4) * t69 + qJD(6) * t103 + t301 + t772, qJD(4) * t7 - t782, qJD(3) * t19 + qJD(4) * t5 - t787; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t882, t883, 0, 0, 0, 0, 0, 0, 0, 0, t157, t158, 0, qJD(4) * t57 - t783; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t514, t727, -t668, t514, t666, t674, t756, t757, 0, 0, t1068, -t107, t377, -t1068, t375, t674, t768 + t921, t767 + t920, -t725, -t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t1010, -pkin(5) * t813, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1010 * t873 + t896, t895 + (-t814 - t813) * pkin(5), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1084, t1083, -t995, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1067, -t80, t298, t1067, t299, t673, qJD(4) * t92 - qJD(5) * t101 + t300 + t770, qJD(4) * t93 - qJD(5) * t103 + t301 + t771, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t158, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1068, -t107, t377, -t1068, t375, t674, t747 + t921, t748 + t920, -t785, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1000 * t690 - t896, pkin(5) * t814 - t895, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t9;
