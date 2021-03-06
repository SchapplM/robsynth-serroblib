% Calculate inertial parameters regressor of coriolis matrix for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PPRRRR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:37
% EndTime: 2019-03-08 19:02:04
% DurationCPUTime: 22.77s
% Computational Cost: add. (21749->727), mult. (58150->1081), div. (0->0), fcn. (71236->14), ass. (0->565)
t927 = qJD(4) + qJD(5);
t1046 = cos(qJ(3));
t1028 = cos(pkin(7));
t1029 = cos(pkin(6));
t713 = sin(pkin(7));
t1027 = sin(pkin(6));
t828 = cos(pkin(13)) * t1027;
t1084 = t1028 * t828 + t1029 * t713;
t717 = sin(qJ(3));
t868 = sin(pkin(13)) * t1027;
t581 = -t1046 * t1084 + t717 * t868;
t718 = cos(qJ(6));
t1045 = cos(qJ(5));
t582 = t1046 * t868 + t1084 * t717;
t716 = sin(qJ(4));
t719 = cos(qJ(4));
t743 = t1028 * t1029 - t713 * t828;
t458 = t582 * t719 + t716 * t743;
t441 = t1045 * t458;
t715 = sin(qJ(5));
t731 = t582 * t716 - t719 * t743;
t729 = t715 * t731;
t1086 = -t729 + t441;
t714 = sin(qJ(6));
t988 = t714 * t1086;
t200 = -t581 * t718 + t988;
t1025 = t200 * t714;
t862 = t1045 * t731 + t715 * t458;
t1101 = t1086 * t718;
t963 = t581 * t714 + t1101;
t865 = t963 * t718;
t1117 = (t1086 - t865 - t1025) * t862;
t1120 = qJD(1) * t1117;
t912 = t713 * t1046;
t989 = t713 * t717;
t664 = t1028 * t716 + t719 * t989;
t644 = t1045 * t664;
t776 = -t1028 * t719 + t716 * t989;
t752 = t715 * t776;
t1085 = -t752 + t644;
t986 = t714 * t1085;
t466 = t718 * t912 + t986;
t1010 = t466 * t714;
t1102 = t1085 * t718;
t805 = -t714 * t912 + t1102;
t785 = t805 * t718;
t861 = t1045 * t776 + t715 * t664;
t1118 = (t1085 - t785 - t1010) * t861;
t1119 = qJD(2) * t1118;
t710 = t718 ^ 2;
t1053 = -t710 / 0.2e1;
t1074 = -t1101 / 0.2e1;
t1067 = -t1102 / 0.2e1;
t908 = t1045 * t716;
t978 = t715 * t719;
t678 = t908 + t978;
t703 = -pkin(4) * t719 - pkin(3);
t952 = qJD(3) * t703;
t1058 = -t678 / 0.2e1;
t666 = t908 / 0.2e1 + t978 / 0.2e1;
t290 = (t1058 + t666) * t581;
t1057 = t678 / 0.2e1;
t849 = t1046 * t1045;
t807 = -t849 / 0.2e1;
t910 = t719 * t1046;
t853 = t715 * t910;
t745 = -t853 / 0.2e1 + t716 * t807;
t488 = (t1046 * t1057 + t745) * t713;
t967 = t290 * qJD(1) + t488 * qJD(2);
t1116 = -t678 * t952 + t967;
t1044 = pkin(4) * t716;
t907 = t1045 * t719;
t979 = t715 * t716;
t676 = -t907 + t979;
t552 = t1044 * t676 + t678 * t703;
t1115 = -qJD(3) * t552 + t967;
t1091 = pkin(10) + pkin(9);
t689 = t1091 * t719;
t860 = t1091 * t908 + t715 * t689;
t1106 = t860 * t714;
t1036 = t678 * pkin(5);
t1037 = t676 * pkin(11);
t600 = t1036 + t1037;
t975 = t718 * t600;
t424 = t975 + t1106;
t1105 = t860 * t718;
t984 = t714 * t600;
t425 = t984 - t1105;
t680 = t1045 * t689;
t855 = t1091 * t979;
t1089 = t680 - t855;
t1100 = t1089 * t714;
t1038 = t676 * pkin(5);
t839 = -t678 * pkin(11) + t1038;
t773 = t703 + t839;
t394 = -t718 * t773 + t1100;
t1099 = t1089 * t718;
t395 = t714 * t773 + t1099;
t789 = (t394 * t718 - t395 * t714) * t676;
t77 = (t424 * t718 + t425 * t714) * t678 + t789;
t1048 = -t718 / 0.2e1;
t1050 = t714 / 0.2e1;
t736 = (t1048 * t466 + t1050 * t805) * t676;
t911 = t716 * t1046;
t616 = (-t715 * t911 + t719 * t849) * t713;
t981 = t714 * t616;
t550 = t718 * t989 - t981;
t913 = t714 * t989;
t972 = t718 * t616;
t551 = t913 + t972;
t794 = t1048 * t551 + t1050 * t550;
t131 = t736 + t794;
t750 = (t1048 * t200 + t1050 * t963) * t676;
t1004 = t582 * t718;
t369 = t581 * t676;
t987 = t714 * t369;
t327 = -t987 + t1004;
t1005 = t582 * t714;
t977 = t718 * t369;
t328 = t977 + t1005;
t798 = t1048 * t328 + t1050 * t327;
t44 = t750 + t798;
t915 = t44 * qJD(1) + t131 * qJD(2);
t1114 = t77 * qJD(3) - t915;
t576 = t600 + t1044;
t976 = t718 * t576;
t416 = t976 + t1106;
t985 = t714 * t576;
t417 = -t1105 + t985;
t72 = (t416 * t718 + t417 * t714) * t678 + t789;
t1113 = t72 * qJD(3) - t915;
t1108 = t860 / 0.2e1;
t1097 = t1108 * t1085;
t1096 = t1108 * t1086;
t708 = t714 ^ 2;
t1054 = t708 / 0.2e1;
t1112 = t1053 - t1054;
t1111 = (-t394 + t1100) * t678;
t1110 = (-t395 + t1099) * t678;
t1109 = -0.2e1 * t678;
t1061 = -t1089 / 0.2e1;
t1107 = t1089 / 0.2e1;
t1041 = t1085 * pkin(5);
t1043 = t1086 * pkin(5);
t1040 = t1089 * pkin(5);
t873 = t1053 - t708 / 0.2e1;
t1104 = t873 * t861;
t1103 = t873 * t862;
t1003 = t1089 * t860;
t1090 = t927 * t678;
t1098 = t676 * t1090;
t1059 = t676 / 0.2e1;
t777 = t907 / 0.2e1 - t979 / 0.2e1;
t291 = (t1059 + t777) * t581;
t844 = t911 / 0.2e1;
t744 = t715 * t844 + t719 * t807;
t925 = -t1046 / 0.2e1;
t489 = (t676 * t925 + t744) * t713;
t966 = t291 * qJD(1) + t489 * qJD(2);
t1095 = t860 * t927 - t966;
t846 = t912 / 0.2e1;
t490 = t676 * t846 + t713 * t744;
t1094 = t490 * qJD(3) + t861 * t927;
t292 = (-t676 / 0.2e1 + t777) * t581;
t1093 = t292 * qJD(3) + t862 * t927;
t1092 = 0.2e1 * t678;
t596 = t927 * t676;
t961 = t708 + t710;
t709 = t716 ^ 2;
t711 = t719 ^ 2;
t1088 = t709 + t711;
t1087 = -t731 * t716 + t582;
t1083 = t1112 * t861;
t1082 = t1112 * t862;
t909 = t1045 * t676;
t923 = t1045 * pkin(4);
t702 = -t923 - pkin(5);
t990 = t702 * t676;
t1035 = t715 * pkin(4);
t701 = pkin(11) + t1035;
t991 = t701 * t678;
t993 = t678 * t715;
t1081 = t1038 / 0.2e1 - t991 / 0.2e1 - t990 / 0.2e1 + (-t909 / 0.2e1 + t993 / 0.2e1) * pkin(4);
t671 = t678 ^ 2;
t926 = -t676 ^ 2 + t671;
t568 = (t1054 + t1053) * t678;
t980 = t714 * t718;
t904 = qJD(3) * t980;
t338 = t568 * t927 + t671 * t904;
t694 = t710 - t708;
t531 = t1109 * t904 + t694 * t927;
t1075 = -t200 / 0.2e1;
t1073 = t862 / 0.2e1;
t1072 = t1086 / 0.2e1;
t1071 = t417 / 0.2e1;
t1070 = -t424 / 0.2e1;
t1069 = t425 / 0.2e1;
t843 = -t441 / 0.2e1;
t1068 = -t466 / 0.2e1;
t1066 = t861 / 0.2e1;
t1065 = t1085 / 0.2e1;
t842 = -t644 / 0.2e1;
t841 = -t680 / 0.2e1;
t1056 = -t702 / 0.2e1;
t1055 = t702 / 0.2e1;
t1051 = -t714 / 0.2e1;
t1049 = t715 / 0.2e1;
t1047 = t718 / 0.2e1;
t368 = t581 * t678;
t1042 = t368 * pkin(5);
t615 = (t716 * t849 + t853) * t713;
t1039 = t615 * pkin(5);
t43 = t750 - t798;
t1034 = t43 * qJD(3);
t878 = t1072 - t1086 / 0.2e1;
t50 = t678 * t878;
t1033 = t50 * qJD(4);
t1032 = t44 * qJD(3);
t1031 = pkin(4) * qJD(5);
t1030 = qJD(4) * pkin(4);
t1026 = t200 * t676;
t1022 = t862 * t368;
t1021 = t368 * t860;
t1020 = t368 * t702;
t1019 = t368 * t714;
t1018 = t368 * t718;
t1017 = t416 * t714;
t1016 = t417 * t718;
t1015 = t424 * t714;
t1014 = t425 * t718;
t1012 = t458 * t719;
t1011 = t466 * t676;
t1009 = t50 * qJD(3);
t1006 = t861 * t615;
t381 = t581 * t716;
t1000 = t860 * t615;
t997 = t615 * t702;
t996 = t615 * t714;
t995 = t615 * t718;
t994 = t678 * t714;
t992 = t678 * t718;
t130 = t736 - t794;
t969 = t130 * qJD(3);
t876 = t1065 - t1085 / 0.2e1;
t135 = t678 * t876;
t968 = t135 * qJD(4);
t962 = t131 * qJD(3);
t486 = t926 * t714;
t960 = qJD(3) * t486;
t487 = t926 * t718;
t959 = qJD(3) * t487;
t553 = t1044 * t678 - t676 * t703;
t957 = qJD(3) * t553;
t956 = qJD(3) * t568;
t579 = t694 * t671;
t955 = qJD(3) * t579;
t954 = qJD(3) * t676;
t953 = qJD(3) * t678;
t951 = qJD(3) * t713;
t950 = qJD(3) * t717;
t949 = qJD(3) * t719;
t948 = qJD(4) * t714;
t947 = qJD(4) * t718;
t946 = qJD(5) * t703;
t945 = qJD(5) * t714;
t944 = qJD(5) * t718;
t943 = qJD(6) * t714;
t942 = qJD(6) * t718;
t940 = t135 * qJD(3);
t935 = t926 * qJD(3);
t934 = t568 * qJD(6);
t569 = t714 * t676;
t933 = t569 * qJD(3);
t572 = t718 * t676;
t563 = t572 * qJD(3);
t932 = t582 * qJD(3);
t931 = t666 * qJD(3);
t695 = t711 - t709;
t930 = t695 * qJD(3);
t929 = t716 * qJD(4);
t928 = t719 * qJD(4);
t924 = -t1045 / 0.2e1;
t922 = pkin(3) * t949;
t921 = pkin(3) * t716 * qJD(3);
t920 = t715 * t1030;
t919 = t715 * t1031;
t917 = pkin(11) * t1051;
t916 = pkin(11) * t1047;
t906 = t676 * t952;
t903 = qJD(6) * t676 * t678;
t902 = t676 * t953;
t901 = t713 * t950;
t698 = t714 * t942;
t900 = t716 * t928;
t287 = t862 * t1050;
t288 = t862 * t1047;
t899 = t1019 / 0.2e1;
t898 = -t1018 / 0.2e1;
t512 = t861 * t1050;
t513 = t861 * t1047;
t897 = -t996 / 0.2e1;
t896 = t995 / 0.2e1;
t895 = -t994 / 0.2e1;
t894 = t992 / 0.2e1;
t889 = t989 / 0.2e1;
t888 = t988 / 0.2e1;
t887 = -t987 / 0.2e1;
t886 = t986 / 0.2e1;
t885 = -t981 / 0.2e1;
t884 = t701 * t1051;
t883 = -t977 / 0.2e1;
t882 = t395 * t1047;
t881 = -t972 / 0.2e1;
t880 = t701 * t1047;
t875 = t1108 - t860 / 0.2e1;
t874 = t1107 + t1061;
t872 = pkin(11) * t961;
t871 = qJD(4) * t1046;
t870 = t1045 * qJD(4);
t869 = t1045 * qJD(5);
t866 = t963 * t676;
t864 = t961 * t701;
t863 = t1088 * t581;
t859 = qJD(4) * t961;
t858 = qJD(5) * t961;
t856 = t927 * t714;
t854 = -t923 / 0.2e1;
t850 = t671 * t698;
t848 = qJD(3) * t912;
t847 = -t912 / 0.2e1;
t845 = t1045 * t1050;
t840 = t927 * t1035;
t832 = t873 * t701;
t831 = t963 * t1058;
t830 = t865 / 0.2e1;
t829 = t582 / 0.2e1 + t862 * t1058;
t826 = t927 * t980;
t824 = t718 * t856;
t758 = t785 / 0.2e1;
t746 = t758 + t1010 / 0.2e1;
t772 = t830 + t1025 / 0.2e1;
t9 = (t1072 - t772) * t861 + (t1065 - t746) * t862;
t823 = t9 * qJD(2) + t1120;
t822 = t9 * qJD(1) + t1119;
t801 = t1066 * t1086 + t1073 * t1085;
t10 = -t746 * t862 - t772 * t861 + t801;
t821 = t10 * qJD(2) + t1120;
t820 = t10 * qJD(1) + t1119;
t799 = -t1066 * t368 + t1073 * t615;
t11 = t328 * t805 / 0.2e1 + t963 * t551 / 0.2e1 + t327 * t1068 + t550 * t1075 + t799;
t19 = -t200 * t327 + t328 * t963 - t1022;
t819 = t19 * qJD(1) + t11 * qJD(2);
t33 = -t1065 * t862 - t1072 * t861 + t801;
t818 = t33 * qJD(2);
t36 = t369 * t1065 + t616 * t1072 + (t582 * t925 + t581 * t717 / 0.2e1) * t713 + t799;
t59 = t1086 * t369 + t581 * t582 - t1022;
t817 = t59 * qJD(1) + t36 * qJD(2);
t816 = t1016 - t1017;
t815 = t1014 - t1015;
t814 = -t990 - t991;
t112 = -t466 * t550 + t551 * t805 + t1006;
t813 = t11 * qJD(1) + t112 * qJD(2);
t111 = (-t1012 + t1087) * t581;
t751 = t716 * t776;
t80 = t846 * t1012 + (-t719 * t664 / 0.2e1 - t751 / 0.2e1 + t889) * t581 + t1087 * t847;
t812 = t111 * qJD(1) + t80 * qJD(2);
t811 = t33 * qJD(1);
t690 = t713 ^ 2 * t717 * t1046;
t230 = t1085 * t616 + t1006 - t690;
t810 = t36 * qJD(1) + t230 * qJD(2);
t439 = t664 * t713 * t910 + t751 * t912 - t690;
t809 = t80 * qJD(1) + t439 * qJD(2);
t808 = t678 * (-qJD(6) - t954);
t806 = t1056 + t854;
t804 = t1037 / 0.2e1 + t1036 / 0.2e1;
t803 = (t1075 + t888) * t678;
t802 = (t1068 + t886) * t678;
t800 = t1107 * t862 + t1096;
t797 = -t1016 / 0.2e1 + t1017 / 0.2e1;
t796 = t1014 / 0.2e1 - t1015 / 0.2e1;
t795 = t1107 * t861 + t1097;
t793 = t1056 * t678 + t1059 * t701;
t792 = t1058 * t861 + t889;
t791 = t1050 * t394 + t882;
t790 = t718 * t808;
t498 = qJD(6) * t666 + t902;
t787 = pkin(5) / 0.2e1 + t806;
t786 = t805 * t676;
t783 = t600 / 0.2e1 + t804;
t782 = t824 * t1092;
t781 = t1049 * t369 - t368 * t924;
t780 = t1049 * t616 + t615 * t924;
t779 = t1088 * t1046;
t778 = t961 * t1045;
t723 = t1068 * t416 + t1071 * t805 - t791 * t861 + t795;
t21 = -t997 / 0.2e1 + t794 * t701 + t723;
t725 = t1071 * t963 + t1075 * t416 - t791 * t862 + t800;
t4 = t1020 / 0.2e1 + t798 * t701 + t725;
t78 = -t394 * t416 + t395 * t417 + t1003;
t775 = t4 * qJD(1) + t21 * qJD(2) + t78 * qJD(3);
t770 = t1107 - t791;
t724 = t1068 * t424 + t1069 * t805 + t770 * t861 + t1097;
t23 = t1039 / 0.2e1 + t794 * pkin(11) + t724;
t726 = t1069 * t963 + t1070 * t200 + t770 * t862 + t1096;
t6 = -t1042 / 0.2e1 + t798 * pkin(11) + t726;
t79 = -t394 * t424 + t395 * t425 + t1003;
t774 = t6 * qJD(1) + t23 * qJD(2) + t79 * qJD(3);
t771 = t576 / 0.2e1 + t793;
t113 = t1111 + (t416 - t1106) * t676;
t757 = t1057 * t1086;
t735 = t1058 * t200 + t714 * t757;
t26 = t898 + t735;
t756 = t1057 * t1085;
t734 = t1058 * t466 + t714 * t756;
t91 = t896 + t734;
t769 = t26 * qJD(1) + t91 * qJD(2) + t113 * qJD(3);
t114 = t1110 + (-t417 - t1105) * t676;
t733 = t718 * t757 + t831;
t31 = t899 + t733;
t759 = t805 * t1058;
t727 = t718 * t756 + t759;
t96 = t897 + t727;
t768 = t31 * qJD(1) + t96 * qJD(2) + t114 * qJD(3);
t115 = t1111 + (t424 - t1106) * t676;
t24 = t898 + t803;
t89 = t896 + t802;
t767 = t24 * qJD(1) + t89 * qJD(2) + t115 * qJD(3);
t116 = t1110 + (-t425 - t1105) * t676;
t755 = t1086 * t894 + t831;
t29 = t899 + t755;
t740 = t1085 * t894 + t759;
t94 = t897 + t740;
t766 = t29 * qJD(1) + t94 * qJD(2) + t116 * qJD(3);
t739 = -t1061 * t862 + t1096 - t800;
t17 = (-t381 / 0.2e1 + t781) * pkin(4) + t739;
t329 = t1044 * t703;
t738 = -t1061 * t861 + t1097 - t795;
t70 = (t713 * t844 + t780) * pkin(4) + t738;
t765 = -t17 * qJD(1) - t70 * qJD(2) + t329 * qJD(3);
t764 = -t50 * qJD(1) - t135 * qJD(2);
t194 = t885 + t786 / 0.2e1 + t792 * t718;
t295 = -t395 * t676 + t860 * t992;
t62 = t887 + t866 / 0.2e1 + t829 * t718;
t763 = qJD(1) * t62 + qJD(2) * t194 - qJD(3) * t295;
t195 = t881 - t1011 / 0.2e1 - t792 * t714;
t294 = t394 * t676 - t860 * t994;
t63 = t883 - t1026 / 0.2e1 - t829 * t714;
t762 = qJD(1) * t63 + qJD(2) * t195 - qJD(3) * t294;
t105 = t1074 + t1101 / 0.2e1;
t730 = (-t701 / 0.2e1 + t1035 / 0.2e1 + pkin(11) / 0.2e1) * t678 + (-pkin(5) / 0.2e1 + t806) * t676;
t151 = t714 * t730 - t718 * t874;
t348 = t1067 + t1102 / 0.2e1;
t761 = qJD(1) * t105 + qJD(2) * t348 + qJD(3) * t151;
t305 = t843 + t441 / 0.2e1;
t527 = t842 + t644 / 0.2e1;
t609 = t841 + t680 / 0.2e1;
t760 = qJD(1) * t305 + qJD(2) * t527 + qJD(3) * t609;
t252 = t783 * t714;
t626 = t787 * t718;
t754 = pkin(5) * t944 - qJD(3) * t252 + qJD(4) * t626;
t254 = t783 * t718;
t625 = t787 * t714;
t753 = pkin(5) * t945 + qJD(3) * t254 + qJD(4) * t625;
t721 = t796 * t701 + (t1045 * t882 + t1049 * t860 + t394 * t845) * pkin(4) + t1089 * t1055;
t52 = t1040 / 0.2e1 + t797 * pkin(11) + t721;
t564 = (t701 * t778 + t702 * t715) * pkin(4);
t720 = (t1045 * t758 + t1049 * t861 + t466 * t845) * pkin(4) + t861 * t832 + t1085 * t1055;
t61 = t1041 / 0.2e1 - pkin(11) * t1104 + t720;
t722 = (t1045 * t830 + t1049 * t862 + t200 * t845) * pkin(4) + t862 * t832 + t1086 * t1055;
t8 = t1043 / 0.2e1 - pkin(11) * t1103 + t722;
t749 = t8 * qJD(1) + t61 * qJD(2) + t52 * qJD(3) + t564 * qJD(4);
t665 = t778 * pkin(4);
t83 = (t1069 - t417 / 0.2e1) * t718 + (t1070 + t416 / 0.2e1) * t714;
t747 = -qJD(3) * t83 - qJD(4) * t665;
t226 = t714 * t875 - t718 * t771;
t742 = -qJD(3) * t226 - t702 * t948;
t224 = t714 * t771 + t718 * t875;
t741 = -qJD(3) * t224 - t702 * t947;
t102 = t878 * t714;
t154 = t714 * t874 + t718 * t730;
t345 = t876 * t714;
t737 = -qJD(1) * t102 - qJD(2) * t345 - qJD(3) * t154 - t714 * t920;
t699 = t716 * t949;
t693 = t714 * t919;
t687 = t694 * qJD(6);
t661 = t665 * qJD(5);
t628 = pkin(5) * t1048 + t1047 * t702 + t718 * t854;
t627 = pkin(5) * t1051 + t1050 * t702 + t714 * t854;
t604 = -t995 / 0.2e1;
t603 = t996 / 0.2e1;
t577 = t927 * t666;
t533 = 0.2e1 * t841 + t855;
t523 = t563 + t942;
t522 = -t933 - t943;
t491 = t678 * t847 + t713 * t745;
t483 = t491 * qJD(3);
t478 = t489 * qJD(3);
t476 = t488 * qJD(3);
t470 = t824 - t956;
t469 = -t826 + t956;
t459 = 0.2e1 * t714 * t790;
t438 = 0.2e1 * t842 + t752;
t413 = t710 * t902 - t934;
t412 = t708 * t902 + t934;
t383 = t581 * t719;
t362 = t1018 / 0.2e1;
t361 = -t1019 / 0.2e1;
t355 = qJD(6) * t572 - t959;
t354 = -qJD(6) * t569 + t960;
t353 = -t1048 * t861 + t513;
t352 = 0.2e1 * t513;
t351 = -t1051 * t861 + t512;
t350 = 0.2e1 * t512;
t349 = 0.2e1 * t1067;
t344 = t1050 * t1085 + t886;
t336 = -t934 + (-t710 * t953 - t826) * t676;
t335 = t934 + (-t708 * t953 + t824) * t676;
t318 = t678 * t856 + t959;
t317 = t1090 * t718 - t960;
t293 = (t1057 + t666) * t581;
t289 = (-qJD(6) + t954) * t980 * t1092 - t694 * t596;
t273 = t293 * qJD(3);
t268 = t291 * qJD(3);
t266 = t290 * qJD(3);
t255 = t1106 + t975 / 0.2e1 - t804 * t718;
t253 = t1105 - t984 / 0.2e1 + t804 * t714;
t227 = t1106 + t976 / 0.2e1 - t793 * t718;
t225 = t1105 - t985 / 0.2e1 + t793 * t714;
t197 = -t786 / 0.2e1 + t861 * t894 + t885 + t718 * t889;
t196 = t1011 / 0.2e1 + t861 * t895 + t881 - t913 / 0.2e1;
t165 = t1104 + t1083;
t156 = 0.2e1 * t843 + t729;
t153 = t1100 - pkin(11) * t992 / 0.2e1 + t1081 * t718;
t152 = pkin(11) * t895 + t1081 * t714 - t1099;
t110 = -t1048 * t862 + t288;
t109 = 0.2e1 * t288;
t108 = -t1051 * t862 + t287;
t107 = 0.2e1 * t287;
t106 = 0.2e1 * t1074;
t101 = t1050 * t1086 + t888;
t95 = t603 + t727;
t93 = t603 + t740;
t92 = t604 + t734;
t90 = t604 + t802;
t82 = t796 - t797;
t71 = pkin(4) * t780 + t1044 * t847 - t738;
t65 = -t866 / 0.2e1 + t862 * t894 + t887 + t1004 / 0.2e1;
t64 = t1026 / 0.2e1 + t862 * t895 + t883 - t1005 / 0.2e1;
t60 = -t1041 / 0.2e1 + t720 + t1083 * pkin(11);
t53 = t1103 + t1082;
t51 = t417 * t916 + t416 * t917 - t1040 / 0.2e1 + t721;
t35 = qJD(3) * t80;
t30 = t361 + t733;
t28 = t361 + t755;
t27 = t362 + t735;
t25 = t362 + t803;
t22 = t551 * t916 + t550 * t917 - t1039 / 0.2e1 + t724;
t20 = t551 * t880 + t550 * t884 + t997 / 0.2e1 + t723;
t18 = -t739 + (t381 / 0.2e1 + t781) * pkin(4);
t7 = -t1043 / 0.2e1 + t722 + t1082 * pkin(11);
t5 = t328 * t916 + t327 * t917 + t1042 / 0.2e1 + t726;
t3 = t328 * t880 + t327 * t884 - t1020 / 0.2e1 + t725;
t2 = qJD(3) * t36 + qJD(4) * t33;
t1 = qJD(3) * t11 + qJD(4) * t10 + qJD(5) * t9;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t111, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t19 + t1117 * t927; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t932, t581 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, t381 * qJD(4) - t719 * t932, qJD(4) * t383 + t716 * t932, -qJD(3) * t863 (-pkin(3) * t582 - pkin(9) * t863) * qJD(3) + t812, 0, 0, 0, 0, 0, 0, t293 * t927 + t676 * t932, t292 * t927 + t678 * t932 (-t368 * t678 - t369 * t676) * qJD(3) - t1033 (t1089 * t369 + t582 * t703 - t1021) * qJD(3) + t18 * qJD(4) + t817, 0, 0, 0, 0, 0, 0 (t327 * t676 - t368 * t994) * qJD(3) + t27 * qJD(4) + t25 * qJD(5) + t65 * qJD(6) (-t328 * t676 - t368 * t992) * qJD(3) + t30 * qJD(4) + t28 * qJD(5) + t64 * qJD(6) (-t327 * t718 - t328 * t714) * t953 + t927 * t43 (-t327 * t394 + t328 * t395 - t1021) * qJD(3) + t3 * qJD(4) + t5 * qJD(5) + t819; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t381 - qJD(4) * t458, t383 * qJD(3) + qJD(4) * t731, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t1086 + qJD(5) * t156 + t273, t1093, -t1009, t18 * qJD(3) + (-t1045 * t1086 - t715 * t862) * t1030 + t818, 0, 0, 0, 0, 0, 0, qJD(3) * t27 + qJD(5) * t106 + qJD(6) * t108 - t1086 * t947, qJD(3) * t30 + qJD(5) * t101 + qJD(6) * t110 + t1086 * t948, t53 * qJD(5) - t859 * t862 + t1034, t3 * qJD(3) + (t1086 * t702 - t862 * t864) * qJD(4) + t7 * qJD(5) + t821; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t156 - qJD(5) * t1086 + t273, t1093, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t25 + qJD(4) * t106 + qJD(6) * t107 - t1086 * t944, qJD(3) * t28 + qJD(4) * t101 + qJD(6) * t109 + t1086 * t945, t53 * qJD(4) - t858 * t862 + t1034, t5 * qJD(3) + t7 * qJD(4) + (-t862 * t872 - t1043) * qJD(5) + t823; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * qJD(3) + t108 * qJD(4) + t107 * qJD(5) - qJD(6) * t963, qJD(3) * t64 + qJD(4) * t110 + qJD(5) * t109 + qJD(6) * t200, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t439, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t230, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t112 + t1118 * t927; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t901, -t848, 0, 0, 0, 0, 0, 0, 0, 0 (-t716 * t871 - t717 * t949) * t713 (t716 * t950 - t719 * t871) * t713, t779 * t951 (-pkin(3) * t717 + pkin(9) * t779) * t951 + t809, 0, 0, 0, 0, 0, 0, t491 * t927 + t676 * t901, t490 * t927 + t678 * t901 (t615 * t678 - t616 * t676) * qJD(3) - t968 (t1089 * t616 + t703 * t989 + t1000) * qJD(3) + t71 * qJD(4) + t810, 0, 0, 0, 0, 0, 0 (t550 * t676 + t615 * t994) * qJD(3) + t92 * qJD(4) + t90 * qJD(5) + t197 * qJD(6) (-t551 * t676 + t615 * t992) * qJD(3) + t95 * qJD(4) + t93 * qJD(5) + t196 * qJD(6) (-t550 * t718 - t551 * t714) * t953 + t927 * t130 (-t394 * t550 + t395 * t551 + t1000) * qJD(3) + t20 * qJD(4) + t22 * qJD(5) + t813; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t664 * qJD(4) - t716 * t848, qJD(4) * t776 - t719 * t848, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t1085 + qJD(5) * t438 + t483, t1094, -t940, t71 * qJD(3) + (-t1045 * t1085 - t715 * t861) * t1030 + t811, 0, 0, 0, 0, 0, 0, qJD(3) * t92 + qJD(5) * t349 + qJD(6) * t351 - t1085 * t947, qJD(3) * t95 + qJD(5) * t344 + qJD(6) * t353 + t1085 * t948, t165 * qJD(5) - t859 * t861 + t969, t20 * qJD(3) + (t1085 * t702 - t861 * t864) * qJD(4) + t60 * qJD(5) + t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t438 - qJD(5) * t1085 + t483, t1094, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t90 + qJD(4) * t349 + qJD(6) * t350 - t1085 * t944, qJD(3) * t93 + qJD(4) * t344 + qJD(6) * t352 + t1085 * t945, t165 * qJD(4) - t858 * t861 + t969, t22 * qJD(3) + t60 * qJD(4) + (-t861 * t872 - t1041) * qJD(5) + t822; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197 * qJD(3) + t351 * qJD(4) + t350 * qJD(5) - qJD(6) * t805, qJD(3) * t196 + qJD(4) * t353 + qJD(5) * t352 + qJD(6) * t466, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t812, 0, 0, 0, 0, 0, 0, -t927 * t290, -t927 * t291, -t1033, -qJD(4) * t17 - t817, 0, 0, 0, 0, 0, 0, qJD(4) * t26 + qJD(5) * t24 - qJD(6) * t62, qJD(4) * t31 + qJD(5) * t29 - qJD(6) * t63, t44 * t927, qJD(4) * t4 + qJD(5) * t6 - t819; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t809, 0, 0, 0, 0, 0, 0, -t927 * t488, -t927 * t489, -t968, -qJD(4) * t70 - t810, 0, 0, 0, 0, 0, 0, qJD(4) * t91 + qJD(5) * t89 - qJD(6) * t194, qJD(4) * t96 + qJD(5) * t94 - qJD(6) * t195, t131 * t927, qJD(4) * t21 + qJD(5) * t23 - t813; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t900, t695 * qJD(4), 0, -t900, 0, 0, -pkin(3) * t929, -pkin(3) * t928, 0, 0, -t1098, -t927 * t926, 0, t1098, 0, 0, qJD(4) * t552 + t678 * t946, qJD(4) * t553 - t676 * t946, 0, qJD(4) * t329, -t1098 * t710 - t850, -qJD(6) * t579 + t676 * t782, t487 * t927 - t714 * t903, -t1098 * t708 + t850, -t486 * t927 - t718 * t903, t1098, qJD(4) * t113 + qJD(5) * t115 + qJD(6) * t295, qJD(4) * t114 + qJD(5) * t116 + qJD(6) * t294, -qJD(4) * t72 - qJD(5) * t77, qJD(4) * t78 + qJD(5) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t699, t930, t928, -t699, -t929, 0, -pkin(9) * t928 - t921, pkin(9) * t929 - t922, 0, 0, -t902, -t935, -t596, t902, -t1090, 0, -qJD(4) * t1089 + qJD(5) * t533 - t1115, t1095 + t957 (t909 - t993) * t1030 + t764 (-t1045 * t1089 - t715 * t860) * t1030 + t765, t336, t289, t318, t335, t317, t498 (t714 * t814 - t1099) * qJD(4) + t152 * qJD(5) + t227 * qJD(6) + t769 (t718 * t814 + t1100) * qJD(4) + t153 * qJD(5) + t225 * qJD(6) + t768, qJD(4) * t816 + t82 * qJD(5) - t1113 (t1089 * t702 + t701 * t816) * qJD(4) + t51 * qJD(5) + t775; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t902, -t935, -t596, t902, -t1090, 0, qJD(4) * t533 - qJD(5) * t1089 - t1116, t1095 - t906, 0, 0, t336, t289, t318, t335, t317, t498, t152 * qJD(4) + (t714 * t839 - t1099) * qJD(5) + t255 * qJD(6) + t767, t153 * qJD(4) + (t718 * t839 + t1100) * qJD(5) + t253 * qJD(6) + t766, t82 * qJD(4) + qJD(5) * t815 - t1114, t51 * qJD(4) + (pkin(11) * t815 - t1040) * qJD(5) + t774; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t338, t1109 * t826 - t955, t714 * t808, t338, t790, t577, qJD(4) * t227 + qJD(5) * t255 - qJD(6) * t395 - t763, qJD(4) * t225 + qJD(5) * t253 + qJD(6) * t394 - t762, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t305 + t266, t268, t1009, qJD(3) * t17 - t818, 0, 0, 0, 0, 0, 0, -qJD(3) * t26 + qJD(5) * t105, -qJD(3) * t31 + qJD(5) * t102, -t1032, -qJD(3) * t4 + qJD(5) * t8 - t821; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5) * t527 + t476, t478, t940, qJD(3) * t70 - t811, 0, 0, 0, 0, 0, 0, -qJD(3) * t91 + qJD(5) * t348, -qJD(3) * t96 + qJD(5) * t345, -t962, -qJD(3) * t21 + qJD(5) * t61 - t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t699, -t930, 0, t699, 0, 0, t921, t922, 0, 0, t902, t935, 0, -t902, 0, 0, qJD(5) * t609 + t1115, -t957 + t966, -t764, -t765, t413, t459, t355, t412, t354, -t498, qJD(5) * t151 + qJD(6) * t226 - t769, qJD(5) * t154 + qJD(6) * t224 - t768, qJD(5) * t83 + t1113, qJD(5) * t52 - t775; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t919, -pkin(4) * t869, 0, 0, t698, t687, 0, -t698, 0, 0, t702 * t943 - t718 * t919, t702 * t942 + t693, t661, t564 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t840 + t760 (-t870 - t869) * pkin(4), 0, 0, t698, t687, 0, -t698, 0, 0, qJD(6) * t627 - t718 * t840 + t761, qJD(6) * t628 + t693 - t737, t661 - t747 (-pkin(5) * t715 + pkin(11) * t778) * t1031 + t749; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t470, t531, t523, t469, t522, -t931, qJD(5) * t627 - t701 * t942 - t742, qJD(5) * t628 + t701 * t943 - t741, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t305 + t266, t268, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t24 - qJD(4) * t105, -qJD(3) * t29 - qJD(4) * t102, -t1032, -qJD(3) * t6 - qJD(4) * t8 - t823; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t527 + t476, t478, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t89 - qJD(4) * t348, -qJD(3) * t94 - qJD(4) * t345, -t962, -qJD(3) * t23 - qJD(4) * t61 - t822; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t902, t935, 0, -t902, 0, 0, -qJD(4) * t609 + t1116, t906 + t966, 0, 0, t413, t459, t355, t412, t354, -t498, -qJD(4) * t151 - qJD(6) * t254 - t767, -qJD(4) * t154 + qJD(6) * t252 - t766, -qJD(4) * t83 + t1114, -qJD(4) * t52 - t774; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t760 + t920, pkin(4) * t870, 0, 0, t698, t687, 0, -t698, 0, 0, -qJD(6) * t625 + t718 * t920 - t761, -qJD(6) * t626 + t737, t747, -t749; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t698, t687, 0, -t698, 0, 0, -pkin(5) * t943, -pkin(5) * t942, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t470, t531, t523, t469, t522, -t931, -pkin(11) * t942 - t753, pkin(11) * t943 - t754, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t62, qJD(3) * t63, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t194, qJD(3) * t195, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t338, t782 + t955, -t572 * t927 + t714 * t902, -t338, t569 * t927 + t718 * t902, t577, -qJD(4) * t226 + qJD(5) * t254 + t763, -qJD(4) * t224 - qJD(5) * t252 + t762, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t469, -t531, -t563, t470, t933, t931, qJD(5) * t625 + t742, qJD(5) * t626 + t741, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t469, -t531, -t563, t470, t933, t931, t753, t754, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t12;
