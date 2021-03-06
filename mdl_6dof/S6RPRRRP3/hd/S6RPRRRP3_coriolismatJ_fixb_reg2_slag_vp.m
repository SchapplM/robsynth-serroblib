% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:21
% EndTime: 2019-03-09 06:05:58
% DurationCPUTime: 25.97s
% Computational Cost: add. (19589->880), mult. (39603->1070), div. (0->0), fcn. (39698->8), ass. (0->657)
t695 = sin(qJ(4));
t1060 = cos(qJ(5));
t1092 = -pkin(9) - pkin(8);
t814 = t1092 * t1060;
t628 = t695 * t814;
t694 = sin(qJ(5));
t697 = cos(qJ(4));
t975 = t694 * t697;
t822 = t1092 * t975;
t522 = -t628 - t822;
t877 = t694 * t1092;
t820 = t697 * t877;
t523 = t628 + t820;
t1148 = t522 + t523;
t624 = t1060 * t695 + t975;
t696 = sin(qJ(3));
t589 = t624 * t696;
t1079 = -t589 / 0.2e1;
t866 = t1060 * t697;
t974 = t695 * t696;
t591 = -t694 * t974 + t696 * t866;
t1077 = t591 / 0.2e1;
t698 = cos(qJ(3));
t590 = t624 * t698;
t973 = t695 * t698;
t650 = t694 * t973;
t593 = t698 * t866 - t650;
t279 = t593 * t589 + t591 * t590;
t1136 = t279 * qJD(1);
t620 = t694 * t695 - t866;
t211 = t589 * t620 - t591 * t624;
t887 = qJD(4) + qJD(5);
t1145 = t887 * t211;
t1147 = t1136 + t1145;
t481 = t593 * t620;
t482 = t624 * t590;
t1146 = t1145 - qJD(3) * (t482 + t481) - t1136;
t1075 = -t620 / 0.2e1;
t1073 = -t624 / 0.2e1;
t939 = t591 * t887;
t1072 = t624 / 0.2e1;
t1074 = t620 / 0.2e1;
t748 = t1072 * t589 + t1074 * t591;
t1144 = t887 * t748;
t890 = t698 * qJD(1);
t862 = t591 * t890;
t1143 = t939 - t862;
t1104 = t696 * t589 - t590 * t698;
t1131 = qJD(1) * t1104;
t892 = t696 * qJD(3);
t600 = t620 * t892;
t1142 = t1131 + t600;
t1051 = t696 * pkin(5);
t1052 = t696 * pkin(4);
t1049 = t698 * pkin(8);
t1053 = t696 * pkin(3);
t648 = -t1049 + t1053;
t630 = t697 * t648;
t669 = sin(pkin(10)) * pkin(1) + pkin(7);
t646 = t696 * t669;
t632 = t695 * t646;
t525 = t632 + t630;
t969 = t697 * t698;
t433 = -pkin(9) * t969 + t1052 + t525;
t869 = t1060 * t433;
t627 = t695 * t648;
t971 = t696 * t697;
t875 = t669 * t971;
t526 = t627 - t875;
t475 = -pkin(9) * t973 + t526;
t976 = t694 * t475;
t948 = t869 / 0.2e1 - t976 / 0.2e1;
t802 = -t1051 / 0.2e1 - t948;
t684 = -pkin(4) * t697 - pkin(3);
t1140 = -t684 / 0.2e1;
t647 = t698 * t669;
t633 = t697 * t647;
t671 = -cos(pkin(10)) * pkin(1) - pkin(2);
t1050 = t698 * pkin(3);
t806 = -t696 * pkin(8) - t1050;
t717 = t806 + t671;
t484 = t695 * t717 + t633;
t773 = pkin(9) * t974 - t484;
t1106 = t1060 * t773;
t1139 = t1106 / 0.2e1;
t1123 = t589 * t887;
t619 = t624 ^ 2;
t388 = t620 ^ 2 - t619;
t1138 = t887 * t388;
t1121 = t887 * t624;
t1135 = t620 * t1121;
t1120 = t887 * t698;
t1094 = t591 ^ 2;
t693 = t698 ^ 2;
t543 = t693 + t1094;
t922 = qJD(3) * t624;
t1134 = -t543 * qJD(1) - t591 * t922 + t1120;
t925 = qJD(1) * t589;
t1133 = (qJD(3) * t620 + t925) * t590 + t1144;
t1132 = -t590 * t925 + t1144;
t1122 = t694 * t773;
t725 = t1092 * t696 + t671;
t709 = (t725 - t1050) * t697;
t989 = t669 * t695;
t834 = -pkin(4) - t989;
t703 = t698 * t834 + t709;
t362 = t1060 * t703;
t241 = -t1122 - t362;
t702 = t694 * t703;
t242 = -t1106 + t702;
t1129 = t242 * t1073 + t1075 * t241;
t860 = t620 * t922;
t1128 = qJD(1) * t748 + t860;
t924 = qJD(1) * t591;
t864 = t589 * t924;
t1127 = qJD(3) * t748 + t864;
t315 = t589 ^ 2 - t1094;
t1126 = t279 * qJD(3) - t315 * t887;
t112 = qJD(1) * t315 + qJD(3) * t211;
t118 = qJD(1) * t211 + qJD(3) * t388;
t1125 = qJD(3) * t1104 - t1120 * t591;
t1124 = 0.2e1 * t695;
t629 = t697 * t814;
t799 = t1092 * t866;
t756 = -t799 / 0.2e1;
t821 = t695 * t877;
t391 = -t629 / 0.2e1 + t756 + t821;
t524 = -t799 + t821;
t1119 = -t391 * qJD(4) - t524 * qJD(5);
t1005 = t524 * t698;
t491 = t1005 / 0.2e1;
t601 = pkin(4) * t974 + t646;
t1057 = t589 * pkin(5);
t999 = t591 * qJ(6);
t801 = -t999 + t1057;
t312 = t801 + t601;
t1055 = t620 * pkin(5);
t991 = t624 * qJ(6);
t800 = -t991 + t1055;
t451 = t684 + t800;
t753 = t312 * t1072 + t451 * t1077;
t1117 = t491 + t753;
t1078 = t589 / 0.2e1;
t521 = -t629 + t821;
t1116 = t1077 * t1148 + t1078 * t521 + t524 * t1079;
t508 = t629 / 0.2e1 + t756;
t900 = t508 * qJD(3);
t701 = -t702 / 0.2e1 + t1139;
t1066 = t694 / 0.2e1;
t844 = t698 * t1066;
t700 = pkin(4) * t844 + t701;
t876 = t695 * t647;
t426 = t709 - t876;
t977 = t694 * t426;
t804 = t1139 - t977 / 0.2e1;
t93 = t700 - t804;
t1115 = -t93 * qJD(1) + t900;
t1061 = -t698 / 0.2e1;
t1054 = t694 * pkin(4);
t670 = qJ(6) + t1054;
t699 = -t701 + (qJ(6) + t670) * t1061;
t77 = -t699 - t804;
t1114 = -t77 * qJD(1) + t900;
t1008 = t522 * t698;
t488 = t1008 / 0.2e1;
t1100 = t1075 * t601 + t1140 * t589;
t389 = t694 * t433;
t450 = t1060 * t475;
t950 = t389 / 0.2e1 + t450 / 0.2e1;
t722 = -t950 - t1100;
t102 = t488 + t722;
t871 = t698 * t1060;
t811 = -t871 / 0.2e1;
t738 = t650 / 0.2e1 + t697 * t811;
t842 = t620 * t1061;
t442 = t842 + t738;
t901 = t442 * qJD(2);
t921 = qJD(3) * t684;
t1113 = qJD(1) * t102 + t620 * t921 + t901;
t1059 = pkin(4) * t695;
t432 = t1059 * t624 - t620 * t684;
t1062 = -t697 / 0.2e1;
t1065 = -t695 / 0.2e1;
t1067 = -t694 / 0.2e1;
t1007 = t523 * t698;
t490 = -t1007 / 0.2e1;
t72 = t490 + (t591 * t1065 + (t1062 * t624 + t1067) * t696) * pkin(4) + t722;
t1112 = qJD(1) * t72 - qJD(3) * t432 + t901;
t1018 = t451 * t620;
t494 = pkin(5) * t624 + qJ(6) * t620;
t230 = -t494 * t624 + t1018;
t487 = -t1008 / 0.2e1;
t1099 = t312 * t1075 + t1079 * t451;
t723 = t950 + t1099;
t1086 = -t494 / 0.2e1;
t1002 = t589 * qJ(6);
t1058 = pkin(5) * t591;
t386 = t1002 + t1058;
t750 = t1073 * t386 + t1086 * t591;
t972 = t696 * qJ(6);
t45 = t487 + t723 - t750 + t972;
t1071 = -t650 / 0.2e1;
t809 = t866 / 0.2e1;
t440 = t1071 + (t809 + t1074) * t698;
t903 = t440 * qJD(2);
t1111 = qJD(1) * t45 - qJD(3) * t230 + t903;
t455 = t494 + t1059;
t216 = -t455 * t624 + t1018;
t489 = t1007 / 0.2e1;
t886 = pkin(4) * t971;
t330 = t386 + t886;
t752 = t1072 * t330 + t1077 * t455;
t1069 = t670 / 0.2e1;
t1091 = qJ(6) / 0.2e1;
t805 = (t1069 + t1091) * t696;
t40 = t489 + t805 + t723 + t752;
t1110 = qJD(1) * t40 - qJD(3) * t216 + t903;
t883 = t1060 * pkin(4);
t810 = t871 / 0.2e1;
t741 = -t362 / 0.2e1 + pkin(4) * t810;
t870 = t1060 * t426;
t95 = t870 / 0.2e1 + t741;
t879 = t95 * qJD(1) - qJD(4) * t883;
t682 = t698 * t696;
t757 = t589 * t590 + t591 * t593 - t682;
t1109 = qJD(2) * t757;
t724 = (t521 - t524) * t624 - t1148 * t620;
t1108 = qJD(3) * t724;
t1107 = qJD(4) * t724;
t1105 = t757 * qJD(3);
t231 = -t1106 + t725 * t975 + (t694 * (-pkin(3) * t697 + t834) - qJ(6)) * t698;
t233 = t698 * pkin(5) + t241;
t1011 = t522 * t589;
t848 = -t1011 / 0.2e1;
t1103 = t231 * t1073 + t233 * t1075 + t848;
t1102 = t848 + t1129;
t1101 = t1073 * t601 + t1140 * t591;
t1098 = (-t481 + t482) * qJD(3);
t1093 = -pkin(5) / 0.2e1;
t1090 = t231 / 0.2e1;
t1089 = t241 / 0.2e1;
t251 = -t1106 + t977;
t1088 = -t251 / 0.2e1;
t252 = t870 + t1122;
t1087 = t252 / 0.2e1;
t1085 = -t521 / 0.2e1;
t1084 = t521 / 0.2e1;
t1083 = t522 / 0.2e1;
t1082 = -t523 / 0.2e1;
t1081 = t523 / 0.2e1;
t1080 = t524 / 0.2e1;
t846 = -t590 / 0.2e1;
t845 = t590 / 0.2e1;
t1076 = t593 / 0.2e1;
t1070 = -t670 / 0.2e1;
t683 = -t883 - pkin(5);
t1068 = t683 / 0.2e1;
t1064 = -t696 / 0.2e1;
t1063 = t696 / 0.2e1;
t1056 = t590 * pkin(5);
t1048 = pkin(4) * qJD(4);
t1047 = pkin(4) * qJD(5);
t94 = t700 + t804;
t1046 = -t251 * qJD(4) + t94 * qJD(5);
t1045 = t94 * qJD(4) - t242 * qJD(5);
t998 = t591 * t524;
t359 = -t998 / 0.2e1;
t360 = t998 / 0.2e1;
t708 = t359 + t1011 / 0.2e1 + t360 - t1129;
t16 = t708 + t1102;
t1044 = qJD(1) * t16;
t1024 = t312 * t591;
t1030 = t251 * t698;
t66 = t330 * t589 + t1024 + t1030;
t1043 = qJD(1) * t66;
t1032 = t242 * t698;
t69 = t386 * t589 + t1024 + t1032;
t1042 = qJD(1) * t69;
t1041 = qJD(3) * t16;
t1036 = t241 * t524;
t156 = -t1036 / 0.2e1;
t24 = t156 + t1036 / 0.2e1;
t1040 = qJD(3) * t24;
t704 = t708 + t1103;
t765 = qJ(6) * t846 + t1093 * t593;
t12 = t704 - t765;
t1039 = t12 * qJD(1);
t839 = t1089 - t233 / 0.2e1;
t841 = t1090 - t242 / 0.2e1;
t20 = t1061 * t386 - t589 * t841 - t591 * t839;
t1038 = t20 * qJD(1);
t837 = t1087 + t233 / 0.2e1;
t840 = t1090 + t1088;
t21 = t1061 * t330 - t589 * t840 + t591 * t837;
t1037 = t21 * qJD(1);
t1034 = t241 * t698;
t1031 = t251 * t522;
t1029 = t252 * t524;
t1028 = t252 * t698;
t28 = (-t231 + t242) * t591 + (-t233 + t241) * t589;
t1027 = t28 * qJD(1);
t836 = t1087 + t1089;
t838 = t242 / 0.2e1 + t1088;
t843 = -t969 / 0.2e1;
t29 = t1052 * t843 - t589 * t838 + t591 * t836;
t1026 = t29 * qJD(1);
t31 = (-t231 + t251) * t591 + (-t233 - t252) * t589;
t1025 = t31 * qJD(1);
t265 = t450 + t389;
t255 = t265 + t972;
t264 = t869 - t976;
t256 = -t264 - t1051;
t33 = -t231 * t590 + t233 * t593 - t255 * t589 + t256 * t591;
t1022 = t33 * qJD(1);
t34 = (-t242 + t251) * t591 + (-t241 - t252) * t589;
t1021 = t34 * qJD(1);
t35 = t241 * t593 - t242 * t590 - t264 * t591 - t265 * t589;
t1020 = t35 * qJD(1);
t1017 = t451 * t624;
t483 = -t697 * t717 + t876;
t1016 = t483 * t698;
t1015 = t484 * t698;
t668 = pkin(4) * t973;
t602 = t647 + t668;
t996 = t593 * qJ(6);
t313 = t602 - t996 + t1056;
t50 = -t233 * t696 + t256 * t698 + t312 * t590 + t313 * t589;
t1014 = t50 * qJD(1);
t1013 = t521 * t591;
t1012 = t521 * t698;
t1009 = t522 * t696;
t1006 = t524 * t696;
t1004 = t525 * t696;
t1003 = t526 * t696;
t997 = t591 * t694;
t60 = -t242 * t696 + t265 * t698 + t602 * t591 + t601 * t593;
t995 = t60 * qJD(1);
t994 = t601 * t589;
t990 = t624 * t694;
t271 = t312 * t589;
t67 = t330 * t591 + t1028 - t271;
t988 = t67 * qJD(1);
t987 = t670 * t589;
t986 = t670 * t591;
t985 = t670 * t624;
t68 = t386 * t591 - t1034 - t271;
t984 = t68 * qJD(1);
t983 = t683 * t589;
t982 = t683 * t591;
t981 = t683 * t620;
t690 = t695 ^ 2;
t978 = t690 * t698;
t691 = t696 ^ 2;
t970 = t697 * t691;
t968 = t698 * t455;
t967 = t698 * t494;
t963 = t95 * qJD(4);
t962 = t95 * qJD(5);
t99 = t231 * t698 + t1024;
t961 = t99 * qJD(1);
t497 = t508 * qJD(4);
t936 = t694 * t843 + t695 * t811;
t439 = t845 + t936;
t904 = t439 * qJD(2);
t946 = t904 + t497;
t945 = -t508 * qJD(5) + t904;
t944 = t887 * t442;
t441 = t1071 + (t809 + t1075) * t698;
t943 = t887 * t441;
t942 = t887 * t439;
t937 = t695 * t810 + t697 * t844;
t437 = t845 + t937;
t941 = t887 * t437;
t662 = t691 - t693;
t692 = t697 ^ 2;
t663 = t692 - t690;
t935 = qJ(6) * qJD(5);
t446 = t601 * t591;
t116 = -t589 * t886 - t1030 - t446;
t934 = qJD(1) * t116;
t117 = t591 * t886 + t1028 - t994;
t933 = qJD(1) * t117;
t138 = -t994 - t1034;
t932 = qJD(1) * t138;
t139 = -t446 - t1032;
t931 = qJD(1) * t139;
t331 = -t691 * t989 - t1016;
t929 = qJD(1) * t331;
t332 = -t669 * t970 - t1015;
t928 = qJD(1) * t332;
t370 = -t591 * t696 + t593 * t698;
t926 = qJD(1) * t370;
t818 = t883 / 0.2e1;
t541 = t591 * t818;
t819 = -t883 / 0.2e1;
t342 = t591 * t819 + t541;
t923 = qJD(2) * t342;
t920 = qJD(3) * t695;
t919 = qJD(3) * t697;
t918 = qJD(4) * t695;
t917 = qJD(4) * t697;
t916 = qJD(4) * t698;
t915 = qJD(5) * t684;
t914 = qJD(6) * t698;
t127 = (t1004 - t1016) * t697 + (t1003 + t1015) * t695;
t913 = t127 * qJD(1);
t261 = t483 * t696 + (t525 - 0.2e1 * t632) * t698;
t911 = t261 * qJD(1);
t262 = t526 * t698 + (-t484 + 0.2e1 * t633) * t696;
t910 = t262 * qJD(1);
t436 = t846 + t937;
t906 = t436 * qJD(3);
t397 = t437 * qJD(1);
t416 = t437 * qJD(3);
t438 = t846 + t936;
t398 = t438 * qJD(1);
t905 = t438 * qJD(3);
t399 = t439 * qJD(3);
t402 = t441 * qJD(1);
t902 = t441 * qJD(3);
t403 = t442 * qJD(3);
t577 = t589 * qJD(6);
t580 = t591 * qJD(6);
t609 = (t690 / 0.2e1 - t692 / 0.2e1) * t696;
t897 = t609 * qJD(4);
t608 = t620 * qJD(6);
t623 = t662 * t695;
t896 = t623 * qJD(1);
t626 = t693 * t697 - t970;
t895 = t626 * qJD(1);
t894 = t662 * qJD(1);
t893 = t696 * qJD(1);
t891 = t696 * qJD(4);
t889 = t698 * qJD(3);
t687 = qJD(5) * t883;
t888 = t687 + qJD(6);
t885 = -t1060 / 0.2e1;
t884 = t1060 / 0.2e1;
t882 = t694 * t1047;
t881 = t1059 / 0.2e1;
t880 = t1054 / 0.2e1;
t874 = -t521 * qJD(4) - t391 * qJD(5) - t904;
t873 = -t904 + t1119;
t260 = -t416 - t939;
t872 = t600 - t941;
t868 = t1060 * t589;
t867 = t1060 * t620;
t863 = t589 * t890;
t861 = t697 * t893;
t859 = t695 * t919;
t858 = t695 * t916;
t857 = t697 * t916;
t599 = t624 * t892;
t856 = t671 * t893;
t855 = t671 * t890;
t854 = t695 * t917;
t853 = t696 * t889;
t852 = t697 * t892;
t851 = t696 * t890;
t850 = -t1029 / 0.2e1;
t849 = t1012 / 0.2e1;
t492 = -t1005 / 0.2e1;
t835 = t1081 + t1083;
t833 = (t690 + t692) * t698;
t832 = t522 * t590 + t593 * t524;
t831 = t522 * t593 - t524 * t590;
t824 = t886 / 0.2e1;
t823 = -qJD(4) + t890;
t817 = t950 - t1099;
t816 = -t950 + t1100;
t815 = t948 - t1101;
t813 = t695 * t852;
t812 = t691 * t854;
t808 = t948 + t1051;
t807 = t890 - qJD(4) / 0.2e1;
t651 = -qJD(5) + t823;
t19 = t1061 * t313 + t1063 * t312 + t1076 * t231 + t1077 * t255 + t1078 * t256 + t233 * t845;
t22 = t231 * t255 + t233 * t256 + t312 * t313;
t796 = t22 * qJD(1) + t19 * qJD(2);
t25 = -t231 * t241 + t233 * t242 + t312 * t386;
t795 = t25 * qJD(1) + t20 * qJD(2);
t26 = t231 * t252 + t233 * t251 + t312 * t330;
t794 = t26 * qJD(1) + t21 * qJD(2);
t793 = qJD(1) * t24;
t27 = t1061 * t602 + t1063 * t601 + t1076 * t242 + t1077 * t265 + t1079 * t264 + t241 * t845;
t39 = -t241 * t264 + t242 * t265 + t601 * t602;
t792 = t39 * qJD(1) + t27 * qJD(2);
t38 = t241 * t251 + t242 * t252 + t601 * t886;
t791 = t38 * qJD(1) + t29 * qJD(2);
t714 = t360 + t835 * t589 - t1013 / 0.2e1;
t747 = t1068 * t593 + t1070 * t590;
t9 = t620 * t837 + t624 * t840 + t714 + t747;
t790 = -t9 * qJD(1) + t1108;
t789 = t521 * t522 + t523 * t524;
t788 = -t525 * t695 + t526 * t697;
t787 = t19 * qJD(1) + t1109;
t59 = t241 * t696 + t264 * t698 - t602 * t589 - t601 * t590;
t786 = t59 * qJD(1);
t49 = -t231 * t696 + t255 * t698 + t312 * t593 + t313 * t591;
t785 = t49 * qJD(1);
t784 = t27 * qJD(1) + t1109;
t718 = (t1067 * t590 + t593 * t885) * pkin(4);
t13 = t620 * t836 + t624 * t838 + t714 + t718;
t783 = -t13 * qJD(1) + t1108;
t215 = t455 * t620 + t1017;
t710 = t1074 * t330 + t1078 * t455 + t753 + t849;
t42 = (t1068 + t1093) * t696 + t710 - t948;
t782 = -qJD(1) * t42 - qJD(3) * t215;
t229 = t494 * t620 + t1017;
t733 = t492 - t753;
t751 = t1074 * t386 + t1078 * t494;
t47 = t733 - t751 + t808;
t780 = qJD(1) * t47 - qJD(3) * t229;
t431 = t1059 * t620 + t624 * t684;
t721 = t948 + t1101;
t73 = -t1012 / 0.2e1 + (t589 * t1065 + (t1062 * t620 + t884) * t696) * pkin(4) + t721;
t778 = qJD(1) * t73 - qJD(3) * t431;
t776 = qJD(4) * t342;
t96 = -t1122 - t870 / 0.2e1 + t741;
t775 = t96 * qJD(4) + t241 * qJD(5);
t774 = -t252 * qJD(4) + t96 * qJD(5);
t772 = t823 * t696;
t105 = (t1003 / 0.2e1 + t1015 / 0.2e1) * t697 + (-t1004 / 0.2e1 + t1016 / 0.2e1) * t695 + (t691 / 0.2e1 - t693 / 0.2e1) * t669;
t140 = t669 ^ 2 * t682 - t483 * t525 + t484 * t526;
t771 = t140 * qJD(1) + t105 * qJD(2);
t530 = t696 * t833 - t682;
t770 = -t105 * qJD(1) - t530 * qJD(2);
t737 = t1093 + t819 - t683 / 0.2e1;
t760 = t1070 + t880 + t1091;
t124 = t589 * t737 + t591 * t760;
t179 = t620 * t737 + t624 * t760;
t769 = qJD(1) * t124 + qJD(3) * t179;
t390 = -t628 - t820 / 0.2e1 - t822 / 0.2e1;
t768 = -qJD(4) * t390 - qJD(5) * t522;
t767 = qJD(4) * t523 - qJD(5) * t390;
t766 = t1049 / 0.2e1 - t1053 / 0.2e1;
t764 = t251 * pkin(5) / 0.2e1 - t252 * qJ(6) / 0.2e1;
t763 = t1091 * t255 + t1093 * t256;
t762 = pkin(5) * t1084 + qJ(6) * t1082;
t761 = -t1056 / 0.2e1 + t996 / 0.2e1;
t735 = t766 * t695;
t501 = t627 / 0.2e1 - t735;
t759 = pkin(3) * t919 - qJD(1) * t501;
t734 = t766 * t697;
t502 = -t630 / 0.2e1 + t734;
t758 = pkin(3) * t920 - qJD(1) * t502;
t755 = t1068 * t256 + t1069 * t255;
t754 = t312 * t1086 - t386 * t451 / 0.2e1;
t749 = t1068 * t590 + t1069 * t593;
t746 = t697 * t772;
t101 = t492 + t721;
t745 = qJD(1) * t101 - t624 * t921;
t290 = t1064 + t748;
t743 = qJD(1) * t290 + t860;
t529 = -qJD(1) * t609 + t859;
t732 = t1066 * t265 + t264 * t884;
t509 = qJD(1) * t695 * t970 + qJD(3) * t609;
t622 = t663 * t691;
t731 = qJD(1) * t622 + 0.2e1 * t813;
t730 = -qJD(3) * t663 + t1124 * t861;
t705 = t231 * t1081 + t233 * t1084 + t1031 / 0.2e1 + t312 * t455 / 0.2e1 + t330 * t451 / 0.2e1;
t3 = t850 - t705 + t755;
t712 = -t835 * t591 + (t1080 + t1085) * t589;
t51 = t968 / 0.2e1 + t712 + t749;
t97 = t451 * t455 + t789;
t729 = -t3 * qJD(1) - t51 * qJD(2) + t97 * qJD(3);
t100 = t451 * t494;
t5 = t522 * t841 + t524 * t839 + t754 + t763;
t53 = t967 / 0.2e1 + t761;
t728 = -t5 * qJD(1) - t53 * qJD(2) + t100 * qJD(3);
t153 = t1059 * t684 + t789;
t719 = (t1066 * t593 + t590 * t885) * pkin(4);
t64 = t668 / 0.2e1 + t719 + t712;
t713 = t241 * t1085 + t242 * t1082 - t1031 / 0.2e1;
t7 = t850 + (t601 * t1065 + t1140 * t971 + t732) * pkin(4) + t713;
t727 = -t7 * qJD(1) - t64 * qJD(2) + t153 * qJD(3);
t720 = (qJD(3) * t590 + t939) * t589;
t62 = t802 + t1117;
t716 = qJD(1) * t62 + qJD(2) * t436 + t451 * t922;
t715 = t252 * t1075 + t523 * t1079 + t359 + t251 * t1072 + t1013 / 0.2e1;
t125 = -t589 * t760 + t591 * t737;
t707 = (t1066 * t233 + t231 * t884) * pkin(4) + t241 * t1070 + t242 * t1068;
t18 = t707 + t764;
t565 = (t1060 * t670 + t683 * t694) * pkin(4);
t706 = (t1066 * t522 + t524 * t884) * pkin(4) + t522 * t1070 + t524 * t1068;
t83 = t706 + t762;
t711 = t18 * qJD(1) - t125 * qJD(2) + t83 * qJD(3) + t565 * qJD(4);
t689 = qJ(6) * qJD(6);
t680 = t692 * t698;
t676 = -t893 / 0.2e1;
t675 = t893 / 0.2e1;
t674 = t892 / 0.2e1;
t666 = t695 * t892;
t661 = t670 * qJD(6);
t618 = t651 * qJ(6);
t617 = t807 * t696;
t585 = (-qJD(5) / 0.2e1 + t807) * t696;
t474 = t624 * t580;
t443 = t842 - t738;
t380 = t632 + t630 / 0.2e1 + t734;
t379 = t875 - t627 / 0.2e1 - t735;
t329 = qJD(3) * t619 + t624 * t924;
t328 = t342 * qJD(5);
t308 = t862 + t416;
t307 = -t863 + t902;
t306 = -t862 + t905;
t291 = t1064 - t748;
t289 = t1121 - t397;
t288 = -t620 * t887 - t402;
t287 = -t1121 - t398;
t259 = -t902 + t1123;
t244 = qJD(3) * t440;
t232 = -t370 * qJD(3) + t1120 * t589;
t194 = -t399 + t1143;
t193 = t863 - t403 - t1123;
t192 = -t906 - t1143;
t178 = -t985 / 0.2e1 - t981 / 0.2e1 - t991 / 0.2e1 + t1055 / 0.2e1 + (-t867 / 0.2e1 + t990 / 0.2e1) * pkin(4);
t170 = t1029 / 0.2e1;
t154 = (qJD(3) * t593 - t1123) * t591;
t144 = t926 - t943;
t141 = qJD(3) * t443 - t1123;
t126 = -t987 / 0.2e1 + t541 + t982 / 0.2e1 + t589 * t880 - t1002 / 0.2e1 - t1058 / 0.2e1;
t123 = -t986 / 0.2e1 - t983 / 0.2e1 - t999 / 0.2e1 + t1057 / 0.2e1 + (-t868 / 0.2e1 + t997 / 0.2e1) * pkin(4);
t122 = t599 - t926 - t944;
t104 = t491 + t815;
t103 = t487 + t816;
t98 = -t593 * t924 - t1144;
t91 = t93 * qJD(4);
t86 = t93 * qJD(5);
t82 = t706 - t762;
t78 = t699 - t804;
t76 = (t922 + t924) * t593 - t1144;
t75 = t1052 * t1067 + t591 * t881 + t624 * t824 + t489 + t816;
t74 = t589 * t881 + t620 * t824 + t696 * t818 + t815 + t849;
t71 = t1048 * t694 + t1115;
t70 = -t1054 * t887 - t1115;
t65 = -t668 / 0.2e1 + t719 + t1116;
t63 = t733 + t802;
t61 = qJD(3) * t105;
t54 = -t967 / 0.2e1 + t761;
t52 = -t968 / 0.2e1 + t749 + t1116;
t48 = t751 + t808 + t1117;
t46 = t488 + t750 + t817 + t972;
t43 = t683 * t1064 + t710 - t802;
t41 = t490 + t805 - t752 + t817;
t23 = t24 * qJD(5);
t17 = t707 - t764;
t15 = t16 * qJD(5);
t14 = t718 + t715 + t1102;
t11 = t704 + t765;
t10 = t715 + t747 + t1103;
t8 = pkin(4) * t732 + t601 * t881 + t684 * t824 + t170 - t713;
t6 = t156 - t231 * t522 / 0.2e1 + t242 * t1083 + t233 * t1080 - t754 + t763;
t4 = t170 + t705 + t755;
t2 = qJD(3) * t27 + qJD(4) * t29;
t1 = qJD(3) * t19 + qJD(4) * t21 + qJD(5) * t20;
t30 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t853, -t662 * qJD(3), 0, -t853, 0, 0, t671 * t892, t671 * t889, 0, 0, t692 * t853 - t812, -t622 * qJD(4) - 0.2e1 * t698 * t813, -t626 * qJD(3) + t696 * t858, t690 * t853 + t812, -t623 * qJD(3) + t696 * t857, -t853, -qJD(3) * t261 - qJD(4) * t332, qJD(3) * t262 + qJD(4) * t331, -qJD(3) * t127, qJD(3) * t140, t154, -t1126, t232, t720, -t1125, -t853, -qJD(3) * t59 - qJD(4) * t116 - qJD(5) * t139, qJD(3) * t60 + qJD(4) * t117 + qJD(5) * t138, qJD(3) * t35 + qJD(4) * t34, qJD(3) * t39 + qJD(4) * t38, t154, t232, t1126, -t853, t1125, t720, qJD(3) * t50 + qJD(4) * t66 + qJD(5) * t69 - t580 * t589, t33 * qJD(3) + t31 * qJD(4) + t28 * qJD(5) + t577 * t698, -qJD(3) * t49 - qJD(4) * t67 - qJD(5) * t68 + qJD(6) * t543, qJD(3) * t22 + qJD(4) * t26 + qJD(5) * t25 - qJD(6) * t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t851, -t894, t889, -t851, -t892, 0, -t669 * t889 + t856, t669 * t892 + t855, 0, 0, -t897 + (t692 * t893 + t859) * t698 (t680 - t978) * qJD(3) + (-qJD(4) - t890) * t971 * t1124, t666 - t895, t897 + (t690 * t893 - t859) * t698, t852 - t896, -t617, -t911 + (t695 * t806 - t633) * qJD(3) + t380 * qJD(4), t910 + (t697 * t806 + t876) * qJD(3) + t379 * qJD(4), qJD(3) * t788 - t913 (-pkin(3) * t647 + pkin(8) * t788) * qJD(3) + t771, t76, t1146, t122, t1133, -t436 * t887 - t1142, -t585 (t590 * t684 + t602 * t620 - t1009) * qJD(3) + t74 * qJD(4) + t104 * qJD(5) - t786, t995 + (t593 * t684 + t602 * t624 - t1006) * qJD(3) + t75 * qJD(4) + t103 * qJD(5), t1020 + (-t264 * t624 - t265 * t620 + t831) * qJD(3) + t14 * qJD(4) + t15 (-t264 * t522 + t265 * t524 + t602 * t684) * qJD(3) + t8 * qJD(4) + t23 + t792, t76, t122, -t1146, -t585, -t942 + t1142, t1133, t1014 + (t313 * t620 + t451 * t590 - t1009) * qJD(3) + t43 * qJD(4) + t48 * qJD(5) + t291 * qJD(6), t1022 + (-t255 * t620 + t256 * t624 + t831) * qJD(3) + t10 * qJD(4) + t11 * qJD(5) - t442 * qJD(6) (-t313 * t624 - t451 * t593 + t1006) * qJD(3) + t41 * qJD(4) + t46 * qJD(5) + t474 - t785 (t255 * t524 + t256 * t522 + t313 * t451) * qJD(3) + t4 * qJD(4) + t6 * qJD(5) + t63 * qJD(6) + t796; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t509, -t731, t695 * t772, t509, t746, t674, qJD(3) * t380 - qJD(4) * t484 - t928, qJD(3) * t379 + qJD(4) * t483 + t929, 0, 0, -t1127, t112, t193, t1127, t192, t674, qJD(3) * t74 + t1046 - t934, qJD(3) * t75 + t774 + t933, t1021 + t14 * qJD(3) + (t868 - t997) * t1048, t8 * qJD(3) + (-t1060 * t251 + t252 * t694) * t1048 + t791, -t1127, t193, -t112, t674, t194, t1127, qJD(3) * t43 + t1043 + t1046, t1025 + t10 * qJD(3) + (-t983 - t986) * qJD(4) + t123 * qJD(5) - t577, t41 * qJD(3) - t774 - t914 - t988, t4 * qJD(3) + (t251 * t683 + t252 * t670) * qJD(4) + t17 * qJD(5) + t78 * qJD(6) + t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1127, t112, t193, t1127, t192, t674, qJD(3) * t104 + t1045 - t931, qJD(3) * t103 + t775 + t932, t1041, t1040, -t1127, t193, -t112, t674, t194, t1127, qJD(3) * t48 + t1042 + t1045, t11 * qJD(3) + t123 * qJD(4) + qJD(5) * t801 + t1027 - t577, t46 * qJD(3) - t775 - t914 - t984, t6 * qJD(3) + t17 * qJD(4) + (-pkin(5) * t242 - qJ(6) * t241) * qJD(5) + t231 * qJD(6) + t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t291 - t864, t193, -t1134, qJD(3) * t63 + qJD(4) * t78 + qJD(5) * t231 - t961; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t530 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1105, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t892, -t889, 0, 0, 0, 0, 0, 0, 0, 0, -t852 - t858, t666 - t857 (t680 + t978) * qJD(3) (pkin(8) * t833 - t1053) * qJD(3) - t770, 0, 0, 0, 0, 0, 0, t872, t599 - t943, t1098 (t684 * t696 + t832) * qJD(3) + t65 * qJD(4) + t784, 0, 0, 0, 0, 0, 0, t872, t1098, t443 * t887 - t599 (t451 * t696 + t832) * qJD(3) + t52 * qJD(4) + t54 * qJD(5) - t438 * qJD(6) + t787; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t695 * t889 - t697 * t891, t695 * t891 - t697 * t889, 0, 0, 0, 0, 0, 0, 0, 0, t260, t259, 0, t1026 + t65 * qJD(3) + (-t1060 * t591 - t589 * t694) * t1048 + t328, 0, 0, 0, 0, 0, 0, t260, 0, t141, t1037 + t52 * qJD(3) + (t982 - t987) * qJD(4) + t126 * qJD(5) + t580; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t259, 0, t776, 0, 0, 0, 0, 0, 0, t260, 0, t141, t54 * qJD(3) + t126 * qJD(4) - qJD(5) * t386 + t1038 + t580; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t905 + t939; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t851, t894, 0, t851, 0, 0, -t856, -t855, 0, 0, -t692 * t851 - t897, t746 * t1124, -t857 + t895, -t690 * t851 + t897, t858 + t896, t617, qJD(4) * t502 + t911, qJD(4) * t501 - t910, t913, -t771, t98, t1147, t144, t1132, -t438 * t887 + t1131, t585, -qJD(4) * t73 - qJD(5) * t101 + t786, -qJD(4) * t72 - qJD(5) * t102 - t995, -qJD(4) * t13 - t1020 + t15, -qJD(4) * t7 + t23 - t792, t98, t144, -t1147, t585, -t1131 - t941, t1132, qJD(4) * t42 - qJD(5) * t47 - qJD(6) * t290 - t1014, -qJD(4) * t9 + qJD(5) * t12 - qJD(6) * t441 - t1022, -qJD(4) * t40 - qJD(5) * t45 + t474 + t785, -qJD(4) * t3 - qJD(5) * t5 - qJD(6) * t62 - t796; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t770, 0, 0, 0, 0, 0, 0, -t942, -t944, 0, -qJD(4) * t64 - t784, 0, 0, 0, 0, 0, 0, -t942, 0, -t440 * t887, -qJD(4) * t51 - qJD(5) * t53 - qJD(6) * t436 - t787; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t854, t663 * qJD(4), 0, -t854, 0, 0, -pkin(3) * t918, -pkin(3) * t917, 0, 0, -t1135, t1138, 0, t1135, 0, 0, qJD(4) * t431 + t624 * t915, qJD(4) * t432 - t620 * t915, t1107, qJD(4) * t153, -t1135, 0, -t1138, 0, 0, t1135, qJD(4) * t215 + qJD(5) * t229 - t608 * t624, t1107, qJD(4) * t216 + qJD(5) * t230 + qJD(6) * t619, qJD(4) * t97 + qJD(5) * t100 - qJD(6) * t1017; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t529, -t730, -t823 * t697, -t529, t823 * t695, t676, -pkin(8) * t917 - t758, pkin(8) * t918 - t759, 0, 0, -t1128, t118, t288, t1128, t287, t676, -t778 + t874, -t767 - t1112 (t867 - t990) * t1048 + t783 (-t1060 * t521 + t523 * t694) * t1048 + t727, -t1128, t288, -t118, t676, t289, t1128, -t782 + t874 (-t981 - t985) * qJD(4) + t178 * qJD(5) + t790 - t608, t767 - t1110 (t521 * t683 + t523 * t670) * qJD(4) + t82 * qJD(5) + t391 * qJD(6) + t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1128, t118, t288, t1128, t287, t676, -t745 + t873, -t768 - t1113, t1044, t793, -t1128, t288, -t118, t676, t289, t1128, -t780 + t873, t178 * qJD(4) + qJD(5) * t800 + t1039 - t608, t768 - t1111, t82 * qJD(4) + (-pkin(5) * t524 - qJ(6) * t522) * qJD(5) + t524 * qJD(6) + t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t743, t288, t329, -t716 - t1119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t509, t731 (-t695 * t893 + t919) * t698, -t509 (-t861 - t920) * t698, t674, -qJD(3) * t502 + t928, -qJD(3) * t501 - t929, 0, 0, t1127, -t112, t307, -t1127, t306, t674, qJD(3) * t73 + t86 + t934, qJD(3) * t72 - t933 + t962, qJD(3) * t13 - t1021, qJD(3) * t7 - t791, t1127, t307, t112, t674, t308, -t1127, -qJD(3) * t42 - t1043 + t86, qJD(3) * t9 + qJD(5) * t124 - t1025, t40 * qJD(3) - t914 - t962 + t988, qJD(3) * t3 + qJD(5) * t18 - qJD(6) * t77 - t794; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t399, t403, 0, qJD(3) * t64 - t1026 + t328, 0, 0, 0, 0, 0, 0, t399, 0, t244, qJD(3) * t51 - qJD(5) * t125 - t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t529, t730, t697 * t890, t529, -t695 * t890, t675, t758, t759, 0, 0, t1128, -t118, t402, -t1128, t398, t675, t778 + t945, t1112, -t783, -t727, t1128, t402, t118, t675, t397, -t1128, t782 + t945, qJD(5) * t179 - t790, t1110, qJD(5) * t83 + qJD(6) * t508 - t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t882, -t687, 0, 0, 0, 0, 0, 0, 0, 0, -t882, 0, t888, qJD(5) * t565 + t661; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t687 + t879, 0, t923, 0, 0, 0, 0, 0, 0, t70, t769, -t879 + t888 (-pkin(5) * t694 + qJ(6) * t1060) * t1047 + t661 + t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t651, t670 * t887 + t1114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1127, -t112, t307, -t1127, t306, t674, qJD(3) * t101 - t91 + t931, qJD(3) * t102 - t932 - t963, -t1041, -t1040, t1127, t307, t112, t674, t308, -t1127, qJD(3) * t47 - t1042 - t91, -qJD(3) * t12 - qJD(4) * t124 - t1027, t45 * qJD(3) - t914 + t963 + t984, -qJ(6) * t914 + t5 * qJD(3) - t18 * qJD(4) - t795; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t399, t403, 0, -t776, 0, 0, 0, 0, 0, 0, t399, 0, t244, qJD(3) * t53 + qJD(4) * t125 - t1038; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1128, -t118, t402, -t1128, t398, t675, t745 + t946, t1113, -t1044, -t793, t1128, t402, t118, t675, t397, -t1128, t780 + t946, -qJD(4) * t179 - t1039, t1111, -qJD(4) * t83 - t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t879, 0, -t923, 0, 0, 0, 0, 0, 0, t71, -t769, qJD(6) + t879, t689 - t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), t689; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t651, -t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t290 + t864, t307, t1134, t62 * qJD(3) + t77 * qJD(4) + t698 * t935 + t961; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t906; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t743, t402, -t329, t716 - t497; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t651, -qJD(4) * t670 - t1114 - t935; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t651, t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t30;
