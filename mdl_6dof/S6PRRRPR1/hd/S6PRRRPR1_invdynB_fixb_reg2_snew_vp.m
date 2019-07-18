% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6PRRRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 07:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6PRRRPR1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR1_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_invdynB_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:01:38
% EndTime: 2019-05-05 07:02:23
% DurationCPUTime: 43.94s
% Computational Cost: add. (376636->1017), mult. (812325->1659), div. (0->0), fcn. (619823->14), ass. (0->730)
t1010 = sin(pkin(12));
t1006 = qJDD(3) + qJDD(4);
t1013 = cos(pkin(12));
t1018 = sin(qJ(4));
t1022 = cos(qJ(4));
t1023 = cos(qJ(3));
t1128 = qJD(2) * t1023;
t1019 = sin(qJ(3));
t1129 = qJD(2) * t1019;
t965 = -t1018 * t1129 + t1022 * t1128;
t1055 = t1018 * t1023 + t1019 * t1022;
t966 = t1055 * qJD(2);
t908 = t1010 * t966 - t1013 * t965;
t910 = t1010 * t965 + t1013 * t966;
t842 = t910 * t908;
t1171 = -t842 + t1006;
t1180 = t1010 * t1171;
t1179 = t1013 * t1171;
t1017 = sin(qJ(6));
t1103 = qJD(3) * t1129;
t1113 = qJDD(2) * t1023;
t1047 = t1103 - t1113;
t1102 = qJD(3) * t1128;
t1115 = qJDD(2) * t1019;
t973 = t1102 + t1115;
t1095 = t1018 * t973 + t1022 * t1047;
t885 = -qJD(4) * t966 - t1095;
t886 = t965 * qJD(4) - t1018 * t1047 + t1022 * t973;
t1097 = t1010 * t886 - t1013 * t885;
t1089 = qJDD(6) + t1097;
t1007 = qJD(3) + qJD(4);
t1021 = cos(qJ(6));
t881 = -t1021 * t1007 + t1017 * t910;
t883 = t1007 * t1017 + t1021 * t910;
t817 = t883 * t881;
t1172 = t1089 - t817;
t1178 = t1017 * t1172;
t926 = t965 * t966;
t1170 = t926 + t1006;
t1177 = t1018 * t1170;
t1176 = t1021 * t1172;
t1175 = t1022 * t1170;
t1012 = sin(pkin(6));
t1015 = cos(pkin(6));
t1011 = sin(pkin(11));
t1014 = cos(pkin(11));
t1101 = g(1) * t1011 - t1014 * g(2);
t1159 = g(3) - qJDD(1);
t1174 = -t1012 * t1159 + t1015 * t1101;
t1158 = t1007 * t910;
t776 = t1097 + t1158;
t813 = t1010 * t885 + t1013 * t886;
t899 = t1007 * t908;
t780 = t813 - t899;
t759 = -t881 * qJD(6) + t1017 * t1006 + t1021 * t813;
t903 = qJD(6) + t908;
t827 = t903 * t881;
t725 = -t827 + t759;
t954 = t1007 * t965;
t867 = t954 - t886;
t1173 = t954 + t886;
t1169 = t1011 * t1159;
t1168 = t1014 * t1159;
t982 = g(1) * t1014 + g(2) * t1011;
t932 = -t1011 * t1101 - t1014 * t982;
t863 = (qJD(4) - t1007) * t966 + t1095;
t1096 = -t1021 * t1006 + t1017 * t813;
t722 = (qJD(6) - t903) * t883 + t1096;
t931 = -t1011 * t982 + t1014 * t1101;
t1165 = qJD(2) ^ 2;
t1020 = sin(qJ(2));
t1024 = cos(qJ(2));
t907 = t1020 * t1174 - t1024 * t982;
t891 = -pkin(2) * t1165 + qJDD(2) * pkin(8) + t907;
t950 = t1012 * t1101 + t1015 * t1159;
t854 = t1019 * t891 + t1023 * t950;
t1167 = -t854 + (-t973 + t1102) * pkin(9);
t878 = t881 ^ 2;
t879 = t883 ^ 2;
t902 = t903 ^ 2;
t904 = t908 ^ 2;
t905 = t910 ^ 2;
t1166 = t965 ^ 2;
t964 = t966 ^ 2;
t1164 = t1007 ^ 2;
t1163 = t1023 ^ 2;
t1162 = 2 * qJD(5);
t906 = -t1020 * t982 - t1024 * t1174;
t832 = t1020 * t906 + t1024 * t907;
t1161 = pkin(7) * t832;
t1160 = pkin(5) * t1010;
t1118 = t1023 * t1165;
t993 = t1019 * t1118;
t983 = qJDD(3) + t993;
t1027 = t983 * pkin(3) + t1167;
t1003 = t1163 * t1165;
t855 = -t1019 * t950 + t1023 * t891;
t986 = qJD(3) * pkin(3) - pkin(9) * t1129;
t1029 = -pkin(9) * t1047 - qJD(3) * t986 + t855;
t1028 = -pkin(3) * t1003 + t1029;
t734 = t1018 * t1028 - t1022 * t1027;
t1026 = pkin(4) * t1170 + qJ(5) * t867 - t734;
t735 = t1018 * t1027 + t1022 * t1028;
t946 = pkin(4) * t1007 - qJ(5) * t966;
t697 = -pkin(4) * t1166 + qJ(5) * t885 - t1007 * t946 + t735;
t599 = -0.2e1 * qJD(5) * t908 + t1010 * t1026 + t1013 * t697;
t837 = pkin(5) * t908 - pkin(10) * t910;
t584 = -pkin(5) * t1164 + pkin(10) * t1006 - t837 * t908 + t599;
t890 = -qJDD(2) * pkin(2) - t1165 * pkin(8) + t906;
t841 = t1047 * pkin(3) - pkin(9) * t1003 + t1129 * t986 + t890;
t754 = -t885 * pkin(4) - t1166 * qJ(5) + t946 * t966 + qJDD(5) + t841;
t651 = pkin(5) * t776 - pkin(10) * t780 + t754;
t550 = t1017 * t651 + t1021 * t584;
t1157 = t1010 * t754;
t835 = t842 + t1006;
t1156 = t1010 * t835;
t1154 = t1013 * t754;
t1153 = t1013 * t835;
t1098 = t1010 * t697 - t1013 * t1026;
t583 = -t1006 * pkin(5) - t1164 * pkin(10) + (t1162 + t837) * t910 + t1098;
t1152 = t1017 * t583;
t745 = t1089 + t817;
t1151 = t1017 * t745;
t1150 = t1017 * t903;
t598 = t1162 * t910 + t1098;
t538 = t1010 * t599 - t1013 * t598;
t1149 = t1018 * t538;
t1148 = t1018 * t841;
t917 = -t926 + t1006;
t1147 = t1018 * t917;
t654 = t1018 * t735 - t1022 * t734;
t1146 = t1019 * t654;
t1145 = t1019 * t890;
t1144 = t1019 * t983;
t984 = qJDD(3) - t993;
t1143 = t1019 * t984;
t1142 = t1020 * t950;
t1141 = t1021 * t583;
t1140 = t1021 * t745;
t1139 = t1021 * t903;
t1138 = t1022 * t538;
t1137 = t1022 * t841;
t1136 = t1022 * t917;
t1135 = t1023 * t654;
t1134 = t1023 * t890;
t974 = -0.2e1 * t1103 + t1113;
t933 = t1023 * t974;
t1133 = t1023 * t984;
t1132 = t1024 * t950;
t1130 = qJD(2) * qJD(3);
t1008 = t1019 ^ 2;
t1127 = t1165 * t1008;
t1126 = t1006 * t1020;
t1125 = t1006 * t1024;
t1124 = t1007 * t1010;
t1123 = t1007 * t1013;
t1122 = t1007 * t1018;
t1121 = t1007 * t1022;
t1120 = t1012 * t1020;
t1119 = t1015 * t1020;
t1116 = qJDD(2) * t1012;
t1114 = qJDD(2) * t1020;
t1112 = qJDD(2) * t1024;
t1111 = t1008 + t1163;
t1110 = t1010 * t817;
t1109 = t1013 * t817;
t1108 = t1020 * t842;
t1107 = t1020 * t926;
t1106 = t1024 * t842;
t1105 = t1024 * t926;
t1104 = -pkin(5) * t1013 - pkin(4);
t539 = t1010 * t598 + t1013 * t599;
t549 = t1017 * t584 - t1021 * t651;
t655 = t1018 * t734 + t1022 * t735;
t773 = t1019 * t854 + t1023 * t855;
t1094 = t1020 * t993;
t1093 = t1024 * t993;
t770 = t1019 * t855 - t1023 * t854;
t975 = t1111 * qJDD(2);
t978 = t1003 + t1127;
t929 = -t1020 * t978 + t1024 * t975;
t1092 = pkin(7) * t929 - t1020 * t770;
t976 = -t1020 * t1165 + t1112;
t1091 = -pkin(7) * t976 - t1142;
t1054 = t1024 * t1165 + t1114;
t1090 = -pkin(7) * t1054 + t1132;
t961 = t1054 * t1015;
t1087 = t1011 * t976 + t1014 * t961;
t921 = t1011 * t961 - t1014 * t976;
t489 = t1017 * t550 - t1021 * t549;
t490 = t1017 * t549 + t1021 * t550;
t468 = t1010 * t490 - t1013 * t583;
t469 = t1010 * t583 + t1013 * t490;
t439 = t1018 * t469 + t1022 * t468;
t440 = -t1018 * t468 + t1022 * t469;
t412 = -t1019 * t439 + t1023 * t440;
t1086 = t1020 * t412 - t1024 * t489;
t484 = t1018 * t539 + t1138;
t485 = t1022 * t539 - t1149;
t450 = -t1019 * t484 + t1023 * t485;
t1085 = t1020 * t450 - t1024 * t754;
t726 = -t827 - t759;
t644 = -t1017 * t726 - t1021 * t722;
t775 = t878 + t879;
t600 = t1010 * t644 + t1013 * t775;
t601 = -t1010 * t775 + t1013 * t644;
t542 = t1018 * t601 + t1022 * t600;
t543 = -t1018 * t600 + t1022 * t601;
t487 = -t1019 * t542 + t1023 * t543;
t642 = -t1017 * t722 + t1021 * t726;
t1084 = t1020 * t487 - t1024 * t642;
t724 = (-qJD(6) - t903) * t883 - t1096;
t645 = -t1017 * t725 + t1021 * t724;
t816 = -t879 + t878;
t613 = t1010 * t645 + t1013 * t816;
t614 = -t1010 * t816 + t1013 * t645;
t553 = t1018 * t614 + t1022 * t613;
t554 = -t1018 * t613 + t1022 * t614;
t497 = -t1019 * t553 + t1023 * t554;
t643 = -t1017 * t724 - t1021 * t725;
t1083 = t1020 * t497 + t1024 * t643;
t788 = -t902 - t878;
t678 = t1021 * t788 - t1178;
t615 = t1010 * t678 + t1013 * t724;
t616 = -t1010 * t724 + t1013 * t678;
t556 = t1018 * t616 + t1022 * t615;
t557 = -t1018 * t615 + t1022 * t616;
t499 = -t1019 * t556 + t1023 * t557;
t677 = t1017 * t788 + t1176;
t1082 = t1020 * t499 - t1024 * t677;
t799 = -t879 - t902;
t687 = -t1017 * t799 - t1140;
t621 = t1010 * t687 - t1013 * t725;
t622 = t1010 * t725 + t1013 * t687;
t563 = t1018 * t622 + t1022 * t621;
t564 = -t1018 * t621 + t1022 * t622;
t503 = -t1019 * t563 + t1023 * t564;
t686 = t1021 * t799 - t1151;
t1081 = t1020 * t503 - t1024 * t686;
t826 = -t879 + t902;
t700 = -t1017 * t826 + t1176;
t633 = t1010 * t700 + t1013 * t726;
t635 = -t1010 * t726 + t1013 * t700;
t569 = t1018 * t635 + t1022 * t633;
t571 = -t1018 * t633 + t1022 * t635;
t507 = -t1019 * t569 + t1023 * t571;
t698 = -t1021 * t826 - t1178;
t1080 = t1020 * t507 + t1024 * t698;
t825 = t878 - t902;
t701 = t1021 * t825 - t1151;
t634 = t1010 * t701 + t1013 * t722;
t636 = -t1010 * t722 + t1013 * t701;
t570 = t1018 * t636 + t1022 * t634;
t572 = -t1018 * t634 + t1022 * t636;
t508 = -t1019 * t570 + t1023 * t572;
t699 = -t1017 * t825 - t1140;
t1079 = t1020 * t508 + t1024 * t699;
t758 = -qJD(6) * t883 - t1096;
t719 = -t1017 * t758 + t1139 * t881;
t668 = t1010 * t719 + t1109;
t670 = t1013 * t719 - t1110;
t586 = t1018 * t670 + t1022 * t668;
t588 = -t1018 * t668 + t1022 * t670;
t534 = -t1019 * t586 + t1023 * t588;
t718 = -t1021 * t758 - t1150 * t881;
t1078 = t1020 * t534 + t1024 * t718;
t721 = t1021 * t759 - t1150 * t883;
t669 = t1010 * t721 - t1109;
t671 = t1013 * t721 + t1110;
t587 = t1018 * t671 + t1022 * t669;
t589 = -t1018 * t669 + t1022 * t671;
t535 = -t1019 * t587 + t1023 * t589;
t720 = -t1017 * t759 - t1139 * t883;
t1077 = t1020 * t535 + t1024 * t720;
t702 = -t1010 * t776 + t1013 * t780;
t704 = -t1010 * t780 - t1013 * t776;
t617 = t1018 * t704 + t1022 * t702;
t619 = -t1018 * t702 + t1022 * t704;
t561 = -t1019 * t617 + t1023 * t619;
t838 = -t905 + t904;
t1076 = t1020 * t561 + t1024 * t838;
t778 = -t1097 + t1158;
t781 = -t813 - t899;
t703 = t1010 * t778 + t1013 * t781;
t705 = -t1010 * t781 + t1013 * t778;
t618 = t1018 * t705 + t1022 * t703;
t620 = -t1018 * t703 + t1022 * t705;
t562 = -t1019 * t618 + t1023 * t620;
t804 = -t904 - t905;
t1075 = t1020 * t562 - t1024 * t804;
t753 = (t1017 * t883 - t1021 * t881) * t903;
t706 = t1010 * t753 - t1013 * t1089;
t707 = t1010 * t1089 + t1013 * t753;
t623 = t1018 * t707 + t1022 * t706;
t624 = -t1018 * t706 + t1022 * t707;
t566 = -t1019 * t623 + t1023 * t624;
t752 = (t1017 * t881 + t1021 * t883) * t903;
t1074 = t1020 * t566 + t1024 * t752;
t580 = t1023 * t655 - t1146;
t1073 = t1020 * t580 - t1024 * t841;
t833 = -t1164 - t904;
t756 = t1010 * t833 + t1179;
t757 = t1013 * t833 - t1180;
t679 = t1018 * t757 + t1022 * t756;
t680 = -t1018 * t756 + t1022 * t757;
t594 = -t1019 * t679 + t1023 * t680;
t1072 = t1020 * t594 - t1024 * t776;
t889 = -t905 - t1164;
t789 = t1013 * t889 - t1156;
t790 = -t1010 * t889 - t1153;
t712 = t1018 * t790 + t1022 * t789;
t713 = -t1018 * t789 + t1022 * t790;
t632 = -t1019 * t712 + t1023 * t713;
t1071 = t1020 * t632 - t1024 * t780;
t895 = -t905 + t1164;
t791 = t1013 * t895 + t1180;
t795 = -t1010 * t895 + t1179;
t714 = t1018 * t795 + t1022 * t791;
t716 = -t1018 * t791 + t1022 * t795;
t640 = -t1019 * t714 + t1023 * t716;
t1070 = t1020 * t640 + t1024 * t781;
t894 = t904 - t1164;
t792 = t1010 * t894 + t1153;
t796 = t1013 * t894 - t1156;
t715 = t1018 * t796 + t1022 * t792;
t717 = -t1018 * t792 + t1022 * t796;
t641 = -t1019 * t715 + t1023 * t717;
t1069 = t1020 * t641 - t1024 * t778;
t862 = (qJD(4) + t1007) * t966 + t1095;
t782 = -t1018 * t862 + t1022 * t1173;
t784 = -t1018 * t1173 - t1022 * t862;
t710 = -t1019 * t782 + t1023 * t784;
t925 = -t964 + t1166;
t1068 = t1020 * t710 + t1024 * t925;
t783 = -t1018 * t863 + t1022 * t867;
t785 = -t1018 * t867 - t1022 * t863;
t711 = -t1019 * t783 + t1023 * t785;
t888 = -t964 - t1166;
t1067 = t1020 * t711 - t1024 * t888;
t911 = -t1164 - t1166;
t839 = t1018 * t911 + t1175;
t840 = t1022 * t911 - t1177;
t761 = -t1019 * t839 + t1023 * t840;
t1066 = t1020 * t761 - t1024 * t862;
t1065 = t1020 * t773 - t1024 * t890;
t943 = -t964 - t1164;
t868 = t1022 * t943 - t1147;
t869 = -t1018 * t943 - t1136;
t787 = -t1019 * t868 + t1023 * t869;
t1064 = t1020 * t787 - t1024 * t1173;
t952 = -t964 + t1164;
t872 = t1022 * t952 + t1177;
t874 = -t1018 * t952 + t1175;
t797 = -t1019 * t872 + t1023 * t874;
t1063 = t1020 * t797 + t1024 * t867;
t951 = -t1164 + t1166;
t873 = t1018 * t951 + t1136;
t875 = t1022 * t951 - t1147;
t798 = -t1019 * t873 + t1023 * t875;
t1062 = t1020 * t798 + t1024 * t863;
t1061 = t1020 * t907 - t1024 * t906;
t972 = 0.2e1 * t1102 + t1115;
t928 = -t1019 * t972 + t933;
t979 = t1003 - t1127;
t1060 = t1020 * t928 + t1024 * t979;
t1025 = qJD(3) ^ 2;
t991 = -t1003 - t1025;
t940 = t1023 * t991 - t1144;
t1059 = t1020 * t940 + t1024 * t974;
t989 = -t1025 - t1127;
t942 = -t1019 * t989 - t1133;
t1058 = t1020 * t942 - t1024 * t972;
t1057 = t1020 * t975 + t1024 * t978;
t970 = t1111 * t1130;
t1056 = -qJDD(3) * t1024 + t1020 * t970;
t767 = -t1013 * t1097 + t1124 * t908;
t768 = t1010 * t1097 + t1123 * t908;
t692 = t1018 * t768 + t1022 * t767;
t694 = -t1018 * t767 + t1022 * t768;
t608 = -t1019 * t692 + t1023 * t694;
t1053 = t1020 * t608 + t1106;
t769 = t1010 * t813 + t1123 * t910;
t771 = t1013 * t813 - t1124 * t910;
t693 = t1018 * t771 + t1022 * t769;
t695 = -t1018 * t769 + t1022 * t771;
t609 = -t1019 * t693 + t1023 * t695;
t1052 = t1020 * t609 - t1106;
t850 = t1022 * t885 - t1122 * t965;
t851 = -t1018 * t885 - t1121 * t965;
t765 = -t1019 * t850 + t1023 * t851;
t1051 = t1020 * t765 - t1105;
t852 = t1018 * t886 + t1121 * t966;
t853 = t1022 * t886 - t1122 * t966;
t766 = -t1019 * t852 + t1023 * t853;
t1050 = t1020 * t766 + t1105;
t971 = t1023 * t983;
t988 = t1025 - t1127;
t941 = -t1019 * t988 + t971;
t1049 = -t1019 * t1112 + t1020 * t941;
t990 = t1003 - t1025;
t939 = t1023 * t990 - t1143;
t1048 = t1020 * t939 - t1023 * t1112;
t413 = qJ(5) * t469 + (-pkin(10) * t1010 + t1104) * t489;
t427 = -qJ(5) * t468 + (-pkin(10) * t1013 + t1160) * t489;
t395 = -pkin(3) * t489 + pkin(9) * t440 + t1018 * t427 + t1022 * t413;
t397 = -pkin(9) * t439 - t1018 * t413 + t1022 * t427;
t411 = t1019 * t440 + t1023 * t439;
t382 = -pkin(8) * t411 - t1019 * t395 + t1023 * t397;
t396 = -pkin(2) * t411 - pkin(3) * t439 - pkin(4) * t468 + pkin(5) * t583 - pkin(10) * t490;
t410 = t1020 * t489 + t1024 * t412;
t1046 = pkin(7) * t410 + t1020 * t382 + t1024 * t396;
t476 = -pkin(10) * t642 - t489;
t457 = qJ(5) * t601 + t1010 * t476 + t1104 * t642;
t459 = -qJ(5) * t600 + t1013 * t476 + t1160 * t642;
t417 = -pkin(3) * t642 + pkin(9) * t543 + t1018 * t459 + t1022 * t457;
t418 = -pkin(9) * t542 - t1018 * t457 + t1022 * t459;
t486 = t1019 * t543 + t1023 * t542;
t400 = -pkin(8) * t486 - t1019 * t417 + t1023 * t418;
t435 = -pkin(2) * t486 - pkin(3) * t542 - pkin(4) * t600 - pkin(5) * t775 - pkin(10) * t644 - t490;
t475 = t1020 * t642 + t1024 * t487;
t1045 = pkin(7) * t475 + t1020 * t400 + t1024 * t435;
t525 = -pkin(5) * t677 + t549;
t567 = -pkin(10) * t677 + t1152;
t464 = -pkin(4) * t677 + qJ(5) * t616 + t1010 * t567 + t1013 * t525;
t467 = -qJ(5) * t615 - t1010 * t525 + t1013 * t567;
t425 = -pkin(3) * t677 + pkin(9) * t557 + t1018 * t467 + t1022 * t464;
t428 = -pkin(9) * t556 - t1018 * t464 + t1022 * t467;
t498 = t1019 * t557 + t1023 * t556;
t403 = -pkin(8) * t498 - t1019 * t425 + t1023 * t428;
t444 = -pkin(2) * t498 - pkin(3) * t556 - pkin(4) * t615 - pkin(5) * t724 - pkin(10) * t678 + t1141;
t480 = t1020 * t677 + t1024 * t499;
t1044 = pkin(7) * t480 + t1020 * t403 + t1024 * t444;
t526 = -pkin(5) * t686 + t550;
t568 = -pkin(10) * t686 + t1141;
t465 = -pkin(4) * t686 + qJ(5) * t622 + t1010 * t568 + t1013 * t526;
t470 = -qJ(5) * t621 - t1010 * t526 + t1013 * t568;
t426 = -pkin(3) * t686 + pkin(9) * t564 + t1018 * t470 + t1022 * t465;
t430 = -pkin(9) * t563 - t1018 * t465 + t1022 * t470;
t502 = t1019 * t564 + t1023 * t563;
t404 = -pkin(8) * t502 - t1019 * t426 + t1023 * t430;
t445 = -pkin(2) * t502 - pkin(3) * t563 - pkin(4) * t621 + pkin(5) * t725 - pkin(10) * t687 - t1152;
t483 = t1020 * t686 + t1024 * t503;
t1043 = pkin(7) * t483 + t1020 * t404 + t1024 * t445;
t530 = -pkin(4) * t754 + qJ(5) * t539;
t443 = -pkin(3) * t754 + pkin(9) * t485 - qJ(5) * t1149 + t1022 * t530;
t447 = -pkin(9) * t484 - qJ(5) * t1138 - t1018 * t530;
t449 = t1019 * t485 + t1023 * t484;
t406 = -pkin(8) * t449 - t1019 * t443 + t1023 * t447;
t421 = -pkin(2) * t449 - pkin(3) * t484 - pkin(4) * t538;
t446 = t1020 * t754 + t1024 * t450;
t1042 = pkin(7) * t446 + t1020 * t406 + t1024 * t421;
t519 = -pkin(4) * t804 + qJ(5) * t705 + t539;
t524 = -qJ(5) * t703 - t538;
t460 = -pkin(3) * t804 + pkin(9) * t620 + t1018 * t524 + t1022 * t519;
t461 = -pkin(9) * t618 - t1018 * t519 + t1022 * t524;
t560 = t1019 * t620 + t1023 * t618;
t422 = -pkin(8) * t560 - t1019 * t460 + t1023 * t461;
t509 = -pkin(2) * t560 - pkin(3) * t618 - pkin(4) * t703;
t544 = t1020 * t804 + t1024 * t562;
t1041 = pkin(7) * t544 + t1020 * t422 + t1024 * t509;
t648 = -pkin(4) * t776 + qJ(5) * t757 - t1154;
t674 = -qJ(5) * t756 + t1157;
t537 = -pkin(3) * t776 + pkin(9) * t680 + t1018 * t674 + t1022 * t648;
t558 = -pkin(9) * t679 - t1018 * t648 + t1022 * t674;
t593 = t1019 * t680 + t1023 * t679;
t471 = -pkin(8) * t593 - t1019 * t537 + t1023 * t558;
t514 = -pkin(2) * t593 - pkin(3) * t679 - pkin(4) * t756 + t598;
t581 = t1020 * t776 + t1024 * t594;
t1040 = pkin(7) * t581 + t1020 * t471 + t1024 * t514;
t652 = -pkin(4) * t780 + qJ(5) * t790 + t1157;
t689 = -qJ(5) * t789 + t1154;
t555 = -pkin(3) * t780 + pkin(9) * t713 + t1018 * t689 + t1022 * t652;
t574 = -pkin(9) * t712 - t1018 * t652 + t1022 * t689;
t631 = t1019 * t713 + t1023 * t712;
t479 = -pkin(8) * t631 - t1019 * t555 + t1023 * t574;
t523 = -pkin(2) * t631 - pkin(3) * t712 - pkin(4) * t789 + t599;
t590 = t1020 * t780 + t1024 * t632;
t1039 = pkin(7) * t590 + t1020 * t479 + t1024 * t523;
t579 = t1019 * t655 + t1135;
t646 = -pkin(3) * t841 + pkin(9) * t655;
t518 = -pkin(8) * t579 - pkin(9) * t1135 - t1019 * t646;
t545 = -pkin(2) * t579 - pkin(3) * t654;
t578 = t1020 * t841 + t1024 * t580;
t1038 = pkin(7) * t578 + t1020 * t518 + t1024 * t545;
t612 = -pkin(3) * t888 + pkin(9) * t785 + t655;
t630 = -pkin(9) * t783 - t654;
t709 = t1019 * t785 + t1023 * t783;
t531 = -pkin(8) * t709 - t1019 * t612 + t1023 * t630;
t656 = -pkin(2) * t709 - pkin(3) * t783;
t683 = t1020 * t888 + t1024 * t711;
t1037 = pkin(7) * t683 + t1020 * t531 + t1024 * t656;
t731 = -pkin(3) * t862 + pkin(9) * t840 - t1137;
t760 = t1019 * t840 + t1023 * t839;
t762 = -pkin(9) * t839 + t1148;
t637 = -pkin(8) * t760 - t1019 * t731 + t1023 * t762;
t653 = -pkin(2) * t760 + t1018 * t1029 - t1022 * t1167 + (-t1022 * qJDD(3) - t1055 * t1118 - t839) * pkin(3);
t729 = t1020 * t862 + t1024 * t761;
t1036 = pkin(7) * t729 + t1020 * t637 + t1024 * t653;
t732 = -pkin(3) * t1173 + pkin(9) * t869 + t1148;
t774 = -pkin(9) * t868 + t1137;
t786 = t1019 * t869 + t1023 * t868;
t647 = -pkin(8) * t786 - t1019 * t732 + t1023 * t774;
t657 = -pkin(2) * t786 - pkin(3) * t868 + t735;
t736 = t1020 * t1173 + t1024 * t787;
t1035 = pkin(7) * t736 + t1020 * t647 + t1024 * t657;
t936 = t1019 * t991 + t971;
t818 = -pkin(2) * t936 + t854;
t848 = -pkin(8) * t936 + t1145;
t896 = -t1020 * t974 + t1024 * t940;
t1034 = pkin(7) * t896 + t1020 * t848 + t1024 * t818;
t938 = t1023 * t989 - t1143;
t819 = -pkin(2) * t938 + t855;
t849 = -pkin(8) * t938 + t1134;
t897 = t1020 * t972 + t1024 * t942;
t1033 = pkin(7) * t897 + t1020 * t849 + t1024 * t819;
t947 = t1019 * t1047 - t1130 * t1163;
t1032 = t1020 * t947 - t1093;
t948 = -t1008 * t1130 + t1023 * t973;
t1031 = t1020 * t948 + t1093;
t743 = t1020 * t890 + t1024 * t773;
t1030 = pkin(7) * t743 + (-pkin(2) * t1024 - pkin(8) * t1020) * t770;
t981 = t1015 * t1125;
t980 = t1012 * t1125;
t962 = t976 * t1015;
t960 = t976 * t1012;
t959 = t1054 * t1012;
t949 = qJDD(3) * t1020 + t1024 * t970;
t937 = t1023 * t988 + t1144;
t935 = t1019 * t990 + t1133;
t934 = (t973 + t1102) * t1019;
t930 = t1056 * t1015;
t927 = t1019 * t974 + t1023 * t972;
t924 = t1057 * t1015;
t923 = t1057 * t1012;
t922 = -t1011 * t962 - t1014 * t1054;
t920 = -t1011 * t1054 + t1014 * t962;
t915 = t1024 * t948 - t1094;
t914 = t1024 * t947 + t1094;
t913 = t1019 * t1114 + t1024 * t941;
t912 = t1020 * t1113 + t1024 * t939;
t893 = (t1018 * t966 + t1022 * t965) * t1007;
t892 = (t1018 * t965 - t1022 * t966) * t1007;
t887 = -t1020 * t979 + t1024 * t928;
t871 = -t1132 + (t1012 * t959 + t1015 * t961) * pkin(7);
t870 = -t1142 + (-t1012 * t960 - t1015 * t962) * pkin(7);
t861 = -t1011 * t924 + t1014 * t929;
t860 = t1011 * t929 + t1014 * t924;
t859 = -t1012 * t934 + t1015 * t1031;
t858 = -t1012 * t933 + t1015 * t1032;
t857 = -t1012 * t937 + t1015 * t1049;
t856 = -t1012 * t935 + t1015 * t1048;
t846 = -t1012 * t938 + t1015 * t1058;
t845 = -t1012 * t936 + t1015 * t1059;
t844 = t1012 * t1058 + t1015 * t938;
t843 = t1012 * t1059 + t1015 * t936;
t830 = -t1012 * t927 + t1015 * t1060;
t829 = pkin(2) * t974 + pkin(8) * t940 - t1134;
t828 = -pkin(2) * t972 + pkin(8) * t942 + t1145;
t824 = t832 * t1015;
t823 = (t1010 * t910 - t1013 * t908) * t1007;
t822 = (-t1010 * t908 - t1013 * t910) * t1007;
t821 = -t1019 * t892 + t1023 * t893;
t820 = t1019 * t893 + t1023 * t892;
t815 = -pkin(1) * t960 + t1012 * t906 + t1015 * t1090;
t814 = pkin(1) * t959 + t1012 * t907 + t1015 * t1091;
t807 = t1024 * t821 + t1126;
t806 = t1012 * t950 + t1015 * t1061;
t805 = t1012 * t1061 - t1015 * t950;
t803 = -t1011 * t846 + t1014 * t897;
t802 = -t1011 * t845 + t1014 * t896;
t801 = t1011 * t897 + t1014 * t846;
t800 = t1011 * t896 + t1014 * t845;
t794 = t1019 * t875 + t1023 * t873;
t793 = t1019 * t874 + t1023 * t872;
t764 = t1019 * t853 + t1023 * t852;
t763 = t1019 * t851 + t1023 * t850;
t755 = pkin(2) * t978 + pkin(8) * t975 + t773;
t751 = -t1018 * t822 + t1022 * t823;
t750 = t1018 * t823 + t1022 * t822;
t749 = t1024 * t766 - t1107;
t748 = t1024 * t765 + t1107;
t747 = -pkin(2) * t890 + pkin(8) * t773;
t742 = -pkin(1) * t805 + t1015 * t1161;
t741 = -t1020 * t863 + t1024 * t798;
t740 = -t1020 * t867 + t1024 * t797;
t739 = -t1011 * t806 + t1014 * t832;
t738 = t1011 * t832 + t1014 * t806;
t737 = -t1012 * t820 + t1119 * t821 - t981;
t730 = -t1024 * t770 + (-t1012 * t923 - t1015 * t924) * pkin(7);
t728 = (-t1012 * t805 - t1015 * t806) * pkin(7);
t708 = t1019 * t784 + t1023 * t782;
t688 = -t1020 * t925 + t1024 * t710;
t685 = -t1020 * t819 + t1024 * t849 + (-t1012 * t844 - t1015 * t846) * pkin(7);
t684 = -t1020 * t818 + t1024 * t848 + (-t1012 * t843 - t1015 * t845) * pkin(7);
t682 = -t1012 * t794 + t1015 * t1062;
t681 = -t1012 * t793 + t1015 * t1063;
t676 = -t1012 * t764 + t1015 * t1050;
t675 = -t1012 * t763 + t1015 * t1051;
t673 = -t1012 * t786 + t1015 * t1064;
t672 = t1012 * t1064 + t1015 * t786;
t667 = -t1012 * t770 + t1015 * t1065;
t666 = t1012 * t1065 + t1015 * t770;
t665 = -t1019 * t750 + t1023 * t751;
t664 = t1019 * t751 + t1023 * t750;
t663 = t1024 * t665 + t1126;
t662 = -pkin(1) * t844 - t1012 * t828 + t1015 * t1033;
t661 = -pkin(1) * t843 - t1012 * t829 + t1015 * t1034;
t660 = -pkin(1) * t923 - t1012 * t755 + t1015 * t1092;
t659 = -t1012 * t760 + t1015 * t1066;
t658 = t1012 * t1066 + t1015 * t760;
t639 = t1019 * t717 + t1023 * t715;
t638 = t1019 * t716 + t1023 * t714;
t629 = -t1011 * t667 + t1014 * t743;
t628 = t1011 * t743 + t1014 * t667;
t627 = -pkin(2) * t1173 + pkin(8) * t787 + t1019 * t774 + t1023 * t732;
t626 = -t1011 * t673 + t1014 * t736;
t625 = t1011 * t736 + t1014 * t673;
t611 = -t1012 * t708 + t1015 * t1068;
t610 = -pkin(2) * t862 + pkin(8) * t761 + t1019 * t762 + t1023 * t731;
t607 = t1019 * t695 + t1023 * t693;
t606 = t1019 * t694 + t1023 * t692;
t605 = -t1012 * t709 + t1015 * t1067;
t604 = t1012 * t1067 + t1015 * t709;
t603 = -t1011 * t659 + t1014 * t729;
t602 = t1011 * t729 + t1014 * t659;
t596 = t1020 * t778 + t1024 * t641;
t595 = -t1020 * t781 + t1024 * t640;
t592 = t1024 * t609 + t1108;
t591 = t1024 * t608 - t1108;
t585 = -t1012 * t664 + t1119 * t665 - t981;
t577 = -t1011 * t605 + t1014 * t683;
t576 = t1011 * t683 + t1014 * t605;
t575 = (pkin(2) * t1020 - pkin(8) * t1024) * t770 + (-t1012 * t666 - t1015 * t667) * pkin(7);
t573 = -pkin(1) * t666 - t1012 * t747 + t1015 * t1030;
t565 = t1019 * t624 + t1023 * t623;
t559 = t1019 * t619 + t1023 * t617;
t552 = -t1012 * t639 + t1015 * t1069;
t551 = -t1012 * t638 + t1015 * t1070;
t548 = -t1012 * t631 + t1015 * t1071;
t547 = t1012 * t1071 + t1015 * t631;
t546 = -t1020 * t838 + t1024 * t561;
t541 = -t1012 * t607 + t1015 * t1052;
t540 = -t1012 * t606 + t1015 * t1053;
t536 = -t1020 * t752 + t1024 * t566;
t533 = t1019 * t589 + t1023 * t587;
t532 = t1019 * t588 + t1023 * t586;
t529 = -pkin(2) * t888 + pkin(8) * t711 + t1019 * t630 + t1023 * t612;
t528 = -t1012 * t593 + t1015 * t1072;
t527 = t1012 * t1072 + t1015 * t593;
t522 = -t1020 * t657 + t1024 * t647 + (-t1012 * t672 - t1015 * t673) * pkin(7);
t521 = -t1020 * t720 + t1024 * t535;
t520 = -t1020 * t718 + t1024 * t534;
t517 = -t1012 * t579 + t1015 * t1073;
t516 = t1012 * t1073 + t1015 * t579;
t515 = -t1020 * t653 + t1024 * t637 + (-t1012 * t658 - t1015 * t659) * pkin(7);
t513 = -pkin(2) * t841 + pkin(8) * t580 - pkin(9) * t1146 + t1023 * t646;
t512 = -t1011 * t548 + t1014 * t590;
t511 = t1011 * t590 + t1014 * t548;
t510 = -pkin(1) * t672 - t1012 * t627 + t1015 * t1035;
t506 = t1019 * t572 + t1023 * t570;
t505 = t1019 * t571 + t1023 * t569;
t504 = -pkin(1) * t658 - t1012 * t610 + t1015 * t1036;
t501 = -t1011 * t528 + t1014 * t581;
t500 = t1011 * t581 + t1014 * t528;
t496 = t1019 * t554 + t1023 * t553;
t495 = -t1012 * t559 + t1015 * t1076;
t494 = -t1012 * t560 + t1015 * t1075;
t493 = t1012 * t1075 + t1015 * t560;
t492 = -t1020 * t699 + t1024 * t508;
t491 = -t1020 * t698 + t1024 * t507;
t488 = -t1012 * t565 + t1015 * t1074;
t482 = -t1011 * t517 + t1014 * t578;
t481 = t1011 * t578 + t1014 * t517;
t478 = -t1020 * t643 + t1024 * t497;
t477 = -pkin(2) * t780 + pkin(8) * t632 + t1019 * t574 + t1023 * t555;
t474 = -t1020 * t656 + t1024 * t531 + (-t1012 * t604 - t1015 * t605) * pkin(7);
t473 = -t1012 * t533 + t1015 * t1077;
t472 = -t1012 * t532 + t1015 * t1078;
t466 = -pkin(2) * t776 + pkin(8) * t594 + t1019 * t558 + t1023 * t537;
t463 = -t1011 * t494 + t1014 * t544;
t462 = t1011 * t544 + t1014 * t494;
t458 = -pkin(1) * t604 - t1012 * t529 + t1015 * t1037;
t456 = -t1012 * t506 + t1015 * t1079;
t455 = -t1012 * t505 + t1015 * t1080;
t454 = -t1012 * t502 + t1015 * t1081;
t453 = t1012 * t1081 + t1015 * t502;
t452 = -t1012 * t498 + t1015 * t1082;
t451 = t1012 * t1082 + t1015 * t498;
t448 = -t1012 * t496 + t1015 * t1083;
t442 = -t1012 * t486 + t1015 * t1084;
t441 = t1012 * t1084 + t1015 * t486;
t438 = -t1020 * t545 + t1024 * t518 + (-t1012 * t516 - t1015 * t517) * pkin(7);
t437 = -t1020 * t523 + t1024 * t479 + (-t1012 * t547 - t1015 * t548) * pkin(7);
t436 = -pkin(1) * t516 - t1012 * t513 + t1015 * t1038;
t434 = -t1011 * t454 + t1014 * t483;
t433 = t1011 * t483 + t1014 * t454;
t432 = -t1011 * t452 + t1014 * t480;
t431 = t1011 * t480 + t1014 * t452;
t429 = -t1020 * t514 + t1024 * t471 + (-t1012 * t527 - t1015 * t528) * pkin(7);
t424 = -t1011 * t442 + t1014 * t475;
t423 = t1011 * t475 + t1014 * t442;
t420 = -pkin(2) * t804 + pkin(8) * t562 + t1019 * t461 + t1023 * t460;
t419 = -pkin(1) * t547 - t1012 * t477 + t1015 * t1039;
t416 = -t1012 * t449 + t1015 * t1085;
t415 = t1012 * t1085 + t1015 * t449;
t414 = -pkin(1) * t527 - t1012 * t466 + t1015 * t1040;
t409 = -t1020 * t509 + t1024 * t422 + (-t1012 * t493 - t1015 * t494) * pkin(7);
t408 = -t1011 * t416 + t1014 * t446;
t407 = t1011 * t446 + t1014 * t416;
t405 = -pkin(2) * t754 + pkin(8) * t450 + t1019 * t447 + t1023 * t443;
t402 = -pkin(2) * t686 + pkin(8) * t503 + t1019 * t430 + t1023 * t426;
t401 = -pkin(2) * t677 + pkin(8) * t499 + t1019 * t428 + t1023 * t425;
t399 = -pkin(2) * t642 + pkin(8) * t487 + t1019 * t418 + t1023 * t417;
t398 = -pkin(1) * t493 - t1012 * t420 + t1015 * t1041;
t394 = -t1012 * t411 + t1015 * t1086;
t393 = t1012 * t1086 + t1015 * t411;
t392 = -t1020 * t445 + t1024 * t404 + (-t1012 * t453 - t1015 * t454) * pkin(7);
t391 = -t1020 * t444 + t1024 * t403 + (-t1012 * t451 - t1015 * t452) * pkin(7);
t390 = -t1020 * t435 + t1024 * t400 + (-t1012 * t441 - t1015 * t442) * pkin(7);
t389 = -t1011 * t394 + t1014 * t410;
t388 = t1011 * t410 + t1014 * t394;
t387 = -t1020 * t421 + t1024 * t406 + (-t1012 * t415 - t1015 * t416) * pkin(7);
t386 = -pkin(1) * t453 - t1012 * t402 + t1015 * t1043;
t385 = -pkin(1) * t451 - t1012 * t401 + t1015 * t1044;
t384 = -pkin(1) * t415 - t1012 * t405 + t1015 * t1042;
t383 = -pkin(1) * t441 - t1012 * t399 + t1015 * t1045;
t381 = -pkin(2) * t489 + pkin(8) * t412 + t1019 * t397 + t1023 * t395;
t380 = -t1020 * t396 + t1024 * t382 + (-t1012 * t393 - t1015 * t394) * pkin(7);
t379 = -pkin(1) * t393 - t1012 * t381 + t1015 * t1046;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t932, 0, 0, 0, 0, 0, 0, t922, t921, 0, t739, 0, 0, 0, 0, 0, 0, t802, t803, t861, t629, 0, 0, 0, 0, 0, 0, t603, t626, t577, t482, 0, 0, 0, 0, 0, 0, t501, t512, t463, t408, 0, 0, 0, 0, 0, 0, t432, t434, t424, t389; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t931, 0, 0, 0, 0, 0, 0, t920, -t1087, 0, t738, 0, 0, 0, 0, 0, 0, t800, t801, t860, t628, 0, 0, 0, 0, 0, 0, t602, t625, t576, t481, 0, 0, 0, 0, 0, 0, t500, t511, t462, t407, 0, 0, 0, 0, 0, 0, t431, t433, t423, t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1159, 0, 0, 0, 0, 0, 0, t960, -t959, 0, t805, 0, 0, 0, 0, 0, 0, t843, t844, t923, t666, 0, 0, 0, 0, 0, 0, t658, t672, t604, t516, 0, 0, 0, 0, 0, 0, t527, t547, t493, t415, 0, 0, 0, 0, 0, 0, t451, t453, t441, t393; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t1169, -t1168, -t931, -qJ(1) * t931, 0, 0, -t921, 0, t922, t1011 * t1116, -qJ(1) * t920 - t1011 * t815 + t1014 * t870, qJ(1) * t1087 - t1011 * t814 + t1014 * t871, -t1011 * t824 - t1014 * t1061, -qJ(1) * t738 - t1011 * t742 + t1014 * t728, -t1011 * t859 + t1014 * t915, -t1011 * t830 + t1014 * t887, -t1011 * t857 + t1014 * t913, -t1011 * t858 + t1014 * t914, -t1011 * t856 + t1014 * t912, -t1011 * t930 + t1014 * t949, -qJ(1) * t800 - t1011 * t661 + t1014 * t684, -qJ(1) * t801 - t1011 * t662 + t1014 * t685, -qJ(1) * t860 - t1011 * t660 + t1014 * t730, -qJ(1) * t628 - t1011 * t573 + t1014 * t575, -t1011 * t676 + t1014 * t749, -t1011 * t611 + t1014 * t688, -t1011 * t681 + t1014 * t740, -t1011 * t675 + t1014 * t748, -t1011 * t682 + t1014 * t741, -t1011 * t737 + t1014 * t807, -qJ(1) * t602 - t1011 * t504 + t1014 * t515, -qJ(1) * t625 - t1011 * t510 + t1014 * t522, -qJ(1) * t576 - t1011 * t458 + t1014 * t474, -qJ(1) * t481 - t1011 * t436 + t1014 * t438, -t1011 * t541 + t1014 * t592, -t1011 * t495 + t1014 * t546, -t1011 * t551 + t1014 * t595, -t1011 * t540 + t1014 * t591, -t1011 * t552 + t1014 * t596, -t1011 * t585 + t1014 * t663, -qJ(1) * t500 - t1011 * t414 + t1014 * t429, -qJ(1) * t511 - t1011 * t419 + t1014 * t437, -qJ(1) * t462 - t1011 * t398 + t1014 * t409, -qJ(1) * t407 - t1011 * t384 + t1014 * t387, -t1011 * t473 + t1014 * t521, -t1011 * t448 + t1014 * t478, -t1011 * t455 + t1014 * t491, -t1011 * t472 + t1014 * t520, -t1011 * t456 + t1014 * t492, -t1011 * t488 + t1014 * t536, -qJ(1) * t431 - t1011 * t385 + t1014 * t391, -qJ(1) * t433 - t1011 * t386 + t1014 * t392, -qJ(1) * t423 - t1011 * t383 + t1014 * t390, -qJ(1) * t388 - t1011 * t379 + t1014 * t380; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t1168, -t1169, t932, qJ(1) * t932, 0, 0, t1087, 0, t920, -t1014 * t1116, qJ(1) * t922 + t1011 * t870 + t1014 * t815, qJ(1) * t921 + t1011 * t871 + t1014 * t814, -t1011 * t1061 + t1014 * t824, qJ(1) * t739 + t1011 * t728 + t1014 * t742, t1011 * t915 + t1014 * t859, t1011 * t887 + t1014 * t830, t1011 * t913 + t1014 * t857, t1011 * t914 + t1014 * t858, t1011 * t912 + t1014 * t856, t1011 * t949 + t1014 * t930, qJ(1) * t802 + t1011 * t684 + t1014 * t661, qJ(1) * t803 + t1011 * t685 + t1014 * t662, qJ(1) * t861 + t1011 * t730 + t1014 * t660, qJ(1) * t629 + t1011 * t575 + t1014 * t573, t1011 * t749 + t1014 * t676, t1011 * t688 + t1014 * t611, t1011 * t740 + t1014 * t681, t1011 * t748 + t1014 * t675, t1011 * t741 + t1014 * t682, t1011 * t807 + t1014 * t737, qJ(1) * t603 + t1011 * t515 + t1014 * t504, qJ(1) * t626 + t1011 * t522 + t1014 * t510, qJ(1) * t577 + t1011 * t474 + t1014 * t458, qJ(1) * t482 + t1011 * t438 + t1014 * t436, t1011 * t592 + t1014 * t541, t1011 * t546 + t1014 * t495, t1011 * t595 + t1014 * t551, t1011 * t591 + t1014 * t540, t1011 * t596 + t1014 * t552, t1011 * t663 + t1014 * t585, qJ(1) * t501 + t1011 * t429 + t1014 * t414, qJ(1) * t512 + t1011 * t437 + t1014 * t419, qJ(1) * t463 + t1011 * t409 + t1014 * t398, qJ(1) * t408 + t1011 * t387 + t1014 * t384, t1011 * t521 + t1014 * t473, t1011 * t478 + t1014 * t448, t1011 * t491 + t1014 * t455, t1011 * t520 + t1014 * t472, t1011 * t492 + t1014 * t456, t1011 * t536 + t1014 * t488, qJ(1) * t432 + t1011 * t391 + t1014 * t385, qJ(1) * t434 + t1011 * t392 + t1014 * t386, qJ(1) * t424 + t1011 * t390 + t1014 * t383, qJ(1) * t389 + t1011 * t380 + t1014 * t379; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t1101, t982, 0, 0, 0, 0, t959, 0, t960, t1015 * qJDD(2), pkin(1) * t962 + t1012 * t1090 - t1015 * t906, -pkin(1) * t961 + t1012 * t1091 - t1015 * t907, t832 * t1012, pkin(1) * t806 + t1012 * t1161, t1012 * t1031 + t1015 * t934, t1012 * t1060 + t1015 * t927, t1012 * t1049 + t1015 * t937, t1012 * t1032 + t1015 * t933, t1012 * t1048 + t1015 * t935, t1056 * t1012, pkin(1) * t845 + t1012 * t1034 + t1015 * t829, pkin(1) * t846 + t1012 * t1033 + t1015 * t828, pkin(1) * t924 + t1012 * t1092 + t1015 * t755, pkin(1) * t667 + t1012 * t1030 + t1015 * t747, t1012 * t1050 + t1015 * t764, t1012 * t1068 + t1015 * t708, t1012 * t1063 + t1015 * t793, t1012 * t1051 + t1015 * t763, t1012 * t1062 + t1015 * t794, t1015 * t820 + t1120 * t821 - t980, pkin(1) * t659 + t1012 * t1036 + t1015 * t610, pkin(1) * t673 + t1012 * t1035 + t1015 * t627, pkin(1) * t605 + t1012 * t1037 + t1015 * t529, pkin(1) * t517 + t1012 * t1038 + t1015 * t513, t1012 * t1052 + t1015 * t607, t1012 * t1076 + t1015 * t559, t1012 * t1070 + t1015 * t638, t1012 * t1053 + t1015 * t606, t1012 * t1069 + t1015 * t639, t1015 * t664 + t1120 * t665 - t980, pkin(1) * t528 + t1012 * t1040 + t1015 * t466, pkin(1) * t548 + t1012 * t1039 + t1015 * t477, pkin(1) * t494 + t1012 * t1041 + t1015 * t420, pkin(1) * t416 + t1012 * t1042 + t1015 * t405, t1012 * t1077 + t1015 * t533, t1012 * t1083 + t1015 * t496, t1012 * t1080 + t1015 * t505, t1012 * t1078 + t1015 * t532, t1012 * t1079 + t1015 * t506, t1012 * t1074 + t1015 * t565, pkin(1) * t452 + t1012 * t1044 + t1015 * t401, pkin(1) * t454 + t1012 * t1043 + t1015 * t402, pkin(1) * t442 + t1012 * t1045 + t1015 * t399, pkin(1) * t394 + t1012 * t1046 + t1015 * t381;];
tauB_reg  = t1;