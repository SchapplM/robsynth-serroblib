% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPPR8_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR8_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR8_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR8_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:39:34
% EndTime: 2019-12-31 19:39:47
% DurationCPUTime: 13.97s
% Computational Cost: add. (53906->648), mult. (124451->785), div. (0->0), fcn. (77484->8), ass. (0->413)
t1053 = cos(qJ(2));
t1129 = qJD(2) * t1053;
t1029 = qJD(1) * t1129;
t1050 = sin(qJ(2));
t1108 = t1050 * qJDD(1);
t1005 = 0.2e1 * t1029 + t1108;
t1051 = sin(qJ(1));
t1054 = cos(qJ(1));
t1055 = qJD(2) ^ 2;
t1044 = t1050 ^ 2;
t1056 = qJD(1) ^ 2;
t1115 = t1044 * t1056;
t1020 = t1055 + t1115;
t1113 = t1053 * t1056;
t1024 = t1050 * t1113;
t1016 = qJDD(2) - t1024;
t1114 = t1053 * t1016;
t966 = -t1020 * t1050 + t1114;
t1197 = pkin(5) * (-t1005 * t1051 + t1054 * t966);
t1196 = pkin(5) * (t1005 * t1054 + t1051 * t966);
t1195 = pkin(6) * t966;
t1165 = pkin(2) + pkin(3);
t1046 = sin(pkin(8));
t1047 = cos(pkin(8));
t1133 = qJD(1) * t1050;
t1014 = -qJD(2) * pkin(3) - qJ(4) * t1133;
t1033 = t1053 * qJDD(1);
t1100 = qJD(2) * t1133;
t1080 = -t1033 + t1100;
t1167 = t1053 ^ 2;
t1112 = t1167 * t1056;
t1128 = qJD(3) * qJD(2);
t1034 = 0.2e1 * t1128;
t1159 = pkin(2) * t1053;
t1082 = -qJ(3) * t1050 - t1159;
t1003 = t1082 * qJD(1);
t1132 = qJD(1) * t1053;
t1018 = g(1) * t1054 + g(2) * t1051;
t994 = -pkin(1) * t1056 + qJDD(1) * pkin(6) - t1018;
t969 = -t1050 * g(3) + t1053 * t994;
t1076 = -t1055 * pkin(2) + qJDD(2) * qJ(3) + t1003 * t1132 + t969;
t905 = t1034 + t1076;
t877 = -pkin(3) * t1112 + qJ(4) * t1080 + qJD(2) * t1014 + t905;
t1006 = t1029 + t1108;
t1158 = g(3) * t1053;
t1073 = -qJDD(2) * pkin(2) - t1055 * qJ(3) + qJDD(3) + t1158;
t1092 = qJD(1) * t1003 + t994;
t879 = -qJDD(2) * pkin(3) + (-t1006 + t1029) * qJ(4) + (-pkin(3) * t1113 + t1092) * t1050 + t1073;
t991 = -t1046 * t1132 + t1047 * t1133;
t812 = 0.2e1 * qJD(4) * t991 + t1046 * t877 - t1047 * t879;
t990 = (t1046 * t1050 + t1047 * t1053) * qJD(1);
t944 = t991 * t990;
t1178 = -qJDD(2) - t944;
t1187 = t1047 * t1178;
t1168 = t990 ^ 2;
t937 = -t1055 - t1168;
t881 = t1046 * t937 + t1187;
t1188 = t1046 * t1178;
t882 = t1047 * t937 - t1188;
t1194 = qJ(3) * t882 - t1165 * t881 + t812;
t1120 = t1016 * t1050;
t960 = t1020 * t1053 + t1120;
t1191 = pkin(1) * t960;
t1190 = pkin(6) * t960;
t1015 = qJDD(2) + t1024;
t1023 = -t1055 - t1112;
t915 = t1050 * t1092 + t1073;
t1189 = pkin(2) * t1015 + qJ(3) * t1023 - t915;
t1049 = sin(qJ(5));
t1039 = qJDD(2) - qJDD(5);
t1052 = cos(qJ(5));
t928 = t1049 * t991 + t1052 * t990;
t930 = -t1049 * t990 + t1052 * t991;
t880 = t930 * t928;
t1179 = -t880 - t1039;
t1186 = t1049 * t1179;
t1185 = t1052 * t1179;
t1022 = -t1055 + t1112;
t964 = -t1022 * t1053 + t1120;
t1184 = t1054 * t1033 + t1051 * t964;
t1013 = (t1044 - t1167) * t1056;
t1127 = t1005 * t1050;
t1007 = t1033 - 0.2e1 * t1100;
t954 = t1007 * t1053;
t949 = -t954 + t1127;
t1183 = t1013 * t1054 + t1051 * t949;
t1182 = -t1051 * t1033 + t1054 * t964;
t1181 = -t1013 * t1051 + t1054 * t949;
t1091 = t1006 * t1046 - t1047 * t1080;
t945 = t1047 * t1006 + t1046 * t1080;
t855 = -t928 * qJD(5) - t1049 * t1091 + t1052 * t945;
t1040 = qJD(2) - qJD(5);
t922 = t928 * t1040;
t1180 = t922 + t855;
t984 = qJD(2) * t990;
t910 = -t984 + t945;
t985 = qJD(2) * t991;
t907 = t1091 + t985;
t1152 = t1046 * t879 + t1047 * t877;
t1149 = qJD(4) * t990;
t980 = -0.2e1 * t1149;
t813 = t980 + t1152;
t769 = t1046 * t813 - t1047 * t812;
t956 = t1022 * t1050 + t1114;
t796 = pkin(4) * t1178 - t910 * pkin(7) - t812;
t974 = -qJD(2) * pkin(4) - pkin(7) * t991;
t799 = -t1168 * pkin(4) - pkin(7) * t1091 + qJD(2) * t974 + t813;
t760 = t1049 * t799 - t1052 * t796;
t1038 = t1040 ^ 2;
t926 = t928 ^ 2;
t871 = -t1038 - t926;
t816 = t1049 * t871 + t1185;
t1085 = -pkin(4) * t816 + t760;
t817 = t1052 * t871 - t1186;
t773 = t1046 * t817 + t1047 * t816;
t774 = -t1046 * t816 + t1047 * t817;
t1176 = qJ(3) * t774 - t1165 * t773 + t1085;
t761 = t1049 * t796 + t1052 * t799;
t874 = -t880 + t1039;
t1139 = t1049 * t874;
t927 = t930 ^ 2;
t914 = -t927 - t1038;
t832 = t1052 * t914 + t1139;
t1106 = -pkin(4) * t832 + t761;
t1136 = t1052 * t874;
t833 = -t1049 * t914 + t1136;
t787 = t1046 * t833 + t1047 * t832;
t788 = -t1046 * t832 + t1047 * t833;
t1175 = qJ(3) * t788 - t1165 * t787 + t1106;
t938 = qJDD(2) - t944;
t1143 = t1046 * t938;
t989 = t991 ^ 2;
t978 = -t989 - t1055;
t888 = t1047 * t978 + t1143;
t1141 = t1047 * t938;
t891 = -t1046 * t978 + t1141;
t1174 = qJ(3) * t891 - t1165 * t888 + t1152;
t1093 = t1049 * t945 + t1052 * t1091;
t825 = (qJD(5) + t1040) * t930 + t1093;
t829 = -t922 + t855;
t783 = -t1049 * t825 - t1052 * t829;
t785 = t1049 * t829 - t1052 * t825;
t745 = t1046 * t785 + t1047 * t783;
t747 = -t1046 * t783 + t1047 * t785;
t781 = pkin(4) * t783;
t1173 = qJ(3) * t747 - t1165 * t745 - t781;
t731 = t1049 * t761 - t1052 * t760;
t1142 = t1047 * t731;
t732 = t1049 * t760 + t1052 * t761;
t710 = t1046 * t732 + t1142;
t1146 = t1046 * t731;
t711 = t1047 * t732 - t1146;
t730 = pkin(4) * t731;
t1172 = qJ(3) * t711 - t1165 * t710 - t730;
t857 = -t1046 * t907 - t1047 * t910;
t859 = t1046 * t910 - t1047 * t907;
t1170 = qJ(3) * t859 - t1165 * t857;
t770 = t1046 * t812 + t1047 * t813;
t1169 = qJ(3) * t770 - t1165 * t769;
t1166 = 0.2e1 * qJD(3);
t998 = t1053 * t1015;
t957 = t1023 * t1050 + t998;
t1164 = pkin(1) * t957;
t1121 = t1015 * t1050;
t963 = t1023 * t1053 - t1121;
t1163 = pkin(5) * (t1007 * t1054 + t1051 * t963);
t1107 = t1044 + t1167;
t1009 = t1107 * qJDD(1);
t1012 = t1107 * t1056;
t1162 = pkin(5) * (t1009 * t1051 + t1012 * t1054);
t1161 = pkin(6) * t957;
t1154 = qJ(4) * t769;
t1153 = qJ(4) * t770;
t1151 = pkin(1) * t1007 + pkin(6) * t963;
t1147 = t1040 * t930;
t1017 = t1051 * g(1) - t1054 * g(2);
t993 = qJDD(1) * pkin(1) + t1056 * pkin(6) + t1017;
t1060 = -pkin(2) * t1100 + t993;
t863 = t1033 * pkin(2) + t1006 * qJ(3) + qJDD(4) - qJ(4) * t1112 - t1080 * pkin(3) + (qJ(3) * t1129 + (-pkin(2) * qJD(2) + t1014 + t1166) * t1050) * qJD(1) + t1060;
t1144 = t1046 * t863;
t862 = t1047 * t863;
t814 = pkin(4) * t1091 - t1168 * pkin(7) + t991 * t974 + t863;
t1140 = t1049 * t814;
t1138 = t1050 * t993;
t1137 = t1052 * t814;
t1135 = t1053 * t993;
t1134 = qJD(1) * qJD(2);
t1131 = qJD(2) * t1046;
t1130 = qJD(2) * t1047;
t1124 = t1007 * t1050;
t1117 = t1040 * t1049;
t1116 = t1040 * t1052;
t1110 = pkin(1) * t1012 + pkin(6) * t1009;
t1109 = qJDD(2) * t1054;
t1105 = t1051 * t880;
t1104 = t1051 * t944;
t1103 = t1054 * t880;
t1102 = t1054 * t944;
t1101 = t990 * t1131;
t1095 = -qJ(4) * t888 + t862;
t968 = t1050 * t994 + t1158;
t904 = t1050 * t968 + t1053 * t969;
t1090 = -t1017 * t1051 - t1054 * t1018;
t1089 = t1051 * t1024;
t1088 = t1054 * t1024;
t1084 = -pkin(2) * t915 + qJ(3) * t905;
t1011 = qJDD(1) * t1054 - t1051 * t1056;
t1083 = -pkin(5) * t1011 - g(3) * t1051;
t1081 = t1006 + t1029;
t1079 = -qJ(4) * t881 + t1144;
t1078 = -qJ(4) * t891 - t1144;
t1077 = -qJ(4) * t882 + t862;
t902 = t1050 * t969 - t1053 * t968;
t946 = t1005 * t1053 + t1124;
t1074 = t1017 * t1054 - t1018 * t1051;
t853 = -t926 - t927;
t717 = -pkin(4) * t853 + pkin(7) * t785 + t732;
t721 = -pkin(7) * t783 - t731;
t1072 = -qJ(4) * t745 - t1046 * t717 + t1047 * t721;
t824 = (qJD(5) - t1040) * t930 + t1093;
t763 = -pkin(4) * t824 + pkin(7) * t817 - t1137;
t772 = -pkin(7) * t816 + t1140;
t1071 = -qJ(4) * t773 - t1046 * t763 + t1047 * t772;
t765 = -pkin(4) * t1180 + pkin(7) * t833 + t1140;
t780 = -pkin(7) * t832 + t1137;
t1070 = -qJ(4) * t787 - t1046 * t765 + t1047 * t780;
t1069 = -qJ(4) * t857 - t769;
t1068 = -qJ(4) * t859 - t770;
t1067 = t1080 * pkin(2);
t1065 = -qJ(4) * t747 - t1046 * t721 - t1047 * t717;
t1064 = -qJ(4) * t774 - t1046 * t772 - t1047 * t763;
t1063 = -qJ(4) * t788 - t1046 * t780 - t1047 * t765;
t728 = -pkin(4) * t814 + pkin(7) * t732;
t1062 = -pkin(7) * t1142 - qJ(4) * t710 - t1046 * t728;
t1061 = pkin(7) * t1146 - qJ(4) * t711 - t1047 * t728;
t1059 = pkin(2) * t1020 + qJ(3) * t1016 + t1076;
t1058 = t1133 * t1166 + t1060;
t1057 = qJ(3) * t1081 + t1058;
t1032 = t1051 * qJDD(2);
t1021 = -t1055 + t1115;
t1010 = qJDD(1) * t1051 + t1054 * t1056;
t1001 = -pkin(2) * t1108 + qJ(3) * t1033;
t997 = t1107 * t1134;
t988 = -pkin(5) * t1010 + g(3) * t1054;
t977 = -t989 + t1055;
t976 = -t1055 + t1168;
t975 = t991 * t1130;
t973 = t1054 * t997 + t1032;
t972 = t1006 * t1053 - t1044 * t1134;
t971 = t1051 * t997 - t1109;
t970 = t1050 * t1080 - t1167 * t1134;
t965 = t1021 * t1050 + t998;
t959 = -t1021 * t1053 + t1121;
t955 = t1081 * t1050;
t950 = pkin(5) * (t1009 * t1054 - t1012 * t1051);
t942 = t989 - t1168;
t936 = t1054 * t972 - t1089;
t935 = t1054 * t970 + t1089;
t934 = t1051 * t972 + t1088;
t933 = t1051 * t970 - t1088;
t932 = t1051 * t1108 + t1054 * t965;
t931 = t1051 * t965 - t1054 * t1108;
t925 = (-t1046 * t991 + t1047 * t990) * qJD(2);
t924 = -t975 - t1101;
t920 = pkin(5) * (-t1007 * t1051 + t1054 * t963);
t919 = -t927 + t1038;
t918 = t926 - t1038;
t917 = -t1135 + t1190;
t916 = -t1138 - t1161;
t913 = t969 + t1191;
t912 = t968 - t1164;
t911 = t984 + t945;
t908 = t1091 - t985;
t906 = -t989 - t1168;
t901 = t1047 * t945 + t991 * t1131;
t900 = t1046 * t945 - t975;
t899 = t1046 * t1091 - t990 * t1130;
t898 = -t1047 * t1091 - t1101;
t897 = t1135 + t1151;
t896 = -pkin(1) * t1005 - t1138 - t1195;
t895 = qJ(3) * t1012 + t915;
t894 = pkin(2) * t1012 + t905;
t893 = -t1067 + t1057;
t892 = -t1046 * t977 + t1187;
t890 = t1047 * t976 + t1143;
t889 = t1047 * t977 + t1188;
t887 = t1046 * t976 - t1141;
t886 = pkin(1) * t993 + pkin(6) * t904;
t885 = (t1007 - t1080) * pkin(2) + t1057;
t884 = -t1067 + (t1005 + t1081) * qJ(3) + t1058;
t883 = t904 + t1110;
t878 = t927 - t926;
t869 = -t1050 * t924 + t1053 * t925;
t868 = t1050 * t925 + t1053 * t924;
t867 = (-t1049 * t930 + t1052 * t928) * t1040;
t866 = (t1049 * t928 + t1052 * t930) * t1040;
t865 = -t1164 - t1189;
t864 = -t1059 - 0.2e1 * t1128 - t1191;
t861 = t1050 * t915 + t1053 * t905;
t860 = t1050 * t905 - t1053 * t915;
t858 = -t1046 * t911 - t1047 * t908;
t856 = -t1046 * t908 + t1047 * t911;
t854 = -qJD(5) * t930 - t1093;
t852 = t1050 * t900 + t1053 * t901;
t851 = t1050 * t898 + t1053 * t899;
t850 = t1050 * t901 - t1053 * t900;
t849 = t1050 * t899 - t1053 * t898;
t848 = -pkin(2) * t1127 + t1053 * t884 - t1190;
t847 = qJ(3) * t954 - t1050 * t885 - t1161;
t846 = -t1050 * t894 + t1053 * t895;
t845 = t1050 * t889 + t1053 * t892;
t844 = t1050 * t888 + t1053 * t891;
t843 = t1050 * t887 + t1053 * t890;
t842 = t1050 * t892 - t1053 * t889;
t841 = t1050 * t891 - t1053 * t888;
t840 = t1050 * t890 - t1053 * t887;
t839 = t1195 + t1050 * t884 + (pkin(1) + t1159) * t1005;
t838 = qJ(3) * t1124 + t1053 * t885 + t1151;
t837 = t1052 * t918 + t1139;
t836 = -t1049 * t919 + t1185;
t835 = t1049 * t918 - t1136;
t834 = t1052 * t919 + t1186;
t830 = t1050 * t895 + t1053 * t894 + t1110;
t823 = t1052 * t855 + t930 * t1117;
t822 = t1049 * t855 - t930 * t1116;
t821 = -t1049 * t854 - t928 * t1116;
t820 = t1052 * t854 - t928 * t1117;
t819 = t1050 * t881 + t1053 * t882;
t818 = t1050 * t882 - t1053 * t881;
t809 = -t1046 * t866 + t1047 * t867;
t808 = t1046 * t867 + t1047 * t866;
t807 = -pkin(1) * t860 - t1084;
t806 = qJ(3) * t911 + t1095;
t805 = t1050 * t857 + t1053 * t859;
t804 = t1050 * t856 + t1053 * t858;
t803 = t1050 * t859 - t1053 * t857;
t802 = t1050 * t858 - t1053 * t856;
t801 = qJ(3) * t908 + t1079;
t800 = -pkin(6) * t860 + (-pkin(2) * t1050 + qJ(3) * t1053) * t893;
t797 = t1165 * t911 + t1078;
t795 = -t1046 * t835 + t1047 * t837;
t794 = -t1046 * t834 + t1047 * t836;
t793 = t1046 * t837 + t1047 * t835;
t792 = t1046 * t836 + t1047 * t834;
t789 = t1165 * t908 + t1077;
t786 = pkin(6) * t861 + (pkin(1) - t1082) * t893;
t784 = -t1049 * t1180 - t1052 * t824;
t782 = -t1049 * t824 + t1052 * t1180;
t779 = -t1046 * t822 + t1047 * t823;
t778 = -t1046 * t820 + t1047 * t821;
t777 = t1046 * t823 + t1047 * t822;
t776 = t1046 * t821 + t1047 * t820;
t768 = t1050 * t808 + t1053 * t809;
t767 = t1050 * t809 - t1053 * t808;
t766 = qJ(3) * t863 - t1154;
t764 = -pkin(1) * t841 + 0.2e1 * t1149 - t1174;
t762 = qJ(3) * t906 + t1069;
t758 = -pkin(1) * t803 - t1170;
t757 = t1050 * t793 + t1053 * t795;
t756 = t1050 * t792 + t1053 * t794;
t755 = t1050 * t795 - t1053 * t793;
t754 = t1050 * t794 - t1053 * t792;
t753 = t1050 * t787 + t1053 * t788;
t752 = t1050 * t788 - t1053 * t787;
t751 = t1165 * t863 - t1153;
t750 = -pkin(1) * t818 - t1194;
t749 = t1165 * t906 + t1068;
t748 = -pkin(6) * t841 - t1050 * t797 + t1053 * t806;
t746 = -t1046 * t782 + t1047 * t784;
t744 = t1046 * t784 + t1047 * t782;
t743 = t1050 * t777 + t1053 * t779;
t742 = t1050 * t776 + t1053 * t778;
t741 = t1050 * t779 - t1053 * t777;
t740 = t1050 * t778 - t1053 * t776;
t739 = pkin(1) * t911 + pkin(6) * t844 + t1050 * t806 + t1053 * t797;
t738 = -pkin(6) * t818 - t1050 * t789 + t1053 * t801;
t737 = t1050 * t773 + t1053 * t774;
t736 = t1050 * t774 - t1053 * t773;
t735 = pkin(1) * t908 + pkin(6) * t819 + t1050 * t801 + t1053 * t789;
t734 = t1050 * t769 + t1053 * t770;
t733 = t1050 * t770 - t1053 * t769;
t729 = qJ(3) * t1180 + t1070;
t727 = t1050 * t745 + t1053 * t747;
t726 = t1050 * t744 + t1053 * t746;
t725 = t1050 * t747 - t1053 * t745;
t724 = t1050 * t746 - t1053 * t744;
t723 = -pkin(6) * t803 - t1050 * t749 + t1053 * t762;
t722 = qJ(3) * t824 + t1071;
t719 = pkin(1) * t906 + pkin(6) * t805 + t1050 * t762 + t1053 * t749;
t718 = t1165 * t1180 + t1063;
t716 = t1165 * t824 + t1064;
t715 = -pkin(6) * t733 - t1050 * t751 + t1053 * t766;
t714 = -pkin(1) * t733 - t1169;
t713 = pkin(1) * t863 + pkin(6) * t734 + t1050 * t766 + t1053 * t751;
t712 = -pkin(1) * t752 - t1175;
t709 = -pkin(1) * t736 - t1176;
t708 = -pkin(6) * t752 - t1050 * t718 + t1053 * t729;
t707 = pkin(1) * t1180 + pkin(6) * t753 + t1050 * t729 + t1053 * t718;
t706 = -pkin(1) * t725 - t1173;
t705 = -pkin(6) * t736 - t1050 * t716 + t1053 * t722;
t704 = qJ(3) * t853 + t1072;
t703 = t1165 * t853 + t1065;
t702 = pkin(1) * t824 + pkin(6) * t737 + t1050 * t722 + t1053 * t716;
t701 = t1050 * t710 + t1053 * t711;
t700 = t1050 * t711 - t1053 * t710;
t699 = qJ(3) * t814 + t1062;
t698 = t1165 * t814 + t1061;
t697 = -pkin(6) * t725 - t1050 * t703 + t1053 * t704;
t696 = pkin(1) * t853 + pkin(6) * t727 + t1050 * t704 + t1053 * t703;
t695 = -pkin(1) * t700 - t1172;
t694 = -pkin(6) * t700 - t1050 * t698 + t1053 * t699;
t693 = pkin(1) * t814 + pkin(6) * t701 + t1050 * t699 + t1053 * t698;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1011, 0, -t1010, 0, t1083, -t988, -t1074, -pkin(5) * t1074, t936, -t1181, t932, t935, -t1182, t973, -t1051 * t912 + t1054 * t916 - t1163, -t1051 * t913 + t1054 * t917 + t1196, -t1054 * t902 - t1162, -pkin(5) * (t1051 * t904 + t1054 * t993) + (pkin(1) * t1051 - pkin(6) * t1054) * t902, t936, t932, t1181, t973, t1182, t935, -t1051 * t865 + t1054 * t847 - t1163, t1001 * t1051 + t1054 * t846 - t1162, -t1051 * t864 + t1054 * t848 - t1196, t1054 * t800 - t1051 * t807 - pkin(5) * (t1051 * t861 + t1054 * t893), t1054 * t852 - t1104, -t1051 * t942 + t1054 * t804, -t1051 * t910 + t1054 * t845, t1054 * t851 + t1104, t1051 * t907 + t1054 * t843, t1054 * t869 + t1032, t1054 * t738 - t1051 * t750 - pkin(5) * (t1051 * t819 + t1054 * t908), t1054 * t748 - t1051 * t764 - pkin(5) * (t1051 * t844 + t1054 * t911), t1054 * t723 - t1051 * t758 - pkin(5) * (t1051 * t805 + t1054 * t906), t1054 * t715 - t1051 * t714 - pkin(5) * (t1051 * t734 + t1054 * t863), t1054 * t743 - t1105, -t1051 * t878 + t1054 * t726, -t1051 * t829 + t1054 * t756, t1054 * t742 + t1105, t1051 * t825 + t1054 * t757, t1039 * t1051 + t1054 * t768, t1054 * t705 - t1051 * t709 - pkin(5) * (t1051 * t737 + t1054 * t824), t1054 * t708 - t1051 * t712 - pkin(5) * (t1051 * t753 + t1054 * t1180), t1054 * t697 - t1051 * t706 - pkin(5) * (t1051 * t727 + t1054 * t853), t1054 * t694 - t1051 * t695 - pkin(5) * (t1051 * t701 + t1054 * t814); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1010, 0, t1011, 0, t988, t1083, t1090, pkin(5) * t1090, t934, -t1183, t931, t933, -t1184, t971, t1051 * t916 + t1054 * t912 + t920, t1051 * t917 + t1054 * t913 - t1197, -t1051 * t902 + t950, pkin(5) * (-t1051 * t993 + t1054 * t904) + (-pkin(1) * t1054 - pkin(6) * t1051) * t902, t934, t931, t1183, t971, t1184, t933, t1051 * t847 + t1054 * t865 + t920, -t1001 * t1054 + t1051 * t846 + t950, t1051 * t848 + t1054 * t864 + t1197, t1051 * t800 + t1054 * t807 + pkin(5) * (-t1051 * t893 + t1054 * t861), t1051 * t852 + t1102, t1051 * t804 + t1054 * t942, t1051 * t845 + t1054 * t910, t1051 * t851 - t1102, t1051 * t843 - t1054 * t907, t1051 * t869 - t1109, t1051 * t738 + t1054 * t750 + pkin(5) * (-t1051 * t908 + t1054 * t819), t1051 * t748 + t1054 * t764 + pkin(5) * (-t1051 * t911 + t1054 * t844), t1051 * t723 + t1054 * t758 + pkin(5) * (-t1051 * t906 + t1054 * t805), t1051 * t715 + t1054 * t714 + pkin(5) * (-t1051 * t863 + t1054 * t734), t1051 * t743 + t1103, t1051 * t726 + t1054 * t878, t1051 * t756 + t1054 * t829, t1051 * t742 - t1103, t1051 * t757 - t1054 * t825, -t1039 * t1054 + t1051 * t768, t1051 * t705 + t1054 * t709 + pkin(5) * (-t1051 * t824 + t1054 * t737), t1051 * t708 + t1054 * t712 + pkin(5) * (-t1051 * t1180 + t1054 * t753), t1051 * t697 + t1054 * t706 + pkin(5) * (-t1051 * t853 + t1054 * t727), t1051 * t694 + t1054 * t695 + pkin(5) * (-t1051 * t814 + t1054 * t701); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1017, t1018, 0, 0, t955, t946, t959, t954, t956, 0, t897, t896, t883, t886, t955, t959, -t946, 0, -t956, t954, t838, t830, t839, t786, t850, t802, t842, t849, t840, t868, t735, t739, t719, t713, t741, t724, t754, t740, t755, t767, t702, t707, t696, t693; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1056, 0, 0, -g(3), -t1017, 0, t972, -t949, t965, t970, -t964, t997, t916, t917, -t902, -pkin(6) * t902, t972, t965, t949, t997, t964, t970, t847, t846, t848, t800, t852, t804, t845, t851, t843, t869, t738, t748, t723, t715, t743, t726, t756, t742, t757, t768, t705, t708, t697, t694; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1056, 0, qJDD(1), 0, g(3), 0, -t1018, 0, t1024, -t1013, -t1108, -t1024, -t1033, -qJDD(2), t912, t913, 0, -pkin(1) * t902, t1024, -t1108, t1013, -qJDD(2), t1033, -t1024, t865, -t1001, t864, t807, t944, t942, t910, -t944, -t907, -qJDD(2), t750, t764, t758, t714, t880, t878, t829, -t880, -t825, -t1039, t709, t712, t706, t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1017, t1018, 0, 0, t955, t946, t959, t954, t956, 0, t897, t896, t883, t886, t955, t959, -t946, 0, -t956, t954, t838, t830, t839, t786, t850, t802, t842, t849, t840, t868, t735, t739, t719, t713, t741, t724, t754, t740, t755, t767, t702, t707, t696, t693; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1006, t1007, t1015, -t1029, t1022, t1029, 0, -t993, t968, 0, t1006, t1015, -t1007, t1029, -t1022, -t1029, qJ(3) * t1007, t895, t884, qJ(3) * t893, t901, t858, t892, t899, t890, t925, t801, t806, t762, t766, t779, t746, t794, t778, t795, t809, t722, t729, t704, t699; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1100, t1005, -t1021, -t1080, t1016, -t1100, t993, 0, t969, 0, t1100, -t1021, -t1005, -t1100, -t1016, -t1080, t885, t894, pkin(2) * t1005, pkin(2) * t893, -t900, -t856, -t889, -t898, -t887, t924, t789, t797, t749, t751, -t777, -t744, -t792, -t776, -t793, -t808, t716, t718, t703, t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1024, t1013, t1108, t1024, t1033, qJDD(2), -t968, -t969, 0, 0, -t1024, t1108, -t1013, qJDD(2), -t1033, t1024, t1189, t1001, t1034 + t1059, t1084, -t944, -t942, -t910, t944, t907, qJDD(2), t1194, t1174 + t980, t1170, t1169, -t880, -t878, -t829, t880, t825, t1039, t1176, t1175, t1173, t1172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1006, t1015, -t1007, t1029, -t1022, -t1029, 0, t915, t893, 0, t901, t858, t892, t899, t890, t925, t1079, t1095, t1069, -t1154, t779, t746, t794, t778, t795, t809, t1071, t1070, t1072, t1062; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1024, t1108, -t1013, qJDD(2), -t1033, t1024, -t915, 0, t905, 0, -t944, -t942, -t910, t944, t907, qJDD(2), -pkin(3) * t881 + t812, -pkin(3) * t888 + t813, -pkin(3) * t857, -pkin(3) * t769, -t880, -t878, -t829, t880, t825, t1039, -pkin(3) * t773 + t1085, -pkin(3) * t787 + t1106, -pkin(3) * t745 - t781, -pkin(3) * t710 - t730; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1100, t1021, t1005, t1100, t1016, t1080, -t893, -t905, 0, 0, t900, t856, t889, t898, t887, -t924, -pkin(3) * t908 - t1077, -pkin(3) * t911 - t1078, -pkin(3) * t906 - t1068, -pkin(3) * t863 + t1153, t777, t744, t792, t776, t793, t808, -pkin(3) * t824 - t1064, -pkin(3) * t1180 - t1063, -pkin(3) * t853 - t1065, -pkin(3) * t814 - t1061; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t945, -t908, t1178, -t984, t976, t984, 0, t863, t812, 0, t823, t784, t836, t821, t837, t867, t772, t780, t721, -pkin(7) * t731; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t985, t911, t977, -t1091, -t938, t985, -t863, 0, t813, 0, t822, t782, t834, t820, t835, t866, t763, t765, t717, t728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t944, t942, t910, -t944, -t907, -qJDD(2), -t812, -t813, 0, 0, t880, t878, t829, -t880, -t825, -t1039, -t1085, -t1106, t781, t730; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t855, -t824, t1179, -t922, t918, t922, 0, t814, t760, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1147, t1180, t919, t854, -t874, t1147, -t814, 0, t761, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t880, t878, t829, -t880, -t825, -t1039, -t760, -t761, 0, 0;];
m_new_reg = t1;
