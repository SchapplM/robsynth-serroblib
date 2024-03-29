% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPRPR4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:23:28
% EndTime: 2022-01-23 09:23:45
% DurationCPUTime: 16.97s
% Computational Cost: add. (92771->642), mult. (203831->922), div. (0->0), fcn. (136613->10), ass. (0->450)
t1143 = sin(pkin(9));
t1145 = cos(pkin(9));
t1149 = sin(qJ(3));
t1152 = cos(qJ(3));
t1093 = (t1143 * t1152 + t1145 * t1149) * qJD(1);
t1150 = sin(qJ(1));
t1153 = cos(qJ(1));
t1121 = g(1) * t1153 + g(2) * t1150;
t1237 = qJD(1) ^ 2;
t1104 = -pkin(1) * t1237 - t1121;
t1144 = sin(pkin(8));
t1146 = cos(pkin(8));
t1120 = g(1) * t1150 - g(2) * t1153;
t1160 = qJDD(1) * pkin(1) + t1120;
t1052 = t1146 * t1104 + t1144 * t1160;
t1033 = -pkin(2) * t1237 + qJDD(1) * pkin(6) + t1052;
t1141 = g(3) - qJDD(2);
t1007 = t1033 * t1149 + t1141 * t1152;
t1214 = qJD(1) * t1152;
t1131 = qJD(3) * t1214;
t1193 = t1149 * qJDD(1);
t1107 = t1131 + t1193;
t1127 = t1152 * t1237 * t1149;
t1118 = qJDD(3) + t1127;
t967 = (-t1107 + t1131) * qJ(4) + t1118 * pkin(3) - t1007;
t1009 = t1033 * t1152 - t1141 * t1149;
t1215 = qJD(1) * t1149;
t1117 = qJD(3) * pkin(3) - qJ(4) * t1215;
t1235 = t1152 ^ 2;
t1135 = t1235 * t1237;
t1184 = qJD(3) * t1215;
t1192 = t1152 * qJDD(1);
t1157 = -t1184 + t1192;
t970 = -pkin(3) * t1135 + qJ(4) * t1157 - qJD(3) * t1117 + t1009;
t1163 = -0.2e1 * qJD(4) * t1093 - t1143 * t970 + t1145 * t967;
t1154 = qJD(3) ^ 2;
t1092 = t1143 * t1215 - t1145 * t1214;
t1236 = t1092 ^ 2;
t1043 = -t1154 - t1236;
t1050 = t1093 * t1092;
t1239 = qJDD(3) - t1050;
t1251 = t1145 * t1239;
t972 = t1043 * t1143 + t1251;
t1255 = pkin(3) * t972 + t1163;
t1051 = t1104 * t1144 - t1146 * t1160;
t1173 = t1051 * t1144 + t1052 * t1146;
t981 = t1051 * t1146 - t1052 * t1144;
t1222 = t1150 * t981;
t1254 = t1153 * t1173 + t1222;
t1218 = t1153 * t981;
t1253 = -t1150 * t1173 + t1218;
t1252 = t1143 * t1239;
t1148 = sin(qJ(5));
t1138 = qJDD(3) + qJDD(5);
t1151 = cos(qJ(5));
t1029 = t1092 * t1151 + t1093 * t1148;
t1031 = -t1092 * t1148 + t1093 * t1151;
t971 = t1031 * t1029;
t1244 = -t971 + t1138;
t1250 = t1148 * t1244;
t1249 = t1151 * t1244;
t1110 = qJDD(1) * t1144 + t1146 * t1237;
t1080 = qJ(2) * t1110 - t1141 * t1146;
t1111 = qJDD(1) * t1146 - t1144 * t1237;
t1162 = -qJ(2) * t1111 - t1141 * t1144;
t1240 = t1110 * t1153 + t1111 * t1150;
t1248 = pkin(5) * t1240 + t1153 * t1080 - t1150 * t1162;
t1056 = -t1110 * t1150 + t1111 * t1153;
t1247 = -pkin(5) * t1056 + t1150 * t1080 + t1153 * t1162;
t953 = t1007 * t1149 + t1009 * t1152;
t1139 = qJD(3) + qJD(5);
t1020 = t1139 * t1029;
t1053 = t1107 * t1145 + t1143 * t1157;
t1171 = t1107 * t1143 - t1145 * t1157;
t945 = -qJD(5) * t1029 + t1053 * t1151 - t1148 * t1171;
t1245 = -t1020 + t945;
t1088 = qJD(3) * t1092;
t1005 = t1088 + t1053;
t1172 = t1053 * t1148 + t1151 * t1171;
t911 = (qJD(5) - t1139) * t1031 + t1172;
t1026 = t1029 ^ 2;
t1027 = t1031 ^ 2;
t1091 = t1093 ^ 2;
t1137 = t1139 ^ 2;
t1211 = qJD(4) * t1092;
t1085 = -0.2e1 * t1211;
t1230 = t1143 * t967 + t1145 * t970;
t888 = t1085 + t1230;
t832 = t1143 * t888 + t1145 * t1163;
t1234 = pkin(3) * t832;
t1212 = qJD(3) * t1093;
t1003 = -t1171 + t1212;
t947 = t1003 * t1143 - t1005 * t1145;
t1233 = pkin(3) * t947;
t859 = pkin(4) * t1239 - pkin(7) * t1005 + t1163;
t1076 = qJD(3) * pkin(4) - pkin(7) * t1093;
t862 = -pkin(4) * t1236 - pkin(7) * t1171 - qJD(3) * t1076 + t888;
t814 = t1148 * t862 - t1151 * t859;
t815 = t1148 * t859 + t1151 * t862;
t770 = t1148 * t815 - t1151 * t814;
t1229 = t1143 * t770;
t1032 = -qJDD(1) * pkin(2) - pkin(6) * t1237 + t1051;
t978 = -pkin(3) * t1157 - qJ(4) * t1135 + t1117 * t1215 + qJDD(4) + t1032;
t1228 = t1143 * t978;
t1227 = t1145 * t770;
t1226 = t1145 * t978;
t896 = pkin(4) * t1171 - pkin(7) * t1236 + t1076 * t1093 + t978;
t1225 = t1148 * t896;
t964 = t971 + t1138;
t1224 = t1148 * t964;
t1223 = t1149 * t832;
t1221 = t1151 * t896;
t1220 = t1151 * t964;
t1219 = t1152 * t832;
t1217 = -pkin(2) * t1032 + pkin(6) * t953;
t1216 = qJD(1) * qJD(3);
t1140 = t1149 ^ 2;
t1213 = t1237 * t1140;
t1045 = qJDD(3) + t1050;
t1209 = t1045 * t1143;
t1208 = t1045 * t1145;
t1207 = t1092 * t1143;
t1206 = t1092 * t1145;
t1205 = t1093 * t1143;
t1204 = t1093 * t1145;
t1108 = -0.2e1 * t1184 + t1192;
t1062 = t1108 * t1152;
t1201 = t1118 * t1149;
t1119 = qJDD(3) - t1127;
t1200 = t1119 * t1149;
t1199 = t1119 * t1152;
t1198 = t1139 * t1031;
t1197 = t1139 * t1148;
t1196 = t1139 * t1151;
t1023 = t1149 * t1032;
t1024 = t1152 * t1032;
t1194 = qJDD(3) * t1146;
t1191 = t1140 + t1235;
t771 = t1148 * t814 + t1151 * t815;
t738 = t1143 * t771 + t1227;
t769 = pkin(4) * t770;
t1190 = pkin(3) * t738 + t769;
t915 = t1020 + t945;
t845 = -t1148 * t911 - t1151 * t915;
t847 = t1148 * t915 - t1151 * t911;
t802 = t1143 * t847 + t1145 * t845;
t843 = pkin(4) * t845;
t1189 = pkin(3) * t802 + t843;
t1188 = t1144 * t971;
t1187 = t1146 * t971;
t1186 = t1144 * t1050;
t1185 = t1146 * t1050;
t1124 = -t1154 - t1213;
t1071 = -t1124 * t1149 - t1199;
t1106 = 0.2e1 * t1131 + t1193;
t1183 = -pkin(2) * t1106 + pkin(6) * t1071 + t1023;
t1126 = -t1135 - t1154;
t1069 = t1126 * t1152 - t1201;
t1182 = pkin(2) * t1108 + pkin(6) * t1069 - t1024;
t833 = -t1143 * t1163 + t1145 * t888;
t739 = t1145 * t771 - t1229;
t761 = -pkin(4) * t896 + pkin(7) * t771;
t718 = -pkin(3) * t896 - pkin(7) * t1229 + qJ(4) * t739 + t1145 * t761;
t722 = -pkin(7) * t1227 - qJ(4) * t738 - t1143 * t761;
t726 = -t1149 * t738 + t1152 * t739;
t1180 = -pkin(2) * t896 + pkin(6) * t726 + t1149 * t722 + t1152 * t718;
t943 = -t1026 - t1027;
t746 = -pkin(4) * t943 + pkin(7) * t847 + t771;
t749 = -pkin(7) * t845 - t770;
t804 = -t1143 * t845 + t1145 * t847;
t730 = -pkin(3) * t943 + qJ(4) * t804 + t1143 * t749 + t1145 * t746;
t733 = -qJ(4) * t802 - t1143 * t746 + t1145 * t749;
t759 = -t1149 * t802 + t1152 * t804;
t1179 = -pkin(2) * t943 + pkin(6) * t759 + t1149 * t733 + t1152 * t730;
t962 = -t1137 - t1026;
t891 = t1151 * t962 - t1250;
t910 = (qJD(5) + t1139) * t1031 + t1172;
t818 = -pkin(4) * t910 + pkin(7) * t891 - t1221;
t890 = t1148 * t962 + t1249;
t836 = -t1143 * t890 + t1145 * t891;
t837 = -pkin(7) * t890 + t1225;
t754 = -pkin(3) * t910 + qJ(4) * t836 + t1143 * t837 + t1145 * t818;
t835 = t1143 * t891 + t1145 * t890;
t765 = -qJ(4) * t835 - t1143 * t818 + t1145 * t837;
t795 = -t1149 * t835 + t1152 * t836;
t1178 = -pkin(2) * t910 + pkin(6) * t795 + t1149 * t765 + t1152 * t754;
t1010 = -t1027 - t1137;
t924 = -t1010 * t1148 - t1220;
t823 = -pkin(4) * t1245 + pkin(7) * t924 + t1225;
t923 = t1010 * t1151 - t1224;
t848 = -pkin(7) * t923 + t1221;
t851 = -t1143 * t923 + t1145 * t924;
t763 = -pkin(3) * t1245 + qJ(4) * t851 + t1143 * t848 + t1145 * t823;
t850 = t1143 * t924 + t1145 * t923;
t773 = -qJ(4) * t850 - t1143 * t823 + t1145 * t848;
t807 = -t1149 * t850 + t1152 * t851;
t1177 = -pkin(2) * t1245 + pkin(6) * t807 + t1149 * t773 + t1152 * t763;
t1000 = -t1091 - t1236;
t949 = t1003 * t1145 + t1005 * t1143;
t816 = -pkin(3) * t1000 + qJ(4) * t949 + t833;
t820 = -qJ(4) * t947 - t832;
t871 = -t1149 * t947 + t1152 * t949;
t1176 = -pkin(2) * t1000 + pkin(6) * t871 + t1149 * t820 + t1152 * t816;
t1001 = t1171 + t1212;
t973 = t1043 * t1145 - t1252;
t873 = -pkin(3) * t1001 + qJ(4) * t973 - t1226;
t894 = -t1149 * t972 + t1152 * t973;
t907 = -qJ(4) * t972 + t1228;
t1175 = -pkin(2) * t1001 + pkin(6) * t894 + t1149 * t907 + t1152 * t873;
t1004 = -t1088 + t1053;
t1083 = -t1091 - t1154;
t990 = -t1083 * t1143 - t1208;
t879 = -pkin(3) * t1004 + qJ(4) * t990 + t1228;
t987 = t1083 * t1145 - t1209;
t919 = -qJ(4) * t987 + t1226;
t935 = -t1149 * t987 + t1152 * t990;
t1174 = -pkin(2) * t1004 + pkin(6) * t935 + t1149 * t919 + t1152 * t879;
t1170 = -t1120 * t1150 - t1121 * t1153;
t1112 = t1191 * qJDD(1);
t1115 = t1135 + t1213;
t1169 = pkin(2) * t1115 + pkin(6) * t1112 + t953;
t1168 = pkin(3) * t987 - t1230;
t1167 = t1144 * t1127;
t1166 = t1146 * t1127;
t1165 = pkin(4) * t890 - t814;
t1114 = qJDD(1) * t1153 - t1150 * t1237;
t1164 = -pkin(5) * t1114 - g(3) * t1150;
t952 = t1007 * t1152 - t1009 * t1149;
t1161 = t1120 * t1153 - t1121 * t1150;
t1159 = pkin(4) * t923 - t815;
t782 = t1152 * t833 - t1223;
t825 = -pkin(3) * t978 + qJ(4) * t833;
t1158 = -pkin(2) * t978 + pkin(6) * t782 - qJ(4) * t1223 + t1152 * t825;
t1156 = pkin(3) * t835 + t1165;
t1155 = pkin(3) * t850 + t1159;
t1133 = t1144 * qJDD(3);
t1125 = t1135 - t1154;
t1123 = t1154 - t1213;
t1116 = -t1135 + t1213;
t1113 = qJDD(1) * t1150 + t1153 * t1237;
t1102 = t1152 * t1118;
t1101 = t1191 * t1216;
t1090 = -pkin(5) * t1113 + g(3) * t1153;
t1082 = -t1091 + t1154;
t1081 = -t1154 + t1236;
t1075 = t1107 * t1152 - t1140 * t1216;
t1074 = -t1149 * t1157 - t1216 * t1235;
t1073 = t1101 * t1146 + t1133;
t1072 = t1101 * t1144 - t1194;
t1070 = -t1123 * t1149 + t1102;
t1068 = t1125 * t1152 - t1200;
t1067 = t1124 * t1152 - t1200;
t1066 = t1123 * t1152 + t1201;
t1065 = t1126 * t1149 + t1102;
t1064 = t1125 * t1149 + t1199;
t1063 = (t1107 + t1131) * t1149;
t1059 = t1112 * t1146 - t1115 * t1144;
t1058 = t1112 * t1144 + t1115 * t1146;
t1055 = -t1106 * t1149 + t1062;
t1054 = t1106 * t1152 + t1108 * t1149;
t1048 = t1091 - t1236;
t1041 = t1075 * t1146 - t1167;
t1040 = t1074 * t1146 + t1167;
t1039 = t1075 * t1144 + t1166;
t1038 = t1074 * t1144 - t1166;
t1037 = t1070 * t1146 + t1144 * t1193;
t1036 = t1068 * t1146 + t1144 * t1192;
t1035 = t1070 * t1144 - t1146 * t1193;
t1034 = t1068 * t1144 - t1146 * t1192;
t1022 = (t1205 - t1206) * qJD(3);
t1021 = (-t1204 - t1207) * qJD(3);
t1018 = t1071 * t1146 + t1106 * t1144;
t1017 = t1069 * t1146 - t1108 * t1144;
t1016 = t1071 * t1144 - t1106 * t1146;
t1015 = t1069 * t1144 + t1108 * t1146;
t1014 = -pkin(1) * t1110 - t1052;
t1013 = pkin(1) * t1111 - t1051;
t1012 = -t1027 + t1137;
t1011 = t1026 - t1137;
t1008 = t1055 * t1146 + t1116 * t1144;
t1006 = t1055 * t1144 - t1116 * t1146;
t994 = -qJD(3) * t1205 + t1053 * t1145;
t993 = qJD(3) * t1204 + t1053 * t1143;
t992 = qJD(3) * t1206 + t1143 * t1171;
t991 = qJD(3) * t1207 - t1145 * t1171;
t989 = -t1082 * t1143 + t1251;
t988 = t1081 * t1145 - t1209;
t986 = t1082 * t1145 + t1252;
t985 = t1081 * t1143 + t1208;
t984 = -pkin(6) * t1067 + t1024;
t983 = -pkin(6) * t1065 + t1023;
t977 = pkin(1) * t981;
t976 = -pkin(2) * t1067 + t1009;
t975 = -pkin(2) * t1065 + t1007;
t969 = pkin(1) * t1141 + qJ(2) * t1173;
t968 = t1027 - t1026;
t959 = -t1021 * t1149 + t1022 * t1152;
t958 = t1021 * t1152 + t1022 * t1149;
t957 = (-t1029 * t1151 + t1031 * t1148) * t1139;
t956 = (-t1029 * t1148 - t1031 * t1151) * t1139;
t955 = t1146 * t959 + t1133;
t954 = t1144 * t959 - t1194;
t948 = -t1001 * t1145 - t1004 * t1143;
t946 = -t1001 * t1143 + t1004 * t1145;
t944 = -qJD(5) * t1031 - t1172;
t942 = -t1149 * t993 + t1152 * t994;
t941 = -t1149 * t991 + t1152 * t992;
t940 = t1149 * t994 + t1152 * t993;
t939 = t1149 * t992 + t1152 * t991;
t937 = pkin(1) * t1015 + t1182;
t936 = pkin(1) * t1016 + t1183;
t934 = -t1149 * t986 + t1152 * t989;
t933 = -t1149 * t985 + t1152 * t988;
t932 = t1149 * t990 + t1152 * t987;
t931 = t1149 * t989 + t1152 * t986;
t930 = t1149 * t988 + t1152 * t985;
t928 = t1011 * t1151 - t1224;
t927 = -t1012 * t1148 + t1249;
t926 = t1011 * t1148 + t1220;
t925 = t1012 * t1151 + t1250;
t922 = -qJ(2) * t1058 + t1146 * t952;
t921 = qJ(2) * t1059 + t1144 * t952;
t917 = t1032 * t1144 + t1146 * t953;
t916 = -t1032 * t1146 + t1144 * t953;
t906 = pkin(1) * t1058 + t1169;
t905 = t1146 * t942 + t1186;
t904 = t1146 * t941 - t1186;
t903 = t1144 * t942 - t1185;
t902 = t1144 * t941 + t1185;
t900 = -t1031 * t1197 + t1151 * t945;
t899 = t1031 * t1196 + t1148 * t945;
t898 = t1029 * t1196 - t1148 * t944;
t897 = t1029 * t1197 + t1151 * t944;
t893 = t1149 * t973 + t1152 * t972;
t886 = t1004 * t1144 + t1146 * t935;
t885 = t1005 * t1144 + t1146 * t934;
t884 = t1003 * t1144 + t1146 * t933;
t883 = -t1004 * t1146 + t1144 * t935;
t882 = -t1005 * t1146 + t1144 * t934;
t881 = -t1003 * t1146 + t1144 * t933;
t878 = -t1143 * t956 + t1145 * t957;
t877 = t1143 * t957 + t1145 * t956;
t875 = -qJ(2) * t1016 - t1144 * t976 + t1146 * t984;
t874 = -qJ(2) * t1015 - t1144 * t975 + t1146 * t983;
t870 = -t1149 * t946 + t1152 * t948;
t869 = t1149 * t949 + t1152 * t947;
t868 = t1149 * t948 + t1152 * t946;
t867 = t1001 * t1144 + t1146 * t894;
t866 = -t1001 * t1146 + t1144 * t894;
t864 = -pkin(1) * t1067 + qJ(2) * t1018 + t1144 * t984 + t1146 * t976;
t863 = -pkin(1) * t1065 + qJ(2) * t1017 + t1144 * t983 + t1146 * t975;
t861 = t1048 * t1144 + t1146 * t870;
t860 = -t1048 * t1146 + t1144 * t870;
t857 = -t1143 * t926 + t1145 * t928;
t856 = -t1143 * t925 + t1145 * t927;
t855 = t1143 * t928 + t1145 * t926;
t854 = t1143 * t927 + t1145 * t925;
t853 = t1000 * t1144 + t1146 * t871;
t852 = -t1000 * t1146 + t1144 * t871;
t849 = pkin(1) * t916 + t1217;
t846 = -t1148 * t1245 - t1151 * t910;
t844 = -t1148 * t910 + t1151 * t1245;
t842 = -t1143 * t899 + t1145 * t900;
t841 = -t1143 * t897 + t1145 * t898;
t840 = t1143 * t900 + t1145 * t899;
t839 = t1143 * t898 + t1145 * t897;
t838 = -pkin(2) * t869 - t1233;
t834 = -qJ(2) * t916 - (pkin(2) * t1144 - pkin(6) * t1146) * t952;
t831 = -pkin(2) * t932 + t1085 - t1168;
t830 = -t1149 * t877 + t1152 * t878;
t829 = t1149 * t878 + t1152 * t877;
t828 = t1138 * t1144 + t1146 * t830;
t827 = -t1138 * t1146 + t1144 * t830;
t826 = -pkin(2) * t893 - t1255;
t822 = -pkin(6) * t932 - t1149 * t879 + t1152 * t919;
t821 = qJ(2) * t917 - (-pkin(2) * t1146 - pkin(6) * t1144 - pkin(1)) * t952;
t817 = -pkin(6) * t893 - t1149 * t873 + t1152 * t907;
t811 = -t1149 * t855 + t1152 * t857;
t810 = -t1149 * t854 + t1152 * t856;
t809 = t1149 * t857 + t1152 * t855;
t808 = t1149 * t856 + t1152 * t854;
t806 = t1149 * t851 + t1152 * t850;
t803 = -t1143 * t844 + t1145 * t846;
t801 = t1143 * t846 + t1145 * t844;
t800 = -t1149 * t840 + t1152 * t842;
t799 = -t1149 * t839 + t1152 * t841;
t798 = t1149 * t842 + t1152 * t840;
t797 = t1149 * t841 + t1152 * t839;
t796 = pkin(1) * t883 + t1174;
t794 = t1149 * t836 + t1152 * t835;
t792 = -t1144 * t911 + t1146 * t811;
t791 = t1144 * t915 + t1146 * t810;
t790 = t1144 * t811 + t1146 * t911;
t789 = t1144 * t810 - t1146 * t915;
t788 = t1144 * t1245 + t1146 * t807;
t787 = t1144 * t807 - t1146 * t1245;
t786 = t1146 * t800 + t1188;
t785 = t1146 * t799 - t1188;
t784 = t1144 * t800 - t1187;
t783 = t1144 * t799 + t1187;
t781 = t1149 * t833 + t1219;
t779 = pkin(1) * t866 + t1175;
t778 = t1144 * t978 + t1146 * t782;
t777 = t1144 * t782 - t1146 * t978;
t776 = t1144 * t910 + t1146 * t795;
t775 = t1144 * t795 - t1146 * t910;
t774 = -qJ(2) * t883 - t1144 * t831 + t1146 * t822;
t768 = -pkin(2) * t781 - t1234;
t767 = -pkin(1) * t932 + qJ(2) * t886 + t1144 * t822 + t1146 * t831;
t766 = -qJ(2) * t866 - t1144 * t826 + t1146 * t817;
t760 = -pkin(6) * t869 - t1149 * t816 + t1152 * t820;
t758 = -t1149 * t801 + t1152 * t803;
t757 = t1149 * t804 + t1152 * t802;
t756 = t1149 * t803 + t1152 * t801;
t752 = -pkin(1) * t893 + qJ(2) * t867 + t1144 * t817 + t1146 * t826;
t751 = t1144 * t968 + t1146 * t758;
t750 = t1144 * t758 - t1146 * t968;
t748 = t1144 * t943 + t1146 * t759;
t747 = t1144 * t759 - t1146 * t943;
t745 = -pkin(2) * t806 - t1155;
t744 = pkin(1) * t852 + t1176;
t743 = -pkin(6) * t781 - qJ(4) * t1219 - t1149 * t825;
t742 = -pkin(2) * t794 - t1156;
t741 = -qJ(2) * t852 - t1144 * t838 + t1146 * t760;
t740 = -pkin(1) * t869 + qJ(2) * t853 + t1144 * t760 + t1146 * t838;
t737 = -pkin(2) * t757 - t1189;
t736 = pkin(1) * t777 + t1158;
t735 = -pkin(6) * t806 - t1149 * t763 + t1152 * t773;
t734 = -pkin(6) * t794 - t1149 * t754 + t1152 * t765;
t731 = -qJ(2) * t777 - t1144 * t768 + t1146 * t743;
t728 = pkin(1) * t787 + t1177;
t727 = pkin(1) * t775 + t1178;
t725 = t1149 * t739 + t1152 * t738;
t723 = -pkin(1) * t781 + qJ(2) * t778 + t1144 * t743 + t1146 * t768;
t720 = t1144 * t896 + t1146 * t726;
t719 = t1144 * t726 - t1146 * t896;
t716 = -qJ(2) * t787 - t1144 * t745 + t1146 * t735;
t715 = -pkin(1) * t806 + qJ(2) * t788 + t1144 * t735 + t1146 * t745;
t714 = -qJ(2) * t775 - t1144 * t742 + t1146 * t734;
t713 = -pkin(1) * t794 + qJ(2) * t776 + t1144 * t734 + t1146 * t742;
t712 = -pkin(2) * t725 - t1190;
t711 = -pkin(6) * t757 - t1149 * t730 + t1152 * t733;
t710 = pkin(1) * t747 + t1179;
t709 = -qJ(2) * t747 - t1144 * t737 + t1146 * t711;
t708 = -pkin(1) * t757 + qJ(2) * t748 + t1144 * t711 + t1146 * t737;
t707 = -pkin(6) * t725 - t1149 * t718 + t1152 * t722;
t706 = pkin(1) * t719 + t1180;
t705 = -qJ(2) * t719 - t1144 * t712 + t1146 * t707;
t704 = -pkin(1) * t725 + qJ(2) * t720 + t1144 * t707 + t1146 * t712;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1114, 0, -t1113, 0, t1164, -t1090, -t1161, -pkin(5) * t1161, 0, 0, t1056, 0, -t1240, 0, t1247, t1248, t1253, pkin(5) * t1253 + qJ(2) * t1218 - t1150 * t969, -t1039 * t1150 + t1041 * t1153, -t1006 * t1150 + t1008 * t1153, -t1035 * t1150 + t1037 * t1153, -t1038 * t1150 + t1040 * t1153, -t1034 * t1150 + t1036 * t1153, -t1072 * t1150 + t1073 * t1153, t1153 * t874 - t1150 * t863 - pkin(5) * (t1015 * t1153 + t1017 * t1150), t1153 * t875 - t1150 * t864 - pkin(5) * (t1016 * t1153 + t1018 * t1150), t1153 * t922 - t1150 * t921 - pkin(5) * (t1058 * t1153 + t1059 * t1150), t1153 * t834 - t1150 * t821 - pkin(5) * (t1150 * t917 + t1153 * t916), -t1150 * t903 + t1153 * t905, -t1150 * t860 + t1153 * t861, -t1150 * t882 + t1153 * t885, -t1150 * t902 + t1153 * t904, -t1150 * t881 + t1153 * t884, -t1150 * t954 + t1153 * t955, t1153 * t766 - t1150 * t752 - pkin(5) * (t1150 * t867 + t1153 * t866), t1153 * t774 - t1150 * t767 - pkin(5) * (t1150 * t886 + t1153 * t883), t1153 * t741 - t1150 * t740 - pkin(5) * (t1150 * t853 + t1153 * t852), t1153 * t731 - t1150 * t723 - pkin(5) * (t1150 * t778 + t1153 * t777), -t1150 * t784 + t1153 * t786, -t1150 * t750 + t1153 * t751, -t1150 * t789 + t1153 * t791, -t1150 * t783 + t1153 * t785, -t1150 * t790 + t1153 * t792, -t1150 * t827 + t1153 * t828, t1153 * t714 - t1150 * t713 - pkin(5) * (t1150 * t776 + t1153 * t775), t1153 * t716 - t1150 * t715 - pkin(5) * (t1150 * t788 + t1153 * t787), t1153 * t709 - t1150 * t708 - pkin(5) * (t1150 * t748 + t1153 * t747), t1153 * t705 - t1150 * t704 - pkin(5) * (t1150 * t720 + t1153 * t719); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1113, 0, t1114, 0, t1090, t1164, t1170, pkin(5) * t1170, 0, 0, t1240, 0, t1056, 0, -t1248, t1247, t1254, pkin(5) * t1254 + qJ(2) * t1222 + t1153 * t969, t1039 * t1153 + t1041 * t1150, t1006 * t1153 + t1008 * t1150, t1035 * t1153 + t1037 * t1150, t1038 * t1153 + t1040 * t1150, t1034 * t1153 + t1036 * t1150, t1072 * t1153 + t1073 * t1150, t1150 * t874 + t1153 * t863 + pkin(5) * (-t1015 * t1150 + t1017 * t1153), t1150 * t875 + t1153 * t864 + pkin(5) * (-t1016 * t1150 + t1018 * t1153), t1150 * t922 + t1153 * t921 + pkin(5) * (-t1058 * t1150 + t1059 * t1153), t1150 * t834 + t1153 * t821 + pkin(5) * (-t1150 * t916 + t1153 * t917), t1150 * t905 + t1153 * t903, t1150 * t861 + t1153 * t860, t1150 * t885 + t1153 * t882, t1150 * t904 + t1153 * t902, t1150 * t884 + t1153 * t881, t1150 * t955 + t1153 * t954, t1150 * t766 + t1153 * t752 + pkin(5) * (-t1150 * t866 + t1153 * t867), t1150 * t774 + t1153 * t767 + pkin(5) * (-t1150 * t883 + t1153 * t886), t1150 * t741 + t1153 * t740 + pkin(5) * (-t1150 * t852 + t1153 * t853), t1150 * t731 + t1153 * t723 + pkin(5) * (-t1150 * t777 + t1153 * t778), t1150 * t786 + t1153 * t784, t1150 * t751 + t1153 * t750, t1150 * t791 + t1153 * t789, t1150 * t785 + t1153 * t783, t1150 * t792 + t1153 * t790, t1150 * t828 + t1153 * t827, t1150 * t714 + t1153 * t713 + pkin(5) * (-t1150 * t775 + t1153 * t776), t1150 * t716 + t1153 * t715 + pkin(5) * (-t1150 * t787 + t1153 * t788), t1150 * t709 + t1153 * t708 + pkin(5) * (-t1150 * t747 + t1153 * t748), t1150 * t705 + t1153 * t704 + pkin(5) * (-t1150 * t719 + t1153 * t720); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1120, t1121, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1013, t1014, 0, -t977, t1063, t1054, t1066, t1062, t1064, 0, t937, t936, t906, t849, t940, t868, t931, t939, t930, t958, t779, t796, t744, t736, t798, t756, t808, t797, t809, t829, t727, t728, t710, t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1237, 0, 0, -g(3), -t1120, 0, 0, 0, t1111, 0, -t1110, 0, t1162, t1080, t981, qJ(2) * t981, t1041, t1008, t1037, t1040, t1036, t1073, t874, t875, t922, t834, t905, t861, t885, t904, t884, t955, t766, t774, t741, t731, t786, t751, t791, t785, t792, t828, t714, t716, t709, t705; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1237, 0, qJDD(1), 0, g(3), 0, -t1121, 0, 0, 0, t1110, 0, t1111, 0, -t1080, t1162, t1173, t969, t1039, t1006, t1035, t1038, t1034, t1072, t863, t864, t921, t821, t903, t860, t882, t902, t881, t954, t752, t767, t740, t723, t784, t750, t789, t783, t790, t827, t713, t715, t708, t704; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1120, t1121, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1013, t1014, 0, -t977, t1063, t1054, t1066, t1062, t1064, 0, t937, t936, t906, t849, t940, t868, t931, t939, t930, t958, t779, t796, t744, t736, t798, t756, t808, t797, t809, t829, t727, t728, t710, t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1237, 0, 0, -t1141, t1051, 0, t1075, t1055, t1070, t1074, t1068, t1101, t983, t984, t952, pkin(6) * t952, t942, t870, t934, t941, t933, t959, t817, t822, t760, t743, t800, t758, t810, t799, t811, t830, t734, t735, t711, t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1237, 0, qJDD(1), 0, t1141, 0, t1052, 0, t1127, -t1116, -t1193, -t1127, -t1192, -qJDD(3), t975, t976, 0, pkin(2) * t952, -t1050, -t1048, -t1005, t1050, -t1003, -qJDD(3), t826, t831, t838, t768, -t971, -t968, -t915, t971, t911, -t1138, t742, t745, t737, t712; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1051, -t1052, 0, 0, t1063, t1054, t1066, t1062, t1064, 0, t1182, t1183, t1169, t1217, t940, t868, t931, t939, t930, t958, t1175, t1174, t1176, t1158, t798, t756, t808, t797, t809, t829, t1178, t1177, t1179, t1180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1107, t1108, t1118, -t1131, t1125, t1131, 0, t1032, t1007, 0, t994, t948, t989, t992, t988, t1022, t907, t919, t820, -qJ(4) * t832, t842, t803, t856, t841, t857, t878, t765, t773, t733, t722; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1184, t1106, t1123, t1157, t1119, -t1184, -t1032, 0, t1009, 0, t993, t946, t986, t991, t985, t1021, t873, t879, t816, t825, t840, t801, t854, t839, t855, t877, t754, t763, t730, t718; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1127, t1116, t1193, t1127, t1192, qJDD(3), -t1007, -t1009, 0, 0, t1050, t1048, t1005, -t1050, t1003, qJDD(3), t1255, t1168 + 0.2e1 * t1211, t1233, t1234, t971, t968, t915, -t971, -t911, t1138, t1156, t1155, t1189, t1190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1053, -t1001, t1239, t1088, t1081, -t1088, 0, t978, -t1163, 0, t900, t846, t927, t898, t928, t957, t837, t848, t749, -pkin(7) * t770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1212, t1004, t1082, -t1171, t1045, -t1212, -t978, 0, t888, 0, t899, t844, t925, t897, t926, t956, t818, t823, t746, t761; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1050, t1048, t1005, -t1050, t1003, qJDD(3), t1163, -t888, 0, 0, t971, t968, t915, -t971, -t911, t1138, t1165, t1159, t843, t769; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t945, -t910, t1244, t1020, t1011, -t1020, 0, t896, t814, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1198, t1245, t1012, t944, t964, -t1198, -t896, 0, t815, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t971, t968, t915, -t971, -t911, t1138, -t814, -t815, 0, 0;];
m_new_reg = t1;
