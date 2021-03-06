% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 20:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRRPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:12:20
% EndTime: 2019-05-07 20:13:18
% DurationCPUTime: 48.11s
% Computational Cost: add. (516646->948), mult. (1045659->1464), div. (0->0), fcn. (784072->12), ass. (0->654)
t1159 = sin(pkin(11));
t1163 = sin(qJ(3));
t1164 = sin(qJ(2));
t1168 = cos(qJ(3));
t1169 = cos(qJ(2));
t1122 = (t1163 * t1169 + t1164 * t1168) * qJD(1);
t1232 = qJD(2) * t1169;
t1192 = qJD(1) * t1232;
t1208 = qJDD(1) * t1164;
t1130 = t1192 + t1208;
t1153 = t1169 * qJDD(1);
t1234 = qJD(1) * t1164;
t1199 = qJD(2) * t1234;
t1131 = t1153 - t1199;
t1188 = t1163 * t1130 - t1168 * t1131;
t1061 = -t1122 * qJD(3) - t1188;
t1058 = qJDD(4) - t1061;
t1162 = sin(qJ(4));
t1167 = cos(qJ(4));
t1212 = qJD(2) + qJD(3);
t1093 = t1162 * t1122 - t1167 * t1212;
t1094 = t1167 * t1122 + t1162 * t1212;
t1160 = cos(pkin(11));
t1043 = t1160 * t1093 + t1094 * t1159;
t1045 = -t1093 * t1159 + t1094 * t1160;
t995 = t1045 * t1043;
t1261 = t1058 - t995;
t1272 = t1159 * t1261;
t1271 = t1160 * t1261;
t1161 = sin(qJ(6));
t1176 = qJDD(6) + t1058;
t1166 = cos(qJ(6));
t991 = t1166 * t1043 + t1045 * t1161;
t993 = -t1043 * t1161 + t1045 * t1166;
t921 = t993 * t991;
t1260 = t1176 - t921;
t1270 = t1161 * t1260;
t1055 = t1094 * t1093;
t1259 = -t1055 + t1058;
t1269 = t1162 * t1259;
t1120 = -t1168 * t1169 * qJD(1) + t1163 * t1234;
t1084 = t1122 * t1120;
t1205 = qJDD(2) + qJDD(3);
t1258 = -t1084 + t1205;
t1268 = t1163 * t1258;
t1267 = t1166 * t1260;
t1266 = t1167 * t1259;
t1265 = t1168 * t1258;
t1062 = -t1120 * qJD(3) + t1168 * t1130 + t1163 * t1131;
t1111 = t1212 * t1120;
t1264 = t1062 - t1111;
t1204 = t1212 ^ 2;
t1009 = -t1093 * qJD(4) + t1167 * t1062 + t1162 * t1205;
t1177 = -t1162 * t1062 + t1167 * t1205;
t1173 = t1094 * qJD(4) - t1177;
t1189 = t1009 * t1159 + t1160 * t1173;
t940 = t1160 * t1009 - t1159 * t1173;
t846 = -qJD(6) * t991 - t1161 * t1189 + t1166 * t940;
t1116 = qJD(4) + t1120;
t1114 = qJD(6) + t1116;
t967 = t1114 * t991;
t1263 = t846 - t967;
t1016 = t1116 * t1043;
t914 = -t1016 - t940;
t1262 = -t1016 + t940;
t1072 = t1116 * t1093;
t987 = -t1072 - t1009;
t985 = -t1072 + t1009;
t1191 = t1161 * t940 + t1166 * t1189;
t808 = (qJD(6) - t1114) * t993 + t1191;
t988 = (-t1116 + qJD(4)) * t1094 - t1177;
t1158 = t1169 ^ 2;
t1165 = sin(qJ(1));
t1170 = cos(qJ(1));
t1140 = t1165 * g(1) - t1170 * g(2);
t1179 = qJDD(1) * pkin(1) + t1140;
t1180 = qJD(2) * pkin(2) - pkin(8) * t1234;
t1257 = qJD(1) ^ 2;
t1064 = t1131 * pkin(2) - t1180 * t1234 + (pkin(8) * t1158 + pkin(7)) * t1257 + t1179;
t989 = t991 ^ 2;
t990 = t993 ^ 2;
t1256 = t1043 ^ 2;
t1042 = t1045 ^ 2;
t1255 = t1093 ^ 2;
t1091 = t1094 ^ 2;
t1112 = t1114 ^ 2;
t1254 = t1116 ^ 2;
t1118 = t1120 ^ 2;
t1119 = t1122 ^ 2;
t1253 = pkin(3) * t1163;
t1252 = t1169 * g(3);
t1190 = t1212 * t1122;
t936 = -t1264 * pkin(9) + (-t1061 + t1190) * pkin(3) - t1064;
t1081 = pkin(3) * t1120 - pkin(9) * t1122;
t1141 = g(1) * t1170 + g(2) * t1165;
t1178 = qJDD(1) * pkin(7) - t1141;
t1174 = -pkin(1) * t1257 + t1178;
t1105 = -t1164 * g(3) + t1169 * t1174;
t1155 = t1158 * t1257;
t1053 = -pkin(2) * t1155 + t1131 * pkin(8) - qJD(2) * t1180 + t1105;
t1172 = -t1164 * t1178 - t1252 - t1130 * pkin(8) + qJDD(2) * pkin(2) + (pkin(8) * t1232 + (pkin(2) * t1169 + pkin(1)) * t1234) * qJD(1);
t999 = t1168 * t1053 + t1163 * t1172;
t946 = -pkin(3) * t1204 + pkin(9) * t1205 - t1120 * t1081 + t999;
t864 = t1162 * t946 - t1167 * t936;
t817 = pkin(4) * t1259 + qJ(5) * t987 - t864;
t1065 = pkin(4) * t1116 - qJ(5) * t1094;
t865 = t1162 * t936 + t1167 * t946;
t828 = -pkin(4) * t1255 - qJ(5) * t1173 - t1116 * t1065 + t865;
t745 = 0.2e1 * qJD(5) * t1045 + t1159 * t828 - t1160 * t817;
t714 = pkin(5) * t1261 + pkin(10) * t914 - t745;
t1010 = pkin(5) * t1116 - pkin(10) * t1045;
t746 = -0.2e1 * qJD(5) * t1043 + t1159 * t817 + t1160 * t828;
t722 = -pkin(5) * t1256 - pkin(10) * t1189 - t1010 * t1116 + t746;
t653 = t1161 * t714 + t1166 * t722;
t652 = t1161 * t722 - t1166 * t714;
t614 = t1161 * t653 - t1166 * t652;
t1251 = t1159 * t614;
t998 = t1053 * t1163 - t1168 * t1172;
t945 = -t1205 * pkin(3) - t1204 * pkin(9) + t1081 * t1122 + t998;
t866 = t1173 * pkin(4) - t1255 * qJ(5) + t1065 * t1094 + qJDD(5) + t945;
t1250 = t1159 * t866;
t950 = t1058 + t995;
t1249 = t1159 * t950;
t1248 = t1160 * t614;
t1247 = t1160 * t866;
t1246 = t1160 * t950;
t786 = pkin(5) * t1189 - pkin(10) * t1256 + t1010 * t1045 + t866;
t1245 = t1161 * t786;
t886 = t1176 + t921;
t1244 = t1161 * t886;
t686 = t1159 * t746 - t1160 * t745;
t1243 = t1162 * t686;
t1242 = t1162 * t945;
t924 = t1163 * t999 - t1168 * t998;
t1241 = t1164 * t924;
t1240 = t1166 * t786;
t1239 = t1166 * t886;
t1238 = t1167 * t686;
t1237 = t1167 * t945;
t1236 = t1169 * t924;
t1235 = qJD(1) * qJD(2);
t1157 = t1164 ^ 2;
t1233 = t1257 * t1157;
t1001 = t1055 + t1058;
t1231 = t1001 * t1162;
t1230 = t1001 * t1167;
t1229 = t1045 * t1116;
t1228 = t1058 * t1163;
t1227 = t1064 * t1163;
t1226 = t1064 * t1168;
t1078 = t1084 + t1205;
t1225 = t1078 * t1163;
t1224 = t1078 * t1168;
t1223 = t1114 * t1161;
t1222 = t1114 * t1166;
t1221 = t1116 * t1159;
t1220 = t1116 * t1160;
t1219 = t1116 * t1162;
t1218 = t1116 * t1167;
t1123 = pkin(7) * t1257 + t1179;
t1217 = t1123 * t1164;
t1216 = t1123 * t1169;
t1147 = t1169 * t1257 * t1164;
t1138 = qJDD(2) + t1147;
t1215 = t1138 * t1164;
t1139 = qJDD(2) - t1147;
t1214 = t1139 * t1164;
t1213 = t1139 * t1169;
t1209 = t1157 + t1158;
t1207 = qJDD(1) * t1165;
t1206 = qJDD(1) * t1170;
t1203 = t1163 * t921;
t1202 = t1168 * t921;
t1201 = -pkin(3) * t1168 - pkin(2);
t1198 = t1163 * t995;
t1197 = t1168 * t995;
t1196 = t1163 * t1055;
t1195 = t1168 * t1055;
t1194 = t1165 * t1084;
t1193 = t1170 * t1084;
t687 = t1159 * t745 + t1160 * t746;
t615 = t1161 * t652 + t1166 * t653;
t925 = t1163 * t998 + t1168 * t999;
t1104 = t1164 * t1174 + t1252;
t1052 = t1104 * t1164 + t1169 * t1105;
t1096 = -t1140 * t1165 - t1170 * t1141;
t1187 = t1165 * t1147;
t1186 = t1170 * t1147;
t1135 = -t1165 * t1257 + t1206;
t1185 = -pkin(6) * t1135 - g(3) * t1165;
t1184 = t1163 * t1111;
t1183 = t1163 * t1190;
t1182 = t1168 * t1111;
t1181 = t1168 * t1190;
t787 = t1162 * t865 - t1167 * t864;
t788 = t1162 * t864 + t1167 * t865;
t1051 = t1104 * t1169 - t1105 * t1164;
t1095 = t1140 * t1170 - t1141 * t1165;
t1027 = qJD(2) * t1122 - t1188;
t915 = t1189 - t1229;
t1171 = qJD(2) ^ 2;
t1145 = -t1155 - t1171;
t1144 = t1155 - t1171;
t1143 = -t1171 - t1233;
t1142 = t1171 - t1233;
t1137 = t1155 - t1233;
t1136 = t1155 + t1233;
t1134 = t1170 * t1257 + t1207;
t1133 = t1209 * qJDD(1);
t1132 = t1153 - 0.2e1 * t1199;
t1129 = 0.2e1 * t1192 + t1208;
t1127 = t1169 * t1138;
t1126 = t1209 * t1235;
t1117 = -pkin(6) * t1134 + g(3) * t1170;
t1109 = t1204 - t1119;
t1108 = t1118 - t1204;
t1107 = t1130 * t1169 - t1157 * t1235;
t1106 = -t1131 * t1164 - t1158 * t1235;
t1103 = -t1119 - t1204;
t1102 = -t1143 * t1164 - t1213;
t1101 = -t1142 * t1164 + t1127;
t1100 = t1145 * t1169 - t1215;
t1099 = t1144 * t1169 - t1214;
t1098 = t1143 * t1169 - t1214;
t1097 = t1145 * t1164 + t1127;
t1090 = t1133 * t1170 - t1136 * t1165;
t1089 = t1133 * t1165 + t1136 * t1170;
t1085 = -t1129 * t1164 + t1132 * t1169;
t1083 = -t1119 + t1118;
t1077 = -t1204 - t1118;
t1076 = t1102 * t1170 + t1129 * t1165;
t1075 = t1100 * t1170 - t1132 * t1165;
t1074 = t1102 * t1165 - t1129 * t1170;
t1073 = t1100 * t1165 + t1132 * t1170;
t1071 = -t1091 + t1254;
t1070 = -t1254 + t1255;
t1069 = -t1182 + t1183;
t1068 = -t1184 - t1181;
t1067 = -pkin(7) * t1098 - t1216;
t1066 = -pkin(7) * t1097 - t1217;
t1063 = -t1118 - t1119;
t1060 = -pkin(1) * t1098 + t1105;
t1059 = -pkin(1) * t1097 + t1104;
t1054 = t1168 * t1058;
t1049 = -t1091 + t1255;
t1041 = -t1091 - t1254;
t1037 = t1108 * t1168 - t1225;
t1036 = -t1109 * t1163 + t1265;
t1035 = t1108 * t1163 + t1224;
t1034 = t1109 * t1168 + t1268;
t1032 = -t1103 * t1163 - t1224;
t1031 = t1103 * t1168 - t1225;
t1029 = -t1062 - t1111;
t1025 = (0.2e1 * qJD(3) + qJD(2)) * t1122 + t1188;
t1024 = -t1254 - t1255;
t1023 = t1168 * t1062 - t1183;
t1022 = t1163 * t1062 + t1181;
t1021 = -t1163 * t1061 + t1182;
t1020 = t1168 * t1061 + t1184;
t1019 = t1052 * t1170 - t1123 * t1165;
t1018 = t1052 * t1165 + t1123 * t1170;
t1015 = t1077 * t1168 - t1268;
t1014 = t1077 * t1163 + t1265;
t1013 = -t1042 + t1254;
t1012 = -t1254 + t1256;
t1011 = t1091 + t1255;
t1006 = (-t1093 * t1167 + t1094 * t1162) * t1116;
t1005 = (t1093 * t1162 + t1094 * t1167) * t1116;
t1004 = -t1042 - t1254;
t1003 = -t1068 * t1164 + t1069 * t1169;
t996 = -pkin(8) * t1031 - t1226;
t994 = -t1042 + t1256;
t983 = (-qJD(4) - t1116) * t1094 + t1177;
t982 = -pkin(8) * t1014 - t1227;
t981 = -t1035 * t1164 + t1037 * t1169;
t980 = -t1034 * t1164 + t1036 * t1169;
t979 = t1009 * t1167 - t1094 * t1219;
t978 = -t1009 * t1162 - t1094 * t1218;
t977 = t1093 * t1218 + t1162 * t1173;
t976 = -t1093 * t1219 + t1167 * t1173;
t975 = -t1254 - t1256;
t974 = -t1031 * t1164 + t1032 * t1169;
t973 = t1031 * t1169 + t1032 * t1164;
t972 = t1027 * t1168 - t1029 * t1163;
t971 = -t1025 * t1168 - t1163 * t1264;
t970 = t1027 * t1163 + t1029 * t1168;
t969 = -t1025 * t1163 + t1168 * t1264;
t966 = -t990 + t1112;
t965 = t989 - t1112;
t964 = t1006 * t1168 + t1228;
t963 = t1006 * t1163 - t1054;
t962 = t1070 * t1167 - t1231;
t961 = -t1071 * t1162 + t1266;
t960 = -t1070 * t1162 - t1230;
t959 = -t1071 * t1167 - t1269;
t958 = -t1022 * t1164 + t1023 * t1169;
t957 = -t1020 * t1164 + t1021 * t1169;
t956 = (-t1043 * t1160 + t1045 * t1159) * t1116;
t955 = (-t1043 * t1159 - t1045 * t1160) * t1116;
t954 = -t1014 * t1164 + t1015 * t1169;
t953 = t1014 * t1169 + t1015 * t1164;
t952 = -t990 - t1112;
t948 = -t1041 * t1162 - t1230;
t947 = t1041 * t1167 - t1231;
t943 = t1024 * t1167 - t1269;
t942 = t1024 * t1162 + t1266;
t941 = -t1042 - t1256;
t935 = -pkin(2) * t1264 + pkin(8) * t1032 - t1227;
t932 = t1168 * t979 + t1196;
t931 = t1168 * t977 - t1196;
t930 = t1163 * t979 - t1195;
t929 = t1163 * t977 + t1195;
t928 = -pkin(2) * t1025 + pkin(8) * t1015 + t1226;
t927 = t1165 * t1264 + t1170 * t974;
t926 = t1165 * t974 - t1170 * t1264;
t923 = t1025 * t1165 + t1170 * t954;
t922 = -t1025 * t1170 + t1165 * t954;
t920 = -t990 + t989;
t919 = t1012 * t1160 - t1249;
t918 = -t1013 * t1159 + t1271;
t917 = t1012 * t1159 + t1246;
t916 = t1013 * t1160 + t1272;
t910 = t1189 + t1229;
t909 = -t1045 * t1221 + t1160 * t940;
t908 = t1045 * t1220 + t1159 * t940;
t907 = t1043 * t1220 + t1159 * t1189;
t906 = t1043 * t1221 - t1160 * t1189;
t905 = -t1162 * t987 - t1167 * t988;
t904 = -t1162 * t985 + t1167 * t983;
t903 = -t1162 * t988 + t1167 * t987;
t902 = -t1162 * t983 - t1167 * t985;
t901 = -t1004 * t1159 - t1246;
t900 = t1004 * t1160 - t1249;
t899 = -t1112 - t989;
t898 = pkin(2) * t1064 + pkin(8) * t925;
t897 = -t1164 * t970 + t1169 * t972;
t896 = -t1164 * t969 + t1169 * t971;
t895 = t1164 * t972 + t1169 * t970;
t894 = -t1163 * t988 + t1168 * t962;
t893 = -t1163 * t987 + t1168 * t961;
t892 = t1163 * t962 + t1168 * t988;
t891 = t1163 * t961 + t1168 * t987;
t890 = (t1161 * t993 - t1166 * t991) * t1114;
t889 = (-t1161 * t991 - t1166 * t993) * t1114;
t888 = -t1164 * t963 + t1169 * t964;
t884 = t1163 * t985 + t1168 * t948;
t883 = t1163 * t948 - t1168 * t985;
t882 = t1160 * t975 - t1272;
t881 = t1159 * t975 + t1271;
t880 = -t1163 * t983 + t1168 * t943;
t879 = t1163 * t943 + t1168 * t983;
t878 = -pkin(1) * t973 - pkin(2) * t1031 + t999;
t877 = -t1049 * t1163 + t1168 * t904;
t876 = t1049 * t1168 + t1163 * t904;
t875 = -t1162 * t955 + t1167 * t956;
t874 = -t1162 * t956 - t1167 * t955;
t873 = t1063 * t1165 + t1170 * t897;
t872 = -t1063 * t1170 + t1165 * t897;
t871 = -pkin(9) * t947 + t1237;
t870 = -t1011 * t1163 + t1168 * t905;
t869 = t1011 * t1168 + t1163 * t905;
t868 = -pkin(1) * t953 - pkin(2) * t1014 + t998;
t867 = -pkin(9) * t942 + t1242;
t863 = -pkin(8) * t970 - t924;
t862 = -t989 - t990;
t861 = t1168 * t875 + t1228;
t860 = t1163 * t875 - t1054;
t859 = -t1164 * t930 + t1169 * t932;
t858 = -t1164 * t929 + t1169 * t931;
t857 = -pkin(2) * t1063 + pkin(8) * t972 + t925;
t856 = -pkin(1) * t895 - pkin(2) * t970;
t855 = -pkin(7) * t973 - t1164 * t935 + t1169 * t996;
t854 = t1166 * t965 - t1244;
t853 = -t1161 * t966 + t1267;
t852 = t1161 * t965 + t1239;
t851 = t1166 * t966 + t1270;
t850 = -t1161 * t952 - t1239;
t849 = t1166 * t952 - t1244;
t848 = t1169 * t925 - t1241;
t847 = t1164 * t925 + t1236;
t845 = -qJD(6) * t993 - t1191;
t844 = -pkin(7) * t953 - t1164 * t928 + t1169 * t982;
t843 = -t1162 * t917 + t1167 * t919;
t842 = -t1162 * t916 + t1167 * t918;
t841 = -t1162 * t919 - t1167 * t917;
t840 = -t1162 * t918 - t1167 * t916;
t839 = -t1159 * t914 - t1160 * t915;
t838 = -t1159 * t1262 - t1160 * t910;
t837 = -t1159 * t915 + t1160 * t914;
t836 = -t1159 * t910 + t1160 * t1262;
t835 = -t1064 * t1165 + t1170 * t848;
t834 = t1064 * t1170 + t1165 * t848;
t833 = -t1162 * t908 + t1167 * t909;
t832 = -t1162 * t906 + t1167 * t907;
t831 = -t1162 * t909 - t1167 * t908;
t830 = -t1162 * t907 - t1167 * t906;
t829 = -pkin(3) * t947 + t865;
t827 = -pkin(3) * t942 + t864;
t826 = -t1162 * t900 + t1167 * t901;
t825 = t1162 * t901 + t1167 * t900;
t823 = t1166 * t899 - t1270;
t822 = t1161 * t899 + t1267;
t821 = -t1164 * t892 + t1169 * t894;
t820 = -t1164 * t891 + t1169 * t893;
t819 = -t1159 * t889 + t1160 * t890;
t818 = t1159 * t890 + t1160 * t889;
t814 = -t1164 * t883 + t1169 * t884;
t813 = t1164 * t884 + t1169 * t883;
t812 = -t1162 * t881 + t1167 * t882;
t811 = t1162 * t882 + t1167 * t881;
t810 = -t1164 * t879 + t1169 * t880;
t809 = t1164 * t880 + t1169 * t879;
t807 = -t846 - t967;
t803 = (qJD(6) + t1114) * t993 + t1191;
t802 = -qJ(5) * t900 + t1247;
t801 = t1166 * t846 - t1223 * t993;
t800 = t1161 * t846 + t1222 * t993;
t799 = -t1161 * t845 + t1222 * t991;
t798 = t1166 * t845 + t1223 * t991;
t797 = -t1164 * t876 + t1169 * t877;
t796 = t1168 * t833 + t1198;
t795 = t1168 * t832 - t1198;
t794 = t1163 * t833 - t1197;
t793 = t1163 * t832 + t1197;
t792 = -qJ(5) * t881 + t1250;
t791 = -t1164 * t869 + t1169 * t870;
t790 = t1164 * t870 + t1169 * t869;
t789 = -pkin(1) * t847 - pkin(2) * t924;
t785 = -t1164 * t860 + t1169 * t861;
t784 = -t1163 * t915 + t1168 * t843;
t783 = -t1163 * t914 + t1168 * t842;
t782 = t1163 * t843 + t1168 * t915;
t781 = t1163 * t842 + t1168 * t914;
t780 = t1165 * t947 + t1170 * t814;
t779 = t1165 * t814 - t1170 * t947;
t778 = t1163 * t1262 + t1168 * t826;
t777 = t1163 * t826 - t1168 * t1262;
t776 = t1165 * t942 + t1170 * t810;
t775 = t1165 * t810 - t1170 * t942;
t774 = t1163 * t910 + t1168 * t812;
t773 = t1163 * t812 - t1168 * t910;
t772 = -pkin(4) * t1262 + qJ(5) * t901 + t1250;
t771 = -t1159 * t852 + t1160 * t854;
t770 = -t1159 * t851 + t1160 * t853;
t769 = t1159 * t854 + t1160 * t852;
t768 = t1159 * t853 + t1160 * t851;
t767 = t1163 * t945 + t1168 * t788;
t766 = t1163 * t788 - t1168 * t945;
t765 = -pkin(4) * t910 + qJ(5) * t882 - t1247;
t764 = -t1159 * t849 + t1160 * t850;
t763 = t1159 * t850 + t1160 * t849;
t762 = t1165 * t903 + t1170 * t791;
t761 = t1165 * t791 - t1170 * t903;
t760 = -pkin(9) * t903 - t787;
t759 = -pkin(7) * t847 - pkin(8) * t1236 - t1164 * t898;
t758 = -t1162 * t837 + t1167 * t839;
t757 = -t1162 * t836 + t1167 * t838;
t756 = t1162 * t839 + t1167 * t837;
t755 = -t1162 * t838 - t1167 * t836;
t754 = -pkin(7) * t895 - t1164 * t857 + t1169 * t863;
t753 = -t1159 * t822 + t1160 * t823;
t752 = t1159 * t823 + t1160 * t822;
t751 = -t1163 * t994 + t1168 * t757;
t750 = t1163 * t757 + t1168 * t994;
t749 = -pkin(8) * t883 - t1163 * t829 + t1168 * t871;
t748 = -t1162 * t818 + t1167 * t819;
t747 = -t1162 * t819 - t1167 * t818;
t743 = -pkin(10) * t849 + t1240;
t742 = -pkin(8) * t879 - t1163 * t827 + t1168 * t867;
t741 = t1163 * t941 + t1168 * t758;
t740 = t1163 * t758 - t1168 * t941;
t739 = t1163 * t1176 + t1168 * t748;
t738 = t1163 * t748 - t1168 * t1176;
t737 = -t1161 * t807 - t1166 * t808;
t736 = -t1161 * t1263 - t1166 * t803;
t735 = -t1161 * t808 + t1166 * t807;
t734 = -t1161 * t803 + t1166 * t1263;
t733 = -t1159 * t800 + t1160 * t801;
t732 = -t1159 * t798 + t1160 * t799;
t731 = t1159 * t801 + t1160 * t800;
t730 = t1159 * t799 + t1160 * t798;
t729 = -t1164 * t794 + t1169 * t796;
t728 = -t1164 * t793 + t1169 * t795;
t727 = -pkin(10) * t822 + t1245;
t726 = -pkin(2) * t947 + pkin(8) * t884 + t1163 * t871 + t1168 * t829;
t725 = -pkin(1) * t813 - pkin(2) * t883 + pkin(3) * t985 - pkin(9) * t948 - t1242;
t724 = -pkin(2) * t942 + pkin(8) * t880 + t1163 * t867 + t1168 * t827;
t723 = -pkin(1) * t809 - pkin(2) * t879 - pkin(3) * t983 - pkin(9) * t943 + t1237;
t720 = -t1164 * t782 + t1169 * t784;
t719 = -t1164 * t781 + t1169 * t783;
t718 = -pkin(3) * t756 - pkin(4) * t837;
t717 = -t1164 * t777 + t1169 * t778;
t716 = t1164 * t778 + t1169 * t777;
t715 = -pkin(8) * t869 + t1168 * t760 + t1253 * t903;
t711 = -t1164 * t773 + t1169 * t774;
t710 = t1164 * t774 + t1169 * t773;
t709 = -t1162 * t769 + t1167 * t771;
t708 = -t1162 * t768 + t1167 * t770;
t707 = -t1162 * t771 - t1167 * t769;
t706 = -t1162 * t770 - t1167 * t768;
t705 = -t1164 * t766 + t1169 * t767;
t704 = t1164 * t767 + t1169 * t766;
t703 = -pkin(5) * t1263 + pkin(10) * t850 + t1245;
t702 = pkin(8) * t870 + t1163 * t760 + t1201 * t903;
t701 = -t1162 * t763 + t1167 * t764;
t700 = t1162 * t764 + t1167 * t763;
t699 = -pkin(5) * t803 + pkin(10) * t823 - t1240;
t698 = -pkin(3) * t825 - pkin(4) * t900 + t746;
t697 = -pkin(9) * t825 - t1162 * t772 + t1167 * t802;
t696 = -pkin(1) * t790 - pkin(2) * t869 - pkin(3) * t1011 - pkin(9) * t905 - t788;
t695 = -pkin(3) * t811 - pkin(4) * t881 + t745;
t694 = t1165 * t825 + t1170 * t717;
t693 = t1165 * t717 - t1170 * t825;
t692 = -pkin(9) * t811 - t1162 * t765 + t1167 * t792;
t691 = -t1162 * t752 + t1167 * t753;
t690 = t1162 * t753 + t1167 * t752;
t689 = -t1164 * t750 + t1169 * t751;
t688 = -pkin(8) * t766 + (-pkin(9) * t1168 + t1253) * t787;
t685 = t1165 * t811 + t1170 * t711;
t684 = t1165 * t711 - t1170 * t811;
t683 = -t1164 * t740 + t1169 * t741;
t682 = t1164 * t741 + t1169 * t740;
t681 = -t1163 * t808 + t1168 * t709;
t680 = -t1163 * t807 + t1168 * t708;
t679 = t1163 * t709 + t1168 * t808;
t678 = t1163 * t708 + t1168 * t807;
t677 = -t1164 * t738 + t1169 * t739;
t676 = t1163 * t1263 + t1168 * t701;
t675 = t1163 * t701 - t1168 * t1263;
t674 = -t1159 * t735 + t1160 * t737;
t673 = -t1159 * t734 + t1160 * t736;
t672 = t1159 * t737 + t1160 * t735;
t671 = t1159 * t736 + t1160 * t734;
t670 = t1165 * t787 + t1170 * t705;
t669 = t1165 * t705 - t1170 * t787;
t668 = -t1162 * t731 + t1167 * t733;
t667 = -t1162 * t730 + t1167 * t732;
t666 = -t1162 * t733 - t1167 * t731;
t665 = -t1162 * t732 - t1167 * t730;
t664 = -pkin(4) * t866 + qJ(5) * t687;
t663 = t1168 * t668 + t1203;
t662 = t1168 * t667 - t1203;
t661 = t1163 * t668 - t1202;
t660 = t1163 * t667 + t1202;
t659 = -qJ(5) * t837 - t686;
t658 = t1163 * t803 + t1168 * t691;
t657 = t1163 * t691 - t1168 * t803;
t656 = -pkin(4) * t941 + qJ(5) * t839 + t687;
t655 = pkin(8) * t767 + (-pkin(9) * t1163 + t1201) * t787;
t654 = -pkin(7) * t813 - t1164 * t726 + t1169 * t749;
t650 = -pkin(7) * t809 - t1164 * t724 + t1169 * t742;
t649 = t1165 * t756 + t1170 * t683;
t648 = t1165 * t683 - t1170 * t756;
t647 = -pkin(1) * t704 - pkin(2) * t766 + pkin(3) * t945 - pkin(9) * t788;
t646 = -qJ(5) * t763 - t1159 * t703 + t1160 * t743;
t645 = -pkin(7) * t790 - t1164 * t702 + t1169 * t715;
t644 = -qJ(5) * t752 - t1159 * t699 + t1160 * t727;
t643 = -pkin(4) * t1263 + qJ(5) * t764 + t1159 * t743 + t1160 * t703;
t642 = -pkin(1) * t716 - pkin(2) * t777 + pkin(3) * t1262 - pkin(9) * t826 - t1162 * t802 - t1167 * t772;
t641 = t1167 * t687 - t1243;
t640 = t1162 * t687 + t1238;
t639 = -pkin(8) * t777 - t1163 * t698 + t1168 * t697;
t638 = -t1164 * t679 + t1169 * t681;
t637 = -t1164 * t678 + t1169 * t680;
t636 = -pkin(4) * t803 + qJ(5) * t753 + t1159 * t727 + t1160 * t699;
t635 = -t1164 * t675 + t1169 * t676;
t634 = t1164 * t676 + t1169 * t675;
t633 = -t1162 * t672 + t1167 * t674;
t632 = -t1162 * t671 + t1167 * t673;
t631 = t1162 * t674 + t1167 * t672;
t630 = -t1162 * t673 - t1167 * t671;
t629 = -pkin(1) * t710 - pkin(2) * t773 + pkin(3) * t910 - pkin(9) * t812 - t1162 * t792 - t1167 * t765;
t628 = t1163 * t866 + t1168 * t641;
t627 = t1163 * t641 - t1168 * t866;
t626 = -pkin(8) * t773 - t1163 * t695 + t1168 * t692;
t625 = -t1163 * t920 + t1168 * t632;
t624 = t1163 * t632 + t1168 * t920;
t623 = -pkin(2) * t825 + pkin(8) * t778 + t1163 * t697 + t1168 * t698;
t622 = t1163 * t862 + t1168 * t633;
t621 = t1163 * t633 - t1168 * t862;
t620 = -t1164 * t661 + t1169 * t663;
t619 = -t1164 * t660 + t1169 * t662;
t618 = -pkin(2) * t811 + pkin(8) * t774 + t1163 * t692 + t1168 * t695;
t617 = -t1164 * t657 + t1169 * t658;
t616 = t1164 * t658 + t1169 * t657;
t613 = -pkin(3) * t700 - pkin(4) * t763 - pkin(5) * t849 + t653;
t612 = t1165 * t700 + t1170 * t635;
t611 = t1165 * t635 - t1170 * t700;
t610 = -pkin(9) * t756 - t1162 * t656 + t1167 * t659;
t609 = -pkin(5) * t786 + pkin(10) * t615;
t608 = -pkin(3) * t690 - pkin(4) * t752 - pkin(5) * t822 + t652;
t607 = -pkin(3) * t640 - pkin(4) * t686;
t606 = -pkin(7) * t704 - t1164 * t655 + t1169 * t688;
t605 = -pkin(10) * t735 - t614;
t604 = -pkin(5) * t862 + pkin(10) * t737 + t615;
t603 = t1165 * t690 + t1170 * t617;
t602 = t1165 * t617 - t1170 * t690;
t601 = -pkin(3) * t631 - pkin(4) * t672 - pkin(5) * t735;
t600 = -pkin(8) * t740 - t1163 * t718 + t1168 * t610;
t599 = -pkin(9) * t640 - qJ(5) * t1238 - t1162 * t664;
t598 = -t1164 * t627 + t1169 * t628;
t597 = t1164 * t628 + t1169 * t627;
t596 = -pkin(9) * t700 - t1162 * t643 + t1167 * t646;
t595 = -t1164 * t624 + t1169 * t625;
t594 = -pkin(2) * t756 + pkin(8) * t741 + t1163 * t610 + t1168 * t718;
t593 = -pkin(1) * t682 - pkin(2) * t740 + pkin(3) * t941 - pkin(9) * t758 - t1162 * t659 - t1167 * t656;
t592 = -t1164 * t621 + t1169 * t622;
t591 = t1164 * t622 + t1169 * t621;
t590 = -pkin(9) * t690 - t1162 * t636 + t1167 * t644;
t589 = -pkin(7) * t716 - t1164 * t623 + t1169 * t639;
t588 = t1160 * t615 - t1251;
t587 = t1159 * t615 + t1248;
t586 = -pkin(7) * t710 - t1164 * t618 + t1169 * t626;
t585 = t1165 * t640 + t1170 * t598;
t584 = t1165 * t598 - t1170 * t640;
t583 = t1165 * t631 + t1170 * t592;
t582 = t1165 * t592 - t1170 * t631;
t581 = -qJ(5) * t672 - t1159 * t604 + t1160 * t605;
t580 = -pkin(4) * t862 + qJ(5) * t674 + t1159 * t605 + t1160 * t604;
t579 = -pkin(1) * t634 - pkin(2) * t675 + pkin(3) * t1263 - pkin(9) * t701 - t1162 * t646 - t1167 * t643;
t578 = -pkin(8) * t675 - t1163 * t613 + t1168 * t596;
t577 = -pkin(1) * t616 - pkin(2) * t657 + pkin(3) * t803 - pkin(9) * t691 - t1162 * t644 - t1167 * t636;
t576 = -pkin(2) * t700 + pkin(8) * t676 + t1163 * t596 + t1168 * t613;
t575 = -pkin(8) * t657 - t1163 * t608 + t1168 * t590;
t574 = -pkin(8) * t627 - t1163 * t607 + t1168 * t599;
t573 = -pkin(7) * t682 - t1164 * t594 + t1169 * t600;
t572 = -pkin(2) * t690 + pkin(8) * t658 + t1163 * t590 + t1168 * t608;
t571 = -t1162 * t587 + t1167 * t588;
t570 = t1162 * t588 + t1167 * t587;
t569 = -pkin(1) * t597 - pkin(2) * t627 + pkin(3) * t866 - pkin(9) * t641 + qJ(5) * t1243 - t1167 * t664;
t568 = -pkin(10) * t1248 - qJ(5) * t587 - t1159 * t609;
t567 = t1163 * t786 + t1168 * t571;
t566 = t1163 * t571 - t1168 * t786;
t565 = -pkin(4) * t786 - pkin(10) * t1251 + qJ(5) * t588 + t1160 * t609;
t564 = -pkin(2) * t640 + pkin(8) * t628 + t1163 * t599 + t1168 * t607;
t563 = -pkin(9) * t631 - t1162 * t580 + t1167 * t581;
t562 = -pkin(3) * t570 - pkin(4) * t587 - pkin(5) * t614;
t561 = -pkin(7) * t634 - t1164 * t576 + t1169 * t578;
t560 = -t1164 * t566 + t1169 * t567;
t559 = t1164 * t567 + t1169 * t566;
t558 = -pkin(7) * t616 - t1164 * t572 + t1169 * t575;
t557 = -pkin(8) * t621 - t1163 * t601 + t1168 * t563;
t556 = -pkin(1) * t591 - pkin(2) * t621 + pkin(3) * t862 - pkin(9) * t633 - t1162 * t581 - t1167 * t580;
t555 = -pkin(7) * t597 - t1164 * t564 + t1169 * t574;
t554 = -pkin(2) * t631 + pkin(8) * t622 + t1163 * t563 + t1168 * t601;
t553 = t1165 * t570 + t1170 * t560;
t552 = t1165 * t560 - t1170 * t570;
t551 = -pkin(9) * t570 - t1162 * t565 + t1167 * t568;
t550 = -pkin(7) * t591 - t1164 * t554 + t1169 * t557;
t549 = -pkin(8) * t566 - t1163 * t562 + t1168 * t551;
t548 = -pkin(1) * t559 - pkin(2) * t566 + pkin(3) * t786 - pkin(9) * t571 - t1162 * t568 - t1167 * t565;
t547 = -pkin(2) * t570 + pkin(8) * t567 + t1163 * t551 + t1168 * t562;
t546 = -pkin(7) * t559 - t1164 * t547 + t1169 * t549;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1134, -t1135, 0, t1096, 0, 0, 0, 0, 0, 0, t1075, t1076, t1090, t1019, 0, 0, 0, 0, 0, 0, t923, t927, t873, t835, 0, 0, 0, 0, 0, 0, t776, t780, t762, t670, 0, 0, 0, 0, 0, 0, t685, t694, t649, t585, 0, 0, 0, 0, 0, 0, t603, t612, t583, t553; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1135, -t1134, 0, t1095, 0, 0, 0, 0, 0, 0, t1073, t1074, t1089, t1018, 0, 0, 0, 0, 0, 0, t922, t926, t872, t834, 0, 0, 0, 0, 0, 0, t775, t779, t761, t669, 0, 0, 0, 0, 0, 0, t684, t693, t648, t584, 0, 0, 0, 0, 0, 0, t602, t611, t582, t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1097, t1098, 0, -t1051, 0, 0, 0, 0, 0, 0, t953, t973, t895, t847, 0, 0, 0, 0, 0, 0, t809, t813, t790, t704, 0, 0, 0, 0, 0, 0, t710, t716, t682, t597, 0, 0, 0, 0, 0, 0, t616, t634, t591, t559; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1135, 0, -t1134, 0, t1185, -t1117, -t1095, -pkin(6) * t1095, t1107 * t1170 - t1187, t1085 * t1170 - t1137 * t1165, t1101 * t1170 + t1164 * t1207, t1106 * t1170 + t1187, t1099 * t1170 + t1153 * t1165, qJDD(2) * t1165 + t1126 * t1170, -pkin(6) * t1073 - t1059 * t1165 + t1066 * t1170, -pkin(6) * t1074 - t1060 * t1165 + t1067 * t1170, -pkin(6) * t1089 + t1051 * t1170, -pkin(6) * t1018 - (pkin(1) * t1165 - pkin(7) * t1170) * t1051, t1170 * t958 + t1194, -t1083 * t1165 + t1170 * t896, -t1029 * t1165 + t1170 * t980, t1170 * t957 - t1194, t1027 * t1165 + t1170 * t981, t1170 * t1003 + t1165 * t1205, -pkin(6) * t922 - t1165 * t868 + t1170 * t844, -pkin(6) * t926 - t1165 * t878 + t1170 * t855, -pkin(6) * t872 - t1165 * t856 + t1170 * t754, -pkin(6) * t834 - t1165 * t789 + t1170 * t759, -t1165 * t978 + t1170 * t859, -t1165 * t902 + t1170 * t797, -t1165 * t959 + t1170 * t820, -t1165 * t976 + t1170 * t858, -t1165 * t960 + t1170 * t821, -t1005 * t1165 + t1170 * t888, -pkin(6) * t775 - t1165 * t723 + t1170 * t650, -pkin(6) * t779 - t1165 * t725 + t1170 * t654, -pkin(6) * t761 - t1165 * t696 + t1170 * t645, -pkin(6) * t669 - t1165 * t647 + t1170 * t606, -t1165 * t831 + t1170 * t729, -t1165 * t755 + t1170 * t689, -t1165 * t840 + t1170 * t719, -t1165 * t830 + t1170 * t728, -t1165 * t841 + t1170 * t720, -t1165 * t874 + t1170 * t785, -pkin(6) * t684 - t1165 * t629 + t1170 * t586, -pkin(6) * t693 - t1165 * t642 + t1170 * t589, -pkin(6) * t648 - t1165 * t593 + t1170 * t573, -pkin(6) * t584 - t1165 * t569 + t1170 * t555, -t1165 * t666 + t1170 * t620, -t1165 * t630 + t1170 * t595, -t1165 * t706 + t1170 * t637, -t1165 * t665 + t1170 * t619, -t1165 * t707 + t1170 * t638, -t1165 * t747 + t1170 * t677, -pkin(6) * t602 - t1165 * t577 + t1170 * t558, -pkin(6) * t611 - t1165 * t579 + t1170 * t561, -pkin(6) * t582 - t1165 * t556 + t1170 * t550, -pkin(6) * t552 - t1165 * t548 + t1170 * t546; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1134, 0, t1135, 0, t1117, t1185, t1096, pkin(6) * t1096, t1107 * t1165 + t1186, t1085 * t1165 + t1137 * t1170, t1101 * t1165 - t1164 * t1206, t1106 * t1165 - t1186, t1099 * t1165 - t1153 * t1170, -qJDD(2) * t1170 + t1126 * t1165, pkin(6) * t1075 + t1059 * t1170 + t1066 * t1165, pkin(6) * t1076 + t1060 * t1170 + t1067 * t1165, pkin(6) * t1090 + t1051 * t1165, pkin(6) * t1019 - (-pkin(1) * t1170 - pkin(7) * t1165) * t1051, t1165 * t958 - t1193, t1083 * t1170 + t1165 * t896, t1029 * t1170 + t1165 * t980, t1165 * t957 + t1193, -t1027 * t1170 + t1165 * t981, t1165 * t1003 - t1170 * t1205, pkin(6) * t923 + t1165 * t844 + t1170 * t868, pkin(6) * t927 + t1165 * t855 + t1170 * t878, pkin(6) * t873 + t1165 * t754 + t1170 * t856, pkin(6) * t835 + t1165 * t759 + t1170 * t789, t1165 * t859 + t1170 * t978, t1165 * t797 + t1170 * t902, t1165 * t820 + t1170 * t959, t1165 * t858 + t1170 * t976, t1165 * t821 + t1170 * t960, t1005 * t1170 + t1165 * t888, pkin(6) * t776 + t1165 * t650 + t1170 * t723, pkin(6) * t780 + t1165 * t654 + t1170 * t725, pkin(6) * t762 + t1165 * t645 + t1170 * t696, pkin(6) * t670 + t1165 * t606 + t1170 * t647, t1165 * t729 + t1170 * t831, t1165 * t689 + t1170 * t755, t1165 * t719 + t1170 * t840, t1165 * t728 + t1170 * t830, t1165 * t720 + t1170 * t841, t1165 * t785 + t1170 * t874, pkin(6) * t685 + t1165 * t586 + t1170 * t629, pkin(6) * t694 + t1165 * t589 + t1170 * t642, pkin(6) * t649 + t1165 * t573 + t1170 * t593, pkin(6) * t585 + t1165 * t555 + t1170 * t569, t1165 * t620 + t1170 * t666, t1165 * t595 + t1170 * t630, t1165 * t637 + t1170 * t706, t1165 * t619 + t1170 * t665, t1165 * t638 + t1170 * t707, t1165 * t677 + t1170 * t747, pkin(6) * t603 + t1165 * t558 + t1170 * t577, pkin(6) * t612 + t1165 * t561 + t1170 * t579, pkin(6) * t583 + t1165 * t550 + t1170 * t556, pkin(6) * t553 + t1165 * t546 + t1170 * t548; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1140, t1141, 0, 0, (t1130 + t1192) * t1164, t1129 * t1169 + t1132 * t1164, t1142 * t1169 + t1215, (t1131 - t1199) * t1169, t1144 * t1164 + t1213, 0, pkin(1) * t1132 + pkin(7) * t1100 + t1216, -pkin(1) * t1129 + pkin(7) * t1102 - t1217, pkin(1) * t1136 + pkin(7) * t1133 + t1052, pkin(1) * t1123 + pkin(7) * t1052, t1022 * t1169 + t1023 * t1164, t1164 * t971 + t1169 * t969, t1034 * t1169 + t1036 * t1164, t1020 * t1169 + t1021 * t1164, t1035 * t1169 + t1037 * t1164, t1068 * t1169 + t1069 * t1164, -pkin(1) * t1025 + pkin(7) * t954 + t1164 * t982 + t1169 * t928, -pkin(1) * t1264 + pkin(7) * t974 + t1164 * t996 + t1169 * t935, -pkin(1) * t1063 + pkin(7) * t897 + t1164 * t863 + t1169 * t857, pkin(1) * t1064 + pkin(7) * t848 - pkin(8) * t1241 + t1169 * t898, t1164 * t932 + t1169 * t930, t1164 * t877 + t1169 * t876, t1164 * t893 + t1169 * t891, t1164 * t931 + t1169 * t929, t1164 * t894 + t1169 * t892, t1164 * t964 + t1169 * t963, -pkin(1) * t942 + pkin(7) * t810 + t1164 * t742 + t1169 * t724, -pkin(1) * t947 + pkin(7) * t814 + t1164 * t749 + t1169 * t726, -pkin(1) * t903 + pkin(7) * t791 + t1164 * t715 + t1169 * t702, -pkin(1) * t787 + pkin(7) * t705 + t1164 * t688 + t1169 * t655, t1164 * t796 + t1169 * t794, t1164 * t751 + t1169 * t750, t1164 * t783 + t1169 * t781, t1164 * t795 + t1169 * t793, t1164 * t784 + t1169 * t782, t1164 * t861 + t1169 * t860, -pkin(1) * t811 + pkin(7) * t711 + t1164 * t626 + t1169 * t618, -pkin(1) * t825 + pkin(7) * t717 + t1164 * t639 + t1169 * t623, -pkin(1) * t756 + pkin(7) * t683 + t1164 * t600 + t1169 * t594, -pkin(1) * t640 + pkin(7) * t598 + t1164 * t574 + t1169 * t564, t1164 * t663 + t1169 * t661, t1164 * t625 + t1169 * t624, t1164 * t680 + t1169 * t678, t1164 * t662 + t1169 * t660, t1164 * t681 + t1169 * t679, t1164 * t739 + t1169 * t738, -pkin(1) * t690 + pkin(7) * t617 + t1164 * t575 + t1169 * t572, -pkin(1) * t700 + pkin(7) * t635 + t1164 * t578 + t1169 * t576, -pkin(1) * t631 + pkin(7) * t592 + t1164 * t557 + t1169 * t554, -pkin(1) * t570 + pkin(7) * t560 + t1164 * t549 + t1169 * t547;];
tauB_reg  = t1;
