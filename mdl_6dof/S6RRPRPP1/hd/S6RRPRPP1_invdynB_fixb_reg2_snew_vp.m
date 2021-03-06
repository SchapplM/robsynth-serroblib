% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRPRPP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRPRPP1_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:19:38
% EndTime: 2019-05-06 12:20:30
% DurationCPUTime: 54.13s
% Computational Cost: add. (165946->818), mult. (386689->1174), div. (0->0), fcn. (281488->10), ass. (0->593)
t1132 = sin(qJ(1));
t1135 = cos(qJ(1));
t1131 = sin(qJ(2));
t1134 = cos(qJ(2));
t1126 = sin(pkin(9));
t1128 = cos(pkin(9));
t1125 = sin(pkin(10));
t1127 = cos(pkin(10));
t1214 = qJD(1) * t1134;
t1215 = qJD(1) * t1131;
t1086 = t1126 * t1214 + t1128 * t1215;
t1130 = sin(qJ(4));
t1133 = cos(qJ(4));
t1061 = qJD(2) * t1130 + t1086 * t1133;
t1170 = qJD(2) * t1214;
t1183 = qJDD(1) * t1131;
t1094 = t1170 + t1183;
t1117 = t1134 * qJDD(1);
t1178 = qJD(2) * t1215;
t1095 = t1117 - t1178;
t1049 = t1094 * t1128 + t1095 * t1126;
t1153 = t1133 * qJDD(2) - t1130 * t1049;
t1144 = qJD(4) * t1061 - t1153;
t1060 = -t1133 * qJD(2) + t1086 * t1130;
t984 = -t1060 * qJD(4) + t1130 * qJDD(2) + t1133 * t1049;
t1139 = -t1125 * t1144 + t1127 * t984;
t1002 = t1127 * t1060 + t1061 * t1125;
t1084 = t1126 * t1215 - t1128 * t1214;
t1080 = qJD(4) + t1084;
t1207 = t1002 * t1080;
t1259 = t1139 - t1207;
t1163 = t1094 * t1126 - t1128 * t1095;
t1046 = qJDD(4) + t1163;
t1004 = -t1060 * t1125 + t1061 * t1127;
t947 = t1004 * t1002;
t1263 = t1046 + t947;
t1234 = t1125 * t1263;
t1000 = t1004 ^ 2;
t1244 = t1080 ^ 2;
t929 = t1244 + t1000;
t826 = t1127 * t929 + t1234;
t1228 = t1127 * t1263;
t847 = t1125 * t929 - t1228;
t774 = t1130 * t826 + t1133 * t847;
t717 = t1126 * t774 - t1128 * t1259;
t719 = t1126 * t1259 + t1128 * t774;
t653 = t1131 * t717 - t1134 * t719;
t772 = t1130 * t847 - t1133 * t826;
t626 = t1132 * t653 + t1135 * t772;
t1396 = pkin(6) * t626;
t628 = -t1132 * t772 + t1135 * t653;
t1395 = pkin(6) * t628;
t660 = t1131 * t719 + t1134 * t717;
t1394 = pkin(7) * t660;
t1393 = pkin(1) * t660 + pkin(2) * t717 - pkin(3) * t1259 + pkin(8) * t774;
t1392 = pkin(1) * t772 + pkin(7) * t653;
t1246 = t1002 ^ 2;
t971 = t1246 - t1244;
t852 = t1125 * t971 + t1228;
t856 = t1127 * t971 - t1234;
t783 = t1130 * t852 - t1133 * t856;
t1167 = t1125 * t984 + t1127 * t1144;
t1201 = t1080 * t1004;
t876 = -t1167 + t1201;
t731 = t1126 * t783 + t1128 * t876;
t735 = -t1126 * t876 + t1128 * t783;
t667 = t1131 * t731 - t1134 * t735;
t779 = t1130 * t856 + t1133 * t852;
t1391 = t1132 * t667 - t1135 * t779;
t1148 = t1167 + t1201;
t796 = -t1125 * t1148 + t1127 * t1259;
t1236 = t1125 * t1259;
t798 = t1127 * t1148 + t1236;
t703 = t1130 * t796 + t1133 * t798;
t940 = t1000 - t1246;
t692 = t1126 * t703 + t1128 * t940;
t694 = -t1126 * t940 + t1128 * t703;
t636 = t1131 * t692 - t1134 * t694;
t702 = t1130 * t798 - t1133 * t796;
t1390 = t1132 * t636 + t1135 * t702;
t1389 = t1132 * t779 + t1135 * t667;
t1388 = -t1132 * t702 + t1135 * t636;
t1387 = qJ(3) * t717;
t1385 = pkin(2) * t772 - qJ(3) * t719;
t1384 = t1131 * t694 + t1134 * t692;
t1383 = t1131 * t735 + t1134 * t731;
t1380 = pkin(8) * t772;
t1257 = -t1207 - t1139;
t1273 = -t1125 * t1257 + t1127 * t876;
t1275 = t1125 * t876 + t1127 * t1257;
t1298 = t1130 * t1273 + t1133 * t1275;
t1299 = -t1130 * t1275 + t1133 * t1273;
t890 = -t1246 - t1000;
t1318 = t1126 * t890 + t1128 * t1299;
t1321 = t1126 * t1299 - t1128 * t890;
t1346 = -t1131 * t1321 + t1134 * t1318;
t1357 = t1132 * t1298 + t1135 * t1346;
t1379 = pkin(6) * t1357;
t1264 = t1046 - t947;
t1233 = t1125 * t1264;
t1256 = -t1244 - t1246;
t1269 = t1127 * t1256 - t1233;
t908 = t1127 * t1264;
t1274 = t1125 * t1256 + t908;
t1300 = t1130 * t1269 + t1133 * t1274;
t1301 = -t1130 * t1274 + t1133 * t1269;
t1320 = t1126 * t1148 + t1128 * t1301;
t1323 = t1126 * t1301 - t1128 * t1148;
t1342 = -t1131 * t1323 + t1134 * t1320;
t1358 = t1132 * t1300 + t1135 * t1342;
t1378 = pkin(6) * t1358;
t1360 = t1132 * t1346 - t1135 * t1298;
t1377 = pkin(6) * t1360;
t1361 = t1132 * t1342 - t1135 * t1300;
t1376 = pkin(6) * t1361;
t1369 = pkin(3) * t772 - pkin(4) * t826;
t1341 = t1131 * t1320 + t1134 * t1323;
t1367 = pkin(7) * t1341;
t1345 = t1131 * t1318 + t1134 * t1321;
t1366 = pkin(7) * t1345;
t1365 = -pkin(1) * t1341 - pkin(2) * t1323 + pkin(3) * t1148 - pkin(8) * t1301;
t1364 = -pkin(1) * t1345 - pkin(2) * t1321 + pkin(3) * t890 - pkin(8) * t1299;
t1363 = -pkin(1) * t1300 + pkin(7) * t1342;
t1362 = -pkin(1) * t1298 + pkin(7) * t1346;
t972 = -t1000 + t1244;
t1288 = t1127 * t972 + t1233;
t1289 = -t1125 * t972 + t908;
t1297 = -t1130 * t1289 - t1133 * t1288;
t1296 = -t1130 * t1288 + t1133 * t1289;
t1319 = -t1126 * t1257 + t1128 * t1296;
t1322 = t1126 * t1296 + t1128 * t1257;
t1344 = -t1131 * t1322 + t1134 * t1319;
t1359 = t1132 * t1344 + t1135 * t1297;
t1356 = -t1132 * t1297 + t1135 * t1344;
t1354 = qJ(5) * t826;
t1353 = qJ(5) * t847;
t1352 = qJ(3) * t1321;
t1351 = qJ(3) * t1323;
t1348 = -pkin(2) * t1298 + qJ(3) * t1318;
t1347 = -pkin(2) * t1300 + qJ(3) * t1320;
t1343 = t1131 * t1319 + t1134 * t1322;
t1336 = pkin(8) * t1298;
t1335 = pkin(8) * t1300;
t670 = -pkin(3) * t1298 - pkin(4) * t1275;
t1326 = -pkin(3) * t1300 - pkin(4) * t1274;
t1145 = (-t1002 * t1125 - t1004 * t1127) * t1080;
t1199 = t1080 * t1127;
t1175 = t1002 * t1199;
t1200 = t1080 * t1125;
t969 = t1004 * t1200;
t1154 = t969 - t1175;
t1253 = -t1130 * t1154 - t1133 * t1145;
t1202 = t1046 * t1126;
t1252 = -t1130 * t1145 + t1133 * t1154;
t1267 = t1128 * t1252 + t1202;
t1036 = t1128 * t1046;
t1271 = t1126 * t1252 - t1036;
t1295 = -t1131 * t1271 + t1134 * t1267;
t1317 = t1132 * t1295 + t1135 * t1253;
t1147 = t1125 * t1167 + t1175;
t1155 = t1002 * t1200 - t1127 * t1167;
t1250 = -t1130 * t1147 - t1133 * t1155;
t1177 = t1126 * t947;
t1251 = -t1130 * t1155 + t1133 * t1147;
t1268 = t1128 * t1251 - t1177;
t1176 = t1128 * t947;
t1270 = t1126 * t1251 + t1176;
t1293 = -t1131 * t1270 + t1134 * t1268;
t1316 = t1132 * t1293 + t1135 * t1250;
t1315 = -t1132 * t1253 + t1135 * t1295;
t1314 = -t1132 * t1250 + t1135 * t1293;
t1310 = qJ(5) * t1269;
t1309 = qJ(5) * t1274;
t1308 = qJ(5) * t1275;
t1307 = qJ(6) * t1259;
t1302 = -pkin(4) * t890 + qJ(5) * t1273;
t1294 = t1131 * t1267 + t1134 * t1271;
t1292 = t1131 * t1268 + t1134 * t1270;
t1291 = 2 * qJD(6);
t1047 = t1086 * t1084;
t1255 = -t1047 + qJDD(2);
t1287 = t1126 * t1255;
t1284 = t1128 * t1255;
t1016 = t1061 * t1060;
t1260 = -t1016 + t1046;
t1281 = t1130 * t1260;
t1278 = t1133 * t1260;
t1212 = qJD(2) * t1086;
t1017 = t1163 + t1212;
t865 = t1004 * t1199 + t1125 * t1139;
t866 = t1127 * t1139 - t969;
t793 = -t1130 * t865 + t1133 * t866;
t1156 = t1126 * t793 - t1176;
t1157 = t1128 * t793 + t1177;
t1248 = -t1131 * t1156 + t1134 * t1157;
t792 = -t1130 * t866 - t1133 * t865;
t1272 = t1132 * t1248 + t1135 * t792;
t1266 = -t1132 * t792 + t1135 * t1248;
t1035 = t1080 * t1060;
t952 = -t1035 - t984;
t953 = -t984 + t1035;
t1037 = pkin(3) * t1084 - pkin(8) * t1086;
t1247 = qJD(2) ^ 2;
t1105 = g(1) * t1135 + g(2) * t1132;
t1136 = qJD(1) ^ 2;
t1143 = -pkin(1) * t1136 + qJDD(1) * pkin(7) - t1105;
t1069 = -t1131 * g(3) + t1134 * t1143;
t1123 = t1134 ^ 2;
t1120 = t1123 * t1136;
t1149 = qJD(2) * pkin(2) - qJ(3) * t1215;
t1009 = -pkin(2) * t1120 + t1095 * qJ(3) - qJD(2) * t1149 + t1069;
t1141 = t1131 * t1143;
t1186 = t1131 * t1136;
t1216 = qJD(1) * qJD(2);
t1137 = -t1141 - t1094 * qJ(3) + qJDD(2) * pkin(2) + (pkin(2) * t1186 + qJ(3) * t1216 - g(3)) * t1134;
t922 = -0.2e1 * qJD(3) * t1084 + t1128 * t1009 + t1126 * t1137;
t884 = -pkin(3) * t1247 + qJDD(2) * pkin(8) - t1037 * t1084 + t922;
t1104 = t1132 * g(1) - t1135 * g(2);
t1152 = qJDD(1) * pkin(1) + t1104;
t1014 = t1095 * pkin(2) + (qJ(3) * t1123 + pkin(7)) * t1136 - t1149 * t1215 - qJDD(3) + t1152;
t1213 = qJD(2) * t1084;
t1164 = -t1049 + t1213;
t907 = pkin(3) * t1017 + t1164 * pkin(8) - t1014;
t812 = t1130 * t884 - t1133 * t907;
t762 = pkin(4) * t1260 + qJ(5) * t952 - t812;
t1024 = pkin(4) * t1080 - qJ(5) * t1061;
t1245 = t1060 ^ 2;
t813 = t1130 * t907 + t1133 * t884;
t768 = -pkin(4) * t1245 - qJ(5) * t1144 - t1080 * t1024 + t813;
t1241 = t1125 * t762 + t1127 * t768;
t938 = pkin(5) * t1002 - qJ(6) * t1004;
t1258 = t1046 * qJ(6) - t1002 * t938 + t1080 * t1291 + t1241;
t1254 = pkin(5) * t1167 - t1307;
t948 = (qJD(4) - t1080) * t1061 - t1153;
t1249 = t1131 * t1157 + t1134 * t1156;
t1059 = t1061 ^ 2;
t1082 = t1084 ^ 2;
t1083 = t1086 ^ 2;
t1243 = pkin(3) * t1126;
t1242 = pkin(5) * t1127;
t1240 = t1125 * t768 - t1127 * t762;
t1239 = qJ(6) * t1127;
t1166 = t1126 * t1009 - t1128 * t1137;
t1151 = -qJDD(2) * pkin(3) - t1247 * pkin(8) + t1166;
t1140 = t1144 * pkin(4) - t1245 * qJ(5) + t1061 * t1024 + qJDD(5) + t1151;
t1165 = (0.2e1 * qJD(3) + t1037) * t1086;
t810 = t1165 + t1140;
t1238 = t1125 * t810;
t1230 = t1127 * t810;
t1209 = qJD(5) * t1004;
t995 = 0.2e1 * t1209;
t683 = t995 + t1240;
t1210 = qJD(5) * t1002;
t993 = -0.2e1 * t1210;
t684 = t993 + t1241;
t624 = t1125 * t684 - t1127 * t683;
t1226 = t1130 * t624;
t883 = t1165 + t1151;
t1225 = t1130 * t883;
t964 = t1016 + t1046;
t1224 = t1130 * t964;
t1211 = qJD(3) * t1086;
t921 = t1166 + 0.2e1 * t1211;
t832 = t1126 * t922 - t1128 * t921;
t1223 = t1131 * t832;
t1221 = t1133 * t624;
t1220 = t1133 * t883;
t1219 = t1133 * t964;
t1218 = t1134 * t832;
t1217 = -t1244 - t890;
t1206 = t1014 * t1126;
t1205 = t1014 * t1128;
t1040 = qJDD(2) + t1047;
t1204 = t1040 * t1126;
t1203 = t1040 * t1128;
t1198 = t1080 * t1130;
t1197 = t1080 * t1133;
t1196 = t1084 * t1126;
t1195 = t1084 * t1128;
t1194 = t1086 * t1126;
t1193 = t1086 * t1128;
t1088 = t1136 * pkin(7) + t1152;
t1192 = t1088 * t1131;
t1191 = t1088 * t1134;
t1113 = t1134 * t1186;
t1102 = qJDD(2) + t1113;
t1190 = t1102 * t1131;
t1103 = qJDD(2) - t1113;
t1189 = t1103 * t1131;
t1188 = t1103 * t1134;
t1122 = t1131 ^ 2;
t1187 = t1122 * t1136;
t1184 = t1122 + t1123;
t1182 = qJDD(1) * t1132;
t1181 = qJDD(1) * t1135;
t1180 = qJDD(2) * t1135;
t1179 = -pkin(3) * t1128 - pkin(2);
t1174 = t1126 * t1016;
t1173 = t1128 * t1016;
t1172 = t1132 * t1047;
t1171 = t1135 * t1047;
t1169 = -qJ(6) * t1125 - pkin(4);
t625 = t1125 * t683 + t1127 * t684;
t833 = t1126 * t921 + t1128 * t922;
t1068 = t1134 * g(3) + t1141;
t1012 = t1068 * t1131 + t1134 * t1069;
t1058 = -t1104 * t1132 - t1135 * t1105;
t1161 = t1004 * t938 + qJDD(6) + t1240;
t1160 = t1132 * t1113;
t1159 = t1135 * t1113;
t1099 = -t1132 * t1136 + t1181;
t1158 = -pkin(6) * t1099 - g(3) * t1132;
t725 = t1130 * t813 - t1133 * t812;
t726 = t1130 * t812 + t1133 * t813;
t1011 = t1068 * t1134 - t1069 * t1131;
t1057 = t1104 * t1135 - t1105 * t1132;
t1150 = t993 + t1258;
t1146 = -t1046 * pkin(5) + t1161;
t1019 = -t1163 + t1212;
t1138 = t1004 * t1291 - t1086 * t1037 - t1140 - 0.2e1 * t1211 - t1254;
t1116 = t1132 * qJDD(2);
t1112 = -t1120 - t1247;
t1111 = t1120 - t1247;
t1110 = -t1187 - t1247;
t1109 = -t1187 + t1247;
t1101 = t1120 - t1187;
t1100 = t1120 + t1187;
t1098 = t1135 * t1136 + t1182;
t1097 = t1184 * qJDD(1);
t1096 = t1117 - 0.2e1 * t1178;
t1093 = 0.2e1 * t1170 + t1183;
t1091 = t1134 * t1102;
t1090 = t1184 * t1216;
t1081 = -pkin(6) * t1098 + g(3) * t1135;
t1075 = -t1083 - t1247;
t1074 = -t1083 + t1247;
t1073 = t1082 - t1247;
t1071 = t1094 * t1134 - t1122 * t1216;
t1070 = -t1095 * t1131 - t1123 * t1216;
t1067 = -t1110 * t1131 - t1188;
t1066 = -t1109 * t1131 + t1091;
t1065 = t1112 * t1134 - t1190;
t1064 = t1111 * t1134 - t1189;
t1063 = t1110 * t1134 - t1189;
t1062 = t1112 * t1131 + t1091;
t1053 = t1097 * t1135 - t1100 * t1132;
t1052 = t1097 * t1132 + t1100 * t1135;
t1050 = -t1093 * t1131 + t1096 * t1134;
t1045 = -t1083 + t1082;
t1038 = -t1082 - t1247;
t1034 = (t1194 - t1195) * qJD(2);
t1033 = (-t1193 - t1196) * qJD(2);
t1032 = t1067 * t1135 + t1093 * t1132;
t1031 = t1065 * t1135 - t1096 * t1132;
t1030 = t1067 * t1132 - t1093 * t1135;
t1029 = t1065 * t1132 + t1096 * t1135;
t1028 = -t1059 + t1244;
t1027 = -t1244 + t1245;
t1026 = -pkin(7) * t1063 - t1191;
t1025 = -pkin(7) * t1062 - t1192;
t1023 = -pkin(1) * t1063 + t1069;
t1022 = -pkin(1) * t1062 + t1068;
t1021 = -t1049 - t1213;
t1015 = -t1082 - t1083;
t1013 = -t1059 + t1245;
t1008 = -qJD(2) * t1194 + t1049 * t1128;
t1007 = qJD(2) * t1193 + t1049 * t1126;
t1006 = qJD(2) * t1195 + t1126 * t1163;
t1005 = qJD(2) * t1196 - t1128 * t1163;
t991 = -t1059 - t1244;
t990 = -t1075 * t1126 - t1203;
t989 = -t1074 * t1126 + t1284;
t988 = t1073 * t1128 - t1204;
t987 = t1075 * t1128 - t1204;
t986 = t1074 * t1128 + t1287;
t985 = t1073 * t1126 + t1203;
t982 = -t1244 - t1245;
t976 = t1059 + t1245;
t975 = t1012 * t1135 - t1088 * t1132;
t974 = t1012 * t1132 + t1088 * t1135;
t968 = t1038 * t1128 - t1287;
t967 = t1038 * t1126 + t1284;
t962 = (t1060 * t1130 + t1061 * t1133) * t1080;
t961 = (-t1060 * t1133 + t1061 * t1130) * t1080;
t960 = -t1033 * t1131 + t1034 * t1134;
t957 = t1019 * t1128 - t1021 * t1126;
t956 = -t1017 * t1128 + t1126 * t1164;
t955 = t1019 * t1126 + t1021 * t1128;
t954 = -t1017 * t1126 - t1128 * t1164;
t949 = (-qJD(4) - t1080) * t1061 + t1153;
t945 = -t1061 * t1198 + t1133 * t984;
t944 = -t1060 * t1198 + t1133 * t1144;
t943 = -t1061 * t1197 - t1130 * t984;
t942 = t1060 * t1197 + t1130 * t1144;
t939 = -qJ(3) * t987 - t1205;
t937 = -t1007 * t1131 + t1008 * t1134;
t936 = -t1005 * t1131 + t1006 * t1134;
t935 = -t1131 * t987 + t1134 * t990;
t934 = -t1131 * t986 + t1134 * t989;
t933 = -t1131 * t985 + t1134 * t988;
t932 = t1131 * t990 + t1134 * t987;
t931 = t1128 * t961 + t1202;
t930 = t1126 * t961 - t1036;
t927 = -t1028 * t1133 - t1281;
t926 = t1027 * t1133 - t1224;
t925 = -t1027 * t1130 - t1219;
t924 = -t1028 * t1130 + t1278;
t923 = -qJ(3) * t967 - t1206;
t910 = -t1130 * t991 - t1219;
t909 = t1133 * t991 - t1224;
t900 = t1133 * t982 - t1281;
t899 = t1130 * t982 + t1278;
t898 = -t1131 * t967 + t1134 * t968;
t897 = t1131 * t968 + t1134 * t967;
t896 = pkin(2) * t1164 + qJ(3) * t990 - t1206;
t895 = t1128 * t945 + t1174;
t894 = t1128 * t942 - t1174;
t893 = t1126 * t945 - t1173;
t892 = t1126 * t942 + t1173;
t887 = -pkin(2) * t1017 + qJ(3) * t968 + t1205;
t886 = -t1132 * t1164 + t1135 * t935;
t885 = t1132 * t935 + t1135 * t1164;
t872 = t1017 * t1132 + t1135 * t898;
t871 = -t1131 * t955 + t1134 * t957;
t870 = -t1131 * t954 + t1134 * t956;
t869 = -t1017 * t1135 + t1132 * t898;
t868 = t1131 * t957 + t1134 * t955;
t860 = -t1130 * t949 + t1133 * t953;
t859 = -t1130 * t952 - t1133 * t948;
t858 = t1130 * t953 + t1133 * t949;
t857 = -t1130 * t948 + t1133 * t952;
t844 = -t1126 * t948 + t1128 * t926;
t843 = -t1126 * t952 + t1128 * t924;
t842 = t1126 * t926 + t1128 * t948;
t841 = t1126 * t924 + t1128 * t952;
t840 = -t1126 * t953 + t1128 * t910;
t839 = t1126 * t910 + t1128 * t953;
t838 = -t1131 * t930 + t1134 * t931;
t837 = -t1126 * t949 + t1128 * t900;
t836 = t1126 * t900 + t1128 * t949;
t835 = t1015 * t1132 + t1135 * t871;
t834 = -t1015 * t1135 + t1132 * t871;
t831 = -t1013 * t1126 + t1128 * t858;
t830 = t1013 * t1128 + t1126 * t858;
t825 = -t1126 * t976 + t1128 * t859;
t824 = t1126 * t859 + t1128 * t976;
t819 = -pkin(1) * t868 - pkin(2) * t955;
t818 = -pkin(1) * t932 - pkin(2) * t987 + t922;
t817 = pkin(2) * t1014 + qJ(3) * t833;
t816 = -pkin(8) * t909 + t1220;
t815 = -t1131 * t893 + t1134 * t895;
t814 = -t1131 * t892 + t1134 * t894;
t811 = -pkin(8) * t899 + t1225;
t805 = -pkin(1) * t897 - pkin(2) * t967 + t921;
t804 = -qJ(3) * t955 - t832;
t803 = -pkin(7) * t932 - t1131 * t896 + t1134 * t939;
t802 = -pkin(2) * t1015 + qJ(3) * t957 + t833;
t776 = -pkin(7) * t897 - t1131 * t887 + t1134 * t923;
t771 = -pkin(3) * t909 + t813;
t770 = -t1131 * t842 + t1134 * t844;
t769 = -t1131 * t841 + t1134 * t843;
t767 = -pkin(3) * t899 + t812;
t764 = -t1131 * t839 + t1134 * t840;
t763 = t1131 * t840 + t1134 * t839;
t759 = -t1131 * t836 + t1134 * t837;
t758 = t1131 * t837 + t1134 * t836;
t757 = t1134 * t833 - t1223;
t756 = t1131 * t833 + t1218;
t747 = -t1131 * t830 + t1134 * t831;
t742 = t1230 + t1354;
t741 = -t1014 * t1132 + t1135 * t757;
t740 = t1014 * t1135 + t1132 * t757;
t739 = -t1131 * t824 + t1134 * t825;
t738 = t1131 * t825 + t1134 * t824;
t737 = t1238 - t1309;
t728 = t1132 * t909 + t1135 * t764;
t727 = t1132 * t764 - t1135 * t909;
t724 = t1132 * t899 + t1135 * t759;
t723 = t1132 * t759 - t1135 * t899;
t712 = -pkin(4) * t1259 + t1238 + t1353;
t711 = -pkin(1) * t756 - pkin(2) * t832;
t710 = (pkin(5) * t1080 - (2 * qJD(6))) * t1004 + t810 + t1254;
t709 = t1132 * t857 + t1135 * t739;
t708 = t1132 * t739 - t1135 * t857;
t707 = -pkin(4) * t1148 - t1230 + t1310;
t698 = t1126 * t883 + t1128 * t726;
t697 = t1126 * t726 - t1128 * t883;
t696 = -pkin(8) * t857 - t725;
t691 = -pkin(7) * t868 - t1131 * t802 + t1134 * t804;
t686 = (-t1148 - t1201) * pkin(5) + t1138;
t685 = -pkin(5) * t1201 + t1138 + t1307;
t681 = -qJ(3) * t839 - t1126 * t771 + t1128 * t816;
t680 = -qJ(3) * t836 - t1126 * t767 + t1128 * t811;
t675 = -pkin(7) * t756 - qJ(3) * t1218 - t1131 * t817;
t674 = -pkin(1) * t763 - pkin(2) * t839 - pkin(3) * t953 - pkin(8) * t910 - t1225;
t673 = -pkin(2) * t909 + qJ(3) * t840 + t1126 * t816 + t1128 * t771;
t672 = -pkin(1) * t758 - pkin(2) * t836 - pkin(3) * t949 - pkin(8) * t900 + t1220;
t671 = -pkin(2) * t899 + qJ(3) * t837 + t1126 * t811 + t1128 * t767;
t669 = qJ(6) * t1244 - t1146 - 0.2e1 * t1209;
t668 = -pkin(5) * t1244 + t1150;
t659 = -qJ(3) * t824 + t1128 * t696 + t1243 * t857;
t658 = qJ(6) * t1217 + t1146 + t995;
t657 = -t1125 * t686 - t1148 * t1239 - t1309;
t656 = pkin(5) * t1217 + t1150;
t655 = -pkin(5) * t1236 + t1127 * t685 - t1354;
t650 = qJ(3) * t825 + t1126 * t696 + t1179 * t857;
t649 = -pkin(5) * t1257 - qJ(6) * t876 + t670;
t648 = t1127 * t686 + t1148 * t1169 + t1310;
t647 = -t1353 + t1125 * t685 + (pkin(4) + t1242) * t1259;
t646 = -t1131 * t697 + t1134 * t698;
t645 = t1131 * t698 + t1134 * t697;
t644 = -t1369 + t684;
t643 = -pkin(1) * t738 - pkin(2) * t824 - pkin(3) * t976 - pkin(8) * t859 - t726;
t642 = -t1130 * t712 + t1133 * t742 - t1380;
t637 = t1326 + t683;
t634 = -t1130 * t707 + t1133 * t737 - t1335;
t623 = -qJ(3) * t697 + (-pkin(8) * t1128 + t1243) * t725;
t622 = t995 + (-t1244 - t1256) * qJ(6) + (-t1046 - t1264) * pkin(5) + t1161 + t1326;
t621 = -qJ(6) * t1263 + 0.2e1 * t1210 + (t1244 - t929) * pkin(5) - t1258 + t1369;
t620 = t1132 * t725 + t1135 * t646;
t619 = t1132 * t646 - t1135 * t725;
t618 = -pkin(4) * t810 + qJ(5) * t625;
t617 = -t624 - t1308;
t616 = -t1125 * t669 + t1127 * t668;
t615 = t1125 * t668 + t1127 * t669;
t614 = t1302 + t625;
t613 = -pkin(7) * t763 - t1131 * t673 + t1134 * t681;
t612 = -pkin(7) * t758 - t1131 * t671 + t1134 * t680;
t607 = qJ(3) * t698 + (-pkin(8) * t1126 + t1179) * t725;
t606 = -t1125 * t656 + t1127 * t658 - t1308;
t605 = -pkin(1) * t645 - pkin(2) * t697 + pkin(3) * t883 - pkin(8) * t726;
t604 = t1125 * t658 + t1127 * t656 + t1302;
t603 = -pkin(7) * t738 - t1131 * t650 + t1134 * t659;
t602 = -t1130 * t648 + t1133 * t657 - t1335;
t601 = -t1130 * t647 + t1133 * t655 + t1380;
t600 = -t1130 * t742 - t1133 * t712 - t1393;
t599 = -qJ(5) * t615 + (pkin(5) * t1125 - t1239) * t710;
t598 = t1133 * t625 - t1226;
t597 = t1130 * t625 + t1221;
t596 = -t1126 * t644 + t1128 * t642 - t1387;
t595 = -t1130 * t737 - t1133 * t707 + t1365;
t594 = t1126 * t810 + t1128 * t598;
t593 = t1126 * t598 - t1128 * t810;
t592 = -t1126 * t637 + t1128 * t634 - t1351;
t591 = t1126 * t642 + t1128 * t644 - t1385;
t590 = qJ(5) * t616 + (t1169 - t1242) * t710;
t589 = -t1130 * t615 + t1133 * t616;
t588 = t1130 * t616 + t1133 * t615;
t587 = t1126 * t634 + t1128 * t637 + t1347;
t586 = -t1130 * t657 - t1133 * t648 + t1365;
t585 = t1126 * t710 + t1128 * t589;
t584 = t1126 * t589 - t1128 * t710;
t583 = -t1130 * t614 + t1133 * t617 - t1336;
t582 = -t1130 * t655 - t1133 * t647 + t1393;
t581 = -pkin(3) * t597 - pkin(4) * t624;
t580 = -t1126 * t622 + t1128 * t602 - t1351;
t579 = -pkin(7) * t645 - t1131 * t607 + t1134 * t623;
t578 = -t1126 * t621 + t1128 * t601 + t1387;
t577 = t1126 * t602 + t1128 * t622 + t1347;
t576 = t1126 * t601 + t1128 * t621 + t1385;
t575 = -t1130 * t604 + t1133 * t606 - t1336;
t574 = -t1126 * t670 + t1128 * t583 - t1352;
t573 = -pkin(8) * t597 - qJ(5) * t1221 - t1130 * t618;
t572 = t1126 * t583 + t1128 * t670 + t1348;
t571 = -t1131 * t593 + t1134 * t594;
t570 = t1131 * t594 + t1134 * t593;
t569 = -t1130 * t617 - t1133 * t614 + t1364;
t568 = -pkin(3) * t588 - pkin(4) * t615 - pkin(5) * t669 - qJ(6) * t668;
t567 = -t1126 * t649 + t1128 * t575 - t1352;
t566 = -t1131 * t591 + t1134 * t596 - t1394;
t565 = -t1130 * t606 - t1133 * t604 + t1364;
t564 = t1126 * t575 + t1128 * t649 + t1348;
t563 = -t1131 * t587 + t1134 * t592 - t1367;
t562 = -t1131 * t584 + t1134 * t585;
t561 = t1131 * t585 + t1134 * t584;
t560 = t1132 * t597 + t1135 * t571;
t559 = t1132 * t571 - t1135 * t597;
t558 = -pkin(8) * t588 - t1130 * t590 + t1133 * t599;
t557 = -t1131 * t577 + t1134 * t580 - t1367;
t556 = -t1131 * t576 + t1134 * t578 + t1394;
t555 = t1132 * t588 + t1135 * t562;
t554 = t1132 * t562 - t1135 * t588;
t553 = -t1131 * t572 + t1134 * t574 - t1366;
t552 = -qJ(3) * t593 - t1126 * t581 + t1128 * t573;
t551 = -pkin(1) * t570 - pkin(2) * t593 + pkin(3) * t810 - pkin(8) * t598 + qJ(5) * t1226 - t1133 * t618;
t550 = -t1131 * t564 + t1134 * t567 - t1366;
t549 = -pkin(2) * t597 + qJ(3) * t594 + t1126 * t573 + t1128 * t581;
t548 = -qJ(3) * t584 - t1126 * t568 + t1128 * t558;
t547 = -pkin(1) * t561 - pkin(2) * t584 + pkin(3) * t710 - pkin(8) * t589 - t1130 * t599 - t1133 * t590;
t546 = -pkin(2) * t588 + qJ(3) * t585 + t1126 * t558 + t1128 * t568;
t545 = -pkin(7) * t570 - t1131 * t549 + t1134 * t552;
t544 = -pkin(7) * t561 - t1131 * t546 + t1134 * t548;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1098, -t1099, 0, t1058, 0, 0, 0, 0, 0, 0, t1031, t1032, t1053, t975, 0, 0, 0, 0, 0, 0, t872, t886, t835, t741, 0, 0, 0, 0, 0, 0, t724, t728, t709, t620, 0, 0, 0, 0, 0, 0, t1358, -t628, t1357, t560, 0, 0, 0, 0, 0, 0, t1358, t1357, t628, t555; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1099, -t1098, 0, t1057, 0, 0, 0, 0, 0, 0, t1029, t1030, t1052, t974, 0, 0, 0, 0, 0, 0, t869, t885, t834, t740, 0, 0, 0, 0, 0, 0, t723, t727, t708, t619, 0, 0, 0, 0, 0, 0, t1361, -t626, t1360, t559, 0, 0, 0, 0, 0, 0, t1361, t1360, t626, t554; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1062, t1063, 0, -t1011, 0, 0, 0, 0, 0, 0, t897, t932, t868, t756, 0, 0, 0, 0, 0, 0, t758, t763, t738, t645, 0, 0, 0, 0, 0, 0, t1341, t660, t1345, t570, 0, 0, 0, 0, 0, 0, t1341, t1345, -t660, t561; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1099, 0, -t1098, 0, t1158, -t1081, -t1057, -pkin(6) * t1057, t1071 * t1135 - t1160, t1050 * t1135 - t1101 * t1132, t1066 * t1135 + t1131 * t1182, t1070 * t1135 + t1160, t1064 * t1135 + t1117 * t1132, t1090 * t1135 + t1116, -pkin(6) * t1029 - t1022 * t1132 + t1025 * t1135, -pkin(6) * t1030 - t1023 * t1132 + t1026 * t1135, -pkin(6) * t1052 + t1011 * t1135, -pkin(6) * t974 - (pkin(1) * t1132 - pkin(7) * t1135) * t1011, t1135 * t937 + t1172, -t1045 * t1132 + t1135 * t870, -t1021 * t1132 + t1135 * t934, t1135 * t936 - t1172, t1019 * t1132 + t1135 * t933, t1135 * t960 + t1116, -pkin(6) * t869 - t1132 * t805 + t1135 * t776, -pkin(6) * t885 - t1132 * t818 + t1135 * t803, -pkin(6) * t834 - t1132 * t819 + t1135 * t691, -pkin(6) * t740 - t1132 * t711 + t1135 * t675, -t1132 * t943 + t1135 * t815, -t1132 * t860 + t1135 * t747, -t1132 * t927 + t1135 * t769, -t1132 * t944 + t1135 * t814, -t1132 * t925 + t1135 * t770, -t1132 * t962 + t1135 * t838, -pkin(6) * t723 - t1132 * t672 + t1135 * t612, -pkin(6) * t727 - t1132 * t674 + t1135 * t613, -pkin(6) * t708 - t1132 * t643 + t1135 * t603, -pkin(6) * t619 - t1132 * t605 + t1135 * t579, t1266, t1388, t1356, t1314, t1389, t1315, -t1132 * t595 + t1135 * t563 - t1376, -t1132 * t600 + t1135 * t566 + t1396, -t1132 * t569 + t1135 * t553 - t1377, -pkin(6) * t559 - t1132 * t551 + t1135 * t545, t1266, t1356, -t1388, t1315, -t1389, t1314, -t1132 * t586 + t1135 * t557 - t1376, -t1132 * t565 + t1135 * t550 - t1377, -t1132 * t582 + t1135 * t556 - t1396, -pkin(6) * t554 - t1132 * t547 + t1135 * t544; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1098, 0, t1099, 0, t1081, t1158, t1058, pkin(6) * t1058, t1071 * t1132 + t1159, t1050 * t1132 + t1101 * t1135, t1066 * t1132 - t1131 * t1181, t1070 * t1132 - t1159, t1064 * t1132 - t1117 * t1135, t1090 * t1132 - t1180, pkin(6) * t1031 + t1022 * t1135 + t1025 * t1132, pkin(6) * t1032 + t1023 * t1135 + t1026 * t1132, pkin(6) * t1053 + t1011 * t1132, pkin(6) * t975 - (-pkin(1) * t1135 - pkin(7) * t1132) * t1011, t1132 * t937 - t1171, t1045 * t1135 + t1132 * t870, t1021 * t1135 + t1132 * t934, t1132 * t936 + t1171, -t1019 * t1135 + t1132 * t933, t1132 * t960 - t1180, pkin(6) * t872 + t1132 * t776 + t1135 * t805, pkin(6) * t886 + t1132 * t803 + t1135 * t818, pkin(6) * t835 + t1132 * t691 + t1135 * t819, pkin(6) * t741 + t1132 * t675 + t1135 * t711, t1132 * t815 + t1135 * t943, t1132 * t747 + t1135 * t860, t1132 * t769 + t1135 * t927, t1132 * t814 + t1135 * t944, t1132 * t770 + t1135 * t925, t1132 * t838 + t1135 * t962, pkin(6) * t724 + t1132 * t612 + t1135 * t672, pkin(6) * t728 + t1132 * t613 + t1135 * t674, pkin(6) * t709 + t1132 * t603 + t1135 * t643, pkin(6) * t620 + t1132 * t579 + t1135 * t605, t1272, t1390, t1359, t1316, t1391, t1317, t1132 * t563 + t1135 * t595 + t1378, t1132 * t566 + t1135 * t600 - t1395, t1132 * t553 + t1135 * t569 + t1379, pkin(6) * t560 + t1132 * t545 + t1135 * t551, t1272, t1359, -t1390, t1317, -t1391, t1316, t1132 * t557 + t1135 * t586 + t1378, t1132 * t550 + t1135 * t565 + t1379, t1132 * t556 + t1135 * t582 + t1395, pkin(6) * t555 + t1132 * t544 + t1135 * t547; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1104, t1105, 0, 0, (t1094 + t1170) * t1131, t1093 * t1134 + t1096 * t1131, t1109 * t1134 + t1190, (t1095 - t1178) * t1134, t1111 * t1131 + t1188, 0, pkin(1) * t1096 + pkin(7) * t1065 + t1191, -pkin(1) * t1093 + pkin(7) * t1067 - t1192, pkin(1) * t1100 + pkin(7) * t1097 + t1012, pkin(1) * t1088 + pkin(7) * t1012, t1007 * t1134 + t1008 * t1131, t1131 * t956 + t1134 * t954, t1131 * t989 + t1134 * t986, t1005 * t1134 + t1006 * t1131, t1131 * t988 + t1134 * t985, t1033 * t1134 + t1034 * t1131, -pkin(1) * t1017 + pkin(7) * t898 + t1131 * t923 + t1134 * t887, pkin(1) * t1164 + pkin(7) * t935 + t1131 * t939 + t1134 * t896, -pkin(1) * t1015 + pkin(7) * t871 + t1131 * t804 + t1134 * t802, pkin(1) * t1014 + pkin(7) * t757 - qJ(3) * t1223 + t1134 * t817, t1131 * t895 + t1134 * t893, t1131 * t831 + t1134 * t830, t1131 * t843 + t1134 * t841, t1131 * t894 + t1134 * t892, t1131 * t844 + t1134 * t842, t1131 * t931 + t1134 * t930, -pkin(1) * t899 + pkin(7) * t759 + t1131 * t680 + t1134 * t671, -pkin(1) * t909 + pkin(7) * t764 + t1131 * t681 + t1134 * t673, -pkin(1) * t857 + pkin(7) * t739 + t1131 * t659 + t1134 * t650, -pkin(1) * t725 + pkin(7) * t646 + t1131 * t623 + t1134 * t607, t1249, -t1384, t1343, t1292, -t1383, t1294, t1131 * t592 + t1134 * t587 + t1363, t1131 * t596 + t1134 * t591 - t1392, t1131 * t574 + t1134 * t572 + t1362, -pkin(1) * t597 + pkin(7) * t571 + t1131 * t552 + t1134 * t549, t1249, t1343, t1384, t1294, t1383, t1292, t1131 * t580 + t1134 * t577 + t1363, t1131 * t567 + t1134 * t564 + t1362, t1131 * t578 + t1134 * t576 + t1392, -pkin(1) * t588 + pkin(7) * t562 + t1131 * t548 + t1134 * t546;];
tauB_reg  = t1;
