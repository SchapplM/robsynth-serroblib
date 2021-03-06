% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 16:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRPRR14_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR14_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 16:44:18
% EndTime: 2019-05-07 16:46:10
% DurationCPUTime: 84.56s
% Computational Cost: add. (330857->974), mult. (711115->1509), div. (0->0), fcn. (563555->12), ass. (0->709)
t1147 = sin(qJ(1));
t1152 = cos(qJ(1));
t1141 = sin(pkin(6));
t1142 = cos(pkin(6));
t1146 = sin(qJ(2));
t1151 = cos(qJ(2));
t1145 = sin(qJ(3));
t1150 = cos(qJ(3));
t1310 = qJD(1) * t1142;
t1261 = qJD(2) + t1310;
t1285 = t1141 * t1146;
t1275 = qJD(1) * t1285;
t1096 = t1145 * t1275 - t1150 * t1261;
t1093 = t1096 ^ 2;
t1284 = t1141 * t1151;
t1126 = qJD(1) * t1284 - qJD(3);
t1122 = t1126 ^ 2;
t1042 = -t1122 - t1093;
t1277 = qJDD(1) * t1141;
t1107 = -qJD(2) * t1275 + t1151 * t1277;
t1101 = -qJDD(3) + t1107;
t1098 = t1145 * t1261 + t1150 * t1275;
t1296 = t1096 * t1098;
t1165 = t1101 + t1296;
t1352 = t1145 * t1165;
t946 = t1042 * t1150 + t1352;
t1276 = t1146 * qJDD(1);
t1309 = qJD(1) * t1151;
t1106 = (qJD(2) * t1309 + t1276) * t1141;
t1133 = qJDD(1) * t1142 + qJDD(2);
t1258 = t1106 * t1145 - t1150 * t1133;
t1190 = qJD(3) * t1098 + t1258;
t1290 = t1126 * t1098;
t992 = t1290 - t1190;
t1223 = t1146 * t946 + t1151 * t992;
t1351 = t1150 * t1165;
t944 = t1042 * t1145 - t1351;
t827 = -t1141 * t944 + t1142 * t1223;
t883 = -t1146 * t992 + t1151 * t946;
t757 = t1147 * t883 + t1152 * t827;
t1431 = pkin(7) * t757;
t760 = t1147 * t827 - t1152 * t883;
t1430 = pkin(7) * t760;
t1337 = t1098 ^ 2;
t1049 = t1337 + t1122;
t1166 = t1101 - t1296;
t1349 = t1166 * t1150;
t955 = t1049 * t1145 + t1349;
t1203 = t1150 * t1106 + t1145 * t1133;
t998 = (qJD(3) - t1126) * t1096 - t1203;
t1221 = t1146 * t955 + t1151 * t998;
t1350 = t1166 * t1145;
t954 = t1049 * t1150 - t1350;
t838 = t1141 * t954 + t1142 * t1221;
t886 = -t1146 * t998 + t1151 * t955;
t762 = t1147 * t886 + t1152 * t838;
t1429 = pkin(7) * t762;
t765 = t1147 * t838 - t1152 * t886;
t1428 = pkin(7) * t765;
t836 = t1141 * t1221 - t1142 * t954;
t1427 = pkin(8) * (t1141 * t836 + t1142 * t838);
t825 = t1141 * t1223 + t1142 * t944;
t1426 = pkin(8) * (t1141 * t825 + t1142 * t827);
t1425 = pkin(1) * t825;
t1424 = pkin(1) * t827;
t1423 = pkin(1) * t836;
t1422 = pkin(1) * t838;
t1066 = t1122 - t1093;
t967 = t1066 * t1150 - t1350;
t990 = t1290 + t1190;
t1216 = t1146 * t967 - t1151 * t990;
t963 = t1066 * t1145 + t1349;
t849 = -t1141 * t963 + t1142 * t1216;
t899 = t1146 * t990 + t1151 * t967;
t1421 = t1147 * t849 - t1152 * t899;
t1420 = t1147 * t899 + t1152 * t849;
t1417 = pkin(8) * t883;
t1416 = pkin(8) * t886;
t1409 = t1141 * t1216 + t1142 * t963;
t1408 = pkin(2) * t944;
t1407 = pkin(2) * t954;
t1406 = pkin(9) * t944;
t1405 = pkin(9) * t946;
t1404 = pkin(9) * t954;
t1403 = pkin(9) * t955;
t1065 = -t1337 + t1122;
t961 = t1065 * t1150 - t1352;
t1400 = t1141 * t961;
t1396 = t1142 * t961;
t964 = t1065 * t1145 + t1351;
t1394 = t1146 * t964;
t1393 = t1151 * t964;
t1336 = -2 * qJD(4);
t1018 = -t1337 - t1093;
t1392 = pkin(2) * t1018;
t1391 = t1145 * t990;
t1390 = t1145 * t992;
t1386 = t1150 * t990;
t1312 = t1150 * t992;
t1382 = t1018 * t1146;
t1381 = t1018 * t1151;
t1074 = t1096 * t1126;
t1169 = qJD(3) * t1096 - t1203;
t997 = t1074 - t1169;
t1333 = pkin(2) * t1151;
t1256 = -pkin(9) * t1146 - t1333;
t1311 = qJD(1) * t1141;
t1105 = t1256 * t1311;
t1257 = t1261 ^ 2;
t1129 = g(1) * t1152 + g(2) * t1147;
t1153 = qJD(1) ^ 2;
t1102 = -pkin(1) * t1153 + pkin(8) * t1277 - t1129;
t1128 = t1147 * g(1) - t1152 * g(2);
t1330 = pkin(8) * t1141;
t1163 = qJDD(1) * pkin(1) + t1153 * t1330 + t1128;
t1161 = t1142 * t1163;
t1259 = t1146 * t1102 - t1151 * t1161;
t974 = (qJD(1) * t1105 * t1146 + g(3) * t1151) * t1141 - t1133 * pkin(2) - t1257 * pkin(9) + t1259;
t1155 = t1190 * pkin(3) - qJ(4) * t997 + t974;
t1380 = t1098 * t1336 + t1155;
t1288 = t1126 * t1150;
t1157 = -t1096 * t1288 + t1145 * t1190;
t1270 = t1146 * t1296;
t1343 = t1151 * t1157 - t1270;
t1289 = t1126 * t1145;
t1253 = -t1096 * t1289 - t1150 * t1190;
t1269 = t1151 * t1296;
t1345 = t1146 * t1157 + t1269;
t1361 = -t1141 * t1253 + t1142 * t1345;
t1379 = -t1147 * t1361 + t1152 * t1343;
t1252 = t1098 * t1289 - t1150 * t1169;
t1341 = t1151 * t1252 + t1270;
t1254 = -t1098 * t1288 - t1145 * t1169;
t1344 = t1146 * t1252 - t1269;
t1362 = -t1141 * t1254 + t1142 * t1344;
t1378 = -t1147 * t1362 + t1152 * t1341;
t1377 = t1147 * t1343 + t1152 * t1361;
t1376 = t1147 * t1341 + t1152 * t1362;
t1143 = sin(qJ(6));
t1036 = qJDD(5) - t1169;
t1034 = qJDD(6) + t1036;
t1144 = sin(qJ(5));
t1149 = cos(qJ(5));
t1058 = -t1149 * t1096 - t1126 * t1144;
t1060 = t1096 * t1144 - t1126 * t1149;
t1148 = cos(qJ(6));
t1002 = t1148 * t1058 + t1060 * t1143;
t1004 = -t1058 * t1143 + t1060 * t1148;
t926 = t1004 * t1002;
t1358 = t1034 - t926;
t1375 = t1143 * t1358;
t1010 = t1060 * t1058;
t1348 = -t1010 + t1036;
t1374 = t1144 * t1348;
t1347 = -t1337 + t1093;
t1373 = t1146 * t1347;
t1370 = t1148 * t1358;
t1369 = t1149 * t1348;
t1368 = t1151 * t1347;
t1365 = -t1169 - t1074;
t1364 = t1141 * t1344 + t1142 * t1254;
t1363 = t1141 * t1345 + t1142 * t1253;
t1188 = (t1096 * t1145 + t1098 * t1150) * t1126;
t1189 = (t1096 * t1150 - t1098 * t1145) * t1126;
t1339 = -t1141 * t1188 + (t1101 * t1151 + t1146 * t1189) * t1142;
t1342 = -t1101 * t1146 + t1151 * t1189;
t1360 = t1147 * t1342 + t1152 * t1339;
t1359 = -t1147 * t1339 + t1152 * t1342;
t1091 = qJD(5) + t1098;
t1301 = t1058 * t1091;
t950 = -t1058 * qJD(5) - t1149 * t1101 + t1144 * t1190;
t1201 = t950 + t1301;
t1138 = t1141 ^ 2;
t1286 = t1138 * t1153;
t1211 = t1261 * qJD(1);
t1357 = t1138 * (-t1142 * t1153 + t1211);
t1046 = pkin(3) * t1096 - qJ(4) * t1098;
t1156 = -g(3) * t1285 + t1146 * t1161;
t975 = t1133 * pkin(9) - t1257 * pkin(2) + (t1105 * t1311 + t1102) * t1151 + t1156;
t1191 = t1146 * t1211;
t1192 = t1151 * t1211;
t1329 = t1142 * g(3);
t976 = -t1107 * pkin(2) - t1106 * pkin(9) - t1329 + (pkin(2) * t1191 - pkin(9) * t1192 - t1163) * t1141;
t1328 = -t1145 * t975 + t1150 * t976;
t1198 = t1101 * pkin(3) - t1122 * qJ(4) + t1098 * t1046 + qJDD(4) - t1328;
t779 = pkin(4) * t1365 + t1165 * pkin(10) + t1198;
t1063 = pkin(4) * t1098 + pkin(10) * t1126;
t1264 = -pkin(3) * t1126 + t1336;
t799 = t1258 * pkin(10) - t1093 * pkin(4) + (pkin(10) * qJD(3) - t1063 + t1264) * t1098 + t1155;
t1262 = -t1144 * t799 + t1149 * t779;
t723 = t1144 * t779 + t1149 * t799;
t656 = -t1144 * t1262 + t1149 * t723;
t655 = t1144 * t723 + t1149 * t1262;
t1340 = t1101 * t1284 + t1142 * t1188 + t1189 * t1285;
t1000 = t1002 ^ 2;
t1001 = t1004 ^ 2;
t1057 = t1058 ^ 2;
t1338 = t1060 ^ 2;
t1086 = qJD(6) + t1091;
t1083 = t1086 ^ 2;
t1088 = t1091 ^ 2;
t1335 = pkin(3) + pkin(10);
t1334 = pkin(2) * t1146;
t1332 = pkin(3) * t1145;
t1331 = pkin(3) * t1150;
t1327 = t1101 * qJ(4);
t1020 = pkin(5) * t1091 - pkin(11) * t1060;
t895 = t1145 * t976 + t1150 * t975;
t1164 = t1122 * pkin(3) + t1096 * t1046 - t895;
t798 = -t1327 - t1190 * pkin(4) - t1093 * pkin(10) + (t1336 - t1063) * t1126 - t1164;
t1260 = t1144 * t1101 + t1149 * t1190;
t949 = -qJD(5) * t1060 + t1260;
t740 = -t949 * pkin(5) - t1057 * pkin(11) + t1060 * t1020 + t798;
t1326 = t1143 * t740;
t876 = t1034 + t926;
t1325 = t1143 * t876;
t688 = pkin(5) * t1348 - pkin(11) * t1201 + t1262;
t692 = -pkin(5) * t1057 + pkin(11) * t949 - t1020 * t1091 + t723;
t633 = t1143 * t692 - t1148 * t688;
t634 = t1143 * t688 + t1148 * t692;
t599 = t1143 * t634 - t1148 * t633;
t1324 = t1144 * t599;
t1322 = t1144 * t798;
t937 = t1010 + t1036;
t1321 = t1144 * t937;
t1320 = t1145 * t974;
t1319 = t1148 * t740;
t1318 = t1148 * t876;
t1317 = t1149 * t599;
t1315 = t1149 * t798;
t1314 = t1149 * t937;
t1313 = t1150 * t974;
t1300 = t1086 * t1143;
t1299 = t1086 * t1148;
t1298 = t1091 * t1144;
t1297 = t1091 * t1149;
t1125 = t1151 * t1146 * t1286;
t1103 = -t1125 + t1133;
t1294 = t1103 * t1146;
t1293 = t1103 * t1151;
t1104 = t1125 + t1133;
t1292 = t1104 * t1146;
t1291 = t1104 * t1151;
t1287 = t1133 * t1141;
t1078 = t1141 * t1163 + t1329;
t1282 = t1146 * t1078;
t1281 = t1151 * t1078;
t1279 = qJD(6) + t1086;
t1139 = t1146 ^ 2;
t1140 = t1151 ^ 2;
t1278 = t1139 + t1140;
t1274 = t1145 * t926;
t1273 = t1150 * t926;
t1272 = t1145 * t1010;
t1271 = t1150 * t1010;
t1268 = t1139 * t1286;
t1267 = t1140 * t1286;
t1266 = -t1088 - t1338;
t1265 = qJ(4) * t1145 + pkin(2);
t600 = t1143 * t633 + t1148 * t634;
t1263 = t1143 * t950 - t1148 * t949;
t809 = -t1145 * t1328 + t1150 * t895;
t1080 = -t1128 * t1147 - t1152 * t1129;
t1124 = qJDD(1) * t1152 - t1147 * t1153;
t1255 = -pkin(7) * t1124 - g(3) * t1147;
t1089 = -t1268 - t1257;
t1045 = -t1089 * t1146 - t1293;
t1251 = pkin(8) * t1045 - t1282;
t1111 = -t1257 - t1267;
t1054 = t1111 * t1151 - t1292;
t1250 = pkin(8) * t1054 + t1281;
t1172 = (-qJD(6) + t1086) * t1004 - t1263;
t1246 = t1143 * t949 + t1148 * t950;
t834 = -qJD(6) * t1002 + t1246;
t959 = t1086 * t1002;
t794 = t834 + t959;
t1249 = t1143 * t1172 - t1148 * t794;
t938 = -t1001 - t1083;
t1248 = -t1148 * t938 + t1325;
t900 = -t1083 - t1000;
t1247 = t1143 * t900 + t1370;
t808 = t1145 * t895 + t1150 * t1328;
t574 = t1144 * t600 + t1317;
t568 = t1145 * t574 + t1150 * t740;
t575 = t1149 * t600 - t1324;
t1245 = t1146 * t568 - t1151 * t575;
t639 = t1145 * t655 + t1150 * t798;
t1244 = t1146 * t639 - t1151 * t656;
t726 = t1143 * t794 + t1148 * t1172;
t659 = t1144 * t726 + t1149 * t1249;
t863 = -t1000 - t1001;
t644 = t1145 * t659 + t1150 * t863;
t661 = -t1144 * t1249 + t1149 * t726;
t1243 = t1146 * t644 - t1151 * t661;
t790 = t1004 * t1279 + t1263;
t793 = t834 - t959;
t724 = -t1143 * t790 + t1148 * t793;
t725 = -t1143 * t793 - t1148 * t790;
t658 = -t1144 * t725 - t1149 * t724;
t925 = t1001 - t1000;
t648 = -t1145 * t658 + t1150 * t925;
t660 = t1144 * t724 - t1149 * t725;
t1242 = t1146 * t648 + t1151 * t660;
t803 = t1148 * t900 - t1375;
t728 = t1144 * t803 + t1149 * t1247;
t682 = t1145 * t728 + t1150 * t790;
t729 = -t1144 * t1247 + t1149 * t803;
t1241 = t1146 * t682 - t1151 * t729;
t823 = -t1143 * t938 - t1318;
t734 = t1144 * t823 - t1149 * t1248;
t795 = -t1002 * t1279 + t1246;
t691 = t1145 * t734 + t1150 * t795;
t735 = t1144 * t1248 + t1149 * t823;
t1240 = t1146 * t691 - t1151 * t735;
t957 = t1000 - t1083;
t830 = t1143 * t957 + t1318;
t832 = t1148 * t957 - t1325;
t743 = -t1144 * t832 - t1149 * t830;
t695 = -t1145 * t743 + t1150 * t1172;
t745 = t1144 * t830 - t1149 * t832;
t1239 = t1146 * t695 + t1151 * t745;
t958 = -t1001 + t1083;
t829 = t1148 * t958 + t1375;
t831 = -t1143 * t958 + t1370;
t742 = -t1144 * t831 - t1149 * t829;
t696 = -t1145 * t742 + t1150 * t794;
t744 = t1144 * t829 - t1149 * t831;
t1238 = t1146 * t696 + t1151 * t744;
t833 = -qJD(6) * t1004 - t1263;
t782 = t1002 * t1300 + t1148 * t833;
t783 = t1002 * t1299 - t1143 * t833;
t717 = -t1144 * t783 - t1149 * t782;
t703 = -t1145 * t717 - t1273;
t719 = t1144 * t782 - t1149 * t783;
t1237 = t1146 * t703 + t1151 * t719;
t784 = t1004 * t1299 + t1143 * t834;
t785 = -t1004 * t1300 + t1148 * t834;
t718 = -t1144 * t785 - t1149 * t784;
t704 = -t1145 * t718 + t1273;
t720 = t1144 * t784 - t1149 * t785;
t1236 = t1146 * t704 + t1151 * t720;
t1159 = 0.2e1 * qJD(4) * t1126 + t1164;
t844 = -t1159 - t1327;
t751 = t1145 * t1198 + t1150 * t844;
t852 = t1098 * t1264 + t1155;
t1235 = t1146 * t751 - t1151 * t852;
t879 = (-t1002 * t1143 - t1004 * t1148) * t1086;
t880 = (-t1002 * t1148 + t1004 * t1143) * t1086;
t780 = -t1144 * t880 - t1149 * t879;
t769 = t1034 * t1150 - t1145 * t780;
t781 = t1144 * t879 - t1149 * t880;
t1234 = t1146 * t769 + t1151 * t781;
t1170 = (-qJD(5) + t1091) * t1060 + t1260;
t818 = t1144 * t1170 - t1149 * t1201;
t951 = -t1057 - t1338;
t773 = t1145 * t818 + t1150 * t951;
t820 = t1144 * t1201 + t1149 * t1170;
t1233 = t1146 * t773 - t1151 * t820;
t1009 = -t1057 + t1338;
t1202 = t950 - t1301;
t908 = (qJD(5) + t1091) * t1060 - t1260;
t817 = t1144 * t908 - t1149 * t1202;
t777 = t1009 * t1150 - t1145 * t817;
t819 = t1144 * t1202 + t1149 * t908;
t1232 = t1146 * t777 + t1151 * t819;
t968 = -t1088 - t1057;
t868 = t1144 * t968 + t1369;
t801 = t1145 * t868 + t1150 * t908;
t869 = t1149 * t968 - t1374;
t1231 = t1146 * t801 - t1151 * t869;
t873 = t1149 * t1266 - t1321;
t807 = t1145 * t873 + t1150 * t1202;
t874 = -t1144 * t1266 - t1314;
t1230 = t1146 * t807 - t1151 * t874;
t1229 = t1146 * t809 - t1151 * t974;
t1021 = t1057 - t1088;
t890 = -t1021 * t1144 - t1314;
t812 = -t1145 * t890 + t1150 * t1170;
t892 = -t1021 * t1149 + t1321;
t1228 = t1146 * t812 + t1151 * t892;
t1022 = t1088 - t1338;
t889 = -t1022 * t1149 - t1374;
t813 = -t1145 * t889 + t1150 * t1201;
t891 = t1022 * t1144 - t1369;
t1227 = t1146 * t813 + t1151 * t891;
t903 = -t1058 * t1298 - t1149 * t949;
t856 = -t1145 * t903 - t1271;
t904 = -t1058 * t1297 + t1144 * t949;
t1226 = t1146 * t856 + t1151 * t904;
t905 = -t1060 * t1297 - t1144 * t950;
t857 = -t1145 * t905 + t1271;
t906 = t1060 * t1298 - t1149 * t950;
t1225 = t1146 * t857 + t1151 * t906;
t941 = (t1058 * t1144 + t1060 * t1149) * t1091;
t902 = t1036 * t1150 - t1145 * t941;
t942 = (t1058 * t1149 - t1060 * t1144) * t1091;
t1224 = t1146 * t902 + t1151 * t942;
t995 = (-qJD(3) - t1126) * t1096 + t1203;
t1219 = t1151 * t995 + t1394;
t1218 = -t1151 * t1365 - t1394;
t917 = t1145 * t995 - t1386;
t1215 = t1146 * t917 - t1381;
t920 = t1145 * t1365 - t1386;
t1214 = t1146 * t920 - t1381;
t918 = -t1145 * t997 + t1312;
t1213 = t1146 * t918 + t1368;
t919 = t1145 * t998 + t1312;
t1212 = t1146 * t919 + t1368;
t1043 = g(3) * t1284 + t1259;
t1044 = t1151 * t1102 + t1156;
t1210 = -t1151 * t1043 + t1146 * t1044;
t952 = t1043 * t1146 + t1044 * t1151;
t1115 = t1141 * t1192;
t1070 = t1115 + t1106;
t1114 = t1141 * t1191;
t1073 = t1107 - t1114;
t1209 = t1070 * t1151 + t1073 * t1146;
t1071 = -t1115 + t1106;
t1072 = t1107 + t1114;
t1208 = -t1071 * t1151 + t1072 * t1146;
t1207 = t1089 * t1151 - t1294;
t1110 = -t1257 + t1267;
t1206 = t1110 * t1146 + t1293;
t1109 = t1257 - t1268;
t1205 = t1109 * t1151 + t1292;
t1204 = t1111 * t1146 + t1291;
t1079 = t1128 * t1152 - t1129 * t1147;
t1193 = t1141 * t1211;
t594 = -pkin(5) * t740 + pkin(11) * t600;
t550 = pkin(4) * t740 + pkin(11) * t1324 - t1149 * t594 - t1335 * t575;
t554 = pkin(4) * t574 + pkin(5) * t599 - qJ(4) * t575;
t567 = t1145 * t740 - t1150 * t574;
t540 = -pkin(9) * t567 - t1145 * t550 + t1150 * t554;
t544 = -pkin(2) * t567 + pkin(11) * t1317 - qJ(4) * t740 + t1144 * t594 + t1335 * t574;
t557 = t1146 * t575 + t1151 * t568;
t1187 = pkin(8) * t557 + t1146 * t540 + t1151 * t544;
t587 = -pkin(5) * t863 + pkin(11) * t726 + t600;
t592 = -pkin(11) * t1249 - t599;
t559 = pkin(4) * t863 - t1144 * t592 - t1149 * t587 - t1335 * t661;
t603 = pkin(4) * t659 + pkin(5) * t1249 - qJ(4) * t661;
t643 = t1145 * t863 - t1150 * t659;
t556 = -pkin(9) * t643 - t1145 * t559 + t1150 * t603;
t558 = -pkin(2) * t643 - qJ(4) * t863 + t1144 * t587 - t1149 * t592 + t1335 * t659;
t609 = t1146 * t661 + t1151 * t644;
t1186 = pkin(8) * t609 + t1146 * t556 + t1151 * t558;
t665 = -pkin(5) * t790 + pkin(11) * t803 - t1319;
t697 = -pkin(11) * t1247 + t1326;
t597 = pkin(4) * t790 - t1144 * t697 - t1149 * t665 - t1335 * t729;
t602 = pkin(4) * t728 + pkin(5) * t1247 - qJ(4) * t729 - t633;
t681 = t1145 * t790 - t1150 * t728;
t563 = -pkin(9) * t681 - t1145 * t597 + t1150 * t602;
t584 = -pkin(2) * t681 - qJ(4) * t790 + t1144 * t665 - t1149 * t697 + t1335 * t728;
t641 = t1146 * t729 + t1151 * t682;
t1185 = pkin(8) * t641 + t1146 * t563 + t1151 * t584;
t670 = -pkin(5) * t795 + pkin(11) * t823 + t1326;
t708 = pkin(11) * t1248 + t1319;
t601 = pkin(4) * t795 - t1144 * t708 - t1149 * t670 - t1335 * t735;
t605 = pkin(4) * t734 - pkin(5) * t1248 - qJ(4) * t735 - t634;
t690 = t1145 * t795 - t1150 * t734;
t569 = -pkin(9) * t690 - t1145 * t601 + t1150 * t605;
t589 = -pkin(2) * t690 - qJ(4) * t795 + t1144 * t670 - t1149 * t708 + t1335 * t734;
t649 = t1146 * t735 + t1151 * t691;
t1184 = pkin(8) * t649 + t1146 * t569 + t1151 * t589;
t607 = pkin(4) * t798 - t1335 * t656;
t611 = pkin(4) * t655 - qJ(4) * t656;
t638 = t1145 * t798 - t1150 * t655;
t570 = -pkin(9) * t638 - t1145 * t607 + t1150 * t611;
t583 = -pkin(2) * t638 - qJ(4) * t798 + t1335 * t655;
t606 = t1146 * t656 + t1151 * t639;
t1183 = pkin(8) * t606 + t1146 * t570 + t1151 * t583;
t628 = pkin(4) * t951 - t1335 * t820 - t656;
t733 = pkin(4) * t818 - qJ(4) * t820;
t772 = t1145 * t951 - t1150 * t818;
t608 = -pkin(9) * t772 - t1145 * t628 + t1150 * t733;
t615 = -pkin(2) * t772 - qJ(4) * t951 + t1335 * t818 + t655;
t727 = t1146 * t820 + t1151 * t773;
t1182 = pkin(8) * t727 + t1146 * t608 + t1151 * t615;
t678 = pkin(4) * t868 - qJ(4) * t869 + t1262;
t700 = pkin(4) * t908 - t1335 * t869 + t1315;
t800 = t1145 * t908 - t1150 * t868;
t624 = -pkin(9) * t800 - t1145 * t700 + t1150 * t678;
t668 = -pkin(2) * t800 - qJ(4) * t908 + t1335 * t868 - t1322;
t747 = t1146 * t869 + t1151 * t801;
t1181 = pkin(8) * t747 + t1146 * t624 + t1151 * t668;
t680 = pkin(4) * t873 - qJ(4) * t874 - t723;
t705 = pkin(4) * t1202 - t1335 * t874 - t1322;
t806 = t1145 * t1202 - t1150 * t873;
t625 = -pkin(9) * t806 - t1145 * t705 + t1150 * t680;
t672 = -pkin(2) * t806 - qJ(4) * t1202 + t1335 * t873 - t1315;
t749 = t1146 * t874 + t1151 * t807;
t1180 = pkin(8) * t749 + t1146 * t625 + t1151 * t672;
t750 = t1145 * t844 - t1150 * t1198;
t685 = -pkin(2) * t750 + pkin(3) * t1198 - qJ(4) * t844;
t686 = -pkin(9) * t750 + (-qJ(4) * t1150 + t1332) * t852;
t714 = t1146 * t852 + t1151 * t751;
t1179 = pkin(8) * t714 + t1146 * t686 + t1151 * t685;
t816 = -pkin(3) * t1018 + t844;
t822 = -qJ(4) * t1018 + t1198;
t913 = -t1150 * t995 - t1391;
t713 = -pkin(9) * t913 - t1145 * t816 + t1150 * t822;
t824 = -pkin(2) * t913 + pkin(3) * t995 + qJ(4) * t990;
t861 = t1151 * t917 + t1382;
t1178 = pkin(8) * t861 + t1146 * t713 + t1151 * t824;
t815 = (-t992 - t1290) * pkin(3) + t1380;
t752 = -qJ(4) * t1312 - t1145 * t815 + t1406;
t761 = -pkin(3) * t1165 + qJ(4) * t1042 - t1198 + t1408;
t1177 = t1146 * t752 + t1151 * t761 - t1417;
t814 = pkin(3) * t1290 - qJ(4) * t998 - t1380;
t755 = t1150 * t814 + t1332 * t998 - t1404;
t766 = -t1407 - pkin(3) * t1049 + (t1166 + t1101) * qJ(4) + t1159;
t1176 = t1146 * t755 + t1151 * t766 - t1416;
t840 = -t1328 - t1408;
t878 = t1320 - t1406;
t1175 = t1146 * t878 + t1151 * t840 + t1417;
t843 = t895 + t1407;
t885 = t1313 + t1404;
t1174 = t1146 * t885 + t1151 * t843 + t1416;
t1011 = t1071 * t1146 + t1072 * t1151;
t1173 = pkin(8) * t1011 + t952;
t916 = -t1150 * t1365 - t1391;
t756 = -pkin(9) * t916 - t808;
t862 = t1151 * t920 + t1382;
t1171 = pkin(8) * t862 + t1146 * t756 - t1333 * t916;
t770 = t1146 * t974 + t1151 * t809;
t1167 = pkin(8) * t770 + t1256 * t808;
t1160 = t1141 * t1286 + t1142 * t1193;
t1123 = qJDD(1) * t1147 + t1152 * t1153;
t1113 = t1278 * t1286;
t1112 = (t1139 - t1140) * t1286;
t1108 = -pkin(7) * t1123 + g(3) * t1152;
t1077 = t1261 * t1278 * t1311;
t1069 = (t1276 + (0.2e1 * qJD(2) + t1310) * t1309) * t1141;
t1062 = t1151 * t1106 - t1139 * t1193;
t1061 = -t1146 * t1107 - t1140 * t1193;
t1053 = t1110 * t1151 - t1294;
t1052 = -t1109 * t1146 + t1291;
t1041 = (t1142 * t1106 + t1151 * t1160) * t1146;
t1040 = (t1142 * t1107 - t1146 * t1160) * t1151;
t1012 = -t1070 * t1146 + t1073 * t1151;
t1008 = t1141 * t1073 + t1142 * t1204;
t1007 = -t1141 * t1072 + t1142 * t1206;
t1006 = -t1141 * t1071 + t1142 * t1205;
t1005 = -t1142 * t1073 + t1141 * t1204;
t988 = -t1141 * t1069 + t1142 * t1207;
t987 = t1142 * t1069 + t1141 * t1207;
t973 = -t1141 * t1112 + t1142 * t1209;
t972 = t1141 * t1113 + t1142 * t1208;
t971 = -t1142 * t1113 + t1141 * t1208;
t940 = -t1008 * t1147 + t1054 * t1152;
t939 = t1008 * t1152 + t1054 * t1147;
t930 = t1045 * t1152 - t1147 * t988;
t929 = t1045 * t1147 + t1152 * t988;
t928 = t1141 * t1078 + t1142 * t1210;
t927 = -t1142 * t1078 + t1141 * t1210;
t922 = t1011 * t1152 - t1147 * t972;
t921 = t1011 * t1147 + t1152 * t972;
t915 = -t1150 * t998 + t1390;
t914 = t1150 * t997 + t1390;
t901 = t1036 * t1145 + t1150 * t941;
t897 = t1146 * t1365 - t1393;
t896 = -t1146 * t995 + t1393;
t888 = -t1282 + (-t1005 * t1141 - t1008 * t1142) * pkin(8);
t882 = -t1281 + (-t1141 * t987 - t1142 * t988) * pkin(8);
t881 = -pkin(1) * t1005 + t1141 * t1043 + t1142 * t1250;
t872 = t1151 * t919 - t1373;
t871 = t1151 * t918 - t1373;
t870 = -pkin(1) * t987 + t1141 * t1044 + t1142 * t1251;
t860 = pkin(8) * t1142 * t952 - pkin(1) * t927;
t859 = -t1147 * t928 + t1152 * t952;
t858 = t1147 * t952 + t1152 * t928;
t855 = t1150 * t905 + t1272;
t854 = t1150 * t903 - t1272;
t853 = -pkin(1) * t971 + t1142 * t1173;
t850 = pkin(2) * t998 + t1320 + t1403;
t847 = t1142 * t1218 - t1400;
t846 = t1142 * t1219 + t1400;
t845 = (-t1141 * t927 - t1142 * t928) * pkin(8);
t842 = pkin(2) * t992 - t1313 + t1405;
t841 = (-t1141 * t971 - t1142 * t972) * pkin(8) - t1210;
t835 = -t1146 * t942 + t1151 * t902;
t811 = t1145 * t1201 + t1150 * t889;
t810 = t1145 * t1170 + t1150 * t890;
t805 = -t1141 * t915 + t1142 * t1212;
t804 = -t1141 * t914 + t1142 * t1213;
t789 = -t1141 * t916 + t1142 * t1214;
t788 = -t1141 * t913 + t1142 * t1215;
t787 = t1141 * t1214 + t1142 * t916;
t786 = t1141 * t1215 + t1142 * t913;
t776 = t1009 * t1145 + t1150 * t817;
t775 = -t1146 * t906 + t1151 * t857;
t774 = -t1146 * t904 + t1151 * t856;
t771 = -pkin(2) * t974 + pkin(9) * t809;
t768 = t1034 * t1145 + t1150 * t780;
t767 = -t1141 * t901 + t1142 * t1224;
t754 = -t1146 * t891 + t1151 * t813;
t753 = -t1146 * t892 + t1151 * t812;
t748 = pkin(9) * t920 - t1392 + t809;
t746 = -t1403 + t1145 * t814 - (pkin(2) + t1331) * t998;
t741 = t1150 * t815 - t1265 * t992 - t1405;
t739 = -t1147 * t789 + t1152 * t862;
t738 = -t1147 * t788 + t1152 * t861;
t737 = t1147 * t862 + t1152 * t789;
t736 = t1147 * t861 + t1152 * t788;
t732 = -t1141 * t855 + t1142 * t1225;
t731 = -t1141 * t854 + t1142 * t1226;
t730 = -t1146 * t819 + t1151 * t777;
t716 = -t1141 * t808 + t1142 * t1229;
t715 = t1141 * t1229 + t1142 * t808;
t712 = -t1146 * t781 + t1151 * t769;
t711 = pkin(9) * t917 + t1145 * t822 + t1150 * t816 - t1392;
t710 = -t1141 * t811 + t1142 * t1227;
t709 = -t1141 * t810 + t1142 * t1228;
t707 = -t1141 * t806 + t1142 * t1230;
t706 = t1141 * t1230 + t1142 * t806;
t702 = t1150 * t718 + t1274;
t701 = t1150 * t717 - t1274;
t699 = -t1141 * t800 + t1142 * t1231;
t698 = t1141 * t1231 + t1142 * t800;
t694 = t1145 * t794 + t1150 * t742;
t693 = t1145 * t1172 + t1150 * t743;
t689 = -t1146 * t843 + t1151 * t885 - t1427;
t684 = -t1146 * t840 + t1151 * t878 - t1426;
t683 = -t1141 * t776 + t1142 * t1232;
t679 = -t1141 * t850 + t1142 * t1174 - t1423;
t677 = -t1141 * t772 + t1142 * t1233;
t676 = t1141 * t1233 + t1142 * t772;
t675 = -t1147 * t716 + t1152 * t770;
t674 = t1147 * t770 + t1152 * t716;
t673 = -t1141 * t842 + t1142 * t1175 - t1425;
t671 = -t1141 * t768 + t1142 * t1234;
t669 = pkin(9) * t751 + (-t1265 - t1331) * t852;
t667 = -t1141 * t750 + t1142 * t1235;
t666 = t1141 * t1235 + t1142 * t750;
t664 = t916 * t1334 + t1151 * t756 + (-t1141 * t787 - t1142 * t789) * pkin(8);
t663 = -t1147 * t707 + t1152 * t749;
t662 = t1147 * t749 + t1152 * t707;
t657 = -t1146 * t766 + t1151 * t755 + t1427;
t654 = -t1147 * t699 + t1152 * t747;
t653 = t1147 * t747 + t1152 * t699;
t652 = -t1146 * t744 + t1151 * t696;
t651 = -t1146 * t745 + t1151 * t695;
t650 = -t1146 * t761 + t1151 * t752 + t1426;
t647 = t1145 * t925 + t1150 * t658;
t646 = -t1146 * t720 + t1151 * t704;
t645 = -t1146 * t719 + t1151 * t703;
t642 = -pkin(1) * t787 - t1141 * t748 + t1142 * t1171;
t640 = -t1146 * t824 + t1151 * t713 + (-t1141 * t786 - t1142 * t788) * pkin(8);
t637 = -t1147 * t677 + t1152 * t727;
t636 = t1147 * t727 + t1152 * t677;
t635 = -t1141 * t746 + t1142 * t1176 + t1423;
t631 = -t1141 * t741 + t1142 * t1177 + t1425;
t630 = -t1147 * t667 + t1152 * t714;
t629 = t1147 * t714 + t1152 * t667;
t627 = (-pkin(9) * t1151 + t1334) * t808 + (-t1141 * t715 - t1142 * t716) * pkin(8);
t626 = -pkin(1) * t715 - t1141 * t771 + t1142 * t1167;
t623 = -t1141 * t694 + t1142 * t1238;
t622 = -t1141 * t693 + t1142 * t1239;
t621 = -pkin(1) * t786 - t1141 * t711 + t1142 * t1178;
t620 = -pkin(2) * t874 + pkin(9) * t807 + t1145 * t680 + t1150 * t705;
t619 = -t1141 * t702 + t1142 * t1236;
t618 = -t1141 * t701 + t1142 * t1237;
t617 = -t1141 * t690 + t1142 * t1240;
t616 = t1141 * t1240 + t1142 * t690;
t614 = -pkin(2) * t869 + pkin(9) * t801 + t1145 * t678 + t1150 * t700;
t613 = -t1141 * t681 + t1142 * t1241;
t612 = t1141 * t1241 + t1142 * t681;
t610 = -t1146 * t660 + t1151 * t648;
t604 = -pkin(2) * t820 + pkin(9) * t773 + t1145 * t733 + t1150 * t628;
t596 = -t1147 * t617 + t1152 * t649;
t595 = t1147 * t649 + t1152 * t617;
t593 = -t1146 * t685 + t1151 * t686 + (-t1141 * t666 - t1142 * t667) * pkin(8);
t591 = -t1147 * t613 + t1152 * t641;
t590 = t1147 * t641 + t1152 * t613;
t588 = -t1141 * t647 + t1142 * t1242;
t586 = -t1141 * t643 + t1142 * t1243;
t585 = t1141 * t1243 + t1142 * t643;
t582 = -t1141 * t638 + t1142 * t1244;
t581 = t1141 * t1244 + t1142 * t638;
t580 = -pkin(1) * t666 - t1141 * t669 + t1142 * t1179;
t579 = -t1146 * t672 + t1151 * t625 + (-t1141 * t706 - t1142 * t707) * pkin(8);
t578 = -t1146 * t668 + t1151 * t624 + (-t1141 * t698 - t1142 * t699) * pkin(8);
t577 = -pkin(1) * t706 - t1141 * t620 + t1142 * t1180;
t576 = -pkin(1) * t698 - t1141 * t614 + t1142 * t1181;
t573 = -t1147 * t586 + t1152 * t609;
t572 = t1147 * t609 + t1152 * t586;
t571 = -t1146 * t615 + t1151 * t608 + (-t1141 * t676 - t1142 * t677) * pkin(8);
t566 = -t1147 * t582 + t1152 * t606;
t565 = t1147 * t606 + t1152 * t582;
t564 = -pkin(2) * t735 + pkin(9) * t691 + t1145 * t605 + t1150 * t601;
t562 = -pkin(2) * t729 + pkin(9) * t682 + t1145 * t602 + t1150 * t597;
t561 = -pkin(2) * t656 + pkin(9) * t639 + t1145 * t611 + t1150 * t607;
t560 = -pkin(1) * t676 - t1141 * t604 + t1142 * t1182;
t555 = -t1146 * t589 + t1151 * t569 + (-t1141 * t616 - t1142 * t617) * pkin(8);
t553 = -pkin(2) * t661 + pkin(9) * t644 + t1145 * t603 + t1150 * t559;
t552 = -t1146 * t584 + t1151 * t563 + (-t1141 * t612 - t1142 * t613) * pkin(8);
t551 = -t1146 * t583 + t1151 * t570 + (-t1141 * t581 - t1142 * t582) * pkin(8);
t549 = -t1141 * t567 + t1142 * t1245;
t548 = t1141 * t1245 + t1142 * t567;
t547 = -pkin(1) * t616 - t1141 * t564 + t1142 * t1184;
t546 = -pkin(1) * t612 - t1141 * t562 + t1142 * t1185;
t545 = -pkin(1) * t581 - t1141 * t561 + t1142 * t1183;
t543 = -t1147 * t549 + t1152 * t557;
t542 = t1147 * t557 + t1152 * t549;
t541 = -t1146 * t558 + t1151 * t556 + (-t1141 * t585 - t1142 * t586) * pkin(8);
t539 = -pkin(1) * t585 - t1141 * t553 + t1142 * t1186;
t538 = -pkin(2) * t575 + pkin(9) * t568 + t1145 * t554 + t1150 * t550;
t537 = -t1146 * t544 + t1151 * t540 + (-t1141 * t548 - t1142 * t549) * pkin(8);
t536 = -pkin(1) * t548 - t1141 * t538 + t1142 * t1187;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1123, -t1124, 0, t1080, 0, 0, 0, 0, 0, 0, t940, t930, t922, t859, 0, 0, 0, 0, 0, 0, -t760, -t765, t739, t675, 0, 0, 0, 0, 0, 0, t738, t760, t765, t630, 0, 0, 0, 0, 0, 0, t654, t663, t637, t566, 0, 0, 0, 0, 0, 0, t591, t596, t573, t543; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1124, -t1123, 0, t1079, 0, 0, 0, 0, 0, 0, t939, t929, t921, t858, 0, 0, 0, 0, 0, 0, t757, t762, t737, t674, 0, 0, 0, 0, 0, 0, t736, -t757, -t762, t629, 0, 0, 0, 0, 0, 0, t653, t662, t636, t565, 0, 0, 0, 0, 0, 0, t590, t595, t572, t542; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1005, t987, t971, t927, 0, 0, 0, 0, 0, 0, t825, t836, t787, t715, 0, 0, 0, 0, 0, 0, t786, -t825, -t836, t666, 0, 0, 0, 0, 0, 0, t698, t706, t676, t581, 0, 0, 0, 0, 0, 0, t612, t616, t585, t548; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1124, 0, -t1123, 0, t1255, -t1108, -t1079, -pkin(7) * t1079, -t1041 * t1147 + t1062 * t1152, t1012 * t1152 - t1147 * t973, -t1006 * t1147 + t1052 * t1152, -t1040 * t1147 + t1061 * t1152, -t1007 * t1147 + t1053 * t1152, t1077 * t1152 + t1147 * t1287, -pkin(7) * t939 - t1147 * t881 + t1152 * t888, -pkin(7) * t929 - t1147 * t870 + t1152 * t882, -pkin(7) * t921 - t1147 * t853 + t1152 * t841, -pkin(7) * t858 - t1147 * t860 + t1152 * t845, t1378, -t1147 * t804 + t1152 * t871, -t1147 * t847 + t1152 * t897, t1379, t1421, t1359, -t1147 * t673 + t1152 * t684 - t1431, -t1147 * t679 + t1152 * t689 - t1429, -pkin(7) * t737 - t1147 * t642 + t1152 * t664, -pkin(7) * t674 - t1147 * t626 + t1152 * t627, t1359, -t1147 * t846 + t1152 * t896, -t1421, t1378, -t1147 * t805 + t1152 * t872, t1379, -pkin(7) * t736 - t1147 * t621 + t1152 * t640, -t1147 * t631 + t1152 * t650 + t1431, -t1147 * t635 + t1152 * t657 + t1429, -pkin(7) * t629 - t1147 * t580 + t1152 * t593, -t1147 * t732 + t1152 * t775, -t1147 * t683 + t1152 * t730, -t1147 * t710 + t1152 * t754, -t1147 * t731 + t1152 * t774, -t1147 * t709 + t1152 * t753, -t1147 * t767 + t1152 * t835, -pkin(7) * t653 - t1147 * t576 + t1152 * t578, -pkin(7) * t662 - t1147 * t577 + t1152 * t579, -pkin(7) * t636 - t1147 * t560 + t1152 * t571, -pkin(7) * t565 - t1147 * t545 + t1152 * t551, -t1147 * t619 + t1152 * t646, -t1147 * t588 + t1152 * t610, -t1147 * t623 + t1152 * t652, -t1147 * t618 + t1152 * t645, -t1147 * t622 + t1152 * t651, -t1147 * t671 + t1152 * t712, -pkin(7) * t590 - t1147 * t546 + t1152 * t552, -pkin(7) * t595 - t1147 * t547 + t1152 * t555, -pkin(7) * t572 - t1147 * t539 + t1152 * t541, -pkin(7) * t542 - t1147 * t536 + t1152 * t537; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1123, 0, t1124, 0, t1108, t1255, t1080, pkin(7) * t1080, t1041 * t1152 + t1062 * t1147, t1012 * t1147 + t1152 * t973, t1006 * t1152 + t1052 * t1147, t1040 * t1152 + t1061 * t1147, t1007 * t1152 + t1053 * t1147, t1077 * t1147 - t1152 * t1287, pkin(7) * t940 + t1147 * t888 + t1152 * t881, pkin(7) * t930 + t1147 * t882 + t1152 * t870, pkin(7) * t922 + t1147 * t841 + t1152 * t853, pkin(7) * t859 + t1147 * t845 + t1152 * t860, t1376, t1147 * t871 + t1152 * t804, t1147 * t897 + t1152 * t847, t1377, -t1420, t1360, t1147 * t684 + t1152 * t673 - t1430, t1147 * t689 + t1152 * t679 - t1428, pkin(7) * t739 + t1147 * t664 + t1152 * t642, pkin(7) * t675 + t1147 * t627 + t1152 * t626, t1360, t1147 * t896 + t1152 * t846, t1420, t1376, t1147 * t872 + t1152 * t805, t1377, pkin(7) * t738 + t1147 * t640 + t1152 * t621, t1147 * t650 + t1152 * t631 + t1430, t1147 * t657 + t1152 * t635 + t1428, pkin(7) * t630 + t1147 * t593 + t1152 * t580, t1147 * t775 + t1152 * t732, t1147 * t730 + t1152 * t683, t1147 * t754 + t1152 * t710, t1147 * t774 + t1152 * t731, t1147 * t753 + t1152 * t709, t1147 * t835 + t1152 * t767, pkin(7) * t654 + t1147 * t578 + t1152 * t576, pkin(7) * t663 + t1147 * t579 + t1152 * t577, pkin(7) * t637 + t1147 * t571 + t1152 * t560, pkin(7) * t566 + t1147 * t551 + t1152 * t545, t1147 * t646 + t1152 * t619, t1147 * t610 + t1152 * t588, t1147 * t652 + t1152 * t623, t1147 * t645 + t1152 * t618, t1147 * t651 + t1152 * t622, t1147 * t712 + t1152 * t671, pkin(7) * t591 + t1147 * t552 + t1152 * t546, pkin(7) * t596 + t1147 * t555 + t1152 * t547, pkin(7) * t573 + t1147 * t541 + t1152 * t539, pkin(7) * t543 + t1147 * t537 + t1152 * t536; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1128, t1129, 0, 0, (t1141 * t1106 + t1151 * t1357) * t1146, t1142 * t1112 + t1141 * t1209, t1142 * t1071 + t1141 * t1205, (t1141 * t1107 - t1146 * t1357) * t1151, t1142 * t1072 + t1141 * t1206, t1142 * t1133, pkin(1) * t1008 - t1142 * t1043 + t1141 * t1250, pkin(1) * t988 - t1142 * t1044 + t1141 * t1251, pkin(1) * t972 + t1141 * t1173, pkin(1) * t928 + t1330 * t952, t1364, t1141 * t1213 + t1142 * t914, t1141 * t1218 + t1396, t1363, -t1409, t1340, t1141 * t1175 + t1142 * t842 + t1424, t1141 * t1174 + t1142 * t850 + t1422, pkin(1) * t789 + t1141 * t1171 + t1142 * t748, pkin(1) * t716 + t1141 * t1167 + t1142 * t771, t1340, t1141 * t1219 - t1396, t1409, t1364, t1141 * t1212 + t1142 * t915, t1363, pkin(1) * t788 + t1141 * t1178 + t1142 * t711, t1141 * t1177 + t1142 * t741 - t1424, t1141 * t1176 + t1142 * t746 - t1422, pkin(1) * t667 + t1141 * t1179 + t1142 * t669, t1141 * t1225 + t1142 * t855, t1141 * t1232 + t1142 * t776, t1141 * t1227 + t1142 * t811, t1141 * t1226 + t1142 * t854, t1141 * t1228 + t1142 * t810, t1141 * t1224 + t1142 * t901, pkin(1) * t699 + t1141 * t1181 + t1142 * t614, pkin(1) * t707 + t1141 * t1180 + t1142 * t620, pkin(1) * t677 + t1141 * t1182 + t1142 * t604, pkin(1) * t582 + t1141 * t1183 + t1142 * t561, t1141 * t1236 + t1142 * t702, t1141 * t1242 + t1142 * t647, t1141 * t1238 + t1142 * t694, t1141 * t1237 + t1142 * t701, t1141 * t1239 + t1142 * t693, t1141 * t1234 + t1142 * t768, pkin(1) * t613 + t1141 * t1185 + t1142 * t562, pkin(1) * t617 + t1141 * t1184 + t1142 * t564, pkin(1) * t586 + t1141 * t1186 + t1142 * t553, pkin(1) * t549 + t1141 * t1187 + t1142 * t538;];
tauB_reg  = t1;
