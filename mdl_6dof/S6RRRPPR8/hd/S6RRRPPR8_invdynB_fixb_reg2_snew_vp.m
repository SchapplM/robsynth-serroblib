% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 06:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRPPR8_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:17:14
% EndTime: 2019-05-07 06:18:22
% DurationCPUTime: 61.51s
% Computational Cost: add. (112215->859), mult. (245011->1254), div. (0->0), fcn. (187802->10), ass. (0->627)
t1083 = sin(qJ(1));
t1087 = cos(qJ(1));
t1079 = cos(pkin(6));
t1086 = cos(qJ(2));
t1081 = sin(qJ(3));
t1085 = cos(qJ(3));
t1233 = qJD(1) * t1079;
t1195 = qJD(2) + t1233;
t1078 = sin(pkin(6));
t1082 = sin(qJ(2));
t1218 = t1078 * t1082;
t1204 = qJD(1) * t1218;
t1031 = t1081 * t1195 + t1085 * t1204;
t1217 = t1078 * t1086;
t1064 = qJD(1) * t1217 - qJD(3);
t1008 = t1031 * t1064;
t1209 = t1082 * qJDD(1);
t1232 = qJD(1) * t1086;
t1041 = (qJD(2) * t1232 + t1209) * t1078;
t1192 = qJDD(1) * t1079 + qJDD(2);
t1193 = t1081 * t1041 - t1085 * t1192;
t961 = qJD(3) * t1031 + t1193;
t1309 = -t1008 + t961;
t1278 = t1064 ^ 2;
t1029 = t1081 * t1204 - t1085 * t1195;
t1280 = t1029 ^ 2;
t1296 = -t1280 - t1278;
t1210 = qJDD(1) * t1078;
t1042 = -qJD(2) * t1204 + t1086 * t1210;
t1036 = -qJDD(3) + t1042;
t975 = t1031 * t1029;
t946 = t1036 + t975;
t1362 = t1081 * t946;
t1351 = t1085 * t1296 + t1362;
t1376 = t1082 * t1351;
t1157 = t1086 * t1309 - t1376;
t1359 = t1085 * t946;
t863 = -t1081 * t1296 + t1359;
t1357 = t1078 * t863;
t737 = t1079 * t1157 - t1357;
t1375 = t1086 * t1351;
t787 = t1082 * t1309 + t1375;
t679 = -t1083 * t787 + t1087 * t737;
t1435 = pkin(7) * t679;
t681 = t1083 * t737 + t1087 * t787;
t1434 = pkin(7) * t681;
t1356 = t1079 * t863;
t734 = t1078 * t1157 + t1356;
t1433 = pkin(8) * (t1078 * t734 + t1079 * t737);
t1432 = pkin(1) * t734;
t1431 = pkin(1) * t737;
t1105 = -t1085 * t1041 - t1081 * t1192;
t1097 = -t1029 * qJD(3) - t1105;
t1222 = t1064 * t1029;
t1299 = t1222 + t1097;
t1279 = t1031 ^ 2;
t1298 = -t1279 - t1278;
t947 = t1036 - t975;
t1358 = t1085 * t947;
t864 = -t1081 * t1298 + t1358;
t1386 = t1086 * t864;
t1405 = -t1082 * t1299 - t1386;
t1361 = t1081 * t947;
t861 = t1085 * t1298 + t1361;
t1389 = t1078 * t861;
t1387 = t1082 * t864;
t1406 = t1086 * t1299 - t1387;
t1410 = t1079 * t1406 + t1389;
t1418 = t1083 * t1405 + t1087 * t1410;
t1429 = pkin(7) * t1418;
t1419 = -t1083 * t1410 + t1087 * t1405;
t1428 = pkin(7) * t1419;
t1294 = t1280 - t1279;
t1256 = t1081 * t1299;
t1213 = qJD(3) - t1064;
t923 = t1031 * t1213 + t1193;
t825 = t1085 * t923 + t1256;
t1165 = t1082 * t825 - t1086 * t1294;
t1344 = -t1081 * t923 + t1085 * t1299;
t706 = t1078 * t1344 + t1079 * t1165;
t775 = t1082 * t1294 + t1086 * t825;
t1421 = t1083 * t706 - t1087 * t775;
t1420 = t1083 * t775 + t1087 * t706;
t1295 = t1280 - t1278;
t1373 = t1085 * t1295 + t1361;
t1212 = qJD(3) + t1064;
t919 = t1031 * t1212 + t1193;
t1381 = -t1082 * t919 + t1086 * t1373;
t1383 = t1082 * t1373 + t1086 * t919;
t872 = -t1081 * t1295 + t1358;
t1394 = t1078 * t872 + t1079 * t1383;
t1427 = t1083 * t1394 - t1087 * t1381;
t1402 = t1083 * t1381 + t1087 * t1394;
t1426 = pkin(8) * t787;
t1425 = pkin(1) * t1410;
t1388 = t1079 * t861;
t1411 = t1078 * t1406 - t1388;
t1424 = pkin(1) * t1411;
t1417 = (-t1078 * t1411 - t1079 * t1410) * pkin(8);
t1416 = pkin(8) * t1405;
t1404 = t1078 * t1165 - t1079 * t1344;
t1393 = t1078 * t1383 - t1079 * t872;
t1392 = pkin(2) * t861;
t1407 = (t1036 + t947) * qJ(4) + t1392;
t1067 = g(1) * t1087 + g(2) * t1083;
t1088 = qJD(1) ^ 2;
t1037 = -pkin(1) * t1088 + pkin(8) * t1210 - t1067;
t1269 = pkin(2) * t1086;
t1189 = -pkin(9) * t1082 - t1269;
t1234 = qJD(1) * t1078;
t1040 = t1189 * t1234;
t1066 = t1083 * g(1) - t1087 * g(2);
t1267 = pkin(8) * t1078;
t1103 = qJDD(1) * pkin(1) + t1088 * t1267 + t1066;
t1100 = t1079 * t1103;
t1094 = -g(3) * t1218 + t1082 * t1100;
t1191 = t1195 ^ 2;
t891 = t1192 * pkin(9) - t1191 * pkin(2) + (t1040 * t1234 + t1037) * t1086 + t1094;
t1148 = t1195 * qJD(1);
t1130 = t1082 * t1148;
t1131 = t1086 * t1148;
t1266 = t1079 * g(3);
t892 = -t1042 * pkin(2) - t1041 * pkin(9) - t1266 + (pkin(2) * t1130 - pkin(9) * t1131 - t1103) * t1078;
t795 = t1081 * t891 - t1085 * t892;
t1127 = t1036 * pkin(3) - qJ(4) * t1278 + qJDD(4) + t795;
t970 = pkin(3) * t1029 - qJ(4) * t1031;
t1102 = t1031 * t970 + t1127;
t1368 = pkin(2) * t863;
t684 = pkin(3) * t946 - qJ(4) * t1296 + t1102 + t1368;
t1300 = t1222 - t1097;
t1297 = -t1279 + t1278;
t1374 = -t1081 * t1297 - t1359;
t1380 = -t1082 * t1300 + t1086 * t1374;
t1372 = t1085 * t1297 - t1362;
t1382 = t1082 * t1374 + t1086 * t1300;
t1395 = -t1078 * t1372 + t1079 * t1382;
t1371 = -t1083 * t1395 + t1087 * t1380;
t1370 = t1083 * t1380 + t1087 * t1395;
t1391 = pkin(9) * t861;
t1367 = pkin(9) * t863;
t1390 = pkin(9) * t864;
t1379 = pkin(9) * t1351;
t1348 = t1078 * t1382 + t1079 * t1372;
t1369 = 2 * qJD(4);
t1293 = t1280 + t1279;
t1343 = pkin(2) * t1293;
t1366 = qJ(4) * t919;
t1365 = qJ(4) * t1293;
t1364 = qJ(4) * t1299;
t1363 = t1081 * t919;
t1360 = t1085 * t919;
t1333 = t1082 * t1293;
t1323 = t1086 * t1293;
t1194 = t1082 * t1037 - t1086 * t1100;
t890 = (qJD(1) * t1040 * t1082 + g(3) * t1086) * t1078 - t1192 * pkin(2) - t1191 * pkin(9) + t1194;
t1092 = t961 * pkin(3) - t1364 + t890;
t1091 = t1031 * t1369 - t1092;
t717 = -(t1309 - t1008) * pkin(3) + t1091;
t716 = pkin(3) * t1008 + t1091 + t1364;
t1220 = t1064 * t1085;
t1201 = t1029 * t1220;
t1133 = t1081 * t961 - t1201;
t1203 = t1082 * t975;
t1287 = t1086 * t1133 - t1203;
t1221 = t1064 * t1081;
t1186 = -t1029 * t1221 - t1085 * t961;
t1202 = t1086 * t975;
t1288 = t1082 * t1133 + t1202;
t1315 = -t1078 * t1186 + t1079 * t1288;
t1346 = -t1083 * t1315 + t1087 * t1287;
t1345 = t1083 * t1287 + t1087 * t1315;
t1080 = sin(qJ(6));
t1084 = cos(qJ(6));
t984 = t1029 * t1080 - t1084 * t1064;
t986 = t1029 * t1084 + t1064 * t1080;
t930 = t986 * t984;
t958 = qJDD(6) + t1097;
t1310 = -t930 + t958;
t1337 = t1080 * t1310;
t991 = t1031 * t1221;
t1264 = t1085 * t1097 + t991;
t1285 = t1086 * t1264 + t1203;
t1329 = t1083 * t1285;
t1327 = t1084 * t1310;
t1319 = t1087 * t1285;
t1316 = t1078 * t1288 + t1079 * t1186;
t1185 = -t991 + t1201;
t1216 = t1079 * t1082;
t1126 = (t1029 * t1081 + t1031 * t1085) * t1064;
t1289 = t1079 * t1086 * t1036 - t1078 * t1126;
t1281 = t1185 * t1216 + t1289;
t1227 = t1036 * t1082;
t1286 = t1086 * t1185 - t1227;
t1312 = t1083 * t1286 + t1087 * t1281;
t1311 = -t1083 * t1281 + t1087 * t1286;
t1075 = t1078 ^ 2;
t1219 = t1075 * t1088;
t1306 = t1075 * (-t1079 * t1088 + t1148);
t1184 = -t1031 * t1220 + t1081 * t1097;
t1305 = t1078 * t1184;
t1303 = t1079 * t1184;
t990 = pkin(4) * t1064 - qJ(5) * t1031;
t1301 = t1031 * t990 + qJDD(5);
t796 = t1081 * t892 + t1085 * t891;
t1205 = pkin(3) * t1278 - t796;
t1129 = t1029 * t970 + t1064 * t1369 + t1205;
t1292 = t1036 * pkin(4) + qJ(5) * t1300;
t1291 = t961 * pkin(4) - t1301;
t1290 = t1036 * t1217 + t1079 * t1126;
t1108 = -t1079 * t1202 + t1216 * t1264 - t1305;
t1284 = t1087 * t1108 + t1329;
t1283 = -t1083 * t1108 + t1319;
t1282 = t1185 * t1218 + t1290;
t982 = t984 ^ 2;
t983 = t986 ^ 2;
t1022 = qJD(6) + t1031;
t1019 = t1022 ^ 2;
t1277 = -2 * qJD(4);
t1276 = -2 * qJD(5);
t1275 = 2 * qJD(5);
t1274 = pkin(3) + pkin(4);
t1273 = -pkin(4) - pkin(10);
t1271 = pkin(5) + qJ(4);
t1270 = pkin(2) * t1082;
t1268 = pkin(3) * t1085;
t1265 = t961 * qJ(5);
t1263 = qJ(4) * t1085;
t1262 = t1022 * t984;
t1261 = t1280 * qJ(5);
t1260 = t1036 * qJ(4);
t1021 = t1280 * pkin(4);
t1095 = -t1260 + t1265 - t1021 + (t1277 - t990) * t1064 - t1205;
t1236 = t1275 - t970;
t971 = pkin(5) * t1031 - pkin(10) * t1029;
t675 = -t1036 * pkin(5) - t1278 * pkin(10) + (t971 + t1236) * t1029 + t1095;
t1259 = t1080 * t675;
t846 = t930 + t958;
t1258 = t1080 * t846;
t1257 = t1081 * t890;
t1255 = t1081 * t1300;
t1245 = t1084 * t675;
t1244 = t1084 * t846;
t1243 = t1085 * t890;
t1235 = pkin(3) - t1273;
t1229 = t1022 * t1080;
t1228 = t1022 * t1084;
t1063 = t1086 * t1082 * t1219;
t1038 = -t1063 + t1192;
t1226 = t1038 * t1082;
t1225 = t1038 * t1086;
t1039 = t1063 + t1192;
t1224 = t1039 * t1082;
t1223 = t1039 * t1086;
t1010 = t1078 * t1103 + t1266;
t1215 = t1082 * t1010;
t1214 = t1086 * t1010;
t1076 = t1082 ^ 2;
t1077 = t1086 ^ 2;
t1211 = t1076 + t1077;
t1208 = -t983 - t1019;
t1207 = t1081 * t930;
t1206 = t1085 * t930;
t858 = -t984 * qJD(6) + t1080 * t1036 + t1084 * t961;
t1200 = t1076 * t1219;
t1199 = t1077 * t1219;
t1198 = -qJ(4) * t1081 - pkin(2);
t1197 = -pkin(4) * t1029 - t970;
t711 = t1081 * t795 + t1085 * t796;
t1196 = -t1084 * t1036 + t1080 * t961;
t1012 = -t1066 * t1083 - t1087 * t1067;
t1190 = t1276 - t1197;
t1062 = qJDD(1) * t1087 - t1083 * t1088;
t1188 = -pkin(7) * t1062 - g(3) * t1083;
t1020 = -t1200 - t1191;
t969 = -t1020 * t1082 - t1225;
t1183 = pkin(8) * t969 - t1215;
t1046 = -t1191 - t1199;
t980 = t1046 * t1086 - t1224;
t1182 = pkin(8) * t980 + t1214;
t661 = t1091 + t1273 * t961 + (t1029 * pkin(5) + (pkin(3) + pkin(10)) * t1031) * t1064 + t1097 * pkin(5) - t1261 + t1301;
t1101 = t1127 + t1292;
t668 = -t1278 * pkin(5) + t1036 * pkin(10) + (t1190 - t971) * t1031 + t1101;
t617 = t1080 * t668 - t1084 * t661;
t618 = t1080 * t661 + t1084 * t668;
t575 = t1080 * t618 - t1084 * t617;
t576 = t1080 * t617 + t1084 * t618;
t710 = t1081 * t796 - t1085 * t795;
t567 = t1081 * t576 + t1085 * t675;
t1181 = t1082 * t567 + t1086 * t575;
t689 = t1031 * t1190 + t1101;
t692 = t1029 * t1236 + t1095;
t632 = t1081 * t689 + t1085 * t692;
t752 = (-pkin(3) * t1064 + t1277) * t1031 + t1092;
t702 = t1261 + t752 + t1291;
t1180 = t1082 * t632 - t1086 * t702;
t742 = -t1260 - t1129;
t660 = t1081 * t1102 + t1085 * t742;
t1179 = t1082 * t660 - t1086 * t752;
t1110 = (-qJD(6) + t1022) * t986 - t1196;
t814 = -t858 - t1262;
t722 = -t1080 * t814 + t1084 * t1110;
t859 = -t982 - t983;
t688 = t1081 * t722 + t1085 * t859;
t720 = t1080 * t1110 + t1084 * t814;
t1178 = t1082 * t688 + t1086 * t720;
t1142 = t858 - t1262;
t810 = (qJD(6) + t1022) * t986 + t1196;
t721 = t1080 * t1142 + t1084 * t810;
t929 = t983 - t982;
t694 = -t1081 * t721 + t1085 * t929;
t719 = -t1080 * t810 + t1084 * t1142;
t1177 = t1082 * t694 + t1086 * t719;
t879 = -t1019 - t982;
t771 = t1084 * t879 - t1337;
t704 = t1081 * t771 + t1085 * t810;
t770 = t1080 * t879 + t1327;
t1176 = t1082 * t704 + t1086 * t770;
t778 = -t1080 * t1208 - t1244;
t709 = t1081 * t778 + t1085 * t1142;
t777 = t1084 * t1208 - t1258;
t1175 = t1082 * t709 + t1086 * t777;
t1174 = t1082 * t711 - t1086 * t890;
t945 = -t983 + t1019;
t792 = t1080 * t945 - t1327;
t714 = -t1081 * t792 - t1085 * t814;
t790 = t1084 * t945 + t1337;
t1173 = t1082 * t714 + t1086 * t790;
t944 = t982 - t1019;
t793 = -t1084 * t944 + t1258;
t715 = -t1081 * t793 + t1085 * t1110;
t791 = t1080 * t944 + t1244;
t1172 = t1082 * t715 + t1086 * t791;
t857 = -qJD(6) * t986 - t1196;
t806 = t1080 * t857 - t1228 * t984;
t756 = -t1081 * t806 - t1206;
t805 = t1084 * t857 + t1229 * t984;
t1171 = t1082 * t756 + t1086 * t805;
t808 = -t1084 * t858 + t1229 * t986;
t757 = -t1081 * t808 + t1206;
t807 = t1080 * t858 + t1228 * t986;
t1170 = t1082 * t757 + t1086 * t807;
t850 = (-t1080 * t986 + t1084 * t984) * t1022;
t804 = -t1081 * t850 + t1085 * t958;
t849 = (-t1080 * t984 - t1084 * t986) * t1022;
t1169 = t1082 * t804 + t1086 * t849;
t912 = -t1008 - t961;
t822 = t1085 * t912 - t1255;
t1168 = t1082 * t822 + t1323;
t914 = -t1029 * t1212 - t1105;
t823 = -t1081 * t914 + t1360;
t1167 = t1082 * t823 - t1323;
t826 = -t1255 - t1360;
t1164 = t1082 * t826 + t1323;
t1160 = -t1086 * t923 + t1376;
t915 = t1029 * t1213 + t1105;
t1159 = t1086 * t915 + t1387;
t967 = g(3) * t1217 + t1194;
t968 = t1086 * t1037 + t1094;
t1150 = t1082 * t968 - t1086 * t967;
t860 = t1082 * t967 + t1086 * t968;
t1049 = t1078 * t1130;
t1003 = t1042 - t1049;
t1050 = t1078 * t1131;
t999 = t1050 + t1041;
t1149 = t1003 * t1082 + t1086 * t999;
t1001 = -t1050 + t1041;
t1002 = t1042 + t1049;
t1147 = -t1001 * t1086 + t1002 * t1082;
t1146 = t1020 * t1086 - t1226;
t1045 = -t1191 + t1199;
t1145 = t1045 * t1082 + t1225;
t1044 = t1191 - t1200;
t1144 = t1044 * t1086 + t1224;
t1143 = t1046 * t1082 + t1223;
t1011 = t1066 * t1087 - t1067 * t1083;
t1140 = t1078 * t1192;
t1134 = t1082 * t1264 - t1202;
t1132 = t1078 * t1148;
t1128 = -t1078 * t1202 + t1218 * t1264 + t1303;
t545 = -qJ(5) * t675 + t1235 * t575;
t546 = -qJ(5) * t576 + t1271 * t575;
t566 = t1081 * t675 - t1085 * t576;
t531 = -pkin(9) * t566 - t1081 * t545 + t1085 * t546;
t539 = -pkin(2) * t566 + t1235 * t576 - t1271 * t675;
t550 = -t1082 * t575 + t1086 * t567;
t1125 = pkin(8) * t550 + t1082 * t531 + t1086 * t539;
t560 = -qJ(5) * t859 + t1235 * t720 + t575;
t620 = -qJ(5) * t722 + t1271 * t720;
t687 = t1081 * t859 - t1085 * t722;
t551 = -pkin(9) * t687 - t1081 * t560 + t1085 * t620;
t552 = -pkin(2) * t687 + t1235 * t722 - t1271 * t859 + t576;
t639 = -t1082 * t720 + t1086 * t688;
t1124 = pkin(8) * t639 + t1082 * t551 + t1086 * t552;
t588 = -qJ(5) * t771 + t1271 * t770 - t617;
t602 = -qJ(5) * t810 + t1235 * t770 - t1259;
t703 = t1081 * t810 - t1085 * t771;
t556 = -pkin(9) * t703 - t1081 * t602 + t1085 * t588;
t586 = -pkin(2) * t703 + t1235 * t771 - t1271 * t810 - t1245;
t654 = -t1082 * t770 + t1086 * t704;
t1123 = pkin(8) * t654 + t1082 * t556 + t1086 * t586;
t589 = -qJ(5) * t778 + t1271 * t777 - t618;
t606 = -qJ(5) * t1142 + t1235 * t777 - t1245;
t708 = t1081 * t1142 - t1085 * t778;
t557 = -pkin(9) * t708 - t1081 * t606 + t1085 * t589;
t587 = -pkin(2) * t708 - t1142 * t1271 + t1235 * t778 + t1259;
t657 = -t1082 * t777 + t1086 * t709;
t1122 = pkin(8) * t657 + t1082 * t557 + t1086 * t587;
t611 = -qJ(5) * t692 - t1274 * t702;
t631 = t1081 * t692 - t1085 * t689;
t638 = -qJ(4) * t702 - qJ(5) * t689;
t564 = -pkin(9) * t631 - t1081 * t611 + t1085 * t638;
t577 = -pkin(2) * t631 - qJ(4) * t692 + t1274 * t689;
t607 = t1082 * t702 + t1086 * t632;
t1121 = pkin(8) * t607 + t1082 * t564 + t1086 * t577;
t1098 = t1029 * t1276 + t1064 * t990 + t1021 + t1129;
t653 = t1260 - t1274 * t1293 + (-t919 - t961) * qJ(5) + t1098;
t1017 = t1031 * t1275;
t658 = qJ(5) * t914 + t1031 * t1197 + t1017 - t1101 - t1365;
t817 = t1085 * t914 + t1363;
t601 = -pkin(9) * t817 - t1081 * t653 + t1085 * t658;
t695 = -pkin(2) * t817 - t1274 * t914 - t1366;
t762 = t1086 * t823 + t1333;
t1120 = pkin(8) * t762 + t1082 * t601 + t1086 * t695;
t659 = t1081 * t742 - t1085 * t1102;
t614 = -pkin(2) * t659 + pkin(3) * t1102 - qJ(4) * t742;
t615 = -pkin(9) * t659 + (pkin(3) * t1081 - t1263) * t752;
t635 = t1082 * t752 + t1086 * t660;
t1119 = pkin(8) * t635 + t1082 * t615 + t1086 * t614;
t667 = (-t1280 - t1298) * qJ(5) - t1291 + t716;
t772 = qJ(5) * t947 + t1274 * t1299;
t630 = -t1081 * t772 + t1085 * t667 + t1391;
t642 = t1274 * t1298 + t1098 - t1265 + t1407;
t1118 = t1082 * t630 + t1086 * t642 + t1416;
t718 = pkin(3) * t1293 + t742;
t723 = t1102 + t1365;
t906 = t1085 * t1300;
t820 = t906 - t1363;
t633 = -pkin(9) * t820 - t1081 * t718 + t1085 * t723;
t724 = -pkin(2) * t820 - pkin(3) * t1300 + t1366;
t763 = t1086 * t826 - t1333;
t1117 = pkin(8) * t763 + t1082 * t633 + t1086 * t724;
t655 = (t1280 + t1296) * qJ(5) + (t1309 + t961) * pkin(4) - t1301 - t717;
t833 = qJ(4) * t1309 + qJ(5) * t946;
t634 = -t1081 * t655 + t1085 * t833 - t1367;
t641 = t1017 + (-t946 - t975) * pkin(4) - t1292 - t684;
t1116 = t1082 * t634 + t1086 * t641 - t1426;
t662 = -pkin(3) * t1256 + t1085 * t716 + t1391;
t676 = pkin(3) * t1298 + t1129 + t1407;
t1115 = t1082 * t662 + t1086 * t676 + t1416;
t665 = -t1081 * t717 - t1263 * t1309 + t1367;
t1114 = t1082 * t665 + t1086 * t684 + t1426;
t738 = t795 + t1368;
t779 = t1257 + t1367;
t784 = t1082 * t923 + t1375;
t1113 = pkin(8) * t784 + t1082 * t779 + t1086 * t738;
t741 = t796 - t1392;
t785 = t1243 - t1391;
t786 = -t1082 * t915 + t1386;
t1112 = pkin(8) * t786 + t1082 * t785 + t1086 * t741;
t931 = t1001 * t1082 + t1002 * t1086;
t1111 = pkin(8) * t931 + t860;
t816 = t1081 * t912 + t906;
t666 = -pkin(9) * t816 - t710;
t761 = t1086 * t822 - t1333;
t1109 = pkin(8) * t761 + t1082 * t666 - t1269 * t816;
t685 = t1082 * t890 + t1086 * t711;
t1106 = pkin(8) * t685 + t1189 * t710;
t1099 = t1078 * t1219 + t1079 * t1132;
t1061 = qJDD(1) * t1083 + t1087 * t1088;
t1048 = t1211 * t1219;
t1047 = (t1076 - t1077) * t1219;
t1043 = -pkin(7) * t1061 + g(3) * t1087;
t1009 = t1195 * t1211 * t1234;
t1000 = (t1209 + (0.2e1 * qJD(2) + t1233) * t1232) * t1078;
t989 = t1086 * t1041 - t1076 * t1132;
t988 = -t1082 * t1042 - t1077 * t1132;
t979 = t1045 * t1086 - t1226;
t978 = -t1044 * t1082 + t1223;
t963 = (t1079 * t1041 + t1086 * t1099) * t1082;
t962 = (t1079 * t1042 - t1082 * t1099) * t1086;
t937 = (t1029 * t1085 - t1031 * t1081) * t1064;
t932 = t1003 * t1086 - t1082 * t999;
t928 = t1078 * t1003 + t1079 * t1143;
t927 = -t1078 * t1002 + t1079 * t1145;
t926 = -t1078 * t1001 + t1079 * t1144;
t925 = -t1079 * t1003 + t1078 * t1143;
t908 = -t1078 * t1000 + t1079 * t1146;
t907 = t1079 * t1000 + t1078 * t1146;
t894 = t1086 * t937 - t1227;
t888 = -t1078 * t1047 + t1079 * t1149;
t887 = t1078 * t1048 + t1079 * t1147;
t886 = -t1079 * t1048 + t1078 * t1147;
t848 = -t1083 * t928 + t1087 * t980;
t847 = t1083 * t980 + t1087 * t928;
t837 = -t1083 * t908 + t1087 * t969;
t836 = t1083 * t969 + t1087 * t908;
t835 = t1078 * t1010 + t1079 * t1150;
t834 = -t1079 * t1010 + t1078 * t1150;
t831 = t1216 * t937 + t1289;
t829 = -t1083 * t887 + t1087 * t931;
t828 = t1083 * t931 + t1087 * t887;
t803 = t1081 * t958 + t1085 * t850;
t789 = -t1215 + (-t1078 * t925 - t1079 * t928) * pkin(8);
t781 = -t1214 + (-t1078 * t907 - t1079 * t908) * pkin(8);
t780 = -pkin(1) * t925 + t1078 * t967 + t1079 * t1182;
t773 = -pkin(1) * t907 + t1078 * t968 + t1079 * t1183;
t767 = t1079 * t1134 - t1305;
t760 = pkin(8) * t1079 * t860 - pkin(1) * t834;
t759 = -t1083 * t835 + t1087 * t860;
t758 = t1083 * t860 + t1087 * t835;
t755 = t1085 * t808 + t1207;
t754 = t1085 * t806 - t1207;
t753 = -pkin(1) * t886 + t1079 * t1111;
t750 = pkin(2) * t915 + t1257 + t1390;
t743 = (-t1078 * t834 - t1079 * t835) * pkin(8);
t740 = -pkin(2) * t923 - t1243 + t1379;
t739 = (-t1078 * t886 - t1079 * t887) * pkin(8) - t1150;
t735 = t1079 * t1159 - t1389;
t732 = t1078 * t1159 + t1388;
t731 = -t1082 * t849 + t1086 * t804;
t730 = t1079 * t1160 + t1357;
t727 = t1078 * t1160 - t1356;
t713 = t1081 * t1110 + t1085 * t793;
t712 = -t1081 * t814 + t1085 * t792;
t701 = -t1078 * t820 + t1079 * t1164;
t700 = -t1078 * t817 + t1079 * t1167;
t699 = -t1078 * t816 + t1079 * t1168;
t698 = t1078 * t1164 + t1079 * t820;
t697 = t1078 * t1167 + t1079 * t817;
t696 = t1078 * t1168 + t1079 * t816;
t693 = t1081 * t929 + t1085 * t721;
t691 = -t1082 * t807 + t1086 * t757;
t690 = -t1082 * t805 + t1086 * t756;
t686 = -pkin(2) * t890 + pkin(9) * t711;
t683 = -t1078 * t803 + t1079 * t1169;
t680 = -t1083 * t735 + t1087 * t786;
t677 = t1083 * t786 + t1087 * t735;
t674 = -t1083 * t730 + t1087 * t784;
t671 = t1083 * t784 + t1087 * t730;
t664 = -t1082 * t791 + t1086 * t715;
t663 = -t1082 * t790 + t1086 * t714;
t656 = pkin(9) * t822 + t1343 + t711;
t652 = t1085 * t717 + t1198 * t1309 + t1379;
t651 = -t1390 + t1081 * t716 + (pkin(2) + t1268) * t1299;
t650 = -t1083 * t701 + t1087 * t763;
t649 = -t1083 * t700 + t1087 * t762;
t648 = -t1083 * t699 + t1087 * t761;
t647 = t1083 * t763 + t1087 * t701;
t646 = t1083 * t762 + t1087 * t700;
t645 = t1083 * t761 + t1087 * t699;
t644 = -t1078 * t755 + t1079 * t1170;
t643 = -t1078 * t754 + t1079 * t1171;
t640 = -t1082 * t719 + t1086 * t694;
t637 = -t1078 * t710 + t1079 * t1174;
t636 = t1078 * t1174 + t1079 * t710;
t629 = pkin(9) * t826 + t1081 * t723 + t1085 * t718 + t1343;
t628 = -t1078 * t713 + t1079 * t1172;
t627 = -t1078 * t712 + t1079 * t1173;
t626 = pkin(2) * t1309 + t1081 * t833 + t1085 * t655 - t1379;
t625 = pkin(2) * t1299 + t1081 * t667 + t1085 * t772 - t1390;
t624 = -t1078 * t708 + t1079 * t1175;
t623 = t1078 * t1175 + t1079 * t708;
t622 = -t1078 * t703 + t1079 * t1176;
t621 = t1078 * t1176 + t1079 * t703;
t619 = -t1082 * t741 + t1086 * t785 + (-t1078 * t732 - t1079 * t735) * pkin(8);
t613 = -t1082 * t738 + t1086 * t779 + (-t1078 * t727 - t1079 * t730) * pkin(8);
t612 = -t1078 * t693 + t1079 * t1177;
t610 = -pkin(1) * t732 - t1078 * t750 + t1079 * t1112;
t609 = -t1078 * t687 + t1079 * t1178;
t608 = t1078 * t1178 + t1079 * t687;
t605 = -t1083 * t637 + t1087 * t685;
t604 = t1083 * t685 + t1087 * t637;
t603 = -pkin(1) * t727 - t1078 * t740 + t1079 * t1113;
t600 = pkin(9) * t660 + (t1198 - t1268) * t752;
t599 = -t1078 * t659 + t1079 * t1179;
t598 = t1078 * t1179 + t1079 * t659;
t597 = pkin(9) * t823 + t1081 * t658 + t1085 * t653 - t1343;
t596 = t816 * t1270 + t1086 * t666 + (-t1078 * t696 - t1079 * t699) * pkin(8);
t595 = -t1083 * t624 + t1087 * t657;
t594 = t1083 * t657 + t1087 * t624;
t593 = -t1082 * t684 + t1086 * t665 + t1433;
t592 = -t1083 * t622 + t1087 * t654;
t591 = t1083 * t654 + t1087 * t622;
t590 = -t1082 * t676 + t1086 * t662 + t1417;
t585 = -pkin(1) * t696 - t1078 * t656 + t1079 * t1109;
t584 = -t1082 * t724 + t1086 * t633 + (-t1078 * t698 - t1079 * t701) * pkin(8);
t583 = -t1083 * t609 + t1087 * t639;
t582 = t1083 * t639 + t1087 * t609;
t581 = -t1078 * t631 + t1079 * t1180;
t580 = t1078 * t1180 + t1079 * t631;
t579 = -t1078 * t652 + t1079 * t1114 + t1432;
t578 = -t1082 * t641 + t1086 * t634 - t1433;
t574 = -t1082 * t642 + t1086 * t630 + t1417;
t573 = -t1078 * t651 + t1079 * t1115 - t1424;
t572 = -t1083 * t599 + t1087 * t635;
t571 = t1083 * t635 + t1087 * t599;
t570 = (-pkin(9) * t1086 + t1270) * t710 + (-t1078 * t636 - t1079 * t637) * pkin(8);
t569 = -pkin(1) * t636 - t1078 * t686 + t1079 * t1106;
t568 = -t1082 * t695 + t1086 * t601 + (-t1078 * t697 - t1079 * t700) * pkin(8);
t565 = -pkin(1) * t698 - t1078 * t629 + t1079 * t1117;
t563 = -t1078 * t626 + t1079 * t1116 - t1432;
t562 = -t1078 * t625 + t1079 * t1118 - t1424;
t561 = -pkin(2) * t702 + pkin(9) * t632 + t1081 * t638 + t1085 * t611;
t559 = -t1083 * t581 + t1087 * t607;
t558 = t1083 * t607 + t1087 * t581;
t555 = pkin(2) * t777 + pkin(9) * t709 + t1081 * t589 + t1085 * t606;
t554 = pkin(2) * t770 + pkin(9) * t704 + t1081 * t588 + t1085 * t602;
t553 = -pkin(1) * t697 - t1078 * t597 + t1079 * t1120;
t549 = -t1082 * t614 + t1086 * t615 + (-t1078 * t598 - t1079 * t599) * pkin(8);
t548 = pkin(2) * t720 + pkin(9) * t688 + t1081 * t620 + t1085 * t560;
t547 = -pkin(1) * t598 - t1078 * t600 + t1079 * t1119;
t544 = -t1082 * t587 + t1086 * t557 + (-t1078 * t623 - t1079 * t624) * pkin(8);
t543 = -t1082 * t586 + t1086 * t556 + (-t1078 * t621 - t1079 * t622) * pkin(8);
t542 = -t1078 * t566 + t1079 * t1181;
t541 = t1078 * t1181 + t1079 * t566;
t540 = -t1082 * t577 + t1086 * t564 + (-t1078 * t580 - t1079 * t581) * pkin(8);
t538 = -pkin(1) * t623 - t1078 * t555 + t1079 * t1122;
t537 = -pkin(1) * t621 - t1078 * t554 + t1079 * t1123;
t536 = -t1082 * t552 + t1086 * t551 + (-t1078 * t608 - t1079 * t609) * pkin(8);
t535 = -pkin(1) * t580 - t1078 * t561 + t1079 * t1121;
t534 = -t1083 * t542 + t1087 * t550;
t533 = t1083 * t550 + t1087 * t542;
t532 = -pkin(1) * t608 - t1078 * t548 + t1079 * t1124;
t530 = pkin(2) * t575 + pkin(9) * t567 + t1081 * t546 + t1085 * t545;
t529 = -t1082 * t539 + t1086 * t531 + (-t1078 * t541 - t1079 * t542) * pkin(8);
t528 = -pkin(1) * t541 - t1078 * t530 + t1079 * t1125;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1061, -t1062, 0, t1012, 0, 0, 0, 0, 0, 0, t848, t837, t829, t759, 0, 0, 0, 0, 0, 0, t674, t680, t648, t605, 0, 0, 0, 0, 0, 0, t681, t650, t1419, t572, 0, 0, 0, 0, 0, 0, t1419, -t681, t649, t559, 0, 0, 0, 0, 0, 0, t592, t595, t583, t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1062, -t1061, 0, t1011, 0, 0, 0, 0, 0, 0, t847, t836, t828, t758, 0, 0, 0, 0, 0, 0, t671, t677, t645, t604, 0, 0, 0, 0, 0, 0, -t679, t647, t1418, t571, 0, 0, 0, 0, 0, 0, t1418, t679, t646, t558, 0, 0, 0, 0, 0, 0, t591, t594, t582, t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t925, t907, t886, t834, 0, 0, 0, 0, 0, 0, t727, t732, t696, t636, 0, 0, 0, 0, 0, 0, -t734, t698, t1411, t598, 0, 0, 0, 0, 0, 0, t1411, t734, t697, t580, 0, 0, 0, 0, 0, 0, t621, t623, t608, t541; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1062, 0, -t1061, 0, t1188, -t1043, -t1011, -pkin(7) * t1011, -t1083 * t963 + t1087 * t989, -t1083 * t888 + t1087 * t932, -t1083 * t926 + t1087 * t978, -t1083 * t962 + t1087 * t988, -t1083 * t927 + t1087 * t979, t1087 * t1009 + t1083 * t1140, -pkin(7) * t847 - t1083 * t780 + t1087 * t789, -pkin(7) * t836 - t1083 * t773 + t1087 * t781, -pkin(7) * t828 - t1083 * t753 + t1087 * t739, -pkin(7) * t758 - t1083 * t760 + t1087 * t743, t1283, t1421, t1371, t1346, -t1427, t1311, -pkin(7) * t671 - t1083 * t603 + t1087 * t613, -pkin(7) * t677 - t1083 * t610 + t1087 * t619, -pkin(7) * t645 - t1083 * t585 + t1087 * t596, -pkin(7) * t604 - t1083 * t569 + t1087 * t570, t1283, t1371, -t1421, t1311, t1427, t1346, -t1083 * t579 + t1087 * t593 + t1435, -pkin(7) * t647 - t1083 * t565 + t1087 * t584, -t1083 * t573 + t1087 * t590 - t1429, -pkin(7) * t571 - t1083 * t547 + t1087 * t549, t1346, t1421, -t1427, -t1083 * t767 + t1319, t1371, -t1083 * t831 + t1087 * t894, -t1083 * t562 + t1087 * t574 - t1429, -t1083 * t563 + t1087 * t578 - t1435, -pkin(7) * t646 - t1083 * t553 + t1087 * t568, -pkin(7) * t558 - t1083 * t535 + t1087 * t540, -t1083 * t644 + t1087 * t691, -t1083 * t612 + t1087 * t640, -t1083 * t627 + t1087 * t663, -t1083 * t643 + t1087 * t690, -t1083 * t628 + t1087 * t664, -t1083 * t683 + t1087 * t731, -pkin(7) * t591 - t1083 * t537 + t1087 * t543, -pkin(7) * t594 - t1083 * t538 + t1087 * t544, -pkin(7) * t582 - t1083 * t532 + t1087 * t536, -pkin(7) * t533 - t1083 * t528 + t1087 * t529; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1061, 0, t1062, 0, t1043, t1188, t1012, pkin(7) * t1012, t1083 * t989 + t1087 * t963, t1083 * t932 + t1087 * t888, t1083 * t978 + t1087 * t926, t1083 * t988 + t1087 * t962, t1083 * t979 + t1087 * t927, t1083 * t1009 - t1087 * t1140, pkin(7) * t848 + t1083 * t789 + t1087 * t780, pkin(7) * t837 + t1083 * t781 + t1087 * t773, pkin(7) * t829 + t1083 * t739 + t1087 * t753, pkin(7) * t759 + t1083 * t743 + t1087 * t760, t1284, -t1420, t1370, t1345, t1402, t1312, pkin(7) * t674 + t1083 * t613 + t1087 * t603, pkin(7) * t680 + t1083 * t619 + t1087 * t610, pkin(7) * t648 + t1083 * t596 + t1087 * t585, pkin(7) * t605 + t1083 * t570 + t1087 * t569, t1284, t1370, t1420, t1312, -t1402, t1345, t1083 * t593 + t1087 * t579 + t1434, pkin(7) * t650 + t1083 * t584 + t1087 * t565, t1083 * t590 + t1087 * t573 + t1428, pkin(7) * t572 + t1083 * t549 + t1087 * t547, t1345, -t1420, t1402, t1087 * t767 + t1329, t1370, t1083 * t894 + t1087 * t831, t1083 * t574 + t1087 * t562 + t1428, t1083 * t578 + t1087 * t563 - t1434, pkin(7) * t649 + t1083 * t568 + t1087 * t553, pkin(7) * t559 + t1083 * t540 + t1087 * t535, t1083 * t691 + t1087 * t644, t1083 * t640 + t1087 * t612, t1083 * t663 + t1087 * t627, t1083 * t690 + t1087 * t643, t1083 * t664 + t1087 * t628, t1083 * t731 + t1087 * t683, pkin(7) * t592 + t1083 * t543 + t1087 * t537, pkin(7) * t595 + t1083 * t544 + t1087 * t538, pkin(7) * t583 + t1083 * t536 + t1087 * t532, pkin(7) * t534 + t1083 * t529 + t1087 * t528; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1066, t1067, 0, 0, (t1078 * t1041 + t1086 * t1306) * t1082, t1079 * t1047 + t1078 * t1149, t1079 * t1001 + t1078 * t1144, (t1078 * t1042 - t1082 * t1306) * t1086, t1079 * t1002 + t1078 * t1145, t1079 * t1192, pkin(1) * t928 + t1078 * t1182 - t1079 * t967, pkin(1) * t908 + t1078 * t1183 - t1079 * t968, pkin(1) * t887 + t1078 * t1111, pkin(1) * t835 + t1267 * t860, t1128, -t1404, t1348, t1316, t1393, t1282, pkin(1) * t730 + t1078 * t1113 + t1079 * t740, pkin(1) * t735 + t1078 * t1112 + t1079 * t750, pkin(1) * t699 + t1078 * t1109 + t1079 * t656, pkin(1) * t637 + t1078 * t1106 + t1079 * t686, t1128, t1348, t1404, t1282, -t1393, t1316, t1078 * t1114 + t1079 * t652 - t1431, pkin(1) * t701 + t1078 * t1117 + t1079 * t629, t1078 * t1115 + t1079 * t651 + t1425, pkin(1) * t599 + t1078 * t1119 + t1079 * t600, t1316, -t1404, t1393, t1078 * t1134 + t1303, t1348, t1218 * t937 + t1290, t1078 * t1118 + t1079 * t625 + t1425, t1078 * t1116 + t1079 * t626 + t1431, pkin(1) * t700 + t1078 * t1120 + t1079 * t597, pkin(1) * t581 + t1078 * t1121 + t1079 * t561, t1078 * t1170 + t1079 * t755, t1078 * t1177 + t1079 * t693, t1078 * t1173 + t1079 * t712, t1078 * t1171 + t1079 * t754, t1078 * t1172 + t1079 * t713, t1078 * t1169 + t1079 * t803, pkin(1) * t622 + t1078 * t1123 + t1079 * t554, pkin(1) * t624 + t1078 * t1122 + t1079 * t555, pkin(1) * t609 + t1078 * t1124 + t1079 * t548, pkin(1) * t542 + t1078 * t1125 + t1079 * t530;];
tauB_reg  = t1;