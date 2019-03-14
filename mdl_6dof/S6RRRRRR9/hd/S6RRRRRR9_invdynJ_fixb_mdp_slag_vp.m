% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR9_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR9_invdynJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:40:26
% EndTime: 2019-03-10 05:41:31
% DurationCPUTime: 46.28s
% Computational Cost: add. (38126->1075), mult. (106826->1475), div. (0->0), fcn. (92002->18), ass. (0->428)
t1114 = sin(pkin(6));
t1115 = cos(pkin(7));
t1124 = cos(qJ(2));
t1347 = cos(qJ(3));
t1252 = t1347 * t1124;
t1119 = sin(qJ(3));
t1120 = sin(qJ(2));
t1275 = t1119 * t1120;
t1157 = -t1115 * t1275 + t1252;
t1029 = t1157 * t1114;
t1017 = qJD(1) * t1029;
t1113 = sin(pkin(7));
t1258 = t1113 * t1347;
t1392 = qJD(3) * t1258 - t1017;
t1340 = pkin(10) * t1113;
t1173 = pkin(2) * t1120 - t1124 * t1340;
t1297 = qJD(1) * t1114;
t1033 = t1173 * t1297;
t1279 = t1115 * t1119;
t1283 = t1113 * t1119;
t1095 = pkin(10) * t1283;
t1256 = t1115 * t1347;
t1354 = pkin(2) * t1256 - t1095;
t1337 = cos(pkin(6));
t1262 = pkin(1) * t1337;
t1101 = t1124 * t1262;
t1091 = qJD(1) * t1101;
t1339 = pkin(10) * t1115;
t1189 = t1114 * (-pkin(9) - t1339);
t1174 = t1120 * t1189;
t994 = qJD(1) * t1174 + t1091;
t1100 = t1120 * t1262;
t1134 = t1124 * t1189 - t1100;
t995 = t1134 * qJD(1);
t1391 = t1354 * qJD(3) - t1033 * t1283 - t995 * t1279 - t1347 * t994;
t1253 = t1347 * t1120;
t1274 = t1119 * t1124;
t1154 = t1115 * t1253 + t1274;
t1028 = t1154 * t1114;
t1016 = qJD(1) * t1028;
t922 = t1115 * t1033 - t1113 * t995;
t1390 = -pkin(3) * t1016 + pkin(11) * t1017 + (pkin(3) * t1119 - pkin(11) * t1347) * t1113 * qJD(3) - t922;
t1282 = t1114 * t1120;
t1240 = t1113 * t1282;
t1201 = qJD(1) * t1240;
t1389 = pkin(11) * t1201 - t1391;
t1118 = sin(qJ(4));
t1123 = cos(qJ(4));
t1045 = -t1123 * t1115 + t1118 * t1283;
t1198 = t1118 * t1240;
t1302 = -qJD(1) * t1198 - qJD(4) * t1045 + t1123 * t1392;
t1046 = t1115 * t1118 + t1123 * t1283;
t1197 = t1123 * t1240;
t1301 = qJD(1) * t1197 + qJD(4) * t1046 + t1118 * t1392;
t1230 = t1337 * qJD(1);
t1184 = t1230 + qJD(2);
t1162 = t1184 * t1113;
t1178 = t1252 * t1297;
t1236 = t1114 * t1275;
t1200 = qJD(1) * t1236;
t1363 = t1115 * t1178 + t1347 * t1162 - t1200;
t960 = qJD(4) - t1363;
t1376 = -qJD(3) * t1283 + t1016;
t1268 = pkin(2) * t1279 + pkin(10) * t1258;
t1388 = t1268 * qJD(3) - t1119 * t994 + t1256 * t995;
t1248 = qJD(1) * t1282;
t1255 = t1347 * t1033;
t1299 = -(-pkin(3) * t1248 - t1255) * t1113 + t1388;
t1038 = (-pkin(3) * t1347 - pkin(11) * t1119 - pkin(2)) * t1113;
t1288 = t1038 * t1123;
t1387 = -qJD(4) * t1288 - t1118 * t1390 + t1389 * t1123;
t1112 = qJ(5) + qJ(6);
t1107 = sin(t1112);
t1108 = cos(t1112);
t1125 = cos(qJ(1));
t1231 = t1125 * t1337;
t1346 = sin(qJ(1));
t1047 = t1120 * t1346 - t1124 * t1231;
t1048 = t1120 * t1231 + t1124 * t1346;
t1280 = t1114 * t1125;
t1238 = t1113 * t1280;
t936 = -t1047 * t1279 + t1048 * t1347 - t1119 * t1238;
t998 = -t1047 * t1113 + t1115 * t1280;
t879 = t1118 * t998 - t1123 * t936;
t935 = t1047 * t1256 + t1048 * t1119 + t1238 * t1347;
t1386 = t1107 * t879 + t1108 * t935;
t1385 = -t1107 * t935 + t1108 * t879;
t1117 = sin(qJ(5));
t1122 = cos(qJ(5));
t1384 = t1117 * t879 + t1122 * t935;
t1383 = -t1117 * t935 + t1122 * t879;
t1037 = pkin(11) * t1115 + t1268;
t1271 = t1123 * t1037 + t1118 * t1038;
t1382 = -qJD(4) * t1271 + t1389 * t1118 + t1123 * t1390;
t1021 = t1118 * t1037;
t1381 = pkin(12) * t1376 + qJD(4) * t1021 + t1387;
t1380 = -t1301 * pkin(4) + pkin(12) * t1302 - t1299;
t1281 = t1114 * t1124;
t1147 = pkin(9) * t1281 + t1100;
t1237 = t1115 * t1281;
t956 = t1147 * qJD(1) + (qJD(1) * t1237 + t1162) * pkin(10);
t1138 = pkin(2) * t1337 + t1174;
t961 = qJD(2) * pkin(2) + qJD(1) * t1138 + t1091;
t1170 = -pkin(2) * t1124 - t1120 * t1340 - pkin(1);
t1025 = t1170 * t1114;
t1012 = qJD(1) * t1025;
t978 = t1012 * t1283;
t848 = t961 * t1279 + t1347 * t956 + t978;
t1378 = -t848 + t960 * (pkin(4) * t1118 - pkin(12) * t1123);
t1277 = t1117 * t1123;
t1156 = t1115 * t1274 + t1253;
t1140 = t1156 * t1114;
t1150 = t1119 * t1162;
t964 = qJD(1) * t1140 + t1150;
t883 = -t1122 * t964 + t1277 * t1363;
t1377 = qJD(4) * t1277 - t883;
t1109 = t1114 ^ 2;
t1375 = 0.2e1 * pkin(1) * t1109;
t1121 = cos(qJ(6));
t1116 = sin(qJ(6));
t1239 = t1113 * t1281;
t1139 = qJD(1) * t1239 - t1115 * t1184 - qJD(3);
t909 = -t1118 * t1139 + t1123 * t964;
t856 = t1117 * t960 + t1122 * t909;
t1309 = t1116 * t856;
t854 = t1117 * t909 - t1122 * t960;
t778 = t1121 * t854 + t1309;
t1011 = t1123 * t1139;
t907 = t1118 * t964 + t1011;
t905 = qJD(5) + t907;
t904 = qJD(6) + t905;
t1374 = t778 * t904;
t1155 = -t1117 * t1046 - t1122 * t1258;
t1321 = qJD(5) * t1155 - t1376 * t1117 + t1122 * t1302;
t1004 = t1122 * t1046 - t1117 * t1258;
t1320 = qJD(5) * t1004 + t1117 * t1302 + t1376 * t1122;
t1180 = t1116 * t854 - t1121 * t856;
t1368 = t1180 * t904;
t1272 = t1122 * t1123;
t884 = t1117 * t964 + t1272 * t1363;
t1367 = qJD(4) * t1272 - t884;
t1355 = qJD(5) + qJD(6);
t1366 = t1355 + t907;
t1276 = t1118 * t1122;
t1364 = qJD(5) * t1276 + t1377;
t1290 = qJD(6) * t1116;
t902 = t1115 * t1012 - t1113 * t961;
t822 = -pkin(3) * t1363 - pkin(11) * t964 + t902;
t830 = -pkin(11) * t1139 + t848;
t761 = t1118 * t822 + t1123 * t830;
t753 = pkin(12) * t960 + t761;
t1207 = t1012 * t1258;
t847 = -t1119 * t956 + t1256 * t961 + t1207;
t829 = pkin(3) * t1139 - t847;
t767 = t907 * pkin(4) - t909 * pkin(12) + t829;
t725 = t1117 * t767 + t1122 * t753;
t719 = -pkin(13) * t854 + t725;
t717 = t719 * t1290;
t760 = -t1118 * t830 + t1123 * t822;
t752 = -pkin(4) * t960 - t760;
t736 = pkin(5) * t854 + t752;
t1202 = t1337 * t1346;
t1142 = t1125 * t1120 + t1124 * t1202;
t1257 = t1114 * t1346;
t1000 = t1113 * t1142 + t1115 * t1257;
t1049 = -t1120 * t1202 + t1124 * t1125;
t1203 = t1113 * t1257;
t940 = t1049 * t1347 + (-t1115 * t1142 + t1203) * t1119;
t881 = t1000 * t1118 + t1123 * t940;
t1131 = t1142 * t1347;
t939 = t1049 * t1119 + t1115 * t1131 - t1203 * t1347;
t818 = t1107 * t939 + t1108 * t881;
t1044 = -t1115 * t1337 + t1239;
t1229 = t1337 * t1113;
t1193 = t1119 * t1229;
t987 = t1193 + t1140;
t930 = -t1044 * t1118 + t1123 * t987;
t1177 = t1347 * t1229;
t1204 = t1115 * t1252;
t1179 = t1114 * t1204;
t986 = -t1177 - t1179 + t1236;
t1362 = t736 * t778 + g(1) * t818 - g(2) * t1385 - g(3) * (-t1107 * t986 - t1108 * t930) + t717;
t817 = -t1107 * t881 + t1108 * t939;
t1361 = t736 * t1180 - g(3) * (-t1107 * t930 + t1108 * t986) - g(2) * t1386 - g(1) * t817;
t1225 = t1337 * qJDD(1);
t1094 = t1225 + qJDD(2);
t1127 = qJD(2) * t1154 + qJD(3) * t1156;
t1265 = qJDD(1) * t1120;
t868 = qJD(3) * t1150 - qJDD(1) * t1179 - t1094 * t1258 + t1114 * (qJD(1) * t1127 + t1119 * t1265);
t1360 = t1118 * t936 + t1123 * t998;
t1264 = qJDD(1) * t1124;
t1227 = t1114 * t1264;
t1228 = t1114 * t1265;
t867 = t1094 * t1283 + t1227 * t1279 + t1347 * t1228 + (-t1115 * t1200 + t1178) * qJD(2) + t1363 * qJD(3);
t1247 = qJD(2) * t1282;
t1196 = qJD(1) * t1247;
t969 = t1094 * t1115 + qJDD(3) + (t1196 - t1227) * t1113;
t1218 = t1118 * t867 - t1123 * t969;
t793 = qJD(4) * t909 + t1218;
t791 = qJDD(5) + t793;
t788 = qJDD(6) + t791;
t1359 = t788 * MDP(36) + (t1180 ^ 2 - t778 ^ 2) * MDP(33) - t778 * MDP(32) * t1180;
t1322 = pkin(4) * t1376 - t1382;
t1358 = t1380 * t1122;
t1357 = t1118 * t960;
t1278 = t1117 * t1118;
t1338 = pkin(11) * qJD(4);
t1356 = t1122 * t1378 + t1278 * t1338;
t1066 = t1116 * t1122 + t1117 * t1121;
t1042 = t1066 * t1118;
t1291 = qJD(5) * t1122;
t1292 = qJD(5) * t1117;
t1036 = t1095 + (-pkin(2) * t1347 - pkin(3)) * t1115;
t931 = t1045 * pkin(4) - t1046 * pkin(12) + t1036;
t934 = -pkin(12) * t1258 + t1271;
t1353 = t1380 * t1117 + t1122 * t1381 - t931 * t1291 + t1292 * t934;
t1298 = qJD(1) * qJD(2);
t1352 = -t1120 * t1298 + t1264;
t1080 = -pkin(4) * t1123 - pkin(12) * t1118 - pkin(3);
t1064 = t1122 * t1080;
t889 = pkin(3) * t964 - pkin(11) * t1363;
t1325 = t1118 * t889 + t1123 * t847;
t776 = pkin(12) * t964 + t1325;
t1351 = -qJD(5) * t1064 - t1117 * t1378 + t1122 * t776;
t1294 = qJD(4) * t1118;
t792 = -qJD(4) * t1011 + t1118 * t969 + t1123 * t867 - t1294 * t964;
t864 = qJDD(4) + t868;
t741 = qJD(5) * t856 + t1117 * t792 - t1122 * t864;
t740 = t1117 * t864 + t1122 * t792 + t960 * t1291 - t1292 * t909;
t1222 = t1116 * t740 + t1121 * t741;
t704 = -qJD(6) * t1180 + t1222;
t1126 = qJD(1) ^ 2;
t1348 = pkin(12) + pkin(13);
t1341 = pkin(5) * t1118;
t1336 = t854 * t905;
t1335 = t856 * t905;
t1334 = t907 * t960;
t1333 = t909 * t960;
t834 = pkin(4) * t909 + pkin(12) * t907;
t1332 = t1117 * t834 + t1122 * t760;
t993 = t1101 + t1138;
t918 = t1115 * t1025 - t1113 * t993;
t843 = pkin(3) * t986 - pkin(11) * t987 + t918;
t1005 = t1025 * t1283;
t979 = (t1229 + t1237) * pkin(10) + t1147;
t1261 = t993 * t1279 + t1347 * t979 + t1005;
t851 = -pkin(11) * t1044 + t1261;
t1324 = t1118 * t843 + t1123 * t851;
t770 = pkin(12) * t986 + t1324;
t1206 = t1025 * t1258;
t1136 = -t1119 * t979 + t1256 * t993 + t1206;
t850 = t1044 * pkin(3) - t1136;
t929 = t1044 * t1123 + t1118 * t987;
t785 = t929 * pkin(4) - t930 * pkin(12) + t850;
t1331 = t1117 * t785 + t1122 * t770;
t912 = t1004 * t1116 - t1121 * t1155;
t1328 = -qJD(6) * t912 - t1116 * t1320 + t1121 * t1321;
t913 = t1004 * t1121 + t1116 * t1155;
t1327 = qJD(6) * t913 + t1116 * t1321 + t1121 * t1320;
t1326 = pkin(5) * t1320 + t1322;
t1065 = t1116 * t1117 - t1121 * t1122;
t1293 = qJD(4) * t1123;
t1319 = -t1042 * t1355 - t1065 * t1293 + t1116 * t883 - t1121 * t884;
t1318 = (t1276 * t1355 + t1377) * t1121 + (-t1278 * t1355 + t1367) * t1116;
t1317 = t1117 * t931 + t1122 * t934;
t1316 = t1366 * t1065;
t1315 = t1366 * t1066;
t1313 = qJD(4) * t960;
t1312 = qJD(5) * t905;
t1311 = t1094 * MDP(8);
t1208 = qJD(3) * t1256;
t1296 = qJD(3) * t1119;
t1081 = t1091 * qJD(2);
t1210 = pkin(1) * t1225;
t1232 = pkin(9) * t1227 + t1120 * t1210 + t1081;
t1141 = -pkin(9) * t1196 + t1232;
t897 = (t1114 * t1115 * t1352 + t1094 * t1113) * pkin(10) + t1141;
t1148 = -t1100 * t1298 + t1124 * t1210;
t1249 = t1124 * t1298;
t1161 = -t1249 - t1265;
t1149 = t1161 * pkin(9);
t903 = t1094 * pkin(2) + (t1161 * t1339 + t1149) * t1114 + t1148;
t1158 = t1173 * qJD(2);
t941 = (qJD(1) * t1158 + qJDD(1) * t1170) * t1114;
t1144 = -qJD(3) * t1207 - t961 * t1208 - t903 * t1279 - t941 * t1283 + t1296 * t956 - t1347 * t897;
t748 = pkin(11) * t969 - t1144;
t837 = -t1113 * t903 + t1115 * t941;
t754 = pkin(3) * t868 - pkin(11) * t867 + t837;
t708 = t1118 * t754 + t1123 * t748 + t822 * t1293 - t1294 * t830;
t706 = pkin(12) * t864 + t708;
t1244 = qJD(3) * t1279;
t1259 = qJD(3) * t1347;
t756 = -qJD(3) * t978 - t1119 * t897 - t961 * t1244 + t903 * t1256 + t941 * t1258 - t956 * t1259;
t749 = -pkin(3) * t969 - t756;
t714 = pkin(4) * t793 - pkin(12) * t792 + t749;
t697 = -qJD(5) * t725 - t1117 * t706 + t1122 * t714;
t694 = pkin(5) * t791 - pkin(13) * t740 + t697;
t1310 = t1116 * t694;
t1308 = t1117 * t740;
t1307 = t1117 * t791;
t1306 = t1117 * t907;
t724 = -t1117 * t753 + t1122 * t767;
t718 = -pkin(13) * t856 + t724;
t711 = pkin(5) * t905 + t718;
t1305 = t1121 * t711;
t1304 = t1121 * t719;
t1303 = t1122 * t791;
t1214 = t1122 * t905;
t841 = t1118 * t847;
t775 = -pkin(4) * t964 - t1123 * t889 + t841;
t1300 = pkin(5) * t1364 + pkin(11) * t1293 - t775;
t1295 = qJD(4) * t1117;
t1289 = qJD(6) * t1121;
t1287 = t1107 * t1123;
t1286 = t1108 * t1123;
t1285 = t1109 * t1126;
t1284 = t1113 * t1118;
t1102 = pkin(11) * t1272;
t1267 = t1117 * t1080 + t1102;
t1110 = t1120 ^ 2;
t1266 = -t1124 ^ 2 + t1110;
t1263 = -t1116 * t741 + t1121 * t740 - t854 * t1289;
t1085 = t1348 * t1117;
t1034 = t1114 * t1158;
t1254 = t1347 * t1034;
t1251 = qJD(4) * t1214;
t1246 = qJD(2) * t1281;
t1241 = t1124 * t1285;
t1165 = -t1117 * t714 - t1122 * t706 - t767 * t1291 + t1292 * t753;
t695 = -pkin(13) * t741 - t1165;
t1224 = qJD(6) * t711 + t695;
t1223 = -t1116 * t695 + t1121 * t694;
t1221 = -t1117 * t770 + t1122 * t785;
t1219 = -t1117 * t934 + t1122 * t931;
t1217 = -t1118 * t851 + t1123 * t843;
t1213 = t1123 * t960;
t997 = t1134 * qJD(2);
t923 = t1115 * t1034 - t1113 * t997;
t709 = -t1118 * t748 + t1123 * t754 - t830 * t1293 - t822 * t1294;
t1211 = -t1021 + t1288;
t1199 = qJD(2) * t1240;
t1195 = -t761 + (t1292 + t1306) * pkin(5);
t1194 = t1114 * t1126 * t1337;
t831 = pkin(13) * t1155 + t1317;
t1192 = -pkin(5) * t1301 + pkin(13) * t1321 + qJD(5) * t1317 + qJD(6) * t831 - t1117 * t1381 + t1358;
t985 = -pkin(13) * t1276 + t1064 + (-pkin(11) * t1117 - pkin(5)) * t1123;
t1191 = -qJD(6) * t985 - (-qJD(4) * t1276 - qJD(5) * t1277) * pkin(11) + t1351 + t1364 * pkin(13);
t811 = pkin(5) * t1045 - pkin(13) * t1004 + t1219;
t1190 = pkin(13) * t1320 - qJD(6) * t811 + t1353;
t1007 = -pkin(13) * t1278 + t1267;
t1187 = -pkin(13) * t884 + qJD(6) * t1007 - t1117 * t776 + t1341 * t1363 - (-pkin(13) * t1272 + t1341) * qJD(4) - (-t1102 + (pkin(13) * t1118 - t1080) * t1117) * qJD(5) - t1356;
t1186 = pkin(13) * t1306 + t1085 * t1355 + t1332;
t1086 = t1348 * t1122;
t833 = t1122 * t834;
t1185 = pkin(13) * t1122 * t907 + pkin(5) * t909 + qJD(6) * t1086 - t1117 * t760 + t1348 * t1291 + t833;
t933 = pkin(4) * t1258 - t1211;
t1183 = 0.2e1 * t1230 + qJD(2);
t702 = t1116 * t711 + t1304;
t876 = t1117 * t986 + t1122 * t930;
t723 = pkin(5) * t929 - pkin(13) * t876 + t1221;
t875 = t1117 * t930 - t986 * t1122;
t726 = -pkin(13) * t875 + t1331;
t1182 = -t1116 * t726 + t1121 * t723;
t1181 = t1116 * t723 + t1121 * t726;
t804 = t1116 * t876 + t1121 * t875;
t805 = -t1116 * t875 + t1121 * t876;
t1092 = qJD(2) * t1101;
t996 = qJD(2) * t1174 + t1092;
t1143 = qJD(3) * t1206 + t1034 * t1283 + t993 * t1208 + t997 * t1279 - t1296 * t979 + t1347 * t996;
t796 = pkin(11) * t1199 + t1143;
t920 = qJD(3) * t1193 + t1114 * t1127;
t921 = qJD(3) * t1177 + ((t1204 - t1275) * qJD(3) + t1157 * qJD(2)) * t1114;
t807 = pkin(3) * t920 - pkin(11) * t921 + t923;
t1176 = -t1118 * t796 + t1123 * t807 - t851 * t1293 - t843 * t1294;
t769 = -pkin(4) * t986 - t1217;
t707 = -pkin(4) * t864 - t709;
t1172 = -pkin(11) * t864 + t829 * t960;
t1171 = -pkin(12) * t791 + t752 * t905;
t880 = t1000 * t1123 - t1118 * t940;
t1169 = g(1) * t880 - g(2) * t1360 - g(3) * t929;
t1168 = g(1) * t939 + g(2) * t935 + g(3) * t986;
t1167 = g(1) * t940 + g(2) * t936 + g(3) * t987;
t1166 = t1118 * t807 + t1123 * t796 + t843 * t1293 - t1294 * t851;
t721 = pkin(12) * t920 + t1166;
t1145 = -qJD(3) * t1005 - t1119 * t996 - t993 * t1244 + t1256 * t997 - t979 * t1259;
t797 = (-pkin(3) * t1247 - t1254) * t1113 - t1145;
t835 = -qJD(2) * t1197 + qJD(4) * t930 + t1118 * t921;
t836 = qJD(2) * t1198 - qJD(4) * t929 + t1123 * t921;
t735 = t835 * pkin(4) - t836 * pkin(12) + t797;
t1164 = t1117 * t735 + t1122 * t721 + t785 * t1291 - t1292 * t770;
t703 = -t1290 * t856 + t1263;
t953 = -t1047 * t1347 - t1048 * t1279;
t955 = -t1049 * t1279 - t1131;
t1160 = -g(1) * t955 - g(2) * t953 - g(3) * t1029;
t1159 = t1114 * (t1225 + t1094);
t1152 = t1168 - t749;
t722 = -pkin(4) * t920 - t1176;
t1146 = g(1) * t1049 + g(2) * t1048 + g(3) * t1282;
t1137 = pkin(11) * t1312 - t1168;
t1135 = pkin(12) * t1312 + t1169 + t707;
t692 = -qJD(6) * t702 + t1223;
t1133 = -qJD(5) * t1331 - t1117 * t721 + t1122 * t735;
t1132 = qJD(3) * t1139;
t1129 = t1139 * t1240;
t1128 = t1184 * t1147;
t1106 = -pkin(5) * t1122 - pkin(4);
t1073 = (pkin(5) * t1117 + pkin(11)) * t1118;
t1043 = t1065 * t1118;
t965 = t1029 * t1123 + t1198;
t954 = t1049 * t1256 - t1119 * t1142;
t952 = -t1047 * t1119 + t1048 * t1256;
t911 = t1049 * t1284 + t1123 * t955;
t910 = t1048 * t1284 + t1123 * t953;
t885 = -pkin(5) * t1155 + t933;
t826 = t1117 * t939 + t1122 * t881;
t825 = -t1117 * t881 + t1122 * t939;
t764 = -qJD(5) * t875 + t1117 * t920 + t1122 * t836;
t763 = qJD(5) * t876 + t1117 * t836 - t920 * t1122;
t747 = pkin(5) * t875 + t769;
t716 = qJD(6) * t805 + t1116 * t764 + t1121 * t763;
t715 = -qJD(6) * t804 - t1116 * t763 + t1121 * t764;
t710 = pkin(5) * t763 + t722;
t701 = -t1116 * t719 + t1305;
t700 = pkin(5) * t741 + t707;
t699 = -pkin(13) * t763 + t1164;
t698 = pkin(5) * t835 - pkin(13) * t764 + t1133;
t691 = t1224 * t1121 + t1310 - t717;
t1 = [(-t867 * t1044 - t1139 * t921 + t1199 * t964 + t987 * t969) * MDP(13) + (g(1) * t1346 - g(2) * t1125) * MDP(2) + (g(1) * t1125 + g(2) * t1346) * MDP(3) + (t1120 * t1159 + t1183 * t1246) * MDP(6) + (t1124 * t1159 - t1183 * t1247) * MDP(7) + (-g(1) * t879 - g(2) * t881 + t1176 * t960 + t1217 * t864 + t709 * t986 + t749 * t929 + t760 * t920 + t850 * t793 + t797 * t907 + t829 * t835) * MDP(23) + (-t1180 * t835 + t703 * t929 + t715 * t904 + t788 * t805) * MDP(34) + (t1180 * t716 - t703 * t804 - t704 * t805 - t715 * t778) * MDP(33) + (-t1180 * t715 + t703 * t805) * MDP(32) + (-g(1) * t935 + g(2) * t939 - t1044 * t1144 + t1139 * t1143 - t1199 * t848 - t1261 * t969 + t837 * t987 + t918 * t867 + t902 * t921 + t923 * t964) * MDP(17) + (-t704 * t929 - t716 * t904 - t778 * t835 - t788 * t804) * MDP(35) + (t788 * t929 + t835 * t904) * MDP(36) + (t791 * t929 + t835 * t905) * MDP(29) + (-t741 * t929 - t763 * t905 - t791 * t875 - t835 * t854) * MDP(28) + (t740 * t929 + t764 * t905 + t791 * t876 + t835 * t856) * MDP(27) + (-t792 * t929 - t793 * t930 - t835 * t909 - t836 * t907) * MDP(19) + (t792 * t930 + t836 * t909) * MDP(18) + (t867 * t987 + t921 * t964) * MDP(11) + (t740 * t876 + t764 * t856) * MDP(25) + (-t740 * t875 - t741 * t876 - t763 * t856 - t764 * t854) * MDP(26) + (-g(1) * t1383 - g(2) * t826 + t1133 * t905 + t1221 * t791 + t697 * t929 + t707 * t875 + t722 * t854 + t724 * t835 + t769 * t741 + t752 * t763) * MDP(30) + (g(1) * t1384 - g(2) * t825 - t1164 * t905 + t1165 * t929 - t1331 * t791 + t707 * t876 + t722 * t856 - t725 * t835 + t769 * t740 + t752 * t764) * MDP(31) + ((-qJD(6) * t1181 - t1116 * t699 + t1121 * t698) * t904 + t1182 * t788 + t692 * t929 + t701 * t835 + t710 * t778 + t747 * t704 + t700 * t804 + t736 * t716 - g(1) * t1385 - g(2) * t818) * MDP(37) + (-(qJD(6) * t1182 + t1116 * t698 + t1121 * t699) * t904 - t1181 * t788 - t691 * t929 - t702 * t835 - t710 * t1180 + t747 * t703 + t700 * t805 + t736 * t715 + g(1) * t1386 - g(2) * t817) * MDP(38) + (-t793 * t986 - t835 * t960 - t864 * t929 - t907 * t920) * MDP(21) + (t792 * t986 + t836 * t960 + t864 * t930 + t909 * t920) * MDP(20) + (t864 * t986 + t920 * t960) * MDP(22) + t1337 * t1311 + (-qJD(2) * t1129 - t969 * t1044) * MDP(15) + qJDD(1) * MDP(1) + (0.2e1 * (t1120 * t1264 - t1266 * t1298) * MDP(5) + (qJDD(1) * t1110 + 0.2e1 * t1120 * t1249) * MDP(4)) * t1109 + (-qJD(2) * t1128 + (-pkin(9) * t1282 + t1101) * t1094 + (t1114 * t1149 + t1148) * t1337 + g(1) * t1048 - g(2) * t1049 + t1352 * t1375) * MDP(9) + (-(-pkin(9) * t1247 + t1092) * t1184 - t1147 * t1094 - t1141 * t1337 - g(1) * t1047 + g(2) * t1142 + t1161 * t1375) * MDP(10) + (-g(1) * t1360 - g(2) * t880 - t1166 * t960 - t1324 * t864 - t708 * t986 + t749 * t930 - t761 * t920 + t850 * t792 + t797 * t909 + t829 * t836) * MDP(24) + (t868 * t1044 + t1139 * t920 + t1199 * t1363 - t986 * t969) * MDP(14) + (-(t1113 * t1254 + t1145) * t1139 + t1136 * t969 - t756 * t1044 + t847 * t1199 - t923 * t1363 + t918 * t868 + t837 * t986 + t902 * t920 + g(1) * t936 - g(2) * t940) * MDP(16) + (t1363 * t921 - t867 * t986 - t868 * t987 - t920 * t964) * MDP(12); (t1354 * t969 + t756 * t1115 - t1113 * pkin(2) * t868 - t837 * t1258 - t847 * t1201 + t922 * t1363 - t1376 * t902 + t1160 + (t1113 * t1255 + t1388) * t1139) * MDP(16) + ((-pkin(9) * t1248 + t1091) * t1230 + pkin(1) * t1241 + t1081 + t1146 - t1232) * MDP(10) + ((-t1116 * t831 + t1121 * t811) * t788 + t692 * t1045 + t885 * t704 + t700 * t912 - g(1) * (t1107 * t954 + t1108 * t911) - g(2) * (t1107 * t952 + t1108 * t910) - g(3) * (t1028 * t1107 + t1108 * t965) + (t1116 * t1190 - t1121 * t1192) * t904 + t1326 * t778 + t1327 * t736 + t1301 * t701) * MDP(37) + (-t1045 * t704 - t1301 * t778 - t1327 * t904 - t788 * t912) * MDP(35) + (t1004 * t791 + t1045 * t740 + t1301 * t856 + t1321 * t905) * MDP(27) + (t1004 * t740 + t1321 * t856) * MDP(25) + t1266 * MDP(5) * t1285 + (-t1268 * t969 + t1144 * t1115 - t922 * t964 - t902 * t1017 + g(1) * t954 + g(2) * t952 + g(3) * t1028 + (-pkin(2) * t867 + t837 * t1119 + t1248 * t848 + t1259 * t902) * t1113 + t1391 * t1139) * MDP(17) + (-t1317 * t791 + t1165 * t1045 + t933 * t740 + t707 * t1004 - g(1) * (-t1117 * t911 + t1122 * t954) - g(2) * (-t1117 * t910 + t1122 * t952) - g(3) * (t1028 * t1122 - t1117 * t965) + t1353 * t905 + t1322 * t856 + t1321 * t752 - t1301 * t725) * MDP(31) + (t1120 * t1194 + t1227) * MDP(7) + (-t1124 * t1194 + t1228) * MDP(6) + (-t1180 * t1328 + t703 * t913) * MDP(32) + (-(t1116 * t811 + t1121 * t831) * t788 - t691 * t1045 + t885 * t703 + t700 * t913 - g(1) * (-t1107 * t911 + t1108 * t954) - g(2) * (-t1107 * t910 + t1108 * t952) - g(3) * (t1028 * t1108 - t1107 * t965) + (t1116 * t1192 + t1121 * t1190) * t904 - t1326 * t1180 + t1328 * t736 - t1301 * t702) * MDP(38) + (t1045 * t703 - t1180 * t1301 + t1328 * t904 + t788 * t913) * MDP(34) + (t1180 * t1327 - t1328 * t778 - t703 * t912 - t704 * t913) * MDP(33) + (-t1045 * t741 + t1155 * t791 - t1301 * t854 - t1320 * t905) * MDP(28) + (-t1004 * t741 + t1155 * t740 - t1320 * t856 - t1321 * t854) * MDP(26) + (t1045 * t791 + t1301 * t905) * MDP(29) + (t1045 * t788 + t1301 * t904) * MDP(36) + (t1046 * t792 + t1302 * t909) * MDP(18) + (-t1045 * t792 - t1046 * t793 - t1301 * t909 - t1302 * t907) * MDP(19) + (t1211 * t864 + t1036 * t793 + t749 * t1045 - t760 * t1016 - g(1) * t911 - g(2) * t910 - g(3) * t965 + t1382 * t960 + t1299 * t907 + t1301 * t829 + (t1296 * t760 - t1347 * t709) * t1113) * MDP(23) + (-t1271 * t864 + t1036 * t792 + t749 * t1046 + t761 * t1016 + t1387 * t960 + t1299 * t909 + t1302 * t829 + (t1037 * t1313 - t1160) * t1118 + (-t1123 * t1146 - t1296 * t761 + t1347 * t708) * t1113) * MDP(24) + t1311 - t1120 * MDP(4) * t1241 + (-t960 * t1016 + (t1296 * t960 - t1347 * t864) * t1113) * MDP(22) + (t907 * t1016 - t1045 * t864 - t1301 * t960 + (-t1296 * t907 + t1347 * t793) * t1113) * MDP(21) + (-t1016 * t909 + t864 * t1046 + t1302 * t960 + (t1296 * t909 - t1347 * t792) * t1113) * MDP(20) + (t1219 * t791 + t697 * t1045 + t933 * t741 - t707 * t1155 - g(1) * (t1117 * t954 + t1122 * t911) - g(2) * (t1117 * t952 + t1122 * t910) - g(3) * (t1028 * t1117 + t1122 * t965) + (-t934 * t1291 + (-qJD(5) * t931 + t1381) * t1117 - t1358) * t905 + t1322 * t854 + t1320 * t752 + t1301 * t724) * MDP(30) + (pkin(1) * t1120 * t1285 + g(1) * t1142 + g(2) * t1047 - g(3) * t1281 + qJD(1) * t1128 + t1148 + (-qJD(1) * t1246 - t1228) * pkin(9)) * MDP(9) + (qJD(1) * t1129 + t969 * t1115) * MDP(15) + (t867 * t1115 + t1017 * t1139 + (t1119 * t969 - t1132 * t1347 - t1248 * t964) * t1113) * MDP(13) + (-t964 * t1017 + (t1119 * t867 + t1259 * t964) * t1113) * MDP(11) + (t964 * t1016 - t1017 * t1363 + (t1347 * t867 - t1119 * t868 + (-t1119 * t964 + t1347 * t1363) * qJD(3)) * t1113) * MDP(12) + (-t868 * t1115 - t1016 * t1139 + (t1119 * t1132 - t1248 * t1363 + t1347 * t969) * t1113) * MDP(14); (-pkin(3) * t793 + t841 * t960 - t848 * t907 + t1172 * t1118 + ((-t889 - t1338) * t960 + t1152) * t1123) * MDP(23) + (-pkin(3) * t792 + t1325 * t960 - t848 * t909 + t1172 * t1123 + (pkin(11) * t1313 - t1152) * t1118) * MDP(24) + (t1118 * t792 + t1213 * t909) * MDP(18) + (t1064 * t791 - t752 * t883 - t775 * t854 + t1356 * t905 + ((-qJD(5) * t1080 + t776) * t905 - t1167) * t1117 + (t752 * t1295 - t697 + (qJD(4) * t854 - t1307) * pkin(11) - t1137 * t1122) * t1123 + (pkin(11) * t741 + t707 * t1117 + t1291 * t752 + t724 * t960) * t1118) * MDP(30) + (-t884 * t905 + (-t740 + t1251) * t1123 + (-t1292 * t905 + t856 * t960 + t1303) * t1118) * MDP(27) + (-t1267 * t791 - t775 * t856 - t752 * t884 + t1351 * t905 - t1167 * t1122 + (-t1165 + (pkin(11) * t856 + t1122 * t752) * qJD(4) + t1137 * t1117) * t1123 + (-t752 * t1292 + t707 * t1122 - t960 * t725 + (t740 + t1251) * pkin(11)) * t1118) * MDP(31) + (-t1043 * t703 - t1180 * t1319) * MDP(32) + (-t1042 * t703 + t1043 * t704 + t1180 * t1318 - t1319 * t778) * MDP(33) + ((-t1007 * t1116 + t1121 * t985) * t788 - t692 * t1123 + t1073 * t704 + t700 * t1042 - g(1) * (t1107 * t940 - t1286 * t939) - g(2) * (t1107 * t936 - t1286 * t935) - g(3) * (t1107 * t987 - t1286 * t986) + (t1116 * t1191 - t1121 * t1187) * t904 + t1300 * t778 + t1318 * t736 + t701 * t1357) * MDP(37) + (-t1123 * t791 + t1357 * t905) * MDP(29) + (-t1123 * t788 + t1357 * t904) * MDP(36) + (-t1043 * t788 - t1123 * t703 - t1180 * t1357 + t1319 * t904) * MDP(34) + (-(t1007 * t1121 + t1116 * t985) * t788 + t691 * t1123 + t1073 * t703 - t700 * t1043 - g(1) * (t1108 * t940 + t1287 * t939) - g(2) * (t1108 * t936 + t1287 * t935) - g(3) * (t1108 * t987 + t1287 * t986) + (t1116 * t1187 + t1121 * t1191) * t904 - t1300 * t1180 + t1319 * t736 - t702 * t1357) * MDP(38) + (-t1042 * t788 + t1123 * t704 - t1318 * t904 - t1357 * t778) * MDP(35) + (t1123 * t864 - t1357 * t960) * MDP(21) + (t1118 * t864 + t1213 * t960) * MDP(20) + ((t792 - t1334) * t1123 + (-t793 - t1333) * t1118) * MDP(19) + t969 * MDP(15) + (t883 * t905 + (-t1295 * t905 + t741) * t1123 + (-t1291 * t905 - t854 * t960 - t1307) * t1118) * MDP(28) + (t854 * t884 + t856 * t883 + (-t1117 * t856 - t1122 * t854) * t1293 + (-t1308 - t1122 * t741 + (t1117 * t854 - t1122 * t856) * qJD(5)) * t1118) * MDP(26) + (-t1139 * t848 + t1168 + t756) * MDP(16) - t868 * MDP(14) + (t740 * t1276 + (-qJD(5) * t1278 + t1367) * t856) * MDP(25) - t1363 ^ 2 * MDP(12) + (-MDP(11) * t1363 + MDP(12) * t964 - MDP(14) * t1139 - MDP(16) * t902 - MDP(20) * t909 + MDP(21) * t907 - MDP(22) * t960 - MDP(23) * t760 + MDP(24) * t761) * t964 + (-t1139 * t847 - t1363 * t902 + t1144 + t1167) * MDP(17) + (t1139 * t1363 + t867) * MDP(13); -t907 ^ 2 * MDP(19) + (t792 + t1334) * MDP(20) + (-t1218 + t1333) * MDP(21) + t864 * MDP(22) + (t761 * t960 - t1169 + t709) * MDP(23) + (g(1) * t881 - g(2) * t879 + g(3) * t930 + t760 * t960 + t829 * t907 - t708) * MDP(24) + (t1214 * t856 + t1308) * MDP(25) + ((t740 - t1336) * t1122 + (-t741 - t1335) * t1117) * MDP(26) + (t1214 * t905 + t1307) * MDP(27) + (-t1117 * t905 ^ 2 + t1303) * MDP(28) + (-pkin(4) * t741 - t761 * t854 - t833 * t905 + (t760 * t905 + t1171) * t1117 - t1135 * t1122) * MDP(30) + (-pkin(4) * t740 + t1117 * t1135 + t1122 * t1171 + t1332 * t905 - t761 * t856) * MDP(31) + (t1066 * t703 + t1180 * t1316) * MDP(32) + (-t1065 * t703 - t1066 * t704 + t1180 * t1315 + t1316 * t778) * MDP(33) + (t1066 * t788 - t1316 * t904) * MDP(34) + (-t1065 * t788 - t1315 * t904) * MDP(35) + ((-t1085 * t1121 - t1086 * t1116) * t788 + t1106 * t704 + t700 * t1065 + (t1116 * t1186 - t1121 * t1185) * t904 + t1195 * t778 + t1315 * t736 - t1169 * t1108) * MDP(37) + (-(-t1085 * t1116 + t1086 * t1121) * t788 + t1106 * t703 + t700 * t1066 + (t1116 * t1185 + t1121 * t1186) * t904 - t1195 * t1180 - t1316 * t736 + t1169 * t1107) * MDP(38) + (MDP(18) * t907 + MDP(19) * t909 - MDP(21) * qJD(4) - MDP(23) * t829 - MDP(27) * t856 + MDP(28) * t854 - MDP(29) * t905 - MDP(30) * t724 + MDP(31) * t725 + MDP(34) * t1180 + MDP(35) * t778 - MDP(36) * t904 - MDP(37) * t701 + MDP(38) * t702) * t909; t856 * t854 * MDP(25) + (-t854 ^ 2 + t856 ^ 2) * MDP(26) + (t740 + t1336) * MDP(27) + (t1335 - t741) * MDP(28) + t791 * MDP(29) + (-g(1) * t825 - g(2) * t1384 + g(3) * t875 + t725 * t905 - t752 * t856 + t697) * MDP(30) + (g(1) * t826 - g(2) * t1383 + g(3) * t876 + t724 * t905 + t752 * t854 + t1165) * MDP(31) + (t703 + t1374) * MDP(34) + (-t704 - t1368) * MDP(35) + (-(-t1116 * t718 - t1304) * t904 + (t1121 * t788 - t1290 * t904 - t778 * t856) * pkin(5) + t692 + t1361) * MDP(37) + ((-t719 * t904 - t694) * t1116 + (t718 * t904 - t1224) * t1121 + (-t1116 * t788 + t1180 * t856 - t1289 * t904) * pkin(5) + t1362) * MDP(38) + t1359; (t1263 + t1374) * MDP(34) + (-t1222 - t1368) * MDP(35) + (t702 * t904 + t1223 + t1361) * MDP(37) + (-t1121 * t695 + t701 * t904 - t1310 + t1362) * MDP(38) + (-MDP(34) * t1309 + MDP(35) * t1180 - MDP(37) * t702 - MDP(38) * t1305) * qJD(6) + t1359;];
tau  = t1;