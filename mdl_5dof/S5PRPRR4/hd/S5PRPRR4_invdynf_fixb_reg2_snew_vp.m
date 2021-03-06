% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5PRPRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:52:13
% EndTime: 2019-12-05 15:52:15
% DurationCPUTime: 2.63s
% Computational Cost: add. (10098->195), mult. (19025->313), div. (0->0), fcn. (14269->12), ass. (0->172)
t1327 = sin(pkin(9));
t1330 = cos(pkin(9));
t1304 = t1327 * g(1) - t1330 * g(2);
t1323 = -g(3) + qJDD(1);
t1328 = sin(pkin(5));
t1331 = cos(pkin(5));
t1373 = t1304 * t1331 + t1323 * t1328;
t1338 = cos(qJ(4));
t1364 = t1338 * qJD(2);
t1313 = -qJD(5) + t1364;
t1372 = qJD(5) - t1313;
t1334 = sin(qJ(5));
t1337 = cos(qJ(5));
t1335 = sin(qJ(4));
t1368 = qJD(2) * t1335;
t1291 = -t1337 * qJD(4) + t1334 * t1368;
t1371 = t1291 ^ 2;
t1293 = t1334 * qJD(4) + t1337 * t1368;
t1370 = t1293 ^ 2;
t1369 = t1313 ^ 2;
t1367 = t1293 * t1291;
t1363 = qJD(5) + t1313;
t1305 = -t1330 * g(1) - t1327 * g(2);
t1336 = sin(qJ(2));
t1339 = cos(qJ(2));
t1264 = t1339 * t1305 + t1373 * t1336;
t1341 = qJD(2) ^ 2;
t1259 = -t1341 * pkin(2) + t1264;
t1326 = sin(pkin(10));
t1329 = cos(pkin(10));
t1263 = -t1336 * t1305 + t1373 * t1339;
t1342 = qJDD(2) * pkin(2) + t1263;
t1239 = t1329 * t1259 + t1326 * t1342;
t1235 = -t1341 * pkin(3) + qJDD(2) * pkin(7) + t1239;
t1282 = -t1328 * t1304 + t1331 * t1323;
t1281 = qJDD(3) + t1282;
t1228 = t1338 * t1235 + t1335 * t1281;
t1321 = t1335 ^ 2;
t1322 = t1338 ^ 2;
t1362 = t1321 + t1322;
t1361 = t1335 * qJDD(2);
t1360 = qJD(4) * t1368;
t1359 = qJD(4) * t1364;
t1238 = -t1326 * t1259 + t1329 * t1342;
t1296 = t1359 + t1361;
t1358 = t1337 * qJDD(4) - t1334 * t1296;
t1319 = t1338 * qJDD(2);
t1357 = -t1319 + 0.2e1 * t1360;
t1294 = (-pkin(4) * t1338 - pkin(8) * t1335) * qJD(2);
t1340 = qJD(4) ^ 2;
t1217 = -t1340 * pkin(4) + qJDD(4) * pkin(8) + t1294 * t1364 + t1228;
t1234 = -qJDD(2) * pkin(3) - t1341 * pkin(7) - t1238;
t1220 = (-t1296 - t1359) * pkin(8) + t1357 * pkin(4) + t1234;
t1202 = -t1334 * t1217 + t1337 * t1220;
t1203 = t1337 * t1217 + t1334 * t1220;
t1191 = -t1334 * t1202 + t1337 * t1203;
t1276 = t1338 * t1281;
t1216 = -t1276 - qJDD(4) * pkin(4) - t1340 * pkin(8) + (qJD(2) * t1294 + t1235) * t1335;
t1186 = t1338 * t1191 + t1335 * t1216;
t1190 = t1337 * t1202 + t1334 * t1203;
t1179 = t1326 * t1186 - t1329 * t1190;
t1180 = t1329 * t1186 + t1326 * t1190;
t1356 = t1179 * t1339 + t1180 * t1336;
t1227 = -t1335 * t1235 + t1276;
t1207 = -t1335 * t1227 + t1338 * t1228;
t1200 = t1326 * t1207 - t1329 * t1234;
t1201 = t1329 * t1207 + t1326 * t1234;
t1355 = t1200 * t1339 + t1201 * t1336;
t1249 = -t1363 * t1293 + t1358;
t1344 = -t1334 * qJDD(4) - t1337 * t1296;
t1251 = t1363 * t1291 + t1344;
t1230 = t1337 * t1249 - t1334 * t1251;
t1260 = -t1370 - t1371;
t1215 = t1338 * t1230 + t1335 * t1260;
t1229 = t1334 * t1249 + t1337 * t1251;
t1204 = t1326 * t1215 - t1329 * t1229;
t1205 = t1329 * t1215 + t1326 * t1229;
t1354 = t1204 * t1339 + t1205 * t1336;
t1343 = -qJDD(5) + t1319 - t1360;
t1262 = -t1343 - t1367;
t1269 = -t1369 - t1371;
t1242 = -t1334 * t1262 + t1337 * t1269;
t1248 = t1293 * t1372 - t1358;
t1219 = t1338 * t1242 + t1335 * t1248;
t1241 = t1337 * t1262 + t1334 * t1269;
t1208 = t1326 * t1219 - t1329 * t1241;
t1209 = t1329 * t1219 + t1326 * t1241;
t1353 = t1208 * t1339 + t1209 * t1336;
t1261 = t1343 - t1367;
t1274 = -t1369 - t1370;
t1244 = t1337 * t1261 - t1334 * t1274;
t1250 = -t1291 * t1372 - t1344;
t1222 = t1338 * t1244 + t1335 * t1250;
t1243 = t1334 * t1261 + t1337 * t1274;
t1210 = t1326 * t1222 - t1329 * t1243;
t1211 = t1329 * t1222 + t1326 * t1243;
t1352 = t1210 * t1339 + t1211 * t1336;
t1212 = t1329 * t1238 + t1326 * t1239;
t1213 = -t1326 * t1238 + t1329 * t1239;
t1351 = t1212 * t1339 + t1213 * t1336;
t1312 = t1335 * t1341 * t1338;
t1306 = qJDD(4) + t1312;
t1311 = -t1322 * t1341 - t1340;
t1279 = -t1335 * t1306 + t1338 * t1311;
t1255 = t1326 * t1279 - t1329 * t1357;
t1257 = t1329 * t1279 + t1326 * t1357;
t1350 = t1255 * t1339 + t1257 * t1336;
t1307 = -qJDD(4) + t1312;
t1310 = -t1321 * t1341 - t1340;
t1280 = t1338 * t1307 - t1335 * t1310;
t1295 = 0.2e1 * t1359 + t1361;
t1256 = t1326 * t1280 - t1329 * t1295;
t1258 = t1329 * t1280 + t1326 * t1295;
t1349 = t1256 * t1339 + t1258 * t1336;
t1348 = t1263 * t1339 + t1264 * t1336;
t1300 = t1362 * qJDD(2);
t1303 = t1362 * t1341;
t1272 = t1326 * t1300 + t1329 * t1303;
t1273 = t1329 * t1300 - t1326 * t1303;
t1347 = t1272 * t1339 + t1273 * t1336;
t1298 = t1329 * qJDD(2) - t1326 * t1341;
t1299 = -t1326 * qJDD(2) - t1329 * t1341;
t1346 = t1339 * t1298 + t1336 * t1299;
t1271 = -t1336 * t1298 + t1339 * t1299;
t1345 = t1339 * qJDD(2) - t1336 * t1341;
t1302 = -t1336 * qJDD(2) - t1339 * t1341;
t1287 = t1345 * t1331;
t1286 = t1302 * t1331;
t1285 = t1345 * t1328;
t1284 = t1302 * t1328;
t1278 = t1335 * t1307 + t1338 * t1310;
t1277 = t1338 * t1306 + t1335 * t1311;
t1268 = t1346 * t1331;
t1267 = t1271 * t1331;
t1266 = t1346 * t1328;
t1265 = t1271 * t1328;
t1247 = -t1336 * t1272 + t1339 * t1273;
t1246 = t1347 * t1331;
t1245 = t1347 * t1328;
t1240 = -t1336 * t1263 + t1339 * t1264;
t1237 = -t1336 * t1256 + t1339 * t1258;
t1236 = -t1336 * t1255 + t1339 * t1257;
t1232 = -t1328 * t1282 + t1348 * t1331;
t1231 = t1331 * t1282 + t1348 * t1328;
t1226 = -t1328 * t1278 + t1349 * t1331;
t1225 = -t1328 * t1277 + t1350 * t1331;
t1224 = t1331 * t1278 + t1349 * t1328;
t1223 = t1331 * t1277 + t1350 * t1328;
t1221 = t1335 * t1244 - t1338 * t1250;
t1218 = t1335 * t1242 - t1338 * t1248;
t1214 = t1335 * t1230 - t1338 * t1260;
t1206 = t1338 * t1227 + t1335 * t1228;
t1199 = -t1336 * t1212 + t1339 * t1213;
t1198 = -t1328 * t1281 + t1351 * t1331;
t1197 = t1331 * t1281 + t1351 * t1328;
t1196 = -t1336 * t1210 + t1339 * t1211;
t1195 = -t1336 * t1208 + t1339 * t1209;
t1194 = -t1336 * t1204 + t1339 * t1205;
t1193 = -t1328 * t1221 + t1352 * t1331;
t1192 = t1331 * t1221 + t1352 * t1328;
t1189 = -t1328 * t1218 + t1353 * t1331;
t1188 = t1331 * t1218 + t1353 * t1328;
t1187 = -t1336 * t1200 + t1339 * t1201;
t1185 = t1335 * t1191 - t1338 * t1216;
t1184 = -t1328 * t1214 + t1354 * t1331;
t1183 = t1331 * t1214 + t1354 * t1328;
t1182 = -t1328 * t1206 + t1355 * t1331;
t1181 = t1331 * t1206 + t1355 * t1328;
t1178 = -t1336 * t1179 + t1339 * t1180;
t1177 = -t1328 * t1185 + t1356 * t1331;
t1176 = t1331 * t1185 + t1356 * t1328;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1327 * t1304 + t1330 * t1305, 0, 0, 0, 0, 0, 0, -t1327 * t1287 + t1330 * t1302, -t1327 * t1286 - t1330 * t1345, 0, -t1327 * t1232 + t1330 * t1240, 0, 0, 0, 0, 0, 0, -t1327 * t1268 + t1330 * t1271, -t1327 * t1267 - t1330 * t1346, 0, -t1327 * t1198 + t1330 * t1199, 0, 0, 0, 0, 0, 0, -t1327 * t1225 + t1330 * t1236, -t1327 * t1226 + t1330 * t1237, -t1327 * t1246 + t1330 * t1247, -t1327 * t1182 + t1330 * t1187, 0, 0, 0, 0, 0, 0, -t1327 * t1189 + t1330 * t1195, -t1327 * t1193 + t1330 * t1196, -t1327 * t1184 + t1330 * t1194, -t1327 * t1177 + t1330 * t1178; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1330 * t1304 + t1327 * t1305, 0, 0, 0, 0, 0, 0, t1330 * t1287 + t1327 * t1302, t1330 * t1286 - t1327 * t1345, 0, t1330 * t1232 + t1327 * t1240, 0, 0, 0, 0, 0, 0, t1330 * t1268 + t1327 * t1271, t1330 * t1267 - t1327 * t1346, 0, t1330 * t1198 + t1327 * t1199, 0, 0, 0, 0, 0, 0, t1330 * t1225 + t1327 * t1236, t1330 * t1226 + t1327 * t1237, t1330 * t1246 + t1327 * t1247, t1330 * t1182 + t1327 * t1187, 0, 0, 0, 0, 0, 0, t1330 * t1189 + t1327 * t1195, t1330 * t1193 + t1327 * t1196, t1330 * t1184 + t1327 * t1194, t1330 * t1177 + t1327 * t1178; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1323, 0, 0, 0, 0, 0, 0, t1285, t1284, 0, t1231, 0, 0, 0, 0, 0, 0, t1266, t1265, 0, t1197, 0, 0, 0, 0, 0, 0, t1223, t1224, t1245, t1181, 0, 0, 0, 0, 0, 0, t1188, t1192, t1183, t1176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1305, 0, 0, 0, 0, 0, 0, t1302, -t1345, 0, t1240, 0, 0, 0, 0, 0, 0, t1271, -t1346, 0, t1199, 0, 0, 0, 0, 0, 0, t1236, t1237, t1247, t1187, 0, 0, 0, 0, 0, 0, t1195, t1196, t1194, t1178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1304, 0, 0, 0, 0, 0, 0, t1287, t1286, 0, t1232, 0, 0, 0, 0, 0, 0, t1268, t1267, 0, t1198, 0, 0, 0, 0, 0, 0, t1225, t1226, t1246, t1182, 0, 0, 0, 0, 0, 0, t1189, t1193, t1184, t1177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1323, 0, 0, 0, 0, 0, 0, t1285, t1284, 0, t1231, 0, 0, 0, 0, 0, 0, t1266, t1265, 0, t1197, 0, 0, 0, 0, 0, 0, t1223, t1224, t1245, t1181, 0, 0, 0, 0, 0, 0, t1188, t1192, t1183, t1176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1341, -qJDD(2), 0, t1264, 0, 0, 0, 0, 0, 0, t1299, -t1298, 0, t1213, 0, 0, 0, 0, 0, 0, t1257, t1258, t1273, t1201, 0, 0, 0, 0, 0, 0, t1209, t1211, t1205, t1180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1341, 0, t1263, 0, 0, 0, 0, 0, 0, t1298, t1299, 0, t1212, 0, 0, 0, 0, 0, 0, t1255, t1256, t1272, t1200, 0, 0, 0, 0, 0, 0, t1208, t1210, t1204, t1179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1282, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1281, 0, 0, 0, 0, 0, 0, t1277, t1278, 0, t1206, 0, 0, 0, 0, 0, 0, t1218, t1221, t1214, t1185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1341, -qJDD(2), 0, t1239, 0, 0, 0, 0, 0, 0, t1279, t1280, t1300, t1207, 0, 0, 0, 0, 0, 0, t1219, t1222, t1215, t1186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1341, 0, t1238, 0, 0, 0, 0, 0, 0, -t1357, -t1295, t1303, -t1234, 0, 0, 0, 0, 0, 0, -t1241, -t1243, -t1229, -t1190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1281, 0, 0, 0, 0, 0, 0, t1277, t1278, 0, t1206, 0, 0, 0, 0, 0, 0, t1218, t1221, t1214, t1185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1311, t1307, t1319, t1228, 0, 0, 0, 0, 0, 0, t1242, t1244, t1230, t1191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1306, t1310, -t1361, t1227, 0, 0, 0, 0, 0, 0, -t1248, -t1250, -t1260, -t1216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1357, t1295, -t1303, t1234, 0, 0, 0, 0, 0, 0, t1241, t1243, t1229, t1190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1269, t1261, t1249, t1203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, t1274, t1251, t1202; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1248, t1250, t1260, t1216;];
f_new_reg = t1;
