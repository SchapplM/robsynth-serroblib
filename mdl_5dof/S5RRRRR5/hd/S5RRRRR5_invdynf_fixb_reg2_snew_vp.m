% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRRRR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:54
% EndTime: 2020-01-03 12:13:56
% DurationCPUTime: 2.41s
% Computational Cost: add. (14328->190), mult. (15217->267), div. (0->0), fcn. (9583->10), ass. (0->155)
t1286 = qJD(1) + qJD(2);
t1282 = qJD(3) + t1286;
t1280 = t1282 ^ 2;
t1284 = qJDD(1) + qJDD(2);
t1281 = qJDD(3) + t1284;
t1291 = sin(qJ(3));
t1296 = cos(qJ(3));
t1255 = t1291 * t1280 - t1296 * t1281;
t1292 = sin(qJ(2));
t1297 = cos(qJ(2));
t1306 = -t1296 * t1280 - t1291 * t1281;
t1230 = t1297 * t1255 - t1292 * t1306;
t1293 = sin(qJ(1));
t1298 = cos(qJ(1));
t1331 = t1292 * t1255 + t1297 * t1306;
t1335 = t1293 * t1230 + t1298 * t1331;
t1334 = -t1298 * t1230 + t1293 * t1331;
t1283 = t1286 ^ 2;
t1262 = t1292 * t1283 - t1297 * t1284;
t1305 = -t1297 * t1283 - t1292 * t1284;
t1330 = t1293 * t1262 + t1298 * t1305;
t1329 = -t1298 * t1262 + t1293 * t1305;
t1285 = qJD(4) + qJD(5);
t1324 = qJD(5) + t1285;
t1289 = sin(qJ(5));
t1294 = cos(qJ(5));
t1295 = cos(qJ(4));
t1317 = t1282 * t1295;
t1290 = sin(qJ(4));
t1318 = t1282 * t1290;
t1244 = t1289 * t1318 - t1294 * t1317;
t1323 = t1244 ^ 2;
t1246 = (t1289 * t1295 + t1290 * t1294) * t1282;
t1322 = t1246 ^ 2;
t1321 = t1285 ^ 2;
t1320 = t1246 * t1244;
t1319 = t1280 * t1290;
t1288 = t1295 ^ 2;
t1316 = t1288 * t1280;
t1275 = -t1298 * g(2) - t1293 * g(3);
t1303 = qJDD(1) * pkin(1) + t1275;
t1274 = -t1293 * g(2) + t1298 * g(3);
t1300 = qJD(1) ^ 2;
t1304 = -t1300 * pkin(1) + t1274;
t1237 = t1292 * t1303 + t1297 * t1304;
t1235 = -t1283 * pkin(2) + t1237;
t1236 = -t1292 * t1304 + t1297 * t1303;
t1301 = t1284 * pkin(2) + t1236;
t1215 = t1296 * t1235 + t1291 * t1301;
t1209 = -t1280 * pkin(3) + t1281 * pkin(8) + t1215;
t1315 = t1290 * t1209;
t1314 = t1290 * t1281;
t1313 = qJD(5) - t1285;
t1287 = t1290 ^ 2;
t1312 = t1287 + t1288;
t1311 = -qJDD(4) - qJDD(5);
t1310 = qJD(4) * t1318;
t1309 = qJD(4) * t1317;
t1202 = -t1290 * g(1) + t1295 * t1209;
t1214 = -t1291 * t1235 + t1296 * t1301;
t1250 = t1309 + t1314;
t1277 = t1295 * t1281;
t1307 = -t1277 + t1310;
t1308 = -t1289 * t1250 - t1294 * t1307;
t1208 = -t1281 * pkin(3) - t1280 * pkin(8) - t1214;
t1302 = -t1294 * t1250 + t1289 * t1307;
t1299 = qJD(4) ^ 2;
t1272 = t1298 * qJDD(1) - t1293 * t1300;
t1271 = -t1293 * qJDD(1) - t1298 * t1300;
t1270 = t1295 * t1319;
t1268 = -t1299 - t1316;
t1267 = -t1287 * t1280 - t1299;
t1266 = qJD(4) * pkin(4) - pkin(9) * t1318;
t1265 = -qJDD(4) + t1270;
t1264 = qJDD(4) + t1270;
t1257 = t1312 * t1280;
t1252 = t1312 * t1281;
t1251 = t1277 - 0.2e1 * t1310;
t1249 = 0.2e1 * t1309 + t1314;
t1242 = -t1321 - t1322;
t1241 = t1295 * t1265 - t1290 * t1267;
t1240 = -t1290 * t1264 + t1295 * t1268;
t1239 = t1290 * t1265 + t1295 * t1267;
t1238 = t1295 * t1264 + t1290 * t1268;
t1229 = t1296 * t1252 - t1291 * t1257;
t1226 = t1291 * t1252 + t1296 * t1257;
t1225 = t1311 - t1320;
t1224 = -t1311 - t1320;
t1223 = -t1321 - t1323;
t1222 = t1296 * t1241 + t1291 * t1249;
t1221 = t1296 * t1240 - t1291 * t1251;
t1220 = t1291 * t1241 - t1296 * t1249;
t1219 = t1291 * t1240 + t1296 * t1251;
t1218 = -t1322 - t1323;
t1217 = -t1292 * t1236 + t1297 * t1237;
t1216 = t1297 * t1236 + t1292 * t1237;
t1213 = t1294 * t1225 - t1289 * t1242;
t1212 = t1289 * t1225 + t1294 * t1242;
t1211 = -t1292 * t1226 + t1297 * t1229;
t1210 = t1297 * t1226 + t1292 * t1229;
t1207 = t1313 * t1244 + t1302;
t1206 = -t1324 * t1244 - t1302;
t1205 = -t1313 * t1246 + t1308;
t1204 = t1324 * t1246 - t1308;
t1201 = -t1295 * g(1) - t1315;
t1200 = t1294 * t1223 - t1289 * t1224;
t1199 = t1289 * t1223 + t1294 * t1224;
t1198 = -t1292 * t1220 + t1297 * t1222;
t1197 = -t1292 * t1219 + t1297 * t1221;
t1196 = t1297 * t1220 + t1292 * t1222;
t1195 = t1297 * t1219 + t1292 * t1221;
t1194 = t1307 * pkin(4) - pkin(9) * t1316 + t1266 * t1318 + t1208;
t1193 = -pkin(4) * t1316 - t1307 * pkin(9) - qJD(4) * t1266 + t1202;
t1192 = qJDD(4) * pkin(4) - t1250 * pkin(9) - t1315 + (qJD(4) * t1282 * pkin(9) + pkin(4) * t1319 - g(1)) * t1295;
t1191 = -t1291 * t1214 + t1296 * t1215;
t1190 = t1296 * t1214 + t1291 * t1215;
t1189 = -t1290 * t1212 + t1295 * t1213;
t1188 = t1295 * t1212 + t1290 * t1213;
t1187 = t1294 * t1205 - t1289 * t1207;
t1186 = t1289 * t1205 + t1294 * t1207;
t1185 = -t1290 * t1201 + t1295 * t1202;
t1184 = t1295 * t1201 + t1290 * t1202;
t1183 = -t1290 * t1199 + t1295 * t1200;
t1182 = t1295 * t1199 + t1290 * t1200;
t1181 = t1296 * t1189 + t1291 * t1206;
t1180 = t1291 * t1189 - t1296 * t1206;
t1179 = t1289 * t1192 + t1294 * t1193;
t1178 = t1294 * t1192 - t1289 * t1193;
t1177 = t1296 * t1185 + t1291 * t1208;
t1176 = t1291 * t1185 - t1296 * t1208;
t1175 = t1296 * t1183 + t1291 * t1204;
t1174 = t1291 * t1183 - t1296 * t1204;
t1173 = -t1292 * t1190 + t1297 * t1191;
t1172 = t1297 * t1190 + t1292 * t1191;
t1171 = -t1290 * t1186 + t1295 * t1187;
t1170 = t1295 * t1186 + t1290 * t1187;
t1169 = t1296 * t1171 + t1291 * t1218;
t1168 = t1291 * t1171 - t1296 * t1218;
t1167 = -t1292 * t1180 + t1297 * t1181;
t1166 = t1297 * t1180 + t1292 * t1181;
t1165 = -t1289 * t1178 + t1294 * t1179;
t1164 = t1294 * t1178 + t1289 * t1179;
t1163 = -t1292 * t1176 + t1297 * t1177;
t1162 = t1297 * t1176 + t1292 * t1177;
t1161 = -t1292 * t1174 + t1297 * t1175;
t1160 = t1297 * t1174 + t1292 * t1175;
t1159 = -t1292 * t1168 + t1297 * t1169;
t1158 = t1297 * t1168 + t1292 * t1169;
t1157 = -t1290 * t1164 + t1295 * t1165;
t1156 = t1295 * t1164 + t1290 * t1165;
t1155 = t1296 * t1157 + t1291 * t1194;
t1154 = t1291 * t1157 - t1296 * t1194;
t1153 = -t1292 * t1154 + t1297 * t1155;
t1152 = t1297 * t1154 + t1292 * t1155;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1238, t1239, 0, t1184, 0, 0, 0, 0, 0, 0, t1182, t1188, t1170, t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1272, t1271, 0, t1293 * t1274 + t1298 * t1275, 0, 0, 0, 0, 0, 0, t1329, t1330, 0, t1298 * t1216 + t1293 * t1217, 0, 0, 0, 0, 0, 0, t1334, t1335, 0, t1298 * t1172 + t1293 * t1173, 0, 0, 0, 0, 0, 0, t1298 * t1195 + t1293 * t1197, t1298 * t1196 + t1293 * t1198, t1298 * t1210 + t1293 * t1211, t1298 * t1162 + t1293 * t1163, 0, 0, 0, 0, 0, 0, t1298 * t1160 + t1293 * t1161, t1298 * t1166 + t1293 * t1167, t1298 * t1158 + t1293 * t1159, t1298 * t1152 + t1293 * t1153; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t1271, t1272, 0, -t1298 * t1274 + t1293 * t1275, 0, 0, 0, 0, 0, 0, -t1330, t1329, 0, t1293 * t1216 - t1298 * t1217, 0, 0, 0, 0, 0, 0, -t1335, t1334, 0, t1293 * t1172 - t1298 * t1173, 0, 0, 0, 0, 0, 0, t1293 * t1195 - t1298 * t1197, t1293 * t1196 - t1298 * t1198, t1293 * t1210 - t1298 * t1211, t1293 * t1162 - t1298 * t1163, 0, 0, 0, 0, 0, 0, t1293 * t1160 - t1298 * t1161, t1293 * t1166 - t1298 * t1167, t1293 * t1158 - t1298 * t1159, t1293 * t1152 - t1298 * t1153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1300, -qJDD(1), 0, t1274, 0, 0, 0, 0, 0, 0, t1305, t1262, 0, t1217, 0, 0, 0, 0, 0, 0, t1331, t1230, 0, t1173, 0, 0, 0, 0, 0, 0, t1197, t1198, t1211, t1163, 0, 0, 0, 0, 0, 0, t1161, t1167, t1159, t1153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1300, 0, t1275, 0, 0, 0, 0, 0, 0, -t1262, t1305, 0, t1216, 0, 0, 0, 0, 0, 0, -t1230, t1331, 0, t1172, 0, 0, 0, 0, 0, 0, t1195, t1196, t1210, t1162, 0, 0, 0, 0, 0, 0, t1160, t1166, t1158, t1152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1238, t1239, 0, t1184, 0, 0, 0, 0, 0, 0, t1182, t1188, t1170, t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1283, -t1284, 0, t1237, 0, 0, 0, 0, 0, 0, t1306, t1255, 0, t1191, 0, 0, 0, 0, 0, 0, t1221, t1222, t1229, t1177, 0, 0, 0, 0, 0, 0, t1175, t1181, t1169, t1155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1284, -t1283, 0, t1236, 0, 0, 0, 0, 0, 0, -t1255, t1306, 0, t1190, 0, 0, 0, 0, 0, 0, t1219, t1220, t1226, t1176, 0, 0, 0, 0, 0, 0, t1174, t1180, t1168, t1154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1238, t1239, 0, t1184, 0, 0, 0, 0, 0, 0, t1182, t1188, t1170, t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1280, -t1281, 0, t1215, 0, 0, 0, 0, 0, 0, t1240, t1241, t1252, t1185, 0, 0, 0, 0, 0, 0, t1183, t1189, t1171, t1157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1281, -t1280, 0, t1214, 0, 0, 0, 0, 0, 0, t1251, -t1249, t1257, -t1208, 0, 0, 0, 0, 0, 0, -t1204, -t1206, -t1218, -t1194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1238, t1239, 0, t1184, 0, 0, 0, 0, 0, 0, t1182, t1188, t1170, t1156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1268, t1265, t1277, t1202, 0, 0, 0, 0, 0, 0, t1200, t1213, t1187, t1165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1264, t1267, -t1314, t1201, 0, 0, 0, 0, 0, 0, t1199, t1212, t1186, t1164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1251, t1249, -t1257, t1208, 0, 0, 0, 0, 0, 0, t1204, t1206, t1218, t1194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1223, t1225, t1205, t1179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1224, t1242, t1207, t1178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1204, t1206, t1218, t1194;];
f_new_reg = t1;