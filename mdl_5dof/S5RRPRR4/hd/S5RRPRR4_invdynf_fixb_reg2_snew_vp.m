% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRPRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:02:35
% EndTime: 2020-01-03 12:02:38
% DurationCPUTime: 2.36s
% Computational Cost: add. (11445->183), mult. (15217->261), div. (0->0), fcn. (9583->10), ass. (0->151)
t1265 = qJD(1) + qJD(2);
t1262 = t1265 ^ 2;
t1263 = qJDD(1) + qJDD(2);
t1269 = sin(pkin(9));
t1270 = cos(pkin(9));
t1235 = t1269 * t1262 - t1270 * t1263;
t1273 = sin(qJ(2));
t1277 = cos(qJ(2));
t1286 = -t1270 * t1262 - t1269 * t1263;
t1211 = t1277 * t1235 - t1273 * t1286;
t1274 = sin(qJ(1));
t1278 = cos(qJ(1));
t1309 = t1273 * t1235 + t1277 * t1286;
t1313 = t1274 * t1211 + t1278 * t1309;
t1312 = -t1278 * t1211 + t1274 * t1309;
t1241 = t1273 * t1262 - t1277 * t1263;
t1285 = -t1277 * t1262 - t1273 * t1263;
t1308 = t1274 * t1241 + t1278 * t1285;
t1307 = -t1278 * t1241 + t1274 * t1285;
t1264 = qJD(4) + qJD(5);
t1302 = qJD(5) + t1264;
t1271 = sin(qJ(5));
t1275 = cos(qJ(5));
t1276 = cos(qJ(4));
t1296 = t1265 * t1276;
t1272 = sin(qJ(4));
t1297 = t1265 * t1272;
t1224 = t1271 * t1297 - t1275 * t1296;
t1301 = t1224 ^ 2;
t1226 = (t1271 * t1276 + t1272 * t1275) * t1265;
t1300 = t1226 ^ 2;
t1299 = t1264 ^ 2;
t1298 = t1226 * t1224;
t1267 = t1276 ^ 2;
t1295 = t1267 * t1262;
t1294 = t1272 * t1263;
t1293 = qJD(5) - t1264;
t1254 = -t1278 * g(2) - t1274 * g(3);
t1283 = qJDD(1) * pkin(1) + t1254;
t1253 = -t1274 * g(2) + t1278 * g(3);
t1280 = qJD(1) ^ 2;
t1284 = -t1280 * pkin(1) + t1253;
t1217 = t1273 * t1283 + t1277 * t1284;
t1215 = -t1262 * pkin(2) + t1217;
t1216 = -t1273 * t1284 + t1277 * t1283;
t1281 = t1263 * pkin(2) + t1216;
t1193 = t1270 * t1215 + t1269 * t1281;
t1185 = -t1262 * pkin(3) + t1263 * pkin(7) + t1193;
t1268 = -g(1) + qJDD(3);
t1180 = t1276 * t1185 + t1272 * t1268;
t1266 = t1272 ^ 2;
t1292 = t1266 + t1267;
t1252 = t1276 * t1262 * t1272;
t1244 = qJDD(4) + t1252;
t1291 = -qJDD(4) - qJDD(5);
t1290 = qJD(4) * t1297;
t1289 = qJD(4) * t1296;
t1179 = -t1272 * t1185 + t1276 * t1268;
t1192 = -t1269 * t1215 + t1270 * t1281;
t1230 = t1289 + t1294;
t1257 = t1276 * t1263;
t1287 = -t1257 + t1290;
t1288 = -t1271 * t1230 - t1275 * t1287;
t1184 = -t1263 * pkin(3) - t1262 * pkin(7) - t1192;
t1282 = -t1275 * t1230 + t1271 * t1287;
t1279 = qJD(4) ^ 2;
t1250 = -t1279 - t1295;
t1249 = -t1266 * t1262 - t1279;
t1248 = t1278 * qJDD(1) - t1274 * t1280;
t1247 = -t1274 * qJDD(1) - t1278 * t1280;
t1246 = qJD(4) * pkin(4) - pkin(8) * t1297;
t1245 = -qJDD(4) + t1252;
t1243 = t1292 * t1262;
t1238 = t1292 * t1263;
t1231 = t1257 - 0.2e1 * t1290;
t1229 = 0.2e1 * t1289 + t1294;
t1222 = -t1299 - t1300;
t1221 = t1276 * t1245 - t1272 * t1249;
t1220 = -t1272 * t1244 + t1276 * t1250;
t1219 = t1272 * t1245 + t1276 * t1249;
t1218 = t1276 * t1244 + t1272 * t1250;
t1214 = t1270 * t1238 - t1269 * t1243;
t1213 = t1269 * t1238 + t1270 * t1243;
t1205 = t1291 - t1298;
t1204 = -t1291 - t1298;
t1203 = -t1299 - t1301;
t1202 = t1270 * t1221 + t1269 * t1229;
t1201 = t1270 * t1220 - t1269 * t1231;
t1200 = t1269 * t1221 - t1270 * t1229;
t1199 = t1269 * t1220 + t1270 * t1231;
t1198 = -t1300 - t1301;
t1197 = -t1273 * t1216 + t1277 * t1217;
t1196 = t1277 * t1216 + t1273 * t1217;
t1195 = t1275 * t1205 - t1271 * t1222;
t1194 = t1271 * t1205 + t1275 * t1222;
t1191 = -t1273 * t1213 + t1277 * t1214;
t1190 = t1277 * t1213 + t1273 * t1214;
t1189 = t1293 * t1224 + t1282;
t1188 = -t1302 * t1224 - t1282;
t1187 = -t1293 * t1226 + t1288;
t1186 = t1302 * t1226 - t1288;
t1182 = t1275 * t1203 - t1271 * t1204;
t1181 = t1271 * t1203 + t1275 * t1204;
t1178 = -t1273 * t1200 + t1277 * t1202;
t1177 = -t1273 * t1199 + t1277 * t1201;
t1176 = t1277 * t1200 + t1273 * t1202;
t1175 = t1277 * t1199 + t1273 * t1201;
t1174 = t1287 * pkin(4) - pkin(8) * t1295 + t1246 * t1297 + t1184;
t1173 = -pkin(4) * t1295 - t1287 * pkin(8) - qJD(4) * t1246 + t1180;
t1172 = (-t1230 + t1289) * pkin(8) + t1244 * pkin(4) + t1179;
t1171 = -t1272 * t1194 + t1276 * t1195;
t1170 = t1276 * t1194 + t1272 * t1195;
t1169 = -t1269 * t1192 + t1270 * t1193;
t1168 = t1270 * t1192 + t1269 * t1193;
t1167 = t1275 * t1187 - t1271 * t1189;
t1166 = t1271 * t1187 + t1275 * t1189;
t1165 = -t1272 * t1181 + t1276 * t1182;
t1164 = t1276 * t1181 + t1272 * t1182;
t1163 = -t1272 * t1179 + t1276 * t1180;
t1162 = t1276 * t1179 + t1272 * t1180;
t1161 = t1270 * t1171 + t1269 * t1188;
t1160 = t1269 * t1171 - t1270 * t1188;
t1159 = t1271 * t1172 + t1275 * t1173;
t1158 = t1270 * t1165 + t1269 * t1186;
t1157 = t1275 * t1172 - t1271 * t1173;
t1156 = t1269 * t1165 - t1270 * t1186;
t1155 = t1270 * t1163 + t1269 * t1184;
t1154 = t1269 * t1163 - t1270 * t1184;
t1153 = -t1273 * t1168 + t1277 * t1169;
t1152 = t1277 * t1168 + t1273 * t1169;
t1151 = -t1272 * t1166 + t1276 * t1167;
t1150 = t1276 * t1166 + t1272 * t1167;
t1149 = t1270 * t1151 + t1269 * t1198;
t1148 = t1269 * t1151 - t1270 * t1198;
t1147 = -t1273 * t1160 + t1277 * t1161;
t1146 = t1277 * t1160 + t1273 * t1161;
t1145 = -t1273 * t1156 + t1277 * t1158;
t1144 = -t1271 * t1157 + t1275 * t1159;
t1143 = t1277 * t1156 + t1273 * t1158;
t1142 = t1275 * t1157 + t1271 * t1159;
t1141 = -t1273 * t1154 + t1277 * t1155;
t1140 = t1277 * t1154 + t1273 * t1155;
t1139 = -t1273 * t1148 + t1277 * t1149;
t1138 = t1277 * t1148 + t1273 * t1149;
t1137 = -t1272 * t1142 + t1276 * t1144;
t1136 = t1276 * t1142 + t1272 * t1144;
t1135 = t1270 * t1137 + t1269 * t1174;
t1134 = t1269 * t1137 - t1270 * t1174;
t1133 = -t1273 * t1134 + t1277 * t1135;
t1132 = t1277 * t1134 + t1273 * t1135;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1268, 0, 0, 0, 0, 0, 0, t1218, t1219, 0, t1162, 0, 0, 0, 0, 0, 0, t1164, t1170, t1150, t1136; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1248, t1247, 0, t1274 * t1253 + t1278 * t1254, 0, 0, 0, 0, 0, 0, t1307, t1308, 0, t1278 * t1196 + t1274 * t1197, 0, 0, 0, 0, 0, 0, t1312, t1313, 0, t1278 * t1152 + t1274 * t1153, 0, 0, 0, 0, 0, 0, t1278 * t1175 + t1274 * t1177, t1278 * t1176 + t1274 * t1178, t1278 * t1190 + t1274 * t1191, t1278 * t1140 + t1274 * t1141, 0, 0, 0, 0, 0, 0, t1278 * t1143 + t1274 * t1145, t1278 * t1146 + t1274 * t1147, t1278 * t1138 + t1274 * t1139, t1278 * t1132 + t1274 * t1133; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, -t1247, t1248, 0, -t1278 * t1253 + t1274 * t1254, 0, 0, 0, 0, 0, 0, -t1308, t1307, 0, t1274 * t1196 - t1278 * t1197, 0, 0, 0, 0, 0, 0, -t1313, t1312, 0, t1274 * t1152 - t1278 * t1153, 0, 0, 0, 0, 0, 0, t1274 * t1175 - t1278 * t1177, t1274 * t1176 - t1278 * t1178, t1274 * t1190 - t1278 * t1191, t1274 * t1140 - t1278 * t1141, 0, 0, 0, 0, 0, 0, t1274 * t1143 - t1278 * t1145, t1274 * t1146 - t1278 * t1147, t1274 * t1138 - t1278 * t1139, t1274 * t1132 - t1278 * t1133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1280, -qJDD(1), 0, t1253, 0, 0, 0, 0, 0, 0, t1285, t1241, 0, t1197, 0, 0, 0, 0, 0, 0, t1309, t1211, 0, t1153, 0, 0, 0, 0, 0, 0, t1177, t1178, t1191, t1141, 0, 0, 0, 0, 0, 0, t1145, t1147, t1139, t1133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1280, 0, t1254, 0, 0, 0, 0, 0, 0, -t1241, t1285, 0, t1196, 0, 0, 0, 0, 0, 0, -t1211, t1309, 0, t1152, 0, 0, 0, 0, 0, 0, t1175, t1176, t1190, t1140, 0, 0, 0, 0, 0, 0, t1143, t1146, t1138, t1132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1268, 0, 0, 0, 0, 0, 0, t1218, t1219, 0, t1162, 0, 0, 0, 0, 0, 0, t1164, t1170, t1150, t1136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1262, -t1263, 0, t1217, 0, 0, 0, 0, 0, 0, t1286, t1235, 0, t1169, 0, 0, 0, 0, 0, 0, t1201, t1202, t1214, t1155, 0, 0, 0, 0, 0, 0, t1158, t1161, t1149, t1135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1263, -t1262, 0, t1216, 0, 0, 0, 0, 0, 0, -t1235, t1286, 0, t1168, 0, 0, 0, 0, 0, 0, t1199, t1200, t1213, t1154, 0, 0, 0, 0, 0, 0, t1156, t1160, t1148, t1134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1268, 0, 0, 0, 0, 0, 0, t1218, t1219, 0, t1162, 0, 0, 0, 0, 0, 0, t1164, t1170, t1150, t1136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1262, -t1263, 0, t1193, 0, 0, 0, 0, 0, 0, t1220, t1221, t1238, t1163, 0, 0, 0, 0, 0, 0, t1165, t1171, t1151, t1137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1263, -t1262, 0, t1192, 0, 0, 0, 0, 0, 0, t1231, -t1229, t1243, -t1184, 0, 0, 0, 0, 0, 0, -t1186, -t1188, -t1198, -t1174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1268, 0, 0, 0, 0, 0, 0, t1218, t1219, 0, t1162, 0, 0, 0, 0, 0, 0, t1164, t1170, t1150, t1136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1250, t1245, t1257, t1180, 0, 0, 0, 0, 0, 0, t1182, t1195, t1167, t1144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1244, t1249, -t1294, t1179, 0, 0, 0, 0, 0, 0, t1181, t1194, t1166, t1142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1231, t1229, -t1243, t1184, 0, 0, 0, 0, 0, 0, t1186, t1188, t1198, t1174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1203, t1205, t1187, t1159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1204, t1222, t1189, t1157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1186, t1188, t1198, t1174;];
f_new_reg = t1;
