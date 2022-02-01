% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPRRR5_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:49:36
% EndTime: 2022-01-20 09:49:39
% DurationCPUTime: 2.93s
% Computational Cost: add. (10884->180), mult. (15217->261), div. (0->0), fcn. (9583->10), ass. (0->151)
t1259 = qJD(1) + qJD(3);
t1256 = t1259 ^ 2;
t1257 = qJDD(1) + qJDD(3);
t1267 = sin(qJ(3));
t1271 = cos(qJ(3));
t1231 = t1267 * t1256 - t1271 * t1257;
t1263 = sin(pkin(9));
t1264 = cos(pkin(9));
t1279 = -t1271 * t1256 - t1267 * t1257;
t1205 = t1264 * t1231 - t1263 * t1279;
t1268 = sin(qJ(1));
t1272 = cos(qJ(1));
t1300 = t1263 * t1231 + t1264 * t1279;
t1304 = t1268 * t1205 + t1272 * t1300;
t1303 = t1272 * t1205 - t1268 * t1300;
t1258 = qJD(4) + qJD(5);
t1297 = qJD(5) + t1258;
t1265 = sin(qJ(5));
t1269 = cos(qJ(5));
t1270 = cos(qJ(4));
t1291 = t1259 * t1270;
t1266 = sin(qJ(4));
t1292 = t1259 * t1266;
t1218 = t1265 * t1292 - t1269 * t1291;
t1296 = t1218 ^ 2;
t1220 = (t1265 * t1270 + t1266 * t1269) * t1259;
t1295 = t1220 ^ 2;
t1294 = t1258 ^ 2;
t1293 = t1220 * t1218;
t1261 = t1270 ^ 2;
t1290 = t1261 * t1256;
t1289 = t1266 * t1257;
t1288 = qJD(5) - t1258;
t1246 = t1268 * g(1) - t1272 * g(2);
t1277 = qJDD(1) * pkin(1) + t1246;
t1247 = -t1272 * g(1) - t1268 * g(2);
t1274 = qJD(1) ^ 2;
t1278 = -t1274 * pkin(1) + t1247;
t1211 = t1263 * t1277 + t1264 * t1278;
t1209 = -t1274 * pkin(2) + t1211;
t1210 = -t1263 * t1278 + t1264 * t1277;
t1275 = qJDD(1) * pkin(2) + t1210;
t1189 = t1271 * t1209 + t1267 * t1275;
t1179 = -t1256 * pkin(3) + t1257 * pkin(7) + t1189;
t1262 = -g(3) + qJDD(2);
t1174 = t1270 * t1179 + t1266 * t1262;
t1260 = t1266 ^ 2;
t1287 = t1260 + t1261;
t1245 = t1270 * t1256 * t1266;
t1235 = qJDD(4) + t1245;
t1286 = -qJDD(4) - qJDD(5);
t1285 = qJD(4) * t1292;
t1284 = qJD(4) * t1291;
t1173 = -t1266 * t1179 + t1270 * t1262;
t1188 = -t1267 * t1209 + t1271 * t1275;
t1224 = t1284 + t1289;
t1250 = t1270 * t1257;
t1281 = -t1250 + t1285;
t1283 = -t1265 * t1224 - t1269 * t1281;
t1238 = -t1263 * qJDD(1) - t1264 * t1274;
t1239 = t1264 * qJDD(1) - t1263 * t1274;
t1282 = t1272 * t1238 - t1268 * t1239;
t1280 = t1268 * t1238 + t1272 * t1239;
t1178 = -t1257 * pkin(3) - t1256 * pkin(7) - t1188;
t1276 = -t1269 * t1224 + t1265 * t1281;
t1273 = qJD(4) ^ 2;
t1243 = -t1273 - t1290;
t1242 = -t1260 * t1256 - t1273;
t1241 = -t1268 * qJDD(1) - t1272 * t1274;
t1240 = t1272 * qJDD(1) - t1268 * t1274;
t1237 = qJD(4) * pkin(4) - pkin(8) * t1292;
t1236 = -qJDD(4) + t1245;
t1233 = t1287 * t1256;
t1228 = t1287 * t1257;
t1225 = t1250 - 0.2e1 * t1285;
t1223 = 0.2e1 * t1284 + t1289;
t1216 = -t1294 - t1295;
t1215 = t1270 * t1236 - t1266 * t1242;
t1214 = -t1266 * t1235 + t1270 * t1243;
t1213 = t1266 * t1236 + t1270 * t1242;
t1212 = t1270 * t1235 + t1266 * t1243;
t1208 = t1271 * t1228 - t1267 * t1233;
t1207 = t1267 * t1228 + t1271 * t1233;
t1199 = t1286 - t1293;
t1198 = -t1286 - t1293;
t1197 = -t1294 - t1296;
t1196 = t1271 * t1215 + t1267 * t1223;
t1195 = t1271 * t1214 - t1267 * t1225;
t1194 = t1267 * t1215 - t1271 * t1223;
t1193 = t1267 * t1214 + t1271 * t1225;
t1192 = -t1295 - t1296;
t1191 = -t1263 * t1210 + t1264 * t1211;
t1190 = t1264 * t1210 + t1263 * t1211;
t1187 = t1269 * t1199 - t1265 * t1216;
t1186 = t1265 * t1199 + t1269 * t1216;
t1185 = -t1263 * t1207 + t1264 * t1208;
t1184 = t1264 * t1207 + t1263 * t1208;
t1183 = t1288 * t1218 + t1276;
t1182 = -t1297 * t1218 - t1276;
t1181 = -t1288 * t1220 + t1283;
t1180 = t1297 * t1220 - t1283;
t1176 = t1269 * t1197 - t1265 * t1198;
t1175 = t1265 * t1197 + t1269 * t1198;
t1172 = -t1263 * t1194 + t1264 * t1196;
t1171 = -t1263 * t1193 + t1264 * t1195;
t1170 = t1264 * t1194 + t1263 * t1196;
t1169 = t1264 * t1193 + t1263 * t1195;
t1168 = t1281 * pkin(4) - pkin(8) * t1290 + t1237 * t1292 + t1178;
t1167 = -pkin(4) * t1290 - t1281 * pkin(8) - qJD(4) * t1237 + t1174;
t1166 = (-t1224 + t1284) * pkin(8) + t1235 * pkin(4) + t1173;
t1165 = -t1267 * t1188 + t1271 * t1189;
t1164 = t1271 * t1188 + t1267 * t1189;
t1163 = -t1266 * t1186 + t1270 * t1187;
t1162 = t1270 * t1186 + t1266 * t1187;
t1161 = t1269 * t1181 - t1265 * t1183;
t1160 = t1265 * t1181 + t1269 * t1183;
t1159 = -t1266 * t1175 + t1270 * t1176;
t1158 = t1270 * t1175 + t1266 * t1176;
t1157 = -t1266 * t1173 + t1270 * t1174;
t1156 = t1270 * t1173 + t1266 * t1174;
t1155 = t1271 * t1163 + t1267 * t1182;
t1154 = t1267 * t1163 - t1271 * t1182;
t1153 = t1265 * t1166 + t1269 * t1167;
t1152 = t1269 * t1166 - t1265 * t1167;
t1151 = t1271 * t1159 + t1267 * t1180;
t1150 = t1267 * t1159 - t1271 * t1180;
t1149 = t1271 * t1157 + t1267 * t1178;
t1148 = t1267 * t1157 - t1271 * t1178;
t1147 = -t1263 * t1164 + t1264 * t1165;
t1146 = t1264 * t1164 + t1263 * t1165;
t1145 = -t1266 * t1160 + t1270 * t1161;
t1144 = t1270 * t1160 + t1266 * t1161;
t1143 = t1271 * t1145 + t1267 * t1192;
t1142 = t1267 * t1145 - t1271 * t1192;
t1141 = -t1263 * t1154 + t1264 * t1155;
t1140 = t1264 * t1154 + t1263 * t1155;
t1139 = -t1265 * t1152 + t1269 * t1153;
t1138 = t1269 * t1152 + t1265 * t1153;
t1137 = -t1263 * t1150 + t1264 * t1151;
t1136 = t1264 * t1150 + t1263 * t1151;
t1135 = -t1263 * t1148 + t1264 * t1149;
t1134 = t1264 * t1148 + t1263 * t1149;
t1133 = -t1263 * t1142 + t1264 * t1143;
t1132 = t1264 * t1142 + t1263 * t1143;
t1131 = -t1266 * t1138 + t1270 * t1139;
t1130 = t1270 * t1138 + t1266 * t1139;
t1129 = t1271 * t1131 + t1267 * t1168;
t1128 = t1267 * t1131 - t1271 * t1168;
t1127 = -t1263 * t1128 + t1264 * t1129;
t1126 = t1264 * t1128 + t1263 * t1129;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1241, -t1240, 0, -t1268 * t1246 + t1272 * t1247, 0, 0, 0, 0, 0, 0, t1282, -t1280, 0, -t1268 * t1190 + t1272 * t1191, 0, 0, 0, 0, 0, 0, t1304, t1303, 0, -t1268 * t1146 + t1272 * t1147, 0, 0, 0, 0, 0, 0, -t1268 * t1169 + t1272 * t1171, -t1268 * t1170 + t1272 * t1172, -t1268 * t1184 + t1272 * t1185, -t1268 * t1134 + t1272 * t1135, 0, 0, 0, 0, 0, 0, -t1268 * t1136 + t1272 * t1137, -t1268 * t1140 + t1272 * t1141, -t1268 * t1132 + t1272 * t1133, -t1268 * t1126 + t1272 * t1127; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1240, t1241, 0, t1272 * t1246 + t1268 * t1247, 0, 0, 0, 0, 0, 0, t1280, t1282, 0, t1272 * t1190 + t1268 * t1191, 0, 0, 0, 0, 0, 0, -t1303, t1304, 0, t1272 * t1146 + t1268 * t1147, 0, 0, 0, 0, 0, 0, t1272 * t1169 + t1268 * t1171, t1272 * t1170 + t1268 * t1172, t1272 * t1184 + t1268 * t1185, t1272 * t1134 + t1268 * t1135, 0, 0, 0, 0, 0, 0, t1272 * t1136 + t1268 * t1137, t1272 * t1140 + t1268 * t1141, t1272 * t1132 + t1268 * t1133, t1272 * t1126 + t1268 * t1127; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, 0, 0, 0, 0, 0, 0, t1212, t1213, 0, t1156, 0, 0, 0, 0, 0, 0, t1158, t1162, t1144, t1130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1274, -qJDD(1), 0, t1247, 0, 0, 0, 0, 0, 0, t1238, -t1239, 0, t1191, 0, 0, 0, 0, 0, 0, t1300, t1205, 0, t1147, 0, 0, 0, 0, 0, 0, t1171, t1172, t1185, t1135, 0, 0, 0, 0, 0, 0, t1137, t1141, t1133, t1127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1274, 0, t1246, 0, 0, 0, 0, 0, 0, t1239, t1238, 0, t1190, 0, 0, 0, 0, 0, 0, -t1205, t1300, 0, t1146, 0, 0, 0, 0, 0, 0, t1169, t1170, t1184, t1134, 0, 0, 0, 0, 0, 0, t1136, t1140, t1132, t1126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, 0, 0, 0, 0, 0, 0, t1212, t1213, 0, t1156, 0, 0, 0, 0, 0, 0, t1158, t1162, t1144, t1130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1274, -qJDD(1), 0, t1211, 0, 0, 0, 0, 0, 0, t1279, t1231, 0, t1165, 0, 0, 0, 0, 0, 0, t1195, t1196, t1208, t1149, 0, 0, 0, 0, 0, 0, t1151, t1155, t1143, t1129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1274, 0, t1210, 0, 0, 0, 0, 0, 0, -t1231, t1279, 0, t1164, 0, 0, 0, 0, 0, 0, t1193, t1194, t1207, t1148, 0, 0, 0, 0, 0, 0, t1150, t1154, t1142, t1128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, 0, 0, 0, 0, 0, 0, t1212, t1213, 0, t1156, 0, 0, 0, 0, 0, 0, t1158, t1162, t1144, t1130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1256, -t1257, 0, t1189, 0, 0, 0, 0, 0, 0, t1214, t1215, t1228, t1157, 0, 0, 0, 0, 0, 0, t1159, t1163, t1145, t1131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1257, -t1256, 0, t1188, 0, 0, 0, 0, 0, 0, t1225, -t1223, t1233, -t1178, 0, 0, 0, 0, 0, 0, -t1180, -t1182, -t1192, -t1168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1262, 0, 0, 0, 0, 0, 0, t1212, t1213, 0, t1156, 0, 0, 0, 0, 0, 0, t1158, t1162, t1144, t1130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1243, t1236, t1250, t1174, 0, 0, 0, 0, 0, 0, t1176, t1187, t1161, t1139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1235, t1242, -t1289, t1173, 0, 0, 0, 0, 0, 0, t1175, t1186, t1160, t1138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1225, t1223, -t1233, t1178, 0, 0, 0, 0, 0, 0, t1180, t1182, t1192, t1168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1197, t1199, t1181, t1153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1198, t1216, t1183, t1152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1180, t1182, t1192, t1168;];
f_new_reg = t1;
