% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PRRPR5
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PRRPR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:28:29
% EndTime: 2019-12-05 16:28:49
% DurationCPUTime: 21.01s
% Computational Cost: add. (114673->757), mult. (248550->1152), div. (0->0), fcn. (180727->12), ass. (0->543)
t1228 = sin(pkin(10));
t1231 = cos(pkin(10));
t1240 = cos(qJ(3));
t1349 = qJD(2) * t1240;
t1237 = sin(qJ(3));
t1350 = qJD(2) * t1237;
t1184 = t1228 * t1350 - t1231 * t1349;
t1186 = t1228 * t1349 + t1231 * t1350;
t1141 = t1186 * t1184;
t1365 = qJDD(3) - t1141;
t1373 = t1228 * t1365;
t1372 = t1231 * t1365;
t1236 = sin(qJ(5));
t1239 = cos(qJ(5));
t1151 = -t1239 * qJD(3) + t1236 * t1186;
t1153 = t1236 * qJD(3) + t1239 * t1186;
t1106 = t1153 * t1151;
t1216 = qJD(3) * t1349;
t1322 = t1237 * qJDD(2);
t1194 = t1216 + t1322;
t1310 = qJD(3) * t1350;
t1320 = t1240 * qJDD(2);
t1261 = -t1310 + t1320;
t1304 = t1228 * t1194 - t1231 * t1261;
t1140 = qJDD(5) + t1304;
t1366 = -t1106 + t1140;
t1371 = t1236 * t1366;
t1370 = t1239 * t1366;
t1230 = sin(pkin(5));
t1233 = cos(pkin(5));
t1229 = sin(pkin(9));
t1232 = cos(pkin(9));
t1308 = t1229 * g(1) - t1232 * g(2);
t1355 = g(3) - qJDD(1);
t1369 = -t1230 * t1355 + t1233 * t1308;
t1343 = t1186 * qJD(3);
t1107 = t1304 + t1343;
t1128 = t1184 * pkin(4) - t1186 * pkin(8);
t1242 = qJD(3) ^ 2;
t1348 = qJD(4) * t1184;
t1172 = -0.2e1 * t1348;
t1168 = t1230 * t1308 + t1233 * t1355;
t1201 = t1232 * g(1) + t1229 * g(2);
t1238 = sin(qJ(2));
t1241 = cos(qJ(2));
t1125 = -t1241 * t1201 + t1369 * t1238;
t1362 = qJD(2) ^ 2;
t1244 = -t1362 * pkin(2) + qJDD(2) * pkin(7) + t1125;
t1078 = -t1237 * t1168 + t1240 * t1244;
t1202 = qJD(3) * pkin(3) - qJ(4) * t1350;
t1361 = t1240 ^ 2;
t1223 = t1361 * t1362;
t1042 = -pkin(3) * t1223 + qJ(4) * t1261 - qJD(3) * t1202 + t1078;
t1077 = t1240 * t1168 + t1237 * t1244;
t1213 = t1237 * t1362 * t1240;
t1203 = qJDD(3) + t1213;
t1243 = -t1077 + (-t1194 + t1216) * qJ(4) + t1203 * pkin(3);
t1325 = t1231 * t1042 + t1228 * t1243;
t960 = t1172 + t1325;
t938 = -t1242 * pkin(4) + qJDD(3) * pkin(8) - t1184 * t1128 + t960;
t1124 = -t1238 * t1201 - t1369 * t1241;
t1115 = -qJDD(2) * pkin(2) - t1362 * pkin(7) + t1124;
t1065 = -t1261 * pkin(3) - qJ(4) * t1223 + t1202 * t1350 + qJDD(4) + t1115;
t1143 = t1231 * t1194 + t1228 * t1261;
t1175 = t1184 * qJD(3);
t1110 = -t1175 + t1143;
t984 = t1107 * pkin(4) - t1110 * pkin(8) + t1065;
t898 = t1236 * t938 - t1239 * t984;
t899 = t1236 * t984 + t1239 * t938;
t856 = t1236 * t898 + t1239 * t899;
t1368 = t1229 * t1355;
t1367 = t1232 * t1355;
t1364 = -t1232 * t1201 - t1229 * t1308;
t1363 = -t1229 * t1201 + t1232 * t1308;
t1149 = t1151 ^ 2;
t1150 = t1153 ^ 2;
t1177 = qJD(5) + t1184;
t1176 = t1177 ^ 2;
t1178 = t1184 ^ 2;
t1179 = t1186 ^ 2;
t1360 = 0.2e1 * qJD(4);
t1306 = t1228 * t1042 - t1231 * t1243;
t959 = t1186 * t1360 + t1306;
t901 = t1228 * t960 - t1231 * t959;
t1359 = pkin(3) * t901;
t1109 = -t1304 + t1343;
t1111 = t1175 + t1143;
t1039 = t1228 * t1109 - t1231 * t1111;
t1358 = pkin(3) * t1039;
t1357 = pkin(4) * t1228;
t1061 = t1238 * t1124 + t1241 * t1125;
t1356 = pkin(6) * t1061;
t937 = -qJDD(3) * pkin(4) - t1242 * pkin(8) + (t1360 + t1128) * t1186 + t1306;
t1354 = -pkin(4) * t937 + pkin(8) * t856;
t934 = t1236 * t937;
t1353 = t1237 * t901;
t935 = t1239 * t937;
t1352 = t1240 * t901;
t1351 = qJD(2) * qJD(3);
t1347 = t1177 * t1236;
t1346 = t1177 * t1239;
t1345 = t1184 * t1228;
t1344 = t1184 * t1231;
t1342 = t1186 * t1228;
t1341 = t1186 * t1231;
t1225 = t1237 ^ 2;
t1340 = t1225 * t1362;
t1339 = t1228 * t1065;
t1133 = qJDD(3) + t1141;
t1338 = t1228 * t1133;
t1336 = t1231 * t1065;
t1335 = t1231 * t1133;
t1063 = t1106 + t1140;
t1334 = t1236 * t1063;
t1333 = t1237 * t1115;
t1332 = t1237 * t1203;
t1204 = qJDD(3) - t1213;
t1331 = t1237 * t1204;
t1330 = t1238 * t1168;
t1329 = t1239 * t1063;
t1328 = t1240 * t1115;
t1195 = -0.2e1 * t1310 + t1320;
t1154 = t1240 * t1195;
t1327 = t1240 * t1204;
t1326 = t1241 * t1168;
t1324 = qJDD(3) * t1241;
t1323 = t1230 * qJDD(2);
t1321 = t1238 * qJDD(2);
t1319 = t1241 * qJDD(2);
t1318 = t1225 + t1361;
t1305 = t1239 * qJDD(3) - t1236 * t1143;
t1089 = -t1153 * qJD(5) + t1305;
t1123 = t1153 * t1177;
t1031 = t1089 - t1123;
t1088 = -t1176 - t1149;
t997 = t1239 * t1088 - t1371;
t1317 = pkin(4) * t1031 + pkin(8) * t997 - t935;
t1316 = -pkin(4) * t1231 - pkin(3);
t1099 = -t1150 - t1176;
t1003 = -t1236 * t1099 - t1329;
t1266 = -t1236 * qJDD(3) - t1239 * t1143;
t1034 = (qJD(5) + t1177) * t1151 + t1266;
t1315 = pkin(4) * t1034 + pkin(8) * t1003 + t934;
t1314 = t1241 * t1141;
t1313 = t1228 * t1106;
t1312 = t1231 * t1106;
t1311 = t1151 * t1347;
t1309 = t1238 * t1141;
t902 = t1228 * t959 + t1231 * t960;
t1000 = t1237 * t1077 + t1240 * t1078;
t1074 = t1149 + t1150;
t1030 = (-qJD(5) + t1177) * t1153 + t1305;
t1090 = -t1151 * qJD(5) - t1266;
t1122 = t1177 * t1151;
t1033 = t1090 + t1122;
t965 = t1239 * t1030 + t1236 * t1033;
t1303 = pkin(4) * t1074 + pkin(8) * t965 + t856;
t841 = t1228 * t856 - t1231 * t937;
t1302 = pkin(3) * t841 + t1354;
t1301 = t1241 * t1213;
t1300 = t1238 * t1213;
t1171 = -t1179 - t1242;
t1095 = t1231 * t1171 - t1338;
t1299 = pkin(3) * t1095 - t1325;
t1196 = t1318 * qJDD(2);
t1199 = t1223 + t1340;
t1146 = t1241 * t1196 - t1238 * t1199;
t999 = t1240 * t1077 - t1237 * t1078;
t1298 = pkin(6) * t1146 + t1238 * t999;
t1197 = -t1238 * t1362 + t1319;
t1296 = -pkin(6) * t1197 - t1330;
t1268 = t1241 * t1362 + t1321;
t1295 = -pkin(6) * t1268 + t1326;
t855 = t1236 * t899 - t1239 * t898;
t842 = t1228 * t937 + t1231 * t856;
t812 = -t1237 * t841 + t1240 * t842;
t1294 = t1238 * t812 - t1241 * t855;
t932 = t1231 * t1074 + t1228 * t965;
t933 = -t1228 * t1074 + t1231 * t965;
t875 = -t1237 * t932 + t1240 * t933;
t963 = t1236 * t1030 - t1239 * t1033;
t1293 = t1238 * t875 - t1241 * t963;
t1104 = t1150 - t1149;
t1032 = t1090 - t1122;
t966 = t1239 * t1031 - t1236 * t1032;
t943 = -t1231 * t1104 + t1228 * t966;
t944 = t1228 * t1104 + t1231 * t966;
t879 = -t1237 * t943 + t1240 * t944;
t964 = t1236 * t1031 + t1239 * t1032;
t1292 = t1238 * t879 - t1241 * t964;
t946 = t1231 * t1031 + t1228 * t997;
t947 = -t1228 * t1031 + t1231 * t997;
t882 = -t1237 * t946 + t1240 * t947;
t996 = t1236 * t1088 + t1370;
t1291 = t1238 * t882 - t1241 * t996;
t1290 = pkin(3) * t946 + t1317;
t1002 = t1239 * t1099 - t1334;
t950 = t1228 * t1003 + t1231 * t1034;
t951 = t1231 * t1003 - t1228 * t1034;
t885 = -t1237 * t950 + t1240 * t951;
t1289 = -t1002 * t1241 + t1238 * t885;
t1117 = -t1150 + t1176;
t1005 = t1239 * t1117 + t1371;
t1007 = -t1236 * t1117 + t1370;
t953 = t1228 * t1007 - t1231 * t1033;
t955 = t1231 * t1007 + t1228 * t1033;
t889 = -t1237 * t953 + t1240 * t955;
t1288 = -t1005 * t1241 + t1238 * t889;
t1116 = t1149 - t1176;
t1006 = t1236 * t1116 + t1329;
t1008 = t1239 * t1116 - t1334;
t1029 = -t1089 - t1123;
t954 = t1228 * t1008 + t1231 * t1029;
t956 = t1231 * t1008 - t1228 * t1029;
t890 = -t1237 * t954 + t1240 * t956;
t1287 = -t1006 * t1241 + t1238 * t890;
t1021 = t1239 * t1089 + t1311;
t1022 = -t1236 * t1089 + t1151 * t1346;
t986 = t1228 * t1022 + t1312;
t988 = t1231 * t1022 - t1313;
t921 = -t1237 * t986 + t1240 * t988;
t1286 = -t1021 * t1241 + t1238 * t921;
t1114 = t1153 * t1346;
t1023 = -t1236 * t1090 - t1114;
t1024 = t1239 * t1090 - t1153 * t1347;
t987 = t1228 * t1024 - t1312;
t989 = t1231 * t1024 + t1313;
t922 = -t1237 * t987 + t1240 * t989;
t1285 = t1023 * t1241 + t1238 * t922;
t1052 = -t1114 - t1311;
t1054 = (-t1151 * t1239 + t1153 * t1236) * t1177;
t1009 = t1228 * t1054 - t1231 * t1140;
t1010 = t1231 * t1054 + t1228 * t1140;
t949 = -t1237 * t1009 + t1240 * t1010;
t1284 = -t1052 * t1241 + t1238 * t949;
t858 = t1240 * t902 - t1353;
t1283 = -t1065 * t1241 + t1238 * t858;
t1105 = -t1178 - t1179;
t1041 = t1231 * t1109 + t1228 * t1111;
t970 = -t1237 * t1039 + t1240 * t1041;
t1282 = -t1105 * t1241 + t1238 * t970;
t1131 = -t1242 - t1178;
t1066 = t1228 * t1131 + t1372;
t1067 = t1231 * t1131 - t1373;
t992 = -t1237 * t1066 + t1240 * t1067;
t1281 = -t1107 * t1241 + t1238 * t992;
t1137 = t1179 - t1178;
t1038 = -t1228 * t1107 + t1231 * t1110;
t1040 = -t1231 * t1107 - t1228 * t1110;
t969 = -t1237 * t1038 + t1240 * t1040;
t1280 = -t1137 * t1241 + t1238 * t969;
t1279 = pkin(3) * t950 + t1315;
t1278 = t1000 * t1238 - t1115 * t1241;
t1169 = t1178 - t1242;
t1093 = t1228 * t1169 + t1335;
t1096 = t1231 * t1169 - t1338;
t1014 = -t1237 * t1093 + t1240 * t1096;
t1277 = t1014 * t1238 - t1109 * t1241;
t1170 = -t1179 + t1242;
t1094 = t1231 * t1170 + t1373;
t1097 = -t1228 * t1170 + t1372;
t1015 = -t1237 * t1094 + t1240 * t1097;
t1276 = t1015 * t1238 - t1111 * t1241;
t1098 = -t1228 * t1171 - t1335;
t1016 = -t1237 * t1095 + t1240 * t1098;
t1275 = t1016 * t1238 - t1110 * t1241;
t1060 = t1241 * t1124 - t1238 * t1125;
t1193 = 0.2e1 * t1216 + t1322;
t1145 = -t1237 * t1193 + t1154;
t1200 = -t1223 + t1340;
t1274 = t1145 * t1238 - t1200 * t1241;
t1212 = -t1223 - t1242;
t1161 = t1240 * t1212 - t1332;
t1273 = t1161 * t1238 + t1195 * t1241;
t1210 = -t1242 - t1340;
t1163 = -t1237 * t1210 - t1327;
t1272 = t1163 * t1238 - t1193 * t1241;
t1182 = t1268 * t1233;
t1271 = t1232 * t1182 + t1229 * t1197;
t1270 = t1229 * t1182 - t1232 * t1197;
t1269 = t1196 * t1238 + t1199 * t1241;
t1120 = (-t1341 - t1345) * qJD(3);
t1121 = (t1342 - t1344) * qJD(3);
t1051 = -t1237 * t1120 + t1240 * t1121;
t1267 = t1051 * t1238 - t1324;
t1191 = t1318 * t1351;
t1265 = t1191 * t1238 - t1324;
t1264 = pkin(3) * t932 + t1303;
t1100 = qJD(3) * t1345 - t1231 * t1304;
t1101 = qJD(3) * t1344 + t1228 * t1304;
t1019 = -t1237 * t1100 + t1240 * t1101;
t1263 = t1019 * t1238 + t1314;
t1102 = qJD(3) * t1341 + t1228 * t1143;
t1103 = -qJD(3) * t1342 + t1231 * t1143;
t1020 = -t1237 * t1102 + t1240 * t1103;
t1262 = t1020 * t1238 - t1314;
t1211 = t1223 - t1242;
t1160 = t1240 * t1211 - t1331;
t1260 = t1160 * t1238 - t1240 * t1319;
t1192 = t1240 * t1203;
t1209 = t1242 - t1340;
t1162 = -t1237 * t1209 + t1192;
t1259 = t1162 * t1238 - t1237 * t1319;
t801 = qJ(4) * t842 + (-pkin(8) * t1228 + t1316) * t855;
t809 = -qJ(4) * t841 + (-pkin(8) * t1231 + t1357) * t855;
t811 = t1237 * t842 + t1240 * t841;
t789 = -pkin(7) * t811 - t1237 * t801 + t1240 * t809;
t798 = -pkin(2) * t811 - t1302;
t806 = t1238 * t855 + t1241 * t812;
t1258 = pkin(6) * t806 + t1238 * t789 + t1241 * t798;
t845 = -pkin(8) * t963 - t855;
t822 = qJ(4) * t933 + t1228 * t845 + t1316 * t963;
t824 = -qJ(4) * t932 + t1231 * t845 + t963 * t1357;
t874 = t1237 * t933 + t1240 * t932;
t800 = -pkin(7) * t874 - t1237 * t822 + t1240 * t824;
t815 = -pkin(2) * t874 - t1264;
t859 = t1238 * t963 + t1241 * t875;
t1257 = pkin(6) * t859 + t1238 * t800 + t1241 * t815;
t871 = -pkin(4) * t996 + t898;
t905 = -pkin(8) * t996 + t934;
t826 = -pkin(3) * t996 + qJ(4) * t947 + t1228 * t905 + t1231 * t871;
t836 = -qJ(4) * t946 - t1228 * t871 + t1231 * t905;
t881 = t1237 * t947 + t1240 * t946;
t807 = -pkin(7) * t881 - t1237 * t826 + t1240 * t836;
t833 = -pkin(2) * t881 - t1290;
t865 = t1238 * t996 + t1241 * t882;
t1256 = pkin(6) * t865 + t1238 * t807 + t1241 * t833;
t873 = -pkin(4) * t1002 + t899;
t906 = -pkin(8) * t1002 + t935;
t827 = -pkin(3) * t1002 + qJ(4) * t951 + t1228 * t906 + t1231 * t873;
t838 = -qJ(4) * t950 - t1228 * t873 + t1231 * t906;
t884 = t1237 * t951 + t1240 * t950;
t808 = -pkin(7) * t884 - t1237 * t827 + t1240 * t838;
t837 = -pkin(2) * t884 - t1279;
t866 = t1238 * t1002 + t1241 * t885;
t1255 = pkin(6) * t866 + t1238 * t808 + t1241 * t837;
t857 = t1237 * t902 + t1352;
t880 = -pkin(3) * t1065 + qJ(4) * t902;
t819 = -pkin(7) * t857 - qJ(4) * t1352 - t1237 * t880;
t832 = -pkin(2) * t857 - t1359;
t853 = t1238 * t1065 + t1241 * t858;
t1254 = pkin(6) * t853 + t1238 * t819 + t1241 * t832;
t872 = -pkin(3) * t1105 + qJ(4) * t1041 + t902;
t877 = -qJ(4) * t1039 - t901;
t968 = t1240 * t1039 + t1237 * t1041;
t831 = -pkin(7) * t968 - t1237 * t872 + t1240 * t877;
t925 = -pkin(2) * t968 - t1358;
t945 = t1238 * t1105 + t1241 * t970;
t1253 = pkin(6) * t945 + t1238 * t831 + t1241 * t925;
t972 = -pkin(3) * t1107 + qJ(4) * t1067 - t1336;
t990 = -qJ(4) * t1066 + t1339;
t991 = t1240 * t1066 + t1237 * t1067;
t883 = -pkin(7) * t991 - t1237 * t972 + t1240 * t990;
t1245 = pkin(3) * t1066 - t959;
t900 = -pkin(2) * t991 - t1245;
t971 = t1238 * t1107 + t1241 * t992;
t1252 = pkin(6) * t971 + t1238 * t883 + t1241 * t900;
t1004 = -qJ(4) * t1095 + t1336;
t1013 = t1240 * t1095 + t1237 * t1098;
t979 = -pkin(3) * t1110 + qJ(4) * t1098 + t1339;
t895 = -pkin(7) * t1013 + t1240 * t1004 - t1237 * t979;
t907 = -pkin(2) * t1013 + t1172 - t1299;
t982 = t1241 * t1016 + t1238 * t1110;
t1251 = pkin(6) * t982 + t1238 * t895 + t1241 * t907;
t1165 = -t1237 * t1261 - t1361 * t1351;
t1250 = t1165 * t1238 - t1301;
t1166 = t1240 * t1194 - t1225 * t1351;
t1249 = t1166 * t1238 + t1301;
t1157 = t1237 * t1212 + t1192;
t1048 = -pkin(2) * t1157 + t1077;
t1075 = -pkin(7) * t1157 + t1333;
t1118 = t1241 * t1161 - t1238 * t1195;
t1248 = pkin(6) * t1118 + t1048 * t1241 + t1075 * t1238;
t1159 = t1240 * t1210 - t1331;
t1049 = -pkin(2) * t1159 + t1078;
t1076 = -pkin(7) * t1159 + t1328;
t1119 = t1241 * t1163 + t1238 * t1193;
t1247 = pkin(6) * t1119 + t1049 * t1241 + t1076 * t1238;
t975 = t1241 * t1000 + t1238 * t1115;
t1246 = pkin(6) * t975 - (-pkin(2) * t1241 - pkin(7) * t1238) * t999;
t1219 = t1238 * qJDD(3);
t1218 = t1233 * qJDD(2);
t1183 = t1197 * t1233;
t1181 = t1197 * t1230;
t1180 = t1268 * t1230;
t1167 = t1241 * t1191 + t1219;
t1158 = t1240 * t1209 + t1332;
t1156 = t1237 * t1211 + t1327;
t1155 = (t1194 + t1216) * t1237;
t1148 = t1265 * t1233;
t1147 = t1265 * t1230;
t1144 = t1240 * t1193 + t1237 * t1195;
t1139 = t1269 * t1233;
t1138 = t1269 * t1230;
t1136 = -t1229 * t1183 - t1232 * t1268;
t1135 = t1232 * t1183 - t1229 * t1268;
t1130 = t1241 * t1166 - t1300;
t1129 = t1241 * t1165 + t1300;
t1127 = t1241 * t1162 + t1237 * t1321;
t1126 = t1241 * t1160 + t1238 * t1320;
t1113 = t1241 * t1145 + t1238 * t1200;
t1092 = -t1326 + (t1180 * t1230 + t1182 * t1233) * pkin(6);
t1091 = -t1330 + (-t1181 * t1230 - t1183 * t1233) * pkin(6);
t1087 = -t1230 * t1155 + t1233 * t1249;
t1086 = -t1230 * t1154 + t1233 * t1250;
t1085 = t1233 * t1155 + t1230 * t1249;
t1084 = t1233 * t1154 + t1230 * t1250;
t1083 = -t1230 * t1158 + t1233 * t1259;
t1082 = -t1230 * t1156 + t1233 * t1260;
t1081 = t1233 * t1158 + t1230 * t1259;
t1080 = t1233 * t1156 + t1230 * t1260;
t1072 = -t1230 * t1159 + t1233 * t1272;
t1071 = -t1230 * t1157 + t1233 * t1273;
t1070 = t1233 * t1159 + t1230 * t1272;
t1069 = t1233 * t1157 + t1230 * t1273;
t1059 = -t1230 * t1144 + t1233 * t1274;
t1058 = t1233 * t1144 + t1230 * t1274;
t1057 = pkin(2) * t1195 + pkin(7) * t1161 - t1328;
t1056 = -pkin(2) * t1193 + pkin(7) * t1163 + t1333;
t1055 = t1061 * t1233;
t1053 = t1061 * t1230;
t1050 = t1240 * t1120 + t1237 * t1121;
t1047 = t1241 * t1051 + t1219;
t1046 = -pkin(1) * t1181 + t1230 * t1124 + t1233 * t1295;
t1045 = pkin(1) * t1180 + t1230 * t1125 + t1233 * t1296;
t1044 = pkin(1) * t1183 - t1233 * t1124 + t1230 * t1295;
t1043 = -pkin(1) * t1182 - t1233 * t1125 + t1230 * t1296;
t1026 = -t1060 * t1233 + t1230 * t1168;
t1025 = -t1060 * t1230 - t1233 * t1168;
t1018 = t1240 * t1102 + t1237 * t1103;
t1017 = t1240 * t1100 + t1237 * t1101;
t1012 = t1240 * t1094 + t1237 * t1097;
t1011 = t1240 * t1093 + t1237 * t1096;
t995 = t1241 * t1020 + t1309;
t994 = t1241 * t1019 - t1309;
t985 = pkin(2) * t1199 + pkin(7) * t1196 + t1000;
t981 = t1241 * t1015 + t1238 * t1111;
t980 = t1241 * t1014 + t1109 * t1238;
t978 = -pkin(2) * t1115 + pkin(7) * t1000;
t977 = -t1230 * t1050 + t1233 * t1267;
t976 = t1233 * t1050 + t1230 * t1267;
t974 = pkin(1) * t1026 + t1230 * t1356;
t973 = -pkin(1) * t1025 + t1233 * t1356;
t967 = t1240 * t1038 + t1237 * t1040;
t961 = t1241 * t999 + (-t1138 * t1230 - t1139 * t1233) * pkin(6);
t957 = (-t1025 * t1230 - t1026 * t1233) * pkin(6);
t952 = t1238 * t1137 + t1241 * t969;
t948 = t1240 * t1009 + t1237 * t1010;
t942 = -t1230 * t1018 + t1233 * t1262;
t941 = -t1230 * t1017 + t1233 * t1263;
t940 = t1233 * t1018 + t1230 * t1262;
t939 = t1233 * t1017 + t1230 * t1263;
t931 = -t1230 * t1013 + t1233 * t1275;
t930 = -t1230 * t1012 + t1233 * t1276;
t929 = -t1230 * t1011 + t1233 * t1277;
t928 = t1233 * t1013 + t1230 * t1275;
t927 = t1233 * t1012 + t1230 * t1276;
t926 = t1233 * t1011 + t1230 * t1277;
t924 = -t1238 * t1049 + t1241 * t1076 + (-t1070 * t1230 - t1072 * t1233) * pkin(6);
t923 = -t1238 * t1048 + t1241 * t1075 + (-t1069 * t1230 - t1071 * t1233) * pkin(6);
t920 = t1237 * t989 + t1240 * t987;
t919 = t1237 * t988 + t1240 * t986;
t918 = t1230 * t999 + t1233 * t1278;
t917 = t1230 * t1278 - t1233 * t999;
t916 = t1238 * t1052 + t1241 * t949;
t915 = -t1230 * t991 + t1233 * t1281;
t914 = t1230 * t1281 + t1233 * t991;
t913 = -pkin(1) * t1070 - t1230 * t1056 + t1233 * t1247;
t912 = -pkin(1) * t1069 - t1230 * t1057 + t1233 * t1248;
t911 = pkin(1) * t1072 + t1233 * t1056 + t1230 * t1247;
t910 = pkin(1) * t1071 + t1233 * t1057 + t1230 * t1248;
t909 = -pkin(1) * t1138 - t1230 * t985 + t1233 * t1298;
t908 = pkin(1) * t1139 + t1230 * t1298 + t1233 * t985;
t904 = -t1238 * t1023 + t1241 * t922;
t903 = t1238 * t1021 + t1241 * t921;
t894 = -t1230 * t967 + t1233 * t1280;
t893 = t1230 * t1280 + t1233 * t967;
t892 = -t1230 * t968 + t1233 * t1282;
t891 = t1230 * t1282 + t1233 * t968;
t888 = t1237 * t956 + t1240 * t954;
t887 = t1237 * t955 + t1240 * t953;
t886 = -pkin(2) * t1110 + pkin(7) * t1016 + t1237 * t1004 + t1240 * t979;
t878 = t1237 * t944 + t1240 * t943;
t876 = -pkin(2) * t1107 + pkin(7) * t992 + t1237 * t990 + t1240 * t972;
t870 = -t1230 * t948 + t1233 * t1284;
t869 = t1230 * t1284 + t1233 * t948;
t868 = t1238 * t1006 + t1241 * t890;
t867 = t1238 * t1005 + t1241 * t889;
t864 = t1238 * t964 + t1241 * t879;
t863 = -t1230 * t920 + t1233 * t1285;
t862 = -t1230 * t919 + t1233 * t1286;
t861 = t1230 * t1285 + t1233 * t920;
t860 = t1230 * t1286 + t1233 * t919;
t852 = -(pkin(2) * t1238 - pkin(7) * t1241) * t999 + (-t1230 * t917 - t1233 * t918) * pkin(6);
t851 = -pkin(1) * t917 - t1230 * t978 + t1233 * t1246;
t850 = pkin(1) * t918 + t1230 * t1246 + t1233 * t978;
t849 = -t1230 * t888 + t1233 * t1287;
t848 = -t1230 * t887 + t1233 * t1288;
t847 = t1230 * t1287 + t1233 * t888;
t846 = t1230 * t1288 + t1233 * t887;
t844 = -t1230 * t884 + t1233 * t1289;
t843 = t1230 * t1289 + t1233 * t884;
t840 = -t1230 * t881 + t1233 * t1291;
t839 = t1230 * t1291 + t1233 * t881;
t835 = -t1230 * t878 + t1233 * t1292;
t834 = t1230 * t1292 + t1233 * t878;
t830 = -t1230 * t874 + t1233 * t1293;
t829 = t1230 * t1293 + t1233 * t874;
t828 = -pkin(2) * t1105 + pkin(7) * t970 + t1237 * t877 + t1240 * t872;
t825 = -t1238 * t907 + t1241 * t895 + (-t1230 * t928 - t1233 * t931) * pkin(6);
t823 = -t1238 * t900 + t1241 * t883 + (-t1230 * t914 - t1233 * t915) * pkin(6);
t821 = -pkin(1) * t928 - t1230 * t886 + t1233 * t1251;
t820 = pkin(1) * t931 + t1230 * t1251 + t1233 * t886;
t818 = -t1230 * t857 + t1233 * t1283;
t817 = t1230 * t1283 + t1233 * t857;
t816 = -pkin(2) * t1065 + pkin(7) * t858 - qJ(4) * t1353 + t1240 * t880;
t814 = -pkin(1) * t914 - t1230 * t876 + t1233 * t1252;
t813 = pkin(1) * t915 + t1230 * t1252 + t1233 * t876;
t810 = -t1238 * t925 + t1241 * t831 + (-t1230 * t891 - t1233 * t892) * pkin(6);
t805 = -pkin(2) * t1002 + pkin(7) * t885 + t1237 * t838 + t1240 * t827;
t804 = -pkin(2) * t996 + pkin(7) * t882 + t1237 * t836 + t1240 * t826;
t803 = -pkin(1) * t891 - t1230 * t828 + t1233 * t1253;
t802 = pkin(1) * t892 + t1230 * t1253 + t1233 * t828;
t799 = -pkin(2) * t963 + pkin(7) * t875 + t1237 * t824 + t1240 * t822;
t797 = -t1230 * t811 + t1233 * t1294;
t796 = t1230 * t1294 + t1233 * t811;
t795 = -t1238 * t837 + t1241 * t808 + (-t1230 * t843 - t1233 * t844) * pkin(6);
t794 = -t1238 * t833 + t1241 * t807 + (-t1230 * t839 - t1233 * t840) * pkin(6);
t793 = -t1238 * t832 + t1241 * t819 + (-t1230 * t817 - t1233 * t818) * pkin(6);
t792 = -pkin(1) * t817 - t1230 * t816 + t1233 * t1254;
t791 = pkin(1) * t818 + t1230 * t1254 + t1233 * t816;
t790 = -t1238 * t815 + t1241 * t800 + (-t1230 * t829 - t1233 * t830) * pkin(6);
t788 = -pkin(1) * t843 - t1230 * t805 + t1233 * t1255;
t787 = pkin(1) * t844 + t1230 * t1255 + t1233 * t805;
t786 = -pkin(1) * t839 - t1230 * t804 + t1233 * t1256;
t785 = pkin(1) * t840 + t1230 * t1256 + t1233 * t804;
t784 = -pkin(2) * t855 + pkin(7) * t812 + t1237 * t809 + t1240 * t801;
t783 = -pkin(1) * t829 - t1230 * t799 + t1233 * t1257;
t782 = pkin(1) * t830 + t1230 * t1257 + t1233 * t799;
t781 = -t1238 * t798 + t1241 * t789 + (-t1230 * t796 - t1233 * t797) * pkin(6);
t780 = -pkin(1) * t796 - t1230 * t784 + t1233 * t1258;
t779 = pkin(1) * t797 + t1230 * t1258 + t1233 * t784;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t1368, -t1367, -t1363, -qJ(1) * t1363, 0, 0, -t1270, 0, t1136, t1229 * t1323, -qJ(1) * t1135 - t1229 * t1046 + t1232 * t1091, qJ(1) * t1271 - t1229 * t1045 + t1232 * t1092, -t1229 * t1055 + t1232 * t1060, t1232 * t957 - t1229 * t973 - qJ(1) * (t1232 * t1026 + t1229 * t1061), -t1229 * t1087 + t1232 * t1130, -t1229 * t1059 + t1232 * t1113, -t1229 * t1083 + t1232 * t1127, -t1229 * t1086 + t1232 * t1129, -t1229 * t1082 + t1232 * t1126, -t1229 * t1148 + t1232 * t1167, t1232 * t923 - t1229 * t912 - qJ(1) * (t1232 * t1071 + t1229 * t1118), t1232 * t924 - t1229 * t913 - qJ(1) * (t1232 * t1072 + t1229 * t1119), t1232 * t961 - t1229 * t909 - qJ(1) * (t1232 * t1139 + t1229 * t1146), t1232 * t852 - t1229 * t851 - qJ(1) * (t1229 * t975 + t1232 * t918), -t1229 * t942 + t1232 * t995, -t1229 * t894 + t1232 * t952, -t1229 * t930 + t1232 * t981, -t1229 * t941 + t1232 * t994, -t1229 * t929 + t1232 * t980, t1232 * t1047 - t1229 * t977, t1232 * t823 - t1229 * t814 - qJ(1) * (t1229 * t971 + t1232 * t915), t1232 * t825 - t1229 * t821 - qJ(1) * (t1229 * t982 + t1232 * t931), t1232 * t810 - t1229 * t803 - qJ(1) * (t1229 * t945 + t1232 * t892), t1232 * t793 - t1229 * t792 - qJ(1) * (t1229 * t853 + t1232 * t818), -t1229 * t863 + t1232 * t904, -t1229 * t835 + t1232 * t864, -t1229 * t848 + t1232 * t867, -t1229 * t862 + t1232 * t903, -t1229 * t849 + t1232 * t868, -t1229 * t870 + t1232 * t916, t1232 * t794 - t1229 * t786 - qJ(1) * (t1229 * t865 + t1232 * t840), t1232 * t795 - t1229 * t788 - qJ(1) * (t1229 * t866 + t1232 * t844), t1232 * t790 - t1229 * t783 - qJ(1) * (t1229 * t859 + t1232 * t830), t1232 * t781 - t1229 * t780 - qJ(1) * (t1229 * t806 + t1232 * t797); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t1367, -t1368, t1364, qJ(1) * t1364, 0, 0, t1271, 0, t1135, -t1232 * t1323, qJ(1) * t1136 + t1232 * t1046 + t1229 * t1091, qJ(1) * t1270 + t1232 * t1045 + t1229 * t1092, t1232 * t1055 + t1229 * t1060, t1229 * t957 + t1232 * t973 + qJ(1) * (-t1229 * t1026 + t1232 * t1061), t1232 * t1087 + t1229 * t1130, t1232 * t1059 + t1229 * t1113, t1232 * t1083 + t1229 * t1127, t1232 * t1086 + t1229 * t1129, t1232 * t1082 + t1229 * t1126, t1232 * t1148 + t1229 * t1167, t1229 * t923 + t1232 * t912 + qJ(1) * (-t1229 * t1071 + t1232 * t1118), t1229 * t924 + t1232 * t913 + qJ(1) * (-t1229 * t1072 + t1232 * t1119), t1229 * t961 + t1232 * t909 + qJ(1) * (-t1229 * t1139 + t1232 * t1146), t1229 * t852 + t1232 * t851 + qJ(1) * (-t1229 * t918 + t1232 * t975), t1229 * t995 + t1232 * t942, t1229 * t952 + t1232 * t894, t1229 * t981 + t1232 * t930, t1229 * t994 + t1232 * t941, t1229 * t980 + t1232 * t929, t1229 * t1047 + t1232 * t977, t1229 * t823 + t1232 * t814 + qJ(1) * (-t1229 * t915 + t1232 * t971), t1229 * t825 + t1232 * t821 + qJ(1) * (-t1229 * t931 + t1232 * t982), t1229 * t810 + t1232 * t803 + qJ(1) * (-t1229 * t892 + t1232 * t945), t1229 * t793 + t1232 * t792 + qJ(1) * (-t1229 * t818 + t1232 * t853), t1229 * t904 + t1232 * t863, t1229 * t864 + t1232 * t835, t1229 * t867 + t1232 * t848, t1229 * t903 + t1232 * t862, t1229 * t868 + t1232 * t849, t1229 * t916 + t1232 * t870, t1229 * t794 + t1232 * t786 + qJ(1) * (-t1229 * t840 + t1232 * t865), t1229 * t795 + t1232 * t788 + qJ(1) * (-t1229 * t844 + t1232 * t866), t1229 * t790 + t1232 * t783 + qJ(1) * (-t1229 * t830 + t1232 * t859), t1229 * t781 + t1232 * t780 + qJ(1) * (-t1229 * t797 + t1232 * t806); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t1308, t1201, 0, 0, 0, 0, t1180, 0, t1181, t1218, t1044, t1043, t1053, t974, t1085, t1058, t1081, t1084, t1080, t1147, t910, t911, t908, t850, t940, t893, t927, t939, t926, t976, t813, t820, t802, t791, t861, t834, t846, t860, t847, t869, t785, t787, t782, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1355, -t1308, 0, 0, 0, t1197, 0, -t1268, 0, t1091, t1092, t1060, t957, t1130, t1113, t1127, t1129, t1126, t1167, t923, t924, t961, t852, t995, t952, t981, t994, t980, t1047, t823, t825, t810, t793, t904, t864, t867, t903, t868, t916, t794, t795, t790, t781; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1355, 0, -t1201, 0, 0, 0, t1182, 0, t1183, -t1323, t1046, t1045, t1055, t973, t1087, t1059, t1083, t1086, t1082, t1148, t912, t913, t909, t851, t942, t894, t930, t941, t929, t977, t814, t821, t803, t792, t863, t835, t848, t862, t849, t870, t786, t788, t783, t780; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1308, t1201, 0, 0, 0, 0, t1180, 0, t1181, t1218, t1044, t1043, t1053, t974, t1085, t1058, t1081, t1084, t1080, t1147, t910, t911, t908, t850, t940, t893, t927, t939, t926, t976, t813, t820, t802, t791, t861, t834, t846, t860, t847, t869, t785, t787, t782, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t1362, 0, 0, -t1168, t1124, 0, t1166, t1145, t1162, t1165, t1160, t1191, t1075, t1076, t999, pkin(7) * t999, t1020, t969, t1015, t1019, t1014, t1051, t883, t895, t831, t819, t922, t879, t889, t921, t890, t949, t807, t808, t800, t789; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1362, 0, qJDD(2), 0, t1168, 0, t1125, 0, t1213, -t1200, -t1322, -t1213, -t1320, -qJDD(3), t1048, t1049, 0, pkin(2) * t999, -t1141, -t1137, -t1111, t1141, -t1109, -qJDD(3), t900, t907, t925, t832, t1023, -t964, -t1005, -t1021, -t1006, -t1052, t833, t837, t815, t798; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1124, -t1125, 0, 0, t1155, t1144, t1158, t1154, t1156, 0, t1057, t1056, t985, t978, t1018, t967, t1012, t1017, t1011, t1050, t876, t886, t828, t816, t920, t878, t887, t919, t888, t948, t804, t805, t799, t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1194, t1195, t1203, -t1216, t1211, t1216, 0, t1115, t1077, 0, t1103, t1040, t1097, t1101, t1096, t1121, t990, t1004, t877, -qJ(4) * t901, t989, t944, t955, t988, t956, t1010, t836, t838, t824, t809; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1310, t1193, t1209, t1261, t1204, -t1310, -t1115, 0, t1078, 0, t1102, t1038, t1094, t1100, t1093, t1120, t972, t979, t872, t880, t987, t943, t953, t986, t954, t1009, t826, t827, t822, t801; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1213, t1200, t1322, t1213, t1320, qJDD(3), -t1077, -t1078, 0, 0, t1141, t1137, t1111, -t1141, t1109, qJDD(3), t1245, t1299 + 0.2e1 * t1348, t1358, t1359, -t1023, t964, t1005, t1021, t1006, t1052, t1290, t1279, t1264, t1302; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1143, -t1107, t1365, t1175, t1169, -t1175, 0, t1065, t959, 0, t1024, t966, t1007, t1022, t1008, t1054, t905, t906, t845, -pkin(8) * t855; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1343, t1110, t1170, -t1304, t1133, -t1343, -t1065, 0, t960, 0, -t1106, -t1104, -t1033, t1106, t1029, -t1140, t871, t873, -pkin(4) * t963, -pkin(4) * t855; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1141, t1137, t1111, -t1141, t1109, qJDD(3), -t959, -t960, 0, 0, -t1023, t964, t1005, t1021, t1006, t1052, t1317, t1315, t1303, t1354; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1090, t1031, t1366, t1122, t1116, -t1122, 0, t937, t898, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1123, t1032, t1117, t1089, t1063, -t1123, -t937, 0, t899, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1106, t1104, t1033, -t1106, -t1029, t1140, -t898, -t899, 0, 0;];
m_new_reg = t1;