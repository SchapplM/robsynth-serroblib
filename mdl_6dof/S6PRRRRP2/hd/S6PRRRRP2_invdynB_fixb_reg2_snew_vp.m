% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 09:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6PRRRRP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:38:17
% EndTime: 2019-05-05 09:39:29
% DurationCPUTime: 67.16s
% Computational Cost: add. (168022->893), mult. (333057->1375), div. (0->0), fcn. (247111->12), ass. (0->663)
t1021 = sin(pkin(11));
t1023 = cos(pkin(11));
t1022 = sin(pkin(6));
t1024 = cos(pkin(6));
t1029 = sin(qJ(2));
t1033 = cos(qJ(2));
t1028 = sin(qJ(3));
t1032 = cos(qJ(3));
t1027 = sin(qJ(4));
t1031 = cos(qJ(4));
t1030 = cos(qJ(5));
t1026 = sin(qJ(5));
t1141 = qJD(3) + qJD(4);
t982 = (t1027 * t1032 + t1028 * t1031) * qJD(2);
t947 = t1026 * t982 - t1030 * t1141;
t949 = t1026 * t1141 + t1030 * t982;
t895 = t949 * t947;
t1014 = t1032 * qJDD(2);
t1146 = qJD(2) * t1028;
t1126 = qJD(3) * t1146;
t1101 = -t1014 + t1126;
t1145 = qJD(2) * t1032;
t1125 = qJD(3) * t1145;
t1138 = qJDD(2) * t1028;
t989 = t1125 + t1138;
t1121 = t1027 * t989 + t1031 * t1101;
t899 = -t982 * qJD(4) - t1121;
t897 = qJDD(5) - t899;
t1214 = t895 + t897;
t1178 = t1026 * t1214;
t980 = t1027 * t1146 - t1031 * t1145;
t972 = qJD(5) + t980;
t1194 = t972 ^ 2;
t1195 = t947 ^ 2;
t909 = t1195 - t1194;
t773 = t1030 * t909 - t1178;
t1135 = qJDD(3) + qJDD(4);
t900 = -t980 * qJD(4) - t1027 * t1101 + t1031 * t989;
t1120 = t1026 * t900 - t1030 * t1135;
t803 = (qJD(5) - t972) * t949 + t1120;
t699 = t1027 * t773 + t1031 * t803;
t703 = -t1027 * t803 + t1031 * t773;
t623 = t1028 * t699 - t1032 * t703;
t1161 = t1030 * t1214;
t768 = t1026 * t909 + t1161;
t1085 = t1029 * t623 + t1033 * t768;
t620 = t1028 * t703 + t1032 * t699;
t547 = t1022 * t620 + t1024 * t1085;
t583 = -t1029 * t768 + t1033 * t623;
t1347 = t1021 * t547 - t1023 * t583;
t1058 = -t1026 * t1135 - t1030 * t900;
t1039 = -t947 * qJD(5) - t1058;
t1189 = t947 * t972;
t1203 = -t1189 + t1039;
t1180 = t1026 * t1203;
t1148 = qJD(5) + t972;
t802 = t1148 * t949 + t1120;
t716 = t1030 * t802 + t1180;
t946 = t949 ^ 2;
t891 = t946 - t1195;
t680 = t1027 * t716 + t1031 * t891;
t682 = -t1027 * t891 + t1031 * t716;
t600 = t1028 * t680 - t1032 * t682;
t712 = -t1026 * t802 + t1030 * t1203;
t1092 = t1029 * t600 - t1033 * t712;
t597 = t1028 * t682 + t1032 * t680;
t529 = t1022 * t597 + t1024 * t1092;
t573 = t1029 * t712 + t1033 * t600;
t1346 = t1021 * t529 - t1023 * t573;
t1345 = t1021 * t583 + t1023 * t547;
t1344 = t1021 * t573 + t1023 * t529;
t1339 = t1022 * t1085 - t1024 * t620;
t1338 = t1022 * t1092 - t1024 * t597;
t1204 = -t1189 - t1039;
t1242 = -t1026 * t803 + t1030 * t1204;
t1211 = t946 + t1195;
t1241 = -t1026 * t1204 - t1030 * t803;
t1262 = -t1027 * t1211 + t1031 * t1241;
t1264 = t1027 * t1241 + t1031 * t1211;
t1284 = -t1028 * t1264 + t1032 * t1262;
t1305 = t1029 * t1242 + t1033 * t1284;
t1283 = t1028 * t1262 + t1032 * t1264;
t1307 = t1029 * t1284 - t1033 * t1242;
t1317 = -t1022 * t1283 + t1024 * t1307;
t1326 = t1021 * t1305 + t1023 * t1317;
t1337 = qJ(1) * t1326;
t1328 = -t1021 * t1317 + t1023 * t1305;
t1336 = qJ(1) * t1328;
t1331 = pkin(1) * t1317;
t1319 = t1022 * t1307 + t1024 * t1283;
t1330 = pkin(1) * t1319;
t1215 = -t895 + t897;
t1177 = t1026 * t1215;
t910 = -t946 + t1194;
t1244 = -t1030 * t910 - t1177;
t1160 = t1030 * t1215;
t1243 = -t1026 * t910 + t1160;
t1261 = -t1027 * t1204 + t1031 * t1243;
t1263 = t1027 * t1243 + t1031 * t1204;
t1286 = -t1028 * t1263 + t1032 * t1261;
t1306 = -t1029 * t1244 + t1033 * t1286;
t1285 = t1028 * t1261 + t1032 * t1263;
t1308 = t1029 * t1286 + t1033 * t1244;
t1316 = -t1022 * t1285 + t1024 * t1308;
t1329 = -t1021 * t1316 + t1023 * t1306;
t1327 = t1021 * t1306 + t1023 * t1316;
t1325 = (-t1022 * t1319 - t1024 * t1317) * pkin(7);
t1324 = pkin(7) * t1305;
t1318 = t1022 * t1308 + t1024 * t1285;
t1315 = pkin(8) * t1283;
t1310 = -pkin(2) * t1283 - pkin(3) * t1264 - pkin(4) * t1211 - pkin(10) * t1241;
t1309 = -pkin(2) * t1242 + pkin(8) * t1284;
t889 = -t946 - t1194;
t755 = t1030 * t889 - t1178;
t1304 = pkin(2) * t755;
t1303 = pkin(3) * t755;
t1302 = pkin(4) * t755;
t752 = t1026 * t889 + t1161;
t1301 = pkin(10) * t752;
t1300 = pkin(10) * t755;
t1299 = pkin(9) * t1262;
t1298 = pkin(9) * t1264;
t1297 = t1027 * t752;
t1296 = t1029 * t755;
t1294 = t1031 * t752;
t1293 = t1033 * t755;
t1175 = t1026 * t972;
t1133 = t947 * t1175;
t1158 = t1030 * t972;
t904 = t949 * t1158;
t1113 = t904 + t1133;
t1130 = t947 * t1158;
t903 = t949 * t1175;
t1111 = t903 - t1130;
t1198 = t1027 * t897 + t1031 * t1111;
t1201 = t1027 * t1111 - t1031 * t897;
t1221 = -t1028 * t1201 + t1032 * t1198;
t1238 = -t1029 * t1113 + t1033 * t1221;
t1220 = t1028 * t1198 + t1032 * t1201;
t1240 = t1029 * t1221 + t1033 * t1113;
t1265 = -t1022 * t1220 + t1024 * t1240;
t1290 = -t1021 * t1265 + t1023 * t1238;
t843 = -qJD(5) * t949 - t1120;
t1064 = -t1030 * t843 - t1133;
t1065 = -t1026 * t843 + t1130;
t1132 = t1027 * t895;
t1199 = t1031 * t1065 - t1132;
t1129 = t1031 * t895;
t1200 = t1027 * t1065 + t1129;
t1219 = -t1028 * t1200 + t1032 * t1199;
t1237 = -t1029 * t1064 + t1033 * t1219;
t1218 = t1028 * t1199 + t1032 * t1200;
t1239 = t1029 * t1219 + t1033 * t1064;
t1266 = -t1022 * t1218 + t1024 * t1239;
t1289 = -t1021 * t1266 + t1023 * t1237;
t1288 = t1021 * t1238 + t1023 * t1265;
t1287 = t1021 * t1237 + t1023 * t1266;
t1281 = pkin(10) * t1242;
t797 = t1030 * t1039 - t903;
t1110 = t1027 * t797 - t1129;
t1112 = t1031 * t797 + t1132;
t1196 = -t1028 * t1110 + t1032 * t1112;
t796 = -t1026 * t1039 - t904;
t1270 = t1029 * t1196 + t1033 * t796;
t1268 = t1022 * t1239 + t1024 * t1218;
t1267 = t1022 * t1240 + t1024 * t1220;
t1202 = -t1194 - t1195;
t1224 = t1026 * t1202 + t1160;
t1260 = pkin(2) * t1224;
t1259 = pkin(3) * t1224;
t1258 = pkin(4) * t1224;
t1223 = t1030 * t1202 - t1177;
t1257 = pkin(10) * t1223;
t1256 = pkin(10) * t1224;
t1249 = t1027 * t1223;
t1247 = t1029 * t1224;
t1246 = t1031 * t1223;
t1245 = t1033 * t1224;
t1197 = t1028 * t1112 + t1032 * t1110;
t1217 = -t1022 * t1197 + t1024 * t1270;
t1222 = -t1029 * t796 + t1033 * t1196;
t1236 = t1021 * t1222 + t1023 * t1217;
t1235 = -t1021 * t1217 + t1023 * t1222;
t1234 = 2 * qJD(6);
t1232 = qJ(6) * t1203;
t936 = t982 * t980;
t1212 = -t936 + t1135;
t1229 = t1027 * t1212;
t1226 = t1031 * t1212;
t1123 = g(1) * t1021 - t1023 * g(2);
t1187 = g(3) - qJDD(1);
t1225 = -t1022 * t1187 + t1024 * t1123;
t1216 = t1022 * t1270 + t1024 * t1197;
t1134 = t1141 ^ 2;
t970 = t1141 * t980;
t1213 = -t970 + t900;
t1210 = t1021 * t1187;
t1209 = t1023 * t1187;
t1193 = qJD(2) ^ 2;
t996 = g(1) * t1023 + g(2) * t1021;
t918 = t1029 * t1225 - t1033 * t996;
t1036 = -pkin(2) * t1193 + qJDD(2) * pkin(8) + t918;
t965 = t1022 * t1123 + t1024 * t1187;
t864 = t1028 * t1036 + t1032 * t965;
t1005 = t1028 * t1193 * t1032;
t997 = qJDD(3) + t1005;
t1035 = -t864 + (t1125 - t989) * pkin(9) + t997 * pkin(3);
t1000 = qJD(3) * pkin(3) - pkin(9) * t1146;
t1019 = t1032 ^ 2;
t1016 = t1019 * t1193;
t865 = -t1028 * t965 + t1032 * t1036;
t821 = -pkin(3) * t1016 - pkin(9) * t1101 - qJD(3) * t1000 + t865;
t730 = t1027 * t1035 + t1031 * t821;
t931 = pkin(4) * t980 - pkin(10) * t982;
t711 = -pkin(4) * t1134 + pkin(10) * t1135 - t980 * t931 + t730;
t917 = -t1029 * t996 - t1033 * t1225;
t905 = -qJDD(2) * pkin(2) - t1193 * pkin(8) + t917;
t852 = t1101 * pkin(3) - pkin(9) * t1016 + t1000 * t1146 + t905;
t728 = -t1213 * pkin(10) + (t1141 * t982 - t899) * pkin(4) + t852;
t641 = t1026 * t728 + t1030 * t711;
t890 = pkin(5) * t947 - qJ(6) * t949;
t1105 = t897 * qJ(6) + t1234 * t972 - t947 * t890 + t641;
t943 = -t1021 * t1123 - t1023 * t996;
t942 = -t1021 * t996 + t1023 * t1123;
t978 = t980 ^ 2;
t979 = t982 ^ 2;
t845 = t1029 * t917 + t1033 * t918;
t1192 = pkin(7) * t845;
t1191 = pkin(4) * t1027;
t1190 = pkin(5) * t1030;
t1188 = t949 * t972;
t640 = t1026 * t711 - t1030 * t728;
t1186 = t1211 - t1194;
t1185 = qJ(6) * t1030;
t729 = t1027 * t821 - t1031 * t1035;
t710 = -t1135 * pkin(4) - t1134 * pkin(10) + t982 * t931 + t729;
t1182 = t1026 * t710;
t1173 = t1027 * t852;
t925 = t936 + t1135;
t1171 = t1027 * t925;
t645 = t1027 * t730 - t1031 * t729;
t1170 = t1028 * t645;
t1169 = t1028 * t905;
t1168 = t1028 * t997;
t998 = qJDD(3) - t1005;
t1167 = t1028 * t998;
t1165 = t1029 * t965;
t1164 = t1030 * t710;
t1156 = t1031 * t852;
t1155 = t1031 * t925;
t1154 = t1032 * t645;
t1153 = t1032 * t905;
t990 = t1014 - 0.2e1 * t1126;
t950 = t1032 * t990;
t1152 = t1032 * t998;
t1150 = t1033 * t965;
t1147 = qJD(2) * qJD(3);
t1018 = t1028 ^ 2;
t1144 = t1193 * t1018;
t1140 = t1018 + t1019;
t1139 = qJDD(2) * t1022;
t1137 = qJDD(2) * t1029;
t1136 = qJDD(2) * t1033;
t1131 = t1029 * t936;
t1128 = t1033 * t936;
t1127 = -pkin(4) * t1031 - pkin(3);
t1124 = qJ(6) * t1026 + pkin(4);
t646 = t1027 * t729 + t1031 * t730;
t776 = t1028 * t864 + t1032 * t865;
t1119 = t1027 * t1141;
t1118 = t1031 * t1141;
t1116 = t949 * t890 + qJDD(6) + t640;
t1115 = t1029 * t1005;
t1114 = t1033 * t1005;
t1109 = t982 * t1119;
t1108 = t980 * t1119;
t1107 = t982 * t1118;
t1106 = t980 * t1118;
t774 = t1028 * t865 - t1032 * t864;
t991 = t1140 * qJDD(2);
t994 = t1016 + t1144;
t940 = -t1029 * t994 + t1033 * t991;
t1104 = pkin(7) * t940 - t1029 * t774;
t992 = -t1029 * t1193 + t1136;
t1103 = -pkin(7) * t992 - t1165;
t1066 = t1033 * t1193 + t1137;
t1102 = -pkin(7) * t1066 + t1150;
t975 = t1066 * t1024;
t1099 = t1021 * t992 + t1023 * t975;
t929 = t1021 * t975 - t1023 * t992;
t560 = t1026 * t641 - t1030 * t640;
t561 = t1026 * t640 + t1030 * t641;
t596 = -pkin(5) * t1194 + t1105;
t1061 = -t897 * pkin(5) + t1116;
t603 = qJ(6) * t1194 - t1061;
t550 = -t1026 * t603 + t1030 * t596;
t1038 = -t843 * pkin(5) - t1232 + t710;
t631 = (pkin(5) * t972 - (2 * qJD(6))) * t949 + t1038;
t515 = t1027 * t550 - t1031 * t631;
t516 = t1027 * t631 + t1031 * t550;
t476 = -t1028 * t515 + t1032 * t516;
t549 = t1026 * t596 + t1030 * t603;
t1098 = t1029 * t476 - t1033 * t549;
t552 = t1027 * t561 - t1031 * t710;
t553 = t1027 * t710 + t1031 * t561;
t493 = -t1028 * t552 + t1032 * t553;
t1097 = t1029 * t493 - t1033 * t560;
t575 = t1032 * t646 - t1170;
t1096 = t1029 * t575 - t1033 * t852;
t684 = t1031 * t1203 + t1297;
t686 = -t1027 * t1203 + t1294;
t608 = -t1028 * t684 + t1032 * t686;
t1091 = t1029 * t608 + t1293;
t685 = -t1031 * t802 + t1249;
t687 = t1027 * t802 + t1246;
t609 = -t1028 * t685 + t1032 * t687;
t1090 = t1029 * t609 - t1245;
t810 = t1148 * t947 + t1058;
t688 = t1031 * t810 - t1297;
t690 = -t1027 * t810 - t1294;
t615 = -t1028 * t688 + t1032 * t690;
t1089 = t1029 * t615 - t1293;
t804 = -t843 + t1188;
t689 = -t1031 * t804 + t1249;
t691 = t1027 * t804 + t1246;
t616 = -t1028 * t689 + t1032 * t691;
t1088 = t1029 * t616 - t1245;
t874 = (0.2e1 * qJD(4) + qJD(3)) * t982 + t1121;
t784 = -t1027 * t874 + t1031 * t1213;
t786 = -t1027 * t1213 - t1031 * t874;
t706 = -t1028 * t784 + t1032 * t786;
t935 = -t979 + t978;
t1079 = t1029 * t706 + t1033 * t935;
t876 = qJD(3) * t982 - t1121;
t879 = -t970 - t900;
t785 = t1027 * t876 + t1031 * t879;
t787 = -t1027 * t879 + t1031 * t876;
t707 = -t1028 * t785 + t1032 * t787;
t902 = -t978 - t979;
t1078 = t1029 * t707 - t1033 * t902;
t919 = -t1134 - t978;
t850 = t1027 * t919 + t1226;
t851 = t1031 * t919 - t1229;
t760 = -t1028 * t850 + t1032 * t851;
t1077 = t1029 * t760 - t1033 * t874;
t1076 = t1029 * t776 - t1033 * t905;
t960 = -t979 - t1134;
t880 = t1031 * t960 - t1171;
t881 = -t1027 * t960 - t1155;
t789 = -t1028 * t880 + t1032 * t881;
t1075 = t1029 * t789 - t1033 * t1213;
t967 = -t979 + t1134;
t884 = t1031 * t967 + t1229;
t886 = -t1027 * t967 + t1226;
t800 = -t1028 * t884 + t1032 * t886;
t1074 = t1029 * t800 + t1033 * t879;
t966 = t978 - t1134;
t885 = t1027 * t966 + t1155;
t887 = t1031 * t966 - t1171;
t801 = -t1028 * t885 + t1032 * t887;
t1073 = t1029 * t801 - t1033 * t876;
t1072 = t1029 * t918 - t1033 * t917;
t988 = 0.2e1 * t1125 + t1138;
t938 = -t1028 * t988 + t950;
t995 = t1016 - t1144;
t1071 = t1029 * t938 + t1033 * t995;
t1034 = qJD(3) ^ 2;
t1004 = -t1016 - t1034;
t957 = t1004 * t1032 - t1168;
t1070 = t1029 * t957 + t1033 * t990;
t1002 = -t1034 - t1144;
t959 = -t1002 * t1028 - t1152;
t1069 = t1029 * t959 - t1033 * t988;
t1068 = t1029 * t991 + t1033 * t994;
t986 = t1140 * t1147;
t1067 = -qJDD(3) * t1033 + t1029 * t986;
t860 = t1031 * t899 + t1108;
t861 = -t1027 * t899 + t1106;
t764 = -t1028 * t860 + t1032 * t861;
t1063 = t1029 * t764 + t1128;
t862 = t1027 * t900 + t1107;
t863 = t1031 * t900 - t1109;
t765 = -t1028 * t862 + t1032 * t863;
t1062 = t1029 * t765 - t1128;
t1001 = t1034 - t1144;
t987 = t1032 * t997;
t958 = -t1001 * t1028 + t987;
t1060 = -t1028 * t1136 + t1029 * t958;
t1003 = t1016 - t1034;
t956 = t1003 * t1032 - t1167;
t1059 = -t1014 * t1033 + t1029 * t956;
t906 = -t1108 - t1107;
t907 = -t1106 + t1109;
t832 = -t1028 * t906 + t1032 * t907;
t1057 = t1029 * t832 - t1033 * t1135;
t491 = -pkin(4) * t549 - pkin(5) * t603 - qJ(6) * t596;
t496 = -pkin(10) * t549 + (pkin(5) * t1026 - t1185) * t631;
t453 = -pkin(3) * t549 + pkin(9) * t516 + t1027 * t496 + t1031 * t491;
t457 = -pkin(9) * t515 - t1027 * t491 + t1031 * t496;
t475 = t1028 * t516 + t1032 * t515;
t433 = -pkin(8) * t475 - t1028 * t453 + t1032 * t457;
t452 = -pkin(2) * t475 - pkin(3) * t515 - pkin(10) * t550 + (t1124 + t1190) * t631;
t468 = t1029 * t549 + t1033 * t476;
t1056 = pkin(7) * t468 + t1029 * t433 + t1033 * t452;
t472 = pkin(9) * t553 + (-pkin(10) * t1027 + t1127) * t560;
t480 = -pkin(9) * t552 + (-pkin(10) * t1031 + t1191) * t560;
t492 = t1028 * t553 + t1032 * t552;
t447 = -pkin(8) * t492 - t1028 * t472 + t1032 * t480;
t467 = -pkin(2) * t492 - pkin(3) * t552 + pkin(4) * t710 - pkin(10) * t561;
t479 = t1029 * t560 + t1033 * t493;
t1055 = pkin(7) * t479 + t1029 * t447 + t1033 * t467;
t580 = pkin(5) * t1186 + t1105;
t585 = qJ(6) * t1186 + t1061;
t521 = -t1026 * t580 + t1030 * t585 - t1281;
t637 = -pkin(4) * t1242 - pkin(5) * t1204 + qJ(6) * t803;
t482 = -pkin(3) * t1242 + t1027 * t521 + t1031 * t637 + t1299;
t490 = -t1027 * t637 + t1031 * t521 - t1298;
t460 = -t1028 * t482 + t1032 * t490 - t1315;
t483 = -t1026 * t585 - t1030 * t580 + t1310;
t1054 = t1029 * t460 + t1033 * t483 + t1324;
t557 = t1302 - qJ(6) * t1214 + (t889 + t1194) * pkin(5) - t1105;
t1037 = t1234 * t949 - t1038;
t592 = -pkin(5) * t1188 + t1037 + t1232;
t558 = -pkin(5) * t1180 + t1030 * t592 + t1300;
t494 = pkin(9) * t686 + t1027 * t558 + t1031 * t557 + t1303;
t498 = -pkin(9) * t684 - t1027 * t557 + t1031 * t558;
t606 = t1028 * t686 + t1032 * t684;
t465 = -pkin(8) * t606 - t1028 * t494 + t1032 * t498;
t511 = -pkin(2) * t606 - pkin(3) * t684 - t1301 - t1026 * t592 + (-pkin(4) - t1190) * t1203;
t576 = t1033 * t608 - t1296;
t1053 = pkin(7) * t576 + t1029 * t465 + t1033 * t511;
t593 = (-t804 - t1188) * pkin(5) + t1037;
t559 = -t1026 * t593 - t1185 * t804 - t1256;
t566 = -t1258 + (-t1202 - t1194) * qJ(6) + (-t1215 - t897) * pkin(5) + t1116;
t495 = pkin(9) * t691 + t1027 * t559 + t1031 * t566 - t1259;
t500 = -pkin(9) * t689 - t1027 * t566 + t1031 * t559;
t614 = t1028 * t691 + t1032 * t689;
t466 = -pkin(8) * t614 - t1028 * t495 + t1032 * t500;
t512 = -pkin(2) * t614 - pkin(3) * t689 - t1030 * t593 + t1124 * t804 - t1257;
t579 = t1033 * t616 + t1247;
t1052 = pkin(7) * t579 + t1029 * t466 + t1033 * t512;
t554 = -t560 - t1281;
t509 = t1027 * t554 + t1127 * t1242 + t1299;
t514 = t1031 * t554 + t1191 * t1242 - t1298;
t470 = -t1028 * t509 + t1032 * t514 - t1315;
t497 = t1310 - t561;
t1051 = t1029 * t470 + t1033 * t497 + t1324;
t594 = t640 - t1258;
t642 = t1182 - t1256;
t527 = pkin(9) * t687 + t1027 * t642 + t1031 * t594 - t1259;
t539 = -pkin(9) * t685 - t1027 * t594 + t1031 * t642;
t607 = t1028 * t687 + t1032 * t685;
t477 = -pkin(8) * t607 - t1028 * t527 + t1032 * t539;
t532 = -pkin(2) * t607 - pkin(3) * t685 + pkin(4) * t802 + t1164 - t1257;
t577 = t1033 * t609 + t1247;
t1050 = pkin(7) * t577 + t1029 * t477 + t1033 * t532;
t595 = t641 - t1302;
t644 = t1164 - t1300;
t530 = pkin(9) * t690 + t1027 * t644 + t1031 * t595 - t1303;
t540 = -pkin(9) * t688 - t1027 * t595 + t1031 * t644;
t613 = t1028 * t690 + t1032 * t688;
t478 = -pkin(8) * t613 - t1028 * t530 + t1032 * t540;
t533 = -pkin(2) * t613 - pkin(3) * t688 - pkin(4) * t810 - t1182 + t1301;
t578 = t1033 * t615 + t1296;
t1049 = pkin(7) * t578 + t1029 * t478 + t1033 * t533;
t574 = t1028 * t646 + t1154;
t638 = -pkin(3) * t852 + pkin(9) * t646;
t520 = -pkin(8) * t574 - pkin(9) * t1154 - t1028 * t638;
t551 = -pkin(2) * t574 - pkin(3) * t645;
t567 = t1029 * t852 + t1033 * t575;
t1048 = pkin(7) * t567 + t1029 * t520 + t1033 * t551;
t612 = -pkin(3) * t902 + pkin(9) * t787 + t646;
t630 = -pkin(9) * t785 - t645;
t705 = t1028 * t787 + t1032 * t785;
t534 = -pkin(8) * t705 - t1028 * t612 + t1032 * t630;
t647 = -pkin(2) * t705 - pkin(3) * t785;
t676 = t1029 * t902 + t1033 * t707;
t1047 = pkin(7) * t676 + t1029 * t534 + t1033 * t647;
t723 = -pkin(3) * t874 + pkin(9) * t851 - t1156;
t759 = t1028 * t851 + t1032 * t850;
t761 = -pkin(9) * t850 + t1173;
t632 = -pkin(8) * t759 - t1028 * t723 + t1032 * t761;
t643 = -pkin(2) * t759 - pkin(3) * t850 + t729;
t721 = t1029 * t874 + t1033 * t760;
t1046 = pkin(7) * t721 + t1029 * t632 + t1033 * t643;
t726 = -pkin(3) * t1213 + pkin(9) * t881 + t1173;
t777 = -pkin(9) * t880 + t1156;
t788 = t1028 * t881 + t1032 * t880;
t639 = -pkin(8) * t788 - t1028 * t726 + t1032 * t777;
t648 = -pkin(2) * t788 - pkin(3) * t880 + t730;
t731 = t1029 * t1213 + t1033 * t789;
t1045 = pkin(7) * t731 + t1029 * t639 + t1033 * t648;
t953 = t1004 * t1028 + t987;
t829 = -pkin(2) * t953 + t864;
t858 = -pkin(8) * t953 + t1169;
t915 = -t1029 * t990 + t1033 * t957;
t1044 = pkin(7) * t915 + t1029 * t858 + t1033 * t829;
t955 = t1002 * t1032 - t1167;
t830 = -pkin(2) * t955 + t865;
t859 = -pkin(8) * t955 + t1153;
t916 = t1029 * t988 + t1033 * t959;
t1043 = pkin(7) * t916 + t1029 * t859 + t1033 * t830;
t962 = -t1019 * t1147 + t1028 * t1101;
t1042 = t1029 * t962 - t1114;
t963 = -t1018 * t1147 + t1032 * t989;
t1041 = t1029 * t963 + t1114;
t738 = t1029 * t905 + t1033 * t776;
t1040 = pkin(7) * t738 + (-pkin(2) * t1033 - pkin(8) * t1029) * t774;
t976 = t992 * t1024;
t974 = t992 * t1022;
t973 = t1066 * t1022;
t964 = qJDD(3) * t1029 + t1033 * t986;
t954 = t1001 * t1032 + t1168;
t952 = t1003 * t1028 + t1152;
t951 = (t989 + t1125) * t1028;
t941 = t1067 * t1024;
t937 = t1028 * t990 + t1032 * t988;
t934 = t1068 * t1024;
t933 = t1068 * t1022;
t930 = -t1021 * t976 - t1023 * t1066;
t928 = -t1021 * t1066 + t1023 * t976;
t923 = t1033 * t963 - t1115;
t922 = t1033 * t962 + t1115;
t921 = t1028 * t1137 + t1033 * t958;
t920 = t1014 * t1029 + t1033 * t956;
t901 = -t1029 * t995 + t1033 * t938;
t883 = -t1150 + (t1022 * t973 + t1024 * t975) * pkin(7);
t882 = -t1165 + (-t1022 * t974 - t1024 * t976) * pkin(7);
t873 = -t1021 * t934 + t1023 * t940;
t872 = t1021 * t940 + t1023 * t934;
t871 = -t1022 * t951 + t1024 * t1041;
t870 = -t1022 * t950 + t1024 * t1042;
t869 = -t1022 * t954 + t1024 * t1060;
t868 = -t1022 * t952 + t1024 * t1059;
t856 = -t1022 * t955 + t1024 * t1069;
t855 = -t1022 * t953 + t1024 * t1070;
t854 = t1022 * t1069 + t1024 * t955;
t853 = t1022 * t1070 + t1024 * t953;
t842 = -t1022 * t937 + t1024 * t1071;
t841 = pkin(2) * t990 + pkin(8) * t957 - t1153;
t840 = -pkin(2) * t988 + pkin(8) * t959 + t1169;
t837 = t845 * t1024;
t831 = t1028 * t907 + t1032 * t906;
t823 = -pkin(1) * t974 + t1022 * t917 + t1024 * t1102;
t822 = pkin(1) * t973 + t1022 * t918 + t1024 * t1103;
t817 = t1029 * t1135 + t1033 * t832;
t816 = t1022 * t965 + t1024 * t1072;
t815 = t1022 * t1072 - t1024 * t965;
t814 = -t1021 * t856 + t1023 * t916;
t813 = -t1021 * t855 + t1023 * t915;
t812 = t1021 * t916 + t1023 * t856;
t811 = t1021 * t915 + t1023 * t855;
t799 = t1028 * t887 + t1032 * t885;
t798 = t1028 * t886 + t1032 * t884;
t763 = t1028 * t863 + t1032 * t862;
t762 = t1028 * t861 + t1032 * t860;
t754 = pkin(2) * t994 + pkin(8) * t991 + t776;
t749 = t1033 * t765 + t1131;
t748 = t1033 * t764 - t1131;
t739 = -pkin(2) * t905 + pkin(8) * t776;
t737 = -pkin(1) * t815 + t1024 * t1192;
t736 = t1029 * t876 + t1033 * t801;
t735 = -t1029 * t879 + t1033 * t800;
t734 = -t1021 * t816 + t1023 * t845;
t733 = t1021 * t845 + t1023 * t816;
t732 = -t1022 * t831 + t1024 * t1057;
t722 = -t1033 * t774 + (-t1022 * t933 - t1024 * t934) * pkin(7);
t720 = (-t1022 * t815 - t1024 * t816) * pkin(7);
t704 = t1028 * t786 + t1032 * t784;
t679 = -t1029 * t935 + t1033 * t706;
t678 = -t1029 * t830 + t1033 * t859 + (-t1022 * t854 - t1024 * t856) * pkin(7);
t677 = -t1029 * t829 + t1033 * t858 + (-t1022 * t853 - t1024 * t855) * pkin(7);
t675 = -t1022 * t799 + t1024 * t1073;
t674 = -t1022 * t798 + t1024 * t1074;
t669 = -t1022 * t763 + t1024 * t1062;
t668 = -t1022 * t762 + t1024 * t1063;
t667 = -t1022 * t788 + t1024 * t1075;
t666 = t1022 * t1075 + t1024 * t788;
t665 = -t1022 * t774 + t1024 * t1076;
t664 = t1022 * t1076 + t1024 * t774;
t655 = -pkin(1) * t854 - t1022 * t840 + t1024 * t1043;
t654 = -pkin(1) * t853 - t1022 * t841 + t1024 * t1044;
t651 = -pkin(1) * t933 - t1022 * t754 + t1024 * t1104;
t650 = -t1022 * t759 + t1024 * t1077;
t649 = t1022 * t1077 + t1024 * t759;
t629 = -t1021 * t665 + t1023 * t738;
t628 = t1021 * t738 + t1023 * t665;
t627 = -pkin(2) * t1213 + pkin(8) * t789 + t1028 * t777 + t1032 * t726;
t626 = -t1021 * t667 + t1023 * t731;
t625 = t1021 * t731 + t1023 * t667;
t611 = -t1022 * t704 + t1024 * t1079;
t610 = -pkin(2) * t874 + pkin(8) * t760 + t1028 * t761 + t1032 * t723;
t605 = -t1022 * t705 + t1024 * t1078;
t604 = t1022 * t1078 + t1024 * t705;
t602 = -t1021 * t650 + t1023 * t721;
t601 = t1021 * t721 + t1023 * t650;
t563 = -t1021 * t605 + t1023 * t676;
t562 = t1021 * t676 + t1023 * t605;
t556 = (pkin(2) * t1029 - pkin(8) * t1033) * t774 + (-t1022 * t664 - t1024 * t665) * pkin(7);
t555 = -pkin(1) * t664 - t1022 * t739 + t1024 * t1040;
t544 = -t1022 * t614 + t1024 * t1088;
t543 = -t1022 * t613 + t1024 * t1089;
t542 = t1022 * t1088 + t1024 * t614;
t541 = t1022 * t1089 + t1024 * t613;
t538 = -t1022 * t607 + t1024 * t1090;
t537 = -t1022 * t606 + t1024 * t1091;
t536 = t1022 * t1090 + t1024 * t607;
t535 = t1022 * t1091 + t1024 * t606;
t531 = -pkin(2) * t902 + pkin(8) * t707 + t1028 * t630 + t1032 * t612;
t526 = -t1029 * t648 + t1033 * t639 + (-t1022 * t666 - t1024 * t667) * pkin(7);
t519 = -t1022 * t574 + t1024 * t1096;
t518 = t1022 * t1096 + t1024 * t574;
t517 = -t1029 * t643 + t1033 * t632 + (-t1022 * t649 - t1024 * t650) * pkin(7);
t513 = -pkin(2) * t852 + pkin(8) * t575 - pkin(9) * t1170 + t1032 * t638;
t510 = -pkin(1) * t666 - t1022 * t627 + t1024 * t1045;
t508 = -t1021 * t544 + t1023 * t579;
t507 = -t1021 * t543 + t1023 * t578;
t506 = t1021 * t579 + t1023 * t544;
t505 = t1021 * t578 + t1023 * t543;
t504 = -t1021 * t538 + t1023 * t577;
t503 = -t1021 * t537 + t1023 * t576;
t502 = t1021 * t577 + t1023 * t538;
t501 = t1021 * t576 + t1023 * t537;
t499 = -pkin(1) * t649 - t1022 * t610 + t1024 * t1046;
t485 = -t1021 * t519 + t1023 * t567;
t484 = t1021 * t567 + t1023 * t519;
t481 = -t1029 * t647 + t1033 * t534 + (-t1022 * t604 - t1024 * t605) * pkin(7);
t474 = pkin(8) * t615 + t1028 * t540 + t1032 * t530 - t1304;
t473 = pkin(8) * t609 + t1028 * t539 + t1032 * t527 - t1260;
t471 = -pkin(1) * t604 - t1022 * t531 + t1024 * t1047;
t469 = t1028 * t514 + t1032 * t509 + t1309;
t464 = pkin(8) * t616 + t1028 * t500 + t1032 * t495 - t1260;
t463 = pkin(8) * t608 + t1028 * t498 + t1032 * t494 + t1304;
t462 = -t1022 * t492 + t1024 * t1097;
t461 = t1022 * t1097 + t1024 * t492;
t459 = t1028 * t490 + t1032 * t482 + t1309;
t458 = -t1029 * t551 + t1033 * t520 + (-t1022 * t518 - t1024 * t519) * pkin(7);
t456 = -pkin(1) * t518 - t1022 * t513 + t1024 * t1048;
t455 = -t1029 * t533 + t1033 * t478 + (-t1022 * t541 - t1024 * t543) * pkin(7);
t454 = -t1029 * t532 + t1033 * t477 + (-t1022 * t536 - t1024 * t538) * pkin(7);
t451 = -t1022 * t475 + t1024 * t1098;
t450 = t1022 * t1098 + t1024 * t475;
t449 = -t1021 * t462 + t1023 * t479;
t448 = t1021 * t479 + t1023 * t462;
t446 = -t1029 * t512 + t1033 * t466 + (-t1022 * t542 - t1024 * t544) * pkin(7);
t445 = -t1029 * t511 + t1033 * t465 + (-t1022 * t535 - t1024 * t537) * pkin(7);
t444 = -t1029 * t497 + t1033 * t470 + t1325;
t443 = -pkin(2) * t560 + pkin(8) * t493 + t1028 * t480 + t1032 * t472;
t442 = -pkin(1) * t541 - t1022 * t474 + t1024 * t1049;
t441 = -pkin(1) * t536 - t1022 * t473 + t1024 * t1050;
t440 = -t1029 * t483 + t1033 * t460 + t1325;
t439 = -t1021 * t451 + t1023 * t468;
t438 = t1021 * t468 + t1023 * t451;
t437 = -t1022 * t469 + t1024 * t1051 - t1330;
t436 = -pkin(1) * t542 - t1022 * t464 + t1024 * t1052;
t435 = -pkin(1) * t535 - t1022 * t463 + t1024 * t1053;
t434 = -t1022 * t459 + t1024 * t1054 - t1330;
t432 = -pkin(2) * t549 + pkin(8) * t476 + t1028 * t457 + t1032 * t453;
t431 = -t1029 * t467 + t1033 * t447 + (-t1022 * t461 - t1024 * t462) * pkin(7);
t430 = -pkin(1) * t461 - t1022 * t443 + t1024 * t1055;
t429 = -t1029 * t452 + t1033 * t433 + (-t1022 * t450 - t1024 * t451) * pkin(7);
t428 = -pkin(1) * t450 - t1022 * t432 + t1024 * t1056;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t943, 0, 0, 0, 0, 0, 0, t930, t929, 0, t734, 0, 0, 0, 0, 0, 0, t813, t814, t873, t629, 0, 0, 0, 0, 0, 0, t602, t626, t563, t485, 0, 0, 0, 0, 0, 0, t504, t507, t1328, t449, 0, 0, 0, 0, 0, 0, t508, t1328, t503, t439; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t942, 0, 0, 0, 0, 0, 0, t928, -t1099, 0, t733, 0, 0, 0, 0, 0, 0, t811, t812, t872, t628, 0, 0, 0, 0, 0, 0, t601, t625, t562, t484, 0, 0, 0, 0, 0, 0, t502, t505, t1326, t448, 0, 0, 0, 0, 0, 0, t506, t1326, t501, t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1187, 0, 0, 0, 0, 0, 0, t974, -t973, 0, t815, 0, 0, 0, 0, 0, 0, t853, t854, t933, t664, 0, 0, 0, 0, 0, 0, t649, t666, t604, t518, 0, 0, 0, 0, 0, 0, t536, t541, t1319, t461, 0, 0, 0, 0, 0, 0, t542, t1319, t535, t450; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t1210, -t1209, -t942, -qJ(1) * t942, 0, 0, -t929, 0, t930, t1021 * t1139, -qJ(1) * t928 - t1021 * t823 + t1023 * t882, qJ(1) * t1099 - t1021 * t822 + t1023 * t883, -t1021 * t837 - t1023 * t1072, -qJ(1) * t733 - t1021 * t737 + t1023 * t720, -t1021 * t871 + t1023 * t923, -t1021 * t842 + t1023 * t901, -t1021 * t869 + t1023 * t921, -t1021 * t870 + t1023 * t922, -t1021 * t868 + t1023 * t920, -t1021 * t941 + t1023 * t964, -qJ(1) * t811 - t1021 * t654 + t1023 * t677, -qJ(1) * t812 - t1021 * t655 + t1023 * t678, -qJ(1) * t872 - t1021 * t651 + t1023 * t722, -qJ(1) * t628 - t1021 * t555 + t1023 * t556, -t1021 * t669 + t1023 * t749, -t1021 * t611 + t1023 * t679, -t1021 * t674 + t1023 * t735, -t1021 * t668 + t1023 * t748, -t1021 * t675 + t1023 * t736, -t1021 * t732 + t1023 * t817, -qJ(1) * t601 - t1021 * t499 + t1023 * t517, -qJ(1) * t625 - t1021 * t510 + t1023 * t526, -qJ(1) * t562 - t1021 * t471 + t1023 * t481, -qJ(1) * t484 - t1021 * t456 + t1023 * t458, t1235, -t1346, t1329, t1289, t1347, t1290, -qJ(1) * t502 - t1021 * t441 + t1023 * t454, -qJ(1) * t505 - t1021 * t442 + t1023 * t455, -t1021 * t437 + t1023 * t444 - t1337, -qJ(1) * t448 - t1021 * t430 + t1023 * t431, t1235, t1329, t1346, t1290, -t1347, t1289, -qJ(1) * t506 - t1021 * t436 + t1023 * t446, -t1021 * t434 + t1023 * t440 - t1337, -qJ(1) * t501 - t1021 * t435 + t1023 * t445, -qJ(1) * t438 - t1021 * t428 + t1023 * t429; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t1209, -t1210, t943, qJ(1) * t943, 0, 0, t1099, 0, t928, -t1023 * t1139, qJ(1) * t930 + t1021 * t882 + t1023 * t823, qJ(1) * t929 + t1021 * t883 + t1023 * t822, -t1021 * t1072 + t1023 * t837, qJ(1) * t734 + t1021 * t720 + t1023 * t737, t1021 * t923 + t1023 * t871, t1021 * t901 + t1023 * t842, t1021 * t921 + t1023 * t869, t1021 * t922 + t1023 * t870, t1021 * t920 + t1023 * t868, t1021 * t964 + t1023 * t941, qJ(1) * t813 + t1021 * t677 + t1023 * t654, qJ(1) * t814 + t1021 * t678 + t1023 * t655, qJ(1) * t873 + t1021 * t722 + t1023 * t651, qJ(1) * t629 + t1021 * t556 + t1023 * t555, t1021 * t749 + t1023 * t669, t1021 * t679 + t1023 * t611, t1021 * t735 + t1023 * t674, t1021 * t748 + t1023 * t668, t1021 * t736 + t1023 * t675, t1021 * t817 + t1023 * t732, qJ(1) * t602 + t1021 * t517 + t1023 * t499, qJ(1) * t626 + t1021 * t526 + t1023 * t510, qJ(1) * t563 + t1021 * t481 + t1023 * t471, qJ(1) * t485 + t1021 * t458 + t1023 * t456, t1236, t1344, t1327, t1287, -t1345, t1288, qJ(1) * t504 + t1021 * t454 + t1023 * t441, qJ(1) * t507 + t1021 * t455 + t1023 * t442, t1021 * t444 + t1023 * t437 + t1336, qJ(1) * t449 + t1021 * t431 + t1023 * t430, t1236, t1327, -t1344, t1288, t1345, t1287, qJ(1) * t508 + t1021 * t446 + t1023 * t436, t1021 * t440 + t1023 * t434 + t1336, qJ(1) * t503 + t1021 * t445 + t1023 * t435, qJ(1) * t439 + t1021 * t429 + t1023 * t428; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t1123, t996, 0, 0, 0, 0, t973, 0, t974, t1024 * qJDD(2), pkin(1) * t976 + t1022 * t1102 - t1024 * t917, -pkin(1) * t975 + t1022 * t1103 - t1024 * t918, t845 * t1022, pkin(1) * t816 + t1022 * t1192, t1022 * t1041 + t1024 * t951, t1022 * t1071 + t1024 * t937, t1022 * t1060 + t1024 * t954, t1022 * t1042 + t1024 * t950, t1022 * t1059 + t1024 * t952, t1067 * t1022, pkin(1) * t855 + t1022 * t1044 + t1024 * t841, pkin(1) * t856 + t1022 * t1043 + t1024 * t840, pkin(1) * t934 + t1022 * t1104 + t1024 * t754, pkin(1) * t665 + t1022 * t1040 + t1024 * t739, t1022 * t1062 + t1024 * t763, t1022 * t1079 + t1024 * t704, t1022 * t1074 + t1024 * t798, t1022 * t1063 + t1024 * t762, t1022 * t1073 + t1024 * t799, t1022 * t1057 + t1024 * t831, pkin(1) * t650 + t1022 * t1046 + t1024 * t610, pkin(1) * t667 + t1022 * t1045 + t1024 * t627, pkin(1) * t605 + t1022 * t1047 + t1024 * t531, pkin(1) * t519 + t1022 * t1048 + t1024 * t513, t1216, t1338, t1318, t1268, -t1339, t1267, pkin(1) * t538 + t1022 * t1050 + t1024 * t473, pkin(1) * t543 + t1022 * t1049 + t1024 * t474, t1022 * t1051 + t1024 * t469 + t1331, pkin(1) * t462 + t1022 * t1055 + t1024 * t443, t1216, t1318, -t1338, t1267, t1339, t1268, pkin(1) * t544 + t1022 * t1052 + t1024 * t464, t1022 * t1054 + t1024 * t459 + t1331, pkin(1) * t537 + t1022 * t1053 + t1024 * t463, pkin(1) * t451 + t1022 * t1056 + t1024 * t432;];
tauB_reg  = t1;
