% Calculate vector of inverse dynamics joint torques for
% S6RRRRRR8
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
%   see S6RRRRRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:07:56
% EndTime: 2019-03-10 05:08:53
% DurationCPUTime: 42.43s
% Computational Cost: add. (39395->1003), mult. (109768->1372), div. (0->0), fcn. (94637->18), ass. (0->409)
t1073 = sin(pkin(7));
t1084 = cos(qJ(3));
t1235 = t1073 * t1084;
t1206 = qJD(3) * t1235;
t1074 = sin(pkin(6));
t1075 = cos(pkin(7));
t1085 = cos(qJ(2));
t1226 = t1084 * t1085;
t1079 = sin(qJ(3));
t1080 = sin(qJ(2));
t1229 = t1079 * t1080;
t1126 = -t1075 * t1229 + t1226;
t999 = t1126 * t1074;
t989 = qJD(1) * t999;
t1347 = t1206 - t989;
t1290 = pkin(10) * t1073;
t1132 = pkin(2) * t1080 - t1085 * t1290;
t1114 = t1132 * t1074;
t1003 = qJD(1) * t1114;
t1231 = t1075 * t1079;
t1236 = t1073 * t1079;
t1055 = pkin(10) * t1236;
t1230 = t1075 * t1084;
t1308 = pkin(2) * t1230 - t1055;
t1288 = cos(pkin(6));
t1215 = pkin(1) * t1288;
t1060 = t1085 * t1215;
t1051 = qJD(1) * t1060;
t1289 = pkin(10) * t1075;
t1152 = t1074 * (-pkin(9) - t1289);
t1133 = t1080 * t1152;
t968 = qJD(1) * t1133 + t1051;
t1059 = t1080 * t1215;
t1096 = t1085 * t1152 - t1059;
t969 = t1096 * qJD(1);
t1346 = t1308 * qJD(3) - t1003 * t1236 - t1084 * t968 - t969 * t1231;
t891 = t1075 * t1003 - t1073 * t969;
t1227 = t1080 * t1084;
t1228 = t1079 * t1085;
t1128 = t1075 * t1227 + t1228;
t998 = t1128 * t1074;
t988 = qJD(1) * t998;
t1345 = -pkin(3) * t988 + pkin(11) * t989 + (pkin(3) * t1079 - pkin(11) * t1084) * t1073 * qJD(3) - t891;
t1234 = t1074 * t1080;
t1202 = t1073 * t1234;
t1169 = qJD(1) * t1202;
t1344 = -pkin(11) * t1169 + t1346;
t1078 = sin(qJ(4));
t1083 = cos(qJ(4));
t1014 = -t1083 * t1075 + t1078 * t1236;
t1166 = t1078 * t1202;
t1268 = -qJD(1) * t1166 - qJD(4) * t1014 + t1083 * t1347;
t1015 = t1075 * t1078 + t1083 * t1236;
t1165 = t1083 * t1202;
t1267 = qJD(1) * t1165 + qJD(4) * t1015 + t1078 * t1347;
t1207 = qJD(3) * t1236;
t1324 = t1207 - t988;
t1081 = cos(qJ(6));
t1243 = qJD(6) * t1081;
t1077 = sin(qJ(5));
t1082 = cos(qJ(5));
t1127 = t1075 * t1228 + t1227;
t1106 = t1127 * t1074;
t1193 = t1288 * qJD(1);
t1149 = t1193 + qJD(2);
t1117 = t1073 * t1149;
t1111 = t1079 * t1117;
t938 = qJD(1) * t1106 + t1111;
t1233 = t1074 * t1085;
t1201 = t1073 * t1233;
t991 = qJD(1) * t1201 - t1075 * t1149 - qJD(3);
t1135 = t1078 * t991 - t1083 * t938;
t878 = -t1078 * t938 - t1083 * t991;
t801 = -t1077 * t1135 - t1082 * t878;
t1325 = t1081 * t801;
t1343 = t1243 + t1325;
t1076 = sin(qJ(6));
t1072 = qJ(4) + qJ(5);
t1067 = sin(t1072);
t1068 = cos(t1072);
t1086 = cos(qJ(1));
t1194 = t1086 * t1288;
t1295 = sin(qJ(1));
t1016 = t1080 * t1295 - t1085 * t1194;
t1017 = t1080 * t1194 + t1085 * t1295;
t1232 = t1074 * t1086;
t1200 = t1073 * t1232;
t902 = t1016 * t1231 - t1017 * t1084 + t1079 * t1200;
t972 = -t1016 * t1073 + t1075 * t1232;
t850 = t1067 * t972 + t1068 * t902;
t1342 = t1076 * t850;
t1341 = t1081 * t850;
t1340 = t1345 * t1083;
t1198 = t1075 * t1233;
t1339 = qJD(1) * t1198 + t1117;
t1008 = (-pkin(3) * t1084 - pkin(11) * t1079 - pkin(2)) * t1073;
t994 = t1083 * t1008;
t1338 = -qJD(4) * t994 - t1078 * t1345 - t1344 * t1083;
t1136 = t1077 * t878 - t1082 * t1135;
t1253 = -qJD(6) - t801;
t1333 = t1076 * t1253;
t1245 = qJD(5) * t1082;
t1246 = qJD(5) * t1077;
t1247 = qJD(4) * t1083;
t1248 = qJD(4) * t1078;
t1188 = t1288 * qJDD(1);
t1054 = t1188 + qJDD(2);
t1222 = qJDD(1) * t1085;
t1190 = t1074 * t1222;
t1156 = t1075 * t1190;
t1197 = t1074 * t1229;
t1168 = qJD(1) * t1197;
t1321 = t1084 * t1339 - t1168;
t1223 = qJDD(1) * t1080;
t1191 = t1074 * t1223;
t1208 = qJD(2) * t1233;
t1323 = qJD(1) * t1208 + t1191;
t844 = -qJD(2) * t1075 * t1168 + qJD(3) * t1321 + t1054 * t1236 + t1079 * t1156 + t1084 * t1323;
t1209 = qJD(2) * t1234;
t1164 = qJD(1) * t1209;
t943 = t1054 * t1075 + qJDD(3) + (t1164 - t1190) * t1073;
t768 = t1078 * t943 + t1083 * t844 - t991 * t1247 - t1248 * t938;
t769 = -qJD(4) * t1135 + t1078 * t844 - t1083 * t943;
t719 = -t1077 * t769 + t1082 * t768 + t1135 * t1246 + t878 * t1245;
t1090 = qJD(2) * t1128 + qJD(3) * t1127;
t845 = qJD(3) * t1111 - t1054 * t1235 + t1074 * (qJD(1) * t1090 + t1079 * t1223) - t1084 * t1156;
t839 = qJDD(4) + t845;
t838 = qJDD(5) + t839;
t932 = qJD(4) - t1321;
t929 = qJD(5) + t932;
t1218 = t1076 * t838 + t1081 * t719 + t929 * t1243;
t1244 = qJD(6) * t1076;
t696 = -t1136 * t1244 + t1218;
t694 = t696 * t1076;
t1186 = t1076 * t719 - t1081 * t838;
t778 = t1076 * t929 + t1081 * t1136;
t697 = qJD(6) * t778 + t1186;
t720 = qJD(5) * t1136 + t1077 * t768 + t1082 * t769;
t718 = qJDD(6) + t720;
t715 = t1081 * t718;
t1261 = t1076 * t1136;
t776 = -t1081 * t929 + t1261;
t1337 = t838 * MDP(29) - t720 * MDP(28) - t801 ^ 2 * MDP(26) + (t801 * t929 + t719) * MDP(27) + (MDP(25) * t801 + MDP(26) * t1136 + MDP(28) * t929 + MDP(36) * t1253) * t1136 + (t1076 * t718 - t1136 * t778 - t1253 * t1343) * MDP(34) + (t1343 * t778 + t694) * MDP(32) + (t1136 * t776 - t1253 * t1333 + t715) * MDP(35) + (-t1076 * t697 + t696 * t1081 + t1333 * t778 - t1343 * t776) * MDP(33);
t1225 = pkin(2) * t1231 + pkin(10) * t1235;
t1007 = pkin(11) * t1075 + t1225;
t1266 = t1083 * t1007 + t1078 * t1008;
t1327 = -t1324 * pkin(4) + pkin(12) * t1268 + qJD(4) * t1266 + t1344 * t1078 - t1340;
t1242 = t1007 * t1078;
t1336 = pkin(12) * t1267 + qJD(4) * t1242 + t1338;
t1296 = pkin(11) + pkin(12);
t1044 = t1296 * t1078;
t1257 = t1078 * t1321;
t1108 = pkin(9) * t1233 + t1059;
t926 = pkin(10) * t1339 + t1108 * qJD(1);
t1101 = pkin(2) * t1288 + t1133;
t935 = qJD(2) * pkin(2) + qJD(1) * t1101 + t1051;
t1129 = -pkin(2) * t1085 - t1080 * t1290 - pkin(1);
t995 = t1129 * t1074;
t982 = qJD(1) * t995;
t820 = -t1079 * t926 + t1084 * (t1073 * t982 + t1075 * t935);
t861 = pkin(3) * t938 - pkin(11) * t1321;
t1278 = t1078 * t861 + t1083 * t820;
t1335 = -pkin(12) * t1257 + qJD(4) * t1044 + t1278;
t860 = t1083 * t861;
t1334 = -pkin(12) * t1083 * t1321 + pkin(4) * t938 - t1078 * t820 + t1296 * t1247 + t860;
t1205 = qJD(3) * t1231;
t1332 = pkin(2) * t1205 + pkin(10) * t1206 - t1079 * t968 + t1230 * t969;
t1069 = t1074 ^ 2;
t1329 = 0.2e1 * pkin(1) * t1069;
t870 = -t1073 * t935 + t1075 * t982;
t792 = -pkin(3) * t1321 - pkin(11) * t938 + t870;
t821 = t1084 * t926 + t935 * t1231 + t982 * t1236;
t798 = -pkin(11) * t991 + t821;
t754 = t1078 * t792 + t1083 * t798;
t738 = pkin(12) * t878 + t754;
t1258 = t1077 * t738;
t753 = -t1078 * t798 + t1083 * t792;
t737 = pkin(12) * t1135 + t753;
t727 = pkin(4) * t932 + t737;
t704 = t1082 * t727 - t1258;
t702 = -pkin(5) * t929 - t704;
t1328 = t702 * t801;
t933 = t1082 * t1014 + t1015 * t1077;
t1275 = -qJD(5) * t933 - t1267 * t1077 + t1082 * t1268;
t934 = -t1014 * t1077 + t1015 * t1082;
t1274 = qJD(5) * t934 + t1077 * t1268 + t1267 * t1082;
t1029 = t1077 * t1078 - t1082 * t1083;
t858 = t1029 * t1321;
t1309 = qJD(4) + qJD(5);
t965 = t1309 * t1029;
t1270 = -t965 + t858;
t1030 = t1077 * t1083 + t1078 * t1082;
t1269 = (t1309 - t1321) * t1030;
t1210 = qJD(1) * t1234;
t1252 = -(-pkin(3) * t1210 - t1003 * t1084) * t1073 + t1332;
t1013 = -t1075 * t1288 + t1201;
t1319 = t1067 * t902 - t1068 * t972;
t1170 = t1288 * t1295;
t1018 = -t1080 * t1170 + t1085 * t1086;
t1103 = t1086 * t1080 + t1085 * t1170;
t1214 = t1074 * t1295;
t1171 = t1073 * t1214;
t904 = t1018 * t1084 + (-t1075 * t1103 + t1171) * t1079;
t974 = t1073 * t1103 + t1075 * t1214;
t851 = -t1067 * t904 + t1068 * t974;
t1192 = t1288 * t1073;
t1158 = t1079 * t1192;
t961 = t1158 + t1106;
t1125 = -g(3) * (-t1013 * t1068 - t1067 * t961) - g(2) * t1319 - g(1) * t851;
t1204 = qJD(3) * t1230;
t1250 = qJD(3) * t1079;
t1042 = t1051 * qJD(2);
t1172 = pkin(1) * t1188;
t1195 = pkin(9) * t1190 + t1080 * t1172 + t1042;
t1102 = -pkin(9) * t1164 + t1195;
t1251 = qJD(1) * qJD(2);
t1306 = -t1080 * t1251 + t1222;
t867 = (t1074 * t1075 * t1306 + t1054 * t1073) * pkin(10) + t1102;
t1109 = -t1059 * t1251 + t1085 * t1172;
t1211 = t1085 * t1251;
t1116 = -t1211 - t1223;
t1110 = t1116 * pkin(9);
t873 = t1054 * pkin(2) + (t1116 * t1289 + t1110) * t1074 + t1109;
t907 = (qJDD(1) * t1129 + t1132 * t1251) * t1074;
t748 = t1084 * t867 + t935 * t1204 + t982 * t1206 + t873 * t1231 + t907 * t1236 - t1250 * t926;
t735 = pkin(11) * t943 + t748;
t812 = -t1073 * t873 + t1075 * t907;
t743 = pkin(3) * t845 - pkin(11) * t844 + t812;
t692 = -qJD(4) * t754 - t1078 * t735 + t1083 * t743;
t685 = pkin(4) * t839 - pkin(12) * t768 + t692;
t1121 = -t1078 * t743 - t1083 * t735 - t792 * t1247 + t1248 * t798;
t689 = -pkin(12) * t769 - t1121;
t1254 = t1082 * t738;
t705 = t1077 * t727 + t1254;
t678 = -qJD(5) * t705 - t1077 * t689 + t1082 * t685;
t676 = -pkin(5) * t838 - t678;
t1112 = t1125 - t676;
t1241 = t1017 * t1079;
t899 = t1016 * t1230 + t1084 * t1200 + t1241;
t1099 = t1103 * t1084;
t903 = t1018 * t1079 + t1075 * t1099 - t1084 * t1171;
t1157 = t1084 * t1192;
t1196 = t1075 * t1226;
t960 = -t1074 * t1196 - t1157 + t1197;
t1124 = g(1) * t903 + g(2) * t899 + g(3) * t960;
t1105 = t1124 * t1068;
t1161 = -t821 + (t1248 - t1257) * pkin(4);
t1066 = -pkin(4) * t1083 - pkin(3);
t951 = pkin(5) * t1029 - pkin(13) * t1030 + t1066;
t1045 = t1296 * t1083;
t985 = -t1044 * t1077 + t1045 * t1082;
t1322 = -t702 * qJD(6) * t1030 - (-pkin(5) * t1269 + pkin(13) * t1270 + qJD(6) * t985 - t1161) * t1253 - t951 * t718 - t1105;
t758 = pkin(5) * t1136 + pkin(13) * t801;
t1174 = -t1077 * t685 - t1082 * t689 - t727 * t1245 + t738 * t1246;
t797 = pkin(3) * t991 - t820;
t770 = -pkin(4) * t878 + t797;
t852 = t1067 * t974 + t1068 * t904;
t890 = -t1013 * t1067 + t1068 * t961;
t1320 = g(1) * t852 - g(2) * t850 + g(3) * t890 + t770 * t801 + t1174;
t1318 = t1078 * t972 + t1083 * t902;
t1317 = t1078 * t902 - t1083 * t972;
t1271 = pkin(4) * t1267 + t1252;
t1064 = pkin(4) * t1077 + pkin(13);
t1313 = (-pkin(4) * t1135 + qJD(6) * t1064 + t758) * t1253;
t1312 = -qJD(5) * t985 + t1077 * t1335 - t1082 * t1334;
t1134 = -t1044 * t1082 - t1045 * t1077;
t1311 = -qJD(5) * t1134 + t1077 * t1334 + t1082 * t1335;
t896 = t1013 * t1083 + t1078 * t961;
t897 = -t1013 * t1078 + t1083 * t961;
t829 = -t1077 * t896 + t1082 * t897;
t952 = t960 * t1081;
t1310 = -t1076 * t829 + t952;
t1175 = t994 - t1242;
t874 = -pkin(4) * t1235 - pkin(12) * t1015 + t1175;
t881 = -pkin(12) * t1014 + t1266;
t1307 = t1077 * t1327 + t1336 * t1082 - t874 * t1245 + t1246 * t881;
t703 = pkin(13) * t929 + t705;
t725 = pkin(5) * t801 - pkin(13) * t1136 + t770;
t1144 = t1076 * t703 - t1081 * t725;
t1305 = t1136 * t1144 + t702 * t1244;
t1303 = -t770 * t1136 + t1125 + t678;
t687 = t1076 * t725 + t1081 * t703;
t1301 = -t1076 * t1112 + t687 * t1136 + t702 * t1243;
t967 = t1060 + t1101;
t884 = -t1073 * t967 + t1075 * t995;
t818 = pkin(3) * t960 - pkin(11) * t961 + t884;
t954 = (t1192 + t1198) * pkin(10) + t1108;
t1217 = t1084 * t954 + t967 * t1231 + t995 * t1236;
t826 = -pkin(11) * t1013 + t1217;
t1182 = -t1078 * t826 + t1083 * t818;
t747 = pkin(4) * t960 - pkin(12) * t897 + t1182;
t1277 = t1078 * t818 + t1083 * t826;
t752 = -pkin(12) * t896 + t1277;
t1139 = t1077 * t747 + t1082 * t752;
t1004 = qJD(2) * t1114;
t1052 = qJD(2) * t1060;
t970 = qJD(2) * t1133 + t1052;
t971 = t1096 * qJD(2);
t1104 = t1004 * t1236 + t1084 * t970 + t967 * t1204 + t995 * t1206 + t971 * t1231 - t1250 * t954;
t1167 = qJD(2) * t1202;
t772 = pkin(11) * t1167 + t1104;
t887 = qJD(3) * t1158 + t1074 * t1090;
t888 = qJD(3) * t1157 + ((t1196 - t1229) * qJD(3) + t1126 * qJD(2)) * t1074;
t892 = t1075 * t1004 - t1073 * t971;
t781 = pkin(3) * t887 - pkin(11) * t888 + t892;
t1095 = -qJD(4) * t1277 - t1078 * t772 + t1083 * t781;
t810 = qJD(2) * t1166 - qJD(4) * t896 + t1083 * t888;
t701 = pkin(4) * t887 - pkin(12) * t810 + t1095;
t1120 = t1078 * t781 + t1083 * t772 + t818 * t1247 - t1248 * t826;
t809 = -qJD(2) * t1165 + qJD(4) * t897 + t1078 * t888;
t707 = -pkin(12) * t809 + t1120;
t1298 = -qJD(5) * t1139 - t1077 * t707 + t1082 * t701;
t1087 = qJD(1) ^ 2;
t1287 = t878 * t932;
t1286 = t1135 * t932;
t1273 = t1077 * t874 + t1082 * t881;
t1283 = -pkin(5) * t1324 + qJD(5) * t1273 - t1077 * t1336 + t1082 * t1327;
t905 = t1076 * t934 + t1081 * t1235;
t1280 = -qJD(6) * t905 + t1076 * t1324 + t1081 * t1275;
t1199 = t1076 * t1235;
t1279 = -qJD(6) * t1199 + t1076 * t1275 - t1081 * t1324 + t1243 * t934;
t1272 = pkin(5) * t938 - t1312;
t1263 = qJD(4) * t932;
t1262 = t1054 * MDP(8);
t1259 = t1076 * t960;
t1249 = qJD(3) * t1084;
t1240 = t1030 * t1076;
t1239 = t1030 * t1081;
t1238 = t1067 * t1073;
t1237 = t1069 * t1087;
t1070 = t1080 ^ 2;
t1224 = -t1085 ^ 2 + t1070;
t1203 = t1085 * t1237;
t675 = pkin(13) * t838 - t1174;
t749 = -t1079 * t867 - t935 * t1205 - t982 * t1207 + t873 * t1230 + t907 * t1235 - t926 * t1249;
t736 = -pkin(3) * t943 - t749;
t716 = pkin(4) * t769 + t736;
t682 = pkin(5) * t720 - pkin(13) * t719 + t716;
t1187 = -t1076 * t675 + t1081 * t682;
t807 = -t1076 * t858 - t1081 * t938;
t1185 = t1076 * t965 + t807;
t808 = t1076 * t938 - t1081 * t858;
t1180 = t1081 * t965 + t808;
t1177 = t1083 * t932;
t1173 = -t1079 * t970 - t967 * t1205 - t995 * t1207 - t954 * t1249;
t708 = t1077 * t737 + t1254;
t1160 = pkin(4) * t1246 - t708;
t1159 = t1074 * t1087 * t1288;
t1006 = t1055 + (-pkin(2) * t1084 - pkin(3)) * t1075;
t942 = pkin(4) * t1014 + t1006;
t827 = pkin(5) * t933 - pkin(13) * t934 + t942;
t1155 = -pkin(13) * t1324 - qJD(6) * t827 + t1307;
t796 = -pkin(13) * t1235 + t1273;
t1153 = -pkin(5) * t1274 + pkin(13) * t1275 + qJD(6) * t796 - t1271;
t1148 = 0.2e1 * t1193 + qJD(2);
t1146 = t1073 * t995 + t1075 * t967;
t1145 = t1076 * t682 + t1081 * t675;
t711 = pkin(13) * t960 + t1139;
t945 = t1079 * t954;
t825 = pkin(3) * t1013 - t1084 * t1146 + t945;
t782 = pkin(4) * t896 + t825;
t828 = t1077 * t897 + t1082 * t896;
t730 = pkin(5) * t828 - pkin(13) * t829 + t782;
t1143 = t1076 * t730 + t1081 * t711;
t1142 = -t1076 * t711 + t1081 * t730;
t788 = t1081 * t829 + t1259;
t1140 = -t1077 * t752 + t1082 * t747;
t1137 = -t1077 * t881 + t1082 * t874;
t1130 = -pkin(11) * t839 + t797 * t932;
t1123 = -g(1) * t904 + g(2) * t902 - g(3) * t961;
t922 = -t1016 * t1084 - t1017 * t1231;
t924 = -t1018 * t1231 - t1099;
t1122 = g(1) * t924 + g(2) * t922 + g(3) * t999;
t1119 = t1077 * t701 + t1082 * t707 + t747 * t1245 - t1246 * t752;
t1115 = t1074 * (t1188 + t1054);
t1113 = -qJD(6) * t1240 - t1180;
t1107 = g(1) * t1018 + g(2) * t1017 + g(3) * t1234;
t1100 = t1007 * t1263 + t1122;
t1097 = pkin(11) * t1263 - t1124 + t736;
t709 = t1082 * t737 - t1258;
t1093 = -t1064 * t718 + t1328 - (-pkin(4) * t1245 + t709) * t1253;
t1091 = t1149 * t1108;
t773 = -t971 * t1230 + (-pkin(3) * t1209 - t1004 * t1084) * t1073 - t1173;
t742 = pkin(4) * t809 + t773;
t1088 = t676 * t1030 - t702 * t965 - t985 * t718 - (pkin(13) * t938 - qJD(6) * t951 + t1311) * t1253 + t1123;
t1065 = -pkin(4) * t1082 - pkin(5);
t925 = t1067 * t1202 + t1068 * t999;
t923 = t1018 * t1230 - t1079 * t1103;
t921 = -t1016 * t1079 + t1017 * t1230;
t906 = t1081 * t934 - t1199;
t901 = -t1241 + (-t1016 * t1075 - t1200) * t1084;
t876 = t1018 * t1238 + t1068 * t924;
t875 = t1017 * t1238 + t1068 * t922;
t856 = t1078 * t974 + t1083 * t904;
t855 = -t1078 * t904 + t1083 * t974;
t795 = pkin(5) * t1235 - t1137;
t791 = t1076 * t903 + t1081 * t852;
t790 = -t1076 * t852 + t1081 * t903;
t732 = qJD(5) * t829 + t1077 * t810 + t1082 * t809;
t731 = -qJD(5) * t828 - t1077 * t809 + t1082 * t810;
t713 = qJD(6) * t788 + t1076 * t731 - t887 * t1081;
t712 = qJD(6) * t1310 + t1076 * t887 + t1081 * t731;
t710 = -pkin(5) * t960 - t1140;
t690 = pkin(5) * t732 - pkin(13) * t731 + t742;
t680 = -pkin(5) * t887 - t1298;
t679 = pkin(13) * t887 + t1119;
t673 = -t687 * qJD(6) + t1187;
t672 = -qJD(6) * t1144 + t1145;
t1 = [(g(1) * t1317 - g(2) * t855 - t1120 * t932 + t1121 * t960 - t773 * t1135 - t1277 * t839 + t736 * t897 - t754 * t887 + t825 * t768 + t797 * t810) * MDP(24) + (-qJD(2) * t1091 + (-pkin(9) * t1234 + t1060) * t1054 + (t1074 * t1110 + t1109) * t1288 + g(1) * t1017 - g(2) * t1018 + t1306 * t1329) * MDP(9) + (-(-pkin(9) * t1209 + t1052) * t1149 - t1108 * t1054 - t1102 * t1288 - g(1) * t1016 + g(2) * t1103 + t1116 * t1329) * MDP(10) + (-g(1) * t1318 - g(2) * t856 + t1095 * t932 + t1182 * t839 + t692 * t960 + t736 * t896 + t753 * t887 + t825 * t769 - t773 * t878 + t797 * t809) * MDP(23) + (-t1253 * t712 + t696 * t828 + t718 * t788 + t732 * t778) * MDP(34) + (-t1253 * t732 + t718 * t828) * MDP(36) + (t1135 * t809 - t768 * t896 - t769 * t897 + t810 * t878) * MDP(19) + (-t1135 * t810 + t768 * t897) * MDP(18) + (-t1135 * t887 + t768 * t960 + t810 * t932 + t839 * t897) * MDP(20) + (t1253 * t713 + t1310 * t718 - t697 * t828 - t732 * t776) * MDP(35) + (t1310 * t696 - t697 * t788 - t712 * t776 - t713 * t778) * MDP(33) + (g(1) * t1319 - g(2) * t851 - t1119 * t929 + t742 * t1136 - t1139 * t838 + t1174 * t960 - t705 * t887 + t716 * t829 + t782 * t719 + t770 * t731) * MDP(31) + (t696 * t788 + t712 * t778) * MDP(32) + (-t720 * t960 - t732 * t929 - t801 * t887 - t828 * t838) * MDP(28) + (t838 * t960 + t887 * t929) * MDP(29) + (-t769 * t960 - t809 * t932 - t839 * t896 + t878 * t887) * MDP(21) + (t839 * t960 + t887 * t932) * MDP(22) + (-g(1) * t850 - g(2) * t852 + t1140 * t838 + t1298 * t929 + t678 * t960 + t704 * t887 + t716 * t828 + t782 * t720 + t770 * t732 + t742 * t801) * MDP(30) + (0.2e1 * (t1080 * t1222 - t1224 * t1251) * MDP(5) + (qJDD(1) * t1070 + 0.2e1 * t1080 * t1211) * MDP(4)) * t1069 + (-(-qJD(6) * t1143 - t1076 * t679 + t1081 * t690) * t1253 + t1142 * t718 + t673 * t828 - t1144 * t732 + t680 * t776 + t710 * t697 - t676 * t1310 + t702 * t713 - g(1) * (t1076 * t901 + t1341) - g(2) * t791) * MDP(37) + (g(1) * t901 + g(2) * t903 + t748 * t1013 + t1104 * t991 - t1167 * t821 - t1217 * t943 + t812 * t961 + t884 * t844 + t870 * t888 + t892 * t938) * MDP(17) + (t1080 * t1115 + t1148 * t1208) * MDP(6) + (t1085 * t1115 - t1148 * t1209) * MDP(7) + (g(1) * t1295 - g(2) * t1086) * MDP(2) + (g(1) * t1086 + g(2) * t1295) * MDP(3) + (-t1136 * t732 - t719 * t828 - t720 * t829 - t731 * t801) * MDP(26) + (t1136 * t731 + t719 * t829) * MDP(25) + (t1136 * t887 + t719 * t960 + t731 * t929 + t829 * t838) * MDP(27) + (t1013 * t845 + t1167 * t1321 + t887 * t991 - t943 * t960) * MDP(14) + (t1321 * t888 - t844 * t960 - t845 * t961 - t887 * t938) * MDP(12) + (-t1173 * t991 - t945 * t943 - t749 * t1013 + t820 * t1167 - t892 * t1321 + t884 * t845 + t812 * t960 + t870 * t887 - g(1) * t902 - g(2) * t904 + (-(t1004 * t1073 + t1075 * t971) * t991 + t1146 * t943) * t1084) * MDP(16) + (-t1013 * t943 - t1167 * t991) * MDP(15) + (-t1013 * t844 + t1167 * t938 - t888 * t991 + t943 * t961) * MDP(13) + qJDD(1) * MDP(1) + ((qJD(6) * t1142 + t1076 * t690 + t1081 * t679) * t1253 - t1143 * t718 - t672 * t828 - t687 * t732 + t680 * t778 + t710 * t696 + t676 * t788 + t702 * t712 - g(1) * (t1081 * t901 - t1342) - g(2) * t790) * MDP(38) + t1288 * t1262 + (t844 * t961 + t888 * t938) * MDP(11); (pkin(1) * t1080 * t1237 - pkin(9) * t1323 + g(1) * t1103 + g(2) * t1016 - g(3) * t1233 + qJD(1) * t1091 + t1109) * MDP(9) + ((-t1076 * t796 + t1081 * t827) * t718 + t673 * t933 + t795 * t697 + t676 * t905 - g(1) * (t1076 * t923 + t1081 * t876) - g(2) * (t1076 * t921 + t1081 * t875) - g(3) * (t1076 * t998 + t1081 * t925) - (t1076 * t1155 - t1081 * t1153) * t1253 + t1283 * t776 + t1279 * t702 - t1274 * t1144) * MDP(37) + (-(t1076 * t827 + t1081 * t796) * t718 - t672 * t933 + t795 * t696 + t676 * t906 - g(1) * (-t1076 * t876 + t1081 * t923) - g(2) * (-t1076 * t875 + t1081 * t921) - g(3) * (-t1076 * t925 + t1081 * t998) - (t1076 * t1153 + t1081 * t1155) * t1253 + t1283 * t778 + t1280 * t702 - t1274 * t687) * MDP(38) + (-t1253 * t1274 + t718 * t933) * MDP(36) + (t1253 * t1279 - t1274 * t776 - t697 * t933 - t718 * t905) * MDP(35) + (-t1253 * t1280 + t1274 * t778 + t696 * t933 + t718 * t906) * MDP(34) + (t1015 * t839 + t1135 * t988 + t1268 * t932 + (-t1084 * t768 - t1135 * t1250) * t1073) * MDP(20) + (t1015 * t768 - t1135 * t1268) * MDP(18) + (-t1014 * t768 - t1015 * t769 + t1135 * t1267 + t1268 * t878) * MDP(19) + (-t1273 * t838 + t942 * t719 + t716 * t934 + t705 * t988 + t1307 * t929 + t1271 * t1136 + t1275 * t770 + t1122 * t1067 + (-t1068 * t1107 - t1084 * t1174 - t1250 * t705) * t1073) * MDP(31) + t1224 * MDP(5) * t1237 + (pkin(1) * t1203 + (-pkin(9) * t1210 + t1051) * t1193 + t1042 + t1107 - t1195) * MDP(10) + (t801 * t988 - t838 * t933 - t1274 * t929 + (t1084 * t720 - t1250 * t801) * t1073) * MDP(28) + (t1075 * t844 + t989 * t991 + (t1079 * t943 - t1210 * t938 - t1249 * t991) * t1073) * MDP(13) + (-t938 * t989 + (t1079 * t844 + t1249 * t938) * t1073) * MDP(11) + (-t932 * t988 + (-t1084 * t839 + t1250 * t932) * t1073) * MDP(22) + (-t929 * t988 + (-t1084 * t838 + t1250 * t929) * t1073) * MDP(29) + (t1075 * t943 + t1169 * t991) * MDP(15) + t1262 + (t1136 * t1275 + t719 * t934) * MDP(25) + (-t1136 * t988 + t838 * t934 + t1275 * t929 + (-t1084 * t719 + t1136 * t1250) * t1073) * MDP(27) + (-t1136 * t1274 - t1275 * t801 - t719 * t933 - t720 * t934) * MDP(26) + (-t1266 * t839 + t1006 * t768 + t736 * t1015 + t754 * t988 + t1338 * t932 - t1252 * t1135 + t1268 * t797 + t1100 * t1078 + (-t1083 * t1107 - t1084 * t1121 - t1250 * t754) * t1073) * MDP(24) + (-t1085 * t1159 + t1191) * MDP(6) + (-t1075 * t845 - t988 * t991 + (t1084 * t943 - t1210 * t1321 + t1250 * t991) * t1073) * MDP(14) + (-t1321 * t989 + t938 * t988 + (-t1079 * t845 + t1084 * t844 + (-t1079 * t938 + t1084 * t1321) * qJD(3)) * t1073) * MDP(12) + (t1308 * t943 + t749 * t1075 + t891 * t1321 - t870 * t988 + t1332 * t991 + (-t820 * t1210 + t870 * t1250 - pkin(2) * t845 + (t1003 * t991 - t812) * t1084) * t1073 - t1122) * MDP(16) + (-t1014 * t839 - t878 * t988 - t1267 * t932 + (t1084 * t769 + t1250 * t878) * t1073) * MDP(21) + (t1280 * t778 + t696 * t906) * MDP(32) + (-t1279 * t778 - t1280 * t776 - t696 * t905 - t697 * t906) * MDP(33) - t1080 * MDP(4) * t1203 + (t1137 * t838 + t942 * t720 + t716 * t933 - t704 * t988 - g(1) * t876 - g(2) * t875 - g(3) * t925 + ((-qJD(5) * t881 - t1327) * t1082 + (-qJD(5) * t874 + t1336) * t1077) * t929 + t1271 * t801 + t1274 * t770 + (-t1084 * t678 + t1250 * t704) * t1073) * MDP(30) + (t1175 * t839 + t1006 * t769 + t736 * t1014 - t753 * t988 + ((-qJD(4) * t1008 - t1344) * t1078 + t1340) * t932 - t1252 * t878 + t1267 * t797 - t1100 * t1083 + (-t1078 * t1107 - t692 * t1084 + t1250 * t753) * t1073) * MDP(23) + (-t1225 * t943 - t748 * t1075 - t891 * t938 - t870 * t989 + g(1) * t923 + g(2) * t921 + g(3) * t998 + t1346 * t991 + (-pkin(2) * t844 + t1079 * t812 + t1210 * t821 + t1249 * t870) * t1073) * MDP(17) + (t1080 * t1159 + t1190) * MDP(7); (t716 * t1029 + t1066 * t720 + t1134 * t838 + t1161 * t801 + t1269 * t770 + t1312 * t929 + t1105) * MDP(30) + (-t718 * t1240 - t1029 * t697 - t1269 * t776 - (-qJD(6) * t1239 + t1185) * t1253) * MDP(35) + (t1029 * t696 - t1113 * t1253 + t1239 * t718 + t1269 * t778) * MDP(34) + (t1029 * t718 - t1253 * t1269) * MDP(36) + (-pkin(3) * t768 + t1078 * t1097 + t1083 * t1130 + t1135 * t821 + t1278 * t932) * MDP(24) + (t1078 * t768 - t1135 * t1177) * MDP(18) + (t716 * t1030 + t1066 * t719 - t1124 * t1067 + t1161 * t1136 + t1270 * t770 + t1311 * t929 - t985 * t838) * MDP(31) - t845 * MDP(14) + (t1185 * t778 + t1180 * t776 + (-t694 - t1081 * t697 + (t1076 * t776 - t1081 * t778) * qJD(6)) * t1030) * MDP(33) + (t1113 * t778 + t1239 * t696) * MDP(32) + (-t672 * t1029 + t1076 * t1322 + t1088 * t1081 - t1134 * t696 - t1269 * t687 + t1272 * t778 - t702 * t808) * MDP(38) + (t673 * t1029 + t1088 * t1076 - t1081 * t1322 - t1134 * t697 - t1269 * t1144 + t1272 * t776 - t702 * t807) * MDP(37) + (t1030 * t719 + t1136 * t1270) * MDP(25) + (-t1029 * t719 - t1030 * t720 - t1136 * t1269 - t1270 * t801) * MDP(26) + (-MDP(11) * t1321 + MDP(12) * t938 - MDP(14) * t991 - MDP(16) * t870 + MDP(20) * t1135 - MDP(21) * t878 - MDP(22) * t932 - MDP(23) * t753 + MDP(24) * t754 - MDP(27) * t1136 + MDP(28) * t801 - MDP(29) * t929 - MDP(30) * t704 + MDP(31) * t705) * t938 + (-t1321 * t870 - t820 * t991 - t1123 - t748) * MDP(17) + (t1321 * t991 + t844) * MDP(13) - t1321 ^ 2 * MDP(12) + (-pkin(3) * t769 + t821 * t878 - t860 * t932 + (t820 * t932 + t1130) * t1078 - t1097 * t1083) * MDP(23) + t943 * MDP(15) + (t1078 * t839 + t1177 * t932) * MDP(20) + (-t1078 * t932 ^ 2 + t1083 * t839) * MDP(21) + (-t821 * t991 + t1124 + t749) * MDP(16) + (t1030 * t838 + t1270 * t929) * MDP(27) + (-t1029 * t838 - t1269 * t929) * MDP(28) + ((t768 + t1287) * t1083 + (-t769 + t1286) * t1078) * MDP(19); t839 * MDP(22) + (t1135 ^ 2 - t878 ^ 2) * MDP(19) + (t708 * t929 + (t1082 * t838 + t1135 * t801 - t1246 * t929) * pkin(4) + t1303) * MDP(30) + (t709 * t929 + (-t1077 * t838 + t1135 * t1136 - t1245 * t929) * pkin(4) + t1320) * MDP(31) + (t768 - t1287) * MDP(20) + (t1065 * t697 + t1160 * t776 + t1093 * t1076 + (t1112 + t1313) * t1081 + t1305) * MDP(37) + (t1065 * t696 - t1076 * t1313 + t1081 * t1093 + t1160 * t778 + t1301) * MDP(38) + t1135 * t878 * MDP(18) + (-t769 - t1286) * MDP(21) + (-g(1) * t855 - g(2) * t1317 + g(3) * t896 + t797 * t1135 + t754 * t932 + t692) * MDP(23) + (g(1) * t856 - g(2) * t1318 + g(3) * t897 + t753 * t932 - t797 * t878 + t1121) * MDP(24) + t1337; (t705 * t929 + t1303) * MDP(30) + (t704 * t929 + t1320) * MDP(31) + (-pkin(5) * t697 - t705 * t776 + (-pkin(13) * t718 - t1253 * t704 + t1328) * t1076 + (-(-pkin(13) * qJD(6) - t758) * t1253 + t1112) * t1081 + t1305) * MDP(37) + (-pkin(5) * t696 - (t1076 * t758 + t1081 * t704) * t1253 - t705 * t778 + t702 * t1325 + (-t1244 * t1253 - t715) * pkin(13) + t1301) * MDP(38) + t1337; t778 * t776 * MDP(32) + (-t776 ^ 2 + t778 ^ 2) * MDP(33) + (-t1253 * t776 + t1218) * MDP(34) + (-t1253 * t778 - t1186) * MDP(35) + t718 * MDP(36) + (-t687 * t1253 - t702 * t778 - g(1) * t790 - g(2) * (t1081 * t899 + t1342) - g(3) * (-t1076 * t890 + t952) + t1187) * MDP(37) + (t1144 * t1253 + t702 * t776 + g(1) * t791 - g(2) * (-t1076 * t899 + t1341) - g(3) * (-t1081 * t890 - t1259) - t1145) * MDP(38) + (-MDP(34) * t1261 - MDP(35) * t778 - MDP(37) * t687 + MDP(38) * t1144) * qJD(6);];
tau  = t1;
