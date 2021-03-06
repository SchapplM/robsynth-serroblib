% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR10_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR10_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR10_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 06:17:55
% EndTime: 2019-03-10 06:19:01
% DurationCPUTime: 37.74s
% Computational Cost: add. (55446->989), mult. (169578->1390), div. (0->0), fcn. (145876->16), ass. (0->394)
t1060 = cos(qJ(3));
t1050 = cos(pkin(7));
t1208 = t1050 * t1060;
t1174 = qJD(3) * t1208;
t1055 = sin(qJ(3));
t1209 = t1050 * t1055;
t1047 = sin(pkin(7));
t1215 = t1047 * t1055;
t1051 = cos(pkin(6));
t1061 = cos(qJ(2));
t1275 = pkin(1) * t1061;
t1193 = t1051 * t1275;
t1031 = qJD(1) * t1193;
t1056 = sin(qJ(2));
t1048 = sin(pkin(6));
t1183 = pkin(11) * t1050 + pkin(10);
t1125 = t1048 * t1183;
t1099 = t1056 * t1125;
t952 = -qJD(1) * t1099 + t1031;
t1276 = pkin(1) * t1056;
t1194 = t1051 * t1276;
t1065 = -t1061 * t1125 - t1194;
t953 = t1065 * qJD(1);
t1274 = pkin(11) * t1047;
t1077 = t1048 * (pkin(2) * t1056 - t1061 * t1274);
t989 = qJD(1) * t1077;
t1302 = -pkin(2) * t1174 + t1060 * t952 + t953 * t1209 + t989 * t1215;
t1049 = cos(pkin(8));
t1046 = sin(pkin(8));
t1213 = t1048 * t1056;
t1164 = t1047 * t1213;
t1136 = t1046 * t1164;
t1202 = t1056 * t1060;
t1203 = t1055 * t1061;
t1087 = t1050 * t1202 + t1203;
t1235 = qJD(1) * t1048;
t977 = t1087 * t1235;
t1069 = qJD(1) * t1136 - t1049 * t977;
t1272 = pkin(12) * t1049;
t1126 = t1047 * (-pkin(11) - t1272);
t1098 = t1055 * t1126;
t1301 = -pkin(12) * t1069 + qJD(3) * t1098 - t1302;
t1214 = t1047 * t1060;
t1127 = -t1055 * t952 + t953 * t1208 + t989 * t1214;
t1141 = qJD(1) * t1164;
t1199 = t1060 * t1061;
t1205 = t1055 * t1056;
t1085 = -t1050 * t1205 + t1199;
t978 = t1085 * t1235;
t813 = pkin(3) * t1141 - t1272 * t978 + t1127;
t1038 = pkin(2) * t1209;
t955 = (t1060 * t1126 - t1038) * qJD(3);
t1300 = t813 - t955;
t1273 = pkin(12) * t1046;
t886 = -t1047 * t953 + t1050 * t989;
t842 = pkin(3) * t977 - t1273 * t978 + t886;
t988 = (pkin(3) * t1055 - t1060 * t1273) * t1047 * qJD(3);
t1299 = t842 - t988;
t1054 = sin(qJ(4));
t1059 = cos(qJ(4));
t1216 = t1046 * t1059;
t1165 = t1050 * t1216;
t1200 = t1059 * t1060;
t1207 = t1054 * t1055;
t1279 = t1049 * t1200 - t1207;
t1253 = -t1054 * t1069 - t1059 * t978 + qJD(4) * t1165 + (t1279 * qJD(4) + (-t1049 * t1207 + t1200) * qJD(3)) * t1047;
t1204 = t1055 * t1059;
t1206 = t1054 * t1060;
t1088 = t1049 * t1206 + t1204;
t1102 = t1059 * t1136;
t1218 = t1046 * t1054;
t1167 = t1050 * t1218;
t1210 = t1049 * t1059;
t1252 = qJD(1) * t1102 - t1054 * t978 - t1210 * t977 + qJD(4) * t1167 + (t1088 * qJD(4) + (t1049 * t1204 + t1206) * qJD(3)) * t1047;
t1176 = qJD(3) * t1215;
t1135 = t1049 * t1164;
t1101 = qJD(1) * t1135;
t900 = -t1046 * t977 - t1101;
t1096 = t1046 * t1176 + t900;
t1251 = t1046 * t1300 - t1049 * t1299;
t1172 = qJD(4) * t1216;
t1211 = t1049 * t1054;
t1160 = t1048 * t1205;
t1134 = qJD(1) * t1160;
t1212 = t1048 * t1061;
t1161 = t1050 * t1212;
t1140 = qJD(1) * t1161;
t1234 = qJD(1) * t1051;
t1033 = qJD(2) + t1234;
t1220 = t1033 * t1047;
t907 = -t1134 + (t1140 + t1220) * t1060;
t1086 = t1050 * t1203 + t1202;
t1075 = t1086 * t1048;
t1169 = t1033 * t1215;
t908 = qJD(1) * t1075 + t1169;
t857 = t1059 * t907 - t1211 * t908;
t1298 = t1172 - t857;
t1034 = pkin(12) * t1218;
t1282 = pkin(3) * t1210 - t1034;
t903 = pkin(11) * t1220 + (t1183 * t1212 + t1194) * qJD(1);
t906 = pkin(2) * t1033 + t952;
t985 = (-pkin(2) * t1061 - t1056 * t1274 - pkin(1)) * t1048;
t971 = qJD(1) * t985;
t829 = -t1055 * t903 + t906 * t1208 + t971 * t1214;
t792 = -t1272 * t908 + t829;
t1244 = t1049 * t907;
t830 = t1055 * (t1047 * t971 + t1050 * t906) + t1060 * t903;
t793 = -pkin(12) * t1244 - t830;
t1246 = t1046 * t907;
t858 = pkin(3) * t908 - pkin(12) * t1246;
t1297 = t1282 * qJD(4) - t1059 * t792 - t793 * t1211 - t858 * t1218;
t1170 = qJD(4) * t1210;
t1229 = qJD(4) * t1054;
t1196 = pkin(11) * t1214 + t1038;
t931 = (t1046 * t1050 + t1049 * t1214) * pkin(12) + t1196;
t1039 = pkin(2) * t1208;
t950 = pkin(3) * t1050 + t1039 + t1098;
t984 = (-pkin(3) * t1060 - t1055 * t1273 - pkin(2)) * t1047;
t1296 = t1059 * t1301 + t950 * t1170 + t984 * t1172 - t1300 * t1211 - t1299 * t1218 - t1229 * t931;
t1245 = t1046 * t908;
t1295 = pkin(13) * t1245 - t1297;
t1294 = -pkin(13) * t1096 - t1296;
t750 = -t1046 * t793 + t1049 * t858;
t856 = t1054 * t907 + t1210 * t908;
t1291 = pkin(4) * t856 - pkin(13) * t857 + t750 - (pkin(4) * t1054 - pkin(13) * t1059) * t1046 * qJD(4);
t1290 = -t1252 * pkin(4) + pkin(13) * t1253 - t1251;
t1053 = sin(qJ(5));
t1058 = cos(qJ(5));
t1001 = t1046 * t1214 - t1049 * t1050;
t941 = t1047 * t1088 + t1167;
t889 = t1058 * t1001 + t1053 * t941;
t1259 = -qJD(5) * t889 + t1096 * t1053 + t1058 * t1253;
t890 = -t1001 * t1053 + t1058 * t941;
t1258 = qJD(5) * t890 + t1053 * t1253 - t1096 * t1058;
t1005 = -t1058 * t1049 + t1053 * t1218;
t1250 = -qJD(5) * t1005 - t1053 * t1245 + t1058 * t1298;
t1217 = t1046 * t1058;
t1006 = t1049 * t1053 + t1054 * t1217;
t1249 = qJD(5) * t1006 + t1053 * t1298 + t1217 * t908;
t1197 = pkin(3) * t1211 + pkin(12) * t1216;
t1289 = t1197 * qJD(4) - t1054 * t792;
t1163 = t1047 * t1212;
t1022 = qJD(1) * t1163;
t1145 = t1033 * t1050 - t1022;
t1122 = qJD(3) + t1145;
t1091 = t1046 * t1122;
t947 = t1059 * t1091;
t839 = t1054 * t908 - t907 * t1210 - t947;
t837 = qJD(5) + t839;
t1066 = t1091 + t1244;
t841 = t1054 * t1066 + t1059 * t908;
t873 = -t1049 * t1122 - qJD(4) + t1246;
t779 = t1053 * t841 + t1058 * t873;
t778 = qJD(6) + t779;
t1173 = qJD(4) * t1218;
t1288 = t1173 - t856;
t1171 = qJD(4) * t1211;
t1228 = qJD(4) * t1059;
t1287 = t1054 * t1301 + t950 * t1171 + t984 * t1173 + t931 * t1228;
t1179 = qJD(1) * t1213;
t1139 = t1060 * t1179;
t1232 = qJD(2) * t1050;
t1177 = qJD(2) * t1212;
t1285 = qJD(1) * t1177 + qJD(3) * t1140;
t1100 = t1139 * t1232 + t1285 * t1055 + (t1139 + t1169) * qJD(3);
t1078 = t1100 * t1049;
t1103 = qJD(2) * t1136;
t1286 = qJD(1) * t1103 - t1078;
t1260 = -t955 * t1210 + (-pkin(4) * t1176 - t1059 * t988) * t1046 - pkin(4) * t900 - (-t1046 * t842 - t1049 * t813) * t1059 + t1287;
t1247 = pkin(4) * t1245 - (-t1046 * t858 - t1049 * t793) * t1059 + t1289;
t1284 = MDP(4) * t1056;
t1283 = MDP(5) * (t1056 ^ 2 - t1061 ^ 2);
t1092 = qJD(3) * t1122;
t1224 = qJD(5) * t1058;
t1226 = qJD(5) * t1053;
t994 = pkin(13) * t1049 + t1197;
t995 = (-pkin(4) * t1059 - pkin(13) * t1054 - pkin(3)) * t1046;
t1281 = t1291 * t1053 + t1058 * t1295 - t995 * t1224 + t1226 * t994;
t878 = -t1046 * t950 + t1049 * t984;
t938 = -t1047 * t1279 - t1165;
t825 = pkin(4) * t938 - pkin(13) * t941 + t878;
t1186 = t1059 * t931 + t950 * t1211 + t984 * t1218;
t834 = -pkin(13) * t1001 + t1186;
t1280 = t1053 * t1290 + t1058 * t1294 - t825 * t1224 + t1226 * t834;
t1052 = sin(qJ(6));
t1057 = cos(qJ(6));
t961 = t1006 * t1057 - t1052 * t1216;
t1108 = t1053 * t873 - t1058 * t841;
t1175 = qJD(3) * t1214;
t1231 = qJD(3) * t1055;
t1027 = qJD(2) * t1031;
t1080 = qJD(2) * t1099;
t928 = -qJD(1) * t1080 + t1027;
t957 = t1065 * qJD(2);
t929 = qJD(1) * t957;
t990 = qJD(2) * t1077;
t980 = qJD(1) * t990;
t1073 = -t1060 * t928 - t906 * t1174 - t971 * t1175 - t929 * t1209 - t980 * t1215 + t1231 * t903;
t727 = pkin(12) * t1286 - t1073;
t1138 = qJD(2) * t1164;
t1105 = qJD(1) * t1138;
t1129 = -t1055 * t928 + t929 * t1208 + t980 * t1214;
t759 = -qJD(3) * t830 + t1129;
t874 = t1033 * t1175 + (-qJD(3) - t1232) * t1134 + t1285 * t1060;
t728 = pkin(3) * t1105 - t1272 * t874 + t759;
t875 = -t1047 * t929 + t1050 * t980;
t769 = pkin(3) * t1100 - t1273 * t874 + t875;
t782 = pkin(12) * t1066 + t830;
t783 = pkin(3) * t1122 + t792;
t871 = -t1047 * t906 + t1050 * t971;
t808 = -pkin(3) * t907 - pkin(12) * t1245 + t871;
t657 = t1059 * t727 + t783 * t1170 + t808 * t1172 + t728 * t1211 + t769 * t1218 - t1229 * t782;
t854 = -qJD(2) * t1101 - t1046 * t1100;
t655 = -pkin(13) * t854 + t657;
t699 = t1059 * t782 + t783 * t1211 + t808 * t1218;
t691 = -pkin(13) * t873 + t699;
t722 = -t1046 * t783 + t1049 * t808;
t693 = pkin(4) * t839 - pkin(13) * t841 + t722;
t665 = t1053 * t693 + t1058 * t691;
t700 = -t1046 * t728 + t1049 * t769;
t740 = qJD(4) * t947 + t1054 * t1286 + t1059 * t874 + t907 * t1170 - t1229 * t908;
t1089 = qJD(2) * t1102;
t741 = -qJD(1) * t1089 + qJD(4) * t841 + t1054 * t874 + t1059 * t1078;
t670 = pkin(4) * t741 - pkin(13) * t740 + t700;
t646 = -qJD(5) * t665 - t1053 * t655 + t1058 * t670;
t644 = -pkin(5) * t741 - t646;
t1278 = (-pkin(5) * t1108 + pkin(14) * t778) * t778 + t644;
t698 = -t1054 * t782 + t1059 * (t1046 * t808 + t1049 * t783);
t696 = -qJD(5) * t1108 + t1053 * t740 + t1058 * t854;
t1002 = -t1050 * t1051 + t1163;
t1157 = t1051 * t1214;
t1158 = t1050 * t1199;
t939 = -t1048 * t1158 - t1157 + t1160;
t1106 = -t1002 * t1046 - t1049 * t939;
t1079 = -pkin(10) * t1212 - t1194;
t932 = (t1047 * t1051 + t1161) * pkin(11) - t1079;
t920 = t1060 * t932;
t951 = (pkin(2) + t1275) * t1051 - t1099;
t1185 = t951 * t1209 + t985 * t1215 + t920;
t804 = pkin(12) * t1106 + t1185;
t1128 = -t1055 * t932 + t951 * t1208 + t985 * t1214;
t1162 = t1051 * t1215;
t940 = t1075 + t1162;
t807 = -pkin(3) * t1002 - t1272 * t940 + t1128;
t879 = -t1047 * t951 + t1050 * t985;
t824 = pkin(3) * t939 - t1273 * t940 + t879;
t1188 = t1059 * t804 + t807 * t1211 + t824 * t1218;
t888 = t1002 * t1049 - t1046 * t939;
t703 = -pkin(13) * t888 + t1188;
t734 = -t1046 * t807 + t1049 * t824;
t1240 = t1054 * t940;
t846 = t1002 * t1216 + t1210 * t939 + t1240;
t847 = t1054 * t1106 + t1059 * t940;
t708 = pkin(4) * t846 - pkin(13) * t847 + t734;
t1111 = t1053 * t708 + t1058 * t703;
t881 = qJD(3) * t1162 + (qJD(2) * t1087 + qJD(3) * t1086) * t1048;
t1068 = -t1049 * t881 + t1103;
t1032 = qJD(2) * t1193;
t956 = t1032 - t1080;
t1072 = t1060 * t956 + t951 * t1174 + t985 * t1175 + t957 * t1209 + t990 * t1215 - t1231 * t932;
t745 = pkin(12) * t1068 + t1072;
t1063 = -t1055 * t956 + t957 * t1208 + t990 * t1214 + (-t920 + (-t1047 * t985 - t1050 * t951) * t1055) * qJD(3);
t882 = qJD(3) * t1157 + ((t1158 - t1205) * qJD(3) + t1085 * qJD(2)) * t1048;
t746 = pkin(3) * t1138 - t1272 * t882 + t1063;
t887 = -t1047 * t957 + t1050 * t990;
t787 = pkin(3) * t881 - t1273 * t882 + t887;
t1071 = t1059 * t745 + t807 * t1170 + t824 * t1172 + t746 * t1211 + t787 * t1218 - t1229 * t804;
t862 = -qJD(2) * t1135 - t1046 * t881;
t662 = -pkin(13) * t862 + t1071;
t713 = -t1046 * t746 + t1049 * t787;
t752 = qJD(4) * t847 + t1054 * t882 + t1210 * t881 - t1089;
t753 = t1059 * t882 + t1068 * t1054 + (t1059 * t1106 - t1240) * qJD(4);
t679 = pkin(4) * t752 - pkin(13) * t753 + t713;
t1277 = -qJD(5) * t1111 - t1053 * t662 + t1058 * t679;
t1271 = pkin(13) * qJD(5);
t1241 = t1052 * t1108;
t729 = -t1057 * t837 - t1241;
t1270 = t729 * t778;
t731 = t1052 * t837 - t1057 * t1108;
t1269 = t731 * t778;
t1268 = t779 * t837;
t1267 = t1108 * t837;
t1257 = t1053 * t825 + t1058 * t834;
t1266 = -pkin(5) * t1252 + qJD(5) * t1257 - t1053 * t1294 + t1290 * t1058;
t757 = pkin(4) * t841 + pkin(13) * t839;
t1265 = t1053 * t757 + t1058 * t698;
t848 = t1052 * t890 - t1057 * t938;
t1262 = -qJD(6) * t848 + t1052 * t1252 + t1057 * t1259;
t849 = t1052 * t938 + t1057 * t890;
t1261 = qJD(6) * t849 + t1052 * t1259 - t1057 * t1252;
t1248 = t1053 * t995 + t1058 * t994;
t1256 = -t1288 * pkin(5) + qJD(5) * t1248 - t1053 * t1295 + t1291 * t1058;
t960 = t1006 * t1052 + t1057 * t1216;
t1255 = -qJD(6) * t960 + t1052 * t1288 + t1057 * t1250;
t1254 = qJD(6) * t961 + t1052 * t1250 - t1057 * t1288;
t1222 = qJD(6) * t1057;
t695 = -t1053 * t854 + t1058 * t740 - t873 * t1224 - t1226 * t841;
t1190 = t1052 * t741 + t1057 * t695 + t837 * t1222;
t1223 = qJD(6) * t1052;
t666 = t1108 * t1223 + t1190;
t1243 = t1052 * t666;
t1242 = t1052 * t696;
t1239 = t1057 * t696;
t1236 = -t699 + t837 * (pkin(5) * t1053 - pkin(14) * t1058);
t1233 = qJD(2) * t1047 ^ 2;
t1230 = qJD(3) * t1060;
t1227 = qJD(5) * t1052;
t1225 = qJD(5) * t1057;
t1043 = t1048 ^ 2;
t1062 = qJD(1) ^ 2;
t1219 = t1043 * t1062;
t1201 = t1057 * t1058;
t1198 = qJD(2) - t1033;
t1191 = t1043 * t1276;
t1182 = t778 * t1227;
t1181 = t778 * t1225;
t1180 = qJD(1) * qJD(2) * t1043;
t1178 = qJD(2) * t1213;
t1168 = t1061 * t1219;
t645 = t1053 * t670 + t1058 * t655 + t693 * t1224 - t1226 * t691;
t643 = pkin(14) * t741 + t645;
t1097 = t1054 * t727 + t783 * t1171 + t808 * t1173 - t728 * t1210 - t769 * t1216 + t782 * t1228;
t656 = pkin(4) * t854 + t1097;
t650 = pkin(5) * t696 - pkin(14) * t695 + t656;
t1155 = -t1052 * t643 + t1057 * t650;
t1154 = t1052 * t695 - t1057 * t741;
t1151 = t1057 * t778;
t1150 = t1058 * t837;
t1026 = -pkin(5) * t1058 - pkin(14) * t1053 - pkin(4);
t1149 = pkin(14) * t841 - qJD(6) * t1026 + t1265;
t1148 = -t1054 * t745 - t807 * t1171 - t824 * t1173 - t804 * t1228;
t1144 = pkin(10) * t1178;
t1143 = t1061 * t1180;
t993 = t1034 + (-pkin(3) * t1059 - pkin(4)) * t1049;
t891 = pkin(5) * t1005 - pkin(14) * t1006 + t993;
t1133 = -pkin(14) * t1288 - qJD(6) * t891 + t1281;
t1118 = t1046 * t984 + t1049 * t950;
t918 = t1054 * t931;
t833 = pkin(4) * t1001 - t1059 * t1118 + t918;
t760 = pkin(5) * t889 - pkin(14) * t890 + t833;
t1132 = -pkin(14) * t1252 - qJD(6) * t760 + t1280;
t893 = -pkin(14) * t1216 + t1248;
t1131 = -pkin(5) * t1249 + pkin(14) * t1250 + qJD(6) * t893 - t1247;
t739 = pkin(14) * t938 + t1257;
t1130 = -pkin(5) * t1258 + pkin(14) * t1259 + qJD(6) * t739 - t1260;
t1121 = t1046 * t787 + t1049 * t746;
t1119 = t1046 * t824 + t1049 * t807;
t1116 = t1052 * t650 + t1057 * t643;
t660 = pkin(14) * t837 + t665;
t690 = pkin(4) * t873 - t698;
t673 = pkin(5) * t779 + pkin(14) * t1108 + t690;
t652 = t1052 * t673 + t1057 * t660;
t1115 = t1052 * t660 - t1057 * t673;
t672 = pkin(14) * t846 + t1111;
t800 = t1054 * t804;
t702 = pkin(4) * t888 - t1059 * t1119 + t800;
t805 = t1053 * t847 + t888 * t1058;
t806 = -t1053 * t888 + t1058 * t847;
t684 = pkin(5) * t805 - pkin(14) * t806 + t702;
t1114 = t1052 * t684 + t1057 * t672;
t1113 = -t1052 * t672 + t1057 * t684;
t748 = t1052 * t846 + t1057 * t806;
t747 = t1052 * t806 - t846 * t1057;
t664 = -t1053 * t691 + t1058 * t693;
t1110 = -t1053 * t703 + t1058 * t708;
t1109 = -t1053 * t834 + t1058 * t825;
t1107 = -t1053 * t994 + t1058 * t995;
t1095 = -pkin(13) * t741 + t690 * t837;
t1094 = -t1222 * t778 - t1242;
t1093 = -t1223 * t778 + t1239;
t1084 = t1053 * t679 + t1058 * t662 + t708 * t1224 - t1226 * t703;
t1081 = -qJD(1) * t1144 + t1027;
t1076 = t1047 * t1092;
t1074 = t1079 * t1033;
t659 = -pkin(5) * t837 - t664;
t1067 = -pkin(14) * t696 + (t659 + t664) * t778;
t663 = pkin(4) * t862 - t1059 * t1121 - t1148;
t892 = pkin(5) * t1216 - t1107;
t755 = t1052 * t841 - t1201 * t839;
t754 = -t1052 * t1058 * t839 - t1057 * t841;
t738 = -pkin(5) * t938 - t1109;
t707 = -qJD(5) * t805 - t1053 * t862 + t1058 * t753;
t706 = qJD(5) * t806 + t1053 * t753 + t862 * t1058;
t682 = -pkin(5) * t841 + t1053 * t698 - t1058 * t757;
t675 = -qJD(6) * t747 + t1052 * t752 + t1057 * t707;
t674 = qJD(6) * t748 + t1052 * t707 - t752 * t1057;
t671 = -pkin(5) * t846 - t1110;
t667 = qJD(6) * t731 + t1154;
t653 = pkin(5) * t706 - pkin(14) * t707 + t663;
t648 = -pkin(5) * t752 - t1277;
t647 = pkin(14) * t752 + t1084;
t642 = -qJD(6) * t652 + t1155;
t641 = -qJD(6) * t1115 + t1116;
t1 = [(t1108 * t706 - t695 * t805 - t696 * t806 - t707 * t779) * MDP(26) + (-t1108 * t707 + t695 * t806) * MDP(25) + (-t1108 * t752 + t695 * t846 + t707 * t837 + t741 * t806) * MDP(27) + (-t1084 * t837 - t1108 * t663 - t1111 * t741 - t645 * t846 + t656 * t806 - t665 * t752 + t690 * t707 + t702 * t695) * MDP(31) + ((-qJD(6) * t1114 - t1052 * t647 + t1057 * t653) * t778 + t1113 * t696 + t642 * t805 - t1115 * t706 + t648 * t729 + t671 * t667 + t644 * t747 + t659 * t674) * MDP(37) + (t713 * t839 + t734 * t741 + t700 * t846 + t722 * t752 - t1148 * t873 + t800 * t854 + t1097 * t888 - t698 * t862 + (-t1119 * t854 - t1121 * t873) * t1059) * MDP(23) + (-t1072 * t1122 - t1073 * t1002 + t887 * t908 + t879 * t874 + t875 * t940 + t871 * t882 + (-qJD(1) * t1185 - t830) * t1138) * MDP(17) + (-t1100 * t940 - t939 * t874 - t881 * t908 + t907 * t882) * MDP(12) + (MDP(6) * t1177 - MDP(7) * t1178) * (t1033 + t1234) + (-qJD(1) * t1002 + t1122) * MDP(15) * t1138 + (t1110 * t741 + t1277 * t837 + t646 * t846 + t656 * t805 + t663 * t779 + t664 * t752 + t690 * t706 + t702 * t696) * MDP(30) + (t696 * t805 + t706 * t778) * MDP(36) + (-t667 * t805 - t674 * t778 - t696 * t747 - t706 * t729) * MDP(35) + (t666 * t805 + t675 * t778 + t696 * t748 + t706 * t731) * MDP(34) + (t1074 + (t1051 * t1079 - 0.2e1 * t1191) * qJD(1)) * qJD(2) * MDP(9) + (-t666 * t747 - t667 * t748 - t674 * t731 - t675 * t729) * MDP(33) + (t666 * t748 + t675 * t731) * MDP(32) + (t874 * t940 + t882 * t908) * MDP(11) + (t854 * t888 + t862 * t873) * MDP(22) + (t741 * t888 + t752 * t873 + t839 * t862 + t846 * t854) * MDP(21) + (-t740 * t888 - t753 * t873 - t841 * t862 - t847 * t854) * MDP(20) + (-t696 * t846 - t706 * t837 - t741 * t805 - t752 * t779) * MDP(28) + (t741 * t846 + t752 * t837) * MDP(29) + (t740 * t847 + t753 * t841) * MDP(18) + (-t740 * t846 - t741 * t847 - t752 * t841 - t753 * t839) * MDP(19) + (-(qJD(6) * t1113 + t1052 * t653 + t1057 * t647) * t778 - t1114 * t696 - t641 * t805 - t652 * t706 + t648 * t731 + t671 * t666 + t644 * t748 + t659 * t675) * MDP(38) + (t1071 * t873 + t1188 * t854 + t657 * t888 + t699 * t862 + t700 * t847 + t713 * t841 + t722 * t753 + t734 * t740) * MDP(24) + (-(t1032 - t1144) * t1033 - t1081 * t1051 - 0.2e1 * pkin(1) * t1143) * MDP(10) - 0.2e1 * t1180 * t1283 + 0.2e1 * t1143 * t1284 + (t1063 * t1122 - t759 * t1002 - t887 * t907 + t879 * t1100 + t875 * t939 + t871 * t881 + (qJD(1) * t1128 + t829) * t1138) * MDP(16) + (-t1002 * t874 + t1122 * t882 + (qJD(1) * t940 + t908) * t1138) * MDP(13) + (t1002 * t1100 - t1122 * t881 + (-qJD(1) * t939 + t907) * t1138) * MDP(14); (t1001 * t854 - t1096 * t873) * MDP(22) + (-t1108 * t1259 + t695 * t890) * MDP(25) + (-t1108 * t1252 + t1259 * t837 + t695 * t938 + t741 * t890) * MDP(27) + (t1108 * t1258 - t1259 * t779 - t695 * t889 - t696 * t890) * MDP(26) + (-t1108 * t1260 - t1252 * t665 - t1257 * t741 + t1259 * t690 + t1280 * t837 - t645 * t938 + t656 * t890 + t833 * t695) * MDP(31) + ((-t1052 * t739 + t1057 * t760) * t696 + t642 * t889 + t738 * t667 + t644 * t848 + (t1052 * t1132 - t1057 * t1130) * t778 + t1266 * t729 + t1261 * t659 - t1258 * t1115) * MDP(37) + (-t907 * t978 + t977 * t908 + (t1060 * t874 - t1100 * t1055 + (-t1055 * t908 + t1060 * t907) * qJD(3)) * t1047) * MDP(12) + t1198 * qJD(1) * MDP(6) * t1212 + (t1050 * t1198 - qJD(3) + t1022) * MDP(15) * t1141 - t1198 * MDP(7) * t1179 + (t1097 * t1001 + t700 * t938 + t878 * t741 + t918 * t854 + t1287 * t873 + t1251 * t839 + t1252 * t722 + t1096 * t698 + (-t1118 * t854 + (t1046 * t1299 + t1049 * t1300) * t873) * t1059) * MDP(23) + (t1073 * t1050 - t886 * t908 - t871 * t978 + (t871 * t1230 - pkin(2) * t874 + (pkin(11) * t1092 + t875) * t1055 + (-qJD(2) * t1196 + t830) * t1179) * t1047 + t1302 * t1122) * MDP(17) + (-t1258 * t729 - t1261 * t778 - t667 * t889 - t696 * t848) * MDP(35) + (t1262 * t731 + t666 * t849) * MDP(32) + (t1258 * t731 + t1262 * t778 + t666 * t889 + t696 * t849) * MDP(34) + (-t1261 * t731 - t1262 * t729 - t666 * t848 - t667 * t849) * MDP(33) + (-(t1052 * t760 + t1057 * t739) * t696 - t641 * t889 + t738 * t666 + t644 * t849 + (t1052 * t1130 + t1057 * t1132) * t778 + t1266 * t731 + t1262 * t659 - t1258 * t652) * MDP(38) + (t1258 * t778 + t696 * t889) * MDP(36) + (-t1252 * t779 - t1258 * t837 - t696 * t938 - t741 * t889) * MDP(28) + (t1252 * t837 + t741 * t938) * MDP(29) + (t1253 * t841 + t740 * t941) * MDP(18) + (-t1001 * t740 + t1096 * t841 - t1253 * t873 - t854 * t941) * MDP(20) + (-t1252 * t841 - t1253 * t839 - t740 * t938 - t741 * t941) * MDP(19) + (t1001 * t741 - t1096 * t839 + t1252 * t873 + t854 * t938) * MDP(21) + (-t1050 * t1100 - t1055 * t1076 + t1122 * t977 + (-t1047 * t907 + t1060 * t1233) * t1179) * MDP(14) + (t1050 * t874 + t1060 * t1076 - t1122 * t978 + (-t1047 * t908 + t1055 * t1233) * t1179) * MDP(13) + (-t908 * t978 + (t1055 * t874 + t1230 * t908) * t1047) * MDP(11) + (-t1196 * t1092 + t759 * t1050 - t1047 * pkin(2) * t1100 - t875 * t1214 - t1127 * t1122 + t886 * t907 + (-t977 + t1176) * t871 + (-t829 + (-pkin(11) * t1215 + t1039) * qJD(2)) * t1141) * MDP(16) + (t657 * t1001 - t1096 * t699 + t1186 * t854 + t1251 * t841 + t1253 * t722 + t1296 * t873 + t700 * t941 + t878 * t740) * MDP(24) + (t1062 * t1191 + (qJD(2) * t1079 - t1074) * qJD(1)) * MDP(9) + ((-pkin(10) * t1179 + t1031) * t1033 + pkin(1) * t1168 - t1081) * MDP(10) + (t1109 * t741 + t646 * t938 + t833 * t696 + t656 * t889 + ((-qJD(5) * t834 - t1290) * t1058 + (-qJD(5) * t825 + t1294) * t1053) * t837 + t1260 * t779 + t1258 * t690 + t1252 * t664) * MDP(30) - t1168 * t1284 + t1219 * t1283; (t830 * t1145 - t871 * t908 + t1129) * MDP(16) + (-t1248 * t741 + t993 * t695 + t656 * t1006 + t665 * t856 + t1281 * t837 - t1247 * t1108 + t1250 * t690 + (t1059 * t645 - t1229 * t665) * t1046) * MDP(31) + ((-t1052 * t893 + t1057 * t891) * t696 + t642 * t1005 + t892 * t667 + t644 * t960 + (t1052 * t1133 - t1057 * t1131) * t778 + t1256 * t729 + t1254 * t659 - t1249 * t1115) * MDP(37) + (t1006 * t695 - t1108 * t1250) * MDP(25) + (t1006 * t741 + t1108 * t856 + t1250 * t837 + (-t1059 * t695 - t1108 * t1229) * t1046) * MDP(27) + (-t1005 * t695 - t1006 * t696 + t1108 * t1249 - t1250 * t779) * MDP(26) + (t839 * t857 + t841 * t856 + (-t1054 * t741 + t1059 * t740 + (-t1054 * t841 - t1059 * t839) * qJD(4)) * t1046) * MDP(19) + MDP(15) * t1105 - t907 * t908 * MDP(11) + (t1197 * t854 - t657 * t1049 - t750 * t841 - t722 * t857 + t1297 * t873 + (-pkin(3) * t740 + t1054 * t700 + t1228 * t722 + t699 * t908) * t1046) * MDP(24) + (-t1122 * t907 + t874) * MDP(13) + (t1122 * t908 - t1100) * MDP(14) + (t1122 * t829 - t871 * t907 + t1073) * MDP(17) + (-t907 ^ 2 + t908 ^ 2) * MDP(12) + (-t1005 * t667 - t1249 * t729 - t1254 * t778 - t696 * t960) * MDP(35) + (t1255 * t731 + t666 * t961) * MDP(32) + (t1005 * t666 + t1249 * t731 + t1255 * t778 + t696 * t961) * MDP(34) + (-t1254 * t731 - t1255 * t729 - t666 * t960 - t667 * t961) * MDP(33) + (-(t1052 * t891 + t1057 * t893) * t696 - t641 * t1005 + t892 * t666 + t644 * t961 + (t1052 * t1131 + t1057 * t1133) * t778 + t1256 * t731 + t1255 * t659 - t1249 * t652) * MDP(38) + (t1005 * t696 + t1249 * t778) * MDP(36) + (-t1005 * t741 + t779 * t856 - t1249 * t837 + (t1059 * t696 - t1229 * t779) * t1046) * MDP(28) + (-t1049 * t854 + t1245 * t873) * MDP(22) + (-t841 * t857 + (t1054 * t740 + t1228 * t841) * t1046) * MDP(18) + (-t1049 * t741 - t856 * t873 + (-t1059 * t854 + t1229 * t873 + t839 * t908) * t1046) * MDP(21) + (-t837 * t856 + (-t1059 * t741 + t1229 * t837) * t1046) * MDP(29) + (t1049 * t740 + t857 * t873 + (-t1054 * t854 - t1228 * t873 - t841 * t908) * t1046) * MDP(20) + (-t1282 * t854 - t1097 * t1049 - t750 * t839 - t722 * t856 + (t1210 * t793 + t1289) * t873 + (t722 * t1229 - pkin(3) * t741 - t698 * t908 + (t858 * t873 - t700) * t1059) * t1046) * MDP(23) + (t1107 * t741 + t993 * t696 + t656 * t1005 - t664 * t856 + ((-qJD(5) * t994 - t1291) * t1058 + (-qJD(5) * t995 + t1295) * t1053) * t837 + t1247 * t779 + t1249 * t690 + (-t1059 * t646 + t1229 * t664) * t1046) * MDP(30); -t839 ^ 2 * MDP(19) + (-t839 * t873 + t740) * MDP(20) - t741 * MDP(21) - t854 * MDP(22) + (-t699 * t873 - t1097) * MDP(23) + (-t698 * t873 + t722 * t839 - t657) * MDP(24) + (t1053 * t695 - t1108 * t1150) * MDP(25) + ((t695 - t1268) * t1058 + (-t696 + t1267) * t1053) * MDP(26) + (t1053 * t741 + t1150 * t837) * MDP(27) + (-t1053 * t837 ^ 2 + t1058 * t741) * MDP(28) + (-pkin(4) * t696 - t699 * t779 + (-t656 + (-t757 - t1271) * t837) * t1058 + (t698 * t837 + t1095) * t1053) * MDP(30) + (-pkin(4) * t695 + t1265 * t837 + t699 * t1108 + (t1271 * t837 + t656) * t1053 + t1095 * t1058) * MDP(31) + (t1053 * t1057 * t666 + (qJD(5) * t1201 - t1053 * t1223 - t755) * t731) * MDP(32) + (t729 * t755 + t731 * t754 + (-t1052 * t731 - t1057 * t729) * t1224 + (-t1243 - t1057 * t667 + (t1052 * t729 - t1057 * t731) * qJD(6)) * t1053) * MDP(33) + (-t755 * t778 + (-t666 + t1181) * t1058 + (t731 * t837 + t1093) * t1053) * MDP(34) + (t754 * t778 + (t667 - t1182) * t1058 + (-t729 * t837 + t1094) * t1053) * MDP(35) + (t1053 * t778 * t837 - t1058 * t696) * MDP(36) + (t1026 * t1239 - t659 * t754 - t682 * t729 + (t1052 * t1149 + t1057 * t1236) * t778 + (t659 * t1227 - t642 + (qJD(5) * t729 + t1094) * pkin(13)) * t1058 + (t659 * t1222 + t644 * t1052 - t837 * t1115 + (t667 + t1182) * pkin(13)) * t1053) * MDP(37) + (-t1026 * t1242 - t659 * t755 - t682 * t731 + (-t1052 * t1236 + t1057 * t1149) * t778 + (t659 * t1225 + t641 + (qJD(5) * t731 - t1093) * pkin(13)) * t1058 + (-t659 * t1223 + t644 * t1057 - t837 * t652 + (t666 + t1181) * pkin(13)) * t1053) * MDP(38) + (MDP(18) * t839 + MDP(19) * t841 - MDP(21) * t873 - MDP(23) * t722 + MDP(27) * t1108 + MDP(28) * t779 - MDP(29) * t837 - MDP(30) * t664 + MDP(31) * t665) * t841; -t779 ^ 2 * MDP(26) + (t695 + t1268) * MDP(27) + (-t696 - t1267) * MDP(28) + t741 * MDP(29) + (t665 * t837 + t646) * MDP(30) + (t664 * t837 + t690 * t779 - t645) * MDP(31) + (t1151 * t731 + t1243) * MDP(32) + ((t666 - t1270) * t1057 + (-t667 - t1269) * t1052) * MDP(33) + (t1151 * t778 + t1242) * MDP(34) + (-t1052 * t778 ^ 2 + t1239) * MDP(35) + (-pkin(5) * t667 + t1067 * t1052 - t1057 * t1278 - t665 * t729) * MDP(37) + (-pkin(5) * t666 + t1052 * t1278 + t1067 * t1057 - t665 * t731) * MDP(38) - (MDP(25) * t779 - MDP(26) * t1108 - MDP(30) * t690 - MDP(34) * t731 + MDP(35) * t729 - MDP(36) * t778 + MDP(37) * t1115 + MDP(38) * t652) * t1108; t731 * t729 * MDP(32) + (-t729 ^ 2 + t731 ^ 2) * MDP(33) + (t1190 + t1270) * MDP(34) + (-t1154 + t1269) * MDP(35) + t696 * MDP(36) + (t652 * t778 - t659 * t731 + t1155) * MDP(37) + (-t1115 * t778 + t659 * t729 - t1116) * MDP(38) + (MDP(34) * t1241 - MDP(35) * t731 - MDP(37) * t652 + MDP(38) * t1115) * qJD(6);];
tauc  = t1;
