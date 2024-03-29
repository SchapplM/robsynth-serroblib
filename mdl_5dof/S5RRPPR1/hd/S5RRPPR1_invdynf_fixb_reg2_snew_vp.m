% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRPPR1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:52:11
% EndTime: 2022-01-20 09:52:13
% DurationCPUTime: 2.85s
% Computational Cost: add. (9585->168), mult. (13584->248), div. (0->0), fcn. (8597->10), ass. (0->138)
t1184 = qJD(1) + qJD(2);
t1180 = t1184 ^ 2;
t1187 = sin(pkin(8));
t1181 = qJDD(1) + qJDD(2);
t1189 = cos(pkin(8));
t1208 = t1189 * t1181;
t1155 = t1187 * t1180 - t1208;
t1191 = sin(qJ(2));
t1194 = cos(qJ(2));
t1210 = t1187 * t1181;
t1205 = -t1189 * t1180 - t1210;
t1133 = t1194 * t1155 - t1191 * t1205;
t1192 = sin(qJ(1));
t1195 = cos(qJ(1));
t1226 = t1191 * t1155 + t1194 * t1205;
t1232 = t1192 * t1133 + t1195 * t1226;
t1231 = t1195 * t1133 - t1192 * t1226;
t1230 = 2 * qJD(5);
t1186 = sin(pkin(9));
t1188 = cos(pkin(9));
t1190 = sin(qJ(5));
t1193 = cos(qJ(5));
t1227 = -t1186 * t1190 + t1188 * t1193;
t1161 = t1191 * t1180 - t1194 * t1181;
t1204 = -t1194 * t1180 - t1191 * t1181;
t1225 = t1192 * t1161 + t1195 * t1204;
t1224 = t1195 * t1161 - t1192 * t1204;
t1203 = t1186 * t1193 + t1188 * t1190;
t1145 = t1203 * t1181;
t1182 = t1186 ^ 2;
t1183 = t1188 ^ 2;
t1206 = t1182 + t1183;
t1158 = t1206 * t1180;
t1146 = t1227 * t1184;
t1219 = t1146 ^ 2;
t1148 = t1203 * t1184;
t1218 = t1148 ^ 2;
t1217 = qJD(4) * t1184;
t1216 = t1148 * t1146;
t1215 = t1180 * t1188;
t1214 = t1182 * t1180;
t1213 = t1183 * t1180;
t1212 = t1186 * t1181;
t1175 = t1188 * t1181;
t1170 = t1192 * g(1) - t1195 * g(2);
t1201 = qJDD(1) * pkin(1) + t1170;
t1171 = -t1195 * g(1) - t1192 * g(2);
t1197 = qJD(1) ^ 2;
t1202 = -t1197 * pkin(1) + t1171;
t1141 = t1191 * t1201 + t1194 * t1202;
t1135 = -t1180 * pkin(2) + t1141;
t1140 = -t1191 * t1202 + t1194 * t1201;
t1200 = t1181 * pkin(2) + t1140;
t1110 = t1189 * t1135 + t1187 * t1200;
t1185 = -g(3) + qJDD(3);
t1207 = t1188 * t1185 - 0.2e1 * t1186 * t1217;
t1106 = -t1180 * pkin(3) + t1181 * qJ(4) + t1110;
t1099 = t1186 * t1185 + (t1106 + 0.2e1 * t1217) * t1188;
t1109 = -t1187 * t1135 + t1189 * t1200;
t1119 = t1227 * t1181;
t1103 = -t1181 * pkin(3) - t1180 * qJ(4) + qJDD(4) - t1109;
t1196 = qJD(5) ^ 2;
t1165 = t1186 * t1215;
t1164 = -t1192 * qJDD(1) - t1195 * t1197;
t1163 = t1195 * qJDD(1) - t1192 * t1197;
t1151 = t1206 * t1181;
t1150 = t1188 * t1158;
t1149 = t1186 * t1158;
t1142 = -t1196 - t1218;
t1139 = -t1189 * t1150 - t1187 * t1175;
t1138 = t1189 * t1149 + t1186 * t1210;
t1137 = -t1187 * t1150 + t1188 * t1208;
t1136 = t1187 * t1149 - t1186 * t1208;
t1130 = t1189 * t1151 - t1187 * t1158;
t1129 = t1187 * t1151 + t1189 * t1158;
t1125 = t1146 * t1230 + t1145;
t1124 = t1148 * t1230 - t1119;
t1123 = -qJDD(5) + t1216;
t1122 = qJDD(5) + t1216;
t1121 = -t1196 - t1219;
t1120 = -t1218 - t1219;
t1118 = -t1191 * t1140 + t1194 * t1141;
t1117 = t1194 * t1140 + t1191 * t1141;
t1116 = -t1191 * t1137 + t1194 * t1139;
t1115 = -t1191 * t1136 + t1194 * t1138;
t1114 = t1194 * t1137 + t1191 * t1139;
t1113 = t1194 * t1136 + t1191 * t1138;
t1112 = t1193 * t1123 - t1190 * t1142;
t1111 = t1190 * t1123 + t1193 * t1142;
t1108 = -t1191 * t1129 + t1194 * t1130;
t1107 = t1194 * t1129 + t1191 * t1130;
t1105 = t1193 * t1119 + t1190 * t1145;
t1104 = t1190 * t1119 - t1193 * t1145;
t1102 = t1193 * t1121 - t1190 * t1122;
t1101 = t1190 * t1121 + t1193 * t1122;
t1098 = -t1186 * t1106 + t1207;
t1097 = -pkin(4) * t1175 + t1103 + (-t1213 - t1214) * pkin(7);
t1096 = -pkin(4) * t1213 + pkin(7) * t1175 + t1099;
t1095 = (pkin(4) * t1215 - pkin(7) * t1181 - t1106) * t1186 + t1207;
t1094 = -t1186 * t1111 + t1188 * t1112;
t1093 = t1188 * t1111 + t1186 * t1112;
t1092 = -t1187 * t1109 + t1189 * t1110;
t1091 = t1189 * t1109 + t1187 * t1110;
t1090 = -t1186 * t1104 + t1188 * t1105;
t1089 = t1188 * t1104 + t1186 * t1105;
t1088 = -t1186 * t1101 + t1188 * t1102;
t1087 = t1188 * t1101 + t1186 * t1102;
t1086 = t1189 * t1094 + t1187 * t1125;
t1085 = t1187 * t1094 - t1189 * t1125;
t1084 = -t1186 * t1098 + t1188 * t1099;
t1083 = t1188 * t1098 + t1186 * t1099;
t1082 = t1189 * t1088 + t1187 * t1124;
t1081 = t1187 * t1088 - t1189 * t1124;
t1080 = t1189 * t1090 + t1187 * t1120;
t1079 = t1187 * t1090 - t1189 * t1120;
t1078 = t1190 * t1095 + t1193 * t1096;
t1077 = t1193 * t1095 - t1190 * t1096;
t1076 = t1189 * t1084 + t1187 * t1103;
t1075 = t1187 * t1084 - t1189 * t1103;
t1074 = -t1191 * t1091 + t1194 * t1092;
t1073 = t1194 * t1091 + t1191 * t1092;
t1072 = -t1191 * t1085 + t1194 * t1086;
t1071 = t1194 * t1085 + t1191 * t1086;
t1070 = -t1191 * t1081 + t1194 * t1082;
t1069 = t1194 * t1081 + t1191 * t1082;
t1068 = -t1191 * t1079 + t1194 * t1080;
t1067 = t1194 * t1079 + t1191 * t1080;
t1066 = -t1190 * t1077 + t1193 * t1078;
t1065 = t1193 * t1077 + t1190 * t1078;
t1064 = -t1191 * t1075 + t1194 * t1076;
t1063 = t1194 * t1075 + t1191 * t1076;
t1062 = -t1186 * t1065 + t1188 * t1066;
t1061 = t1188 * t1065 + t1186 * t1066;
t1060 = t1189 * t1062 + t1187 * t1097;
t1059 = t1187 * t1062 - t1189 * t1097;
t1058 = -t1191 * t1059 + t1194 * t1060;
t1057 = t1194 * t1059 + t1191 * t1060;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1164, -t1163, 0, -t1192 * t1170 + t1195 * t1171, 0, 0, 0, 0, 0, 0, t1225, t1224, 0, -t1192 * t1117 + t1195 * t1118, 0, 0, 0, 0, 0, 0, t1232, t1231, 0, -t1192 * t1073 + t1195 * t1074, 0, 0, 0, 0, 0, 0, -t1192 * t1114 + t1195 * t1116, -t1192 * t1113 + t1195 * t1115, -t1192 * t1107 + t1195 * t1108, -t1192 * t1063 + t1195 * t1064, 0, 0, 0, 0, 0, 0, -t1192 * t1069 + t1195 * t1070, -t1192 * t1071 + t1195 * t1072, -t1192 * t1067 + t1195 * t1068, -t1192 * t1057 + t1195 * t1058; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1163, t1164, 0, t1195 * t1170 + t1192 * t1171, 0, 0, 0, 0, 0, 0, -t1224, t1225, 0, t1195 * t1117 + t1192 * t1118, 0, 0, 0, 0, 0, 0, -t1231, t1232, 0, t1195 * t1073 + t1192 * t1074, 0, 0, 0, 0, 0, 0, t1195 * t1114 + t1192 * t1116, t1195 * t1113 + t1192 * t1115, t1195 * t1107 + t1192 * t1108, t1195 * t1063 + t1192 * t1064, 0, 0, 0, 0, 0, 0, t1195 * t1069 + t1192 * t1070, t1195 * t1071 + t1192 * t1072, t1195 * t1067 + t1192 * t1068, t1195 * t1057 + t1192 * t1058; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1083, 0, 0, 0, 0, 0, 0, t1087, t1093, t1089, t1061; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1197, -qJDD(1), 0, t1171, 0, 0, 0, 0, 0, 0, t1204, t1161, 0, t1118, 0, 0, 0, 0, 0, 0, t1226, t1133, 0, t1074, 0, 0, 0, 0, 0, 0, t1116, t1115, t1108, t1064, 0, 0, 0, 0, 0, 0, t1070, t1072, t1068, t1058; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1197, 0, t1170, 0, 0, 0, 0, 0, 0, -t1161, t1204, 0, t1117, 0, 0, 0, 0, 0, 0, -t1133, t1226, 0, t1073, 0, 0, 0, 0, 0, 0, t1114, t1113, t1107, t1063, 0, 0, 0, 0, 0, 0, t1069, t1071, t1067, t1057; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1083, 0, 0, 0, 0, 0, 0, t1087, t1093, t1089, t1061; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1180, -t1181, 0, t1141, 0, 0, 0, 0, 0, 0, t1205, t1155, 0, t1092, 0, 0, 0, 0, 0, 0, t1139, t1138, t1130, t1076, 0, 0, 0, 0, 0, 0, t1082, t1086, t1080, t1060; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1181, -t1180, 0, t1140, 0, 0, 0, 0, 0, 0, -t1155, t1205, 0, t1091, 0, 0, 0, 0, 0, 0, t1137, t1136, t1129, t1075, 0, 0, 0, 0, 0, 0, t1081, t1085, t1079, t1059; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1083, 0, 0, 0, 0, 0, 0, t1087, t1093, t1089, t1061; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1180, -t1181, 0, t1110, 0, 0, 0, 0, 0, 0, -t1150, t1149, t1151, t1084, 0, 0, 0, 0, 0, 0, t1088, t1094, t1090, t1062; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1181, -t1180, 0, t1109, 0, 0, 0, 0, 0, 0, t1175, -t1212, t1158, -t1103, 0, 0, 0, 0, 0, 0, -t1124, -t1125, -t1120, -t1097; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1083, 0, 0, 0, 0, 0, 0, t1087, t1093, t1089, t1061; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1213, t1165, t1175, t1099, 0, 0, 0, 0, 0, 0, t1102, t1112, t1105, t1066; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1165, -t1214, -t1212, t1098, 0, 0, 0, 0, 0, 0, t1101, t1111, t1104, t1065; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1175, t1212, -t1158, t1103, 0, 0, 0, 0, 0, 0, t1124, t1125, t1120, t1097; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1121, t1123, t1119, t1078; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1122, t1142, -t1145, t1077; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1124, t1125, t1120, t1097;];
f_new_reg = t1;
