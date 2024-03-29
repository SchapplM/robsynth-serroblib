% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RRPPR2
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RRPPR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:06:18
% EndTime: 2022-01-20 10:06:20
% DurationCPUTime: 2.35s
% Computational Cost: add. (8661->175), mult. (12134->252), div. (0->0), fcn. (7317->10), ass. (0->142)
t1161 = qJD(1) + qJD(2);
t1157 = t1161 ^ 2;
t1166 = sin(pkin(8));
t1158 = qJDD(1) + qJDD(2);
t1168 = cos(pkin(8));
t1186 = t1168 * t1158;
t1131 = t1166 * t1157 - t1186;
t1170 = sin(qJ(2));
t1173 = cos(qJ(2));
t1188 = t1166 * t1158;
t1180 = -t1168 * t1157 - t1188;
t1109 = t1173 * t1131 - t1170 * t1180;
t1171 = sin(qJ(1));
t1174 = cos(qJ(1));
t1203 = t1170 * t1131 + t1173 * t1180;
t1207 = t1171 * t1109 + t1174 * t1203;
t1206 = t1174 * t1109 - t1171 * t1203;
t1137 = t1170 * t1157 - t1173 * t1158;
t1179 = -t1173 * t1157 - t1170 * t1158;
t1202 = t1171 * t1137 + t1174 * t1179;
t1201 = t1174 * t1137 - t1171 * t1179;
t1165 = sin(pkin(9));
t1159 = t1165 ^ 2;
t1167 = cos(pkin(9));
t1160 = t1167 ^ 2;
t1185 = t1159 + t1160;
t1134 = t1185 * t1157;
t1187 = t1167 * t1161;
t1145 = -qJD(5) + t1187;
t1196 = t1145 ^ 2;
t1195 = 2 * qJD(4);
t1169 = sin(qJ(5));
t1194 = t1158 * t1169;
t1172 = cos(qJ(5));
t1193 = t1158 * t1172;
t1192 = t1159 * t1157;
t1191 = t1161 * t1169;
t1190 = t1161 * t1172;
t1189 = t1165 * t1158;
t1151 = t1167 * t1158;
t1149 = t1171 * g(1) - t1174 * g(2);
t1139 = qJDD(1) * pkin(1) + t1149;
t1150 = -t1174 * g(1) - t1171 * g(2);
t1175 = qJD(1) ^ 2;
t1140 = -t1175 * pkin(1) + t1150;
t1117 = t1170 * t1139 + t1173 * t1140;
t1111 = -t1157 * pkin(2) + t1117;
t1116 = t1173 * t1139 - t1170 * t1140;
t1178 = t1158 * pkin(2) + t1116;
t1087 = t1168 * t1111 + t1166 * t1178;
t1184 = t1151 - qJDD(5);
t1183 = t1145 * t1191;
t1081 = -t1157 * pkin(3) + t1158 * qJ(4) + t1087;
t1164 = -g(3) + qJDD(3);
t1073 = t1167 * t1081 + t1165 * t1164 + t1187 * t1195;
t1086 = -t1166 * t1111 + t1168 * t1178;
t1182 = t1169 * t1172 * t1192;
t1181 = -pkin(4) * t1167 - pkin(7) * t1165;
t1080 = -t1158 * pkin(3) - t1157 * qJ(4) + qJDD(4) - t1086;
t1163 = t1172 ^ 2;
t1162 = t1169 ^ 2;
t1153 = t1167 * t1164;
t1144 = t1167 * t1157 * t1165;
t1143 = -t1171 * qJDD(1) - t1174 * t1175;
t1142 = t1174 * qJDD(1) - t1171 * t1175;
t1141 = t1165 * qJD(5) * t1191;
t1126 = t1185 * t1158;
t1125 = t1181 * t1161;
t1124 = t1167 * t1134;
t1123 = t1165 * t1134;
t1122 = (t1162 + t1163) * t1192;
t1121 = -t1162 * t1192 - t1196;
t1120 = -t1182 - t1184;
t1119 = -t1182 + t1184;
t1118 = -t1163 * t1192 - t1196;
t1115 = -t1168 * t1124 - t1166 * t1151;
t1114 = t1168 * t1123 + t1165 * t1188;
t1113 = -t1166 * t1124 + t1167 * t1186;
t1112 = t1166 * t1123 - t1165 * t1186;
t1106 = t1168 * t1126 - t1166 * t1134;
t1105 = t1166 * t1126 + t1168 * t1134;
t1101 = t1141 + (t1183 - t1193) * t1165;
t1100 = -t1141 + (t1183 + t1193) * t1165;
t1099 = (-t1194 + (-qJD(5) - t1145) * t1190) * t1165;
t1098 = (t1194 + (qJD(5) - t1145) * t1190) * t1165;
t1097 = -t1169 * t1120 + t1172 * t1121;
t1096 = t1172 * t1120 + t1169 * t1121;
t1095 = -t1169 * t1118 + t1172 * t1119;
t1094 = t1172 * t1118 + t1169 * t1119;
t1093 = -t1170 * t1116 + t1173 * t1117;
t1092 = t1173 * t1116 + t1170 * t1117;
t1091 = -t1170 * t1113 + t1173 * t1115;
t1090 = -t1170 * t1112 + t1173 * t1114;
t1089 = t1173 * t1113 + t1170 * t1115;
t1088 = t1173 * t1112 + t1170 * t1114;
t1085 = -t1170 * t1105 + t1173 * t1106;
t1084 = t1173 * t1105 + t1170 * t1106;
t1083 = t1172 * t1099 - t1169 * t1101;
t1082 = t1169 * t1099 + t1172 * t1101;
t1078 = t1167 * t1097 + t1165 * t1098;
t1077 = t1165 * t1097 - t1167 * t1098;
t1076 = t1167 * t1095 + t1165 * t1100;
t1075 = t1165 * t1095 - t1167 * t1100;
t1074 = t1181 * t1158 + t1080;
t1072 = t1153 + (-0.2e1 * qJD(4) * t1161 - t1081) * t1165;
t1071 = t1167 * t1083 - t1165 * t1122;
t1070 = t1165 * t1083 + t1167 * t1122;
t1069 = t1125 * t1187 + t1073;
t1068 = -t1153 + (t1081 + (t1195 + t1125) * t1161) * t1165;
t1067 = t1168 * t1078 + t1166 * t1096;
t1066 = t1166 * t1078 - t1168 * t1096;
t1065 = -t1166 * t1086 + t1168 * t1087;
t1064 = t1168 * t1086 + t1166 * t1087;
t1063 = t1168 * t1076 + t1166 * t1094;
t1062 = t1166 * t1076 - t1168 * t1094;
t1061 = t1168 * t1071 + t1166 * t1082;
t1060 = t1166 * t1071 - t1168 * t1082;
t1059 = -t1165 * t1072 + t1167 * t1073;
t1058 = t1167 * t1072 + t1165 * t1073;
t1057 = t1172 * t1069 + t1169 * t1074;
t1056 = -t1169 * t1069 + t1172 * t1074;
t1055 = -t1170 * t1066 + t1173 * t1067;
t1054 = t1173 * t1066 + t1170 * t1067;
t1053 = t1168 * t1059 + t1166 * t1080;
t1052 = t1166 * t1059 - t1168 * t1080;
t1051 = -t1170 * t1064 + t1173 * t1065;
t1050 = t1173 * t1064 + t1170 * t1065;
t1049 = -t1170 * t1062 + t1173 * t1063;
t1048 = t1173 * t1062 + t1170 * t1063;
t1047 = -t1170 * t1060 + t1173 * t1061;
t1046 = t1173 * t1060 + t1170 * t1061;
t1045 = -t1169 * t1056 + t1172 * t1057;
t1044 = t1172 * t1056 + t1169 * t1057;
t1043 = t1167 * t1045 + t1165 * t1068;
t1042 = t1165 * t1045 - t1167 * t1068;
t1041 = -t1170 * t1052 + t1173 * t1053;
t1040 = t1173 * t1052 + t1170 * t1053;
t1039 = t1168 * t1043 + t1166 * t1044;
t1038 = t1166 * t1043 - t1168 * t1044;
t1037 = -t1170 * t1038 + t1173 * t1039;
t1036 = t1173 * t1038 + t1170 * t1039;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1143, -t1142, 0, -t1171 * t1149 + t1174 * t1150, 0, 0, 0, 0, 0, 0, t1202, t1201, 0, -t1171 * t1092 + t1174 * t1093, 0, 0, 0, 0, 0, 0, t1207, t1206, 0, -t1171 * t1050 + t1174 * t1051, 0, 0, 0, 0, 0, 0, -t1171 * t1089 + t1174 * t1091, -t1171 * t1088 + t1174 * t1090, -t1171 * t1084 + t1174 * t1085, -t1171 * t1040 + t1174 * t1041, 0, 0, 0, 0, 0, 0, -t1171 * t1054 + t1174 * t1055, -t1171 * t1048 + t1174 * t1049, -t1171 * t1046 + t1174 * t1047, -t1171 * t1036 + t1174 * t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1142, t1143, 0, t1174 * t1149 + t1171 * t1150, 0, 0, 0, 0, 0, 0, -t1201, t1202, 0, t1174 * t1092 + t1171 * t1093, 0, 0, 0, 0, 0, 0, -t1206, t1207, 0, t1174 * t1050 + t1171 * t1051, 0, 0, 0, 0, 0, 0, t1174 * t1089 + t1171 * t1091, t1174 * t1088 + t1171 * t1090, t1174 * t1084 + t1171 * t1085, t1174 * t1040 + t1171 * t1041, 0, 0, 0, 0, 0, 0, t1174 * t1054 + t1171 * t1055, t1174 * t1048 + t1171 * t1049, t1174 * t1046 + t1171 * t1047, t1174 * t1036 + t1171 * t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1164, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1058, 0, 0, 0, 0, 0, 0, t1077, t1075, t1070, t1042; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1175, -qJDD(1), 0, t1150, 0, 0, 0, 0, 0, 0, t1179, t1137, 0, t1093, 0, 0, 0, 0, 0, 0, t1203, t1109, 0, t1051, 0, 0, 0, 0, 0, 0, t1091, t1090, t1085, t1041, 0, 0, 0, 0, 0, 0, t1055, t1049, t1047, t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1175, 0, t1149, 0, 0, 0, 0, 0, 0, -t1137, t1179, 0, t1092, 0, 0, 0, 0, 0, 0, -t1109, t1203, 0, t1050, 0, 0, 0, 0, 0, 0, t1089, t1088, t1084, t1040, 0, 0, 0, 0, 0, 0, t1054, t1048, t1046, t1036; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1164, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1058, 0, 0, 0, 0, 0, 0, t1077, t1075, t1070, t1042; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1157, -t1158, 0, t1117, 0, 0, 0, 0, 0, 0, t1180, t1131, 0, t1065, 0, 0, 0, 0, 0, 0, t1115, t1114, t1106, t1053, 0, 0, 0, 0, 0, 0, t1067, t1063, t1061, t1039; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1158, -t1157, 0, t1116, 0, 0, 0, 0, 0, 0, -t1131, t1180, 0, t1064, 0, 0, 0, 0, 0, 0, t1113, t1112, t1105, t1052, 0, 0, 0, 0, 0, 0, t1066, t1062, t1060, t1038; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1164, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1058, 0, 0, 0, 0, 0, 0, t1077, t1075, t1070, t1042; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1157, -t1158, 0, t1087, 0, 0, 0, 0, 0, 0, -t1124, t1123, t1126, t1059, 0, 0, 0, 0, 0, 0, t1078, t1076, t1071, t1043; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1158, -t1157, 0, t1086, 0, 0, 0, 0, 0, 0, t1151, -t1189, t1134, -t1080, 0, 0, 0, 0, 0, 0, -t1096, -t1094, -t1082, -t1044; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1164, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1058, 0, 0, 0, 0, 0, 0, t1077, t1075, t1070, t1042; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1160 * t1157, t1144, t1151, t1073, 0, 0, 0, 0, 0, 0, t1097, t1095, t1083, t1045; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1144, -t1192, -t1189, t1072, 0, 0, 0, 0, 0, 0, -t1098, -t1100, t1122, -t1068; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1151, t1189, -t1134, t1080, 0, 0, 0, 0, 0, 0, t1096, t1094, t1082, t1044; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1121, t1119, t1099, t1057; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1120, t1118, t1101, t1056; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1098, t1100, -t1122, t1068;];
f_new_reg = t1;
