% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5RPPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5RPPRP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynf_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:13:02
% EndTime: 2022-01-23 09:13:03
% DurationCPUTime: 1.73s
% Computational Cost: add. (4179->168), mult. (9290->230), div. (0->0), fcn. (5585->8), ass. (0->131)
t1184 = -2 * qJD(5);
t1155 = qJD(1) ^ 2;
t1147 = sin(pkin(8));
t1141 = t1147 ^ 2;
t1149 = cos(pkin(8));
t1142 = t1149 ^ 2;
t1174 = t1141 + t1142;
t1124 = t1174 * t1155;
t1176 = t1149 * qJD(1);
t1132 = -qJD(4) + t1176;
t1183 = t1132 ^ 2;
t1182 = 2 * qJD(3);
t1181 = -g(3) + qJDD(2);
t1180 = qJD(1) * t1147;
t1151 = sin(qJ(4));
t1179 = qJD(1) * t1151;
t1153 = cos(qJ(4));
t1178 = qJD(1) * t1153;
t1177 = t1141 * t1155;
t1161 = -pkin(3) * t1149 - pkin(6) * t1147;
t1119 = t1161 * qJD(1);
t1175 = t1182 + t1119;
t1152 = sin(qJ(1));
t1154 = cos(qJ(1));
t1128 = -t1154 * g(1) - t1152 * g(2);
t1120 = -t1155 * pkin(1) + t1128;
t1148 = sin(pkin(7));
t1150 = cos(pkin(7));
t1127 = t1152 * g(1) - t1154 * g(2);
t1159 = qJDD(1) * pkin(1) + t1127;
t1099 = t1150 * t1120 + t1148 * t1159;
t1093 = -t1155 * pkin(2) + qJDD(1) * qJ(3) + t1099;
t1083 = t1149 * t1093 + t1147 * t1181 + t1176 * t1182;
t1077 = t1119 * t1176 + t1083;
t1098 = -t1148 * t1120 + t1150 * t1159;
t1092 = -qJDD(1) * pkin(2) - t1155 * qJ(3) + qJDD(3) - t1098;
t1086 = t1161 * qJDD(1) + t1092;
t1063 = t1153 * t1077 + t1151 * t1086;
t1173 = qJDD(1) * t1151;
t1172 = qJDD(1) * t1153;
t1171 = t1147 * qJDD(1);
t1170 = t1148 * qJDD(1);
t1138 = t1149 * qJDD(1);
t1169 = t1150 * qJDD(1);
t1168 = -t1138 + qJDD(4);
t1167 = t1132 * t1179;
t1166 = t1147 * t1178;
t1165 = t1153 * t1177;
t1144 = t1151 ^ 2;
t1164 = t1144 * t1177;
t1122 = -t1150 * t1155 - t1170;
t1123 = -t1148 * t1155 + t1169;
t1163 = t1154 * t1122 - t1152 * t1123;
t1162 = t1151 * t1165;
t1160 = t1152 * t1122 + t1154 * t1123;
t1158 = -qJD(4) * t1178 - t1173;
t1145 = t1153 ^ 2;
t1137 = t1149 * t1181;
t1131 = t1149 * t1155 * t1147;
t1129 = t1147 * qJD(4) * t1179;
t1126 = -t1152 * qJDD(1) - t1154 * t1155;
t1125 = t1154 * qJDD(1) - t1152 * t1155;
t1121 = t1174 * qJDD(1);
t1118 = t1149 * t1124;
t1117 = t1147 * t1124;
t1112 = (t1144 + t1145) * t1177;
t1111 = -t1164 - t1183;
t1110 = -t1162 + t1168;
t1109 = -t1162 - t1168;
t1108 = -t1132 * pkin(4) - qJ(5) * t1166;
t1106 = -t1145 * t1177 - t1183;
t1105 = -t1150 * t1118 - t1148 * t1138;
t1104 = t1150 * t1117 + t1147 * t1170;
t1103 = -t1148 * t1118 + t1149 * t1169;
t1102 = t1148 * t1117 - t1147 * t1169;
t1101 = t1150 * t1121 - t1148 * t1124;
t1100 = t1148 * t1121 + t1150 * t1124;
t1097 = t1129 + (t1167 - t1172) * t1147;
t1096 = -t1129 + (t1167 + t1172) * t1147;
t1095 = (-t1173 + (-qJD(4) - t1132) * t1178) * t1147;
t1094 = (t1173 + (qJD(4) - t1132) * t1178) * t1147;
t1091 = -t1151 * t1110 + t1153 * t1111;
t1090 = t1153 * t1110 + t1151 * t1111;
t1088 = -t1151 * t1106 + t1153 * t1109;
t1087 = t1153 * t1106 + t1151 * t1109;
t1085 = t1153 * t1086;
t1082 = -0.2e1 * qJD(3) * t1180 - t1147 * t1093 + t1137;
t1081 = t1153 * t1095 - t1151 * t1097;
t1080 = t1151 * t1095 + t1153 * t1097;
t1079 = -t1148 * t1098 + t1150 * t1099;
t1078 = t1150 * t1098 + t1148 * t1099;
t1076 = -t1137 + (t1175 * qJD(1) + t1093) * t1147;
t1074 = t1149 * t1091 + t1147 * t1094;
t1073 = t1147 * t1091 - t1149 * t1094;
t1072 = t1149 * t1088 + t1147 * t1096;
t1071 = t1147 * t1088 - t1149 * t1096;
t1070 = t1149 * t1081 - t1147 * t1112;
t1069 = t1147 * t1081 + t1149 * t1112;
t1068 = -t1147 * t1082 + t1149 * t1083;
t1067 = t1149 * t1082 + t1147 * t1083;
t1066 = -qJ(5) * t1164 + qJDD(5) - t1137 + (-t1158 * pkin(4) + t1093 + (t1108 * t1153 + t1175) * qJD(1)) * t1147;
t1065 = t1150 * t1074 + t1148 * t1090;
t1064 = t1148 * t1074 - t1150 * t1090;
t1062 = -t1151 * t1077 + t1085;
t1061 = t1150 * t1072 + t1148 * t1087;
t1060 = t1148 * t1072 - t1150 * t1087;
t1059 = t1150 * t1068 + t1148 * t1092;
t1058 = t1148 * t1068 - t1150 * t1092;
t1057 = t1150 * t1070 + t1148 * t1080;
t1056 = t1148 * t1070 - t1150 * t1080;
t1055 = -pkin(4) * t1164 + t1132 * t1108 + (t1158 * qJ(5) + t1179 * t1184) * t1147 + t1063;
t1054 = t1085 - (t1153 * t1171 - t1129) * qJ(5) + t1168 * pkin(4) + t1166 * t1184 + (t1132 * qJ(5) * t1180 - pkin(4) * t1165 - t1077) * t1151;
t1053 = -t1152 * t1064 + t1154 * t1065;
t1052 = t1154 * t1064 + t1152 * t1065;
t1051 = -t1151 * t1062 + t1153 * t1063;
t1050 = t1153 * t1062 + t1151 * t1063;
t1049 = -t1152 * t1060 + t1154 * t1061;
t1048 = t1154 * t1060 + t1152 * t1061;
t1047 = t1149 * t1051 + t1147 * t1076;
t1046 = t1147 * t1051 - t1149 * t1076;
t1045 = -t1152 * t1056 + t1154 * t1057;
t1044 = t1154 * t1056 + t1152 * t1057;
t1043 = -t1151 * t1054 + t1153 * t1055;
t1042 = t1153 * t1054 + t1151 * t1055;
t1041 = t1149 * t1043 + t1147 * t1066;
t1040 = t1147 * t1043 - t1149 * t1066;
t1039 = t1150 * t1047 + t1148 * t1050;
t1038 = t1148 * t1047 - t1150 * t1050;
t1037 = t1150 * t1041 + t1148 * t1042;
t1036 = t1148 * t1041 - t1150 * t1042;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t1126, -t1125, 0, -t1152 * t1127 + t1154 * t1128, 0, 0, 0, 0, 0, 0, t1163, -t1160, 0, -t1152 * t1078 + t1154 * t1079, 0, 0, 0, 0, 0, 0, -t1152 * t1103 + t1154 * t1105, -t1152 * t1102 + t1154 * t1104, -t1152 * t1100 + t1154 * t1101, -t1152 * t1058 + t1154 * t1059, 0, 0, 0, 0, 0, 0, t1053, t1049, t1045, -t1152 * t1038 + t1154 * t1039, 0, 0, 0, 0, 0, 0, t1053, t1049, t1045, -t1152 * t1036 + t1154 * t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1125, t1126, 0, t1154 * t1127 + t1152 * t1128, 0, 0, 0, 0, 0, 0, t1160, t1163, 0, t1154 * t1078 + t1152 * t1079, 0, 0, 0, 0, 0, 0, t1154 * t1103 + t1152 * t1105, t1154 * t1102 + t1152 * t1104, t1154 * t1100 + t1152 * t1101, t1154 * t1058 + t1152 * t1059, 0, 0, 0, 0, 0, 0, t1052, t1048, t1044, t1154 * t1038 + t1152 * t1039, 0, 0, 0, 0, 0, 0, t1052, t1048, t1044, t1154 * t1036 + t1152 * t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1181, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1067, 0, 0, 0, 0, 0, 0, t1073, t1071, t1069, t1046, 0, 0, 0, 0, 0, 0, t1073, t1071, t1069, t1040; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1155, -qJDD(1), 0, t1128, 0, 0, 0, 0, 0, 0, t1122, -t1123, 0, t1079, 0, 0, 0, 0, 0, 0, t1105, t1104, t1101, t1059, 0, 0, 0, 0, 0, 0, t1065, t1061, t1057, t1039, 0, 0, 0, 0, 0, 0, t1065, t1061, t1057, t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1155, 0, t1127, 0, 0, 0, 0, 0, 0, t1123, t1122, 0, t1078, 0, 0, 0, 0, 0, 0, t1103, t1102, t1100, t1058, 0, 0, 0, 0, 0, 0, t1064, t1060, t1056, t1038, 0, 0, 0, 0, 0, 0, t1064, t1060, t1056, t1036; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1181, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1067, 0, 0, 0, 0, 0, 0, t1073, t1071, t1069, t1046, 0, 0, 0, 0, 0, 0, t1073, t1071, t1069, t1040; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1155, -qJDD(1), 0, t1099, 0, 0, 0, 0, 0, 0, -t1118, t1117, t1121, t1068, 0, 0, 0, 0, 0, 0, t1074, t1072, t1070, t1047, 0, 0, 0, 0, 0, 0, t1074, t1072, t1070, t1041; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1155, 0, t1098, 0, 0, 0, 0, 0, 0, t1138, -t1171, t1124, -t1092, 0, 0, 0, 0, 0, 0, -t1090, -t1087, -t1080, -t1050, 0, 0, 0, 0, 0, 0, -t1090, -t1087, -t1080, -t1042; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1181, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1067, 0, 0, 0, 0, 0, 0, t1073, t1071, t1069, t1046, 0, 0, 0, 0, 0, 0, t1073, t1071, t1069, t1040; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1142 * t1155, t1131, t1138, t1083, 0, 0, 0, 0, 0, 0, t1091, t1088, t1081, t1051, 0, 0, 0, 0, 0, 0, t1091, t1088, t1081, t1043; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1131, -t1177, -t1171, t1082, 0, 0, 0, 0, 0, 0, -t1094, -t1096, t1112, -t1076, 0, 0, 0, 0, 0, 0, -t1094, -t1096, t1112, -t1066; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1138, t1171, -t1124, t1092, 0, 0, 0, 0, 0, 0, t1090, t1087, t1080, t1050, 0, 0, 0, 0, 0, 0, t1090, t1087, t1080, t1042; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1111, t1109, t1095, t1063, 0, 0, 0, 0, 0, 0, t1111, t1109, t1095, t1055; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1110, t1106, t1097, t1062, 0, 0, 0, 0, 0, 0, t1110, t1106, t1097, t1054; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1094, t1096, -t1112, t1076, 0, 0, 0, 0, 0, 0, t1094, t1096, -t1112, t1066; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1111, t1109, t1095, t1055; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1110, t1106, t1097, t1054; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1094, t1096, -t1112, t1066;];
f_new_reg = t1;
