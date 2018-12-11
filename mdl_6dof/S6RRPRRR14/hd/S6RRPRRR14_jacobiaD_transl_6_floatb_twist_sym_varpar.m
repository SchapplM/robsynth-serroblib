% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobiaD_transl_6_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:38
% EndTime: 2018-12-10 18:38:44
% DurationCPUTime: 6.57s
% Computational Cost: add. (17923->354), mult. (21062->534), div. (0->0), fcn. (18069->30), ass. (0->209)
t1087 = sin(qJ(6));
t1090 = cos(qJ(6));
t1180 = pkin(6) + qJ(2);
t1159 = cos(t1180) / 0.2e1;
t1181 = pkin(6) - qJ(2);
t1172 = cos(t1181);
t1139 = t1172 / 0.2e1 + t1159;
t1215 = sin(qJ(2));
t1218 = cos(qJ(1));
t1164 = t1218 * t1215;
t1216 = sin(qJ(1));
t1114 = t1139 * t1216 + t1164;
t1169 = sin(t1180);
t1156 = t1169 / 0.2e1;
t1170 = sin(t1181);
t1137 = t1156 - t1170 / 0.2e1;
t1131 = qJD(2) * t1137;
t1217 = cos(qJ(2));
t1163 = t1216 * t1217;
t1046 = qJD(1) * t1114 + qJD(2) * t1163 + t1131 * t1218;
t1165 = t1218 * t1217;
t1115 = -t1216 * t1137 + t1165;
t1066 = t1139 * qJD(2);
t1082 = t1216 * t1215;
t1151 = qJD(2) * t1082 - t1066 * t1218;
t1047 = t1115 * qJD(1) - t1151;
t1176 = pkin(7) + pkin(14);
t1152 = sin(t1176) / 0.2e1;
t1177 = pkin(7) - pkin(14);
t1166 = sin(t1177);
t1068 = t1152 - t1166 / 0.2e1;
t1085 = cos(pkin(14));
t1153 = cos(t1177) / 0.2e1;
t1167 = cos(t1176);
t1069 = t1153 - t1167 / 0.2e1;
t1212 = sin(pkin(6));
t1161 = t1216 * t1212;
t1146 = t1069 * t1161;
t1003 = qJD(1) * t1146 - t1046 * t1068 + t1047 * t1085;
t1133 = t1152 + t1166 / 0.2e1;
t1122 = t1133 * t1212;
t1118 = t1216 * t1122;
t1134 = t1153 + t1167 / 0.2e1;
t1211 = sin(pkin(14));
t1004 = -qJD(1) * t1118 + t1046 * t1134 + t1047 * t1211;
t1089 = sin(qJ(4));
t1213 = cos(pkin(7));
t1154 = t1213 * t1212;
t1142 = t1216 * t1154;
t1084 = sin(pkin(7));
t1187 = t1046 * t1084;
t1126 = qJD(1) * t1142 + t1187;
t1178 = pkin(8) + qJ(4);
t1155 = sin(t1178) / 0.2e1;
t1179 = pkin(8) - qJ(4);
t1168 = sin(t1179);
t1070 = t1155 - t1168 / 0.2e1;
t1129 = qJD(4) * t1070;
t1158 = cos(t1178) / 0.2e1;
t1171 = cos(t1179);
t1071 = t1158 - t1171 / 0.2e1;
t1130 = qJD(4) * t1071;
t1135 = t1155 + t1168 / 0.2e1;
t1138 = t1171 / 0.2e1 + t1158;
t1113 = t1137 * t1218 + t1163;
t1119 = t1218 * t1122;
t1121 = -t1218 * t1139 + t1082;
t1220 = t1113 * t1211 + t1121 * t1134 + t1119;
t1143 = t1218 * t1154;
t1221 = -t1084 * t1121 + t1143;
t1162 = t1218 * t1212;
t1147 = t1069 * t1162;
t1026 = t1068 * t1121 - t1085 * t1113 + t1147;
t1092 = cos(qJ(4));
t1228 = t1026 * t1092;
t957 = qJD(4) * t1228 - t1003 * t1089 - t1004 * t1138 + t1126 * t1135 + t1129 * t1220 - t1130 * t1221;
t1206 = t1090 * t957;
t1088 = sin(qJ(5));
t1091 = cos(qJ(5));
t1063 = t1135 * qJD(4);
t1064 = t1138 * qJD(4);
t1229 = t1026 * t1089;
t1095 = qJD(4) * t1229 + t1003 * t1092 - t1004 * t1070 - t1063 * t1221 - t1064 * t1220 - t1071 * t1126;
t1083 = sin(pkin(8));
t1086 = cos(pkin(8));
t1008 = t1083 * t1220 - t1086 * t1221;
t983 = t1070 * t1220 - t1071 * t1221 + t1228;
t1150 = t1008 * t1091 + t1088 * t983;
t988 = t1004 * t1083 + t1086 * t1126;
t944 = -qJD(5) * t1150 - t1088 * t988 - t1091 * t1095;
t1243 = t1087 * t944 - t1206;
t1209 = t1087 * t957;
t1242 = t1090 * t944 + t1209;
t967 = -t1008 * t1088 + t1091 * t983;
t1241 = qJD(5) * t967 - t1088 * t1095 + t1091 * t988;
t982 = -t1135 * t1221 - t1138 * t1220 + t1229;
t1240 = -t1087 * t967 + t1090 * t982;
t1239 = t1087 * t982 + t1090 * t967;
t1219 = r_i_i_C(3) + pkin(13);
t1107 = t1121 * qJD(1) - qJD(2) * t1165 + t1216 * t1131;
t1132 = -qJD(2) * t1164 - t1066 * t1216;
t1109 = -qJD(1) * t1113 + t1132;
t1097 = -qJD(1) * t1119 - t1107 * t1134 + t1109 * t1211;
t1105 = -qJD(1) * t1143 + t1107 * t1084;
t1094 = t1097 * t1083 - t1105 * t1086;
t1223 = t1114 * t1084 + t1142;
t1072 = t1159 - t1172 / 0.2e1;
t1157 = t1170 / 0.2e1;
t1145 = t1157 - t1169 / 0.2e1;
t1120 = -t1145 * t1218 + t1163;
t1045 = qJD(1) * t1120 - t1132;
t1011 = t1045 * t1134 - t1107 * t1211;
t1182 = t1084 * t1086;
t991 = -t1011 * t1083 - t1045 * t1182;
t1101 = t1114 * t1134 + t1115 * t1211 - t1118;
t1010 = t1083 * t1101 + t1086 * t1223;
t1027 = -t1068 * t1114 + t1085 * t1115 + t1146;
t1023 = t1027 * t1092;
t1098 = -t1070 * t1101 - t1071 * t1223 + t1023;
t1222 = t1010 * t1091 - t1098 * t1088;
t1214 = cos(pkin(6));
t1002 = qJD(1) * t1147 + t1068 * t1107 + t1085 * t1109;
t951 = qJD(4) * t1023 + t1002 * t1089 + t1097 * t1138 - t1101 * t1129 + t1105 * t1135 - t1130 * t1223;
t1210 = t1087 * t951;
t1136 = t1156 + t1157;
t1042 = t1068 * t1136 + t1069 * t1214 - t1072 * t1085;
t1041 = t1042 * t1092;
t1065 = t1136 * qJD(2);
t1067 = t1072 * qJD(2);
t1049 = t1065 * t1085 + t1067 * t1068;
t1108 = t1072 * t1211 + t1133 * t1214 + t1134 * t1136;
t1116 = t1084 * t1136 - t1213 * t1214;
t1117 = -t1065 * t1211 + t1067 * t1134;
t1127 = t1084 * t1135;
t970 = qJD(4) * t1041 + t1049 * t1089 + t1067 * t1127 + t1108 * t1129 + t1116 * t1130 - t1117 * t1138;
t1208 = t1087 * t970;
t1207 = t1090 * t951;
t1205 = t1090 * t970;
t1203 = qJD(6) * t1087;
t1202 = qJD(6) * t1090;
t1060 = -t1145 * t1216 - t1165;
t1048 = qJD(1) * t1060 + t1151;
t1013 = t1046 * t1211 + t1048 * t1134;
t1197 = t1013 * t1083;
t1193 = t1027 * t1089;
t1038 = -t1068 * t1120 - t1085 * t1121;
t1192 = t1038 * t1089;
t1191 = t1038 * t1092;
t1040 = t1060 * t1068 - t1085 * t1114;
t1190 = t1040 * t1089;
t1189 = t1040 * t1092;
t1188 = t1042 * t1089;
t1050 = -t1065 * t1134 - t1067 * t1211;
t1186 = t1050 * t1083;
t1055 = t1072 * t1068 + t1085 * t1136;
t1185 = t1055 * t1089;
t1184 = t1055 * t1092;
t1183 = t1084 * t1071;
t1175 = qJD(5) * t1219;
t1173 = -pkin(11) * t1086 - qJ(3);
t969 = t1010 * t1088 + t1091 * t1098;
t1037 = -t1120 * t1134 + t1121 * t1211;
t1017 = -t1037 * t1083 + t1120 * t1182;
t994 = t1037 * t1070 - t1120 * t1183 + t1191;
t975 = t1017 * t1088 + t1091 * t994;
t1039 = t1060 * t1134 + t1114 * t1211;
t1018 = -t1039 * t1083 - t1060 * t1182;
t996 = t1039 * t1070 + t1060 * t1183 + t1189;
t976 = t1018 * t1088 + t1091 * t996;
t1021 = -t1083 * t1108 - t1086 * t1116;
t1103 = t1070 * t1108 + t1071 * t1116 + t1041;
t1149 = t1021 * t1091 - t1088 * t1103;
t978 = t1021 * t1088 + t1091 * t1103;
t1148 = r_i_i_C(1) * t1090 - r_i_i_C(2) * t1087 + pkin(5);
t1054 = t1072 * t1134 - t1136 * t1211;
t1007 = t1054 * t1070 + t1072 * t1183 + t1184;
t1035 = -t1054 * t1083 - t1072 * t1182;
t990 = t1007 * t1091 + t1035 * t1088;
t1144 = qJD(6) * (-t1087 * r_i_i_C(1) - t1090 * r_i_i_C(2));
t1141 = qJD(5) * t1148;
t1123 = t1084 * t1130;
t1100 = -qJD(4) * t1188 + t1049 * t1092 - t1063 * t1116 + t1064 * t1108 + t1067 * t1183 + t1070 * t1117;
t1093 = -qJD(4) * t1193 + t1002 * t1092 + t1063 * t1223 - t1064 * t1101 - t1070 * t1097 + t1071 * t1105;
t1051 = -t1065 * t1068 + t1067 * t1085;
t1031 = t1065 * t1182 - t1186;
t1030 = -t1067 * t1182 - t1083 * t1117;
t1014 = -t1046 * t1085 + t1048 * t1068;
t1012 = t1045 * t1068 + t1085 * t1107;
t1006 = -t1054 * t1138 + t1072 * t1127 + t1185;
t997 = -t1108 * t1138 + t1116 * t1135 + t1188;
t995 = -t1039 * t1138 + t1060 * t1127 + t1190;
t993 = -t1037 * t1138 - t1120 * t1127 + t1192;
t992 = -t1048 * t1182 - t1197;
t984 = t1101 * t1138 - t1135 * t1223 + t1193;
t974 = -qJD(4) * t1185 + t1050 * t1070 + t1051 * t1092 + t1054 * t1064 + (-t1063 * t1072 - t1065 * t1071) * t1084;
t973 = qJD(4) * t1184 - t1050 * t1138 + t1051 * t1089 + t1054 * t1129 - t1065 * t1127 + t1072 * t1123;
t964 = -qJD(4) * t1192 + t1013 * t1070 + t1014 * t1092 + t1037 * t1064 + (t1048 * t1071 + t1063 * t1120) * t1084;
t963 = qJD(4) * t1191 - t1013 * t1138 + t1014 * t1089 + t1037 * t1129 + t1048 * t1127 - t1120 * t1123;
t962 = -qJD(4) * t1190 + t1011 * t1070 + t1012 * t1092 + t1039 * t1064 + (t1045 * t1071 - t1060 * t1063) * t1084;
t961 = qJD(4) * t1189 - t1011 * t1138 + t1012 * t1089 + t1039 * t1129 + t1045 * t1127 + t1060 * t1123;
t960 = t1031 * t1088 + t1091 * t974 + (-t1007 * t1088 + t1035 * t1091) * qJD(5);
t950 = qJD(5) * t1149 + t1030 * t1088 + t1091 * t1100;
t948 = t1088 * t992 + t1091 * t964 + (t1017 * t1091 - t1088 * t994) * qJD(5);
t946 = t1088 * t991 + t1091 * t962 + (t1018 * t1091 - t1088 * t996) * qJD(5);
t940 = t1222 * qJD(5) + t1088 * t1094 + t1093 * t1091;
t939 = qJD(5) * t969 + t1088 * t1093 - t1091 * t1094;
t938 = t1210 + t1090 * t940 + (-t1087 * t969 + t1090 * t984) * qJD(6);
t937 = -t1087 * t940 + t1207 + (-t1087 * t984 - t1090 * t969) * qJD(6);
t1 = [t1242 * r_i_i_C(1) - t1243 * r_i_i_C(2) + t944 * pkin(5) - t1095 * pkin(4) + t957 * pkin(12) - t1003 * pkin(3) - t1047 * pkin(2) - qJ(3) * t1187 + t1219 * t1241 + (r_i_i_C(1) * t1240 - r_i_i_C(2) * t1239) * qJD(6) + t1221 * qJD(3) - t988 * pkin(11) + (-pkin(1) * t1218 - pkin(10) * t1161 - qJ(3) * t1142) * qJD(1) (t1087 * t961 + t1090 * t946 + (-t1087 * t976 + t1090 * t995) * qJD(6)) * r_i_i_C(1) + (-t1087 * t946 + t1090 * t961 + (-t1087 * t995 - t1090 * t976) * qJD(6)) * r_i_i_C(2) + t946 * pkin(5) + t962 * pkin(4) + t961 * pkin(12) + t1012 * pkin(3) + t1107 * pkin(2) + t1219 * (qJD(5) * t976 + t1088 * t962 - t1091 * t991) + (-qJ(3) * t1045 - qJD(3) * t1060) * t1084 + t991 * pkin(11), -t1105 (t1087 * t1093 + t1098 * t1202) * r_i_i_C(1) + (t1090 * t1093 - t1098 * t1203) * r_i_i_C(2) - t951 * pkin(4) + t1093 * pkin(12) + (t984 * t1141 - t1219 * t951) * t1088 + ((t1203 * t984 - t1207) * r_i_i_C(1) + (t1202 * t984 + t1210) * r_i_i_C(2) - t951 * pkin(5) - t984 * t1175) * t1091, t1222 * t1144 - t1148 * t939 + t1219 * t940, r_i_i_C(1) * t937 - r_i_i_C(2) * t938; t1109 * pkin(2) + t1002 * pkin(3) + t1093 * pkin(4) + t940 * pkin(5) + t951 * pkin(12) + t938 * r_i_i_C(1) + t937 * r_i_i_C(2) + t1219 * t939 + t1223 * qJD(3) + (-pkin(1) * t1216 + pkin(10) * t1162) * qJD(1) - t1105 * qJ(3) + t1094 * pkin(11) (t1087 * t963 + t1090 * t948) * r_i_i_C(1) + (-t1087 * t948 + t1090 * t963) * r_i_i_C(2) + t948 * pkin(5) + t964 * pkin(4) + t963 * pkin(12) + t1014 * pkin(3) - pkin(11) * t1197 - t1046 * pkin(2) + t1219 * (qJD(5) * t975 + t1088 * t964 - t1091 * t992) + ((-t1087 * t975 + t1090 * t993) * r_i_i_C(1) + (-t1087 * t993 - t1090 * t975) * r_i_i_C(2)) * qJD(6) + (qJD(3) * t1120 + t1048 * t1173) * t1084, t1126 (t1087 * t1095 - t1202 * t983) * r_i_i_C(1) + (t1090 * t1095 + t1203 * t983) * r_i_i_C(2) + t957 * pkin(4) + t1095 * pkin(12) + (-t1141 * t982 + t1219 * t957) * t1088 + ((-t1203 * t982 + t1206) * r_i_i_C(1) + (-t1202 * t982 - t1209) * r_i_i_C(2) + t957 * pkin(5) + t982 * t1175) * t1091, t1150 * t1144 + t1148 * t1241 - t1219 * t944, t1243 * r_i_i_C(1) + t1242 * r_i_i_C(2) + (r_i_i_C(1) * t1239 + r_i_i_C(2) * t1240) * qJD(6); 0 (t1087 * t973 + t1090 * t960) * r_i_i_C(1) + (-t1087 * t960 + t1090 * t973) * r_i_i_C(2) + t960 * pkin(5) + t974 * pkin(4) + t973 * pkin(12) + t1051 * pkin(3) - pkin(11) * t1186 + t1067 * pkin(2) + t1219 * (qJD(5) * t990 - t1031 * t1091 + t1088 * t974) + ((t1006 * t1090 - t1087 * t990) * r_i_i_C(1) + (-t1006 * t1087 - t1090 * t990) * r_i_i_C(2)) * qJD(6) + (-t1072 * qJD(3) - t1065 * t1173) * t1084, -t1067 * t1084 (t1087 * t1100 + t1103 * t1202) * r_i_i_C(1) + (t1090 * t1100 - t1103 * t1203) * r_i_i_C(2) - t970 * pkin(4) + t1100 * pkin(12) + (t997 * t1141 - t1219 * t970) * t1088 + ((t1203 * t997 - t1205) * r_i_i_C(1) + (t1202 * t997 + t1208) * r_i_i_C(2) - t970 * pkin(5) - t997 * t1175) * t1091, t1219 * t950 + t1149 * t1144 + t1148 * (-qJD(5) * t978 + t1030 * t1091 - t1088 * t1100) (-t1087 * t950 + t1205) * r_i_i_C(1) + (-t1090 * t950 - t1208) * r_i_i_C(2) + ((-t1087 * t997 - t1090 * t978) * r_i_i_C(1) + (t1087 * t978 - t1090 * t997) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
