% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S5PRPRR3
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
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
%
% Output:
% f_new_reg [(3*6)x(6*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S5PRPRR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_invdynf_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_invdynf_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:48:09
% EndTime: 2019-12-05 15:48:11
% DurationCPUTime: 2.03s
% Computational Cost: add. (6449->175), mult. (12027->250), div. (0->0), fcn. (8297->10), ass. (0->143)
t1130 = qJD(4) + qJD(5);
t1167 = qJD(5) + t1130;
t1140 = sin(qJ(5));
t1143 = cos(qJ(5));
t1144 = cos(qJ(4));
t1162 = qJD(2) * t1144;
t1141 = sin(qJ(4));
t1163 = qJD(2) * t1141;
t1094 = t1140 * t1163 - t1143 * t1162;
t1166 = t1094 ^ 2;
t1096 = (t1140 * t1144 + t1141 * t1143) * qJD(2);
t1165 = t1096 ^ 2;
t1164 = t1130 ^ 2;
t1161 = t1096 * t1094;
t1132 = t1144 ^ 2;
t1147 = qJD(2) ^ 2;
t1160 = t1132 * t1147;
t1136 = sin(pkin(8));
t1138 = cos(pkin(8));
t1114 = t1136 * g(1) - t1138 * g(2);
t1159 = t1136 * t1114;
t1158 = qJD(5) - t1130;
t1115 = -t1138 * g(1) - t1136 * g(2);
t1133 = -g(3) + qJDD(1);
t1142 = sin(qJ(2));
t1145 = cos(qJ(2));
t1092 = t1145 * t1115 + t1142 * t1133;
t1090 = -t1147 * pkin(2) + t1092;
t1135 = sin(pkin(9));
t1137 = cos(pkin(9));
t1091 = -t1142 * t1115 + t1145 * t1133;
t1149 = qJDD(2) * pkin(2) + t1091;
t1065 = t1137 * t1090 + t1135 * t1149;
t1063 = -t1147 * pkin(3) + qJDD(2) * pkin(6) + t1065;
t1107 = -qJDD(3) + t1114;
t1050 = t1144 * t1063 - t1141 * t1107;
t1131 = t1141 ^ 2;
t1157 = t1131 + t1132;
t1156 = t1141 * qJDD(2);
t1121 = t1141 * t1147 * t1144;
t1116 = qJDD(4) + t1121;
t1155 = -qJDD(4) - qJDD(5);
t1154 = qJD(4) * t1163;
t1153 = qJD(4) * t1162;
t1049 = -t1141 * t1063 - t1144 * t1107;
t1064 = -t1135 * t1090 + t1137 * t1149;
t1105 = t1153 + t1156;
t1127 = t1144 * qJDD(2);
t1150 = -t1127 + t1154;
t1152 = -t1140 * t1105 - t1143 * t1150;
t1108 = t1137 * qJDD(2) - t1135 * t1147;
t1109 = -t1135 * qJDD(2) - t1137 * t1147;
t1151 = -t1142 * t1108 + t1145 * t1109;
t1077 = t1145 * t1108 + t1142 * t1109;
t1062 = -qJDD(2) * pkin(3) - t1147 * pkin(6) - t1064;
t1148 = -t1143 * t1105 + t1140 * t1150;
t1146 = qJD(4) ^ 2;
t1120 = -t1146 - t1160;
t1119 = -t1131 * t1147 - t1146;
t1118 = qJD(4) * pkin(4) - pkin(7) * t1163;
t1117 = -qJDD(4) + t1121;
t1113 = t1157 * t1147;
t1112 = t1145 * qJDD(2) - t1142 * t1147;
t1111 = -t1142 * qJDD(2) - t1145 * t1147;
t1110 = t1157 * qJDD(2);
t1106 = t1127 - 0.2e1 * t1154;
t1104 = 0.2e1 * t1153 + t1156;
t1101 = t1138 * t1114;
t1089 = -t1164 - t1165;
t1088 = t1144 * t1117 - t1141 * t1119;
t1087 = -t1141 * t1116 + t1144 * t1120;
t1086 = t1141 * t1117 + t1144 * t1119;
t1085 = t1144 * t1116 + t1141 * t1120;
t1081 = t1137 * t1110 - t1135 * t1113;
t1080 = t1135 * t1110 + t1137 * t1113;
t1075 = t1155 - t1161;
t1074 = -t1155 - t1161;
t1073 = -t1164 - t1166;
t1072 = t1137 * t1088 + t1135 * t1104;
t1071 = t1137 * t1087 - t1135 * t1106;
t1070 = t1135 * t1088 - t1137 * t1104;
t1069 = t1135 * t1087 + t1137 * t1106;
t1068 = -t1142 * t1091 + t1145 * t1092;
t1067 = t1145 * t1091 + t1142 * t1092;
t1066 = -t1165 - t1166;
t1060 = -t1142 * t1080 + t1145 * t1081;
t1059 = t1145 * t1080 + t1142 * t1081;
t1058 = t1143 * t1075 - t1140 * t1089;
t1057 = t1140 * t1075 + t1143 * t1089;
t1056 = t1158 * t1094 + t1148;
t1055 = -t1167 * t1094 - t1148;
t1054 = -t1158 * t1096 + t1152;
t1053 = t1167 * t1096 - t1152;
t1052 = t1143 * t1073 - t1140 * t1074;
t1051 = t1140 * t1073 + t1143 * t1074;
t1048 = t1150 * pkin(4) - pkin(7) * t1160 + t1118 * t1163 + t1062;
t1047 = -t1142 * t1070 + t1145 * t1072;
t1046 = -t1142 * t1069 + t1145 * t1071;
t1045 = t1145 * t1070 + t1142 * t1072;
t1044 = t1145 * t1069 + t1142 * t1071;
t1043 = -pkin(4) * t1160 - t1150 * pkin(7) - qJD(4) * t1118 + t1050;
t1042 = (-t1105 + t1153) * pkin(7) + t1116 * pkin(4) + t1049;
t1041 = -t1135 * t1064 + t1137 * t1065;
t1040 = t1137 * t1064 + t1135 * t1065;
t1039 = -t1141 * t1057 + t1144 * t1058;
t1038 = t1144 * t1057 + t1141 * t1058;
t1037 = t1143 * t1054 - t1140 * t1056;
t1036 = t1140 * t1054 + t1143 * t1056;
t1035 = -t1141 * t1051 + t1144 * t1052;
t1034 = t1144 * t1051 + t1141 * t1052;
t1033 = -t1141 * t1049 + t1144 * t1050;
t1032 = t1144 * t1049 + t1141 * t1050;
t1031 = t1137 * t1039 + t1135 * t1055;
t1030 = t1135 * t1039 - t1137 * t1055;
t1029 = t1140 * t1042 + t1143 * t1043;
t1028 = t1143 * t1042 - t1140 * t1043;
t1027 = t1137 * t1033 + t1135 * t1062;
t1026 = t1135 * t1033 - t1137 * t1062;
t1025 = t1137 * t1035 + t1135 * t1053;
t1024 = t1135 * t1035 - t1137 * t1053;
t1023 = -t1142 * t1040 + t1145 * t1041;
t1022 = t1145 * t1040 + t1142 * t1041;
t1021 = -t1141 * t1036 + t1144 * t1037;
t1020 = t1144 * t1036 + t1141 * t1037;
t1019 = t1137 * t1021 + t1135 * t1066;
t1018 = t1135 * t1021 - t1137 * t1066;
t1017 = -t1142 * t1030 + t1145 * t1031;
t1016 = t1145 * t1030 + t1142 * t1031;
t1015 = -t1140 * t1028 + t1143 * t1029;
t1014 = t1143 * t1028 + t1140 * t1029;
t1013 = -t1142 * t1026 + t1145 * t1027;
t1012 = t1145 * t1026 + t1142 * t1027;
t1011 = -t1142 * t1024 + t1145 * t1025;
t1010 = t1145 * t1024 + t1142 * t1025;
t1009 = -t1142 * t1018 + t1145 * t1019;
t1008 = t1145 * t1018 + t1142 * t1019;
t1007 = -t1141 * t1014 + t1144 * t1015;
t1006 = t1144 * t1014 + t1141 * t1015;
t1005 = t1137 * t1007 + t1135 * t1048;
t1004 = t1135 * t1007 - t1137 * t1048;
t1003 = -t1142 * t1004 + t1145 * t1005;
t1002 = t1145 * t1004 + t1142 * t1005;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1138 * t1115 - t1159, 0, 0, 0, 0, 0, 0, t1138 * t1111, -t1138 * t1112, 0, t1138 * t1068 - t1159, 0, 0, 0, 0, 0, 0, t1138 * t1151, -t1138 * t1077, 0, t1138 * t1023 - t1136 * t1107, 0, 0, 0, 0, 0, 0, t1138 * t1046 + t1136 * t1085, t1138 * t1047 + t1136 * t1086, t1138 * t1060, t1138 * t1013 + t1136 * t1032, 0, 0, 0, 0, 0, 0, t1138 * t1011 + t1136 * t1034, t1138 * t1017 + t1136 * t1038, t1138 * t1009 + t1136 * t1020, t1138 * t1003 + t1136 * t1006; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1136 * t1115 + t1101, 0, 0, 0, 0, 0, 0, t1136 * t1111, -t1136 * t1112, 0, t1136 * t1068 + t1101, 0, 0, 0, 0, 0, 0, t1136 * t1151, -t1136 * t1077, 0, t1136 * t1023 + t1138 * t1107, 0, 0, 0, 0, 0, 0, t1136 * t1046 - t1138 * t1085, t1136 * t1047 - t1138 * t1086, t1136 * t1060, t1136 * t1013 - t1138 * t1032, 0, 0, 0, 0, 0, 0, t1136 * t1011 - t1138 * t1034, t1136 * t1017 - t1138 * t1038, t1136 * t1009 - t1138 * t1020, t1136 * t1003 - t1138 * t1006; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t1133, 0, 0, 0, 0, 0, 0, t1112, t1111, 0, t1067, 0, 0, 0, 0, 0, 0, t1077, t1151, 0, t1022, 0, 0, 0, 0, 0, 0, t1044, t1045, t1059, t1012, 0, 0, 0, 0, 0, 0, t1010, t1016, t1008, t1002; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1115, 0, 0, 0, 0, 0, 0, t1111, -t1112, 0, t1068, 0, 0, 0, 0, 0, 0, t1151, -t1077, 0, t1023, 0, 0, 0, 0, 0, 0, t1046, t1047, t1060, t1013, 0, 0, 0, 0, 0, 0, t1011, t1017, t1009, t1003; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1114, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1114, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1107, 0, 0, 0, 0, 0, 0, -t1085, -t1086, 0, -t1032, 0, 0, 0, 0, 0, 0, -t1034, -t1038, -t1020, -t1006; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1133, 0, 0, 0, 0, 0, 0, t1112, t1111, 0, t1067, 0, 0, 0, 0, 0, 0, t1077, t1151, 0, t1022, 0, 0, 0, 0, 0, 0, t1044, t1045, t1059, t1012, 0, 0, 0, 0, 0, 0, t1010, t1016, t1008, t1002; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1147, -qJDD(2), 0, t1092, 0, 0, 0, 0, 0, 0, t1109, -t1108, 0, t1041, 0, 0, 0, 0, 0, 0, t1071, t1072, t1081, t1027, 0, 0, 0, 0, 0, 0, t1025, t1031, t1019, t1005; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1147, 0, t1091, 0, 0, 0, 0, 0, 0, t1108, t1109, 0, t1040, 0, 0, 0, 0, 0, 0, t1069, t1070, t1080, t1026, 0, 0, 0, 0, 0, 0, t1024, t1030, t1018, t1004; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1114, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1107, 0, 0, 0, 0, 0, 0, t1085, t1086, 0, t1032, 0, 0, 0, 0, 0, 0, t1034, t1038, t1020, t1006; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1147, -qJDD(2), 0, t1065, 0, 0, 0, 0, 0, 0, t1087, t1088, t1110, t1033, 0, 0, 0, 0, 0, 0, t1035, t1039, t1021, t1007; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1147, 0, t1064, 0, 0, 0, 0, 0, 0, t1106, -t1104, t1113, -t1062, 0, 0, 0, 0, 0, 0, -t1053, -t1055, -t1066, -t1048; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1107, 0, 0, 0, 0, 0, 0, t1085, t1086, 0, t1032, 0, 0, 0, 0, 0, 0, t1034, t1038, t1020, t1006; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1120, t1117, t1127, t1050, 0, 0, 0, 0, 0, 0, t1052, t1058, t1037, t1015; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1116, t1119, -t1156, t1049, 0, 0, 0, 0, 0, 0, t1051, t1057, t1036, t1014; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1106, t1104, -t1113, t1062, 0, 0, 0, 0, 0, 0, t1053, t1055, t1066, t1048; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1073, t1075, t1054, t1029; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1074, t1089, t1056, t1028; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1053, t1055, t1066, t1048;];
f_new_reg = t1;