% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPRR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:03:21
% EndTime: 2022-01-20 11:03:39
% DurationCPUTime: 18.06s
% Computational Cost: add. (133512->604), mult. (187064->863), div. (0->0), fcn. (132391->10), ass. (0->418)
t1145 = cos(qJ(2));
t1195 = qJD(1) + qJD(2);
t1193 = t1195 ^ 2;
t1135 = qJDD(1) + qJDD(2);
t1141 = sin(qJ(2));
t1197 = t1141 * t1135;
t1105 = t1145 * t1193 + t1197;
t1086 = pkin(6) * t1105 - g(3) * t1145;
t1142 = sin(qJ(1));
t1146 = cos(qJ(1));
t1196 = t1145 * t1135;
t1108 = t1141 * t1193 - t1196;
t1162 = t1105 * t1146 - t1108 * t1142;
t1237 = pkin(6) * t1108 - g(3) * t1141;
t1249 = pkin(5) * t1162 + t1146 * t1086 - t1142 * t1237;
t1234 = t1105 * t1142 + t1108 * t1146;
t1248 = pkin(5) * t1234 + t1142 * t1086 + t1146 * t1237;
t1120 = g(1) * t1146 + g(2) * t1142;
t1148 = qJD(1) ^ 2;
t1111 = -pkin(1) * t1148 - t1120;
t1119 = g(1) * t1142 - g(2) * t1146;
t1158 = qJDD(1) * pkin(1) + t1119;
t1070 = t1111 * t1141 - t1145 * t1158;
t1071 = t1145 * t1111 + t1141 * t1158;
t1174 = t1070 * t1141 + t1071 * t1145;
t1014 = t1070 * t1145 - t1071 * t1141;
t1211 = t1014 * t1142;
t1246 = t1146 * t1174 + t1211;
t1210 = t1014 * t1146;
t1245 = -t1142 * t1174 + t1210;
t1139 = sin(qJ(5));
t1134 = qJDD(4) + qJDD(5);
t1137 = sin(pkin(9));
t1138 = cos(pkin(9));
t1144 = cos(qJ(4));
t1181 = t1144 * t1195;
t1140 = sin(qJ(4));
t1182 = t1140 * t1195;
t1090 = t1137 * t1182 - t1138 * t1181;
t1091 = t1137 * t1181 + t1138 * t1182;
t1143 = cos(qJ(5));
t1036 = t1090 * t1143 + t1091 * t1139;
t1038 = -t1090 * t1139 + t1091 * t1143;
t986 = t1038 * t1036;
t1232 = -t986 + t1134;
t1242 = t1139 * t1232;
t1056 = t1091 * t1090;
t1229 = qJDD(4) - t1056;
t1241 = t1140 * t1229;
t1240 = t1143 * t1232;
t1239 = t1144 * t1229;
t1165 = -pkin(2) * t1193 + qJ(3) * t1135 + 0.2e1 * qJD(3) * t1195 + t1071;
t1238 = t1135 * pkin(7) + t1165;
t1151 = t1138 ^ 2;
t1129 = t1151 * t1193;
t1149 = t1137 ^ 2;
t1180 = t1149 * t1193;
t1103 = t1180 + t1129;
t1095 = t1103 * t1138;
t1066 = -t1095 * t1141 + t1138 * t1196;
t1068 = t1095 * t1145 + t1138 * t1197;
t1236 = t1066 * t1146 - t1068 * t1142;
t1235 = t1066 * t1142 + t1068 * t1146;
t1136 = qJD(4) + qJD(5);
t1029 = t1136 * t1036;
t1081 = t1090 * qJD(4);
t1088 = (t1137 * t1144 + t1138 * t1140) * t1135;
t1053 = -t1081 + t1088;
t1126 = t1137 * t1135;
t1128 = t1138 * t1135;
t1087 = t1126 * t1140 - t1128 * t1144;
t1203 = t1091 * qJD(4);
t1154 = t1087 + t1203;
t968 = -qJD(5) * t1036 + t1053 * t1143 - t1139 * t1154;
t1233 = -t1029 + t968;
t1221 = t1138 * g(3);
t1023 = t1137 * t1165 + t1221;
t1222 = t1137 * g(3);
t1024 = t1138 * t1165 - t1222;
t972 = t1023 * t1137 + t1024 * t1138;
t1175 = t1053 * t1139 + t1143 * t1154;
t916 = (qJD(5) - t1136) * t1038 + t1175;
t1032 = t1036 ^ 2;
t1033 = t1038 ^ 2;
t1228 = t1090 ^ 2;
t1089 = t1091 ^ 2;
t1133 = t1136 ^ 2;
t1183 = t1138 * t1193;
t996 = -t1221 + (pkin(3) * t1183 - t1238) * t1137;
t997 = -pkin(3) * t1129 + t1138 * t1238 - t1222;
t942 = t1140 * t997 - t1144 * t996;
t943 = t1140 * t996 + t1144 * t997;
t874 = t1140 * t943 - t1144 * t942;
t1227 = pkin(3) * t874;
t989 = -t1087 * t1140 - t1088 * t1144;
t1226 = pkin(3) * t989;
t1220 = t1137 * t874;
t1219 = t1138 * t874;
t1046 = -pkin(2) * t1135 - qJ(3) * t1193 + qJDD(3) + t1070;
t1017 = -pkin(3) * t1128 - pkin(7) * t1103 + t1046;
t1074 = qJD(4) * pkin(4) - pkin(8) * t1091;
t931 = pkin(4) * t1154 - pkin(8) * t1228 + t1074 * t1091 + t1017;
t1218 = t1139 * t931;
t983 = t986 + t1134;
t1217 = t1139 * t983;
t882 = (-t1053 - t1081) * pkin(8) + t1229 * pkin(4) - t942;
t885 = -pkin(4) * t1228 - pkin(8) * t1154 - qJD(4) * t1074 + t943;
t842 = t1139 * t885 - t1143 * t882;
t843 = t1139 * t882 + t1143 * t885;
t797 = t1139 * t843 - t1143 * t842;
t1216 = t1140 * t797;
t1215 = t1143 * t931;
t1214 = t1143 * t983;
t1213 = t1144 * t797;
t1212 = -pkin(2) * t1046 + qJ(3) * t972;
t1209 = t1017 * t1140;
t1208 = t1017 * t1144;
t1049 = qJDD(4) + t1056;
t1207 = t1049 * t1140;
t1206 = t1049 * t1144;
t1205 = t1090 * t1140;
t1204 = t1090 * t1144;
t1202 = t1091 * t1140;
t1201 = t1091 * t1144;
t1200 = t1136 * t1038;
t1199 = t1136 * t1139;
t1198 = t1136 * t1143;
t1039 = t1137 * t1046;
t1040 = t1138 * t1046;
t798 = t1139 * t842 + t1143 * t843;
t763 = t1140 * t798 + t1213;
t796 = pkin(4) * t797;
t1192 = pkin(3) * t763 + t796;
t920 = t1029 + t968;
t865 = -t1139 * t916 - t1143 * t920;
t867 = t1139 * t920 - t1143 * t916;
t821 = t1140 * t867 + t1144 * t865;
t863 = pkin(4) * t865;
t1191 = pkin(3) * t821 + t863;
t1190 = t1141 * t986;
t1189 = t1145 * t986;
t1188 = t1141 * t1056;
t1187 = t1145 * t1056;
t1186 = t1137 * t1128;
t1185 = pkin(2) * t1128 - qJ(3) * t1095 - t1040;
t875 = t1140 * t942 + t1144 * t943;
t764 = t1144 * t798 - t1216;
t791 = -pkin(4) * t931 + pkin(8) * t798;
t744 = -pkin(3) * t931 + pkin(7) * t764 - pkin(8) * t1216 + t1144 * t791;
t748 = -pkin(7) * t763 - pkin(8) * t1213 - t1140 * t791;
t752 = -t1137 * t763 + t1138 * t764;
t1179 = -pkin(2) * t931 + qJ(3) * t752 + t1137 * t748 + t1138 * t744;
t966 = -t1032 - t1033;
t773 = -pkin(4) * t966 + pkin(8) * t867 + t798;
t782 = -pkin(8) * t865 - t797;
t823 = -t1140 * t865 + t1144 * t867;
t754 = -pkin(3) * t966 + pkin(7) * t823 + t1140 * t782 + t1144 * t773;
t756 = -pkin(7) * t821 - t1140 * t773 + t1144 * t782;
t780 = -t1137 * t821 + t1138 * t823;
t1178 = -pkin(2) * t966 + qJ(3) * t780 + t1137 * t756 + t1138 * t754;
t981 = -t1133 - t1032;
t896 = t1143 * t981 - t1242;
t915 = (qJD(5) + t1136) * t1038 + t1175;
t837 = -pkin(4) * t915 + pkin(8) * t896 - t1215;
t895 = t1139 * t981 + t1240;
t856 = -t1140 * t895 + t1144 * t896;
t869 = -pkin(8) * t895 + t1218;
t776 = -pkin(3) * t915 + pkin(7) * t856 + t1140 * t869 + t1144 * t837;
t855 = t1140 * t896 + t1144 * t895;
t788 = -pkin(7) * t855 - t1140 * t837 + t1144 * t869;
t815 = -t1137 * t855 + t1138 * t856;
t1177 = -pkin(2) * t915 + qJ(3) * t815 + t1137 * t788 + t1138 * t776;
t1025 = -t1033 - t1133;
t941 = -t1025 * t1139 - t1214;
t840 = -pkin(4) * t1233 + pkin(8) * t941 + t1218;
t940 = t1025 * t1143 - t1217;
t871 = -pkin(8) * t940 + t1215;
t873 = -t1140 * t940 + t1144 * t941;
t785 = -pkin(3) * t1233 + pkin(7) * t873 + t1140 * t871 + t1144 * t840;
t872 = t1140 * t941 + t1144 * t940;
t793 = -pkin(7) * t872 - t1140 * t840 + t1144 * t871;
t829 = -t1137 * t872 + t1138 * t873;
t1176 = -pkin(2) * t1233 + qJ(3) * t829 + t1137 * t793 + t1138 * t785;
t1173 = -t1119 * t1142 - t1120 * t1146;
t1022 = -t1089 - t1228;
t991 = -t1087 * t1144 + t1088 * t1140;
t846 = -pkin(3) * t1022 + pkin(7) * t991 + t875;
t854 = -pkin(7) * t989 - t874;
t923 = -t1137 * t989 + t1138 * t991;
t1172 = -pkin(2) * t1022 + qJ(3) * t923 + t1137 * t854 + t1138 * t846;
t1051 = t1087 + 0.2e1 * t1203;
t1147 = qJD(4) ^ 2;
t1047 = -t1147 - t1228;
t988 = t1047 * t1144 - t1241;
t893 = -pkin(3) * t1051 + pkin(7) * t988 - t1208;
t987 = t1047 * t1140 + t1239;
t911 = -t1137 * t987 + t1138 * t988;
t937 = -pkin(7) * t987 + t1209;
t1171 = -pkin(2) * t1051 + qJ(3) * t911 + t1137 * t937 + t1138 * t893;
t1052 = -0.2e1 * t1081 + t1088;
t1077 = -t1089 - t1147;
t1004 = -t1077 * t1140 - t1206;
t925 = -pkin(3) * t1052 + pkin(7) * t1004 + t1209;
t1001 = t1077 * t1144 - t1207;
t954 = -t1001 * t1137 + t1004 * t1138;
t960 = -pkin(7) * t1001 + t1208;
t1170 = -pkin(2) * t1052 + qJ(3) * t954 + t1137 * t960 + t1138 * t925;
t1169 = pkin(3) * t1001 - t943;
t1125 = t1149 * t1135;
t1127 = t1151 * t1135;
t1100 = t1127 + t1125;
t1168 = pkin(2) * t1103 + qJ(3) * t1100 + t972;
t1167 = pkin(4) * t895 - t842;
t1114 = qJDD(1) * t1146 - t1142 * t1148;
t1166 = -pkin(5) * t1114 - g(3) * t1142;
t1115 = t1137 * t1183;
t971 = t1023 * t1138 - t1024 * t1137;
t1072 = t1105 * t1138 * t1137;
t1073 = -t1115 * t1141 + t1145 * t1186;
t1164 = t1072 * t1146 + t1073 * t1142;
t1163 = t1072 * t1142 - t1073 * t1146;
t1161 = t1119 * t1146 - t1120 * t1142;
t1094 = t1103 * t1137;
t1160 = -pkin(2) * t1126 + qJ(3) * t1094 + t1039;
t1159 = pkin(3) * t987 - t942;
t1157 = pkin(4) * t940 - t843;
t832 = t1138 * t875 - t1220;
t862 = -pkin(3) * t1017 + pkin(7) * t875;
t1156 = -pkin(2) * t1017 - pkin(7) * t1220 + qJ(3) * t832 + t1138 * t862;
t1155 = pkin(3) * t855 + t1167;
t1153 = pkin(3) * t872 + t1157;
t1113 = qJDD(1) * t1142 + t1146 * t1148;
t1112 = 0.2e1 * t1186;
t1104 = t1180 - t1129;
t1101 = t1127 - t1125;
t1098 = -pkin(5) * t1113 + g(3) * t1146;
t1076 = -t1089 + t1147;
t1075 = -t1147 + t1228;
t1067 = t1094 * t1145 + t1137 * t1197;
t1064 = t1094 * t1141 - t1137 * t1196;
t1060 = t1101 * t1145 + t1104 * t1141;
t1059 = t1100 * t1145 - t1103 * t1141;
t1058 = t1101 * t1141 - t1104 * t1145;
t1057 = t1100 * t1141 + t1103 * t1145;
t1054 = t1089 - t1228;
t1035 = -pkin(1) * t1105 - t1071;
t1034 = -pkin(1) * t1108 - t1070;
t1031 = (t1202 - t1204) * qJD(4);
t1030 = (-t1201 - t1205) * qJD(4);
t1027 = -t1033 + t1133;
t1026 = t1032 - t1133;
t1011 = -t1064 * t1142 + t1067 * t1146;
t1010 = t1064 * t1146 + t1067 * t1142;
t1009 = pkin(1) * t1014;
t1008 = -qJD(4) * t1202 + t1053 * t1144;
t1007 = qJD(4) * t1201 + t1053 * t1140;
t1006 = qJD(4) * t1204 + t1140 * t1154;
t1005 = qJD(4) * t1205 - t1144 * t1154;
t1003 = -t1076 * t1140 + t1239;
t1002 = t1075 * t1144 - t1207;
t1000 = t1076 * t1144 + t1241;
t999 = t1075 * t1140 + t1206;
t998 = pkin(1) * g(3) + pkin(6) * t1174;
t992 = -t1051 * t1144 - t1052 * t1140;
t990 = -t1051 * t1140 + t1052 * t1144;
t985 = t1033 - t1032;
t980 = -t1030 * t1137 + t1031 * t1138;
t979 = t1030 * t1138 + t1031 * t1137;
t978 = (-t1036 * t1143 + t1038 * t1139) * t1136;
t977 = (-t1036 * t1139 - t1038 * t1143) * t1136;
t976 = pkin(1) * t1066 + t1185;
t975 = pkin(1) * t1064 + t1160;
t974 = qJDD(4) * t1141 + t1145 * t980;
t973 = -qJDD(4) * t1145 + t1141 * t980;
t967 = -qJD(5) * t1038 - t1175;
t964 = -pkin(6) * t1064 - t1024 * t1141 + t1040 * t1145;
t963 = -pkin(6) * t1066 - t1023 * t1141 + t1039 * t1145;
t962 = pkin(6) * t1067 + t1024 * t1145 + t1040 * t1141;
t961 = -pkin(6) * t1068 + t1023 * t1145 + t1039 * t1141;
t958 = -t1007 * t1137 + t1008 * t1138;
t957 = -t1005 * t1137 + t1006 * t1138;
t956 = t1007 * t1138 + t1008 * t1137;
t955 = t1005 * t1138 + t1006 * t1137;
t953 = -t1000 * t1137 + t1003 * t1138;
t952 = t1002 * t1138 - t1137 * t999;
t951 = t1001 * t1138 + t1004 * t1137;
t950 = t1000 * t1138 + t1003 * t1137;
t949 = t1002 * t1137 + t1138 * t999;
t947 = t1026 * t1143 - t1217;
t946 = -t1027 * t1139 + t1240;
t945 = t1026 * t1139 + t1214;
t944 = t1027 * t1143 + t1242;
t936 = -pkin(6) * t1057 + t1145 * t971;
t935 = pkin(6) * t1059 + t1141 * t971;
t933 = t1046 * t1141 + t1145 * t972;
t932 = -t1046 * t1145 + t1141 * t972;
t929 = t1088 * t1141 + t1145 * t953;
t928 = -t1087 * t1141 + t1145 * t952;
t927 = -t1088 * t1145 + t1141 * t953;
t926 = t1087 * t1145 + t1141 * t952;
t924 = -t1137 * t990 + t1138 * t992;
t922 = t1137 * t992 + t1138 * t990;
t921 = t1137 * t991 + t1138 * t989;
t910 = t1137 * t988 + t1138 * t987;
t907 = t1145 * t958 + t1188;
t906 = t1145 * t957 - t1188;
t905 = t1141 * t958 - t1187;
t904 = t1141 * t957 + t1187;
t903 = pkin(1) * t1057 + t1168;
t902 = -t1038 * t1199 + t1143 * t968;
t901 = t1038 * t1198 + t1139 * t968;
t900 = t1036 * t1198 - t1139 * t967;
t899 = t1036 * t1199 + t1143 * t967;
t898 = t1052 * t1141 + t1145 * t954;
t897 = -t1052 * t1145 + t1141 * t954;
t891 = -t1140 * t977 + t1144 * t978;
t890 = t1140 * t978 + t1144 * t977;
t889 = t1054 * t1141 + t1145 * t924;
t888 = -t1054 * t1145 + t1141 * t924;
t887 = t1051 * t1141 + t1145 * t911;
t886 = -t1051 * t1145 + t1141 * t911;
t884 = t1022 * t1141 + t1145 * t923;
t883 = -t1022 * t1145 + t1141 * t923;
t880 = -pkin(2) * t921 - t1226;
t879 = -t1140 * t945 + t1144 * t947;
t878 = -t1140 * t944 + t1144 * t946;
t877 = t1140 * t947 + t1144 * t945;
t876 = t1140 * t946 + t1144 * t944;
t870 = pkin(1) * t932 + t1212;
t868 = -pkin(2) * t951 - t1169;
t866 = -t1139 * t1233 - t1143 * t915;
t864 = -t1139 * t915 + t1143 * t1233;
t860 = -t1140 * t901 + t1144 * t902;
t859 = -t1140 * t899 + t1144 * t900;
t858 = t1140 * t902 + t1144 * t901;
t857 = t1140 * t900 + t1144 * t899;
t852 = -pkin(6) * t932 - (pkin(2) * t1141 - qJ(3) * t1145) * t971;
t851 = -pkin(2) * t910 - t1159;
t850 = -t1137 * t890 + t1138 * t891;
t849 = t1137 * t891 + t1138 * t890;
t848 = t1134 * t1141 + t1145 * t850;
t847 = -t1134 * t1145 + t1141 * t850;
t844 = -qJ(3) * t951 - t1137 * t925 + t1138 * t960;
t839 = pkin(6) * t933 - (-pkin(2) * t1145 - qJ(3) * t1141 - pkin(1)) * t971;
t838 = -qJ(3) * t910 - t1137 * t893 + t1138 * t937;
t836 = -t1137 * t877 + t1138 * t879;
t835 = -t1137 * t876 + t1138 * t878;
t834 = t1137 * t879 + t1138 * t877;
t833 = t1137 * t878 + t1138 * t876;
t831 = t1137 * t875 + t1219;
t828 = t1137 * t873 + t1138 * t872;
t826 = t1017 * t1141 + t1145 * t832;
t825 = -t1017 * t1145 + t1141 * t832;
t824 = pkin(1) * t897 + t1170;
t822 = -t1140 * t864 + t1144 * t866;
t820 = t1140 * t866 + t1144 * t864;
t819 = -t1137 * t858 + t1138 * t860;
t818 = -t1137 * t857 + t1138 * t859;
t817 = t1137 * t860 + t1138 * t858;
t816 = t1137 * t859 + t1138 * t857;
t814 = t1137 * t856 + t1138 * t855;
t812 = -t1141 * t916 + t1145 * t836;
t811 = t1141 * t920 + t1145 * t835;
t810 = t1141 * t836 + t1145 * t916;
t809 = t1141 * t835 - t1145 * t920;
t808 = t1141 * t1233 + t1145 * t829;
t807 = t1141 * t829 - t1145 * t1233;
t806 = pkin(1) * t886 + t1171;
t805 = t1145 * t819 + t1190;
t804 = t1145 * t818 - t1190;
t803 = t1141 * t819 - t1189;
t802 = t1141 * t818 + t1189;
t801 = -pkin(2) * t831 - t1227;
t800 = t1141 * t915 + t1145 * t815;
t799 = t1141 * t815 - t1145 * t915;
t795 = -pkin(6) * t897 - t1141 * t868 + t1145 * t844;
t794 = -qJ(3) * t921 - t1137 * t846 + t1138 * t854;
t790 = -pkin(1) * t951 + pkin(6) * t898 + t1141 * t844 + t1145 * t868;
t789 = -pkin(6) * t886 - t1141 * t851 + t1145 * t838;
t787 = -pkin(7) * t1219 - qJ(3) * t831 - t1137 * t862;
t783 = pkin(1) * t883 + t1172;
t781 = -pkin(1) * t910 + pkin(6) * t887 + t1141 * t838 + t1145 * t851;
t779 = -t1137 * t820 + t1138 * t822;
t778 = t1137 * t823 + t1138 * t821;
t777 = t1137 * t822 + t1138 * t820;
t772 = -pkin(6) * t883 - t1141 * t880 + t1145 * t794;
t771 = t1141 * t985 + t1145 * t779;
t770 = t1141 * t779 - t1145 * t985;
t769 = -pkin(2) * t828 - t1153;
t768 = t1141 * t966 + t1145 * t780;
t767 = t1141 * t780 - t1145 * t966;
t766 = -pkin(1) * t921 + pkin(6) * t884 + t1141 * t794 + t1145 * t880;
t765 = -pkin(2) * t814 - t1155;
t762 = pkin(1) * t825 + t1156;
t761 = -pkin(2) * t778 - t1191;
t760 = -pkin(6) * t825 - t1141 * t801 + t1145 * t787;
t759 = -qJ(3) * t828 - t1137 * t785 + t1138 * t793;
t758 = -pkin(1) * t831 + pkin(6) * t826 + t1141 * t787 + t1145 * t801;
t757 = -qJ(3) * t814 - t1137 * t776 + t1138 * t788;
t751 = t1137 * t764 + t1138 * t763;
t749 = pkin(1) * t807 + t1176;
t746 = t1141 * t931 + t1145 * t752;
t745 = t1141 * t752 - t1145 * t931;
t742 = pkin(1) * t799 + t1177;
t741 = -pkin(6) * t807 - t1141 * t769 + t1145 * t759;
t740 = -pkin(1) * t828 + pkin(6) * t808 + t1141 * t759 + t1145 * t769;
t739 = -pkin(6) * t799 - t1141 * t765 + t1145 * t757;
t738 = -pkin(2) * t751 - t1192;
t737 = -pkin(1) * t814 + pkin(6) * t800 + t1141 * t757 + t1145 * t765;
t736 = -qJ(3) * t778 - t1137 * t754 + t1138 * t756;
t735 = pkin(1) * t767 + t1178;
t734 = -qJ(3) * t751 - t1137 * t744 + t1138 * t748;
t733 = -pkin(6) * t767 - t1141 * t761 + t1145 * t736;
t732 = -pkin(1) * t778 + pkin(6) * t768 + t1141 * t736 + t1145 * t761;
t731 = pkin(1) * t745 + t1179;
t730 = -pkin(6) * t745 - t1141 * t738 + t1145 * t734;
t729 = -pkin(1) * t751 + pkin(6) * t746 + t1141 * t734 + t1145 * t738;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1114, 0, -t1113, 0, t1166, -t1098, -t1161, -pkin(5) * t1161, 0, 0, -t1234, 0, -t1162, 0, t1248, t1249, t1245, pkin(5) * t1245 + pkin(6) * t1210 - t1142 * t998, -t1163, -t1058 * t1142 + t1060 * t1146, t1011, t1163, t1235, 0, -pkin(5) * t1236 - t1142 * t961 + t1146 * t963, -pkin(5) * t1010 - t1142 * t962 + t1146 * t964, t1146 * t936 - t1142 * t935 - pkin(5) * (t1057 * t1146 + t1059 * t1142), t1146 * t852 - t1142 * t839 - pkin(5) * (t1142 * t933 + t1146 * t932), -t1142 * t905 + t1146 * t907, -t1142 * t888 + t1146 * t889, -t1142 * t927 + t1146 * t929, -t1142 * t904 + t1146 * t906, -t1142 * t926 + t1146 * t928, -t1142 * t973 + t1146 * t974, t1146 * t789 - t1142 * t781 - pkin(5) * (t1142 * t887 + t1146 * t886), t1146 * t795 - t1142 * t790 - pkin(5) * (t1142 * t898 + t1146 * t897), t1146 * t772 - t1142 * t766 - pkin(5) * (t1142 * t884 + t1146 * t883), t1146 * t760 - t1142 * t758 - pkin(5) * (t1142 * t826 + t1146 * t825), -t1142 * t803 + t1146 * t805, -t1142 * t770 + t1146 * t771, -t1142 * t809 + t1146 * t811, -t1142 * t802 + t1146 * t804, -t1142 * t810 + t1146 * t812, -t1142 * t847 + t1146 * t848, t1146 * t739 - t1142 * t737 - pkin(5) * (t1142 * t800 + t1146 * t799), t1146 * t741 - t1142 * t740 - pkin(5) * (t1142 * t808 + t1146 * t807), t1146 * t733 - t1142 * t732 - pkin(5) * (t1142 * t768 + t1146 * t767), t1146 * t730 - t1142 * t729 - pkin(5) * (t1142 * t746 + t1146 * t745); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1113, 0, t1114, 0, t1098, t1166, t1173, pkin(5) * t1173, 0, 0, t1162, 0, -t1234, 0, -t1249, t1248, t1246, pkin(5) * t1246 + pkin(6) * t1211 + t1146 * t998, t1164, t1058 * t1146 + t1060 * t1142, t1010, -t1164, -t1236, 0, -pkin(5) * t1235 + t1142 * t963 + t1146 * t961, pkin(5) * t1011 + t1142 * t964 + t1146 * t962, t1142 * t936 + t1146 * t935 + pkin(5) * (-t1057 * t1142 + t1059 * t1146), t1142 * t852 + t1146 * t839 + pkin(5) * (-t1142 * t932 + t1146 * t933), t1142 * t907 + t1146 * t905, t1142 * t889 + t1146 * t888, t1142 * t929 + t1146 * t927, t1142 * t906 + t1146 * t904, t1142 * t928 + t1146 * t926, t1142 * t974 + t1146 * t973, t1142 * t789 + t1146 * t781 + pkin(5) * (-t1142 * t886 + t1146 * t887), t1142 * t795 + t1146 * t790 + pkin(5) * (-t1142 * t897 + t1146 * t898), t1142 * t772 + t1146 * t766 + pkin(5) * (-t1142 * t883 + t1146 * t884), t1142 * t760 + t1146 * t758 + pkin(5) * (-t1142 * t825 + t1146 * t826), t1142 * t805 + t1146 * t803, t1142 * t771 + t1146 * t770, t1142 * t811 + t1146 * t809, t1142 * t804 + t1146 * t802, t1142 * t812 + t1146 * t810, t1142 * t848 + t1146 * t847, t1142 * t739 + t1146 * t737 + pkin(5) * (-t1142 * t799 + t1146 * t800), t1142 * t741 + t1146 * t740 + pkin(5) * (-t1142 * t807 + t1146 * t808), t1142 * t733 + t1146 * t732 + pkin(5) * (-t1142 * t767 + t1146 * t768), t1142 * t730 + t1146 * t729 + pkin(5) * (-t1142 * t745 + t1146 * t746); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1119, t1120, 0, 0, 0, 0, 0, 0, 0, t1135, t1034, t1035, 0, -t1009, t1125, t1112, 0, t1127, 0, 0, t976, t975, t903, t870, t956, t922, t950, t955, t949, t979, t806, t824, t783, t762, t817, t777, t833, t816, t834, t849, t742, t749, t735, t731; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1148, 0, 0, -g(3), -t1119, 0, 0, 0, -t1108, 0, -t1105, 0, t1237, t1086, t1014, pkin(6) * t1014, t1073, t1060, t1067, -t1073, t1068, 0, t963, t964, t936, t852, t907, t889, t929, t906, t928, t974, t789, t795, t772, t760, t805, t771, t811, t804, t812, t848, t739, t741, t733, t730; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1148, 0, qJDD(1), 0, g(3), 0, -t1120, 0, 0, 0, t1105, 0, -t1108, 0, -t1086, t1237, t1174, t998, t1072, t1058, t1064, -t1072, -t1066, 0, t961, t962, t935, t839, t905, t888, t927, t904, t926, t973, t781, t790, t766, t758, t803, t770, t809, t802, t810, t847, t737, t740, t732, t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1119, t1120, 0, 0, 0, 0, 0, 0, 0, t1135, t1034, t1035, 0, -t1009, t1125, t1112, 0, t1127, 0, 0, t976, t975, t903, t870, t956, t922, t950, t955, t949, t979, t806, t824, t783, t762, t817, t777, t833, t816, t834, t849, t742, t749, t735, t731; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1135, 0, -t1193, 0, 0, -g(3), t1070, 0, t1186, t1101, t1094, -t1186, t1095, 0, t1039, t1040, t971, qJ(3) * t971, t958, t924, t953, t957, t952, t980, t838, t844, t794, t787, t819, t779, t835, t818, t836, t850, t757, t759, t736, t734; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1193, 0, t1135, 0, g(3), 0, t1071, 0, t1115, -t1104, -t1126, -t1115, -t1128, 0, t1023, t1024, 0, pkin(2) * t971, -t1056, -t1054, -t1088, t1056, t1087, -qJDD(4), t851, t868, t880, t801, -t986, -t985, -t920, t986, t916, -t1134, t765, t769, t761, t738; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1135, -t1070, -t1071, 0, 0, t1125, t1112, 0, t1127, 0, 0, t1185, t1160, t1168, t1212, t956, t922, t950, t955, t949, t979, t1171, t1170, t1172, t1156, t817, t777, t833, t816, t834, t849, t1177, t1176, t1178, t1179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1126, t1128, t1115, 0, t1129, 0, 0, t1046, t1023, 0, t1008, t992, t1003, t1006, t1002, t1031, t937, t960, t854, -pkin(7) * t874, t860, t822, t878, t859, t879, t891, t788, t793, t756, t748; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1126, -t1180, t1128, -t1115, 0, -t1046, 0, t1024, 0, t1007, t990, t1000, t1005, t999, t1030, t893, t925, t846, t862, t858, t820, t876, t857, t877, t890, t776, t785, t754, t744; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1115, t1104, t1126, t1115, t1128, 0, -t1023, -t1024, 0, 0, t1056, t1054, t1088, -t1056, -t1087, qJDD(4), t1159, t1169, t1226, t1227, t986, t985, t920, -t986, -t916, t1134, t1155, t1153, t1191, t1192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1053, -t1051, t1229, t1081, t1075, -t1081, 0, t1017, t942, 0, t902, t866, t946, t900, t947, t978, t869, t871, t782, -pkin(8) * t797; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1203, t1052, t1076, -t1154, t1049, -t1203, -t1017, 0, t943, 0, t901, t864, t944, t899, t945, t977, t837, t840, t773, t791; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1056, t1054, t1088, -t1056, -t1087, qJDD(4), -t942, -t943, 0, 0, t986, t985, t920, -t986, -t916, t1134, t1167, t1157, t863, t796; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t968, -t915, t1232, t1029, t1026, -t1029, 0, t931, t842, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1200, t1233, t1027, t967, t983, -t1200, -t931, 0, t843, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t986, t985, t920, -t986, -t916, t1134, -t842, -t843, 0, 0;];
m_new_reg = t1;
