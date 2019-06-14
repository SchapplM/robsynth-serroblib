% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR13
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR13_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:55:35
% EndTime: 2019-05-05 20:56:23
% DurationCPUTime: 47.15s
% Computational Cost: add. (668696->385), mult. (2118812->508), div. (0->0), fcn. (1780803->14), ass. (0->182)
t1119 = cos(pkin(12));
t1117 = sin(pkin(7));
t1183 = cos(pkin(6));
t1158 = t1183 * t1117;
t1118 = sin(pkin(6));
t1182 = cos(pkin(7));
t1161 = t1118 * t1182;
t1195 = t1119 * t1161 + t1158;
t1194 = Ifges(4,1) + Ifges(5,2);
t1187 = Ifges(4,5) - Ifges(5,4);
t1193 = -Ifges(4,2) - Ifges(5,3);
t1186 = Ifges(4,6) - Ifges(5,5);
t1185 = -Ifges(5,6) - Ifges(4,4);
t1192 = Ifges(4,3) + Ifges(5,1);
t1095 = t1195 * qJD(1) * pkin(9);
t1123 = sin(qJ(1));
t1126 = cos(qJ(1));
t1112 = -g(1) * t1126 - g(2) * t1123;
t1127 = qJD(1) ^ 2;
t1180 = qJ(2) * t1118;
t1100 = -pkin(1) * t1127 + qJDD(1) * t1180 + t1112;
t1116 = sin(pkin(12));
t1184 = pkin(9) * t1117;
t1145 = -pkin(2) * t1119 - t1116 * t1184;
t1169 = t1182 * pkin(9);
t1179 = qJD(1) * t1118;
t1140 = qJD(1) * t1145 * t1179 + qJDD(1) * t1169;
t1111 = t1123 * g(1) - g(2) * t1126;
t1099 = qJDD(1) * pkin(1) + t1127 * t1180 + t1111;
t1160 = t1119 * t1183;
t1167 = qJD(2) * t1179;
t1174 = t1118 * t1119;
t1146 = -g(3) * t1174 + t1099 * t1160 - 0.2e1 * t1116 * t1167;
t1157 = qJDD(1) * t1183;
t1163 = qJD(1) * t1183;
t1033 = pkin(2) * t1157 + t1095 * t1163 + (-t1118 * t1140 - t1100) * t1116 + t1146;
t1101 = (-pkin(9) * t1116 * t1161 + pkin(2) * t1183) * qJD(1);
t1162 = t1116 * t1183;
t1164 = t1099 * t1162 + (t1100 + 0.2e1 * t1167) * t1119;
t1034 = t1157 * t1184 - t1101 * t1163 + (-g(3) * t1116 + t1119 * t1140) * t1118 + t1164;
t1150 = -g(3) * t1183 + qJDD(2);
t1044 = (-t1099 + t1145 * qJDD(1) + (-t1095 * t1119 + t1101 * t1116) * qJD(1)) * t1118 + t1150;
t1122 = sin(qJ(3));
t1159 = t1122 * t1182;
t1175 = t1117 * t1122;
t1189 = cos(qJ(3));
t1007 = t1033 * t1159 + t1034 * t1189 + t1044 * t1175;
t1176 = t1116 * t1118;
t1190 = t1122 * t1176 - t1189 * t1195;
t1080 = t1190 * qJD(1);
t1129 = t1122 * t1158 + (t1116 * t1189 + t1119 * t1159) * t1118;
t1081 = t1129 * qJD(1);
t1060 = pkin(3) * t1080 - qJ(4) * t1081;
t1191 = -t1117 * t1174 + t1183 * t1182;
t1096 = -qJD(1) * t1191 - qJD(3);
t1092 = t1096 ^ 2;
t1093 = qJDD(1) * t1191 + qJDD(3);
t1002 = pkin(3) * t1092 - t1093 * qJ(4) + 0.2e1 * qJD(4) * t1096 + t1080 * t1060 - t1007;
t1188 = mrSges(4,1) - mrSges(5,2);
t1066 = -t1116 * t1100 + t1146;
t1149 = -mrSges(3,1) * t1119 + mrSges(3,2) * t1116;
t1098 = t1149 * t1179;
t1141 = -mrSges(3,2) * t1183 + mrSges(3,3) * t1174;
t1103 = t1141 * qJD(1);
t1142 = mrSges(3,1) * t1183 - mrSges(3,3) * t1176;
t1153 = t1182 * t1189;
t1168 = t1117 * t1189;
t1006 = t1033 * t1153 - t1122 * t1034 + t1044 * t1168;
t1061 = mrSges(4,1) * t1080 + mrSges(4,2) * t1081;
t1064 = -t1080 * qJD(3) + qJDD(1) * t1129;
t1004 = -t1093 * pkin(3) - t1092 * qJ(4) + t1081 * t1060 + qJDD(4) - t1006;
t1062 = -mrSges(5,2) * t1080 - mrSges(5,3) * t1081;
t1121 = sin(qJ(5));
t1125 = cos(qJ(5));
t1063 = qJD(3) * t1081 + qJDD(1) * t1190;
t1069 = t1080 * t1121 - t1096 * t1125;
t1028 = -qJD(5) * t1069 + t1063 * t1125 - t1093 * t1121;
t1068 = t1080 * t1125 + t1096 * t1121;
t1035 = -mrSges(6,1) * t1068 + mrSges(6,2) * t1069;
t1078 = qJD(5) + t1081;
t1046 = mrSges(6,1) * t1078 - mrSges(6,3) * t1069;
t1059 = qJDD(5) + t1064;
t1120 = sin(qJ(6));
t1124 = cos(qJ(6));
t1029 = qJD(5) * t1068 + t1063 * t1121 + t1093 * t1125;
t1042 = -t1069 * t1120 + t1078 * t1124;
t1010 = qJD(6) * t1042 + t1029 * t1124 + t1059 * t1120;
t1043 = t1069 * t1124 + t1078 * t1120;
t1017 = -mrSges(7,1) * t1042 + mrSges(7,2) * t1043;
t1065 = qJD(6) - t1068;
t1018 = -mrSges(7,2) * t1065 + mrSges(7,3) * t1042;
t1026 = qJDD(6) - t1028;
t1036 = -pkin(5) * t1068 - pkin(11) * t1069;
t1077 = t1078 ^ 2;
t1074 = pkin(4) * t1081 + pkin(10) * t1096;
t1079 = t1080 ^ 2;
t1015 = -t1033 * t1117 + t1182 * t1044;
t1177 = t1080 * t1096;
t1132 = (-t1064 - t1177) * qJ(4) + t1015 + (-pkin(3) * t1096 - 0.2e1 * qJD(4)) * t1081;
t1001 = -pkin(4) * t1079 - t1074 * t1081 + (pkin(3) + pkin(10)) * t1063 + t1132;
t997 = (t1080 * t1081 - t1093) * pkin(10) + (t1064 - t1177) * pkin(4) + t1004;
t994 = t1125 * t1001 + t1121 * t997;
t991 = -pkin(5) * t1077 + pkin(11) * t1059 + t1036 * t1068 + t994;
t999 = -pkin(4) * t1063 - pkin(10) * t1079 - t1096 * t1074 - t1002;
t995 = (-t1068 * t1078 - t1029) * pkin(11) + (t1069 * t1078 - t1028) * pkin(5) + t999;
t987 = -t1120 * t991 + t1124 * t995;
t985 = m(7) * t987 + mrSges(7,1) * t1026 - mrSges(7,3) * t1010 - t1017 * t1043 + t1018 * t1065;
t1009 = -qJD(6) * t1043 - t1029 * t1120 + t1059 * t1124;
t1019 = mrSges(7,1) * t1065 - mrSges(7,3) * t1043;
t988 = t1120 * t995 + t1124 * t991;
t986 = m(7) * t988 - mrSges(7,2) * t1026 + mrSges(7,3) * t1009 + t1017 * t1042 - t1019 * t1065;
t1156 = -t1120 * t985 + t1124 * t986;
t974 = m(6) * t994 - mrSges(6,2) * t1059 + mrSges(6,3) * t1028 + t1035 * t1068 - t1046 * t1078 + t1156;
t1045 = -mrSges(6,2) * t1078 + mrSges(6,3) * t1068;
t993 = -t1001 * t1121 + t1125 * t997;
t990 = -pkin(5) * t1059 - pkin(11) * t1077 + t1036 * t1069 - t993;
t1137 = -m(7) * t990 + t1009 * mrSges(7,1) - mrSges(7,2) * t1010 + t1042 * t1018 - t1019 * t1043;
t981 = m(6) * t993 + mrSges(6,1) * t1059 - mrSges(6,3) * t1029 - t1035 * t1069 + t1045 * t1078 + t1137;
t968 = t1121 * t974 + t1125 * t981;
t1136 = -m(5) * t1004 - t1064 * mrSges(5,1) - t1081 * t1062 - t968;
t1071 = mrSges(5,1) * t1080 + mrSges(5,3) * t1096;
t1170 = mrSges(4,2) * t1096 - mrSges(4,3) * t1080 - t1071;
t963 = m(4) * t1006 - mrSges(4,3) * t1064 - t1061 * t1081 + t1093 * t1188 - t1096 * t1170 + t1136;
t1073 = -mrSges(4,1) * t1096 - mrSges(4,3) * t1081;
t1005 = pkin(3) * t1063 + t1132;
t1072 = mrSges(5,1) * t1081 - mrSges(5,2) * t1096;
t1155 = -t1121 * t981 + t1125 * t974;
t1143 = m(5) * t1005 - t1064 * mrSges(5,3) - t1081 * t1072 + t1155;
t965 = m(4) * t1015 + mrSges(4,2) * t1064 + t1063 * t1188 + t1073 * t1081 + t1080 * t1170 + t1143;
t976 = t1120 * t986 + t1124 * t985;
t1135 = -m(6) * t999 + mrSges(6,1) * t1028 - t1029 * mrSges(6,2) + t1045 * t1068 - t1069 * t1046 - t976;
t1130 = -m(5) * t1002 + t1093 * mrSges(5,3) - t1096 * t1072 - t1135;
t972 = t1130 + (-mrSges(4,3) - mrSges(5,1)) * t1063 + (-t1061 - t1062) * t1080 + m(4) * t1007 - mrSges(4,2) * t1093 + t1073 * t1096;
t953 = -t1117 * t965 + t963 * t1153 + t972 * t1159;
t949 = m(3) * t1066 + t1142 * qJDD(1) + (-t1098 * t1176 + t1103 * t1183) * qJD(1) + t953;
t1082 = -t1118 * t1099 + t1150;
t1102 = t1142 * qJD(1);
t952 = t963 * t1168 + t972 * t1175 + t1182 * t965;
t951 = m(3) * t1082 + (t1149 * qJDD(1) + (t1102 * t1116 - t1103 * t1119) * qJD(1)) * t1118 + t952;
t1067 = -g(3) * t1176 + t1164;
t958 = -t1122 * t963 + t1189 * t972;
t957 = m(3) * t1067 + t1141 * qJDD(1) + (t1098 * t1174 - t1102 * t1183) * qJD(1) + t958;
t938 = -t1118 * t951 + t949 * t1160 + t957 * t1162;
t935 = m(2) * t1111 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1127 + t938;
t945 = -t1116 * t949 + t1119 * t957;
t943 = m(2) * t1112 - mrSges(2,1) * t1127 - qJDD(1) * mrSges(2,2) + t945;
t1181 = t1123 * t943 + t1126 * t935;
t1173 = t1080 * t1186 - t1081 * t1187 + t1096 * t1192;
t1172 = t1080 * t1193 - t1081 * t1185 - t1096 * t1186;
t1171 = -t1185 * t1080 - t1081 * t1194 + t1187 * t1096;
t937 = t949 * t1174 + t957 * t1176 + t1183 * t951;
t1154 = -t1123 * t935 + t1126 * t943;
t1148 = Ifges(3,5) * t1116 + Ifges(3,6) * t1119;
t1133 = t1183 * Ifges(3,6) + (Ifges(3,4) * t1116 + Ifges(3,2) * t1119) * t1118;
t1086 = t1133 * qJD(1);
t1134 = t1183 * Ifges(3,5) + (Ifges(3,1) * t1116 + Ifges(3,4) * t1119) * t1118;
t1087 = t1134 * qJD(1);
t1020 = Ifges(6,5) * t1069 + Ifges(6,6) * t1068 + Ifges(6,3) * t1078;
t1021 = Ifges(6,4) * t1069 + Ifges(6,2) * t1068 + Ifges(6,6) * t1078;
t1011 = Ifges(7,5) * t1043 + Ifges(7,6) * t1042 + Ifges(7,3) * t1065;
t1013 = Ifges(7,1) * t1043 + Ifges(7,4) * t1042 + Ifges(7,5) * t1065;
t979 = -mrSges(7,1) * t990 + mrSges(7,3) * t988 + Ifges(7,4) * t1010 + Ifges(7,2) * t1009 + Ifges(7,6) * t1026 - t1011 * t1043 + t1013 * t1065;
t1012 = Ifges(7,4) * t1043 + Ifges(7,2) * t1042 + Ifges(7,6) * t1065;
t980 = mrSges(7,2) * t990 - mrSges(7,3) * t987 + Ifges(7,1) * t1010 + Ifges(7,4) * t1009 + Ifges(7,5) * t1026 + t1011 * t1042 - t1012 * t1065;
t959 = mrSges(6,2) * t999 - mrSges(6,3) * t993 + Ifges(6,1) * t1029 + Ifges(6,4) * t1028 + Ifges(6,5) * t1059 - pkin(11) * t976 + t1020 * t1068 - t1021 * t1078 - t1120 * t979 + t1124 * t980;
t1022 = Ifges(6,1) * t1069 + Ifges(6,4) * t1068 + Ifges(6,5) * t1078;
t1128 = mrSges(7,1) * t987 - mrSges(7,2) * t988 + Ifges(7,5) * t1010 + Ifges(7,6) * t1009 + Ifges(7,3) * t1026 + t1012 * t1043 - t1013 * t1042;
t960 = -mrSges(6,1) * t999 + mrSges(6,3) * t994 + Ifges(6,4) * t1029 + Ifges(6,2) * t1028 + Ifges(6,6) * t1059 - pkin(5) * t976 - t1020 * t1069 + t1022 * t1078 - t1128;
t967 = mrSges(5,2) * t1093 - t1071 * t1096 - t1136;
t939 = mrSges(4,1) * t1006 - mrSges(4,2) * t1007 + mrSges(5,2) * t1004 - mrSges(5,3) * t1002 + t1125 * t959 - t1121 * t960 - pkin(10) * t968 - pkin(3) * t967 + qJ(4) * t1130 + t1192 * t1093 + t1172 * t1081 + (-qJ(4) * t1062 - t1171) * t1080 + t1187 * t1064 + (-qJ(4) * mrSges(5,1) - t1186) * t1063;
t966 = -mrSges(5,2) * t1063 - t1071 * t1080 + t1143;
t940 = -mrSges(4,1) * t1015 - mrSges(5,1) * t1002 + mrSges(5,2) * t1005 + mrSges(4,3) * t1007 - pkin(3) * t966 - pkin(4) * t1135 - pkin(10) * t1155 + t1063 * t1193 - t1185 * t1064 + t1173 * t1081 + t1186 * t1093 + t1171 * t1096 - t1121 * t959 - t1125 * t960;
t1131 = mrSges(6,1) * t993 - mrSges(6,2) * t994 + Ifges(6,5) * t1029 + Ifges(6,6) * t1028 + Ifges(6,3) * t1059 + pkin(5) * t1137 + pkin(11) * t1156 + t1069 * t1021 - t1068 * t1022 + t1120 * t980 + t1124 * t979;
t946 = mrSges(5,1) * t1004 + mrSges(4,2) * t1015 - mrSges(4,3) * t1006 - mrSges(5,3) * t1005 + pkin(4) * t968 - qJ(4) * t966 + t1185 * t1063 + t1064 * t1194 + t1173 * t1080 + t1187 * t1093 + t1172 * t1096 + t1131;
t929 = Ifges(3,3) * t1157 + mrSges(3,1) * t1066 - mrSges(3,2) * t1067 + t1182 * t939 + pkin(2) * t953 + (pkin(9) * t958 + t1122 * t946 + t1189 * t940) * t1117 + (t1148 * qJDD(1) + (t1086 * t1116 - t1087 * t1119) * qJD(1)) * t1118;
t1085 = (Ifges(3,3) * t1183 + t1118 * t1148) * qJD(1);
t931 = t940 * t1153 + t958 * t1169 + t946 * t1159 - mrSges(3,1) * t1082 + mrSges(3,3) * t1067 - pkin(2) * t952 - t1117 * t939 + (-t1085 * t1176 + t1087 * t1183) * qJD(1) + t1133 * qJDD(1);
t933 = mrSges(3,2) * t1082 - mrSges(3,3) * t1066 + t1189 * t946 - t1122 * t940 + (t1085 * t1174 - t1086 * t1183) * qJD(1) + (-t1117 * t952 - t1182 * t953) * pkin(9) + t1134 * qJDD(1);
t1138 = mrSges(2,1) * t1111 - mrSges(2,2) * t1112 + Ifges(2,3) * qJDD(1) + pkin(1) * t938 + t931 * t1174 + t933 * t1176 + t945 * t1180 + t1183 * t929;
t927 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1111 + Ifges(2,5) * qJDD(1) - t1127 * Ifges(2,6) - t1116 * t931 + t1119 * t933 + (-t1118 * t937 - t1183 * t938) * qJ(2);
t926 = qJ(2) * t1183 * t945 + mrSges(2,1) * g(3) + mrSges(2,3) * t1112 + t1127 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t937 - t1118 * t929 + t1160 * t931 + t1162 * t933;
t1 = [-m(1) * g(1) + t1154; -m(1) * g(2) + t1181; (-m(1) - m(2)) * g(3) + t937; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1181 - t1123 * t926 + t1126 * t927; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1154 + t1123 * t927 + t1126 * t926; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1138; t1138; t951; t939; t967; t1131; t1128;];
tauJB  = t1;
