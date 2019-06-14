% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-06 00:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:43:45
% EndTime: 2019-05-06 00:44:39
% DurationCPUTime: 53.32s
% Computational Cost: add. (754109->389), mult. (2352425->513), div. (0->0), fcn. (1993669->14), ass. (0->182)
t1157 = Ifges(5,1) + Ifges(6,2);
t1150 = Ifges(5,4) + Ifges(6,6);
t1149 = Ifges(5,5) - Ifges(6,4);
t1156 = -Ifges(5,2) - Ifges(6,3);
t1148 = Ifges(5,6) - Ifges(6,5);
t1155 = Ifges(5,3) + Ifges(6,1);
t1081 = sin(pkin(12));
t1083 = sin(pkin(6));
t1084 = cos(pkin(12));
t1086 = cos(pkin(6));
t1092 = cos(qJ(3));
t1085 = cos(pkin(7));
t1089 = sin(qJ(3));
t1127 = t1085 * t1089;
t1082 = sin(pkin(7));
t1132 = t1082 * t1089;
t1098 = t1086 * t1132 + (t1081 * t1092 + t1084 * t1127) * t1083;
t1051 = t1098 * qJD(1);
t1126 = t1085 * t1092;
t1123 = t1084 * t1126;
t1131 = t1082 * t1092;
t1124 = t1086 * t1131;
t1037 = -t1051 * qJD(3) + (t1083 * (-t1081 * t1089 + t1123) + t1124) * qJDD(1);
t1129 = t1083 * t1085;
t1064 = (t1082 * t1086 + t1084 * t1129) * qJD(1) * pkin(9);
t1090 = sin(qJ(1));
t1093 = cos(qJ(1));
t1078 = -g(1) * t1093 - g(2) * t1090;
t1094 = qJD(1) ^ 2;
t1139 = qJ(2) * t1083;
t1068 = -pkin(1) * t1094 + qJDD(1) * t1139 + t1078;
t1146 = pkin(9) * t1081;
t1114 = -pkin(2) * t1084 - t1082 * t1146;
t1137 = qJD(1) * t1083;
t1144 = pkin(9) * qJDD(1);
t1111 = qJD(1) * t1114 * t1137 + t1085 * t1144;
t1077 = t1090 * g(1) - g(2) * t1093;
t1067 = qJDD(1) * pkin(1) + t1094 * t1139 + t1077;
t1134 = t1081 * t1083;
t1125 = qJD(1) * t1134;
t1128 = t1084 * t1086;
t1130 = t1083 * t1084;
t1115 = -g(3) * t1130 - 0.2e1 * qJD(2) * t1125 + t1067 * t1128;
t1013 = (pkin(2) * qJDD(1) + qJD(1) * t1064) * t1086 + (-t1083 * t1111 - t1068) * t1081 + t1115;
t1069 = (pkin(2) * t1086 - t1129 * t1146) * qJD(1);
t1133 = t1081 * t1086;
t1122 = 0.2e1 * qJD(2) * qJD(1) * t1130 + t1067 * t1133 + t1084 * t1068;
t1014 = (-qJD(1) * t1069 + t1082 * t1144) * t1086 + (-g(3) * t1081 + t1084 * t1111) * t1083 + t1122;
t1121 = -t1086 * g(3) + qJDD(2);
t1023 = (-t1067 + t1114 * qJDD(1) + (-t1064 * t1084 + t1069 * t1081) * qJD(1)) * t1083 + t1121;
t976 = -t1089 * t1014 + (t1013 * t1085 + t1023 * t1082) * t1092;
t1050 = qJD(1) * t1124 - t1089 * t1125 + t1123 * t1137;
t1038 = t1050 * qJD(3) + qJDD(1) * t1098;
t1108 = -t1082 * t1130 + t1085 * t1086;
t1065 = qJD(1) * t1108 + qJD(3);
t1088 = sin(qJ(4));
t1151 = cos(qJ(4));
t1045 = t1051 * t1151 + t1088 * t1065;
t1062 = qJDD(1) * t1108 + qJDD(3);
t1007 = t1045 * qJD(4) + t1088 * t1038 - t1062 * t1151;
t1044 = t1088 * t1051 - t1065 * t1151;
t1008 = -t1044 * qJD(4) + t1038 * t1151 + t1088 * t1062;
t1017 = -mrSges(6,2) * t1044 - mrSges(6,3) * t1045;
t1034 = qJDD(4) - t1037;
t1087 = sin(qJ(6));
t1091 = cos(qJ(6));
t1049 = qJD(4) - t1050;
t1025 = mrSges(6,1) * t1045 + mrSges(6,2) * t1049;
t1021 = t1044 * t1091 - t1049 * t1087;
t1022 = t1044 * t1087 + t1049 * t1091;
t1028 = pkin(5) * t1045 - pkin(11) * t1049;
t1043 = t1044 ^ 2;
t1015 = pkin(4) * t1044 - qJ(5) * t1045;
t1048 = t1049 ^ 2;
t1036 = -pkin(3) * t1050 - pkin(10) * t1051;
t1061 = t1065 ^ 2;
t977 = t1013 * t1127 + t1092 * t1014 + t1023 * t1132;
t973 = -pkin(3) * t1061 + pkin(10) * t1062 + t1036 * t1050 + t977;
t988 = -t1082 * t1013 + t1085 * t1023;
t975 = (-t1050 * t1065 - t1038) * pkin(10) + (t1051 * t1065 - t1037) * pkin(3) + t988;
t969 = t1088 * t975 + t1151 * t973;
t1101 = -t1048 * pkin(4) + t1034 * qJ(5) - t1044 * t1015 + t969;
t963 = -t1007 * pkin(5) - t1043 * pkin(11) + ((2 * qJD(5)) + t1028) * t1049 + t1101;
t982 = -qJD(6) * t1022 + t1007 * t1091 - t1034 * t1087;
t983 = qJD(6) * t1021 + t1007 * t1087 + t1034 * t1091;
t1040 = qJD(6) + t1045;
t992 = -mrSges(7,2) * t1040 + mrSges(7,3) * t1021;
t993 = mrSges(7,1) * t1040 - mrSges(7,3) * t1022;
t1106 = -m(7) * t963 + t982 * mrSges(7,1) - t983 * mrSges(7,2) + t1021 * t992 - t1022 * t993;
t1152 = -2 * qJD(5);
t965 = t1049 * t1152 - t1101;
t1099 = -m(6) * t965 + t1034 * mrSges(6,3) + t1049 * t1025 - t1106;
t1138 = t1150 * t1044 - t1045 * t1157 - t1149 * t1049;
t1140 = t1044 * t1156 + t1045 * t1150 + t1049 * t1148;
t1024 = mrSges(6,1) * t1044 - mrSges(6,3) * t1049;
t1004 = qJDD(6) + t1008;
t1136 = t1044 * t1049;
t968 = -t1088 * t973 + t1151 * t975;
t966 = -t1034 * pkin(4) - t1048 * qJ(5) + t1045 * t1015 + qJDD(5) - t968;
t961 = (t1044 * t1045 - t1034) * pkin(11) + (t1008 + t1136) * pkin(5) + t966;
t972 = -t1062 * pkin(3) - t1061 * pkin(10) + t1051 * t1036 - t976;
t1095 = (-t1008 + t1136) * qJ(5) + t972 + (t1049 * pkin(4) + t1152) * t1045;
t964 = -t1043 * pkin(5) - t1045 * t1028 + (pkin(4) + pkin(11)) * t1007 + t1095;
t959 = -t1087 * t964 + t1091 * t961;
t989 = -mrSges(7,1) * t1021 + mrSges(7,2) * t1022;
t956 = m(7) * t959 + mrSges(7,1) * t1004 - mrSges(7,3) * t983 - t1022 * t989 + t1040 * t992;
t960 = t1087 * t961 + t1091 * t964;
t957 = m(7) * t960 - mrSges(7,2) * t1004 + mrSges(7,3) * t982 + t1021 * t989 - t1040 * t993;
t947 = t1087 * t957 + t1091 * t956;
t1105 = -m(6) * t966 - t1008 * mrSges(6,1) - t1045 * t1017 - t947;
t945 = t1034 * mrSges(6,2) + t1049 * t1024 - t1105;
t984 = Ifges(7,5) * t1022 + Ifges(7,6) * t1021 + Ifges(7,3) * t1040;
t986 = Ifges(7,1) * t1022 + Ifges(7,4) * t1021 + Ifges(7,5) * t1040;
t948 = -mrSges(7,1) * t963 + mrSges(7,3) * t960 + Ifges(7,4) * t983 + Ifges(7,2) * t982 + Ifges(7,6) * t1004 - t1022 * t984 + t1040 * t986;
t985 = Ifges(7,4) * t1022 + Ifges(7,2) * t1021 + Ifges(7,6) * t1040;
t949 = mrSges(7,2) * t963 - mrSges(7,3) * t959 + Ifges(7,1) * t983 + Ifges(7,4) * t982 + Ifges(7,5) * t1004 + t1021 * t984 - t1040 * t985;
t1153 = -t1007 * t1148 + t1008 * t1149 + t1155 * t1034 - t1044 * t1138 + t1045 * t1140 + mrSges(5,1) * t968 - mrSges(5,2) * t969 + mrSges(6,2) * t966 - mrSges(6,3) * t965 - pkin(4) * t945 - pkin(11) * t947 + qJ(5) * (-t1007 * mrSges(6,1) - t1044 * t1017 + t1099) - t1087 * t948 + t1091 * t949;
t1145 = Ifges(3,3) * t1086;
t1041 = -t1081 * t1068 + t1115;
t1118 = -mrSges(3,1) * t1084 + mrSges(3,2) * t1081;
t1066 = t1118 * t1137;
t1112 = -mrSges(3,2) * t1086 + mrSges(3,3) * t1130;
t1071 = t1112 * qJD(1);
t1113 = mrSges(3,1) * t1086 - mrSges(3,3) * t1134;
t1035 = -mrSges(4,1) * t1050 + mrSges(4,2) * t1051;
t1047 = mrSges(4,1) * t1065 - mrSges(4,3) * t1051;
t1016 = mrSges(5,1) * t1044 + mrSges(5,2) * t1045;
t1026 = -mrSges(5,2) * t1049 - mrSges(5,3) * t1044;
t944 = m(5) * t968 - t1008 * mrSges(5,3) - t1045 * t1016 + (-t1024 + t1026) * t1049 + (mrSges(5,1) - mrSges(6,2)) * t1034 + t1105;
t1027 = mrSges(5,1) * t1049 - mrSges(5,3) * t1045;
t952 = m(5) * t969 - t1034 * mrSges(5,2) - t1049 * t1027 + (-t1016 - t1017) * t1044 + (-mrSges(5,3) - mrSges(6,1)) * t1007 + t1099;
t1120 = -t1088 * t944 + t1151 * t952;
t936 = m(4) * t977 - mrSges(4,2) * t1062 + mrSges(4,3) * t1037 + t1035 * t1050 - t1047 * t1065 + t1120;
t1046 = -mrSges(4,2) * t1065 + mrSges(4,3) * t1050;
t939 = t1088 * t952 + t1151 * t944;
t938 = m(4) * t988 - mrSges(4,1) * t1037 + mrSges(4,2) * t1038 - t1046 * t1050 + t1047 * t1051 + t939;
t1142 = -t1087 * t956 + t1091 * t957;
t967 = t1007 * pkin(4) + t1095;
t1110 = -m(6) * t967 + t1007 * mrSges(6,2) + t1044 * t1024 - t1142;
t1097 = -m(5) * t972 - t1007 * mrSges(5,1) - t1044 * t1026 + (t1025 - t1027) * t1045 + (-mrSges(5,2) + mrSges(6,3)) * t1008 + t1110;
t942 = m(4) * t976 + t1062 * mrSges(4,1) - t1038 * mrSges(4,3) - t1051 * t1035 + t1065 * t1046 + t1097;
t925 = -t1082 * t938 + t942 * t1126 + t936 * t1127;
t921 = m(3) * t1041 + t1113 * qJDD(1) + (-t1066 * t1134 + t1071 * t1086) * qJD(1) + t925;
t1052 = -t1083 * t1067 + t1121;
t1070 = t1113 * qJD(1);
t924 = t1085 * t938 + t942 * t1131 + t936 * t1132;
t923 = m(3) * t1052 + (t1118 * qJDD(1) + (t1070 * t1081 - t1071 * t1084) * qJD(1)) * t1083 + t924;
t1042 = -g(3) * t1134 + t1122;
t932 = -t1089 * t942 + t1092 * t936;
t931 = m(3) * t1042 + t1112 * qJDD(1) + (t1066 * t1130 - t1070 * t1086) * qJD(1) + t932;
t910 = -t1083 * t923 + t921 * t1128 + t931 * t1133;
t907 = m(2) * t1077 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1094 + t910;
t917 = -t1081 * t921 + t1084 * t931;
t915 = m(2) * t1078 - mrSges(2,1) * t1094 - qJDD(1) * mrSges(2,2) + t917;
t1143 = t1090 * t915 + t1093 * t907;
t1141 = t1044 * t1148 - t1045 * t1149 - t1049 * t1155;
t909 = t1086 * t923 + t921 * t1130 + t931 * t1134;
t1119 = -t1090 * t907 + t1093 * t915;
t1117 = Ifges(3,5) * t1081 + Ifges(3,6) * t1084;
t1030 = Ifges(4,5) * t1051 + Ifges(4,6) * t1050 + Ifges(4,3) * t1065;
t1031 = Ifges(4,4) * t1051 + Ifges(4,2) * t1050 + Ifges(4,6) * t1065;
t946 = -t1008 * mrSges(6,3) - t1045 * t1025 - t1110;
t926 = -mrSges(5,1) * t972 - mrSges(6,1) * t965 + mrSges(6,2) * t967 + mrSges(5,3) * t969 - pkin(4) * t946 - pkin(5) * t1106 - pkin(11) * t1142 + t1007 * t1156 + t1150 * t1008 + t1148 * t1034 + t1141 * t1045 - t1138 * t1049 - t1087 * t949 - t1091 * t948;
t1100 = mrSges(7,1) * t959 - mrSges(7,2) * t960 + Ifges(7,5) * t983 + Ifges(7,6) * t982 + Ifges(7,3) * t1004 - t1021 * t986 + t1022 * t985;
t927 = mrSges(6,1) * t966 + mrSges(5,2) * t972 - mrSges(5,3) * t968 - mrSges(6,3) * t967 + pkin(5) * t947 - qJ(5) * t946 - t1150 * t1007 + t1008 * t1157 + t1149 * t1034 + t1141 * t1044 - t1140 * t1049 + t1100;
t912 = mrSges(4,2) * t988 - mrSges(4,3) * t976 + Ifges(4,1) * t1038 + Ifges(4,4) * t1037 + Ifges(4,5) * t1062 - pkin(10) * t939 + t1050 * t1030 - t1065 * t1031 - t1088 * t926 + t1151 * t927;
t1032 = Ifges(4,1) * t1051 + Ifges(4,4) * t1050 + Ifges(4,5) * t1065;
t918 = -mrSges(4,1) * t988 + mrSges(4,3) * t977 + Ifges(4,4) * t1038 + Ifges(4,2) * t1037 + Ifges(4,6) * t1062 - pkin(3) * t939 - t1051 * t1030 + t1065 * t1032 - t1153;
t1107 = pkin(9) * t932 + t1089 * t912 + t1092 * t918;
t1102 = Ifges(3,6) * t1086 + (Ifges(3,4) * t1081 + Ifges(3,2) * t1084) * t1083;
t1056 = t1102 * qJD(1);
t1103 = Ifges(3,5) * t1086 + (Ifges(3,1) * t1081 + Ifges(3,4) * t1084) * t1083;
t1057 = t1103 * qJD(1);
t911 = mrSges(4,1) * t976 - mrSges(4,2) * t977 + Ifges(4,5) * t1038 + Ifges(4,6) * t1037 + Ifges(4,3) * t1062 + pkin(3) * t1097 + pkin(10) * t1120 + t1051 * t1031 - t1050 * t1032 + t1088 * t927 + t1151 * t926;
t901 = qJDD(1) * t1145 + mrSges(3,1) * t1041 - mrSges(3,2) * t1042 + pkin(2) * t925 + t1085 * t911 + t1107 * t1082 + (t1117 * qJDD(1) + (t1056 * t1081 - t1057 * t1084) * qJD(1)) * t1083;
t1055 = (t1083 * t1117 + t1145) * qJD(1);
t903 = -mrSges(3,1) * t1052 + mrSges(3,3) * t1042 - pkin(2) * t924 - t1082 * t911 + (-t1055 * t1134 + t1057 * t1086) * qJD(1) + t1107 * t1085 + t1102 * qJDD(1);
t905 = mrSges(3,2) * t1052 - mrSges(3,3) * t1041 - t1089 * t918 + t1092 * t912 + (t1055 * t1130 - t1056 * t1086) * qJD(1) + (-t1082 * t924 - t1085 * t925) * pkin(9) + t1103 * qJDD(1);
t1104 = mrSges(2,1) * t1077 - mrSges(2,2) * t1078 + Ifges(2,3) * qJDD(1) + pkin(1) * t910 + t1086 * t901 + t903 * t1130 + t905 * t1134 + t917 * t1139;
t899 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1077 + Ifges(2,5) * qJDD(1) - t1094 * Ifges(2,6) - t1081 * t903 + t1084 * t905 + (-t1083 * t909 - t1086 * t910) * qJ(2);
t898 = mrSges(2,1) * g(3) + mrSges(2,3) * t1078 + t1094 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t909 - t1083 * t901 + (qJ(2) * t917 + t1081 * t905 + t1084 * t903) * t1086;
t1 = [-m(1) * g(1) + t1119; -m(1) * g(2) + t1143; (-m(1) - m(2)) * g(3) + t909; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1143 - t1090 * t898 + t1093 * t899; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1119 + t1090 * t899 + t1093 * t898; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1104; t1104; t923; t911; t1153; t945; t1100;];
tauJB  = t1;
