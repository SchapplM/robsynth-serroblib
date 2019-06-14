% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 11:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:32:27
% EndTime: 2019-05-06 11:32:38
% DurationCPUTime: 7.87s
% Computational Cost: add. (107641->368), mult. (251291->448), div. (0->0), fcn. (167889->10), ass. (0->161)
t1100 = cos(pkin(6));
t1092 = qJD(1) * t1100 + qJD(2);
t1167 = -0.2e1 * t1092;
t1166 = Ifges(3,1) + Ifges(5,3) + Ifges(4,2);
t1140 = Ifges(3,4) + Ifges(4,6) - Ifges(5,6);
t1139 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t1165 = Ifges(3,2) + Ifges(4,3) + Ifges(5,2);
t1138 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t1164 = Ifges(3,3) + Ifges(4,1) + Ifges(5,1);
t1099 = sin(pkin(6));
t1103 = sin(qJ(2));
t1107 = cos(qJ(2));
t1068 = (qJD(1) * qJD(2) * t1107 + qJDD(1) * t1103) * t1099;
t1091 = qJDD(1) * t1100 + qJDD(2);
t1147 = t1099 * t1107;
t1135 = qJD(1) * t1147;
t1126 = t1092 * t1135;
t1109 = qJD(1) ^ 2;
t1150 = t1099 ^ 2 * t1109;
t1133 = t1107 * t1150;
t1104 = sin(qJ(1));
t1108 = cos(qJ(1));
t1085 = t1104 * g(1) - g(2) * t1108;
t1160 = pkin(8) * t1099;
t1061 = qJDD(1) * pkin(1) + t1109 * t1160 + t1085;
t1086 = -g(1) * t1108 - g(2) * t1104;
t1141 = qJDD(1) * t1099;
t1062 = -pkin(1) * t1109 + pkin(8) * t1141 + t1086;
t1145 = t1100 * t1107;
t1012 = -g(3) * t1147 + t1061 * t1145 - t1103 * t1062;
t1155 = qJD(1) * t1099;
t1063 = (-pkin(2) * t1107 - qJ(3) * t1103) * t1155;
t1090 = t1092 ^ 2;
t1148 = t1099 * t1103;
t1136 = qJD(1) * t1148;
t995 = -t1091 * pkin(2) - t1090 * qJ(3) + t1063 * t1136 + qJDD(3) - t1012;
t982 = t995 - (t1103 * t1133 + t1091) * qJ(4) - (-t1068 + t1126) * pkin(3) + qJD(4) * t1167;
t1146 = t1100 * t1103;
t1144 = t1061 * t1146 + t1107 * t1062;
t1163 = t1090 * pkin(2) - t1091 * qJ(3) + qJD(3) * t1167 - t1063 * t1135 - t1144;
t1162 = mrSges(3,1) - mrSges(4,2);
t1161 = -mrSges(3,3) - mrSges(5,1);
t1159 = t1100 * g(3);
t1158 = t1068 * mrSges(5,1);
t1149 = t1099 * t1061;
t1029 = -t1149 - t1159;
t1054 = mrSges(3,1) * t1092 - mrSges(3,3) * t1136;
t1055 = -mrSges(3,2) * t1092 + mrSges(3,3) * t1135;
t1059 = mrSges(5,1) * t1135 + mrSges(5,2) * t1092;
t1069 = -qJD(2) * t1136 + t1107 * t1141;
t1058 = -mrSges(4,1) * t1135 - mrSges(4,3) * t1092;
t1102 = sin(qJ(5));
t1106 = cos(qJ(5));
t1044 = t1092 * t1106 + t1102 * t1136;
t1010 = -qJD(5) * t1044 + t1068 * t1106 - t1091 * t1102;
t1043 = -t1092 * t1102 + t1106 * t1136;
t1014 = -mrSges(6,1) * t1043 + mrSges(6,2) * t1044;
t1077 = qJD(5) + t1135;
t1019 = mrSges(6,1) * t1077 - mrSges(6,3) * t1044;
t1052 = qJDD(5) + t1069;
t1101 = sin(qJ(6));
t1105 = cos(qJ(6));
t1007 = qJDD(6) - t1010;
t1017 = t1044 * t1105 + t1077 * t1101;
t1042 = qJD(6) - t1043;
t1015 = -pkin(5) * t1043 - pkin(10) * t1044;
t1074 = t1077 ^ 2;
t1067 = pkin(4) * t1135 - pkin(9) * t1092;
t1097 = t1103 ^ 2;
t1098 = t1107 ^ 2;
t1122 = -0.2e1 * qJD(3) * t1136 - t1159 + (t1092 * t1136 - t1069) * pkin(2);
t1118 = -t1069 * qJ(4) - 0.2e1 * qJD(4) * t1135 + t1122;
t1056 = pkin(3) * t1136 - qJ(4) * t1092;
t1152 = t1056 * t1103;
t1156 = qJ(3) * t1092;
t975 = (-pkin(3) * t1098 - pkin(4) * t1097) * t1150 + (pkin(9) - qJ(3)) * t1068 + (-t1061 + (-t1152 + (-t1067 - t1156) * t1107) * qJD(1)) * t1099 + t1118;
t1134 = t1098 * t1150;
t1112 = -qJ(4) * t1134 + t1092 * t1056 + qJDD(4) - t1163;
t979 = -t1091 * pkin(9) + (pkin(3) + pkin(4)) * t1069 + (pkin(9) * t1133 + (pkin(4) * qJD(1) * t1092 - g(3)) * t1099) * t1103 + t1112;
t972 = t1102 * t979 + t1106 * t975;
t969 = -pkin(5) * t1074 + pkin(10) * t1052 + t1015 * t1043 + t972;
t1011 = qJD(5) * t1043 + t1068 * t1102 + t1091 * t1106;
t978 = -pkin(9) * t1097 * t1150 - t1068 * pkin(4) + t1092 * t1067 - t982;
t973 = t978 + (-t1043 * t1077 - t1011) * pkin(10) + (t1044 * t1077 - t1010) * pkin(5);
t966 = -t1101 * t969 + t1105 * t973;
t1016 = -t1044 * t1101 + t1077 * t1105;
t987 = qJD(6) * t1016 + t1011 * t1105 + t1052 * t1101;
t996 = -mrSges(7,1) * t1016 + mrSges(7,2) * t1017;
t998 = -mrSges(7,2) * t1042 + mrSges(7,3) * t1016;
t963 = m(7) * t966 + mrSges(7,1) * t1007 - mrSges(7,3) * t987 - t1017 * t996 + t1042 * t998;
t967 = t1101 * t973 + t1105 * t969;
t986 = -qJD(6) * t1017 - t1011 * t1101 + t1052 * t1105;
t999 = mrSges(7,1) * t1042 - mrSges(7,3) * t1017;
t964 = m(7) * t967 - mrSges(7,2) * t1007 + mrSges(7,3) * t986 + t1016 * t996 - t1042 * t999;
t1129 = -t1101 * t963 + t1105 * t964;
t951 = m(6) * t972 - mrSges(6,2) * t1052 + mrSges(6,3) * t1010 + t1014 * t1043 - t1019 * t1077 + t1129;
t1018 = -mrSges(6,2) * t1077 + mrSges(6,3) * t1043;
t971 = -t1102 * t975 + t1106 * t979;
t968 = -pkin(5) * t1052 - pkin(10) * t1074 + t1015 * t1044 - t971;
t1119 = -m(7) * t968 + t986 * mrSges(7,1) - mrSges(7,2) * t987 + t1016 * t998 - t1017 * t999;
t959 = m(6) * t971 + mrSges(6,1) * t1052 - mrSges(6,3) * t1011 - t1014 * t1044 + t1018 * t1077 + t1119;
t1128 = -t1102 * t959 + t1106 * t951;
t981 = -pkin(3) * t1134 - t1068 * qJ(3) + (-t1061 + (-t1107 * t1156 - t1152) * qJD(1)) * t1099 + t1118;
t1125 = m(5) * t981 - t1069 * mrSges(5,3) + t1128;
t990 = -t1149 + (-t1068 - t1126) * qJ(3) + t1122;
t1120 = m(4) * t990 - t1068 * mrSges(4,3) + t1058 * t1135 + t1125;
t1057 = mrSges(5,1) * t1136 - mrSges(5,3) * t1092;
t1060 = mrSges(4,1) * t1136 + mrSges(4,2) * t1092;
t1143 = -t1057 - t1060;
t936 = m(3) * t1029 - t1162 * t1069 + (mrSges(3,2) - mrSges(5,2)) * t1068 + ((-t1055 - t1059) * t1107 + (t1054 + t1143) * t1103) * t1155 + t1120;
t1137 = g(3) * t1148;
t1013 = -t1137 + t1144;
t1065 = (-mrSges(3,1) * t1107 + mrSges(3,2) * t1103) * t1155;
t1064 = (mrSges(4,2) * t1107 - mrSges(4,3) * t1103) * t1155;
t1066 = (-mrSges(5,2) * t1103 - mrSges(5,3) * t1107) * t1155;
t943 = t1102 * t951 + t1106 * t959;
t984 = t1069 * pkin(3) + t1112 - t1137;
t1124 = m(5) * t984 + t1091 * mrSges(5,2) + t1092 * t1057 + t1066 * t1135 + t943;
t989 = t1137 + t1163;
t1116 = -m(4) * t989 + t1091 * mrSges(4,3) + t1092 * t1060 + t1064 * t1135 + t1124;
t940 = t1065 * t1135 + t1116 - t1091 * mrSges(3,2) - t1092 * t1054 + m(3) * t1013 + (mrSges(4,1) - t1161) * t1069;
t953 = t1101 * t964 + t1105 * t963;
t1121 = -m(6) * t978 + t1010 * mrSges(6,1) - t1011 * mrSges(6,2) + t1043 * t1018 - t1044 * t1019 - t953;
t1115 = -m(5) * t982 + t1091 * mrSges(5,3) + t1092 * t1059 - t1121;
t1111 = m(4) * t995 + t1068 * mrSges(4,1) - t1115;
t1142 = t1064 + t1066;
t946 = -t1111 + (t1055 - t1058) * t1092 + t1162 * t1091 + t1161 * t1068 + m(3) * t1012 + (-t1065 - t1142) * t1136;
t927 = -t1099 * t936 + t946 * t1145 + t940 * t1146;
t924 = m(2) * t1085 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1109 + t927;
t932 = -t1103 * t946 + t1107 * t940;
t930 = m(2) * t1086 - mrSges(2,1) * t1109 - qJDD(1) * mrSges(2,2) + t932;
t1157 = t1104 * t930 + t1108 * t924;
t1151 = t1059 * t1107;
t926 = t1100 * t936 + t946 * t1147 + t940 * t1148;
t1132 = (t1103 * t1139 + t1107 * t1138) * t1155 + t1164 * t1092;
t1131 = (t1103 * t1140 + t1107 * t1165) * t1155 + t1138 * t1092;
t1130 = (t1103 * t1166 + t1107 * t1140) * t1155 + t1139 * t1092;
t1127 = -t1104 * t924 + t1108 * t930;
t1000 = Ifges(6,5) * t1044 + Ifges(6,6) * t1043 + Ifges(6,3) * t1077;
t1001 = Ifges(6,4) * t1044 + Ifges(6,2) * t1043 + Ifges(6,6) * t1077;
t991 = Ifges(7,5) * t1017 + Ifges(7,6) * t1016 + Ifges(7,3) * t1042;
t993 = Ifges(7,1) * t1017 + Ifges(7,4) * t1016 + Ifges(7,5) * t1042;
t956 = -mrSges(7,1) * t968 + mrSges(7,3) * t967 + Ifges(7,4) * t987 + Ifges(7,2) * t986 + Ifges(7,6) * t1007 - t1017 * t991 + t1042 * t993;
t992 = Ifges(7,4) * t1017 + Ifges(7,2) * t1016 + Ifges(7,6) * t1042;
t957 = mrSges(7,2) * t968 - mrSges(7,3) * t966 + Ifges(7,1) * t987 + Ifges(7,4) * t986 + Ifges(7,5) * t1007 + t1016 * t991 - t1042 * t992;
t933 = mrSges(6,2) * t978 - mrSges(6,3) * t971 + Ifges(6,1) * t1011 + Ifges(6,4) * t1010 + Ifges(6,5) * t1052 - pkin(10) * t953 + t1000 * t1043 - t1001 * t1077 - t1101 * t956 + t1105 * t957;
t1002 = Ifges(6,1) * t1044 + Ifges(6,4) * t1043 + Ifges(6,5) * t1077;
t1110 = mrSges(7,1) * t966 - mrSges(7,2) * t967 + Ifges(7,5) * t987 + Ifges(7,6) * t986 + Ifges(7,3) * t1007 - t1016 * t993 + t1017 * t992;
t934 = -mrSges(6,1) * t978 + mrSges(6,3) * t972 + Ifges(6,4) * t1011 + Ifges(6,2) * t1010 + Ifges(6,6) * t1052 - pkin(5) * t953 - t1000 * t1044 + t1002 * t1077 - t1110;
t947 = t1091 * mrSges(4,2) + t1092 * t1058 + t1136 * t1142 + t1111 + t1158;
t948 = t1066 * t1136 - t1115 + t1158;
t918 = t1106 * t933 - t1102 * t934 + qJ(3) * t1116 + mrSges(3,1) * t1012 - mrSges(3,2) * t1013 + mrSges(4,2) * t995 - mrSges(5,3) * t982 + mrSges(5,2) * t984 - mrSges(4,3) * t989 - qJ(4) * t948 - pkin(2) * t947 - pkin(9) * t943 + t1164 * t1091 + t1139 * t1068 + (qJ(3) * (mrSges(4,1) + mrSges(5,1)) + t1138) * t1069 + (t1103 * t1131 - t1107 * t1130) * t1155;
t941 = t1069 * mrSges(4,2) - t1068 * mrSges(5,2) + (t1103 * t1143 - t1151) * t1155 + t1120;
t920 = t1106 * t934 + pkin(9) * t1128 + t1102 * t933 + mrSges(3,2) * t1029 + pkin(4) * t1121 - mrSges(3,3) * t1012 + mrSges(4,1) * t995 + mrSges(5,1) * t982 - mrSges(4,3) * t990 - mrSges(5,2) * t981 + pkin(3) * t948 - qJ(3) * t941 - t1131 * t1092 + t1139 * t1091 + t1140 * t1069 + t1166 * t1068 + t1132 * t1135;
t1113 = mrSges(6,1) * t971 - mrSges(6,2) * t972 + Ifges(6,5) * t1011 + Ifges(6,6) * t1010 + Ifges(6,3) * t1052 + pkin(5) * t1119 + pkin(10) * t1129 + t1044 * t1001 - t1043 * t1002 + t1101 * t957 + t1105 * t956;
t942 = mrSges(5,1) * t1069 + t1124;
t922 = t1113 + t1130 * t1092 + t1138 * t1091 + t1165 * t1069 + (mrSges(5,2) * qJ(4) + t1140) * t1068 + (qJ(4) * t1151 + (qJ(4) * t1057 - t1132) * t1103) * t1155 - mrSges(3,1) * t1029 + mrSges(3,3) * t1013 + mrSges(5,1) * t984 - mrSges(4,1) * t989 + mrSges(4,2) * t990 - mrSges(5,3) * t981 + pkin(4) * t943 + pkin(3) * t942 - pkin(2) * t941 - qJ(4) * t1125;
t1117 = mrSges(2,1) * t1085 - mrSges(2,2) * t1086 + Ifges(2,3) * qJDD(1) + pkin(1) * t927 + t1100 * t918 + t922 * t1147 + t920 * t1148 + t932 * t1160;
t916 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1085 + Ifges(2,5) * qJDD(1) - t1109 * Ifges(2,6) - t1103 * t922 + t1107 * t920 + (-t1099 * t926 - t1100 * t927) * pkin(8);
t915 = mrSges(2,1) * g(3) + mrSges(2,3) * t1086 + t1109 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t926 - t1099 * t918 + (pkin(8) * t932 + t1103 * t920 + t1107 * t922) * t1100;
t1 = [-m(1) * g(1) + t1127; -m(1) * g(2) + t1157; (-m(1) - m(2)) * g(3) + t926; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1157 - t1104 * t915 + t1108 * t916; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1127 + t1104 * t916 + t1108 * t915; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1117; t1117; t918; t947; t942; t1113; t1110;];
tauJB  = t1;
