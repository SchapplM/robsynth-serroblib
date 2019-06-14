% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR10
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 04:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:47:09
% EndTime: 2019-05-06 04:49:22
% DurationCPUTime: 137.98s
% Computational Cost: add. (2132773->405), mult. (6613941->548), div. (0->0), fcn. (5718045->16), ass. (0->187)
t1090 = sin(qJ(3));
t1095 = cos(qJ(3));
t1082 = sin(pkin(7));
t1084 = cos(pkin(13));
t1086 = cos(pkin(6));
t1083 = sin(pkin(6));
t1085 = cos(pkin(7));
t1132 = t1083 * t1085;
t1112 = t1082 * t1086 + t1084 * t1132;
t1081 = sin(pkin(13));
t1136 = t1081 * t1083;
t1145 = -t1090 * t1136 + t1112 * t1095;
t1052 = t1145 * qJD(1);
t1130 = t1085 * t1090;
t1134 = t1082 * t1090;
t1102 = t1086 * t1134 + (t1081 * t1095 + t1084 * t1130) * t1083;
t1053 = t1102 * qJD(1);
t1039 = -t1053 * qJD(3) + t1145 * qJDD(1);
t1064 = t1112 * qJD(1) * pkin(9);
t1091 = sin(qJ(1));
t1096 = cos(qJ(1));
t1078 = -t1096 * g(1) - t1091 * g(2);
t1097 = qJD(1) ^ 2;
t1139 = qJ(2) * t1083;
t1068 = -t1097 * pkin(1) + qJDD(1) * t1139 + t1078;
t1143 = pkin(9) * t1081;
t1116 = -pkin(2) * t1084 - t1082 * t1143;
t1137 = qJD(1) * t1083;
t1141 = pkin(9) * qJDD(1);
t1113 = qJD(1) * t1116 * t1137 + t1085 * t1141;
t1077 = t1091 * g(1) - t1096 * g(2);
t1067 = qJDD(1) * pkin(1) + t1097 * t1139 + t1077;
t1128 = qJD(2) * t1137;
t1131 = t1084 * t1086;
t1133 = t1083 * t1084;
t1117 = -g(3) * t1133 + t1067 * t1131 - 0.2e1 * t1081 * t1128;
t1018 = (pkin(2) * qJDD(1) + qJD(1) * t1064) * t1086 + (-t1083 * t1113 - t1068) * t1081 + t1117;
t1069 = (pkin(2) * t1086 - t1132 * t1143) * qJD(1);
t1135 = t1081 * t1086;
t1126 = t1067 * t1135 + (t1068 + 0.2e1 * t1128) * t1084;
t1019 = (-qJD(1) * t1069 + t1082 * t1141) * t1086 + (-g(3) * t1081 + t1084 * t1113) * t1083 + t1126;
t1125 = -t1086 * g(3) + qJDD(2);
t1027 = (-t1067 + t1116 * qJDD(1) + (-t1064 * t1084 + t1069 * t1081) * qJD(1)) * t1083 + t1125;
t991 = -t1090 * t1019 + (t1018 * t1085 + t1027 * t1082) * t1095;
t1142 = Ifges(3,3) * t1086;
t1041 = -t1081 * t1068 + t1117;
t1120 = -mrSges(3,1) * t1084 + mrSges(3,2) * t1081;
t1066 = t1120 * t1137;
t1114 = -mrSges(3,2) * t1086 + mrSges(3,3) * t1133;
t1071 = t1114 * qJD(1);
t1115 = mrSges(3,1) * t1086 - mrSges(3,3) * t1136;
t1037 = -t1052 * mrSges(4,1) + t1053 * mrSges(4,2);
t1040 = t1052 * qJD(3) + qJDD(1) * t1102;
t1111 = -t1082 * t1133 + t1085 * t1086;
t1065 = qJD(1) * t1111 + qJD(3);
t1046 = -t1065 * mrSges(4,2) + t1052 * mrSges(4,3);
t1062 = qJDD(1) * t1111 + qJDD(3);
t1089 = sin(qJ(4));
t1094 = cos(qJ(4));
t1045 = t1094 * t1053 + t1089 * t1065;
t1013 = -t1045 * qJD(4) - t1089 * t1040 + t1094 * t1062;
t1044 = -t1089 * t1053 + t1094 * t1065;
t1014 = t1044 * qJD(4) + t1094 * t1040 + t1089 * t1062;
t1051 = qJD(4) - t1052;
t1028 = -t1051 * mrSges(5,2) + t1044 * mrSges(5,3);
t1029 = t1051 * mrSges(5,1) - t1045 * mrSges(5,3);
t1088 = sin(qJ(5));
t1093 = cos(qJ(5));
t1021 = t1093 * t1044 - t1088 * t1045;
t1049 = qJD(5) + t1051;
t1005 = -t1049 * mrSges(6,2) + t1021 * mrSges(6,3);
t1022 = t1088 * t1044 + t1093 * t1045;
t1006 = t1049 * mrSges(6,1) - t1022 * mrSges(6,3);
t1087 = sin(qJ(6));
t1092 = cos(qJ(6));
t1004 = t1092 * t1022 + t1087 * t1049;
t1020 = qJD(6) - t1021;
t1001 = -t1021 * pkin(5) - t1022 * pkin(12);
t1036 = qJDD(4) - t1039;
t1035 = qJDD(5) + t1036;
t1048 = t1049 ^ 2;
t1038 = -t1052 * pkin(3) - t1053 * pkin(10);
t1061 = t1065 ^ 2;
t992 = t1018 * t1130 + t1095 * t1019 + t1027 * t1134;
t976 = -t1061 * pkin(3) + t1062 * pkin(10) + t1052 * t1038 + t992;
t1002 = -t1082 * t1018 + t1085 * t1027;
t979 = (-t1052 * t1065 - t1040) * pkin(10) + (t1053 * t1065 - t1039) * pkin(3) + t1002;
t968 = -t1089 * t976 + t1094 * t979;
t965 = (t1044 * t1051 - t1014) * pkin(11) + (t1044 * t1045 + t1036) * pkin(4) + t968;
t1030 = t1051 * pkin(4) - t1045 * pkin(11);
t1043 = t1044 ^ 2;
t969 = t1089 * t979 + t1094 * t976;
t967 = -t1043 * pkin(4) + t1013 * pkin(11) - t1051 * t1030 + t969;
t962 = t1088 * t965 + t1093 * t967;
t959 = -t1048 * pkin(5) + t1035 * pkin(12) + t1021 * t1001 + t962;
t975 = -t1062 * pkin(3) - t1061 * pkin(10) + t1053 * t1038 - t991;
t970 = -t1013 * pkin(4) - t1043 * pkin(11) + t1045 * t1030 + t975;
t988 = -t1022 * qJD(5) + t1093 * t1013 - t1088 * t1014;
t989 = t1021 * qJD(5) + t1088 * t1013 + t1093 * t1014;
t963 = (-t1021 * t1049 - t989) * pkin(12) + (t1022 * t1049 - t988) * pkin(5) + t970;
t956 = -t1087 * t959 + t1092 * t963;
t1003 = -t1087 * t1022 + t1092 * t1049;
t973 = t1003 * qJD(6) + t1087 * t1035 + t1092 * t989;
t987 = qJDD(6) - t988;
t993 = -t1003 * mrSges(7,1) + t1004 * mrSges(7,2);
t994 = -t1020 * mrSges(7,2) + t1003 * mrSges(7,3);
t952 = m(7) * t956 + t987 * mrSges(7,1) - t973 * mrSges(7,3) - t1004 * t993 + t1020 * t994;
t957 = t1087 * t963 + t1092 * t959;
t972 = -t1004 * qJD(6) + t1092 * t1035 - t1087 * t989;
t995 = t1020 * mrSges(7,1) - t1004 * mrSges(7,3);
t953 = m(7) * t957 - t987 * mrSges(7,2) + t972 * mrSges(7,3) + t1003 * t993 - t1020 * t995;
t941 = t1087 * t953 + t1092 * t952;
t1103 = m(6) * t970 - t988 * mrSges(6,1) + t989 * mrSges(6,2) - t1021 * t1005 + t1022 * t1006 + t941;
t1099 = -m(5) * t975 + t1013 * mrSges(5,1) - t1014 * mrSges(5,2) + t1044 * t1028 - t1045 * t1029 - t1103;
t936 = m(4) * t991 + t1062 * mrSges(4,1) - t1040 * mrSges(4,3) - t1053 * t1037 + t1065 * t1046 + t1099;
t1138 = t1095 * t936;
t1047 = t1065 * mrSges(4,1) - t1053 * mrSges(4,3);
t1023 = -t1044 * mrSges(5,1) + t1045 * mrSges(5,2);
t1000 = -t1021 * mrSges(6,1) + t1022 * mrSges(6,2);
t1124 = -t1087 * t952 + t1092 * t953;
t939 = m(6) * t962 - t1035 * mrSges(6,2) + t988 * mrSges(6,3) + t1021 * t1000 - t1049 * t1006 + t1124;
t961 = -t1088 * t967 + t1093 * t965;
t958 = -t1035 * pkin(5) - t1048 * pkin(12) + t1022 * t1001 - t961;
t1107 = -m(7) * t958 + t972 * mrSges(7,1) - t973 * mrSges(7,2) + t1003 * t994 - t1004 * t995;
t948 = m(6) * t961 + t1035 * mrSges(6,1) - t989 * mrSges(6,3) - t1022 * t1000 + t1049 * t1005 + t1107;
t933 = t1088 * t939 + t1093 * t948;
t931 = m(5) * t968 + t1036 * mrSges(5,1) - t1014 * mrSges(5,3) - t1045 * t1023 + t1051 * t1028 + t933;
t1123 = -t1088 * t948 + t1093 * t939;
t932 = m(5) * t969 - t1036 * mrSges(5,2) + t1013 * mrSges(5,3) + t1044 * t1023 - t1051 * t1029 + t1123;
t1122 = -t1089 * t931 + t1094 * t932;
t922 = m(4) * t992 - t1062 * mrSges(4,2) + t1039 * mrSges(4,3) + t1052 * t1037 - t1065 * t1047 + t1122;
t925 = t1089 * t932 + t1094 * t931;
t924 = m(4) * t1002 - t1039 * mrSges(4,1) + t1040 * mrSges(4,2) - t1052 * t1046 + t1053 * t1047 + t925;
t911 = -t1082 * t924 + t1085 * t1138 + t922 * t1130;
t907 = m(3) * t1041 + t1115 * qJDD(1) + (-t1066 * t1136 + t1071 * t1086) * qJD(1) + t911;
t1054 = -t1083 * t1067 + t1125;
t1070 = t1115 * qJD(1);
t910 = t1082 * t1138 + t1085 * t924 + t922 * t1134;
t909 = m(3) * t1054 + (t1120 * qJDD(1) + (t1070 * t1081 - t1071 * t1084) * qJD(1)) * t1083 + t910;
t1042 = -g(3) * t1136 + t1126;
t918 = -t1090 * t936 + t1095 * t922;
t917 = m(3) * t1042 + t1114 * qJDD(1) + (t1066 * t1133 - t1070 * t1086) * qJD(1) + t918;
t896 = -t1083 * t909 + t907 * t1131 + t917 * t1135;
t893 = m(2) * t1077 + qJDD(1) * mrSges(2,1) - t1097 * mrSges(2,2) + t896;
t903 = -t1081 * t907 + t1084 * t917;
t901 = m(2) * t1078 - t1097 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t903;
t1140 = t1091 * t901 + t1096 * t893;
t895 = t1086 * t909 + t907 * t1133 + t917 * t1136;
t1121 = -t1091 * t893 + t1096 * t901;
t1119 = Ifges(3,5) * t1081 + Ifges(3,6) * t1084;
t1031 = Ifges(4,5) * t1053 + Ifges(4,6) * t1052 + Ifges(4,3) * t1065;
t1032 = Ifges(4,4) * t1053 + Ifges(4,2) * t1052 + Ifges(4,6) * t1065;
t1007 = Ifges(5,5) * t1045 + Ifges(5,6) * t1044 + Ifges(5,3) * t1051;
t1009 = Ifges(5,1) * t1045 + Ifges(5,4) * t1044 + Ifges(5,5) * t1051;
t980 = Ifges(7,5) * t1004 + Ifges(7,6) * t1003 + Ifges(7,3) * t1020;
t982 = Ifges(7,1) * t1004 + Ifges(7,4) * t1003 + Ifges(7,5) * t1020;
t945 = -mrSges(7,1) * t958 + mrSges(7,3) * t957 + Ifges(7,4) * t973 + Ifges(7,2) * t972 + Ifges(7,6) * t987 - t1004 * t980 + t1020 * t982;
t981 = Ifges(7,4) * t1004 + Ifges(7,2) * t1003 + Ifges(7,6) * t1020;
t946 = mrSges(7,2) * t958 - mrSges(7,3) * t956 + Ifges(7,1) * t973 + Ifges(7,4) * t972 + Ifges(7,5) * t987 + t1003 * t980 - t1020 * t981;
t996 = Ifges(6,5) * t1022 + Ifges(6,6) * t1021 + Ifges(6,3) * t1049;
t997 = Ifges(6,4) * t1022 + Ifges(6,2) * t1021 + Ifges(6,6) * t1049;
t926 = mrSges(6,2) * t970 - mrSges(6,3) * t961 + Ifges(6,1) * t989 + Ifges(6,4) * t988 + Ifges(6,5) * t1035 - pkin(12) * t941 + t1021 * t996 - t1049 * t997 - t1087 * t945 + t1092 * t946;
t1100 = mrSges(7,1) * t956 - mrSges(7,2) * t957 + Ifges(7,5) * t973 + Ifges(7,6) * t972 + Ifges(7,3) * t987 - t1003 * t982 + t1004 * t981;
t998 = Ifges(6,1) * t1022 + Ifges(6,4) * t1021 + Ifges(6,5) * t1049;
t927 = -mrSges(6,1) * t970 + mrSges(6,3) * t962 + Ifges(6,4) * t989 + Ifges(6,2) * t988 + Ifges(6,6) * t1035 - pkin(5) * t941 - t1022 * t996 + t1049 * t998 - t1100;
t912 = -mrSges(5,1) * t975 + mrSges(5,3) * t969 + Ifges(5,4) * t1014 + Ifges(5,2) * t1013 + Ifges(5,6) * t1036 - pkin(4) * t1103 + pkin(11) * t1123 - t1045 * t1007 + t1051 * t1009 + t1088 * t926 + t1093 * t927;
t1008 = Ifges(5,4) * t1045 + Ifges(5,2) * t1044 + Ifges(5,6) * t1051;
t913 = mrSges(5,2) * t975 - mrSges(5,3) * t968 + Ifges(5,1) * t1014 + Ifges(5,4) * t1013 + Ifges(5,5) * t1036 - pkin(11) * t933 + t1044 * t1007 - t1051 * t1008 - t1088 * t927 + t1093 * t926;
t898 = mrSges(4,2) * t1002 - mrSges(4,3) * t991 + Ifges(4,1) * t1040 + Ifges(4,4) * t1039 + Ifges(4,5) * t1062 - pkin(10) * t925 + t1052 * t1031 - t1065 * t1032 - t1089 * t912 + t1094 * t913;
t1033 = Ifges(4,1) * t1053 + Ifges(4,4) * t1052 + Ifges(4,5) * t1065;
t1101 = -mrSges(6,1) * t961 + mrSges(6,2) * t962 - Ifges(6,5) * t989 - Ifges(6,6) * t988 - Ifges(6,3) * t1035 - pkin(5) * t1107 - pkin(12) * t1124 + t1021 * t998 - t1022 * t997 - t1087 * t946 - t1092 * t945;
t1098 = mrSges(5,1) * t968 - mrSges(5,2) * t969 + Ifges(5,5) * t1014 + Ifges(5,6) * t1013 + Ifges(5,3) * t1036 + pkin(4) * t933 + t1045 * t1008 - t1044 * t1009 - t1101;
t904 = -mrSges(4,1) * t1002 + mrSges(4,3) * t992 + Ifges(4,4) * t1040 + Ifges(4,2) * t1039 + Ifges(4,6) * t1062 - pkin(3) * t925 - t1053 * t1031 + t1065 * t1033 - t1098;
t1110 = pkin(9) * t918 + t1090 * t898 + t1095 * t904;
t1104 = Ifges(3,6) * t1086 + (Ifges(3,4) * t1081 + Ifges(3,2) * t1084) * t1083;
t1058 = t1104 * qJD(1);
t1105 = Ifges(3,5) * t1086 + (Ifges(3,1) * t1081 + Ifges(3,4) * t1084) * t1083;
t1059 = t1105 * qJD(1);
t897 = mrSges(4,1) * t991 - mrSges(4,2) * t992 + Ifges(4,5) * t1040 + Ifges(4,6) * t1039 + Ifges(4,3) * t1062 + pkin(3) * t1099 + pkin(10) * t1122 + t1053 * t1032 - t1052 * t1033 + t1089 * t913 + t1094 * t912;
t887 = qJDD(1) * t1142 + mrSges(3,1) * t1041 - mrSges(3,2) * t1042 + pkin(2) * t911 + t1085 * t897 + t1110 * t1082 + (t1119 * qJDD(1) + (t1058 * t1081 - t1059 * t1084) * qJD(1)) * t1083;
t1057 = (t1083 * t1119 + t1142) * qJD(1);
t889 = -mrSges(3,1) * t1054 + mrSges(3,3) * t1042 - pkin(2) * t910 - t1082 * t897 + (-t1057 * t1136 + t1059 * t1086) * qJD(1) + t1110 * t1085 + t1104 * qJDD(1);
t891 = mrSges(3,2) * t1054 - mrSges(3,3) * t1041 - t1090 * t904 + t1095 * t898 + (t1057 * t1133 - t1058 * t1086) * qJD(1) + (-t1082 * t910 - t1085 * t911) * pkin(9) + t1105 * qJDD(1);
t1106 = mrSges(2,1) * t1077 - mrSges(2,2) * t1078 + Ifges(2,3) * qJDD(1) + pkin(1) * t896 + t1086 * t887 + t889 * t1133 + t891 * t1136 + t903 * t1139;
t885 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1077 + Ifges(2,5) * qJDD(1) - t1097 * Ifges(2,6) - t1081 * t889 + t1084 * t891 + (-t1083 * t895 - t1086 * t896) * qJ(2);
t884 = mrSges(2,1) * g(3) + mrSges(2,3) * t1078 + t1097 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t895 - t1083 * t887 + (qJ(2) * t903 + t1081 * t891 + t1084 * t889) * t1086;
t1 = [-m(1) * g(1) + t1121; -m(1) * g(2) + t1140; (-m(1) - m(2)) * g(3) + t895; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1140 - t1091 * t884 + t1096 * t885; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1121 + t1091 * t885 + t1096 * t884; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1106; t1106; t909; t897; t1098; -t1101; t1100;];
tauJB  = t1;
