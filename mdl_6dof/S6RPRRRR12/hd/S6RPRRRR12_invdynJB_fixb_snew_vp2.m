% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 07:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 06:37:41
% EndTime: 2019-05-06 06:43:23
% DurationCPUTime: 352.04s
% Computational Cost: add. (5357154->420), mult. (17472500->582), div. (0->0), fcn. (15362513->18), ass. (0->200)
t1092 = sin(pkin(7));
t1094 = cos(pkin(14));
t1097 = cos(pkin(6));
t1093 = sin(pkin(6));
t1096 = cos(pkin(7));
t1141 = t1093 * t1096;
t1075 = (t1092 * t1097 + t1094 * t1141) * qJD(1) * pkin(10);
t1102 = sin(qJ(1));
t1107 = cos(qJ(1));
t1087 = -g(1) * t1107 - g(2) * t1102;
t1108 = qJD(1) ^ 2;
t1150 = qJ(2) * t1093;
t1079 = -pkin(1) * t1108 + qJDD(1) * t1150 + t1087;
t1090 = sin(pkin(14));
t1156 = pkin(10) * t1090;
t1123 = -pkin(2) * t1094 - t1092 * t1156;
t1148 = qJD(1) * t1093;
t1152 = pkin(10) * qJDD(1);
t1120 = qJD(1) * t1123 * t1148 + t1096 * t1152;
t1086 = t1102 * g(1) - g(2) * t1107;
t1078 = qJDD(1) * pkin(1) + t1108 * t1150 + t1086;
t1135 = qJD(2) * t1148;
t1140 = t1094 * t1097;
t1142 = t1093 * t1094;
t1124 = -g(3) * t1142 + t1078 * t1140 - 0.2e1 * t1090 * t1135;
t1038 = (pkin(2) * qJDD(1) + qJD(1) * t1075) * t1097 + (-t1093 * t1120 - t1079) * t1090 + t1124;
t1080 = (pkin(2) * t1097 - t1141 * t1156) * qJD(1);
t1133 = -t1097 * g(3) + qJDD(2);
t1047 = (-t1078 + t1123 * qJDD(1) + (-t1075 * t1094 + t1080 * t1090) * qJD(1)) * t1093 + t1133;
t1027 = -t1038 * t1092 + t1096 * t1047;
t1106 = cos(qJ(3));
t1101 = sin(qJ(3));
t1138 = t1096 * t1101;
t1144 = t1092 * t1101;
t1113 = t1097 * t1144 + (t1090 * t1106 + t1094 * t1138) * t1093;
t1065 = t1113 * qJD(1);
t1137 = t1096 * t1106;
t1143 = t1092 * t1106;
t1112 = t1097 * t1143 + (-t1090 * t1101 + t1094 * t1137) * t1093;
t1054 = -t1065 * qJD(3) + qJDD(1) * t1112;
t1064 = t1112 * qJD(1);
t1055 = t1064 * qJD(3) + qJDD(1) * t1113;
t1119 = -t1092 * t1142 + t1096 * t1097;
t1076 = qJD(1) * t1119 + qJD(3);
t1091 = sin(pkin(8));
t1095 = cos(pkin(8));
t1125 = t1064 * t1095 + t1076 * t1091;
t1057 = t1125 * pkin(11);
t1154 = pkin(11) * t1095;
t1061 = pkin(3) * t1076 - t1065 * t1154;
t1155 = pkin(11) * t1091;
t1000 = -pkin(3) * t1054 - t1055 * t1155 - t1057 * t1064 + t1061 * t1065 + t1027;
t1100 = sin(qJ(4));
t1105 = cos(qJ(4));
t1146 = t1090 * t1097;
t1134 = t1078 * t1146 + (t1079 + 0.2e1 * t1135) * t1094;
t1039 = (-qJD(1) * t1080 + t1092 * t1152) * t1097 + (-g(3) * t1090 + t1094 * t1120) * t1093 + t1134;
t1012 = t1038 * t1137 - t1039 * t1101 + t1047 * t1143;
t1052 = -pkin(3) * t1064 - t1065 * t1155;
t1073 = qJDD(1) * t1119 + qJDD(3);
t996 = pkin(3) * t1073 - t1052 * t1065 - t1055 * t1154 + t1057 * t1076 + t1012;
t1013 = t1038 * t1138 + t1106 * t1039 + t1047 * t1144;
t1126 = t1054 * t1095 + t1073 * t1091;
t997 = pkin(11) * t1126 + t1064 * t1052 - t1076 * t1061 + t1013;
t983 = -t1100 * t997 + (t1000 * t1091 + t1095 * t996) * t1105;
t1045 = -t1100 * t1065 + t1105 * t1125;
t1046 = t1105 * t1065 + t1100 * t1125;
t1020 = -t1046 * qJD(4) - t1100 * t1055 + t1105 * t1126;
t1153 = Ifges(3,3) * t1097;
t1059 = -t1090 * t1079 + t1124;
t1130 = -mrSges(3,1) * t1094 + mrSges(3,2) * t1090;
t1077 = t1130 * t1148;
t1121 = -mrSges(3,2) * t1097 + mrSges(3,3) * t1142;
t1082 = t1121 * qJD(1);
t1147 = t1090 * t1093;
t1122 = mrSges(3,1) * t1097 - mrSges(3,3) * t1147;
t1053 = -mrSges(4,1) * t1064 + mrSges(4,2) * t1065;
t1062 = -mrSges(4,2) * t1076 + mrSges(4,3) * t1064;
t1139 = t1095 * t1100;
t1021 = t1045 * qJD(4) + t1105 * t1055 + t1100 * t1126;
t1028 = -mrSges(5,1) * t1045 + mrSges(5,2) * t1046;
t1058 = -t1064 * t1091 + t1076 * t1095 + qJD(4);
t1033 = -mrSges(5,2) * t1058 + mrSges(5,3) * t1045;
t1048 = -t1054 * t1091 + t1073 * t1095 + qJDD(4);
t1099 = sin(qJ(5));
t1104 = cos(qJ(5));
t1032 = t1046 * t1104 + t1058 * t1099;
t1004 = -qJD(5) * t1032 - t1021 * t1099 + t1048 * t1104;
t1031 = -t1046 * t1099 + t1058 * t1104;
t1005 = qJD(5) * t1031 + t1021 * t1104 + t1048 * t1099;
t1043 = qJD(5) - t1045;
t1022 = -mrSges(6,2) * t1043 + mrSges(6,3) * t1031;
t1023 = mrSges(6,1) * t1043 - mrSges(6,3) * t1032;
t1098 = sin(qJ(6));
t1103 = cos(qJ(6));
t1017 = -t1032 * t1098 + t1043 * t1103;
t1018 = t1032 * t1103 + t1043 * t1098;
t1001 = -mrSges(7,1) * t1017 + mrSges(7,2) * t1018;
t1003 = qJDD(6) - t1004;
t1030 = qJD(6) - t1031;
t1006 = -mrSges(7,2) * t1030 + mrSges(7,3) * t1017;
t1015 = -pkin(5) * t1031 - pkin(13) * t1032;
t1019 = qJDD(5) - t1020;
t1042 = t1043 ^ 2;
t1029 = -pkin(4) * t1045 - pkin(12) * t1046;
t1056 = t1058 ^ 2;
t1145 = t1091 * t1100;
t984 = t1000 * t1145 + t1105 * t997 + t996 * t1139;
t980 = -pkin(4) * t1056 + pkin(12) * t1048 + t1029 * t1045 + t984;
t985 = t1095 * t1000 - t1091 * t996;
t982 = (-t1045 * t1058 - t1021) * pkin(12) + (t1046 * t1058 - t1020) * pkin(4) + t985;
t976 = t1099 * t982 + t1104 * t980;
t974 = -pkin(5) * t1042 + pkin(13) * t1019 + t1015 * t1031 + t976;
t979 = -t1048 * pkin(4) - t1056 * pkin(12) + t1046 * t1029 - t983;
t977 = (-t1031 * t1043 - t1005) * pkin(13) + (t1032 * t1043 - t1004) * pkin(5) + t979;
t970 = -t1098 * t974 + t1103 * t977;
t988 = qJD(6) * t1017 + t1005 * t1103 + t1019 * t1098;
t968 = m(7) * t970 + mrSges(7,1) * t1003 - mrSges(7,3) * t988 - t1001 * t1018 + t1006 * t1030;
t1007 = mrSges(7,1) * t1030 - mrSges(7,3) * t1018;
t971 = t1098 * t977 + t1103 * t974;
t987 = -qJD(6) * t1018 - t1005 * t1098 + t1019 * t1103;
t969 = m(7) * t971 - mrSges(7,2) * t1003 + mrSges(7,3) * t987 + t1001 * t1017 - t1007 * t1030;
t961 = t1098 * t969 + t1103 * t968;
t1111 = -m(6) * t979 + t1004 * mrSges(6,1) - mrSges(6,2) * t1005 + t1031 * t1022 - t1023 * t1032 - t961;
t957 = m(5) * t983 + mrSges(5,1) * t1048 - mrSges(5,3) * t1021 - t1028 * t1046 + t1033 * t1058 + t1111;
t1149 = t1105 * t957;
t1034 = mrSges(5,1) * t1058 - mrSges(5,3) * t1046;
t1014 = -mrSges(6,1) * t1031 + mrSges(6,2) * t1032;
t962 = -t1098 * t968 + t1103 * t969;
t960 = m(6) * t976 - mrSges(6,2) * t1019 + mrSges(6,3) * t1004 + t1014 * t1031 - t1023 * t1043 + t962;
t975 = -t1099 * t980 + t1104 * t982;
t973 = -pkin(5) * t1019 - pkin(13) * t1042 + t1015 * t1032 - t975;
t972 = -m(7) * t973 + t987 * mrSges(7,1) - mrSges(7,2) * t988 + t1017 * t1006 - t1007 * t1018;
t966 = m(6) * t975 + mrSges(6,1) * t1019 - mrSges(6,3) * t1005 - t1014 * t1032 + t1022 * t1043 + t972;
t1132 = -t1099 * t966 + t1104 * t960;
t951 = m(5) * t984 - mrSges(5,2) * t1048 + mrSges(5,3) * t1020 + t1028 * t1045 - t1034 * t1058 + t1132;
t954 = t1099 * t960 + t1104 * t966;
t953 = m(5) * t985 - mrSges(5,1) * t1020 + mrSges(5,2) * t1021 - t1033 * t1045 + t1034 * t1046 + t954;
t940 = -t1091 * t953 + t1095 * t1149 + t951 * t1139;
t936 = m(4) * t1012 + mrSges(4,1) * t1073 - mrSges(4,3) * t1055 - t1053 * t1065 + t1062 * t1076 + t940;
t1063 = mrSges(4,1) * t1076 - mrSges(4,3) * t1065;
t939 = t1091 * t1149 + t1095 * t953 + t951 * t1145;
t938 = m(4) * t1027 - mrSges(4,1) * t1054 + mrSges(4,2) * t1055 - t1062 * t1064 + t1063 * t1065 + t939;
t945 = -t1100 * t957 + t1105 * t951;
t944 = m(4) * t1013 - mrSges(4,2) * t1073 + mrSges(4,3) * t1054 + t1053 * t1064 - t1063 * t1076 + t945;
t925 = -t1092 * t938 + t936 * t1137 + t944 * t1138;
t921 = m(3) * t1059 + t1122 * qJDD(1) + (-t1077 * t1147 + t1082 * t1097) * qJD(1) + t925;
t1066 = -t1093 * t1078 + t1133;
t1081 = t1122 * qJD(1);
t924 = t1096 * t938 + t936 * t1143 + t944 * t1144;
t923 = m(3) * t1066 + (t1130 * qJDD(1) + (t1081 * t1090 - t1082 * t1094) * qJD(1)) * t1093 + t924;
t1060 = -g(3) * t1147 + t1134;
t930 = -t1101 * t936 + t1106 * t944;
t929 = m(3) * t1060 + t1121 * qJDD(1) + (t1077 * t1142 - t1081 * t1097) * qJD(1) + t930;
t910 = -t1093 * t923 + t921 * t1140 + t929 * t1146;
t907 = m(2) * t1086 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1108 + t910;
t918 = -t1090 * t921 + t1094 * t929;
t916 = m(2) * t1087 - mrSges(2,1) * t1108 - qJDD(1) * mrSges(2,2) + t918;
t1151 = t1102 * t916 + t1107 * t907;
t909 = t1097 * t923 + t921 * t1142 + t929 * t1147;
t1131 = -t1102 * t907 + t1107 * t916;
t1129 = Ifges(3,5) * t1090 + Ifges(3,6) * t1094;
t1049 = Ifges(4,5) * t1065 + Ifges(4,6) * t1064 + Ifges(4,3) * t1076;
t1051 = Ifges(4,1) * t1065 + Ifges(4,4) * t1064 + Ifges(4,5) * t1076;
t1024 = Ifges(5,5) * t1046 + Ifges(5,6) * t1045 + Ifges(5,3) * t1058;
t1025 = Ifges(5,4) * t1046 + Ifges(5,2) * t1045 + Ifges(5,6) * t1058;
t1008 = Ifges(6,5) * t1032 + Ifges(6,6) * t1031 + Ifges(6,3) * t1043;
t1009 = Ifges(6,4) * t1032 + Ifges(6,2) * t1031 + Ifges(6,6) * t1043;
t989 = Ifges(7,5) * t1018 + Ifges(7,6) * t1017 + Ifges(7,3) * t1030;
t991 = Ifges(7,1) * t1018 + Ifges(7,4) * t1017 + Ifges(7,5) * t1030;
t963 = -mrSges(7,1) * t973 + mrSges(7,3) * t971 + Ifges(7,4) * t988 + Ifges(7,2) * t987 + Ifges(7,6) * t1003 - t1018 * t989 + t1030 * t991;
t990 = Ifges(7,4) * t1018 + Ifges(7,2) * t1017 + Ifges(7,6) * t1030;
t964 = mrSges(7,2) * t973 - mrSges(7,3) * t970 + Ifges(7,1) * t988 + Ifges(7,4) * t987 + Ifges(7,5) * t1003 + t1017 * t989 - t1030 * t990;
t946 = mrSges(6,2) * t979 - mrSges(6,3) * t975 + Ifges(6,1) * t1005 + Ifges(6,4) * t1004 + Ifges(6,5) * t1019 - pkin(13) * t961 + t1008 * t1031 - t1009 * t1043 - t1098 * t963 + t1103 * t964;
t1010 = Ifges(6,1) * t1032 + Ifges(6,4) * t1031 + Ifges(6,5) * t1043;
t1110 = mrSges(7,1) * t970 - mrSges(7,2) * t971 + Ifges(7,5) * t988 + Ifges(7,6) * t987 + Ifges(7,3) * t1003 - t1017 * t991 + t1018 * t990;
t947 = -mrSges(6,1) * t979 + mrSges(6,3) * t976 + Ifges(6,4) * t1005 + Ifges(6,2) * t1004 + Ifges(6,6) * t1019 - pkin(5) * t961 - t1008 * t1032 + t1010 * t1043 - t1110;
t932 = mrSges(5,2) * t985 - mrSges(5,3) * t983 + Ifges(5,1) * t1021 + Ifges(5,4) * t1020 + Ifges(5,5) * t1048 - pkin(12) * t954 + t1024 * t1045 - t1025 * t1058 - t1099 * t947 + t1104 * t946;
t1026 = Ifges(5,1) * t1046 + Ifges(5,4) * t1045 + Ifges(5,5) * t1058;
t1109 = mrSges(6,1) * t975 - mrSges(6,2) * t976 + Ifges(6,5) * t1005 + Ifges(6,6) * t1004 + Ifges(6,3) * t1019 + pkin(5) * t972 + pkin(13) * t962 + t1032 * t1009 - t1031 * t1010 + t1098 * t964 + t1103 * t963;
t933 = -mrSges(5,1) * t985 + mrSges(5,3) * t984 + Ifges(5,4) * t1021 + Ifges(5,2) * t1020 + Ifges(5,6) * t1048 - pkin(4) * t954 - t1046 * t1024 + t1058 * t1026 - t1109;
t1117 = pkin(11) * t945 + t1100 * t932 + t1105 * t933;
t931 = mrSges(5,1) * t983 - mrSges(5,2) * t984 + Ifges(5,5) * t1021 + Ifges(5,6) * t1020 + Ifges(5,3) * t1048 + pkin(4) * t1111 + pkin(12) * t1132 + t1046 * t1025 - t1045 * t1026 + t1099 * t946 + t1104 * t947;
t912 = -mrSges(4,1) * t1027 + mrSges(4,3) * t1013 + Ifges(4,4) * t1055 + Ifges(4,2) * t1054 + Ifges(4,6) * t1073 - pkin(3) * t939 - t1065 * t1049 + t1076 * t1051 - t1091 * t931 + t1095 * t1117;
t1050 = Ifges(4,4) * t1065 + Ifges(4,2) * t1064 + Ifges(4,6) * t1076;
t913 = mrSges(4,2) * t1027 - mrSges(4,3) * t1012 + Ifges(4,1) * t1055 + Ifges(4,4) * t1054 + Ifges(4,5) * t1073 + t1064 * t1049 - t1076 * t1050 - t1100 * t933 + t1105 * t932 + (-t1091 * t939 - t1095 * t940) * pkin(11);
t1118 = pkin(10) * t930 + t1101 * t913 + t1106 * t912;
t1114 = Ifges(3,6) * t1097 + (Ifges(3,4) * t1090 + Ifges(3,2) * t1094) * t1093;
t1070 = t1114 * qJD(1);
t1115 = Ifges(3,5) * t1097 + (Ifges(3,1) * t1090 + Ifges(3,4) * t1094) * t1093;
t1071 = t1115 * qJD(1);
t911 = mrSges(4,1) * t1012 - mrSges(4,2) * t1013 + Ifges(4,5) * t1055 + Ifges(4,6) * t1054 + Ifges(4,3) * t1073 + pkin(3) * t940 + t1065 * t1050 - t1064 * t1051 + t1091 * t1117 + t1095 * t931;
t901 = qJDD(1) * t1153 + mrSges(3,1) * t1059 - mrSges(3,2) * t1060 + pkin(2) * t925 + t1096 * t911 + t1118 * t1092 + (t1129 * qJDD(1) + (t1070 * t1090 - t1071 * t1094) * qJD(1)) * t1093;
t1069 = (t1093 * t1129 + t1153) * qJD(1);
t903 = -mrSges(3,1) * t1066 + mrSges(3,3) * t1060 - pkin(2) * t924 - t1092 * t911 + (-t1069 * t1147 + t1071 * t1097) * qJD(1) + t1118 * t1096 + t1114 * qJDD(1);
t905 = mrSges(3,2) * t1066 - mrSges(3,3) * t1059 - t1101 * t912 + t1106 * t913 + (t1069 * t1142 - t1070 * t1097) * qJD(1) + (-t1092 * t924 - t1096 * t925) * pkin(10) + t1115 * qJDD(1);
t1116 = mrSges(2,1) * t1086 - mrSges(2,2) * t1087 + Ifges(2,3) * qJDD(1) + pkin(1) * t910 + t1097 * t901 + t903 * t1142 + t905 * t1147 + t918 * t1150;
t899 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1086 + Ifges(2,5) * qJDD(1) - t1108 * Ifges(2,6) - t1090 * t903 + t1094 * t905 + (-t1093 * t909 - t1097 * t910) * qJ(2);
t898 = mrSges(2,1) * g(3) + mrSges(2,3) * t1087 + t1108 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t909 - t1093 * t901 + (qJ(2) * t918 + t1090 * t905 + t1094 * t903) * t1097;
t1 = [-m(1) * g(1) + t1131; -m(1) * g(2) + t1151; (-m(1) - m(2)) * g(3) + t909; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(9) * t1151 - t1102 * t898 + t1107 * t899; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(9) * t1131 + t1102 * t899 + t1107 * t898; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1116; t1116; t923; t911; t931; t1109; t1110;];
tauJB  = t1;
