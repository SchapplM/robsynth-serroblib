% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-09 12:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-09 11:35:59
% EndTime: 2019-05-09 11:42:41
% DurationCPUTime: 409.22s
% Computational Cost: add. (6620790->433), mult. (17956716->588), div. (0->0), fcn. (15481279->18), ass. (0->197)
t1089 = sin(pkin(8));
t1093 = cos(pkin(8));
t1098 = sin(qJ(4));
t1103 = cos(qJ(4));
t1095 = cos(pkin(6));
t1085 = qJD(1) * t1095 + qJD(2);
t1090 = sin(pkin(7));
t1094 = cos(pkin(7));
t1091 = sin(pkin(6));
t1104 = cos(qJ(2));
t1130 = t1091 * t1104;
t1121 = qJD(1) * t1130;
t1067 = (t1085 * t1090 + t1094 * t1121) * qJ(3);
t1099 = sin(qJ(2));
t1138 = qJD(1) * t1091;
t1141 = qJ(3) * t1090;
t1073 = (-pkin(2) * t1104 - t1099 * t1141) * t1138;
t1136 = qJD(1) * t1104;
t1079 = (qJD(2) * t1136 + qJDD(1) * t1099) * t1091;
t1084 = qJDD(1) * t1095 + qJDD(2);
t1100 = sin(qJ(1));
t1105 = cos(qJ(1));
t1082 = t1100 * g(1) - g(2) * t1105;
t1106 = qJD(1) ^ 2;
t1146 = pkin(10) * t1091;
t1076 = qJDD(1) * pkin(1) + t1106 * t1146 + t1082;
t1083 = -g(1) * t1105 - g(2) * t1100;
t1077 = -pkin(1) * t1106 + qJDD(1) * t1146 + t1083;
t1125 = t1095 * t1104;
t1118 = t1076 * t1125 - t1099 * t1077;
t1137 = qJD(1) * t1099;
t1140 = qJ(3) * t1094;
t1032 = -t1079 * t1140 + t1084 * pkin(2) + t1085 * t1067 + (-g(3) * t1104 - t1073 * t1137) * t1091 + t1118;
t1131 = t1091 * t1099;
t1122 = qJD(1) * t1131;
t1072 = pkin(2) * t1085 - t1122 * t1140;
t1080 = (-qJD(2) * t1137 + qJDD(1) * t1104) * t1091;
t1113 = t1080 * t1094 + t1084 * t1090;
t1126 = t1095 * t1099;
t1123 = t1076 * t1126 + t1104 * t1077;
t1033 = -t1085 * t1072 + (-g(3) * t1099 + t1073 * t1136) * t1091 + t1113 * qJ(3) + t1123;
t1143 = t1095 * g(3);
t1037 = -t1079 * t1141 - t1080 * pkin(2) - t1143 + (-t1076 + (-t1067 * t1104 + t1072 * t1099) * qJD(1)) * t1091;
t1088 = sin(pkin(14));
t1092 = cos(pkin(14));
t1127 = t1094 * t1104;
t1135 = t1088 * t1090;
t1061 = t1085 * t1135 + (t1088 * t1127 + t1092 * t1099) * t1138;
t1129 = t1092 * t1094;
t1132 = t1090 * t1092;
t1006 = -0.2e1 * qJD(3) * t1061 + t1032 * t1129 - t1033 * t1088 + t1037 * t1132;
t1060 = t1085 * t1132 + (-t1088 * t1099 + t1092 * t1127) * t1138;
t1145 = pkin(11) * t1089;
t1046 = -pkin(3) * t1060 - t1061 * t1145;
t1070 = t1085 * t1094 - t1090 * t1121;
t1114 = t1060 * t1093 + t1070 * t1089;
t1049 = t1114 * pkin(11);
t1055 = t1092 * t1079 + t1113 * t1088;
t1062 = -t1080 * t1090 + t1084 * t1094;
t1144 = pkin(11) * t1093;
t990 = pkin(3) * t1062 - t1046 * t1061 + t1049 * t1070 - t1055 * t1144 + t1006;
t1134 = t1088 * t1094;
t1007 = 0.2e1 * qJD(3) * t1060 + t1032 * t1134 + t1092 * t1033 + t1037 * t1135;
t1051 = pkin(3) * t1070 - t1061 * t1144;
t1054 = -t1088 * t1079 + t1113 * t1092;
t1115 = t1054 * t1093 + t1062 * t1089;
t991 = t1115 * pkin(11) + t1060 * t1046 - t1070 * t1051 + t1007;
t1021 = -t1032 * t1090 + t1094 * t1037 + qJDD(3);
t995 = -pkin(3) * t1054 - t1049 * t1060 + t1051 * t1061 - t1055 * t1145 + t1021;
t977 = -t1098 * t991 + (t1089 * t995 + t1093 * t990) * t1103;
t1040 = -t1098 * t1061 + t1114 * t1103;
t1041 = t1103 * t1061 + t1114 * t1098;
t1019 = -t1041 * qJD(4) - t1098 * t1055 + t1115 * t1103;
t1056 = -g(3) * t1130 + t1118;
t1075 = -mrSges(3,2) * t1085 + mrSges(3,3) * t1121;
t1078 = (-mrSges(3,1) * t1104 + mrSges(3,2) * t1099) * t1138;
t1047 = -mrSges(4,1) * t1060 + mrSges(4,2) * t1061;
t1052 = -mrSges(4,2) * t1070 + mrSges(4,3) * t1060;
t1128 = t1093 * t1098;
t1020 = t1040 * qJD(4) + t1103 * t1055 + t1115 * t1098;
t1022 = -mrSges(5,1) * t1040 + mrSges(5,2) * t1041;
t1050 = -t1060 * t1089 + t1070 * t1093 + qJD(4);
t1027 = -mrSges(5,2) * t1050 + mrSges(5,3) * t1040;
t1045 = -t1054 * t1089 + t1062 * t1093 + qJDD(4);
t1097 = sin(qJ(5));
t1102 = cos(qJ(5));
t1025 = -t1041 * t1097 + t1050 * t1102;
t1039 = qJD(5) - t1040;
t1013 = -mrSges(6,2) * t1039 + mrSges(6,3) * t1025;
t1026 = t1041 * t1102 + t1050 * t1097;
t1014 = mrSges(6,1) * t1039 - mrSges(6,3) * t1026;
t1096 = sin(qJ(6));
t1101 = cos(qJ(6));
t1011 = -t1026 * t1096 + t1039 * t1101;
t1024 = qJD(6) - t1025;
t1000 = -mrSges(7,2) * t1024 + mrSges(7,3) * t1011;
t1012 = t1026 * t1101 + t1039 * t1096;
t1009 = -pkin(5) * t1025 - pkin(13) * t1026;
t1018 = qJDD(5) - t1019;
t1038 = t1039 ^ 2;
t1023 = -pkin(4) * t1040 - pkin(12) * t1041;
t1048 = t1050 ^ 2;
t1133 = t1089 * t1098;
t978 = t1103 * t991 + t990 * t1128 + t995 * t1133;
t974 = -pkin(4) * t1048 + pkin(12) * t1045 + t1023 * t1040 + t978;
t979 = -t1089 * t990 + t1093 * t995;
t976 = (-t1040 * t1050 - t1020) * pkin(12) + (t1041 * t1050 - t1019) * pkin(4) + t979;
t970 = t1097 * t976 + t1102 * t974;
t968 = -pkin(5) * t1038 + pkin(13) * t1018 + t1009 * t1025 + t970;
t973 = -t1045 * pkin(4) - t1048 * pkin(12) + t1041 * t1023 - t977;
t998 = -qJD(5) * t1026 - t1020 * t1097 + t1045 * t1102;
t999 = qJD(5) * t1025 + t1020 * t1102 + t1045 * t1097;
t971 = (-t1025 * t1039 - t999) * pkin(13) + (t1026 * t1039 - t998) * pkin(5) + t973;
t964 = -t1096 * t968 + t1101 * t971;
t982 = qJD(6) * t1011 + t1018 * t1096 + t1101 * t999;
t992 = -mrSges(7,1) * t1011 + mrSges(7,2) * t1012;
t997 = qJDD(6) - t998;
t962 = m(7) * t964 + mrSges(7,1) * t997 - mrSges(7,3) * t982 + t1000 * t1024 - t1012 * t992;
t1001 = mrSges(7,1) * t1024 - mrSges(7,3) * t1012;
t965 = t1096 * t971 + t1101 * t968;
t981 = -qJD(6) * t1012 + t1018 * t1101 - t1096 * t999;
t963 = m(7) * t965 - mrSges(7,2) * t997 + mrSges(7,3) * t981 - t1001 * t1024 + t1011 * t992;
t955 = t1096 * t963 + t1101 * t962;
t1109 = -m(6) * t973 + t998 * mrSges(6,1) - mrSges(6,2) * t999 + t1025 * t1013 - t1014 * t1026 - t955;
t951 = m(5) * t977 + mrSges(5,1) * t1045 - mrSges(5,3) * t1020 - t1022 * t1041 + t1027 * t1050 + t1109;
t1139 = t1103 * t951;
t1028 = mrSges(5,1) * t1050 - mrSges(5,3) * t1041;
t1008 = -mrSges(6,1) * t1025 + mrSges(6,2) * t1026;
t956 = -t1096 * t962 + t1101 * t963;
t954 = m(6) * t970 - mrSges(6,2) * t1018 + mrSges(6,3) * t998 + t1008 * t1025 - t1014 * t1039 + t956;
t969 = -t1097 * t974 + t1102 * t976;
t967 = -pkin(5) * t1018 - pkin(13) * t1038 + t1009 * t1026 - t969;
t966 = -m(7) * t967 + t981 * mrSges(7,1) - mrSges(7,2) * t982 + t1011 * t1000 - t1001 * t1012;
t960 = m(6) * t969 + mrSges(6,1) * t1018 - mrSges(6,3) * t999 - t1008 * t1026 + t1013 * t1039 + t966;
t1120 = -t1097 * t960 + t1102 * t954;
t945 = m(5) * t978 - mrSges(5,2) * t1045 + mrSges(5,3) * t1019 + t1022 * t1040 - t1028 * t1050 + t1120;
t948 = t1097 * t954 + t1102 * t960;
t947 = m(5) * t979 - mrSges(5,1) * t1019 + mrSges(5,2) * t1020 - t1027 * t1040 + t1028 * t1041 + t948;
t934 = -t1089 * t947 + t1093 * t1139 + t945 * t1128;
t930 = m(4) * t1006 + mrSges(4,1) * t1062 - mrSges(4,3) * t1055 - t1047 * t1061 + t1052 * t1070 + t934;
t1053 = mrSges(4,1) * t1070 - mrSges(4,3) * t1061;
t933 = t1089 * t1139 + t1093 * t947 + t945 * t1133;
t932 = m(4) * t1021 - mrSges(4,1) * t1054 + mrSges(4,2) * t1055 - t1052 * t1060 + t1053 * t1061 + t933;
t939 = -t1098 * t951 + t1103 * t945;
t938 = m(4) * t1007 - mrSges(4,2) * t1062 + mrSges(4,3) * t1054 + t1047 * t1060 - t1053 * t1070 + t939;
t919 = -t1090 * t932 + t930 * t1129 + t938 * t1134;
t915 = m(3) * t1056 + mrSges(3,1) * t1084 - mrSges(3,3) * t1079 + t1075 * t1085 - t1078 * t1122 + t919;
t1066 = -t1091 * t1076 - t1143;
t1074 = mrSges(3,1) * t1085 - mrSges(3,3) * t1122;
t918 = t1094 * t932 + t930 * t1132 + t938 * t1135;
t917 = m(3) * t1066 - t1080 * mrSges(3,1) + t1079 * mrSges(3,2) + (t1074 * t1099 - t1075 * t1104) * t1138 + t918;
t1057 = -g(3) * t1131 + t1123;
t924 = -t1088 * t930 + t1092 * t938;
t923 = m(3) * t1057 - mrSges(3,2) * t1084 + mrSges(3,3) * t1080 - t1074 * t1085 + t1078 * t1121 + t924;
t904 = -t1091 * t917 + t915 * t1125 + t923 * t1126;
t901 = m(2) * t1082 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1106 + t904;
t912 = -t1099 * t915 + t1104 * t923;
t910 = m(2) * t1083 - mrSges(2,1) * t1106 - qJDD(1) * mrSges(2,2) + t912;
t1142 = t1100 * t910 + t1105 * t901;
t903 = t1095 * t917 + t915 * t1130 + t923 * t1131;
t1119 = -t1100 * t901 + t1105 * t910;
t1015 = Ifges(5,5) * t1041 + Ifges(5,6) * t1040 + Ifges(5,3) * t1050;
t1016 = Ifges(5,4) * t1041 + Ifges(5,2) * t1040 + Ifges(5,6) * t1050;
t1002 = Ifges(6,5) * t1026 + Ifges(6,6) * t1025 + Ifges(6,3) * t1039;
t1003 = Ifges(6,4) * t1026 + Ifges(6,2) * t1025 + Ifges(6,6) * t1039;
t983 = Ifges(7,5) * t1012 + Ifges(7,6) * t1011 + Ifges(7,3) * t1024;
t985 = Ifges(7,1) * t1012 + Ifges(7,4) * t1011 + Ifges(7,5) * t1024;
t957 = -mrSges(7,1) * t967 + mrSges(7,3) * t965 + Ifges(7,4) * t982 + Ifges(7,2) * t981 + Ifges(7,6) * t997 - t1012 * t983 + t1024 * t985;
t984 = Ifges(7,4) * t1012 + Ifges(7,2) * t1011 + Ifges(7,6) * t1024;
t958 = mrSges(7,2) * t967 - mrSges(7,3) * t964 + Ifges(7,1) * t982 + Ifges(7,4) * t981 + Ifges(7,5) * t997 + t1011 * t983 - t1024 * t984;
t940 = mrSges(6,2) * t973 - mrSges(6,3) * t969 + Ifges(6,1) * t999 + Ifges(6,4) * t998 + Ifges(6,5) * t1018 - pkin(13) * t955 + t1002 * t1025 - t1003 * t1039 - t1096 * t957 + t1101 * t958;
t1004 = Ifges(6,1) * t1026 + Ifges(6,4) * t1025 + Ifges(6,5) * t1039;
t1108 = mrSges(7,1) * t964 - mrSges(7,2) * t965 + Ifges(7,5) * t982 + Ifges(7,6) * t981 + Ifges(7,3) * t997 - t1011 * t985 + t1012 * t984;
t941 = -mrSges(6,1) * t973 + mrSges(6,3) * t970 + Ifges(6,4) * t999 + Ifges(6,2) * t998 + Ifges(6,6) * t1018 - pkin(5) * t955 - t1002 * t1026 + t1004 * t1039 - t1108;
t926 = mrSges(5,2) * t979 - mrSges(5,3) * t977 + Ifges(5,1) * t1020 + Ifges(5,4) * t1019 + Ifges(5,5) * t1045 - pkin(12) * t948 + t1015 * t1040 - t1016 * t1050 - t1097 * t941 + t1102 * t940;
t1017 = Ifges(5,1) * t1041 + Ifges(5,4) * t1040 + Ifges(5,5) * t1050;
t1107 = mrSges(6,1) * t969 - mrSges(6,2) * t970 + Ifges(6,5) * t999 + Ifges(6,6) * t998 + Ifges(6,3) * t1018 + pkin(5) * t966 + pkin(13) * t956 + t1026 * t1003 - t1025 * t1004 + t1096 * t958 + t1101 * t957;
t927 = -mrSges(5,1) * t979 + mrSges(5,3) * t978 + Ifges(5,4) * t1020 + Ifges(5,2) * t1019 + Ifges(5,6) * t1045 - pkin(4) * t948 - t1041 * t1015 + t1050 * t1017 - t1107;
t1112 = pkin(11) * t939 + t1098 * t926 + t1103 * t927;
t1042 = Ifges(4,5) * t1061 + Ifges(4,6) * t1060 + Ifges(4,3) * t1070;
t1044 = Ifges(4,1) * t1061 + Ifges(4,4) * t1060 + Ifges(4,5) * t1070;
t925 = mrSges(5,1) * t977 - mrSges(5,2) * t978 + Ifges(5,5) * t1020 + Ifges(5,6) * t1019 + Ifges(5,3) * t1045 + pkin(4) * t1109 + pkin(12) * t1120 + t1041 * t1016 - t1040 * t1017 + t1097 * t940 + t1102 * t941;
t906 = -mrSges(4,1) * t1021 + mrSges(4,3) * t1007 + Ifges(4,4) * t1055 + Ifges(4,2) * t1054 + Ifges(4,6) * t1062 - pkin(3) * t933 - t1061 * t1042 + t1070 * t1044 - t1089 * t925 + t1112 * t1093;
t1043 = Ifges(4,4) * t1061 + Ifges(4,2) * t1060 + Ifges(4,6) * t1070;
t907 = mrSges(4,2) * t1021 - mrSges(4,3) * t1006 + Ifges(4,1) * t1055 + Ifges(4,4) * t1054 + Ifges(4,5) * t1062 + t1060 * t1042 - t1070 * t1043 - t1098 * t927 + t1103 * t926 + (-t1089 * t933 - t1093 * t934) * pkin(11);
t1111 = qJ(3) * t924 + t1088 * t907 + t1092 * t906;
t1064 = Ifges(3,6) * t1085 + (Ifges(3,4) * t1099 + Ifges(3,2) * t1104) * t1138;
t1065 = Ifges(3,5) * t1085 + (Ifges(3,1) * t1099 + Ifges(3,4) * t1104) * t1138;
t905 = mrSges(4,1) * t1006 - mrSges(4,2) * t1007 + Ifges(4,5) * t1055 + Ifges(4,6) * t1054 + Ifges(4,3) * t1062 + pkin(3) * t934 + t1061 * t1043 - t1060 * t1044 + t1112 * t1089 + t1093 * t925;
t895 = mrSges(3,1) * t1056 - mrSges(3,2) * t1057 + Ifges(3,5) * t1079 + Ifges(3,6) * t1080 + Ifges(3,3) * t1084 + pkin(2) * t919 + t1094 * t905 + (t1064 * t1099 - t1065 * t1104) * t1138 + t1111 * t1090;
t1063 = Ifges(3,3) * t1085 + (Ifges(3,5) * t1099 + Ifges(3,6) * t1104) * t1138;
t897 = -mrSges(3,1) * t1066 + mrSges(3,3) * t1057 + Ifges(3,4) * t1079 + Ifges(3,2) * t1080 + Ifges(3,6) * t1084 - pkin(2) * t918 - t1063 * t1122 + t1085 * t1065 - t1090 * t905 + t1111 * t1094;
t899 = t1063 * t1121 + mrSges(3,2) * t1066 - mrSges(3,3) * t1056 + Ifges(3,1) * t1079 + Ifges(3,4) * t1080 + Ifges(3,5) * t1084 - t1085 * t1064 - t1088 * t906 + t1092 * t907 + (-t1090 * t918 - t1094 * t919) * qJ(3);
t1110 = mrSges(2,1) * t1082 - mrSges(2,2) * t1083 + Ifges(2,3) * qJDD(1) + pkin(1) * t904 + t1095 * t895 + t897 * t1130 + t899 * t1131 + t912 * t1146;
t893 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1082 + Ifges(2,5) * qJDD(1) - t1106 * Ifges(2,6) - t1099 * t897 + t1104 * t899 + (-t1091 * t903 - t1095 * t904) * pkin(10);
t892 = mrSges(2,1) * g(3) + mrSges(2,3) * t1083 + t1106 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t903 - t1091 * t895 + (pkin(10) * t912 + t1099 * t899 + t1104 * t897) * t1095;
t1 = [-m(1) * g(1) + t1119; -m(1) * g(2) + t1142; (-m(1) - m(2)) * g(3) + t903; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(9) * t1142 - t1100 * t892 + t1105 * t893; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(9) * t1119 + t1100 * t893 + t1105 * t892; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1110; t1110; t895; t932; t925; t1107; t1108;];
tauJB  = t1;
