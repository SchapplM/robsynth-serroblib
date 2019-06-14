% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR5
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
% Datum: 2019-05-06 10:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:39:13
% EndTime: 2019-05-06 10:39:30
% DurationCPUTime: 8.52s
% Computational Cost: add. (107004->367), mult. (246617->448), div. (0->0), fcn. (166871->10), ass. (0->155)
t1145 = Ifges(3,1) + Ifges(4,1) + Ifges(5,1);
t1120 = Ifges(3,4) - Ifges(4,5) - Ifges(5,4);
t1144 = Ifges(5,5) - Ifges(3,5) - Ifges(4,4);
t1118 = Ifges(3,6) + Ifges(5,6) - Ifges(4,6);
t1143 = Ifges(4,3) + Ifges(5,2) + Ifges(3,2);
t1142 = Ifges(5,3) + Ifges(3,3) + Ifges(4,2);
t1081 = cos(pkin(6));
t1072 = qJD(1) * t1081 + qJD(2);
t1080 = sin(pkin(6));
t1088 = cos(qJ(2));
t1127 = t1080 * t1088;
t1113 = qJD(1) * t1127;
t1041 = mrSges(4,2) * t1113 + mrSges(4,3) * t1072;
t1071 = qJDD(1) * t1081 + qJDD(2);
t1084 = sin(qJ(2));
t1131 = qJD(1) * t1080;
t1044 = (-pkin(2) * t1088 - qJ(3) * t1084) * t1131;
t1128 = t1080 * t1084;
t1114 = qJD(1) * t1128;
t1140 = t1072 ^ 2;
t1085 = sin(qJ(1));
t1089 = cos(qJ(1));
t1063 = t1085 * g(1) - g(2) * t1089;
t1090 = qJD(1) ^ 2;
t1134 = pkin(8) * t1080;
t1042 = qJDD(1) * pkin(1) + t1090 * t1134 + t1063;
t1064 = -g(1) * t1089 - g(2) * t1085;
t1121 = qJDD(1) * t1080;
t1043 = -pkin(1) * t1090 + pkin(8) * t1121 + t1064;
t1125 = t1081 * t1088;
t990 = -g(3) * t1127 + t1042 * t1125 - t1084 * t1043;
t973 = -t1071 * pkin(2) - qJ(3) * t1140 + t1044 * t1114 + qJDD(3) - t990;
t1141 = -m(4) * t973 + t1071 * mrSges(4,1) + t1072 * t1041;
t1139 = -2 * qJD(4);
t1137 = -mrSges(5,2) - mrSges(4,3);
t1136 = mrSges(3,3) + mrSges(4,2);
t1135 = pkin(2) * t1072;
t1126 = t1081 * t1084;
t1046 = (mrSges(5,1) * t1088 + mrSges(5,2) * t1084) * t1131;
t1047 = (-mrSges(3,1) * t1088 + mrSges(3,2) * t1084) * t1131;
t1050 = -qJD(2) * t1114 + t1088 * t1121;
t1038 = -mrSges(4,1) * t1072 + mrSges(4,2) * t1114;
t1045 = (-mrSges(4,1) * t1088 - mrSges(4,3) * t1084) * t1131;
t1083 = sin(qJ(5));
t1087 = cos(qJ(5));
t1021 = -t1072 * t1087 - t1083 * t1114;
t1033 = qJDD(5) + t1050;
t1056 = qJD(5) + t1113;
t1082 = sin(qJ(6));
t1086 = cos(qJ(6));
t1019 = qJD(6) - t1021;
t1054 = t1056 ^ 2;
t1130 = qJD(1) * t1088;
t1049 = (qJD(2) * t1130 + qJDD(1) * t1084) * t1080;
t1035 = -pkin(3) * t1072 - qJ(4) * t1114;
t1007 = -t1081 * g(3) - t1080 * t1042;
t1106 = t1072 * t1113;
t1105 = -t1050 * pkin(2) + t1007 + (-t1049 - t1106) * qJ(3);
t1107 = 0.2e1 * t1114;
t1129 = t1080 ^ 2 * t1090;
t1112 = t1088 ^ 2 * t1129;
t1094 = -qJ(4) * t1112 + qJD(3) * t1107 + t1035 * t1114 + qJDD(4) - t1105;
t953 = -t1049 * pkin(9) + (pkin(3) + pkin(4)) * t1050 + (-pkin(9) * t1088 + (-pkin(2) - pkin(4)) * t1084) * t1072 * t1131 + t1094;
t1048 = (pkin(4) * t1088 - pkin(9) * t1084) * t1131;
t1124 = t1042 * t1126 + t1088 * t1043;
t1102 = -pkin(2) * t1140 + t1071 * qJ(3) + 0.2e1 * qJD(3) * t1072 + t1044 * t1113 + t1124;
t1092 = -pkin(3) * t1112 - t1050 * qJ(4) + t1072 * t1035 + t1113 * t1139 + t1102;
t955 = -t1140 * pkin(4) - t1071 * pkin(9) + (-g(3) * t1084 - t1048 * t1130) * t1080 + t1092;
t950 = t1083 * t953 + t1087 * t955;
t1022 = -t1072 * t1083 + t1087 * t1114;
t993 = -pkin(5) * t1021 - pkin(10) * t1022;
t947 = -pkin(5) * t1054 + pkin(10) * t1033 + t1021 * t993 + t950;
t1096 = -t1049 * qJ(4) + t973 + (-t1084 * t1088 * t1129 - t1071) * pkin(3);
t957 = t1071 * pkin(4) - pkin(9) * t1140 - qJ(4) * t1106 + qJD(4) * t1107 + t1048 * t1114 - t1096;
t988 = -qJD(5) * t1022 - t1049 * t1083 - t1071 * t1087;
t989 = qJD(5) * t1021 + t1049 * t1087 - t1071 * t1083;
t951 = t957 + (t1022 * t1056 - t988) * pkin(5) + (-t1021 * t1056 - t989) * pkin(10);
t944 = -t1082 * t947 + t1086 * t951;
t994 = -t1022 * t1082 + t1056 * t1086;
t964 = t994 * qJD(6) + t1033 * t1082 + t1086 * t989;
t995 = t1022 * t1086 + t1056 * t1082;
t974 = -mrSges(7,1) * t994 + mrSges(7,2) * t995;
t976 = -mrSges(7,2) * t1019 + t994 * mrSges(7,3);
t985 = qJDD(6) - t988;
t941 = m(7) * t944 + t985 * mrSges(7,1) - t964 * mrSges(7,3) + t1019 * t976 - t995 * t974;
t945 = t1082 * t951 + t1086 * t947;
t963 = -t995 * qJD(6) + t1033 * t1086 - t1082 * t989;
t977 = mrSges(7,1) * t1019 - t995 * mrSges(7,3);
t942 = m(7) * t945 - t985 * mrSges(7,2) + t963 * mrSges(7,3) - t1019 * t977 + t994 * t974;
t1110 = -t1082 * t941 + t1086 * t942;
t992 = -mrSges(6,1) * t1021 + mrSges(6,2) * t1022;
t997 = mrSges(6,1) * t1056 - mrSges(6,3) * t1022;
t930 = m(6) * t950 - mrSges(6,2) * t1033 + t988 * mrSges(6,3) + t1021 * t992 - t1056 * t997 + t1110;
t949 = -t1083 * t955 + t1087 * t953;
t946 = -pkin(5) * t1033 - pkin(10) * t1054 + t1022 * t993 - t949;
t1100 = -m(7) * t946 + t963 * mrSges(7,1) - t964 * mrSges(7,2) + t994 * t976 - t995 * t977;
t996 = -mrSges(6,2) * t1056 + mrSges(6,3) * t1021;
t937 = m(6) * t949 + mrSges(6,1) * t1033 - t989 * mrSges(6,3) - t1022 * t992 + t1056 * t996 + t1100;
t1109 = -t1083 * t937 + t1087 * t930;
t1117 = g(3) * t1128;
t959 = t1092 - t1117;
t1104 = m(5) * t959 - t1050 * mrSges(5,3) + t1109;
t967 = t1102 - t1117;
t1098 = m(4) * t967 + t1071 * mrSges(4,3) + t1072 * t1038 + t1045 * t1113 + t1104;
t1036 = -mrSges(5,1) * t1072 - mrSges(5,3) * t1114;
t1123 = -mrSges(3,1) * t1072 + mrSges(3,3) * t1114 + t1036;
t991 = -t1117 + t1124;
t918 = m(3) * t991 + t1123 * t1072 + (-mrSges(3,2) + mrSges(5,2)) * t1071 + t1136 * t1050 + (-t1046 + t1047) * t1113 + t1098;
t923 = t1083 * t930 + t1087 * t937;
t961 = t1050 * pkin(3) - t1114 * t1135 + t1094;
t1101 = m(5) * t961 + t1050 * mrSges(5,1) + t923;
t968 = (-0.2e1 * qJD(3) + t1135) * t1114 + t1105;
t1097 = m(4) * t968 - t1050 * mrSges(4,1) - t1101;
t1039 = mrSges(5,2) * t1072 - mrSges(5,3) * t1113;
t1122 = -mrSges(3,2) * t1072 + mrSges(3,3) * t1113 + t1039;
t920 = m(3) * t1007 - t1050 * mrSges(3,1) + (mrSges(3,2) + t1137) * t1049 + ((-t1041 - t1122) * t1088 + (-t1038 - t1123) * t1084) * t1131 + t1097;
t932 = t1082 * t942 + t1086 * t941;
t1103 = -m(6) * t957 + t988 * mrSges(6,1) - t989 * mrSges(6,2) + t1021 * t996 - t1022 * t997 - t932;
t1132 = qJ(4) * t1088;
t960 = (t1072 * t1132 + t1084 * t1139) * t1131 + t1096;
t1095 = -m(5) * t960 + t1049 * mrSges(5,3) + t1046 * t1114 - t1103;
t926 = (mrSges(3,1) + mrSges(5,1)) * t1071 - t1136 * t1049 + m(3) * t990 + t1122 * t1072 + t1095 + (-t1045 - t1047) * t1114 + t1141;
t907 = -t1080 * t920 + t926 * t1125 + t918 * t1126;
t904 = m(2) * t1063 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1090 + t907;
t912 = -t1084 * t926 + t1088 * t918;
t910 = m(2) * t1064 - mrSges(2,1) * t1090 - qJDD(1) * mrSges(2,2) + t912;
t1133 = t1085 * t910 + t1089 * t904;
t906 = t1081 * t920 + t926 * t1127 + t918 * t1128;
t1116 = (-t1144 * t1084 + t1118 * t1088) * t1131 + t1142 * t1072;
t1115 = (-t1120 * t1084 - t1143 * t1088) * t1131 - t1118 * t1072;
t1111 = (-t1145 * t1084 - t1120 * t1088) * t1131 + t1144 * t1072;
t1108 = -t1085 * t904 + t1089 * t910;
t969 = Ifges(7,5) * t995 + Ifges(7,6) * t994 + Ifges(7,3) * t1019;
t971 = Ifges(7,1) * t995 + Ifges(7,4) * t994 + Ifges(7,5) * t1019;
t935 = -mrSges(7,1) * t946 + mrSges(7,3) * t945 + Ifges(7,4) * t964 + Ifges(7,2) * t963 + Ifges(7,6) * t985 + t1019 * t971 - t995 * t969;
t970 = Ifges(7,4) * t995 + Ifges(7,2) * t994 + Ifges(7,6) * t1019;
t936 = mrSges(7,2) * t946 - mrSges(7,3) * t944 + Ifges(7,1) * t964 + Ifges(7,4) * t963 + Ifges(7,5) * t985 - t1019 * t970 + t994 * t969;
t978 = Ifges(6,5) * t1022 + Ifges(6,6) * t1021 + Ifges(6,3) * t1056;
t979 = Ifges(6,4) * t1022 + Ifges(6,2) * t1021 + Ifges(6,6) * t1056;
t913 = mrSges(6,2) * t957 - mrSges(6,3) * t949 + Ifges(6,1) * t989 + Ifges(6,4) * t988 + Ifges(6,5) * t1033 - pkin(10) * t932 + t1021 * t978 - t1056 * t979 - t1082 * t935 + t1086 * t936;
t1091 = mrSges(7,1) * t944 - mrSges(7,2) * t945 + Ifges(7,5) * t964 + Ifges(7,6) * t963 + Ifges(7,3) * t985 + t995 * t970 - t994 * t971;
t980 = Ifges(6,1) * t1022 + Ifges(6,4) * t1021 + Ifges(6,5) * t1056;
t914 = -mrSges(6,1) * t957 + mrSges(6,3) * t950 + Ifges(6,4) * t989 + Ifges(6,2) * t988 + Ifges(6,6) * t1033 - pkin(5) * t932 - t1022 * t978 + t1056 * t980 - t1091;
t928 = -t1071 * mrSges(5,1) - t1072 * t1039 - t1095;
t927 = t1049 * mrSges(4,2) + t1045 * t1114 - t1141 + t928;
t898 = -pkin(3) * t928 - pkin(2) * t927 - mrSges(4,1) * t973 - t1087 * t914 + qJ(3) * (t1072 * t1036 + t1098) - t1083 * t913 - pkin(9) * t1109 - pkin(4) * t1103 + mrSges(3,1) * t990 - mrSges(3,2) * t991 + mrSges(5,2) * t959 - mrSges(5,1) * t960 + mrSges(4,3) * t967 + (mrSges(5,2) * qJ(3) + t1142) * t1071 + (mrSges(4,2) * qJ(3) + t1118) * t1050 - t1144 * t1049 + (-t1115 * t1084 + (-qJ(3) * t1046 + t1111) * t1088) * t1131;
t1093 = mrSges(6,1) * t949 - mrSges(6,2) * t950 + Ifges(6,5) * t989 + Ifges(6,6) * t988 + Ifges(6,3) * t1033 + pkin(5) * t1100 + pkin(10) * t1110 - t1021 * t980 + t1022 * t979 + t1082 * t936 + t1086 * t935;
t921 = t1137 * t1049 + ((-t1039 - t1041) * t1088 + (-t1036 - t1038) * t1084) * t1131 + t1097;
t922 = t1049 * mrSges(5,2) + (t1036 * t1084 + t1039 * t1088) * t1131 + t1101;
t900 = t1093 + pkin(4) * t923 + pkin(3) * t922 - mrSges(3,1) * t1007 + mrSges(3,3) * t991 - pkin(2) * t921 - mrSges(5,3) * t959 + mrSges(5,1) * t961 + mrSges(4,2) * t967 - mrSges(4,1) * t968 - qJ(4) * t1104 + (t1046 * t1132 - t1084 * t1116) * t1131 + (-qJ(4) * t1036 - t1111) * t1072 + (-mrSges(5,2) * qJ(4) + t1118) * t1071 + t1143 * t1050 + t1120 * t1049;
t902 = mrSges(3,2) * t1007 + mrSges(4,2) * t973 + mrSges(5,2) * t961 - mrSges(3,3) * t990 - mrSges(4,3) * t968 - mrSges(5,3) * t960 - pkin(9) * t923 - qJ(3) * t921 - qJ(4) * t928 - t1083 * t914 + t1087 * t913 + t1115 * t1072 - t1144 * t1071 + t1120 * t1050 + t1145 * t1049 + t1116 * t1113;
t1099 = mrSges(2,1) * t1063 - mrSges(2,2) * t1064 + Ifges(2,3) * qJDD(1) + pkin(1) * t907 + t1081 * t898 + t900 * t1127 + t902 * t1128 + t912 * t1134;
t896 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1063 + Ifges(2,5) * qJDD(1) - t1090 * Ifges(2,6) - t1084 * t900 + t1088 * t902 + (-t1080 * t906 - t1081 * t907) * pkin(8);
t895 = mrSges(2,1) * g(3) + mrSges(2,3) * t1064 + t1090 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t906 - t1080 * t898 + (pkin(8) * t912 + t1084 * t902 + t1088 * t900) * t1081;
t1 = [-m(1) * g(1) + t1108; -m(1) * g(2) + t1133; (-m(1) - m(2)) * g(3) + t906; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1133 - t1085 * t895 + t1089 * t896; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1108 + t1085 * t896 + t1089 * t895; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1099; t1099; t898; t927; t922; t1093; t1091;];
tauJB  = t1;
