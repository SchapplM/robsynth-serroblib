% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-05-07 06:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:10:18
% EndTime: 2019-05-07 06:10:37
% DurationCPUTime: 9.88s
% Computational Cost: add. (141326->362), mult. (305109->436), div. (0->0), fcn. (227439->10), ass. (0->154)
t1066 = cos(pkin(6));
t1061 = qJD(1) * t1066 + qJD(2);
t1068 = sin(qJ(3));
t1065 = sin(pkin(6));
t1069 = sin(qJ(2));
t1106 = t1065 * t1069;
t1092 = qJD(1) * t1106;
t1115 = cos(qJ(3));
t1030 = t1068 * t1061 + t1115 * t1092;
t1072 = cos(qJ(2));
t1105 = t1065 * t1072;
t1091 = qJD(1) * t1105;
t1052 = -qJD(3) + t1091;
t1012 = pkin(4) * t1052 - qJ(5) * t1030;
t1125 = (2 * qJD(4)) + t1012;
t1124 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t1099 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t1098 = Ifges(5,5) - Ifges(6,4) - Ifges(4,4);
t1097 = Ifges(5,6) - Ifges(6,5) - Ifges(4,6);
t1123 = Ifges(5,3) + Ifges(6,1) + Ifges(4,2);
t1122 = Ifges(6,3) + Ifges(4,3) + Ifges(5,2);
t1029 = -t1115 * t1061 + t1068 * t1092;
t1001 = mrSges(5,1) * t1029 - mrSges(5,3) * t1030;
t1015 = -mrSges(6,2) * t1052 - mrSges(6,3) * t1030;
t1100 = qJDD(1) * t1065;
t1045 = -qJD(2) * t1092 + t1072 * t1100;
t1037 = -qJDD(3) + t1045;
t1067 = sin(qJ(6));
t1071 = cos(qJ(6));
t1014 = mrSges(5,1) * t1052 + mrSges(5,2) * t1030;
t1008 = -t1029 * t1067 + t1052 * t1071;
t1009 = t1029 * t1071 + t1052 * t1067;
t1003 = pkin(5) * t1030 - pkin(10) * t1029;
t1119 = -2 * qJD(4);
t1046 = t1052 * t1119;
t1028 = t1029 ^ 2;
t1000 = pkin(3) * t1029 - qJ(4) * t1030;
t1120 = t1052 ^ 2;
t1109 = qJD(1) * t1065;
t1043 = (-pkin(2) * t1072 - pkin(9) * t1069) * t1109;
t1059 = t1061 ^ 2;
t1060 = qJDD(1) * t1066 + qJDD(2);
t1070 = sin(qJ(1));
t1073 = cos(qJ(1));
t1054 = t1070 * g(1) - g(2) * t1073;
t1074 = qJD(1) ^ 2;
t1112 = pkin(8) * t1065;
t1040 = qJDD(1) * pkin(1) + t1074 * t1112 + t1054;
t1055 = -g(1) * t1073 - g(2) * t1070;
t1041 = -pkin(1) * t1074 + pkin(8) * t1100 + t1055;
t1104 = t1066 * t1069;
t1101 = t1040 * t1104 + t1072 * t1041;
t1108 = qJD(1) * t1072;
t964 = -t1059 * pkin(2) + t1060 * pkin(9) + (-g(3) * t1069 + t1043 * t1108) * t1065 + t1101;
t1044 = (qJD(2) * t1108 + qJDD(1) * t1069) * t1065;
t1111 = t1066 * g(3);
t965 = -t1045 * pkin(2) - t1044 * pkin(9) - t1111 + (-t1040 + (pkin(2) * t1069 - pkin(9) * t1072) * t1061 * qJD(1)) * t1065;
t949 = t1068 * t965 + t1115 * t964;
t1088 = pkin(3) * t1120 + t1037 * qJ(4) + t1029 * t1000 - t949;
t995 = t1030 * qJD(3) + t1068 * t1044 - t1115 * t1060;
t1082 = t1028 * pkin(4) - t995 * qJ(5) + t1088;
t939 = -t1037 * pkin(5) - t1120 * pkin(10) - t1052 * t1012 + t1046 + ((2 * qJD(5)) + t1003) * t1029 - t1082;
t954 = -qJD(6) * t1009 + t1037 * t1071 - t1067 * t995;
t955 = qJD(6) * t1008 + t1037 * t1067 + t1071 * t995;
t1027 = qJD(6) + t1030;
t970 = -mrSges(7,2) * t1027 + mrSges(7,3) * t1008;
t971 = mrSges(7,1) * t1027 - mrSges(7,3) * t1009;
t935 = -m(7) * t939 + t954 * mrSges(7,1) - t955 * mrSges(7,2) + t1008 * t970 - t1009 * t971;
t1117 = -2 * qJD(5);
t941 = t1029 * t1117 + t1125 * t1052 + t1082;
t999 = mrSges(6,1) * t1030 + mrSges(6,2) * t1029;
t1081 = -m(6) * t941 + t995 * mrSges(6,3) + t1029 * t999 - t935;
t945 = t1046 - t1088;
t1079 = m(5) * t945 - t1037 * mrSges(5,3) - t1052 * t1014 + t1081;
t1094 = -t1098 * t1029 - t1124 * t1030 + t1099 * t1052;
t1095 = t1123 * t1029 + t1098 * t1030 - t1097 * t1052;
t1016 = -mrSges(5,2) * t1029 - mrSges(5,3) * t1052;
t1010 = mrSges(6,1) * t1052 - mrSges(6,3) * t1029;
t1107 = t1029 * t1052;
t1103 = t1066 * t1072;
t997 = -g(3) * t1105 + t1040 * t1103 - t1069 * t1041;
t963 = -t1060 * pkin(2) - t1059 * pkin(9) + t1043 * t1092 - t997;
t996 = -t1029 * qJD(3) + t1115 * t1044 + t1068 * t1060;
t1084 = t995 * pkin(3) + t963 + (-t1107 - t996) * qJ(4);
t1078 = -t1028 * qJ(5) + t1125 * t1030 + qJDD(5) - t1084;
t1116 = -pkin(4) - pkin(10);
t936 = t1078 + (pkin(5) * t1029 + (pkin(3) + pkin(10)) * t1030) * t1052 + t1116 * t995 + t996 * pkin(5);
t948 = -t1068 * t964 + t1115 * t965;
t1086 = -qJ(4) * t1120 + t1030 * t1000 + qJDD(4) - t948;
t1077 = (-t996 + t1107) * qJ(5) + t1086 + (t1029 * pkin(4) + t1117) * t1030;
t937 = -t1120 * pkin(5) - t1030 * t1003 + (pkin(3) - t1116) * t1037 + t1077;
t933 = -t1067 * t937 + t1071 * t936;
t966 = -mrSges(7,1) * t1008 + mrSges(7,2) * t1009;
t988 = qJDD(6) + t996;
t930 = m(7) * t933 + mrSges(7,1) * t988 - t955 * mrSges(7,3) - t1009 * t966 + t1027 * t970;
t934 = t1067 * t936 + t1071 * t937;
t931 = m(7) * t934 - mrSges(7,2) * t988 + t954 * mrSges(7,3) + t1008 * t966 - t1027 * t971;
t920 = -t1067 * t930 + t1071 * t931;
t940 = (pkin(3) + pkin(4)) * t1037 + t1077;
t1087 = m(6) * t940 - t1037 * mrSges(6,2) - t1052 * t1010 + t920;
t946 = t1037 * pkin(3) + t1086;
t1080 = m(5) * t946 + t1037 * mrSges(5,1) + t1052 * t1016 + t1087;
t916 = (mrSges(5,2) - mrSges(6,3)) * t996 + (t1001 - t999) * t1030 + t1080;
t917 = -t996 * mrSges(6,3) - t1030 * t999 + t1087;
t956 = Ifges(7,5) * t1009 + Ifges(7,6) * t1008 + Ifges(7,3) * t1027;
t958 = Ifges(7,1) * t1009 + Ifges(7,4) * t1008 + Ifges(7,5) * t1027;
t924 = -mrSges(7,1) * t939 + mrSges(7,3) * t934 + Ifges(7,4) * t955 + Ifges(7,2) * t954 + Ifges(7,6) * t988 - t1009 * t956 + t1027 * t958;
t957 = Ifges(7,4) * t1009 + Ifges(7,2) * t1008 + Ifges(7,6) * t1027;
t925 = mrSges(7,2) * t939 - mrSges(7,3) * t933 + Ifges(7,1) * t955 + Ifges(7,4) * t954 + Ifges(7,5) * t988 + t1008 * t956 - t1027 * t957;
t1121 = -t1094 * t1029 - t1095 * t1030 - t1122 * t1037 + t1097 * t995 + t1099 * t996 + mrSges(4,1) * t948 - mrSges(5,1) * t946 - mrSges(6,1) * t941 - mrSges(4,2) * t949 + mrSges(6,2) * t940 + mrSges(5,3) * t945 - pkin(3) * t916 - pkin(4) * t917 - pkin(5) * t935 - pkin(10) * t920 + qJ(4) * (-t1037 * mrSges(6,1) - t995 * mrSges(5,2) - t1029 * t1001 - t1052 * t1015 + t1079) - t1067 * t925 - t1071 * t924;
t1114 = -mrSges(4,3) - mrSges(5,2);
t1113 = pkin(3) * t1052;
t1038 = mrSges(3,1) * t1061 - mrSges(3,3) * t1092;
t1042 = (-mrSges(3,1) * t1072 + mrSges(3,2) * t1069) * t1109;
t1011 = mrSges(4,2) * t1052 - mrSges(4,3) * t1029;
t1102 = -mrSges(4,1) * t1029 - mrSges(4,2) * t1030 - t1001;
t914 = m(4) * t948 - t1037 * mrSges(4,1) - t1052 * t1011 + (mrSges(6,3) + t1114) * t996 + (t999 + t1102) * t1030 - t1080;
t1013 = -mrSges(4,1) * t1052 - mrSges(4,3) * t1030;
t923 = t1079 + t1114 * t995 + (t1013 - t1015) * t1052 + (mrSges(4,2) - mrSges(6,1)) * t1037 + t1102 * t1029 + m(4) * t949;
t1090 = -t1068 * t914 + t1115 * t923;
t998 = -g(3) * t1106 + t1101;
t906 = m(3) * t998 - mrSges(3,2) * t1060 + mrSges(3,3) * t1045 - t1038 * t1061 + t1042 * t1091 + t1090;
t1021 = -t1065 * t1040 - t1111;
t1039 = -mrSges(3,2) * t1061 + mrSges(3,3) * t1091;
t909 = t1068 * t923 + t1115 * t914;
t908 = m(3) * t1021 - t1045 * mrSges(3,1) + t1044 * mrSges(3,2) + (t1038 * t1069 - t1039 * t1072) * t1109 + t909;
t919 = t1067 * t931 + t1071 * t930;
t943 = -t995 * pkin(4) + t1030 * t1113 + t1078;
t918 = m(6) * t943 + t996 * mrSges(6,1) + t995 * mrSges(6,2) + t1029 * t1010 + t1030 * t1015 + t919;
t947 = (t1119 - t1113) * t1030 + t1084;
t915 = m(5) * t947 + t995 * mrSges(5,1) - t996 * mrSges(5,3) - t1030 * t1014 + t1029 * t1016 - t918;
t1076 = -m(4) * t963 - t995 * mrSges(4,1) - t996 * mrSges(4,2) - t1029 * t1011 - t1030 * t1013 - t915;
t912 = m(3) * t997 + t1060 * mrSges(3,1) - t1044 * mrSges(3,3) + t1061 * t1039 - t1042 * t1092 + t1076;
t896 = -t1065 * t908 + t912 * t1103 + t906 * t1104;
t893 = m(2) * t1054 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1074 + t896;
t902 = -t1069 * t912 + t1072 * t906;
t900 = m(2) * t1055 - mrSges(2,1) * t1074 - qJDD(1) * mrSges(2,2) + t902;
t1110 = t1070 * t900 + t1073 * t893;
t895 = t1066 * t908 + t912 * t1105 + t906 * t1106;
t1093 = -t1097 * t1029 - t1099 * t1030 + t1122 * t1052;
t1089 = -t1070 * t893 + t1073 * t900;
t1019 = Ifges(3,6) * t1061 + (Ifges(3,4) * t1069 + Ifges(3,2) * t1072) * t1109;
t1020 = Ifges(3,5) * t1061 + (Ifges(3,1) * t1069 + Ifges(3,4) * t1072) * t1109;
t891 = -mrSges(4,1) * t963 + t1067 * t924 - t1071 * t925 + pkin(10) * t919 + pkin(4) * t918 + mrSges(6,3) * t941 - mrSges(6,2) * t943 - pkin(3) * t915 - qJ(5) * t1081 + mrSges(5,2) * t945 - mrSges(5,1) * t947 + mrSges(4,3) * t949 - t1098 * t996 - t1123 * t995 + (qJ(5) * t1015 + t1094) * t1052 + (mrSges(6,1) * qJ(5) + t1097) * t1037 + t1093 * t1030;
t1083 = mrSges(7,1) * t933 - mrSges(7,2) * t934 + Ifges(7,5) * t955 + Ifges(7,6) * t954 + Ifges(7,3) * t988 - t1008 * t958 + t1009 * t957;
t897 = mrSges(6,1) * t943 + mrSges(4,2) * t963 + mrSges(5,2) * t946 - mrSges(4,3) * t948 - mrSges(5,3) * t947 - mrSges(6,3) * t940 + pkin(5) * t919 - qJ(4) * t915 - qJ(5) * t917 + t1093 * t1029 - t1099 * t1037 - t1095 * t1052 + t1098 * t995 + t1124 * t996 + t1083;
t886 = Ifges(3,5) * t1044 + Ifges(3,6) * t1045 + Ifges(3,3) * t1060 + mrSges(3,1) * t997 - mrSges(3,2) * t998 + t1068 * t897 + t1115 * t891 + pkin(2) * t1076 + pkin(9) * t1090 + (t1019 * t1069 - t1020 * t1072) * t1109;
t1018 = Ifges(3,3) * t1061 + (Ifges(3,5) * t1069 + Ifges(3,6) * t1072) * t1109;
t888 = mrSges(3,2) * t1021 - mrSges(3,3) * t997 + Ifges(3,1) * t1044 + Ifges(3,4) * t1045 + Ifges(3,5) * t1060 - pkin(9) * t909 + t1018 * t1091 - t1061 * t1019 - t1068 * t891 + t1115 * t897;
t890 = -mrSges(3,1) * t1021 + mrSges(3,3) * t998 + Ifges(3,4) * t1044 + Ifges(3,2) * t1045 + Ifges(3,6) * t1060 - pkin(2) * t909 - t1018 * t1092 + t1061 * t1020 - t1121;
t1085 = mrSges(2,1) * t1054 - mrSges(2,2) * t1055 + Ifges(2,3) * qJDD(1) + pkin(1) * t896 + t1066 * t886 + t890 * t1105 + t888 * t1106 + t902 * t1112;
t884 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1054 + Ifges(2,5) * qJDD(1) - t1074 * Ifges(2,6) - t1069 * t890 + t1072 * t888 + (-t1065 * t895 - t1066 * t896) * pkin(8);
t883 = mrSges(2,1) * g(3) + mrSges(2,3) * t1055 + t1074 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t895 - t1065 * t886 + (pkin(8) * t902 + t1069 * t888 + t1072 * t890) * t1066;
t1 = [-m(1) * g(1) + t1089; -m(1) * g(2) + t1110; (-m(1) - m(2)) * g(3) + t895; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1110 - t1070 * t883 + t1073 * t884; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1089 + t1070 * t884 + t1073 * t883; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1085; t1085; t886; t1121; t916; t918; t1083;];
tauJB  = t1;
