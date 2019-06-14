% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
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
% Datum: 2019-05-06 17:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR14_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:55:33
% EndTime: 2019-05-06 16:55:48
% DurationCPUTime: 9.73s
% Computational Cost: add. (126339->369), mult. (283980->439), div. (0->0), fcn. (200652->10), ass. (0->163)
t1107 = Ifges(5,4) + Ifges(6,6);
t1125 = -Ifges(5,2) - Ifges(6,3);
t1119 = Ifges(5,6) - Ifges(6,5);
t1124 = -2 * qJD(3);
t1123 = Ifges(3,1) + Ifges(4,2);
t1122 = Ifges(5,1) + Ifges(6,2);
t1108 = Ifges(3,4) + Ifges(4,6);
t1106 = Ifges(3,5) - Ifges(4,4);
t1121 = Ifges(5,5) - Ifges(6,4);
t1120 = Ifges(3,2) + Ifges(4,3);
t1105 = Ifges(3,6) - Ifges(4,5);
t1118 = Ifges(3,3) + Ifges(4,1);
t1117 = Ifges(5,3) + Ifges(6,1);
t1052 = cos(pkin(6));
t1045 = qJD(1) * t1052 + qJD(2);
t1054 = sin(qJ(4));
t1058 = cos(qJ(4));
t1051 = sin(pkin(6));
t1059 = cos(qJ(2));
t1088 = t1051 * t1059;
t1082 = qJD(1) * t1088;
t1008 = t1045 * t1054 + t1058 * t1082;
t1009 = t1045 * t1058 - t1054 * t1082;
t1055 = sin(qJ(2));
t1089 = t1051 * t1055;
t1083 = qJD(1) * t1089;
t1036 = qJD(4) + t1083;
t1116 = t1125 * t1008 + t1107 * t1009 + t1119 * t1036;
t1115 = (pkin(2) * t1045 + t1124) * t1083;
t1093 = qJD(1) * t1051;
t1025 = (-pkin(2) * t1059 - qJ(3) * t1055) * t1093;
t1043 = t1045 ^ 2;
t1044 = qJDD(1) * t1052 + qJDD(2);
t1056 = sin(qJ(1));
t1060 = cos(qJ(1));
t1040 = t1056 * g(1) - g(2) * t1060;
t1061 = qJD(1) ^ 2;
t1104 = pkin(8) * t1051;
t1023 = qJDD(1) * pkin(1) + t1061 * t1104 + t1040;
t1041 = -g(1) * t1060 - g(2) * t1056;
t1084 = qJDD(1) * t1051;
t1024 = -pkin(1) * t1061 + pkin(8) * t1084 + t1041;
t1087 = t1052 * t1055;
t978 = -g(3) * t1089 + t1023 * t1087 + t1059 * t1024;
t944 = t1043 * pkin(2) - t1044 * qJ(3) - t1025 * t1082 + t1045 * t1124 - t978;
t1030 = -qJD(2) * t1083 + t1059 * t1084;
t975 = qJD(4) * t1009 + t1058 * t1030 + t1044 * t1054;
t985 = mrSges(6,1) * t1008 - mrSges(6,3) * t1036;
t1114 = -t975 * mrSges(6,2) - t1008 * t985;
t987 = -mrSges(5,2) * t1036 - mrSges(5,3) * t1008;
t1098 = t985 - t987;
t1109 = mrSges(5,1) - mrSges(6,2);
t1113 = -t1098 * t1008 + t1109 * t975;
t1112 = -2 * qJD(5);
t1111 = -pkin(2) - pkin(9);
t1110 = mrSges(3,1) - mrSges(4,2);
t1103 = t1052 * g(3);
t1086 = t1052 * t1059;
t1019 = mrSges(3,1) * t1045 - mrSges(3,3) * t1083;
t1020 = -mrSges(3,2) * t1045 + mrSges(3,3) * t1082;
t1022 = mrSges(4,1) * t1083 + mrSges(4,2) * t1045;
t1029 = (qJD(1) * qJD(2) * t1059 + qJDD(1) * t1055) * t1051;
t1021 = -mrSges(4,1) * t1082 - mrSges(4,3) * t1045;
t1018 = qJDD(4) + t1029;
t1053 = sin(qJ(6));
t1057 = cos(qJ(6));
t1006 = qJD(6) + t1009;
t1091 = t1008 * t1036;
t1033 = t1036 ^ 2;
t1028 = pkin(3) * t1083 - pkin(9) * t1045;
t1090 = t1051 ^ 2 * t1061;
t1081 = t1059 ^ 2 * t1090;
t934 = -pkin(3) * t1081 - t1103 - t1029 * qJ(3) + t1111 * t1030 + (-t1023 + (-qJ(3) * t1045 * t1059 - t1028 * t1055) * qJD(1)) * t1051 + t1115;
t1085 = g(3) * t1088 + t1055 * t1024;
t1074 = -t1043 * qJ(3) + t1025 * t1083 + qJDD(3) + t1085;
t936 = t1029 * pkin(3) + t1111 * t1044 + (-pkin(3) * t1045 * t1093 - pkin(9) * t1055 * t1090 - t1023 * t1052) * t1059 + t1074;
t928 = -t1054 * t934 + t1058 * t936;
t979 = pkin(4) * t1008 - qJ(5) * t1009;
t924 = -t1018 * pkin(4) - t1033 * qJ(5) + t1009 * t979 + qJDD(5) - t928;
t976 = -qJD(4) * t1008 - t1030 * t1054 + t1044 * t1058;
t918 = (t1008 * t1009 - t1018) * pkin(10) + (t976 + t1091) * pkin(5) + t924;
t1007 = t1008 ^ 2;
t933 = t1030 * pkin(3) - pkin(9) * t1081 + t1045 * t1028 - t944;
t1062 = (-t976 + t1091) * qJ(5) + t933 + (t1036 * pkin(4) + t1112) * t1009;
t989 = pkin(5) * t1009 - pkin(10) * t1036;
t921 = t1062 + (pkin(4) + pkin(10)) * t975 - t1009 * t989 - t1007 * pkin(5);
t916 = -t1053 * t921 + t1057 * t918;
t983 = t1008 * t1057 - t1036 * t1053;
t942 = qJD(6) * t983 + t1018 * t1057 + t1053 * t975;
t984 = t1008 * t1053 + t1036 * t1057;
t951 = -mrSges(7,1) * t983 + mrSges(7,2) * t984;
t956 = -mrSges(7,2) * t1006 + mrSges(7,3) * t983;
t972 = qJDD(6) + t976;
t913 = m(7) * t916 + mrSges(7,1) * t972 - t942 * mrSges(7,3) + t1006 * t956 - t951 * t984;
t917 = t1053 * t918 + t1057 * t921;
t941 = -qJD(6) * t984 - t1018 * t1053 + t1057 * t975;
t957 = mrSges(7,1) * t1006 - mrSges(7,3) * t984;
t914 = m(7) * t917 - mrSges(7,2) * t972 + t941 * mrSges(7,3) - t1006 * t957 + t951 * t983;
t905 = t1053 * t914 + t1057 * t913;
t981 = -mrSges(6,2) * t1008 - mrSges(6,3) * t1009;
t1071 = -m(6) * t924 - t976 * mrSges(6,1) - t1009 * t981 - t905;
t980 = mrSges(5,1) * t1008 + mrSges(5,2) * t1009;
t901 = m(5) * t928 - t976 * mrSges(5,3) - t1009 * t980 + t1109 * t1018 - t1098 * t1036 + t1071;
t929 = t1054 * t936 + t1058 * t934;
t1069 = -t1033 * pkin(4) + t1018 * qJ(5) - t1008 * t979 + t929;
t920 = -t975 * pkin(5) - t1007 * pkin(10) + ((2 * qJD(5)) + t989) * t1036 + t1069;
t1073 = -m(7) * t920 + t941 * mrSges(7,1) - t942 * mrSges(7,2) + t983 * t956 - t984 * t957;
t922 = t1036 * t1112 - t1069;
t986 = mrSges(6,1) * t1009 + mrSges(6,2) * t1036;
t1065 = -m(6) * t922 + t1018 * mrSges(6,3) + t1036 * t986 - t1073;
t988 = mrSges(5,1) * t1036 - mrSges(5,3) * t1009;
t910 = m(5) * t929 - t1018 * mrSges(5,2) - t1036 * t988 + (-mrSges(5,3) - mrSges(6,1)) * t975 + (-t980 - t981) * t1008 + t1065;
t1078 = -t1054 * t901 + t1058 * t910;
t997 = -t1051 * t1023 - t1103;
t945 = -t1030 * pkin(2) + (-t1045 * t1082 - t1029) * qJ(3) + t997 + t1115;
t1075 = m(4) * t945 - t1029 * mrSges(4,3) + t1021 * t1082 + t1078;
t892 = m(3) * t997 + t1029 * mrSges(3,2) - t1110 * t1030 + (-t1020 * t1059 + (t1019 - t1022) * t1055) * t1093 + t1075;
t1026 = (mrSges(4,2) * t1059 - mrSges(4,3) * t1055) * t1093;
t1027 = (-mrSges(3,1) * t1059 + mrSges(3,2) * t1055) * t1093;
t896 = t1054 * t910 + t1058 * t901;
t1080 = t1023 * t1086;
t950 = -t1044 * pkin(2) + t1074 - t1080;
t1072 = -m(4) * t950 - t1029 * mrSges(4,1) - t896;
t977 = t1080 - t1085;
t893 = m(3) * t977 - t1029 * mrSges(3,3) + (t1020 - t1021) * t1045 + t1110 * t1044 + (-t1026 - t1027) * t1083 + t1072;
t1079 = -t1053 * t913 + t1057 * t914;
t926 = t975 * pkin(4) + t1062;
t1076 = -m(6) * t926 + t976 * mrSges(6,3) + t1009 * t986 - t1079;
t1070 = -m(5) * t933 - t976 * mrSges(5,2) - t1009 * t988 + t1076;
t1064 = -m(4) * t944 + t1044 * mrSges(4,3) + t1045 * t1022 + t1026 * t1082 - t1070;
t900 = t1027 * t1082 + t1064 - t1044 * mrSges(3,2) - t1045 * t1019 + m(3) * t978 + (mrSges(3,3) + mrSges(4,1)) * t1030 + t1113;
t881 = -t1051 * t892 + t893 * t1086 + t900 * t1087;
t878 = m(2) * t1040 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1061 + t881;
t888 = -t1055 * t893 + t1059 * t900;
t886 = m(2) * t1041 - mrSges(2,1) * t1061 - qJDD(1) * mrSges(2,2) + t888;
t1101 = t1056 * t886 + t1060 * t878;
t1100 = t1119 * t1008 - t1121 * t1009 - t1117 * t1036;
t1099 = -t1107 * t1008 + t1122 * t1009 + t1121 * t1036;
t1097 = (t1106 * t1055 + t1105 * t1059) * t1093 + t1118 * t1045;
t1096 = (t1108 * t1055 + t1120 * t1059) * t1093 + t1105 * t1045;
t1095 = (t1123 * t1055 + t1108 * t1059) * t1093 + t1106 * t1045;
t880 = t1052 * t892 + t893 * t1088 + t900 * t1089;
t1077 = -t1056 * t878 + t1060 * t886;
t904 = -t1076 + t1114;
t946 = Ifges(7,5) * t984 + Ifges(7,6) * t983 + Ifges(7,3) * t1006;
t948 = Ifges(7,1) * t984 + Ifges(7,4) * t983 + Ifges(7,5) * t1006;
t907 = -mrSges(7,1) * t920 + mrSges(7,3) * t917 + Ifges(7,4) * t942 + Ifges(7,2) * t941 + Ifges(7,6) * t972 + t1006 * t948 - t946 * t984;
t947 = Ifges(7,4) * t984 + Ifges(7,2) * t983 + Ifges(7,6) * t1006;
t908 = mrSges(7,2) * t920 - mrSges(7,3) * t916 + Ifges(7,1) * t942 + Ifges(7,4) * t941 + Ifges(7,5) * t972 - t1006 * t947 + t946 * t983;
t882 = -mrSges(5,1) * t933 - mrSges(6,1) * t922 + mrSges(6,2) * t926 + mrSges(5,3) * t929 - pkin(4) * t904 - pkin(5) * t1073 - pkin(10) * t1079 + t1100 * t1009 + t1119 * t1018 + t1099 * t1036 - t1053 * t908 - t1057 * t907 + t1107 * t976 + t1125 * t975;
t1067 = mrSges(7,1) * t916 - mrSges(7,2) * t917 + Ifges(7,5) * t942 + Ifges(7,6) * t941 + Ifges(7,3) * t972 + t984 * t947 - t983 * t948;
t883 = mrSges(6,1) * t924 + mrSges(5,2) * t933 - mrSges(5,3) * t928 - mrSges(6,3) * t926 + pkin(5) * t905 - qJ(5) * t904 + t1100 * t1008 + t1121 * t1018 - t1116 * t1036 - t1107 * t975 + t1122 * t976 + t1067;
t895 = t1044 * mrSges(4,2) + t1045 * t1021 + t1026 * t1083 - t1072;
t872 = mrSges(3,1) * t977 - mrSges(3,2) * t978 + mrSges(4,2) * t950 - mrSges(4,3) * t944 + t1058 * t883 - t1054 * t882 - pkin(9) * t896 - pkin(2) * t895 + qJ(3) * (t975 * mrSges(5,1) + t1008 * t987 + t1064 + t1114) + t1118 * t1044 + (mrSges(4,1) * qJ(3) + t1105) * t1030 + t1106 * t1029 + (t1096 * t1055 - t1095 * t1059) * t1093;
t894 = t1030 * mrSges(4,2) - t1022 * t1083 + t1075;
t874 = -mrSges(3,1) * t997 + mrSges(3,3) * t978 - mrSges(4,1) * t944 + mrSges(4,2) * t945 - t1054 * t883 - t1058 * t882 - pkin(3) * (t1070 - t1113) - pkin(9) * t1078 - pkin(2) * t894 + t1095 * t1045 + t1105 * t1044 + t1120 * t1030 + t1108 * t1029 - t1097 * t1083;
t903 = t1018 * mrSges(6,2) + t1036 * t985 - t1071;
t1063 = -mrSges(5,2) * t929 - mrSges(6,3) * t922 - pkin(4) * t903 - pkin(10) * t905 + t1099 * t1008 - t1053 * t907 + t1057 * t908 + qJ(5) * (-t1008 * t981 + t1065) + mrSges(6,2) * t924 + mrSges(5,1) * t928 + t1121 * t976 + (-qJ(5) * mrSges(6,1) - t1119) * t975 + t1117 * t1018 + t1116 * t1009;
t876 = mrSges(4,1) * t950 + mrSges(3,2) * t997 - mrSges(3,3) * t977 - mrSges(4,3) * t945 + pkin(3) * t896 - qJ(3) * t894 + t1123 * t1029 + t1108 * t1030 + t1106 * t1044 - t1096 * t1045 + t1097 * t1082 + t1063;
t1068 = mrSges(2,1) * t1040 - mrSges(2,2) * t1041 + Ifges(2,3) * qJDD(1) + pkin(1) * t881 + t1052 * t872 + t874 * t1088 + t876 * t1089 + t888 * t1104;
t870 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1040 + Ifges(2,5) * qJDD(1) - t1061 * Ifges(2,6) - t1055 * t874 + t1059 * t876 + (-t1051 * t880 - t1052 * t881) * pkin(8);
t869 = mrSges(2,1) * g(3) + mrSges(2,3) * t1041 + t1061 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t880 - t1051 * t872 + (pkin(8) * t888 + t1055 * t876 + t1059 * t874) * t1052;
t1 = [-m(1) * g(1) + t1077; -m(1) * g(2) + t1101; (-m(1) - m(2)) * g(3) + t880; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1101 - t1056 * t869 + t1060 * t870; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1077 + t1056 * t870 + t1060 * t869; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1068; t1068; t872; t895; t1063; t903; t1067;];
tauJB  = t1;
