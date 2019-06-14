% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR7
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
% Datum: 2019-05-06 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:03:35
% EndTime: 2019-05-06 11:03:46
% DurationCPUTime: 8.16s
% Computational Cost: add. (108157->369), mult. (249188->448), div. (0->0), fcn. (169145->10), ass. (0->160)
t1125 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t1101 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t1100 = Ifges(5,5) + Ifges(3,6) - Ifges(4,6);
t1099 = Ifges(5,6) + Ifges(3,5) + Ifges(4,4);
t1124 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t1123 = -Ifges(5,3) - Ifges(3,3) - Ifges(4,2);
t1067 = cos(qJ(2));
t1059 = sin(pkin(6));
t1063 = sin(qJ(2));
t1108 = t1059 * t1063;
t1094 = qJD(1) * t1108;
t1102 = qJDD(1) * t1059;
t1031 = -qJD(2) * t1094 + t1067 * t1102;
t1111 = qJD(1) * t1067;
t1030 = (qJD(2) * t1111 + qJDD(1) * t1063) * t1059;
t1060 = cos(pkin(6));
t1051 = qJD(1) * t1060 + qJD(2);
t1107 = t1059 * t1067;
t1093 = qJD(1) * t1107;
t1087 = t1051 * t1093;
t1064 = sin(qJ(1));
t1068 = cos(qJ(1));
t1045 = t1064 * g(1) - g(2) * t1068;
t1069 = qJD(1) ^ 2;
t1117 = pkin(8) * t1059;
t1023 = qJDD(1) * pkin(1) + t1069 * t1117 + t1045;
t990 = -t1060 * g(3) - t1059 * t1023;
t1086 = -t1031 * pkin(2) + t990 + (-t1030 - t1087) * qJ(3);
t1118 = pkin(2) * t1051;
t953 = (-(2 * qJD(3)) + t1118) * t1094 + t1086;
t1122 = m(4) * t953 - t1031 * mrSges(4,1);
t1121 = t1051 ^ 2;
t1120 = pkin(3) + pkin(9);
t1119 = mrSges(3,3) + mrSges(4,2);
t1116 = g(3) * t1063;
t1050 = qJDD(1) * t1060 + qJDD(2);
t1115 = t1050 * pkin(2);
t1105 = t1060 * t1067;
t1106 = t1060 * t1063;
t1019 = -mrSges(4,1) * t1051 + mrSges(4,2) * t1094;
t1021 = -mrSges(3,2) * t1051 + mrSges(3,3) * t1093;
t1022 = mrSges(4,2) * t1093 + mrSges(4,3) * t1051;
t1020 = -mrSges(5,1) * t1051 + mrSges(5,3) * t1093;
t1062 = sin(qJ(5));
t1066 = cos(qJ(5));
t1004 = -t1051 * t1066 + t1062 * t1093;
t1014 = qJDD(5) + t1030;
t1037 = qJD(5) + t1094;
t1061 = sin(qJ(6));
t1065 = cos(qJ(6));
t1002 = qJD(6) - t1004;
t1035 = t1037 ^ 2;
t1016 = -pkin(3) * t1051 - qJ(4) * t1094;
t1109 = t1059 ^ 2 * t1069;
t1092 = t1067 ^ 2 * t1109;
t1074 = -qJ(4) * t1092 + qJDD(4) - t1086 + ((2 * qJD(3)) + t1016) * t1094;
t1112 = qJD(1) * t1059;
t939 = t1030 * pkin(4) + t1120 * t1031 + (pkin(4) * t1067 + (-pkin(2) - pkin(9)) * t1063) * t1051 * t1112 + t1074;
t1029 = (pkin(4) * t1063 + pkin(9) * t1067) * t1112;
t1025 = (-pkin(2) * t1067 - qJ(3) * t1063) * t1112;
t1046 = -g(1) * t1068 - g(2) * t1064;
t1024 = -pkin(1) * t1069 + pkin(8) * t1102 + t1046;
t973 = -g(3) * t1107 + t1023 * t1105 - t1063 * t1024;
t1079 = -qJ(3) * t1121 + t1025 * t1094 + qJDD(3) - t973;
t1088 = -0.2e1 * qJD(4) * t1112;
t1072 = t1063 * t1088 + t1079 + (-t1030 + t1087) * qJ(4);
t1091 = t1067 * t1109;
t943 = -t1121 * pkin(4) + (-pkin(3) * t1091 - t1029 * t1112) * t1063 + (-pkin(2) - t1120) * t1050 + t1072;
t936 = t1062 * t939 + t1066 * t943;
t1005 = -t1051 * t1062 - t1066 * t1093;
t976 = -pkin(5) * t1004 - pkin(10) * t1005;
t933 = -pkin(5) * t1035 + pkin(10) * t1014 + t1004 * t976 + t936;
t1110 = qJD(3) * t1051;
t1036 = 0.2e1 * t1110;
t1104 = t1023 * t1106 + t1067 * t1024;
t1085 = pkin(2) * t1121 - t1050 * qJ(3) - t1025 * t1093 - t1104;
t1075 = pkin(3) * t1092 + t1031 * qJ(4) - t1051 * t1016 + t1085;
t941 = t1050 * pkin(4) - t1121 * pkin(9) + t1036 + t1067 * t1088 + (-t1029 * t1111 - t1116) * t1059 - t1075;
t971 = -qJD(5) * t1005 + t1031 * t1062 - t1050 * t1066;
t972 = qJD(5) * t1004 - t1031 * t1066 - t1050 * t1062;
t937 = t941 + (t1005 * t1037 - t971) * pkin(5) + (-t1004 * t1037 - t972) * pkin(10);
t930 = -t1061 * t933 + t1065 * t937;
t977 = -t1005 * t1061 + t1037 * t1065;
t949 = qJD(6) * t977 + t1014 * t1061 + t1065 * t972;
t978 = t1005 * t1065 + t1037 * t1061;
t959 = -mrSges(7,1) * t977 + mrSges(7,2) * t978;
t961 = -mrSges(7,2) * t1002 + mrSges(7,3) * t977;
t969 = qJDD(6) - t971;
t927 = m(7) * t930 + mrSges(7,1) * t969 - mrSges(7,3) * t949 + t1002 * t961 - t959 * t978;
t931 = t1061 * t937 + t1065 * t933;
t948 = -qJD(6) * t978 + t1014 * t1065 - t1061 * t972;
t962 = mrSges(7,1) * t1002 - mrSges(7,3) * t978;
t928 = m(7) * t931 - mrSges(7,2) * t969 + mrSges(7,3) * t948 - t1002 * t962 + t959 * t977;
t1090 = -t1061 * t927 + t1065 * t928;
t975 = -mrSges(6,1) * t1004 + mrSges(6,2) * t1005;
t980 = mrSges(6,1) * t1037 - mrSges(6,3) * t1005;
t915 = m(6) * t936 - mrSges(6,2) * t1014 + mrSges(6,3) * t971 + t1004 * t975 - t1037 * t980 + t1090;
t935 = -t1062 * t943 + t1066 * t939;
t932 = -pkin(5) * t1014 - pkin(10) * t1035 + t1005 * t976 - t935;
t1083 = -m(7) * t932 + t948 * mrSges(7,1) - mrSges(7,2) * t949 + t977 * t961 - t962 * t978;
t979 = -mrSges(6,2) * t1037 + mrSges(6,3) * t1004;
t923 = m(6) * t935 + mrSges(6,1) * t1014 - mrSges(6,3) * t972 - t1005 * t975 + t1037 * t979 + t1083;
t909 = t1062 * t915 + t1066 * t923;
t946 = t1031 * pkin(3) - t1094 * t1118 + t1074;
t1082 = -m(5) * t946 - t1030 * mrSges(5,1) + t1020 * t1093 - t909;
t1017 = mrSges(5,2) * t1051 - mrSges(5,3) * t1094;
t1103 = -mrSges(3,1) * t1051 + mrSges(3,3) * t1094 + t1017;
t901 = m(3) * t990 + (-mrSges(3,1) + mrSges(5,2)) * t1031 + (mrSges(3,2) - mrSges(4,3)) * t1030 + ((-t1021 - t1022) * t1067 + (-t1019 - t1103) * t1063) * t1112 + t1082 + t1122;
t1026 = (-mrSges(4,1) * t1067 - mrSges(4,3) * t1063) * t1112;
t1027 = (-mrSges(3,1) * t1067 + mrSges(3,2) * t1063) * t1112;
t1028 = (mrSges(5,1) * t1063 - mrSges(5,2) * t1067) * t1112;
t1113 = -t1062 * t923 + t1066 * t915;
t945 = -t1115 + (-t1063 * t1091 - t1050) * pkin(3) + t1072;
t1084 = m(5) * t945 + t1050 * mrSges(5,2) + t1051 * t1020 - t1028 * t1094 + t1113;
t958 = t1079 - t1115;
t1078 = m(4) * t958 - t1050 * mrSges(4,1) - t1051 * t1022 + t1084;
t904 = m(3) * t973 + t1050 * mrSges(3,1) + t1051 * t1021 + (-t1026 - t1027) * t1094 + (mrSges(5,3) - t1119) * t1030 - t1078;
t917 = t1061 * t928 + t1065 * t927;
t1080 = -m(6) * t941 + t971 * mrSges(6,1) - t972 * mrSges(6,2) + t1004 * t979 - t1005 * t980 - t917;
t944 = -0.2e1 * t1110 + (0.2e1 * qJD(4) * t1111 + t1116) * t1059 + t1075;
t1076 = -m(5) * t944 - t1031 * mrSges(5,3) - t1080;
t1098 = g(3) * t1108;
t952 = t1036 - t1085 - t1098;
t1070 = m(4) * t952 + t1050 * mrSges(4,3) + t1051 * t1019 + t1026 * t1093 + t1076;
t974 = -t1098 + t1104;
t913 = t1070 + (t1027 - t1028) * t1093 + m(3) * t974 + t1119 * t1031 + (-mrSges(3,2) + mrSges(5,1)) * t1050 + t1103 * t1051;
t892 = -t1059 * t901 + t904 * t1105 + t913 * t1106;
t889 = m(2) * t1045 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1069 + t892;
t897 = -t1063 * t904 + t1067 * t913;
t895 = m(2) * t1046 - mrSges(2,1) * t1069 - qJDD(1) * mrSges(2,2) + t897;
t1114 = t1064 * t895 + t1068 * t889;
t891 = t1060 * t901 + t904 * t1107 + t913 * t1108;
t1097 = (-t1063 * t1099 - t1067 * t1100) * t1112 + t1123 * t1051;
t1096 = (-t1063 * t1101 - t1067 * t1124) * t1112 - t1100 * t1051;
t1095 = (-t1063 * t1125 - t1101 * t1067) * t1112 - t1099 * t1051;
t1089 = -t1064 * t889 + t1068 * t895;
t954 = Ifges(7,5) * t978 + Ifges(7,6) * t977 + Ifges(7,3) * t1002;
t956 = Ifges(7,1) * t978 + Ifges(7,4) * t977 + Ifges(7,5) * t1002;
t920 = -mrSges(7,1) * t932 + mrSges(7,3) * t931 + Ifges(7,4) * t949 + Ifges(7,2) * t948 + Ifges(7,6) * t969 + t1002 * t956 - t954 * t978;
t955 = Ifges(7,4) * t978 + Ifges(7,2) * t977 + Ifges(7,6) * t1002;
t921 = mrSges(7,2) * t932 - mrSges(7,3) * t930 + Ifges(7,1) * t949 + Ifges(7,4) * t948 + Ifges(7,5) * t969 - t1002 * t955 + t954 * t977;
t963 = Ifges(6,5) * t1005 + Ifges(6,6) * t1004 + Ifges(6,3) * t1037;
t964 = Ifges(6,4) * t1005 + Ifges(6,2) * t1004 + Ifges(6,6) * t1037;
t898 = mrSges(6,2) * t941 - mrSges(6,3) * t935 + Ifges(6,1) * t972 + Ifges(6,4) * t971 + Ifges(6,5) * t1014 - pkin(10) * t917 + t1004 * t963 - t1037 * t964 - t1061 * t920 + t1065 * t921;
t1071 = mrSges(7,1) * t930 - mrSges(7,2) * t931 + Ifges(7,5) * t949 + Ifges(7,6) * t948 + Ifges(7,3) * t969 + t955 * t978 - t956 * t977;
t965 = Ifges(6,1) * t1005 + Ifges(6,4) * t1004 + Ifges(6,5) * t1037;
t899 = -mrSges(6,1) * t941 + mrSges(6,3) * t936 + Ifges(6,4) * t972 + Ifges(6,2) * t971 + Ifges(6,6) * t1014 - pkin(5) * t917 - t1005 * t963 + t1037 * t965 - t1071;
t906 = t1026 * t1094 + (mrSges(4,2) - mrSges(5,3)) * t1030 + t1078;
t908 = -t1030 * mrSges(5,3) + t1084;
t883 = -pkin(3) * t908 - mrSges(5,1) * t944 + mrSges(5,2) * t945 - pkin(9) * t1113 + mrSges(4,3) * t952 - mrSges(4,1) * t958 + mrSges(3,1) * t973 - mrSges(3,2) * t974 - pkin(4) * t1080 - t1062 * t898 - t1066 * t899 + qJ(3) * (t1051 * t1017 + t1070) - pkin(2) * t906 + (mrSges(5,1) * qJ(3) - t1123) * t1050 + (mrSges(4,2) * qJ(3) + t1100) * t1031 + t1099 * t1030 + (-t1096 * t1063 + (-qJ(3) * t1028 + t1095) * t1067) * t1112;
t1077 = t1031 * mrSges(5,2) + t1082;
t905 = -t1030 * mrSges(4,3) + (-t1022 * t1067 + (-t1017 - t1019) * t1063) * t1112 + t1077 + t1122;
t907 = t1017 * t1094 - t1077;
t885 = pkin(3) * t907 - pkin(2) * t905 + mrSges(5,3) * t944 - mrSges(5,2) * t946 + pkin(9) * t909 + mrSges(4,2) * t952 - mrSges(4,1) * t953 + mrSges(3,3) * t974 - mrSges(3,1) * t990 + t1062 * t899 - t1066 * t898 - qJ(4) * t1076 + (-qJ(4) * t1017 - t1095) * t1051 + (-mrSges(5,1) * qJ(4) + t1100) * t1050 + t1124 * t1031 + t1101 * t1030 + (qJ(4) * t1028 * t1067 + t1063 * t1097) * t1112;
t1073 = mrSges(6,1) * t935 - mrSges(6,2) * t936 + Ifges(6,5) * t972 + Ifges(6,6) * t971 + Ifges(6,3) * t1014 + pkin(5) * t1083 + pkin(10) * t1090 - t1004 * t965 + t1005 * t964 + t1061 * t921 + t1065 * t920;
t887 = mrSges(5,1) * t946 + mrSges(3,2) * t990 + mrSges(4,2) * t958 - mrSges(3,3) * t973 - mrSges(4,3) * t953 - mrSges(5,3) * t945 + pkin(4) * t909 - qJ(3) * t905 - qJ(4) * t908 + t1030 * t1125 + t1101 * t1031 + t1099 * t1050 + t1096 * t1051 - t1097 * t1093 + t1073;
t1081 = mrSges(2,1) * t1045 - mrSges(2,2) * t1046 + Ifges(2,3) * qJDD(1) + pkin(1) * t892 + t1060 * t883 + t885 * t1107 + t887 * t1108 + t897 * t1117;
t881 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1045 + Ifges(2,5) * qJDD(1) - t1069 * Ifges(2,6) - t1063 * t885 + t1067 * t887 + (-t1059 * t891 - t1060 * t892) * pkin(8);
t880 = mrSges(2,1) * g(3) + mrSges(2,3) * t1046 + t1069 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t891 - t1059 * t883 + (pkin(8) * t897 + t1063 * t887 + t1067 * t885) * t1060;
t1 = [-m(1) * g(1) + t1089; -m(1) * g(2) + t1114; (-m(1) - m(2)) * g(3) + t891; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1114 - t1064 * t880 + t1068 * t881; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1089 + t1064 * t881 + t1068 * t880; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1081; t1081; t883; t906; t907; t1073; t1071;];
tauJB  = t1;
