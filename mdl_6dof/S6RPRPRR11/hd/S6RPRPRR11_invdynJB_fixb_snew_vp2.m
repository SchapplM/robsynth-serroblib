% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 20:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:16:35
% EndTime: 2019-05-05 20:18:42
% DurationCPUTime: 131.16s
% Computational Cost: add. (1965721->404), mult. (6238273->551), div. (0->0), fcn. (5358305->16), ass. (0->185)
t1072 = cos(pkin(6));
t1075 = sin(qJ(3));
t1070 = cos(pkin(12));
t1068 = sin(pkin(6));
t1071 = cos(pkin(7));
t1111 = t1068 * t1071;
t1104 = t1070 * t1111;
t1067 = sin(pkin(7));
t1122 = cos(qJ(3));
t1108 = t1067 * t1122;
t1066 = sin(pkin(12));
t1115 = t1066 * t1068;
t1123 = -t1072 * t1108 + t1075 * t1115 - t1122 * t1104;
t1121 = pkin(9) * t1066;
t1120 = Ifges(3,3) * t1072;
t1119 = pkin(9) * qJDD(1);
t1076 = sin(qJ(1));
t1079 = cos(qJ(1));
t1061 = t1076 * g(1) - g(2) * t1079;
t1080 = qJD(1) ^ 2;
t1110 = t1070 * t1072;
t1114 = t1066 * t1072;
t1062 = -g(1) * t1079 - g(2) * t1076;
t1117 = qJ(2) * t1068;
t1050 = -pkin(1) * t1080 + qJDD(1) * t1117 + t1062;
t1049 = qJDD(1) * pkin(1) + t1080 * t1117 + t1061;
t1116 = qJD(1) * t1068;
t1106 = qJD(2) * t1116;
t1112 = t1068 * t1070;
t1094 = -g(3) * t1112 + t1049 * t1110 - 0.2e1 * t1066 * t1106;
t1026 = -t1066 * t1050 + t1094;
t1097 = -mrSges(3,1) * t1070 + mrSges(3,2) * t1066;
t1048 = t1097 * t1116;
t1091 = -mrSges(3,2) * t1072 + mrSges(3,3) * t1112;
t1053 = t1091 * qJD(1);
t1092 = mrSges(3,1) * t1072 - mrSges(3,3) * t1115;
t1107 = t1071 * t1122;
t1109 = t1071 * t1075;
t1034 = t1123 * qJD(1);
t1113 = t1067 * t1075;
t1083 = t1072 * t1113 + (t1066 * t1122 + t1070 * t1109) * t1068;
t1035 = t1083 * qJD(1);
t1020 = mrSges(4,1) * t1034 + mrSges(4,2) * t1035;
t1021 = qJD(3) * t1035 + t1123 * qJDD(1);
t1089 = -t1067 * t1112 + t1071 * t1072;
t1047 = qJD(1) * t1089 + qJD(3);
t1031 = mrSges(4,1) * t1047 - mrSges(4,3) * t1035;
t1044 = qJDD(1) * t1089 + qJDD(3);
t1065 = sin(pkin(13));
t1069 = cos(pkin(13));
t1028 = -t1035 * t1065 + t1047 * t1069;
t1029 = t1035 * t1069 + t1047 * t1065;
t1004 = -mrSges(5,1) * t1028 + mrSges(5,2) * t1029;
t1010 = -mrSges(5,2) * t1034 + mrSges(5,3) * t1028;
t1022 = -t1034 * qJD(3) + qJDD(1) * t1083;
t1014 = t1022 * t1069 + t1044 * t1065;
t1074 = sin(qJ(5));
t1078 = cos(qJ(5));
t1002 = t1028 * t1078 - t1029 * t1074;
t1018 = qJDD(5) + t1021;
t1033 = qJD(5) + t1034;
t1073 = sin(qJ(6));
t1077 = cos(qJ(6));
t1001 = qJD(6) - t1002;
t1032 = t1033 ^ 2;
t1019 = pkin(3) * t1034 - qJ(4) * t1035;
t1043 = t1047 ^ 2;
t1051 = (pkin(2) * t1072 - t1111 * t1121) * qJD(1);
t1093 = -pkin(2) * t1070 - t1067 * t1121;
t1090 = qJD(1) * t1093 * t1116 + t1071 * t1119;
t1103 = t1049 * t1114 + (t1050 + 0.2e1 * t1106) * t1070;
t1000 = (-qJD(1) * t1051 + t1067 * t1119) * t1072 + (-g(3) * t1066 + t1070 * t1090) * t1068 + t1103;
t1046 = (t1067 * t1072 + t1104) * qJD(1) * pkin(9);
t1102 = -t1072 * g(3) + qJDD(2);
t1008 = (-t1049 + t1093 * qJDD(1) + (-t1046 * t1070 + t1051 * t1066) * qJD(1)) * t1068 + t1102;
t999 = (pkin(2) * qJDD(1) + qJD(1) * t1046) * t1072 + (-t1068 * t1090 - t1050) * t1066 + t1094;
t972 = t1122 * t1000 + t1008 * t1113 + t999 * t1109;
t963 = -pkin(3) * t1043 + qJ(4) * t1044 - t1019 * t1034 + t972;
t986 = t1071 * t1008 - t1067 * t999;
t966 = (t1034 * t1047 - t1022) * qJ(4) + (t1035 * t1047 + t1021) * pkin(3) + t986;
t955 = -0.2e1 * qJD(4) * t1029 - t1065 * t963 + t1069 * t966;
t952 = (t1028 * t1034 - t1014) * pkin(10) + (t1028 * t1029 + t1021) * pkin(4) + t955;
t1012 = pkin(4) * t1034 - pkin(10) * t1029;
t1013 = -t1022 * t1065 + t1044 * t1069;
t1025 = t1028 ^ 2;
t956 = 0.2e1 * qJD(4) * t1028 + t1065 * t966 + t1069 * t963;
t954 = -pkin(4) * t1025 + pkin(10) * t1013 - t1012 * t1034 + t956;
t949 = t1074 * t952 + t1078 * t954;
t1003 = t1028 * t1074 + t1029 * t1078;
t985 = -pkin(5) * t1002 - pkin(11) * t1003;
t947 = -pkin(5) * t1032 + pkin(11) * t1018 + t1002 * t985 + t949;
t971 = -t1075 * t1000 + t1008 * t1108 + t1107 * t999;
t962 = -t1044 * pkin(3) - t1043 * qJ(4) + t1035 * t1019 + qJDD(4) - t971;
t957 = -t1013 * pkin(4) - t1025 * pkin(10) + t1029 * t1012 + t962;
t976 = -qJD(5) * t1003 + t1013 * t1078 - t1014 * t1074;
t977 = qJD(5) * t1002 + t1013 * t1074 + t1014 * t1078;
t950 = (-t1002 * t1033 - t977) * pkin(11) + (t1003 * t1033 - t976) * pkin(5) + t957;
t944 = -t1073 * t947 + t1077 * t950;
t987 = -t1003 * t1073 + t1033 * t1077;
t960 = qJD(6) * t987 + t1018 * t1073 + t1077 * t977;
t988 = t1003 * t1077 + t1033 * t1073;
t973 = -mrSges(7,1) * t987 + mrSges(7,2) * t988;
t975 = qJDD(6) - t976;
t978 = -mrSges(7,2) * t1001 + mrSges(7,3) * t987;
t941 = m(7) * t944 + mrSges(7,1) * t975 - mrSges(7,3) * t960 + t1001 * t978 - t973 * t988;
t945 = t1073 * t950 + t1077 * t947;
t959 = -qJD(6) * t988 + t1018 * t1077 - t1073 * t977;
t979 = mrSges(7,1) * t1001 - mrSges(7,3) * t988;
t942 = m(7) * t945 - mrSges(7,2) * t975 + mrSges(7,3) * t959 - t1001 * t979 + t973 * t987;
t933 = -t1073 * t941 + t1077 * t942;
t984 = -mrSges(6,1) * t1002 + mrSges(6,2) * t1003;
t990 = mrSges(6,1) * t1033 - mrSges(6,3) * t1003;
t927 = m(6) * t949 - mrSges(6,2) * t1018 + mrSges(6,3) * t976 + t1002 * t984 - t1033 * t990 + t933;
t948 = -t1074 * t954 + t1078 * t952;
t946 = -pkin(5) * t1018 - pkin(11) * t1032 + t1003 * t985 - t948;
t943 = -m(7) * t946 + t959 * mrSges(7,1) - mrSges(7,2) * t960 + t987 * t978 - t979 * t988;
t989 = -mrSges(6,2) * t1033 + mrSges(6,3) * t1002;
t937 = m(6) * t948 + mrSges(6,1) * t1018 - mrSges(6,3) * t977 - t1003 * t984 + t1033 * t989 + t943;
t924 = t1074 * t927 + t1078 * t937;
t922 = m(5) * t955 + mrSges(5,1) * t1021 - mrSges(5,3) * t1014 - t1004 * t1029 + t1010 * t1034 + t924;
t1011 = mrSges(5,1) * t1034 - mrSges(5,3) * t1029;
t1100 = -t1074 * t937 + t1078 * t927;
t923 = m(5) * t956 - mrSges(5,2) * t1021 + mrSges(5,3) * t1013 + t1004 * t1028 - t1011 * t1034 + t1100;
t1101 = -t1065 * t922 + t1069 * t923;
t913 = m(4) * t972 - mrSges(4,2) * t1044 - mrSges(4,3) * t1021 - t1020 * t1034 - t1031 * t1047 + t1101;
t1030 = -mrSges(4,2) * t1047 - mrSges(4,3) * t1034;
t916 = t1065 * t923 + t1069 * t922;
t915 = m(4) * t986 + mrSges(4,1) * t1021 + mrSges(4,2) * t1022 + t1030 * t1034 + t1031 * t1035 + t916;
t932 = t1073 * t942 + t1077 * t941;
t1084 = m(6) * t957 - t976 * mrSges(6,1) + mrSges(6,2) * t977 - t1002 * t989 + t1003 * t990 + t932;
t931 = m(5) * t962 - t1013 * mrSges(5,1) + mrSges(5,2) * t1014 - t1028 * t1010 + t1011 * t1029 + t1084;
t930 = m(4) * t971 + mrSges(4,1) * t1044 - mrSges(4,3) * t1022 - t1020 * t1035 + t1030 * t1047 - t931;
t902 = -t1067 * t915 + t930 * t1107 + t913 * t1109;
t898 = m(3) * t1026 + t1092 * qJDD(1) + (-t1048 * t1115 + t1053 * t1072) * qJD(1) + t902;
t1036 = -t1068 * t1049 + t1102;
t1052 = t1092 * qJD(1);
t901 = t1071 * t915 + t930 * t1108 + t913 * t1113;
t900 = m(3) * t1036 + (t1097 * qJDD(1) + (t1052 * t1066 - t1053 * t1070) * qJD(1)) * t1068 + t901;
t1027 = -g(3) * t1115 + t1103;
t909 = -t1075 * t930 + t1122 * t913;
t908 = m(3) * t1027 + t1091 * qJDD(1) + (t1048 * t1112 - t1052 * t1072) * qJD(1) + t909;
t887 = -t1068 * t900 + t898 * t1110 + t908 * t1114;
t884 = m(2) * t1061 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1080 + t887;
t894 = -t1066 * t898 + t1070 * t908;
t892 = m(2) * t1062 - mrSges(2,1) * t1080 - qJDD(1) * mrSges(2,2) + t894;
t1118 = t1076 * t892 + t1079 * t884;
t886 = t1072 * t900 + t898 * t1112 + t908 * t1115;
t1099 = -t1076 * t884 + t1079 * t892;
t1096 = Ifges(3,5) * t1066 + Ifges(3,6) * t1070;
t1015 = Ifges(4,5) * t1035 - Ifges(4,6) * t1034 + Ifges(4,3) * t1047;
t1016 = Ifges(4,4) * t1035 - Ifges(4,2) * t1034 + Ifges(4,6) * t1047;
t967 = Ifges(7,5) * t988 + Ifges(7,6) * t987 + Ifges(7,3) * t1001;
t969 = Ifges(7,1) * t988 + Ifges(7,4) * t987 + Ifges(7,5) * t1001;
t934 = -mrSges(7,1) * t946 + mrSges(7,3) * t945 + Ifges(7,4) * t960 + Ifges(7,2) * t959 + Ifges(7,6) * t975 + t1001 * t969 - t967 * t988;
t968 = Ifges(7,4) * t988 + Ifges(7,2) * t987 + Ifges(7,6) * t1001;
t935 = mrSges(7,2) * t946 - mrSges(7,3) * t944 + Ifges(7,1) * t960 + Ifges(7,4) * t959 + Ifges(7,5) * t975 - t1001 * t968 + t967 * t987;
t980 = Ifges(6,5) * t1003 + Ifges(6,6) * t1002 + Ifges(6,3) * t1033;
t981 = Ifges(6,4) * t1003 + Ifges(6,2) * t1002 + Ifges(6,6) * t1033;
t917 = mrSges(6,2) * t957 - mrSges(6,3) * t948 + Ifges(6,1) * t977 + Ifges(6,4) * t976 + Ifges(6,5) * t1018 - pkin(11) * t932 + t1002 * t980 - t1033 * t981 - t1073 * t934 + t1077 * t935;
t1082 = mrSges(7,1) * t944 - mrSges(7,2) * t945 + Ifges(7,5) * t960 + Ifges(7,6) * t959 + Ifges(7,3) * t975 + t968 * t988 - t969 * t987;
t982 = Ifges(6,1) * t1003 + Ifges(6,4) * t1002 + Ifges(6,5) * t1033;
t918 = -mrSges(6,1) * t957 + mrSges(6,3) * t949 + Ifges(6,4) * t977 + Ifges(6,2) * t976 + Ifges(6,6) * t1018 - pkin(5) * t932 - t1003 * t980 + t1033 * t982 - t1082;
t991 = Ifges(5,5) * t1029 + Ifges(5,6) * t1028 + Ifges(5,3) * t1034;
t993 = Ifges(5,1) * t1029 + Ifges(5,4) * t1028 + Ifges(5,5) * t1034;
t903 = -mrSges(5,1) * t962 + mrSges(5,3) * t956 + Ifges(5,4) * t1014 + Ifges(5,2) * t1013 + Ifges(5,6) * t1021 - pkin(4) * t1084 + pkin(10) * t1100 - t1029 * t991 + t1034 * t993 + t1074 * t917 + t1078 * t918;
t992 = Ifges(5,4) * t1029 + Ifges(5,2) * t1028 + Ifges(5,6) * t1034;
t904 = mrSges(5,2) * t962 - mrSges(5,3) * t955 + Ifges(5,1) * t1014 + Ifges(5,4) * t1013 + Ifges(5,5) * t1021 - pkin(10) * t924 + t1028 * t991 - t1034 * t992 - t1074 * t918 + t1078 * t917;
t889 = mrSges(4,2) * t986 - mrSges(4,3) * t971 + Ifges(4,1) * t1022 - Ifges(4,4) * t1021 + Ifges(4,5) * t1044 - qJ(4) * t916 - t1015 * t1034 - t1016 * t1047 - t1065 * t903 + t1069 * t904;
t1017 = Ifges(4,1) * t1035 - Ifges(4,4) * t1034 + Ifges(4,5) * t1047;
t1081 = mrSges(6,1) * t948 - mrSges(6,2) * t949 + Ifges(6,5) * t977 + Ifges(6,6) * t976 + Ifges(6,3) * t1018 + pkin(5) * t943 + pkin(11) * t933 - t1002 * t982 + t1003 * t981 + t1073 * t935 + t1077 * t934;
t895 = -pkin(4) * t924 - pkin(3) * t916 - t1035 * t1015 + Ifges(4,6) * t1044 + t1047 * t1017 + t1028 * t993 - t1029 * t992 + Ifges(4,4) * t1022 - Ifges(5,6) * t1013 - Ifges(5,5) * t1014 - mrSges(4,1) * t986 + mrSges(4,3) * t972 + (-Ifges(5,3) - Ifges(4,2)) * t1021 + mrSges(5,2) * t956 - mrSges(5,1) * t955 - t1081;
t1088 = pkin(9) * t909 + t1075 * t889 + t1122 * t895;
t1085 = Ifges(3,6) * t1072 + (Ifges(3,4) * t1066 + Ifges(3,2) * t1070) * t1068;
t1040 = t1085 * qJD(1);
t1086 = Ifges(3,5) * t1072 + (Ifges(3,1) * t1066 + Ifges(3,4) * t1070) * t1068;
t1041 = t1086 * qJD(1);
t888 = mrSges(4,1) * t971 - mrSges(4,2) * t972 + Ifges(4,5) * t1022 - Ifges(4,6) * t1021 + Ifges(4,3) * t1044 - pkin(3) * t931 + qJ(4) * t1101 + t1035 * t1016 + t1034 * t1017 + t1065 * t904 + t1069 * t903;
t878 = qJDD(1) * t1120 + mrSges(3,1) * t1026 - mrSges(3,2) * t1027 + pkin(2) * t902 + t1071 * t888 + t1088 * t1067 + (t1096 * qJDD(1) + (t1040 * t1066 - t1041 * t1070) * qJD(1)) * t1068;
t1039 = (t1068 * t1096 + t1120) * qJD(1);
t880 = -mrSges(3,1) * t1036 + mrSges(3,3) * t1027 - pkin(2) * t901 - t1067 * t888 + (-t1039 * t1115 + t1041 * t1072) * qJD(1) + t1088 * t1071 + t1085 * qJDD(1);
t882 = mrSges(3,2) * t1036 - mrSges(3,3) * t1026 + t1122 * t889 - t1075 * t895 + (t1039 * t1112 - t1040 * t1072) * qJD(1) + (-t1067 * t901 - t1071 * t902) * pkin(9) + t1086 * qJDD(1);
t1087 = mrSges(2,1) * t1061 - mrSges(2,2) * t1062 + Ifges(2,3) * qJDD(1) + pkin(1) * t887 + t1072 * t878 + t880 * t1112 + t882 * t1115 + t894 * t1117;
t876 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1061 + Ifges(2,5) * qJDD(1) - t1080 * Ifges(2,6) - t1066 * t880 + t1070 * t882 + (-t1068 * t886 - t1072 * t887) * qJ(2);
t875 = mrSges(2,1) * g(3) + mrSges(2,3) * t1062 + t1080 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t886 - t1068 * t878 + (qJ(2) * t894 + t1066 * t882 + t1070 * t880) * t1072;
t1 = [-m(1) * g(1) + t1099; -m(1) * g(2) + t1118; (-m(1) - m(2)) * g(3) + t886; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1118 - t1076 * t875 + t1079 * t876; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1099 + t1076 * t876 + t1079 * t875; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1087; t1087; t900; t888; t931; t1081; t1082;];
tauJB  = t1;
