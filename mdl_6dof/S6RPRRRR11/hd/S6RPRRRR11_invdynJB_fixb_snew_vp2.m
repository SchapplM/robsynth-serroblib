% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR11
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
% Datum: 2019-05-06 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 05:36:41
% EndTime: 2019-05-06 05:38:58
% DurationCPUTime: 141.31s
% Computational Cost: add. (2211447->406), mult. (6857487->552), div. (0->0), fcn. (5913505->16), ass. (0->187)
t1062 = sin(pkin(13));
t1064 = sin(pkin(6));
t1065 = cos(pkin(13));
t1067 = cos(pkin(6));
t1071 = sin(qJ(3));
t1066 = cos(pkin(7));
t1076 = cos(qJ(3));
t1108 = t1066 * t1076;
t1063 = sin(pkin(7));
t1113 = t1063 * t1076;
t1082 = t1064 * (-t1062 * t1071 + t1065 * t1108) + t1067 * t1113;
t1035 = t1082 * qJD(1);
t1109 = t1066 * t1071;
t1114 = t1063 * t1071;
t1084 = t1067 * t1114 + (t1062 * t1076 + t1065 * t1109) * t1064;
t1036 = t1084 * qJD(1);
t1021 = -t1036 * qJD(3) + qJDD(1) * t1082;
t1111 = t1064 * t1066;
t1047 = (t1063 * t1067 + t1065 * t1111) * qJD(1) * pkin(9);
t1072 = sin(qJ(1));
t1077 = cos(qJ(1));
t1059 = -g(1) * t1077 - g(2) * t1072;
t1078 = qJD(1) ^ 2;
t1119 = qJ(2) * t1064;
t1051 = -pkin(1) * t1078 + qJDD(1) * t1119 + t1059;
t1123 = pkin(9) * t1062;
t1096 = -pkin(2) * t1065 - t1063 * t1123;
t1118 = qJD(1) * t1064;
t1121 = pkin(9) * qJDD(1);
t1093 = qJD(1) * t1096 * t1118 + t1066 * t1121;
t1058 = t1072 * g(1) - g(2) * t1077;
t1050 = qJDD(1) * pkin(1) + t1078 * t1119 + t1058;
t1107 = qJD(2) * t1118;
t1110 = t1065 * t1067;
t1112 = t1064 * t1065;
t1097 = -g(3) * t1112 + t1050 * t1110 - 0.2e1 * t1062 * t1107;
t1001 = (pkin(2) * qJDD(1) + qJD(1) * t1047) * t1067 + (-t1064 * t1093 - t1051) * t1062 + t1097;
t1052 = (pkin(2) * t1067 - t1111 * t1123) * qJD(1);
t1115 = t1062 * t1067;
t1105 = t1050 * t1115 + (t1051 + 0.2e1 * t1107) * t1065;
t1002 = (-qJD(1) * t1052 + t1063 * t1121) * t1067 + (-g(3) * t1062 + t1065 * t1093) * t1064 + t1105;
t1104 = -t1067 * g(3) + qJDD(2);
t1011 = (-t1050 + t1096 * qJDD(1) + (-t1047 * t1065 + t1052 * t1062) * qJD(1)) * t1064 + t1104;
t969 = -t1071 * t1002 + (t1001 * t1066 + t1011 * t1063) * t1076;
t1122 = Ifges(3,3) * t1067;
t1026 = -t1062 * t1051 + t1097;
t1100 = -mrSges(3,1) * t1065 + mrSges(3,2) * t1062;
t1049 = t1100 * t1118;
t1094 = -mrSges(3,2) * t1067 + mrSges(3,3) * t1112;
t1054 = t1094 * qJD(1);
t1116 = t1062 * t1064;
t1095 = mrSges(3,1) * t1067 - mrSges(3,3) * t1116;
t1019 = -mrSges(4,1) * t1035 + mrSges(4,2) * t1036;
t1091 = -t1063 * t1112 + t1066 * t1067;
t1048 = qJD(1) * t1091 + qJD(3);
t1031 = mrSges(4,1) * t1048 - mrSges(4,3) * t1036;
t1045 = qJDD(1) * t1091 + qJDD(3);
t1070 = sin(qJ(4));
t1075 = cos(qJ(4));
t1028 = -t1070 * t1036 + t1048 * t1075;
t1029 = t1036 * t1075 + t1048 * t1070;
t1003 = -mrSges(5,1) * t1028 + mrSges(5,2) * t1029;
t1034 = qJD(4) - t1035;
t1013 = mrSges(5,1) * t1034 - mrSges(5,3) * t1029;
t1018 = qJDD(4) - t1021;
t1069 = sin(qJ(5));
t1074 = cos(qJ(5));
t1010 = t1029 * t1074 + t1034 * t1069;
t1025 = qJD(5) - t1028;
t1068 = sin(qJ(6));
t1073 = cos(qJ(6));
t1023 = qJD(6) + t1025;
t1009 = -t1029 * t1069 + t1034 * t1074;
t1004 = -pkin(4) * t1028 - pkin(11) * t1029;
t1033 = t1034 ^ 2;
t1020 = -pkin(3) * t1035 - pkin(10) * t1036;
t1044 = t1048 ^ 2;
t970 = t1001 * t1109 + t1076 * t1002 + t1011 * t1114;
t961 = -pkin(3) * t1044 + pkin(10) * t1045 + t1020 * t1035 + t970;
t1022 = t1035 * qJD(3) + qJDD(1) * t1084;
t980 = -t1063 * t1001 + t1066 * t1011;
t966 = (-t1035 * t1048 - t1022) * pkin(10) + (t1036 * t1048 - t1021) * pkin(3) + t980;
t951 = t1070 * t966 + t1075 * t961;
t946 = -pkin(4) * t1033 + pkin(11) * t1018 + t1004 * t1028 + t951;
t960 = -t1045 * pkin(3) - t1044 * pkin(10) + t1036 * t1020 - t969;
t996 = -t1029 * qJD(4) - t1070 * t1022 + t1045 * t1075;
t997 = qJD(4) * t1028 + t1022 * t1075 + t1045 * t1070;
t949 = (-t1028 * t1034 - t997) * pkin(11) + (t1029 * t1034 - t996) * pkin(4) + t960;
t941 = -t1069 * t946 + t1074 * t949;
t973 = qJD(5) * t1009 + t1018 * t1069 + t1074 * t997;
t994 = qJDD(5) - t996;
t939 = (t1009 * t1025 - t973) * pkin(12) + (t1009 * t1010 + t994) * pkin(5) + t941;
t1008 = t1009 ^ 2;
t942 = t1069 * t949 + t1074 * t946;
t972 = -qJD(5) * t1010 + t1018 * t1074 - t1069 * t997;
t987 = pkin(5) * t1025 - pkin(12) * t1010;
t940 = -pkin(5) * t1008 + pkin(12) * t972 - t1025 * t987 + t942;
t937 = -t1068 * t940 + t1073 * t939;
t981 = t1009 * t1073 - t1010 * t1068;
t956 = qJD(6) * t981 + t1068 * t972 + t1073 * t973;
t982 = t1009 * t1068 + t1010 * t1073;
t968 = -mrSges(7,1) * t981 + mrSges(7,2) * t982;
t974 = -mrSges(7,2) * t1023 + mrSges(7,3) * t981;
t989 = qJDD(6) + t994;
t933 = m(7) * t937 + mrSges(7,1) * t989 - mrSges(7,3) * t956 + t1023 * t974 - t968 * t982;
t938 = t1068 * t939 + t1073 * t940;
t955 = -qJD(6) * t982 - t1068 * t973 + t1073 * t972;
t975 = mrSges(7,1) * t1023 - mrSges(7,3) * t982;
t934 = m(7) * t938 - mrSges(7,2) * t989 + mrSges(7,3) * t955 - t1023 * t975 + t968 * t981;
t925 = t1068 * t934 + t1073 * t933;
t983 = -mrSges(6,1) * t1009 + mrSges(6,2) * t1010;
t985 = -mrSges(6,2) * t1025 + mrSges(6,3) * t1009;
t923 = m(6) * t941 + mrSges(6,1) * t994 - mrSges(6,3) * t973 - t1010 * t983 + t1025 * t985 + t925;
t1103 = -t1068 * t933 + t1073 * t934;
t986 = mrSges(6,1) * t1025 - mrSges(6,3) * t1010;
t924 = m(6) * t942 - mrSges(6,2) * t994 + mrSges(6,3) * t972 + t1009 * t983 - t1025 * t986 + t1103;
t921 = -t1069 * t923 + t1074 * t924;
t919 = m(5) * t951 - mrSges(5,2) * t1018 + mrSges(5,3) * t996 + t1003 * t1028 - t1013 * t1034 + t921;
t1012 = -mrSges(5,2) * t1034 + mrSges(5,3) * t1028;
t950 = -t1070 * t961 + t1075 * t966;
t945 = -pkin(4) * t1018 - pkin(11) * t1033 + t1029 * t1004 - t950;
t943 = -pkin(5) * t972 - pkin(12) * t1008 + t1010 * t987 + t945;
t1089 = m(7) * t943 - t955 * mrSges(7,1) + mrSges(7,2) * t956 - t981 * t974 + t975 * t982;
t935 = -m(6) * t945 + t972 * mrSges(6,1) - mrSges(6,2) * t973 + t1009 * t985 - t1010 * t986 - t1089;
t929 = m(5) * t950 + mrSges(5,1) * t1018 - mrSges(5,3) * t997 - t1003 * t1029 + t1012 * t1034 + t935;
t1102 = -t1070 * t929 + t1075 * t919;
t908 = m(4) * t970 - mrSges(4,2) * t1045 + mrSges(4,3) * t1021 + t1019 * t1035 - t1031 * t1048 + t1102;
t1030 = -mrSges(4,2) * t1048 + mrSges(4,3) * t1035;
t911 = t1070 * t919 + t1075 * t929;
t910 = m(4) * t980 - mrSges(4,1) * t1021 + mrSges(4,2) * t1022 - t1030 * t1035 + t1031 * t1036 + t911;
t920 = t1069 * t924 + t1074 * t923;
t1081 = -m(5) * t960 + t996 * mrSges(5,1) - mrSges(5,2) * t997 + t1028 * t1012 - t1013 * t1029 - t920;
t916 = m(4) * t969 + mrSges(4,1) * t1045 - mrSges(4,3) * t1022 - t1019 * t1036 + t1030 * t1048 + t1081;
t897 = -t1063 * t910 + t916 * t1108 + t908 * t1109;
t893 = m(3) * t1026 + t1095 * qJDD(1) + (-t1049 * t1116 + t1054 * t1067) * qJD(1) + t897;
t1037 = -t1064 * t1050 + t1104;
t1053 = t1095 * qJD(1);
t896 = t1066 * t910 + t916 * t1113 + t908 * t1114;
t895 = m(3) * t1037 + (t1100 * qJDD(1) + (t1053 * t1062 - t1054 * t1065) * qJD(1)) * t1064 + t896;
t1027 = -g(3) * t1116 + t1105;
t903 = -t1071 * t916 + t1076 * t908;
t902 = m(3) * t1027 + t1094 * qJDD(1) + (t1049 * t1112 - t1053 * t1067) * qJD(1) + t903;
t882 = -t1064 * t895 + t893 * t1110 + t902 * t1115;
t879 = m(2) * t1058 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1078 + t882;
t889 = -t1062 * t893 + t1065 * t902;
t887 = m(2) * t1059 - mrSges(2,1) * t1078 - qJDD(1) * mrSges(2,2) + t889;
t1120 = t1072 * t887 + t1077 * t879;
t881 = t1067 * t895 + t893 * t1112 + t902 * t1116;
t1101 = -t1072 * t879 + t1077 * t887;
t1099 = Ifges(3,5) * t1062 + Ifges(3,6) * t1065;
t1014 = Ifges(4,5) * t1036 + Ifges(4,6) * t1035 + Ifges(4,3) * t1048;
t1015 = Ifges(4,4) * t1036 + Ifges(4,2) * t1035 + Ifges(4,6) * t1048;
t963 = Ifges(7,5) * t982 + Ifges(7,6) * t981 + Ifges(7,3) * t1023;
t965 = Ifges(7,1) * t982 + Ifges(7,4) * t981 + Ifges(7,5) * t1023;
t926 = -mrSges(7,1) * t943 + mrSges(7,3) * t938 + Ifges(7,4) * t956 + Ifges(7,2) * t955 + Ifges(7,6) * t989 + t1023 * t965 - t963 * t982;
t964 = Ifges(7,4) * t982 + Ifges(7,2) * t981 + Ifges(7,6) * t1023;
t927 = mrSges(7,2) * t943 - mrSges(7,3) * t937 + Ifges(7,1) * t956 + Ifges(7,4) * t955 + Ifges(7,5) * t989 - t1023 * t964 + t963 * t981;
t976 = Ifges(6,5) * t1010 + Ifges(6,6) * t1009 + Ifges(6,3) * t1025;
t978 = Ifges(6,1) * t1010 + Ifges(6,4) * t1009 + Ifges(6,5) * t1025;
t912 = -mrSges(6,1) * t945 + mrSges(6,3) * t942 + Ifges(6,4) * t973 + Ifges(6,2) * t972 + Ifges(6,6) * t994 - pkin(5) * t1089 + pkin(12) * t1103 - t1010 * t976 + t1025 * t978 + t1068 * t927 + t1073 * t926;
t977 = Ifges(6,4) * t1010 + Ifges(6,2) * t1009 + Ifges(6,6) * t1025;
t913 = mrSges(6,2) * t945 - mrSges(6,3) * t941 + Ifges(6,1) * t973 + Ifges(6,4) * t972 + Ifges(6,5) * t994 - pkin(12) * t925 + t1009 * t976 - t1025 * t977 - t1068 * t926 + t1073 * t927;
t990 = Ifges(5,5) * t1029 + Ifges(5,6) * t1028 + Ifges(5,3) * t1034;
t991 = Ifges(5,4) * t1029 + Ifges(5,2) * t1028 + Ifges(5,6) * t1034;
t898 = mrSges(5,2) * t960 - mrSges(5,3) * t950 + Ifges(5,1) * t997 + Ifges(5,4) * t996 + Ifges(5,5) * t1018 - pkin(11) * t920 + t1028 * t990 - t1034 * t991 - t1069 * t912 + t1074 * t913;
t1085 = -mrSges(7,1) * t937 + mrSges(7,2) * t938 - Ifges(7,5) * t956 - Ifges(7,6) * t955 - Ifges(7,3) * t989 - t982 * t964 + t981 * t965;
t1079 = mrSges(6,1) * t941 - mrSges(6,2) * t942 + Ifges(6,5) * t973 + Ifges(6,6) * t972 + Ifges(6,3) * t994 + pkin(5) * t925 - t1009 * t978 + t1010 * t977 - t1085;
t992 = Ifges(5,1) * t1029 + Ifges(5,4) * t1028 + Ifges(5,5) * t1034;
t904 = -mrSges(5,1) * t960 + mrSges(5,3) * t951 + Ifges(5,4) * t997 + Ifges(5,2) * t996 + Ifges(5,6) * t1018 - pkin(4) * t920 - t1029 * t990 + t1034 * t992 - t1079;
t884 = mrSges(4,2) * t980 - mrSges(4,3) * t969 + Ifges(4,1) * t1022 + Ifges(4,4) * t1021 + Ifges(4,5) * t1045 - pkin(10) * t911 + t1014 * t1035 - t1015 * t1048 - t1070 * t904 + t1075 * t898;
t1016 = Ifges(4,1) * t1036 + Ifges(4,4) * t1035 + Ifges(4,5) * t1048;
t1080 = mrSges(5,1) * t950 - mrSges(5,2) * t951 + Ifges(5,5) * t997 + Ifges(5,6) * t996 + Ifges(5,3) * t1018 + pkin(4) * t935 + pkin(11) * t921 - t1028 * t992 + t1029 * t991 + t1069 * t913 + t1074 * t912;
t890 = -mrSges(4,1) * t980 + mrSges(4,3) * t970 + Ifges(4,4) * t1022 + Ifges(4,2) * t1021 + Ifges(4,6) * t1045 - pkin(3) * t911 - t1036 * t1014 + t1048 * t1016 - t1080;
t1090 = pkin(9) * t903 + t1071 * t884 + t1076 * t890;
t1086 = Ifges(3,6) * t1067 + (Ifges(3,4) * t1062 + Ifges(3,2) * t1065) * t1064;
t1041 = t1086 * qJD(1);
t1087 = Ifges(3,5) * t1067 + (Ifges(3,1) * t1062 + Ifges(3,4) * t1065) * t1064;
t1042 = t1087 * qJD(1);
t883 = mrSges(4,1) * t969 - mrSges(4,2) * t970 + Ifges(4,5) * t1022 + Ifges(4,6) * t1021 + Ifges(4,3) * t1045 + pkin(3) * t1081 + pkin(10) * t1102 + t1036 * t1015 - t1035 * t1016 + t1070 * t898 + t1075 * t904;
t873 = qJDD(1) * t1122 + mrSges(3,1) * t1026 - mrSges(3,2) * t1027 + pkin(2) * t897 + t1066 * t883 + t1090 * t1063 + (t1099 * qJDD(1) + (t1041 * t1062 - t1042 * t1065) * qJD(1)) * t1064;
t1040 = (t1064 * t1099 + t1122) * qJD(1);
t875 = -mrSges(3,1) * t1037 + mrSges(3,3) * t1027 - pkin(2) * t896 - t1063 * t883 + (-t1040 * t1116 + t1042 * t1067) * qJD(1) + t1090 * t1066 + t1086 * qJDD(1);
t877 = mrSges(3,2) * t1037 - mrSges(3,3) * t1026 - t1071 * t890 + t1076 * t884 + (t1040 * t1112 - t1041 * t1067) * qJD(1) + (-t1063 * t896 - t1066 * t897) * pkin(9) + t1087 * qJDD(1);
t1088 = mrSges(2,1) * t1058 - mrSges(2,2) * t1059 + Ifges(2,3) * qJDD(1) + pkin(1) * t882 + t1067 * t873 + t875 * t1112 + t877 * t1116 + t889 * t1119;
t871 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1058 + Ifges(2,5) * qJDD(1) - t1078 * Ifges(2,6) - t1062 * t875 + t1065 * t877 + (-t1064 * t881 - t1067 * t882) * qJ(2);
t870 = mrSges(2,1) * g(3) + mrSges(2,3) * t1059 + t1078 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t881 - t1064 * t873 + (qJ(2) * t889 + t1062 * t877 + t1065 * t875) * t1067;
t1 = [-m(1) * g(1) + t1101; -m(1) * g(2) + t1120; (-m(1) - m(2)) * g(3) + t881; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1120 - t1072 * t870 + t1077 * t871; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1101 + t1072 * t871 + t1077 * t870; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1088; t1088; t895; t883; t1080; t1079; -t1085;];
tauJB  = t1;
