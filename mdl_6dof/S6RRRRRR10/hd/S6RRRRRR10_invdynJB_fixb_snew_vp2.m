% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 18:28:55
% EndTime: 2019-05-08 18:37:44
% DurationCPUTime: 472.95s
% Computational Cost: add. (7623451->434), mult. (19545963->586), div. (0->0), fcn. (16939649->18), ass. (0->197)
t1077 = sin(pkin(8));
t1080 = cos(pkin(8));
t1085 = sin(qJ(4));
t1091 = cos(qJ(4));
t1082 = cos(pkin(6));
t1074 = qJD(1) * t1082 + qJD(2);
t1086 = sin(qJ(3));
t1087 = sin(qJ(2));
t1092 = cos(qJ(3));
t1081 = cos(pkin(7));
t1093 = cos(qJ(2));
t1116 = t1081 * t1093;
t1078 = sin(pkin(7));
t1122 = t1078 * t1092;
t1079 = sin(pkin(6));
t1127 = qJD(1) * t1079;
t1049 = t1074 * t1122 + (-t1086 * t1087 + t1092 * t1116) * t1127;
t1125 = qJD(1) * t1093;
t1068 = (qJD(2) * t1125 + qJDD(1) * t1087) * t1079;
t1126 = qJD(1) * t1087;
t1069 = (-qJD(2) * t1126 + qJDD(1) * t1093) * t1079;
t1073 = qJDD(1) * t1082 + qJDD(2);
t1102 = t1069 * t1081 + t1073 * t1078;
t1038 = t1049 * qJD(3) + t1092 * t1068 + t1102 * t1086;
t1123 = t1078 * t1086;
t1050 = t1074 * t1123 + (t1086 * t1116 + t1087 * t1092) * t1127;
t1132 = pkin(12) * t1077;
t1039 = -pkin(3) * t1049 - t1050 * t1132;
t1120 = t1079 * t1093;
t1110 = qJD(1) * t1120;
t1059 = t1074 * t1081 - t1078 * t1110 + qJD(3);
t1103 = t1049 * t1080 + t1059 * t1077;
t1042 = t1103 * pkin(12);
t1051 = -t1069 * t1078 + t1073 * t1081 + qJDD(3);
t1131 = pkin(12) * t1080;
t1058 = (t1074 * t1078 + t1081 * t1110) * pkin(11);
t1134 = pkin(11) * t1078;
t1062 = (-pkin(2) * t1093 - t1087 * t1134) * t1127;
t1088 = sin(qJ(1));
t1094 = cos(qJ(1));
t1071 = t1088 * g(1) - g(2) * t1094;
t1095 = qJD(1) ^ 2;
t1135 = pkin(10) * t1079;
t1065 = qJDD(1) * pkin(1) + t1095 * t1135 + t1071;
t1072 = -g(1) * t1094 - g(2) * t1088;
t1066 = -pkin(1) * t1095 + qJDD(1) * t1135 + t1072;
t1114 = t1082 * t1093;
t1107 = t1065 * t1114 - t1087 * t1066;
t1133 = pkin(11) * t1081;
t1023 = -t1068 * t1133 + t1073 * pkin(2) + t1074 * t1058 + (-g(3) * t1093 - t1062 * t1126) * t1079 + t1107;
t1121 = t1079 * t1087;
t1111 = qJD(1) * t1121;
t1061 = pkin(2) * t1074 - t1111 * t1133;
t1115 = t1082 * t1087;
t1112 = t1065 * t1115 + t1093 * t1066;
t1024 = -t1074 * t1061 + (-g(3) * t1087 + t1062 * t1125) * t1079 + t1102 * pkin(11) + t1112;
t1130 = t1082 * g(3);
t1029 = -t1068 * t1134 - t1069 * pkin(2) - t1130 + (-t1065 + (-t1058 * t1093 + t1061 * t1087) * qJD(1)) * t1079;
t1117 = t1081 * t1092;
t997 = t1023 * t1117 - t1024 * t1086 + t1029 * t1122;
t981 = pkin(3) * t1051 - t1038 * t1131 - t1039 * t1050 + t1042 * t1059 + t997;
t1044 = pkin(3) * t1059 - t1050 * t1131;
t1037 = -t1050 * qJD(3) - t1086 * t1068 + t1102 * t1092;
t1104 = t1037 * t1080 + t1051 * t1077;
t1118 = t1081 * t1086;
t998 = t1023 * t1118 + t1092 * t1024 + t1029 * t1123;
t982 = t1104 * pkin(12) + t1049 * t1039 - t1059 * t1044 + t998;
t1012 = -t1023 * t1078 + t1081 * t1029;
t989 = -pkin(3) * t1037 - t1038 * t1132 - t1042 * t1049 + t1044 * t1050 + t1012;
t968 = -t1085 * t982 + (t1077 * t989 + t1080 * t981) * t1091;
t1032 = -t1085 * t1050 + t1103 * t1091;
t1033 = t1091 * t1050 + t1103 * t1085;
t1002 = -t1033 * qJD(4) - t1085 * t1038 + t1104 * t1091;
t1047 = -g(3) * t1120 + t1107;
t1064 = -mrSges(3,2) * t1074 + mrSges(3,3) * t1110;
t1067 = (-mrSges(3,1) * t1093 + mrSges(3,2) * t1087) * t1127;
t1040 = -mrSges(4,1) * t1049 + mrSges(4,2) * t1050;
t1045 = -mrSges(4,2) * t1059 + mrSges(4,3) * t1049;
t1119 = t1080 * t1085;
t1003 = t1032 * qJD(4) + t1091 * t1038 + t1104 * t1085;
t1013 = -mrSges(5,1) * t1032 + mrSges(5,2) * t1033;
t1043 = -t1049 * t1077 + t1059 * t1080 + qJD(4);
t1018 = -mrSges(5,2) * t1043 + mrSges(5,3) * t1032;
t1025 = -t1037 * t1077 + t1051 * t1080 + qJDD(4);
t1084 = sin(qJ(5));
t1090 = cos(qJ(5));
t1016 = -t1033 * t1084 + t1043 * t1090;
t1031 = qJD(5) - t1032;
t1007 = -mrSges(6,2) * t1031 + mrSges(6,3) * t1016;
t1017 = t1033 * t1090 + t1043 * t1084;
t1008 = mrSges(6,1) * t1031 - mrSges(6,3) * t1017;
t1083 = sin(qJ(6));
t1089 = cos(qJ(6));
t1006 = t1017 * t1089 + t1031 * t1083;
t1015 = qJD(6) - t1016;
t1000 = -pkin(5) * t1016 - pkin(14) * t1017;
t1001 = qJDD(5) - t1002;
t1030 = t1031 ^ 2;
t1014 = -pkin(4) * t1032 - pkin(13) * t1033;
t1041 = t1043 ^ 2;
t1124 = t1077 * t1085;
t969 = t1091 * t982 + t981 * t1119 + t989 * t1124;
t965 = -pkin(4) * t1041 + pkin(13) * t1025 + t1014 * t1032 + t969;
t970 = -t1077 * t981 + t1080 * t989;
t967 = (-t1032 * t1043 - t1003) * pkin(13) + (t1033 * t1043 - t1002) * pkin(4) + t970;
t961 = t1084 * t967 + t1090 * t965;
t959 = -pkin(5) * t1030 + pkin(14) * t1001 + t1000 * t1016 + t961;
t964 = -t1025 * pkin(4) - t1041 * pkin(13) + t1033 * t1014 - t968;
t985 = -qJD(5) * t1017 - t1003 * t1084 + t1025 * t1090;
t986 = qJD(5) * t1016 + t1003 * t1090 + t1025 * t1084;
t962 = (-t1016 * t1031 - t986) * pkin(14) + (t1017 * t1031 - t985) * pkin(5) + t964;
t955 = -t1083 * t959 + t1089 * t962;
t1005 = -t1017 * t1083 + t1031 * t1089;
t973 = qJD(6) * t1005 + t1001 * t1083 + t1089 * t986;
t984 = qJDD(6) - t985;
t990 = -mrSges(7,1) * t1005 + mrSges(7,2) * t1006;
t991 = -mrSges(7,2) * t1015 + mrSges(7,3) * t1005;
t953 = m(7) * t955 + mrSges(7,1) * t984 - mrSges(7,3) * t973 - t1006 * t990 + t1015 * t991;
t956 = t1083 * t962 + t1089 * t959;
t972 = -qJD(6) * t1006 + t1001 * t1089 - t1083 * t986;
t992 = mrSges(7,1) * t1015 - mrSges(7,3) * t1006;
t954 = m(7) * t956 - mrSges(7,2) * t984 + mrSges(7,3) * t972 + t1005 * t990 - t1015 * t992;
t946 = t1083 * t954 + t1089 * t953;
t1098 = -m(6) * t964 + t985 * mrSges(6,1) - mrSges(6,2) * t986 + t1016 * t1007 - t1008 * t1017 - t946;
t942 = m(5) * t968 + mrSges(5,1) * t1025 - mrSges(5,3) * t1003 - t1013 * t1033 + t1018 * t1043 + t1098;
t1128 = t1091 * t942;
t1019 = mrSges(5,1) * t1043 - mrSges(5,3) * t1033;
t947 = -t1083 * t953 + t1089 * t954;
t999 = -mrSges(6,1) * t1016 + mrSges(6,2) * t1017;
t945 = m(6) * t961 - mrSges(6,2) * t1001 + mrSges(6,3) * t985 - t1008 * t1031 + t1016 * t999 + t947;
t960 = -t1084 * t965 + t1090 * t967;
t958 = -pkin(5) * t1001 - pkin(14) * t1030 + t1000 * t1017 - t960;
t957 = -m(7) * t958 + t972 * mrSges(7,1) - mrSges(7,2) * t973 + t1005 * t991 - t1006 * t992;
t951 = m(6) * t960 + mrSges(6,1) * t1001 - mrSges(6,3) * t986 + t1007 * t1031 - t1017 * t999 + t957;
t1109 = -t1084 * t951 + t1090 * t945;
t936 = m(5) * t969 - mrSges(5,2) * t1025 + mrSges(5,3) * t1002 + t1013 * t1032 - t1019 * t1043 + t1109;
t939 = t1084 * t945 + t1090 * t951;
t938 = m(5) * t970 - mrSges(5,1) * t1002 + mrSges(5,2) * t1003 - t1018 * t1032 + t1019 * t1033 + t939;
t925 = -t1077 * t938 + t1080 * t1128 + t936 * t1119;
t921 = m(4) * t997 + mrSges(4,1) * t1051 - mrSges(4,3) * t1038 - t1040 * t1050 + t1045 * t1059 + t925;
t1046 = mrSges(4,1) * t1059 - mrSges(4,3) * t1050;
t924 = t1077 * t1128 + t1080 * t938 + t936 * t1124;
t923 = m(4) * t1012 - mrSges(4,1) * t1037 + mrSges(4,2) * t1038 - t1045 * t1049 + t1046 * t1050 + t924;
t930 = -t1085 * t942 + t1091 * t936;
t929 = m(4) * t998 - mrSges(4,2) * t1051 + mrSges(4,3) * t1037 + t1040 * t1049 - t1046 * t1059 + t930;
t910 = -t1078 * t923 + t921 * t1117 + t929 * t1118;
t906 = m(3) * t1047 + mrSges(3,1) * t1073 - mrSges(3,3) * t1068 + t1064 * t1074 - t1067 * t1111 + t910;
t1055 = -t1079 * t1065 - t1130;
t1063 = mrSges(3,1) * t1074 - mrSges(3,3) * t1111;
t909 = t1081 * t923 + t921 * t1122 + t929 * t1123;
t908 = m(3) * t1055 - t1069 * mrSges(3,1) + t1068 * mrSges(3,2) + (t1063 * t1087 - t1064 * t1093) * t1127 + t909;
t1048 = -g(3) * t1121 + t1112;
t915 = -t1086 * t921 + t1092 * t929;
t914 = m(3) * t1048 - mrSges(3,2) * t1073 + mrSges(3,3) * t1069 - t1063 * t1074 + t1067 * t1110 + t915;
t895 = -t1079 * t908 + t906 * t1114 + t914 * t1115;
t892 = m(2) * t1071 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1095 + t895;
t903 = -t1087 * t906 + t1093 * t914;
t901 = m(2) * t1072 - mrSges(2,1) * t1095 - qJDD(1) * mrSges(2,2) + t903;
t1129 = t1088 * t901 + t1094 * t892;
t894 = t1082 * t908 + t906 * t1120 + t914 * t1121;
t1108 = -t1088 * t892 + t1094 * t901;
t1034 = Ifges(4,5) * t1050 + Ifges(4,6) * t1049 + Ifges(4,3) * t1059;
t1036 = Ifges(4,1) * t1050 + Ifges(4,4) * t1049 + Ifges(4,5) * t1059;
t1009 = Ifges(5,5) * t1033 + Ifges(5,6) * t1032 + Ifges(5,3) * t1043;
t1010 = Ifges(5,4) * t1033 + Ifges(5,2) * t1032 + Ifges(5,6) * t1043;
t974 = Ifges(7,5) * t1006 + Ifges(7,6) * t1005 + Ifges(7,3) * t1015;
t976 = Ifges(7,1) * t1006 + Ifges(7,4) * t1005 + Ifges(7,5) * t1015;
t948 = -mrSges(7,1) * t958 + mrSges(7,3) * t956 + Ifges(7,4) * t973 + Ifges(7,2) * t972 + Ifges(7,6) * t984 - t1006 * t974 + t1015 * t976;
t975 = Ifges(7,4) * t1006 + Ifges(7,2) * t1005 + Ifges(7,6) * t1015;
t949 = mrSges(7,2) * t958 - mrSges(7,3) * t955 + Ifges(7,1) * t973 + Ifges(7,4) * t972 + Ifges(7,5) * t984 + t1005 * t974 - t1015 * t975;
t993 = Ifges(6,5) * t1017 + Ifges(6,6) * t1016 + Ifges(6,3) * t1031;
t994 = Ifges(6,4) * t1017 + Ifges(6,2) * t1016 + Ifges(6,6) * t1031;
t931 = mrSges(6,2) * t964 - mrSges(6,3) * t960 + Ifges(6,1) * t986 + Ifges(6,4) * t985 + Ifges(6,5) * t1001 - pkin(14) * t946 + t1016 * t993 - t1031 * t994 - t1083 * t948 + t1089 * t949;
t1097 = mrSges(7,1) * t955 - mrSges(7,2) * t956 + Ifges(7,5) * t973 + Ifges(7,6) * t972 + Ifges(7,3) * t984 - t1005 * t976 + t1006 * t975;
t995 = Ifges(6,1) * t1017 + Ifges(6,4) * t1016 + Ifges(6,5) * t1031;
t932 = -mrSges(6,1) * t964 + mrSges(6,3) * t961 + Ifges(6,4) * t986 + Ifges(6,2) * t985 + Ifges(6,6) * t1001 - pkin(5) * t946 - t1017 * t993 + t1031 * t995 - t1097;
t917 = mrSges(5,2) * t970 - mrSges(5,3) * t968 + Ifges(5,1) * t1003 + Ifges(5,4) * t1002 + Ifges(5,5) * t1025 - pkin(13) * t939 + t1009 * t1032 - t1010 * t1043 - t1084 * t932 + t1090 * t931;
t1011 = Ifges(5,1) * t1033 + Ifges(5,4) * t1032 + Ifges(5,5) * t1043;
t1096 = mrSges(6,1) * t960 - mrSges(6,2) * t961 + Ifges(6,5) * t986 + Ifges(6,6) * t985 + Ifges(6,3) * t1001 + pkin(5) * t957 + pkin(14) * t947 - t1016 * t995 + t1017 * t994 + t1083 * t949 + t1089 * t948;
t918 = -mrSges(5,1) * t970 + mrSges(5,3) * t969 + Ifges(5,4) * t1003 + Ifges(5,2) * t1002 + Ifges(5,6) * t1025 - pkin(4) * t939 - t1033 * t1009 + t1043 * t1011 - t1096;
t1100 = pkin(12) * t930 + t1085 * t917 + t1091 * t918;
t916 = mrSges(5,1) * t968 - mrSges(5,2) * t969 + Ifges(5,5) * t1003 + Ifges(5,6) * t1002 + Ifges(5,3) * t1025 + pkin(4) * t1098 + pkin(13) * t1109 + t1033 * t1010 - t1032 * t1011 + t1084 * t931 + t1090 * t932;
t897 = -mrSges(4,1) * t1012 + mrSges(4,3) * t998 + Ifges(4,4) * t1038 + Ifges(4,2) * t1037 + Ifges(4,6) * t1051 - pkin(3) * t924 - t1050 * t1034 + t1059 * t1036 - t1077 * t916 + t1100 * t1080;
t1035 = Ifges(4,4) * t1050 + Ifges(4,2) * t1049 + Ifges(4,6) * t1059;
t898 = mrSges(4,2) * t1012 - mrSges(4,3) * t997 + Ifges(4,1) * t1038 + Ifges(4,4) * t1037 + Ifges(4,5) * t1051 + t1049 * t1034 - t1059 * t1035 - t1085 * t918 + t1091 * t917 + (-t1077 * t924 - t1080 * t925) * pkin(12);
t1101 = pkin(11) * t915 + t1086 * t898 + t1092 * t897;
t1053 = Ifges(3,6) * t1074 + (Ifges(3,4) * t1087 + Ifges(3,2) * t1093) * t1127;
t1054 = Ifges(3,5) * t1074 + (Ifges(3,1) * t1087 + Ifges(3,4) * t1093) * t1127;
t896 = mrSges(4,1) * t997 - mrSges(4,2) * t998 + Ifges(4,5) * t1038 + Ifges(4,6) * t1037 + Ifges(4,3) * t1051 + pkin(3) * t925 + t1050 * t1035 - t1049 * t1036 + t1100 * t1077 + t1080 * t916;
t886 = mrSges(3,1) * t1047 - mrSges(3,2) * t1048 + Ifges(3,5) * t1068 + Ifges(3,6) * t1069 + Ifges(3,3) * t1073 + pkin(2) * t910 + t1081 * t896 + (t1053 * t1087 - t1054 * t1093) * t1127 + t1101 * t1078;
t1052 = Ifges(3,3) * t1074 + (Ifges(3,5) * t1087 + Ifges(3,6) * t1093) * t1127;
t888 = -mrSges(3,1) * t1055 + mrSges(3,3) * t1048 + Ifges(3,4) * t1068 + Ifges(3,2) * t1069 + Ifges(3,6) * t1073 - pkin(2) * t909 - t1052 * t1111 + t1074 * t1054 - t1078 * t896 + t1101 * t1081;
t890 = t1052 * t1110 + mrSges(3,2) * t1055 - mrSges(3,3) * t1047 + Ifges(3,1) * t1068 + Ifges(3,4) * t1069 + Ifges(3,5) * t1073 - t1074 * t1053 - t1086 * t897 + t1092 * t898 + (-t1078 * t909 - t1081 * t910) * pkin(11);
t1099 = mrSges(2,1) * t1071 - mrSges(2,2) * t1072 + Ifges(2,3) * qJDD(1) + pkin(1) * t895 + t1082 * t886 + t888 * t1120 + t890 * t1121 + t903 * t1135;
t884 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1071 + Ifges(2,5) * qJDD(1) - t1095 * Ifges(2,6) - t1087 * t888 + t1093 * t890 + (-t1079 * t894 - t1082 * t895) * pkin(10);
t883 = mrSges(2,1) * g(3) + mrSges(2,3) * t1072 + t1095 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t894 - t1079 * t886 + (pkin(10) * t903 + t1087 * t890 + t1093 * t888) * t1082;
t1 = [-m(1) * g(1) + t1108; -m(1) * g(2) + t1129; (-m(1) - m(2)) * g(3) + t894; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(9) * t1129 - t1088 * t883 + t1094 * t884; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(9) * t1108 + t1088 * t884 + t1094 * t883; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1099; t1099; t886; t896; t916; t1096; t1097;];
tauJB  = t1;
