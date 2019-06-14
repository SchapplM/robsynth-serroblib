% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 17:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR15_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:09:02
% EndTime: 2019-05-07 17:10:04
% DurationCPUTime: 58.05s
% Computational Cost: add. (948555->403), mult. (2375290->520), div. (0->0), fcn. (1972460->14), ass. (0->178)
t1146 = Ifges(4,5) - Ifges(5,4);
t1152 = -Ifges(4,2) - Ifges(5,3);
t1151 = Ifges(5,2) + Ifges(4,1);
t1145 = Ifges(4,6) - Ifges(5,5);
t1144 = -Ifges(5,6) - Ifges(4,4);
t1150 = Ifges(4,3) + Ifges(5,1);
t1090 = sin(qJ(3));
t1087 = sin(pkin(6));
t1091 = sin(qJ(2));
t1132 = t1087 * t1091;
t1122 = qJD(1) * t1132;
t1148 = cos(qJ(3));
t1095 = cos(qJ(2));
t1086 = sin(pkin(7));
t1141 = cos(pkin(6));
t1111 = qJD(1) * t1141 + qJD(2);
t1109 = t1111 * t1086;
t1138 = qJD(1) * t1087;
t1140 = cos(pkin(7));
t1112 = t1140 * t1138;
t1149 = t1095 * t1112 + t1109;
t1049 = t1090 * t1122 - t1148 * t1149;
t1120 = t1090 * t1140;
t1050 = t1090 * t1109 + (t1091 * t1148 + t1095 * t1120) * t1138;
t1031 = pkin(3) * t1049 - qJ(4) * t1050;
t1137 = qJD(1) * t1091;
t1075 = (-qJD(2) * t1137 + qJDD(1) * t1095) * t1087;
t1083 = qJDD(1) * t1141 + qJDD(2);
t1053 = -t1086 * t1075 + t1083 * t1140 + qJDD(3);
t1131 = t1087 * t1095;
t1121 = qJD(1) * t1131;
t1064 = t1086 * t1121 - t1111 * t1140 - qJD(3);
t1061 = t1064 ^ 2;
t1062 = t1149 * pkin(10);
t1142 = pkin(10) * t1091;
t1068 = (-pkin(2) * t1095 - t1086 * t1142) * t1138;
t1136 = qJD(1) * t1095;
t1074 = (qJD(2) * t1136 + qJDD(1) * t1091) * t1087;
t1092 = sin(qJ(1));
t1096 = cos(qJ(1));
t1081 = t1092 * g(1) - g(2) * t1096;
t1097 = qJD(1) ^ 2;
t1143 = pkin(9) * t1087;
t1071 = qJDD(1) * pkin(1) + t1097 * t1143 + t1081;
t1082 = -g(1) * t1096 - g(2) * t1092;
t1072 = -pkin(1) * t1097 + qJDD(1) * t1143 + t1082;
t1118 = t1095 * t1141;
t1114 = t1071 * t1118 - t1091 * t1072;
t1124 = t1140 * pkin(10);
t1002 = -t1074 * t1124 + t1083 * pkin(2) + t1111 * t1062 + (-g(3) * t1095 - t1068 * t1137) * t1087 + t1114;
t1067 = pkin(2) * t1111 - t1112 * t1142;
t1108 = t1075 * t1140 + t1083 * t1086;
t1119 = t1091 * t1141;
t1126 = t1071 * t1119 + t1095 * t1072;
t1003 = -t1111 * t1067 + (-g(3) * t1091 + t1068 * t1136) * t1087 + t1108 * pkin(10) + t1126;
t1125 = t1141 * g(3);
t1008 = -t1074 * t1086 * pkin(10) - t1125 - t1075 * pkin(2) + (-t1071 + (-t1062 * t1095 + t1067 * t1091) * qJD(1)) * t1087;
t1133 = t1086 * t1090;
t979 = t1002 * t1120 + t1148 * t1003 + t1008 * t1133;
t973 = t1061 * pkin(3) - t1053 * qJ(4) + 0.2e1 * qJD(4) * t1064 + t1049 * t1031 - t979;
t1147 = mrSges(4,1) - mrSges(5,2);
t1043 = -g(3) * t1131 + t1114;
t1070 = -mrSges(3,2) * t1111 + mrSges(3,3) * t1121;
t1073 = (-mrSges(3,1) * t1095 + mrSges(3,2) * t1091) * t1138;
t1113 = t1140 * t1148;
t1123 = t1086 * t1148;
t1027 = t1050 * qJD(3) + t1090 * t1074 - t1075 * t1113 - t1083 * t1123;
t1028 = -t1049 * qJD(3) + t1074 * t1148 + t1108 * t1090;
t1040 = -mrSges(4,1) * t1064 - mrSges(4,3) * t1050;
t1039 = mrSges(5,1) * t1050 - mrSges(5,2) * t1064;
t1089 = sin(qJ(5));
t1094 = cos(qJ(5));
t1035 = t1049 * t1094 + t1064 * t1089;
t1036 = t1049 * t1089 - t1064 * t1094;
t1004 = -mrSges(6,1) * t1035 + mrSges(6,2) * t1036;
t1047 = qJD(5) + t1050;
t1015 = mrSges(6,1) * t1047 - mrSges(6,3) * t1036;
t1026 = qJDD(5) + t1028;
t1088 = sin(qJ(6));
t1093 = cos(qJ(6));
t1013 = t1036 * t1093 + t1047 * t1088;
t1034 = qJD(6) - t1035;
t1005 = -pkin(5) * t1035 - pkin(12) * t1036;
t1046 = t1047 ^ 2;
t1134 = t1049 * t1064;
t978 = t1002 * t1113 - t1090 * t1003 + t1008 * t1123;
t974 = -t1053 * pkin(3) - t1061 * qJ(4) + t1050 * t1031 + qJDD(4) - t978;
t967 = (t1049 * t1050 - t1053) * pkin(11) + (t1028 - t1134) * pkin(4) + t974;
t1041 = pkin(4) * t1050 + pkin(11) * t1064;
t1048 = t1049 ^ 2;
t984 = -t1086 * t1002 + t1140 * t1008;
t1101 = (-t1028 - t1134) * qJ(4) + t984 + (-pkin(3) * t1064 - 0.2e1 * qJD(4)) * t1050;
t968 = -t1048 * pkin(4) - t1050 * t1041 + (pkin(3) + pkin(11)) * t1027 + t1101;
t963 = t1089 * t967 + t1094 * t968;
t960 = -pkin(5) * t1046 + pkin(12) * t1026 + t1005 * t1035 + t963;
t970 = -t1027 * pkin(4) - t1048 * pkin(11) - t1064 * t1041 - t973;
t991 = -qJD(5) * t1036 + t1027 * t1094 - t1053 * t1089;
t992 = qJD(5) * t1035 + t1027 * t1089 + t1053 * t1094;
t964 = (-t1035 * t1047 - t992) * pkin(12) + (t1036 * t1047 - t991) * pkin(5) + t970;
t956 = -t1088 * t960 + t1093 * t964;
t1012 = -t1036 * t1088 + t1047 * t1093;
t977 = qJD(6) * t1012 + t1026 * t1088 + t1093 * t992;
t986 = -mrSges(7,1) * t1012 + mrSges(7,2) * t1013;
t990 = qJDD(6) - t991;
t993 = -mrSges(7,2) * t1034 + mrSges(7,3) * t1012;
t954 = m(7) * t956 + mrSges(7,1) * t990 - mrSges(7,3) * t977 - t1013 * t986 + t1034 * t993;
t957 = t1088 * t964 + t1093 * t960;
t976 = -qJD(6) * t1013 + t1026 * t1093 - t1088 * t992;
t994 = mrSges(7,1) * t1034 - mrSges(7,3) * t1013;
t955 = m(7) * t957 - mrSges(7,2) * t990 + mrSges(7,3) * t976 + t1012 * t986 - t1034 * t994;
t1117 = -t1088 * t954 + t1093 * t955;
t943 = m(6) * t963 - mrSges(6,2) * t1026 + mrSges(6,3) * t991 + t1004 * t1035 - t1015 * t1047 + t1117;
t1014 = -mrSges(6,2) * t1047 + mrSges(6,3) * t1035;
t962 = -t1089 * t968 + t1094 * t967;
t959 = -pkin(5) * t1026 - pkin(12) * t1046 + t1005 * t1036 - t962;
t1105 = -m(7) * t959 + t976 * mrSges(7,1) - mrSges(7,2) * t977 + t1012 * t993 - t1013 * t994;
t950 = m(6) * t962 + mrSges(6,1) * t1026 - mrSges(6,3) * t992 - t1004 * t1036 + t1014 * t1047 + t1105;
t1116 = -t1089 * t950 + t1094 * t943;
t972 = t1027 * pkin(3) + t1101;
t1107 = m(5) * t972 - t1028 * mrSges(5,3) - t1050 * t1039 + t1116;
t1038 = mrSges(5,1) * t1049 + mrSges(5,3) * t1064;
t1127 = mrSges(4,2) * t1064 - mrSges(4,3) * t1049 - t1038;
t931 = m(4) * t984 + t1028 * mrSges(4,2) + t1027 * t1147 + t1050 * t1040 + t1049 * t1127 + t1107;
t1032 = mrSges(4,1) * t1049 + mrSges(4,2) * t1050;
t1033 = -mrSges(5,2) * t1049 - mrSges(5,3) * t1050;
t937 = t1089 * t943 + t1094 * t950;
t1103 = -m(5) * t974 - t1028 * mrSges(5,1) - t1050 * t1033 - t937;
t934 = m(4) * t978 - t1028 * mrSges(4,3) - t1050 * t1032 + t1053 * t1147 - t1064 * t1127 + t1103;
t945 = t1088 * t955 + t1093 * t954;
t1102 = -m(6) * t970 + t991 * mrSges(6,1) - t992 * mrSges(6,2) + t1035 * t1014 - t1036 * t1015 - t945;
t1099 = -m(5) * t973 + t1053 * mrSges(5,3) - t1064 * t1039 - t1102;
t941 = t1099 + m(4) * t979 - t1053 * mrSges(4,2) + (-mrSges(4,3) - mrSges(5,1)) * t1027 + (-t1032 - t1033) * t1049 + t1064 * t1040;
t922 = -t1086 * t931 + t934 * t1113 + t941 * t1120;
t918 = m(3) * t1043 + t1083 * mrSges(3,1) - t1074 * mrSges(3,3) + t1070 * t1111 - t1073 * t1122 + t922;
t1057 = -t1087 * t1071 - t1125;
t1069 = mrSges(3,1) * t1111 - mrSges(3,3) * t1122;
t921 = t934 * t1123 + t941 * t1133 + t1140 * t931;
t920 = m(3) * t1057 - t1075 * mrSges(3,1) + t1074 * mrSges(3,2) + (t1069 * t1091 - t1070 * t1095) * t1138 + t921;
t1044 = -g(3) * t1132 + t1126;
t927 = -t1090 * t934 + t1148 * t941;
t926 = m(3) * t1044 - t1083 * mrSges(3,2) + t1075 * mrSges(3,3) - t1069 * t1111 + t1073 * t1121 + t927;
t907 = -t1087 * t920 + t918 * t1118 + t926 * t1119;
t904 = m(2) * t1081 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1097 + t907;
t914 = -t1091 * t918 + t1095 * t926;
t912 = m(2) * t1082 - mrSges(2,1) * t1097 - qJDD(1) * mrSges(2,2) + t914;
t1139 = t1092 * t912 + t1096 * t904;
t1130 = t1145 * t1049 - t1146 * t1050 + t1150 * t1064;
t1129 = t1152 * t1049 - t1144 * t1050 - t1145 * t1064;
t1128 = -t1144 * t1049 - t1151 * t1050 + t1146 * t1064;
t906 = t918 * t1131 + t926 * t1132 + t1141 * t920;
t1115 = -t1092 * t904 + t1096 * t912;
t1055 = Ifges(3,6) * qJD(2) + (Ifges(3,6) * t1141 + (Ifges(3,4) * t1091 + Ifges(3,2) * t1095) * t1087) * qJD(1);
t1056 = Ifges(3,5) * qJD(2) + (Ifges(3,5) * t1141 + (Ifges(3,1) * t1091 + Ifges(3,4) * t1095) * t1087) * qJD(1);
t980 = Ifges(7,5) * t1013 + Ifges(7,6) * t1012 + Ifges(7,3) * t1034;
t982 = Ifges(7,1) * t1013 + Ifges(7,4) * t1012 + Ifges(7,5) * t1034;
t948 = -mrSges(7,1) * t959 + mrSges(7,3) * t957 + Ifges(7,4) * t977 + Ifges(7,2) * t976 + Ifges(7,6) * t990 - t1013 * t980 + t1034 * t982;
t981 = Ifges(7,4) * t1013 + Ifges(7,2) * t1012 + Ifges(7,6) * t1034;
t949 = mrSges(7,2) * t959 - mrSges(7,3) * t956 + Ifges(7,1) * t977 + Ifges(7,4) * t976 + Ifges(7,5) * t990 + t1012 * t980 - t1034 * t981;
t995 = Ifges(6,5) * t1036 + Ifges(6,6) * t1035 + Ifges(6,3) * t1047;
t996 = Ifges(6,4) * t1036 + Ifges(6,2) * t1035 + Ifges(6,6) * t1047;
t928 = mrSges(6,2) * t970 - mrSges(6,3) * t962 + Ifges(6,1) * t992 + Ifges(6,4) * t991 + Ifges(6,5) * t1026 - pkin(12) * t945 + t1035 * t995 - t1047 * t996 - t1088 * t948 + t1093 * t949;
t1098 = mrSges(7,1) * t956 - mrSges(7,2) * t957 + Ifges(7,5) * t977 + Ifges(7,6) * t976 + Ifges(7,3) * t990 - t1012 * t982 + t1013 * t981;
t997 = Ifges(6,1) * t1036 + Ifges(6,4) * t1035 + Ifges(6,5) * t1047;
t929 = -mrSges(6,1) * t970 + mrSges(6,3) * t963 + Ifges(6,4) * t992 + Ifges(6,2) * t991 + Ifges(6,6) * t1026 - pkin(5) * t945 - t1036 * t995 + t1047 * t997 - t1098;
t936 = t1053 * mrSges(5,2) - t1064 * t1038 - t1103;
t908 = mrSges(4,1) * t978 - mrSges(4,2) * t979 + mrSges(5,2) * t974 - mrSges(5,3) * t973 + t1094 * t928 - t1089 * t929 - pkin(11) * t937 - pkin(3) * t936 + qJ(4) * t1099 + t1150 * t1053 + t1129 * t1050 + (-qJ(4) * t1033 - t1128) * t1049 + t1146 * t1028 + (-mrSges(5,1) * qJ(4) - t1145) * t1027;
t935 = -t1027 * mrSges(5,2) - t1049 * t1038 + t1107;
t909 = -mrSges(4,1) * t984 - mrSges(5,1) * t973 + mrSges(5,2) * t972 + mrSges(4,3) * t979 - pkin(3) * t935 - pkin(4) * t1102 - pkin(11) * t1116 + t1152 * t1027 - t1144 * t1028 + t1130 * t1050 + t1145 * t1053 + t1128 * t1064 - t1089 * t928 - t1094 * t929;
t1100 = mrSges(6,1) * t962 - mrSges(6,2) * t963 + Ifges(6,5) * t992 + Ifges(6,6) * t991 + Ifges(6,3) * t1026 + pkin(5) * t1105 + pkin(12) * t1117 - t1035 * t997 + t1036 * t996 + t1088 * t949 + t1093 * t948;
t915 = mrSges(5,1) * t974 + mrSges(4,2) * t984 - mrSges(4,3) * t978 - mrSges(5,3) * t972 + pkin(4) * t937 - qJ(4) * t935 + t1144 * t1027 + t1151 * t1028 + t1130 * t1049 + t1146 * t1053 + t1129 * t1064 + t1100;
t898 = mrSges(3,1) * t1043 - mrSges(3,2) * t1044 + t1140 * t908 + Ifges(3,5) * t1074 + Ifges(3,6) * t1075 + Ifges(3,3) * t1083 + pkin(2) * t922 + (t1055 * t1091 - t1056 * t1095) * t1138 + (pkin(10) * t927 + t1090 * t915 + t1148 * t909) * t1086;
t1054 = Ifges(3,3) * qJD(2) + (Ifges(3,3) * t1141 + (Ifges(3,5) * t1091 + Ifges(3,6) * t1095) * t1087) * qJD(1);
t900 = -mrSges(3,1) * t1057 + mrSges(3,3) * t1044 + Ifges(3,4) * t1074 + Ifges(3,2) * t1075 + Ifges(3,6) * t1083 - pkin(2) * t921 - t1054 * t1122 + t1056 * t1111 - t1086 * t908 + t1113 * t909 + t1120 * t915 + t1124 * t927;
t902 = Ifges(3,1) * t1074 + Ifges(3,4) * t1075 + Ifges(3,5) * t1083 + t1054 * t1121 - t1111 * t1055 + mrSges(3,2) * t1057 - mrSges(3,3) * t1043 + t1148 * t915 - t1090 * t909 + (-t1086 * t921 - t1140 * t922) * pkin(10);
t1104 = mrSges(2,1) * t1081 - mrSges(2,2) * t1082 + Ifges(2,3) * qJDD(1) + pkin(1) * t907 + t900 * t1131 + t902 * t1132 + t1141 * t898 + t914 * t1143;
t896 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1081 + Ifges(2,5) * qJDD(1) - t1097 * Ifges(2,6) - t1091 * t900 + t1095 * t902 + (-t1087 * t906 - t1141 * t907) * pkin(9);
t895 = pkin(9) * t1141 * t914 + mrSges(2,1) * g(3) + mrSges(2,3) * t1082 + t1097 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t906 - t1087 * t898 + t1118 * t900 + t1119 * t902;
t1 = [-m(1) * g(1) + t1115; -m(1) * g(2) + t1139; (-m(1) - m(2)) * g(3) + t906; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1139 - t1092 * t895 + t1096 * t896; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1115 + t1092 * t896 + t1096 * t895; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1104; t1104; t898; t908; t936; t1100; t1098;];
tauJB  = t1;
