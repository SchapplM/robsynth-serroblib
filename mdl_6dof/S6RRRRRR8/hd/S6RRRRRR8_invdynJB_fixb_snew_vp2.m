% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 13:44:05
% EndTime: 2019-05-08 13:47:09
% DurationCPUTime: 181.21s
% Computational Cost: add. (3012777->419), mult. (7386433->552), div. (0->0), fcn. (6302888->16), ass. (0->183)
t1061 = cos(pkin(6));
t1055 = qJD(1) * t1061 + qJD(2);
t1058 = sin(pkin(7));
t1060 = cos(pkin(7));
t1059 = sin(pkin(6));
t1072 = cos(qJ(2));
t1098 = t1059 * t1072;
t1091 = qJD(1) * t1098;
t1082 = t1055 * t1058 + t1060 * t1091;
t1037 = t1082 * pkin(10);
t1066 = sin(qJ(2));
t1099 = t1059 * t1066;
t1092 = qJD(1) * t1099;
t1107 = pkin(10) * t1060;
t1040 = pkin(2) * t1055 - t1092 * t1107;
t1067 = sin(qJ(1));
t1073 = cos(qJ(1));
t1052 = t1067 * g(1) - g(2) * t1073;
t1074 = qJD(1) ^ 2;
t1109 = pkin(9) * t1059;
t1045 = qJDD(1) * pkin(1) + t1074 * t1109 + t1052;
t1101 = qJD(1) * t1072;
t1048 = (qJD(2) * t1101 + qJDD(1) * t1066) * t1059;
t1102 = qJD(1) * t1066;
t1049 = (-qJD(2) * t1102 + qJDD(1) * t1072) * t1059;
t1106 = t1061 * g(3);
t1108 = pkin(10) * t1058;
t1001 = -t1048 * t1108 - t1049 * pkin(2) - t1106 + (-t1045 + (-t1037 * t1072 + t1040 * t1066) * qJD(1)) * t1059;
t1065 = sin(qJ(3));
t1071 = cos(qJ(3));
t1103 = qJD(1) * t1059;
t1042 = (-pkin(2) * t1072 - t1066 * t1108) * t1103;
t1054 = qJDD(1) * t1061 + qJDD(2);
t1053 = -g(1) * t1073 - g(2) * t1067;
t1046 = -pkin(1) * t1074 + qJDD(1) * t1109 + t1053;
t1095 = t1061 * t1072;
t1086 = t1045 * t1095 - t1066 * t1046;
t994 = -t1048 * t1107 + t1054 * pkin(2) + t1055 * t1037 + (-g(3) * t1072 - t1042 * t1102) * t1059 + t1086;
t1084 = t1049 * t1060 + t1054 * t1058;
t1096 = t1061 * t1066;
t1093 = t1045 * t1096 + t1072 * t1046;
t995 = -t1055 * t1040 + (-g(3) * t1066 + t1042 * t1101) * t1059 + t1084 * pkin(10) + t1093;
t966 = -t1065 * t995 + (t1001 * t1058 + t1060 * t994) * t1071;
t1027 = -t1065 * t1092 + t1082 * t1071;
t1097 = t1060 * t1065;
t1100 = t1058 * t1065;
t1028 = t1055 * t1100 + (t1066 * t1071 + t1072 * t1097) * t1103;
t1012 = -t1028 * qJD(3) - t1065 * t1048 + t1084 * t1071;
t1021 = -g(3) * t1098 + t1086;
t1044 = -mrSges(3,2) * t1055 + mrSges(3,3) * t1091;
t1047 = (-mrSges(3,1) * t1072 + mrSges(3,2) * t1066) * t1103;
t1013 = t1027 * qJD(3) + t1071 * t1048 + t1084 * t1065;
t1014 = -mrSges(4,1) * t1027 + mrSges(4,2) * t1028;
t1038 = t1055 * t1060 - t1058 * t1091 + qJD(3);
t1019 = -mrSges(4,2) * t1038 + mrSges(4,3) * t1027;
t1029 = -t1049 * t1058 + t1054 * t1060 + qJDD(3);
t1064 = sin(qJ(4));
t1070 = cos(qJ(4));
t1017 = -t1028 * t1064 + t1038 * t1070;
t1026 = qJD(4) - t1027;
t1003 = -mrSges(5,2) * t1026 + mrSges(5,3) * t1017;
t1018 = t1028 * t1070 + t1038 * t1064;
t1004 = mrSges(5,1) * t1026 - mrSges(5,3) * t1018;
t1062 = sin(qJ(6));
t1068 = cos(qJ(6));
t1011 = qJDD(4) - t1012;
t1010 = qJDD(5) + t1011;
t1024 = qJD(5) + t1026;
t1023 = t1024 ^ 2;
t1063 = sin(qJ(5));
t1069 = cos(qJ(5));
t1015 = -pkin(3) * t1027 - pkin(11) * t1028;
t1036 = t1038 ^ 2;
t967 = t1001 * t1100 + t1071 * t995 + t994 * t1097;
t951 = -pkin(3) * t1036 + pkin(11) * t1029 + t1015 * t1027 + t967;
t977 = t1060 * t1001 - t1058 * t994;
t954 = (-t1027 * t1038 - t1013) * pkin(11) + (t1028 * t1038 - t1012) * pkin(3) + t977;
t943 = -t1064 * t951 + t1070 * t954;
t980 = qJD(4) * t1017 + t1013 * t1070 + t1029 * t1064;
t940 = (t1017 * t1026 - t980) * pkin(12) + (t1017 * t1018 + t1011) * pkin(4) + t943;
t1005 = pkin(4) * t1026 - pkin(12) * t1018;
t1016 = t1017 ^ 2;
t944 = t1064 * t954 + t1070 * t951;
t979 = -qJD(4) * t1018 - t1013 * t1064 + t1029 * t1070;
t942 = -pkin(4) * t1016 + pkin(12) * t979 - t1005 * t1026 + t944;
t937 = t1063 * t940 + t1069 * t942;
t996 = t1017 * t1069 - t1018 * t1063;
t997 = t1017 * t1063 + t1018 * t1069;
t976 = -pkin(5) * t996 - pkin(13) * t997;
t934 = -pkin(5) * t1023 + pkin(13) * t1010 + t976 * t996 + t937;
t950 = -t1029 * pkin(3) - t1036 * pkin(11) + t1028 * t1015 - t966;
t945 = -t979 * pkin(4) - t1016 * pkin(12) + t1018 * t1005 + t950;
t959 = -qJD(5) * t997 - t1063 * t980 + t1069 * t979;
t960 = qJD(5) * t996 + t1063 * t979 + t1069 * t980;
t938 = (-t1024 * t996 - t960) * pkin(13) + (t1024 * t997 - t959) * pkin(5) + t945;
t931 = -t1062 * t934 + t1068 * t938;
t981 = t1024 * t1068 - t1062 * t997;
t948 = qJD(6) * t981 + t1010 * t1062 + t1068 * t960;
t958 = qJDD(6) - t959;
t982 = t1024 * t1062 + t1068 * t997;
t968 = -mrSges(7,1) * t981 + mrSges(7,2) * t982;
t993 = qJD(6) - t996;
t969 = -mrSges(7,2) * t993 + mrSges(7,3) * t981;
t927 = m(7) * t931 + mrSges(7,1) * t958 - mrSges(7,3) * t948 - t968 * t982 + t969 * t993;
t932 = t1062 * t938 + t1068 * t934;
t947 = -qJD(6) * t982 + t1010 * t1068 - t1062 * t960;
t970 = mrSges(7,1) * t993 - mrSges(7,3) * t982;
t928 = m(7) * t932 - mrSges(7,2) * t958 + mrSges(7,3) * t947 + t968 * t981 - t970 * t993;
t916 = t1062 * t928 + t1068 * t927;
t983 = -mrSges(6,2) * t1024 + mrSges(6,3) * t996;
t984 = mrSges(6,1) * t1024 - mrSges(6,3) * t997;
t1079 = m(6) * t945 - t959 * mrSges(6,1) + mrSges(6,2) * t960 - t996 * t983 + t984 * t997 + t916;
t1076 = -m(5) * t950 + t979 * mrSges(5,1) - mrSges(5,2) * t980 + t1017 * t1003 - t1004 * t1018 - t1079;
t911 = m(4) * t966 + mrSges(4,1) * t1029 - mrSges(4,3) * t1013 - t1014 * t1028 + t1019 * t1038 + t1076;
t1104 = t1071 * t911;
t1020 = mrSges(4,1) * t1038 - mrSges(4,3) * t1028;
t1090 = -t1062 * t927 + t1068 * t928;
t975 = -mrSges(6,1) * t996 + mrSges(6,2) * t997;
t914 = m(6) * t937 - mrSges(6,2) * t1010 + mrSges(6,3) * t959 - t1024 * t984 + t975 * t996 + t1090;
t936 = -t1063 * t942 + t1069 * t940;
t933 = -pkin(5) * t1010 - pkin(13) * t1023 + t976 * t997 - t936;
t1081 = -m(7) * t933 + t947 * mrSges(7,1) - mrSges(7,2) * t948 + t981 * t969 - t970 * t982;
t923 = m(6) * t936 + mrSges(6,1) * t1010 - mrSges(6,3) * t960 + t1024 * t983 - t975 * t997 + t1081;
t908 = t1063 * t914 + t1069 * t923;
t998 = -mrSges(5,1) * t1017 + mrSges(5,2) * t1018;
t906 = m(5) * t943 + mrSges(5,1) * t1011 - mrSges(5,3) * t980 + t1003 * t1026 - t1018 * t998 + t908;
t1089 = -t1063 * t923 + t1069 * t914;
t907 = m(5) * t944 - mrSges(5,2) * t1011 + mrSges(5,3) * t979 - t1004 * t1026 + t1017 * t998 + t1089;
t1088 = -t1064 * t906 + t1070 * t907;
t897 = m(4) * t967 - mrSges(4,2) * t1029 + mrSges(4,3) * t1012 + t1014 * t1027 - t1020 * t1038 + t1088;
t900 = t1064 * t907 + t1070 * t906;
t899 = m(4) * t977 - mrSges(4,1) * t1012 + mrSges(4,2) * t1013 - t1019 * t1027 + t1020 * t1028 + t900;
t886 = -t1058 * t899 + t1060 * t1104 + t897 * t1097;
t882 = m(3) * t1021 + mrSges(3,1) * t1054 - mrSges(3,3) * t1048 + t1044 * t1055 - t1047 * t1092 + t886;
t1033 = -t1059 * t1045 - t1106;
t1043 = mrSges(3,1) * t1055 - mrSges(3,3) * t1092;
t885 = t1058 * t1104 + t1060 * t899 + t897 * t1100;
t884 = m(3) * t1033 - t1049 * mrSges(3,1) + t1048 * mrSges(3,2) + (t1043 * t1066 - t1044 * t1072) * t1103 + t885;
t1022 = -g(3) * t1099 + t1093;
t893 = -t1065 * t911 + t1071 * t897;
t892 = m(3) * t1022 - mrSges(3,2) * t1054 + mrSges(3,3) * t1049 - t1043 * t1055 + t1047 * t1091 + t893;
t871 = -t1059 * t884 + t882 * t1095 + t892 * t1096;
t868 = m(2) * t1052 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1074 + t871;
t878 = -t1066 * t882 + t1072 * t892;
t876 = m(2) * t1053 - mrSges(2,1) * t1074 - qJDD(1) * mrSges(2,2) + t878;
t1105 = t1067 * t876 + t1073 * t868;
t870 = t1061 * t884 + t882 * t1098 + t892 * t1099;
t1087 = -t1067 * t868 + t1073 * t876;
t1006 = Ifges(4,5) * t1028 + Ifges(4,6) * t1027 + Ifges(4,3) * t1038;
t1007 = Ifges(4,4) * t1028 + Ifges(4,2) * t1027 + Ifges(4,6) * t1038;
t961 = Ifges(7,5) * t982 + Ifges(7,6) * t981 + Ifges(7,3) * t993;
t963 = Ifges(7,1) * t982 + Ifges(7,4) * t981 + Ifges(7,5) * t993;
t920 = -mrSges(7,1) * t933 + mrSges(7,3) * t932 + Ifges(7,4) * t948 + Ifges(7,2) * t947 + Ifges(7,6) * t958 - t961 * t982 + t963 * t993;
t962 = Ifges(7,4) * t982 + Ifges(7,2) * t981 + Ifges(7,6) * t993;
t921 = mrSges(7,2) * t933 - mrSges(7,3) * t931 + Ifges(7,1) * t948 + Ifges(7,4) * t947 + Ifges(7,5) * t958 + t961 * t981 - t962 * t993;
t971 = Ifges(6,5) * t997 + Ifges(6,6) * t996 + Ifges(6,3) * t1024;
t972 = Ifges(6,4) * t997 + Ifges(6,2) * t996 + Ifges(6,6) * t1024;
t901 = mrSges(6,2) * t945 - mrSges(6,3) * t936 + Ifges(6,1) * t960 + Ifges(6,4) * t959 + Ifges(6,5) * t1010 - pkin(13) * t916 - t1024 * t972 - t1062 * t920 + t1068 * t921 + t971 * t996;
t1077 = mrSges(7,1) * t931 - mrSges(7,2) * t932 + Ifges(7,5) * t948 + Ifges(7,6) * t947 + Ifges(7,3) * t958 + t962 * t982 - t963 * t981;
t973 = Ifges(6,1) * t997 + Ifges(6,4) * t996 + Ifges(6,5) * t1024;
t902 = -mrSges(6,1) * t945 + mrSges(6,3) * t937 + Ifges(6,4) * t960 + Ifges(6,2) * t959 + Ifges(6,6) * t1010 - pkin(5) * t916 + t1024 * t973 - t971 * t997 - t1077;
t985 = Ifges(5,5) * t1018 + Ifges(5,6) * t1017 + Ifges(5,3) * t1026;
t987 = Ifges(5,1) * t1018 + Ifges(5,4) * t1017 + Ifges(5,5) * t1026;
t887 = -mrSges(5,1) * t950 + mrSges(5,3) * t944 + Ifges(5,4) * t980 + Ifges(5,2) * t979 + Ifges(5,6) * t1011 - pkin(4) * t1079 + pkin(12) * t1089 - t1018 * t985 + t1026 * t987 + t1063 * t901 + t1069 * t902;
t986 = Ifges(5,4) * t1018 + Ifges(5,2) * t1017 + Ifges(5,6) * t1026;
t888 = mrSges(5,2) * t950 - mrSges(5,3) * t943 + Ifges(5,1) * t980 + Ifges(5,4) * t979 + Ifges(5,5) * t1011 - pkin(12) * t908 + t1017 * t985 - t1026 * t986 - t1063 * t902 + t1069 * t901;
t873 = mrSges(4,2) * t977 - mrSges(4,3) * t966 + Ifges(4,1) * t1013 + Ifges(4,4) * t1012 + Ifges(4,5) * t1029 - pkin(11) * t900 + t1006 * t1027 - t1007 * t1038 - t1064 * t887 + t1070 * t888;
t1008 = Ifges(4,1) * t1028 + Ifges(4,4) * t1027 + Ifges(4,5) * t1038;
t1078 = -mrSges(6,1) * t936 + mrSges(6,2) * t937 - Ifges(6,5) * t960 - Ifges(6,6) * t959 - Ifges(6,3) * t1010 - pkin(5) * t1081 - pkin(13) * t1090 - t1062 * t921 - t1068 * t920 - t997 * t972 + t996 * t973;
t1075 = mrSges(5,1) * t943 - mrSges(5,2) * t944 + Ifges(5,5) * t980 + Ifges(5,6) * t979 + Ifges(5,3) * t1011 + pkin(4) * t908 - t1017 * t987 + t1018 * t986 - t1078;
t879 = -mrSges(4,1) * t977 + mrSges(4,3) * t967 + Ifges(4,4) * t1013 + Ifges(4,2) * t1012 + Ifges(4,6) * t1029 - pkin(3) * t900 - t1028 * t1006 + t1038 * t1008 - t1075;
t1083 = pkin(10) * t893 + t1065 * t873 + t1071 * t879;
t1031 = Ifges(3,6) * t1055 + (Ifges(3,4) * t1066 + Ifges(3,2) * t1072) * t1103;
t1032 = Ifges(3,5) * t1055 + (Ifges(3,1) * t1066 + Ifges(3,4) * t1072) * t1103;
t872 = mrSges(4,1) * t966 - mrSges(4,2) * t967 + Ifges(4,5) * t1013 + Ifges(4,6) * t1012 + Ifges(4,3) * t1029 + pkin(3) * t1076 + pkin(11) * t1088 + t1028 * t1007 - t1027 * t1008 + t1064 * t888 + t1070 * t887;
t862 = mrSges(3,1) * t1021 - mrSges(3,2) * t1022 + Ifges(3,5) * t1048 + Ifges(3,6) * t1049 + Ifges(3,3) * t1054 + pkin(2) * t886 + t1060 * t872 + (t1031 * t1066 - t1032 * t1072) * t1103 + t1083 * t1058;
t1030 = Ifges(3,3) * t1055 + (Ifges(3,5) * t1066 + Ifges(3,6) * t1072) * t1103;
t864 = -mrSges(3,1) * t1033 + mrSges(3,3) * t1022 + Ifges(3,4) * t1048 + Ifges(3,2) * t1049 + Ifges(3,6) * t1054 - pkin(2) * t885 - t1030 * t1092 + t1055 * t1032 - t1058 * t872 + t1083 * t1060;
t866 = t1030 * t1091 + mrSges(3,2) * t1033 - mrSges(3,3) * t1021 + Ifges(3,1) * t1048 + Ifges(3,4) * t1049 + Ifges(3,5) * t1054 - t1055 * t1031 - t1065 * t879 + t1071 * t873 + (-t1058 * t885 - t1060 * t886) * pkin(10);
t1080 = mrSges(2,1) * t1052 - mrSges(2,2) * t1053 + Ifges(2,3) * qJDD(1) + pkin(1) * t871 + t1061 * t862 + t864 * t1098 + t866 * t1099 + t878 * t1109;
t860 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1052 + Ifges(2,5) * qJDD(1) - t1074 * Ifges(2,6) - t1066 * t864 + t1072 * t866 + (-t1059 * t870 - t1061 * t871) * pkin(9);
t859 = mrSges(2,1) * g(3) + mrSges(2,3) * t1053 + t1074 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t870 - t1059 * t862 + (pkin(9) * t878 + t1066 * t866 + t1072 * t864) * t1061;
t1 = [-m(1) * g(1) + t1087; -m(1) * g(2) + t1105; (-m(1) - m(2)) * g(3) + t870; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1105 - t1067 * t859 + t1073 * t860; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1087 + t1067 * t860 + t1073 * t859; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1080; t1080; t862; t872; t1075; -t1078; t1077;];
tauJB  = t1;
