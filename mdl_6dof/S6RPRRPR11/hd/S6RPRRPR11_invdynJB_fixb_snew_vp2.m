% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-06 00:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:11:18
% EndTime: 2019-05-06 00:13:32
% DurationCPUTime: 138.92s
% Computational Cost: add. (2131083->404), mult. (6670359->553), div. (0->0), fcn. (5729246->16), ass. (0->184)
t1046 = sin(pkin(12));
t1048 = sin(pkin(6));
t1050 = cos(pkin(12));
t1052 = cos(pkin(6));
t1055 = sin(qJ(3));
t1051 = cos(pkin(7));
t1058 = cos(qJ(3));
t1090 = t1051 * t1058;
t1047 = sin(pkin(7));
t1095 = t1047 * t1058;
t1064 = t1048 * (-t1046 * t1055 + t1050 * t1090) + t1052 * t1095;
t1016 = t1064 * qJD(1);
t1091 = t1051 * t1055;
t1096 = t1047 * t1055;
t1066 = t1052 * t1096 + (t1046 * t1058 + t1050 * t1091) * t1048;
t1017 = t1066 * qJD(1);
t1005 = -t1017 * qJD(3) + qJDD(1) * t1064;
t1093 = t1048 * t1051;
t1030 = (t1047 * t1052 + t1050 * t1093) * qJD(1) * pkin(9);
t1056 = sin(qJ(1));
t1059 = cos(qJ(1));
t1042 = -g(1) * t1059 - g(2) * t1056;
t1060 = qJD(1) ^ 2;
t1101 = qJ(2) * t1048;
t1034 = -pkin(1) * t1060 + qJDD(1) * t1101 + t1042;
t1105 = pkin(9) * t1046;
t1077 = -pkin(2) * t1050 - t1047 * t1105;
t1100 = qJD(1) * t1048;
t1103 = pkin(9) * qJDD(1);
t1074 = qJD(1) * t1077 * t1100 + t1051 * t1103;
t1041 = t1056 * g(1) - g(2) * t1059;
t1033 = qJDD(1) * pkin(1) + t1060 * t1101 + t1041;
t1089 = qJD(2) * t1100;
t1092 = t1050 * t1052;
t1094 = t1048 * t1050;
t1078 = -g(3) * t1094 + t1033 * t1092 - 0.2e1 * t1046 * t1089;
t984 = (pkin(2) * qJDD(1) + qJD(1) * t1030) * t1052 + (-t1048 * t1074 - t1034) * t1046 + t1078;
t1035 = (pkin(2) * t1052 - t1093 * t1105) * qJD(1);
t1097 = t1046 * t1052;
t1087 = t1033 * t1097 + (t1034 + 0.2e1 * t1089) * t1050;
t985 = (-qJD(1) * t1035 + t1047 * t1103) * t1052 + (-g(3) * t1046 + t1050 * t1074) * t1048 + t1087;
t1086 = -t1052 * g(3) + qJDD(2);
t996 = (-t1033 + t1077 * qJDD(1) + (-t1030 * t1050 + t1035 * t1046) * qJD(1)) * t1048 + t1086;
t955 = -t1055 * t985 + (t1047 * t996 + t1051 * t984) * t1058;
t1106 = cos(qJ(4));
t1104 = Ifges(3,3) * t1052;
t1008 = -t1046 * t1034 + t1078;
t1081 = -mrSges(3,1) * t1050 + mrSges(3,2) * t1046;
t1032 = t1081 * t1100;
t1075 = -mrSges(3,2) * t1052 + mrSges(3,3) * t1094;
t1037 = t1075 * qJD(1);
t1098 = t1046 * t1048;
t1076 = mrSges(3,1) * t1052 - mrSges(3,3) * t1098;
t1003 = -mrSges(4,1) * t1016 + mrSges(4,2) * t1017;
t1072 = -t1047 * t1094 + t1051 * t1052;
t1031 = qJD(1) * t1072 + qJD(3);
t1013 = mrSges(4,1) * t1031 - mrSges(4,3) * t1017;
t1028 = qJDD(1) * t1072 + qJDD(3);
t1054 = sin(qJ(4));
t1002 = qJDD(4) - t1005;
t1010 = t1017 * t1054 - t1106 * t1031;
t1015 = qJD(4) - t1016;
t1045 = sin(pkin(13));
t1049 = cos(pkin(13));
t1011 = t1106 * t1017 + t1054 * t1031;
t994 = -t1011 * t1045 + t1015 * t1049;
t1082 = -mrSges(6,2) * t1010 + mrSges(6,3) * t994;
t1053 = sin(qJ(6));
t1057 = cos(qJ(6));
t1007 = qJD(6) + t1010;
t1014 = t1015 ^ 2;
t1004 = -pkin(3) * t1016 - pkin(10) * t1017;
t1027 = t1031 ^ 2;
t956 = t1058 * t985 + t984 * t1091 + t996 * t1096;
t947 = -pkin(3) * t1027 + pkin(10) * t1028 + t1004 * t1016 + t956;
t1006 = t1016 * qJD(3) + qJDD(1) * t1066;
t964 = -t1047 * t984 + t1051 * t996;
t953 = (-t1016 * t1031 - t1006) * pkin(10) + (t1017 * t1031 - t1005) * pkin(3) + t964;
t940 = t1054 * t953 + t1106 * t947;
t986 = pkin(4) * t1010 - qJ(5) * t1011;
t935 = -pkin(4) * t1014 + qJ(5) * t1002 - t1010 * t986 + t940;
t946 = -t1028 * pkin(3) - t1027 * pkin(10) + t1017 * t1004 - t955;
t979 = qJD(4) * t1011 + t1006 * t1054 - t1106 * t1028;
t980 = -t1010 * qJD(4) + t1106 * t1006 + t1054 * t1028;
t938 = (t1010 * t1015 - t980) * qJ(5) + (t1011 * t1015 + t979) * pkin(4) + t946;
t995 = t1011 * t1049 + t1015 * t1045;
t930 = -0.2e1 * qJD(5) * t995 - t1045 * t935 + t1049 * t938;
t966 = t1002 * t1045 + t1049 * t980;
t928 = (t1010 * t994 - t966) * pkin(11) + (t994 * t995 + t979) * pkin(5) + t930;
t931 = 0.2e1 * qJD(5) * t994 + t1045 * t938 + t1049 * t935;
t965 = t1002 * t1049 - t1045 * t980;
t972 = pkin(5) * t1010 - pkin(11) * t995;
t993 = t994 ^ 2;
t929 = -pkin(5) * t993 + pkin(11) * t965 - t1010 * t972 + t931;
t926 = -t1053 * t929 + t1057 * t928;
t967 = -t1053 * t995 + t1057 * t994;
t943 = qJD(6) * t967 + t1053 * t965 + t1057 * t966;
t968 = t1053 * t994 + t1057 * t995;
t954 = -mrSges(7,1) * t967 + mrSges(7,2) * t968;
t957 = -mrSges(7,2) * t1007 + mrSges(7,3) * t967;
t977 = qJDD(6) + t979;
t923 = m(7) * t926 + mrSges(7,1) * t977 - mrSges(7,3) * t943 + t1007 * t957 - t954 * t968;
t927 = t1053 * t928 + t1057 * t929;
t942 = -qJD(6) * t968 - t1053 * t966 + t1057 * t965;
t958 = mrSges(7,1) * t1007 - mrSges(7,3) * t968;
t924 = m(7) * t927 - mrSges(7,2) * t977 + mrSges(7,3) * t942 - t1007 * t958 + t954 * t967;
t915 = t1053 * t924 + t1057 * t923;
t969 = -mrSges(6,1) * t994 + mrSges(6,2) * t995;
t913 = m(6) * t930 + t979 * mrSges(6,1) - t966 * mrSges(6,3) + t1010 * t1082 - t995 * t969 + t915;
t1085 = -t1053 * t923 + t1057 * t924;
t971 = mrSges(6,1) * t1010 - mrSges(6,3) * t995;
t914 = m(6) * t931 - mrSges(6,2) * t979 + mrSges(6,3) * t965 - t1010 * t971 + t969 * t994 + t1085;
t911 = -t1045 * t913 + t1049 * t914;
t987 = mrSges(5,1) * t1010 + mrSges(5,2) * t1011;
t998 = mrSges(5,1) * t1015 - mrSges(5,3) * t1011;
t909 = m(5) * t940 - mrSges(5,2) * t1002 - mrSges(5,3) * t979 - t1010 * t987 - t1015 * t998 + t911;
t939 = -t1054 * t947 + t1106 * t953;
t934 = -t1002 * pkin(4) - t1014 * qJ(5) + t1011 * t986 + qJDD(5) - t939;
t932 = -t965 * pkin(5) - t993 * pkin(11) + t995 * t972 + t934;
t1070 = m(7) * t932 - t942 * mrSges(7,1) + mrSges(7,2) * t943 - t967 * t957 + t958 * t968;
t925 = m(6) * t934 - t965 * mrSges(6,1) + mrSges(6,2) * t966 - t994 * t1082 + t971 * t995 + t1070;
t997 = -mrSges(5,2) * t1015 - mrSges(5,3) * t1010;
t919 = m(5) * t939 + mrSges(5,1) * t1002 - mrSges(5,3) * t980 - t1011 * t987 + t1015 * t997 - t925;
t1084 = -t1054 * t919 + t1106 * t909;
t898 = m(4) * t956 - mrSges(4,2) * t1028 + mrSges(4,3) * t1005 + t1003 * t1016 - t1013 * t1031 + t1084;
t1012 = -mrSges(4,2) * t1031 + mrSges(4,3) * t1016;
t901 = t1054 * t909 + t1106 * t919;
t900 = m(4) * t964 - mrSges(4,1) * t1005 + mrSges(4,2) * t1006 - t1012 * t1016 + t1013 * t1017 + t901;
t910 = t1045 * t914 + t1049 * t913;
t1063 = -m(5) * t946 - t979 * mrSges(5,1) - mrSges(5,2) * t980 - t1010 * t997 - t1011 * t998 - t910;
t906 = m(4) * t955 + mrSges(4,1) * t1028 - mrSges(4,3) * t1006 - t1003 * t1017 + t1012 * t1031 + t1063;
t887 = -t1047 * t900 + t906 * t1090 + t898 * t1091;
t883 = m(3) * t1008 + t1076 * qJDD(1) + (-t1032 * t1098 + t1037 * t1052) * qJD(1) + t887;
t1018 = -t1048 * t1033 + t1086;
t1036 = t1076 * qJD(1);
t886 = t1051 * t900 + t906 * t1095 + t898 * t1096;
t885 = m(3) * t1018 + (t1081 * qJDD(1) + (t1036 * t1046 - t1037 * t1050) * qJD(1)) * t1048 + t886;
t1009 = -g(3) * t1098 + t1087;
t893 = -t1055 * t906 + t1058 * t898;
t892 = m(3) * t1009 + t1075 * qJDD(1) + (t1032 * t1094 - t1036 * t1052) * qJD(1) + t893;
t872 = -t1048 * t885 + t883 * t1092 + t892 * t1097;
t869 = m(2) * t1041 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1060 + t872;
t879 = -t1046 * t883 + t1050 * t892;
t877 = m(2) * t1042 - mrSges(2,1) * t1060 - qJDD(1) * mrSges(2,2) + t879;
t1102 = t1056 * t877 + t1059 * t869;
t871 = t1052 * t885 + t883 * t1094 + t892 * t1098;
t1083 = -t1056 * t869 + t1059 * t877;
t1080 = Ifges(3,5) * t1046 + Ifges(3,6) * t1050;
t1000 = Ifges(4,4) * t1017 + Ifges(4,2) * t1016 + Ifges(4,6) * t1031;
t949 = Ifges(7,5) * t968 + Ifges(7,6) * t967 + Ifges(7,3) * t1007;
t951 = Ifges(7,1) * t968 + Ifges(7,4) * t967 + Ifges(7,5) * t1007;
t916 = -mrSges(7,1) * t932 + mrSges(7,3) * t927 + Ifges(7,4) * t943 + Ifges(7,2) * t942 + Ifges(7,6) * t977 + t1007 * t951 - t949 * t968;
t950 = Ifges(7,4) * t968 + Ifges(7,2) * t967 + Ifges(7,6) * t1007;
t917 = mrSges(7,2) * t932 - mrSges(7,3) * t926 + Ifges(7,1) * t943 + Ifges(7,4) * t942 + Ifges(7,5) * t977 - t1007 * t950 + t949 * t967;
t959 = Ifges(6,5) * t995 + Ifges(6,6) * t994 + Ifges(6,3) * t1010;
t961 = Ifges(6,1) * t995 + Ifges(6,4) * t994 + Ifges(6,5) * t1010;
t902 = -mrSges(6,1) * t934 + mrSges(6,3) * t931 + Ifges(6,4) * t966 + Ifges(6,2) * t965 + Ifges(6,6) * t979 - pkin(5) * t1070 + pkin(11) * t1085 + t1010 * t961 + t1053 * t917 + t1057 * t916 - t995 * t959;
t960 = Ifges(6,4) * t995 + Ifges(6,2) * t994 + Ifges(6,6) * t1010;
t903 = mrSges(6,2) * t934 - mrSges(6,3) * t930 + Ifges(6,1) * t966 + Ifges(6,4) * t965 + Ifges(6,5) * t979 - pkin(11) * t915 - t1010 * t960 - t1053 * t916 + t1057 * t917 + t959 * t994;
t973 = Ifges(5,5) * t1011 - Ifges(5,6) * t1010 + Ifges(5,3) * t1015;
t974 = Ifges(5,4) * t1011 - Ifges(5,2) * t1010 + Ifges(5,6) * t1015;
t888 = mrSges(5,2) * t946 - mrSges(5,3) * t939 + Ifges(5,1) * t980 - Ifges(5,4) * t979 + Ifges(5,5) * t1002 - qJ(5) * t910 - t1010 * t973 - t1015 * t974 - t1045 * t902 + t1049 * t903;
t1062 = mrSges(7,1) * t926 - mrSges(7,2) * t927 + Ifges(7,5) * t943 + Ifges(7,6) * t942 + Ifges(7,3) * t977 + t968 * t950 - t967 * t951;
t975 = Ifges(5,1) * t1011 - Ifges(5,4) * t1010 + Ifges(5,5) * t1015;
t894 = -pkin(4) * t910 + (-Ifges(5,2) - Ifges(6,3)) * t979 - t1011 * t973 + t1015 * t975 + Ifges(5,6) * t1002 + t994 * t961 - t995 * t960 + Ifges(5,4) * t980 - Ifges(6,6) * t965 - Ifges(6,5) * t966 - mrSges(5,1) * t946 + mrSges(5,3) * t940 + mrSges(6,2) * t931 - mrSges(6,1) * t930 - t1062 - pkin(5) * t915;
t999 = Ifges(4,5) * t1017 + Ifges(4,6) * t1016 + Ifges(4,3) * t1031;
t874 = mrSges(4,2) * t964 - mrSges(4,3) * t955 + Ifges(4,1) * t1006 + Ifges(4,4) * t1005 + Ifges(4,5) * t1028 - pkin(10) * t901 - t1031 * t1000 + t1016 * t999 - t1054 * t894 + t1106 * t888;
t1001 = Ifges(4,1) * t1017 + Ifges(4,4) * t1016 + Ifges(4,5) * t1031;
t1061 = mrSges(5,1) * t939 - mrSges(5,2) * t940 + Ifges(5,5) * t980 - Ifges(5,6) * t979 + Ifges(5,3) * t1002 - pkin(4) * t925 + qJ(5) * t911 + t1010 * t975 + t1011 * t974 + t1045 * t903 + t1049 * t902;
t880 = -mrSges(4,1) * t964 + mrSges(4,3) * t956 + Ifges(4,4) * t1006 + Ifges(4,2) * t1005 + Ifges(4,6) * t1028 - pkin(3) * t901 + t1031 * t1001 - t1017 * t999 - t1061;
t1071 = pkin(9) * t893 + t1055 * t874 + t1058 * t880;
t1067 = Ifges(3,6) * t1052 + (Ifges(3,4) * t1046 + Ifges(3,2) * t1050) * t1048;
t1022 = t1067 * qJD(1);
t1068 = Ifges(3,5) * t1052 + (Ifges(3,1) * t1046 + Ifges(3,4) * t1050) * t1048;
t1023 = t1068 * qJD(1);
t873 = mrSges(4,1) * t955 - mrSges(4,2) * t956 + Ifges(4,5) * t1006 + Ifges(4,6) * t1005 + Ifges(4,3) * t1028 + pkin(3) * t1063 + pkin(10) * t1084 + t1017 * t1000 - t1016 * t1001 + t1054 * t888 + t1106 * t894;
t863 = qJDD(1) * t1104 + mrSges(3,1) * t1008 - mrSges(3,2) * t1009 + pkin(2) * t887 + t1051 * t873 + t1071 * t1047 + (t1080 * qJDD(1) + (t1022 * t1046 - t1023 * t1050) * qJD(1)) * t1048;
t1021 = (t1048 * t1080 + t1104) * qJD(1);
t865 = -mrSges(3,1) * t1018 + mrSges(3,3) * t1009 - pkin(2) * t886 - t1047 * t873 + (-t1021 * t1098 + t1023 * t1052) * qJD(1) + t1071 * t1051 + t1067 * qJDD(1);
t867 = mrSges(3,2) * t1018 - mrSges(3,3) * t1008 - t1055 * t880 + t1058 * t874 + (t1021 * t1094 - t1022 * t1052) * qJD(1) + (-t1047 * t886 - t1051 * t887) * pkin(9) + t1068 * qJDD(1);
t1069 = mrSges(2,1) * t1041 - mrSges(2,2) * t1042 + Ifges(2,3) * qJDD(1) + pkin(1) * t872 + t1052 * t863 + t865 * t1094 + t867 * t1098 + t879 * t1101;
t861 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1041 + Ifges(2,5) * qJDD(1) - t1060 * Ifges(2,6) - t1046 * t865 + t1050 * t867 + (-t1048 * t871 - t1052 * t872) * qJ(2);
t860 = mrSges(2,1) * g(3) + mrSges(2,3) * t1042 + t1060 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t871 - t1048 * t863 + (qJ(2) * t879 + t1046 * t867 + t1050 * t865) * t1052;
t1 = [-m(1) * g(1) + t1083; -m(1) * g(2) + t1102; (-m(1) - m(2)) * g(3) + t871; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1102 - t1056 * t860 + t1059 * t861; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1083 + t1056 * t861 + t1059 * t860; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1069; t1069; t885; t873; t1061; t925; t1062;];
tauJB  = t1;
