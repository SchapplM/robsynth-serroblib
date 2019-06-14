% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 03:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR15_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 03:24:05
% EndTime: 2019-05-08 03:25:33
% DurationCPUTime: 64.97s
% Computational Cost: add. (1064077->401), mult. (2630440->511), div. (0->0), fcn. (2202642->14), ass. (0->177)
t1112 = Ifges(5,1) + Ifges(6,2);
t1106 = Ifges(5,4) + Ifges(6,6);
t1105 = Ifges(5,5) - Ifges(6,4);
t1111 = -Ifges(5,2) - Ifges(6,3);
t1104 = Ifges(5,6) - Ifges(6,5);
t1110 = Ifges(5,3) + Ifges(6,1);
t1048 = sin(pkin(7));
t1050 = cos(pkin(7));
t1054 = sin(qJ(3));
t1058 = cos(qJ(3));
t1051 = cos(pkin(6));
t1045 = qJD(1) * t1051 + qJD(2);
t1049 = sin(pkin(6));
t1059 = cos(qJ(2));
t1086 = t1049 * t1059;
t1080 = qJD(1) * t1086;
t1076 = t1050 * t1080;
t1027 = (t1045 * t1048 + t1076) * pkin(10);
t1055 = sin(qJ(2));
t1093 = qJD(1) * t1049;
t1101 = pkin(10) * t1048;
t1031 = (-pkin(2) * t1059 - t1055 * t1101) * t1093;
t1091 = qJD(1) * t1059;
t1037 = (qJD(2) * t1091 + qJDD(1) * t1055) * t1049;
t1044 = qJDD(1) * t1051 + qJDD(2);
t1056 = sin(qJ(1));
t1060 = cos(qJ(1));
t1042 = t1056 * g(1) - g(2) * t1060;
t1061 = qJD(1) ^ 2;
t1102 = pkin(9) * t1049;
t1034 = qJDD(1) * pkin(1) + t1061 * t1102 + t1042;
t1043 = -g(1) * t1060 - g(2) * t1056;
t1035 = -pkin(1) * t1061 + qJDD(1) * t1102 + t1043;
t1083 = t1051 * t1059;
t1077 = t1034 * t1083 - t1055 * t1035;
t1092 = qJD(1) * t1055;
t1100 = pkin(10) * t1050;
t978 = -t1037 * t1100 + t1044 * pkin(2) + t1045 * t1027 + (-g(3) * t1059 - t1031 * t1092) * t1049 + t1077;
t1087 = t1049 * t1055;
t1081 = qJD(1) * t1087;
t1030 = pkin(2) * t1045 - t1081 * t1100;
t1038 = (-qJD(2) * t1092 + qJDD(1) * t1059) * t1049;
t1073 = t1038 * t1050 + t1044 * t1048;
t1084 = t1051 * t1055;
t1082 = t1034 * t1084 + t1059 * t1035;
t979 = -t1045 * t1030 + (-g(3) * t1055 + t1031 * t1091) * t1049 + t1073 * pkin(10) + t1082;
t1099 = t1051 * g(3);
t985 = -t1037 * t1101 - t1038 * pkin(2) - t1099 + (-t1034 + (-t1027 * t1059 + t1030 * t1055) * qJD(1)) * t1049;
t946 = -t1054 * t979 + (t1048 * t985 + t1050 * t978) * t1058;
t1085 = t1050 * t1054;
t1089 = t1048 * t1054;
t1017 = t1045 * t1089 + (t1055 * t1058 + t1059 * t1085) * t1093;
t1000 = -t1017 * qJD(3) - t1054 * t1037 + t1058 * t1073;
t1028 = t1045 * t1050 - t1048 * t1080 + qJD(3);
t1053 = sin(qJ(4));
t1107 = cos(qJ(4));
t1007 = t1017 * t1053 - t1028 * t1107;
t1008 = t1017 * t1107 + t1053 * t1028;
t1052 = sin(qJ(6));
t1057 = cos(qJ(6));
t1088 = t1048 * t1058;
t1016 = t1045 * t1088 - t1054 * t1081 + t1058 * t1076;
t1014 = qJD(4) - t1016;
t1006 = t1007 ^ 2;
t1013 = t1014 ^ 2;
t1003 = -pkin(3) * t1016 - pkin(11) * t1017;
t1018 = -t1038 * t1048 + t1044 * t1050 + qJDD(3);
t1026 = t1028 ^ 2;
t947 = t1058 * t979 + t978 * t1085 + t985 * t1089;
t938 = -pkin(3) * t1026 + pkin(11) * t1018 + t1003 * t1016 + t947;
t1001 = t1016 * qJD(3) + t1058 * t1037 + t1054 * t1073;
t953 = -t1048 * t978 + t1050 * t985;
t940 = (-t1016 * t1028 - t1001) * pkin(11) + (t1017 * t1028 - t1000) * pkin(3) + t953;
t934 = t1053 * t940 + t1107 * t938;
t980 = pkin(4) * t1007 - qJ(5) * t1008;
t999 = qJDD(4) - t1000;
t1068 = -t1013 * pkin(4) + t999 * qJ(5) - t1007 * t980 + t934;
t960 = qJD(4) * t1008 + t1001 * t1053 - t1018 * t1107;
t993 = pkin(5) * t1008 - pkin(12) * t1014;
t928 = -t960 * pkin(5) - t1006 * pkin(12) + ((2 * qJD(5)) + t993) * t1014 + t1068;
t988 = t1007 * t1052 + t1014 * t1057;
t944 = -qJD(6) * t988 - t1052 * t999 + t1057 * t960;
t987 = t1007 * t1057 - t1014 * t1052;
t945 = qJD(6) * t987 + t1052 * t960 + t1057 * t999;
t1005 = qJD(6) + t1008;
t964 = -mrSges(7,2) * t1005 + mrSges(7,3) * t987;
t965 = mrSges(7,1) * t1005 - mrSges(7,3) * t988;
t1070 = -m(7) * t928 + mrSges(7,1) * t944 - t945 * mrSges(7,2) + t964 * t987 - t988 * t965;
t1108 = -2 * qJD(5);
t930 = t1014 * t1108 - t1068;
t990 = mrSges(6,1) * t1008 + mrSges(6,2) * t1014;
t1065 = -m(6) * t930 + t999 * mrSges(6,3) + t1014 * t990 - t1070;
t1094 = -t1106 * t1007 + t1008 * t1112 + t1105 * t1014;
t1095 = t1007 * t1111 + t1008 * t1106 + t1014 * t1104;
t1090 = t1007 * t1014;
t933 = -t1053 * t938 + t1107 * t940;
t931 = -t999 * pkin(4) - t1013 * qJ(5) + t1008 * t980 + qJDD(5) - t933;
t961 = -t1007 * qJD(4) + t1001 * t1107 + t1053 * t1018;
t926 = (t1007 * t1008 - t999) * pkin(12) + (t961 + t1090) * pkin(5) + t931;
t937 = -t1018 * pkin(3) - t1026 * pkin(11) + t1017 * t1003 - t946;
t1062 = (-t961 + t1090) * qJ(5) + t937 + (pkin(4) * t1014 + t1108) * t1008;
t929 = -pkin(5) * t1006 - t1008 * t993 + (pkin(4) + pkin(12)) * t960 + t1062;
t924 = -t1052 * t929 + t1057 * t926;
t954 = -mrSges(7,1) * t987 + mrSges(7,2) * t988;
t959 = qJDD(6) + t961;
t921 = m(7) * t924 + mrSges(7,1) * t959 - mrSges(7,3) * t945 + t1005 * t964 - t954 * t988;
t925 = t1052 * t926 + t1057 * t929;
t922 = m(7) * t925 - mrSges(7,2) * t959 + mrSges(7,3) * t944 - t1005 * t965 + t954 * t987;
t912 = t1052 * t922 + t1057 * t921;
t982 = -mrSges(6,2) * t1007 - mrSges(6,3) * t1008;
t1069 = -m(6) * t931 - t961 * mrSges(6,1) - t1008 * t982 - t912;
t989 = mrSges(6,1) * t1007 - mrSges(6,3) * t1014;
t910 = mrSges(6,2) * t999 + t1014 * t989 - t1069;
t949 = Ifges(7,5) * t988 + Ifges(7,6) * t987 + Ifges(7,3) * t1005;
t951 = Ifges(7,1) * t988 + Ifges(7,4) * t987 + Ifges(7,5) * t1005;
t913 = -mrSges(7,1) * t928 + mrSges(7,3) * t925 + Ifges(7,4) * t945 + Ifges(7,2) * t944 + Ifges(7,6) * t959 + t1005 * t951 - t949 * t988;
t950 = Ifges(7,4) * t988 + Ifges(7,2) * t987 + Ifges(7,6) * t1005;
t914 = mrSges(7,2) * t928 - mrSges(7,3) * t924 + Ifges(7,1) * t945 + Ifges(7,4) * t944 + Ifges(7,5) * t959 - t1005 * t950 + t949 * t987;
t1109 = t1007 * t1094 + t1008 * t1095 + t1110 * t999 - t1104 * t960 + t1105 * t961 + mrSges(5,1) * t933 - mrSges(5,2) * t934 + mrSges(6,2) * t931 - mrSges(6,3) * t930 - pkin(4) * t910 - pkin(12) * t912 + qJ(5) * (-mrSges(6,1) * t960 - t1007 * t982 + t1065) - t1052 * t913 + t1057 * t914;
t1011 = -g(3) * t1086 + t1077;
t1033 = -mrSges(3,2) * t1045 + mrSges(3,3) * t1080;
t1036 = (-mrSges(3,1) * t1059 + mrSges(3,2) * t1055) * t1093;
t1002 = -mrSges(4,1) * t1016 + mrSges(4,2) * t1017;
t1010 = mrSges(4,1) * t1028 - mrSges(4,3) * t1017;
t981 = mrSges(5,1) * t1007 + mrSges(5,2) * t1008;
t991 = -mrSges(5,2) * t1014 - mrSges(5,3) * t1007;
t909 = m(5) * t933 - mrSges(5,3) * t961 - t1008 * t981 + (mrSges(5,1) - mrSges(6,2)) * t999 + (-t989 + t991) * t1014 + t1069;
t992 = mrSges(5,1) * t1014 - mrSges(5,3) * t1008;
t917 = m(5) * t934 - mrSges(5,2) * t999 - t1014 * t992 + (-mrSges(5,3) - mrSges(6,1)) * t960 + (-t981 - t982) * t1007 + t1065;
t1079 = -t1053 * t909 + t1107 * t917;
t901 = m(4) * t947 - mrSges(4,2) * t1018 + mrSges(4,3) * t1000 + t1002 * t1016 - t1010 * t1028 + t1079;
t1009 = -mrSges(4,2) * t1028 + mrSges(4,3) * t1016;
t904 = t1053 * t917 + t1107 * t909;
t903 = m(4) * t953 - mrSges(4,1) * t1000 + mrSges(4,2) * t1001 - t1009 * t1016 + t1010 * t1017 + t904;
t1097 = -t1052 * t921 + t1057 * t922;
t932 = pkin(4) * t960 + t1062;
t1072 = -m(6) * t932 + t960 * mrSges(6,2) + t1007 * t989 - t1097;
t1064 = -m(5) * t937 - t960 * mrSges(5,1) - t1007 * t991 + (-mrSges(5,2) + mrSges(6,3)) * t961 + (t990 - t992) * t1008 + t1072;
t907 = m(4) * t946 + mrSges(4,1) * t1018 - mrSges(4,3) * t1001 - t1002 * t1017 + t1009 * t1028 + t1064;
t890 = t1050 * t1058 * t907 - t1048 * t903 + t901 * t1085;
t886 = m(3) * t1011 + mrSges(3,1) * t1044 - mrSges(3,3) * t1037 + t1033 * t1045 - t1036 * t1081 + t890;
t1022 = -t1049 * t1034 - t1099;
t1032 = mrSges(3,1) * t1045 - mrSges(3,3) * t1081;
t889 = t1050 * t903 + t907 * t1088 + t901 * t1089;
t888 = m(3) * t1022 - t1038 * mrSges(3,1) + t1037 * mrSges(3,2) + (t1032 * t1055 - t1033 * t1059) * t1093 + t889;
t1012 = -g(3) * t1087 + t1082;
t897 = -t1054 * t907 + t1058 * t901;
t896 = m(3) * t1012 - mrSges(3,2) * t1044 + mrSges(3,3) * t1038 - t1032 * t1045 + t1036 * t1080 + t897;
t875 = -t1049 * t888 + t886 * t1083 + t896 * t1084;
t872 = m(2) * t1042 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1061 + t875;
t882 = -t1055 * t886 + t1059 * t896;
t880 = m(2) * t1043 - mrSges(2,1) * t1061 - qJDD(1) * mrSges(2,2) + t882;
t1098 = t1056 * t880 + t1060 * t872;
t1096 = t1007 * t1104 - t1008 * t1105 - t1014 * t1110;
t874 = t1051 * t888 + t886 * t1086 + t896 * t1087;
t1078 = -t1056 * t872 + t1060 * t880;
t911 = -mrSges(6,3) * t961 - t1008 * t990 - t1072;
t891 = -mrSges(5,1) * t937 - mrSges(6,1) * t930 + mrSges(6,2) * t932 + mrSges(5,3) * t934 - pkin(4) * t911 - pkin(5) * t1070 - pkin(12) * t1097 + t1096 * t1008 + t1094 * t1014 - t1052 * t914 - t1057 * t913 + t1104 * t999 + t1106 * t961 + t1111 * t960;
t1066 = mrSges(7,1) * t924 - mrSges(7,2) * t925 + Ifges(7,5) * t945 + Ifges(7,6) * t944 + Ifges(7,3) * t959 + t988 * t950 - t987 * t951;
t892 = mrSges(6,1) * t931 + mrSges(5,2) * t937 - mrSges(5,3) * t933 - mrSges(6,3) * t932 + pkin(5) * t912 - qJ(5) * t911 + t1096 * t1007 - t1095 * t1014 + t1105 * t999 - t1106 * t960 + t1112 * t961 + t1066;
t995 = Ifges(4,5) * t1017 + Ifges(4,6) * t1016 + Ifges(4,3) * t1028;
t996 = Ifges(4,4) * t1017 + Ifges(4,2) * t1016 + Ifges(4,6) * t1028;
t877 = mrSges(4,2) * t953 - mrSges(4,3) * t946 + Ifges(4,1) * t1001 + Ifges(4,4) * t1000 + Ifges(4,5) * t1018 - pkin(11) * t904 + t1016 * t995 - t1028 * t996 - t1053 * t891 + t1107 * t892;
t997 = Ifges(4,1) * t1017 + Ifges(4,4) * t1016 + Ifges(4,5) * t1028;
t883 = -mrSges(4,1) * t953 + mrSges(4,3) * t947 + Ifges(4,4) * t1001 + Ifges(4,2) * t1000 + Ifges(4,6) * t1018 - pkin(3) * t904 - t1017 * t995 + t1028 * t997 - t1109;
t1071 = pkin(10) * t897 + t1054 * t877 + t1058 * t883;
t1020 = Ifges(3,6) * t1045 + (Ifges(3,4) * t1055 + Ifges(3,2) * t1059) * t1093;
t1021 = Ifges(3,5) * t1045 + (Ifges(3,1) * t1055 + Ifges(3,4) * t1059) * t1093;
t876 = mrSges(4,1) * t946 - mrSges(4,2) * t947 + Ifges(4,5) * t1001 + Ifges(4,6) * t1000 + Ifges(4,3) * t1018 + pkin(3) * t1064 + pkin(11) * t1079 - t1016 * t997 + t1017 * t996 + t1053 * t892 + t1107 * t891;
t866 = mrSges(3,1) * t1011 - mrSges(3,2) * t1012 + Ifges(3,5) * t1037 + Ifges(3,6) * t1038 + Ifges(3,3) * t1044 + pkin(2) * t890 + t1050 * t876 + (t1020 * t1055 - t1021 * t1059) * t1093 + t1071 * t1048;
t1019 = Ifges(3,3) * t1045 + (Ifges(3,5) * t1055 + Ifges(3,6) * t1059) * t1093;
t868 = -mrSges(3,1) * t1022 + mrSges(3,3) * t1012 + Ifges(3,4) * t1037 + Ifges(3,2) * t1038 + Ifges(3,6) * t1044 - pkin(2) * t889 - t1019 * t1081 + t1021 * t1045 - t1048 * t876 + t1050 * t1071;
t870 = t1019 * t1080 + mrSges(3,2) * t1022 - mrSges(3,3) * t1011 + Ifges(3,1) * t1037 + Ifges(3,4) * t1038 + Ifges(3,5) * t1044 - t1020 * t1045 - t1054 * t883 + t1058 * t877 + (-t1048 * t889 - t1050 * t890) * pkin(10);
t1067 = mrSges(2,1) * t1042 - mrSges(2,2) * t1043 + Ifges(2,3) * qJDD(1) + pkin(1) * t875 + t1051 * t866 + t868 * t1086 + t870 * t1087 + t882 * t1102;
t864 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1042 + Ifges(2,5) * qJDD(1) - t1061 * Ifges(2,6) - t1055 * t868 + t1059 * t870 + (-t1049 * t874 - t1051 * t875) * pkin(9);
t863 = mrSges(2,1) * g(3) + mrSges(2,3) * t1043 + t1061 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t874 - t1049 * t866 + (pkin(9) * t882 + t1055 * t870 + t1059 * t868) * t1051;
t1 = [-m(1) * g(1) + t1078; -m(1) * g(2) + t1098; (-m(1) - m(2)) * g(3) + t874; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1098 - t1056 * t863 + t1060 * t864; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1078 + t1056 * t864 + t1060 * t863; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1067; t1067; t866; t876; t1109; t910; t1066;];
tauJB  = t1;
