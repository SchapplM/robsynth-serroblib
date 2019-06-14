% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:04:41
% EndTime: 2019-05-06 02:05:45
% DurationCPUTime: 65.69s
% Computational Cost: add. (989496->383), mult. (3073257->510), div. (0->0), fcn. (2617472->14), ass. (0->178)
t1040 = sin(pkin(12));
t1042 = sin(pkin(6));
t1043 = cos(pkin(12));
t1045 = cos(pkin(6));
t1048 = sin(qJ(3));
t1044 = cos(pkin(7));
t1052 = cos(qJ(3));
t1084 = t1044 * t1052;
t1041 = sin(pkin(7));
t1089 = t1041 * t1052;
t1057 = t1042 * (-t1040 * t1048 + t1043 * t1084) + t1045 * t1089;
t1013 = t1057 * qJD(1);
t1085 = t1044 * t1048;
t1090 = t1041 * t1048;
t1059 = t1045 * t1090 + (t1040 * t1052 + t1043 * t1085) * t1042;
t1014 = t1059 * qJD(1);
t1002 = -t1014 * qJD(3) + qJDD(1) * t1057;
t1112 = Ifges(6,1) + Ifges(7,1);
t1107 = Ifges(6,4) + Ifges(7,4);
t1106 = Ifges(6,5) + Ifges(7,5);
t1111 = Ifges(6,2) + Ifges(7,2);
t1105 = Ifges(6,6) + Ifges(7,6);
t1110 = Ifges(6,3) + Ifges(7,3);
t1087 = t1042 * t1044;
t1025 = (t1041 * t1045 + t1043 * t1087) * qJD(1) * pkin(9);
t1049 = sin(qJ(1));
t1053 = cos(qJ(1));
t1037 = -g(1) * t1053 - g(2) * t1049;
t1054 = qJD(1) ^ 2;
t1095 = qJ(2) * t1042;
t1029 = -pkin(1) * t1054 + qJDD(1) * t1095 + t1037;
t1103 = pkin(9) * t1040;
t1070 = -pkin(2) * t1043 - t1041 * t1103;
t1094 = qJD(1) * t1042;
t1101 = pkin(9) * qJDD(1);
t1067 = qJD(1) * t1070 * t1094 + t1044 * t1101;
t1036 = t1049 * g(1) - g(2) * t1053;
t1028 = qJDD(1) * pkin(1) + t1054 * t1095 + t1036;
t1081 = qJD(2) * t1094;
t1086 = t1043 * t1045;
t1088 = t1042 * t1043;
t1071 = -g(3) * t1088 + t1028 * t1086 - 0.2e1 * t1040 * t1081;
t983 = (pkin(2) * qJDD(1) + qJD(1) * t1025) * t1045 + (-t1042 * t1067 - t1029) * t1040 + t1071;
t1030 = (pkin(2) * t1045 - t1087 * t1103) * qJD(1);
t1091 = t1040 * t1045;
t1079 = t1028 * t1091 + (t1029 + 0.2e1 * t1081) * t1043;
t984 = (-qJD(1) * t1030 + t1041 * t1101) * t1045 + (-g(3) * t1040 + t1043 * t1067) * t1042 + t1079;
t1078 = -t1045 * g(3) + qJDD(2);
t993 = (-t1028 + t1070 * qJDD(1) + (-t1025 * t1043 + t1030 * t1040) * qJD(1)) * t1042 + t1078;
t945 = -t1048 * t984 + (t1041 * t993 + t1044 * t983) * t1052;
t1065 = -t1041 * t1088 + t1044 * t1045;
t1026 = qJD(1) * t1065 + qJD(3);
t1047 = sin(qJ(4));
t1051 = cos(qJ(4));
t1007 = -t1014 * t1047 + t1026 * t1051;
t1004 = qJD(5) - t1007;
t1008 = t1014 * t1051 + t1026 * t1047;
t1012 = qJD(4) - t1013;
t1046 = sin(qJ(5));
t1050 = cos(qJ(5));
t991 = -t1008 * t1046 + t1012 * t1050;
t992 = t1008 * t1050 + t1012 * t1046;
t1097 = t1106 * t1004 + t1107 * t991 + t1112 * t992;
t1098 = -t1004 * t1105 - t1107 * t992 - t1111 * t991;
t1011 = t1012 ^ 2;
t1001 = -pkin(3) * t1013 - pkin(10) * t1014;
t1022 = t1026 ^ 2;
t1023 = qJDD(1) * t1065 + qJDD(3);
t946 = t1052 * t984 + t983 * t1085 + t993 * t1090;
t942 = -pkin(3) * t1022 + pkin(10) * t1023 + t1001 * t1013 + t946;
t1003 = t1013 * qJD(3) + qJDD(1) * t1059;
t961 = -t1041 * t983 + t1044 * t993;
t944 = (-t1013 * t1026 - t1003) * pkin(10) + (t1014 * t1026 - t1002) * pkin(3) + t961;
t938 = t1047 * t944 + t1051 * t942;
t986 = -pkin(4) * t1007 - pkin(11) * t1008;
t999 = qJDD(4) - t1002;
t933 = -pkin(4) * t1011 + pkin(11) * t999 + t1007 * t986 + t938;
t941 = -t1023 * pkin(3) - t1022 * pkin(10) + t1014 * t1001 - t945;
t978 = -qJD(4) * t1008 - t1003 * t1047 + t1023 * t1051;
t979 = qJD(4) * t1007 + t1003 * t1051 + t1023 * t1047;
t936 = (-t1007 * t1012 - t979) * pkin(11) + (t1008 * t1012 - t978) * pkin(4) + t941;
t928 = -t1046 * t933 + t1050 * t936;
t951 = qJD(5) * t991 + t1046 * t999 + t1050 * t979;
t976 = qJDD(5) - t978;
t924 = -0.2e1 * qJD(6) * t992 + (t1004 * t991 - t951) * qJ(6) + (t991 * t992 + t976) * pkin(5) + t928;
t966 = -mrSges(7,2) * t1004 + mrSges(7,3) * t991;
t1083 = m(7) * t924 + t976 * mrSges(7,1) + t1004 * t966;
t963 = -mrSges(7,1) * t991 + mrSges(7,2) * t992;
t922 = -t951 * mrSges(7,3) - t992 * t963 + t1083;
t929 = t1046 * t936 + t1050 * t933;
t950 = -qJD(5) * t992 - t1046 * t979 + t1050 * t999;
t968 = pkin(5) * t1004 - qJ(6) * t992;
t990 = t991 ^ 2;
t927 = -pkin(5) * t990 + qJ(6) * t950 + 0.2e1 * qJD(6) * t991 - t1004 * t968 + t929;
t1109 = mrSges(6,1) * t928 + mrSges(7,1) * t924 - mrSges(6,2) * t929 - mrSges(7,2) * t927 + pkin(5) * t922 - t1097 * t991 - t1098 * t992 + t1105 * t950 + t1106 * t951 + t1110 * t976;
t1108 = -mrSges(6,2) - mrSges(7,2);
t1102 = Ifges(3,3) * t1045;
t1005 = -t1040 * t1029 + t1071;
t1074 = -mrSges(3,1) * t1043 + mrSges(3,2) * t1040;
t1027 = t1074 * t1094;
t1068 = -mrSges(3,2) * t1045 + mrSges(3,3) * t1088;
t1032 = t1068 * qJD(1);
t1092 = t1040 * t1042;
t1069 = mrSges(3,1) * t1045 - mrSges(3,3) * t1092;
t1000 = -mrSges(4,1) * t1013 + mrSges(4,2) * t1014;
t1010 = mrSges(4,1) * t1026 - mrSges(4,3) * t1014;
t964 = -mrSges(6,1) * t991 + mrSges(6,2) * t992;
t967 = -mrSges(6,2) * t1004 + mrSges(6,3) * t991;
t916 = m(6) * t928 + t976 * mrSges(6,1) + t1004 * t967 + (-t963 - t964) * t992 + (-mrSges(6,3) - mrSges(7,3)) * t951 + t1083;
t1082 = m(7) * t927 + t950 * mrSges(7,3) + t991 * t963;
t969 = mrSges(7,1) * t1004 - mrSges(7,3) * t992;
t1096 = -mrSges(6,1) * t1004 + mrSges(6,3) * t992 - t969;
t918 = m(6) * t929 + t950 * mrSges(6,3) + t1096 * t1004 + t1108 * t976 + t991 * t964 + t1082;
t915 = -t1046 * t916 + t1050 * t918;
t985 = -mrSges(5,1) * t1007 + mrSges(5,2) * t1008;
t995 = mrSges(5,1) * t1012 - mrSges(5,3) * t1008;
t912 = m(5) * t938 - mrSges(5,2) * t999 + mrSges(5,3) * t978 + t1007 * t985 - t1012 * t995 + t915;
t937 = -t1047 * t942 + t1051 * t944;
t932 = -pkin(4) * t999 - pkin(11) * t1011 + t1008 * t986 - t937;
t930 = -pkin(5) * t950 - qJ(6) * t990 + t968 * t992 + qJDD(6) + t932;
t1075 = -m(7) * t930 + t950 * mrSges(7,1) + t991 * t966;
t921 = -m(6) * t932 + t950 * mrSges(6,1) + t1096 * t992 + t1108 * t951 + t991 * t967 + t1075;
t994 = -mrSges(5,2) * t1012 + mrSges(5,3) * t1007;
t920 = m(5) * t937 + t999 * mrSges(5,1) - t979 * mrSges(5,3) - t1008 * t985 + t1012 * t994 + t921;
t1077 = -t1047 * t920 + t1051 * t912;
t902 = m(4) * t946 - mrSges(4,2) * t1023 + mrSges(4,3) * t1002 + t1000 * t1013 - t1010 * t1026 + t1077;
t1009 = -mrSges(4,2) * t1026 + mrSges(4,3) * t1013;
t905 = t1047 * t912 + t1051 * t920;
t904 = m(4) * t961 - mrSges(4,1) * t1002 + mrSges(4,2) * t1003 - t1009 * t1013 + t1010 * t1014 + t905;
t914 = t1046 * t918 + t1050 * t916;
t1056 = -m(5) * t941 + t978 * mrSges(5,1) - mrSges(5,2) * t979 + t1007 * t994 - t1008 * t995 - t914;
t909 = m(4) * t945 + mrSges(4,1) * t1023 - mrSges(4,3) * t1003 - t1000 * t1014 + t1009 * t1026 + t1056;
t891 = -t1041 * t904 + t909 * t1084 + t902 * t1085;
t887 = m(3) * t1005 + t1069 * qJDD(1) + (-t1027 * t1092 + t1032 * t1045) * qJD(1) + t891;
t1015 = -t1042 * t1028 + t1078;
t1031 = t1069 * qJD(1);
t890 = t1044 * t904 + t909 * t1089 + t902 * t1090;
t889 = m(3) * t1015 + (t1074 * qJDD(1) + (t1031 * t1040 - t1032 * t1043) * qJD(1)) * t1042 + t890;
t1006 = -g(3) * t1092 + t1079;
t897 = -t1048 * t909 + t1052 * t902;
t896 = m(3) * t1006 + t1068 * qJDD(1) + (t1027 * t1088 - t1031 * t1045) * qJD(1) + t897;
t876 = -t1042 * t889 + t887 * t1086 + t896 * t1091;
t873 = m(2) * t1036 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1054 + t876;
t883 = -t1040 * t887 + t1043 * t896;
t881 = m(2) * t1037 - mrSges(2,1) * t1054 - qJDD(1) * mrSges(2,2) + t883;
t1100 = t1049 * t881 + t1053 * t873;
t1099 = -t1004 * t1110 - t1105 * t991 - t1106 * t992;
t875 = t1045 * t889 + t887 * t1088 + t896 * t1092;
t1076 = -t1049 * t873 + t1053 * t881;
t1073 = Ifges(3,5) * t1040 + Ifges(3,6) * t1043;
t925 = t951 * mrSges(7,2) + t992 * t969 - t1075;
t906 = -mrSges(6,1) * t932 + mrSges(6,3) * t929 - mrSges(7,1) * t930 + mrSges(7,3) * t927 - pkin(5) * t925 + qJ(6) * t1082 + t1099 * t992 + (-mrSges(7,2) * qJ(6) + t1105) * t976 + t1107 * t951 + t1111 * t950 + (-qJ(6) * t969 + t1097) * t1004;
t913 = mrSges(6,2) * t932 + mrSges(7,2) * t930 - mrSges(6,3) * t928 - mrSges(7,3) * t924 - qJ(6) * t922 + t1098 * t1004 - t1099 * t991 + t1106 * t976 + t1107 * t950 + t1112 * t951;
t972 = Ifges(5,5) * t1008 + Ifges(5,6) * t1007 + Ifges(5,3) * t1012;
t973 = Ifges(5,4) * t1008 + Ifges(5,2) * t1007 + Ifges(5,6) * t1012;
t892 = mrSges(5,2) * t941 - mrSges(5,3) * t937 + Ifges(5,1) * t979 + Ifges(5,4) * t978 + Ifges(5,5) * t999 - pkin(11) * t914 + t1007 * t972 - t1012 * t973 - t1046 * t906 + t1050 * t913;
t974 = Ifges(5,1) * t1008 + Ifges(5,4) * t1007 + Ifges(5,5) * t1012;
t898 = -mrSges(5,1) * t941 + mrSges(5,3) * t938 + Ifges(5,4) * t979 + Ifges(5,2) * t978 + Ifges(5,6) * t999 - pkin(4) * t914 - t1008 * t972 + t1012 * t974 - t1109;
t996 = Ifges(4,5) * t1014 + Ifges(4,6) * t1013 + Ifges(4,3) * t1026;
t997 = Ifges(4,4) * t1014 + Ifges(4,2) * t1013 + Ifges(4,6) * t1026;
t878 = mrSges(4,2) * t961 - mrSges(4,3) * t945 + Ifges(4,1) * t1003 + Ifges(4,4) * t1002 + Ifges(4,5) * t1023 - pkin(10) * t905 + t1013 * t996 - t1026 * t997 - t1047 * t898 + t1051 * t892;
t1055 = mrSges(5,1) * t937 - mrSges(5,2) * t938 + Ifges(5,5) * t979 + Ifges(5,6) * t978 + Ifges(5,3) * t999 + pkin(4) * t921 + pkin(11) * t915 - t1007 * t974 + t1008 * t973 + t1046 * t913 + t1050 * t906;
t998 = Ifges(4,1) * t1014 + Ifges(4,4) * t1013 + Ifges(4,5) * t1026;
t884 = -mrSges(4,1) * t961 + mrSges(4,3) * t946 + Ifges(4,4) * t1003 + Ifges(4,2) * t1002 + Ifges(4,6) * t1023 - pkin(3) * t905 - t1014 * t996 + t1026 * t998 - t1055;
t1064 = pkin(9) * t897 + t1048 * t878 + t1052 * t884;
t1061 = Ifges(3,6) * t1045 + (Ifges(3,4) * t1040 + Ifges(3,2) * t1043) * t1042;
t1019 = t1061 * qJD(1);
t1062 = Ifges(3,5) * t1045 + (Ifges(3,1) * t1040 + Ifges(3,4) * t1043) * t1042;
t1020 = t1062 * qJD(1);
t877 = mrSges(4,1) * t945 - mrSges(4,2) * t946 + Ifges(4,5) * t1003 + Ifges(4,6) * t1002 + Ifges(4,3) * t1023 + pkin(3) * t1056 + pkin(10) * t1077 - t1013 * t998 + t1014 * t997 + t1047 * t892 + t1051 * t898;
t867 = qJDD(1) * t1102 + mrSges(3,1) * t1005 - mrSges(3,2) * t1006 + pkin(2) * t891 + t1044 * t877 + t1064 * t1041 + (t1073 * qJDD(1) + (t1019 * t1040 - t1020 * t1043) * qJD(1)) * t1042;
t1018 = (t1042 * t1073 + t1102) * qJD(1);
t869 = -mrSges(3,1) * t1015 + mrSges(3,3) * t1006 - pkin(2) * t890 - t1041 * t877 + (-t1018 * t1092 + t1020 * t1045) * qJD(1) + t1064 * t1044 + t1061 * qJDD(1);
t871 = mrSges(3,2) * t1015 - mrSges(3,3) * t1005 - t1048 * t884 + t1052 * t878 + (t1018 * t1088 - t1019 * t1045) * qJD(1) + (-t1041 * t890 - t1044 * t891) * pkin(9) + t1062 * qJDD(1);
t1063 = mrSges(2,1) * t1036 - mrSges(2,2) * t1037 + Ifges(2,3) * qJDD(1) + pkin(1) * t876 + t1045 * t867 + t869 * t1088 + t871 * t1092 + t883 * t1095;
t865 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1036 + Ifges(2,5) * qJDD(1) - t1054 * Ifges(2,6) - t1040 * t869 + t1043 * t871 + (-t1042 * t875 - t1045 * t876) * qJ(2);
t864 = mrSges(2,1) * g(3) + mrSges(2,3) * t1037 + t1054 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t875 - t1042 * t867 + (qJ(2) * t883 + t1040 * t871 + t1043 * t869) * t1045;
t1 = [-m(1) * g(1) + t1076; -m(1) * g(2) + t1100; (-m(1) - m(2)) * g(3) + t875; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1100 - t1049 * t864 + t1053 * t865; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1076 + t1049 * t865 + t1053 * t864; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1063; t1063; t889; t877; t1055; t1109; t925;];
tauJB  = t1;
