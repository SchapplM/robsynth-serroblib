% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR9
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
% Datum: 2019-05-05 19:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:34:40
% EndTime: 2019-05-05 19:36:43
% DurationCPUTime: 126.48s
% Computational Cost: add. (1856984->404), mult. (6190045->552), div. (0->0), fcn. (5325550->16), ass. (0->184)
t1038 = sin(pkin(12));
t1040 = sin(pkin(6));
t1042 = cos(pkin(12));
t1044 = cos(pkin(6));
t1047 = sin(qJ(3));
t1043 = cos(pkin(7));
t1051 = cos(qJ(3));
t1078 = t1043 * t1051;
t1039 = sin(pkin(7));
t1083 = t1039 * t1051;
t1057 = t1044 * t1083 + (-t1038 * t1047 + t1042 * t1078) * t1040;
t1010 = t1057 * qJD(1);
t1079 = t1043 * t1047;
t1084 = t1039 * t1047;
t1058 = t1044 * t1084 + (t1038 * t1051 + t1042 * t1079) * t1040;
t1011 = t1058 * qJD(1);
t1037 = sin(pkin(13));
t1041 = cos(pkin(13));
t1000 = t1010 * t1037 + t1011 * t1041;
t1003 = t1010 * qJD(3) + t1058 * qJDD(1);
t1082 = t1040 * t1042;
t1063 = -t1039 * t1082 + t1043 * t1044;
t1020 = t1063 * qJDD(1) + qJDD(3);
t1023 = t1063 * qJD(1) + qJD(3);
t1081 = t1040 * t1043;
t1022 = (t1039 * t1044 + t1042 * t1081) * qJD(1) * pkin(9);
t1048 = sin(qJ(1));
t1052 = cos(qJ(1));
t1034 = -g(1) * t1052 - g(2) * t1048;
t1053 = qJD(1) ^ 2;
t1089 = qJ(2) * t1040;
t1026 = -pkin(1) * t1053 + qJDD(1) * t1089 + t1034;
t1093 = pkin(9) * t1038;
t1067 = -pkin(2) * t1042 - t1039 * t1093;
t1088 = qJD(1) * t1040;
t1091 = pkin(9) * qJDD(1);
t1064 = qJD(1) * t1067 * t1088 + t1043 * t1091;
t1033 = t1048 * g(1) - g(2) * t1052;
t1025 = qJDD(1) * pkin(1) + t1053 * t1089 + t1033;
t1077 = qJD(2) * t1088;
t1080 = t1042 * t1044;
t1068 = -g(3) * t1082 + t1025 * t1080 - 0.2e1 * t1038 * t1077;
t980 = (pkin(2) * qJDD(1) + qJD(1) * t1022) * t1044 + (-t1064 * t1040 - t1026) * t1038 + t1068;
t1027 = (pkin(2) * t1044 - t1081 * t1093) * qJD(1);
t1085 = t1038 * t1044;
t1076 = t1025 * t1085 + (t1026 + 0.2e1 * t1077) * t1042;
t981 = (-qJD(1) * t1027 + t1039 * t1091) * t1044 + (-g(3) * t1038 + t1064 * t1042) * t1040 + t1076;
t1075 = -t1044 * g(3) + qJDD(2);
t988 = (-t1025 + t1067 * qJDD(1) + (-t1022 * t1042 + t1027 * t1038) * qJD(1)) * t1040 + t1075;
t946 = -t1047 * t981 + t980 * t1078 + t988 * t1083;
t936 = (t1010 * t1023 - t1003) * qJ(4) + (t1010 * t1011 + t1020) * pkin(3) + t946;
t1002 = -t1011 * qJD(3) + t1057 * qJDD(1);
t1007 = pkin(3) * t1023 - qJ(4) * t1011;
t1009 = t1010 ^ 2;
t947 = t1051 * t981 + t980 * t1079 + t988 * t1084;
t939 = -pkin(3) * t1009 + qJ(4) * t1002 - t1007 * t1023 + t947;
t928 = -0.2e1 * qJD(4) * t1000 - t1037 * t939 + t1041 * t936;
t1092 = Ifges(3,3) * t1044;
t1004 = -t1038 * t1026 + t1068;
t1070 = -mrSges(3,1) * t1042 + mrSges(3,2) * t1038;
t1024 = t1070 * t1088;
t1065 = -mrSges(3,2) * t1044 + mrSges(3,3) * t1082;
t1029 = t1065 * qJD(1);
t1086 = t1038 * t1040;
t1066 = mrSges(3,1) * t1044 - mrSges(3,3) * t1086;
t1001 = -mrSges(4,1) * t1010 + mrSges(4,2) * t1011;
t1006 = -mrSges(4,2) * t1023 + mrSges(4,3) * t1010;
t1046 = sin(qJ(5));
t1050 = cos(qJ(5));
t1045 = sin(qJ(6));
t1049 = cos(qJ(6));
t1019 = t1023 ^ 2;
t999 = t1010 * t1041 - t1011 * t1037;
t929 = 0.2e1 * qJD(4) * t999 + t1037 * t936 + t1041 * t939;
t972 = -pkin(4) * t999 - pkin(10) * t1000;
t927 = -pkin(4) * t1019 + pkin(10) * t1020 + t972 * t999 + t929;
t959 = -t1039 * t980 + t1043 * t988;
t945 = -t1002 * pkin(3) - t1009 * qJ(4) + t1011 * t1007 + qJDD(4) + t959;
t975 = t1002 * t1041 - t1003 * t1037;
t976 = t1002 * t1037 + t1003 * t1041;
t931 = (-t1023 * t999 - t976) * pkin(10) + (t1000 * t1023 - t975) * pkin(4) + t945;
t923 = t1046 * t931 + t1050 * t927;
t986 = -t1000 * t1046 + t1023 * t1050;
t987 = t1000 * t1050 + t1023 * t1046;
t962 = -pkin(5) * t986 - pkin(11) * t987;
t974 = qJDD(5) - t975;
t998 = qJD(5) - t999;
t997 = t998 ^ 2;
t921 = -pkin(5) * t997 + pkin(11) * t974 + t962 * t986 + t923;
t926 = -t1020 * pkin(4) - t1019 * pkin(10) + t1000 * t972 - t928;
t957 = -qJD(5) * t987 + t1020 * t1050 - t1046 * t976;
t958 = qJD(5) * t986 + t1020 * t1046 + t1050 * t976;
t924 = (-t986 * t998 - t958) * pkin(11) + (t987 * t998 - t957) * pkin(5) + t926;
t917 = -t1045 * t921 + t1049 * t924;
t963 = -t1045 * t987 + t1049 * t998;
t934 = qJD(6) * t963 + t1045 * t974 + t1049 * t958;
t964 = t1045 * t998 + t1049 * t987;
t948 = -mrSges(7,1) * t963 + mrSges(7,2) * t964;
t985 = qJD(6) - t986;
t949 = -mrSges(7,2) * t985 + mrSges(7,3) * t963;
t956 = qJDD(6) - t957;
t915 = m(7) * t917 + mrSges(7,1) * t956 - mrSges(7,3) * t934 - t948 * t964 + t949 * t985;
t918 = t1045 * t924 + t1049 * t921;
t933 = -qJD(6) * t964 - t1045 * t958 + t1049 * t974;
t950 = mrSges(7,1) * t985 - mrSges(7,3) * t964;
t916 = m(7) * t918 - mrSges(7,2) * t956 + mrSges(7,3) * t933 + t948 * t963 - t950 * t985;
t909 = -t1045 * t915 + t1049 * t916;
t961 = -mrSges(6,1) * t986 + mrSges(6,2) * t987;
t966 = mrSges(6,1) * t998 - mrSges(6,3) * t987;
t907 = m(6) * t923 - mrSges(6,2) * t974 + mrSges(6,3) * t957 + t961 * t986 - t966 * t998 + t909;
t922 = -t1046 * t927 + t1050 * t931;
t920 = -pkin(5) * t974 - pkin(11) * t997 + t962 * t987 - t922;
t919 = -m(7) * t920 + t933 * mrSges(7,1) - mrSges(7,2) * t934 + t963 * t949 - t950 * t964;
t965 = -mrSges(6,2) * t998 + mrSges(6,3) * t986;
t913 = m(6) * t922 + mrSges(6,1) * t974 - mrSges(6,3) * t958 - t961 * t987 + t965 * t998 + t919;
t1072 = -t1046 * t913 + t1050 * t907;
t971 = -mrSges(5,1) * t999 + mrSges(5,2) * t1000;
t990 = mrSges(5,1) * t1023 - mrSges(5,3) * t1000;
t898 = m(5) * t929 - mrSges(5,2) * t1020 + mrSges(5,3) * t975 - t1023 * t990 + t971 * t999 + t1072;
t908 = t1045 * t916 + t1049 * t915;
t1056 = -m(6) * t926 + t957 * mrSges(6,1) - mrSges(6,2) * t958 + t986 * t965 - t966 * t987 - t908;
t989 = -mrSges(5,2) * t1023 + mrSges(5,3) * t999;
t904 = m(5) * t928 + mrSges(5,1) * t1020 - mrSges(5,3) * t976 - t1000 * t971 + t1023 * t989 + t1056;
t893 = t1037 * t898 + t1041 * t904;
t891 = m(4) * t946 + mrSges(4,1) * t1020 - mrSges(4,3) * t1003 - t1001 * t1011 + t1006 * t1023 + t893;
t1008 = mrSges(4,1) * t1023 - mrSges(4,3) * t1011;
t1074 = -t1037 * t904 + t1041 * t898;
t892 = m(4) * t947 - mrSges(4,2) * t1020 + mrSges(4,3) * t1002 + t1001 * t1010 - t1008 * t1023 + t1074;
t902 = t1046 * t907 + t1050 * t913;
t901 = m(5) * t945 - mrSges(5,1) * t975 + t976 * mrSges(5,2) + t1000 * t990 - t989 * t999 + t902;
t900 = m(4) * t959 - mrSges(4,1) * t1002 + mrSges(4,2) * t1003 - t1006 * t1010 + t1008 * t1011 + t901;
t878 = -t1039 * t900 + t891 * t1078 + t892 * t1079;
t874 = m(3) * t1004 + t1066 * qJDD(1) + (-t1024 * t1086 + t1029 * t1044) * qJD(1) + t878;
t1012 = -t1040 * t1025 + t1075;
t1028 = t1066 * qJD(1);
t877 = t1043 * t900 + t891 * t1083 + t892 * t1084;
t876 = m(3) * t1012 + (t1070 * qJDD(1) + (t1028 * t1038 - t1029 * t1042) * qJD(1)) * t1040 + t877;
t1005 = -g(3) * t1086 + t1076;
t884 = -t1047 * t891 + t1051 * t892;
t883 = m(3) * t1005 + t1065 * qJDD(1) + (t1024 * t1082 - t1028 * t1044) * qJD(1) + t884;
t863 = -t1040 * t876 + t874 * t1080 + t883 * t1085;
t860 = m(2) * t1033 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1053 + t863;
t870 = -t1038 * t874 + t1042 * t883;
t868 = m(2) * t1034 - mrSges(2,1) * t1053 - qJDD(1) * mrSges(2,2) + t870;
t1090 = t1048 * t868 + t1052 * t860;
t862 = t1044 * t876 + t874 * t1082 + t883 * t1086;
t1071 = -t1048 * t860 + t1052 * t868;
t1069 = Ifges(3,5) * t1038 + Ifges(3,6) * t1042;
t940 = Ifges(7,5) * t964 + Ifges(7,6) * t963 + Ifges(7,3) * t985;
t942 = Ifges(7,1) * t964 + Ifges(7,4) * t963 + Ifges(7,5) * t985;
t910 = -mrSges(7,1) * t920 + mrSges(7,3) * t918 + Ifges(7,4) * t934 + Ifges(7,2) * t933 + Ifges(7,6) * t956 - t940 * t964 + t942 * t985;
t941 = Ifges(7,4) * t964 + Ifges(7,2) * t963 + Ifges(7,6) * t985;
t911 = mrSges(7,2) * t920 - mrSges(7,3) * t917 + Ifges(7,1) * t934 + Ifges(7,4) * t933 + Ifges(7,5) * t956 + t940 * t963 - t941 * t985;
t951 = Ifges(6,5) * t987 + Ifges(6,6) * t986 + Ifges(6,3) * t998;
t952 = Ifges(6,4) * t987 + Ifges(6,2) * t986 + Ifges(6,6) * t998;
t894 = mrSges(6,2) * t926 - mrSges(6,3) * t922 + Ifges(6,1) * t958 + Ifges(6,4) * t957 + Ifges(6,5) * t974 - pkin(11) * t908 - t1045 * t910 + t1049 * t911 + t951 * t986 - t952 * t998;
t1055 = mrSges(7,1) * t917 - mrSges(7,2) * t918 + Ifges(7,5) * t934 + Ifges(7,6) * t933 + Ifges(7,3) * t956 + t941 * t964 - t942 * t963;
t953 = Ifges(6,1) * t987 + Ifges(6,4) * t986 + Ifges(6,5) * t998;
t895 = -mrSges(6,1) * t926 + mrSges(6,3) * t923 + Ifges(6,4) * t958 + Ifges(6,2) * t957 + Ifges(6,6) * t974 - pkin(5) * t908 - t951 * t987 + t953 * t998 - t1055;
t967 = Ifges(5,5) * t1000 + Ifges(5,6) * t999 + Ifges(5,3) * t1023;
t968 = Ifges(5,4) * t1000 + Ifges(5,2) * t999 + Ifges(5,6) * t1023;
t879 = mrSges(5,2) * t945 - mrSges(5,3) * t928 + Ifges(5,1) * t976 + Ifges(5,4) * t975 + Ifges(5,5) * t1020 - pkin(10) * t902 - t1023 * t968 - t1046 * t895 + t1050 * t894 + t967 * t999;
t1054 = mrSges(6,1) * t922 - mrSges(6,2) * t923 + Ifges(6,5) * t958 + Ifges(6,6) * t957 + Ifges(6,3) * t974 + pkin(5) * t919 + pkin(11) * t909 + t1045 * t911 + t1049 * t910 + t987 * t952 - t986 * t953;
t969 = Ifges(5,1) * t1000 + Ifges(5,4) * t999 + Ifges(5,5) * t1023;
t885 = -mrSges(5,1) * t945 + mrSges(5,3) * t929 + Ifges(5,4) * t976 + Ifges(5,2) * t975 + Ifges(5,6) * t1020 - pkin(4) * t902 - t1000 * t967 + t1023 * t969 - t1054;
t991 = Ifges(4,5) * t1011 + Ifges(4,6) * t1010 + Ifges(4,3) * t1023;
t993 = Ifges(4,1) * t1011 + Ifges(4,4) * t1010 + Ifges(4,5) * t1023;
t864 = -mrSges(4,1) * t959 + mrSges(4,3) * t947 + Ifges(4,4) * t1003 + Ifges(4,2) * t1002 + Ifges(4,6) * t1020 - pkin(3) * t901 + qJ(4) * t1074 - t1011 * t991 + t1023 * t993 + t1037 * t879 + t1041 * t885;
t992 = Ifges(4,4) * t1011 + Ifges(4,2) * t1010 + Ifges(4,6) * t1023;
t865 = mrSges(4,2) * t959 - mrSges(4,3) * t946 + Ifges(4,1) * t1003 + Ifges(4,4) * t1002 + Ifges(4,5) * t1020 - qJ(4) * t893 + t1010 * t991 - t1023 * t992 - t1037 * t885 + t1041 * t879;
t1062 = pkin(9) * t884 + t1047 * t865 + t1051 * t864;
t1059 = Ifges(3,6) * t1044 + (Ifges(3,4) * t1038 + Ifges(3,2) * t1042) * t1040;
t1016 = t1059 * qJD(1);
t1060 = Ifges(3,5) * t1044 + (Ifges(3,1) * t1038 + Ifges(3,4) * t1042) * t1040;
t1017 = t1060 * qJD(1);
t871 = Ifges(4,5) * t1003 + Ifges(4,6) * t1002 + t1011 * t992 - t1010 * t993 + mrSges(4,1) * t946 - mrSges(4,2) * t947 + Ifges(5,5) * t976 + Ifges(5,6) * t975 + t1000 * t968 - t999 * t969 + mrSges(5,1) * t928 - mrSges(5,2) * t929 + t1046 * t894 + t1050 * t895 + pkin(4) * t1056 + pkin(10) * t1072 + pkin(3) * t893 + (Ifges(4,3) + Ifges(5,3)) * t1020;
t854 = qJDD(1) * t1092 + mrSges(3,1) * t1004 - mrSges(3,2) * t1005 + pkin(2) * t878 + t1043 * t871 + t1062 * t1039 + (t1069 * qJDD(1) + (t1016 * t1038 - t1017 * t1042) * qJD(1)) * t1040;
t1015 = (t1069 * t1040 + t1092) * qJD(1);
t856 = -mrSges(3,1) * t1012 + mrSges(3,3) * t1005 - pkin(2) * t877 - t1039 * t871 + (-t1015 * t1086 + t1017 * t1044) * qJD(1) + t1062 * t1043 + t1059 * qJDD(1);
t858 = mrSges(3,2) * t1012 - mrSges(3,3) * t1004 - t1047 * t864 + t1051 * t865 + (t1015 * t1082 - t1016 * t1044) * qJD(1) + (-t1039 * t877 - t1043 * t878) * pkin(9) + t1060 * qJDD(1);
t1061 = mrSges(2,1) * t1033 - mrSges(2,2) * t1034 + Ifges(2,3) * qJDD(1) + pkin(1) * t863 + t1044 * t854 + t856 * t1082 + t858 * t1086 + t870 * t1089;
t852 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1033 + Ifges(2,5) * qJDD(1) - t1053 * Ifges(2,6) - t1038 * t856 + t1042 * t858 + (-t1040 * t862 - t1044 * t863) * qJ(2);
t851 = mrSges(2,1) * g(3) + mrSges(2,3) * t1034 + t1053 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t862 - t1040 * t854 + (qJ(2) * t870 + t1038 * t858 + t1042 * t856) * t1044;
t1 = [-m(1) * g(1) + t1071; -m(1) * g(2) + t1090; (-m(1) - m(2)) * g(3) + t862; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1090 - t1048 * t851 + t1052 * t852; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1071 + t1048 * t852 + t1052 * t851; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1061; t1061; t876; t871; t901; t1054; t1055;];
tauJB  = t1;
