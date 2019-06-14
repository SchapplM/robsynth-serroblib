% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-05-07 01:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR13_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 01:10:15
% EndTime: 2019-05-07 01:10:45
% DurationCPUTime: 22.60s
% Computational Cost: add. (374414->386), mult. (830088->480), div. (0->0), fcn. (616638->12), ass. (0->165)
t1101 = -2 * qJD(3);
t1100 = Ifges(3,1) + Ifges(4,2);
t1094 = Ifges(3,4) + Ifges(4,6);
t1093 = Ifges(3,5) - Ifges(4,4);
t1099 = Ifges(3,2) + Ifges(4,3);
t1092 = Ifges(3,6) - Ifges(4,5);
t1098 = Ifges(3,3) + Ifges(4,1);
t1045 = cos(pkin(6));
t1038 = t1045 * qJD(1) + qJD(2);
t1044 = sin(pkin(6));
t1049 = sin(qJ(2));
t1082 = t1044 * t1049;
t1074 = qJD(1) * t1082;
t1097 = (pkin(2) * t1038 + t1101) * t1074;
t1054 = cos(qJ(2));
t1085 = qJD(1) * t1044;
t1016 = (-pkin(2) * t1054 - qJ(3) * t1049) * t1085;
t1036 = t1038 ^ 2;
t1037 = t1045 * qJDD(1) + qJDD(2);
t1081 = t1044 * t1054;
t1076 = qJD(1) * t1081;
t1050 = sin(qJ(1));
t1055 = cos(qJ(1));
t1033 = t1050 * g(1) - t1055 * g(2);
t1056 = qJD(1) ^ 2;
t1091 = pkin(8) * t1044;
t1014 = qJDD(1) * pkin(1) + t1056 * t1091 + t1033;
t1034 = -t1055 * g(1) - t1050 * g(2);
t1077 = qJDD(1) * t1044;
t1015 = -t1056 * pkin(1) + pkin(8) * t1077 + t1034;
t1080 = t1045 * t1049;
t977 = -g(3) * t1082 + t1014 * t1080 + t1054 * t1015;
t949 = t1036 * pkin(2) - t1037 * qJ(3) - t1016 * t1076 + t1038 * t1101 - t977;
t1096 = -pkin(2) - pkin(9);
t1095 = mrSges(3,1) - mrSges(4,2);
t1090 = t1045 * g(3);
t1079 = t1045 * t1054;
t1010 = mrSges(3,1) * t1038 - mrSges(3,3) * t1074;
t1011 = -mrSges(3,2) * t1038 + mrSges(3,3) * t1076;
t1013 = mrSges(4,1) * t1074 + mrSges(4,2) * t1038;
t1020 = (qJD(1) * qJD(2) * t1054 + qJDD(1) * t1049) * t1044;
t1021 = -qJD(2) * t1074 + t1054 * t1077;
t1012 = -mrSges(4,1) * t1076 - mrSges(4,3) * t1038;
t1048 = sin(qJ(4));
t1053 = cos(qJ(4));
t1002 = -t1048 * t1038 - t1053 * t1076;
t1009 = qJDD(4) + t1020;
t1029 = qJD(4) + t1074;
t1047 = sin(qJ(5));
t1052 = cos(qJ(5));
t1001 = qJD(5) - t1002;
t1046 = sin(qJ(6));
t1051 = cos(qJ(6));
t1026 = t1029 ^ 2;
t1019 = pkin(3) * t1074 - t1038 * pkin(9);
t1083 = t1044 ^ 2 * t1056;
t1075 = t1054 ^ 2 * t1083;
t940 = -pkin(3) * t1075 - t1090 - t1020 * qJ(3) + t1096 * t1021 + (-t1014 + (-qJ(3) * t1038 * t1054 - t1019 * t1049) * qJD(1)) * t1044 + t1097;
t1078 = g(3) * t1081 + t1049 * t1015;
t1067 = -t1036 * qJ(3) + t1016 * t1074 + qJDD(3) + t1078;
t942 = t1020 * pkin(3) + t1096 * t1037 + (-pkin(3) * t1038 * t1085 - pkin(9) * t1049 * t1083 - t1014 * t1045) * t1054 + t1067;
t929 = t1048 * t942 + t1053 * t940;
t1003 = t1038 * t1053 - t1048 * t1076;
t979 = -pkin(4) * t1002 - pkin(10) * t1003;
t918 = -pkin(4) * t1026 + pkin(10) * t1009 + t1002 * t979 + t929;
t939 = t1021 * pkin(3) - pkin(9) * t1075 + t1038 * t1019 - t949;
t974 = -t1003 * qJD(4) - t1021 * t1053 - t1048 * t1037;
t975 = qJD(4) * t1002 - t1021 * t1048 + t1037 * t1053;
t924 = (-t1002 * t1029 - t975) * pkin(10) + (t1003 * t1029 - t974) * pkin(4) + t939;
t913 = -t1047 * t918 + t1052 * t924;
t981 = -t1003 * t1047 + t1029 * t1052;
t945 = qJD(5) * t981 + t1009 * t1047 + t1052 * t975;
t972 = qJDD(5) - t974;
t982 = t1003 * t1052 + t1029 * t1047;
t911 = (t1001 * t981 - t945) * pkin(11) + (t981 * t982 + t972) * pkin(5) + t913;
t914 = t1047 * t924 + t1052 * t918;
t944 = -qJD(5) * t982 + t1009 * t1052 - t1047 * t975;
t963 = pkin(5) * t1001 - pkin(11) * t982;
t980 = t981 ^ 2;
t912 = -pkin(5) * t980 + t944 * pkin(11) - t1001 * t963 + t914;
t909 = -t1046 * t912 + t1051 * t911;
t956 = -t1046 * t982 + t1051 * t981;
t926 = qJD(6) * t956 + t1046 * t944 + t1051 * t945;
t957 = t1046 * t981 + t1051 * t982;
t935 = -mrSges(7,1) * t956 + mrSges(7,2) * t957;
t999 = qJD(6) + t1001;
t946 = -mrSges(7,2) * t999 + mrSges(7,3) * t956;
t965 = qJDD(6) + t972;
t905 = m(7) * t909 + mrSges(7,1) * t965 - t926 * mrSges(7,3) - t935 * t957 + t946 * t999;
t910 = t1046 * t911 + t1051 * t912;
t925 = -qJD(6) * t957 - t1046 * t945 + t1051 * t944;
t947 = mrSges(7,1) * t999 - mrSges(7,3) * t957;
t906 = m(7) * t910 - mrSges(7,2) * t965 + t925 * mrSges(7,3) + t935 * t956 - t947 * t999;
t898 = t1046 * t906 + t1051 * t905;
t958 = -mrSges(6,1) * t981 + mrSges(6,2) * t982;
t961 = -mrSges(6,2) * t1001 + mrSges(6,3) * t981;
t896 = m(6) * t913 + mrSges(6,1) * t972 - t945 * mrSges(6,3) + t1001 * t961 - t958 * t982 + t898;
t1072 = -t1046 * t905 + t1051 * t906;
t962 = mrSges(6,1) * t1001 - mrSges(6,3) * t982;
t897 = m(6) * t914 - mrSges(6,2) * t972 + t944 * mrSges(6,3) - t1001 * t962 + t958 * t981 + t1072;
t1071 = -t1047 * t896 + t1052 * t897;
t978 = -mrSges(5,1) * t1002 + mrSges(5,2) * t1003;
t984 = mrSges(5,1) * t1029 - mrSges(5,3) * t1003;
t890 = m(5) * t929 - mrSges(5,2) * t1009 + mrSges(5,3) * t974 + t1002 * t978 - t1029 * t984 + t1071;
t928 = -t1048 * t940 + t1053 * t942;
t917 = -pkin(4) * t1009 - pkin(10) * t1026 + t1003 * t979 - t928;
t915 = -t944 * pkin(5) - pkin(11) * t980 + t963 * t982 + t917;
t1065 = m(7) * t915 - t925 * mrSges(7,1) + t926 * mrSges(7,2) - t956 * t946 + t947 * t957;
t1058 = -m(6) * t917 + t944 * mrSges(6,1) - t945 * mrSges(6,2) + t981 * t961 - t962 * t982 - t1065;
t983 = -mrSges(5,2) * t1029 + mrSges(5,3) * t1002;
t901 = m(5) * t928 + mrSges(5,1) * t1009 - mrSges(5,3) * t975 - t1003 * t978 + t1029 * t983 + t1058;
t1070 = -t1048 * t901 + t1053 * t890;
t991 = -t1044 * t1014 - t1090;
t950 = -t1021 * pkin(2) + (-t1038 * t1076 - t1020) * qJ(3) + t991 + t1097;
t1068 = m(4) * t950 - t1020 * mrSges(4,3) + t1012 * t1076 + t1070;
t876 = m(3) * t991 + t1020 * mrSges(3,2) - t1095 * t1021 + (-t1011 * t1054 + (t1010 - t1013) * t1049) * t1085 + t1068;
t1017 = (mrSges(4,2) * t1054 - mrSges(4,3) * t1049) * t1085;
t1018 = (-mrSges(3,1) * t1054 + mrSges(3,2) * t1049) * t1085;
t880 = t1048 * t890 + t1053 * t901;
t1073 = t1014 * t1079;
t955 = -t1037 * pkin(2) + t1067 - t1073;
t1066 = -m(4) * t955 - t1020 * mrSges(4,1) - t880;
t976 = t1073 - t1078;
t877 = m(3) * t976 - t1020 * mrSges(3,3) + (t1011 - t1012) * t1038 + t1095 * t1037 + (-t1017 - t1018) * t1074 + t1066;
t892 = t1047 * t897 + t1052 * t896;
t1062 = -m(5) * t939 + t974 * mrSges(5,1) - t975 * mrSges(5,2) + t1002 * t983 - t1003 * t984 - t892;
t1059 = -m(4) * t949 + t1037 * mrSges(4,3) + t1038 * t1013 + t1017 * t1076 - t1062;
t888 = t1018 * t1076 + t1059 - t1037 * mrSges(3,2) - t1038 * t1010 + m(3) * t977 + (mrSges(3,3) + mrSges(4,1)) * t1021;
t865 = -t1044 * t876 + t877 * t1079 + t888 * t1080;
t862 = m(2) * t1033 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1056 + t865;
t871 = -t1049 * t877 + t1054 * t888;
t869 = m(2) * t1034 - mrSges(2,1) * t1056 - qJDD(1) * mrSges(2,2) + t871;
t1089 = t1050 * t869 + t1055 * t862;
t1088 = (t1049 * t1093 + t1054 * t1092) * t1085 + t1098 * t1038;
t1087 = (t1049 * t1094 + t1054 * t1099) * t1085 + t1092 * t1038;
t1086 = (t1049 * t1100 + t1054 * t1094) * t1085 + t1093 * t1038;
t864 = t1045 * t876 + t877 * t1081 + t888 * t1082;
t1069 = -t1050 * t862 + t1055 * t869;
t931 = Ifges(7,5) * t957 + Ifges(7,6) * t956 + Ifges(7,3) * t999;
t933 = Ifges(7,1) * t957 + Ifges(7,4) * t956 + Ifges(7,5) * t999;
t899 = -mrSges(7,1) * t915 + mrSges(7,3) * t910 + Ifges(7,4) * t926 + Ifges(7,2) * t925 + Ifges(7,6) * t965 - t931 * t957 + t933 * t999;
t932 = Ifges(7,4) * t957 + Ifges(7,2) * t956 + Ifges(7,6) * t999;
t900 = mrSges(7,2) * t915 - mrSges(7,3) * t909 + Ifges(7,1) * t926 + Ifges(7,4) * t925 + Ifges(7,5) * t965 + t931 * t956 - t932 * t999;
t951 = Ifges(6,5) * t982 + Ifges(6,6) * t981 + Ifges(6,3) * t1001;
t953 = Ifges(6,1) * t982 + Ifges(6,4) * t981 + Ifges(6,5) * t1001;
t882 = -mrSges(6,1) * t917 + mrSges(6,3) * t914 + Ifges(6,4) * t945 + Ifges(6,2) * t944 + Ifges(6,6) * t972 - pkin(5) * t1065 + pkin(11) * t1072 + t1001 * t953 + t1046 * t900 + t1051 * t899 - t982 * t951;
t952 = Ifges(6,4) * t982 + Ifges(6,2) * t981 + Ifges(6,6) * t1001;
t884 = mrSges(6,2) * t917 - mrSges(6,3) * t913 + Ifges(6,1) * t945 + Ifges(6,4) * t944 + Ifges(6,5) * t972 - pkin(11) * t898 - t1001 * t952 - t1046 * t899 + t1051 * t900 + t951 * t981;
t966 = Ifges(5,5) * t1003 + Ifges(5,6) * t1002 + Ifges(5,3) * t1029;
t967 = Ifges(5,4) * t1003 + Ifges(5,2) * t1002 + Ifges(5,6) * t1029;
t866 = mrSges(5,2) * t939 - mrSges(5,3) * t928 + Ifges(5,1) * t975 + Ifges(5,4) * t974 + Ifges(5,5) * t1009 - pkin(10) * t892 + t1002 * t966 - t1029 * t967 - t1047 * t882 + t1052 * t884;
t1063 = -mrSges(7,1) * t909 + mrSges(7,2) * t910 - Ifges(7,5) * t926 - Ifges(7,6) * t925 - Ifges(7,3) * t965 - t957 * t932 + t956 * t933;
t1057 = mrSges(6,1) * t913 - mrSges(6,2) * t914 + Ifges(6,5) * t945 + Ifges(6,6) * t944 + Ifges(6,3) * t972 + pkin(5) * t898 + t982 * t952 - t981 * t953 - t1063;
t968 = Ifges(5,1) * t1003 + Ifges(5,4) * t1002 + Ifges(5,5) * t1029;
t872 = -mrSges(5,1) * t939 + mrSges(5,3) * t929 + Ifges(5,4) * t975 + Ifges(5,2) * t974 + Ifges(5,6) * t1009 - pkin(4) * t892 - t1003 * t966 + t1029 * t968 - t1057;
t879 = t1037 * mrSges(4,2) + t1038 * t1012 + t1017 * t1074 - t1066;
t856 = mrSges(3,1) * t976 - mrSges(3,2) * t977 + mrSges(4,2) * t955 - mrSges(4,3) * t949 + t1053 * t866 - t1048 * t872 - pkin(9) * t880 - pkin(2) * t879 + qJ(3) * t1059 + t1098 * t1037 + (mrSges(4,1) * qJ(3) + t1092) * t1021 + t1093 * t1020 + (t1049 * t1087 - t1054 * t1086) * t1085;
t878 = t1021 * mrSges(4,2) - t1013 * t1074 + t1068;
t858 = -mrSges(3,1) * t991 - mrSges(4,1) * t949 + mrSges(4,2) * t950 + mrSges(3,3) * t977 - pkin(2) * t878 - pkin(3) * t1062 - pkin(9) * t1070 + t1094 * t1020 + t1021 * t1099 + t1092 * t1037 + t1086 * t1038 - t1048 * t866 - t1053 * t872 - t1088 * t1074;
t1060 = mrSges(5,1) * t928 - mrSges(5,2) * t929 + Ifges(5,5) * t975 + Ifges(5,6) * t974 + Ifges(5,3) * t1009 + pkin(4) * t1058 + pkin(10) * t1071 - t1002 * t968 + t1003 * t967 + t1047 * t884 + t1052 * t882;
t860 = mrSges(4,1) * t955 + mrSges(3,2) * t991 - mrSges(3,3) * t976 - mrSges(4,3) * t950 + pkin(3) * t880 - qJ(3) * t878 + t1020 * t1100 + t1094 * t1021 + t1093 * t1037 - t1087 * t1038 + t1088 * t1076 + t1060;
t1064 = mrSges(2,1) * t1033 - mrSges(2,2) * t1034 + Ifges(2,3) * qJDD(1) + pkin(1) * t865 + t1045 * t856 + t858 * t1081 + t860 * t1082 + t871 * t1091;
t854 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1033 + Ifges(2,5) * qJDD(1) - t1056 * Ifges(2,6) - t1049 * t858 + t1054 * t860 + (-t1044 * t864 - t1045 * t865) * pkin(8);
t853 = mrSges(2,1) * g(3) + mrSges(2,3) * t1034 + t1056 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t864 - t1044 * t856 + (pkin(8) * t871 + t1049 * t860 + t1054 * t858) * t1045;
t1 = [-m(1) * g(1) + t1069; -m(1) * g(2) + t1089; (-m(1) - m(2)) * g(3) + t864; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1089 - t1050 * t853 + t1055 * t854; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1069 + t1050 * t854 + t1055 * t853; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1064; t1064; t856; t879; t1060; t1057; -t1063;];
tauJB  = t1;
