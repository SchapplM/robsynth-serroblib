% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 14:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:23:51
% EndTime: 2019-05-07 14:24:14
% DurationCPUTime: 22.19s
% Computational Cost: add. (369139->380), mult. (791041->474), div. (0->0), fcn. (618389->12), ass. (0->163)
t1097 = Ifges(4,1) + Ifges(5,1);
t1090 = Ifges(4,4) - Ifges(5,5);
t1089 = Ifges(4,5) + Ifges(5,4);
t1096 = -Ifges(4,2) - Ifges(5,3);
t1088 = Ifges(5,6) - Ifges(4,6);
t1095 = Ifges(4,3) + Ifges(5,2);
t1043 = cos(pkin(6));
t1038 = qJD(1) * t1043 + qJD(2);
t1046 = sin(qJ(3));
t1042 = sin(pkin(6));
t1047 = sin(qJ(2));
t1075 = t1042 * t1047;
t1069 = qJD(1) * t1075;
t1092 = cos(qJ(3));
t1006 = -t1092 * t1038 + t1046 * t1069;
t1007 = t1046 * t1038 + t1092 * t1069;
t1051 = cos(qJ(2));
t1070 = qJDD(1) * t1042;
t1022 = -qJD(2) * t1069 + t1051 * t1070;
t1014 = qJDD(3) - t1022;
t1013 = qJDD(5) - t1014;
t1044 = sin(qJ(6));
t1049 = cos(qJ(6));
t1074 = t1042 * t1051;
t1068 = qJD(1) * t1074;
t1030 = qJD(3) - t1068;
t1025 = qJD(5) - t1030;
t1024 = t1025 ^ 2;
t1045 = sin(qJ(5));
t1050 = cos(qJ(5));
t1076 = t1006 * t1030;
t1029 = t1030 ^ 2;
t1078 = qJD(1) * t1042;
t1020 = (-pkin(2) * t1051 - pkin(9) * t1047) * t1078;
t1036 = t1038 ^ 2;
t1037 = qJDD(1) * t1043 + qJDD(2);
t1048 = sin(qJ(1));
t1052 = cos(qJ(1));
t1031 = t1048 * g(1) - g(2) * t1052;
t1053 = qJD(1) ^ 2;
t1085 = pkin(8) * t1042;
t1017 = qJDD(1) * pkin(1) + t1053 * t1085 + t1031;
t1032 = -g(1) * t1052 - g(2) * t1048;
t1018 = -pkin(1) * t1053 + pkin(8) * t1070 + t1032;
t1073 = t1043 * t1047;
t1071 = t1017 * t1073 + t1051 * t1018;
t1077 = qJD(1) * t1051;
t957 = -t1036 * pkin(2) + t1037 * pkin(9) + (-g(3) * t1047 + t1020 * t1077) * t1042 + t1071;
t1021 = (qJD(2) * t1077 + qJDD(1) * t1047) * t1042;
t1084 = t1043 * g(3);
t958 = -t1022 * pkin(2) - t1021 * pkin(9) - t1084 + (-t1017 + (pkin(2) * t1047 - pkin(9) * t1051) * t1038 * qJD(1)) * t1042;
t929 = -t1046 * t957 + t1092 * t958;
t985 = pkin(3) * t1006 - qJ(4) * t1007;
t924 = -t1014 * pkin(3) - t1029 * qJ(4) + t1007 * t985 + qJDD(4) - t929;
t979 = -t1006 * qJD(3) + t1092 * t1021 + t1046 * t1037;
t917 = (-t979 - t1076) * pkin(10) + (t1006 * t1007 - t1014) * pkin(4) + t924;
t1005 = t1006 ^ 2;
t1093 = 2 * qJD(4);
t930 = t1046 * t958 + t1092 * t957;
t923 = -pkin(3) * t1029 + t1014 * qJ(4) - t1006 * t985 + t1030 * t1093 + t930;
t978 = qJD(3) * t1007 + t1021 * t1046 - t1092 * t1037;
t994 = -pkin(4) * t1030 - pkin(10) * t1007;
t920 = -pkin(4) * t1005 + pkin(10) * t978 + t1030 * t994 + t923;
t914 = -t1045 * t920 + t1050 * t917;
t983 = t1006 * t1050 - t1007 * t1045;
t984 = t1006 * t1045 + t1007 * t1050;
t952 = -pkin(5) * t983 - pkin(11) * t984;
t910 = -pkin(5) * t1013 - pkin(11) * t1024 + t952 * t984 - t914;
t941 = qJD(5) * t983 + t1045 * t978 + t1050 * t979;
t962 = t1025 * t1044 + t1049 * t984;
t927 = -qJD(6) * t962 + t1013 * t1049 - t1044 * t941;
t961 = t1025 * t1049 - t1044 * t984;
t928 = qJD(6) * t961 + t1013 * t1044 + t1049 * t941;
t982 = qJD(6) - t983;
t944 = -mrSges(7,2) * t982 + mrSges(7,3) * t961;
t945 = mrSges(7,1) * t982 - mrSges(7,3) * t962;
t1061 = -m(7) * t910 + t927 * mrSges(7,1) - mrSges(7,2) * t928 + t961 * t944 - t945 * t962;
t915 = t1045 * t917 + t1050 * t920;
t911 = -pkin(5) * t1024 + pkin(11) * t1013 + t952 * t983 + t915;
t1072 = t1043 * t1051;
t980 = -g(3) * t1074 + t1017 * t1072 - t1047 * t1018;
t956 = -t1037 * pkin(2) - t1036 * pkin(9) + t1020 * t1069 - t980;
t1060 = t978 * pkin(3) + t956 + (t1076 - t979) * qJ(4);
t1086 = pkin(3) * t1030;
t921 = -t978 * pkin(4) - t1005 * pkin(10) - t1060 + (-t1086 + t1093 + t994) * t1007;
t940 = -qJD(5) * t984 - t1045 * t979 + t1050 * t978;
t912 = t921 + (-t1025 * t983 - t941) * pkin(11) + (t1025 * t984 - t940) * pkin(5);
t908 = -t1044 * t911 + t1049 * t912;
t939 = qJDD(6) - t940;
t943 = -mrSges(7,1) * t961 + mrSges(7,2) * t962;
t904 = m(7) * t908 + mrSges(7,1) * t939 - mrSges(7,3) * t928 - t943 * t962 + t944 * t982;
t909 = t1044 * t912 + t1049 * t911;
t905 = m(7) * t909 - mrSges(7,2) * t939 + mrSges(7,3) * t927 + t943 * t961 - t945 * t982;
t1067 = -t1044 * t904 + t1049 * t905;
t931 = Ifges(7,5) * t962 + Ifges(7,6) * t961 + Ifges(7,3) * t982;
t933 = Ifges(7,1) * t962 + Ifges(7,4) * t961 + Ifges(7,5) * t982;
t898 = -mrSges(7,1) * t910 + mrSges(7,3) * t909 + Ifges(7,4) * t928 + Ifges(7,2) * t927 + Ifges(7,6) * t939 - t931 * t962 + t933 * t982;
t932 = Ifges(7,4) * t962 + Ifges(7,2) * t961 + Ifges(7,6) * t982;
t899 = mrSges(7,2) * t910 - mrSges(7,3) * t908 + Ifges(7,1) * t928 + Ifges(7,4) * t927 + Ifges(7,5) * t939 + t931 * t961 - t932 * t982;
t947 = Ifges(6,4) * t984 + Ifges(6,2) * t983 + Ifges(6,6) * t1025;
t948 = Ifges(6,1) * t984 + Ifges(6,4) * t983 + Ifges(6,5) * t1025;
t1057 = -mrSges(6,1) * t914 + mrSges(6,2) * t915 - Ifges(6,5) * t941 - Ifges(6,6) * t940 - Ifges(6,3) * t1013 - pkin(5) * t1061 - pkin(11) * t1067 - t1044 * t899 - t1049 * t898 - t984 * t947 + t983 * t948;
t951 = -mrSges(6,1) * t983 + mrSges(6,2) * t984;
t964 = mrSges(6,1) * t1025 - mrSges(6,3) * t984;
t892 = m(6) * t915 - mrSges(6,2) * t1013 + mrSges(6,3) * t940 - t1025 * t964 + t951 * t983 + t1067;
t963 = -mrSges(6,2) * t1025 + mrSges(6,3) * t983;
t900 = m(6) * t914 + mrSges(6,1) * t1013 - mrSges(6,3) * t941 + t1025 * t963 - t951 * t984 + t1061;
t1066 = -t1045 * t900 + t1050 * t892;
t992 = -mrSges(5,1) * t1030 + mrSges(5,2) * t1007;
t1063 = m(5) * t923 + t1014 * mrSges(5,3) + t1030 * t992 + t1066;
t1080 = -t1090 * t1006 + t1007 * t1097 + t1089 * t1030;
t1081 = t1006 * t1096 + t1007 * t1090 - t1030 * t1088;
t886 = t1045 * t892 + t1050 * t900;
t993 = -mrSges(5,2) * t1006 + mrSges(5,3) * t1030;
t1059 = -m(5) * t924 + t1014 * mrSges(5,1) + t1030 * t993 - t886;
t986 = mrSges(5,1) * t1006 - mrSges(5,3) * t1007;
t885 = t979 * mrSges(5,2) + t1007 * t986 - t1059;
t1094 = t1080 * t1006 + t1081 * t1007 + t1095 * t1014 + t1088 * t978 + t1089 * t979 + mrSges(4,1) * t929 - mrSges(5,1) * t924 - mrSges(4,2) * t930 + mrSges(5,3) * t923 - pkin(3) * t885 - pkin(4) * t886 + qJ(4) * (-t978 * mrSges(5,2) - t1006 * t986 + t1063) + t1057;
t1091 = -mrSges(4,3) - mrSges(5,2);
t1015 = mrSges(3,1) * t1038 - mrSges(3,3) * t1069;
t1019 = (-mrSges(3,1) * t1051 + mrSges(3,2) * t1047) * t1078;
t1079 = -mrSges(4,1) * t1006 - mrSges(4,2) * t1007 - t986;
t991 = mrSges(4,1) * t1030 - mrSges(4,3) * t1007;
t881 = m(4) * t930 - t1014 * mrSges(4,2) + t1079 * t1006 - t1030 * t991 + t1091 * t978 + t1063;
t990 = -mrSges(4,2) * t1030 - mrSges(4,3) * t1006;
t883 = m(4) * t929 + t1014 * mrSges(4,1) + t1079 * t1007 + t1030 * t990 + t1091 * t979 + t1059;
t1065 = -t1046 * t883 + t1092 * t881;
t981 = -g(3) * t1075 + t1071;
t873 = m(3) * t981 - mrSges(3,2) * t1037 + mrSges(3,3) * t1022 - t1015 * t1038 + t1019 * t1068 + t1065;
t1016 = -mrSges(3,2) * t1038 + mrSges(3,3) * t1068;
t876 = t1046 * t881 + t1092 * t883;
t999 = -t1042 * t1017 - t1084;
t875 = m(3) * t999 - t1022 * mrSges(3,1) + t1021 * mrSges(3,2) + (t1015 * t1047 - t1016 * t1051) * t1078 + t876;
t894 = t1044 * t905 + t1049 * t904;
t1062 = -m(6) * t921 + t940 * mrSges(6,1) - t941 * mrSges(6,2) + t983 * t963 - t984 * t964 - t894;
t925 = (-(2 * qJD(4)) + t1086) * t1007 + t1060;
t890 = m(5) * t925 + t978 * mrSges(5,1) - t979 * mrSges(5,3) + t1006 * t993 - t1007 * t992 + t1062;
t1055 = -m(4) * t956 - t978 * mrSges(4,1) - t979 * mrSges(4,2) - t1006 * t990 - t1007 * t991 - t890;
t889 = m(3) * t980 + t1037 * mrSges(3,1) - t1021 * mrSges(3,3) + t1038 * t1016 - t1019 * t1069 + t1055;
t863 = -t1042 * t875 + t889 * t1072 + t873 * t1073;
t860 = m(2) * t1031 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1053 + t863;
t869 = -t1047 * t889 + t1051 * t873;
t867 = m(2) * t1032 - mrSges(2,1) * t1053 - qJDD(1) * mrSges(2,2) + t869;
t1083 = t1048 * t867 + t1052 * t860;
t1082 = -t1006 * t1088 - t1007 * t1089 - t1030 * t1095;
t862 = t1043 * t875 + t889 * t1074 + t873 * t1075;
t1064 = -t1048 * t860 + t1052 * t867;
t946 = Ifges(6,5) * t984 + Ifges(6,6) * t983 + Ifges(6,3) * t1025;
t877 = mrSges(6,2) * t921 - mrSges(6,3) * t914 + Ifges(6,1) * t941 + Ifges(6,4) * t940 + Ifges(6,5) * t1013 - pkin(11) * t894 - t1025 * t947 - t1044 * t898 + t1049 * t899 + t946 * t983;
t1056 = mrSges(7,1) * t908 - mrSges(7,2) * t909 + Ifges(7,5) * t928 + Ifges(7,6) * t927 + Ifges(7,3) * t939 + t932 * t962 - t933 * t961;
t878 = -mrSges(6,1) * t921 + mrSges(6,3) * t915 + Ifges(6,4) * t941 + Ifges(6,2) * t940 + Ifges(6,6) * t1013 - pkin(5) * t894 + t1025 * t948 - t946 * t984 - t1056;
t858 = -mrSges(4,1) * t956 - mrSges(5,1) * t925 + mrSges(5,2) * t923 + mrSges(4,3) * t930 - pkin(3) * t890 - pkin(4) * t1062 - pkin(10) * t1066 + t1082 * t1007 - t1088 * t1014 + t1080 * t1030 - t1045 * t877 - t1050 * t878 + t1090 * t979 + t1096 * t978;
t864 = mrSges(4,2) * t956 + mrSges(5,2) * t924 - mrSges(4,3) * t929 - mrSges(5,3) * t925 - pkin(10) * t886 - qJ(4) * t890 + t1082 * t1006 + t1089 * t1014 - t1081 * t1030 - t1045 * t878 + t1050 * t877 - t1090 * t978 + t1097 * t979;
t997 = Ifges(3,6) * t1038 + (Ifges(3,4) * t1047 + Ifges(3,2) * t1051) * t1078;
t998 = Ifges(3,5) * t1038 + (Ifges(3,1) * t1047 + Ifges(3,4) * t1051) * t1078;
t853 = Ifges(3,5) * t1021 + Ifges(3,6) * t1022 + Ifges(3,3) * t1037 + mrSges(3,1) * t980 - mrSges(3,2) * t981 + t1046 * t864 + t1092 * t858 + pkin(2) * t1055 + pkin(9) * t1065 + (t1047 * t997 - t1051 * t998) * t1078;
t996 = Ifges(3,3) * t1038 + (Ifges(3,5) * t1047 + Ifges(3,6) * t1051) * t1078;
t855 = mrSges(3,2) * t999 - mrSges(3,3) * t980 + Ifges(3,1) * t1021 + Ifges(3,4) * t1022 + Ifges(3,5) * t1037 - pkin(9) * t876 - t1038 * t997 - t1046 * t858 + t996 * t1068 + t1092 * t864;
t857 = -mrSges(3,1) * t999 + mrSges(3,3) * t981 + Ifges(3,4) * t1021 + Ifges(3,2) * t1022 + Ifges(3,6) * t1037 - pkin(2) * t876 + t1038 * t998 - t996 * t1069 - t1094;
t1058 = mrSges(2,1) * t1031 - mrSges(2,2) * t1032 + Ifges(2,3) * qJDD(1) + pkin(1) * t863 + t1043 * t853 + t857 * t1074 + t855 * t1075 + t869 * t1085;
t851 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1031 + Ifges(2,5) * qJDD(1) - t1053 * Ifges(2,6) - t1047 * t857 + t1051 * t855 + (-t1042 * t862 - t1043 * t863) * pkin(8);
t850 = mrSges(2,1) * g(3) + mrSges(2,3) * t1032 + t1053 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t862 - t1042 * t853 + (pkin(8) * t869 + t1047 * t855 + t1051 * t857) * t1043;
t1 = [-m(1) * g(1) + t1064; -m(1) * g(2) + t1083; (-m(1) - m(2)) * g(3) + t862; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1083 - t1048 * t850 + t1052 * t851; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1064 + t1048 * t851 + t1052 * t850; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1058; t1058; t853; t1094; t885; -t1057; t1056;];
tauJB  = t1;
