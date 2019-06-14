% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR9
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
% Datum: 2019-05-08 16:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 15:55:10
% EndTime: 2019-05-08 15:58:58
% DurationCPUTime: 186.67s
% Computational Cost: add. (3110771->420), mult. (7651729->556), div. (0->0), fcn. (6514976->16), ass. (0->183)
t1045 = sin(pkin(7));
t1047 = cos(pkin(7));
t1052 = sin(qJ(3));
t1058 = cos(qJ(3));
t1048 = cos(pkin(6));
t1042 = t1048 * qJD(1) + qJD(2);
t1046 = sin(pkin(6));
t1059 = cos(qJ(2));
t1085 = t1046 * t1059;
t1078 = qJD(1) * t1085;
t1026 = (t1042 * t1045 + t1047 * t1078) * pkin(10);
t1053 = sin(qJ(2));
t1091 = qJD(1) * t1046;
t1095 = pkin(10) * t1045;
t1030 = (-pkin(2) * t1059 - t1053 * t1095) * t1091;
t1089 = qJD(1) * t1059;
t1036 = (qJD(2) * t1089 + qJDD(1) * t1053) * t1046;
t1041 = t1048 * qJDD(1) + qJDD(2);
t1054 = sin(qJ(1));
t1060 = cos(qJ(1));
t1039 = t1054 * g(1) - t1060 * g(2);
t1061 = qJD(1) ^ 2;
t1096 = pkin(9) * t1046;
t1033 = qJDD(1) * pkin(1) + t1061 * t1096 + t1039;
t1040 = -t1060 * g(1) - t1054 * g(2);
t1034 = -t1061 * pkin(1) + qJDD(1) * t1096 + t1040;
t1081 = t1048 * t1059;
t1073 = t1033 * t1081 - t1053 * t1034;
t1090 = qJD(1) * t1053;
t1094 = pkin(10) * t1047;
t982 = -t1036 * t1094 + t1041 * pkin(2) + t1042 * t1026 + (-g(3) * t1059 - t1030 * t1090) * t1046 + t1073;
t1086 = t1046 * t1053;
t1079 = qJD(1) * t1086;
t1029 = t1042 * pkin(2) - t1079 * t1094;
t1037 = (-qJD(2) * t1090 + qJDD(1) * t1059) * t1046;
t1070 = t1037 * t1047 + t1041 * t1045;
t1082 = t1048 * t1053;
t1080 = t1033 * t1082 + t1059 * t1034;
t983 = -t1042 * t1029 + (-g(3) * t1053 + t1030 * t1089) * t1046 + t1070 * pkin(10) + t1080;
t1093 = t1048 * g(3);
t988 = -t1036 * t1095 - t1037 * pkin(2) - t1093 + (-t1033 + (-t1026 * t1059 + t1029 * t1053) * qJD(1)) * t1046;
t953 = -t1052 * t983 + (t1045 * t988 + t1047 * t982) * t1058;
t1083 = t1047 * t1059;
t1088 = t1045 * t1052;
t1017 = t1042 * t1088 + (t1052 * t1083 + t1053 * t1058) * t1091;
t1000 = -t1017 * qJD(3) - t1052 * t1036 + t1058 * t1070;
t1087 = t1045 * t1058;
t1016 = (-t1052 * t1053 + t1058 * t1083) * t1091 + t1042 * t1087;
t1011 = -g(3) * t1085 + t1073;
t1032 = -t1042 * mrSges(3,2) + mrSges(3,3) * t1078;
t1035 = (-mrSges(3,1) * t1059 + mrSges(3,2) * t1053) * t1091;
t1084 = t1047 * t1052;
t1002 = -t1016 * mrSges(4,1) + t1017 * mrSges(4,2);
t1027 = t1047 * t1042 - t1045 * t1078 + qJD(3);
t1010 = t1027 * mrSges(4,1) - t1017 * mrSges(4,3);
t1018 = -t1045 * t1037 + t1047 * t1041 + qJDD(3);
t1051 = sin(qJ(4));
t1057 = cos(qJ(4));
t1007 = -t1051 * t1017 + t1057 * t1027;
t1015 = qJD(4) - t1016;
t1050 = sin(qJ(5));
t1056 = cos(qJ(5));
t1006 = qJD(5) - t1007;
t1049 = sin(qJ(6));
t1055 = cos(qJ(6));
t1004 = qJD(6) + t1006;
t1014 = t1015 ^ 2;
t1003 = -t1016 * pkin(3) - t1017 * pkin(11);
t1025 = t1027 ^ 2;
t954 = t1058 * t983 + t982 * t1084 + t988 * t1088;
t942 = -t1025 * pkin(3) + t1018 * pkin(11) + t1016 * t1003 + t954;
t1001 = t1016 * qJD(3) + t1058 * t1036 + t1052 * t1070;
t961 = -t1045 * t982 + t1047 * t988;
t944 = (-t1016 * t1027 - t1001) * pkin(11) + (t1017 * t1027 - t1000) * pkin(3) + t961;
t932 = t1051 * t944 + t1057 * t942;
t1008 = t1057 * t1017 + t1051 * t1027;
t985 = -t1007 * pkin(4) - t1008 * pkin(12);
t999 = qJDD(4) - t1000;
t927 = -t1014 * pkin(4) + t999 * pkin(12) + t1007 * t985 + t932;
t941 = -t1018 * pkin(3) - t1025 * pkin(11) + t1017 * t1003 - t953;
t969 = -t1008 * qJD(4) - t1051 * t1001 + t1057 * t1018;
t970 = t1007 * qJD(4) + t1057 * t1001 + t1051 * t1018;
t930 = (-t1007 * t1015 - t970) * pkin(12) + (t1008 * t1015 - t969) * pkin(4) + t941;
t922 = -t1050 * t927 + t1056 * t930;
t991 = -t1050 * t1008 + t1056 * t1015;
t952 = t991 * qJD(5) + t1050 * t999 + t1056 * t970;
t968 = qJDD(5) - t969;
t992 = t1056 * t1008 + t1050 * t1015;
t920 = (t1006 * t991 - t952) * pkin(13) + (t991 * t992 + t968) * pkin(5) + t922;
t923 = t1050 * t930 + t1056 * t927;
t951 = -t992 * qJD(5) - t1050 * t970 + t1056 * t999;
t974 = t1006 * pkin(5) - t992 * pkin(13);
t990 = t991 ^ 2;
t921 = -t990 * pkin(5) + t951 * pkin(13) - t1006 * t974 + t923;
t918 = -t1049 * t921 + t1055 * t920;
t962 = -t1049 * t992 + t1055 * t991;
t937 = t962 * qJD(6) + t1049 * t951 + t1055 * t952;
t963 = t1049 * t991 + t1055 * t992;
t949 = -t962 * mrSges(7,1) + t963 * mrSges(7,2);
t955 = -t1004 * mrSges(7,2) + t962 * mrSges(7,3);
t966 = qJDD(6) + t968;
t914 = m(7) * t918 + t966 * mrSges(7,1) - t937 * mrSges(7,3) + t1004 * t955 - t963 * t949;
t919 = t1049 * t920 + t1055 * t921;
t936 = -t963 * qJD(6) - t1049 * t952 + t1055 * t951;
t956 = t1004 * mrSges(7,1) - t963 * mrSges(7,3);
t915 = m(7) * t919 - t966 * mrSges(7,2) + t936 * mrSges(7,3) - t1004 * t956 + t962 * t949;
t906 = t1049 * t915 + t1055 * t914;
t964 = -t991 * mrSges(6,1) + t992 * mrSges(6,2);
t972 = -t1006 * mrSges(6,2) + t991 * mrSges(6,3);
t904 = m(6) * t922 + t968 * mrSges(6,1) - t952 * mrSges(6,3) + t1006 * t972 - t992 * t964 + t906;
t1076 = -t1049 * t914 + t1055 * t915;
t973 = t1006 * mrSges(6,1) - t992 * mrSges(6,3);
t905 = m(6) * t923 - t968 * mrSges(6,2) + t951 * mrSges(6,3) - t1006 * t973 + t991 * t964 + t1076;
t902 = -t1050 * t904 + t1056 * t905;
t984 = -t1007 * mrSges(5,1) + t1008 * mrSges(5,2);
t994 = t1015 * mrSges(5,1) - t1008 * mrSges(5,3);
t900 = m(5) * t932 - t999 * mrSges(5,2) + t969 * mrSges(5,3) + t1007 * t984 - t1015 * t994 + t902;
t931 = -t1051 * t942 + t1057 * t944;
t926 = -t999 * pkin(4) - t1014 * pkin(12) + t1008 * t985 - t931;
t924 = -t951 * pkin(5) - t990 * pkin(13) + t992 * t974 + t926;
t1067 = m(7) * t924 - t936 * mrSges(7,1) + t937 * mrSges(7,2) - t962 * t955 + t963 * t956;
t916 = -m(6) * t926 + t951 * mrSges(6,1) - t952 * mrSges(6,2) + t991 * t972 - t992 * t973 - t1067;
t993 = -t1015 * mrSges(5,2) + t1007 * mrSges(5,3);
t910 = m(5) * t931 + t999 * mrSges(5,1) - t970 * mrSges(5,3) - t1008 * t984 + t1015 * t993 + t916;
t1075 = -t1051 * t910 + t1057 * t900;
t889 = m(4) * t954 - t1018 * mrSges(4,2) + t1000 * mrSges(4,3) + t1016 * t1002 - t1027 * t1010 + t1075;
t1009 = -t1027 * mrSges(4,2) + t1016 * mrSges(4,3);
t892 = t1051 * t900 + t1057 * t910;
t891 = m(4) * t961 - t1000 * mrSges(4,1) + t1001 * mrSges(4,2) - t1016 * t1009 + t1017 * t1010 + t892;
t901 = t1050 * t905 + t1056 * t904;
t1064 = -m(5) * t941 + t969 * mrSges(5,1) - t970 * mrSges(5,2) + t1007 * t993 - t1008 * t994 - t901;
t897 = m(4) * t953 + t1018 * mrSges(4,1) - t1001 * mrSges(4,3) - t1017 * t1002 + t1027 * t1009 + t1064;
t878 = t1047 * t1058 * t897 - t1045 * t891 + t889 * t1084;
t874 = m(3) * t1011 + t1041 * mrSges(3,1) - t1036 * mrSges(3,3) + t1042 * t1032 - t1035 * t1079 + t878;
t1022 = -t1046 * t1033 - t1093;
t1031 = t1042 * mrSges(3,1) - mrSges(3,3) * t1079;
t877 = t1047 * t891 + t897 * t1087 + t889 * t1088;
t876 = m(3) * t1022 - t1037 * mrSges(3,1) + t1036 * mrSges(3,2) + (t1031 * t1053 - t1032 * t1059) * t1091 + t877;
t1012 = -g(3) * t1086 + t1080;
t884 = -t1052 * t897 + t1058 * t889;
t883 = m(3) * t1012 - t1041 * mrSges(3,2) + t1037 * mrSges(3,3) - t1042 * t1031 + t1035 * t1078 + t884;
t863 = -t1046 * t876 + t874 * t1081 + t883 * t1082;
t860 = m(2) * t1039 + qJDD(1) * mrSges(2,1) - t1061 * mrSges(2,2) + t863;
t870 = -t1053 * t874 + t1059 * t883;
t868 = m(2) * t1040 - t1061 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t870;
t1092 = t1054 * t868 + t1060 * t860;
t862 = t1048 * t876 + t874 * t1085 + t883 * t1086;
t1074 = -t1054 * t860 + t1060 * t868;
t945 = Ifges(7,5) * t963 + Ifges(7,6) * t962 + Ifges(7,3) * t1004;
t947 = Ifges(7,1) * t963 + Ifges(7,4) * t962 + Ifges(7,5) * t1004;
t907 = -mrSges(7,1) * t924 + mrSges(7,3) * t919 + Ifges(7,4) * t937 + Ifges(7,2) * t936 + Ifges(7,6) * t966 + t1004 * t947 - t963 * t945;
t946 = Ifges(7,4) * t963 + Ifges(7,2) * t962 + Ifges(7,6) * t1004;
t908 = mrSges(7,2) * t924 - mrSges(7,3) * t918 + Ifges(7,1) * t937 + Ifges(7,4) * t936 + Ifges(7,5) * t966 - t1004 * t946 + t962 * t945;
t957 = Ifges(6,5) * t992 + Ifges(6,6) * t991 + Ifges(6,3) * t1006;
t959 = Ifges(6,1) * t992 + Ifges(6,4) * t991 + Ifges(6,5) * t1006;
t893 = -mrSges(6,1) * t926 + mrSges(6,3) * t923 + Ifges(6,4) * t952 + Ifges(6,2) * t951 + Ifges(6,6) * t968 - pkin(5) * t1067 + pkin(13) * t1076 + t1006 * t959 + t1049 * t908 + t1055 * t907 - t992 * t957;
t958 = Ifges(6,4) * t992 + Ifges(6,2) * t991 + Ifges(6,6) * t1006;
t894 = mrSges(6,2) * t926 - mrSges(6,3) * t922 + Ifges(6,1) * t952 + Ifges(6,4) * t951 + Ifges(6,5) * t968 - pkin(13) * t906 - t1006 * t958 - t1049 * t907 + t1055 * t908 + t991 * t957;
t975 = Ifges(5,5) * t1008 + Ifges(5,6) * t1007 + Ifges(5,3) * t1015;
t976 = Ifges(5,4) * t1008 + Ifges(5,2) * t1007 + Ifges(5,6) * t1015;
t879 = mrSges(5,2) * t941 - mrSges(5,3) * t931 + Ifges(5,1) * t970 + Ifges(5,4) * t969 + Ifges(5,5) * t999 - pkin(12) * t901 + t1007 * t975 - t1015 * t976 - t1050 * t893 + t1056 * t894;
t1065 = -mrSges(7,1) * t918 + mrSges(7,2) * t919 - Ifges(7,5) * t937 - Ifges(7,6) * t936 - Ifges(7,3) * t966 - t963 * t946 + t962 * t947;
t1062 = mrSges(6,1) * t922 - mrSges(6,2) * t923 + Ifges(6,5) * t952 + Ifges(6,6) * t951 + Ifges(6,3) * t968 + pkin(5) * t906 + t992 * t958 - t991 * t959 - t1065;
t977 = Ifges(5,1) * t1008 + Ifges(5,4) * t1007 + Ifges(5,5) * t1015;
t885 = -mrSges(5,1) * t941 + mrSges(5,3) * t932 + Ifges(5,4) * t970 + Ifges(5,2) * t969 + Ifges(5,6) * t999 - pkin(4) * t901 - t1008 * t975 + t1015 * t977 - t1062;
t995 = Ifges(4,5) * t1017 + Ifges(4,6) * t1016 + Ifges(4,3) * t1027;
t996 = Ifges(4,4) * t1017 + Ifges(4,2) * t1016 + Ifges(4,6) * t1027;
t865 = mrSges(4,2) * t961 - mrSges(4,3) * t953 + Ifges(4,1) * t1001 + Ifges(4,4) * t1000 + Ifges(4,5) * t1018 - pkin(11) * t892 + t1016 * t995 - t1027 * t996 - t1051 * t885 + t1057 * t879;
t1063 = mrSges(5,1) * t931 - mrSges(5,2) * t932 + Ifges(5,5) * t970 + Ifges(5,6) * t969 + Ifges(5,3) * t999 + pkin(4) * t916 + pkin(12) * t902 - t1007 * t977 + t1008 * t976 + t1050 * t894 + t1056 * t893;
t997 = Ifges(4,1) * t1017 + Ifges(4,4) * t1016 + Ifges(4,5) * t1027;
t871 = -mrSges(4,1) * t961 + mrSges(4,3) * t954 + Ifges(4,4) * t1001 + Ifges(4,2) * t1000 + Ifges(4,6) * t1018 - pkin(3) * t892 - t1017 * t995 + t1027 * t997 - t1063;
t1068 = pkin(10) * t884 + t1052 * t865 + t1058 * t871;
t1020 = Ifges(3,6) * t1042 + (Ifges(3,4) * t1053 + Ifges(3,2) * t1059) * t1091;
t1021 = Ifges(3,5) * t1042 + (Ifges(3,1) * t1053 + Ifges(3,4) * t1059) * t1091;
t864 = mrSges(4,1) * t953 - mrSges(4,2) * t954 + Ifges(4,5) * t1001 + Ifges(4,6) * t1000 + Ifges(4,3) * t1018 + pkin(3) * t1064 + pkin(11) * t1075 - t1016 * t997 + t1017 * t996 + t1051 * t879 + t1057 * t885;
t854 = mrSges(3,1) * t1011 - mrSges(3,2) * t1012 + Ifges(3,5) * t1036 + Ifges(3,6) * t1037 + Ifges(3,3) * t1041 + pkin(2) * t878 + t1047 * t864 + (t1020 * t1053 - t1021 * t1059) * t1091 + t1068 * t1045;
t1019 = Ifges(3,3) * t1042 + (Ifges(3,5) * t1053 + Ifges(3,6) * t1059) * t1091;
t856 = -mrSges(3,1) * t1022 + mrSges(3,3) * t1012 + Ifges(3,4) * t1036 + Ifges(3,2) * t1037 + Ifges(3,6) * t1041 - pkin(2) * t877 - t1019 * t1079 + t1042 * t1021 - t1045 * t864 + t1047 * t1068;
t858 = t1019 * t1078 + mrSges(3,2) * t1022 - mrSges(3,3) * t1011 + Ifges(3,1) * t1036 + Ifges(3,4) * t1037 + Ifges(3,5) * t1041 - t1042 * t1020 - t1052 * t871 + t1058 * t865 + (-t1045 * t877 - t1047 * t878) * pkin(10);
t1066 = mrSges(2,1) * t1039 - mrSges(2,2) * t1040 + Ifges(2,3) * qJDD(1) + pkin(1) * t863 + t1048 * t854 + t856 * t1085 + t858 * t1086 + t1096 * t870;
t852 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1039 + Ifges(2,5) * qJDD(1) - t1061 * Ifges(2,6) - t1053 * t856 + t1059 * t858 + (-t1046 * t862 - t1048 * t863) * pkin(9);
t851 = mrSges(2,1) * g(3) + mrSges(2,3) * t1040 + t1061 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t862 - t1046 * t854 + (pkin(9) * t870 + t1053 * t858 + t1059 * t856) * t1048;
t1 = [-m(1) * g(1) + t1074; -m(1) * g(2) + t1092; (-m(1) - m(2)) * g(3) + t862; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1092 - t1054 * t851 + t1060 * t852; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1074 + t1054 * t852 + t1060 * t851; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1066; t1066; t854; t864; t1063; t1062; -t1065;];
tauJB  = t1;
