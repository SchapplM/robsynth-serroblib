% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
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
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:31:23
% EndTime: 2019-05-07 04:31:32
% DurationCPUTime: 7.60s
% Computational Cost: add. (67563->349), mult. (141365->406), div. (0->0), fcn. (92168->8), ass. (0->139)
t1101 = Ifges(5,1) + Ifges(4,1) + Ifges(6,2);
t1092 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t1100 = -Ifges(6,4) + Ifges(5,5) - Ifges(4,4);
t1042 = sin(qJ(3));
t1046 = cos(qJ(2));
t1075 = qJD(1) * t1046;
t1043 = sin(qJ(2));
t1076 = qJD(1) * t1043;
t1082 = cos(qJ(3));
t1003 = t1042 * t1076 - t1075 * t1082;
t1004 = (t1042 * t1046 + t1043 * t1082) * qJD(1);
t1035 = qJD(2) + qJD(3);
t1099 = t1100 * t1003 + t1101 * t1004 + t1092 * t1035;
t1098 = Ifges(6,1) + Ifges(5,3) + Ifges(4,2);
t1091 = Ifges(6,5) + Ifges(4,6) - Ifges(5,6);
t1034 = qJDD(2) + qJDD(3);
t1041 = sin(qJ(6));
t1045 = cos(qJ(6));
t1085 = 2 * qJD(4);
t1022 = t1035 * t1085;
t1087 = t1035 ^ 2;
t1013 = qJD(2) * t1075 + qJDD(1) * t1043;
t1048 = qJD(1) ^ 2;
t1044 = sin(qJ(1));
t1047 = cos(qJ(1));
t1020 = -g(1) * t1047 - g(2) * t1044;
t1006 = -pkin(1) * t1048 + qJDD(1) * pkin(7) + t1020;
t1072 = t1043 * t1006;
t935 = qJDD(2) * pkin(2) - t1013 * pkin(8) - t1072 + (pkin(2) * t1043 * t1048 + pkin(8) * qJD(1) * qJD(2) - g(3)) * t1046;
t1014 = -qJD(2) * t1076 + qJDD(1) * t1046;
t1018 = qJD(2) * pkin(2) - pkin(8) * t1076;
t1073 = t1046 ^ 2 * t1048;
t986 = -g(3) * t1043 + t1046 * t1006;
t936 = -pkin(2) * t1073 + pkin(8) * t1014 - qJD(2) * t1018 + t986;
t917 = t1042 * t935 + t1082 * t936;
t975 = pkin(3) * t1003 - qJ(4) * t1004;
t1064 = pkin(3) * t1087 - t1034 * qJ(4) + t1003 * t975 - t917;
t954 = qJD(3) * t1004 + t1013 * t1042 - t1014 * t1082;
t999 = t1003 ^ 2;
t1058 = t999 * pkin(4) - t954 * qJ(5) + t1064;
t978 = pkin(5) * t1004 - pkin(9) * t1003;
t989 = -pkin(4) * t1035 - qJ(5) * t1004;
t903 = t1034 * pkin(5) - t1087 * pkin(9) + t1035 * t989 + t1022 + ((2 * qJD(5)) + t978) * t1003 - t1058;
t984 = t1003 * t1045 - t1035 * t1041;
t922 = -qJD(6) * t984 - t1034 * t1045 - t1041 * t954;
t983 = -t1003 * t1041 - t1035 * t1045;
t923 = qJD(6) * t983 - t1034 * t1041 + t1045 * t954;
t998 = qJD(6) + t1004;
t957 = -mrSges(7,2) * t998 + mrSges(7,3) * t983;
t958 = mrSges(7,1) * t998 - mrSges(7,3) * t984;
t899 = -m(7) * t903 + t922 * mrSges(7,1) - t923 * mrSges(7,2) + t983 * t957 - t984 * t958;
t1084 = -2 * qJD(5);
t1086 = -2 * qJD(4);
t908 = t1003 * t1084 + (t1086 - t989) * t1035 + t1058;
t974 = mrSges(6,1) * t1004 + mrSges(6,2) * t1003;
t1057 = -m(6) * t908 + t954 * mrSges(6,3) + t1003 * t974 - t899;
t913 = t1022 - t1064;
t991 = -mrSges(5,1) * t1035 + mrSges(5,2) * t1004;
t1055 = m(5) * t913 + t1034 * mrSges(5,3) + t1035 * t991 + t1057;
t1089 = t1098 * t1003 + t1004 * t1100 - t1091 * t1035;
t1090 = -Ifges(5,2) - Ifges(4,3) - Ifges(6,3);
t1074 = t1003 * t1035;
t955 = -t1003 * qJD(3) + t1013 * t1082 + t1042 * t1014;
t1019 = t1044 * g(1) - t1047 * g(2);
t1005 = -qJDD(1) * pkin(1) - t1048 * pkin(7) - t1019;
t956 = -t1014 * pkin(2) - pkin(8) * t1073 + t1018 * t1076 + t1005;
t1059 = t954 * pkin(3) + t956 + (t1074 - t955) * qJ(4);
t1054 = -t999 * qJ(5) + qJDD(5) - t1059 + (t1085 + t989) * t1004;
t1083 = -pkin(4) - pkin(9);
t900 = t1054 + (-pkin(5) * t1003 + (-pkin(3) - pkin(9)) * t1004) * t1035 + t1083 * t954 + t955 * pkin(5);
t916 = -t1042 * t936 + t1082 * t935;
t1061 = -qJ(4) * t1087 + t1004 * t975 + qJDD(4) - t916;
t1053 = (-t955 - t1074) * qJ(5) + t1061 + (t1003 * pkin(4) + t1084) * t1004;
t901 = -t1087 * pkin(5) - t1004 * t978 + (-pkin(3) + t1083) * t1034 + t1053;
t897 = -t1041 * t901 + t1045 * t900;
t931 = -mrSges(7,1) * t983 + mrSges(7,2) * t984;
t952 = qJDD(6) + t955;
t894 = m(7) * t897 + mrSges(7,1) * t952 - mrSges(7,3) * t923 - t931 * t984 + t957 * t998;
t898 = t1041 * t900 + t1045 * t901;
t895 = m(7) * t898 - mrSges(7,2) * t952 + mrSges(7,3) * t922 + t931 * t983 - t958 * t998;
t883 = -t1041 * t894 + t1045 * t895;
t907 = (-pkin(3) - pkin(4)) * t1034 + t1053;
t987 = -mrSges(6,1) * t1035 - mrSges(6,3) * t1003;
t1063 = m(6) * t907 + t1034 * mrSges(6,2) + t1035 * t987 + t883;
t914 = -t1034 * pkin(3) + t1061;
t993 = -mrSges(5,2) * t1003 + mrSges(5,3) * t1035;
t1056 = m(5) * t914 - t1034 * mrSges(5,1) - t1035 * t993 + t1063;
t976 = mrSges(5,1) * t1003 - mrSges(5,3) * t1004;
t879 = (mrSges(5,2) - mrSges(6,3)) * t955 + (-t974 + t976) * t1004 + t1056;
t881 = -t955 * mrSges(6,3) - t1004 * t974 + t1063;
t924 = Ifges(7,5) * t984 + Ifges(7,6) * t983 + Ifges(7,3) * t998;
t926 = Ifges(7,1) * t984 + Ifges(7,4) * t983 + Ifges(7,5) * t998;
t887 = -mrSges(7,1) * t903 + mrSges(7,3) * t898 + Ifges(7,4) * t923 + Ifges(7,2) * t922 + Ifges(7,6) * t952 - t924 * t984 + t926 * t998;
t925 = Ifges(7,4) * t984 + Ifges(7,2) * t983 + Ifges(7,6) * t998;
t888 = mrSges(7,2) * t903 - mrSges(7,3) * t897 + Ifges(7,1) * t923 + Ifges(7,4) * t922 + Ifges(7,5) * t952 + t924 * t983 - t925 * t998;
t992 = mrSges(6,2) * t1035 - mrSges(6,3) * t1004;
t1097 = qJ(4) * (t1035 * t992 + t1055) - t1045 * t887 - t1041 * t888 + mrSges(4,1) * t916 - mrSges(4,2) * t917 - mrSges(6,1) * t908 + mrSges(5,3) * t913 - mrSges(5,1) * t914 + mrSges(6,2) * t907 - pkin(5) * t899 - pkin(9) * t883 - pkin(3) * t879 - pkin(4) * t881 + t1092 * t955 + (-qJ(4) * mrSges(5,2) - t1091) * t954 + (qJ(4) * mrSges(6,1) - t1090) * t1034 - t1089 * t1004 + (-qJ(4) * t976 + t1099) * t1003;
t1001 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t1043 + Ifges(3,2) * t1046) * qJD(1);
t1002 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t1043 + Ifges(3,4) * t1046) * qJD(1);
t1077 = -mrSges(4,1) * t1003 - mrSges(4,2) * t1004 - t976;
t1081 = -mrSges(4,3) - mrSges(5,2);
t988 = -mrSges(4,2) * t1035 - mrSges(4,3) * t1003;
t876 = m(4) * t916 + t1034 * mrSges(4,1) + t1035 * t988 + (mrSges(6,3) + t1081) * t955 + (t974 + t1077) * t1004 - t1056;
t990 = mrSges(4,1) * t1035 - mrSges(4,3) * t1004;
t886 = t1055 + t1081 * t954 + (-t990 + t992) * t1035 + (-mrSges(4,2) + mrSges(6,1)) * t1034 + t1077 * t1003 + m(4) * t917;
t871 = t1042 * t886 + t1082 * t876;
t985 = -t1046 * g(3) - t1072;
t1093 = mrSges(3,1) * t985 - mrSges(3,2) * t986 + Ifges(3,5) * t1013 + Ifges(3,6) * t1014 + Ifges(3,3) * qJDD(2) + pkin(2) * t871 + (t1001 * t1043 - t1002 * t1046) * qJD(1) + t1097;
t1080 = pkin(3) * t1035;
t1012 = (-mrSges(3,1) * t1046 + mrSges(3,2) * t1043) * qJD(1);
t1017 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1075;
t869 = m(3) * t985 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t1013 + qJD(2) * t1017 - t1012 * t1076 + t871;
t1016 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1076;
t1068 = -t1042 * t876 + t1082 * t886;
t870 = m(3) * t986 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t1014 - qJD(2) * t1016 + t1012 * t1075 + t1068;
t1067 = -t1043 * t869 + t1046 * t870;
t863 = m(2) * t1020 - mrSges(2,1) * t1048 - qJDD(1) * mrSges(2,2) + t1067;
t882 = t1041 * t895 + t1045 * t894;
t905 = -t954 * pkin(4) - t1004 * t1080 + t1054;
t880 = m(6) * t905 + t955 * mrSges(6,1) + t954 * mrSges(6,2) + t1003 * t987 + t1004 * t992 + t882;
t910 = (t1086 + t1080) * t1004 + t1059;
t877 = m(5) * t910 + t954 * mrSges(5,1) - t955 * mrSges(5,3) + t1003 * t993 - t1004 * t991 - t880;
t1052 = m(4) * t956 + t954 * mrSges(4,1) + t955 * mrSges(4,2) + t1003 * t988 + t1004 * t990 + t877;
t1051 = -m(3) * t1005 + t1014 * mrSges(3,1) - t1013 * mrSges(3,2) - t1016 * t1076 + t1017 * t1075 - t1052;
t873 = m(2) * t1019 + qJDD(1) * mrSges(2,1) - t1048 * mrSges(2,2) + t1051;
t1079 = t1044 * t863 + t1047 * t873;
t865 = t1043 * t870 + t1046 * t869;
t1070 = t1091 * t1003 - t1092 * t1004 + t1090 * t1035;
t1066 = -t1044 * t873 + t1047 * t863;
t1000 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t1043 + Ifges(3,6) * t1046) * qJD(1);
t859 = -t1045 * t888 + t1041 * t887 - qJ(5) * t1057 - mrSges(4,1) * t956 + mrSges(4,3) * t917 + mrSges(6,3) * t908 - mrSges(5,1) * t910 + mrSges(5,2) * t913 - mrSges(6,2) * t905 - pkin(3) * t877 + pkin(4) * t880 + pkin(9) * t882 - t1100 * t955 - t1098 * t954 + (-qJ(5) * t992 + t1099) * t1035 + (-mrSges(6,1) * qJ(5) + t1091) * t1034 + t1070 * t1004;
t1060 = mrSges(7,1) * t897 - mrSges(7,2) * t898 + Ifges(7,5) * t923 + Ifges(7,6) * t922 + Ifges(7,3) * t952 + t984 * t925 - t983 * t926;
t860 = mrSges(6,1) * t905 + mrSges(4,2) * t956 + mrSges(5,2) * t914 - mrSges(4,3) * t916 - mrSges(5,3) * t910 - mrSges(6,3) * t907 + pkin(5) * t882 - qJ(4) * t877 - qJ(5) * t881 + t1070 * t1003 + t1092 * t1034 + t1089 * t1035 + t1100 * t954 + t1101 * t955 + t1060;
t855 = -mrSges(3,1) * t1005 + mrSges(3,3) * t986 + Ifges(3,4) * t1013 + Ifges(3,2) * t1014 + Ifges(3,6) * qJDD(2) - pkin(2) * t1052 + pkin(8) * t1068 + qJD(2) * t1002 - t1000 * t1076 + t1042 * t860 + t1082 * t859;
t857 = mrSges(3,2) * t1005 - mrSges(3,3) * t985 + Ifges(3,1) * t1013 + Ifges(3,4) * t1014 + Ifges(3,5) * qJDD(2) - pkin(8) * t871 - qJD(2) * t1001 + t1000 * t1075 - t1042 * t859 + t1082 * t860;
t1062 = mrSges(2,1) * t1019 - mrSges(2,2) * t1020 + Ifges(2,3) * qJDD(1) + pkin(1) * t1051 + pkin(7) * t1067 + t1043 * t857 + t1046 * t855;
t858 = mrSges(2,1) * g(3) + mrSges(2,3) * t1020 + t1048 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t865 - t1093;
t853 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1019 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t1048 - pkin(7) * t865 - t1043 * t855 + t1046 * t857;
t1 = [-m(1) * g(1) + t1066; -m(1) * g(2) + t1079; (-m(1) - m(2)) * g(3) + t865; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1079 - t1044 * t858 + t1047 * t853; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1066 + t1044 * t853 + t1047 * t858; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1062; t1062; t1093; t1097; t879; t880; t1060;];
tauJB  = t1;
