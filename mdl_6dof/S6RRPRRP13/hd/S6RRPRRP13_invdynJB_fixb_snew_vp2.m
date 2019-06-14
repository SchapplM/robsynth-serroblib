% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-05-06 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP13_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:07:35
% EndTime: 2019-05-06 19:07:57
% DurationCPUTime: 11.06s
% Computational Cost: add. (166197->363), mult. (370509->438), div. (0->0), fcn. (267280->10), ass. (0->156)
t1084 = -2 * qJD(3);
t1083 = Ifges(3,1) + Ifges(4,2);
t1082 = Ifges(6,1) + Ifges(7,1);
t1072 = Ifges(3,4) + Ifges(4,6);
t1071 = Ifges(6,4) + Ifges(7,4);
t1070 = Ifges(3,5) - Ifges(4,4);
t1069 = Ifges(6,5) + Ifges(7,5);
t1081 = Ifges(3,2) + Ifges(4,3);
t1080 = Ifges(6,2) + Ifges(7,2);
t1068 = Ifges(3,6) - Ifges(4,5);
t1067 = Ifges(6,6) + Ifges(7,6);
t1079 = Ifges(3,3) + Ifges(4,1);
t1078 = Ifges(6,3) + Ifges(7,3);
t1017 = cos(pkin(6));
t1010 = qJD(1) * t1017 + qJD(2);
t1016 = sin(pkin(6));
t1020 = sin(qJ(2));
t1051 = t1016 * t1020;
t1041 = qJD(1) * t1051;
t1077 = (pkin(2) * t1010 + t1084) * t1041;
t1008 = t1010 ^ 2;
t1009 = qJDD(1) * t1017 + qJDD(2);
t1024 = cos(qJ(2));
t1050 = t1016 * t1024;
t1043 = qJD(1) * t1050;
t1049 = t1017 * t1020;
t1021 = sin(qJ(1));
t1025 = cos(qJ(1));
t1005 = g(1) * t1021 - g(2) * t1025;
t1026 = qJD(1) ^ 2;
t1065 = pkin(8) * t1016;
t989 = qJDD(1) * pkin(1) + t1026 * t1065 + t1005;
t1006 = -g(1) * t1025 - g(2) * t1021;
t1047 = qJDD(1) * t1016;
t990 = -pkin(1) * t1026 + pkin(8) * t1047 + t1006;
t954 = -g(3) * t1051 + t1024 * t990 + t1049 * t989;
t1054 = qJD(1) * t1016;
t991 = (-pkin(2) * t1024 - qJ(3) * t1020) * t1054;
t921 = t1008 * pkin(2) - qJ(3) * t1009 + t1010 * t1084 - t1043 * t991 - t954;
t1001 = qJD(4) + t1041;
t1018 = sin(qJ(5));
t1022 = cos(qJ(5));
t1019 = sin(qJ(4));
t1023 = cos(qJ(4));
t978 = t1010 * t1023 - t1019 * t1043;
t958 = t1001 * t1022 - t1018 * t978;
t959 = t1001 * t1018 + t1022 * t978;
t977 = -t1010 * t1019 - t1023 * t1043;
t976 = qJD(5) - t977;
t1060 = t1069 * t976 + t1071 * t958 + t1082 * t959;
t1061 = -t1067 * t976 - t1071 * t959 - t1080 * t958;
t1052 = t1016 ^ 2 * t1026;
t1042 = t1024 ^ 2 * t1052;
t1064 = t1017 * g(3);
t1075 = -pkin(2) - pkin(9);
t994 = pkin(3) * t1041 - pkin(9) * t1010;
t995 = (qJD(1) * qJD(2) * t1024 + qJDD(1) * t1020) * t1016;
t996 = -qJD(2) * t1041 + t1024 * t1047;
t911 = -pkin(3) * t1042 - t1064 - t995 * qJ(3) + t1075 * t996 + (-t989 + (-qJ(3) * t1010 * t1024 - t1020 * t994) * qJD(1)) * t1016 + t1077;
t1055 = g(3) * t1050 + t1020 * t990;
t1035 = -t1008 * qJ(3) + t1041 * t991 + qJDD(3) + t1055;
t913 = t995 * pkin(3) + t1075 * t1009 + (-pkin(3) * t1010 * t1054 - pkin(9) * t1020 * t1052 - t1017 * t989) * t1024 + t1035;
t906 = t1019 * t913 + t1023 * t911;
t956 = -pkin(4) * t977 - pkin(10) * t978;
t984 = qJDD(4) + t995;
t999 = t1001 ^ 2;
t900 = -pkin(4) * t999 + pkin(10) * t984 + t956 * t977 + t906;
t910 = t996 * pkin(3) - pkin(9) * t1042 + t1010 * t994 - t921;
t951 = -qJD(4) * t978 - t1009 * t1019 - t1023 * t996;
t952 = qJD(4) * t977 + t1009 * t1023 - t1019 * t996;
t903 = (-t1001 * t977 - t952) * pkin(10) + (t1001 * t978 - t951) * pkin(4) + t910;
t895 = -t1018 * t900 + t1022 * t903;
t918 = qJD(5) * t958 + t1018 * t984 + t1022 * t952;
t949 = qJDD(5) - t951;
t892 = -0.2e1 * qJD(6) * t959 + (t958 * t976 - t918) * qJ(6) + (t958 * t959 + t949) * pkin(5) + t895;
t937 = -mrSges(7,2) * t976 + mrSges(7,3) * t958;
t1046 = m(7) * t892 + mrSges(7,1) * t949 + t937 * t976;
t933 = -mrSges(7,1) * t958 + mrSges(7,2) * t959;
t889 = -t918 * mrSges(7,3) - t959 * t933 + t1046;
t896 = t1018 * t903 + t1022 * t900;
t917 = -qJD(5) * t959 - t1018 * t952 + t1022 * t984;
t939 = pkin(5) * t976 - qJ(6) * t959;
t957 = t958 ^ 2;
t894 = -pkin(5) * t957 + qJ(6) * t917 + 0.2e1 * qJD(6) * t958 - t939 * t976 + t896;
t1076 = mrSges(6,1) * t895 + mrSges(7,1) * t892 - mrSges(6,2) * t896 - mrSges(7,2) * t894 + pkin(5) * t889 - t1060 * t958 - t1061 * t959 + t1067 * t917 + t1069 * t918 + t1078 * t949;
t1074 = mrSges(3,1) - mrSges(4,2);
t1073 = -mrSges(6,2) - mrSges(7,2);
t1048 = t1017 * t1024;
t934 = -mrSges(6,1) * t958 + mrSges(6,2) * t959;
t938 = -mrSges(6,2) * t976 + mrSges(6,3) * t958;
t883 = m(6) * t895 + t949 * mrSges(6,1) + t976 * t938 + (-t933 - t934) * t959 + (-mrSges(6,3) - mrSges(7,3)) * t918 + t1046;
t1045 = m(7) * t894 + mrSges(7,3) * t917 + t933 * t958;
t940 = mrSges(7,1) * t976 - mrSges(7,3) * t959;
t1059 = -mrSges(6,1) * t976 + mrSges(6,3) * t959 - t940;
t886 = m(6) * t896 + t917 * mrSges(6,3) + t1059 * t976 + t1073 * t949 + t958 * t934 + t1045;
t1040 = -t1018 * t883 + t1022 * t886;
t955 = -mrSges(5,1) * t977 + mrSges(5,2) * t978;
t961 = mrSges(5,1) * t1001 - mrSges(5,3) * t978;
t877 = m(5) * t906 - mrSges(5,2) * t984 + mrSges(5,3) * t951 - t1001 * t961 + t955 * t977 + t1040;
t905 = -t1019 * t911 + t1023 * t913;
t899 = -pkin(4) * t984 - pkin(10) * t999 + t956 * t978 - t905;
t897 = -pkin(5) * t917 - qJ(6) * t957 + t939 * t959 + qJDD(6) + t899;
t1037 = -m(7) * t897 + mrSges(7,1) * t917 + t937 * t958;
t1027 = -m(6) * t899 + mrSges(6,1) * t917 + t1059 * t959 + t1073 * t918 + t938 * t958 + t1037;
t960 = -mrSges(5,2) * t1001 + mrSges(5,3) * t977;
t887 = m(5) * t905 + t984 * mrSges(5,1) - t952 * mrSges(5,3) + t1001 * t960 - t978 * t955 + t1027;
t1039 = -t1019 * t887 + t1023 * t877;
t968 = -t1016 * t989 - t1064;
t922 = -t996 * pkin(2) + (-t1010 * t1043 - t995) * qJ(3) + t968 + t1077;
t987 = -mrSges(4,1) * t1043 - mrSges(4,3) * t1010;
t1036 = m(4) * t922 - mrSges(4,3) * t995 + t1043 * t987 + t1039;
t985 = mrSges(3,1) * t1010 - mrSges(3,3) * t1041;
t986 = -mrSges(3,2) * t1010 + mrSges(3,3) * t1043;
t988 = mrSges(4,1) * t1041 + mrSges(4,2) * t1010;
t865 = m(3) * t968 + t995 * mrSges(3,2) - t1074 * t996 + (-t1024 * t986 + (t985 - t988) * t1020) * t1054 + t1036;
t869 = t1019 * t877 + t1023 * t887;
t1044 = t989 * t1048;
t931 = -t1009 * pkin(2) + t1035 - t1044;
t1034 = -m(4) * t931 - mrSges(4,1) * t995 - t869;
t953 = t1044 - t1055;
t992 = (mrSges(4,2) * t1024 - mrSges(4,3) * t1020) * t1054;
t993 = (-mrSges(3,1) * t1024 + mrSges(3,2) * t1020) * t1054;
t866 = m(3) * t953 - t995 * mrSges(3,3) + (t986 - t987) * t1010 + t1074 * t1009 + (-t992 - t993) * t1041 + t1034;
t881 = t1018 * t886 + t1022 * t883;
t1032 = -m(5) * t910 + t951 * mrSges(5,1) - mrSges(5,2) * t952 + t977 * t960 - t961 * t978 - t881;
t1028 = -m(4) * t921 + t1009 * mrSges(4,3) + t1010 * t988 + t992 * t1043 - t1032;
t875 = t993 * t1043 + t1028 - t1010 * t985 - t1009 * mrSges(3,2) + m(3) * t954 + (mrSges(3,3) + mrSges(4,1)) * t996;
t854 = -t1016 * t865 + t1048 * t866 + t1049 * t875;
t851 = m(2) * t1005 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1026 + t854;
t860 = -t1020 * t866 + t1024 * t875;
t858 = m(2) * t1006 - mrSges(2,1) * t1026 - qJDD(1) * mrSges(2,2) + t860;
t1063 = t1021 * t858 + t1025 * t851;
t1062 = -t1067 * t958 - t1069 * t959 - t1078 * t976;
t1058 = (t1020 * t1070 + t1024 * t1068) * t1054 + t1079 * t1010;
t1057 = (t1020 * t1072 + t1024 * t1081) * t1054 + t1068 * t1010;
t1056 = (t1020 * t1083 + t1024 * t1072) * t1054 + t1070 * t1010;
t853 = t1017 * t865 + t1050 * t866 + t1051 * t875;
t1038 = -t1021 * t851 + t1025 * t858;
t890 = t918 * mrSges(7,2) + t959 * t940 - t1037;
t871 = -mrSges(6,1) * t899 + mrSges(6,3) * t896 - mrSges(7,1) * t897 + mrSges(7,3) * t894 - pkin(5) * t890 + qJ(6) * t1045 + (-qJ(6) * t940 + t1060) * t976 + t1062 * t959 + (-mrSges(7,2) * qJ(6) + t1067) * t949 + t1071 * t918 + t1080 * t917;
t879 = mrSges(6,2) * t899 + mrSges(7,2) * t897 - mrSges(6,3) * t895 - mrSges(7,3) * t892 - qJ(6) * t889 + t1061 * t976 - t1062 * t958 + t1069 * t949 + t1071 * t917 + t1082 * t918;
t943 = Ifges(5,5) * t978 + Ifges(5,6) * t977 + Ifges(5,3) * t1001;
t944 = Ifges(5,4) * t978 + Ifges(5,2) * t977 + Ifges(5,6) * t1001;
t855 = mrSges(5,2) * t910 - mrSges(5,3) * t905 + Ifges(5,1) * t952 + Ifges(5,4) * t951 + Ifges(5,5) * t984 - pkin(10) * t881 - t1001 * t944 - t1018 * t871 + t1022 * t879 + t943 * t977;
t945 = Ifges(5,1) * t978 + Ifges(5,4) * t977 + Ifges(5,5) * t1001;
t861 = -mrSges(5,1) * t910 + mrSges(5,3) * t906 + Ifges(5,4) * t952 + Ifges(5,2) * t951 + Ifges(5,6) * t984 - pkin(4) * t881 + t1001 * t945 - t978 * t943 - t1076;
t868 = t1009 * mrSges(4,2) + t1010 * t987 + t1041 * t992 - t1034;
t845 = mrSges(3,1) * t953 - mrSges(3,2) * t954 + mrSges(4,2) * t931 - mrSges(4,3) * t921 + t1023 * t855 - t1019 * t861 - pkin(9) * t869 - pkin(2) * t868 + qJ(3) * t1028 + (mrSges(4,1) * qJ(3) + t1068) * t996 + t1070 * t995 + t1079 * t1009 + (t1020 * t1057 - t1024 * t1056) * t1054;
t867 = t996 * mrSges(4,2) - t1041 * t988 + t1036;
t847 = -mrSges(3,1) * t968 - mrSges(4,1) * t921 + mrSges(4,2) * t922 + mrSges(3,3) * t954 - pkin(2) * t867 - pkin(3) * t1032 - pkin(9) * t1039 + t1009 * t1068 + t1010 * t1056 - t1019 * t855 - t1023 * t861 - t1041 * t1058 + t1072 * t995 + t1081 * t996;
t1029 = mrSges(5,1) * t905 - mrSges(5,2) * t906 + Ifges(5,5) * t952 + Ifges(5,6) * t951 + Ifges(5,3) * t984 + pkin(4) * t1027 + pkin(10) * t1040 + t1018 * t879 + t1022 * t871 + t944 * t978 - t977 * t945;
t849 = mrSges(4,1) * t931 + mrSges(3,2) * t968 - mrSges(3,3) * t953 - mrSges(4,3) * t922 + pkin(3) * t869 - qJ(3) * t867 + t1009 * t1070 - t1010 * t1057 + t1043 * t1058 + t1072 * t996 + t1083 * t995 + t1029;
t1033 = mrSges(2,1) * t1005 - mrSges(2,2) * t1006 + Ifges(2,3) * qJDD(1) + pkin(1) * t854 + t1017 * t845 + t1050 * t847 + t1051 * t849 + t1065 * t860;
t843 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1005 + Ifges(2,5) * qJDD(1) - t1026 * Ifges(2,6) - t1020 * t847 + t1024 * t849 + (-t1016 * t853 - t1017 * t854) * pkin(8);
t842 = mrSges(2,1) * g(3) + mrSges(2,3) * t1006 + t1026 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t853 - t1016 * t845 + (pkin(8) * t860 + t1020 * t849 + t1024 * t847) * t1017;
t1 = [-m(1) * g(1) + t1038; -m(1) * g(2) + t1063; (-m(1) - m(2)) * g(3) + t853; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1063 - t1021 * t842 + t1025 * t843; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1038 + t1021 * t843 + t1025 * t842; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1033; t1033; t845; t868; t1029; t1076; t890;];
tauJB  = t1;
