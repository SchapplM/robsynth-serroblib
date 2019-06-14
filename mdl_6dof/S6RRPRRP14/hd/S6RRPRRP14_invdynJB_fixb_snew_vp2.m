% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP14
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
% Datum: 2019-05-06 19:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP14_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP14_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:18:31
% EndTime: 2019-05-06 19:18:48
% DurationCPUTime: 11.00s
% Computational Cost: add. (162417->361), mult. (361572->438), div. (0->0), fcn. (260092->10), ass. (0->155)
t1081 = -2 * qJD(3);
t1080 = Ifges(3,1) + Ifges(4,2);
t1079 = Ifges(6,1) + Ifges(7,1);
t1068 = Ifges(3,4) + Ifges(4,6);
t1067 = Ifges(6,4) - Ifges(7,5);
t1066 = Ifges(3,5) - Ifges(4,4);
t1065 = -Ifges(6,5) - Ifges(7,4);
t1078 = Ifges(3,2) + Ifges(4,3);
t1077 = Ifges(6,2) + Ifges(7,3);
t1064 = Ifges(3,6) - Ifges(4,5);
t1063 = Ifges(6,6) - Ifges(7,6);
t1076 = Ifges(3,3) + Ifges(4,1);
t1075 = -Ifges(6,3) - Ifges(7,2);
t1015 = cos(pkin(6));
t1008 = qJD(1) * t1015 + qJD(2);
t1014 = sin(pkin(6));
t1018 = sin(qJ(2));
t1047 = t1014 * t1018;
t1040 = qJD(1) * t1047;
t1074 = (pkin(2) * t1008 + t1081) * t1040;
t1006 = t1008 ^ 2;
t1007 = qJDD(1) * t1015 + qJDD(2);
t1021 = cos(qJ(2));
t1046 = t1014 * t1021;
t1039 = qJD(1) * t1046;
t1045 = t1015 * t1018;
t1019 = sin(qJ(1));
t1022 = cos(qJ(1));
t1003 = t1019 * g(1) - g(2) * t1022;
t1023 = qJD(1) ^ 2;
t1061 = pkin(8) * t1014;
t986 = qJDD(1) * pkin(1) + t1023 * t1061 + t1003;
t1004 = -g(1) * t1022 - g(2) * t1019;
t1043 = qJDD(1) * t1014;
t987 = -pkin(1) * t1023 + pkin(8) * t1043 + t1004;
t950 = -g(3) * t1047 + t1021 * t987 + t986 * t1045;
t1050 = qJD(1) * t1014;
t988 = (-pkin(2) * t1021 - qJ(3) * t1018) * t1050;
t915 = t1006 * pkin(2) - t1007 * qJ(3) + t1008 * t1081 - t988 * t1039 - t950;
t1016 = sin(qJ(5));
t1071 = cos(qJ(5));
t1017 = sin(qJ(4));
t1020 = cos(qJ(4));
t1048 = t1014 ^ 2 * t1023;
t1038 = t1021 ^ 2 * t1048;
t1060 = t1015 * g(3);
t1072 = -pkin(2) - pkin(9);
t991 = pkin(3) * t1040 - pkin(9) * t1008;
t992 = (qJD(1) * qJD(2) * t1021 + qJDD(1) * t1018) * t1014;
t993 = -qJD(2) * t1040 + t1021 * t1043;
t907 = -pkin(3) * t1038 - t1060 - t992 * qJ(3) + t1072 * t993 + (-t986 + (-qJ(3) * t1008 * t1021 - t1018 * t991) * qJD(1)) * t1014 + t1074;
t1051 = g(3) * t1046 + t1018 * t987;
t1032 = -t1006 * qJ(3) + t988 * t1040 + qJDD(3) + t1051;
t909 = t992 * pkin(3) + t1072 * t1007 + (-pkin(3) * t1008 * t1050 - pkin(9) * t1018 * t1048 - t1015 * t986) * t1021 + t1032;
t902 = t1017 * t909 + t1020 * t907;
t974 = -t1008 * t1017 - t1020 * t1039;
t975 = t1008 * t1020 - t1017 * t1039;
t952 = -pkin(4) * t974 - pkin(10) * t975;
t981 = qJDD(4) + t992;
t999 = qJD(4) + t1040;
t997 = t999 ^ 2;
t897 = -pkin(4) * t997 + pkin(10) * t981 + t952 * t974 + t902;
t906 = t993 * pkin(3) - pkin(9) * t1038 + t1008 * t991 - t915;
t947 = -qJD(4) * t975 - t1007 * t1017 - t1020 * t993;
t948 = qJD(4) * t974 + t1007 * t1020 - t1017 * t993;
t899 = (-t974 * t999 - t948) * pkin(10) + (t975 * t999 - t947) * pkin(4) + t906;
t894 = t1016 * t899 + t1071 * t897;
t953 = t1016 * t975 - t1071 * t999;
t954 = t1016 * t999 + t1071 * t975;
t928 = pkin(5) * t953 - qJ(6) * t954;
t945 = qJDD(5) - t947;
t972 = qJD(5) - t974;
t971 = t972 ^ 2;
t890 = -pkin(5) * t971 + qJ(6) * t945 + 0.2e1 * qJD(6) * t972 - t928 * t953 + t894;
t936 = -mrSges(7,1) * t972 + mrSges(7,2) * t954;
t1042 = m(7) * t890 + t945 * mrSges(7,3) + t972 * t936;
t1056 = t1065 * t972 + t1067 * t953 - t1079 * t954;
t1057 = -t1063 * t972 - t1067 * t954 + t1077 * t953;
t893 = -t1016 * t897 + t1071 * t899;
t891 = -t945 * pkin(5) - t971 * qJ(6) + t954 * t928 + qJDD(6) - t893;
t933 = -mrSges(7,2) * t953 + mrSges(7,3) * t972;
t1034 = -m(7) * t891 + t945 * mrSges(7,1) + t972 * t933;
t913 = -t953 * qJD(5) + t1016 * t981 + t1071 * t948;
t929 = mrSges(7,1) * t953 - mrSges(7,3) * t954;
t887 = t913 * mrSges(7,2) + t954 * t929 - t1034;
t912 = qJD(5) * t954 + t1016 * t948 - t1071 * t981;
t1073 = -t1056 * t953 - t1057 * t954 - t1075 * t945 - t1063 * t912 - t1065 * t913 + mrSges(6,1) * t893 - mrSges(7,1) * t891 - mrSges(6,2) * t894 + mrSges(7,3) * t890 - pkin(5) * t887 + qJ(6) * (-t912 * mrSges(7,2) - t953 * t929 + t1042);
t1070 = mrSges(3,1) - mrSges(4,2);
t1069 = -mrSges(6,3) - mrSges(7,2);
t1044 = t1015 * t1021;
t1055 = -mrSges(6,1) * t953 - mrSges(6,2) * t954 - t929;
t935 = mrSges(6,1) * t972 - mrSges(6,3) * t954;
t882 = m(6) * t894 - t945 * mrSges(6,2) + t1055 * t953 + t1069 * t912 - t972 * t935 + t1042;
t934 = -mrSges(6,2) * t972 - mrSges(6,3) * t953;
t884 = m(6) * t893 + t945 * mrSges(6,1) + t1055 * t954 + t1069 * t913 + t972 * t934 + t1034;
t1037 = -t1016 * t884 + t1071 * t882;
t951 = -mrSges(5,1) * t974 + mrSges(5,2) * t975;
t956 = mrSges(5,1) * t999 - mrSges(5,3) * t975;
t872 = m(5) * t902 - mrSges(5,2) * t981 + mrSges(5,3) * t947 + t951 * t974 - t956 * t999 + t1037;
t901 = -t1017 * t907 + t1020 * t909;
t896 = -t981 * pkin(4) - t997 * pkin(10) + t975 * t952 - t901;
t892 = -0.2e1 * qJD(6) * t954 + (t953 * t972 - t913) * qJ(6) + (t954 * t972 + t912) * pkin(5) + t896;
t888 = m(7) * t892 + mrSges(7,1) * t912 - t913 * mrSges(7,3) + t933 * t953 - t954 * t936;
t1024 = -m(6) * t896 - t912 * mrSges(6,1) - mrSges(6,2) * t913 - t953 * t934 - t935 * t954 - t888;
t955 = -mrSges(5,2) * t999 + mrSges(5,3) * t974;
t879 = m(5) * t901 + mrSges(5,1) * t981 - mrSges(5,3) * t948 - t951 * t975 + t955 * t999 + t1024;
t1036 = -t1017 * t879 + t1020 * t872;
t963 = -t1014 * t986 - t1060;
t916 = -t993 * pkin(2) + (-t1008 * t1039 - t992) * qJ(3) + t963 + t1074;
t984 = -mrSges(4,1) * t1039 - mrSges(4,3) * t1008;
t1033 = m(4) * t916 - t992 * mrSges(4,3) + t984 * t1039 + t1036;
t982 = mrSges(3,1) * t1008 - mrSges(3,3) * t1040;
t983 = -mrSges(3,2) * t1008 + mrSges(3,3) * t1039;
t985 = mrSges(4,1) * t1040 + mrSges(4,2) * t1008;
t862 = m(3) * t963 + t992 * mrSges(3,2) - t1070 * t993 + (-t1021 * t983 + (t982 - t985) * t1018) * t1050 + t1033;
t866 = t1017 * t872 + t1020 * t879;
t1041 = t986 * t1044;
t925 = -t1007 * pkin(2) + t1032 - t1041;
t1031 = -m(4) * t925 - t992 * mrSges(4,1) - t866;
t949 = t1041 - t1051;
t989 = (mrSges(4,2) * t1021 - mrSges(4,3) * t1018) * t1050;
t990 = (-mrSges(3,1) * t1021 + mrSges(3,2) * t1018) * t1050;
t863 = m(3) * t949 - t992 * mrSges(3,3) + (t983 - t984) * t1008 + t1070 * t1007 + (-t989 - t990) * t1040 + t1031;
t878 = t1016 * t882 + t1071 * t884;
t1029 = -m(5) * t906 + t947 * mrSges(5,1) - t948 * mrSges(5,2) + t974 * t955 - t975 * t956 - t878;
t1025 = -m(4) * t915 + t1007 * mrSges(4,3) + t1008 * t985 + t989 * t1039 - t1029;
t870 = t1025 - t1008 * t982 - t1007 * mrSges(3,2) + m(3) * t950 + (mrSges(3,3) + mrSges(4,1)) * t993 + t990 * t1039;
t851 = -t1014 * t862 + t863 * t1044 + t870 * t1045;
t848 = m(2) * t1003 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1023 + t851;
t857 = -t1018 * t863 + t1021 * t870;
t855 = m(2) * t1004 - mrSges(2,1) * t1023 - qJDD(1) * mrSges(2,2) + t857;
t1059 = t1019 * t855 + t1022 * t848;
t1058 = t1063 * t953 + t1065 * t954 + t1075 * t972;
t1054 = (t1018 * t1066 + t1021 * t1064) * t1050 + t1076 * t1008;
t1053 = (t1018 * t1068 + t1021 * t1078) * t1050 + t1064 * t1008;
t1052 = (t1018 * t1080 + t1021 * t1068) * t1050 + t1066 * t1008;
t850 = t1015 * t862 + t863 * t1046 + t870 * t1047;
t1035 = -t1019 * t848 + t1022 * t855;
t874 = -mrSges(6,1) * t896 - mrSges(7,1) * t892 + mrSges(7,2) * t890 + mrSges(6,3) * t894 - pkin(5) * t888 - t1056 * t972 + t1058 * t954 + t1063 * t945 + t1067 * t913 - t1077 * t912;
t876 = mrSges(6,2) * t896 + mrSges(7,2) * t891 - mrSges(6,3) * t893 - mrSges(7,3) * t892 - qJ(6) * t888 + t1057 * t972 + t1058 * t953 - t1065 * t945 - t1067 * t912 + t1079 * t913;
t939 = Ifges(5,5) * t975 + Ifges(5,6) * t974 + Ifges(5,3) * t999;
t940 = Ifges(5,4) * t975 + Ifges(5,2) * t974 + Ifges(5,6) * t999;
t852 = mrSges(5,2) * t906 - mrSges(5,3) * t901 + Ifges(5,1) * t948 + Ifges(5,4) * t947 + Ifges(5,5) * t981 - pkin(10) * t878 - t1016 * t874 + t1071 * t876 + t974 * t939 - t999 * t940;
t941 = Ifges(5,1) * t975 + Ifges(5,4) * t974 + Ifges(5,5) * t999;
t858 = -mrSges(5,1) * t906 + mrSges(5,3) * t902 + Ifges(5,4) * t948 + Ifges(5,2) * t947 + Ifges(5,6) * t981 - pkin(4) * t878 - t975 * t939 + t999 * t941 - t1073;
t865 = t1007 * mrSges(4,2) + t1008 * t984 + t989 * t1040 - t1031;
t842 = mrSges(3,1) * t949 - mrSges(3,2) * t950 + mrSges(4,2) * t925 - mrSges(4,3) * t915 + t1020 * t852 - t1017 * t858 - pkin(9) * t866 - pkin(2) * t865 + qJ(3) * t1025 + (mrSges(4,1) * qJ(3) + t1064) * t993 + t1066 * t992 + t1076 * t1007 + (t1053 * t1018 - t1052 * t1021) * t1050;
t864 = t993 * mrSges(4,2) - t985 * t1040 + t1033;
t844 = -mrSges(3,1) * t963 - mrSges(4,1) * t915 + mrSges(4,2) * t916 + mrSges(3,3) * t950 - pkin(2) * t864 - pkin(3) * t1029 - pkin(9) * t1036 + t1064 * t1007 + t1052 * t1008 - t1017 * t852 - t1020 * t858 - t1054 * t1040 + t1068 * t992 + t1078 * t993;
t1027 = mrSges(5,1) * t901 - mrSges(5,2) * t902 + Ifges(5,5) * t948 + Ifges(5,6) * t947 + Ifges(5,3) * t981 + pkin(4) * t1024 + pkin(10) * t1037 + t1016 * t876 + t1071 * t874 + t975 * t940 - t974 * t941;
t846 = mrSges(4,1) * t925 + mrSges(3,2) * t963 - mrSges(3,3) * t949 - mrSges(4,3) * t916 + pkin(3) * t866 - qJ(3) * t864 + t1066 * t1007 - t1053 * t1008 + t1054 * t1039 + t1068 * t993 + t1080 * t992 + t1027;
t1030 = mrSges(2,1) * t1003 - mrSges(2,2) * t1004 + Ifges(2,3) * qJDD(1) + pkin(1) * t851 + t1015 * t842 + t844 * t1046 + t846 * t1047 + t857 * t1061;
t840 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1003 + Ifges(2,5) * qJDD(1) - t1023 * Ifges(2,6) - t1018 * t844 + t1021 * t846 + (-t1014 * t850 - t1015 * t851) * pkin(8);
t839 = mrSges(2,1) * g(3) + mrSges(2,3) * t1004 + t1023 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t850 - t1014 * t842 + (pkin(8) * t857 + t1018 * t846 + t1021 * t844) * t1015;
t1 = [-m(1) * g(1) + t1035; -m(1) * g(2) + t1059; (-m(1) - m(2)) * g(3) + t850; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1059 - t1019 * t839 + t1022 * t840; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1035 + t1019 * t840 + t1022 * t839; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1030; t1030; t842; t865; t1027; t1073; t887;];
tauJB  = t1;
