% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:15:26
% EndTime: 2019-05-06 16:15:50
% DurationCPUTime: 21.59s
% Computational Cost: add. (343275->385), mult. (788428->481), div. (0->0), fcn. (584798->12), ass. (0->164)
t1086 = -2 * qJD(3);
t1085 = Ifges(3,1) + Ifges(4,2);
t1078 = Ifges(3,4) + Ifges(4,6);
t1077 = Ifges(3,5) - Ifges(4,4);
t1084 = Ifges(3,2) + Ifges(4,3);
t1076 = Ifges(3,6) - Ifges(4,5);
t1083 = Ifges(3,3) + Ifges(4,1);
t1032 = cos(pkin(6));
t1023 = t1032 * qJD(1) + qJD(2);
t1030 = sin(pkin(6));
t1035 = sin(qJ(2));
t1065 = t1030 * t1035;
t1058 = qJD(1) * t1065;
t1082 = (pkin(2) * t1023 + t1086) * t1058;
t1039 = cos(qJ(2));
t1068 = qJD(1) * t1030;
t1004 = (-pkin(2) * t1039 - qJ(3) * t1035) * t1068;
t1021 = t1023 ^ 2;
t1022 = t1032 * qJDD(1) + qJDD(2);
t1064 = t1030 * t1039;
t1060 = qJD(1) * t1064;
t1036 = sin(qJ(1));
t1040 = cos(qJ(1));
t1018 = t1036 * g(1) - t1040 * g(2);
t1041 = qJD(1) ^ 2;
t1075 = pkin(8) * t1030;
t1002 = qJDD(1) * pkin(1) + t1041 * t1075 + t1018;
t1019 = -t1040 * g(1) - t1036 * g(2);
t1061 = qJDD(1) * t1030;
t1003 = -t1041 * pkin(1) + pkin(8) * t1061 + t1019;
t1063 = t1032 * t1035;
t965 = -g(3) * t1065 + t1002 * t1063 + t1039 * t1003;
t944 = t1021 * pkin(2) - t1022 * qJ(3) - t1004 * t1060 + t1023 * t1086 - t965;
t1081 = 2 * qJD(5);
t1080 = -pkin(2) - pkin(9);
t1079 = mrSges(3,1) - mrSges(4,2);
t1074 = t1032 * g(3);
t1062 = t1032 * t1039;
t1001 = mrSges(4,1) * t1058 + t1023 * mrSges(4,2);
t1008 = (qJD(1) * qJD(2) * t1039 + qJDD(1) * t1035) * t1030;
t1009 = -qJD(2) * t1058 + t1039 * t1061;
t1000 = -mrSges(4,1) * t1060 - t1023 * mrSges(4,3);
t1034 = sin(qJ(4));
t1038 = cos(qJ(4));
t1014 = qJD(4) + t1058;
t1029 = sin(pkin(11));
t1031 = cos(pkin(11));
t1033 = sin(qJ(6));
t1037 = cos(qJ(6));
t1012 = t1014 ^ 2;
t1007 = pkin(3) * t1058 - t1023 * pkin(9);
t1066 = t1030 ^ 2 * t1041;
t1059 = t1039 ^ 2 * t1066;
t926 = -pkin(3) * t1059 - t1074 - t1008 * qJ(3) + t1080 * t1009 + (-t1002 + (-qJ(3) * t1023 * t1039 - t1007 * t1035) * qJD(1)) * t1030 + t1082;
t1069 = g(3) * t1064 + t1035 * t1003;
t1050 = -t1021 * qJ(3) + t1004 * t1058 + qJDD(3) + t1069;
t929 = t1008 * pkin(3) + t1080 * t1022 + (-pkin(3) * t1023 * t1068 - pkin(9) * t1035 * t1066 - t1002 * t1032) * t1039 + t1050;
t911 = -t1034 * t926 + t1038 * t929;
t989 = -t1034 * t1023 - t1038 * t1060;
t963 = t989 * qJD(4) - t1034 * t1009 + t1038 * t1022;
t990 = t1038 * t1023 - t1034 * t1060;
t997 = qJDD(4) + t1008;
t907 = (t1014 * t989 - t963) * qJ(5) + (t989 * t990 + t997) * pkin(4) + t911;
t912 = t1034 * t929 + t1038 * t926;
t962 = -t990 * qJD(4) - t1038 * t1009 - t1034 * t1022;
t972 = t1014 * pkin(4) - t990 * qJ(5);
t988 = t989 ^ 2;
t909 = -t988 * pkin(4) + t962 * qJ(5) - t1014 * t972 + t912;
t968 = -t1029 * t990 + t1031 * t989;
t904 = t1029 * t907 + t1031 * t909 + t1081 * t968;
t969 = t1029 * t989 + t1031 * t990;
t947 = -t968 * pkin(5) - t969 * pkin(10);
t901 = -t1012 * pkin(5) + t997 * pkin(10) + t968 * t947 + t904;
t925 = t1009 * pkin(3) - pkin(9) * t1059 + t1023 * t1007 - t944;
t914 = -t962 * pkin(4) - t988 * qJ(5) + t990 * t972 + qJDD(5) + t925;
t937 = -t1029 * t963 + t1031 * t962;
t938 = t1029 * t962 + t1031 * t963;
t905 = (t1014 * t969 - t937) * pkin(5) + (-t1014 * t968 - t938) * pkin(10) + t914;
t898 = -t1033 * t901 + t1037 * t905;
t950 = t1037 * t1014 - t1033 * t969;
t917 = t950 * qJD(6) + t1033 * t997 + t1037 * t938;
t951 = t1033 * t1014 + t1037 * t969;
t930 = -t950 * mrSges(7,1) + t951 * mrSges(7,2);
t967 = qJD(6) - t968;
t931 = -t967 * mrSges(7,2) + t950 * mrSges(7,3);
t936 = qJDD(6) - t937;
t895 = m(7) * t898 + t936 * mrSges(7,1) - t917 * mrSges(7,3) - t951 * t930 + t967 * t931;
t899 = t1033 * t905 + t1037 * t901;
t916 = -t951 * qJD(6) - t1033 * t938 + t1037 * t997;
t932 = t967 * mrSges(7,1) - t951 * mrSges(7,3);
t896 = m(7) * t899 - t936 * mrSges(7,2) + t916 * mrSges(7,3) + t950 * t930 - t967 * t932;
t1055 = -t1033 * t895 + t1037 * t896;
t946 = -t968 * mrSges(6,1) + t969 * mrSges(6,2);
t953 = t1014 * mrSges(6,1) - t969 * mrSges(6,3);
t882 = m(6) * t904 - t997 * mrSges(6,2) + t937 * mrSges(6,3) - t1014 * t953 + t968 * t946 + t1055;
t1052 = t1029 * t909 - t1031 * t907;
t900 = -t997 * pkin(5) - t1012 * pkin(10) + (t1081 + t947) * t969 + t1052;
t1048 = -m(7) * t900 + t916 * mrSges(7,1) - t917 * mrSges(7,2) + t950 * t931 - t951 * t932;
t903 = -0.2e1 * qJD(5) * t969 - t1052;
t952 = -t1014 * mrSges(6,2) + t968 * mrSges(6,3);
t891 = m(6) * t903 + t997 * mrSges(6,1) - t938 * mrSges(6,3) + t1014 * t952 - t969 * t946 + t1048;
t875 = t1029 * t882 + t1031 * t891;
t970 = -t989 * mrSges(5,1) + t990 * mrSges(5,2);
t971 = -t1014 * mrSges(5,2) + t989 * mrSges(5,3);
t872 = m(5) * t911 + t997 * mrSges(5,1) - t963 * mrSges(5,3) + t1014 * t971 - t990 * t970 + t875;
t1056 = -t1029 * t891 + t1031 * t882;
t973 = t1014 * mrSges(5,1) - t990 * mrSges(5,3);
t873 = m(5) * t912 - t997 * mrSges(5,2) + t962 * mrSges(5,3) - t1014 * t973 + t989 * t970 + t1056;
t1054 = -t1034 * t872 + t1038 * t873;
t980 = -t1030 * t1002 - t1074;
t945 = -t1009 * pkin(2) + (-t1023 * t1060 - t1008) * qJ(3) + t980 + t1082;
t1051 = m(4) * t945 - t1008 * mrSges(4,3) + t1000 * t1060 + t1054;
t998 = t1023 * mrSges(3,1) - mrSges(3,3) * t1058;
t999 = -t1023 * mrSges(3,2) + mrSges(3,3) * t1060;
t864 = m(3) * t980 + t1008 * mrSges(3,2) - t1079 * t1009 + (-t1039 * t999 + (-t1001 + t998) * t1035) * t1068 + t1051;
t1005 = (mrSges(4,2) * t1039 - mrSges(4,3) * t1035) * t1068;
t1006 = (-mrSges(3,1) * t1039 + mrSges(3,2) * t1035) * t1068;
t868 = t1034 * t873 + t1038 * t872;
t1057 = t1002 * t1062;
t948 = -t1022 * pkin(2) + t1050 - t1057;
t1049 = -m(4) * t948 - t1008 * mrSges(4,1) - t868;
t964 = t1057 - t1069;
t865 = m(3) * t964 - t1008 * mrSges(3,3) + (-t1000 + t999) * t1023 + t1079 * t1022 + (-t1005 - t1006) * t1058 + t1049;
t885 = t1033 * t896 + t1037 * t895;
t883 = m(6) * t914 - t937 * mrSges(6,1) + t938 * mrSges(6,2) - t968 * t952 + t969 * t953 + t885;
t1044 = -m(5) * t925 + t962 * mrSges(5,1) - t963 * mrSges(5,2) + t989 * t971 - t990 * t973 - t883;
t1043 = -m(4) * t944 + t1022 * mrSges(4,3) + t1023 * t1001 + t1005 * t1060 - t1044;
t879 = t1043 + t1006 * t1060 + (mrSges(3,3) + mrSges(4,1)) * t1009 - t1022 * mrSges(3,2) - t1023 * t998 + m(3) * t965;
t853 = -t1030 * t864 + t865 * t1062 + t879 * t1063;
t850 = m(2) * t1018 + qJDD(1) * mrSges(2,1) - t1041 * mrSges(2,2) + t853;
t860 = -t1035 * t865 + t1039 * t879;
t858 = m(2) * t1019 - t1041 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t860;
t1073 = t1036 * t858 + t1040 * t850;
t1072 = (t1035 * t1077 + t1039 * t1076) * t1068 + t1083 * t1023;
t1071 = (t1035 * t1078 + t1039 * t1084) * t1068 + t1076 * t1023;
t1070 = (t1035 * t1085 + t1039 * t1078) * t1068 + t1077 * t1023;
t852 = t1032 * t864 + t865 * t1064 + t879 * t1065;
t1053 = -t1036 * t850 + t1040 * t858;
t918 = Ifges(7,5) * t951 + Ifges(7,6) * t950 + Ifges(7,3) * t967;
t920 = Ifges(7,1) * t951 + Ifges(7,4) * t950 + Ifges(7,5) * t967;
t888 = -mrSges(7,1) * t900 + mrSges(7,3) * t899 + Ifges(7,4) * t917 + Ifges(7,2) * t916 + Ifges(7,6) * t936 - t951 * t918 + t967 * t920;
t919 = Ifges(7,4) * t951 + Ifges(7,2) * t950 + Ifges(7,6) * t967;
t889 = mrSges(7,2) * t900 - mrSges(7,3) * t898 + Ifges(7,1) * t917 + Ifges(7,4) * t916 + Ifges(7,5) * t936 + t950 * t918 - t967 * t919;
t939 = Ifges(6,5) * t969 + Ifges(6,6) * t968 + Ifges(6,3) * t1014;
t940 = Ifges(6,4) * t969 + Ifges(6,2) * t968 + Ifges(6,6) * t1014;
t869 = mrSges(6,2) * t914 - mrSges(6,3) * t903 + Ifges(6,1) * t938 + Ifges(6,4) * t937 + Ifges(6,5) * t997 - pkin(10) * t885 - t1014 * t940 - t1033 * t888 + t1037 * t889 + t968 * t939;
t1045 = mrSges(7,1) * t898 - mrSges(7,2) * t899 + Ifges(7,5) * t917 + Ifges(7,6) * t916 + Ifges(7,3) * t936 + t951 * t919 - t950 * t920;
t941 = Ifges(6,1) * t969 + Ifges(6,4) * t968 + Ifges(6,5) * t1014;
t870 = -mrSges(6,1) * t914 + mrSges(6,3) * t904 + Ifges(6,4) * t938 + Ifges(6,2) * t937 + Ifges(6,6) * t997 - pkin(5) * t885 + t1014 * t941 - t969 * t939 - t1045;
t954 = Ifges(5,5) * t990 + Ifges(5,6) * t989 + Ifges(5,3) * t1014;
t956 = Ifges(5,1) * t990 + Ifges(5,4) * t989 + Ifges(5,5) * t1014;
t854 = -mrSges(5,1) * t925 + mrSges(5,3) * t912 + Ifges(5,4) * t963 + Ifges(5,2) * t962 + Ifges(5,6) * t997 - pkin(4) * t883 + qJ(5) * t1056 + t1014 * t956 + t1029 * t869 + t1031 * t870 - t990 * t954;
t955 = Ifges(5,4) * t990 + Ifges(5,2) * t989 + Ifges(5,6) * t1014;
t855 = mrSges(5,2) * t925 - mrSges(5,3) * t911 + Ifges(5,1) * t963 + Ifges(5,4) * t962 + Ifges(5,5) * t997 - qJ(5) * t875 - t1014 * t955 - t1029 * t870 + t1031 * t869 + t989 * t954;
t867 = t1022 * mrSges(4,2) + t1023 * t1000 + t1005 * t1058 - t1049;
t844 = mrSges(3,1) * t964 - mrSges(3,2) * t965 + mrSges(4,2) * t948 - mrSges(4,3) * t944 + t1038 * t855 - t1034 * t854 - pkin(9) * t868 - pkin(2) * t867 + qJ(3) * t1043 + t1083 * t1022 + (qJ(3) * mrSges(4,1) + t1076) * t1009 + t1077 * t1008 + (t1035 * t1071 - t1039 * t1070) * t1068;
t866 = t1009 * mrSges(4,2) - t1001 * t1058 + t1051;
t846 = -mrSges(3,1) * t980 - mrSges(4,1) * t944 + mrSges(4,2) * t945 + mrSges(3,3) * t965 - pkin(2) * t866 - pkin(3) * t1044 - pkin(9) * t1054 + t1078 * t1008 + t1009 * t1084 + t1076 * t1022 + t1070 * t1023 - t1034 * t855 - t1038 * t854 - t1072 * t1058;
t1042 = mrSges(5,1) * t911 + mrSges(6,1) * t903 - mrSges(5,2) * t912 - mrSges(6,2) * t904 + pkin(4) * t875 + pkin(5) * t1048 + pkin(10) * t1055 + t1033 * t889 + t1037 * t888 + t969 * t940 - t968 * t941 - t989 * t956 + Ifges(6,6) * t937 + Ifges(6,5) * t938 + t990 * t955 + Ifges(5,6) * t962 + Ifges(5,5) * t963 + (Ifges(6,3) + Ifges(5,3)) * t997;
t848 = mrSges(4,1) * t948 + mrSges(3,2) * t980 - mrSges(3,3) * t964 - mrSges(4,3) * t945 + pkin(3) * t868 - qJ(3) * t866 + t1008 * t1085 + t1078 * t1009 + t1077 * t1022 - t1071 * t1023 + t1072 * t1060 + t1042;
t1047 = mrSges(2,1) * t1018 - mrSges(2,2) * t1019 + Ifges(2,3) * qJDD(1) + pkin(1) * t853 + t1032 * t844 + t846 * t1064 + t848 * t1065 + t860 * t1075;
t842 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1018 + Ifges(2,5) * qJDD(1) - t1041 * Ifges(2,6) - t1035 * t846 + t1039 * t848 + (-t1030 * t852 - t1032 * t853) * pkin(8);
t841 = mrSges(2,1) * g(3) + mrSges(2,3) * t1019 + t1041 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t852 - t1030 * t844 + (pkin(8) * t860 + t1035 * t848 + t1039 * t846) * t1032;
t1 = [-m(1) * g(1) + t1053; -m(1) * g(2) + t1073; (-m(1) - m(2)) * g(3) + t852; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1073 - t1036 * t841 + t1040 * t842; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1053 + t1036 * t842 + t1040 * t841; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1047; t1047; t844; t867; t1042; t883; t1045;];
tauJB  = t1;
