% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR13
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
% Datum: 2019-05-06 16:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:33:33
% EndTime: 2019-05-06 16:33:56
% DurationCPUTime: 21.84s
% Computational Cost: add. (355311->384), mult. (805690->481), div. (0->0), fcn. (593656->12), ass. (0->162)
t1075 = -2 * qJD(3);
t1074 = Ifges(3,1) + Ifges(4,2);
t1068 = Ifges(3,4) + Ifges(4,6);
t1067 = Ifges(3,5) - Ifges(4,4);
t1073 = Ifges(3,2) + Ifges(4,3);
t1066 = Ifges(3,6) - Ifges(4,5);
t1072 = Ifges(3,3) + Ifges(4,1);
t1023 = cos(pkin(6));
t1014 = qJD(1) * t1023 + qJD(2);
t1021 = sin(pkin(6));
t1026 = sin(qJ(2));
t1055 = t1021 * t1026;
t1049 = qJD(1) * t1055;
t1071 = (pkin(2) * t1014 + t1075) * t1049;
t1012 = t1014 ^ 2;
t1013 = qJDD(1) * t1023 + qJDD(2);
t1030 = cos(qJ(2));
t1054 = t1021 * t1030;
t1048 = qJD(1) * t1054;
t1053 = t1023 * t1026;
t1027 = sin(qJ(1));
t1031 = cos(qJ(1));
t1009 = t1027 * g(1) - g(2) * t1031;
t1032 = qJD(1) ^ 2;
t1065 = pkin(8) * t1021;
t993 = qJDD(1) * pkin(1) + t1032 * t1065 + t1009;
t1010 = -g(1) * t1031 - g(2) * t1027;
t1051 = qJDD(1) * t1021;
t994 = -pkin(1) * t1032 + pkin(8) * t1051 + t1010;
t955 = -g(3) * t1055 + t1030 * t994 + t993 * t1053;
t1058 = qJD(1) * t1021;
t995 = (-pkin(2) * t1030 - qJ(3) * t1026) * t1058;
t926 = pkin(2) * t1012 - t1013 * qJ(3) + t1014 * t1075 - t995 * t1048 - t955;
t1070 = -pkin(2) - pkin(9);
t1069 = mrSges(3,1) - mrSges(4,2);
t1064 = g(3) * t1023;
t1052 = t1023 * t1030;
t1000 = -qJD(2) * t1049 + t1030 * t1051;
t1025 = sin(qJ(4));
t1029 = cos(qJ(4));
t1005 = qJD(4) + t1049;
t1020 = sin(pkin(11));
t1022 = cos(pkin(11));
t1024 = sin(qJ(6));
t1028 = cos(qJ(6));
t1003 = t1005 ^ 2;
t1056 = t1021 ^ 2 * t1032;
t1047 = t1030 ^ 2 * t1056;
t998 = pkin(3) * t1049 - pkin(9) * t1014;
t999 = (qJD(1) * qJD(2) * t1030 + qJDD(1) * t1026) * t1021;
t920 = -pkin(3) * t1047 - t1064 - qJ(3) * t999 + t1070 * t1000 + (-t993 + (-qJ(3) * t1014 * t1030 - t1026 * t998) * qJD(1)) * t1021 + t1071;
t1059 = g(3) * t1054 + t1026 * t994;
t1041 = -qJ(3) * t1012 + t995 * t1049 + qJDD(3) + t1059;
t922 = pkin(3) * t999 + t1070 * t1013 + (-pkin(3) * t1014 * t1058 - pkin(9) * t1026 * t1056 - t1023 * t993) * t1030 + t1041;
t907 = t1025 * t922 + t1029 * t920;
t980 = t1014 * t1025 + t1029 * t1048;
t981 = t1014 * t1029 - t1025 * t1048;
t956 = pkin(4) * t980 - qJ(5) * t981;
t988 = qJDD(4) + t999;
t901 = -pkin(4) * t1003 + qJ(5) * t988 - t956 * t980 + t907;
t919 = pkin(3) * t1000 - pkin(9) * t1047 + t1014 * t998 - t926;
t952 = qJD(4) * t981 + t1029 * t1000 + t1013 * t1025;
t953 = -qJD(4) * t980 - t1000 * t1025 + t1013 * t1029;
t904 = (t1005 * t980 - t953) * qJ(5) + (t1005 * t981 + t952) * pkin(4) + t919;
t962 = t1005 * t1020 + t1022 * t981;
t896 = -0.2e1 * qJD(5) * t962 - t1020 * t901 + t1022 * t904;
t937 = t1020 * t988 + t1022 * t953;
t961 = t1005 * t1022 - t1020 * t981;
t894 = (t961 * t980 - t937) * pkin(10) + (t961 * t962 + t952) * pkin(5) + t896;
t897 = 0.2e1 * qJD(5) * t961 + t1020 * t904 + t1022 * t901;
t936 = -t1020 * t953 + t1022 * t988;
t943 = pkin(5) * t980 - pkin(10) * t962;
t960 = t961 ^ 2;
t895 = -pkin(5) * t960 + pkin(10) * t936 - t943 * t980 + t897;
t892 = -t1024 * t895 + t1028 * t894;
t933 = -t1024 * t962 + t1028 * t961;
t910 = qJD(6) * t933 + t1024 * t936 + t1028 * t937;
t934 = t1024 * t961 + t1028 * t962;
t915 = -mrSges(7,1) * t933 + mrSges(7,2) * t934;
t979 = qJD(6) + t980;
t923 = -mrSges(7,2) * t979 + mrSges(7,3) * t933;
t950 = qJDD(6) + t952;
t888 = m(7) * t892 + mrSges(7,1) * t950 - mrSges(7,3) * t910 - t915 * t934 + t923 * t979;
t893 = t1024 * t894 + t1028 * t895;
t909 = -qJD(6) * t934 - t1024 * t937 + t1028 * t936;
t924 = mrSges(7,1) * t979 - mrSges(7,3) * t934;
t889 = m(7) * t893 - mrSges(7,2) * t950 + mrSges(7,3) * t909 + t915 * t933 - t924 * t979;
t881 = t1024 * t889 + t1028 * t888;
t938 = -mrSges(6,1) * t961 + mrSges(6,2) * t962;
t941 = -mrSges(6,2) * t980 + mrSges(6,3) * t961;
t879 = m(6) * t896 + mrSges(6,1) * t952 - mrSges(6,3) * t937 - t938 * t962 + t941 * t980 + t881;
t1045 = -t1024 * t888 + t1028 * t889;
t942 = mrSges(6,1) * t980 - mrSges(6,3) * t962;
t880 = m(6) * t897 - mrSges(6,2) * t952 + mrSges(6,3) * t936 + t938 * t961 - t942 * t980 + t1045;
t1046 = -t1020 * t879 + t1022 * t880;
t957 = mrSges(5,1) * t980 + mrSges(5,2) * t981;
t964 = mrSges(5,1) * t1005 - mrSges(5,3) * t981;
t873 = m(5) * t907 - mrSges(5,2) * t988 - mrSges(5,3) * t952 - t1005 * t964 - t957 * t980 + t1046;
t906 = -t1025 * t920 + t1029 * t922;
t900 = -pkin(4) * t988 - qJ(5) * t1003 + t981 * t956 + qJDD(5) - t906;
t898 = -pkin(5) * t936 - pkin(10) * t960 + t943 * t962 + t900;
t1039 = m(7) * t898 - t909 * mrSges(7,1) + mrSges(7,2) * t910 - t933 * t923 + t924 * t934;
t891 = m(6) * t900 - t936 * mrSges(6,1) + mrSges(6,2) * t937 - t961 * t941 + t942 * t962 + t1039;
t963 = -mrSges(5,2) * t1005 - mrSges(5,3) * t980;
t884 = m(5) * t906 + mrSges(5,1) * t988 - mrSges(5,3) * t953 + t1005 * t963 - t957 * t981 - t891;
t1044 = -t1025 * t884 + t1029 * t873;
t971 = -t1021 * t993 - t1064;
t927 = -pkin(2) * t1000 + (-t1014 * t1048 - t999) * qJ(3) + t971 + t1071;
t991 = -mrSges(4,1) * t1048 - mrSges(4,3) * t1014;
t1042 = m(4) * t927 - t999 * mrSges(4,3) + t991 * t1048 + t1044;
t989 = mrSges(3,1) * t1014 - mrSges(3,3) * t1049;
t990 = -mrSges(3,2) * t1014 + mrSges(3,3) * t1048;
t992 = mrSges(4,1) * t1049 + mrSges(4,2) * t1014;
t859 = m(3) * t971 + mrSges(3,2) * t999 - t1069 * t1000 + (-t1030 * t990 + (t989 - t992) * t1026) * t1058 + t1042;
t863 = t1025 * t873 + t1029 * t884;
t1050 = t993 * t1052;
t932 = -pkin(2) * t1013 + t1041 - t1050;
t1040 = -m(4) * t932 - t999 * mrSges(4,1) - t863;
t954 = t1050 - t1059;
t996 = (mrSges(4,2) * t1030 - mrSges(4,3) * t1026) * t1058;
t997 = (-mrSges(3,1) * t1030 + mrSges(3,2) * t1026) * t1058;
t860 = m(3) * t954 - mrSges(3,3) * t999 + (t990 - t991) * t1014 + t1069 * t1013 + (-t996 - t997) * t1049 + t1040;
t875 = t1020 * t880 + t1022 * t879;
t1037 = m(5) * t919 + mrSges(5,1) * t952 + t953 * mrSges(5,2) + t963 * t980 + t981 * t964 + t875;
t1034 = -m(4) * t926 + t1013 * mrSges(4,3) + t1014 * t992 + t996 * t1048 + t1037;
t871 = t1034 + m(3) * t955 - mrSges(3,2) * t1013 - t1014 * t989 + (mrSges(3,3) + mrSges(4,1)) * t1000 + t997 * t1048;
t848 = -t1021 * t859 + t860 * t1052 + t871 * t1053;
t845 = m(2) * t1009 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1032 + t848;
t854 = -t1026 * t860 + t1030 * t871;
t852 = m(2) * t1010 - mrSges(2,1) * t1032 - qJDD(1) * mrSges(2,2) + t854;
t1063 = t1027 * t852 + t1031 * t845;
t1062 = (t1067 * t1026 + t1066 * t1030) * t1058 + t1072 * t1014;
t1061 = (t1068 * t1026 + t1073 * t1030) * t1058 + t1066 * t1014;
t1060 = (t1074 * t1026 + t1068 * t1030) * t1058 + t1067 * t1014;
t847 = t1023 * t859 + t860 * t1054 + t871 * t1055;
t1043 = -t1027 * t845 + t1031 * t852;
t911 = Ifges(7,5) * t934 + Ifges(7,6) * t933 + Ifges(7,3) * t979;
t913 = Ifges(7,1) * t934 + Ifges(7,4) * t933 + Ifges(7,5) * t979;
t882 = -mrSges(7,1) * t898 + mrSges(7,3) * t893 + Ifges(7,4) * t910 + Ifges(7,2) * t909 + Ifges(7,6) * t950 - t911 * t934 + t913 * t979;
t912 = Ifges(7,4) * t934 + Ifges(7,2) * t933 + Ifges(7,6) * t979;
t883 = mrSges(7,2) * t898 - mrSges(7,3) * t892 + Ifges(7,1) * t910 + Ifges(7,4) * t909 + Ifges(7,5) * t950 + t911 * t933 - t912 * t979;
t928 = Ifges(6,5) * t962 + Ifges(6,6) * t961 + Ifges(6,3) * t980;
t930 = Ifges(6,1) * t962 + Ifges(6,4) * t961 + Ifges(6,5) * t980;
t865 = -mrSges(6,1) * t900 + mrSges(6,3) * t897 + Ifges(6,4) * t937 + Ifges(6,2) * t936 + Ifges(6,6) * t952 - pkin(5) * t1039 + pkin(10) * t1045 + t1024 * t883 + t1028 * t882 - t962 * t928 + t980 * t930;
t929 = Ifges(6,4) * t962 + Ifges(6,2) * t961 + Ifges(6,6) * t980;
t867 = mrSges(6,2) * t900 - mrSges(6,3) * t896 + Ifges(6,1) * t937 + Ifges(6,4) * t936 + Ifges(6,5) * t952 - pkin(10) * t881 - t1024 * t882 + t1028 * t883 + t928 * t961 - t929 * t980;
t944 = Ifges(5,5) * t981 - Ifges(5,6) * t980 + Ifges(5,3) * t1005;
t945 = Ifges(5,4) * t981 - Ifges(5,2) * t980 + Ifges(5,6) * t1005;
t849 = mrSges(5,2) * t919 - mrSges(5,3) * t906 + Ifges(5,1) * t953 - Ifges(5,4) * t952 + Ifges(5,5) * t988 - qJ(5) * t875 - t1005 * t945 - t1020 * t865 + t1022 * t867 - t944 * t980;
t1033 = mrSges(7,1) * t892 - mrSges(7,2) * t893 + Ifges(7,5) * t910 + Ifges(7,6) * t909 + Ifges(7,3) * t950 + t934 * t912 - t933 * t913;
t946 = Ifges(5,1) * t981 - Ifges(5,4) * t980 + Ifges(5,5) * t1005;
t855 = -t1033 + t1005 * t946 + Ifges(5,6) * t988 - t981 * t944 + t961 * t930 - t962 * t929 + Ifges(5,4) * t953 - Ifges(6,6) * t936 - Ifges(6,5) * t937 - mrSges(5,1) * t919 + mrSges(5,3) * t907 - mrSges(6,1) * t896 + mrSges(6,2) * t897 - pkin(5) * t881 - pkin(4) * t875 + (-Ifges(5,2) - Ifges(6,3)) * t952;
t862 = mrSges(4,2) * t1013 + t1014 * t991 + t996 * t1049 - t1040;
t839 = mrSges(3,1) * t954 - mrSges(3,2) * t955 + mrSges(4,2) * t932 - mrSges(4,3) * t926 + t1029 * t849 - t1025 * t855 - pkin(9) * t863 - pkin(2) * t862 + qJ(3) * t1034 + t1067 * t999 + t1072 * t1013 + (qJ(3) * mrSges(4,1) + t1066) * t1000 + (t1061 * t1026 - t1060 * t1030) * t1058;
t861 = mrSges(4,2) * t1000 - t992 * t1049 + t1042;
t841 = -mrSges(3,1) * t971 - mrSges(4,1) * t926 + mrSges(4,2) * t927 + mrSges(3,3) * t955 - pkin(2) * t861 + pkin(3) * t1037 - pkin(9) * t1044 + t1073 * t1000 + t1066 * t1013 + t1060 * t1014 - t1025 * t849 - t1029 * t855 - t1062 * t1049 + t1068 * t999;
t1035 = mrSges(5,1) * t906 - mrSges(5,2) * t907 + Ifges(5,5) * t953 - Ifges(5,6) * t952 + Ifges(5,3) * t988 - pkin(4) * t891 + qJ(5) * t1046 + t1020 * t867 + t1022 * t865 + t981 * t945 + t980 * t946;
t843 = mrSges(4,1) * t932 + mrSges(3,2) * t971 - mrSges(3,3) * t954 - mrSges(4,3) * t927 + pkin(3) * t863 - qJ(3) * t861 + t1068 * t1000 + t1067 * t1013 - t1061 * t1014 + t1062 * t1048 + t1074 * t999 + t1035;
t1038 = mrSges(2,1) * t1009 - mrSges(2,2) * t1010 + Ifges(2,3) * qJDD(1) + pkin(1) * t848 + t1023 * t839 + t841 * t1054 + t843 * t1055 + t854 * t1065;
t837 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1009 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t1032 - t1026 * t841 + t1030 * t843 + (-t1021 * t847 - t1023 * t848) * pkin(8);
t836 = mrSges(2,1) * g(3) + mrSges(2,3) * t1010 + Ifges(2,5) * t1032 + Ifges(2,6) * qJDD(1) - pkin(1) * t847 - t1021 * t839 + (pkin(8) * t854 + t1026 * t843 + t1030 * t841) * t1023;
t1 = [-m(1) * g(1) + t1043; -m(1) * g(2) + t1063; (-m(1) - m(2)) * g(3) + t847; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1063 - t1027 * t836 + t1031 * t837; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1043 + t1027 * t837 + t1031 * t836; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1038; t1038; t839; t862; t1035; t891; t1033;];
tauJB  = t1;
