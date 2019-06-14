% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR12
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
% Datum: 2019-05-07 00:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 00:38:48
% EndTime: 2019-05-07 00:39:11
% DurationCPUTime: 22.26s
% Computational Cost: add. (366664->385), mult. (818836->480), div. (0->0), fcn. (613392->12), ass. (0->165)
t1084 = -2 * qJD(3);
t1083 = Ifges(3,1) + Ifges(4,2);
t1077 = Ifges(3,4) + Ifges(4,6);
t1076 = Ifges(3,5) - Ifges(4,4);
t1082 = Ifges(3,2) + Ifges(4,3);
t1075 = Ifges(3,6) - Ifges(4,5);
t1081 = Ifges(3,3) + Ifges(4,1);
t1028 = sin(pkin(6));
t1033 = sin(qJ(2));
t1064 = t1028 * t1033;
t1020 = qJD(1) * t1064;
t1029 = cos(pkin(6));
t1022 = qJD(1) * t1029 + qJD(2);
t1080 = (pkin(2) * t1022 + t1084) * t1020;
t1038 = cos(qJ(2));
t1067 = qJD(1) * t1028;
t1001 = (-pkin(2) * t1038 - qJ(3) * t1033) * t1067;
t1019 = t1022 ^ 2;
t1021 = qJDD(1) * t1029 + qJDD(2);
t1063 = t1028 * t1038;
t1058 = qJD(1) * t1063;
t1034 = sin(qJ(1));
t1039 = cos(qJ(1));
t1017 = -g(1) * t1039 - g(2) * t1034;
t1040 = qJD(1) ^ 2;
t1060 = qJDD(1) * t1028;
t1000 = -pkin(1) * t1040 + pkin(8) * t1060 + t1017;
t1062 = t1029 * t1033;
t1016 = t1034 * g(1) - g(2) * t1039;
t1074 = pkin(8) * t1028;
t999 = qJDD(1) * pkin(1) + t1040 * t1074 + t1016;
t962 = -g(3) * t1064 + t1038 * t1000 + t999 * t1062;
t941 = t1019 * pkin(2) - t1021 * qJ(3) - t1001 * t1058 + t1022 * t1084 - t962;
t1079 = -pkin(2) - pkin(9);
t1078 = mrSges(3,1) - mrSges(4,2);
t1073 = t1029 * g(3);
t1061 = t1029 * t1038;
t1005 = (qJD(1) * qJD(2) * t1038 + qJDD(1) * t1033) * t1028;
t1006 = -qJD(2) * t1020 + t1038 * t1060;
t1032 = sin(qJ(4));
t1037 = cos(qJ(4));
t1012 = t1020 + qJD(4);
t1031 = sin(qJ(5));
t1036 = cos(qJ(5));
t1009 = qJD(5) + t1012;
t1030 = sin(qJ(6));
t1035 = cos(qJ(6));
t1008 = t1009 ^ 2;
t1004 = pkin(3) * t1020 - pkin(9) * t1022;
t1065 = t1028 ^ 2 * t1040;
t1057 = t1038 ^ 2 * t1065;
t929 = -pkin(3) * t1057 - t1073 - t1005 * qJ(3) + t1079 * t1006 + (-t999 + (-qJ(3) * t1022 * t1038 - t1004 * t1033) * qJD(1)) * t1028 + t1080;
t1068 = g(3) * t1063 + t1033 * t1000;
t1051 = -t1019 * qJ(3) + t1001 * t1020 + qJDD(3) + t1068;
t932 = t1005 * pkin(3) + t1079 * t1021 + (-pkin(3) * t1022 * t1067 - pkin(9) * t1033 * t1065 - t1029 * t999) * t1038 + t1051;
t908 = -t1032 * t929 + t1037 * t932;
t986 = -t1022 * t1032 - t1037 * t1058;
t960 = qJD(4) * t986 - t1006 * t1032 + t1021 * t1037;
t987 = t1022 * t1037 - t1032 * t1058;
t994 = qJDD(4) + t1005;
t904 = (t1012 * t986 - t960) * pkin(10) + (t986 * t987 + t994) * pkin(4) + t908;
t909 = t1032 * t932 + t1037 * t929;
t959 = -qJD(4) * t987 - t1006 * t1037 - t1021 * t1032;
t969 = pkin(4) * t1012 - pkin(10) * t987;
t985 = t986 ^ 2;
t906 = -pkin(4) * t985 + pkin(10) * t959 - t1012 * t969 + t909;
t901 = t1031 * t904 + t1036 * t906;
t964 = -t1031 * t987 + t1036 * t986;
t965 = t1031 * t986 + t1036 * t987;
t944 = -pkin(5) * t964 - pkin(11) * t965;
t991 = qJDD(5) + t994;
t898 = -pkin(5) * t1008 + pkin(11) * t991 + t944 * t964 + t901;
t928 = t1006 * pkin(3) - pkin(9) * t1057 + t1022 * t1004 - t941;
t911 = -t959 * pkin(4) - t985 * pkin(10) + t987 * t969 + t928;
t923 = -qJD(5) * t965 - t1031 * t960 + t1036 * t959;
t924 = qJD(5) * t964 + t1031 * t959 + t1036 * t960;
t902 = t911 + (t1009 * t965 - t923) * pkin(5) + (-t1009 * t964 - t924) * pkin(11);
t895 = -t1030 * t898 + t1035 * t902;
t947 = t1009 * t1035 - t1030 * t965;
t914 = qJD(6) * t947 + t1030 * t991 + t1035 * t924;
t922 = qJDD(6) - t923;
t948 = t1009 * t1030 + t1035 * t965;
t933 = -mrSges(7,1) * t947 + mrSges(7,2) * t948;
t963 = qJD(6) - t964;
t934 = -mrSges(7,2) * t963 + mrSges(7,3) * t947;
t892 = m(7) * t895 + mrSges(7,1) * t922 - mrSges(7,3) * t914 - t933 * t948 + t934 * t963;
t896 = t1030 * t902 + t1035 * t898;
t913 = -qJD(6) * t948 - t1030 * t924 + t1035 * t991;
t935 = mrSges(7,1) * t963 - mrSges(7,3) * t948;
t893 = m(7) * t896 - mrSges(7,2) * t922 + mrSges(7,3) * t913 + t933 * t947 - t935 * t963;
t1056 = -t1030 * t892 + t1035 * t893;
t943 = -mrSges(6,1) * t964 + mrSges(6,2) * t965;
t950 = mrSges(6,1) * t1009 - mrSges(6,3) * t965;
t880 = m(6) * t901 - mrSges(6,2) * t991 + mrSges(6,3) * t923 - t1009 * t950 + t943 * t964 + t1056;
t900 = -t1031 * t906 + t1036 * t904;
t897 = -pkin(5) * t991 - pkin(11) * t1008 + t944 * t965 - t900;
t1049 = -m(7) * t897 + t913 * mrSges(7,1) - mrSges(7,2) * t914 + t947 * t934 - t935 * t948;
t949 = -mrSges(6,2) * t1009 + mrSges(6,3) * t964;
t888 = m(6) * t900 + mrSges(6,1) * t991 - mrSges(6,3) * t924 + t1009 * t949 - t943 * t965 + t1049;
t873 = t1031 * t880 + t1036 * t888;
t966 = -mrSges(5,1) * t986 + mrSges(5,2) * t987;
t967 = -mrSges(5,2) * t1012 + mrSges(5,3) * t986;
t870 = m(5) * t908 + mrSges(5,1) * t994 - mrSges(5,3) * t960 + t1012 * t967 - t966 * t987 + t873;
t1055 = -t1031 * t888 + t1036 * t880;
t968 = mrSges(5,1) * t1012 - mrSges(5,3) * t987;
t871 = m(5) * t909 - mrSges(5,2) * t994 + mrSges(5,3) * t959 - t1012 * t968 + t966 * t986 + t1055;
t1054 = -t1032 * t870 + t1037 * t871;
t976 = -t1028 * t999 - t1073;
t942 = -t1006 * pkin(2) + (-t1022 * t1058 - t1005) * qJ(3) + t976 + t1080;
t997 = -mrSges(4,1) * t1058 - mrSges(4,3) * t1022;
t1052 = m(4) * t942 - t1005 * mrSges(4,3) + t997 * t1058 + t1054;
t995 = mrSges(3,1) * t1022 - mrSges(3,3) * t1020;
t996 = -mrSges(3,2) * t1022 + mrSges(3,3) * t1058;
t998 = mrSges(4,1) * t1020 + mrSges(4,2) * t1022;
t862 = m(3) * t976 + t1005 * mrSges(3,2) - t1078 * t1006 + (-t1038 * t996 + (t995 - t998) * t1033) * t1067 + t1052;
t1002 = (mrSges(4,2) * t1038 - mrSges(4,3) * t1033) * t1067;
t1003 = (-mrSges(3,1) * t1038 + mrSges(3,2) * t1033) * t1067;
t866 = t1032 * t871 + t1037 * t870;
t1059 = t999 * t1061;
t945 = -t1021 * pkin(2) + t1051 - t1059;
t1050 = -m(4) * t945 - t1005 * mrSges(4,1) - t866;
t961 = t1059 - t1068;
t863 = m(3) * t961 - t1005 * mrSges(3,3) + (t996 - t997) * t1022 + t1078 * t1021 + (-t1002 - t1003) * t1020 + t1050;
t882 = t1030 * t893 + t1035 * t892;
t1047 = m(6) * t911 - t923 * mrSges(6,1) + t924 * mrSges(6,2) - t964 * t949 + t965 * t950 + t882;
t1043 = -m(5) * t928 + t959 * mrSges(5,1) - t960 * mrSges(5,2) + t986 * t967 - t987 * t968 - t1047;
t1042 = -m(4) * t941 + t1021 * mrSges(4,3) + t1002 * t1058 + t1022 * t998 - t1043;
t877 = t1042 + t1003 * t1058 + (mrSges(3,3) + mrSges(4,1)) * t1006 - t1021 * mrSges(3,2) - t1022 * t995 + m(3) * t962;
t851 = -t1028 * t862 + t863 * t1061 + t877 * t1062;
t848 = m(2) * t1016 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1040 + t851;
t858 = -t1033 * t863 + t1038 * t877;
t856 = m(2) * t1017 - mrSges(2,1) * t1040 - qJDD(1) * mrSges(2,2) + t858;
t1072 = t1034 * t856 + t1039 * t848;
t1071 = (t1076 * t1033 + t1075 * t1038) * t1067 + t1081 * t1022;
t1070 = (t1077 * t1033 + t1082 * t1038) * t1067 + t1075 * t1022;
t1069 = (t1083 * t1033 + t1077 * t1038) * t1067 + t1076 * t1022;
t850 = t1029 * t862 + t863 * t1063 + t877 * t1064;
t1053 = -t1034 * t848 + t1039 * t856;
t915 = Ifges(7,5) * t948 + Ifges(7,6) * t947 + Ifges(7,3) * t963;
t917 = Ifges(7,1) * t948 + Ifges(7,4) * t947 + Ifges(7,5) * t963;
t885 = -mrSges(7,1) * t897 + mrSges(7,3) * t896 + Ifges(7,4) * t914 + Ifges(7,2) * t913 + Ifges(7,6) * t922 - t915 * t948 + t917 * t963;
t916 = Ifges(7,4) * t948 + Ifges(7,2) * t947 + Ifges(7,6) * t963;
t886 = mrSges(7,2) * t897 - mrSges(7,3) * t895 + Ifges(7,1) * t914 + Ifges(7,4) * t913 + Ifges(7,5) * t922 + t915 * t947 - t916 * t963;
t936 = Ifges(6,5) * t965 + Ifges(6,6) * t964 + Ifges(6,3) * t1009;
t937 = Ifges(6,4) * t965 + Ifges(6,2) * t964 + Ifges(6,6) * t1009;
t867 = mrSges(6,2) * t911 - mrSges(6,3) * t900 + Ifges(6,1) * t924 + Ifges(6,4) * t923 + Ifges(6,5) * t991 - pkin(11) * t882 - t1009 * t937 - t1030 * t885 + t1035 * t886 + t936 * t964;
t1044 = mrSges(7,1) * t895 - mrSges(7,2) * t896 + Ifges(7,5) * t914 + Ifges(7,6) * t913 + Ifges(7,3) * t922 + t916 * t948 - t917 * t947;
t938 = Ifges(6,1) * t965 + Ifges(6,4) * t964 + Ifges(6,5) * t1009;
t868 = -mrSges(6,1) * t911 + mrSges(6,3) * t901 + Ifges(6,4) * t924 + Ifges(6,2) * t923 + Ifges(6,6) * t991 - pkin(5) * t882 + t1009 * t938 - t936 * t965 - t1044;
t951 = Ifges(5,5) * t987 + Ifges(5,6) * t986 + Ifges(5,3) * t1012;
t953 = Ifges(5,1) * t987 + Ifges(5,4) * t986 + Ifges(5,5) * t1012;
t852 = -mrSges(5,1) * t928 + mrSges(5,3) * t909 + Ifges(5,4) * t960 + Ifges(5,2) * t959 + Ifges(5,6) * t994 - pkin(4) * t1047 + pkin(10) * t1055 + t1012 * t953 + t1031 * t867 + t1036 * t868 - t987 * t951;
t952 = Ifges(5,4) * t987 + Ifges(5,2) * t986 + Ifges(5,6) * t1012;
t853 = mrSges(5,2) * t928 - mrSges(5,3) * t908 + Ifges(5,1) * t960 + Ifges(5,4) * t959 + Ifges(5,5) * t994 - pkin(10) * t873 - t1012 * t952 - t1031 * t868 + t1036 * t867 + t951 * t986;
t865 = t1021 * mrSges(4,2) + t1002 * t1020 + t1022 * t997 - t1050;
t842 = mrSges(3,1) * t961 - mrSges(3,2) * t962 + mrSges(4,2) * t945 - mrSges(4,3) * t941 + t1037 * t853 - t1032 * t852 - pkin(9) * t866 - pkin(2) * t865 + qJ(3) * t1042 + t1081 * t1021 + (mrSges(4,1) * qJ(3) + t1075) * t1006 + t1076 * t1005 + (t1070 * t1033 - t1069 * t1038) * t1067;
t864 = t1006 * mrSges(4,2) - t998 * t1020 + t1052;
t844 = -mrSges(3,1) * t976 - mrSges(4,1) * t941 + mrSges(4,2) * t942 + mrSges(3,3) * t962 - pkin(2) * t864 - pkin(3) * t1043 - pkin(9) * t1054 + t1077 * t1005 + t1082 * t1006 - t1071 * t1020 + t1075 * t1021 + t1069 * t1022 - t1032 * t853 - t1037 * t852;
t1045 = mrSges(6,1) * t900 - mrSges(6,2) * t901 + Ifges(6,5) * t924 + Ifges(6,6) * t923 + Ifges(6,3) * t991 + pkin(5) * t1049 + pkin(11) * t1056 + t1030 * t886 + t1035 * t885 + t965 * t937 - t964 * t938;
t1041 = mrSges(5,1) * t908 - mrSges(5,2) * t909 + Ifges(5,5) * t960 + Ifges(5,6) * t959 + Ifges(5,3) * t994 + pkin(4) * t873 + t987 * t952 - t986 * t953 + t1045;
t846 = mrSges(4,1) * t945 + mrSges(3,2) * t976 - mrSges(3,3) * t961 - mrSges(4,3) * t942 + pkin(3) * t866 - qJ(3) * t864 + t1083 * t1005 + t1077 * t1006 + t1076 * t1021 - t1070 * t1022 + t1071 * t1058 + t1041;
t1048 = mrSges(2,1) * t1016 - mrSges(2,2) * t1017 + Ifges(2,3) * qJDD(1) + pkin(1) * t851 + t1029 * t842 + t844 * t1063 + t846 * t1064 + t858 * t1074;
t840 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1016 + Ifges(2,5) * qJDD(1) - t1040 * Ifges(2,6) - t1033 * t844 + t1038 * t846 + (-t1028 * t850 - t1029 * t851) * pkin(8);
t839 = mrSges(2,1) * g(3) + mrSges(2,3) * t1017 + t1040 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t850 - t1028 * t842 + (pkin(8) * t858 + t1033 * t846 + t1038 * t844) * t1029;
t1 = [-m(1) * g(1) + t1053; -m(1) * g(2) + t1072; (-m(1) - m(2)) * g(3) + t850; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1072 - t1034 * t839 + t1039 * t840; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1053 + t1034 * t840 + t1039 * t839; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1048; t1048; t842; t865; t1041; t1045; t1044;];
tauJB  = t1;
