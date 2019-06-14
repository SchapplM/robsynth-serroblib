% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR13_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:40:45
% EndTime: 2019-05-07 15:43:35
% DurationCPUTime: 166.22s
% Computational Cost: add. (2794982->419), mult. (6994726->555), div. (0->0), fcn. (5928108->16), ass. (0->181)
t1087 = cos(qJ(3));
t1043 = sin(pkin(6));
t1086 = pkin(9) * t1043;
t1042 = sin(pkin(7));
t1085 = pkin(10) * t1042;
t1045 = cos(pkin(7));
t1084 = pkin(10) * t1045;
t1046 = cos(pkin(6));
t1083 = t1046 * g(3);
t1051 = sin(qJ(1));
t1055 = cos(qJ(1));
t1035 = t1051 * g(1) - g(2) * t1055;
t1056 = qJD(1) ^ 2;
t1054 = cos(qJ(2));
t1073 = t1046 * t1054;
t1050 = sin(qJ(2));
t1074 = t1046 * t1050;
t1026 = qJDD(1) * pkin(1) + t1056 * t1086 + t1035;
t1036 = -g(1) * t1055 - g(2) * t1051;
t1027 = -pkin(1) * t1056 + qJDD(1) * t1086 + t1036;
t1064 = t1026 * t1073 - t1050 * t1027;
t1076 = t1043 * t1054;
t1004 = -g(3) * t1076 + t1064;
t1038 = qJD(1) * t1046 + qJD(2);
t1068 = qJD(1) * t1076;
t1025 = -mrSges(3,2) * t1038 + mrSges(3,3) * t1068;
t1081 = qJD(1) * t1043;
t1028 = (-mrSges(3,1) * t1054 + mrSges(3,2) * t1050) * t1081;
t1079 = qJD(1) * t1054;
t1029 = (qJD(2) * t1079 + qJDD(1) * t1050) * t1043;
t1037 = qJDD(1) * t1046 + qJDD(2);
t1077 = t1043 * t1050;
t1069 = qJD(1) * t1077;
t1070 = t1045 * t1087;
t1049 = sin(qJ(3));
t1075 = t1045 * t1049;
t1078 = t1042 * t1049;
t1009 = t1038 * t1078 + (t1087 * t1050 + t1054 * t1075) * t1081;
t1020 = t1038 * t1045 - t1042 * t1068 + qJD(3);
t1003 = mrSges(4,1) * t1020 - mrSges(4,3) * t1009;
t1063 = t1045 * t1068;
t1071 = t1042 * t1087;
t1008 = -t1038 * t1071 + t1049 * t1069 - t1087 * t1063;
t1080 = qJD(1) * t1050;
t1030 = (-qJD(2) * t1080 + qJDD(1) * t1054) * t1043;
t1010 = -t1030 * t1042 + t1037 * t1045 + qJDD(3);
t1041 = sin(pkin(13));
t1044 = cos(pkin(13));
t1001 = t1009 * t1044 + t1020 * t1041;
t1048 = sin(qJ(5));
t1053 = cos(qJ(5));
t1007 = qJD(5) + t1008;
t1047 = sin(qJ(6));
t1052 = cos(qJ(6));
t1006 = t1007 ^ 2;
t1000 = -t1009 * t1041 + t1020 * t1044;
t1017 = t1020 ^ 2;
t1018 = (t1038 * t1042 + t1063) * pkin(10);
t1023 = (-pkin(2) * t1054 - t1050 * t1085) * t1081;
t976 = -t1029 * t1084 + t1037 * pkin(2) + t1038 * t1018 + (-g(3) * t1054 - t1023 * t1080) * t1043 + t1064;
t1022 = pkin(2) * t1038 - t1069 * t1084;
t1062 = t1030 * t1045 + t1037 * t1042;
t1072 = t1026 * t1074 + t1054 * t1027;
t977 = -t1038 * t1022 + (-g(3) * t1050 + t1023 * t1079) * t1043 + t1062 * pkin(10) + t1072;
t984 = -t1029 * t1085 - t1030 * pkin(2) - t1083 + (-t1026 + (-t1018 * t1054 + t1022 * t1050) * qJD(1)) * t1043;
t950 = t976 * t1075 + t984 * t1078 + t1087 * t977;
t995 = pkin(3) * t1008 - qJ(4) * t1009;
t937 = -pkin(3) * t1017 + qJ(4) * t1010 - t1008 * t995 + t950;
t960 = -t1042 * t976 + t1045 * t984;
t993 = qJD(3) * t1009 + t1029 * t1049 - t1030 * t1070 - t1037 * t1071;
t994 = -t1008 * qJD(3) + t1087 * t1029 + t1062 * t1049;
t940 = (t1008 * t1020 - t994) * qJ(4) + (t1009 * t1020 + t993) * pkin(3) + t960;
t929 = -0.2e1 * qJD(4) * t1001 - t1041 * t937 + t1044 * t940;
t981 = t1010 * t1041 + t1044 * t994;
t926 = (t1000 * t1008 - t981) * pkin(11) + (t1000 * t1001 + t993) * pkin(4) + t929;
t930 = 0.2e1 * qJD(4) * t1000 + t1041 * t940 + t1044 * t937;
t980 = t1010 * t1044 - t1041 * t994;
t988 = pkin(4) * t1008 - pkin(11) * t1001;
t999 = t1000 ^ 2;
t928 = -pkin(4) * t999 + pkin(11) * t980 - t1008 * t988 + t930;
t923 = t1048 * t926 + t1053 * t928;
t974 = t1000 * t1053 - t1001 * t1048;
t975 = t1000 * t1048 + t1001 * t1053;
t959 = -pkin(5) * t974 - pkin(12) * t975;
t992 = qJDD(5) + t993;
t921 = -pkin(5) * t1006 + pkin(12) * t992 + t959 * t974 + t923;
t949 = -t1049 * t977 + t976 * t1070 + t984 * t1071;
t936 = -t1010 * pkin(3) - t1017 * qJ(4) + t1009 * t995 + qJDD(4) - t949;
t931 = -t980 * pkin(4) - t999 * pkin(11) + t1001 * t988 + t936;
t947 = -qJD(5) * t975 - t1048 * t981 + t1053 * t980;
t948 = qJD(5) * t974 + t1048 * t980 + t1053 * t981;
t924 = (-t1007 * t974 - t948) * pkin(12) + (t1007 * t975 - t947) * pkin(5) + t931;
t918 = -t1047 * t921 + t1052 * t924;
t961 = t1007 * t1052 - t1047 * t975;
t934 = qJD(6) * t961 + t1047 * t992 + t1052 * t948;
t946 = qJDD(6) - t947;
t962 = t1007 * t1047 + t1052 * t975;
t951 = -mrSges(7,1) * t961 + mrSges(7,2) * t962;
t973 = qJD(6) - t974;
t952 = -mrSges(7,2) * t973 + mrSges(7,3) * t961;
t915 = m(7) * t918 + mrSges(7,1) * t946 - mrSges(7,3) * t934 - t951 * t962 + t952 * t973;
t919 = t1047 * t924 + t1052 * t921;
t933 = -qJD(6) * t962 - t1047 * t948 + t1052 * t992;
t953 = mrSges(7,1) * t973 - mrSges(7,3) * t962;
t916 = m(7) * t919 - mrSges(7,2) * t946 + mrSges(7,3) * t933 + t951 * t961 - t953 * t973;
t907 = -t1047 * t915 + t1052 * t916;
t958 = -mrSges(6,1) * t974 + mrSges(6,2) * t975;
t964 = mrSges(6,1) * t1007 - mrSges(6,3) * t975;
t901 = m(6) * t923 - mrSges(6,2) * t992 + mrSges(6,3) * t947 - t1007 * t964 + t958 * t974 + t907;
t922 = -t1048 * t928 + t1053 * t926;
t920 = -pkin(5) * t992 - pkin(12) * t1006 + t959 * t975 - t922;
t917 = -m(7) * t920 + t933 * mrSges(7,1) - mrSges(7,2) * t934 + t961 * t952 - t953 * t962;
t963 = -mrSges(6,2) * t1007 + mrSges(6,3) * t974;
t911 = m(6) * t922 + mrSges(6,1) * t992 - mrSges(6,3) * t948 + t1007 * t963 - t958 * t975 + t917;
t898 = t1048 * t901 + t1053 * t911;
t979 = -mrSges(5,1) * t1000 + mrSges(5,2) * t1001;
t986 = -mrSges(5,2) * t1008 + mrSges(5,3) * t1000;
t896 = m(5) * t929 + mrSges(5,1) * t993 - mrSges(5,3) * t981 - t1001 * t979 + t1008 * t986 + t898;
t1066 = -t1048 * t911 + t1053 * t901;
t987 = mrSges(5,1) * t1008 - mrSges(5,3) * t1001;
t897 = m(5) * t930 - mrSges(5,2) * t993 + mrSges(5,3) * t980 + t1000 * t979 - t1008 * t987 + t1066;
t1067 = -t1041 * t896 + t1044 * t897;
t996 = mrSges(4,1) * t1008 + mrSges(4,2) * t1009;
t887 = m(4) * t950 - mrSges(4,2) * t1010 - mrSges(4,3) * t993 - t1003 * t1020 - t1008 * t996 + t1067;
t1002 = -mrSges(4,2) * t1020 - mrSges(4,3) * t1008;
t890 = t1041 * t897 + t1044 * t896;
t889 = m(4) * t960 + mrSges(4,1) * t993 + mrSges(4,2) * t994 + t1002 * t1008 + t1003 * t1009 + t890;
t906 = t1047 * t916 + t1052 * t915;
t1059 = m(6) * t931 - t947 * mrSges(6,1) + mrSges(6,2) * t948 - t974 * t963 + t964 * t975 + t906;
t905 = m(5) * t936 - t980 * mrSges(5,1) + mrSges(5,2) * t981 - t1000 * t986 + t1001 * t987 + t1059;
t904 = m(4) * t949 + mrSges(4,1) * t1010 - mrSges(4,3) * t994 + t1002 * t1020 - t1009 * t996 - t905;
t876 = -t1042 * t889 + t904 * t1070 + t887 * t1075;
t872 = m(3) * t1004 + mrSges(3,1) * t1037 - mrSges(3,3) * t1029 + t1025 * t1038 - t1028 * t1069 + t876;
t1014 = -t1043 * t1026 - t1083;
t1024 = mrSges(3,1) * t1038 - mrSges(3,3) * t1069;
t875 = t1045 * t889 + t904 * t1071 + t887 * t1078;
t874 = m(3) * t1014 - t1030 * mrSges(3,1) + t1029 * mrSges(3,2) + (t1024 * t1050 - t1025 * t1054) * t1081 + t875;
t1005 = -g(3) * t1077 + t1072;
t883 = -t1049 * t904 + t1087 * t887;
t882 = m(3) * t1005 - mrSges(3,2) * t1037 + mrSges(3,3) * t1030 - t1024 * t1038 + t1028 * t1068 + t883;
t861 = -t1043 * t874 + t872 * t1073 + t882 * t1074;
t858 = m(2) * t1035 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1056 + t861;
t868 = -t1050 * t872 + t1054 * t882;
t866 = m(2) * t1036 - mrSges(2,1) * t1056 - qJDD(1) * mrSges(2,2) + t868;
t1082 = t1051 * t866 + t1055 * t858;
t860 = t1046 * t874 + t872 * t1076 + t882 * t1077;
t1065 = -t1051 * t858 + t1055 * t866;
t941 = Ifges(7,5) * t962 + Ifges(7,6) * t961 + Ifges(7,3) * t973;
t943 = Ifges(7,1) * t962 + Ifges(7,4) * t961 + Ifges(7,5) * t973;
t908 = -mrSges(7,1) * t920 + mrSges(7,3) * t919 + Ifges(7,4) * t934 + Ifges(7,2) * t933 + Ifges(7,6) * t946 - t941 * t962 + t943 * t973;
t942 = Ifges(7,4) * t962 + Ifges(7,2) * t961 + Ifges(7,6) * t973;
t909 = mrSges(7,2) * t920 - mrSges(7,3) * t918 + Ifges(7,1) * t934 + Ifges(7,4) * t933 + Ifges(7,5) * t946 + t941 * t961 - t942 * t973;
t954 = Ifges(6,5) * t975 + Ifges(6,6) * t974 + Ifges(6,3) * t1007;
t955 = Ifges(6,4) * t975 + Ifges(6,2) * t974 + Ifges(6,6) * t1007;
t891 = mrSges(6,2) * t931 - mrSges(6,3) * t922 + Ifges(6,1) * t948 + Ifges(6,4) * t947 + Ifges(6,5) * t992 - pkin(12) * t906 - t1007 * t955 - t1047 * t908 + t1052 * t909 + t954 * t974;
t1058 = mrSges(7,1) * t918 - mrSges(7,2) * t919 + Ifges(7,5) * t934 + Ifges(7,6) * t933 + Ifges(7,3) * t946 + t942 * t962 - t943 * t961;
t956 = Ifges(6,1) * t975 + Ifges(6,4) * t974 + Ifges(6,5) * t1007;
t892 = -mrSges(6,1) * t931 + mrSges(6,3) * t923 + Ifges(6,4) * t948 + Ifges(6,2) * t947 + Ifges(6,6) * t992 - pkin(5) * t906 + t1007 * t956 - t954 * t975 - t1058;
t965 = Ifges(5,5) * t1001 + Ifges(5,6) * t1000 + Ifges(5,3) * t1008;
t967 = Ifges(5,1) * t1001 + Ifges(5,4) * t1000 + Ifges(5,5) * t1008;
t877 = -mrSges(5,1) * t936 + mrSges(5,3) * t930 + Ifges(5,4) * t981 + Ifges(5,2) * t980 + Ifges(5,6) * t993 - pkin(4) * t1059 + pkin(11) * t1066 - t1001 * t965 + t1008 * t967 + t1048 * t891 + t1053 * t892;
t966 = Ifges(5,4) * t1001 + Ifges(5,2) * t1000 + Ifges(5,6) * t1008;
t878 = mrSges(5,2) * t936 - mrSges(5,3) * t929 + Ifges(5,1) * t981 + Ifges(5,4) * t980 + Ifges(5,5) * t993 - pkin(11) * t898 + t1000 * t965 - t1008 * t966 - t1048 * t892 + t1053 * t891;
t989 = Ifges(4,5) * t1009 - Ifges(4,6) * t1008 + Ifges(4,3) * t1020;
t990 = Ifges(4,4) * t1009 - Ifges(4,2) * t1008 + Ifges(4,6) * t1020;
t863 = mrSges(4,2) * t960 - mrSges(4,3) * t949 + Ifges(4,1) * t994 - Ifges(4,4) * t993 + Ifges(4,5) * t1010 - qJ(4) * t890 - t1008 * t989 - t1020 * t990 - t1041 * t877 + t1044 * t878;
t1057 = mrSges(6,1) * t922 - mrSges(6,2) * t923 + Ifges(6,5) * t948 + Ifges(6,6) * t947 + Ifges(6,3) * t992 + pkin(5) * t917 + pkin(12) * t907 + t1047 * t909 + t1052 * t908 + t975 * t955 - t974 * t956;
t991 = Ifges(4,1) * t1009 - Ifges(4,4) * t1008 + Ifges(4,5) * t1020;
t869 = -t1057 + Ifges(4,6) * t1010 + t1020 * t991 - t1009 * t989 + t1000 * t967 - t1001 * t966 + Ifges(4,4) * t994 - Ifges(5,6) * t980 - Ifges(5,5) * t981 - mrSges(4,1) * t960 + mrSges(4,3) * t950 - mrSges(5,1) * t929 + mrSges(5,2) * t930 + (-Ifges(5,3) - Ifges(4,2)) * t993 - pkin(3) * t890 - pkin(4) * t898;
t1061 = pkin(10) * t883 + t1049 * t863 + t1087 * t869;
t1012 = Ifges(3,6) * t1038 + (Ifges(3,4) * t1050 + Ifges(3,2) * t1054) * t1081;
t1013 = Ifges(3,5) * t1038 + (Ifges(3,1) * t1050 + Ifges(3,4) * t1054) * t1081;
t862 = mrSges(4,1) * t949 - mrSges(4,2) * t950 + Ifges(4,5) * t994 - Ifges(4,6) * t993 + Ifges(4,3) * t1010 - pkin(3) * t905 + qJ(4) * t1067 + t1008 * t991 + t1009 * t990 + t1041 * t878 + t1044 * t877;
t852 = mrSges(3,1) * t1004 - mrSges(3,2) * t1005 + Ifges(3,5) * t1029 + Ifges(3,6) * t1030 + Ifges(3,3) * t1037 + pkin(2) * t876 + t1045 * t862 + (t1012 * t1050 - t1013 * t1054) * t1081 + t1061 * t1042;
t1011 = Ifges(3,3) * t1038 + (Ifges(3,5) * t1050 + Ifges(3,6) * t1054) * t1081;
t854 = -mrSges(3,1) * t1014 + mrSges(3,3) * t1005 + Ifges(3,4) * t1029 + Ifges(3,2) * t1030 + Ifges(3,6) * t1037 - pkin(2) * t875 - t1011 * t1069 + t1038 * t1013 - t1042 * t862 + t1061 * t1045;
t856 = t1011 * t1068 + mrSges(3,2) * t1014 - mrSges(3,3) * t1004 + t1087 * t863 + Ifges(3,1) * t1029 + Ifges(3,4) * t1030 + Ifges(3,5) * t1037 - t1038 * t1012 - t1049 * t869 + (-t1042 * t875 - t1045 * t876) * pkin(10);
t1060 = mrSges(2,1) * t1035 - mrSges(2,2) * t1036 + Ifges(2,3) * qJDD(1) + pkin(1) * t861 + t1046 * t852 + t854 * t1076 + t856 * t1077 + t868 * t1086;
t850 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1035 + Ifges(2,5) * qJDD(1) - t1056 * Ifges(2,6) - t1050 * t854 + t1054 * t856 + (-t1043 * t860 - t1046 * t861) * pkin(9);
t849 = mrSges(2,1) * g(3) + mrSges(2,3) * t1036 + t1056 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t860 - t1043 * t852 + (pkin(9) * t868 + t1050 * t856 + t1054 * t854) * t1046;
t1 = [-m(1) * g(1) + t1065; -m(1) * g(2) + t1082; (-m(1) - m(2)) * g(3) + t860; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1082 - t1051 * t849 + t1055 * t850; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1065 + t1051 * t850 + t1055 * t849; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1060; t1060; t852; t862; t905; t1057; t1058;];
tauJB  = t1;
