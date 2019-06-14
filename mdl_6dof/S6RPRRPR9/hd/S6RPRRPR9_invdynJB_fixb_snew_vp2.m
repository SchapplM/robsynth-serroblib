% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 23:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:24:38
% EndTime: 2019-05-05 23:26:48
% DurationCPUTime: 134.91s
% Computational Cost: add. (2037381->405), mult. (6394801->553), div. (0->0), fcn. (5502847->16), ass. (0->186)
t1054 = sin(pkin(12));
t1056 = sin(pkin(6));
t1058 = cos(pkin(12));
t1060 = cos(pkin(6));
t1063 = sin(qJ(3));
t1059 = cos(pkin(7));
t1067 = cos(qJ(3));
t1098 = t1059 * t1067;
t1055 = sin(pkin(7));
t1103 = t1055 * t1067;
t1073 = (-t1054 * t1063 + t1058 * t1098) * t1056 + t1060 * t1103;
t1026 = t1073 * qJD(1);
t1099 = t1059 * t1063;
t1104 = t1055 * t1063;
t1075 = t1060 * t1104 + (t1054 * t1067 + t1058 * t1099) * t1056;
t1027 = t1075 * qJD(1);
t1015 = -t1027 * qJD(3) + t1073 * qJDD(1);
t1101 = t1056 * t1059;
t1038 = (t1055 * t1060 + t1058 * t1101) * qJD(1) * pkin(9);
t1064 = sin(qJ(1));
t1068 = cos(qJ(1));
t1049 = t1064 * g(1) - g(2) * t1068;
t1069 = qJD(1) ^ 2;
t1109 = qJ(2) * t1056;
t1041 = qJDD(1) * pkin(1) + t1069 * t1109 + t1049;
t1113 = pkin(9) * t1054;
t1043 = (pkin(2) * t1060 - t1101 * t1113) * qJD(1);
t1085 = -pkin(2) * t1058 - t1055 * t1113;
t1094 = -t1060 * g(3) + qJDD(2);
t1005 = (-t1041 + t1085 * qJDD(1) + (-t1038 * t1058 + t1043 * t1054) * qJD(1)) * t1056 + t1094;
t1050 = -g(1) * t1068 - g(2) * t1064;
t1042 = -pkin(1) * t1069 + qJDD(1) * t1109 + t1050;
t1108 = qJD(1) * t1056;
t1111 = pkin(9) * qJDD(1);
t1082 = qJD(1) * t1085 * t1108 + t1059 * t1111;
t1097 = qJD(2) * t1108;
t1100 = t1058 * t1060;
t1102 = t1056 * t1058;
t1086 = -g(3) * t1102 + t1041 * t1100 - 0.2e1 * t1054 * t1097;
t995 = (pkin(2) * qJDD(1) + qJD(1) * t1038) * t1060 + (-t1082 * t1056 - t1042) * t1054 + t1086;
t1105 = t1054 * t1060;
t1095 = t1041 * t1105 + (t1042 + 0.2e1 * t1097) * t1058;
t996 = (-qJD(1) * t1043 + t1055 * t1111) * t1060 + (-g(3) * t1054 + t1082 * t1058) * t1056 + t1095;
t964 = -t1063 * t996 + (t1005 * t1055 + t1059 * t995) * t1067;
t1080 = -t1055 * t1102 + t1059 * t1060;
t1039 = t1080 * qJD(1) + qJD(3);
t1062 = sin(qJ(4));
t1066 = cos(qJ(4));
t1020 = -t1027 * t1062 + t1039 * t1066;
t1021 = t1027 * t1066 + t1039 * t1062;
t1053 = sin(pkin(13));
t1057 = cos(pkin(13));
t1000 = t1020 * t1053 + t1021 * t1057;
t1012 = qJDD(4) - t1015;
t1061 = sin(qJ(6));
t1065 = cos(qJ(6));
t1025 = qJD(4) - t1026;
t1024 = t1025 ^ 2;
t1115 = 2 * qJD(5);
t1014 = -pkin(3) * t1026 - pkin(10) * t1027;
t1035 = t1039 ^ 2;
t1036 = t1080 * qJDD(1) + qJDD(3);
t965 = t1005 * t1104 + t1067 * t996 + t995 * t1099;
t956 = -pkin(3) * t1035 + pkin(10) * t1036 + t1014 * t1026 + t965;
t1016 = t1026 * qJD(3) + t1075 * qJDD(1);
t979 = t1059 * t1005 - t1055 * t995;
t959 = (-t1026 * t1039 - t1016) * pkin(10) + (t1027 * t1039 - t1015) * pkin(3) + t979;
t948 = -t1062 * t956 + t1066 * t959;
t991 = qJD(4) * t1020 + t1016 * t1066 + t1036 * t1062;
t945 = (t1020 * t1025 - t991) * qJ(5) + (t1020 * t1021 + t1012) * pkin(4) + t948;
t1007 = pkin(4) * t1025 - qJ(5) * t1021;
t1019 = t1020 ^ 2;
t949 = t1062 * t959 + t1066 * t956;
t990 = -qJD(4) * t1021 - t1016 * t1062 + t1036 * t1066;
t947 = -pkin(4) * t1019 + qJ(5) * t990 - t1007 * t1025 + t949;
t999 = t1020 * t1057 - t1021 * t1053;
t942 = t1053 * t945 + t1057 * t947 + t999 * t1115;
t978 = -pkin(5) * t999 - pkin(11) * t1000;
t940 = -pkin(5) * t1024 + pkin(11) * t1012 + t978 * t999 + t942;
t955 = -t1036 * pkin(3) - t1035 * pkin(10) + t1027 * t1014 - t964;
t950 = -t990 * pkin(4) - t1019 * qJ(5) + t1021 * t1007 + qJDD(5) + t955;
t969 = -t1053 * t991 + t1057 * t990;
t970 = t1053 * t990 + t1057 * t991;
t943 = (-t1025 * t999 - t970) * pkin(11) + (t1000 * t1025 - t969) * pkin(5) + t950;
t937 = -t1061 * t940 + t1065 * t943;
t980 = -t1000 * t1061 + t1025 * t1065;
t953 = qJD(6) * t980 + t1012 * t1061 + t1065 * t970;
t981 = t1000 * t1065 + t1025 * t1061;
t966 = -mrSges(7,1) * t980 + mrSges(7,2) * t981;
t968 = qJDD(6) - t969;
t998 = qJD(6) - t999;
t971 = -mrSges(7,2) * t998 + mrSges(7,3) * t980;
t934 = m(7) * t937 + mrSges(7,1) * t968 - mrSges(7,3) * t953 - t966 * t981 + t971 * t998;
t938 = t1061 * t943 + t1065 * t940;
t952 = -qJD(6) * t981 + t1012 * t1065 - t1061 * t970;
t972 = mrSges(7,1) * t998 - mrSges(7,3) * t981;
t935 = m(7) * t938 - mrSges(7,2) * t968 + mrSges(7,3) * t952 + t966 * t980 - t972 * t998;
t926 = -t1061 * t934 + t1065 * t935;
t977 = -mrSges(6,1) * t999 + mrSges(6,2) * t1000;
t983 = mrSges(6,1) * t1025 - mrSges(6,3) * t1000;
t922 = m(6) * t942 - mrSges(6,2) * t1012 + mrSges(6,3) * t969 - t1025 * t983 + t977 * t999 + t926;
t1088 = t1053 * t947 - t1057 * t945;
t939 = -t1012 * pkin(5) - t1024 * pkin(11) + (t1115 + t978) * t1000 + t1088;
t936 = -m(7) * t939 + t952 * mrSges(7,1) - mrSges(7,2) * t953 + t980 * t971 - t972 * t981;
t941 = -0.2e1 * qJD(5) * t1000 - t1088;
t982 = -mrSges(6,2) * t1025 + mrSges(6,3) * t999;
t930 = m(6) * t941 + mrSges(6,1) * t1012 - mrSges(6,3) * t970 - t1000 * t977 + t1025 * t982 + t936;
t917 = t1053 * t922 + t1057 * t930;
t960 = Ifges(7,5) * t981 + Ifges(7,6) * t980 + Ifges(7,3) * t998;
t962 = Ifges(7,1) * t981 + Ifges(7,4) * t980 + Ifges(7,5) * t998;
t927 = -mrSges(7,1) * t939 + mrSges(7,3) * t938 + Ifges(7,4) * t953 + Ifges(7,2) * t952 + Ifges(7,6) * t968 - t960 * t981 + t962 * t998;
t961 = Ifges(7,4) * t981 + Ifges(7,2) * t980 + Ifges(7,6) * t998;
t928 = mrSges(7,2) * t939 - mrSges(7,3) * t937 + Ifges(7,1) * t953 + Ifges(7,4) * t952 + Ifges(7,5) * t968 + t960 * t980 - t961 * t998;
t974 = Ifges(6,4) * t1000 + Ifges(6,2) * t999 + Ifges(6,6) * t1025;
t975 = Ifges(6,1) * t1000 + Ifges(6,4) * t999 + Ifges(6,5) * t1025;
t985 = Ifges(5,4) * t1021 + Ifges(5,2) * t1020 + Ifges(5,6) * t1025;
t986 = Ifges(5,1) * t1021 + Ifges(5,4) * t1020 + Ifges(5,5) * t1025;
t1116 = Ifges(5,5) * t991 + Ifges(5,6) * t990 + t1021 * t985 - t1020 * t986 + mrSges(5,1) * t948 - mrSges(5,2) * t949 + Ifges(6,5) * t970 + Ifges(6,6) * t969 + t1000 * t974 - t999 * t975 + mrSges(6,1) * t941 - mrSges(6,2) * t942 + t1061 * t928 + t1065 * t927 + pkin(5) * t936 + pkin(11) * t926 + pkin(4) * t917 + (Ifges(5,3) + Ifges(6,3)) * t1012;
t1112 = Ifges(3,3) * t1060;
t1017 = -t1054 * t1042 + t1086;
t1090 = -mrSges(3,1) * t1058 + mrSges(3,2) * t1054;
t1040 = t1090 * t1108;
t1083 = -mrSges(3,2) * t1060 + mrSges(3,3) * t1102;
t1045 = t1083 * qJD(1);
t1106 = t1054 * t1056;
t1084 = mrSges(3,1) * t1060 - mrSges(3,3) * t1106;
t1013 = -mrSges(4,1) * t1026 + mrSges(4,2) * t1027;
t1023 = mrSges(4,1) * t1039 - mrSges(4,3) * t1027;
t1001 = -mrSges(5,1) * t1020 + mrSges(5,2) * t1021;
t1006 = -mrSges(5,2) * t1025 + mrSges(5,3) * t1020;
t915 = m(5) * t948 + mrSges(5,1) * t1012 - mrSges(5,3) * t991 - t1001 * t1021 + t1006 * t1025 + t917;
t1008 = mrSges(5,1) * t1025 - mrSges(5,3) * t1021;
t1093 = -t1053 * t930 + t1057 * t922;
t916 = m(5) * t949 - mrSges(5,2) * t1012 + mrSges(5,3) * t990 + t1001 * t1020 - t1008 * t1025 + t1093;
t1092 = -t1062 * t915 + t1066 * t916;
t906 = m(4) * t965 - mrSges(4,2) * t1036 + mrSges(4,3) * t1015 + t1013 * t1026 - t1023 * t1039 + t1092;
t1022 = -mrSges(4,2) * t1039 + mrSges(4,3) * t1026;
t909 = t1062 * t916 + t1066 * t915;
t908 = m(4) * t979 - mrSges(4,1) * t1015 + mrSges(4,2) * t1016 - t1022 * t1026 + t1023 * t1027 + t909;
t925 = t1061 * t935 + t1065 * t934;
t924 = m(6) * t950 - t969 * mrSges(6,1) + mrSges(6,2) * t970 + t1000 * t983 - t999 * t982 + t925;
t1071 = -m(5) * t955 + t990 * mrSges(5,1) - mrSges(5,2) * t991 + t1020 * t1006 - t1008 * t1021 - t924;
t923 = m(4) * t964 + mrSges(4,1) * t1036 - mrSges(4,3) * t1016 - t1013 * t1027 + t1022 * t1039 + t1071;
t895 = -t1055 * t908 + t923 * t1098 + t906 * t1099;
t891 = m(3) * t1017 + t1084 * qJDD(1) + (-t1040 * t1106 + t1045 * t1060) * qJD(1) + t895;
t1028 = -t1056 * t1041 + t1094;
t1044 = t1084 * qJD(1);
t894 = t1059 * t908 + t923 * t1103 + t906 * t1104;
t893 = m(3) * t1028 + (t1090 * qJDD(1) + (t1044 * t1054 - t1045 * t1058) * qJD(1)) * t1056 + t894;
t1018 = -g(3) * t1106 + t1095;
t902 = -t1063 * t923 + t1067 * t906;
t901 = m(3) * t1018 + t1083 * qJDD(1) + (t1040 * t1102 - t1044 * t1060) * qJD(1) + t902;
t880 = -t1056 * t893 + t891 * t1100 + t901 * t1105;
t877 = m(2) * t1049 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1069 + t880;
t887 = -t1054 * t891 + t1058 * t901;
t885 = m(2) * t1050 - mrSges(2,1) * t1069 - qJDD(1) * mrSges(2,2) + t887;
t1110 = t1064 * t885 + t1068 * t877;
t879 = t1060 * t893 + t891 * t1102 + t901 * t1106;
t1091 = -t1064 * t877 + t1068 * t885;
t1089 = Ifges(3,5) * t1054 + Ifges(3,6) * t1058;
t1009 = Ifges(4,5) * t1027 + Ifges(4,6) * t1026 + Ifges(4,3) * t1039;
t1010 = Ifges(4,4) * t1027 + Ifges(4,2) * t1026 + Ifges(4,6) * t1039;
t973 = Ifges(6,5) * t1000 + Ifges(6,6) * t999 + Ifges(6,3) * t1025;
t910 = mrSges(6,2) * t950 - mrSges(6,3) * t941 + Ifges(6,1) * t970 + Ifges(6,4) * t969 + Ifges(6,5) * t1012 - pkin(11) * t925 - t1025 * t974 - t1061 * t927 + t1065 * t928 + t973 * t999;
t1072 = mrSges(7,1) * t937 - mrSges(7,2) * t938 + Ifges(7,5) * t953 + Ifges(7,6) * t952 + Ifges(7,3) * t968 + t961 * t981 - t962 * t980;
t911 = -mrSges(6,1) * t950 + mrSges(6,3) * t942 + Ifges(6,4) * t970 + Ifges(6,2) * t969 + Ifges(6,6) * t1012 - pkin(5) * t925 - t1000 * t973 + t1025 * t975 - t1072;
t984 = Ifges(5,5) * t1021 + Ifges(5,6) * t1020 + Ifges(5,3) * t1025;
t896 = -mrSges(5,1) * t955 + mrSges(5,3) * t949 + Ifges(5,4) * t991 + Ifges(5,2) * t990 + Ifges(5,6) * t1012 - pkin(4) * t924 + qJ(5) * t1093 - t1021 * t984 + t1025 * t986 + t1053 * t910 + t1057 * t911;
t897 = mrSges(5,2) * t955 - mrSges(5,3) * t948 + Ifges(5,1) * t991 + Ifges(5,4) * t990 + Ifges(5,5) * t1012 - qJ(5) * t917 + t1020 * t984 - t1025 * t985 - t1053 * t911 + t1057 * t910;
t882 = mrSges(4,2) * t979 - mrSges(4,3) * t964 + Ifges(4,1) * t1016 + Ifges(4,4) * t1015 + Ifges(4,5) * t1036 - pkin(10) * t909 + t1009 * t1026 - t1010 * t1039 - t1062 * t896 + t1066 * t897;
t1011 = Ifges(4,1) * t1027 + Ifges(4,4) * t1026 + Ifges(4,5) * t1039;
t888 = -mrSges(4,1) * t979 + mrSges(4,3) * t965 + Ifges(4,4) * t1016 + Ifges(4,2) * t1015 + Ifges(4,6) * t1036 - pkin(3) * t909 - t1027 * t1009 + t1039 * t1011 - t1116;
t1079 = pkin(9) * t902 + t1063 * t882 + t1067 * t888;
t1076 = Ifges(3,6) * t1060 + (Ifges(3,4) * t1054 + Ifges(3,2) * t1058) * t1056;
t1032 = t1076 * qJD(1);
t1077 = Ifges(3,5) * t1060 + (Ifges(3,1) * t1054 + Ifges(3,4) * t1058) * t1056;
t1033 = t1077 * qJD(1);
t881 = mrSges(4,1) * t964 - mrSges(4,2) * t965 + Ifges(4,5) * t1016 + Ifges(4,6) * t1015 + Ifges(4,3) * t1036 + pkin(3) * t1071 + pkin(10) * t1092 + t1027 * t1010 - t1026 * t1011 + t1062 * t897 + t1066 * t896;
t871 = qJDD(1) * t1112 + mrSges(3,1) * t1017 - mrSges(3,2) * t1018 + pkin(2) * t895 + t1059 * t881 + t1079 * t1055 + (t1089 * qJDD(1) + (t1032 * t1054 - t1033 * t1058) * qJD(1)) * t1056;
t1031 = (t1089 * t1056 + t1112) * qJD(1);
t873 = -mrSges(3,1) * t1028 + mrSges(3,3) * t1018 - pkin(2) * t894 - t1055 * t881 + (-t1031 * t1106 + t1033 * t1060) * qJD(1) + t1079 * t1059 + t1076 * qJDD(1);
t875 = mrSges(3,2) * t1028 - mrSges(3,3) * t1017 - t1063 * t888 + t1067 * t882 + (t1031 * t1102 - t1032 * t1060) * qJD(1) + (-t1055 * t894 - t1059 * t895) * pkin(9) + t1077 * qJDD(1);
t1078 = mrSges(2,1) * t1049 - mrSges(2,2) * t1050 + Ifges(2,3) * qJDD(1) + pkin(1) * t880 + t1060 * t871 + t873 * t1102 + t875 * t1106 + t887 * t1109;
t869 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1049 + Ifges(2,5) * qJDD(1) - t1069 * Ifges(2,6) - t1054 * t873 + t1058 * t875 + (-t1056 * t879 - t1060 * t880) * qJ(2);
t868 = mrSges(2,1) * g(3) + mrSges(2,3) * t1050 + t1069 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t879 - t1056 * t871 + (qJ(2) * t887 + t1054 * t875 + t1058 * t873) * t1060;
t1 = [-m(1) * g(1) + t1091; -m(1) * g(2) + t1110; (-m(1) - m(2)) * g(3) + t879; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1110 - t1064 * t868 + t1068 * t869; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1091 + t1064 * t869 + t1068 * t868; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1078; t1078; t893; t881; t1116; t924; t1072;];
tauJB  = t1;
