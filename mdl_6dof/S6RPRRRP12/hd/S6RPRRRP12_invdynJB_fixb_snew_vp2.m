% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 02:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRP12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:18:21
% EndTime: 2019-05-06 02:19:25
% DurationCPUTime: 64.99s
% Computational Cost: add. (966331->381), mult. (2998350->510), div. (0->0), fcn. (2549420->14), ass. (0->177)
t1038 = sin(pkin(12));
t1040 = sin(pkin(6));
t1041 = cos(pkin(12));
t1043 = cos(pkin(6));
t1046 = sin(qJ(3));
t1042 = cos(pkin(7));
t1049 = cos(qJ(3));
t1080 = t1042 * t1049;
t1039 = sin(pkin(7));
t1085 = t1039 * t1049;
t1055 = t1040 * (-t1038 * t1046 + t1041 * t1080) + t1043 * t1085;
t1011 = t1055 * qJD(1);
t1081 = t1042 * t1046;
t1086 = t1039 * t1046;
t1057 = t1043 * t1086 + (t1038 * t1049 + t1041 * t1081) * t1040;
t1012 = t1057 * qJD(1);
t998 = -t1012 * qJD(3) + qJDD(1) * t1055;
t1109 = Ifges(6,1) + Ifges(7,1);
t1103 = Ifges(6,4) - Ifges(7,5);
t1102 = -Ifges(6,5) - Ifges(7,4);
t1108 = Ifges(6,2) + Ifges(7,3);
t1101 = Ifges(6,6) - Ifges(7,6);
t1107 = -Ifges(6,3) - Ifges(7,2);
t1083 = t1040 * t1042;
t1023 = (t1039 * t1043 + t1041 * t1083) * qJD(1) * pkin(9);
t1047 = sin(qJ(1));
t1050 = cos(qJ(1));
t1035 = -t1050 * g(1) - t1047 * g(2);
t1051 = qJD(1) ^ 2;
t1091 = qJ(2) * t1040;
t1027 = -t1051 * pkin(1) + qJDD(1) * t1091 + t1035;
t1099 = pkin(9) * t1038;
t1067 = -pkin(2) * t1041 - t1039 * t1099;
t1090 = qJD(1) * t1040;
t1097 = pkin(9) * qJDD(1);
t1064 = qJD(1) * t1067 * t1090 + t1042 * t1097;
t1034 = t1047 * g(1) - t1050 * g(2);
t1026 = qJDD(1) * pkin(1) + t1051 * t1091 + t1034;
t1078 = qJD(2) * t1090;
t1082 = t1041 * t1043;
t1084 = t1040 * t1041;
t1068 = -g(3) * t1084 + t1026 * t1082 - 0.2e1 * t1038 * t1078;
t979 = (pkin(2) * qJDD(1) + qJD(1) * t1023) * t1043 + (-t1040 * t1064 - t1027) * t1038 + t1068;
t1028 = (pkin(2) * t1043 - t1083 * t1099) * qJD(1);
t1087 = t1038 * t1043;
t1076 = t1026 * t1087 + (t1027 + 0.2e1 * t1078) * t1041;
t980 = (-qJD(1) * t1028 + t1039 * t1097) * t1043 + (-g(3) * t1038 + t1041 * t1064) * t1040 + t1076;
t1075 = -t1043 * g(3) + qJDD(2);
t988 = (-t1026 + t1067 * qJDD(1) + (-t1023 * t1041 + t1028 * t1038) * qJD(1)) * t1040 + t1075;
t941 = -t1046 * t980 + (t1039 * t988 + t1042 * t979) * t1049;
t1062 = -t1039 * t1084 + t1042 * t1043;
t1024 = qJD(1) * t1062 + qJD(3);
t1045 = sin(qJ(4));
t1048 = cos(qJ(4));
t1004 = -t1045 * t1012 + t1048 * t1024;
t1001 = qJD(5) - t1004;
t1000 = t1001 ^ 2;
t1044 = sin(qJ(5));
t1105 = cos(qJ(5));
t1010 = qJD(4) - t1011;
t1009 = t1010 ^ 2;
t1020 = t1024 ^ 2;
t1021 = qJDD(1) * t1062 + qJDD(3);
t942 = t1049 * t980 + t979 * t1081 + t988 * t1086;
t997 = -t1011 * pkin(3) - t1012 * pkin(10);
t938 = -t1020 * pkin(3) + t1021 * pkin(10) + t1011 * t997 + t942;
t955 = -t1039 * t979 + t1042 * t988;
t999 = t1011 * qJD(3) + qJDD(1) * t1057;
t940 = (-t1011 * t1024 - t999) * pkin(10) + (t1012 * t1024 - t998) * pkin(3) + t955;
t934 = t1045 * t940 + t1048 * t938;
t1005 = t1048 * t1012 + t1045 * t1024;
t982 = -t1004 * pkin(4) - t1005 * pkin(11);
t995 = qJDD(4) - t998;
t930 = -t1009 * pkin(4) + t995 * pkin(11) + t1004 * t982 + t934;
t937 = -t1021 * pkin(3) - t1020 * pkin(10) + t1012 * t997 - t941;
t974 = -t1005 * qJD(4) + t1048 * t1021 - t1045 * t999;
t975 = t1004 * qJD(4) + t1045 * t1021 + t1048 * t999;
t932 = (-t1004 * t1010 - t975) * pkin(11) + (t1005 * t1010 - t974) * pkin(4) + t937;
t927 = t1044 * t932 + t1105 * t930;
t986 = t1044 * t1005 - t1010 * t1105;
t987 = t1005 * t1105 + t1044 * t1010;
t958 = t986 * pkin(5) - t987 * qJ(6);
t972 = qJDD(5) - t974;
t923 = -t1000 * pkin(5) + t972 * qJ(6) + 0.2e1 * qJD(6) * t1001 - t986 * t958 + t927;
t965 = -t1001 * mrSges(7,1) + t987 * mrSges(7,2);
t1079 = m(7) * t923 + t972 * mrSges(7,3) + t1001 * t965;
t1093 = t1102 * t1001 + t1103 * t986 - t1109 * t987;
t1094 = -t1001 * t1101 - t1103 * t987 + t1108 * t986;
t926 = -t1044 * t930 + t1105 * t932;
t924 = -t972 * pkin(5) - t1000 * qJ(6) + t987 * t958 + qJDD(6) - t926;
t962 = -t986 * mrSges(7,2) + t1001 * mrSges(7,3);
t1072 = -m(7) * t924 + t972 * mrSges(7,1) + t1001 * t962;
t946 = -t986 * qJD(5) + t1044 * t995 + t1105 * t975;
t959 = t986 * mrSges(7,1) - t987 * mrSges(7,3);
t920 = t946 * mrSges(7,2) + t987 * t959 - t1072;
t945 = t987 * qJD(5) + t1044 * t975 - t1105 * t995;
t1106 = -t1093 * t986 - t1094 * t987 - t1107 * t972 - t1101 * t945 - t1102 * t946 + mrSges(6,1) * t926 - mrSges(7,1) * t924 - mrSges(6,2) * t927 + mrSges(7,3) * t923 - pkin(5) * t920 + qJ(6) * (-t945 * mrSges(7,2) - t986 * t959 + t1079);
t1104 = -mrSges(6,3) - mrSges(7,2);
t1098 = Ifges(3,3) * t1043;
t1002 = -t1038 * t1027 + t1068;
t1071 = -mrSges(3,1) * t1041 + mrSges(3,2) * t1038;
t1025 = t1071 * t1090;
t1065 = -mrSges(3,2) * t1043 + mrSges(3,3) * t1084;
t1030 = t1065 * qJD(1);
t1088 = t1038 * t1040;
t1066 = mrSges(3,1) * t1043 - mrSges(3,3) * t1088;
t1007 = t1024 * mrSges(4,1) - t1012 * mrSges(4,3);
t1092 = -t986 * mrSges(6,1) - t987 * mrSges(6,2) - t959;
t964 = t1001 * mrSges(6,1) - t987 * mrSges(6,3);
t916 = m(6) * t927 - t972 * mrSges(6,2) - t1001 * t964 + t1092 * t986 + t1104 * t945 + t1079;
t963 = -t1001 * mrSges(6,2) - t986 * mrSges(6,3);
t917 = m(6) * t926 + t972 * mrSges(6,1) + t1001 * t963 + t1092 * t987 + t1104 * t946 + t1072;
t912 = -t1044 * t917 + t1105 * t916;
t981 = -t1004 * mrSges(5,1) + t1005 * mrSges(5,2);
t990 = t1010 * mrSges(5,1) - t1005 * mrSges(5,3);
t908 = m(5) * t934 - t995 * mrSges(5,2) + t974 * mrSges(5,3) + t1004 * t981 - t1010 * t990 + t912;
t933 = -t1045 * t938 + t1048 * t940;
t929 = -t995 * pkin(4) - t1009 * pkin(11) + t1005 * t982 - t933;
t925 = -0.2e1 * qJD(6) * t987 + (t1001 * t986 - t946) * qJ(6) + (t1001 * t987 + t945) * pkin(5) + t929;
t921 = m(7) * t925 + t945 * mrSges(7,1) - t946 * mrSges(7,3) + t986 * t962 - t987 * t965;
t918 = -m(6) * t929 - t945 * mrSges(6,1) - t946 * mrSges(6,2) - t986 * t963 - t987 * t964 - t921;
t989 = -t1010 * mrSges(5,2) + t1004 * mrSges(5,3);
t914 = m(5) * t933 + t995 * mrSges(5,1) - t975 * mrSges(5,3) - t1005 * t981 + t1010 * t989 + t918;
t1074 = -t1045 * t914 + t1048 * t908;
t996 = -t1011 * mrSges(4,1) + t1012 * mrSges(4,2);
t899 = m(4) * t942 - t1021 * mrSges(4,2) + t998 * mrSges(4,3) - t1024 * t1007 + t1011 * t996 + t1074;
t1006 = -t1024 * mrSges(4,2) + t1011 * mrSges(4,3);
t902 = t1045 * t908 + t1048 * t914;
t901 = m(4) * t955 - t998 * mrSges(4,1) + t999 * mrSges(4,2) - t1011 * t1006 + t1012 * t1007 + t902;
t911 = t1044 * t916 + t1105 * t917;
t1053 = -m(5) * t937 + t974 * mrSges(5,1) - t975 * mrSges(5,2) + t1004 * t989 - t1005 * t990 - t911;
t905 = m(4) * t941 + t1021 * mrSges(4,1) - t999 * mrSges(4,3) + t1024 * t1006 - t1012 * t996 + t1053;
t888 = -t1039 * t901 + t905 * t1080 + t899 * t1081;
t884 = m(3) * t1002 + t1066 * qJDD(1) + (-t1025 * t1088 + t1030 * t1043) * qJD(1) + t888;
t1013 = -t1040 * t1026 + t1075;
t1029 = t1066 * qJD(1);
t887 = t1042 * t901 + t905 * t1085 + t899 * t1086;
t886 = m(3) * t1013 + (t1071 * qJDD(1) + (t1029 * t1038 - t1030 * t1041) * qJD(1)) * t1040 + t887;
t1003 = -g(3) * t1088 + t1076;
t894 = -t1046 * t905 + t1049 * t899;
t893 = m(3) * t1003 + t1065 * qJDD(1) + (t1025 * t1084 - t1029 * t1043) * qJD(1) + t894;
t873 = -t1040 * t886 + t884 * t1082 + t893 * t1087;
t870 = m(2) * t1034 + qJDD(1) * mrSges(2,1) - t1051 * mrSges(2,2) + t873;
t880 = -t1038 * t884 + t1041 * t893;
t878 = m(2) * t1035 - t1051 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t880;
t1096 = t1047 * t878 + t1050 * t870;
t1095 = t1001 * t1107 + t1101 * t986 + t1102 * t987;
t872 = t1043 * t886 + t884 * t1084 + t893 * t1088;
t1073 = -t1047 * t870 + t1050 * t878;
t1070 = Ifges(3,5) * t1038 + Ifges(3,6) * t1041;
t909 = -mrSges(6,1) * t929 - mrSges(7,1) * t925 + mrSges(7,2) * t923 + mrSges(6,3) * t927 - pkin(5) * t921 - t1093 * t1001 + t1095 * t987 + t1101 * t972 + t1103 * t946 - t1108 * t945;
t910 = mrSges(6,2) * t929 + mrSges(7,2) * t924 - mrSges(6,3) * t926 - mrSges(7,3) * t925 - qJ(6) * t921 + t1094 * t1001 + t1095 * t986 - t1102 * t972 - t1103 * t945 + t1109 * t946;
t968 = Ifges(5,5) * t1005 + Ifges(5,6) * t1004 + Ifges(5,3) * t1010;
t969 = Ifges(5,4) * t1005 + Ifges(5,2) * t1004 + Ifges(5,6) * t1010;
t889 = mrSges(5,2) * t937 - mrSges(5,3) * t933 + Ifges(5,1) * t975 + Ifges(5,4) * t974 + Ifges(5,5) * t995 - pkin(11) * t911 + t1004 * t968 - t1010 * t969 - t1044 * t909 + t1105 * t910;
t970 = Ifges(5,1) * t1005 + Ifges(5,4) * t1004 + Ifges(5,5) * t1010;
t895 = -mrSges(5,1) * t937 + mrSges(5,3) * t934 + Ifges(5,4) * t975 + Ifges(5,2) * t974 + Ifges(5,6) * t995 - pkin(4) * t911 - t1005 * t968 + t1010 * t970 - t1106;
t991 = Ifges(4,5) * t1012 + Ifges(4,6) * t1011 + Ifges(4,3) * t1024;
t992 = Ifges(4,4) * t1012 + Ifges(4,2) * t1011 + Ifges(4,6) * t1024;
t875 = mrSges(4,2) * t955 - mrSges(4,3) * t941 + Ifges(4,1) * t999 + Ifges(4,4) * t998 + Ifges(4,5) * t1021 - pkin(10) * t902 + t1011 * t991 - t1024 * t992 - t1045 * t895 + t1048 * t889;
t1052 = mrSges(5,1) * t933 - mrSges(5,2) * t934 + Ifges(5,5) * t975 + Ifges(5,6) * t974 + Ifges(5,3) * t995 + pkin(4) * t918 + pkin(11) * t912 - t1004 * t970 + t1005 * t969 + t1044 * t910 + t1105 * t909;
t993 = Ifges(4,1) * t1012 + Ifges(4,4) * t1011 + Ifges(4,5) * t1024;
t881 = -mrSges(4,1) * t955 + mrSges(4,3) * t942 + Ifges(4,4) * t999 + Ifges(4,2) * t998 + Ifges(4,6) * t1021 - pkin(3) * t902 - t1012 * t991 + t1024 * t993 - t1052;
t1061 = pkin(9) * t894 + t1046 * t875 + t1049 * t881;
t1058 = Ifges(3,6) * t1043 + (Ifges(3,4) * t1038 + Ifges(3,2) * t1041) * t1040;
t1017 = t1058 * qJD(1);
t1059 = Ifges(3,5) * t1043 + (Ifges(3,1) * t1038 + Ifges(3,4) * t1041) * t1040;
t1018 = t1059 * qJD(1);
t874 = mrSges(4,1) * t941 - mrSges(4,2) * t942 + Ifges(4,5) * t999 + Ifges(4,6) * t998 + Ifges(4,3) * t1021 + pkin(3) * t1053 + pkin(10) * t1074 - t1011 * t993 + t1012 * t992 + t1045 * t889 + t1048 * t895;
t864 = qJDD(1) * t1098 + mrSges(3,1) * t1002 - mrSges(3,2) * t1003 + pkin(2) * t888 + t1042 * t874 + t1061 * t1039 + (t1070 * qJDD(1) + (t1017 * t1038 - t1018 * t1041) * qJD(1)) * t1040;
t1016 = (t1040 * t1070 + t1098) * qJD(1);
t866 = -mrSges(3,1) * t1013 + mrSges(3,3) * t1003 - pkin(2) * t887 - t1039 * t874 + (-t1016 * t1088 + t1018 * t1043) * qJD(1) + t1061 * t1042 + t1058 * qJDD(1);
t868 = mrSges(3,2) * t1013 - mrSges(3,3) * t1002 - t1046 * t881 + t1049 * t875 + (t1016 * t1084 - t1017 * t1043) * qJD(1) + (-t1039 * t887 - t1042 * t888) * pkin(9) + t1059 * qJDD(1);
t1060 = mrSges(2,1) * t1034 - mrSges(2,2) * t1035 + Ifges(2,3) * qJDD(1) + pkin(1) * t873 + t1043 * t864 + t866 * t1084 + t868 * t1088 + t880 * t1091;
t862 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1034 + Ifges(2,5) * qJDD(1) - t1051 * Ifges(2,6) - t1038 * t866 + t1041 * t868 + (-t1040 * t872 - t1043 * t873) * qJ(2);
t861 = mrSges(2,1) * g(3) + mrSges(2,3) * t1035 + t1051 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t872 - t1040 * t864 + (qJ(2) * t880 + t1038 * t868 + t1041 * t866) * t1043;
t1 = [-m(1) * g(1) + t1073; -m(1) * g(2) + t1096; (-m(1) - m(2)) * g(3) + t872; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1096 - t1047 * t861 + t1050 * t862; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1073 + t1047 * t862 + t1050 * t861; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1060; t1060; t886; t874; t1052; t1106; t920;];
tauJB  = t1;
