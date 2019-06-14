% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-05-06 09:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:14:45
% EndTime: 2019-05-06 09:14:52
% DurationCPUTime: 3.87s
% Computational Cost: add. (25869->331), mult. (53405->370), div. (0->0), fcn. (26316->6), ass. (0->134)
t1103 = Ifges(6,4) + Ifges(7,4);
t1127 = Ifges(6,2) + Ifges(7,2);
t1124 = Ifges(6,6) + Ifges(7,6);
t1126 = Ifges(6,1) + Ifges(7,1);
t1125 = Ifges(6,5) + Ifges(7,5);
t1123 = Ifges(6,3) + Ifges(7,3);
t1053 = sin(qJ(5));
t1056 = cos(qJ(5));
t1057 = cos(qJ(2));
t1094 = qJD(1) * t1057;
t1010 = -t1056 * qJD(2) + t1053 * t1094;
t1011 = t1053 * qJD(2) + t1056 * t1094;
t1054 = sin(qJ(2));
t1091 = t1054 * qJD(1);
t1035 = qJD(5) + t1091;
t1122 = -t1127 * t1010 + t1103 * t1011 - t1124 * t1035;
t1121 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t1088 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t1087 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t1086 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t1120 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t1119 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t1076 = qJD(2) * t1094;
t1017 = t1054 * qJDD(1) + t1076;
t1077 = qJD(2) * t1091;
t1018 = t1057 * qJDD(1) - t1077;
t1013 = (-mrSges(4,1) * t1057 - mrSges(4,3) * t1054) * qJD(1);
t1027 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t1091;
t1025 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t1091;
t1007 = t1010 ^ 2;
t1016 = (pkin(4) * t1054 + pkin(8) * t1057) * qJD(1);
t1093 = qJD(3) * qJD(2);
t1040 = 0.2e1 * t1093;
t1059 = qJD(2) ^ 2;
t1108 = -pkin(2) - pkin(8);
t1024 = -qJD(2) * pkin(3) - qJ(4) * t1091;
t1060 = qJD(1) ^ 2;
t1092 = t1057 ^ 2 * t1060;
t1012 = (-pkin(2) * t1057 - qJ(3) * t1054) * qJD(1);
t1055 = sin(qJ(1));
t1058 = cos(qJ(1));
t1032 = -t1058 * g(1) - t1055 * g(2);
t990 = -t1060 * pkin(1) + qJDD(1) * pkin(7) + t1032;
t973 = -t1054 * g(3) + t1057 * t990;
t1115 = qJDD(2) * qJ(3) + t1012 * t1094 + t973;
t1117 = pkin(3) * t1092 + t1018 * qJ(4) - qJD(2) * t1024 + 0.2e1 * qJD(4) * t1094 - t1115;
t934 = qJDD(2) * pkin(4) - t1016 * t1094 + t1108 * t1059 + t1040 - t1117;
t962 = t1011 * qJD(5) - t1056 * qJDD(2) + t1053 * t1018;
t969 = t1035 * pkin(5) + t1011 * qJ(6);
t928 = -t962 * pkin(5) - t1007 * qJ(6) - t1011 * t969 + qJDD(6) + t934;
t963 = t1010 * qJD(5) - t1053 * qJDD(2) - t1056 * t1018;
t970 = t1035 * mrSges(7,1) + t1011 * mrSges(7,3);
t1082 = -m(7) * t928 - t963 * mrSges(7,2) + t1011 * t970;
t971 = t1035 * mrSges(6,1) + t1011 * mrSges(6,3);
t1072 = -m(6) * t934 - t963 * mrSges(6,2) + t1011 * t971 + t1082;
t1102 = t1059 * pkin(2);
t938 = -0.2e1 * t1093 + t1102 + t1117;
t1067 = -m(5) * t938 + qJDD(2) * mrSges(5,1) - t1018 * mrSges(5,3) + qJD(2) * t1025 - t1072;
t967 = -t1035 * mrSges(7,2) + t1010 * mrSges(7,3);
t968 = -t1035 * mrSges(6,2) + t1010 * mrSges(6,3);
t1112 = -(-t968 - t967) * t1010 + (mrSges(6,1) + mrSges(7,1)) * t962;
t944 = t1040 - t1102 + t1115;
t1063 = m(4) * t944 + qJDD(2) * mrSges(4,3) + qJD(2) * t1027 + t1013 * t1094 + t1067 - t1112;
t1079 = -t1087 * qJD(2) + (-t1121 * t1054 - t1088 * t1057) * qJD(1);
t1080 = -t1086 * qJD(2) + (-t1088 * t1054 - t1119 * t1057) * qJD(1);
t1015 = (mrSges(5,1) * t1054 - mrSges(5,2) * t1057) * qJD(1);
t1090 = t1057 * t1015;
t1008 = qJDD(5) + t1017;
t1109 = 2 * qJD(6);
t1031 = t1055 * g(1) - t1058 * g(2);
t989 = -qJDD(1) * pkin(1) - t1060 * pkin(7) - t1031;
t1070 = -t1018 * pkin(2) + t989 + (-t1017 - t1076) * qJ(3);
t1065 = -qJ(4) * t1092 + qJDD(4) - t1070 + (0.2e1 * qJD(3) + t1024) * t1091;
t1107 = pkin(3) + pkin(8);
t931 = t1065 + (pkin(4) * t1057 + t1108 * t1054) * qJD(2) * qJD(1) + t1107 * t1018 + t1017 * pkin(4);
t972 = -t1057 * g(3) - t1054 * t990;
t1073 = t1012 * t1091 + qJDD(3) - t972;
t1089 = t1057 * t1060;
t1114 = -0.2e1 * qJD(4) * t1091 + (-t1017 + t1076) * qJ(4);
t935 = (-pkin(4) - qJ(3)) * t1059 + (-pkin(3) * t1089 - qJD(1) * t1016) * t1054 + (-pkin(2) - t1107) * qJDD(2) + t1073 + t1114;
t926 = t1053 * t931 + t1056 * t935;
t922 = -t1007 * pkin(5) + t962 * qJ(6) + t1010 * t1109 - t1035 * t969 + t926;
t965 = -t1010 * mrSges(7,1) - t1011 * mrSges(7,2);
t1083 = m(7) * t922 + t962 * mrSges(7,3) + t1010 * t965;
t1098 = -t1103 * t1010 + t1126 * t1011 - t1125 * t1035;
t1099 = t1124 * t1010 - t1125 * t1011 + t1123 * t1035;
t1116 = -t962 * mrSges(7,1) - t1010 * t967;
t923 = -t1082 + t1116;
t893 = -mrSges(6,1) * t934 + mrSges(6,3) * t926 - mrSges(7,1) * t928 + mrSges(7,3) * t922 - pkin(5) * t923 + qJ(6) * t1083 + t1103 * t963 + t1127 * t962 + (-qJ(6) * t970 - t1098) * t1035 + t1099 * t1011 + (-qJ(6) * mrSges(7,2) + t1124) * t1008;
t1028 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t1094;
t925 = -t1053 * t935 + t1056 * t931;
t920 = t1011 * t1109 + (t1010 * t1035 - t963) * qJ(6) + (-t1010 * t1011 + t1008) * pkin(5) + t925;
t1084 = m(7) * t920 + t1008 * mrSges(7,1) + t1035 * t967;
t966 = -t1010 * mrSges(6,1) - t1011 * mrSges(6,2);
t911 = m(6) * t925 + t1008 * mrSges(6,1) + t1035 * t968 + (-mrSges(6,3) - mrSges(7,3)) * t963 + (t965 + t966) * t1011 + t1084;
t913 = m(6) * t926 + t962 * mrSges(6,3) + t1010 * t966 + (-t970 - t971) * t1035 + (-mrSges(6,2) - mrSges(7,2)) * t1008 + t1083;
t905 = -t1053 * t911 + t1056 * t913;
t945 = -qJDD(2) * pkin(2) - t1059 * qJ(3) + t1073;
t939 = (-t1054 * t1089 - qJDD(2)) * pkin(3) + t945 + t1114;
t1071 = -m(5) * t939 + t1015 * t1091 - t905;
t1066 = qJDD(2) * mrSges(5,2) + qJD(2) * t1028 - t1071;
t1030 = mrSges(4,2) * t1094 + qJD(2) * mrSges(4,3);
t1113 = m(4) * t945 - qJDD(2) * mrSges(4,1) - qJD(2) * t1030;
t900 = t1013 * t1091 + (mrSges(4,2) - mrSges(5,3)) * t1017 + t1066 + t1113;
t902 = -t1017 * mrSges(5,3) + t1066;
t917 = -t963 * mrSges(7,3) + t1011 * t965 + t1084;
t903 = mrSges(6,2) * t934 + mrSges(7,2) * t928 - mrSges(6,3) * t925 - mrSges(7,3) * t920 - qJ(6) * t917 + t1125 * t1008 + t1099 * t1010 + t1122 * t1035 + t1103 * t962 + t1126 * t963;
t1118 = -(t1080 * t1054 - t1079 * t1057) * qJD(1) + t1120 * qJDD(2) + t1087 * t1017 + t1086 * t1018 + mrSges(3,1) * t972 - mrSges(4,1) * t945 - mrSges(5,1) * t938 - mrSges(3,2) * t973 + mrSges(5,2) * t939 + mrSges(4,3) * t944 - pkin(2) * t900 - pkin(3) * t902 - pkin(4) * (t1072 + t1112) - pkin(8) * t905 + qJ(3) * (t1018 * mrSges(4,2) - qJD(1) * t1090 + t1063) - t1053 * t903 - t1056 * t893;
t1104 = mrSges(3,3) + mrSges(4,2);
t1014 = (-mrSges(3,1) * t1057 + mrSges(3,2) * t1054) * qJD(1);
t1029 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t1094;
t898 = m(3) * t972 + (mrSges(3,1) - mrSges(5,2)) * qJDD(2) + (-t1028 + t1029) * qJD(2) + (-t1013 - t1014) * t1091 + (mrSges(5,3) - t1104) * t1017 + t1071 - t1113;
t1026 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t1091;
t908 = (t1014 - t1015) * t1094 + t1063 - qJDD(2) * mrSges(3,2) + t1104 * t1018 - qJD(2) * t1026 + m(3) * t973;
t1075 = -t1054 * t898 + t1057 * t908;
t890 = m(2) * t1032 - t1060 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t1075;
t904 = t1053 * t913 + t1056 * t911;
t937 = -pkin(2) * t1077 + t1018 * pkin(3) + t1065;
t901 = m(5) * t937 + t1017 * mrSges(5,1) - t1018 * mrSges(5,2) + t1025 * t1091 - t1028 * t1094 + t904;
t940 = (pkin(2) * qJD(2) - 0.2e1 * qJD(3)) * t1091 + t1070;
t899 = m(4) * t940 - t1018 * mrSges(4,1) - t1017 * mrSges(4,3) - t1027 * t1091 - t1030 * t1094 - t901;
t1062 = -m(3) * t989 + t1018 * mrSges(3,1) - t1017 * mrSges(3,2) - t1026 * t1091 + t1029 * t1094 - t899;
t895 = m(2) * t1031 + qJDD(1) * mrSges(2,1) - t1060 * mrSges(2,2) + t1062;
t1100 = t1055 * t890 + t1058 * t895;
t892 = t1054 * t908 + t1057 * t898;
t1081 = -t1120 * qJD(2) + (-t1087 * t1054 - t1086 * t1057) * qJD(1);
t1074 = -t1055 * t895 + t1058 * t890;
t885 = t1053 * t893 - t1056 * t903 - qJ(4) * (-t962 * mrSges(6,1) - t1010 * t968 + t1067 + t1116) - mrSges(3,1) * t989 + mrSges(3,3) * t973 - mrSges(5,2) * t937 + mrSges(5,3) * t938 - mrSges(4,1) * t940 + mrSges(4,2) * t944 + pkin(8) * t904 + pkin(3) * t901 - pkin(2) * t899 + t1119 * t1018 + t1088 * t1017 + t1086 * qJDD(2) - t1079 * qJD(2) + (qJ(4) * t1090 + t1081 * t1054) * qJD(1);
t1064 = mrSges(6,1) * t925 + mrSges(7,1) * t920 - mrSges(6,2) * t926 - mrSges(7,2) * t922 + pkin(5) * t917 + t1123 * t1008 + t1098 * t1010 + t1122 * t1011 + t1124 * t962 + t1125 * t963;
t887 = mrSges(5,1) * t937 + mrSges(3,2) * t989 + mrSges(4,2) * t945 - mrSges(3,3) * t972 - mrSges(4,3) * t940 - mrSges(5,3) * t939 + pkin(4) * t904 - qJ(3) * t899 - qJ(4) * t902 + t1080 * qJD(2) + t1087 * qJDD(2) + t1121 * t1017 + t1088 * t1018 - t1081 * t1094 + t1064;
t1069 = mrSges(2,1) * t1031 - mrSges(2,2) * t1032 + Ifges(2,3) * qJDD(1) + pkin(1) * t1062 + pkin(7) * t1075 + t1054 * t887 + t1057 * t885;
t883 = mrSges(2,1) * g(3) + mrSges(2,3) * t1032 + t1060 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t892 - t1118;
t882 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1031 + Ifges(2,5) * qJDD(1) - t1060 * Ifges(2,6) - pkin(7) * t892 - t1054 * t885 + t1057 * t887;
t1 = [-m(1) * g(1) + t1074; -m(1) * g(2) + t1100; (-m(1) - m(2)) * g(3) + t892; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t1100 - t1055 * t883 + t1058 * t882; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t1074 + t1055 * t882 + t1058 * t883; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1069; t1069; t1118; t900; t901; t1064; t923;];
tauJB  = t1;
