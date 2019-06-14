% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-05-07 09:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:08:53
% EndTime: 2019-05-07 09:09:11
% DurationCPUTime: 12.57s
% Computational Cost: add. (187080->362), mult. (401193->433), div. (0->0), fcn. (305203->10), ass. (0->150)
t1070 = Ifges(6,4) + Ifges(7,4);
t1083 = Ifges(6,2) + Ifges(7,2);
t1078 = Ifges(6,6) + Ifges(7,6);
t1082 = Ifges(4,1) + Ifges(5,2);
t1081 = Ifges(6,1) + Ifges(7,1);
t1069 = Ifges(4,5) - Ifges(5,4);
t1080 = Ifges(6,5) + Ifges(7,5);
t1079 = -Ifges(4,2) - Ifges(5,3);
t1068 = Ifges(4,6) - Ifges(5,5);
t1067 = -Ifges(5,6) - Ifges(4,4);
t1077 = Ifges(4,3) + Ifges(5,1);
t1076 = Ifges(6,3) + Ifges(7,3);
t1019 = sin(pkin(6));
t1026 = cos(qJ(2));
t1049 = t1019 * t1026;
t1042 = qJD(1) * t1049;
t1009 = -qJD(3) + t1042;
t1021 = sin(qJ(5));
t1025 = cos(qJ(5));
t1020 = cos(pkin(6));
t1015 = qJD(1) * t1020 + qJD(2);
t1022 = sin(qJ(3));
t1023 = sin(qJ(2));
t1050 = t1019 * t1023;
t1043 = qJD(1) * t1050;
t1072 = cos(qJ(3));
t989 = -t1015 * t1072 + t1022 * t1043;
t972 = t1009 * t1021 + t1025 * t989;
t973 = -t1009 * t1025 + t1021 * t989;
t990 = t1022 * t1015 + t1043 * t1072;
t987 = qJD(5) + t990;
t1075 = t1070 * t973 + t1078 * t987 + t1083 * t972;
t1003 = (qJD(1) * qJD(2) * t1023 - qJDD(1) * t1026) * t1019;
t1008 = t1009 ^ 2;
t1053 = qJD(1) * t1019;
t1001 = (-pkin(2) * t1026 - pkin(9) * t1023) * t1053;
t1013 = t1015 ^ 2;
t1014 = qJDD(1) * t1020 + qJDD(2);
t1052 = qJD(1) * t1026;
t1048 = t1020 * t1023;
t1024 = sin(qJ(1));
t1027 = cos(qJ(1));
t1011 = t1024 * g(1) - g(2) * t1027;
t1028 = qJD(1) ^ 2;
t1065 = pkin(8) * t1019;
t998 = qJDD(1) * pkin(1) + t1028 * t1065 + t1011;
t1012 = -g(1) * t1027 - g(2) * t1024;
t999 = -pkin(1) * t1028 + qJDD(1) * t1065 + t1012;
t1055 = t1026 * t999 + t998 * t1048;
t934 = -pkin(2) * t1013 + pkin(9) * t1014 + (-g(3) * t1023 + t1001 * t1052) * t1019 + t1055;
t1002 = (qJD(2) * t1052 + qJDD(1) * t1023) * t1019;
t1064 = g(3) * t1020;
t935 = pkin(2) * t1003 - pkin(9) * t1002 - t1064 + (-t998 + (pkin(2) * t1023 - pkin(9) * t1026) * t1015 * qJD(1)) * t1019;
t909 = t1022 * t935 + t1072 * t934;
t966 = pkin(3) * t989 - qJ(4) * t990;
t995 = qJDD(3) + t1003;
t905 = pkin(3) * t1008 - t995 * qJ(4) + 0.2e1 * qJD(4) * t1009 + t989 * t966 - t909;
t962 = qJD(3) * t990 + t1002 * t1022 - t1014 * t1072;
t978 = pkin(4) * t990 + pkin(10) * t1009;
t988 = t989 ^ 2;
t903 = -pkin(4) * t962 - pkin(10) * t988 - t1009 * t978 - t905;
t919 = -qJD(5) * t973 - t1021 * t995 + t1025 * t962;
t944 = pkin(5) * t987 - qJ(6) * t973;
t971 = t972 ^ 2;
t897 = -t919 * pkin(5) - qJ(6) * t971 + t944 * t973 + qJDD(6) + t903;
t920 = qJD(5) * t972 + t1021 * t962 + t1025 * t995;
t945 = mrSges(7,1) * t987 - mrSges(7,3) * t973;
t1044 = m(7) * t897 + t920 * mrSges(7,2) + t973 * t945;
t942 = -mrSges(7,2) * t987 + mrSges(7,3) * t972;
t943 = -mrSges(6,2) * t987 + mrSges(6,3) * t972;
t946 = mrSges(6,1) * t987 - mrSges(6,3) * t973;
t1074 = -m(6) * t903 - t920 * mrSges(6,2) + (t942 + t943) * t972 + (mrSges(6,1) + mrSges(7,1)) * t919 - t973 * t946 - t1044;
t977 = mrSges(5,1) * t990 - mrSges(5,2) * t1009;
t1033 = -m(5) * t905 + t995 * mrSges(5,3) - t1009 * t977 - t1074;
t1056 = -t1069 * t1009 + t1067 * t989 + t1082 * t990;
t1057 = -t1068 * t1009 - t1067 * t990 + t1079 * t989;
t1054 = t1009 * t989;
t908 = -t1022 * t934 + t1072 * t935;
t906 = -t995 * pkin(3) - t1008 * qJ(4) + t990 * t966 + qJDD(4) - t908;
t963 = -t989 * qJD(3) + t1002 * t1072 + t1022 * t1014;
t900 = (t989 * t990 - t995) * pkin(10) + (t963 - t1054) * pkin(4) + t906;
t1047 = t1020 * t1026;
t964 = -g(3) * t1049 - t1023 * t999 + t1047 * t998;
t933 = -pkin(2) * t1014 - pkin(9) * t1013 + t1001 * t1043 - t964;
t1031 = (-t963 - t1054) * qJ(4) + t933 + (-t1009 * pkin(3) - 0.2e1 * qJD(4)) * t990;
t904 = -pkin(4) * t988 - t978 * t990 + (pkin(3) + pkin(10)) * t962 + t1031;
t895 = t1021 * t900 + t1025 * t904;
t891 = -pkin(5) * t971 + t919 * qJ(6) + 0.2e1 * qJD(6) * t972 - t944 * t987 + t895;
t937 = -mrSges(7,1) * t972 + mrSges(7,2) * t973;
t1045 = m(7) * t891 + t919 * mrSges(7,3) + t972 * t937;
t1060 = -t1070 * t972 - t1080 * t987 - t1081 * t973;
t1061 = -t1076 * t987 - t1078 * t972 - t1080 * t973;
t892 = -t919 * mrSges(7,1) - t942 * t972 + t1044;
t959 = qJDD(5) + t963;
t867 = -mrSges(6,1) * t903 + mrSges(6,3) * t895 - mrSges(7,1) * t897 + mrSges(7,3) * t891 - pkin(5) * t892 + qJ(6) * t1045 + (-qJ(6) * t945 - t1060) * t987 + t1061 * t973 + (-qJ(6) * mrSges(7,2) + t1078) * t959 + t1070 * t920 + t1083 * t919;
t894 = -t1021 * t904 + t1025 * t900;
t889 = -0.2e1 * qJD(6) * t973 + (t972 * t987 - t920) * qJ(6) + (t972 * t973 + t959) * pkin(5) + t894;
t1046 = m(7) * t889 + t959 * mrSges(7,1) + t987 * t942;
t938 = -mrSges(6,1) * t972 + mrSges(6,2) * t973;
t881 = m(6) * t894 + mrSges(6,1) * t959 + t943 * t987 + (-t937 - t938) * t973 + (-mrSges(6,3) - mrSges(7,3)) * t920 + t1046;
t883 = m(6) * t895 + t919 * mrSges(6,3) + t938 * t972 + (-t945 - t946) * t987 + (-mrSges(6,2) - mrSges(7,2)) * t959 + t1045;
t876 = t1021 * t883 + t1025 * t881;
t968 = -mrSges(5,2) * t989 - mrSges(5,3) * t990;
t1035 = -m(5) * t906 - t963 * mrSges(5,1) - t990 * t968 - t876;
t976 = mrSges(5,1) * t989 + mrSges(5,3) * t1009;
t874 = mrSges(5,2) * t995 - t1009 * t976 - t1035;
t886 = -t920 * mrSges(7,3) - t937 * t973 + t1046;
t875 = mrSges(6,2) * t903 + mrSges(7,2) * t897 - mrSges(6,3) * t894 - mrSges(7,3) * t889 - qJ(6) * t886 - t1061 * t972 + t1070 * t919 - t1075 * t987 + t1080 * t959 + t1081 * t920;
t1073 = t1056 * t989 + t1057 * t990 + t1077 * t995 - t1068 * t962 + t1069 * t963 + mrSges(4,1) * t908 - mrSges(4,2) * t909 + mrSges(5,2) * t906 - mrSges(5,3) * t905 - pkin(3) * t874 - pkin(10) * t876 + qJ(4) * (-mrSges(5,1) * t962 - t968 * t989 + t1033) - t1021 * t867 + t1025 * t875;
t1000 = (-mrSges(3,1) * t1026 + mrSges(3,2) * t1023) * t1053;
t967 = mrSges(4,1) * t989 + mrSges(4,2) * t990;
t974 = mrSges(4,2) * t1009 - mrSges(4,3) * t989;
t872 = m(4) * t908 - mrSges(4,3) * t963 - t967 * t990 + (mrSges(4,1) - mrSges(5,2)) * t995 + (-t974 + t976) * t1009 + t1035;
t975 = -mrSges(4,1) * t1009 - mrSges(4,3) * t990;
t879 = t1033 + m(4) * t909 + (-t967 - t968) * t989 + (-mrSges(4,3) - mrSges(5,1)) * t962 + t1009 * t975 - mrSges(4,2) * t995;
t1041 = -t1022 * t872 + t1072 * t879;
t965 = -g(3) * t1050 + t1055;
t996 = mrSges(3,1) * t1015 - mrSges(3,3) * t1043;
t863 = m(3) * t965 - mrSges(3,2) * t1014 - mrSges(3,3) * t1003 + t1000 * t1042 - t1015 * t996 + t1041;
t866 = t1022 * t879 + t1072 * t872;
t983 = -t1019 * t998 - t1064;
t997 = -mrSges(3,2) * t1015 + mrSges(3,3) * t1042;
t865 = m(3) * t983 + mrSges(3,1) * t1003 + mrSges(3,2) * t1002 + (t1023 * t996 - t1026 * t997) * t1053 + t866;
t1062 = -t1021 * t881 + t1025 * t883;
t907 = pkin(3) * t962 + t1031;
t1038 = -m(5) * t907 + t962 * mrSges(5,2) + t989 * t976 - t1062;
t1032 = -m(4) * t933 - t962 * mrSges(4,1) - t989 * t974 + (-t975 + t977) * t990 + (-mrSges(4,2) + mrSges(5,3)) * t963 + t1038;
t870 = m(3) * t964 + mrSges(3,1) * t1014 - mrSges(3,3) * t1002 - t1000 * t1043 + t1015 * t997 + t1032;
t853 = -t1019 * t865 + t870 * t1047 + t863 * t1048;
t850 = m(2) * t1011 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1028 + t853;
t859 = -t1023 * t870 + t1026 * t863;
t857 = m(2) * t1012 - mrSges(2,1) * t1028 - qJDD(1) * mrSges(2,2) + t859;
t1063 = t1024 * t857 + t1027 * t850;
t1058 = t1077 * t1009 + t1068 * t989 - t1069 * t990;
t852 = t1020 * t865 + t870 * t1049 + t863 * t1050;
t1040 = -t1024 * t850 + t1027 * t857;
t873 = -mrSges(5,3) * t963 - t977 * t990 - t1038;
t848 = -mrSges(4,1) * t933 - mrSges(5,1) * t905 + mrSges(5,2) * t907 + mrSges(4,3) * t909 - pkin(3) * t873 - pkin(4) * t1074 - pkin(10) * t1062 - t1056 * t1009 - t1021 * t875 - t1025 * t867 + t1058 * t990 - t1067 * t963 + t1068 * t995 + t1079 * t962;
t1030 = mrSges(6,1) * t894 + mrSges(7,1) * t889 - mrSges(6,2) * t895 - mrSges(7,2) * t891 + pkin(5) * t886 + t1060 * t972 + t1075 * t973 + t1076 * t959 + t1078 * t919 + t1080 * t920;
t854 = mrSges(5,1) * t906 + mrSges(4,2) * t933 - mrSges(4,3) * t908 - mrSges(5,3) * t907 + pkin(4) * t876 - qJ(4) * t873 + t1057 * t1009 + t1058 * t989 + t1067 * t962 + t1069 * t995 + t1082 * t963 + t1030;
t981 = Ifges(3,6) * t1015 + (Ifges(3,4) * t1023 + Ifges(3,2) * t1026) * t1053;
t982 = Ifges(3,5) * t1015 + (Ifges(3,1) * t1023 + Ifges(3,4) * t1026) * t1053;
t843 = Ifges(3,5) * t1002 - Ifges(3,6) * t1003 + Ifges(3,3) * t1014 + mrSges(3,1) * t964 - mrSges(3,2) * t965 + t1022 * t854 + t1072 * t848 + pkin(2) * t1032 + pkin(9) * t1041 + (t1023 * t981 - t1026 * t982) * t1053;
t980 = Ifges(3,3) * t1015 + (Ifges(3,5) * t1023 + Ifges(3,6) * t1026) * t1053;
t845 = mrSges(3,2) * t983 - mrSges(3,3) * t964 + Ifges(3,1) * t1002 - Ifges(3,4) * t1003 + Ifges(3,5) * t1014 - pkin(9) * t866 - t1015 * t981 - t1022 * t848 + t980 * t1042 + t1072 * t854;
t847 = -mrSges(3,1) * t983 + mrSges(3,3) * t965 + Ifges(3,4) * t1002 - Ifges(3,2) * t1003 + Ifges(3,6) * t1014 - pkin(2) * t866 + t1015 * t982 - t1043 * t980 - t1073;
t1034 = mrSges(2,1) * t1011 - mrSges(2,2) * t1012 + Ifges(2,3) * qJDD(1) + pkin(1) * t853 + t1020 * t843 + t847 * t1049 + t845 * t1050 + t859 * t1065;
t841 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1011 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t1028 - t1023 * t847 + t1026 * t845 + (-t1019 * t852 - t1020 * t853) * pkin(8);
t840 = mrSges(2,1) * g(3) + mrSges(2,3) * t1012 + Ifges(2,5) * t1028 + Ifges(2,6) * qJDD(1) - pkin(1) * t852 - t1019 * t843 + (pkin(8) * t859 + t1023 * t845 + t1026 * t847) * t1020;
t1 = [-m(1) * g(1) + t1040; -m(1) * g(2) + t1063; (-m(1) - m(2)) * g(3) + t852; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1063 - t1024 * t840 + t1027 * t841; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1040 + t1024 * t841 + t1027 * t840; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1034; t1034; t843; t1073; t874; t1030; t892;];
tauJB  = t1;
