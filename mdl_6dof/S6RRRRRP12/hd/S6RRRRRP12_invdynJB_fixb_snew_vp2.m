% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-05-08 07:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRP12_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP12_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 07:23:52
% EndTime: 2019-05-08 07:25:44
% DurationCPUTime: 81.25s
% Computational Cost: add. (1355950->395), mult. (3342734->514), div. (0->0), fcn. (2808994->14), ass. (0->173)
t1080 = Ifges(6,1) + Ifges(7,1);
t1074 = Ifges(6,4) - Ifges(7,5);
t1073 = -Ifges(6,5) - Ifges(7,4);
t1079 = Ifges(6,2) + Ifges(7,3);
t1072 = Ifges(6,6) - Ifges(7,6);
t1078 = -Ifges(6,3) - Ifges(7,2);
t1019 = sin(pkin(7));
t1021 = cos(pkin(7));
t1025 = sin(qJ(3));
t1029 = cos(qJ(3));
t1022 = cos(pkin(6));
t1016 = qJD(1) * t1022 + qJD(2);
t1020 = sin(pkin(6));
t1030 = cos(qJ(2));
t1054 = t1020 * t1030;
t1047 = qJD(1) * t1054;
t1000 = (t1016 * t1019 + t1021 * t1047) * pkin(10);
t1026 = sin(qJ(2));
t1060 = qJD(1) * t1020;
t1069 = pkin(10) * t1019;
t1004 = (-pkin(2) * t1030 - t1026 * t1069) * t1060;
t1058 = qJD(1) * t1030;
t1010 = (qJD(2) * t1058 + qJDD(1) * t1026) * t1020;
t1015 = qJDD(1) * t1022 + qJDD(2);
t1027 = sin(qJ(1));
t1031 = cos(qJ(1));
t1013 = g(1) * t1027 - g(2) * t1031;
t1032 = qJD(1) ^ 2;
t1070 = pkin(9) * t1020;
t1007 = qJDD(1) * pkin(1) + t1032 * t1070 + t1013;
t1014 = -g(1) * t1031 - g(2) * t1027;
t1008 = -pkin(1) * t1032 + qJDD(1) * t1070 + t1014;
t1050 = t1022 * t1030;
t1043 = t1007 * t1050 - t1026 * t1008;
t1059 = qJD(1) * t1026;
t1068 = pkin(10) * t1021;
t958 = -t1010 * t1068 + t1015 * pkin(2) + t1016 * t1000 + (-g(3) * t1030 - t1004 * t1059) * t1020 + t1043;
t1055 = t1020 * t1026;
t1048 = qJD(1) * t1055;
t1003 = pkin(2) * t1016 - t1048 * t1068;
t1011 = (-qJD(2) * t1059 + qJDD(1) * t1030) * t1020;
t1039 = t1011 * t1021 + t1015 * t1019;
t1051 = t1022 * t1026;
t1061 = t1007 * t1051 + t1008 * t1030;
t959 = -t1016 * t1003 + (-g(3) * t1026 + t1004 * t1058) * t1020 + t1039 * pkin(10) + t1061;
t1067 = t1022 * g(3);
t964 = -t1010 * t1069 - t1011 * pkin(2) - t1067 + (-t1007 + (-t1000 * t1030 + t1003 * t1026) * qJD(1)) * t1020;
t924 = -t1025 * t959 + (t1019 * t964 + t1021 * t958) * t1029;
t1052 = t1021 * t1030;
t1057 = t1019 * t1025;
t991 = t1016 * t1057 + (t1025 * t1052 + t1026 * t1029) * t1060;
t975 = -t991 * qJD(3) - t1025 * t1010 + t1029 * t1039;
t1056 = t1019 * t1029;
t990 = (-t1025 * t1026 + t1029 * t1052) * t1060 + t1016 * t1056;
t1023 = sin(qJ(5));
t1076 = cos(qJ(5));
t1024 = sin(qJ(4));
t1028 = cos(qJ(4));
t1053 = t1021 * t1025;
t925 = t1029 * t959 + t1053 * t958 + t1057 * t964;
t978 = -pkin(3) * t990 - pkin(11) * t991;
t992 = -t1011 * t1019 + t1015 * t1021 + qJDD(3);
t1001 = t1016 * t1021 - t1019 * t1047 + qJD(3);
t999 = t1001 ^ 2;
t917 = -pkin(3) * t999 + pkin(11) * t992 + t978 * t990 + t925;
t934 = -t1019 * t958 + t1021 * t964;
t976 = t990 * qJD(3) + t1029 * t1010 + t1025 * t1039;
t919 = (-t1001 * t990 - t976) * pkin(11) + (t1001 * t991 - t975) * pkin(3) + t934;
t913 = t1024 * t919 + t1028 * t917;
t981 = t1001 * t1028 - t1024 * t991;
t982 = t1001 * t1024 + t1028 * t991;
t961 = -pkin(4) * t981 - pkin(12) * t982;
t974 = qJDD(4) - t975;
t989 = qJD(4) - t990;
t988 = t989 ^ 2;
t909 = -pkin(4) * t988 + pkin(12) * t974 + t961 * t981 + t913;
t916 = -t992 * pkin(3) - t999 * pkin(11) + t978 * t991 - t924;
t944 = -qJD(4) * t982 - t1024 * t976 + t1028 * t992;
t945 = qJD(4) * t981 + t1024 * t992 + t1028 * t976;
t911 = (-t981 * t989 - t945) * pkin(12) + (t982 * t989 - t944) * pkin(4) + t916;
t906 = t1023 * t911 + t1076 * t909;
t966 = t1023 * t982 - t1076 * t989;
t967 = t1023 * t989 + t1076 * t982;
t937 = pkin(5) * t966 - qJ(6) * t967;
t943 = qJDD(5) - t944;
t980 = qJD(5) - t981;
t979 = t980 ^ 2;
t902 = -pkin(5) * t979 + qJ(6) * t943 + 0.2e1 * qJD(6) * t980 - t937 * t966 + t906;
t950 = -mrSges(7,1) * t980 + mrSges(7,2) * t967;
t1049 = m(7) * t902 + mrSges(7,3) * t943 + t950 * t980;
t1063 = t1073 * t980 + t1074 * t966 - t1080 * t967;
t1064 = -t1072 * t980 - t1074 * t967 + t1079 * t966;
t905 = -t1023 * t909 + t1076 * t911;
t903 = -pkin(5) * t943 - qJ(6) * t979 + t937 * t967 + qJDD(6) - t905;
t947 = -mrSges(7,2) * t966 + mrSges(7,3) * t980;
t1042 = -m(7) * t903 + mrSges(7,1) * t943 + t947 * t980;
t923 = -qJD(5) * t966 + t1023 * t974 + t1076 * t945;
t938 = mrSges(7,1) * t966 - mrSges(7,3) * t967;
t899 = t923 * mrSges(7,2) + t967 * t938 - t1042;
t922 = qJD(5) * t967 + t1023 * t945 - t1076 * t974;
t1077 = -t1063 * t966 - t1064 * t967 - t1078 * t943 - t1072 * t922 - t1073 * t923 + mrSges(6,1) * t905 - mrSges(7,1) * t903 - mrSges(6,2) * t906 + mrSges(7,3) * t902 - pkin(5) * t899 + qJ(6) * (-t922 * mrSges(7,2) - t966 * t938 + t1049);
t1075 = -mrSges(6,3) - mrSges(7,2);
t1006 = -mrSges(3,2) * t1016 + mrSges(3,3) * t1047;
t1009 = (-mrSges(3,1) * t1030 + mrSges(3,2) * t1026) * t1060;
t1062 = -mrSges(6,1) * t966 - mrSges(6,2) * t967 - t938;
t949 = mrSges(6,1) * t980 - mrSges(6,3) * t967;
t895 = m(6) * t906 - t943 * mrSges(6,2) + t1062 * t966 + t1075 * t922 - t980 * t949 + t1049;
t948 = -mrSges(6,2) * t980 - mrSges(6,3) * t966;
t896 = m(6) * t905 + t943 * mrSges(6,1) + t1062 * t967 + t1075 * t923 + t980 * t948 + t1042;
t891 = -t1023 * t896 + t1076 * t895;
t960 = -mrSges(5,1) * t981 + mrSges(5,2) * t982;
t969 = mrSges(5,1) * t989 - mrSges(5,3) * t982;
t887 = m(5) * t913 - mrSges(5,2) * t974 + mrSges(5,3) * t944 + t960 * t981 - t969 * t989 + t891;
t912 = -t1024 * t917 + t1028 * t919;
t908 = -t974 * pkin(4) - t988 * pkin(12) + t961 * t982 - t912;
t904 = -0.2e1 * qJD(6) * t967 + (t966 * t980 - t923) * qJ(6) + (t967 * t980 + t922) * pkin(5) + t908;
t900 = m(7) * t904 + mrSges(7,1) * t922 - mrSges(7,3) * t923 + t947 * t966 - t950 * t967;
t897 = -m(6) * t908 - mrSges(6,1) * t922 - mrSges(6,2) * t923 - t948 * t966 - t949 * t967 - t900;
t968 = -mrSges(5,2) * t989 + mrSges(5,3) * t981;
t893 = m(5) * t912 + mrSges(5,1) * t974 - mrSges(5,3) * t945 - t960 * t982 + t968 * t989 + t897;
t1045 = -t1024 * t893 + t1028 * t887;
t977 = -mrSges(4,1) * t990 + mrSges(4,2) * t991;
t984 = mrSges(4,1) * t1001 - mrSges(4,3) * t991;
t878 = m(4) * t925 - mrSges(4,2) * t992 + mrSges(4,3) * t975 - t1001 * t984 + t977 * t990 + t1045;
t881 = t1024 * t887 + t1028 * t893;
t983 = -mrSges(4,2) * t1001 + mrSges(4,3) * t990;
t880 = m(4) * t934 - mrSges(4,1) * t975 + mrSges(4,2) * t976 - t983 * t990 + t984 * t991 + t881;
t890 = t1023 * t895 + t1076 * t896;
t1034 = -m(5) * t916 + mrSges(5,1) * t944 - mrSges(5,2) * t945 + t968 * t981 - t969 * t982 - t890;
t884 = m(4) * t924 + mrSges(4,1) * t992 - mrSges(4,3) * t976 + t1001 * t983 - t977 * t991 + t1034;
t867 = t1021 * t1029 * t884 - t1019 * t880 + t1053 * t878;
t985 = -g(3) * t1054 + t1043;
t863 = m(3) * t985 + mrSges(3,1) * t1015 - mrSges(3,3) * t1010 + t1006 * t1016 - t1009 * t1048 + t867;
t1005 = mrSges(3,1) * t1016 - mrSges(3,3) * t1048;
t866 = t1021 * t880 + t1056 * t884 + t1057 * t878;
t996 = -t1020 * t1007 - t1067;
t865 = m(3) * t996 - t1011 * mrSges(3,1) + t1010 * mrSges(3,2) + (t1005 * t1026 - t1006 * t1030) * t1060 + t866;
t873 = -t1025 * t884 + t1029 * t878;
t986 = -g(3) * t1055 + t1061;
t872 = m(3) * t986 - mrSges(3,2) * t1015 + mrSges(3,3) * t1011 - t1005 * t1016 + t1009 * t1047 + t873;
t852 = -t1020 * t865 + t1050 * t863 + t1051 * t872;
t849 = m(2) * t1013 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1032 + t852;
t859 = -t1026 * t863 + t1030 * t872;
t857 = m(2) * t1014 - mrSges(2,1) * t1032 - qJDD(1) * mrSges(2,2) + t859;
t1066 = t1027 * t857 + t1031 * t849;
t1065 = t1072 * t966 + t1073 * t967 + t1078 * t980;
t851 = t1022 * t865 + t1054 * t863 + t1055 * t872;
t1044 = -t1027 * t849 + t1031 * t857;
t888 = -mrSges(6,1) * t908 - mrSges(7,1) * t904 + mrSges(7,2) * t902 + mrSges(6,3) * t906 - pkin(5) * t900 - t1063 * t980 + t1065 * t967 + t1072 * t943 + t1074 * t923 - t1079 * t922;
t889 = mrSges(6,2) * t908 + mrSges(7,2) * t903 - mrSges(6,3) * t905 - mrSges(7,3) * t904 - qJ(6) * t900 + t1064 * t980 + t1065 * t966 - t1073 * t943 - t1074 * t922 + t1080 * t923;
t951 = Ifges(5,5) * t982 + Ifges(5,6) * t981 + Ifges(5,3) * t989;
t952 = Ifges(5,4) * t982 + Ifges(5,2) * t981 + Ifges(5,6) * t989;
t868 = mrSges(5,2) * t916 - mrSges(5,3) * t912 + Ifges(5,1) * t945 + Ifges(5,4) * t944 + Ifges(5,5) * t974 - pkin(12) * t890 - t1023 * t888 + t1076 * t889 + t951 * t981 - t952 * t989;
t953 = Ifges(5,1) * t982 + Ifges(5,4) * t981 + Ifges(5,5) * t989;
t874 = -mrSges(5,1) * t916 + mrSges(5,3) * t913 + Ifges(5,4) * t945 + Ifges(5,2) * t944 + Ifges(5,6) * t974 - pkin(4) * t890 - t982 * t951 + t989 * t953 - t1077;
t970 = Ifges(4,5) * t991 + Ifges(4,6) * t990 + Ifges(4,3) * t1001;
t971 = Ifges(4,4) * t991 + Ifges(4,2) * t990 + Ifges(4,6) * t1001;
t854 = mrSges(4,2) * t934 - mrSges(4,3) * t924 + Ifges(4,1) * t976 + Ifges(4,4) * t975 + Ifges(4,5) * t992 - pkin(11) * t881 - t1001 * t971 - t1024 * t874 + t1028 * t868 + t970 * t990;
t1033 = mrSges(5,1) * t912 - mrSges(5,2) * t913 + Ifges(5,5) * t945 + Ifges(5,6) * t944 + Ifges(5,3) * t974 + pkin(4) * t897 + pkin(12) * t891 + t1023 * t889 + t1076 * t888 + t952 * t982 - t953 * t981;
t972 = Ifges(4,1) * t991 + Ifges(4,4) * t990 + Ifges(4,5) * t1001;
t860 = -mrSges(4,1) * t934 + mrSges(4,3) * t925 + Ifges(4,4) * t976 + Ifges(4,2) * t975 + Ifges(4,6) * t992 - pkin(3) * t881 + t1001 * t972 - t970 * t991 - t1033;
t1037 = pkin(10) * t873 + t1025 * t854 + t1029 * t860;
t853 = mrSges(4,1) * t924 - mrSges(4,2) * t925 + Ifges(4,5) * t976 + Ifges(4,6) * t975 + Ifges(4,3) * t992 + pkin(3) * t1034 + pkin(11) * t1045 + t1024 * t868 + t1028 * t874 + t991 * t971 - t990 * t972;
t994 = Ifges(3,6) * t1016 + (Ifges(3,4) * t1026 + Ifges(3,2) * t1030) * t1060;
t995 = Ifges(3,5) * t1016 + (Ifges(3,1) * t1026 + Ifges(3,4) * t1030) * t1060;
t843 = mrSges(3,1) * t985 - mrSges(3,2) * t986 + Ifges(3,5) * t1010 + Ifges(3,6) * t1011 + Ifges(3,3) * t1015 + pkin(2) * t867 + t1021 * t853 + (t1026 * t994 - t1030 * t995) * t1060 + t1037 * t1019;
t993 = Ifges(3,3) * t1016 + (Ifges(3,5) * t1026 + Ifges(3,6) * t1030) * t1060;
t845 = -mrSges(3,1) * t996 + mrSges(3,3) * t986 + Ifges(3,4) * t1010 + Ifges(3,2) * t1011 + Ifges(3,6) * t1015 - pkin(2) * t866 + t1016 * t995 - t1019 * t853 + t1021 * t1037 - t1048 * t993;
t847 = t993 * t1047 + mrSges(3,2) * t996 - mrSges(3,3) * t985 + Ifges(3,1) * t1010 + Ifges(3,4) * t1011 + Ifges(3,5) * t1015 - t1016 * t994 - t1025 * t860 + t1029 * t854 + (-t1019 * t866 - t1021 * t867) * pkin(10);
t1036 = mrSges(2,1) * t1013 - mrSges(2,2) * t1014 + Ifges(2,3) * qJDD(1) + pkin(1) * t852 + t1022 * t843 + t1054 * t845 + t1055 * t847 + t1070 * t859;
t841 = -mrSges(2,2) * g(3) - mrSges(2,3) * t1013 + Ifges(2,5) * qJDD(1) - t1032 * Ifges(2,6) - t1026 * t845 + t1030 * t847 + (-t1020 * t851 - t1022 * t852) * pkin(9);
t840 = mrSges(2,1) * g(3) + mrSges(2,3) * t1014 + t1032 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t851 - t1020 * t843 + (pkin(9) * t859 + t1026 * t847 + t1030 * t845) * t1022;
t1 = [-m(1) * g(1) + t1044; -m(1) * g(2) + t1066; (-m(1) - m(2)) * g(3) + t851; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t1066 - t1027 * t840 + t1031 * t841; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t1044 + t1027 * t841 + t1031 * t840; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1036; t1036; t843; t853; t1033; t1077; t899;];
tauJB  = t1;
