% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 11:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 11:12:30
% EndTime: 2019-05-08 11:14:16
% DurationCPUTime: 75.92s
% Computational Cost: add. (1284449->403), mult. (2714599->517), div. (0->0), fcn. (2229659->14), ass. (0->169)
t1007 = sin(pkin(6));
t1048 = pkin(8) * t1007;
t1008 = cos(pkin(6));
t1047 = t1008 * g(3);
t1014 = sin(qJ(1));
t1020 = cos(qJ(1));
t1021 = qJD(1) ^ 2;
t1019 = cos(qJ(2));
t1039 = t1008 * t1019;
t1013 = sin(qJ(2));
t1040 = t1008 * t1013;
t1002 = qJDD(1) * t1008 + qJDD(2);
t1003 = qJD(1) * t1008 + qJD(2);
t1012 = sin(qJ(3));
t1018 = cos(qJ(3));
t1011 = sin(qJ(4));
t1017 = cos(qJ(4));
t1010 = sin(qJ(5));
t1016 = cos(qJ(5));
t1009 = sin(qJ(6));
t1015 = cos(qJ(6));
t1001 = t1003 ^ 2;
t1043 = qJD(1) * t1019;
t998 = t1014 * g(1) - g(2) * t1020;
t988 = qJDD(1) * pkin(1) + t1021 * t1048 + t998;
t1038 = qJDD(1) * t1007;
t999 = -g(1) * t1020 - g(2) * t1014;
t989 = -pkin(1) * t1021 + pkin(8) * t1038 + t999;
t1045 = t1019 * t989 + t988 * t1040;
t1044 = qJD(1) * t1007;
t991 = (-pkin(2) * t1019 - pkin(9) * t1013) * t1044;
t944 = -t1001 * pkin(2) + t1002 * pkin(9) + (-g(3) * t1013 + t991 * t1043) * t1007 + t1045;
t992 = (qJD(2) * t1043 + qJDD(1) * t1013) * t1007;
t1042 = t1007 * t1013;
t1037 = qJD(1) * t1042;
t993 = -qJD(2) * t1037 + t1019 * t1038;
t945 = -t993 * pkin(2) - t992 * pkin(9) - t1047 + (-t988 + (pkin(2) * t1013 - pkin(9) * t1019) * t1003 * qJD(1)) * t1007;
t912 = -t1012 * t944 + t1018 * t945;
t980 = t1003 * t1018 - t1012 * t1037;
t959 = qJD(3) * t980 + t1002 * t1012 + t1018 * t992;
t981 = t1003 * t1012 + t1018 * t1037;
t985 = qJDD(3) - t993;
t1041 = t1007 * t1019;
t1036 = qJD(1) * t1041;
t997 = qJD(3) - t1036;
t897 = (t980 * t997 - t959) * pkin(10) + (t980 * t981 + t985) * pkin(3) + t912;
t913 = t1012 * t945 + t1018 * t944;
t958 = -qJD(3) * t981 + t1002 * t1018 - t1012 * t992;
t970 = pkin(3) * t997 - pkin(10) * t981;
t979 = t980 ^ 2;
t904 = -pkin(3) * t979 + pkin(10) * t958 - t970 * t997 + t913;
t894 = t1011 * t897 + t1017 * t904;
t965 = -t1011 * t981 + t1017 * t980;
t966 = t1011 * t980 + t1017 * t981;
t939 = -pkin(4) * t965 - pkin(11) * t966;
t984 = qJDD(4) + t985;
t995 = qJD(4) + t997;
t994 = t995 ^ 2;
t883 = -pkin(4) * t994 + pkin(11) * t984 + t939 * t965 + t894;
t960 = -g(3) * t1041 - t1013 * t989 + t988 * t1039;
t943 = -t1002 * pkin(2) - t1001 * pkin(9) + t991 * t1037 - t960;
t909 = -t958 * pkin(3) - t979 * pkin(10) + t981 * t970 + t943;
t924 = -t966 * qJD(4) - t1011 * t959 + t1017 * t958;
t925 = qJD(4) * t965 + t1011 * t958 + t1017 * t959;
t886 = (-t965 * t995 - t925) * pkin(11) + (t966 * t995 - t924) * pkin(4) + t909;
t878 = -t1010 * t883 + t1016 * t886;
t947 = -t1010 * t966 + t1016 * t995;
t908 = qJD(5) * t947 + t1010 * t984 + t1016 * t925;
t923 = qJDD(5) - t924;
t948 = t1010 * t995 + t1016 * t966;
t964 = qJD(5) - t965;
t876 = (t947 * t964 - t908) * pkin(12) + (t947 * t948 + t923) * pkin(5) + t878;
t879 = t1010 * t886 + t1016 * t883;
t907 = -qJD(5) * t948 - t1010 * t925 + t1016 * t984;
t933 = pkin(5) * t964 - pkin(12) * t948;
t946 = t947 ^ 2;
t877 = -pkin(5) * t946 + pkin(12) * t907 - t933 * t964 + t879;
t874 = -t1009 * t877 + t1015 * t876;
t928 = -t1009 * t948 + t1015 * t947;
t891 = qJD(6) * t928 + t1009 * t907 + t1015 * t908;
t929 = t1009 * t947 + t1015 * t948;
t905 = -mrSges(7,1) * t928 + mrSges(7,2) * t929;
t962 = qJD(6) + t964;
t910 = -mrSges(7,2) * t962 + mrSges(7,3) * t928;
t919 = qJDD(6) + t923;
t869 = m(7) * t874 + mrSges(7,1) * t919 - mrSges(7,3) * t891 - t905 * t929 + t910 * t962;
t875 = t1009 * t876 + t1015 * t877;
t890 = -qJD(6) * t929 - t1009 * t908 + t1015 * t907;
t911 = mrSges(7,1) * t962 - mrSges(7,3) * t929;
t870 = m(7) * t875 - mrSges(7,2) * t919 + mrSges(7,3) * t890 + t905 * t928 - t911 * t962;
t861 = t1009 * t870 + t1015 * t869;
t930 = -mrSges(6,1) * t947 + mrSges(6,2) * t948;
t931 = -mrSges(6,2) * t964 + mrSges(6,3) * t947;
t859 = m(6) * t878 + mrSges(6,1) * t923 - mrSges(6,3) * t908 - t930 * t948 + t931 * t964 + t861;
t1035 = -t1009 * t869 + t1015 * t870;
t932 = mrSges(6,1) * t964 - mrSges(6,3) * t948;
t860 = m(6) * t879 - mrSges(6,2) * t923 + mrSges(6,3) * t907 + t930 * t947 - t932 * t964 + t1035;
t1034 = -t1010 * t859 + t1016 * t860;
t938 = -mrSges(5,1) * t965 + mrSges(5,2) * t966;
t950 = mrSges(5,1) * t995 - mrSges(5,3) * t966;
t852 = m(5) * t894 - mrSges(5,2) * t984 + mrSges(5,3) * t924 + t938 * t965 - t950 * t995 + t1034;
t893 = -t1011 * t904 + t1017 * t897;
t882 = -pkin(4) * t984 - pkin(11) * t994 + t966 * t939 - t893;
t880 = -pkin(5) * t907 - pkin(12) * t946 + t933 * t948 + t882;
t1030 = m(7) * t880 - t890 * mrSges(7,1) + mrSges(7,2) * t891 - t928 * t910 + t911 * t929;
t1025 = -m(6) * t882 + t907 * mrSges(6,1) - mrSges(6,2) * t908 + t947 * t931 - t932 * t948 - t1030;
t949 = -mrSges(5,2) * t995 + mrSges(5,3) * t965;
t865 = m(5) * t893 + mrSges(5,1) * t984 - mrSges(5,3) * t925 - t938 * t966 + t949 * t995 + t1025;
t842 = t1011 * t852 + t1017 * t865;
t967 = -mrSges(4,1) * t980 + mrSges(4,2) * t981;
t968 = -mrSges(4,2) * t997 + mrSges(4,3) * t980;
t840 = m(4) * t912 + mrSges(4,1) * t985 - mrSges(4,3) * t959 - t967 * t981 + t968 * t997 + t842;
t1033 = -t1011 * t865 + t1017 * t852;
t969 = mrSges(4,1) * t997 - mrSges(4,3) * t981;
t841 = m(4) * t913 - mrSges(4,2) * t985 + mrSges(4,3) * t958 + t967 * t980 - t969 * t997 + t1033;
t1032 = -t1012 * t840 + t1018 * t841;
t961 = -g(3) * t1042 + t1045;
t986 = mrSges(3,1) * t1003 - mrSges(3,3) * t1037;
t990 = (-mrSges(3,1) * t1019 + mrSges(3,2) * t1013) * t1044;
t832 = m(3) * t961 - mrSges(3,2) * t1002 + mrSges(3,3) * t993 - t1003 * t986 + t990 * t1036 + t1032;
t835 = t1012 * t841 + t1018 * t840;
t974 = -t1007 * t988 - t1047;
t987 = -mrSges(3,2) * t1003 + mrSges(3,3) * t1036;
t834 = m(3) * t974 - t993 * mrSges(3,1) + t992 * mrSges(3,2) + (t1013 * t986 - t1019 * t987) * t1044 + t835;
t854 = t1010 * t860 + t1016 * t859;
t1027 = m(5) * t909 - t924 * mrSges(5,1) + mrSges(5,2) * t925 - t965 * t949 + t950 * t966 + t854;
t1024 = -m(4) * t943 + t958 * mrSges(4,1) - mrSges(4,2) * t959 + t980 * t968 - t969 * t981 - t1027;
t849 = m(3) * t960 + mrSges(3,1) * t1002 - mrSges(3,3) * t992 + t1003 * t987 - t990 * t1037 + t1024;
t822 = -t1007 * t834 + t849 * t1039 + t832 * t1040;
t819 = m(2) * t998 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1021 + t822;
t827 = -t1013 * t849 + t1019 * t832;
t825 = m(2) * t999 - mrSges(2,1) * t1021 - qJDD(1) * mrSges(2,2) + t827;
t1046 = t1014 * t825 + t1020 * t819;
t821 = t1008 * t834 + t849 * t1041 + t832 * t1042;
t1031 = -t1014 * t819 + t1020 * t825;
t898 = Ifges(7,5) * t929 + Ifges(7,6) * t928 + Ifges(7,3) * t962;
t900 = Ifges(7,1) * t929 + Ifges(7,4) * t928 + Ifges(7,5) * t962;
t862 = -mrSges(7,1) * t880 + mrSges(7,3) * t875 + Ifges(7,4) * t891 + Ifges(7,2) * t890 + Ifges(7,6) * t919 - t898 * t929 + t900 * t962;
t899 = Ifges(7,4) * t929 + Ifges(7,2) * t928 + Ifges(7,6) * t962;
t863 = mrSges(7,2) * t880 - mrSges(7,3) * t874 + Ifges(7,1) * t891 + Ifges(7,4) * t890 + Ifges(7,5) * t919 + t898 * t928 - t899 * t962;
t914 = Ifges(6,5) * t948 + Ifges(6,6) * t947 + Ifges(6,3) * t964;
t916 = Ifges(6,1) * t948 + Ifges(6,4) * t947 + Ifges(6,5) * t964;
t844 = -mrSges(6,1) * t882 + mrSges(6,3) * t879 + Ifges(6,4) * t908 + Ifges(6,2) * t907 + Ifges(6,6) * t923 - pkin(5) * t1030 + pkin(12) * t1035 + t1009 * t863 + t1015 * t862 - t948 * t914 + t964 * t916;
t915 = Ifges(6,4) * t948 + Ifges(6,2) * t947 + Ifges(6,6) * t964;
t846 = mrSges(6,2) * t882 - mrSges(6,3) * t878 + Ifges(6,1) * t908 + Ifges(6,4) * t907 + Ifges(6,5) * t923 - pkin(12) * t861 - t1009 * t862 + t1015 * t863 + t914 * t947 - t915 * t964;
t934 = Ifges(5,5) * t966 + Ifges(5,6) * t965 + Ifges(5,3) * t995;
t935 = Ifges(5,4) * t966 + Ifges(5,2) * t965 + Ifges(5,6) * t995;
t828 = mrSges(5,2) * t909 - mrSges(5,3) * t893 + Ifges(5,1) * t925 + Ifges(5,4) * t924 + Ifges(5,5) * t984 - pkin(11) * t854 - t1010 * t844 + t1016 * t846 + t934 * t965 - t935 * t995;
t1028 = -mrSges(7,1) * t874 + mrSges(7,2) * t875 - Ifges(7,5) * t891 - Ifges(7,6) * t890 - Ifges(7,3) * t919 - t929 * t899 + t928 * t900;
t1023 = mrSges(6,1) * t878 - mrSges(6,2) * t879 + Ifges(6,5) * t908 + Ifges(6,6) * t907 + Ifges(6,3) * t923 + pkin(5) * t861 + t948 * t915 - t947 * t916 - t1028;
t936 = Ifges(5,1) * t966 + Ifges(5,4) * t965 + Ifges(5,5) * t995;
t836 = -mrSges(5,1) * t909 + mrSges(5,3) * t894 + Ifges(5,4) * t925 + Ifges(5,2) * t924 + Ifges(5,6) * t984 - pkin(4) * t854 - t966 * t934 + t995 * t936 - t1023;
t952 = Ifges(4,5) * t981 + Ifges(4,6) * t980 + Ifges(4,3) * t997;
t954 = Ifges(4,1) * t981 + Ifges(4,4) * t980 + Ifges(4,5) * t997;
t814 = -mrSges(4,1) * t943 + mrSges(4,3) * t913 + Ifges(4,4) * t959 + Ifges(4,2) * t958 + Ifges(4,6) * t985 - pkin(3) * t1027 + pkin(10) * t1033 + t1011 * t828 + t1017 * t836 - t981 * t952 + t997 * t954;
t953 = Ifges(4,4) * t981 + Ifges(4,2) * t980 + Ifges(4,6) * t997;
t817 = mrSges(4,2) * t943 - mrSges(4,3) * t912 + Ifges(4,1) * t959 + Ifges(4,4) * t958 + Ifges(4,5) * t985 - pkin(10) * t842 - t1011 * t836 + t1017 * t828 + t952 * t980 - t953 * t997;
t972 = Ifges(3,6) * t1003 + (Ifges(3,4) * t1013 + Ifges(3,2) * t1019) * t1044;
t973 = Ifges(3,5) * t1003 + (Ifges(3,1) * t1013 + Ifges(3,4) * t1019) * t1044;
t811 = Ifges(3,5) * t992 + Ifges(3,6) * t993 + Ifges(3,3) * t1002 + mrSges(3,1) * t960 - mrSges(3,2) * t961 + t1012 * t817 + t1018 * t814 + pkin(2) * t1024 + pkin(9) * t1032 + (t1013 * t972 - t1019 * t973) * t1044;
t971 = Ifges(3,3) * t1003 + (Ifges(3,5) * t1013 + Ifges(3,6) * t1019) * t1044;
t813 = mrSges(3,2) * t974 - mrSges(3,3) * t960 + Ifges(3,1) * t992 + Ifges(3,4) * t993 + Ifges(3,5) * t1002 - pkin(9) * t835 - t1003 * t972 - t1012 * t814 + t1018 * t817 + t971 * t1036;
t1026 = -mrSges(5,1) * t893 + mrSges(5,2) * t894 - Ifges(5,5) * t925 - Ifges(5,6) * t924 - Ifges(5,3) * t984 - pkin(4) * t1025 - pkin(11) * t1034 - t1010 * t846 - t1016 * t844 - t966 * t935 + t965 * t936;
t1022 = mrSges(4,1) * t912 - mrSges(4,2) * t913 + Ifges(4,5) * t959 + Ifges(4,6) * t958 + Ifges(4,3) * t985 + pkin(3) * t842 + t981 * t953 - t980 * t954 - t1026;
t816 = -mrSges(3,1) * t974 + mrSges(3,3) * t961 + Ifges(3,4) * t992 + Ifges(3,2) * t993 + Ifges(3,6) * t1002 - pkin(2) * t835 + t1003 * t973 - t971 * t1037 - t1022;
t1029 = mrSges(2,1) * t998 - mrSges(2,2) * t999 + Ifges(2,3) * qJDD(1) + pkin(1) * t822 + t1008 * t811 + t816 * t1041 + t813 * t1042 + t827 * t1048;
t809 = -mrSges(2,2) * g(3) - mrSges(2,3) * t998 + Ifges(2,5) * qJDD(1) - t1021 * Ifges(2,6) - t1013 * t816 + t1019 * t813 + (-t1007 * t821 - t1008 * t822) * pkin(8);
t808 = mrSges(2,1) * g(3) + mrSges(2,3) * t999 + t1021 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t821 - t1007 * t811 + (pkin(8) * t827 + t1013 * t813 + t1019 * t816) * t1008;
t1 = [-m(1) * g(1) + t1031; -m(1) * g(2) + t1046; (-m(1) - m(2)) * g(3) + t821; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1046 - t1014 * t808 + t1020 * t809; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1031 + t1014 * t809 + t1020 * t808; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1029; t1029; t811; t1022; -t1026; t1023; -t1028;];
tauJB  = t1;
