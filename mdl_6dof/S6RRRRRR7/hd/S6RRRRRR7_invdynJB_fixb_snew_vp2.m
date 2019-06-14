% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRRR7
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
% Datum: 2019-05-08 12:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 12:22:35
% EndTime: 2019-05-08 12:24:30
% DurationCPUTime: 85.73s
% Computational Cost: add. (1480549->403), mult. (3154708->518), div. (0->0), fcn. (2569485->14), ass. (0->169)
t1009 = cos(qJ(2));
t1003 = sin(qJ(2));
t1028 = qJD(1) * t1003;
t997 = sin(pkin(6));
t984 = (qJD(2) * t1028 - qJDD(1) * t1009) * t997;
t1037 = pkin(8) * t997;
t998 = cos(pkin(6));
t1036 = t998 * g(3);
t1004 = sin(qJ(1));
t1010 = cos(qJ(1));
t1011 = qJD(1) ^ 2;
t1029 = t1009 * t998;
t1031 = t1003 * t998;
t1002 = sin(qJ(3));
t1008 = cos(qJ(3));
t1001 = sin(qJ(4));
t1007 = cos(qJ(4));
t1000 = sin(qJ(5));
t1006 = cos(qJ(5));
t1005 = cos(qJ(6));
t1027 = qJD(1) * t1009;
t989 = t1004 * g(1) - t1010 * g(2);
t979 = qJDD(1) * pkin(1) + t1011 * t1037 + t989;
t990 = -t1010 * g(1) - t1004 * g(2);
t980 = -t1011 * pkin(1) + qJDD(1) * t1037 + t990;
t1034 = t1009 * t980 + t979 * t1031;
t1033 = qJD(1) * t997;
t982 = (-pkin(2) * t1009 - pkin(9) * t1003) * t1033;
t993 = t998 * qJD(1) + qJD(2);
t991 = t993 ^ 2;
t992 = t998 * qJDD(1) + qJDD(2);
t930 = -t991 * pkin(2) + t992 * pkin(9) + (-g(3) * t1003 + t1027 * t982) * t997 + t1034;
t983 = (qJD(2) * t1027 + qJDD(1) * t1003) * t997;
t931 = t984 * pkin(2) - t983 * pkin(9) - t1036 + (-t979 + (pkin(2) * t1003 - pkin(9) * t1009) * t993 * qJD(1)) * t997;
t907 = t1002 * t931 + t1008 * t930;
t1026 = t997 * t1028;
t971 = -t1002 * t1026 + t1008 * t993;
t972 = t1002 * t993 + t1008 * t1026;
t955 = -t971 * pkin(3) - t972 * pkin(10);
t976 = qJDD(3) + t984;
t1025 = t997 * t1027;
t988 = qJD(3) - t1025;
t986 = t988 ^ 2;
t899 = -t986 * pkin(3) + t976 * pkin(10) + t971 * t955 + t907;
t1030 = t1009 * t997;
t952 = -g(3) * t1030 - t1003 * t980 + t1029 * t979;
t929 = -t992 * pkin(2) - t991 * pkin(9) + t982 * t1026 - t952;
t950 = -t972 * qJD(3) - t1002 * t983 + t1008 * t992;
t951 = t971 * qJD(3) + t1002 * t992 + t1008 * t983;
t903 = (-t971 * t988 - t951) * pkin(10) + (t972 * t988 - t950) * pkin(3) + t929;
t882 = -t1001 * t899 + t1007 * t903;
t957 = -t1001 * t972 + t1007 * t988;
t917 = t957 * qJD(4) + t1001 * t976 + t1007 * t951;
t948 = qJDD(4) - t950;
t958 = t1001 * t988 + t1007 * t972;
t970 = qJD(4) - t971;
t873 = (t957 * t970 - t917) * pkin(11) + (t957 * t958 + t948) * pkin(4) + t882;
t883 = t1001 * t903 + t1007 * t899;
t916 = -t958 * qJD(4) - t1001 * t951 + t1007 * t976;
t939 = t970 * pkin(4) - t958 * pkin(11);
t956 = t957 ^ 2;
t881 = -t956 * pkin(4) + t916 * pkin(11) - t970 * t939 + t883;
t867 = -t1000 * t881 + t1006 * t873;
t933 = -t1000 * t958 + t1006 * t957;
t896 = t933 * qJD(5) + t1000 * t916 + t1006 * t917;
t934 = t1000 * t957 + t1006 * t958;
t943 = qJDD(5) + t948;
t968 = qJD(5) + t970;
t864 = (t933 * t968 - t896) * pkin(12) + (t933 * t934 + t943) * pkin(5) + t867;
t868 = t1000 * t873 + t1006 * t881;
t895 = -t934 * qJD(5) - t1000 * t917 + t1006 * t916;
t920 = t968 * pkin(5) - t934 * pkin(12);
t932 = t933 ^ 2;
t865 = -t932 * pkin(5) + t895 * pkin(12) - t968 * t920 + t868;
t999 = sin(qJ(6));
t862 = t1005 * t864 - t999 * t865;
t912 = t1005 * t933 - t999 * t934;
t879 = t912 * qJD(6) + t1005 * t896 + t999 * t895;
t913 = t1005 * t934 + t999 * t933;
t891 = -t912 * mrSges(7,1) + t913 * mrSges(7,2);
t965 = qJD(6) + t968;
t904 = -t965 * mrSges(7,2) + t912 * mrSges(7,3);
t942 = qJDD(6) + t943;
t858 = m(7) * t862 + t942 * mrSges(7,1) - t879 * mrSges(7,3) - t913 * t891 + t965 * t904;
t863 = t1005 * t865 + t999 * t864;
t878 = -t913 * qJD(6) + t1005 * t895 - t999 * t896;
t905 = t965 * mrSges(7,1) - t913 * mrSges(7,3);
t859 = m(7) * t863 - t942 * mrSges(7,2) + t878 * mrSges(7,3) + t912 * t891 - t965 * t905;
t850 = t1005 * t858 + t999 * t859;
t914 = -t933 * mrSges(6,1) + t934 * mrSges(6,2);
t918 = -t968 * mrSges(6,2) + t933 * mrSges(6,3);
t847 = m(6) * t867 + t943 * mrSges(6,1) - t896 * mrSges(6,3) - t934 * t914 + t968 * t918 + t850;
t1024 = t1005 * t859 - t999 * t858;
t919 = t968 * mrSges(6,1) - t934 * mrSges(6,3);
t848 = m(6) * t868 - t943 * mrSges(6,2) + t895 * mrSges(6,3) + t933 * t914 - t968 * t919 + t1024;
t843 = t1000 * t848 + t1006 * t847;
t935 = -t957 * mrSges(5,1) + t958 * mrSges(5,2);
t937 = -t970 * mrSges(5,2) + t957 * mrSges(5,3);
t841 = m(5) * t882 + t948 * mrSges(5,1) - t917 * mrSges(5,3) - t958 * t935 + t970 * t937 + t843;
t1023 = -t1000 * t847 + t1006 * t848;
t938 = t970 * mrSges(5,1) - t958 * mrSges(5,3);
t842 = m(5) * t883 - t948 * mrSges(5,2) + t916 * mrSges(5,3) + t957 * t935 - t970 * t938 + t1023;
t837 = -t1001 * t841 + t1007 * t842;
t954 = -t971 * mrSges(4,1) + t972 * mrSges(4,2);
t960 = t988 * mrSges(4,1) - t972 * mrSges(4,3);
t835 = m(4) * t907 - t976 * mrSges(4,2) + t950 * mrSges(4,3) + t971 * t954 - t988 * t960 + t837;
t906 = -t1002 * t930 + t1008 * t931;
t898 = -t976 * pkin(3) - t986 * pkin(10) + t972 * t955 - t906;
t884 = -t916 * pkin(4) - t956 * pkin(11) + t958 * t939 + t898;
t870 = -t895 * pkin(5) - t932 * pkin(12) + t934 * t920 + t884;
t1020 = m(7) * t870 - t878 * mrSges(7,1) + t879 * mrSges(7,2) - t912 * t904 + t913 * t905;
t1016 = m(6) * t884 - t895 * mrSges(6,1) + t896 * mrSges(6,2) - t933 * t918 + t934 * t919 + t1020;
t860 = -m(5) * t898 + t916 * mrSges(5,1) - t917 * mrSges(5,2) + t957 * t937 - t958 * t938 - t1016;
t959 = -t988 * mrSges(4,2) + t971 * mrSges(4,3);
t854 = m(4) * t906 + t976 * mrSges(4,1) - t951 * mrSges(4,3) - t972 * t954 + t988 * t959 + t860;
t1022 = -t1002 * t854 + t1008 * t835;
t1032 = t1003 * t997;
t953 = -g(3) * t1032 + t1034;
t977 = t993 * mrSges(3,1) - mrSges(3,3) * t1026;
t981 = (-mrSges(3,1) * t1009 + mrSges(3,2) * t1003) * t1033;
t826 = m(3) * t953 - t992 * mrSges(3,2) - t984 * mrSges(3,3) + t1025 * t981 - t993 * t977 + t1022;
t829 = t1002 * t835 + t1008 * t854;
t964 = -t997 * t979 - t1036;
t978 = -t993 * mrSges(3,2) + mrSges(3,3) * t1025;
t828 = m(3) * t964 + t984 * mrSges(3,1) + t983 * mrSges(3,2) + (t1003 * t977 - t1009 * t978) * t1033 + t829;
t836 = t1001 * t842 + t1007 * t841;
t1015 = -m(4) * t929 + t950 * mrSges(4,1) - t951 * mrSges(4,2) + t971 * t959 - t972 * t960 - t836;
t832 = m(3) * t952 + t992 * mrSges(3,1) - t983 * mrSges(3,3) - t1026 * t981 + t993 * t978 + t1015;
t814 = t832 * t1029 + t826 * t1031 - t997 * t828;
t811 = m(2) * t989 + qJDD(1) * mrSges(2,1) - t1011 * mrSges(2,2) + t814;
t819 = -t1003 * t832 + t1009 * t826;
t817 = m(2) * t990 - t1011 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t819;
t1035 = t1004 * t817 + t1010 * t811;
t813 = t832 * t1030 + t826 * t1032 + t998 * t828;
t1021 = -t1004 * t811 + t1010 * t817;
t886 = Ifges(7,5) * t913 + Ifges(7,6) * t912 + Ifges(7,3) * t965;
t888 = Ifges(7,1) * t913 + Ifges(7,4) * t912 + Ifges(7,5) * t965;
t851 = -mrSges(7,1) * t870 + mrSges(7,3) * t863 + Ifges(7,4) * t879 + Ifges(7,2) * t878 + Ifges(7,6) * t942 - t913 * t886 + t965 * t888;
t887 = Ifges(7,4) * t913 + Ifges(7,2) * t912 + Ifges(7,6) * t965;
t852 = mrSges(7,2) * t870 - mrSges(7,3) * t862 + Ifges(7,1) * t879 + Ifges(7,4) * t878 + Ifges(7,5) * t942 + t912 * t886 - t965 * t887;
t908 = Ifges(6,5) * t934 + Ifges(6,6) * t933 + Ifges(6,3) * t968;
t910 = Ifges(6,1) * t934 + Ifges(6,4) * t933 + Ifges(6,5) * t968;
t838 = -mrSges(6,1) * t884 + mrSges(6,3) * t868 + Ifges(6,4) * t896 + Ifges(6,2) * t895 + Ifges(6,6) * t943 - pkin(5) * t1020 + pkin(12) * t1024 + t1005 * t851 + t999 * t852 - t934 * t908 + t968 * t910;
t909 = Ifges(6,4) * t934 + Ifges(6,2) * t933 + Ifges(6,6) * t968;
t839 = mrSges(6,2) * t884 - mrSges(6,3) * t867 + Ifges(6,1) * t896 + Ifges(6,4) * t895 + Ifges(6,5) * t943 - pkin(12) * t850 + t1005 * t852 - t999 * t851 + t933 * t908 - t968 * t909;
t921 = Ifges(5,5) * t958 + Ifges(5,6) * t957 + Ifges(5,3) * t970;
t923 = Ifges(5,1) * t958 + Ifges(5,4) * t957 + Ifges(5,5) * t970;
t821 = -mrSges(5,1) * t898 + mrSges(5,3) * t883 + Ifges(5,4) * t917 + Ifges(5,2) * t916 + Ifges(5,6) * t948 - pkin(4) * t1016 + pkin(11) * t1023 + t1000 * t839 + t1006 * t838 - t958 * t921 + t970 * t923;
t922 = Ifges(5,4) * t958 + Ifges(5,2) * t957 + Ifges(5,6) * t970;
t822 = mrSges(5,2) * t898 - mrSges(5,3) * t882 + Ifges(5,1) * t917 + Ifges(5,4) * t916 + Ifges(5,5) * t948 - pkin(11) * t843 - t1000 * t838 + t1006 * t839 + t957 * t921 - t970 * t922;
t944 = Ifges(4,5) * t972 + Ifges(4,6) * t971 + Ifges(4,3) * t988;
t945 = Ifges(4,4) * t972 + Ifges(4,2) * t971 + Ifges(4,6) * t988;
t809 = mrSges(4,2) * t929 - mrSges(4,3) * t906 + Ifges(4,1) * t951 + Ifges(4,4) * t950 + Ifges(4,5) * t976 - pkin(10) * t836 - t1001 * t821 + t1007 * t822 + t971 * t944 - t988 * t945;
t1017 = -mrSges(7,1) * t862 + mrSges(7,2) * t863 - Ifges(7,5) * t879 - Ifges(7,6) * t878 - Ifges(7,3) * t942 - t913 * t887 + t912 * t888;
t1014 = -mrSges(6,1) * t867 + mrSges(6,2) * t868 - Ifges(6,5) * t896 - Ifges(6,6) * t895 - Ifges(6,3) * t943 - pkin(5) * t850 - t934 * t909 + t933 * t910 + t1017;
t1012 = mrSges(5,1) * t882 - mrSges(5,2) * t883 + Ifges(5,5) * t917 + Ifges(5,6) * t916 + Ifges(5,3) * t948 + pkin(4) * t843 + t958 * t922 - t957 * t923 - t1014;
t946 = Ifges(4,1) * t972 + Ifges(4,4) * t971 + Ifges(4,5) * t988;
t820 = -mrSges(4,1) * t929 + mrSges(4,3) * t907 + Ifges(4,4) * t951 + Ifges(4,2) * t950 + Ifges(4,6) * t976 - pkin(3) * t836 - t972 * t944 + t988 * t946 - t1012;
t962 = Ifges(3,6) * t993 + (Ifges(3,4) * t1003 + Ifges(3,2) * t1009) * t1033;
t963 = Ifges(3,5) * t993 + (Ifges(3,1) * t1003 + Ifges(3,4) * t1009) * t1033;
t804 = Ifges(3,5) * t983 - Ifges(3,6) * t984 + Ifges(3,3) * t992 + mrSges(3,1) * t952 - mrSges(3,2) * t953 + t1002 * t809 + t1008 * t820 + pkin(2) * t1015 + pkin(9) * t1022 + (t1003 * t962 - t1009 * t963) * t1033;
t961 = Ifges(3,3) * t993 + (Ifges(3,5) * t1003 + Ifges(3,6) * t1009) * t1033;
t806 = mrSges(3,2) * t964 - mrSges(3,3) * t952 + Ifges(3,1) * t983 - Ifges(3,4) * t984 + Ifges(3,5) * t992 - pkin(9) * t829 - t1002 * t820 + t1008 * t809 + t1025 * t961 - t993 * t962;
t1013 = mrSges(4,1) * t906 - mrSges(4,2) * t907 + Ifges(4,5) * t951 + Ifges(4,6) * t950 + Ifges(4,3) * t976 + pkin(3) * t860 + pkin(10) * t837 + t1001 * t822 + t1007 * t821 + t972 * t945 - t971 * t946;
t808 = -mrSges(3,1) * t964 + mrSges(3,3) * t953 + Ifges(3,4) * t983 - Ifges(3,2) * t984 + Ifges(3,6) * t992 - pkin(2) * t829 - t1026 * t961 + t993 * t963 - t1013;
t1018 = mrSges(2,1) * t989 - mrSges(2,2) * t990 + Ifges(2,3) * qJDD(1) + pkin(1) * t814 + t808 * t1030 + t806 * t1032 + t1037 * t819 + t998 * t804;
t802 = -mrSges(2,2) * g(3) - mrSges(2,3) * t989 + Ifges(2,5) * qJDD(1) - t1011 * Ifges(2,6) - t1003 * t808 + t1009 * t806 + (-t813 * t997 - t814 * t998) * pkin(8);
t801 = mrSges(2,1) * g(3) + mrSges(2,3) * t990 + t1011 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t813 - t997 * t804 + (pkin(8) * t819 + t1003 * t806 + t1009 * t808) * t998;
t1 = [-m(1) * g(1) + t1021; -m(1) * g(2) + t1035; (-m(1) - m(2)) * g(3) + t813; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1035 - t1004 * t801 + t1010 * t802; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1021 + t1004 * t802 + t1010 * t801; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1018; t1018; t804; t1013; t1012; -t1014; -t1017;];
tauJB  = t1;
