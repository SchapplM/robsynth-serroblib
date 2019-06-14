% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 08:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:01:58
% EndTime: 2019-05-07 08:02:31
% DurationCPUTime: 31.00s
% Computational Cost: add. (523671->378), mult. (1155080->477), div. (0->0), fcn. (923342->12), ass. (0->157)
t1028 = Ifges(6,1) + Ifges(7,1);
t1020 = Ifges(6,4) + Ifges(7,4);
t1019 = Ifges(6,5) + Ifges(7,5);
t1027 = Ifges(6,2) + Ifges(7,2);
t1018 = Ifges(6,6) + Ifges(7,6);
t1026 = Ifges(6,3) + Ifges(7,3);
t1002 = qJD(1) * qJD(2);
t976 = sin(pkin(6));
t981 = sin(qJ(2));
t985 = cos(qJ(2));
t963 = (-qJDD(1) * t985 + t1002 * t981) * t976;
t1004 = qJD(1) * t985;
t978 = cos(pkin(6));
t1013 = t978 * t981;
t1023 = pkin(8) * t976;
t982 = sin(qJ(1));
t986 = cos(qJ(1));
t967 = t982 * g(1) - g(2) * t986;
t987 = qJD(1) ^ 2;
t958 = qJDD(1) * pkin(1) + t1023 * t987 + t967;
t968 = -g(1) * t986 - g(2) * t982;
t959 = -pkin(1) * t987 + qJDD(1) * t1023 + t968;
t1006 = t958 * t1013 + t985 * t959;
t1005 = qJD(1) * t976;
t961 = (-pkin(2) * t985 - pkin(9) * t981) * t1005;
t971 = qJD(1) * t978 + qJD(2);
t969 = t971 ^ 2;
t970 = qJDD(1) * t978 + qJDD(2);
t918 = -t969 * pkin(2) + t970 * pkin(9) + (-g(3) * t981 + t1004 * t961) * t976 + t1006;
t1022 = t978 * g(3);
t962 = (qJDD(1) * t981 + t1002 * t985) * t976;
t919 = t963 * pkin(2) - t962 * pkin(9) - t1022 + (-t958 + (pkin(2) * t981 - pkin(9) * t985) * t971 * qJD(1)) * t976;
t980 = sin(qJ(3));
t984 = cos(qJ(3));
t884 = -t980 * t918 + t984 * t919;
t999 = t981 * t1005;
t951 = t971 * t984 - t980 * t999;
t932 = qJD(3) * t951 + t962 * t984 + t970 * t980;
t952 = t971 * t980 + t984 * t999;
t955 = qJDD(3) + t963;
t998 = t976 * t1004;
t966 = qJD(3) - t998;
t873 = (t951 * t966 - t932) * qJ(4) + (t951 * t952 + t955) * pkin(3) + t884;
t885 = t984 * t918 + t980 * t919;
t931 = -qJD(3) * t952 - t962 * t980 + t970 * t984;
t942 = pkin(3) * t966 - qJ(4) * t952;
t950 = t951 ^ 2;
t876 = -pkin(3) * t950 + qJ(4) * t931 - t942 * t966 + t885;
t975 = sin(pkin(11));
t977 = cos(pkin(11));
t939 = t951 * t975 + t952 * t977;
t867 = -0.2e1 * qJD(4) * t939 + t873 * t977 - t975 * t876;
t938 = t951 * t977 - t952 * t975;
t868 = 0.2e1 * qJD(4) * t938 + t975 * t873 + t977 * t876;
t913 = -pkin(4) * t938 - pkin(10) * t939;
t965 = t966 ^ 2;
t866 = -pkin(4) * t965 + pkin(10) * t955 + t913 * t938 + t868;
t1012 = t978 * t985;
t1014 = t976 * t985;
t933 = -g(3) * t1014 + t1012 * t958 - t981 * t959;
t917 = -t970 * pkin(2) - t969 * pkin(9) + t961 * t999 - t933;
t877 = -t931 * pkin(3) - t950 * qJ(4) + t952 * t942 + qJDD(4) + t917;
t906 = t931 * t977 - t932 * t975;
t907 = t931 * t975 + t932 * t977;
t871 = (-t938 * t966 - t907) * pkin(10) + (t939 * t966 - t906) * pkin(4) + t877;
t979 = sin(qJ(5));
t983 = cos(qJ(5));
t861 = -t979 * t866 + t983 * t871;
t921 = -t939 * t979 + t966 * t983;
t882 = qJD(5) * t921 + t907 * t983 + t955 * t979;
t905 = qJDD(5) - t906;
t922 = t939 * t983 + t966 * t979;
t937 = qJD(5) - t938;
t858 = -0.2e1 * qJD(6) * t922 + (t921 * t937 - t882) * qJ(6) + (t921 * t922 + t905) * pkin(5) + t861;
t898 = -mrSges(7,2) * t937 + mrSges(7,3) * t921;
t1001 = m(7) * t858 + t905 * mrSges(7,1) + t937 * t898;
t896 = -mrSges(7,1) * t921 + mrSges(7,2) * t922;
t897 = -mrSges(6,1) * t921 + mrSges(6,2) * t922;
t899 = -mrSges(6,2) * t937 + mrSges(6,3) * t921;
t848 = m(6) * t861 + t905 * mrSges(6,1) + t937 * t899 + (-t896 - t897) * t922 + (-mrSges(6,3) - mrSges(7,3)) * t882 + t1001;
t862 = t983 * t866 + t979 * t871;
t881 = -qJD(5) * t922 - t907 * t979 + t955 * t983;
t900 = pkin(5) * t937 - qJ(6) * t922;
t920 = t921 ^ 2;
t860 = -pkin(5) * t920 + qJ(6) * t881 + 0.2e1 * qJD(6) * t921 - t900 * t937 + t862;
t1000 = m(7) * t860 + t881 * mrSges(7,3) + t921 * t896;
t901 = mrSges(7,1) * t937 - mrSges(7,3) * t922;
t1007 = -mrSges(6,1) * t937 + mrSges(6,3) * t922 - t901;
t1021 = -mrSges(6,2) - mrSges(7,2);
t853 = m(6) * t862 + t881 * mrSges(6,3) + t1007 * t937 + t1021 * t905 + t921 * t897 + t1000;
t846 = -t848 * t979 + t983 * t853;
t912 = -mrSges(5,1) * t938 + mrSges(5,2) * t939;
t924 = mrSges(5,1) * t966 - mrSges(5,3) * t939;
t842 = m(5) * t868 - mrSges(5,2) * t955 + mrSges(5,3) * t906 + t912 * t938 - t924 * t966 + t846;
t865 = -pkin(4) * t955 - pkin(10) * t965 + t939 * t913 - t867;
t863 = -pkin(5) * t881 - qJ(6) * t920 + t900 * t922 + qJDD(6) + t865;
t993 = -m(7) * t863 + t881 * mrSges(7,1) + t921 * t898;
t854 = -m(6) * t865 + t881 * mrSges(6,1) + t1007 * t922 + t1021 * t882 + t921 * t899 + t993;
t923 = -mrSges(5,2) * t966 + mrSges(5,3) * t938;
t850 = m(5) * t867 + t955 * mrSges(5,1) - t907 * mrSges(5,3) - t939 * t912 + t966 * t923 + t854;
t835 = t975 * t842 + t977 * t850;
t1008 = t1019 * t937 + t1020 * t921 + t1028 * t922;
t1010 = -t1018 * t921 - t1019 * t922 - t1026 * t937;
t856 = t882 * mrSges(7,2) + t922 * t901 - t993;
t836 = -mrSges(6,1) * t865 + mrSges(6,3) * t862 - mrSges(7,1) * t863 + mrSges(7,3) * t860 - pkin(5) * t856 + qJ(6) * t1000 + (-qJ(6) * t901 + t1008) * t937 + t1010 * t922 + (-mrSges(7,2) * qJ(6) + t1018) * t905 + t1020 * t882 + t1027 * t881;
t1009 = -t1018 * t937 - t1020 * t922 - t1027 * t921;
t855 = -t882 * mrSges(7,3) - t922 * t896 + t1001;
t843 = mrSges(6,2) * t865 + mrSges(7,2) * t863 - mrSges(6,3) * t861 - mrSges(7,3) * t858 - qJ(6) * t855 + t1009 * t937 - t1010 * t921 + t1019 * t905 + t1020 * t881 + t1028 * t882;
t909 = Ifges(5,4) * t939 + Ifges(5,2) * t938 + Ifges(5,6) * t966;
t910 = Ifges(5,1) * t939 + Ifges(5,4) * t938 + Ifges(5,5) * t966;
t926 = Ifges(4,4) * t952 + Ifges(4,2) * t951 + Ifges(4,6) * t966;
t927 = Ifges(4,1) * t952 + Ifges(4,4) * t951 + Ifges(4,5) * t966;
t1025 = Ifges(4,5) * t932 + Ifges(4,6) * t931 + t952 * t926 - t951 * t927 + mrSges(4,1) * t884 - mrSges(4,2) * t885 + Ifges(5,5) * t907 + Ifges(5,6) * t906 + t939 * t909 - t938 * t910 + mrSges(5,1) * t867 - mrSges(5,2) * t868 + t979 * t843 + t983 * t836 + pkin(4) * t854 + pkin(10) * t846 + pkin(3) * t835 + (Ifges(4,3) + Ifges(5,3)) * t955;
t1024 = mrSges(6,1) * t861 + mrSges(7,1) * t858 - mrSges(6,2) * t862 - mrSges(7,2) * t860 + pkin(5) * t855 - t1008 * t921 - t1009 * t922 + t1018 * t881 + t1019 * t882 + t1026 * t905;
t1015 = t976 * t981;
t934 = -g(3) * t1015 + t1006;
t956 = mrSges(3,1) * t971 - mrSges(3,3) * t999;
t960 = (-mrSges(3,1) * t985 + mrSges(3,2) * t981) * t1005;
t940 = -mrSges(4,1) * t951 + mrSges(4,2) * t952;
t941 = -mrSges(4,2) * t966 + mrSges(4,3) * t951;
t833 = m(4) * t884 + mrSges(4,1) * t955 - mrSges(4,3) * t932 - t940 * t952 + t941 * t966 + t835;
t943 = mrSges(4,1) * t966 - mrSges(4,3) * t952;
t995 = t977 * t842 - t850 * t975;
t834 = m(4) * t885 - mrSges(4,2) * t955 + mrSges(4,3) * t931 + t940 * t951 - t943 * t966 + t995;
t996 = -t833 * t980 + t984 * t834;
t825 = m(3) * t934 - mrSges(3,2) * t970 - mrSges(3,3) * t963 - t956 * t971 + t960 * t998 + t996;
t828 = t984 * t833 + t980 * t834;
t947 = -t976 * t958 - t1022;
t957 = -mrSges(3,2) * t971 + mrSges(3,3) * t998;
t827 = m(3) * t947 + t963 * mrSges(3,1) + t962 * mrSges(3,2) + (t956 * t981 - t957 * t985) * t1005 + t828;
t845 = t983 * t848 + t979 * t853;
t844 = m(5) * t877 - t906 * mrSges(5,1) + mrSges(5,2) * t907 - t938 * t923 + t924 * t939 + t845;
t989 = -m(4) * t917 + t931 * mrSges(4,1) - mrSges(4,2) * t932 + t951 * t941 - t943 * t952 - t844;
t839 = m(3) * t933 + mrSges(3,1) * t970 - mrSges(3,3) * t962 + t957 * t971 - t960 * t999 + t989;
t814 = t839 * t1012 + t825 * t1013 - t827 * t976;
t811 = m(2) * t967 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t987 + t814;
t820 = t985 * t825 - t839 * t981;
t818 = m(2) * t968 - mrSges(2,1) * t987 - qJDD(1) * mrSges(2,2) + t820;
t1011 = t986 * t811 + t982 * t818;
t813 = t839 * t1014 + t825 * t1015 + t978 * t827;
t997 = -t811 * t982 + t986 * t818;
t908 = Ifges(5,5) * t939 + Ifges(5,6) * t938 + Ifges(5,3) * t966;
t821 = mrSges(5,2) * t877 - mrSges(5,3) * t867 + Ifges(5,1) * t907 + Ifges(5,4) * t906 + Ifges(5,5) * t955 - pkin(10) * t845 - t836 * t979 + t843 * t983 + t908 * t938 - t909 * t966;
t829 = -mrSges(5,1) * t877 + mrSges(5,3) * t868 + Ifges(5,4) * t907 + Ifges(5,2) * t906 + Ifges(5,6) * t955 - pkin(4) * t845 - t939 * t908 + t966 * t910 - t1024;
t925 = Ifges(4,5) * t952 + Ifges(4,6) * t951 + Ifges(4,3) * t966;
t809 = -mrSges(4,1) * t917 + mrSges(4,3) * t885 + Ifges(4,4) * t932 + Ifges(4,2) * t931 + Ifges(4,6) * t955 - pkin(3) * t844 + qJ(4) * t995 + t975 * t821 + t977 * t829 - t952 * t925 + t966 * t927;
t815 = mrSges(4,2) * t917 - mrSges(4,3) * t884 + Ifges(4,1) * t932 + Ifges(4,4) * t931 + Ifges(4,5) * t955 - qJ(4) * t835 + t821 * t977 - t829 * t975 + t925 * t951 - t926 * t966;
t945 = Ifges(3,6) * t971 + (Ifges(3,4) * t981 + Ifges(3,2) * t985) * t1005;
t946 = Ifges(3,5) * t971 + (Ifges(3,1) * t981 + Ifges(3,4) * t985) * t1005;
t804 = Ifges(3,5) * t962 - Ifges(3,6) * t963 + Ifges(3,3) * t970 + mrSges(3,1) * t933 - mrSges(3,2) * t934 + t980 * t815 + t984 * t809 + pkin(2) * t989 + pkin(9) * t996 + (t945 * t981 - t946 * t985) * t1005;
t944 = Ifges(3,3) * t971 + (Ifges(3,5) * t981 + Ifges(3,6) * t985) * t1005;
t806 = mrSges(3,2) * t947 - mrSges(3,3) * t933 + Ifges(3,1) * t962 - Ifges(3,4) * t963 + Ifges(3,5) * t970 - pkin(9) * t828 - t809 * t980 + t815 * t984 + t944 * t998 - t945 * t971;
t808 = -mrSges(3,1) * t947 + mrSges(3,3) * t934 + Ifges(3,4) * t962 - Ifges(3,2) * t963 + Ifges(3,6) * t970 - pkin(2) * t828 - t944 * t999 + t971 * t946 - t1025;
t991 = mrSges(2,1) * t967 - mrSges(2,2) * t968 + Ifges(2,3) * qJDD(1) + pkin(1) * t814 + t808 * t1014 + t806 * t1015 + t820 * t1023 + t978 * t804;
t802 = -mrSges(2,2) * g(3) - mrSges(2,3) * t967 + Ifges(2,5) * qJDD(1) - t987 * Ifges(2,6) + t985 * t806 - t981 * t808 + (-t813 * t976 - t814 * t978) * pkin(8);
t801 = mrSges(2,1) * g(3) + mrSges(2,3) * t968 + t987 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t813 - t976 * t804 + (pkin(8) * t820 + t806 * t981 + t808 * t985) * t978;
t1 = [-m(1) * g(1) + t997; -m(1) * g(2) + t1011; (-m(1) - m(2)) * g(3) + t813; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1011 - t982 * t801 + t986 * t802; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t997 + t986 * t801 + t982 * t802; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t991; t991; t804; t1025; t844; t1024; t856;];
tauJB  = t1;
