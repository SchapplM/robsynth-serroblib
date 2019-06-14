% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:51
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRR10_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR10_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:40:56
% EndTime: 2019-05-06 23:42:35
% DurationCPUTime: 66.17s
% Computational Cost: add. (1097794->401), mult. (2481085->518), div. (0->0), fcn. (2029015->14), ass. (0->166)
t988 = sin(pkin(6));
t1024 = pkin(8) * t988;
t990 = cos(pkin(6));
t1023 = t990 * g(3);
t994 = sin(qJ(2));
t1022 = t988 * t994;
t999 = cos(qJ(2));
t1021 = t988 * t999;
t1020 = t990 * t994;
t1019 = t990 * t999;
t1000 = cos(qJ(1));
t1001 = qJD(1) ^ 2;
t1015 = qJD(1) * t999;
t1012 = t988 * t1015;
t995 = sin(qJ(1));
t978 = t995 * g(1) - g(2) * t1000;
t969 = qJDD(1) * pkin(1) + t1001 * t1024 + t978;
t1014 = qJDD(1) * t988;
t979 = -g(1) * t1000 - g(2) * t995;
t970 = -pkin(1) * t1001 + pkin(8) * t1014 + t979;
t1017 = t969 * t1020 + t999 * t970;
t1016 = qJD(1) * t988;
t971 = (-pkin(2) * t999 - qJ(3) * t994) * t1016;
t983 = qJD(1) * t990 + qJD(2);
t981 = t983 ^ 2;
t982 = qJDD(1) * t990 + qJDD(2);
t925 = -t981 * pkin(2) + t982 * qJ(3) + (-g(3) * t994 + t971 * t1015) * t988 + t1017;
t973 = (qJD(2) * t1015 + qJDD(1) * t994) * t988;
t1013 = t994 * t1016;
t974 = -qJD(2) * t1013 + t999 * t1014;
t926 = -t974 * pkin(2) - t1023 - t973 * qJ(3) + (-t969 + (pkin(2) * t994 - qJ(3) * t999) * t983 * qJD(1)) * t988;
t987 = sin(pkin(12));
t989 = cos(pkin(12));
t963 = t989 * t1013 + t983 * t987;
t890 = -0.2e1 * qJD(3) * t963 - t987 * t925 + t989 * t926;
t950 = t973 * t989 + t982 * t987;
t962 = -t987 * t1013 + t983 * t989;
t881 = (-t962 * t1012 - t950) * pkin(9) + (t962 * t963 - t974) * pkin(3) + t890;
t891 = 0.2e1 * qJD(3) * t962 + t989 * t925 + t987 * t926;
t949 = -t973 * t987 + t982 * t989;
t951 = -pkin(3) * t1012 - pkin(9) * t963;
t961 = t962 ^ 2;
t888 = -pkin(3) * t961 + pkin(9) * t949 + t951 * t1012 + t891;
t993 = sin(qJ(4));
t998 = cos(qJ(4));
t870 = t993 * t881 + t998 * t888;
t942 = t962 * t998 - t993 * t963;
t943 = t962 * t993 + t963 * t998;
t920 = -pkin(4) * t942 - pkin(10) * t943;
t966 = qJDD(4) - t974;
t977 = qJD(4) - t1012;
t976 = t977 ^ 2;
t868 = -pkin(4) * t976 + pkin(10) * t966 + t920 * t942 + t870;
t938 = -g(3) * t1021 + t969 * t1019 - t994 * t970;
t924 = -t982 * pkin(2) - t981 * qJ(3) + t971 * t1013 + qJDD(3) - t938;
t895 = -t949 * pkin(3) - t961 * pkin(9) + t963 * t951 + t924;
t910 = -t943 * qJD(4) + t949 * t998 - t993 * t950;
t911 = qJD(4) * t942 + t949 * t993 + t950 * t998;
t878 = (-t942 * t977 - t911) * pkin(10) + (t943 * t977 - t910) * pkin(4) + t895;
t992 = sin(qJ(5));
t997 = cos(qJ(5));
t863 = -t992 * t868 + t997 * t878;
t928 = -t943 * t992 + t977 * t997;
t894 = qJD(5) * t928 + t911 * t997 + t966 * t992;
t909 = qJDD(5) - t910;
t929 = t943 * t997 + t977 * t992;
t941 = qJD(5) - t942;
t861 = (t928 * t941 - t894) * pkin(11) + (t928 * t929 + t909) * pkin(5) + t863;
t864 = t997 * t868 + t992 * t878;
t893 = -qJD(5) * t929 - t911 * t992 + t966 * t997;
t914 = pkin(5) * t941 - pkin(11) * t929;
t927 = t928 ^ 2;
t862 = -pkin(5) * t927 + pkin(11) * t893 - t914 * t941 + t864;
t991 = sin(qJ(6));
t996 = cos(qJ(6));
t859 = t861 * t996 - t862 * t991;
t903 = t928 * t996 - t929 * t991;
t875 = qJD(6) * t903 + t893 * t991 + t894 * t996;
t904 = t928 * t991 + t929 * t996;
t889 = -mrSges(7,1) * t903 + mrSges(7,2) * t904;
t937 = qJD(6) + t941;
t896 = -mrSges(7,2) * t937 + mrSges(7,3) * t903;
t907 = qJDD(6) + t909;
t855 = m(7) * t859 + mrSges(7,1) * t907 - mrSges(7,3) * t875 - t889 * t904 + t896 * t937;
t860 = t861 * t991 + t862 * t996;
t874 = -qJD(6) * t904 + t893 * t996 - t894 * t991;
t897 = mrSges(7,1) * t937 - mrSges(7,3) * t904;
t856 = m(7) * t860 - mrSges(7,2) * t907 + t874 * mrSges(7,3) + t889 * t903 - t897 * t937;
t847 = t996 * t855 + t991 * t856;
t905 = -mrSges(6,1) * t928 + mrSges(6,2) * t929;
t912 = -mrSges(6,2) * t941 + mrSges(6,3) * t928;
t845 = m(6) * t863 + mrSges(6,1) * t909 - mrSges(6,3) * t894 - t905 * t929 + t912 * t941 + t847;
t1008 = -t855 * t991 + t996 * t856;
t913 = mrSges(6,1) * t941 - mrSges(6,3) * t929;
t846 = m(6) * t864 - mrSges(6,2) * t909 + mrSges(6,3) * t893 + t905 * t928 - t913 * t941 + t1008;
t841 = -t845 * t992 + t997 * t846;
t919 = -mrSges(5,1) * t942 + mrSges(5,2) * t943;
t931 = mrSges(5,1) * t977 - mrSges(5,3) * t943;
t838 = m(5) * t870 - mrSges(5,2) * t966 + mrSges(5,3) * t910 + t919 * t942 - t931 * t977 + t841;
t869 = t881 * t998 - t993 * t888;
t867 = -pkin(4) * t966 - pkin(10) * t976 + t943 * t920 - t869;
t865 = -pkin(5) * t893 - pkin(11) * t927 + t914 * t929 + t867;
t1007 = m(7) * t865 - t874 * mrSges(7,1) + mrSges(7,2) * t875 - t903 * t896 + t897 * t904;
t857 = -m(6) * t867 + t893 * mrSges(6,1) - mrSges(6,2) * t894 + t928 * t912 - t913 * t929 - t1007;
t930 = -mrSges(5,2) * t977 + mrSges(5,3) * t942;
t851 = m(5) * t869 + mrSges(5,1) * t966 - mrSges(5,3) * t911 - t919 * t943 + t930 * t977 + t857;
t830 = t993 * t838 + t998 * t851;
t944 = -mrSges(4,1) * t962 + mrSges(4,2) * t963;
t947 = mrSges(4,2) * t1012 + mrSges(4,3) * t962;
t828 = m(4) * t890 - mrSges(4,1) * t974 - mrSges(4,3) * t950 - t947 * t1012 - t944 * t963 + t830;
t1009 = t998 * t838 - t851 * t993;
t948 = -mrSges(4,1) * t1012 - mrSges(4,3) * t963;
t829 = m(4) * t891 + mrSges(4,2) * t974 + mrSges(4,3) * t949 + t948 * t1012 + t944 * t962 + t1009;
t1010 = -t828 * t987 + t989 * t829;
t939 = -g(3) * t1022 + t1017;
t967 = mrSges(3,1) * t983 - mrSges(3,3) * t1013;
t972 = (-mrSges(3,1) * t999 + mrSges(3,2) * t994) * t1016;
t820 = m(3) * t939 - mrSges(3,2) * t982 + mrSges(3,3) * t974 + t972 * t1012 - t967 * t983 + t1010;
t823 = t989 * t828 + t987 * t829;
t955 = -t988 * t969 - t1023;
t968 = -mrSges(3,2) * t983 + mrSges(3,3) * t1012;
t822 = m(3) * t955 - t974 * mrSges(3,1) + t973 * mrSges(3,2) + (t967 * t994 - t968 * t999) * t1016 + t823;
t840 = t997 * t845 + t992 * t846;
t1004 = m(5) * t895 - t910 * mrSges(5,1) + mrSges(5,2) * t911 - t942 * t930 + t931 * t943 + t840;
t839 = m(4) * t924 - t949 * mrSges(4,1) + mrSges(4,2) * t950 - t962 * t947 + t948 * t963 + t1004;
t835 = m(3) * t938 + mrSges(3,1) * t982 - mrSges(3,3) * t973 - t972 * t1013 + t968 * t983 - t839;
t810 = t835 * t1019 + t820 * t1020 - t822 * t988;
t807 = m(2) * t978 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1001 + t810;
t815 = t999 * t820 - t835 * t994;
t813 = m(2) * t979 - mrSges(2,1) * t1001 - qJDD(1) * mrSges(2,2) + t815;
t1018 = t1000 * t807 + t995 * t813;
t809 = t835 * t1021 + t820 * t1022 + t990 * t822;
t1011 = t1000 * t813 - t807 * t995;
t882 = Ifges(7,5) * t904 + Ifges(7,6) * t903 + Ifges(7,3) * t937;
t884 = Ifges(7,1) * t904 + Ifges(7,4) * t903 + Ifges(7,5) * t937;
t848 = -mrSges(7,1) * t865 + mrSges(7,3) * t860 + Ifges(7,4) * t875 + Ifges(7,2) * t874 + Ifges(7,6) * t907 - t882 * t904 + t884 * t937;
t883 = Ifges(7,4) * t904 + Ifges(7,2) * t903 + Ifges(7,6) * t937;
t849 = mrSges(7,2) * t865 - mrSges(7,3) * t859 + Ifges(7,1) * t875 + Ifges(7,4) * t874 + Ifges(7,5) * t907 + t882 * t903 - t883 * t937;
t898 = Ifges(6,5) * t929 + Ifges(6,6) * t928 + Ifges(6,3) * t941;
t900 = Ifges(6,1) * t929 + Ifges(6,4) * t928 + Ifges(6,5) * t941;
t831 = -mrSges(6,1) * t867 + mrSges(6,3) * t864 + Ifges(6,4) * t894 + Ifges(6,2) * t893 + Ifges(6,6) * t909 - pkin(5) * t1007 + pkin(11) * t1008 + t996 * t848 + t991 * t849 - t929 * t898 + t941 * t900;
t899 = Ifges(6,4) * t929 + Ifges(6,2) * t928 + Ifges(6,6) * t941;
t832 = mrSges(6,2) * t867 - mrSges(6,3) * t863 + Ifges(6,1) * t894 + Ifges(6,4) * t893 + Ifges(6,5) * t909 - pkin(11) * t847 - t848 * t991 + t849 * t996 + t898 * t928 - t899 * t941;
t915 = Ifges(5,5) * t943 + Ifges(5,6) * t942 + Ifges(5,3) * t977;
t916 = Ifges(5,4) * t943 + Ifges(5,2) * t942 + Ifges(5,6) * t977;
t816 = mrSges(5,2) * t895 - mrSges(5,3) * t869 + Ifges(5,1) * t911 + Ifges(5,4) * t910 + Ifges(5,5) * t966 - pkin(10) * t840 - t831 * t992 + t832 * t997 + t915 * t942 - t916 * t977;
t1005 = -mrSges(7,1) * t859 + mrSges(7,2) * t860 - Ifges(7,5) * t875 - Ifges(7,6) * t874 - Ifges(7,3) * t907 - t904 * t883 + t903 * t884;
t1002 = mrSges(6,1) * t863 - mrSges(6,2) * t864 + Ifges(6,5) * t894 + Ifges(6,6) * t893 + Ifges(6,3) * t909 + pkin(5) * t847 + t929 * t899 - t928 * t900 - t1005;
t917 = Ifges(5,1) * t943 + Ifges(5,4) * t942 + Ifges(5,5) * t977;
t824 = -mrSges(5,1) * t895 + mrSges(5,3) * t870 + Ifges(5,4) * t911 + Ifges(5,2) * t910 + Ifges(5,6) * t966 - pkin(4) * t840 - t943 * t915 + t977 * t917 - t1002;
t932 = Ifges(4,5) * t963 + Ifges(4,6) * t962 - Ifges(4,3) * t1012;
t934 = Ifges(4,1) * t963 + Ifges(4,4) * t962 - Ifges(4,5) * t1012;
t802 = -mrSges(4,1) * t924 + mrSges(4,3) * t891 + Ifges(4,4) * t950 + Ifges(4,2) * t949 - Ifges(4,6) * t974 - pkin(3) * t1004 + pkin(9) * t1009 - t934 * t1012 + t993 * t816 + t998 * t824 - t963 * t932;
t933 = Ifges(4,4) * t963 + Ifges(4,2) * t962 - Ifges(4,6) * t1012;
t805 = mrSges(4,2) * t924 - mrSges(4,3) * t890 + Ifges(4,1) * t950 + Ifges(4,4) * t949 - Ifges(4,5) * t974 - pkin(9) * t830 + t933 * t1012 + t816 * t998 - t824 * t993 + t932 * t962;
t953 = Ifges(3,6) * t983 + (Ifges(3,4) * t994 + Ifges(3,2) * t999) * t1016;
t954 = Ifges(3,5) * t983 + (Ifges(3,1) * t994 + Ifges(3,4) * t999) * t1016;
t799 = Ifges(3,5) * t973 + Ifges(3,6) * t974 + Ifges(3,3) * t982 + mrSges(3,1) * t938 - mrSges(3,2) * t939 + t987 * t805 + t989 * t802 - pkin(2) * t839 + qJ(3) * t1010 + (t953 * t994 - t954 * t999) * t1016;
t952 = Ifges(3,3) * t983 + (Ifges(3,5) * t994 + Ifges(3,6) * t999) * t1016;
t801 = mrSges(3,2) * t955 - mrSges(3,3) * t938 + Ifges(3,1) * t973 + Ifges(3,4) * t974 + Ifges(3,5) * t982 - qJ(3) * t823 + t952 * t1012 - t802 * t987 + t805 * t989 - t953 * t983;
t1003 = mrSges(5,1) * t869 - mrSges(5,2) * t870 + Ifges(5,5) * t911 + Ifges(5,6) * t910 + Ifges(5,3) * t966 + pkin(4) * t857 + pkin(10) * t841 + t997 * t831 + t992 * t832 + t943 * t916 - t942 * t917;
t804 = -t1003 - mrSges(4,1) * t890 + mrSges(4,2) * t891 - pkin(2) * t823 + (Ifges(3,2) + Ifges(4,3)) * t974 + t983 * t954 - t952 * t1013 + Ifges(3,6) * t982 + Ifges(3,4) * t973 + t962 * t934 - t963 * t933 - Ifges(4,6) * t949 - Ifges(4,5) * t950 - mrSges(3,1) * t955 + mrSges(3,3) * t939 - pkin(3) * t830;
t1006 = mrSges(2,1) * t978 - mrSges(2,2) * t979 + Ifges(2,3) * qJDD(1) + pkin(1) * t810 + t804 * t1021 + t801 * t1022 + t815 * t1024 + t990 * t799;
t797 = -mrSges(2,2) * g(3) - mrSges(2,3) * t978 + Ifges(2,5) * qJDD(1) - t1001 * Ifges(2,6) + t999 * t801 - t994 * t804 + (-t809 * t988 - t810 * t990) * pkin(8);
t796 = mrSges(2,1) * g(3) + mrSges(2,3) * t979 + t1001 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t809 - t988 * t799 + (pkin(8) * t815 + t801 * t994 + t804 * t999) * t990;
t1 = [-m(1) * g(1) + t1011; -m(1) * g(2) + t1018; (-m(1) - m(2)) * g(3) + t809; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1018 + t1000 * t797 - t995 * t796; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1011 + t1000 * t796 + t995 * t797; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1006; t1006; t799; t839; t1003; t1002; -t1005;];
tauJB  = t1;
