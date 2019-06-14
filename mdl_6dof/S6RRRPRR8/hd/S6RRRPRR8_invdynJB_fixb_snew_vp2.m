% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 12:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:18:10
% EndTime: 2019-05-07 12:19:49
% DurationCPUTime: 70.48s
% Computational Cost: add. (1152128->401), mult. (2537565->519), div. (0->0), fcn. (2064714->14), ass. (0->166)
t1005 = cos(qJ(2));
t1000 = sin(qJ(2));
t1023 = qJD(1) * t1000;
t994 = sin(pkin(6));
t981 = (qJD(2) * t1023 - qJDD(1) * t1005) * t994;
t1004 = cos(qJ(3));
t1022 = qJD(1) * t1005;
t996 = cos(pkin(6));
t1026 = t1000 * t996;
t1007 = qJD(1) ^ 2;
t1034 = pkin(8) * t994;
t1001 = sin(qJ(1));
t1006 = cos(qJ(1));
t985 = t1001 * g(1) - g(2) * t1006;
t976 = qJDD(1) * pkin(1) + t1007 * t1034 + t985;
t986 = -g(1) * t1006 - g(2) * t1001;
t977 = -pkin(1) * t1007 + qJDD(1) * t1034 + t986;
t1030 = t1005 * t977 + t976 * t1026;
t1029 = qJD(1) * t994;
t979 = (-pkin(2) * t1005 - pkin(9) * t1000) * t1029;
t989 = qJD(1) * t996 + qJD(2);
t987 = t989 ^ 2;
t988 = qJDD(1) * t996 + qJDD(2);
t933 = -t987 * pkin(2) + t988 * pkin(9) + (-g(3) * t1000 + t1022 * t979) * t994 + t1030;
t1033 = t996 * g(3);
t980 = (qJD(2) * t1022 + qJDD(1) * t1000) * t994;
t934 = t981 * pkin(2) - t980 * pkin(9) - t1033 + (-t976 + (pkin(2) * t1000 - pkin(9) * t1005) * t989 * qJD(1)) * t994;
t999 = sin(qJ(3));
t904 = t1004 * t934 - t999 * t933;
t1021 = t994 * t1023;
t969 = t1004 * t989 - t1021 * t999;
t948 = qJD(3) * t969 + t1004 * t980 + t988 * t999;
t970 = t1004 * t1021 + t989 * t999;
t973 = qJDD(3) + t981;
t1020 = t994 * t1022;
t984 = qJD(3) - t1020;
t889 = (t969 * t984 - t948) * qJ(4) + (t969 * t970 + t973) * pkin(3) + t904;
t905 = t1004 * t933 + t999 * t934;
t947 = -qJD(3) * t970 + t1004 * t988 - t980 * t999;
t959 = pkin(3) * t984 - qJ(4) * t970;
t968 = t969 ^ 2;
t896 = -pkin(3) * t968 + qJ(4) * t947 - t959 * t984 + t905;
t993 = sin(pkin(12));
t995 = cos(pkin(12));
t956 = t969 * t993 + t970 * t995;
t877 = -0.2e1 * qJD(4) * t956 + t889 * t995 - t993 * t896;
t1003 = cos(qJ(5));
t1002 = cos(qJ(6));
t955 = t969 * t995 - t993 * t970;
t878 = 0.2e1 * qJD(4) * t955 + t993 * t889 + t995 * t896;
t928 = -pkin(4) * t955 - pkin(10) * t956;
t983 = t984 ^ 2;
t876 = -pkin(4) * t983 + pkin(10) * t973 + t928 * t955 + t878;
t1024 = t1005 * t996;
t1025 = t1005 * t994;
t950 = -g(3) * t1025 - t1000 * t977 + t1024 * t976;
t932 = -t988 * pkin(2) - t987 * pkin(9) + t979 * t1021 - t950;
t898 = -t947 * pkin(3) - t968 * qJ(4) + t970 * t959 + qJDD(4) + t932;
t921 = t947 * t995 - t993 * t948;
t922 = t947 * t993 + t948 * t995;
t881 = (-t955 * t984 - t922) * pkin(10) + (t956 * t984 - t921) * pkin(4) + t898;
t998 = sin(qJ(5));
t871 = t1003 * t881 - t998 * t876;
t936 = t1003 * t984 - t956 * t998;
t901 = qJD(5) * t936 + t1003 * t922 + t973 * t998;
t920 = qJDD(5) - t921;
t937 = t1003 * t956 + t984 * t998;
t954 = qJD(5) - t955;
t869 = (t936 * t954 - t901) * pkin(11) + (t936 * t937 + t920) * pkin(5) + t871;
t872 = t1003 * t876 + t998 * t881;
t900 = -qJD(5) * t937 + t1003 * t973 - t922 * t998;
t917 = pkin(5) * t954 - pkin(11) * t937;
t935 = t936 ^ 2;
t870 = -pkin(5) * t935 + pkin(11) * t900 - t917 * t954 + t872;
t997 = sin(qJ(6));
t867 = t1002 * t869 - t870 * t997;
t911 = t1002 * t936 - t937 * t997;
t886 = qJD(6) * t911 + t1002 * t901 + t900 * t997;
t912 = t1002 * t937 + t936 * t997;
t897 = -mrSges(7,1) * t911 + mrSges(7,2) * t912;
t949 = qJD(6) + t954;
t902 = -mrSges(7,2) * t949 + mrSges(7,3) * t911;
t918 = qJDD(6) + t920;
t863 = m(7) * t867 + mrSges(7,1) * t918 - mrSges(7,3) * t886 - t897 * t912 + t902 * t949;
t868 = t1002 * t870 + t869 * t997;
t885 = -qJD(6) * t912 + t1002 * t900 - t901 * t997;
t903 = mrSges(7,1) * t949 - mrSges(7,3) * t912;
t864 = m(7) * t868 - mrSges(7,2) * t918 + mrSges(7,3) * t885 + t897 * t911 - t903 * t949;
t855 = t1002 * t863 + t997 * t864;
t913 = -mrSges(6,1) * t936 + mrSges(6,2) * t937;
t915 = -mrSges(6,2) * t954 + mrSges(6,3) * t936;
t853 = m(6) * t871 + mrSges(6,1) * t920 - mrSges(6,3) * t901 - t913 * t937 + t915 * t954 + t855;
t1017 = t1002 * t864 - t863 * t997;
t916 = mrSges(6,1) * t954 - mrSges(6,3) * t937;
t854 = m(6) * t872 - mrSges(6,2) * t920 + mrSges(6,3) * t900 + t913 * t936 - t916 * t954 + t1017;
t849 = t1003 * t854 - t853 * t998;
t927 = -mrSges(5,1) * t955 + mrSges(5,2) * t956;
t939 = mrSges(5,1) * t984 - mrSges(5,3) * t956;
t846 = m(5) * t878 - mrSges(5,2) * t973 + mrSges(5,3) * t921 + t927 * t955 - t939 * t984 + t849;
t875 = -pkin(4) * t973 - pkin(10) * t983 + t956 * t928 - t877;
t873 = -pkin(5) * t900 - pkin(11) * t935 + t917 * t937 + t875;
t1013 = m(7) * t873 - t885 * mrSges(7,1) + mrSges(7,2) * t886 - t911 * t902 + t903 * t912;
t865 = -m(6) * t875 + t900 * mrSges(6,1) - mrSges(6,2) * t901 + t936 * t915 - t916 * t937 - t1013;
t938 = -mrSges(5,2) * t984 + mrSges(5,3) * t955;
t859 = m(5) * t877 + mrSges(5,1) * t973 - mrSges(5,3) * t922 - t927 * t956 + t938 * t984 + t865;
t838 = t993 * t846 + t995 * t859;
t890 = Ifges(7,5) * t912 + Ifges(7,6) * t911 + Ifges(7,3) * t949;
t892 = Ifges(7,1) * t912 + Ifges(7,4) * t911 + Ifges(7,5) * t949;
t856 = -mrSges(7,1) * t873 + mrSges(7,3) * t868 + Ifges(7,4) * t886 + Ifges(7,2) * t885 + Ifges(7,6) * t918 - t890 * t912 + t892 * t949;
t891 = Ifges(7,4) * t912 + Ifges(7,2) * t911 + Ifges(7,6) * t949;
t857 = mrSges(7,2) * t873 - mrSges(7,3) * t867 + Ifges(7,1) * t886 + Ifges(7,4) * t885 + Ifges(7,5) * t918 + t890 * t911 - t891 * t949;
t906 = Ifges(6,5) * t937 + Ifges(6,6) * t936 + Ifges(6,3) * t954;
t908 = Ifges(6,1) * t937 + Ifges(6,4) * t936 + Ifges(6,5) * t954;
t839 = -mrSges(6,1) * t875 + mrSges(6,3) * t872 + Ifges(6,4) * t901 + Ifges(6,2) * t900 + Ifges(6,6) * t920 - pkin(5) * t1013 + pkin(11) * t1017 + t1002 * t856 + t997 * t857 - t937 * t906 + t954 * t908;
t907 = Ifges(6,4) * t937 + Ifges(6,2) * t936 + Ifges(6,6) * t954;
t840 = mrSges(6,2) * t875 - mrSges(6,3) * t871 + Ifges(6,1) * t901 + Ifges(6,4) * t900 + Ifges(6,5) * t920 - pkin(11) * t855 + t1002 * t857 - t856 * t997 + t906 * t936 - t907 * t954;
t924 = Ifges(5,4) * t956 + Ifges(5,2) * t955 + Ifges(5,6) * t984;
t925 = Ifges(5,1) * t956 + Ifges(5,4) * t955 + Ifges(5,5) * t984;
t942 = Ifges(4,4) * t970 + Ifges(4,2) * t969 + Ifges(4,6) * t984;
t943 = Ifges(4,1) * t970 + Ifges(4,4) * t969 + Ifges(4,5) * t984;
t1035 = Ifges(4,5) * t948 + Ifges(4,6) * t947 + t970 * t942 - t969 * t943 + mrSges(4,1) * t904 - mrSges(4,2) * t905 + Ifges(5,5) * t922 + Ifges(5,6) * t921 + t956 * t924 - t955 * t925 + mrSges(5,1) * t877 - mrSges(5,2) * t878 + t998 * t840 + t1003 * t839 + pkin(4) * t865 + pkin(10) * t849 + pkin(3) * t838 + (Ifges(4,3) + Ifges(5,3)) * t973;
t957 = -mrSges(4,1) * t969 + mrSges(4,2) * t970;
t958 = -mrSges(4,2) * t984 + mrSges(4,3) * t969;
t836 = m(4) * t904 + mrSges(4,1) * t973 - mrSges(4,3) * t948 - t957 * t970 + t958 * t984 + t838;
t1018 = t995 * t846 - t859 * t993;
t960 = mrSges(4,1) * t984 - mrSges(4,3) * t970;
t837 = m(4) * t905 - mrSges(4,2) * t973 + mrSges(4,3) * t947 + t957 * t969 - t960 * t984 + t1018;
t1019 = t1004 * t837 - t836 * t999;
t1027 = t1000 * t994;
t951 = -g(3) * t1027 + t1030;
t974 = mrSges(3,1) * t989 - mrSges(3,3) * t1021;
t978 = (-mrSges(3,1) * t1005 + mrSges(3,2) * t1000) * t1029;
t828 = m(3) * t951 - mrSges(3,2) * t988 - mrSges(3,3) * t981 + t1020 * t978 - t974 * t989 + t1019;
t831 = t1004 * t836 + t999 * t837;
t964 = -t994 * t976 - t1033;
t975 = -mrSges(3,2) * t989 + mrSges(3,3) * t1020;
t830 = m(3) * t964 + t981 * mrSges(3,1) + t980 * mrSges(3,2) + (t1000 * t974 - t1005 * t975) * t1029 + t831;
t848 = t1003 * t853 + t998 * t854;
t847 = m(5) * t898 - t921 * mrSges(5,1) + mrSges(5,2) * t922 - t955 * t938 + t939 * t956 + t848;
t1010 = -m(4) * t932 + t947 * mrSges(4,1) - mrSges(4,2) * t948 + t969 * t958 - t960 * t970 - t847;
t843 = m(3) * t950 + mrSges(3,1) * t988 - mrSges(3,3) * t980 - t1021 * t978 + t975 * t989 + t1010;
t818 = t843 * t1024 + t828 * t1026 - t830 * t994;
t815 = m(2) * t985 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1007 + t818;
t823 = -t1000 * t843 + t1005 * t828;
t821 = m(2) * t986 - mrSges(2,1) * t1007 - qJDD(1) * mrSges(2,2) + t823;
t1031 = t1001 * t821 + t1006 * t815;
t817 = t843 * t1025 + t828 * t1027 + t996 * t830;
t1015 = -t1001 * t815 + t1006 * t821;
t923 = Ifges(5,5) * t956 + Ifges(5,6) * t955 + Ifges(5,3) * t984;
t824 = mrSges(5,2) * t898 - mrSges(5,3) * t877 + Ifges(5,1) * t922 + Ifges(5,4) * t921 + Ifges(5,5) * t973 - pkin(10) * t848 + t1003 * t840 - t839 * t998 + t923 * t955 - t924 * t984;
t1011 = -mrSges(7,1) * t867 + mrSges(7,2) * t868 - Ifges(7,5) * t886 - Ifges(7,6) * t885 - Ifges(7,3) * t918 - t912 * t891 + t911 * t892;
t1009 = mrSges(6,1) * t871 - mrSges(6,2) * t872 + Ifges(6,5) * t901 + Ifges(6,6) * t900 + Ifges(6,3) * t920 + pkin(5) * t855 + t937 * t907 - t936 * t908 - t1011;
t832 = -mrSges(5,1) * t898 + mrSges(5,3) * t878 + Ifges(5,4) * t922 + Ifges(5,2) * t921 + Ifges(5,6) * t973 - pkin(4) * t848 - t956 * t923 + t984 * t925 - t1009;
t941 = Ifges(4,5) * t970 + Ifges(4,6) * t969 + Ifges(4,3) * t984;
t810 = -mrSges(4,1) * t932 + mrSges(4,3) * t905 + Ifges(4,4) * t948 + Ifges(4,2) * t947 + Ifges(4,6) * t973 - pkin(3) * t847 + qJ(4) * t1018 + t993 * t824 + t995 * t832 - t970 * t941 + t984 * t943;
t813 = mrSges(4,2) * t932 - mrSges(4,3) * t904 + Ifges(4,1) * t948 + Ifges(4,4) * t947 + Ifges(4,5) * t973 - qJ(4) * t838 + t824 * t995 - t832 * t993 + t941 * t969 - t942 * t984;
t962 = Ifges(3,6) * t989 + (Ifges(3,4) * t1000 + Ifges(3,2) * t1005) * t1029;
t963 = Ifges(3,5) * t989 + (Ifges(3,1) * t1000 + Ifges(3,4) * t1005) * t1029;
t807 = Ifges(3,5) * t980 - Ifges(3,6) * t981 + Ifges(3,3) * t988 + mrSges(3,1) * t950 - mrSges(3,2) * t951 + t999 * t813 + t1004 * t810 + pkin(2) * t1010 + pkin(9) * t1019 + (t1000 * t962 - t1005 * t963) * t1029;
t961 = Ifges(3,3) * t989 + (Ifges(3,5) * t1000 + Ifges(3,6) * t1005) * t1029;
t809 = mrSges(3,2) * t964 - mrSges(3,3) * t950 + Ifges(3,1) * t980 - Ifges(3,4) * t981 + Ifges(3,5) * t988 - pkin(9) * t831 + t1004 * t813 + t1020 * t961 - t810 * t999 - t962 * t989;
t812 = -mrSges(3,1) * t964 + mrSges(3,3) * t951 + Ifges(3,4) * t980 - Ifges(3,2) * t981 + Ifges(3,6) * t988 - pkin(2) * t831 - t1021 * t961 + t989 * t963 - t1035;
t1012 = mrSges(2,1) * t985 - mrSges(2,2) * t986 + Ifges(2,3) * qJDD(1) + pkin(1) * t818 + t812 * t1025 + t809 * t1027 + t823 * t1034 + t996 * t807;
t805 = -mrSges(2,2) * g(3) - mrSges(2,3) * t985 + Ifges(2,5) * qJDD(1) - t1007 * Ifges(2,6) - t1000 * t812 + t1005 * t809 + (-t817 * t994 - t818 * t996) * pkin(8);
t804 = mrSges(2,1) * g(3) + mrSges(2,3) * t986 + t1007 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t817 - t994 * t807 + (pkin(8) * t823 + t1000 * t809 + t1005 * t812) * t996;
t1 = [-m(1) * g(1) + t1015; -m(1) * g(2) + t1031; (-m(1) - m(2)) * g(3) + t817; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1031 - t1001 * t804 + t1006 * t805; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1015 + t1001 * t805 + t1006 * t804; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1012; t1012; t807; t1035; t847; t1009; -t1011;];
tauJB  = t1;
