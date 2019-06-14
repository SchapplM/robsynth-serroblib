% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 23:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPR11_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 23:27:51
% EndTime: 2019-05-07 23:29:44
% DurationCPUTime: 83.10s
% Computational Cost: add. (1425752->401), mult. (3067810->519), div. (0->0), fcn. (2484560->14), ass. (0->166)
t1015 = qJD(1) * qJD(2);
t988 = sin(pkin(6));
t994 = sin(qJ(2));
t999 = cos(qJ(2));
t974 = (-qJDD(1) * t999 + t1015 * t994) * t988;
t1016 = qJD(1) * t999;
t990 = cos(pkin(6));
t1021 = t990 * t994;
t1001 = qJD(1) ^ 2;
t1026 = pkin(8) * t988;
t1000 = cos(qJ(1));
t995 = sin(qJ(1));
t979 = t995 * g(1) - g(2) * t1000;
t969 = qJDD(1) * pkin(1) + t1001 * t1026 + t979;
t980 = -g(1) * t1000 - g(2) * t995;
t970 = -pkin(1) * t1001 + qJDD(1) * t1026 + t980;
t1018 = t969 * t1021 + t999 * t970;
t1017 = qJD(1) * t988;
t972 = (-pkin(2) * t999 - pkin(9) * t994) * t1017;
t983 = qJD(1) * t990 + qJD(2);
t981 = t983 ^ 2;
t982 = qJDD(1) * t990 + qJDD(2);
t921 = -t981 * pkin(2) + t982 * pkin(9) + (-g(3) * t994 + t1016 * t972) * t988 + t1018;
t1025 = t990 * g(3);
t973 = (qJDD(1) * t994 + t1015 * t999) * t988;
t922 = t974 * pkin(2) - t973 * pkin(9) - t1025 + (-t969 + (pkin(2) * t994 - pkin(9) * t999) * t983 * qJD(1)) * t988;
t993 = sin(qJ(3));
t998 = cos(qJ(3));
t898 = t998 * t921 + t993 * t922;
t1014 = t994 * t1017;
t961 = -t993 * t1014 + t983 * t998;
t962 = t1014 * t998 + t983 * t993;
t946 = -pkin(3) * t961 - pkin(10) * t962;
t966 = qJDD(3) + t974;
t1013 = t988 * t1016;
t978 = qJD(3) - t1013;
t976 = t978 ^ 2;
t888 = -pkin(3) * t976 + pkin(10) * t966 + t946 * t961 + t898;
t1020 = t990 * t999;
t1022 = t988 * t999;
t943 = -g(3) * t1022 + t1020 * t969 - t994 * t970;
t920 = -t982 * pkin(2) - t981 * pkin(9) + t972 * t1014 - t943;
t941 = -t962 * qJD(3) - t993 * t973 + t982 * t998;
t942 = qJD(3) * t961 + t973 * t998 + t982 * t993;
t891 = (-t961 * t978 - t942) * pkin(10) + (t962 * t978 - t941) * pkin(3) + t920;
t992 = sin(qJ(4));
t997 = cos(qJ(4));
t876 = -t992 * t888 + t997 * t891;
t948 = -t962 * t992 + t978 * t997;
t908 = qJD(4) * t948 + t942 * t997 + t966 * t992;
t939 = qJDD(4) - t941;
t949 = t962 * t997 + t978 * t992;
t960 = qJD(4) - t961;
t867 = (t948 * t960 - t908) * qJ(5) + (t948 * t949 + t939) * pkin(4) + t876;
t877 = t997 * t888 + t992 * t891;
t907 = -qJD(4) * t949 - t942 * t992 + t966 * t997;
t931 = pkin(4) * t960 - qJ(5) * t949;
t947 = t948 ^ 2;
t869 = -pkin(4) * t947 + qJ(5) * t907 - t931 * t960 + t877;
t987 = sin(pkin(12));
t989 = cos(pkin(12));
t927 = t948 * t987 + t949 * t989;
t861 = -0.2e1 * qJD(5) * t927 + t989 * t867 - t987 * t869;
t894 = t907 * t987 + t908 * t989;
t926 = t948 * t989 - t949 * t987;
t859 = (t926 * t960 - t894) * pkin(11) + (t926 * t927 + t939) * pkin(5) + t861;
t862 = 0.2e1 * qJD(5) * t926 + t987 * t867 + t989 * t869;
t893 = t907 * t989 - t908 * t987;
t911 = pkin(5) * t960 - pkin(11) * t927;
t925 = t926 ^ 2;
t860 = -pkin(5) * t925 + pkin(11) * t893 - t911 * t960 + t862;
t991 = sin(qJ(6));
t996 = cos(qJ(6));
t856 = t859 * t996 - t860 * t991;
t857 = t859 * t991 + t860 * t996;
t904 = t926 * t991 + t927 * t996;
t874 = -qJD(6) * t904 + t893 * t996 - t894 * t991;
t903 = t926 * t996 - t927 * t991;
t875 = qJD(6) * t903 + t893 * t991 + t894 * t996;
t958 = qJD(6) + t960;
t881 = Ifges(7,4) * t904 + Ifges(7,2) * t903 + Ifges(7,6) * t958;
t882 = Ifges(7,1) * t904 + Ifges(7,4) * t903 + Ifges(7,5) * t958;
t934 = qJDD(6) + t939;
t1005 = -mrSges(7,1) * t856 + mrSges(7,2) * t857 - Ifges(7,5) * t875 - Ifges(7,6) * t874 - Ifges(7,3) * t934 - t904 * t881 + t903 * t882;
t885 = -mrSges(7,1) * t903 + mrSges(7,2) * t904;
t895 = -mrSges(7,2) * t958 + mrSges(7,3) * t903;
t850 = m(7) * t856 + mrSges(7,1) * t934 - mrSges(7,3) * t875 - t885 * t904 + t895 * t958;
t896 = mrSges(7,1) * t958 - mrSges(7,3) * t904;
t851 = m(7) * t857 - mrSges(7,2) * t934 + mrSges(7,3) * t874 + t885 * t903 - t896 * t958;
t844 = t996 * t850 + t991 * t851;
t905 = -mrSges(6,1) * t926 + mrSges(6,2) * t927;
t909 = -mrSges(6,2) * t960 + mrSges(6,3) * t926;
t842 = m(6) * t861 + mrSges(6,1) * t939 - mrSges(6,3) * t894 - t905 * t927 + t909 * t960 + t844;
t1009 = -t850 * t991 + t996 * t851;
t910 = mrSges(6,1) * t960 - mrSges(6,3) * t927;
t843 = m(6) * t862 - mrSges(6,2) * t939 + mrSges(6,3) * t893 + t905 * t926 - t910 * t960 + t1009;
t838 = t989 * t842 + t987 * t843;
t900 = Ifges(6,4) * t927 + Ifges(6,2) * t926 + Ifges(6,6) * t960;
t901 = Ifges(6,1) * t927 + Ifges(6,4) * t926 + Ifges(6,5) * t960;
t913 = Ifges(5,4) * t949 + Ifges(5,2) * t948 + Ifges(5,6) * t960;
t914 = Ifges(5,1) * t949 + Ifges(5,4) * t948 + Ifges(5,5) * t960;
t1027 = mrSges(5,1) * t876 + mrSges(6,1) * t861 - mrSges(5,2) * t877 - mrSges(6,2) * t862 + Ifges(5,5) * t908 + Ifges(6,5) * t894 + Ifges(5,6) * t907 + Ifges(6,6) * t893 + pkin(4) * t838 + pkin(5) * t844 + (Ifges(5,3) + Ifges(6,3)) * t939 + t927 * t900 - t926 * t901 + t949 * t913 - t948 * t914 - t1005;
t1023 = t988 * t994;
t928 = -mrSges(5,1) * t948 + mrSges(5,2) * t949;
t930 = -mrSges(5,2) * t960 + mrSges(5,3) * t948;
t836 = m(5) * t876 + mrSges(5,1) * t939 - mrSges(5,3) * t908 - t928 * t949 + t930 * t960 + t838;
t1010 = -t842 * t987 + t989 * t843;
t932 = mrSges(5,1) * t960 - mrSges(5,3) * t949;
t837 = m(5) * t877 - mrSges(5,2) * t939 + mrSges(5,3) * t907 + t928 * t948 - t932 * t960 + t1010;
t832 = -t836 * t992 + t997 * t837;
t945 = -mrSges(4,1) * t961 + mrSges(4,2) * t962;
t951 = mrSges(4,1) * t978 - mrSges(4,3) * t962;
t830 = m(4) * t898 - mrSges(4,2) * t966 + mrSges(4,3) * t941 + t945 * t961 - t951 * t978 + t832;
t897 = -t993 * t921 + t922 * t998;
t887 = -pkin(3) * t966 - pkin(10) * t976 + t962 * t946 - t897;
t878 = -pkin(4) * t907 - qJ(5) * t947 + t949 * t931 + qJDD(5) + t887;
t864 = -pkin(5) * t893 - pkin(11) * t925 + t911 * t927 + t878;
t1008 = m(7) * t864 - t874 * mrSges(7,1) + t875 * mrSges(7,2) - t903 * t895 + t904 * t896;
t858 = m(6) * t878 - t893 * mrSges(6,1) + mrSges(6,2) * t894 - t926 * t909 + t910 * t927 + t1008;
t854 = -m(5) * t887 + t907 * mrSges(5,1) - mrSges(5,2) * t908 + t948 * t930 - t932 * t949 - t858;
t950 = -mrSges(4,2) * t978 + mrSges(4,3) * t961;
t853 = m(4) * t897 + mrSges(4,1) * t966 - mrSges(4,3) * t942 - t945 * t962 + t950 * t978 + t854;
t1011 = t998 * t830 - t853 * t993;
t944 = -g(3) * t1023 + t1018;
t967 = mrSges(3,1) * t983 - mrSges(3,3) * t1014;
t971 = (-mrSges(3,1) * t999 + mrSges(3,2) * t994) * t1017;
t821 = m(3) * t944 - mrSges(3,2) * t982 - mrSges(3,3) * t974 + t1013 * t971 - t967 * t983 + t1011;
t824 = t993 * t830 + t998 * t853;
t955 = -t988 * t969 - t1025;
t968 = -mrSges(3,2) * t983 + mrSges(3,3) * t1013;
t823 = m(3) * t955 + t974 * mrSges(3,1) + t973 * mrSges(3,2) + (t967 * t994 - t968 * t999) * t1017 + t824;
t831 = t836 * t997 + t837 * t992;
t1004 = -m(4) * t920 + t941 * mrSges(4,1) - mrSges(4,2) * t942 + t961 * t950 - t951 * t962 - t831;
t827 = m(3) * t943 + mrSges(3,1) * t982 - mrSges(3,3) * t973 - t1014 * t971 + t968 * t983 + t1004;
t809 = t827 * t1020 + t821 * t1021 - t823 * t988;
t806 = m(2) * t979 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t1001 + t809;
t814 = t999 * t821 - t827 * t994;
t812 = m(2) * t980 - mrSges(2,1) * t1001 - qJDD(1) * mrSges(2,2) + t814;
t1019 = t1000 * t806 + t995 * t812;
t808 = t827 * t1022 + t821 * t1023 + t990 * t823;
t1012 = t1000 * t812 - t806 * t995;
t880 = Ifges(7,5) * t904 + Ifges(7,6) * t903 + Ifges(7,3) * t958;
t845 = -mrSges(7,1) * t864 + mrSges(7,3) * t857 + Ifges(7,4) * t875 + Ifges(7,2) * t874 + Ifges(7,6) * t934 - t880 * t904 + t882 * t958;
t846 = mrSges(7,2) * t864 - mrSges(7,3) * t856 + Ifges(7,1) * t875 + Ifges(7,4) * t874 + Ifges(7,5) * t934 + t880 * t903 - t881 * t958;
t899 = Ifges(6,5) * t927 + Ifges(6,6) * t926 + Ifges(6,3) * t960;
t833 = -mrSges(6,1) * t878 + mrSges(6,3) * t862 + Ifges(6,4) * t894 + Ifges(6,2) * t893 + Ifges(6,6) * t939 - pkin(5) * t1008 + pkin(11) * t1009 + t996 * t845 + t991 * t846 - t927 * t899 + t960 * t901;
t834 = mrSges(6,2) * t878 - mrSges(6,3) * t861 + Ifges(6,1) * t894 + Ifges(6,4) * t893 + Ifges(6,5) * t939 - pkin(11) * t844 - t845 * t991 + t846 * t996 + t899 * t926 - t900 * t960;
t912 = Ifges(5,5) * t949 + Ifges(5,6) * t948 + Ifges(5,3) * t960;
t816 = -mrSges(5,1) * t887 + mrSges(5,3) * t877 + Ifges(5,4) * t908 + Ifges(5,2) * t907 + Ifges(5,6) * t939 - pkin(4) * t858 + qJ(5) * t1010 + t989 * t833 + t987 * t834 - t949 * t912 + t960 * t914;
t817 = mrSges(5,2) * t887 - mrSges(5,3) * t876 + Ifges(5,1) * t908 + Ifges(5,4) * t907 + Ifges(5,5) * t939 - qJ(5) * t838 - t833 * t987 + t834 * t989 + t912 * t948 - t913 * t960;
t935 = Ifges(4,5) * t962 + Ifges(4,6) * t961 + Ifges(4,3) * t978;
t936 = Ifges(4,4) * t962 + Ifges(4,2) * t961 + Ifges(4,6) * t978;
t804 = mrSges(4,2) * t920 - mrSges(4,3) * t897 + Ifges(4,1) * t942 + Ifges(4,4) * t941 + Ifges(4,5) * t966 - pkin(10) * t831 - t816 * t992 + t817 * t997 + t935 * t961 - t936 * t978;
t937 = Ifges(4,1) * t962 + Ifges(4,4) * t961 + Ifges(4,5) * t978;
t815 = -mrSges(4,1) * t920 + mrSges(4,3) * t898 + Ifges(4,4) * t942 + Ifges(4,2) * t941 + Ifges(4,6) * t966 - pkin(3) * t831 - t962 * t935 + t978 * t937 - t1027;
t953 = Ifges(3,6) * t983 + (Ifges(3,4) * t994 + Ifges(3,2) * t999) * t1017;
t954 = Ifges(3,5) * t983 + (Ifges(3,1) * t994 + Ifges(3,4) * t999) * t1017;
t799 = Ifges(3,5) * t973 - Ifges(3,6) * t974 + Ifges(3,3) * t982 + mrSges(3,1) * t943 - mrSges(3,2) * t944 + t993 * t804 + t998 * t815 + pkin(2) * t1004 + pkin(9) * t1011 + (t953 * t994 - t954 * t999) * t1017;
t952 = Ifges(3,3) * t983 + (Ifges(3,5) * t994 + Ifges(3,6) * t999) * t1017;
t801 = mrSges(3,2) * t955 - mrSges(3,3) * t943 + Ifges(3,1) * t973 - Ifges(3,4) * t974 + Ifges(3,5) * t982 - pkin(9) * t824 + t1013 * t952 + t804 * t998 - t815 * t993 - t953 * t983;
t1003 = mrSges(4,1) * t897 - mrSges(4,2) * t898 + Ifges(4,5) * t942 + Ifges(4,6) * t941 + Ifges(4,3) * t966 + pkin(3) * t854 + pkin(10) * t832 + t997 * t816 + t992 * t817 + t962 * t936 - t961 * t937;
t803 = -mrSges(3,1) * t955 + mrSges(3,3) * t944 + Ifges(3,4) * t973 - Ifges(3,2) * t974 + Ifges(3,6) * t982 - pkin(2) * t824 - t1014 * t952 + t983 * t954 - t1003;
t1006 = mrSges(2,1) * t979 - mrSges(2,2) * t980 + Ifges(2,3) * qJDD(1) + pkin(1) * t809 + t803 * t1022 + t801 * t1023 + t814 * t1026 + t990 * t799;
t797 = -mrSges(2,2) * g(3) - mrSges(2,3) * t979 + Ifges(2,5) * qJDD(1) - t1001 * Ifges(2,6) + t999 * t801 - t994 * t803 + (-t808 * t988 - t809 * t990) * pkin(8);
t796 = mrSges(2,1) * g(3) + mrSges(2,3) * t980 + t1001 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t808 - t988 * t799 + (pkin(8) * t814 + t801 * t994 + t803 * t999) * t990;
t1 = [-m(1) * g(1) + t1012; -m(1) * g(2) + t1019; (-m(1) - m(2)) * g(3) + t808; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t1019 + t1000 * t797 - t995 * t796; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t1012 + t1000 * t796 + t995 * t797; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t1006; t1006; t799; t1003; t1027; t858; -t1005;];
tauJB  = t1;
