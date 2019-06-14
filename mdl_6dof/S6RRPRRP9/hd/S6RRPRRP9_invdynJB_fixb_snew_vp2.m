% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP9_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:29:27
% EndTime: 2019-05-06 18:29:59
% DurationCPUTime: 30.51s
% Computational Cost: add. (496945->378), mult. (1126828->476), div. (0->0), fcn. (904152->12), ass. (0->157)
t1012 = Ifges(6,1) + Ifges(7,1);
t1005 = Ifges(6,4) + Ifges(7,4);
t1004 = Ifges(6,5) + Ifges(7,5);
t1011 = Ifges(6,2) + Ifges(7,2);
t1003 = Ifges(6,6) + Ifges(7,6);
t1010 = Ifges(6,3) + Ifges(7,3);
t967 = cos(pkin(6));
t960 = qJD(1) * t967 + qJD(2);
t964 = sin(pkin(11));
t966 = cos(pkin(11));
t970 = sin(qJ(2));
t965 = sin(pkin(6));
t991 = qJD(1) * t965;
t986 = t970 * t991;
t939 = t960 * t966 - t964 * t986;
t940 = t960 * t964 + t966 * t986;
t969 = sin(qJ(4));
t973 = cos(qJ(4));
t921 = t939 * t973 - t940 * t969;
t974 = cos(qJ(2));
t990 = qJD(1) * t974;
t950 = (qJD(2) * t990 + qJDD(1) * t970) * t965;
t959 = qJDD(1) * t967 + qJDD(2);
t927 = -t950 * t964 + t959 * t966;
t928 = t950 * t966 + t959 * t964;
t890 = qJD(4) * t921 + t927 * t969 + t928 * t973;
t922 = t939 * t969 + t940 * t973;
t985 = t965 * t990;
t954 = qJD(4) - t985;
t968 = sin(qJ(5));
t972 = cos(qJ(5));
t909 = -t922 * t968 + t954 * t972;
t989 = qJDD(1) * t965;
t951 = -qJD(2) * t986 + t974 * t989;
t943 = qJDD(4) - t951;
t871 = qJD(5) * t909 + t890 * t972 + t943 * t968;
t910 = t922 * t972 + t954 * t968;
t884 = -mrSges(7,1) * t909 + mrSges(7,2) * t910;
t948 = (-pkin(2) * t974 - qJ(3) * t970) * t991;
t958 = t960 ^ 2;
t1008 = pkin(8) * t965;
t971 = sin(qJ(1));
t975 = cos(qJ(1));
t955 = t971 * g(1) - g(2) * t975;
t976 = qJD(1) ^ 2;
t946 = qJDD(1) * pkin(1) + t976 * t1008 + t955;
t956 = -g(1) * t975 - g(2) * t971;
t947 = -pkin(1) * t976 + pkin(8) * t989 + t956;
t999 = t967 * t970;
t992 = t946 * t999 + t974 * t947;
t906 = -t958 * pkin(2) + t959 * qJ(3) + (-g(3) * t970 + t948 * t990) * t965 + t992;
t1007 = t967 * g(3);
t907 = -t951 * pkin(2) - t1007 - t950 * qJ(3) + (-t946 + (pkin(2) * t970 - qJ(3) * t974) * t960 * qJD(1)) * t965;
t865 = -0.2e1 * qJD(3) * t940 - t964 * t906 + t966 * t907;
t861 = (-t939 * t985 - t928) * pkin(9) + (t939 * t940 - t951) * pkin(3) + t865;
t866 = 0.2e1 * qJD(3) * t939 + t966 * t906 + t964 * t907;
t929 = -pkin(3) * t985 - pkin(9) * t940;
t938 = t939 ^ 2;
t864 = -pkin(3) * t938 + pkin(9) * t927 + t929 * t985 + t866;
t856 = t969 * t861 + t973 * t864;
t901 = -pkin(4) * t921 - pkin(10) * t922;
t953 = t954 ^ 2;
t854 = -pkin(4) * t953 + pkin(10) * t943 + t901 * t921 + t856;
t1000 = t965 * t974;
t998 = t967 * t974;
t918 = -g(3) * t1000 + t946 * t998 - t970 * t947;
t905 = -t959 * pkin(2) - t958 * qJ(3) + t948 * t986 + qJDD(3) - t918;
t872 = -t927 * pkin(3) - t938 * pkin(9) + t940 * t929 + t905;
t889 = -qJD(4) * t922 + t927 * t973 - t928 * t969;
t859 = (-t921 * t954 - t890) * pkin(10) + (t922 * t954 - t889) * pkin(4) + t872;
t849 = -t968 * t854 + t972 * t859;
t888 = qJDD(5) - t889;
t920 = qJD(5) - t921;
t846 = -0.2e1 * qJD(6) * t910 + (t909 * t920 - t871) * qJ(6) + (t909 * t910 + t888) * pkin(5) + t849;
t891 = -mrSges(7,2) * t920 + mrSges(7,3) * t909;
t988 = m(7) * t846 + t888 * mrSges(7,1) + t920 * t891;
t843 = -t871 * mrSges(7,3) - t910 * t884 + t988;
t850 = t972 * t854 + t968 * t859;
t870 = -qJD(5) * t910 - t890 * t968 + t943 * t972;
t893 = pkin(5) * t920 - qJ(6) * t910;
t908 = t909 ^ 2;
t848 = -pkin(5) * t908 + qJ(6) * t870 + 0.2e1 * qJD(6) * t909 - t893 * t920 + t850;
t994 = t1004 * t920 + t1005 * t909 + t1012 * t910;
t995 = -t1003 * t920 - t1005 * t910 - t1011 * t909;
t1009 = mrSges(6,1) * t849 + mrSges(7,1) * t846 - mrSges(6,2) * t850 - mrSges(7,2) * t848 + pkin(5) * t843 + t1003 * t870 + t1004 * t871 + t1010 * t888 - t994 * t909 - t995 * t910;
t1006 = -mrSges(6,2) - mrSges(7,2);
t1001 = t965 * t970;
t919 = -g(3) * t1001 + t992;
t944 = mrSges(3,1) * t960 - mrSges(3,3) * t986;
t949 = (-mrSges(3,1) * t974 + mrSges(3,2) * t970) * t991;
t885 = -mrSges(6,1) * t909 + mrSges(6,2) * t910;
t892 = -mrSges(6,2) * t920 + mrSges(6,3) * t909;
t836 = m(6) * t849 + t888 * mrSges(6,1) + t920 * t892 + (-t884 - t885) * t910 + (-mrSges(6,3) - mrSges(7,3)) * t871 + t988;
t987 = m(7) * t848 + t870 * mrSges(7,3) + t909 * t884;
t894 = mrSges(7,1) * t920 - mrSges(7,3) * t910;
t993 = -mrSges(6,1) * t920 + mrSges(6,3) * t910 - t894;
t841 = m(6) * t850 + t870 * mrSges(6,3) + t1006 * t888 + t909 * t885 + t993 * t920 + t987;
t834 = -t836 * t968 + t972 * t841;
t900 = -mrSges(5,1) * t921 + mrSges(5,2) * t922;
t912 = mrSges(5,1) * t954 - mrSges(5,3) * t922;
t830 = m(5) * t856 - mrSges(5,2) * t943 + mrSges(5,3) * t889 + t900 * t921 - t912 * t954 + t834;
t855 = t861 * t973 - t969 * t864;
t853 = -pkin(4) * t943 - pkin(10) * t953 + t922 * t901 - t855;
t851 = -pkin(5) * t870 - qJ(6) * t908 + t893 * t910 + qJDD(6) + t853;
t981 = -m(7) * t851 + t870 * mrSges(7,1) + t909 * t891;
t842 = -m(6) * t853 + t870 * mrSges(6,1) + t1006 * t871 + t909 * t892 + t993 * t910 + t981;
t911 = -mrSges(5,2) * t954 + mrSges(5,3) * t921;
t838 = m(5) * t855 + t943 * mrSges(5,1) - t890 * mrSges(5,3) - t922 * t900 + t954 * t911 + t842;
t823 = t969 * t830 + t973 * t838;
t923 = -mrSges(4,1) * t939 + mrSges(4,2) * t940;
t925 = mrSges(4,2) * t985 + mrSges(4,3) * t939;
t821 = m(4) * t865 - mrSges(4,1) * t951 - mrSges(4,3) * t928 - t923 * t940 - t925 * t985 + t823;
t926 = -mrSges(4,1) * t985 - mrSges(4,3) * t940;
t982 = t973 * t830 - t838 * t969;
t822 = m(4) * t866 + mrSges(4,2) * t951 + mrSges(4,3) * t927 + t923 * t939 + t926 * t985 + t982;
t983 = -t821 * t964 + t966 * t822;
t813 = m(3) * t919 - mrSges(3,2) * t959 + mrSges(3,3) * t951 - t944 * t960 + t949 * t985 + t983;
t816 = t966 * t821 + t964 * t822;
t933 = -t965 * t946 - t1007;
t945 = -mrSges(3,2) * t960 + mrSges(3,3) * t985;
t815 = m(3) * t933 - t951 * mrSges(3,1) + t950 * mrSges(3,2) + (t944 * t970 - t945 * t974) * t991 + t816;
t833 = t972 * t836 + t968 * t841;
t979 = m(5) * t872 - t889 * mrSges(5,1) + mrSges(5,2) * t890 - t921 * t911 + t912 * t922 + t833;
t831 = m(4) * t905 - t927 * mrSges(4,1) + mrSges(4,2) * t928 - t939 * t925 + t926 * t940 + t979;
t827 = m(3) * t918 + mrSges(3,1) * t959 - mrSges(3,3) * t950 + t945 * t960 - t949 * t986 - t831;
t802 = t813 * t999 - t815 * t965 + t827 * t998;
t799 = m(2) * t955 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t976 + t802;
t808 = t974 * t813 - t827 * t970;
t806 = m(2) * t956 - mrSges(2,1) * t976 - qJDD(1) * mrSges(2,2) + t808;
t997 = t975 * t799 + t971 * t806;
t996 = -t1003 * t909 - t1004 * t910 - t1010 * t920;
t801 = t827 * t1000 + t813 * t1001 + t967 * t815;
t984 = -t799 * t971 + t975 * t806;
t844 = t871 * mrSges(7,2) + t910 * t894 - t981;
t824 = -mrSges(6,1) * t853 + mrSges(6,3) * t850 - mrSges(7,1) * t851 + mrSges(7,3) * t848 - pkin(5) * t844 + qJ(6) * t987 + (-qJ(6) * t894 + t994) * t920 + t996 * t910 + (-mrSges(7,2) * qJ(6) + t1003) * t888 + t1005 * t871 + t1011 * t870;
t832 = mrSges(6,2) * t853 + mrSges(7,2) * t851 - mrSges(6,3) * t849 - mrSges(7,3) * t846 - qJ(6) * t843 + t1004 * t888 + t1005 * t870 + t1012 * t871 - t996 * t909 + t995 * t920;
t896 = Ifges(5,5) * t922 + Ifges(5,6) * t921 + Ifges(5,3) * t954;
t897 = Ifges(5,4) * t922 + Ifges(5,2) * t921 + Ifges(5,6) * t954;
t809 = mrSges(5,2) * t872 - mrSges(5,3) * t855 + Ifges(5,1) * t890 + Ifges(5,4) * t889 + Ifges(5,5) * t943 - pkin(10) * t833 - t824 * t968 + t832 * t972 + t896 * t921 - t897 * t954;
t898 = Ifges(5,1) * t922 + Ifges(5,4) * t921 + Ifges(5,5) * t954;
t817 = -mrSges(5,1) * t872 + mrSges(5,3) * t856 + Ifges(5,4) * t890 + Ifges(5,2) * t889 + Ifges(5,6) * t943 - pkin(4) * t833 - t922 * t896 + t954 * t898 - t1009;
t913 = Ifges(4,5) * t940 + Ifges(4,6) * t939 - Ifges(4,3) * t985;
t915 = Ifges(4,1) * t940 + Ifges(4,4) * t939 - Ifges(4,5) * t985;
t797 = -mrSges(4,1) * t905 + mrSges(4,3) * t866 + Ifges(4,4) * t928 + Ifges(4,2) * t927 - Ifges(4,6) * t951 - pkin(3) * t979 + pkin(9) * t982 + t969 * t809 + t973 * t817 - t940 * t913 - t915 * t985;
t914 = Ifges(4,4) * t940 + Ifges(4,2) * t939 - Ifges(4,6) * t985;
t803 = mrSges(4,2) * t905 - mrSges(4,3) * t865 + Ifges(4,1) * t928 + Ifges(4,4) * t927 - Ifges(4,5) * t951 - pkin(9) * t823 + t809 * t973 - t817 * t969 + t913 * t939 + t914 * t985;
t931 = Ifges(3,6) * t960 + (Ifges(3,4) * t970 + Ifges(3,2) * t974) * t991;
t932 = Ifges(3,5) * t960 + (Ifges(3,1) * t970 + Ifges(3,4) * t974) * t991;
t792 = Ifges(3,5) * t950 + Ifges(3,6) * t951 + Ifges(3,3) * t959 + mrSges(3,1) * t918 - mrSges(3,2) * t919 + t964 * t803 + t966 * t797 - pkin(2) * t831 + qJ(3) * t983 + (t931 * t970 - t932 * t974) * t991;
t930 = Ifges(3,3) * t960 + (Ifges(3,5) * t970 + Ifges(3,6) * t974) * t991;
t794 = mrSges(3,2) * t933 - mrSges(3,3) * t918 + Ifges(3,1) * t950 + Ifges(3,4) * t951 + Ifges(3,5) * t959 - qJ(3) * t816 - t797 * t964 + t803 * t966 + t930 * t985 - t931 * t960;
t977 = mrSges(5,1) * t855 - mrSges(5,2) * t856 + Ifges(5,5) * t890 + Ifges(5,6) * t889 + Ifges(5,3) * t943 + pkin(4) * t842 + pkin(10) * t834 + t972 * t824 + t968 * t832 + t922 * t897 - t921 * t898;
t796 = -t977 + mrSges(3,3) * t919 - pkin(2) * t816 + Ifges(3,4) * t950 - pkin(3) * t823 + Ifges(3,6) * t959 + t960 * t932 - Ifges(4,6) * t927 - Ifges(4,5) * t928 - mrSges(3,1) * t933 - mrSges(4,1) * t865 + mrSges(4,2) * t866 + t939 * t915 - t940 * t914 + (Ifges(3,2) + Ifges(4,3)) * t951 - t930 * t986;
t980 = mrSges(2,1) * t955 - mrSges(2,2) * t956 + Ifges(2,3) * qJDD(1) + pkin(1) * t802 + t796 * t1000 + t794 * t1001 + t808 * t1008 + t967 * t792;
t790 = -mrSges(2,2) * g(3) - mrSges(2,3) * t955 + Ifges(2,5) * qJDD(1) - t976 * Ifges(2,6) + t974 * t794 - t970 * t796 + (-t801 * t965 - t802 * t967) * pkin(8);
t789 = mrSges(2,1) * g(3) + mrSges(2,3) * t956 + t976 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t801 - t965 * t792 + (pkin(8) * t808 + t794 * t970 + t796 * t974) * t967;
t1 = [-m(1) * g(1) + t984; -m(1) * g(2) + t997; (-m(1) - m(2)) * g(3) + t801; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t997 - t971 * t789 + t975 * t790; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t984 + t975 * t789 + t971 * t790; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t980; t980; t792; t831; t977; t1009; t844;];
tauJB  = t1;
