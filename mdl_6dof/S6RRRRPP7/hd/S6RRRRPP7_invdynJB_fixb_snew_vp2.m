% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRRPP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-05-07 18:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRRPP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:42:24
% EndTime: 2019-05-07 18:42:58
% DurationCPUTime: 32.47s
% Computational Cost: add. (563672->377), mult. (1210212->477), div. (0->0), fcn. (962704->12), ass. (0->158)
t1007 = Ifges(6,1) + Ifges(7,1);
t999 = Ifges(6,4) - Ifges(7,5);
t998 = Ifges(6,5) + Ifges(7,4);
t1006 = -Ifges(6,2) - Ifges(7,3);
t1005 = -Ifges(7,2) - Ifges(6,3);
t997 = Ifges(6,6) - Ifges(7,6);
t958 = sin(pkin(6));
t962 = sin(qJ(2));
t966 = cos(qJ(2));
t983 = qJD(1) * qJD(2);
t945 = (-qJDD(1) * t966 + t962 * t983) * t958;
t1000 = -mrSges(6,3) - mrSges(7,2);
t1003 = -2 * qJD(5);
t985 = qJD(1) * t958;
t943 = (-pkin(2) * t966 - pkin(9) * t962) * t985;
t959 = cos(pkin(6));
t953 = qJD(1) * t959 + qJD(2);
t951 = t953 ^ 2;
t952 = qJDD(1) * t959 + qJDD(2);
t984 = qJD(1) * t966;
t1002 = pkin(8) * t958;
t963 = sin(qJ(1));
t967 = cos(qJ(1));
t949 = t963 * g(1) - g(2) * t967;
t968 = qJD(1) ^ 2;
t940 = qJDD(1) * pkin(1) + t1002 * t968 + t949;
t950 = -g(1) * t967 - g(2) * t963;
t941 = -pkin(1) * t968 + qJDD(1) * t1002 + t950;
t993 = t959 * t962;
t986 = t940 * t993 + t966 * t941;
t895 = -t951 * pkin(2) + t952 * pkin(9) + (-g(3) * t962 + t943 * t984) * t958 + t986;
t1001 = t959 * g(3);
t944 = (qJDD(1) * t962 + t966 * t983) * t958;
t896 = t945 * pkin(2) - t944 * pkin(9) - t1001 + (-t940 + (pkin(2) * t962 - pkin(9) * t966) * t953 * qJD(1)) * t958;
t961 = sin(qJ(3));
t965 = cos(qJ(3));
t863 = t965 * t895 + t961 * t896;
t980 = t962 * t985;
t933 = t953 * t965 - t961 * t980;
t934 = t953 * t961 + t965 * t980;
t918 = -pkin(3) * t933 - pkin(10) * t934;
t937 = qJDD(3) + t945;
t979 = t958 * t984;
t948 = qJD(3) - t979;
t947 = t948 ^ 2;
t853 = -pkin(3) * t947 + pkin(10) * t937 + t918 * t933 + t863;
t992 = t959 * t966;
t994 = t958 * t966;
t915 = -g(3) * t994 + t940 * t992 - t962 * t941;
t894 = -t952 * pkin(2) - t951 * pkin(9) + t943 * t980 - t915;
t913 = -qJD(3) * t934 - t944 * t961 + t952 * t965;
t914 = qJD(3) * t933 + t944 * t965 + t952 * t961;
t856 = (-t933 * t948 - t914) * pkin(10) + (t934 * t948 - t913) * pkin(3) + t894;
t960 = sin(qJ(4));
t964 = cos(qJ(4));
t848 = -t960 * t853 + t964 * t856;
t921 = -t934 * t960 + t948 * t964;
t881 = qJD(4) * t921 + t914 * t964 + t937 * t960;
t911 = qJDD(4) - t913;
t922 = t934 * t964 + t948 * t960;
t932 = qJD(4) - t933;
t845 = (t921 * t932 - t881) * qJ(5) + (t921 * t922 + t911) * pkin(4) + t848;
t849 = t964 * t853 + t960 * t856;
t880 = -qJD(4) * t922 - t914 * t960 + t937 * t964;
t903 = pkin(4) * t932 - qJ(5) * t922;
t920 = t921 ^ 2;
t847 = -pkin(4) * t920 + qJ(5) * t880 - t903 * t932 + t849;
t957 = sin(pkin(11));
t996 = cos(pkin(11));
t898 = -t921 * t996 + t957 * t922;
t841 = t898 * t1003 + t957 * t845 + t996 * t847;
t860 = -t880 * t996 + t957 * t881;
t899 = t957 * t921 + t922 * t996;
t884 = mrSges(6,1) * t932 - mrSges(6,3) * t899;
t873 = pkin(5) * t898 - qJ(6) * t899;
t931 = t932 ^ 2;
t838 = -pkin(5) * t931 + qJ(6) * t911 + 0.2e1 * qJD(6) * t932 - t873 * t898 + t841;
t885 = -mrSges(7,1) * t932 + mrSges(7,2) * t899;
t981 = m(7) * t838 + t911 * mrSges(7,3) + t932 * t885;
t874 = mrSges(7,1) * t898 - mrSges(7,3) * t899;
t987 = -mrSges(6,1) * t898 - mrSges(6,2) * t899 - t874;
t827 = m(6) * t841 - t911 * mrSges(6,2) + t1000 * t860 - t932 * t884 + t898 * t987 + t981;
t974 = t845 * t996 - t957 * t847;
t840 = t1003 * t899 + t974;
t861 = t957 * t880 + t881 * t996;
t883 = -mrSges(6,2) * t932 - mrSges(6,3) * t898;
t839 = -t911 * pkin(5) - t931 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t873) * t899 - t974;
t882 = -mrSges(7,2) * t898 + mrSges(7,3) * t932;
t975 = -m(7) * t839 + t911 * mrSges(7,1) + t932 * t882;
t829 = m(6) * t840 + t911 * mrSges(6,1) + t1000 * t861 + t932 * t883 + t899 * t987 + t975;
t822 = t957 * t827 + t996 * t829;
t835 = t861 * mrSges(7,2) + t899 * t874 - t975;
t887 = Ifges(5,4) * t922 + Ifges(5,2) * t921 + Ifges(5,6) * t932;
t888 = Ifges(5,1) * t922 + Ifges(5,4) * t921 + Ifges(5,5) * t932;
t988 = t1007 * t899 - t999 * t898 + t998 * t932;
t989 = t1006 * t898 + t899 * t999 + t932 * t997;
t1004 = -t860 * t997 + t861 * t998 + t898 * t988 + t899 * t989 + (Ifges(5,3) - t1005) * t911 + mrSges(5,1) * t848 + mrSges(6,1) * t840 - mrSges(7,1) * t839 - mrSges(5,2) * t849 - mrSges(6,2) * t841 + mrSges(7,3) * t838 + Ifges(5,5) * t881 + Ifges(5,6) * t880 + pkin(4) * t822 - pkin(5) * t835 + qJ(6) * (-t860 * mrSges(7,2) - t898 * t874 + t981) + t922 * t887 - t921 * t888;
t995 = t958 * t962;
t916 = -g(3) * t995 + t986;
t938 = mrSges(3,1) * t953 - mrSges(3,3) * t980;
t942 = (-mrSges(3,1) * t966 + mrSges(3,2) * t962) * t985;
t900 = -mrSges(5,1) * t921 + mrSges(5,2) * t922;
t902 = -mrSges(5,2) * t932 + mrSges(5,3) * t921;
t820 = m(5) * t848 + mrSges(5,1) * t911 - mrSges(5,3) * t881 - t900 * t922 + t902 * t932 + t822;
t904 = mrSges(5,1) * t932 - mrSges(5,3) * t922;
t976 = t996 * t827 - t829 * t957;
t821 = m(5) * t849 - mrSges(5,2) * t911 + mrSges(5,3) * t880 + t900 * t921 - t904 * t932 + t976;
t818 = -t820 * t960 + t964 * t821;
t917 = -mrSges(4,1) * t933 + mrSges(4,2) * t934;
t924 = mrSges(4,1) * t948 - mrSges(4,3) * t934;
t816 = m(4) * t863 - mrSges(4,2) * t937 + mrSges(4,3) * t913 + t917 * t933 - t924 * t948 + t818;
t862 = -t961 * t895 + t965 * t896;
t852 = -t937 * pkin(3) - t947 * pkin(10) + t934 * t918 - t862;
t850 = -t880 * pkin(4) - t920 * qJ(5) + t922 * t903 + qJDD(5) + t852;
t843 = -0.2e1 * qJD(6) * t899 + (t898 * t932 - t861) * qJ(6) + (t899 * t932 + t860) * pkin(5) + t850;
t836 = m(7) * t843 + t860 * mrSges(7,1) - t861 * mrSges(7,3) + t898 * t882 - t899 * t885;
t833 = m(6) * t850 + t860 * mrSges(6,1) + mrSges(6,2) * t861 + t898 * t883 + t884 * t899 + t836;
t832 = -m(5) * t852 + t880 * mrSges(5,1) - mrSges(5,2) * t881 + t921 * t902 - t904 * t922 - t833;
t923 = -mrSges(4,2) * t948 + mrSges(4,3) * t933;
t831 = m(4) * t862 + mrSges(4,1) * t937 - mrSges(4,3) * t914 - t917 * t934 + t923 * t948 + t832;
t977 = t965 * t816 - t831 * t961;
t807 = m(3) * t916 - mrSges(3,2) * t952 - mrSges(3,3) * t945 - t938 * t953 + t942 * t979 + t977;
t810 = t961 * t816 + t965 * t831;
t928 = -t958 * t940 - t1001;
t939 = -mrSges(3,2) * t953 + mrSges(3,3) * t979;
t809 = m(3) * t928 + t945 * mrSges(3,1) + t944 * mrSges(3,2) + (t938 * t962 - t939 * t966) * t985 + t810;
t817 = t820 * t964 + t821 * t960;
t971 = -m(4) * t894 + t913 * mrSges(4,1) - mrSges(4,2) * t914 + t933 * t923 - t924 * t934 - t817;
t813 = m(3) * t915 + mrSges(3,1) * t952 - mrSges(3,3) * t944 + t939 * t953 - t942 * t980 + t971;
t795 = t807 * t993 - t809 * t958 + t813 * t992;
t792 = m(2) * t949 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t968 + t795;
t800 = t966 * t807 - t813 * t962;
t798 = m(2) * t950 - mrSges(2,1) * t968 - qJDD(1) * mrSges(2,2) + t800;
t991 = t967 * t792 + t963 * t798;
t990 = t1005 * t932 + t898 * t997 - t899 * t998;
t794 = t807 * t995 + t959 * t809 + t813 * t994;
t978 = -t792 * t963 + t967 * t798;
t823 = -mrSges(6,1) * t850 - mrSges(7,1) * t843 + mrSges(7,2) * t838 + mrSges(6,3) * t841 - pkin(5) * t836 + t1006 * t860 + t999 * t861 + t990 * t899 + t997 * t911 + t988 * t932;
t824 = mrSges(6,2) * t850 + mrSges(7,2) * t839 - mrSges(6,3) * t840 - mrSges(7,3) * t843 - qJ(6) * t836 + t1007 * t861 - t999 * t860 + t990 * t898 + t998 * t911 - t989 * t932;
t886 = Ifges(5,5) * t922 + Ifges(5,6) * t921 + Ifges(5,3) * t932;
t802 = -mrSges(5,1) * t852 + mrSges(5,3) * t849 + Ifges(5,4) * t881 + Ifges(5,2) * t880 + Ifges(5,6) * t911 - pkin(4) * t833 + qJ(5) * t976 + t823 * t996 + t957 * t824 - t922 * t886 + t932 * t888;
t803 = mrSges(5,2) * t852 - mrSges(5,3) * t848 + Ifges(5,1) * t881 + Ifges(5,4) * t880 + Ifges(5,5) * t911 - qJ(5) * t822 - t957 * t823 + t824 * t996 + t921 * t886 - t932 * t887;
t907 = Ifges(4,5) * t934 + Ifges(4,6) * t933 + Ifges(4,3) * t948;
t908 = Ifges(4,4) * t934 + Ifges(4,2) * t933 + Ifges(4,6) * t948;
t790 = mrSges(4,2) * t894 - mrSges(4,3) * t862 + Ifges(4,1) * t914 + Ifges(4,4) * t913 + Ifges(4,5) * t937 - pkin(10) * t817 - t802 * t960 + t803 * t964 + t907 * t933 - t908 * t948;
t909 = Ifges(4,1) * t934 + Ifges(4,4) * t933 + Ifges(4,5) * t948;
t801 = -mrSges(4,1) * t894 + mrSges(4,3) * t863 + Ifges(4,4) * t914 + Ifges(4,2) * t913 + Ifges(4,6) * t937 - pkin(3) * t817 - t934 * t907 + t948 * t909 - t1004;
t926 = Ifges(3,6) * t953 + (Ifges(3,4) * t962 + Ifges(3,2) * t966) * t985;
t927 = Ifges(3,5) * t953 + (Ifges(3,1) * t962 + Ifges(3,4) * t966) * t985;
t785 = Ifges(3,5) * t944 - Ifges(3,6) * t945 + Ifges(3,3) * t952 + mrSges(3,1) * t915 - mrSges(3,2) * t916 + t961 * t790 + t965 * t801 + pkin(2) * t971 + pkin(9) * t977 + (t926 * t962 - t927 * t966) * t985;
t925 = Ifges(3,3) * t953 + (Ifges(3,5) * t962 + Ifges(3,6) * t966) * t985;
t787 = mrSges(3,2) * t928 - mrSges(3,3) * t915 + Ifges(3,1) * t944 - Ifges(3,4) * t945 + Ifges(3,5) * t952 - pkin(9) * t810 + t790 * t965 - t801 * t961 + t925 * t979 - t926 * t953;
t970 = mrSges(4,1) * t862 - mrSges(4,2) * t863 + Ifges(4,5) * t914 + Ifges(4,6) * t913 + Ifges(4,3) * t937 + pkin(3) * t832 + pkin(10) * t818 + t964 * t802 + t960 * t803 + t934 * t908 - t933 * t909;
t789 = -mrSges(3,1) * t928 + mrSges(3,3) * t916 + Ifges(3,4) * t944 - Ifges(3,2) * t945 + Ifges(3,6) * t952 - pkin(2) * t810 - t925 * t980 + t953 * t927 - t970;
t972 = mrSges(2,1) * t949 - mrSges(2,2) * t950 + Ifges(2,3) * qJDD(1) + pkin(1) * t795 + t800 * t1002 + t959 * t785 + t787 * t995 + t789 * t994;
t783 = -mrSges(2,2) * g(3) - mrSges(2,3) * t949 + Ifges(2,5) * qJDD(1) - t968 * Ifges(2,6) + t966 * t787 - t962 * t789 + (-t794 * t958 - t795 * t959) * pkin(8);
t782 = mrSges(2,1) * g(3) + mrSges(2,3) * t950 + t968 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t794 - t958 * t785 + (pkin(8) * t800 + t787 * t962 + t789 * t966) * t959;
t1 = [-m(1) * g(1) + t978; -m(1) * g(2) + t991; (-m(1) - m(2)) * g(3) + t794; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t991 - t963 * t782 + t967 * t783; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t978 + t967 * t782 + t963 * t783; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t972; t972; t785; t970; t1004; t833; t835;];
tauJB  = t1;
