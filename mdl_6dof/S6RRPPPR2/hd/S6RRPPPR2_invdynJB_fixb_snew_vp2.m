% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-05-06 08:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPPR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:22:57
% EndTime: 2019-05-06 08:23:08
% DurationCPUTime: 10.56s
% Computational Cost: add. (150633->367), mult. (358875->447), div. (0->0), fcn. (243964->10), ass. (0->147)
t1011 = -2 * qJD(4);
t1010 = -Ifges(5,1) - Ifges(4,3);
t1004 = Ifges(4,5) - Ifges(5,4);
t1009 = Ifges(4,2) + Ifges(5,3);
t1008 = Ifges(5,2) + Ifges(4,1);
t1003 = Ifges(4,6) - Ifges(5,5);
t1002 = -Ifges(5,6) - Ifges(4,4);
t1001 = cos(pkin(9));
t960 = sin(pkin(9));
t966 = cos(qJ(2));
t992 = qJD(1) * t966;
t963 = sin(qJ(2));
t993 = qJD(1) * t963;
t932 = -t1001 * t992 + t960 * t993;
t933 = (t1001 * t963 + t960 * t966) * qJD(1);
t900 = t932 * pkin(3) - t933 * qJ(4);
t968 = qJD(2) ^ 2;
t964 = sin(qJ(1));
t967 = cos(qJ(1));
t951 = -t967 * g(1) - t964 * g(2);
t969 = qJD(1) ^ 2;
t939 = -t969 * pkin(1) + qJDD(1) * pkin(7) + t951;
t1000 = t963 * t939;
t1005 = pkin(2) * t969;
t988 = qJD(1) * qJD(2);
t944 = t963 * qJDD(1) + t966 * t988;
t878 = qJDD(2) * pkin(2) - t944 * qJ(3) - t1000 + (qJ(3) * t988 + t963 * t1005 - g(3)) * t966;
t919 = -t963 * g(3) + t966 * t939;
t945 = t966 * qJDD(1) - t963 * t988;
t947 = qJD(2) * pkin(2) - qJ(3) * t993;
t958 = t966 ^ 2;
t880 = t945 * qJ(3) - qJD(2) * t947 - t958 * t1005 + t919;
t997 = t1001 * t880 + t960 * t878;
t1007 = t968 * pkin(3) - qJDD(2) * qJ(4) + qJD(2) * t1011 + t932 * t900 - t997;
t864 = -0.2e1 * qJD(3) * t933 + t1001 * t878 - t960 * t880;
t901 = t932 * mrSges(4,1) + t933 * mrSges(4,2);
t910 = t1001 * t944 + t960 * t945;
t920 = -qJD(2) * mrSges(4,2) - t932 * mrSges(4,3);
t923 = t932 * mrSges(5,1) - qJD(2) * mrSges(5,3);
t852 = -qJDD(2) * pkin(3) - t968 * qJ(4) + t933 * t900 + qJDD(4) - t864;
t991 = qJD(2) * t932;
t845 = (t932 * t933 - qJDD(2)) * qJ(5) + (t910 + t991) * pkin(4) + t852;
t909 = -t1001 * t945 + t960 * t944;
t922 = t933 * pkin(4) - qJD(2) * qJ(5);
t931 = t932 ^ 2;
t950 = t964 * g(1) - t967 * g(2);
t979 = -qJDD(1) * pkin(1) - t950;
t884 = -t945 * pkin(2) + qJDD(3) + t947 * t993 + (-qJ(3) * t958 - pkin(7)) * t969 + t979;
t972 = (-t910 + t991) * qJ(4) + t884 + (pkin(3) * qJD(2) + t1011) * t933;
t849 = -t931 * pkin(4) - t933 * t922 + (pkin(3) + qJ(5)) * t909 + t972;
t959 = sin(pkin(10));
t961 = cos(pkin(10));
t915 = t961 * qJD(2) + t959 * t932;
t839 = -0.2e1 * qJD(5) * t915 + t961 * t845 - t959 * t849;
t890 = t961 * qJDD(2) + t959 * t909;
t914 = -t959 * qJD(2) + t961 * t932;
t837 = (t914 * t933 - t890) * pkin(8) + (t914 * t915 + t910) * pkin(5) + t839;
t840 = 0.2e1 * qJD(5) * t914 + t959 * t845 + t961 * t849;
t887 = t933 * pkin(5) - t915 * pkin(8);
t889 = -t959 * qJDD(2) + t961 * t909;
t913 = t914 ^ 2;
t838 = -t913 * pkin(5) + t889 * pkin(8) - t933 * t887 + t840;
t962 = sin(qJ(6));
t965 = cos(qJ(6));
t835 = t965 * t837 - t962 * t838;
t874 = t965 * t914 - t962 * t915;
t859 = t874 * qJD(6) + t962 * t889 + t965 * t890;
t875 = t962 * t914 + t965 * t915;
t866 = -t874 * mrSges(7,1) + t875 * mrSges(7,2);
t930 = qJD(6) + t933;
t867 = -t930 * mrSges(7,2) + t874 * mrSges(7,3);
t908 = qJDD(6) + t910;
t830 = m(7) * t835 + t908 * mrSges(7,1) - t859 * mrSges(7,3) - t875 * t866 + t930 * t867;
t836 = t962 * t837 + t965 * t838;
t858 = -t875 * qJD(6) + t965 * t889 - t962 * t890;
t868 = t930 * mrSges(7,1) - t875 * mrSges(7,3);
t831 = m(7) * t836 - t908 * mrSges(7,2) + t858 * mrSges(7,3) + t874 * t866 - t930 * t868;
t821 = t965 * t830 + t962 * t831;
t879 = -t914 * mrSges(6,1) + t915 * mrSges(6,2);
t885 = -t933 * mrSges(6,2) + t914 * mrSges(6,3);
t819 = m(6) * t839 + t910 * mrSges(6,1) - t890 * mrSges(6,3) - t915 * t879 + t933 * t885 + t821;
t886 = t933 * mrSges(6,1) - t915 * mrSges(6,3);
t982 = -t962 * t830 + t965 * t831;
t820 = m(6) * t840 - t910 * mrSges(6,2) + t889 * mrSges(6,3) + t914 * t879 - t933 * t886 + t982;
t816 = t961 * t819 + t959 * t820;
t902 = -t932 * mrSges(5,2) - t933 * mrSges(5,3);
t975 = -m(5) * t852 - t910 * mrSges(5,1) - t933 * t902 - t816;
t812 = m(4) * t864 - t910 * mrSges(4,3) - t933 * t901 + (mrSges(4,1) - mrSges(5,2)) * qJDD(2) + (t920 - t923) * qJD(2) + t975;
t990 = qJD(3) * t932;
t927 = -0.2e1 * t990;
t865 = t927 + t997;
t921 = qJD(2) * mrSges(4,1) - t933 * mrSges(4,3);
t847 = -t909 * pkin(4) - t931 * qJ(5) + qJD(2) * t922 + qJDD(5) - t1007 + t927;
t842 = -t889 * pkin(5) - t913 * pkin(8) + t915 * t887 + t847;
t977 = m(7) * t842 - t858 * mrSges(7,1) + t859 * mrSges(7,2) - t874 * t867 + t875 * t868;
t833 = m(6) * t847 - t889 * mrSges(6,1) + t890 * mrSges(6,2) - t914 * t885 + t915 * t886 + t977;
t850 = 0.2e1 * t990 + t1007;
t924 = t933 * mrSges(5,1) + qJD(2) * mrSges(5,2);
t973 = -m(5) * t850 + qJDD(2) * mrSges(5,3) + qJD(2) * t924 + t833;
t826 = t973 + (-mrSges(4,3) - mrSges(5,1)) * t909 + (-t901 - t902) * t932 - qJD(2) * t921 - qJDD(2) * mrSges(4,2) + m(4) * t865;
t805 = t1001 * t812 + t960 * t826;
t860 = Ifges(7,5) * t875 + Ifges(7,6) * t874 + Ifges(7,3) * t930;
t862 = Ifges(7,1) * t875 + Ifges(7,4) * t874 + Ifges(7,5) * t930;
t822 = -mrSges(7,1) * t842 + mrSges(7,3) * t836 + Ifges(7,4) * t859 + Ifges(7,2) * t858 + Ifges(7,6) * t908 - t875 * t860 + t930 * t862;
t861 = Ifges(7,4) * t875 + Ifges(7,2) * t874 + Ifges(7,6) * t930;
t823 = mrSges(7,2) * t842 - mrSges(7,3) * t835 + Ifges(7,1) * t859 + Ifges(7,4) * t858 + Ifges(7,5) * t908 + t874 * t860 - t930 * t861;
t869 = Ifges(6,5) * t915 + Ifges(6,6) * t914 + Ifges(6,3) * t933;
t871 = Ifges(6,1) * t915 + Ifges(6,4) * t914 + Ifges(6,5) * t933;
t806 = -mrSges(6,1) * t847 + mrSges(6,3) * t840 + Ifges(6,4) * t890 + Ifges(6,2) * t889 + Ifges(6,6) * t910 - pkin(5) * t977 + pkin(8) * t982 + t965 * t822 + t962 * t823 - t915 * t869 + t933 * t871;
t870 = Ifges(6,4) * t915 + Ifges(6,2) * t914 + Ifges(6,6) * t933;
t807 = mrSges(6,2) * t847 - mrSges(6,3) * t839 + Ifges(6,1) * t890 + Ifges(6,4) * t889 + Ifges(6,5) * t910 - pkin(8) * t821 - t962 * t822 + t965 * t823 + t914 * t869 - t933 * t870;
t815 = qJDD(2) * mrSges(5,2) + qJD(2) * t923 - t975;
t918 = -t966 * g(3) - t1000;
t935 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t963 + Ifges(3,2) * t966) * qJD(1);
t936 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t963 + Ifges(3,4) * t966) * qJD(1);
t994 = t1004 * qJD(2) + t1002 * t932 + t1008 * t933;
t996 = -t1003 * qJD(2) + t1002 * t933 + t1009 * t932;
t1006 = (t963 * t935 - t966 * t936) * qJD(1) + (Ifges(3,3) - t1010) * qJDD(2) - t1003 * t909 + t1004 * t910 + t994 * t932 - t996 * t933 + mrSges(3,1) * t918 + mrSges(4,1) * t864 - mrSges(3,2) * t919 - mrSges(4,2) * t865 + mrSges(5,2) * t852 - mrSges(5,3) * t850 + Ifges(3,5) * t944 + Ifges(3,6) * t945 + pkin(2) * t805 - pkin(3) * t815 + qJ(4) * (-t909 * mrSges(5,1) - t932 * t902 + t973) - qJ(5) * t816 - t959 * t806 + t961 * t807;
t943 = (-mrSges(3,1) * t966 + mrSges(3,2) * t963) * qJD(1);
t949 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t992;
t803 = m(3) * t918 + qJDD(2) * mrSges(3,1) - t944 * mrSges(3,3) + qJD(2) * t949 - t943 * t993 + t805;
t948 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t993;
t983 = t1001 * t826 - t960 * t812;
t804 = m(3) * t919 - qJDD(2) * mrSges(3,2) + t945 * mrSges(3,3) - qJD(2) * t948 + t943 * t992 + t983;
t984 = -t963 * t803 + t966 * t804;
t797 = m(2) * t951 - t969 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t984;
t854 = t909 * pkin(3) + t972;
t998 = -t959 * t819 + t961 * t820;
t814 = m(5) * t854 - t909 * mrSges(5,2) - t910 * mrSges(5,3) - t932 * t923 - t933 * t924 + t998;
t813 = m(4) * t884 + t909 * mrSges(4,1) + t910 * mrSges(4,2) + t932 * t920 + t933 * t921 + t814;
t938 = -t969 * pkin(7) + t979;
t971 = -m(3) * t938 + t945 * mrSges(3,1) - t944 * mrSges(3,2) - t948 * t993 + t949 * t992 - t813;
t809 = m(2) * t950 + qJDD(1) * mrSges(2,1) - t969 * mrSges(2,2) + t971;
t999 = t964 * t797 + t967 * t809;
t799 = t966 * t803 + t963 * t804;
t995 = t1010 * qJD(2) + t1003 * t932 - t1004 * t933;
t985 = t967 * t797 - t964 * t809;
t793 = -mrSges(4,1) * t884 - mrSges(5,1) * t850 + mrSges(5,2) * t854 + mrSges(4,3) * t865 - pkin(3) * t814 + pkin(4) * t833 - qJ(5) * t998 + t994 * qJD(2) + t1003 * qJDD(2) - t1002 * t910 - t1009 * t909 - t961 * t806 - t959 * t807 + t995 * t933;
t974 = mrSges(7,1) * t835 - mrSges(7,2) * t836 + Ifges(7,5) * t859 + Ifges(7,6) * t858 + Ifges(7,3) * t908 + t875 * t861 - t874 * t862;
t794 = -t914 * t871 + t915 * t870 + mrSges(4,2) * t884 + Ifges(6,6) * t889 + Ifges(6,5) * t890 - mrSges(4,3) * t864 + mrSges(5,1) * t852 - mrSges(5,3) * t854 - mrSges(6,2) * t840 + mrSges(6,1) * t839 + t974 + pkin(5) * t821 + pkin(4) * t816 - qJ(4) * t814 + t995 * t932 + (Ifges(6,3) + t1008) * t910 + t1002 * t909 + t996 * qJD(2) + t1004 * qJDD(2);
t934 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t963 + Ifges(3,6) * t966) * qJD(1);
t789 = -mrSges(3,1) * t938 + mrSges(3,3) * t919 + Ifges(3,4) * t944 + Ifges(3,2) * t945 + Ifges(3,6) * qJDD(2) - pkin(2) * t813 + qJ(3) * t983 + qJD(2) * t936 + t1001 * t793 + t960 * t794 - t934 * t993;
t791 = mrSges(3,2) * t938 - mrSges(3,3) * t918 + Ifges(3,1) * t944 + Ifges(3,4) * t945 + Ifges(3,5) * qJDD(2) - qJ(3) * t805 - qJD(2) * t935 + t1001 * t794 - t960 * t793 + t934 * t992;
t976 = mrSges(2,1) * t950 - mrSges(2,2) * t951 + Ifges(2,3) * qJDD(1) + pkin(1) * t971 + pkin(7) * t984 + t966 * t789 + t963 * t791;
t792 = mrSges(2,1) * g(3) + mrSges(2,3) * t951 + t969 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t799 - t1006;
t787 = -mrSges(2,2) * g(3) - mrSges(2,3) * t950 + Ifges(2,5) * qJDD(1) - t969 * Ifges(2,6) - pkin(7) * t799 - t963 * t789 + t966 * t791;
t1 = [-m(1) * g(1) + t985; -m(1) * g(2) + t999; (-m(1) - m(2)) * g(3) + t799; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t999 + t967 * t787 - t964 * t792; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t985 + t964 * t787 + t967 * t792; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t976; t976; t1006; t813; t815; t833; t974;];
tauJB  = t1;
