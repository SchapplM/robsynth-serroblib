% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 04:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:08:29
% EndTime: 2019-05-07 04:09:01
% DurationCPUTime: 32.31s
% Computational Cost: add. (544814->387), mult. (1242887->490), div. (0->0), fcn. (927001->12), ass. (0->152)
t994 = -2 * qJD(4);
t961 = sin(qJ(2));
t965 = cos(qJ(2));
t985 = qJD(1) * qJD(2);
t940 = qJDD(1) * t961 + t965 * t985;
t962 = sin(qJ(1));
t966 = cos(qJ(1));
t947 = -g(1) * t966 - g(2) * t962;
t967 = qJD(1) ^ 2;
t935 = -pkin(1) * t967 + qJDD(1) * pkin(7) + t947;
t990 = t961 * t935;
t992 = pkin(2) * t967;
t894 = qJDD(2) * pkin(2) - t940 * pkin(8) - t990 + (pkin(8) * t985 + t961 * t992 - g(3)) * t965;
t922 = -g(3) * t961 + t965 * t935;
t941 = qJDD(1) * t965 - t961 * t985;
t988 = qJD(1) * t961;
t945 = qJD(2) * pkin(2) - pkin(8) * t988;
t955 = t965 ^ 2;
t895 = pkin(8) * t941 - qJD(2) * t945 - t955 * t992 + t922;
t960 = sin(qJ(3));
t964 = cos(qJ(3));
t869 = t964 * t894 - t960 * t895;
t932 = (-t960 * t961 + t964 * t965) * qJD(1);
t906 = qJD(3) * t932 + t940 * t964 + t941 * t960;
t933 = (t960 * t965 + t961 * t964) * qJD(1);
t952 = qJDD(2) + qJDD(3);
t953 = qJD(2) + qJD(3);
t850 = (t932 * t953 - t906) * qJ(4) + (t932 * t933 + t952) * pkin(3) + t869;
t870 = t960 * t894 + t964 * t895;
t905 = -qJD(3) * t933 - t940 * t960 + t941 * t964;
t924 = pkin(3) * t953 - qJ(4) * t933;
t928 = t932 ^ 2;
t854 = -pkin(3) * t928 + qJ(4) * t905 - t924 * t953 + t870;
t957 = sin(pkin(10));
t991 = cos(pkin(10));
t919 = t957 * t932 + t933 * t991;
t837 = t850 * t991 - t957 * t854 + t919 * t994;
t918 = -t991 * t932 + t933 * t957;
t838 = t957 * t850 + t991 * t854 + t918 * t994;
t879 = -t991 * t905 + t906 * t957;
t889 = mrSges(5,1) * t918 + mrSges(5,2) * t919;
t909 = mrSges(5,1) * t953 - mrSges(5,3) * t919;
t888 = pkin(4) * t918 - qJ(5) * t919;
t951 = t953 ^ 2;
t835 = -pkin(4) * t951 + qJ(5) * t952 - t888 * t918 + t838;
t946 = t962 * g(1) - t966 * g(2);
t975 = -qJDD(1) * pkin(1) - t946;
t907 = -t941 * pkin(2) + t945 * t988 + (-pkin(8) * t955 - pkin(7)) * t967 + t975;
t856 = -t905 * pkin(3) - t928 * qJ(4) + t933 * t924 + qJDD(4) + t907;
t880 = t957 * t905 + t906 * t991;
t841 = (t918 * t953 - t880) * qJ(5) + (t919 * t953 + t879) * pkin(4) + t856;
t956 = sin(pkin(11));
t958 = cos(pkin(11));
t901 = t919 * t958 + t953 * t956;
t830 = -0.2e1 * qJD(5) * t901 - t956 * t835 + t958 * t841;
t868 = t880 * t958 + t952 * t956;
t900 = -t919 * t956 + t953 * t958;
t828 = (t900 * t918 - t868) * pkin(9) + (t900 * t901 + t879) * pkin(5) + t830;
t831 = 0.2e1 * qJD(5) * t900 + t958 * t835 + t956 * t841;
t867 = -t880 * t956 + t952 * t958;
t882 = pkin(5) * t918 - pkin(9) * t901;
t899 = t900 ^ 2;
t829 = -pkin(5) * t899 + pkin(9) * t867 - t882 * t918 + t831;
t959 = sin(qJ(6));
t963 = cos(qJ(6));
t826 = t828 * t963 - t829 * t959;
t871 = t900 * t963 - t901 * t959;
t844 = qJD(6) * t871 + t867 * t959 + t868 * t963;
t872 = t900 * t959 + t901 * t963;
t851 = -mrSges(7,1) * t871 + mrSges(7,2) * t872;
t912 = qJD(6) + t918;
t857 = -mrSges(7,2) * t912 + mrSges(7,3) * t871;
t877 = qJDD(6) + t879;
t822 = m(7) * t826 + mrSges(7,1) * t877 - t844 * mrSges(7,3) - t851 * t872 + t857 * t912;
t827 = t828 * t959 + t829 * t963;
t843 = -qJD(6) * t872 + t867 * t963 - t868 * t959;
t858 = mrSges(7,1) * t912 - mrSges(7,3) * t872;
t823 = m(7) * t827 - mrSges(7,2) * t877 + t843 * mrSges(7,3) + t851 * t871 - t858 * t912;
t814 = t963 * t822 + t959 * t823;
t878 = -mrSges(6,1) * t900 + mrSges(6,2) * t901;
t978 = -mrSges(6,2) * t918 + mrSges(6,3) * t900;
t812 = m(6) * t830 + t879 * mrSges(6,1) - t868 * mrSges(6,3) - t901 * t878 + t918 * t978 + t814;
t881 = mrSges(6,1) * t918 - mrSges(6,3) * t901;
t979 = -t822 * t959 + t963 * t823;
t813 = m(6) * t831 - mrSges(6,2) * t879 + mrSges(6,3) * t867 + t878 * t900 - t881 * t918 + t979;
t980 = -t812 * t956 + t958 * t813;
t804 = m(5) * t838 - mrSges(5,2) * t952 - mrSges(5,3) * t879 - t889 * t918 - t909 * t953 + t980;
t834 = -t952 * pkin(4) - t951 * qJ(5) + t919 * t888 + qJDD(5) - t837;
t832 = -t867 * pkin(5) - t899 * pkin(9) + t901 * t882 + t834;
t973 = m(7) * t832 - t843 * mrSges(7,1) + t844 * mrSges(7,2) - t871 * t857 + t858 * t872;
t825 = m(6) * t834 - t867 * mrSges(6,1) + mrSges(6,2) * t868 + t881 * t901 - t900 * t978 + t973;
t908 = -mrSges(5,2) * t953 - mrSges(5,3) * t918;
t818 = m(5) * t837 + mrSges(5,1) * t952 - mrSges(5,3) * t880 - t889 * t919 + t908 * t953 - t825;
t794 = t957 * t804 + t991 * t818;
t920 = -mrSges(4,1) * t932 + mrSges(4,2) * t933;
t923 = -mrSges(4,2) * t953 + mrSges(4,3) * t932;
t791 = m(4) * t869 + mrSges(4,1) * t952 - mrSges(4,3) * t906 - t920 * t933 + t923 * t953 + t794;
t925 = mrSges(4,1) * t953 - mrSges(4,3) * t933;
t981 = t991 * t804 - t818 * t957;
t792 = m(4) * t870 - mrSges(4,2) * t952 + mrSges(4,3) * t905 + t920 * t932 - t925 * t953 + t981;
t786 = t964 * t791 + t960 * t792;
t921 = -t965 * g(3) - t990;
t930 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t961 + Ifges(3,2) * t965) * qJD(1);
t931 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t961 + Ifges(3,4) * t965) * qJD(1);
t845 = Ifges(7,5) * t872 + Ifges(7,6) * t871 + Ifges(7,3) * t912;
t847 = Ifges(7,1) * t872 + Ifges(7,4) * t871 + Ifges(7,5) * t912;
t815 = -mrSges(7,1) * t832 + mrSges(7,3) * t827 + Ifges(7,4) * t844 + Ifges(7,2) * t843 + Ifges(7,6) * t877 - t845 * t872 + t847 * t912;
t846 = Ifges(7,4) * t872 + Ifges(7,2) * t871 + Ifges(7,6) * t912;
t816 = mrSges(7,2) * t832 - mrSges(7,3) * t826 + Ifges(7,1) * t844 + Ifges(7,4) * t843 + Ifges(7,5) * t877 + t845 * t871 - t846 * t912;
t859 = Ifges(6,5) * t901 + Ifges(6,6) * t900 + Ifges(6,3) * t918;
t861 = Ifges(6,1) * t901 + Ifges(6,4) * t900 + Ifges(6,5) * t918;
t796 = -mrSges(6,1) * t834 + mrSges(6,3) * t831 + Ifges(6,4) * t868 + Ifges(6,2) * t867 + Ifges(6,6) * t879 - pkin(5) * t973 + pkin(9) * t979 + t963 * t815 + t959 * t816 - t901 * t859 + t918 * t861;
t860 = Ifges(6,4) * t901 + Ifges(6,2) * t900 + Ifges(6,6) * t918;
t798 = mrSges(6,2) * t834 - mrSges(6,3) * t830 + Ifges(6,1) * t868 + Ifges(6,4) * t867 + Ifges(6,5) * t879 - pkin(9) * t814 - t815 * t959 + t816 * t963 + t859 * t900 - t860 * t918;
t884 = Ifges(5,4) * t919 - Ifges(5,2) * t918 + Ifges(5,6) * t953;
t885 = Ifges(5,1) * t919 - Ifges(5,4) * t918 + Ifges(5,5) * t953;
t914 = Ifges(4,4) * t933 + Ifges(4,2) * t932 + Ifges(4,6) * t953;
t915 = Ifges(4,1) * t933 + Ifges(4,4) * t932 + Ifges(4,5) * t953;
t970 = -mrSges(4,1) * t869 - mrSges(5,1) * t837 + mrSges(4,2) * t870 + mrSges(5,2) * t838 - pkin(3) * t794 + pkin(4) * t825 - qJ(5) * t980 - t958 * t796 - t956 * t798 - t919 * t884 - t918 * t885 + t932 * t915 + Ifges(5,6) * t879 - Ifges(5,5) * t880 - t933 * t914 - Ifges(4,6) * t905 - Ifges(4,5) * t906 + (-Ifges(4,3) - Ifges(5,3)) * t952;
t993 = mrSges(3,1) * t921 - mrSges(3,2) * t922 + Ifges(3,5) * t940 + Ifges(3,6) * t941 + Ifges(3,3) * qJDD(2) + pkin(2) * t786 + (t930 * t961 - t931 * t965) * qJD(1) - t970;
t939 = (-mrSges(3,1) * t965 + mrSges(3,2) * t961) * qJD(1);
t987 = qJD(1) * t965;
t944 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t987;
t784 = m(3) * t921 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t940 + qJD(2) * t944 - t939 * t988 + t786;
t943 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t988;
t982 = -t791 * t960 + t964 * t792;
t785 = m(3) * t922 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t941 - qJD(2) * t943 + t939 * t987 + t982;
t983 = -t784 * t961 + t965 * t785;
t777 = m(2) * t947 - mrSges(2,1) * t967 - qJDD(1) * mrSges(2,2) + t983;
t934 = -t967 * pkin(7) + t975;
t807 = t958 * t812 + t956 * t813;
t805 = m(5) * t856 + t879 * mrSges(5,1) + t880 * mrSges(5,2) + t918 * t908 + t919 * t909 + t807;
t972 = m(4) * t907 - t905 * mrSges(4,1) + mrSges(4,2) * t906 - t932 * t923 + t925 * t933 + t805;
t969 = -m(3) * t934 + t941 * mrSges(3,1) - mrSges(3,2) * t940 - t943 * t988 + t944 * t987 - t972;
t800 = m(2) * t946 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t967 + t969;
t989 = t962 * t777 + t966 * t800;
t779 = t965 * t784 + t961 * t785;
t984 = t966 * t777 - t800 * t962;
t883 = Ifges(5,5) * t919 - Ifges(5,6) * t918 + Ifges(5,3) * t953;
t780 = mrSges(5,2) * t856 - mrSges(5,3) * t837 + Ifges(5,1) * t880 - Ifges(5,4) * t879 + Ifges(5,5) * t952 - qJ(5) * t807 - t796 * t956 + t798 * t958 - t883 * t918 - t884 * t953;
t971 = mrSges(7,1) * t826 - mrSges(7,2) * t827 + Ifges(7,5) * t844 + Ifges(7,6) * t843 + Ifges(7,3) * t877 + t872 * t846 - t871 * t847;
t787 = -t971 + (-Ifges(5,2) - Ifges(6,3)) * t879 + t953 * t885 + Ifges(5,6) * t952 - t919 * t883 + t900 * t861 - t901 * t860 + Ifges(5,4) * t880 - Ifges(6,6) * t867 - Ifges(6,5) * t868 - mrSges(5,1) * t856 + mrSges(5,3) * t838 - mrSges(6,1) * t830 + mrSges(6,2) * t831 - pkin(5) * t814 - pkin(4) * t807;
t913 = Ifges(4,5) * t933 + Ifges(4,6) * t932 + Ifges(4,3) * t953;
t773 = -mrSges(4,1) * t907 + mrSges(4,3) * t870 + Ifges(4,4) * t906 + Ifges(4,2) * t905 + Ifges(4,6) * t952 - pkin(3) * t805 + qJ(4) * t981 + t957 * t780 + t787 * t991 - t933 * t913 + t953 * t915;
t774 = mrSges(4,2) * t907 - mrSges(4,3) * t869 + Ifges(4,1) * t906 + Ifges(4,4) * t905 + Ifges(4,5) * t952 - qJ(4) * t794 + t780 * t991 - t957 * t787 + t932 * t913 - t953 * t914;
t929 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t961 + Ifges(3,6) * t965) * qJD(1);
t769 = -mrSges(3,1) * t934 + mrSges(3,3) * t922 + Ifges(3,4) * t940 + Ifges(3,2) * t941 + Ifges(3,6) * qJDD(2) - pkin(2) * t972 + pkin(8) * t982 + qJD(2) * t931 + t964 * t773 + t960 * t774 - t929 * t988;
t771 = mrSges(3,2) * t934 - mrSges(3,3) * t921 + Ifges(3,1) * t940 + Ifges(3,4) * t941 + Ifges(3,5) * qJDD(2) - pkin(8) * t786 - qJD(2) * t930 - t773 * t960 + t774 * t964 + t929 * t987;
t974 = mrSges(2,1) * t946 - mrSges(2,2) * t947 + Ifges(2,3) * qJDD(1) + pkin(1) * t969 + pkin(7) * t983 + t965 * t769 + t961 * t771;
t772 = mrSges(2,1) * g(3) + mrSges(2,3) * t947 + t967 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t779 - t993;
t767 = -mrSges(2,2) * g(3) - mrSges(2,3) * t946 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t967 - pkin(7) * t779 - t769 * t961 + t771 * t965;
t1 = [-m(1) * g(1) + t984; -m(1) * g(2) + t989; (-m(1) - m(2)) * g(3) + t779; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t989 + t966 * t767 - t962 * t772; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t984 + t962 * t767 + t966 * t772; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t974; t974; t993; -t970; t805; t825; t971;];
tauJB  = t1;
