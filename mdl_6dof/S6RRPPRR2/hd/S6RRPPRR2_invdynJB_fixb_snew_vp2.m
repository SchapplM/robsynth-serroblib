% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-06 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPPRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:45:53
% EndTime: 2019-05-06 09:46:24
% DurationCPUTime: 32.12s
% Computational Cost: add. (508893->387), mult. (1209016->489), div. (0->0), fcn. (881458->12), ass. (0->152)
t973 = -2 * qJD(3);
t939 = sin(qJ(2));
t943 = cos(qJ(2));
t963 = qJD(1) * qJD(2);
t922 = t939 * qJDD(1) + t943 * t963;
t940 = sin(qJ(1));
t944 = cos(qJ(1));
t929 = -t944 * g(1) - t940 * g(2);
t946 = qJD(1) ^ 2;
t917 = -t946 * pkin(1) + qJDD(1) * pkin(7) + t929;
t968 = t939 * t917;
t971 = pkin(2) * t946;
t874 = qJDD(2) * pkin(2) - t922 * qJ(3) - t968 + (qJ(3) * t963 + t939 * t971 - g(3)) * t943;
t902 = -t939 * g(3) + t943 * t917;
t923 = t943 * qJDD(1) - t939 * t963;
t966 = qJD(1) * t939;
t925 = qJD(2) * pkin(2) - qJ(3) * t966;
t933 = t943 ^ 2;
t876 = t923 * qJ(3) - qJD(2) * t925 - t933 * t971 + t902;
t935 = sin(pkin(10));
t969 = cos(pkin(10));
t911 = (t935 * t943 + t939 * t969) * qJD(1);
t855 = t874 * t969 - t935 * t876 + t911 * t973;
t965 = qJD(1) * t943;
t910 = t935 * t966 - t969 * t965;
t856 = t935 * t874 + t969 * t876 + t910 * t973;
t888 = t910 * pkin(3) - t911 * qJ(4);
t945 = qJD(2) ^ 2;
t842 = -t945 * pkin(3) + qJDD(2) * qJ(4) - t910 * t888 + t856;
t928 = t940 * g(1) - t944 * g(2);
t954 = -qJDD(1) * pkin(1) - t928;
t878 = -t923 * pkin(2) + qJDD(3) + t925 * t966 + (-qJ(3) * t933 - pkin(7)) * t946 + t954;
t894 = t935 * t922 - t969 * t923;
t895 = t922 * t969 + t935 * t923;
t845 = (qJD(2) * t910 - t895) * qJ(4) + (qJD(2) * t911 + t894) * pkin(3) + t878;
t934 = sin(pkin(11));
t936 = cos(pkin(11));
t900 = t934 * qJD(2) + t936 * t911;
t832 = -0.2e1 * qJD(4) * t900 - t934 * t842 + t936 * t845;
t883 = t934 * qJDD(2) + t936 * t895;
t899 = t936 * qJD(2) - t934 * t911;
t823 = (t910 * t899 - t883) * pkin(8) + (t899 * t900 + t894) * pkin(4) + t832;
t833 = 0.2e1 * qJD(4) * t899 + t936 * t842 + t934 * t845;
t880 = t910 * pkin(4) - t900 * pkin(8);
t882 = t936 * qJDD(2) - t934 * t895;
t898 = t899 ^ 2;
t825 = -t898 * pkin(4) + t882 * pkin(8) - t910 * t880 + t833;
t938 = sin(qJ(5));
t942 = cos(qJ(5));
t818 = t942 * t823 - t938 * t825;
t870 = t942 * t899 - t938 * t900;
t848 = t870 * qJD(5) + t938 * t882 + t942 * t883;
t871 = t938 * t899 + t942 * t900;
t893 = qJDD(5) + t894;
t909 = qJD(5) + t910;
t816 = (t870 * t909 - t848) * pkin(9) + (t870 * t871 + t893) * pkin(5) + t818;
t819 = t938 * t823 + t942 * t825;
t847 = -t871 * qJD(5) + t942 * t882 - t938 * t883;
t862 = t909 * pkin(5) - t871 * pkin(9);
t869 = t870 ^ 2;
t817 = -t869 * pkin(5) + t847 * pkin(9) - t909 * t862 + t819;
t937 = sin(qJ(6));
t941 = cos(qJ(6));
t815 = t937 * t816 + t941 * t817;
t841 = -qJDD(2) * pkin(3) - t945 * qJ(4) + t911 * t888 + qJDD(4) - t855;
t834 = -t882 * pkin(4) - t898 * pkin(8) + t900 * t880 + t841;
t820 = -t847 * pkin(5) - t869 * pkin(9) + t871 * t862 + t834;
t858 = t937 * t870 + t941 * t871;
t829 = -t858 * qJD(6) + t941 * t847 - t937 * t848;
t857 = t941 * t870 - t937 * t871;
t830 = t857 * qJD(6) + t937 * t847 + t941 * t848;
t905 = qJD(6) + t909;
t835 = Ifges(7,5) * t858 + Ifges(7,6) * t857 + Ifges(7,3) * t905;
t837 = Ifges(7,1) * t858 + Ifges(7,4) * t857 + Ifges(7,5) * t905;
t891 = qJDD(6) + t893;
t803 = -mrSges(7,1) * t820 + mrSges(7,3) * t815 + Ifges(7,4) * t830 + Ifges(7,2) * t829 + Ifges(7,6) * t891 - t858 * t835 + t905 * t837;
t814 = t941 * t816 - t937 * t817;
t836 = Ifges(7,4) * t858 + Ifges(7,2) * t857 + Ifges(7,6) * t905;
t804 = mrSges(7,2) * t820 - mrSges(7,3) * t814 + Ifges(7,1) * t830 + Ifges(7,4) * t829 + Ifges(7,5) * t891 + t857 * t835 - t905 * t836;
t851 = Ifges(6,5) * t871 + Ifges(6,6) * t870 + Ifges(6,3) * t909;
t853 = Ifges(6,1) * t871 + Ifges(6,4) * t870 + Ifges(6,5) * t909;
t849 = -t905 * mrSges(7,2) + t857 * mrSges(7,3);
t850 = t905 * mrSges(7,1) - t858 * mrSges(7,3);
t952 = m(7) * t820 - t829 * mrSges(7,1) + t830 * mrSges(7,2) - t857 * t849 + t858 * t850;
t839 = -t857 * mrSges(7,1) + t858 * mrSges(7,2);
t808 = m(7) * t814 + t891 * mrSges(7,1) - t830 * mrSges(7,3) - t858 * t839 + t905 * t849;
t809 = m(7) * t815 - t891 * mrSges(7,2) + t829 * mrSges(7,3) + t857 * t839 - t905 * t850;
t958 = -t937 * t808 + t941 * t809;
t789 = -mrSges(6,1) * t834 + mrSges(6,3) * t819 + Ifges(6,4) * t848 + Ifges(6,2) * t847 + Ifges(6,6) * t893 - pkin(5) * t952 + pkin(9) * t958 + t941 * t803 + t937 * t804 - t871 * t851 + t909 * t853;
t802 = t941 * t808 + t937 * t809;
t852 = Ifges(6,4) * t871 + Ifges(6,2) * t870 + Ifges(6,6) * t909;
t790 = mrSges(6,2) * t834 - mrSges(6,3) * t818 + Ifges(6,1) * t848 + Ifges(6,4) * t847 + Ifges(6,5) * t893 - pkin(9) * t802 - t937 * t803 + t941 * t804 + t870 * t851 - t909 * t852;
t863 = Ifges(5,5) * t900 + Ifges(5,6) * t899 + Ifges(5,3) * t910;
t865 = Ifges(5,1) * t900 + Ifges(5,4) * t899 + Ifges(5,5) * t910;
t860 = -t909 * mrSges(6,2) + t870 * mrSges(6,3);
t861 = t909 * mrSges(6,1) - t871 * mrSges(6,3);
t950 = m(6) * t834 - t847 * mrSges(6,1) + t848 * mrSges(6,2) - t870 * t860 + t871 * t861 + t952;
t859 = -t870 * mrSges(6,1) + t871 * mrSges(6,2);
t800 = m(6) * t818 + t893 * mrSges(6,1) - t848 * mrSges(6,3) - t871 * t859 + t909 * t860 + t802;
t801 = m(6) * t819 - t893 * mrSges(6,2) + t847 * mrSges(6,3) + t870 * t859 - t909 * t861 + t958;
t959 = -t938 * t800 + t942 * t801;
t772 = -mrSges(5,1) * t841 + mrSges(5,3) * t833 + Ifges(5,4) * t883 + Ifges(5,2) * t882 + Ifges(5,6) * t894 - pkin(4) * t950 + pkin(8) * t959 + t942 * t789 + t938 * t790 - t900 * t863 + t910 * t865;
t796 = t942 * t800 + t938 * t801;
t864 = Ifges(5,4) * t900 + Ifges(5,2) * t899 + Ifges(5,6) * t910;
t773 = mrSges(5,2) * t841 - mrSges(5,3) * t832 + Ifges(5,1) * t883 + Ifges(5,4) * t882 + Ifges(5,5) * t894 - pkin(8) * t796 - t938 * t789 + t942 * t790 + t899 * t863 - t910 * t864;
t875 = -t899 * mrSges(5,1) + t900 * mrSges(5,2);
t957 = -t910 * mrSges(5,2) + t899 * mrSges(5,3);
t794 = m(5) * t832 + t894 * mrSges(5,1) - t883 * mrSges(5,3) - t900 * t875 + t910 * t957 + t796;
t879 = t910 * mrSges(5,1) - t900 * mrSges(5,3);
t795 = m(5) * t833 - t894 * mrSges(5,2) + t882 * mrSges(5,3) + t899 * t875 - t910 * t879 + t959;
t788 = -t934 * t794 + t936 * t795;
t889 = t910 * mrSges(4,1) + t911 * mrSges(4,2);
t904 = qJD(2) * mrSges(4,1) - t911 * mrSges(4,3);
t785 = m(4) * t856 - qJDD(2) * mrSges(4,2) - t894 * mrSges(4,3) - qJD(2) * t904 - t910 * t889 + t788;
t812 = m(5) * t841 - t882 * mrSges(5,1) + t883 * mrSges(5,2) + t900 * t879 - t899 * t957 + t950;
t903 = -qJD(2) * mrSges(4,2) - t910 * mrSges(4,3);
t811 = m(4) * t855 + qJDD(2) * mrSges(4,1) - t895 * mrSges(4,3) + qJD(2) * t903 - t911 * t889 - t812;
t779 = t935 * t785 + t969 * t811;
t885 = Ifges(4,4) * t911 - Ifges(4,2) * t910 + Ifges(4,6) * qJD(2);
t886 = Ifges(4,1) * t911 - Ifges(4,4) * t910 + Ifges(4,5) * qJD(2);
t901 = -t943 * g(3) - t968;
t913 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t939 + Ifges(3,2) * t943) * qJD(1);
t914 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t939 + Ifges(3,4) * t943) * qJD(1);
t972 = (t939 * t913 - t943 * t914) * qJD(1) - (-Ifges(4,3) - Ifges(3,3)) * qJDD(2) + mrSges(3,1) * t901 + mrSges(4,1) * t855 - mrSges(3,2) * t902 - mrSges(4,2) * t856 + Ifges(3,5) * t922 + Ifges(4,5) * t895 + Ifges(3,6) * t923 - Ifges(4,6) * t894 + pkin(2) * t779 - pkin(3) * t812 + qJ(4) * t788 + t936 * t772 + t934 * t773 + t911 * t885 + t910 * t886;
t921 = (-mrSges(3,1) * t943 + mrSges(3,2) * t939) * qJD(1);
t927 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t965;
t777 = m(3) * t901 + qJDD(2) * mrSges(3,1) - t922 * mrSges(3,3) + qJD(2) * t927 - t921 * t966 + t779;
t926 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t966;
t960 = t969 * t785 - t935 * t811;
t778 = m(3) * t902 - qJDD(2) * mrSges(3,2) + t923 * mrSges(3,3) - qJD(2) * t926 + t921 * t965 + t960;
t961 = -t939 * t777 + t943 * t778;
t768 = m(2) * t929 - t946 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t961;
t787 = t936 * t794 + t934 * t795;
t786 = m(4) * t878 + t894 * mrSges(4,1) + t895 * mrSges(4,2) + t910 * t903 + t911 * t904 + t787;
t916 = -t946 * pkin(7) + t954;
t949 = -m(3) * t916 + t923 * mrSges(3,1) - t922 * mrSges(3,2) - t926 * t966 + t927 * t965 - t786;
t781 = m(2) * t928 + qJDD(1) * mrSges(2,1) - t946 * mrSges(2,2) + t949;
t967 = t940 * t768 + t944 * t781;
t770 = t943 * t777 + t939 * t778;
t962 = t944 * t768 - t940 * t781;
t884 = Ifges(4,5) * t911 - Ifges(4,6) * t910 + Ifges(4,3) * qJD(2);
t765 = mrSges(4,2) * t878 - mrSges(4,3) * t855 + Ifges(4,1) * t895 - Ifges(4,4) * t894 + Ifges(4,5) * qJDD(2) - qJ(4) * t787 - qJD(2) * t885 - t934 * t772 + t936 * t773 - t910 * t884;
t951 = -mrSges(7,1) * t814 + mrSges(7,2) * t815 - Ifges(7,5) * t830 - Ifges(7,6) * t829 - Ifges(7,3) * t891 - t858 * t836 + t857 * t837;
t948 = mrSges(6,1) * t818 - mrSges(6,2) * t819 + Ifges(6,5) * t848 + Ifges(6,6) * t847 + Ifges(6,3) * t893 + pkin(5) * t802 + t871 * t852 - t870 * t853 - t951;
t771 = (-Ifges(5,3) - Ifges(4,2)) * t894 - t911 * t884 + Ifges(4,4) * t895 + t899 * t865 - t900 * t864 - Ifges(5,6) * t882 - Ifges(5,5) * t883 + qJD(2) * t886 - mrSges(4,1) * t878 + mrSges(4,3) * t856 - mrSges(5,1) * t832 + mrSges(5,2) * t833 - t948 - pkin(4) * t796 - pkin(3) * t787 + Ifges(4,6) * qJDD(2);
t912 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t939 + Ifges(3,6) * t943) * qJD(1);
t761 = -mrSges(3,1) * t916 + mrSges(3,3) * t902 + Ifges(3,4) * t922 + Ifges(3,2) * t923 + Ifges(3,6) * qJDD(2) - pkin(2) * t786 + qJ(3) * t960 + qJD(2) * t914 + t935 * t765 + t771 * t969 - t912 * t966;
t764 = mrSges(3,2) * t916 - mrSges(3,3) * t901 + Ifges(3,1) * t922 + Ifges(3,4) * t923 + Ifges(3,5) * qJDD(2) - qJ(3) * t779 - qJD(2) * t913 + t765 * t969 - t935 * t771 + t912 * t965;
t953 = mrSges(2,1) * t928 - mrSges(2,2) * t929 + Ifges(2,3) * qJDD(1) + pkin(1) * t949 + pkin(7) * t961 + t943 * t761 + t939 * t764;
t762 = mrSges(2,1) * g(3) + mrSges(2,3) * t929 + t946 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t770 - t972;
t759 = -mrSges(2,2) * g(3) - mrSges(2,3) * t928 + Ifges(2,5) * qJDD(1) - t946 * Ifges(2,6) - pkin(7) * t770 - t939 * t761 + t943 * t764;
t1 = [-m(1) * g(1) + t962; -m(1) * g(2) + t967; (-m(1) - m(2)) * g(3) + t770; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t967 + t944 * t759 - t940 * t762; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t962 + t940 * t759 + t944 * t762; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t953; t953; t972; t786; t812; t948; -t951;];
tauJB  = t1;
