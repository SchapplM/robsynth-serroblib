% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 15:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRPR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:55:45
% EndTime: 2019-05-06 14:55:58
% DurationCPUTime: 12.72s
% Computational Cost: add. (189353->366), mult. (409152->442), div. (0->0), fcn. (288028->10), ass. (0->146)
t989 = Ifges(5,1) + Ifges(6,1);
t980 = Ifges(5,4) - Ifges(6,5);
t979 = Ifges(5,5) + Ifges(6,4);
t988 = -Ifges(5,2) - Ifges(6,3);
t978 = Ifges(5,6) - Ifges(6,6);
t987 = Ifges(5,3) + Ifges(6,2);
t938 = sin(pkin(10));
t939 = cos(pkin(10));
t942 = sin(qJ(2));
t970 = qJD(1) * t942;
t917 = qJD(2) * t939 - t938 * t970;
t918 = qJD(2) * t938 + t939 * t970;
t941 = sin(qJ(4));
t982 = cos(qJ(4));
t888 = -t917 * t982 + t918 * t941;
t945 = cos(qJ(2));
t968 = qJD(1) * qJD(2);
t967 = t945 * t968;
t923 = qJDD(1) * t942 + t967;
t896 = qJDD(2) * t939 - t923 * t938;
t897 = qJDD(2) * t938 + t923 * t939;
t855 = -qJD(4) * t888 + t896 * t941 + t897 * t982;
t969 = qJD(1) * t945;
t898 = -pkin(3) * t969 - pkin(8) * t918;
t916 = t917 ^ 2;
t943 = sin(qJ(1));
t946 = cos(qJ(1));
t930 = -g(1) * t946 - g(2) * t943;
t948 = qJD(1) ^ 2;
t910 = -pkin(1) * t948 + qJDD(1) * pkin(7) + t930;
t892 = -g(3) * t945 - t910 * t942;
t921 = (-pkin(2) * t945 - qJ(3) * t942) * qJD(1);
t947 = qJD(2) ^ 2;
t955 = qJDD(2) * pkin(2) + t947 * qJ(3) - t921 * t970 - qJDD(3) + t892;
t953 = t896 * pkin(3) + t916 * pkin(8) - t898 * t918 + t955;
t933 = qJD(4) - t969;
t976 = t888 * t933;
t986 = (-t855 + t976) * qJ(5) - t953;
t929 = t943 * g(1) - g(2) * t946;
t909 = -qJDD(1) * pkin(1) - t948 * pkin(7) - t929;
t934 = t942 * t968;
t924 = qJDD(1) * t945 - t934;
t871 = (-t923 - t967) * qJ(3) + (-t924 + t934) * pkin(2) + t909;
t893 = -g(3) * t942 + t910 * t945;
t876 = -pkin(2) * t947 + qJDD(2) * qJ(3) + t921 * t969 + t893;
t845 = -0.2e1 * qJD(3) * t918 + t871 * t939 - t938 * t876;
t835 = (-t917 * t969 - t897) * pkin(8) + (t917 * t918 - t924) * pkin(3) + t845;
t846 = 0.2e1 * qJD(3) * t917 + t871 * t938 + t876 * t939;
t838 = -pkin(3) * t916 + pkin(8) * t896 + t898 * t969 + t846;
t825 = t835 * t982 - t838 * t941;
t889 = t917 * t941 + t918 * t982;
t866 = pkin(4) * t888 - qJ(5) * t889;
t920 = qJDD(4) - t924;
t932 = t933 ^ 2;
t823 = -t920 * pkin(4) - t932 * qJ(5) + t866 * t889 + qJDD(5) - t825;
t817 = (-t855 - t976) * pkin(9) + (t888 * t889 - t920) * pkin(5) + t823;
t826 = t835 * t941 + t838 * t982;
t983 = 2 * qJD(5);
t822 = -pkin(4) * t932 + qJ(5) * t920 - t866 * t888 + t933 * t983 + t826;
t854 = qJD(4) * t889 - t896 * t982 + t897 * t941;
t881 = -pkin(5) * t933 - pkin(9) * t889;
t887 = t888 ^ 2;
t818 = -pkin(5) * t887 + pkin(9) * t854 + t881 * t933 + t822;
t940 = sin(qJ(6));
t944 = cos(qJ(6));
t816 = t817 * t940 + t818 * t944;
t820 = -t887 * pkin(9) + (-pkin(4) - pkin(5)) * t854 + (-pkin(4) * t933 + t881 + t983) * t889 - t986;
t865 = t888 * t940 + t889 * t944;
t831 = -qJD(6) * t865 + t854 * t944 - t855 * t940;
t864 = t888 * t944 - t889 * t940;
t832 = qJD(6) * t864 + t854 * t940 + t855 * t944;
t931 = qJD(6) - t933;
t839 = Ifges(7,5) * t865 + Ifges(7,6) * t864 + Ifges(7,3) * t931;
t841 = Ifges(7,1) * t865 + Ifges(7,4) * t864 + Ifges(7,5) * t931;
t914 = qJDD(6) - t920;
t805 = -mrSges(7,1) * t820 + mrSges(7,3) * t816 + Ifges(7,4) * t832 + Ifges(7,2) * t831 + Ifges(7,6) * t914 - t839 * t865 + t841 * t931;
t815 = t817 * t944 - t818 * t940;
t840 = Ifges(7,4) * t865 + Ifges(7,2) * t864 + Ifges(7,6) * t931;
t806 = mrSges(7,2) * t820 - mrSges(7,3) * t815 + Ifges(7,1) * t832 + Ifges(7,4) * t831 + Ifges(7,5) * t914 + t839 * t864 - t840 * t931;
t824 = -0.2e1 * qJD(5) * t889 + (t889 * t933 + t854) * pkin(4) + t986;
t879 = -mrSges(6,1) * t933 + mrSges(6,2) * t889;
t880 = -mrSges(6,2) * t888 + mrSges(6,3) * t933;
t850 = -mrSges(7,2) * t931 + mrSges(7,3) * t864;
t851 = mrSges(7,1) * t931 - mrSges(7,3) * t865;
t961 = -m(7) * t820 + mrSges(7,1) * t831 - mrSges(7,2) * t832 + t850 * t864 - t851 * t865;
t810 = m(6) * t824 + mrSges(6,1) * t854 - mrSges(6,3) * t855 - t879 * t889 + t880 * t888 + t961;
t844 = -mrSges(7,1) * t864 + mrSges(7,2) * t865;
t812 = m(7) * t815 + mrSges(7,1) * t914 - mrSges(7,3) * t832 - t844 * t865 + t850 * t931;
t813 = m(7) * t816 - mrSges(7,2) * t914 + mrSges(7,3) * t831 + t844 * t864 - t851 * t931;
t963 = -t940 * t812 + t813 * t944;
t972 = -t888 * t980 + t889 * t989 + t933 * t979;
t974 = t888 * t978 - t889 * t979 - t933 * t987;
t791 = mrSges(5,1) * t953 - mrSges(6,1) * t824 + mrSges(6,2) * t822 + mrSges(5,3) * t826 - pkin(4) * t810 - pkin(5) * t961 - pkin(9) * t963 - t944 * t805 - t940 * t806 + t854 * t988 + t855 * t980 + t889 * t974 + t920 * t978 + t933 * t972;
t804 = t944 * t812 + t940 * t813;
t973 = t888 * t988 + t889 * t980 + t933 * t978;
t792 = -mrSges(5,2) * t953 + mrSges(6,2) * t823 - mrSges(5,3) * t825 - mrSges(6,3) * t824 - pkin(9) * t804 - qJ(5) * t810 - t940 * t805 + t944 * t806 - t854 * t980 + t855 * t989 + t888 * t974 + t920 * t979 - t933 * t973;
t882 = Ifges(4,5) * t918 + Ifges(4,6) * t917 - Ifges(4,3) * t969;
t884 = Ifges(4,1) * t918 + Ifges(4,4) * t917 - Ifges(4,5) * t969;
t877 = -mrSges(5,2) * t933 - mrSges(5,3) * t888;
t878 = mrSges(5,1) * t933 - mrSges(5,3) * t889;
t950 = -m(5) * t953 + mrSges(5,1) * t854 + mrSges(5,2) * t855 + t877 * t888 + t878 * t889 + t810;
t958 = m(6) * t822 + mrSges(6,3) * t920 + t879 * t933 + t963;
t867 = mrSges(6,1) * t888 - mrSges(6,3) * t889;
t971 = -mrSges(5,1) * t888 - mrSges(5,2) * t889 - t867;
t981 = -mrSges(5,3) - mrSges(6,2);
t799 = m(5) * t826 - t920 * mrSges(5,2) + t854 * t981 - t933 * t878 + t888 * t971 + t958;
t956 = -m(6) * t823 + mrSges(6,1) * t920 + t880 * t933 - t804;
t801 = m(5) * t825 + t920 * mrSges(5,1) + t855 * t981 + t933 * t877 + t889 * t971 + t956;
t964 = t799 * t982 - t801 * t941;
t776 = mrSges(4,1) * t955 + mrSges(4,3) * t846 + Ifges(4,4) * t897 + Ifges(4,2) * t896 - Ifges(4,6) * t924 - pkin(3) * t950 + pkin(8) * t964 + t791 * t982 + t941 * t792 - t918 * t882 - t884 * t969;
t796 = t799 * t941 + t801 * t982;
t883 = Ifges(4,4) * t918 + Ifges(4,2) * t917 - Ifges(4,6) * t969;
t777 = -mrSges(4,2) * t955 - mrSges(4,3) * t845 + Ifges(4,1) * t897 + Ifges(4,4) * t896 - Ifges(4,5) * t924 - pkin(8) * t796 - t791 * t941 + t792 * t982 + t882 * t917 + t883 * t969;
t890 = -mrSges(4,1) * t917 + mrSges(4,2) * t918;
t959 = mrSges(4,2) * t969 + mrSges(4,3) * t917;
t794 = m(4) * t845 - t924 * mrSges(4,1) - t897 * mrSges(4,3) - t918 * t890 - t959 * t969 + t796;
t895 = -mrSges(4,1) * t969 - mrSges(4,3) * t918;
t795 = m(4) * t846 + mrSges(4,2) * t924 + mrSges(4,3) * t896 + t890 * t917 + t895 * t969 + t964;
t790 = -t794 * t938 + t795 * t939;
t809 = -m(4) * t955 - mrSges(4,1) * t896 + mrSges(4,2) * t897 + t895 * t918 - t917 * t959 + t950;
t907 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t942 + Ifges(3,2) * t945) * qJD(1);
t908 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t942 + Ifges(3,4) * t945) * qJD(1);
t985 = mrSges(3,1) * t892 - mrSges(3,2) * t893 + Ifges(3,5) * t923 + Ifges(3,6) * t924 + Ifges(3,3) * qJDD(2) - pkin(2) * t809 + qJ(3) * t790 + t939 * t776 + t938 * t777 + (t907 * t942 - t908 * t945) * qJD(1);
t803 = t855 * mrSges(6,2) + t889 * t867 - t956;
t954 = -mrSges(7,1) * t815 + mrSges(7,2) * t816 - Ifges(7,5) * t832 - Ifges(7,6) * t831 - Ifges(7,3) * t914 - t840 * t865 + t864 * t841;
t984 = t978 * t854 - t979 * t855 - t972 * t888 - t973 * t889 - t987 * t920 - mrSges(5,1) * t825 + mrSges(6,1) * t823 + mrSges(5,2) * t826 - mrSges(6,3) * t822 + pkin(4) * t803 + pkin(5) * t804 - qJ(5) * (-t854 * mrSges(6,2) - t888 * t867 + t958) - t954;
t922 = (-mrSges(3,1) * t945 + mrSges(3,2) * t942) * qJD(1);
t926 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t970;
t788 = m(3) * t893 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t924 - qJD(2) * t926 + t922 * t969 + t790;
t927 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t969;
t808 = m(3) * t892 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t923 + qJD(2) * t927 - t922 * t970 - t809;
t965 = t788 * t945 - t808 * t942;
t780 = m(2) * t930 - mrSges(2,1) * t948 - qJDD(1) * mrSges(2,2) + t965;
t789 = t794 * t939 + t795 * t938;
t952 = -m(3) * t909 + mrSges(3,1) * t924 - mrSges(3,2) * t923 - t926 * t970 + t927 * t969 - t789;
t784 = m(2) * t929 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t948 + t952;
t975 = t780 * t943 + t784 * t946;
t782 = t788 * t942 + t808 * t945;
t966 = t780 * t946 - t784 * t943;
t906 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t942 + Ifges(3,6) * t945) * qJD(1);
t773 = mrSges(3,2) * t909 - mrSges(3,3) * t892 + Ifges(3,1) * t923 + Ifges(3,4) * t924 + Ifges(3,5) * qJDD(2) - qJ(3) * t789 - qJD(2) * t907 - t776 * t938 + t777 * t939 + t906 * t969;
t775 = t984 - t906 * t970 + Ifges(3,6) * qJDD(2) + (Ifges(3,2) + Ifges(4,3)) * t924 - t918 * t883 + Ifges(3,4) * t923 + qJD(2) * t908 - mrSges(3,1) * t909 + t917 * t884 + mrSges(3,3) * t893 - Ifges(4,6) * t896 - Ifges(4,5) * t897 - mrSges(4,1) * t845 + mrSges(4,2) * t846 - pkin(3) * t796 - pkin(2) * t789;
t957 = mrSges(2,1) * t929 - mrSges(2,2) * t930 + Ifges(2,3) * qJDD(1) + pkin(1) * t952 + pkin(7) * t965 + t773 * t942 + t775 * t945;
t771 = mrSges(2,1) * g(3) + mrSges(2,3) * t930 + t948 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t782 - t985;
t770 = -mrSges(2,2) * g(3) - mrSges(2,3) * t929 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t948 - pkin(7) * t782 + t773 * t945 - t775 * t942;
t1 = [-m(1) * g(1) + t966; -m(1) * g(2) + t975; (-m(1) - m(2)) * g(3) + t782; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t975 + t770 * t946 - t771 * t943; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t966 + t943 * t770 + t946 * t771; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t957; t957; t985; t809; -t984; t803; -t954;];
tauJB  = t1;
