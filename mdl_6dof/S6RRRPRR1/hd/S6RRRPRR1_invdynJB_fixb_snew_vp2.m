% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 09:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:43:11
% EndTime: 2019-05-07 09:43:55
% DurationCPUTime: 43.71s
% Computational Cost: add. (641406->388), mult. (1500162->490), div. (0->0), fcn. (1142144->12), ass. (0->154)
t961 = sin(qJ(2));
t966 = cos(qJ(2));
t986 = qJD(1) * qJD(2);
t937 = qJDD(1) * t961 + t966 * t986;
t962 = sin(qJ(1));
t967 = cos(qJ(1));
t944 = -g(1) * t967 - g(2) * t962;
t968 = qJD(1) ^ 2;
t932 = -pkin(1) * t968 + qJDD(1) * pkin(7) + t944;
t990 = t961 * t932;
t991 = pkin(2) * t968;
t897 = qJDD(2) * pkin(2) - t937 * pkin(8) - t990 + (pkin(8) * t986 + t961 * t991 - g(3)) * t966;
t920 = -g(3) * t961 + t966 * t932;
t938 = qJDD(1) * t966 - t961 * t986;
t988 = qJD(1) * t961;
t942 = qJD(2) * pkin(2) - pkin(8) * t988;
t955 = t966 ^ 2;
t898 = pkin(8) * t938 - qJD(2) * t942 - t955 * t991 + t920;
t960 = sin(qJ(3));
t965 = cos(qJ(3));
t871 = t965 * t897 - t960 * t898;
t929 = (-t960 * t961 + t965 * t966) * qJD(1);
t903 = qJD(3) * t929 + t937 * t965 + t938 * t960;
t930 = (t960 * t966 + t961 * t965) * qJD(1);
t952 = qJDD(2) + qJDD(3);
t953 = qJD(2) + qJD(3);
t854 = (t929 * t953 - t903) * qJ(4) + (t929 * t930 + t952) * pkin(3) + t871;
t872 = t960 * t897 + t965 * t898;
t902 = -qJD(3) * t930 - t937 * t960 + t938 * t965;
t922 = pkin(3) * t953 - qJ(4) * t930;
t925 = t929 ^ 2;
t856 = -pkin(3) * t925 + qJ(4) * t902 - t922 * t953 + t872;
t956 = sin(pkin(11));
t957 = cos(pkin(11));
t917 = t929 * t956 + t930 * t957;
t833 = -0.2e1 * qJD(4) * t917 + t957 * t854 - t956 * t856;
t878 = t902 * t956 + t903 * t957;
t916 = t929 * t957 - t930 * t956;
t829 = (t916 * t953 - t878) * pkin(9) + (t916 * t917 + t952) * pkin(4) + t833;
t834 = 0.2e1 * qJD(4) * t916 + t956 * t854 + t957 * t856;
t877 = t902 * t957 - t903 * t956;
t907 = pkin(4) * t953 - pkin(9) * t917;
t913 = t916 ^ 2;
t831 = -pkin(4) * t913 + pkin(9) * t877 - t907 * t953 + t834;
t959 = sin(qJ(5));
t964 = cos(qJ(5));
t826 = t959 * t829 + t964 * t831;
t891 = t916 * t959 + t917 * t964;
t845 = -qJD(5) * t891 + t877 * t964 - t878 * t959;
t890 = t916 * t964 - t917 * t959;
t865 = -mrSges(6,1) * t890 + mrSges(6,2) * t891;
t950 = qJD(5) + t953;
t882 = mrSges(6,1) * t950 - mrSges(6,3) * t891;
t949 = qJDD(5) + t952;
t866 = -pkin(5) * t890 - pkin(10) * t891;
t948 = t950 ^ 2;
t823 = -pkin(5) * t948 + pkin(10) * t949 + t866 * t890 + t826;
t943 = t962 * g(1) - t967 * g(2);
t978 = -qJDD(1) * pkin(1) - t943;
t904 = -t938 * pkin(2) + t942 * t988 + (-pkin(8) * t955 - pkin(7)) * t968 + t978;
t868 = -t902 * pkin(3) - t925 * qJ(4) + t930 * t922 + qJDD(4) + t904;
t839 = -t877 * pkin(4) - t913 * pkin(9) + t917 * t907 + t868;
t846 = qJD(5) * t890 + t877 * t959 + t878 * t964;
t827 = (t891 * t950 - t845) * pkin(5) + (-t890 * t950 - t846) * pkin(10) + t839;
t958 = sin(qJ(6));
t963 = cos(qJ(6));
t820 = -t823 * t958 + t827 * t963;
t879 = -t891 * t958 + t950 * t963;
t837 = qJD(6) * t879 + t846 * t963 + t949 * t958;
t844 = qJDD(6) - t845;
t880 = t891 * t963 + t950 * t958;
t857 = -mrSges(7,1) * t879 + mrSges(7,2) * t880;
t886 = qJD(6) - t890;
t858 = -mrSges(7,2) * t886 + mrSges(7,3) * t879;
t816 = m(7) * t820 + mrSges(7,1) * t844 - mrSges(7,3) * t837 - t857 * t880 + t858 * t886;
t821 = t823 * t963 + t827 * t958;
t836 = -qJD(6) * t880 - t846 * t958 + t949 * t963;
t859 = mrSges(7,1) * t886 - mrSges(7,3) * t880;
t817 = m(7) * t821 - mrSges(7,2) * t844 + mrSges(7,3) * t836 + t857 * t879 - t859 * t886;
t980 = -t816 * t958 + t963 * t817;
t802 = m(6) * t826 - mrSges(6,2) * t949 + mrSges(6,3) * t845 + t865 * t890 - t882 * t950 + t980;
t825 = t829 * t964 - t831 * t959;
t881 = -mrSges(6,2) * t950 + mrSges(6,3) * t890;
t822 = -pkin(5) * t949 - pkin(10) * t948 + t866 * t891 - t825;
t975 = -m(7) * t822 + t836 * mrSges(7,1) - mrSges(7,2) * t837 + t879 * t858 - t859 * t880;
t812 = m(6) * t825 + mrSges(6,1) * t949 - mrSges(6,3) * t846 - t865 * t891 + t881 * t950 + t975;
t796 = t959 * t802 + t964 * t812;
t892 = -mrSges(5,1) * t916 + mrSges(5,2) * t917;
t905 = -mrSges(5,2) * t953 + mrSges(5,3) * t916;
t793 = m(5) * t833 + mrSges(5,1) * t952 - mrSges(5,3) * t878 - t892 * t917 + t905 * t953 + t796;
t906 = mrSges(5,1) * t953 - mrSges(5,3) * t917;
t981 = t964 * t802 - t812 * t959;
t794 = m(5) * t834 - mrSges(5,2) * t952 + mrSges(5,3) * t877 + t892 * t916 - t906 * t953 + t981;
t787 = t957 * t793 + t956 * t794;
t918 = -mrSges(4,1) * t929 + mrSges(4,2) * t930;
t921 = -mrSges(4,2) * t953 + mrSges(4,3) * t929;
t784 = m(4) * t871 + mrSges(4,1) * t952 - mrSges(4,3) * t903 - t918 * t930 + t921 * t953 + t787;
t923 = mrSges(4,1) * t953 - mrSges(4,3) * t930;
t982 = -t793 * t956 + t957 * t794;
t785 = m(4) * t872 - mrSges(4,2) * t952 + mrSges(4,3) * t902 + t918 * t929 - t923 * t953 + t982;
t779 = t965 * t784 + t960 * t785;
t919 = -t966 * g(3) - t990;
t927 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t961 + Ifges(3,2) * t966) * qJD(1);
t928 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t961 + Ifges(3,4) * t966) * qJD(1);
t884 = Ifges(5,4) * t917 + Ifges(5,2) * t916 + Ifges(5,6) * t953;
t885 = Ifges(5,1) * t917 + Ifges(5,4) * t916 + Ifges(5,5) * t953;
t911 = Ifges(4,4) * t930 + Ifges(4,2) * t929 + Ifges(4,6) * t953;
t912 = Ifges(4,1) * t930 + Ifges(4,4) * t929 + Ifges(4,5) * t953;
t847 = Ifges(7,5) * t880 + Ifges(7,6) * t879 + Ifges(7,3) * t886;
t849 = Ifges(7,1) * t880 + Ifges(7,4) * t879 + Ifges(7,5) * t886;
t809 = -mrSges(7,1) * t822 + mrSges(7,3) * t821 + Ifges(7,4) * t837 + Ifges(7,2) * t836 + Ifges(7,6) * t844 - t847 * t880 + t849 * t886;
t848 = Ifges(7,4) * t880 + Ifges(7,2) * t879 + Ifges(7,6) * t886;
t810 = mrSges(7,2) * t822 - mrSges(7,3) * t820 + Ifges(7,1) * t837 + Ifges(7,4) * t836 + Ifges(7,5) * t844 + t847 * t879 - t848 * t886;
t861 = Ifges(6,4) * t891 + Ifges(6,2) * t890 + Ifges(6,6) * t950;
t862 = Ifges(6,1) * t891 + Ifges(6,4) * t890 + Ifges(6,5) * t950;
t974 = -mrSges(6,1) * t825 + mrSges(6,2) * t826 - Ifges(6,5) * t846 - Ifges(6,6) * t845 - Ifges(6,3) * t949 - pkin(5) * t975 - pkin(10) * t980 - t963 * t809 - t958 * t810 - t891 * t861 + t890 * t862;
t970 = mrSges(4,1) * t871 + mrSges(5,1) * t833 - mrSges(4,2) * t872 - mrSges(5,2) * t834 + Ifges(4,5) * t903 + Ifges(5,5) * t878 + Ifges(4,6) * t902 + Ifges(5,6) * t877 + pkin(3) * t787 + pkin(4) * t796 + t917 * t884 - t916 * t885 + t930 * t911 - t929 * t912 - t974 + (Ifges(5,3) + Ifges(4,3)) * t952;
t992 = mrSges(3,1) * t919 - mrSges(3,2) * t920 + Ifges(3,5) * t937 + Ifges(3,6) * t938 + Ifges(3,3) * qJDD(2) + pkin(2) * t779 + (t927 * t961 - t928 * t966) * qJD(1) + t970;
t936 = (-mrSges(3,1) * t966 + mrSges(3,2) * t961) * qJD(1);
t987 = qJD(1) * t966;
t941 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t987;
t777 = m(3) * t919 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t937 + qJD(2) * t941 - t936 * t988 + t779;
t940 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t988;
t983 = -t784 * t960 + t965 * t785;
t778 = m(3) * t920 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t938 - qJD(2) * t940 + t936 * t987 + t983;
t984 = -t777 * t961 + t966 * t778;
t770 = m(2) * t944 - mrSges(2,1) * t968 - qJDD(1) * mrSges(2,2) + t984;
t931 = -t968 * pkin(7) + t978;
t805 = t963 * t816 + t958 * t817;
t977 = m(6) * t839 - t845 * mrSges(6,1) + t846 * mrSges(6,2) - t890 * t881 + t891 * t882 + t805;
t803 = m(5) * t868 - t877 * mrSges(5,1) + t878 * mrSges(5,2) - t916 * t905 + t917 * t906 + t977;
t972 = m(4) * t904 - t902 * mrSges(4,1) + t903 * mrSges(4,2) - t929 * t921 + t930 * t923 + t803;
t971 = -m(3) * t931 + t938 * mrSges(3,1) - t937 * mrSges(3,2) - t940 * t988 + t941 * t987 - t972;
t798 = m(2) * t943 + qJDD(1) * mrSges(2,1) - t968 * mrSges(2,2) + t971;
t989 = t962 * t770 + t967 * t798;
t772 = t966 * t777 + t961 * t778;
t985 = t967 * t770 - t798 * t962;
t860 = Ifges(6,5) * t891 + Ifges(6,6) * t890 + Ifges(6,3) * t950;
t788 = mrSges(6,2) * t839 - mrSges(6,3) * t825 + Ifges(6,1) * t846 + Ifges(6,4) * t845 + Ifges(6,5) * t949 - pkin(10) * t805 - t809 * t958 + t810 * t963 + t860 * t890 - t861 * t950;
t973 = mrSges(7,1) * t820 - mrSges(7,2) * t821 + Ifges(7,5) * t837 + Ifges(7,6) * t836 + Ifges(7,3) * t844 + t848 * t880 - t849 * t879;
t789 = -mrSges(6,1) * t839 + mrSges(6,3) * t826 + Ifges(6,4) * t846 + Ifges(6,2) * t845 + Ifges(6,6) * t949 - pkin(5) * t805 - t860 * t891 + t862 * t950 - t973;
t883 = Ifges(5,5) * t917 + Ifges(5,6) * t916 + Ifges(5,3) * t953;
t773 = -mrSges(5,1) * t868 + mrSges(5,3) * t834 + Ifges(5,4) * t878 + Ifges(5,2) * t877 + Ifges(5,6) * t952 - pkin(4) * t977 + pkin(9) * t981 + t959 * t788 + t964 * t789 - t917 * t883 + t953 * t885;
t780 = mrSges(5,2) * t868 - mrSges(5,3) * t833 + Ifges(5,1) * t878 + Ifges(5,4) * t877 + Ifges(5,5) * t952 - pkin(9) * t796 + t788 * t964 - t789 * t959 + t883 * t916 - t884 * t953;
t910 = Ifges(4,5) * t930 + Ifges(4,6) * t929 + Ifges(4,3) * t953;
t766 = -mrSges(4,1) * t904 + mrSges(4,3) * t872 + Ifges(4,4) * t903 + Ifges(4,2) * t902 + Ifges(4,6) * t952 - pkin(3) * t803 + qJ(4) * t982 + t957 * t773 + t956 * t780 - t930 * t910 + t953 * t912;
t767 = mrSges(4,2) * t904 - mrSges(4,3) * t871 + Ifges(4,1) * t903 + Ifges(4,4) * t902 + Ifges(4,5) * t952 - qJ(4) * t787 - t773 * t956 + t780 * t957 + t910 * t929 - t911 * t953;
t926 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t961 + Ifges(3,6) * t966) * qJD(1);
t762 = -mrSges(3,1) * t931 + mrSges(3,3) * t920 + Ifges(3,4) * t937 + Ifges(3,2) * t938 + Ifges(3,6) * qJDD(2) - pkin(2) * t972 + pkin(8) * t983 + qJD(2) * t928 + t965 * t766 + t960 * t767 - t926 * t988;
t764 = mrSges(3,2) * t931 - mrSges(3,3) * t919 + Ifges(3,1) * t937 + Ifges(3,4) * t938 + Ifges(3,5) * qJDD(2) - pkin(8) * t779 - qJD(2) * t927 - t766 * t960 + t767 * t965 + t926 * t987;
t976 = mrSges(2,1) * t943 - mrSges(2,2) * t944 + Ifges(2,3) * qJDD(1) + pkin(1) * t971 + pkin(7) * t984 + t966 * t762 + t961 * t764;
t765 = mrSges(2,1) * g(3) + mrSges(2,3) * t944 + t968 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t772 - t992;
t760 = -mrSges(2,2) * g(3) - mrSges(2,3) * t943 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t968 - pkin(7) * t772 - t762 * t961 + t764 * t966;
t1 = [-m(1) * g(1) + t985; -m(1) * g(2) + t989; (-m(1) - m(2)) * g(3) + t772; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t989 + t967 * t760 - t962 * t765; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t985 + t962 * t760 + t967 * t765; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t976; t976; t992; t970; t803; -t974; t973;];
tauJB  = t1;
