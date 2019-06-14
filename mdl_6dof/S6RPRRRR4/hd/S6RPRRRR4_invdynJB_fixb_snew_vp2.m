% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:14:17
% EndTime: 2019-05-06 03:14:51
% DurationCPUTime: 35.31s
% Computational Cost: add. (559289->367), mult. (1398177->456), div. (0->0), fcn. (1132509->12), ass. (0->156)
t956 = qJD(1) ^ 2;
t945 = cos(pkin(11));
t988 = pkin(2) * t945;
t944 = sin(pkin(11));
t987 = t944 * mrSges(3,2);
t941 = t945 ^ 2;
t986 = t941 * t956;
t950 = sin(qJ(1));
t955 = cos(qJ(1));
t927 = -t955 * g(1) - t950 * g(2);
t922 = -t956 * pkin(1) + qJDD(1) * qJ(2) + t927;
t982 = qJD(1) * qJD(2);
t980 = -t945 * g(3) - 0.2e1 * t944 * t982;
t897 = (-pkin(7) * qJDD(1) + t956 * t988 - t922) * t944 + t980;
t913 = -t944 * g(3) + (t922 + 0.2e1 * t982) * t945;
t981 = qJDD(1) * t945;
t898 = -pkin(2) * t986 + pkin(7) * t981 + t913;
t949 = sin(qJ(3));
t954 = cos(qJ(3));
t877 = t954 * t897 - t949 * t898;
t969 = t944 * t954 + t945 * t949;
t968 = -t944 * t949 + t945 * t954;
t920 = t968 * qJD(1);
t983 = t920 * qJD(3);
t911 = qJDD(1) * t969 + t983;
t921 = t969 * qJD(1);
t856 = (-t911 + t983) * pkin(8) + (t920 * t921 + qJDD(3)) * pkin(3) + t877;
t878 = t949 * t897 + t954 * t898;
t910 = -t921 * qJD(3) + qJDD(1) * t968;
t916 = qJD(3) * pkin(3) - t921 * pkin(8);
t919 = t920 ^ 2;
t865 = -t919 * pkin(3) + t910 * pkin(8) - qJD(3) * t916 + t878;
t948 = sin(qJ(4));
t953 = cos(qJ(4));
t837 = t953 * t856 - t948 * t865;
t903 = t953 * t920 - t948 * t921;
t874 = t903 * qJD(4) + t948 * t910 + t953 * t911;
t904 = t948 * t920 + t953 * t921;
t939 = qJDD(3) + qJDD(4);
t942 = qJD(3) + qJD(4);
t828 = (t903 * t942 - t874) * pkin(9) + (t903 * t904 + t939) * pkin(4) + t837;
t838 = t948 * t856 + t953 * t865;
t873 = -t904 * qJD(4) + t953 * t910 - t948 * t911;
t895 = t942 * pkin(4) - t904 * pkin(9);
t899 = t903 ^ 2;
t830 = -t899 * pkin(4) + t873 * pkin(9) - t942 * t895 + t838;
t947 = sin(qJ(5));
t952 = cos(qJ(5));
t826 = t947 * t828 + t952 * t830;
t889 = t947 * t903 + t952 * t904;
t844 = -t889 * qJD(5) + t952 * t873 - t947 * t874;
t888 = t952 * t903 - t947 * t904;
t862 = -t888 * mrSges(6,1) + t889 * mrSges(6,2);
t937 = qJD(5) + t942;
t880 = t937 * mrSges(6,1) - t889 * mrSges(6,3);
t936 = qJDD(5) + t939;
t864 = -t888 * pkin(5) - t889 * pkin(10);
t935 = t937 ^ 2;
t822 = -t935 * pkin(5) + t936 * pkin(10) + t888 * t864 + t826;
t940 = t944 ^ 2;
t926 = t950 * g(1) - t955 * g(2);
t973 = qJDD(2) - t926;
t909 = (-pkin(1) - t988) * qJDD(1) + (-qJ(2) + (-t940 - t941) * pkin(7)) * t956 + t973;
t868 = -t910 * pkin(3) - t919 * pkin(8) + t921 * t916 + t909;
t835 = -t873 * pkin(4) - t899 * pkin(9) + t904 * t895 + t868;
t845 = t888 * qJD(5) + t947 * t873 + t952 * t874;
t823 = t835 + (t889 * t937 - t844) * pkin(5) + (-t888 * t937 - t845) * pkin(10);
t946 = sin(qJ(6));
t951 = cos(qJ(6));
t819 = -t946 * t822 + t951 * t823;
t875 = -t946 * t889 + t951 * t937;
t833 = t875 * qJD(6) + t951 * t845 + t946 * t936;
t843 = qJDD(6) - t844;
t876 = t951 * t889 + t946 * t937;
t851 = -t875 * mrSges(7,1) + t876 * mrSges(7,2);
t884 = qJD(6) - t888;
t854 = -t884 * mrSges(7,2) + t875 * mrSges(7,3);
t815 = m(7) * t819 + t843 * mrSges(7,1) - t833 * mrSges(7,3) - t876 * t851 + t884 * t854;
t820 = t951 * t822 + t946 * t823;
t832 = -t876 * qJD(6) - t946 * t845 + t951 * t936;
t855 = t884 * mrSges(7,1) - t876 * mrSges(7,3);
t816 = m(7) * t820 - t843 * mrSges(7,2) + t832 * mrSges(7,3) + t875 * t851 - t884 * t855;
t974 = -t946 * t815 + t951 * t816;
t802 = m(6) * t826 - t936 * mrSges(6,2) + t844 * mrSges(6,3) + t888 * t862 - t937 * t880 + t974;
t825 = t952 * t828 - t947 * t830;
t879 = -t937 * mrSges(6,2) + t888 * mrSges(6,3);
t821 = -t936 * pkin(5) - t935 * pkin(10) + t889 * t864 - t825;
t964 = -m(7) * t821 + t832 * mrSges(7,1) - t833 * mrSges(7,2) + t875 * t854 - t876 * t855;
t811 = m(6) * t825 + t936 * mrSges(6,1) - t845 * mrSges(6,3) - t889 * t862 + t937 * t879 + t964;
t795 = t947 * t802 + t952 * t811;
t890 = -t903 * mrSges(5,1) + t904 * mrSges(5,2);
t893 = -t942 * mrSges(5,2) + t903 * mrSges(5,3);
t792 = m(5) * t837 + t939 * mrSges(5,1) - t874 * mrSges(5,3) - t904 * t890 + t942 * t893 + t795;
t894 = t942 * mrSges(5,1) - t904 * mrSges(5,3);
t975 = t952 * t802 - t947 * t811;
t793 = m(5) * t838 - t939 * mrSges(5,2) + t873 * mrSges(5,3) + t903 * t890 - t942 * t894 + t975;
t786 = t953 * t792 + t948 * t793;
t907 = -t920 * mrSges(4,1) + t921 * mrSges(4,2);
t914 = -qJD(3) * mrSges(4,2) + t920 * mrSges(4,3);
t784 = m(4) * t877 + qJDD(3) * mrSges(4,1) - t911 * mrSges(4,3) + qJD(3) * t914 - t921 * t907 + t786;
t915 = qJD(3) * mrSges(4,1) - t921 * mrSges(4,3);
t976 = -t948 * t792 + t953 * t793;
t785 = m(4) * t878 - qJDD(3) * mrSges(4,2) + t910 * mrSges(4,3) - qJD(3) * t915 + t920 * t907 + t976;
t779 = t954 * t784 + t949 * t785;
t912 = -t944 * t922 + t980;
t967 = mrSges(3,3) * qJDD(1) + t956 * (-mrSges(3,1) * t945 + t987);
t777 = m(3) * t912 - t944 * t967 + t779;
t977 = -t949 * t784 + t954 * t785;
t778 = m(3) * t913 + t945 * t967 + t977;
t978 = -t944 * t777 + t945 * t778;
t770 = m(2) * t927 - t956 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t978;
t918 = -qJDD(1) * pkin(1) - t956 * qJ(2) + t973;
t804 = t951 * t815 + t946 * t816;
t966 = m(6) * t835 - t844 * mrSges(6,1) + t845 * mrSges(6,2) - t888 * t879 + t889 * t880 + t804;
t963 = m(5) * t868 - t873 * mrSges(5,1) + t874 * mrSges(5,2) - t903 * t893 + t904 * t894 + t966;
t960 = m(4) * t909 - t910 * mrSges(4,1) + t911 * mrSges(4,2) - t920 * t914 + t921 * t915 + t963;
t958 = -m(3) * t918 + mrSges(3,1) * t981 - t960 + (t940 * t956 + t986) * mrSges(3,3);
t797 = t958 + (mrSges(2,1) - t987) * qJDD(1) - t956 * mrSges(2,2) + m(2) * t926;
t985 = t950 * t770 + t955 * t797;
t772 = t945 * t777 + t944 * t778;
t970 = Ifges(3,5) * t944 + Ifges(3,6) * t945;
t984 = t956 * t970;
t979 = t955 * t770 - t950 * t797;
t972 = Ifges(3,1) * t944 + Ifges(3,4) * t945;
t971 = Ifges(3,4) * t944 + Ifges(3,2) * t945;
t846 = Ifges(7,5) * t876 + Ifges(7,6) * t875 + Ifges(7,3) * t884;
t848 = Ifges(7,1) * t876 + Ifges(7,4) * t875 + Ifges(7,5) * t884;
t808 = -mrSges(7,1) * t821 + mrSges(7,3) * t820 + Ifges(7,4) * t833 + Ifges(7,2) * t832 + Ifges(7,6) * t843 - t876 * t846 + t884 * t848;
t847 = Ifges(7,4) * t876 + Ifges(7,2) * t875 + Ifges(7,6) * t884;
t809 = mrSges(7,2) * t821 - mrSges(7,3) * t819 + Ifges(7,1) * t833 + Ifges(7,4) * t832 + Ifges(7,5) * t843 + t875 * t846 - t884 * t847;
t857 = Ifges(6,5) * t889 + Ifges(6,6) * t888 + Ifges(6,3) * t937;
t858 = Ifges(6,4) * t889 + Ifges(6,2) * t888 + Ifges(6,6) * t937;
t787 = mrSges(6,2) * t835 - mrSges(6,3) * t825 + Ifges(6,1) * t845 + Ifges(6,4) * t844 + Ifges(6,5) * t936 - pkin(10) * t804 - t946 * t808 + t951 * t809 + t888 * t857 - t937 * t858;
t859 = Ifges(6,1) * t889 + Ifges(6,4) * t888 + Ifges(6,5) * t937;
t961 = mrSges(7,1) * t819 - mrSges(7,2) * t820 + Ifges(7,5) * t833 + Ifges(7,6) * t832 + Ifges(7,3) * t843 + t876 * t847 - t875 * t848;
t788 = -mrSges(6,1) * t835 + mrSges(6,3) * t826 + Ifges(6,4) * t845 + Ifges(6,2) * t844 + Ifges(6,6) * t936 - pkin(5) * t804 - t889 * t857 + t937 * t859 - t961;
t881 = Ifges(5,5) * t904 + Ifges(5,6) * t903 + Ifges(5,3) * t942;
t883 = Ifges(5,1) * t904 + Ifges(5,4) * t903 + Ifges(5,5) * t942;
t773 = -mrSges(5,1) * t868 + mrSges(5,3) * t838 + Ifges(5,4) * t874 + Ifges(5,2) * t873 + Ifges(5,6) * t939 - pkin(4) * t966 + pkin(9) * t975 + t947 * t787 + t952 * t788 - t904 * t881 + t942 * t883;
t882 = Ifges(5,4) * t904 + Ifges(5,2) * t903 + Ifges(5,6) * t942;
t780 = mrSges(5,2) * t868 - mrSges(5,3) * t837 + Ifges(5,1) * t874 + Ifges(5,4) * t873 + Ifges(5,5) * t939 - pkin(9) * t795 + t952 * t787 - t947 * t788 + t903 * t881 - t942 * t882;
t900 = Ifges(4,5) * t921 + Ifges(4,6) * t920 + Ifges(4,3) * qJD(3);
t902 = Ifges(4,1) * t921 + Ifges(4,4) * t920 + Ifges(4,5) * qJD(3);
t766 = -mrSges(4,1) * t909 + mrSges(4,3) * t878 + Ifges(4,4) * t911 + Ifges(4,2) * t910 + Ifges(4,6) * qJDD(3) - pkin(3) * t963 + pkin(8) * t976 + qJD(3) * t902 + t953 * t773 + t948 * t780 - t921 * t900;
t901 = Ifges(4,4) * t921 + Ifges(4,2) * t920 + Ifges(4,6) * qJD(3);
t767 = mrSges(4,2) * t909 - mrSges(4,3) * t877 + Ifges(4,1) * t911 + Ifges(4,4) * t910 + Ifges(4,5) * qJDD(3) - pkin(8) * t786 - qJD(3) * t901 - t948 * t773 + t953 * t780 + t920 * t900;
t762 = -mrSges(3,1) * t918 + mrSges(3,3) * t913 - pkin(2) * t960 + pkin(7) * t977 + qJDD(1) * t971 + t954 * t766 + t949 * t767 - t944 * t984;
t764 = mrSges(3,2) * t918 - mrSges(3,3) * t912 - pkin(7) * t779 + qJDD(1) * t972 - t949 * t766 + t954 * t767 + t945 * t984;
t799 = qJDD(1) * t987 - t958;
t965 = mrSges(2,1) * t926 - mrSges(2,2) * t927 + Ifges(2,3) * qJDD(1) - pkin(1) * t799 + qJ(2) * t978 + t945 * t762 + t944 * t764;
t962 = -mrSges(6,1) * t825 + mrSges(6,2) * t826 - Ifges(6,5) * t845 - Ifges(6,6) * t844 - Ifges(6,3) * t936 - pkin(5) * t964 - pkin(10) * t974 - t951 * t808 - t946 * t809 - t889 * t858 + t888 * t859;
t959 = -mrSges(5,1) * t837 + mrSges(5,2) * t838 - Ifges(5,5) * t874 - Ifges(5,6) * t873 - Ifges(5,3) * t939 - pkin(4) * t795 - t904 * t882 + t903 * t883 + t962;
t957 = mrSges(4,1) * t877 - mrSges(4,2) * t878 + Ifges(4,5) * t911 + Ifges(4,6) * t910 + Ifges(4,3) * qJDD(3) + pkin(3) * t786 + t921 * t901 - t920 * t902 - t959;
t765 = -t957 - pkin(2) * t779 + mrSges(2,1) * g(3) + (Ifges(2,6) - t970) * qJDD(1) + mrSges(2,3) * t927 - mrSges(3,1) * t912 + mrSges(3,2) * t913 - pkin(1) * t772 + (-t944 * t971 + t945 * t972 + Ifges(2,5)) * t956;
t760 = -mrSges(2,2) * g(3) - mrSges(2,3) * t926 + Ifges(2,5) * qJDD(1) - t956 * Ifges(2,6) - qJ(2) * t772 - t944 * t762 + t945 * t764;
t1 = [-m(1) * g(1) + t979; -m(1) * g(2) + t985; (-m(1) - m(2)) * g(3) + t772; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t985 + t955 * t760 - t950 * t765; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t979 + t950 * t760 + t955 * t765; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t965; t965; t799; t957; -t959; -t962; t961;];
tauJB  = t1;
