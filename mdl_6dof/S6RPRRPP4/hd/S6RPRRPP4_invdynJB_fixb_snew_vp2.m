% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-05-05 21:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:31:02
% EndTime: 2019-05-05 21:31:15
% DurationCPUTime: 13.11s
% Computational Cost: add. (195167->341), mult. (467598->414), div. (0->0), fcn. (350622->10), ass. (0->146)
t963 = Ifges(6,1) + Ifges(7,1);
t956 = Ifges(6,4) - Ifges(7,5);
t955 = Ifges(6,5) + Ifges(7,4);
t962 = -Ifges(6,2) - Ifges(7,3);
t961 = -Ifges(7,2) - Ifges(6,3);
t954 = Ifges(6,6) - Ifges(7,6);
t919 = qJD(1) ^ 2;
t910 = sin(pkin(9));
t911 = cos(pkin(9));
t913 = sin(qJ(3));
t916 = cos(qJ(3));
t927 = t910 * t913 - t911 * t916;
t892 = t927 * qJD(1);
t928 = t910 * t916 + t911 * t913;
t893 = t928 * qJD(1);
t943 = t893 * qJD(3);
t878 = -t927 * qJDD(1) - t943;
t914 = sin(qJ(1));
t917 = cos(qJ(1));
t899 = -t917 * g(1) - t914 * g(2);
t894 = -t919 * pkin(1) + qJDD(1) * qJ(2) + t899;
t942 = qJD(1) * qJD(2);
t938 = -t911 * g(3) - 0.2e1 * t910 * t942;
t958 = pkin(2) * t911;
t865 = (-pkin(7) * qJDD(1) + t919 * t958 - t894) * t910 + t938;
t881 = -t910 * g(3) + (t894 + 0.2e1 * t942) * t911;
t941 = qJDD(1) * t911;
t907 = t911 ^ 2;
t951 = t907 * t919;
t866 = -pkin(2) * t951 + pkin(7) * t941 + t881;
t840 = t913 * t865 + t916 * t866;
t876 = t892 * pkin(3) - t893 * pkin(8);
t918 = qJD(3) ^ 2;
t816 = -t918 * pkin(3) + qJDD(3) * pkin(8) - t892 * t876 + t840;
t906 = t910 ^ 2;
t898 = t914 * g(1) - t917 * g(2);
t932 = qJDD(2) - t898;
t877 = (-pkin(1) - t958) * qJDD(1) + (-qJ(2) + (-t906 - t907) * pkin(7)) * t919 + t932;
t944 = t892 * qJD(3);
t879 = t928 * qJDD(1) - t944;
t824 = (-t879 + t944) * pkin(8) + (-t878 + t943) * pkin(3) + t877;
t912 = sin(qJ(4));
t915 = cos(qJ(4));
t812 = -t912 * t816 + t915 * t824;
t884 = t915 * qJD(3) - t912 * t893;
t853 = t884 * qJD(4) + t912 * qJDD(3) + t915 * t879;
t875 = qJDD(4) - t878;
t885 = t912 * qJD(3) + t915 * t893;
t890 = qJD(4) + t892;
t808 = (t884 * t890 - t853) * qJ(5) + (t884 * t885 + t875) * pkin(4) + t812;
t813 = t915 * t816 + t912 * t824;
t852 = -t885 * qJD(4) + t915 * qJDD(3) - t912 * t879;
t861 = t890 * pkin(4) - t885 * qJ(5);
t883 = t884 ^ 2;
t810 = -t883 * pkin(4) + t852 * qJ(5) - t890 * t861 + t813;
t909 = sin(pkin(10));
t952 = cos(pkin(10));
t855 = -t952 * t884 + t909 * t885;
t959 = -2 * qJD(5);
t804 = t909 * t808 + t952 * t810 + t855 * t959;
t820 = -t952 * t852 + t909 * t853;
t856 = t909 * t884 + t952 * t885;
t843 = t890 * mrSges(6,1) - t856 * mrSges(6,3);
t834 = t855 * pkin(5) - t856 * qJ(6);
t889 = t890 ^ 2;
t801 = -t889 * pkin(5) + t875 * qJ(6) + 0.2e1 * qJD(6) * t890 - t855 * t834 + t804;
t844 = -t890 * mrSges(7,1) + t856 * mrSges(7,2);
t939 = m(7) * t801 + t875 * mrSges(7,3) + t890 * t844;
t835 = t855 * mrSges(7,1) - t856 * mrSges(7,3);
t946 = -t855 * mrSges(6,1) - t856 * mrSges(6,2) - t835;
t957 = -mrSges(6,3) - mrSges(7,2);
t790 = m(6) * t804 - t875 * mrSges(6,2) + t957 * t820 - t890 * t843 + t946 * t855 + t939;
t925 = t952 * t808 - t909 * t810;
t803 = t856 * t959 + t925;
t821 = t909 * t852 + t952 * t853;
t842 = -t890 * mrSges(6,2) - t855 * mrSges(6,3);
t802 = -t875 * pkin(5) - t889 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t834) * t856 - t925;
t841 = -t855 * mrSges(7,2) + t890 * mrSges(7,3);
t933 = -m(7) * t802 + t875 * mrSges(7,1) + t890 * t841;
t792 = m(6) * t803 + t875 * mrSges(6,1) + t957 * t821 + t890 * t842 + t946 * t856 + t933;
t785 = t909 * t790 + t952 * t792;
t798 = t821 * mrSges(7,2) + t856 * t835 - t933;
t847 = Ifges(5,4) * t885 + Ifges(5,2) * t884 + Ifges(5,6) * t890;
t848 = Ifges(5,1) * t885 + Ifges(5,4) * t884 + Ifges(5,5) * t890;
t947 = -t956 * t855 + t963 * t856 + t955 * t890;
t948 = t962 * t855 + t956 * t856 + t954 * t890;
t960 = -t954 * t820 + t955 * t821 + t947 * t855 + t948 * t856 - (-Ifges(5,3) + t961) * t875 + mrSges(5,1) * t812 + mrSges(6,1) * t803 - mrSges(7,1) * t802 - mrSges(5,2) * t813 - mrSges(6,2) * t804 + mrSges(7,3) * t801 + Ifges(5,5) * t853 + Ifges(5,6) * t852 + pkin(4) * t785 - pkin(5) * t798 + qJ(6) * (-t820 * mrSges(7,2) - t855 * t835 + t939) + t885 * t847 - t884 * t848;
t953 = mrSges(3,2) * t910;
t857 = -t884 * mrSges(5,1) + t885 * mrSges(5,2);
t860 = -t890 * mrSges(5,2) + t884 * mrSges(5,3);
t783 = m(5) * t812 + t875 * mrSges(5,1) - t853 * mrSges(5,3) - t885 * t857 + t890 * t860 + t785;
t862 = t890 * mrSges(5,1) - t885 * mrSges(5,3);
t934 = t952 * t790 - t909 * t792;
t784 = m(5) * t813 - t875 * mrSges(5,2) + t852 * mrSges(5,3) + t884 * t857 - t890 * t862 + t934;
t779 = -t912 * t783 + t915 * t784;
t873 = t892 * mrSges(4,1) + t893 * mrSges(4,2);
t887 = qJD(3) * mrSges(4,1) - t893 * mrSges(4,3);
t777 = m(4) * t840 - qJDD(3) * mrSges(4,2) + t878 * mrSges(4,3) - qJD(3) * t887 - t892 * t873 + t779;
t839 = t916 * t865 - t913 * t866;
t815 = -qJDD(3) * pkin(3) - t918 * pkin(8) + t893 * t876 - t839;
t811 = -t852 * pkin(4) - t883 * qJ(5) + t885 * t861 + qJDD(5) + t815;
t806 = -0.2e1 * qJD(6) * t856 + (t855 * t890 - t821) * qJ(6) + (t856 * t890 + t820) * pkin(5) + t811;
t799 = m(7) * t806 + t820 * mrSges(7,1) - t821 * mrSges(7,3) + t855 * t841 - t856 * t844;
t796 = m(6) * t811 + t820 * mrSges(6,1) + t821 * mrSges(6,2) + t855 * t842 + t856 * t843 + t799;
t795 = -m(5) * t815 + t852 * mrSges(5,1) - t853 * mrSges(5,2) + t884 * t860 - t885 * t862 - t796;
t886 = -qJD(3) * mrSges(4,2) - t892 * mrSges(4,3);
t794 = m(4) * t839 + qJDD(3) * mrSges(4,1) - t879 * mrSges(4,3) + qJD(3) * t886 - t893 * t873 + t795;
t770 = t913 * t777 + t916 * t794;
t880 = -t910 * t894 + t938;
t926 = mrSges(3,3) * qJDD(1) + t919 * (-mrSges(3,1) * t911 + t953);
t768 = m(3) * t880 - t926 * t910 + t770;
t935 = t916 * t777 - t913 * t794;
t769 = m(3) * t881 + t926 * t911 + t935;
t936 = -t910 * t768 + t911 * t769;
t759 = m(2) * t899 - t919 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t936;
t891 = -qJDD(1) * pkin(1) - t919 * qJ(2) + t932;
t778 = t915 * t783 + t912 * t784;
t923 = m(4) * t877 - t878 * mrSges(4,1) + t879 * mrSges(4,2) + t892 * t886 + t893 * t887 + t778;
t922 = -m(3) * t891 + mrSges(3,1) * t941 - t923 + (t906 * t919 + t951) * mrSges(3,3);
t772 = t922 + (mrSges(2,1) - t953) * qJDD(1) + m(2) * t898 - t919 * mrSges(2,2);
t950 = t914 * t759 + t917 * t772;
t761 = t911 * t768 + t910 * t769;
t949 = t954 * t855 - t955 * t856 + t961 * t890;
t929 = Ifges(3,5) * t910 + Ifges(3,6) * t911;
t945 = t919 * t929;
t937 = t917 * t759 - t914 * t772;
t931 = Ifges(3,1) * t910 + Ifges(3,4) * t911;
t930 = Ifges(3,4) * t910 + Ifges(3,2) * t911;
t786 = -mrSges(6,1) * t811 - mrSges(7,1) * t806 + mrSges(7,2) * t801 + mrSges(6,3) * t804 - pkin(5) * t799 + t962 * t820 + t956 * t821 + t949 * t856 + t954 * t875 + t947 * t890;
t787 = mrSges(6,2) * t811 + mrSges(7,2) * t802 - mrSges(6,3) * t803 - mrSges(7,3) * t806 - qJ(6) * t799 - t956 * t820 + t963 * t821 + t949 * t855 + t955 * t875 - t948 * t890;
t846 = Ifges(5,5) * t885 + Ifges(5,6) * t884 + Ifges(5,3) * t890;
t763 = -mrSges(5,1) * t815 + mrSges(5,3) * t813 + Ifges(5,4) * t853 + Ifges(5,2) * t852 + Ifges(5,6) * t875 - pkin(4) * t796 + qJ(5) * t934 + t952 * t786 + t909 * t787 - t885 * t846 + t890 * t848;
t764 = mrSges(5,2) * t815 - mrSges(5,3) * t812 + Ifges(5,1) * t853 + Ifges(5,4) * t852 + Ifges(5,5) * t875 - qJ(5) * t785 - t909 * t786 + t952 * t787 + t884 * t846 - t890 * t847;
t867 = Ifges(4,5) * t893 - Ifges(4,6) * t892 + Ifges(4,3) * qJD(3);
t868 = Ifges(4,4) * t893 - Ifges(4,2) * t892 + Ifges(4,6) * qJD(3);
t756 = mrSges(4,2) * t877 - mrSges(4,3) * t839 + Ifges(4,1) * t879 + Ifges(4,4) * t878 + Ifges(4,5) * qJDD(3) - pkin(8) * t778 - qJD(3) * t868 - t912 * t763 + t915 * t764 - t892 * t867;
t869 = Ifges(4,1) * t893 - Ifges(4,4) * t892 + Ifges(4,5) * qJD(3);
t762 = -mrSges(4,1) * t877 + mrSges(4,3) * t840 + Ifges(4,4) * t879 + Ifges(4,2) * t878 + Ifges(4,6) * qJDD(3) - pkin(3) * t778 + qJD(3) * t869 - t893 * t867 - t960;
t752 = -mrSges(3,1) * t891 + mrSges(3,3) * t881 - pkin(2) * t923 + pkin(7) * t935 + t930 * qJDD(1) + t913 * t756 + t916 * t762 - t910 * t945;
t755 = mrSges(3,2) * t891 - mrSges(3,3) * t880 - pkin(7) * t770 + t931 * qJDD(1) + t916 * t756 - t913 * t762 + t911 * t945;
t774 = qJDD(1) * t953 - t922;
t924 = mrSges(2,1) * t898 - mrSges(2,2) * t899 + Ifges(2,3) * qJDD(1) - pkin(1) * t774 + qJ(2) * t936 + t911 * t752 + t910 * t755;
t921 = mrSges(4,1) * t839 - mrSges(4,2) * t840 + Ifges(4,5) * t879 + Ifges(4,6) * t878 + Ifges(4,3) * qJDD(3) + pkin(3) * t795 + pkin(8) * t779 + t915 * t763 + t912 * t764 + t893 * t868 + t892 * t869;
t753 = -t921 - pkin(2) * t770 + mrSges(2,1) * g(3) - pkin(1) * t761 - mrSges(3,1) * t880 + mrSges(3,2) * t881 + mrSges(2,3) * t899 + (Ifges(2,6) - t929) * qJDD(1) + (-t910 * t930 + t911 * t931 + Ifges(2,5)) * t919;
t750 = -mrSges(2,2) * g(3) - mrSges(2,3) * t898 + Ifges(2,5) * qJDD(1) - t919 * Ifges(2,6) - qJ(2) * t761 - t910 * t752 + t911 * t755;
t1 = [-m(1) * g(1) + t937; -m(1) * g(2) + t950; (-m(1) - m(2)) * g(3) + t761; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t950 + t917 * t750 - t914 * t753; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t937 + t914 * t750 + t917 * t753; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t924; t924; t774; t921; t960; t796; t798;];
tauJB  = t1;
