% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-05-07 07:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:25:40
% EndTime: 2019-05-07 07:25:57
% DurationCPUTime: 16.70s
% Computational Cost: add. (257603->366), mult. (578683->448), div. (0->0), fcn. (424186->10), ass. (0->145)
t996 = Ifges(6,1) + Ifges(7,1);
t989 = Ifges(6,4) + Ifges(7,4);
t988 = Ifges(6,5) + Ifges(7,5);
t995 = Ifges(6,2) + Ifges(7,2);
t987 = Ifges(6,6) + Ifges(7,6);
t994 = Ifges(6,3) + Ifges(7,3);
t951 = sin(qJ(2));
t955 = cos(qJ(2));
t976 = qJD(1) * qJD(2);
t931 = qJDD(1) * t951 + t955 * t976;
t952 = sin(qJ(1));
t956 = cos(qJ(1));
t938 = -g(1) * t956 - g(2) * t952;
t957 = qJD(1) ^ 2;
t926 = -pkin(1) * t957 + qJDD(1) * pkin(7) + t938;
t985 = t951 * t926;
t991 = pkin(2) * t957;
t889 = qJDD(2) * pkin(2) - t931 * pkin(8) - t985 + (pkin(8) * t976 + t951 * t991 - g(3)) * t955;
t914 = -g(3) * t951 + t955 * t926;
t932 = qJDD(1) * t955 - t951 * t976;
t979 = qJD(1) * t951;
t936 = qJD(2) * pkin(2) - pkin(8) * t979;
t946 = t955 ^ 2;
t890 = pkin(8) * t932 - qJD(2) * t936 - t946 * t991 + t914;
t950 = sin(qJ(3));
t954 = cos(qJ(3));
t861 = t954 * t889 - t950 * t890;
t923 = (-t950 * t951 + t954 * t955) * qJD(1);
t898 = qJD(3) * t923 + t931 * t954 + t932 * t950;
t924 = (t950 * t955 + t951 * t954) * qJD(1);
t943 = qJDD(2) + qJDD(3);
t944 = qJD(2) + qJD(3);
t837 = (t923 * t944 - t898) * qJ(4) + (t923 * t924 + t943) * pkin(3) + t861;
t862 = t950 * t889 + t954 * t890;
t897 = -qJD(3) * t924 - t931 * t950 + t932 * t954;
t916 = pkin(3) * t944 - qJ(4) * t924;
t919 = t923 ^ 2;
t840 = -pkin(3) * t919 + qJ(4) * t897 - t916 * t944 + t862;
t947 = sin(pkin(10));
t948 = cos(pkin(10));
t911 = t923 * t947 + t924 * t948;
t831 = -0.2e1 * qJD(4) * t911 + t837 * t948 - t947 * t840;
t910 = t923 * t948 - t924 * t947;
t832 = 0.2e1 * qJD(4) * t910 + t947 * t837 + t948 * t840;
t869 = t897 * t948 - t898 * t947;
t883 = -mrSges(5,1) * t910 + mrSges(5,2) * t911;
t901 = mrSges(5,1) * t944 - mrSges(5,3) * t911;
t884 = -pkin(4) * t910 - pkin(9) * t911;
t942 = t944 ^ 2;
t829 = -pkin(4) * t942 + pkin(9) * t943 + t884 * t910 + t832;
t937 = t952 * g(1) - t956 * g(2);
t965 = -qJDD(1) * pkin(1) - t937;
t899 = -t932 * pkin(2) + t936 * t979 + (-pkin(8) * t946 - pkin(7)) * t957 + t965;
t847 = -t897 * pkin(3) - t919 * qJ(4) + t924 * t916 + qJDD(4) + t899;
t870 = t897 * t947 + t898 * t948;
t835 = (-t910 * t944 - t870) * pkin(9) + (t911 * t944 - t869) * pkin(4) + t847;
t949 = sin(qJ(5));
t953 = cos(qJ(5));
t824 = -t949 * t829 + t953 * t835;
t895 = -t911 * t949 + t944 * t953;
t845 = qJD(5) * t895 + t870 * t953 + t943 * t949;
t868 = qJDD(5) - t869;
t896 = t911 * t953 + t944 * t949;
t871 = -mrSges(7,1) * t895 + mrSges(7,2) * t896;
t872 = -mrSges(6,1) * t895 + mrSges(6,2) * t896;
t904 = qJD(5) - t910;
t874 = -mrSges(6,2) * t904 + mrSges(6,3) * t895;
t821 = -0.2e1 * qJD(6) * t896 + (t895 * t904 - t845) * qJ(6) + (t895 * t896 + t868) * pkin(5) + t824;
t873 = -mrSges(7,2) * t904 + mrSges(7,3) * t895;
t975 = m(7) * t821 + t868 * mrSges(7,1) + t904 * t873;
t810 = m(6) * t824 + t868 * mrSges(6,1) + t904 * t874 + (-t871 - t872) * t896 + (-mrSges(6,3) - mrSges(7,3)) * t845 + t975;
t825 = t953 * t829 + t949 * t835;
t844 = -qJD(5) * t896 - t870 * t949 + t943 * t953;
t875 = pkin(5) * t904 - qJ(6) * t896;
t891 = t895 ^ 2;
t823 = -pkin(5) * t891 + qJ(6) * t844 + 0.2e1 * qJD(6) * t895 - t875 * t904 + t825;
t974 = m(7) * t823 + t844 * mrSges(7,3) + t895 * t871;
t876 = mrSges(7,1) * t904 - mrSges(7,3) * t896;
t980 = -mrSges(6,1) * t904 + mrSges(6,3) * t896 - t876;
t990 = -mrSges(6,2) - mrSges(7,2);
t813 = m(6) * t825 + t844 * mrSges(6,3) + t868 * t990 + t895 * t872 + t904 * t980 + t974;
t969 = -t810 * t949 + t953 * t813;
t802 = m(5) * t832 - mrSges(5,2) * t943 + mrSges(5,3) * t869 + t883 * t910 - t901 * t944 + t969;
t900 = -mrSges(5,2) * t944 + mrSges(5,3) * t910;
t828 = -pkin(4) * t943 - pkin(9) * t942 + t911 * t884 - t831;
t826 = -pkin(5) * t844 - qJ(6) * t891 + t875 * t896 + qJDD(6) + t828;
t967 = -m(7) * t826 + t844 * mrSges(7,1) + t895 * t873;
t961 = -m(6) * t828 + t844 * mrSges(6,1) + t845 * t990 + t895 * t874 + t896 * t980 + t967;
t815 = m(5) * t831 + t943 * mrSges(5,1) - t870 * mrSges(5,3) - t911 * t883 + t944 * t900 + t961;
t794 = t947 * t802 + t948 * t815;
t912 = -mrSges(4,1) * t923 + mrSges(4,2) * t924;
t915 = -mrSges(4,2) * t944 + mrSges(4,3) * t923;
t791 = m(4) * t861 + mrSges(4,1) * t943 - mrSges(4,3) * t898 - t912 * t924 + t915 * t944 + t794;
t917 = mrSges(4,1) * t944 - mrSges(4,3) * t924;
t970 = t948 * t802 - t815 * t947;
t792 = m(4) * t862 - mrSges(4,2) * t943 + mrSges(4,3) * t897 + t912 * t923 - t917 * t944 + t970;
t786 = t954 * t791 + t950 * t792;
t913 = -t955 * g(3) - t985;
t921 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t951 + Ifges(3,2) * t955) * qJD(1);
t922 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t951 + Ifges(3,4) * t955) * qJD(1);
t819 = t845 * mrSges(7,2) + t896 * t876 - t967;
t981 = t989 * t895 + t996 * t896 + t988 * t904;
t983 = -t987 * t895 - t988 * t896 - t994 * t904;
t796 = -mrSges(6,1) * t828 + mrSges(6,3) * t825 - mrSges(7,1) * t826 + mrSges(7,3) * t823 - pkin(5) * t819 + qJ(6) * t974 + (-qJ(6) * t876 + t981) * t904 + t983 * t896 + (-mrSges(7,2) * qJ(6) + t987) * t868 + t989 * t845 + t995 * t844;
t818 = -t845 * mrSges(7,3) - t896 * t871 + t975;
t982 = -t995 * t895 - t989 * t896 - t987 * t904;
t804 = mrSges(6,2) * t828 + mrSges(7,2) * t826 - mrSges(6,3) * t824 - mrSges(7,3) * t821 - qJ(6) * t818 + t989 * t844 + t996 * t845 + t988 * t868 - t983 * t895 + t982 * t904;
t879 = Ifges(5,4) * t911 + Ifges(5,2) * t910 + Ifges(5,6) * t944;
t880 = Ifges(5,1) * t911 + Ifges(5,4) * t910 + Ifges(5,5) * t944;
t906 = Ifges(4,4) * t924 + Ifges(4,2) * t923 + Ifges(4,6) * t944;
t907 = Ifges(4,1) * t924 + Ifges(4,4) * t923 + Ifges(4,5) * t944;
t960 = -mrSges(4,1) * t861 - mrSges(5,1) * t831 + mrSges(4,2) * t862 + mrSges(5,2) * t832 - pkin(3) * t794 - pkin(4) * t961 - pkin(9) * t969 - t953 * t796 - t949 * t804 - t911 * t879 + t910 * t880 + t923 * t907 - Ifges(5,6) * t869 - Ifges(5,5) * t870 - t924 * t906 - Ifges(4,6) * t897 - Ifges(4,5) * t898 + (-Ifges(4,3) - Ifges(5,3)) * t943;
t993 = mrSges(3,1) * t913 - mrSges(3,2) * t914 + Ifges(3,5) * t931 + Ifges(3,6) * t932 + Ifges(3,3) * qJDD(2) + pkin(2) * t786 + (t921 * t951 - t922 * t955) * qJD(1) - t960;
t992 = mrSges(6,1) * t824 + mrSges(7,1) * t821 - mrSges(6,2) * t825 - mrSges(7,2) * t823 + pkin(5) * t818 + t844 * t987 + t845 * t988 + t994 * t868 - t895 * t981 - t896 * t982;
t930 = (-mrSges(3,1) * t955 + mrSges(3,2) * t951) * qJD(1);
t978 = qJD(1) * t955;
t935 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t978;
t784 = m(3) * t913 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t931 + qJD(2) * t935 - t930 * t979 + t786;
t934 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t979;
t971 = -t791 * t950 + t954 * t792;
t785 = m(3) * t914 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t932 - qJD(2) * t934 + t930 * t978 + t971;
t972 = -t784 * t951 + t955 * t785;
t777 = m(2) * t938 - mrSges(2,1) * t957 - qJDD(1) * mrSges(2,2) + t972;
t925 = -t957 * pkin(7) + t965;
t807 = t953 * t810 + t949 * t813;
t805 = m(5) * t847 - t869 * mrSges(5,1) + t870 * mrSges(5,2) - t910 * t900 + t911 * t901 + t807;
t962 = m(4) * t899 - t897 * mrSges(4,1) + mrSges(4,2) * t898 - t923 * t915 + t917 * t924 + t805;
t959 = -m(3) * t925 + t932 * mrSges(3,1) - mrSges(3,2) * t931 - t934 * t979 + t935 * t978 - t962;
t798 = m(2) * t937 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t957 + t959;
t984 = t952 * t777 + t956 * t798;
t779 = t955 * t784 + t951 * t785;
t973 = t956 * t777 - t798 * t952;
t878 = Ifges(5,5) * t911 + Ifges(5,6) * t910 + Ifges(5,3) * t944;
t780 = mrSges(5,2) * t847 - mrSges(5,3) * t831 + Ifges(5,1) * t870 + Ifges(5,4) * t869 + Ifges(5,5) * t943 - pkin(9) * t807 - t796 * t949 + t804 * t953 + t878 * t910 - t879 * t944;
t787 = -mrSges(5,1) * t847 + mrSges(5,3) * t832 + Ifges(5,4) * t870 + Ifges(5,2) * t869 + Ifges(5,6) * t943 - pkin(4) * t807 - t911 * t878 + t944 * t880 - t992;
t905 = Ifges(4,5) * t924 + Ifges(4,6) * t923 + Ifges(4,3) * t944;
t773 = -mrSges(4,1) * t899 + mrSges(4,3) * t862 + Ifges(4,4) * t898 + Ifges(4,2) * t897 + Ifges(4,6) * t943 - pkin(3) * t805 + qJ(4) * t970 + t947 * t780 + t948 * t787 - t924 * t905 + t944 * t907;
t774 = mrSges(4,2) * t899 - mrSges(4,3) * t861 + Ifges(4,1) * t898 + Ifges(4,4) * t897 + Ifges(4,5) * t943 - qJ(4) * t794 + t780 * t948 - t787 * t947 + t905 * t923 - t906 * t944;
t920 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t951 + Ifges(3,6) * t955) * qJD(1);
t769 = -mrSges(3,1) * t925 + mrSges(3,3) * t914 + Ifges(3,4) * t931 + Ifges(3,2) * t932 + Ifges(3,6) * qJDD(2) - pkin(2) * t962 + pkin(8) * t971 + qJD(2) * t922 + t954 * t773 + t950 * t774 - t920 * t979;
t771 = mrSges(3,2) * t925 - mrSges(3,3) * t913 + Ifges(3,1) * t931 + Ifges(3,4) * t932 + Ifges(3,5) * qJDD(2) - pkin(8) * t786 - qJD(2) * t921 - t773 * t950 + t774 * t954 + t920 * t978;
t964 = mrSges(2,1) * t937 - mrSges(2,2) * t938 + Ifges(2,3) * qJDD(1) + pkin(1) * t959 + pkin(7) * t972 + t955 * t769 + t951 * t771;
t772 = mrSges(2,1) * g(3) + mrSges(2,3) * t938 + t957 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t779 - t993;
t767 = -mrSges(2,2) * g(3) - mrSges(2,3) * t937 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t957 - pkin(7) * t779 - t769 * t951 + t771 * t955;
t1 = [-m(1) * g(1) + t973; -m(1) * g(2) + t984; (-m(1) - m(2)) * g(3) + t779; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t984 + t956 * t767 - t952 * t772; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t973 + t952 * t767 + t956 * t772; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t964; t964; t993; -t960; t805; t992; t819;];
tauJB  = t1;
