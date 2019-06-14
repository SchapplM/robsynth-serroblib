% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR4
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
% Datum: 2019-05-07 10:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:26:45
% EndTime: 2019-05-07 10:27:22
% DurationCPUTime: 37.12s
% Computational Cost: add. (655864->389), mult. (1352907->489), div. (0->0), fcn. (998222->12), ass. (0->154)
t946 = sin(qJ(2));
t950 = cos(qJ(2));
t971 = qJD(1) * qJD(2);
t925 = t946 * qJDD(1) + t950 * t971;
t947 = sin(qJ(1));
t951 = cos(qJ(1));
t932 = -t951 * g(1) - t947 * g(2);
t952 = qJD(1) ^ 2;
t919 = -t952 * pkin(1) + qJDD(1) * pkin(7) + t932;
t975 = t946 * t919;
t976 = pkin(2) * t952;
t881 = qJDD(2) * pkin(2) - t925 * pkin(8) - t975 + (pkin(8) * t971 + t946 * t976 - g(3)) * t950;
t907 = -t946 * g(3) + t950 * t919;
t926 = t950 * qJDD(1) - t946 * t971;
t973 = qJD(1) * t946;
t930 = qJD(2) * pkin(2) - pkin(8) * t973;
t940 = t950 ^ 2;
t882 = t926 * pkin(8) - qJD(2) * t930 - t940 * t976 + t907;
t945 = sin(qJ(3));
t977 = cos(qJ(3));
t861 = t945 * t881 + t977 * t882;
t917 = (t945 * t950 + t977 * t946) * qJD(1);
t890 = t917 * qJD(3) + t945 * t925 - t977 * t926;
t972 = qJD(1) * t950;
t916 = t945 * t973 - t977 * t972;
t900 = t916 * mrSges(4,1) + t917 * mrSges(4,2);
t938 = qJD(2) + qJD(3);
t909 = t938 * mrSges(4,1) - t917 * mrSges(4,3);
t937 = qJDD(2) + qJDD(3);
t891 = -t916 * qJD(3) + t977 * t925 + t945 * t926;
t931 = t947 * g(1) - t951 * g(2);
t962 = -qJDD(1) * pkin(1) - t931;
t892 = -t926 * pkin(2) + t930 * t973 + (-pkin(8) * t940 - pkin(7)) * t952 + t962;
t846 = (t916 * t938 - t891) * qJ(4) + (t917 * t938 + t890) * pkin(3) + t892;
t899 = t916 * pkin(3) - t917 * qJ(4);
t936 = t938 ^ 2;
t849 = -t936 * pkin(3) + t937 * qJ(4) - t916 * t899 + t861;
t941 = sin(pkin(11));
t942 = cos(pkin(11));
t905 = t942 * t917 + t941 * t938;
t833 = -0.2e1 * qJD(4) * t905 + t942 * t846 - t941 * t849;
t875 = t942 * t891 + t941 * t937;
t904 = -t941 * t917 + t942 * t938;
t824 = (t916 * t904 - t875) * pkin(9) + (t904 * t905 + t890) * pkin(4) + t833;
t834 = 0.2e1 * qJD(4) * t904 + t941 * t846 + t942 * t849;
t874 = -t941 * t891 + t942 * t937;
t894 = t916 * pkin(4) - t905 * pkin(9);
t903 = t904 ^ 2;
t826 = -t903 * pkin(4) + t874 * pkin(9) - t916 * t894 + t834;
t944 = sin(qJ(5));
t949 = cos(qJ(5));
t819 = t949 * t824 - t944 * t826;
t872 = t949 * t904 - t944 * t905;
t845 = t872 * qJD(5) + t944 * t874 + t949 * t875;
t873 = t944 * t904 + t949 * t905;
t889 = qJDD(5) + t890;
t912 = qJD(5) + t916;
t817 = (t872 * t912 - t845) * pkin(10) + (t872 * t873 + t889) * pkin(5) + t819;
t820 = t944 * t824 + t949 * t826;
t844 = -t873 * qJD(5) + t949 * t874 - t944 * t875;
t864 = t912 * pkin(5) - t873 * pkin(10);
t871 = t872 ^ 2;
t818 = -t871 * pkin(5) + t844 * pkin(10) - t912 * t864 + t820;
t943 = sin(qJ(6));
t948 = cos(qJ(6));
t815 = t948 * t817 - t943 * t818;
t856 = t948 * t872 - t943 * t873;
t831 = t856 * qJD(6) + t943 * t844 + t948 * t845;
t857 = t943 * t872 + t948 * t873;
t840 = -t856 * mrSges(7,1) + t857 * mrSges(7,2);
t911 = qJD(6) + t912;
t850 = -t911 * mrSges(7,2) + t856 * mrSges(7,3);
t884 = qJDD(6) + t889;
t808 = m(7) * t815 + t884 * mrSges(7,1) - t831 * mrSges(7,3) - t857 * t840 + t911 * t850;
t816 = t943 * t817 + t948 * t818;
t830 = -t857 * qJD(6) + t948 * t844 - t943 * t845;
t851 = t911 * mrSges(7,1) - t857 * mrSges(7,3);
t809 = m(7) * t816 - t884 * mrSges(7,2) + t830 * mrSges(7,3) + t856 * t840 - t911 * t851;
t802 = t948 * t808 + t943 * t809;
t858 = -t872 * mrSges(6,1) + t873 * mrSges(6,2);
t862 = -t912 * mrSges(6,2) + t872 * mrSges(6,3);
t800 = m(6) * t819 + t889 * mrSges(6,1) - t845 * mrSges(6,3) - t873 * t858 + t912 * t862 + t802;
t863 = t912 * mrSges(6,1) - t873 * mrSges(6,3);
t965 = -t943 * t808 + t948 * t809;
t801 = m(6) * t820 - t889 * mrSges(6,2) + t844 * mrSges(6,3) + t872 * t858 - t912 * t863 + t965;
t796 = t949 * t800 + t944 * t801;
t877 = -t904 * mrSges(5,1) + t905 * mrSges(5,2);
t964 = -t916 * mrSges(5,2) + t904 * mrSges(5,3);
t794 = m(5) * t833 + t890 * mrSges(5,1) - t875 * mrSges(5,3) - t905 * t877 + t916 * t964 + t796;
t893 = t916 * mrSges(5,1) - t905 * mrSges(5,3);
t966 = -t944 * t800 + t949 * t801;
t795 = m(5) * t834 - t890 * mrSges(5,2) + t874 * mrSges(5,3) + t904 * t877 - t916 * t893 + t966;
t967 = -t941 * t794 + t942 * t795;
t785 = m(4) * t861 - t937 * mrSges(4,2) - t890 * mrSges(4,3) - t916 * t900 - t938 * t909 + t967;
t860 = t977 * t881 - t945 * t882;
t848 = -t937 * pkin(3) - t936 * qJ(4) + t917 * t899 + qJDD(4) - t860;
t835 = -t874 * pkin(4) - t903 * pkin(9) + t905 * t894 + t848;
t821 = -t844 * pkin(5) - t871 * pkin(10) + t873 * t864 + t835;
t960 = m(7) * t821 - t830 * mrSges(7,1) + t831 * mrSges(7,2) - t856 * t850 + t857 * t851;
t956 = m(6) * t835 - t844 * mrSges(6,1) + t845 * mrSges(6,2) - t872 * t862 + t873 * t863 + t960;
t813 = m(5) * t848 - t874 * mrSges(5,1) + t875 * mrSges(5,2) + t905 * t893 - t904 * t964 + t956;
t908 = -t938 * mrSges(4,2) - t916 * mrSges(4,3);
t811 = m(4) * t860 + t937 * mrSges(4,1) - t891 * mrSges(4,3) - t917 * t900 + t938 * t908 - t813;
t779 = t945 * t785 + t977 * t811;
t906 = -t950 * g(3) - t975;
t914 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t946 + Ifges(3,2) * t950) * qJD(1);
t915 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t946 + Ifges(3,4) * t950) * qJD(1);
t836 = Ifges(7,5) * t857 + Ifges(7,6) * t856 + Ifges(7,3) * t911;
t838 = Ifges(7,1) * t857 + Ifges(7,4) * t856 + Ifges(7,5) * t911;
t803 = -mrSges(7,1) * t821 + mrSges(7,3) * t816 + Ifges(7,4) * t831 + Ifges(7,2) * t830 + Ifges(7,6) * t884 - t857 * t836 + t911 * t838;
t837 = Ifges(7,4) * t857 + Ifges(7,2) * t856 + Ifges(7,6) * t911;
t804 = mrSges(7,2) * t821 - mrSges(7,3) * t815 + Ifges(7,1) * t831 + Ifges(7,4) * t830 + Ifges(7,5) * t884 + t856 * t836 - t911 * t837;
t852 = Ifges(6,5) * t873 + Ifges(6,6) * t872 + Ifges(6,3) * t912;
t854 = Ifges(6,1) * t873 + Ifges(6,4) * t872 + Ifges(6,5) * t912;
t789 = -mrSges(6,1) * t835 + mrSges(6,3) * t820 + Ifges(6,4) * t845 + Ifges(6,2) * t844 + Ifges(6,6) * t889 - pkin(5) * t960 + pkin(10) * t965 + t948 * t803 + t943 * t804 - t873 * t852 + t912 * t854;
t853 = Ifges(6,4) * t873 + Ifges(6,2) * t872 + Ifges(6,6) * t912;
t790 = mrSges(6,2) * t835 - mrSges(6,3) * t819 + Ifges(6,1) * t845 + Ifges(6,4) * t844 + Ifges(6,5) * t889 - pkin(10) * t802 - t943 * t803 + t948 * t804 + t872 * t852 - t912 * t853;
t865 = Ifges(5,5) * t905 + Ifges(5,6) * t904 + Ifges(5,3) * t916;
t867 = Ifges(5,1) * t905 + Ifges(5,4) * t904 + Ifges(5,5) * t916;
t771 = -mrSges(5,1) * t848 + mrSges(5,3) * t834 + Ifges(5,4) * t875 + Ifges(5,2) * t874 + Ifges(5,6) * t890 - pkin(4) * t956 + pkin(9) * t966 + t949 * t789 + t944 * t790 - t905 * t865 + t916 * t867;
t866 = Ifges(5,4) * t905 + Ifges(5,2) * t904 + Ifges(5,6) * t916;
t773 = mrSges(5,2) * t848 - mrSges(5,3) * t833 + Ifges(5,1) * t875 + Ifges(5,4) * t874 + Ifges(5,5) * t890 - pkin(9) * t796 - t944 * t789 + t949 * t790 + t904 * t865 - t916 * t866;
t896 = Ifges(4,4) * t917 - Ifges(4,2) * t916 + Ifges(4,6) * t938;
t897 = Ifges(4,1) * t917 - Ifges(4,4) * t916 + Ifges(4,5) * t938;
t957 = -mrSges(4,1) * t860 + mrSges(4,2) * t861 - Ifges(4,5) * t891 + Ifges(4,6) * t890 - Ifges(4,3) * t937 + pkin(3) * t813 - qJ(4) * t967 - t942 * t771 - t941 * t773 - t917 * t896 - t916 * t897;
t978 = mrSges(3,1) * t906 - mrSges(3,2) * t907 + Ifges(3,5) * t925 + Ifges(3,6) * t926 + Ifges(3,3) * qJDD(2) + pkin(2) * t779 + (t946 * t914 - t950 * t915) * qJD(1) - t957;
t924 = (-mrSges(3,1) * t950 + mrSges(3,2) * t946) * qJD(1);
t929 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t972;
t777 = m(3) * t906 + qJDD(2) * mrSges(3,1) - t925 * mrSges(3,3) + qJD(2) * t929 - t924 * t973 + t779;
t928 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t973;
t968 = t977 * t785 - t945 * t811;
t778 = m(3) * t907 - qJDD(2) * mrSges(3,2) + t926 * mrSges(3,3) - qJD(2) * t928 + t924 * t972 + t968;
t969 = -t946 * t777 + t950 * t778;
t766 = m(2) * t932 - t952 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t969;
t918 = -t952 * pkin(7) + t962;
t787 = t942 * t794 + t941 * t795;
t958 = m(4) * t892 + t890 * mrSges(4,1) + t891 * mrSges(4,2) + t916 * t908 + t917 * t909 + t787;
t955 = -m(3) * t918 + t926 * mrSges(3,1) - t925 * mrSges(3,2) - t928 * t973 + t929 * t972 - t958;
t781 = m(2) * t931 + qJDD(1) * mrSges(2,1) - t952 * mrSges(2,2) + t955;
t974 = t947 * t766 + t951 * t781;
t768 = t950 * t777 + t946 * t778;
t970 = t951 * t766 - t947 * t781;
t895 = Ifges(4,5) * t917 - Ifges(4,6) * t916 + Ifges(4,3) * t938;
t763 = mrSges(4,2) * t892 - mrSges(4,3) * t860 + Ifges(4,1) * t891 - Ifges(4,4) * t890 + Ifges(4,5) * t937 - qJ(4) * t787 - t941 * t771 + t942 * t773 - t916 * t895 - t938 * t896;
t959 = -mrSges(7,1) * t815 + mrSges(7,2) * t816 - Ifges(7,5) * t831 - Ifges(7,6) * t830 - Ifges(7,3) * t884 - t857 * t837 + t856 * t838;
t953 = mrSges(6,1) * t819 - mrSges(6,2) * t820 + Ifges(6,5) * t845 + Ifges(6,6) * t844 + Ifges(6,3) * t889 + pkin(5) * t802 + t873 * t853 - t872 * t854 - t959;
t769 = (-Ifges(5,3) - Ifges(4,2)) * t890 + Ifges(4,6) * t937 + t938 * t897 - t917 * t895 - t905 * t866 - mrSges(4,1) * t892 + t904 * t867 + Ifges(4,4) * t891 - Ifges(5,6) * t874 - Ifges(5,5) * t875 + mrSges(4,3) * t861 - mrSges(5,1) * t833 + mrSges(5,2) * t834 - pkin(4) * t796 - t953 - pkin(3) * t787;
t913 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t946 + Ifges(3,6) * t950) * qJD(1);
t759 = -mrSges(3,1) * t918 + mrSges(3,3) * t907 + Ifges(3,4) * t925 + Ifges(3,2) * t926 + Ifges(3,6) * qJDD(2) - pkin(2) * t958 + pkin(8) * t968 + qJD(2) * t915 + t945 * t763 + t977 * t769 - t913 * t973;
t762 = mrSges(3,2) * t918 - mrSges(3,3) * t906 + Ifges(3,1) * t925 + Ifges(3,4) * t926 + Ifges(3,5) * qJDD(2) - pkin(8) * t779 - qJD(2) * t914 + t977 * t763 - t945 * t769 + t913 * t972;
t961 = mrSges(2,1) * t931 - mrSges(2,2) * t932 + Ifges(2,3) * qJDD(1) + pkin(1) * t955 + pkin(7) * t969 + t950 * t759 + t946 * t762;
t760 = mrSges(2,1) * g(3) + mrSges(2,3) * t932 + t952 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t768 - t978;
t757 = -mrSges(2,2) * g(3) - mrSges(2,3) * t931 + Ifges(2,5) * qJDD(1) - t952 * Ifges(2,6) - pkin(7) * t768 - t946 * t759 + t950 * t762;
t1 = [-m(1) * g(1) + t970; -m(1) * g(2) + t974; (-m(1) - m(2)) * g(3) + t768; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t974 + t951 * t757 - t947 * t760; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t970 + t947 * t757 + t951 * t760; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t961; t961; t978; -t957; t813; t953; -t959;];
tauJB  = t1;
