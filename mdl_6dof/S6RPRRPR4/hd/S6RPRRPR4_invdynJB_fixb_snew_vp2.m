% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:28:18
% EndTime: 2019-05-05 22:28:49
% DurationCPUTime: 30.69s
% Computational Cost: add. (476469->365), mult. (1160977->457), div. (0->0), fcn. (922148->12), ass. (0->153)
t944 = qJD(1) ^ 2;
t976 = cos(qJ(4));
t936 = cos(pkin(10));
t975 = pkin(2) * t936;
t934 = sin(pkin(10));
t974 = mrSges(3,2) * t934;
t930 = t936 ^ 2;
t973 = t930 * t944;
t940 = sin(qJ(1));
t943 = cos(qJ(1));
t919 = -g(1) * t943 - g(2) * t940;
t914 = -pkin(1) * t944 + qJDD(1) * qJ(2) + t919;
t969 = qJD(1) * qJD(2);
t967 = -t936 * g(3) - 0.2e1 * t934 * t969;
t887 = (-pkin(7) * qJDD(1) + t944 * t975 - t914) * t934 + t967;
t904 = -g(3) * t934 + (t914 + 0.2e1 * t969) * t936;
t968 = qJDD(1) * t936;
t888 = -pkin(2) * t973 + pkin(7) * t968 + t904;
t939 = sin(qJ(3));
t942 = cos(qJ(3));
t866 = t942 * t887 - t939 * t888;
t955 = t934 * t942 + t936 * t939;
t954 = -t934 * t939 + t936 * t942;
t912 = t954 * qJD(1);
t970 = t912 * qJD(3);
t902 = t955 * qJDD(1) + t970;
t913 = t955 * qJD(1);
t839 = (-t902 + t970) * pkin(8) + (t912 * t913 + qJDD(3)) * pkin(3) + t866;
t867 = t939 * t887 + t942 * t888;
t901 = -t913 * qJD(3) + t954 * qJDD(1);
t907 = qJD(3) * pkin(3) - pkin(8) * t913;
t911 = t912 ^ 2;
t842 = -pkin(3) * t911 + pkin(8) * t901 - qJD(3) * t907 + t867;
t938 = sin(qJ(4));
t832 = t938 * t839 + t976 * t842;
t894 = t938 * t912 + t976 * t913;
t861 = qJD(4) * t894 - t976 * t901 + t902 * t938;
t893 = -t976 * t912 + t913 * t938;
t876 = mrSges(5,1) * t893 + mrSges(5,2) * t894;
t931 = qJD(3) + qJD(4);
t885 = mrSges(5,1) * t931 - mrSges(5,3) * t894;
t928 = qJDD(3) + qJDD(4);
t875 = pkin(4) * t893 - qJ(5) * t894;
t927 = t931 ^ 2;
t823 = -pkin(4) * t927 + qJ(5) * t928 - t875 * t893 + t832;
t929 = t934 ^ 2;
t918 = t940 * g(1) - t943 * g(2);
t960 = qJDD(2) - t918;
t900 = (-pkin(1) - t975) * qJDD(1) + (-qJ(2) + (-t929 - t930) * pkin(7)) * t944 + t960;
t854 = -t901 * pkin(3) - t911 * pkin(8) + t913 * t907 + t900;
t862 = -t893 * qJD(4) + t938 * t901 + t976 * t902;
t826 = (t893 * t931 - t862) * qJ(5) + (t894 * t931 + t861) * pkin(4) + t854;
t933 = sin(pkin(11));
t935 = cos(pkin(11));
t881 = t894 * t935 + t931 * t933;
t818 = -0.2e1 * qJD(5) * t881 - t933 * t823 + t935 * t826;
t851 = t862 * t935 + t928 * t933;
t880 = -t894 * t933 + t931 * t935;
t816 = (t880 * t893 - t851) * pkin(9) + (t880 * t881 + t861) * pkin(5) + t818;
t819 = 0.2e1 * qJD(5) * t880 + t935 * t823 + t933 * t826;
t850 = -t862 * t933 + t928 * t935;
t869 = pkin(5) * t893 - pkin(9) * t881;
t879 = t880 ^ 2;
t817 = -pkin(5) * t879 + pkin(9) * t850 - t869 * t893 + t819;
t937 = sin(qJ(6));
t941 = cos(qJ(6));
t814 = t816 * t941 - t817 * t937;
t863 = t880 * t941 - t881 * t937;
t829 = qJD(6) * t863 + t850 * t937 + t851 * t941;
t864 = t880 * t937 + t881 * t941;
t837 = -mrSges(7,1) * t863 + mrSges(7,2) * t864;
t889 = qJD(6) + t893;
t843 = -mrSges(7,2) * t889 + mrSges(7,3) * t863;
t860 = qJDD(6) + t861;
t810 = m(7) * t814 + mrSges(7,1) * t860 - mrSges(7,3) * t829 - t837 * t864 + t843 * t889;
t815 = t816 * t937 + t817 * t941;
t828 = -qJD(6) * t864 + t850 * t941 - t851 * t937;
t844 = mrSges(7,1) * t889 - mrSges(7,3) * t864;
t811 = m(7) * t815 - mrSges(7,2) * t860 + mrSges(7,3) * t828 + t837 * t863 - t844 * t889;
t802 = t941 * t810 + t937 * t811;
t865 = -mrSges(6,1) * t880 + mrSges(6,2) * t881;
t959 = -mrSges(6,2) * t893 + mrSges(6,3) * t880;
t800 = m(6) * t818 + t861 * mrSges(6,1) - t851 * mrSges(6,3) - t881 * t865 + t893 * t959 + t802;
t868 = mrSges(6,1) * t893 - mrSges(6,3) * t881;
t961 = -t810 * t937 + t941 * t811;
t801 = m(6) * t819 - mrSges(6,2) * t861 + mrSges(6,3) * t850 + t865 * t880 - t868 * t893 + t961;
t962 = -t800 * t933 + t935 * t801;
t793 = m(5) * t832 - mrSges(5,2) * t928 - mrSges(5,3) * t861 - t876 * t893 - t885 * t931 + t962;
t831 = t976 * t839 - t938 * t842;
t822 = -t928 * pkin(4) - t927 * qJ(5) + t894 * t875 + qJDD(5) - t831;
t820 = -t850 * pkin(5) - t879 * pkin(9) + t881 * t869 + t822;
t950 = m(7) * t820 - t828 * mrSges(7,1) + mrSges(7,2) * t829 - t863 * t843 + t844 * t864;
t813 = m(6) * t822 - t850 * mrSges(6,1) + mrSges(6,2) * t851 + t868 * t881 - t880 * t959 + t950;
t884 = -mrSges(5,2) * t931 - mrSges(5,3) * t893;
t806 = m(5) * t831 + mrSges(5,1) * t928 - mrSges(5,3) * t862 - t876 * t894 + t884 * t931 - t813;
t782 = t938 * t793 + t976 * t806;
t898 = -mrSges(4,1) * t912 + mrSges(4,2) * t913;
t905 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t912;
t780 = m(4) * t866 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t902 + qJD(3) * t905 - t898 * t913 + t782;
t906 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t913;
t963 = t976 * t793 - t806 * t938;
t781 = m(4) * t867 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t901 - qJD(3) * t906 + t898 * t912 + t963;
t775 = t942 * t780 + t939 * t781;
t903 = -t934 * t914 + t967;
t953 = mrSges(3,3) * qJDD(1) + t944 * (-mrSges(3,1) * t936 + t974);
t773 = m(3) * t903 - t953 * t934 + t775;
t964 = -t939 * t780 + t942 * t781;
t774 = m(3) * t904 + t953 * t936 + t964;
t965 = -t773 * t934 + t936 * t774;
t766 = m(2) * t919 - mrSges(2,1) * t944 - qJDD(1) * mrSges(2,2) + t965;
t910 = -qJDD(1) * pkin(1) - t944 * qJ(2) + t960;
t795 = t935 * t800 + t933 * t801;
t952 = m(5) * t854 + t861 * mrSges(5,1) + t862 * mrSges(5,2) + t893 * t884 + t894 * t885 + t795;
t948 = m(4) * t900 - t901 * mrSges(4,1) + t902 * mrSges(4,2) - t912 * t905 + t913 * t906 + t952;
t946 = -m(3) * t910 + mrSges(3,1) * t968 - t948 + (t929 * t944 + t973) * mrSges(3,3);
t788 = t946 + (mrSges(2,1) - t974) * qJDD(1) - t944 * mrSges(2,2) + m(2) * t918;
t972 = t940 * t766 + t943 * t788;
t768 = t936 * t773 + t934 * t774;
t956 = Ifges(3,5) * t934 + Ifges(3,6) * t936;
t971 = t944 * t956;
t966 = t943 * t766 - t788 * t940;
t958 = Ifges(3,1) * t934 + Ifges(3,4) * t936;
t957 = Ifges(3,4) * t934 + Ifges(3,2) * t936;
t833 = Ifges(7,5) * t864 + Ifges(7,6) * t863 + Ifges(7,3) * t889;
t835 = Ifges(7,1) * t864 + Ifges(7,4) * t863 + Ifges(7,5) * t889;
t803 = -mrSges(7,1) * t820 + mrSges(7,3) * t815 + Ifges(7,4) * t829 + Ifges(7,2) * t828 + Ifges(7,6) * t860 - t833 * t864 + t835 * t889;
t834 = Ifges(7,4) * t864 + Ifges(7,2) * t863 + Ifges(7,6) * t889;
t804 = mrSges(7,2) * t820 - mrSges(7,3) * t814 + Ifges(7,1) * t829 + Ifges(7,4) * t828 + Ifges(7,5) * t860 + t833 * t863 - t834 * t889;
t845 = Ifges(6,5) * t881 + Ifges(6,6) * t880 + Ifges(6,3) * t893;
t847 = Ifges(6,1) * t881 + Ifges(6,4) * t880 + Ifges(6,5) * t893;
t784 = -mrSges(6,1) * t822 + mrSges(6,3) * t819 + Ifges(6,4) * t851 + Ifges(6,2) * t850 + Ifges(6,6) * t861 - pkin(5) * t950 + pkin(9) * t961 + t941 * t803 + t937 * t804 - t881 * t845 + t893 * t847;
t846 = Ifges(6,4) * t881 + Ifges(6,2) * t880 + Ifges(6,6) * t893;
t786 = mrSges(6,2) * t822 - mrSges(6,3) * t818 + Ifges(6,1) * t851 + Ifges(6,4) * t850 + Ifges(6,5) * t861 - pkin(9) * t802 - t803 * t937 + t804 * t941 + t845 * t880 - t846 * t893;
t870 = Ifges(5,5) * t894 - Ifges(5,6) * t893 + Ifges(5,3) * t931;
t871 = Ifges(5,4) * t894 - Ifges(5,2) * t893 + Ifges(5,6) * t931;
t769 = mrSges(5,2) * t854 - mrSges(5,3) * t831 + Ifges(5,1) * t862 - Ifges(5,4) * t861 + Ifges(5,5) * t928 - qJ(5) * t795 - t784 * t933 + t786 * t935 - t870 * t893 - t871 * t931;
t872 = Ifges(5,1) * t894 - Ifges(5,4) * t893 + Ifges(5,5) * t931;
t947 = mrSges(7,1) * t814 - mrSges(7,2) * t815 + Ifges(7,5) * t829 + Ifges(7,6) * t828 + Ifges(7,3) * t860 + t864 * t834 - t863 * t835;
t776 = -t947 + (-Ifges(5,2) - Ifges(6,3)) * t861 + t931 * t872 + Ifges(5,6) * t928 - t894 * t870 - t881 * t846 + t880 * t847 + Ifges(5,4) * t862 - Ifges(6,6) * t850 - Ifges(6,5) * t851 - mrSges(5,1) * t854 + mrSges(5,3) * t832 + mrSges(6,2) * t819 - mrSges(6,1) * t818 - pkin(5) * t802 - pkin(4) * t795;
t890 = Ifges(4,5) * t913 + Ifges(4,6) * t912 + Ifges(4,3) * qJD(3);
t892 = Ifges(4,1) * t913 + Ifges(4,4) * t912 + Ifges(4,5) * qJD(3);
t762 = -mrSges(4,1) * t900 + mrSges(4,3) * t867 + Ifges(4,4) * t902 + Ifges(4,2) * t901 + Ifges(4,6) * qJDD(3) - pkin(3) * t952 + pkin(8) * t963 + qJD(3) * t892 + t938 * t769 + t976 * t776 - t913 * t890;
t891 = Ifges(4,4) * t913 + Ifges(4,2) * t912 + Ifges(4,6) * qJD(3);
t763 = mrSges(4,2) * t900 - mrSges(4,3) * t866 + Ifges(4,1) * t902 + Ifges(4,4) * t901 + Ifges(4,5) * qJDD(3) - pkin(8) * t782 - qJD(3) * t891 + t976 * t769 - t938 * t776 + t912 * t890;
t758 = -mrSges(3,1) * t910 + mrSges(3,3) * t904 - pkin(2) * t948 + pkin(7) * t964 + t957 * qJDD(1) + t942 * t762 + t939 * t763 - t934 * t971;
t760 = mrSges(3,2) * t910 - mrSges(3,3) * t903 - pkin(7) * t775 + t958 * qJDD(1) - t939 * t762 + t942 * t763 + t936 * t971;
t790 = qJDD(1) * t974 - t946;
t951 = mrSges(2,1) * t918 - mrSges(2,2) * t919 + Ifges(2,3) * qJDD(1) - pkin(1) * t790 + qJ(2) * t965 + t936 * t758 + t934 * t760;
t949 = -mrSges(5,1) * t831 + mrSges(5,2) * t832 - Ifges(5,5) * t862 + Ifges(5,6) * t861 - Ifges(5,3) * t928 + pkin(4) * t813 - qJ(5) * t962 - t935 * t784 - t933 * t786 - t894 * t871 - t893 * t872;
t945 = mrSges(4,1) * t866 - mrSges(4,2) * t867 + Ifges(4,5) * t902 + Ifges(4,6) * t901 + Ifges(4,3) * qJDD(3) + pkin(3) * t782 + t913 * t891 - t912 * t892 - t949;
t761 = mrSges(2,1) * g(3) + (Ifges(2,6) - t956) * qJDD(1) - t945 - pkin(2) * t775 + mrSges(2,3) * t919 + mrSges(3,2) * t904 - mrSges(3,1) * t903 - pkin(1) * t768 + (-t934 * t957 + t936 * t958 + Ifges(2,5)) * t944;
t756 = -mrSges(2,2) * g(3) - mrSges(2,3) * t918 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t944 - qJ(2) * t768 - t758 * t934 + t760 * t936;
t1 = [-m(1) * g(1) + t966; -m(1) * g(2) + t972; (-m(1) - m(2)) * g(3) + t768; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t972 + t943 * t756 - t940 * t761; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t966 + t940 * t756 + t943 * t761; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t951; t951; t790; t945; -t949; t813; t947;];
tauJB  = t1;
