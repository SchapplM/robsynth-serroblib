% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRRPR6
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
% Datum: 2019-05-05 22:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:49:49
% EndTime: 2019-05-05 22:50:21
% DurationCPUTime: 31.62s
% Computational Cost: add. (491843->366), mult. (1186529->459), div. (0->0), fcn. (914610->12), ass. (0->156)
t950 = qJD(1) ^ 2;
t944 = sin(qJ(1));
t948 = cos(qJ(1));
t925 = -g(1) * t948 - g(2) * t944;
t920 = -pkin(1) * t950 + qJDD(1) * qJ(2) + t925;
t938 = sin(pkin(10));
t940 = cos(pkin(10));
t972 = qJD(1) * qJD(2);
t969 = -t940 * g(3) - 0.2e1 * t938 * t972;
t980 = pkin(2) * t940;
t891 = (-pkin(7) * qJDD(1) + t950 * t980 - t920) * t938 + t969;
t907 = -g(3) * t938 + (t920 + 0.2e1 * t972) * t940;
t970 = qJDD(1) * t940;
t935 = t940 ^ 2;
t977 = t935 * t950;
t892 = -pkin(2) * t977 + pkin(7) * t970 + t907;
t943 = sin(qJ(3));
t947 = cos(qJ(3));
t866 = t943 * t891 + t947 * t892;
t974 = qJD(1) * t940;
t975 = qJD(1) * t938;
t918 = -t943 * t975 + t947 * t974;
t958 = t938 * t947 + t940 * t943;
t919 = t958 * qJD(1);
t902 = -pkin(3) * t918 - pkin(8) * t919;
t949 = qJD(3) ^ 2;
t849 = -pkin(3) * t949 + qJDD(3) * pkin(8) + t902 * t918 + t866;
t934 = t938 ^ 2;
t924 = t944 * g(1) - t948 * g(2);
t963 = qJDD(2) - t924;
t903 = (-pkin(1) - t980) * qJDD(1) + (-qJ(2) + (-t934 - t935) * pkin(7)) * t950 + t963;
t915 = t919 * qJD(3);
t971 = qJDD(1) * t938;
t904 = -t943 * t971 + t947 * t970 - t915;
t973 = t918 * qJD(3);
t905 = t958 * qJDD(1) + t973;
t857 = (-t905 - t973) * pkin(8) + (-t904 + t915) * pkin(3) + t903;
t942 = sin(qJ(4));
t946 = cos(qJ(4));
t839 = -t942 * t849 + t946 * t857;
t909 = qJD(3) * t946 - t919 * t942;
t877 = qJD(4) * t909 + qJDD(3) * t942 + t905 * t946;
t901 = qJDD(4) - t904;
t910 = qJD(3) * t942 + t919 * t946;
t916 = qJD(4) - t918;
t828 = (t909 * t916 - t877) * qJ(5) + (t909 * t910 + t901) * pkin(4) + t839;
t840 = t946 * t849 + t942 * t857;
t876 = -qJD(4) * t910 + qJDD(3) * t946 - t905 * t942;
t887 = pkin(4) * t916 - qJ(5) * t910;
t908 = t909 ^ 2;
t830 = -pkin(4) * t908 + qJ(5) * t876 - t887 * t916 + t840;
t937 = sin(pkin(11));
t939 = cos(pkin(11));
t882 = t909 * t937 + t910 * t939;
t822 = -0.2e1 * qJD(5) * t882 + t939 * t828 - t937 * t830;
t852 = t876 * t937 + t877 * t939;
t881 = t909 * t939 - t910 * t937;
t820 = (t881 * t916 - t852) * pkin(9) + (t881 * t882 + t901) * pkin(5) + t822;
t823 = 0.2e1 * qJD(5) * t881 + t937 * t828 + t939 * t830;
t851 = t876 * t939 - t877 * t937;
t869 = pkin(5) * t916 - pkin(9) * t882;
t880 = t881 ^ 2;
t821 = -pkin(5) * t880 + pkin(9) * t851 - t869 * t916 + t823;
t941 = sin(qJ(6));
t945 = cos(qJ(6));
t818 = t820 * t945 - t821 * t941;
t862 = t881 * t945 - t882 * t941;
t836 = qJD(6) * t862 + t851 * t941 + t852 * t945;
t863 = t881 * t941 + t882 * t945;
t846 = -mrSges(7,1) * t862 + mrSges(7,2) * t863;
t914 = qJD(6) + t916;
t855 = -mrSges(7,2) * t914 + mrSges(7,3) * t862;
t898 = qJDD(6) + t901;
t811 = m(7) * t818 + mrSges(7,1) * t898 - mrSges(7,3) * t836 - t846 * t863 + t855 * t914;
t819 = t820 * t941 + t821 * t945;
t835 = -qJD(6) * t863 + t851 * t945 - t852 * t941;
t856 = mrSges(7,1) * t914 - mrSges(7,3) * t863;
t812 = m(7) * t819 - mrSges(7,2) * t898 + mrSges(7,3) * t835 + t846 * t862 - t856 * t914;
t805 = t945 * t811 + t941 * t812;
t864 = -mrSges(6,1) * t881 + mrSges(6,2) * t882;
t867 = -mrSges(6,2) * t916 + mrSges(6,3) * t881;
t803 = m(6) * t822 + mrSges(6,1) * t901 - mrSges(6,3) * t852 - t864 * t882 + t867 * t916 + t805;
t868 = mrSges(6,1) * t916 - mrSges(6,3) * t882;
t964 = -t811 * t941 + t945 * t812;
t804 = m(6) * t823 - mrSges(6,2) * t901 + mrSges(6,3) * t851 + t864 * t881 - t868 * t916 + t964;
t799 = t939 * t803 + t937 * t804;
t859 = Ifges(6,4) * t882 + Ifges(6,2) * t881 + Ifges(6,6) * t916;
t860 = Ifges(6,1) * t882 + Ifges(6,4) * t881 + Ifges(6,5) * t916;
t871 = Ifges(5,4) * t910 + Ifges(5,2) * t909 + Ifges(5,6) * t916;
t872 = Ifges(5,1) * t910 + Ifges(5,4) * t909 + Ifges(5,5) * t916;
t842 = Ifges(7,4) * t863 + Ifges(7,2) * t862 + Ifges(7,6) * t914;
t843 = Ifges(7,1) * t863 + Ifges(7,4) * t862 + Ifges(7,5) * t914;
t955 = -mrSges(7,1) * t818 + mrSges(7,2) * t819 - Ifges(7,5) * t836 - Ifges(7,6) * t835 - Ifges(7,3) * t898 - t863 * t842 + t862 * t843;
t981 = mrSges(5,1) * t839 + mrSges(6,1) * t822 - mrSges(5,2) * t840 - mrSges(6,2) * t823 + Ifges(5,5) * t877 + Ifges(6,5) * t852 + Ifges(5,6) * t876 + Ifges(6,6) * t851 + pkin(4) * t799 + pkin(5) * t805 + t882 * t859 - t881 * t860 + t910 * t871 - t909 * t872 + (Ifges(5,3) + Ifges(6,3)) * t901 - t955;
t978 = mrSges(3,2) * t938;
t883 = -mrSges(5,1) * t909 + mrSges(5,2) * t910;
t886 = -mrSges(5,2) * t916 + mrSges(5,3) * t909;
t797 = m(5) * t839 + mrSges(5,1) * t901 - mrSges(5,3) * t877 - t883 * t910 + t886 * t916 + t799;
t888 = mrSges(5,1) * t916 - mrSges(5,3) * t910;
t965 = -t803 * t937 + t939 * t804;
t798 = m(5) * t840 - mrSges(5,2) * t901 + mrSges(5,3) * t876 + t883 * t909 - t888 * t916 + t965;
t791 = -t797 * t942 + t946 * t798;
t899 = -mrSges(4,1) * t918 + mrSges(4,2) * t919;
t912 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t919;
t789 = m(4) * t866 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t904 - qJD(3) * t912 + t899 * t918 + t791;
t865 = t891 * t947 - t943 * t892;
t848 = -qJDD(3) * pkin(3) - pkin(8) * t949 + t919 * t902 - t865;
t838 = -pkin(4) * t876 - qJ(5) * t908 + t910 * t887 + qJDD(5) + t848;
t825 = -pkin(5) * t851 - pkin(9) * t880 + t869 * t882 + t838;
t959 = m(7) * t825 - t835 * mrSges(7,1) + t836 * mrSges(7,2) - t862 * t855 + t863 * t856;
t816 = m(6) * t838 - t851 * mrSges(6,1) + mrSges(6,2) * t852 - t881 * t867 + t868 * t882 + t959;
t815 = -m(5) * t848 + t876 * mrSges(5,1) - mrSges(5,2) * t877 + t909 * t886 - t888 * t910 - t816;
t911 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t918;
t814 = m(4) * t865 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t905 + qJD(3) * t911 - t899 * t919 + t815;
t782 = t943 * t789 + t947 * t814;
t906 = -t938 * t920 + t969;
t957 = mrSges(3,3) * qJDD(1) + t950 * (-mrSges(3,1) * t940 + t978);
t780 = m(3) * t906 - t957 * t938 + t782;
t966 = t947 * t789 - t943 * t814;
t781 = m(3) * t907 + t957 * t940 + t966;
t967 = -t780 * t938 + t940 * t781;
t771 = m(2) * t925 - mrSges(2,1) * t950 - qJDD(1) * mrSges(2,2) + t967;
t917 = -qJDD(1) * pkin(1) - t950 * qJ(2) + t963;
t790 = t946 * t797 + t942 * t798;
t954 = m(4) * t903 - t904 * mrSges(4,1) + t905 * mrSges(4,2) - t918 * t911 + t919 * t912 + t790;
t953 = -m(3) * t917 + mrSges(3,1) * t970 - t954 + (t934 * t950 + t977) * mrSges(3,3);
t784 = t953 + (mrSges(2,1) - t978) * qJDD(1) + m(2) * t924 - t950 * mrSges(2,2);
t976 = t944 * t771 + t948 * t784;
t773 = t940 * t780 + t938 * t781;
t968 = t948 * t771 - t784 * t944;
t962 = Ifges(3,1) * t938 + Ifges(3,4) * t940;
t961 = Ifges(3,4) * t938 + Ifges(3,2) * t940;
t960 = Ifges(3,5) * t938 + Ifges(3,6) * t940;
t841 = Ifges(7,5) * t863 + Ifges(7,6) * t862 + Ifges(7,3) * t914;
t806 = -mrSges(7,1) * t825 + mrSges(7,3) * t819 + Ifges(7,4) * t836 + Ifges(7,2) * t835 + Ifges(7,6) * t898 - t841 * t863 + t843 * t914;
t807 = mrSges(7,2) * t825 - mrSges(7,3) * t818 + Ifges(7,1) * t836 + Ifges(7,4) * t835 + Ifges(7,5) * t898 + t841 * t862 - t842 * t914;
t858 = Ifges(6,5) * t882 + Ifges(6,6) * t881 + Ifges(6,3) * t916;
t792 = -mrSges(6,1) * t838 + mrSges(6,3) * t823 + Ifges(6,4) * t852 + Ifges(6,2) * t851 + Ifges(6,6) * t901 - pkin(5) * t959 + pkin(9) * t964 + t945 * t806 + t941 * t807 - t882 * t858 + t916 * t860;
t793 = mrSges(6,2) * t838 - mrSges(6,3) * t822 + Ifges(6,1) * t852 + Ifges(6,4) * t851 + Ifges(6,5) * t901 - pkin(9) * t805 - t806 * t941 + t807 * t945 + t858 * t881 - t859 * t916;
t870 = Ifges(5,5) * t910 + Ifges(5,6) * t909 + Ifges(5,3) * t916;
t775 = -mrSges(5,1) * t848 + mrSges(5,3) * t840 + Ifges(5,4) * t877 + Ifges(5,2) * t876 + Ifges(5,6) * t901 - pkin(4) * t816 + qJ(5) * t965 + t939 * t792 + t937 * t793 - t910 * t870 + t916 * t872;
t776 = mrSges(5,2) * t848 - mrSges(5,3) * t839 + Ifges(5,1) * t877 + Ifges(5,4) * t876 + Ifges(5,5) * t901 - qJ(5) * t799 - t792 * t937 + t793 * t939 + t870 * t909 - t871 * t916;
t893 = Ifges(4,5) * t919 + Ifges(4,6) * t918 + Ifges(4,3) * qJD(3);
t894 = Ifges(4,4) * t919 + Ifges(4,2) * t918 + Ifges(4,6) * qJD(3);
t768 = mrSges(4,2) * t903 - mrSges(4,3) * t865 + Ifges(4,1) * t905 + Ifges(4,4) * t904 + Ifges(4,5) * qJDD(3) - pkin(8) * t790 - qJD(3) * t894 - t775 * t942 + t776 * t946 + t893 * t918;
t895 = Ifges(4,1) * t919 + Ifges(4,4) * t918 + Ifges(4,5) * qJD(3);
t774 = -mrSges(4,1) * t903 + mrSges(4,3) * t866 + Ifges(4,4) * t905 + Ifges(4,2) * t904 + Ifges(4,6) * qJDD(3) - pkin(3) * t790 + qJD(3) * t895 - t919 * t893 - t981;
t922 = t960 * qJD(1);
t764 = -mrSges(3,1) * t917 + mrSges(3,3) * t907 - pkin(2) * t954 + pkin(7) * t966 + t961 * qJDD(1) + t943 * t768 + t947 * t774 - t922 * t975;
t767 = mrSges(3,2) * t917 - mrSges(3,3) * t906 - pkin(7) * t782 + t962 * qJDD(1) + t947 * t768 - t943 * t774 + t922 * t974;
t786 = mrSges(3,2) * t971 - t953;
t956 = mrSges(2,1) * t924 - mrSges(2,2) * t925 + Ifges(2,3) * qJDD(1) - pkin(1) * t786 + qJ(2) * t967 + t940 * t764 + t938 * t767;
t952 = mrSges(4,1) * t865 - mrSges(4,2) * t866 + Ifges(4,5) * t905 + Ifges(4,6) * t904 + Ifges(4,3) * qJDD(3) + pkin(3) * t815 + pkin(8) * t791 + t946 * t775 + t942 * t776 + t919 * t894 - t918 * t895;
t765 = mrSges(2,1) * g(3) - pkin(2) * t782 + (Ifges(2,6) - t960) * qJDD(1) - t952 - pkin(1) * t773 - mrSges(3,1) * t906 + mrSges(3,2) * t907 + mrSges(2,3) * t925 + (-t938 * t961 + t940 * t962 + Ifges(2,5)) * t950;
t762 = -mrSges(2,2) * g(3) - mrSges(2,3) * t924 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t950 - qJ(2) * t773 - t764 * t938 + t767 * t940;
t1 = [-m(1) * g(1) + t968; -m(1) * g(2) + t976; (-m(1) - m(2)) * g(3) + t773; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t976 + t948 * t762 - t944 * t765; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t968 + t944 * t762 + t948 * t765; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t956; t956; t786; t952; t981; t816; -t955;];
tauJB  = t1;
