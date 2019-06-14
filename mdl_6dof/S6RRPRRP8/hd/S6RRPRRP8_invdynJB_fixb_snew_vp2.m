% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 18:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRPRRP8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:20:31
% EndTime: 2019-05-06 18:20:50
% DurationCPUTime: 18.09s
% Computational Cost: add. (279751->362), mult. (608478->442), div. (0->0), fcn. (438566->10), ass. (0->143)
t969 = Ifges(6,1) + Ifges(7,1);
t959 = Ifges(6,4) - Ifges(7,5);
t967 = Ifges(7,4) + Ifges(6,5);
t968 = Ifges(6,2) + Ifges(7,3);
t966 = Ifges(6,6) - Ifges(7,6);
t965 = -Ifges(6,3) - Ifges(7,2);
t926 = sin(pkin(10));
t927 = cos(pkin(10));
t930 = sin(qJ(2));
t955 = qJD(1) * t930;
t906 = t927 * qJD(2) - t926 * t955;
t907 = t926 * qJD(2) + t927 * t955;
t929 = sin(qJ(4));
t932 = cos(qJ(4));
t879 = t932 * t906 - t929 * t907;
t880 = t929 * t906 + t932 * t907;
t928 = sin(qJ(5));
t961 = cos(qJ(5));
t860 = -t879 * t961 + t928 * t880;
t861 = t928 * t879 + t880 * t961;
t933 = cos(qJ(2));
t954 = t933 * qJD(1);
t921 = qJD(4) - t954;
t920 = qJD(5) + t921;
t964 = t860 * t968 - t861 * t959 - t920 * t966;
t963 = -t959 * t860 + t861 * t969 + t967 * t920;
t931 = sin(qJ(1));
t934 = cos(qJ(1));
t918 = -t934 * g(1) - t931 * g(2);
t936 = qJD(1) ^ 2;
t902 = -t936 * pkin(1) + qJDD(1) * pkin(7) + t918;
t882 = -t933 * g(3) - t930 * t902;
t910 = (-pkin(2) * t933 - qJ(3) * t930) * qJD(1);
t935 = qJD(2) ^ 2;
t867 = -qJDD(2) * pkin(2) - t935 * qJ(3) + t910 * t955 + qJDD(3) - t882;
t953 = qJD(1) * qJD(2);
t951 = t933 * t953;
t912 = t930 * qJDD(1) + t951;
t886 = t927 * qJDD(2) - t926 * t912;
t888 = -pkin(3) * t954 - t907 * pkin(8);
t905 = t906 ^ 2;
t844 = -t886 * pkin(3) - t905 * pkin(8) + t907 * t888 + t867;
t887 = t926 * qJDD(2) + t927 * t912;
t853 = -t880 * qJD(4) + t932 * t886 - t929 * t887;
t871 = t921 * pkin(4) - t880 * pkin(9);
t878 = t879 ^ 2;
t813 = -t853 * pkin(4) - t878 * pkin(9) + t880 * t871 + t844;
t854 = t879 * qJD(4) + t929 * t886 + t932 * t887;
t821 = t861 * qJD(5) - t853 * t961 + t928 * t854;
t822 = -t860 * qJD(5) + t928 * t853 + t854 * t961;
t806 = t813 + (t860 * t920 - t822) * qJ(6) + (t861 * t920 + t821) * pkin(5) - 0.2e1 * qJD(6) * t861;
t848 = -t860 * mrSges(7,2) + t920 * mrSges(7,3);
t851 = -t920 * mrSges(7,1) + t861 * mrSges(7,2);
t797 = m(7) * t806 + t821 * mrSges(7,1) - t822 * mrSges(7,3) + t860 * t848 - t861 * t851;
t917 = t931 * g(1) - t934 * g(2);
t901 = -qJDD(1) * pkin(1) - t936 * pkin(7) - t917;
t922 = t930 * t953;
t913 = t933 * qJDD(1) - t922;
t865 = (-t912 - t951) * qJ(3) + (-t913 + t922) * pkin(2) + t901;
t883 = -t930 * g(3) + t933 * t902;
t868 = -t935 * pkin(2) + qJDD(2) * qJ(3) + t910 * t954 + t883;
t842 = -0.2e1 * qJD(3) * t907 + t927 * t865 - t926 * t868;
t827 = (-t906 * t954 - t887) * pkin(8) + (t906 * t907 - t913) * pkin(3) + t842;
t843 = 0.2e1 * qJD(3) * t906 + t926 * t865 + t927 * t868;
t829 = -t905 * pkin(3) + t886 * pkin(8) + t888 * t954 + t843;
t811 = t932 * t827 - t929 * t829;
t909 = qJDD(4) - t913;
t808 = (t879 * t921 - t854) * pkin(9) + (t879 * t880 + t909) * pkin(4) + t811;
t812 = t929 * t827 + t932 * t829;
t810 = -t878 * pkin(4) + t853 * pkin(9) - t921 * t871 + t812;
t804 = t928 * t808 + t810 * t961;
t839 = t860 * pkin(5) - t861 * qJ(6);
t903 = qJDD(5) + t909;
t919 = t920 ^ 2;
t800 = -pkin(5) * t919 + qJ(6) * t903 + 0.2e1 * qJD(6) * t920 - t839 * t860 + t804;
t957 = t966 * t860 - t967 * t861 + t965 * t920;
t782 = -mrSges(6,1) * t813 - mrSges(7,1) * t806 + mrSges(7,2) * t800 + mrSges(6,3) * t804 - pkin(5) * t797 - t821 * t968 + t959 * t822 + t957 * t861 + t966 * t903 + t963 * t920;
t803 = t808 * t961 - t928 * t810;
t801 = -t903 * pkin(5) - t919 * qJ(6) + t861 * t839 + qJDD(6) - t803;
t783 = mrSges(6,2) * t813 + mrSges(7,2) * t801 - mrSges(6,3) * t803 - mrSges(7,3) * t806 - qJ(6) * t797 - t959 * t821 + t822 * t969 + t957 * t860 + t967 * t903 + t964 * t920;
t855 = Ifges(5,5) * t880 + Ifges(5,6) * t879 + Ifges(5,3) * t921;
t857 = Ifges(5,1) * t880 + Ifges(5,4) * t879 + Ifges(5,5) * t921;
t849 = -t920 * mrSges(6,2) - t860 * mrSges(6,3);
t850 = t920 * mrSges(6,1) - t861 * mrSges(6,3);
t942 = m(6) * t813 + t821 * mrSges(6,1) + t822 * mrSges(6,2) + t860 * t849 + t861 * t850 + t797;
t952 = m(7) * t800 + t903 * mrSges(7,3) + t920 * t851;
t840 = t860 * mrSges(7,1) - t861 * mrSges(7,3);
t956 = -t860 * mrSges(6,1) - t861 * mrSges(6,2) - t840;
t960 = -mrSges(6,3) - mrSges(7,2);
t787 = m(6) * t804 - t903 * mrSges(6,2) + t821 * t960 - t920 * t850 + t860 * t956 + t952;
t946 = -m(7) * t801 + t903 * mrSges(7,1) + t920 * t848;
t789 = m(6) * t803 + t903 * mrSges(6,1) + t822 * t960 + t920 * t849 + t861 * t956 + t946;
t947 = t787 * t961 - t789 * t928;
t771 = -mrSges(5,1) * t844 + mrSges(5,3) * t812 + Ifges(5,4) * t854 + Ifges(5,2) * t853 + Ifges(5,6) * t909 - pkin(4) * t942 + pkin(9) * t947 + t782 * t961 + t928 * t783 - t880 * t855 + t921 * t857;
t784 = t928 * t787 + t789 * t961;
t856 = Ifges(5,4) * t880 + Ifges(5,2) * t879 + Ifges(5,6) * t921;
t772 = mrSges(5,2) * t844 - mrSges(5,3) * t811 + Ifges(5,1) * t854 + Ifges(5,4) * t853 + Ifges(5,5) * t909 - pkin(9) * t784 - t928 * t782 + t783 * t961 + t879 * t855 - t921 * t856;
t873 = Ifges(4,5) * t907 + Ifges(4,6) * t906 - Ifges(4,3) * t954;
t875 = Ifges(4,1) * t907 + Ifges(4,4) * t906 - Ifges(4,5) * t954;
t869 = -t921 * mrSges(5,2) + t879 * mrSges(5,3);
t870 = t921 * mrSges(5,1) - t880 * mrSges(5,3);
t938 = m(5) * t844 - t853 * mrSges(5,1) + t854 * mrSges(5,2) - t879 * t869 + t880 * t870 + t942;
t862 = -t879 * mrSges(5,1) + t880 * mrSges(5,2);
t780 = m(5) * t811 + mrSges(5,1) * t909 - mrSges(5,3) * t854 - t862 * t880 + t869 * t921 + t784;
t781 = m(5) * t812 - mrSges(5,2) * t909 + mrSges(5,3) * t853 + t862 * t879 - t870 * t921 + t947;
t948 = -t780 * t929 + t932 * t781;
t754 = -mrSges(4,1) * t867 + mrSges(4,3) * t843 + Ifges(4,4) * t887 + Ifges(4,2) * t886 - Ifges(4,6) * t913 - pkin(3) * t938 + pkin(8) * t948 + t932 * t771 + t929 * t772 - t907 * t873 - t875 * t954;
t776 = t932 * t780 + t929 * t781;
t874 = Ifges(4,4) * t907 + Ifges(4,2) * t906 - Ifges(4,6) * t954;
t757 = mrSges(4,2) * t867 - mrSges(4,3) * t842 + Ifges(4,1) * t887 + Ifges(4,4) * t886 - Ifges(4,5) * t913 - pkin(8) * t776 - t771 * t929 + t772 * t932 + t873 * t906 + t874 * t954;
t881 = -t906 * mrSges(4,1) + t907 * mrSges(4,2);
t944 = mrSges(4,2) * t954 + t906 * mrSges(4,3);
t774 = m(4) * t842 - t913 * mrSges(4,1) - t887 * mrSges(4,3) - t907 * t881 - t944 * t954 + t776;
t885 = -mrSges(4,1) * t954 - t907 * mrSges(4,3);
t775 = m(4) * t843 + mrSges(4,2) * t913 + mrSges(4,3) * t886 + t881 * t906 + t885 * t954 + t948;
t770 = -t774 * t926 + t927 * t775;
t792 = m(4) * t867 - t886 * mrSges(4,1) + t887 * mrSges(4,2) + t907 * t885 - t906 * t944 + t938;
t899 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t930 + Ifges(3,2) * t933) * qJD(1);
t900 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t930 + Ifges(3,4) * t933) * qJD(1);
t962 = mrSges(3,1) * t882 - mrSges(3,2) * t883 + Ifges(3,5) * t912 + Ifges(3,6) * t913 + Ifges(3,3) * qJDD(2) - pkin(2) * t792 + qJ(3) * t770 + t927 * t754 + t926 * t757 + (t899 * t930 - t900 * t933) * qJD(1);
t911 = (-mrSges(3,1) * t933 + mrSges(3,2) * t930) * qJD(1);
t915 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t955;
t768 = m(3) * t883 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t913 - qJD(2) * t915 + t911 * t954 + t770;
t916 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t954;
t791 = m(3) * t882 + qJDD(2) * mrSges(3,1) - t912 * mrSges(3,3) + qJD(2) * t916 - t911 * t955 - t792;
t949 = t933 * t768 - t791 * t930;
t760 = m(2) * t918 - mrSges(2,1) * t936 - qJDD(1) * mrSges(2,2) + t949;
t769 = t774 * t927 + t775 * t926;
t941 = -m(3) * t901 + t913 * mrSges(3,1) - mrSges(3,2) * t912 - t915 * t955 + t916 * t954 - t769;
t764 = m(2) * t917 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t936 + t941;
t958 = t931 * t760 + t934 * t764;
t762 = t930 * t768 + t933 * t791;
t950 = t934 * t760 - t764 * t931;
t898 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t930 + Ifges(3,6) * t933) * qJD(1);
t753 = mrSges(3,2) * t901 - mrSges(3,3) * t882 + Ifges(3,1) * t912 + Ifges(3,4) * t913 + Ifges(3,5) * qJDD(2) - qJ(3) * t769 - qJD(2) * t899 - t754 * t926 + t757 * t927 + t898 * t954;
t796 = t822 * mrSges(7,2) + t861 * t840 - t946;
t940 = -mrSges(6,1) * t803 + mrSges(7,1) * t801 + mrSges(6,2) * t804 - mrSges(7,3) * t800 + pkin(5) * t796 - qJ(6) * t952 + t965 * t903 + t964 * t861 + (qJ(6) * t840 - t963) * t860 - t967 * t822 + (qJ(6) * mrSges(7,2) + t966) * t821;
t937 = mrSges(5,1) * t811 - mrSges(5,2) * t812 + Ifges(5,5) * t854 + Ifges(5,6) * t853 + Ifges(5,3) * t909 + pkin(4) * t784 + t880 * t856 - t879 * t857 - t940;
t756 = -t937 + Ifges(3,6) * qJDD(2) - t898 * t955 + (Ifges(3,2) + Ifges(4,3)) * t913 - pkin(2) * t769 + Ifges(3,4) * t912 + t906 * t875 - t907 * t874 - Ifges(4,6) * t886 - Ifges(4,5) * t887 + qJD(2) * t900 - mrSges(3,1) * t901 + mrSges(3,3) * t883 - mrSges(4,1) * t842 + mrSges(4,2) * t843 - pkin(3) * t776;
t943 = mrSges(2,1) * t917 - mrSges(2,2) * t918 + Ifges(2,3) * qJDD(1) + pkin(1) * t941 + pkin(7) * t949 + t930 * t753 + t933 * t756;
t751 = mrSges(2,1) * g(3) + mrSges(2,3) * t918 + t936 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t762 - t962;
t750 = -mrSges(2,2) * g(3) - mrSges(2,3) * t917 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t936 - pkin(7) * t762 + t753 * t933 - t756 * t930;
t1 = [-m(1) * g(1) + t950; -m(1) * g(2) + t958; (-m(1) - m(2)) * g(3) + t762; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t958 + t934 * t750 - t931 * t751; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t950 + t931 * t750 + t934 * t751; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t943; t943; t962; t792; t937; -t940; t796;];
tauJB  = t1;
