% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RRRPRR6
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
% Datum: 2019-05-07 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RRRPRR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:01:57
% EndTime: 2019-05-07 11:03:41
% DurationCPUTime: 48.28s
% Computational Cost: add. (790132->387), mult. (1676050->484), div. (0->0), fcn. (1239434->12), ass. (0->153)
t944 = sin(qJ(1));
t949 = cos(qJ(1));
t930 = t944 * g(1) - t949 * g(2);
t951 = qJD(1) ^ 2;
t913 = -qJDD(1) * pkin(1) - t951 * pkin(7) - t930;
t943 = sin(qJ(2));
t948 = cos(qJ(2));
t967 = qJD(1) * qJD(2);
t966 = t948 * t967;
t924 = qJDD(1) * t943 + t966;
t934 = t943 * t967;
t925 = qJDD(1) * t948 - t934;
t878 = (-t924 - t966) * pkin(8) + (-t925 + t934) * pkin(2) + t913;
t931 = -g(1) * t949 - g(2) * t944;
t914 = -pkin(1) * t951 + qJDD(1) * pkin(7) + t931;
t903 = -g(3) * t943 + t948 * t914;
t923 = (-pkin(2) * t948 - pkin(8) * t943) * qJD(1);
t950 = qJD(2) ^ 2;
t968 = qJD(1) * t948;
t881 = -pkin(2) * t950 + qJDD(2) * pkin(8) + t923 * t968 + t903;
t942 = sin(qJ(3));
t947 = cos(qJ(3));
t859 = t947 * t878 - t942 * t881;
t969 = qJD(1) * t943;
t920 = qJD(2) * t947 - t942 * t969;
t894 = qJD(3) * t920 + qJDD(2) * t942 + t924 * t947;
t919 = qJDD(3) - t925;
t921 = qJD(2) * t942 + t947 * t969;
t933 = qJD(3) - t968;
t845 = (t920 * t933 - t894) * qJ(4) + (t920 * t921 + t919) * pkin(3) + t859;
t860 = t942 * t878 + t947 * t881;
t893 = -qJD(3) * t921 + qJDD(2) * t947 - t924 * t942;
t900 = pkin(3) * t933 - qJ(4) * t921;
t918 = t920 ^ 2;
t847 = -pkin(3) * t918 + qJ(4) * t893 - t900 * t933 + t860;
t938 = sin(pkin(11));
t939 = cos(pkin(11));
t897 = t920 * t938 + t921 * t939;
t826 = -0.2e1 * qJD(4) * t897 + t939 * t845 - t938 * t847;
t869 = t893 * t938 + t894 * t939;
t896 = t920 * t939 - t921 * t938;
t816 = (t896 * t933 - t869) * pkin(9) + (t896 * t897 + t919) * pkin(4) + t826;
t827 = 0.2e1 * qJD(4) * t896 + t938 * t845 + t939 * t847;
t868 = t893 * t939 - t894 * t938;
t884 = pkin(4) * t933 - pkin(9) * t897;
t895 = t896 ^ 2;
t824 = -pkin(4) * t895 + pkin(9) * t868 - t884 * t933 + t827;
t941 = sin(qJ(5));
t946 = cos(qJ(5));
t810 = t946 * t816 - t941 * t824;
t873 = t896 * t946 - t897 * t941;
t841 = qJD(5) * t873 + t868 * t941 + t869 * t946;
t874 = t896 * t941 + t897 * t946;
t915 = qJDD(5) + t919;
t932 = qJD(5) + t933;
t807 = (t873 * t932 - t841) * pkin(10) + (t873 * t874 + t915) * pkin(5) + t810;
t811 = t941 * t816 + t946 * t824;
t840 = -qJD(5) * t874 + t868 * t946 - t869 * t941;
t863 = pkin(5) * t932 - pkin(10) * t874;
t872 = t873 ^ 2;
t808 = -pkin(5) * t872 + pkin(10) * t840 - t863 * t932 + t811;
t940 = sin(qJ(6));
t945 = cos(qJ(6));
t805 = t807 * t945 - t808 * t940;
t855 = t873 * t945 - t874 * t940;
t822 = qJD(6) * t855 + t840 * t940 + t841 * t945;
t856 = t873 * t940 + t874 * t945;
t835 = -mrSges(7,1) * t855 + mrSges(7,2) * t856;
t927 = qJD(6) + t932;
t848 = -mrSges(7,2) * t927 + mrSges(7,3) * t855;
t909 = qJDD(6) + t915;
t798 = m(7) * t805 + mrSges(7,1) * t909 - mrSges(7,3) * t822 - t835 * t856 + t848 * t927;
t806 = t807 * t940 + t808 * t945;
t821 = -qJD(6) * t856 + t840 * t945 - t841 * t940;
t849 = mrSges(7,1) * t927 - mrSges(7,3) * t856;
t799 = m(7) * t806 - mrSges(7,2) * t909 + mrSges(7,3) * t821 + t835 * t855 - t849 * t927;
t792 = t945 * t798 + t940 * t799;
t857 = -mrSges(6,1) * t873 + mrSges(6,2) * t874;
t861 = -mrSges(6,2) * t932 + mrSges(6,3) * t873;
t789 = m(6) * t810 + mrSges(6,1) * t915 - mrSges(6,3) * t841 - t857 * t874 + t861 * t932 + t792;
t862 = mrSges(6,1) * t932 - mrSges(6,3) * t874;
t961 = -t798 * t940 + t945 * t799;
t790 = m(6) * t811 - mrSges(6,2) * t915 + mrSges(6,3) * t840 + t857 * t873 - t862 * t932 + t961;
t785 = t946 * t789 + t941 * t790;
t875 = -mrSges(5,1) * t896 + mrSges(5,2) * t897;
t882 = -mrSges(5,2) * t933 + mrSges(5,3) * t896;
t783 = m(5) * t826 + mrSges(5,1) * t919 - mrSges(5,3) * t869 - t875 * t897 + t882 * t933 + t785;
t883 = mrSges(5,1) * t933 - mrSges(5,3) * t897;
t962 = -t789 * t941 + t946 * t790;
t784 = m(5) * t827 - mrSges(5,2) * t919 + mrSges(5,3) * t868 + t875 * t896 - t883 * t933 + t962;
t777 = t939 * t783 + t938 * t784;
t866 = Ifges(5,4) * t897 + Ifges(5,2) * t896 + Ifges(5,6) * t933;
t867 = Ifges(5,1) * t897 + Ifges(5,4) * t896 + Ifges(5,5) * t933;
t886 = Ifges(4,4) * t921 + Ifges(4,2) * t920 + Ifges(4,6) * t933;
t887 = Ifges(4,1) * t921 + Ifges(4,4) * t920 + Ifges(4,5) * t933;
t851 = Ifges(6,4) * t874 + Ifges(6,2) * t873 + Ifges(6,6) * t932;
t852 = Ifges(6,1) * t874 + Ifges(6,4) * t873 + Ifges(6,5) * t932;
t831 = Ifges(7,4) * t856 + Ifges(7,2) * t855 + Ifges(7,6) * t927;
t832 = Ifges(7,1) * t856 + Ifges(7,4) * t855 + Ifges(7,5) * t927;
t956 = -mrSges(7,1) * t805 + mrSges(7,2) * t806 - Ifges(7,5) * t822 - Ifges(7,6) * t821 - Ifges(7,3) * t909 - t856 * t831 + t855 * t832;
t954 = -mrSges(6,1) * t810 + mrSges(6,2) * t811 - Ifges(6,5) * t841 - Ifges(6,6) * t840 - Ifges(6,3) * t915 - pkin(5) * t792 - t874 * t851 + t873 * t852 + t956;
t973 = mrSges(4,1) * t859 + mrSges(5,1) * t826 - mrSges(4,2) * t860 - mrSges(5,2) * t827 + Ifges(4,5) * t894 + Ifges(5,5) * t869 + Ifges(4,6) * t893 + Ifges(5,6) * t868 + pkin(3) * t777 + pkin(4) * t785 + t897 * t866 - t896 * t867 + t921 * t886 - t920 * t887 + (Ifges(4,3) + Ifges(5,3)) * t919 - t954;
t902 = -t948 * g(3) - t943 * t914;
t880 = -qJDD(2) * pkin(2) - pkin(8) * t950 + t923 * t969 - t902;
t858 = -pkin(3) * t893 - qJ(4) * t918 + t921 * t900 + qJDD(4) + t880;
t829 = -pkin(4) * t868 - pkin(9) * t895 + t897 * t884 + t858;
t813 = -pkin(5) * t840 - pkin(10) * t872 + t863 * t874 + t829;
t830 = Ifges(7,5) * t856 + Ifges(7,6) * t855 + Ifges(7,3) * t927;
t793 = -mrSges(7,1) * t813 + mrSges(7,3) * t806 + Ifges(7,4) * t822 + Ifges(7,2) * t821 + Ifges(7,6) * t909 - t830 * t856 + t832 * t927;
t794 = mrSges(7,2) * t813 - mrSges(7,3) * t805 + Ifges(7,1) * t822 + Ifges(7,4) * t821 + Ifges(7,5) * t909 + t830 * t855 - t831 * t927;
t850 = Ifges(6,5) * t874 + Ifges(6,6) * t873 + Ifges(6,3) * t932;
t960 = m(7) * t813 - t821 * mrSges(7,1) + t822 * mrSges(7,2) - t855 * t848 + t856 * t849;
t778 = -mrSges(6,1) * t829 + mrSges(6,3) * t811 + Ifges(6,4) * t841 + Ifges(6,2) * t840 + Ifges(6,6) * t915 - pkin(5) * t960 + pkin(10) * t961 + t945 * t793 + t940 * t794 - t874 * t850 + t932 * t852;
t779 = mrSges(6,2) * t829 - mrSges(6,3) * t810 + Ifges(6,1) * t841 + Ifges(6,4) * t840 + Ifges(6,5) * t915 - pkin(10) * t792 - t793 * t940 + t794 * t945 + t850 * t873 - t851 * t932;
t865 = Ifges(5,5) * t897 + Ifges(5,6) * t896 + Ifges(5,3) * t933;
t957 = m(6) * t829 - t840 * mrSges(6,1) + t841 * mrSges(6,2) - t873 * t861 + t874 * t862 + t960;
t772 = -mrSges(5,1) * t858 + mrSges(5,3) * t827 + Ifges(5,4) * t869 + Ifges(5,2) * t868 + Ifges(5,6) * t919 - pkin(4) * t957 + pkin(9) * t962 + t946 * t778 + t941 * t779 - t897 * t865 + t933 * t867;
t773 = mrSges(5,2) * t858 - mrSges(5,3) * t826 + Ifges(5,1) * t869 + Ifges(5,4) * t868 + Ifges(5,5) * t919 - pkin(9) * t785 - t778 * t941 + t779 * t946 + t865 * t896 - t866 * t933;
t803 = m(5) * t858 - t868 * mrSges(5,1) + t869 * mrSges(5,2) - t896 * t882 + t897 * t883 + t957;
t885 = Ifges(4,5) * t921 + Ifges(4,6) * t920 + Ifges(4,3) * t933;
t963 = -t783 * t938 + t939 * t784;
t755 = -mrSges(4,1) * t880 + mrSges(4,3) * t860 + Ifges(4,4) * t894 + Ifges(4,2) * t893 + Ifges(4,6) * t919 - pkin(3) * t803 + qJ(4) * t963 + t939 * t772 + t938 * t773 - t921 * t885 + t933 * t887;
t756 = mrSges(4,2) * t880 - mrSges(4,3) * t859 + Ifges(4,1) * t894 + Ifges(4,4) * t893 + Ifges(4,5) * t919 - qJ(4) * t777 - t772 * t938 + t773 * t939 + t885 * t920 - t886 * t933;
t898 = -mrSges(4,1) * t920 + mrSges(4,2) * t921;
t899 = -mrSges(4,2) * t933 + mrSges(4,3) * t920;
t775 = m(4) * t859 + mrSges(4,1) * t919 - mrSges(4,3) * t894 - t898 * t921 + t899 * t933 + t777;
t901 = mrSges(4,1) * t933 - mrSges(4,3) * t921;
t776 = m(4) * t860 - mrSges(4,2) * t919 + mrSges(4,3) * t893 + t898 * t920 - t901 * t933 + t963;
t771 = -t775 * t942 + t947 * t776;
t802 = -m(4) * t880 + t893 * mrSges(4,1) - t894 * mrSges(4,2) + t920 * t899 - t921 * t901 - t803;
t911 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t943 + Ifges(3,2) * t948) * qJD(1);
t912 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t943 + Ifges(3,4) * t948) * qJD(1);
t972 = mrSges(3,1) * t902 - mrSges(3,2) * t903 + Ifges(3,5) * t924 + Ifges(3,6) * t925 + Ifges(3,3) * qJDD(2) + pkin(2) * t802 + pkin(8) * t771 + t947 * t755 + t942 * t756 + (t911 * t943 - t912 * t948) * qJD(1);
t922 = (-mrSges(3,1) * t948 + mrSges(3,2) * t943) * qJD(1);
t928 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t969;
t769 = m(3) * t903 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t925 - qJD(2) * t928 + t922 * t968 + t771;
t929 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t968;
t801 = m(3) * t902 + qJDD(2) * mrSges(3,1) - t924 * mrSges(3,3) + qJD(2) * t929 - t922 * t969 + t802;
t964 = t948 * t769 - t801 * t943;
t761 = m(2) * t931 - mrSges(2,1) * t951 - qJDD(1) * mrSges(2,2) + t964;
t770 = t775 * t947 + t776 * t942;
t955 = -m(3) * t913 + t925 * mrSges(3,1) - mrSges(3,2) * t924 - t928 * t969 + t929 * t968 - t770;
t765 = m(2) * t930 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t951 + t955;
t970 = t944 * t761 + t949 * t765;
t763 = t943 * t769 + t948 * t801;
t965 = t949 * t761 - t765 * t944;
t910 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t943 + Ifges(3,6) * t948) * qJD(1);
t754 = mrSges(3,2) * t913 - mrSges(3,3) * t902 + Ifges(3,1) * t924 + Ifges(3,4) * t925 + Ifges(3,5) * qJDD(2) - pkin(8) * t770 - qJD(2) * t911 - t755 * t942 + t756 * t947 + t910 * t968;
t758 = -mrSges(3,1) * t913 + mrSges(3,3) * t903 + Ifges(3,4) * t924 + Ifges(3,2) * t925 + Ifges(3,6) * qJDD(2) - pkin(2) * t770 + qJD(2) * t912 - t910 * t969 - t973;
t958 = mrSges(2,1) * t930 - mrSges(2,2) * t931 + Ifges(2,3) * qJDD(1) + pkin(1) * t955 + pkin(7) * t964 + t943 * t754 + t948 * t758;
t752 = mrSges(2,1) * g(3) + mrSges(2,3) * t931 + t951 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t763 - t972;
t751 = -mrSges(2,2) * g(3) - mrSges(2,3) * t930 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t951 - pkin(7) * t763 + t754 * t948 - t758 * t943;
t1 = [-m(1) * g(1) + t965; -m(1) * g(2) + t970; (-m(1) - m(2)) * g(3) + t763; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t970 + t949 * t751 - t944 * t752; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t965 + t944 * t751 + t949 * t752; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t958; t958; t972; t973; t803; -t954; -t956;];
tauJB  = t1;
