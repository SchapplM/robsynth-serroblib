% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:37:11
% EndTime: 2019-05-05 03:37:23
% DurationCPUTime: 10.07s
% Computational Cost: add. (160990->319), mult. (344971->400), div. (0->0), fcn. (242943->12), ass. (0->137)
t917 = -2 * qJD(4);
t916 = Ifges(6,1) + Ifges(7,1);
t910 = Ifges(6,4) + Ifges(7,4);
t909 = Ifges(6,5) + Ifges(7,5);
t915 = Ifges(6,2) + Ifges(7,2);
t908 = Ifges(6,6) + Ifges(7,6);
t914 = Ifges(6,3) + Ifges(7,3);
t864 = sin(pkin(10));
t867 = cos(pkin(10));
t854 = g(1) * t864 - g(2) * t867;
t855 = -g(1) * t867 - g(2) * t864;
t862 = -g(3) + qJDD(1);
t874 = cos(qJ(2));
t868 = cos(pkin(6));
t871 = sin(qJ(2));
t903 = t868 * t871;
t865 = sin(pkin(6));
t904 = t865 * t871;
t820 = t854 * t903 + t874 * t855 + t862 * t904;
t876 = qJD(2) ^ 2;
t815 = -pkin(2) * t876 + qJDD(2) * pkin(8) + t820;
t836 = -t854 * t865 + t862 * t868;
t870 = sin(qJ(3));
t873 = cos(qJ(3));
t789 = -t870 * t815 + t873 * t836;
t894 = qJD(2) * qJD(3);
t891 = t873 * t894;
t852 = qJDD(2) * t870 + t891;
t785 = (-t852 + t891) * qJ(4) + (t870 * t873 * t876 + qJDD(3)) * pkin(3) + t789;
t790 = t873 * t815 + t870 * t836;
t853 = qJDD(2) * t873 - t870 * t894;
t897 = qJD(2) * t870;
t856 = qJD(3) * pkin(3) - qJ(4) * t897;
t861 = t873 ^ 2;
t786 = -pkin(3) * t861 * t876 + qJ(4) * t853 - qJD(3) * t856 + t790;
t863 = sin(pkin(11));
t866 = cos(pkin(11));
t841 = (t863 * t873 + t866 * t870) * qJD(2);
t777 = t785 * t866 - t863 * t786 + t841 * t917;
t840 = (t863 * t870 - t866 * t873) * qJD(2);
t819 = -t871 * t855 + (t854 * t868 + t862 * t865) * t874;
t778 = t863 * t785 + t866 * t786 + t840 * t917;
t823 = pkin(4) * t840 - pkin(9) * t841;
t875 = qJD(3) ^ 2;
t776 = -pkin(4) * t875 + qJDD(3) * pkin(9) - t823 * t840 + t778;
t880 = -qJDD(2) * pkin(2) - t819;
t787 = -t853 * pkin(3) + qJDD(4) + t856 * t897 + (-qJ(4) * t861 - pkin(8)) * t876 + t880;
t827 = -t852 * t863 + t853 * t866;
t828 = t852 * t866 + t853 * t863;
t781 = (qJD(3) * t840 - t828) * pkin(9) + (qJD(3) * t841 - t827) * pkin(4) + t787;
t869 = sin(qJ(5));
t872 = cos(qJ(5));
t771 = -t869 * t776 + t872 * t781;
t830 = qJD(3) * t872 - t841 * t869;
t803 = qJD(5) * t830 + qJDD(3) * t869 + t828 * t872;
t831 = qJD(3) * t869 + t841 * t872;
t805 = -mrSges(7,1) * t830 + mrSges(7,2) * t831;
t806 = -mrSges(6,1) * t830 + mrSges(6,2) * t831;
t839 = qJD(5) + t840;
t810 = -mrSges(6,2) * t839 + mrSges(6,3) * t830;
t826 = qJDD(5) - t827;
t768 = -0.2e1 * qJD(6) * t831 + (t830 * t839 - t803) * qJ(6) + (t830 * t831 + t826) * pkin(5) + t771;
t809 = -mrSges(7,2) * t839 + mrSges(7,3) * t830;
t893 = m(7) * t768 + t826 * mrSges(7,1) + t839 * t809;
t758 = m(6) * t771 + t826 * mrSges(6,1) + t839 * t810 + (-t805 - t806) * t831 + (-mrSges(6,3) - mrSges(7,3)) * t803 + t893;
t772 = t872 * t776 + t869 * t781;
t802 = -qJD(5) * t831 + qJDD(3) * t872 - t828 * t869;
t811 = pkin(5) * t839 - qJ(6) * t831;
t829 = t830 ^ 2;
t770 = -pkin(5) * t829 + qJ(6) * t802 + 0.2e1 * qJD(6) * t830 - t811 * t839 + t772;
t892 = m(7) * t770 + t802 * mrSges(7,3) + t830 * t805;
t812 = mrSges(7,1) * t839 - mrSges(7,3) * t831;
t898 = -mrSges(6,1) * t839 + mrSges(6,3) * t831 - t812;
t911 = -mrSges(6,2) - mrSges(7,2);
t763 = m(6) * t772 + t802 * mrSges(6,3) + t830 * t806 + t911 * t826 + t898 * t839 + t892;
t756 = -t758 * t869 + t872 * t763;
t822 = mrSges(5,1) * t840 + mrSges(5,2) * t841;
t835 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t841;
t752 = m(5) * t778 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t827 - qJD(3) * t835 - t822 * t840 + t756;
t775 = -qJDD(3) * pkin(4) - pkin(9) * t875 + t841 * t823 - t777;
t773 = -pkin(5) * t802 - qJ(6) * t829 + t811 * t831 + qJDD(6) + t775;
t885 = -m(7) * t773 + t802 * mrSges(7,1) + t830 * t809;
t764 = -m(6) * t775 + t802 * mrSges(6,1) + t911 * t803 + t830 * t810 + t898 * t831 + t885;
t834 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t840;
t760 = m(5) * t777 + qJDD(3) * mrSges(5,1) - t828 * mrSges(5,3) + qJD(3) * t834 - t841 * t822 + t764;
t745 = t863 * t752 + t866 * t760;
t766 = t803 * mrSges(7,2) + t831 * t812 - t885;
t899 = t830 * t910 + t831 * t916 + t839 * t909;
t901 = -t830 * t908 - t831 * t909 - t839 * t914;
t746 = -mrSges(6,1) * t775 + mrSges(6,3) * t772 - mrSges(7,1) * t773 + mrSges(7,3) * t770 - pkin(5) * t766 + qJ(6) * t892 + (-qJ(6) * t812 + t899) * t839 + t901 * t831 + (-mrSges(7,2) * qJ(6) + t908) * t826 + t910 * t803 + t915 * t802;
t765 = -t803 * mrSges(7,3) - t831 * t805 + t893;
t900 = -t830 * t915 - t831 * t910 - t839 * t908;
t753 = mrSges(6,2) * t775 + mrSges(7,2) * t773 - mrSges(6,3) * t771 - mrSges(7,3) * t768 - qJ(6) * t765 + t910 * t802 + t803 * t916 + t909 * t826 - t901 * t830 + t900 * t839;
t817 = Ifges(5,4) * t841 - Ifges(5,2) * t840 + Ifges(5,6) * qJD(3);
t818 = Ifges(5,1) * t841 - Ifges(5,4) * t840 + Ifges(5,5) * qJD(3);
t844 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t870 + Ifges(4,2) * t873) * qJD(2);
t845 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t870 + Ifges(4,4) * t873) * qJD(2);
t913 = (t844 * t870 - t845 * t873) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + mrSges(4,1) * t789 + mrSges(5,1) * t777 - mrSges(4,2) * t790 - mrSges(5,2) * t778 + Ifges(4,5) * t852 + Ifges(5,5) * t828 + Ifges(4,6) * t853 + Ifges(5,6) * t827 + pkin(3) * t745 + pkin(4) * t764 + pkin(9) * t756 + t872 * t746 + t869 * t753 + t841 * t817 + t840 * t818;
t912 = mrSges(6,1) * t771 + mrSges(7,1) * t768 - mrSges(6,2) * t772 - mrSges(7,2) * t770 + pkin(5) * t765 + t908 * t802 + t909 * t803 + t826 * t914 - t899 * t830 - t900 * t831;
t755 = t872 * t758 + t869 * t763;
t754 = m(5) * t787 - t827 * mrSges(5,1) + mrSges(5,2) * t828 + t840 * t834 + t835 * t841 + t755;
t814 = -t876 * pkin(8) + t880;
t857 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t897;
t896 = qJD(2) * t873;
t858 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t896;
t878 = -m(4) * t814 + t853 * mrSges(4,1) - mrSges(4,2) * t852 - t857 * t897 + t858 * t896 - t754;
t749 = m(3) * t819 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t876 + t878;
t905 = t749 * t874;
t851 = (-mrSges(4,1) * t873 + mrSges(4,2) * t870) * qJD(2);
t743 = m(4) * t789 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t852 + qJD(3) * t858 - t851 * t897 + t745;
t888 = t866 * t752 - t760 * t863;
t744 = m(4) * t790 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t853 - qJD(3) * t857 + t851 * t896 + t888;
t889 = -t743 * t870 + t873 * t744;
t735 = m(3) * t820 - mrSges(3,1) * t876 - qJDD(2) * mrSges(3,2) + t889;
t738 = t873 * t743 + t870 * t744;
t737 = m(3) * t836 + t738;
t725 = t735 * t903 - t737 * t865 + t868 * t905;
t723 = m(2) * t854 + t725;
t730 = t874 * t735 - t749 * t871;
t729 = m(2) * t855 + t730;
t902 = t867 * t723 + t864 * t729;
t724 = t735 * t904 + t868 * t737 + t865 * t905;
t890 = -t723 * t864 + t867 * t729;
t886 = m(2) * t862 + t724;
t816 = Ifges(5,5) * t841 - Ifges(5,6) * t840 + Ifges(5,3) * qJD(3);
t731 = mrSges(5,2) * t787 - mrSges(5,3) * t777 + Ifges(5,1) * t828 + Ifges(5,4) * t827 + Ifges(5,5) * qJDD(3) - pkin(9) * t755 - qJD(3) * t817 - t746 * t869 + t753 * t872 - t816 * t840;
t739 = -mrSges(5,1) * t787 + mrSges(5,3) * t778 + Ifges(5,4) * t828 + Ifges(5,2) * t827 + Ifges(5,6) * qJDD(3) - pkin(4) * t755 + qJD(3) * t818 - t841 * t816 - t912;
t843 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t870 + Ifges(4,6) * t873) * qJD(2);
t721 = -mrSges(4,1) * t814 + mrSges(4,3) * t790 + Ifges(4,4) * t852 + Ifges(4,2) * t853 + Ifges(4,6) * qJDD(3) - pkin(3) * t754 + qJ(4) * t888 + qJD(3) * t845 + t863 * t731 + t866 * t739 - t843 * t897;
t726 = mrSges(4,2) * t814 - mrSges(4,3) * t789 + Ifges(4,1) * t852 + Ifges(4,4) * t853 + Ifges(4,5) * qJDD(3) - qJ(4) * t745 - qJD(3) * t844 + t731 * t866 - t739 * t863 + t843 * t896;
t719 = mrSges(3,2) * t836 - mrSges(3,3) * t819 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t876 - pkin(8) * t738 - t721 * t870 + t726 * t873;
t720 = -mrSges(3,1) * t836 + mrSges(3,3) * t820 + t876 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t738 - t913;
t881 = pkin(7) * t730 + t719 * t871 + t720 * t874;
t718 = mrSges(3,1) * t819 - mrSges(3,2) * t820 + Ifges(3,3) * qJDD(2) + pkin(2) * t878 + pkin(8) * t889 + t873 * t721 + t870 * t726;
t717 = mrSges(2,2) * t862 - mrSges(2,3) * t854 + t874 * t719 - t871 * t720 + (-t724 * t865 - t725 * t868) * pkin(7);
t716 = -mrSges(2,1) * t862 + mrSges(2,3) * t855 - pkin(1) * t724 - t865 * t718 + t881 * t868;
t1 = [-m(1) * g(1) + t890; -m(1) * g(2) + t902; -m(1) * g(3) + t886; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t902 - t864 * t716 + t867 * t717; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t890 + t867 * t716 + t864 * t717; -mrSges(1,1) * g(2) + mrSges(2,1) * t854 + mrSges(1,2) * g(1) - mrSges(2,2) * t855 + pkin(1) * t725 + t868 * t718 + t881 * t865; t886; t718; t913; t754; t912; t766;];
tauJB  = t1;
