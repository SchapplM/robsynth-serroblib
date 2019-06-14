% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRRPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRRPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:15:37
% EndTime: 2019-05-05 04:16:02
% DurationCPUTime: 23.61s
% Computational Cost: add. (389772->342), mult. (855372->443), div. (0->0), fcn. (637037->14), ass. (0->145)
t874 = sin(pkin(11));
t877 = cos(pkin(11));
t860 = g(1) * t874 - g(2) * t877;
t861 = -g(1) * t877 - g(2) * t874;
t872 = -g(3) + qJDD(1);
t875 = sin(pkin(6));
t878 = cos(pkin(6));
t882 = sin(qJ(2));
t886 = cos(qJ(2));
t829 = -t882 * t861 + (t860 * t878 + t872 * t875) * t886;
t909 = t878 * t882;
t910 = t875 * t882;
t830 = t860 * t909 + t886 * t861 + t872 * t910;
t887 = qJD(2) ^ 2;
t822 = -pkin(2) * t887 + qJDD(2) * pkin(8) + t830;
t842 = -t860 * t875 + t872 * t878;
t881 = sin(qJ(3));
t885 = cos(qJ(3));
t814 = -t881 * t822 + t885 * t842;
t905 = qJD(2) * qJD(3);
t904 = t885 * t905;
t858 = qJDD(2) * t881 + t904;
t801 = (-t858 + t904) * qJ(4) + (t881 * t885 * t887 + qJDD(3)) * pkin(3) + t814;
t815 = t885 * t822 + t881 * t842;
t859 = qJDD(2) * t885 - t881 * t905;
t907 = qJD(2) * t881;
t862 = qJD(3) * pkin(3) - qJ(4) * t907;
t871 = t885 ^ 2;
t803 = -pkin(3) * t871 * t887 + qJ(4) * t859 - qJD(3) * t862 + t815;
t873 = sin(pkin(12));
t876 = cos(pkin(12));
t847 = (t873 * t885 + t876 * t881) * qJD(2);
t779 = -0.2e1 * qJD(4) * t847 + t876 * t801 - t873 * t803;
t836 = t858 * t876 + t859 * t873;
t846 = (-t873 * t881 + t876 * t885) * qJD(2);
t776 = (qJD(3) * t846 - t836) * pkin(9) + (t846 * t847 + qJDD(3)) * pkin(4) + t779;
t780 = 0.2e1 * qJD(4) * t846 + t873 * t801 + t876 * t803;
t835 = -t858 * t873 + t859 * t876;
t841 = qJD(3) * pkin(4) - pkin(9) * t847;
t845 = t846 ^ 2;
t778 = -pkin(4) * t845 + pkin(9) * t835 - qJD(3) * t841 + t780;
t880 = sin(qJ(5));
t884 = cos(qJ(5));
t773 = t880 * t776 + t884 * t778;
t828 = t846 * t880 + t847 * t884;
t796 = -qJD(5) * t828 + t835 * t884 - t836 * t880;
t827 = t846 * t884 - t847 * t880;
t811 = -mrSges(6,1) * t827 + mrSges(6,2) * t828;
t870 = qJD(3) + qJD(5);
t820 = mrSges(6,1) * t870 - mrSges(6,3) * t828;
t869 = qJDD(3) + qJDD(5);
t812 = -pkin(5) * t827 - pkin(10) * t828;
t868 = t870 ^ 2;
t770 = -pkin(5) * t868 + pkin(10) * t869 + t812 * t827 + t773;
t892 = -qJDD(2) * pkin(2) - t829;
t813 = -t859 * pkin(3) + qJDD(4) + t862 * t907 + (-qJ(4) * t871 - pkin(8)) * t887 + t892;
t785 = -t835 * pkin(4) - t845 * pkin(9) + t847 * t841 + t813;
t797 = qJD(5) * t827 + t835 * t880 + t836 * t884;
t774 = (-t827 * t870 - t797) * pkin(10) + (t828 * t870 - t796) * pkin(5) + t785;
t879 = sin(qJ(6));
t883 = cos(qJ(6));
t767 = -t770 * t879 + t774 * t883;
t816 = -t828 * t879 + t870 * t883;
t783 = qJD(6) * t816 + t797 * t883 + t869 * t879;
t795 = qJDD(6) - t796;
t817 = t828 * t883 + t870 * t879;
t802 = -mrSges(7,1) * t816 + mrSges(7,2) * t817;
t823 = qJD(6) - t827;
t804 = -mrSges(7,2) * t823 + mrSges(7,3) * t816;
t763 = m(7) * t767 + mrSges(7,1) * t795 - mrSges(7,3) * t783 - t802 * t817 + t804 * t823;
t768 = t770 * t883 + t774 * t879;
t782 = -qJD(6) * t817 - t797 * t879 + t869 * t883;
t805 = mrSges(7,1) * t823 - mrSges(7,3) * t817;
t764 = m(7) * t768 - mrSges(7,2) * t795 + mrSges(7,3) * t782 + t802 * t816 - t805 * t823;
t899 = -t763 * t879 + t883 * t764;
t749 = m(6) * t773 - mrSges(6,2) * t869 + mrSges(6,3) * t796 + t811 * t827 - t820 * t870 + t899;
t772 = t776 * t884 - t778 * t880;
t819 = -mrSges(6,2) * t870 + mrSges(6,3) * t827;
t769 = -pkin(5) * t869 - pkin(10) * t868 + t812 * t828 - t772;
t893 = -m(7) * t769 + t782 * mrSges(7,1) - mrSges(7,2) * t783 + t816 * t804 - t805 * t817;
t759 = m(6) * t772 + mrSges(6,1) * t869 - mrSges(6,3) * t797 - t811 * t828 + t819 * t870 + t893;
t743 = t880 * t749 + t884 * t759;
t833 = -mrSges(5,1) * t846 + mrSges(5,2) * t847;
t839 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t846;
t741 = m(5) * t779 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t836 + qJD(3) * t839 - t833 * t847 + t743;
t840 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t847;
t900 = t884 * t749 - t759 * t880;
t742 = m(5) * t780 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t835 - qJD(3) * t840 + t833 * t846 + t900;
t735 = t876 * t741 + t873 * t742;
t825 = Ifges(5,4) * t847 + Ifges(5,2) * t846 + Ifges(5,6) * qJD(3);
t826 = Ifges(5,1) * t847 + Ifges(5,4) * t846 + Ifges(5,5) * qJD(3);
t850 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t881 + Ifges(4,2) * t885) * qJD(2);
t851 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t881 + Ifges(4,4) * t885) * qJD(2);
t786 = Ifges(7,5) * t817 + Ifges(7,6) * t816 + Ifges(7,3) * t823;
t788 = Ifges(7,1) * t817 + Ifges(7,4) * t816 + Ifges(7,5) * t823;
t756 = -mrSges(7,1) * t769 + mrSges(7,3) * t768 + Ifges(7,4) * t783 + Ifges(7,2) * t782 + Ifges(7,6) * t795 - t786 * t817 + t788 * t823;
t787 = Ifges(7,4) * t817 + Ifges(7,2) * t816 + Ifges(7,6) * t823;
t757 = mrSges(7,2) * t769 - mrSges(7,3) * t767 + Ifges(7,1) * t783 + Ifges(7,4) * t782 + Ifges(7,5) * t795 + t786 * t816 - t787 * t823;
t807 = Ifges(6,4) * t828 + Ifges(6,2) * t827 + Ifges(6,6) * t870;
t808 = Ifges(6,1) * t828 + Ifges(6,4) * t827 + Ifges(6,5) * t870;
t891 = -mrSges(6,1) * t772 + mrSges(6,2) * t773 - Ifges(6,5) * t797 - Ifges(6,6) * t796 - Ifges(6,3) * t869 - pkin(5) * t893 - pkin(10) * t899 - t883 * t756 - t879 * t757 - t828 * t807 + t827 * t808;
t913 = mrSges(4,1) * t814 + mrSges(5,1) * t779 - mrSges(4,2) * t815 - mrSges(5,2) * t780 + Ifges(4,5) * t858 + Ifges(5,5) * t836 + Ifges(4,6) * t859 + Ifges(5,6) * t835 + pkin(3) * t735 + pkin(4) * t743 + (t850 * t881 - t851 * t885) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t847 * t825 - t846 * t826 - t891;
t752 = t883 * t763 + t879 * t764;
t895 = m(6) * t785 - t796 * mrSges(6,1) + t797 * mrSges(6,2) - t827 * t819 + t828 * t820 + t752;
t750 = m(5) * t813 - t835 * mrSges(5,1) + mrSges(5,2) * t836 - t846 * t839 + t840 * t847 + t895;
t821 = -t887 * pkin(8) + t892;
t863 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t907;
t906 = qJD(2) * t885;
t864 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t906;
t889 = -m(4) * t821 + t859 * mrSges(4,1) - mrSges(4,2) * t858 - t863 * t907 + t864 * t906 - t750;
t746 = m(3) * t829 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t887 + t889;
t911 = t746 * t886;
t857 = (-mrSges(4,1) * t885 + mrSges(4,2) * t881) * qJD(2);
t733 = m(4) * t814 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t858 + qJD(3) * t864 - t857 * t907 + t735;
t901 = -t741 * t873 + t876 * t742;
t734 = m(4) * t815 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t859 - qJD(3) * t863 + t857 * t906 + t901;
t902 = -t733 * t881 + t885 * t734;
t725 = m(3) * t830 - mrSges(3,1) * t887 - qJDD(2) * mrSges(3,2) + t902;
t728 = t885 * t733 + t881 * t734;
t727 = m(3) * t842 + t728;
t716 = t725 * t909 - t727 * t875 + t878 * t911;
t714 = m(2) * t860 + t716;
t720 = t886 * t725 - t746 * t882;
t719 = m(2) * t861 + t720;
t908 = t877 * t714 + t874 * t719;
t715 = t725 * t910 + t878 * t727 + t875 * t911;
t903 = -t714 * t874 + t877 * t719;
t898 = m(2) * t872 + t715;
t806 = Ifges(6,5) * t828 + Ifges(6,6) * t827 + Ifges(6,3) * t870;
t736 = mrSges(6,2) * t785 - mrSges(6,3) * t772 + Ifges(6,1) * t797 + Ifges(6,4) * t796 + Ifges(6,5) * t869 - pkin(10) * t752 - t756 * t879 + t757 * t883 + t806 * t827 - t807 * t870;
t890 = mrSges(7,1) * t767 - mrSges(7,2) * t768 + Ifges(7,5) * t783 + Ifges(7,6) * t782 + Ifges(7,3) * t795 + t787 * t817 - t788 * t816;
t737 = -mrSges(6,1) * t785 + mrSges(6,3) * t773 + Ifges(6,4) * t797 + Ifges(6,2) * t796 + Ifges(6,6) * t869 - pkin(5) * t752 - t806 * t828 + t808 * t870 - t890;
t824 = Ifges(5,5) * t847 + Ifges(5,6) * t846 + Ifges(5,3) * qJD(3);
t721 = -mrSges(5,1) * t813 + mrSges(5,3) * t780 + Ifges(5,4) * t836 + Ifges(5,2) * t835 + Ifges(5,6) * qJDD(3) - pkin(4) * t895 + pkin(9) * t900 + qJD(3) * t826 + t880 * t736 + t884 * t737 - t847 * t824;
t729 = mrSges(5,2) * t813 - mrSges(5,3) * t779 + Ifges(5,1) * t836 + Ifges(5,4) * t835 + Ifges(5,5) * qJDD(3) - pkin(9) * t743 - qJD(3) * t825 + t736 * t884 - t737 * t880 + t824 * t846;
t849 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t881 + Ifges(4,6) * t885) * qJD(2);
t710 = -mrSges(4,1) * t821 + mrSges(4,3) * t815 + Ifges(4,4) * t858 + Ifges(4,2) * t859 + Ifges(4,6) * qJDD(3) - pkin(3) * t750 + qJ(4) * t901 + qJD(3) * t851 + t876 * t721 + t873 * t729 - t849 * t907;
t711 = mrSges(4,2) * t821 - mrSges(4,3) * t814 + Ifges(4,1) * t858 + Ifges(4,4) * t859 + Ifges(4,5) * qJDD(3) - qJ(4) * t735 - qJD(3) * t850 - t721 * t873 + t729 * t876 + t849 * t906;
t709 = mrSges(3,2) * t842 - mrSges(3,3) * t829 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t887 - pkin(8) * t728 - t710 * t881 + t711 * t885;
t712 = -mrSges(3,1) * t842 + mrSges(3,3) * t830 + t887 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t728 - t913;
t894 = pkin(7) * t720 + t709 * t882 + t712 * t886;
t708 = mrSges(3,1) * t829 - mrSges(3,2) * t830 + Ifges(3,3) * qJDD(2) + pkin(2) * t889 + pkin(8) * t902 + t885 * t710 + t881 * t711;
t707 = mrSges(2,2) * t872 - mrSges(2,3) * t860 + t886 * t709 - t882 * t712 + (-t715 * t875 - t716 * t878) * pkin(7);
t706 = -mrSges(2,1) * t872 + mrSges(2,3) * t861 - pkin(1) * t715 - t875 * t708 + t878 * t894;
t1 = [-m(1) * g(1) + t903; -m(1) * g(2) + t908; -m(1) * g(3) + t898; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t908 - t874 * t706 + t877 * t707; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t903 + t877 * t706 + t874 * t707; -mrSges(1,1) * g(2) + mrSges(2,1) * t860 + mrSges(1,2) * g(1) - mrSges(2,2) * t861 + pkin(1) * t716 + t878 * t708 + t875 * t894; t898; t708; t913; t750; -t891; t890;];
tauJB  = t1;
