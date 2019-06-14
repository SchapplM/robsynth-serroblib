% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPPRRR2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:21:05
% EndTime: 2019-05-05 15:21:18
% DurationCPUTime: 13.06s
% Computational Cost: add. (206740->324), mult. (462602->403), div. (0->0), fcn. (331636->12), ass. (0->144)
t865 = qJD(1) ^ 2;
t854 = cos(pkin(11));
t894 = pkin(3) * t854;
t852 = sin(pkin(11));
t893 = mrSges(4,2) * t852;
t848 = t854 ^ 2;
t892 = t848 * t865;
t859 = sin(qJ(1));
t863 = cos(qJ(1));
t834 = t859 * g(1) - g(2) * t863;
t831 = qJDD(1) * pkin(1) + t834;
t835 = -g(1) * t863 - g(2) * t859;
t832 = -pkin(1) * t865 + t835;
t853 = sin(pkin(10));
t855 = cos(pkin(10));
t815 = t853 * t831 + t855 * t832;
t804 = -pkin(2) * t865 + qJDD(1) * qJ(3) + t815;
t851 = -g(3) + qJDD(2);
t886 = qJD(1) * qJD(3);
t890 = t854 * t851 - 0.2e1 * t852 * t886;
t784 = (-pkin(7) * qJDD(1) + t865 * t894 - t804) * t852 + t890;
t792 = t852 * t851 + (t804 + 0.2e1 * t886) * t854;
t884 = qJDD(1) * t854;
t787 = -pkin(3) * t892 + pkin(7) * t884 + t792;
t858 = sin(qJ(4));
t862 = cos(qJ(4));
t768 = t858 * t784 + t862 * t787;
t888 = qJD(1) * t854;
t889 = qJD(1) * t852;
t824 = -t858 * t889 + t862 * t888;
t874 = t852 * t862 + t854 * t858;
t825 = t874 * qJD(1);
t810 = -pkin(4) * t824 - pkin(8) * t825;
t864 = qJD(4) ^ 2;
t762 = -pkin(4) * t864 + qJDD(4) * pkin(8) + t810 * t824 + t768;
t847 = t852 ^ 2;
t814 = t831 * t855 - t853 * t832;
t875 = qJDD(3) - t814;
t788 = (-pkin(2) - t894) * qJDD(1) + (-qJ(3) + (-t847 - t848) * pkin(7)) * t865 + t875;
t822 = t825 * qJD(4);
t885 = qJDD(1) * t852;
t811 = -t858 * t885 + t862 * t884 - t822;
t887 = qJD(4) * t824;
t812 = t874 * qJDD(1) + t887;
t766 = (-t812 - t887) * pkin(8) + (-t811 + t822) * pkin(4) + t788;
t857 = sin(qJ(5));
t861 = cos(qJ(5));
t752 = -t762 * t857 + t861 * t766;
t817 = qJD(4) * t861 - t825 * t857;
t786 = qJD(5) * t817 + qJDD(4) * t857 + t812 * t861;
t809 = qJDD(5) - t811;
t818 = qJD(4) * t857 + t825 * t861;
t823 = qJD(5) - t824;
t750 = (t817 * t823 - t786) * pkin(9) + (t817 * t818 + t809) * pkin(5) + t752;
t753 = t861 * t762 + t857 * t766;
t785 = -qJD(5) * t818 + qJDD(4) * t861 - t812 * t857;
t797 = pkin(5) * t823 - pkin(9) * t818;
t816 = t817 ^ 2;
t751 = -pkin(5) * t816 + pkin(9) * t785 - t797 * t823 + t753;
t856 = sin(qJ(6));
t860 = cos(qJ(6));
t748 = t750 * t860 - t751 * t856;
t789 = t817 * t860 - t818 * t856;
t760 = qJD(6) * t789 + t785 * t856 + t786 * t860;
t790 = t817 * t856 + t818 * t860;
t773 = -mrSges(7,1) * t789 + mrSges(7,2) * t790;
t821 = qJD(6) + t823;
t774 = -mrSges(7,2) * t821 + mrSges(7,3) * t789;
t806 = qJDD(6) + t809;
t744 = m(7) * t748 + mrSges(7,1) * t806 - mrSges(7,3) * t760 - t773 * t790 + t774 * t821;
t749 = t750 * t856 + t751 * t860;
t759 = -qJD(6) * t790 + t785 * t860 - t786 * t856;
t775 = mrSges(7,1) * t821 - mrSges(7,3) * t790;
t745 = m(7) * t749 - mrSges(7,2) * t806 + mrSges(7,3) * t759 + t773 * t789 - t775 * t821;
t736 = t860 * t744 + t856 * t745;
t793 = -mrSges(6,1) * t817 + mrSges(6,2) * t818;
t795 = -mrSges(6,2) * t823 + mrSges(6,3) * t817;
t734 = m(6) * t752 + mrSges(6,1) * t809 - mrSges(6,3) * t786 - t793 * t818 + t795 * t823 + t736;
t796 = mrSges(6,1) * t823 - mrSges(6,3) * t818;
t879 = -t744 * t856 + t860 * t745;
t735 = m(6) * t753 - mrSges(6,2) * t809 + mrSges(6,3) * t785 + t793 * t817 - t796 * t823 + t879;
t730 = -t734 * t857 + t861 * t735;
t807 = -mrSges(5,1) * t824 + mrSges(5,2) * t825;
t820 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t825;
t728 = m(5) * t768 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t811 - qJD(4) * t820 + t807 * t824 + t730;
t767 = t784 * t862 - t858 * t787;
t761 = -qJDD(4) * pkin(4) - pkin(8) * t864 + t825 * t810 - t767;
t754 = -pkin(5) * t785 - pkin(9) * t816 + t797 * t818 + t761;
t872 = m(7) * t754 - t759 * mrSges(7,1) + mrSges(7,2) * t760 - t789 * t774 + t775 * t790;
t746 = -m(6) * t761 + t785 * mrSges(6,1) - mrSges(6,2) * t786 + t817 * t795 - t796 * t818 - t872;
t819 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t824;
t740 = m(5) * t767 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t812 + qJD(4) * t819 - t807 * t825 + t746;
t719 = t858 * t728 + t862 * t740;
t791 = -t804 * t852 + t890;
t873 = mrSges(4,3) * qJDD(1) + t865 * (-mrSges(4,1) * t854 + t893);
t717 = m(4) * t791 - t873 * t852 + t719;
t880 = t862 * t728 - t858 * t740;
t718 = m(4) * t792 + t873 * t854 + t880;
t881 = -t717 * t852 + t854 * t718;
t709 = m(3) * t815 - mrSges(3,1) * t865 - qJDD(1) * mrSges(3,2) + t881;
t799 = -qJDD(1) * pkin(2) - qJ(3) * t865 + t875;
t729 = t861 * t734 + t857 * t735;
t870 = m(5) * t788 - t811 * mrSges(5,1) + t812 * mrSges(5,2) - t824 * t819 + t825 * t820 + t729;
t868 = -m(4) * t799 + mrSges(4,1) * t884 - t870 + (t847 * t865 + t892) * mrSges(4,3);
t723 = t868 + m(3) * t814 - mrSges(3,2) * t865 + (mrSges(3,1) - t893) * qJDD(1);
t705 = t853 * t709 + t855 * t723;
t702 = m(2) * t834 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t865 + t705;
t882 = t855 * t709 - t723 * t853;
t703 = m(2) * t835 - mrSges(2,1) * t865 - qJDD(1) * mrSges(2,2) + t882;
t891 = t863 * t702 + t859 * t703;
t712 = t854 * t717 + t852 * t718;
t710 = m(3) * t851 + t712;
t883 = -t702 * t859 + t863 * t703;
t878 = Ifges(4,1) * t852 + Ifges(4,4) * t854;
t877 = Ifges(4,4) * t852 + Ifges(4,2) * t854;
t876 = Ifges(4,5) * t852 + Ifges(4,6) * t854;
t770 = Ifges(7,4) * t790 + Ifges(7,2) * t789 + Ifges(7,6) * t821;
t771 = Ifges(7,1) * t790 + Ifges(7,4) * t789 + Ifges(7,5) * t821;
t871 = -mrSges(7,1) * t748 + mrSges(7,2) * t749 - Ifges(7,5) * t760 - Ifges(7,6) * t759 - Ifges(7,3) * t806 - t790 * t770 + t789 * t771;
t769 = Ifges(7,5) * t790 + Ifges(7,6) * t789 + Ifges(7,3) * t821;
t737 = -mrSges(7,1) * t754 + mrSges(7,3) * t749 + Ifges(7,4) * t760 + Ifges(7,2) * t759 + Ifges(7,6) * t806 - t769 * t790 + t771 * t821;
t738 = mrSges(7,2) * t754 - mrSges(7,3) * t748 + Ifges(7,1) * t760 + Ifges(7,4) * t759 + Ifges(7,5) * t806 + t769 * t789 - t770 * t821;
t777 = Ifges(6,5) * t818 + Ifges(6,6) * t817 + Ifges(6,3) * t823;
t779 = Ifges(6,1) * t818 + Ifges(6,4) * t817 + Ifges(6,5) * t823;
t720 = -mrSges(6,1) * t761 + mrSges(6,3) * t753 + Ifges(6,4) * t786 + Ifges(6,2) * t785 + Ifges(6,6) * t809 - pkin(5) * t872 + pkin(9) * t879 + t860 * t737 + t856 * t738 - t818 * t777 + t823 * t779;
t778 = Ifges(6,4) * t818 + Ifges(6,2) * t817 + Ifges(6,6) * t823;
t721 = mrSges(6,2) * t761 - mrSges(6,3) * t752 + Ifges(6,1) * t786 + Ifges(6,4) * t785 + Ifges(6,5) * t809 - pkin(9) * t736 - t737 * t856 + t738 * t860 + t777 * t817 - t778 * t823;
t800 = Ifges(5,5) * t825 + Ifges(5,6) * t824 + Ifges(5,3) * qJD(4);
t801 = Ifges(5,4) * t825 + Ifges(5,2) * t824 + Ifges(5,6) * qJD(4);
t706 = mrSges(5,2) * t788 - mrSges(5,3) * t767 + Ifges(5,1) * t812 + Ifges(5,4) * t811 + Ifges(5,5) * qJDD(4) - pkin(8) * t729 - qJD(4) * t801 - t720 * t857 + t721 * t861 + t800 * t824;
t802 = Ifges(5,1) * t825 + Ifges(5,4) * t824 + Ifges(5,5) * qJD(4);
t866 = mrSges(6,1) * t752 - mrSges(6,2) * t753 + Ifges(6,5) * t786 + Ifges(6,6) * t785 + Ifges(6,3) * t809 + pkin(5) * t736 + t818 * t778 - t817 * t779 - t871;
t713 = -mrSges(5,1) * t788 + mrSges(5,3) * t768 + Ifges(5,4) * t812 + Ifges(5,2) * t811 + Ifges(5,6) * qJDD(4) - pkin(4) * t729 + qJD(4) * t802 - t825 * t800 - t866;
t830 = t876 * qJD(1);
t695 = -mrSges(4,1) * t799 + mrSges(4,3) * t792 - pkin(3) * t870 + pkin(7) * t880 + t877 * qJDD(1) + t858 * t706 + t862 * t713 - t830 * t889;
t698 = mrSges(4,2) * t799 - mrSges(4,3) * t791 - pkin(7) * t719 + t878 * qJDD(1) + t706 * t862 - t713 * t858 + t830 * t888;
t725 = mrSges(4,2) * t885 - t868;
t869 = mrSges(2,1) * t834 + mrSges(3,1) * t814 - mrSges(2,2) * t835 - mrSges(3,2) * t815 + pkin(1) * t705 - pkin(2) * t725 + qJ(3) * t881 + t854 * t695 + t852 * t698 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t867 = mrSges(5,1) * t767 - mrSges(5,2) * t768 + Ifges(5,5) * t812 + Ifges(5,6) * t811 + Ifges(5,3) * qJDD(4) + pkin(4) * t746 + pkin(8) * t730 + t861 * t720 + t857 * t721 + t825 * t801 - t824 * t802;
t696 = -pkin(2) * t712 - t867 - pkin(3) * t719 + (Ifges(3,6) - t876) * qJDD(1) - mrSges(3,1) * t851 + mrSges(3,3) * t815 - mrSges(4,1) * t791 + mrSges(4,2) * t792 + (-t852 * t877 + t854 * t878 + Ifges(3,5)) * t865;
t693 = mrSges(3,2) * t851 - mrSges(3,3) * t814 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t865 - qJ(3) * t712 - t695 * t852 + t698 * t854;
t692 = -mrSges(2,2) * g(3) - mrSges(2,3) * t834 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t865 - qJ(2) * t705 + t693 * t855 - t696 * t853;
t691 = mrSges(2,1) * g(3) + mrSges(2,3) * t835 + t865 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t710 + qJ(2) * t882 + t853 * t693 + t855 * t696;
t1 = [-m(1) * g(1) + t883; -m(1) * g(2) + t891; (-m(1) - m(2)) * g(3) + t710; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t891 - t859 * t691 + t863 * t692; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t883 + t863 * t691 + t859 * t692; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t869; t869; t710; t725; t867; t866; -t871;];
tauJB  = t1;
