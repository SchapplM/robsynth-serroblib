% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRP7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 17:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRP7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP7_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:57:13
% EndTime: 2019-05-05 17:57:21
% DurationCPUTime: 5.53s
% Computational Cost: add. (59264->316), mult. (128611->373), div. (0->0), fcn. (82393->8), ass. (0->127)
t898 = -2 * qJD(4);
t897 = Ifges(6,1) + Ifges(7,1);
t890 = Ifges(6,4) + Ifges(7,4);
t888 = Ifges(6,5) + Ifges(7,5);
t896 = Ifges(6,2) + Ifges(7,2);
t886 = Ifges(6,6) + Ifges(7,6);
t895 = Ifges(6,3) + Ifges(7,3);
t849 = sin(qJ(1));
t852 = cos(qJ(1));
t830 = t849 * g(1) - t852 * g(2);
t854 = qJD(1) ^ 2;
t862 = -t854 * qJ(2) + qJDD(2) - t830;
t893 = -pkin(1) - pkin(7);
t803 = t893 * qJDD(1) + t862;
t848 = sin(qJ(3));
t851 = cos(qJ(3));
t792 = t848 * g(3) + t851 * t803;
t876 = qJD(1) * qJD(3);
t872 = t848 * t876;
t825 = t851 * qJDD(1) - t872;
t765 = (-t825 - t872) * qJ(4) + (-t848 * t851 * t854 + qJDD(3)) * pkin(3) + t792;
t793 = -t851 * g(3) + t848 * t803;
t824 = -t848 * qJDD(1) - t851 * t876;
t878 = qJD(1) * t851;
t828 = qJD(3) * pkin(3) - qJ(4) * t878;
t842 = t848 ^ 2;
t766 = -t842 * t854 * pkin(3) + t824 * qJ(4) - qJD(3) * t828 + t793;
t845 = sin(pkin(9));
t846 = cos(pkin(9));
t879 = qJD(1) * t848;
t813 = -t845 * t879 + t846 * t878;
t745 = t846 * t765 - t845 * t766 + t813 * t898;
t831 = -t852 * g(1) - t849 * g(2);
t863 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t831;
t812 = (t845 * t851 + t846 * t848) * qJD(1);
t790 = t845 * t824 + t846 * t825;
t847 = sin(qJ(5));
t850 = cos(qJ(5));
t795 = t850 * qJD(3) - t847 * t813;
t760 = t795 * qJD(5) + t847 * qJDD(3) + t850 * t790;
t796 = t847 * qJD(3) + t850 * t813;
t770 = -t795 * mrSges(7,1) + t796 * mrSges(7,2);
t746 = t845 * t765 + t846 * t766 + t812 * t898;
t783 = t812 * pkin(4) - t813 * pkin(8);
t853 = qJD(3) ^ 2;
t740 = -t853 * pkin(4) + qJDD(3) * pkin(8) - t812 * t783 + t746;
t768 = -t824 * pkin(3) + qJDD(4) + t828 * t878 + (-qJ(4) * t842 + t893) * t854 + t863;
t789 = t846 * t824 - t845 * t825;
t743 = (qJD(3) * t812 - t790) * pkin(8) + (qJD(3) * t813 - t789) * pkin(4) + t768;
t736 = -t847 * t740 + t850 * t743;
t788 = qJDD(5) - t789;
t810 = qJD(5) + t812;
t732 = -0.2e1 * qJD(6) * t796 + (t795 * t810 - t760) * qJ(6) + (t795 * t796 + t788) * pkin(5) + t736;
t773 = -t810 * mrSges(7,2) + t795 * mrSges(7,3);
t874 = m(7) * t732 + t788 * mrSges(7,1) + t810 * t773;
t729 = -t760 * mrSges(7,3) - t796 * t770 + t874;
t737 = t850 * t740 + t847 * t743;
t759 = -t796 * qJD(5) + t850 * qJDD(3) - t847 * t790;
t775 = t810 * pkin(5) - t796 * qJ(6);
t794 = t795 ^ 2;
t734 = -t794 * pkin(5) + t759 * qJ(6) + 0.2e1 * qJD(6) * t795 - t810 * t775 + t737;
t881 = t890 * t795 + t897 * t796 + t888 * t810;
t882 = -t896 * t795 - t890 * t796 - t886 * t810;
t894 = mrSges(6,1) * t736 + mrSges(7,1) * t732 - mrSges(6,2) * t737 - mrSges(7,2) * t734 + pkin(5) * t729 + t886 * t759 + t888 * t760 + t895 * t788 - t881 * t795 - t882 * t796;
t892 = mrSges(2,1) - mrSges(3,2);
t891 = -mrSges(6,2) - mrSges(7,2);
t889 = Ifges(2,5) - Ifges(3,4);
t887 = -Ifges(2,6) + Ifges(3,5);
t782 = t812 * mrSges(5,1) + t813 * mrSges(5,2);
t802 = qJD(3) * mrSges(5,1) - t813 * mrSges(5,3);
t771 = -t795 * mrSges(6,1) + t796 * mrSges(6,2);
t774 = -t810 * mrSges(6,2) + t795 * mrSges(6,3);
t722 = m(6) * t736 + t788 * mrSges(6,1) + t810 * t774 + (-t770 - t771) * t796 + (-mrSges(6,3) - mrSges(7,3)) * t760 + t874;
t873 = m(7) * t734 + t759 * mrSges(7,3) + t795 * t770;
t776 = t810 * mrSges(7,1) - t796 * mrSges(7,3);
t880 = -t810 * mrSges(6,1) + t796 * mrSges(6,3) - t776;
t725 = m(6) * t737 + t759 * mrSges(6,3) + t795 * t771 + t891 * t788 + t880 * t810 + t873;
t868 = -t847 * t722 + t850 * t725;
t715 = m(5) * t746 - qJDD(3) * mrSges(5,2) + t789 * mrSges(5,3) - qJD(3) * t802 - t812 * t782 + t868;
t801 = -qJD(3) * mrSges(5,2) - t812 * mrSges(5,3);
t739 = -qJDD(3) * pkin(4) - t853 * pkin(8) + t813 * t783 - t745;
t735 = -t759 * pkin(5) - t794 * qJ(6) + t796 * t775 + qJDD(6) + t739;
t866 = -m(7) * t735 + t759 * mrSges(7,1) + t795 * t773;
t858 = -m(6) * t739 + t759 * mrSges(6,1) + t891 * t760 + t795 * t774 + t880 * t796 + t866;
t727 = m(5) * t745 + qJDD(3) * mrSges(5,1) - t790 * mrSges(5,3) + qJD(3) * t801 - t813 * t782 + t858;
t706 = t845 * t715 + t846 * t727;
t823 = (mrSges(4,1) * t848 + mrSges(4,2) * t851) * qJD(1);
t827 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t879;
t703 = m(4) * t792 + qJDD(3) * mrSges(4,1) - t825 * mrSges(4,3) + qJD(3) * t827 - t823 * t878 + t706;
t829 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t878;
t869 = t846 * t715 - t845 * t727;
t704 = m(4) * t793 - qJDD(3) * mrSges(4,2) + t824 * mrSges(4,3) - qJD(3) * t829 - t823 * t879 + t869;
t700 = t851 * t703 + t848 * t704;
t811 = -qJDD(1) * pkin(1) + t862;
t861 = -m(3) * t811 + t854 * mrSges(3,3) - t700;
t696 = m(2) * t830 - t854 * mrSges(2,2) + t892 * qJDD(1) + t861;
t806 = t854 * pkin(1) - t863;
t720 = t850 * t722 + t847 * t725;
t718 = m(5) * t768 - t789 * mrSges(5,1) + t790 * mrSges(5,2) + t812 * t801 + t813 * t802 + t720;
t800 = t893 * t854 + t863;
t859 = -m(4) * t800 + t824 * mrSges(4,1) - t825 * mrSges(4,2) - t827 * t879 - t829 * t878 - t718;
t856 = -m(3) * t806 + t854 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t859;
t711 = m(2) * t831 - t854 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t856;
t884 = t852 * t696 + t849 * t711;
t883 = -t886 * t795 - t888 * t796 - t895 * t810;
t871 = -t849 * t696 + t852 * t711;
t870 = -t848 * t703 + t851 * t704;
t730 = t760 * mrSges(7,2) + t796 * t776 - t866;
t708 = -mrSges(6,1) * t739 + mrSges(6,3) * t737 - mrSges(7,1) * t735 + mrSges(7,3) * t734 - pkin(5) * t730 + qJ(6) * t873 + (-qJ(6) * t776 + t881) * t810 + t883 * t796 + (-qJ(6) * mrSges(7,2) + t886) * t788 + t890 * t760 + t896 * t759;
t717 = mrSges(6,2) * t739 + mrSges(7,2) * t735 - mrSges(6,3) * t736 - mrSges(7,3) * t732 - qJ(6) * t729 + t890 * t759 + t897 * t760 + t888 * t788 - t883 * t795 + t882 * t810;
t778 = Ifges(5,5) * t813 - Ifges(5,6) * t812 + Ifges(5,3) * qJD(3);
t779 = Ifges(5,4) * t813 - Ifges(5,2) * t812 + Ifges(5,6) * qJD(3);
t694 = mrSges(5,2) * t768 - mrSges(5,3) * t745 + Ifges(5,1) * t790 + Ifges(5,4) * t789 + Ifges(5,5) * qJDD(3) - pkin(8) * t720 - qJD(3) * t779 - t847 * t708 + t850 * t717 - t812 * t778;
t780 = Ifges(5,1) * t813 - Ifges(5,4) * t812 + Ifges(5,5) * qJD(3);
t701 = -mrSges(5,1) * t768 + mrSges(5,3) * t746 + Ifges(5,4) * t790 + Ifges(5,2) * t789 + Ifges(5,6) * qJDD(3) - pkin(4) * t720 + qJD(3) * t780 - t813 * t778 - t894;
t814 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t851 - Ifges(4,6) * t848) * qJD(1);
t816 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t851 - Ifges(4,4) * t848) * qJD(1);
t691 = -mrSges(4,1) * t800 + mrSges(4,3) * t793 + Ifges(4,4) * t825 + Ifges(4,2) * t824 + Ifges(4,6) * qJDD(3) - pkin(3) * t718 + qJ(4) * t869 + qJD(3) * t816 + t845 * t694 + t846 * t701 - t814 * t878;
t815 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t851 - Ifges(4,2) * t848) * qJD(1);
t693 = mrSges(4,2) * t800 - mrSges(4,3) * t792 + Ifges(4,1) * t825 + Ifges(4,4) * t824 + Ifges(4,5) * qJDD(3) - qJ(4) * t706 - qJD(3) * t815 + t846 * t694 - t845 * t701 - t814 * t879;
t698 = qJDD(1) * mrSges(3,2) - t861;
t857 = mrSges(2,1) * t830 - mrSges(2,2) * t831 + mrSges(3,2) * t811 - mrSges(3,3) * t806 - pkin(1) * t698 - pkin(7) * t700 + qJ(2) * t856 - t848 * t691 + t851 * t693 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t855 = mrSges(5,1) * t745 - mrSges(4,2) * t793 - mrSges(5,2) * t746 + Ifges(5,5) * t790 + Ifges(5,6) * t789 + pkin(3) * t706 + pkin(4) * t858 + pkin(8) * t868 + t850 * t708 + t847 * t717 + t813 * t779 + t812 * t780 + mrSges(4,1) * t792 + t816 * t879 + t815 * t878 + Ifges(4,6) * t824 + Ifges(4,5) * t825 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3);
t699 = -m(3) * g(3) + t870;
t690 = t855 - mrSges(2,3) * t830 + mrSges(3,1) * t811 + pkin(2) * t700 + (-mrSges(2,2) + mrSges(3,3)) * g(3) - qJ(2) * t699 + t889 * qJDD(1) + t887 * t854;
t689 = -mrSges(3,1) * t806 + mrSges(2,3) * t831 - pkin(1) * t699 - pkin(2) * t859 - pkin(7) * t870 + t892 * g(3) - t887 * qJDD(1) - t851 * t691 - t848 * t693 + t889 * t854;
t1 = [-m(1) * g(1) + t871; -m(1) * g(2) + t884; (-m(1) - m(2) - m(3)) * g(3) + t870; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t884 - t849 * t689 + t852 * t690; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t871 + t852 * t689 + t849 * t690; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t857; t857; t698; t855; t718; t894; t730;];
tauJB  = t1;
