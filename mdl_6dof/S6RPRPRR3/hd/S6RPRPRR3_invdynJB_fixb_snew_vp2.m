% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6RPRPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:33:57
% EndTime: 2019-05-05 18:34:15
% DurationCPUTime: 17.57s
% Computational Cost: add. (288019->343), mult. (602105->429), div. (0->0), fcn. (409795->12), ass. (0->138)
t830 = sin(qJ(1));
t834 = cos(qJ(1));
t812 = t830 * g(1) - g(2) * t834;
t803 = qJDD(1) * pkin(1) + t812;
t813 = -g(1) * t834 - g(2) * t830;
t836 = qJD(1) ^ 2;
t806 = -pkin(1) * t836 + t813;
t824 = sin(pkin(10));
t826 = cos(pkin(10));
t778 = t826 * t803 - t824 * t806;
t767 = -qJDD(1) * pkin(2) - t836 * pkin(7) - t778;
t829 = sin(qJ(3));
t833 = cos(qJ(3));
t852 = qJD(1) * qJD(3);
t851 = t833 * t852;
t807 = qJDD(1) * t829 + t851;
t817 = t829 * t852;
t808 = qJDD(1) * t833 - t817;
t753 = (-t807 - t851) * qJ(4) + (-t808 + t817) * pkin(3) + t767;
t779 = t824 * t803 + t826 * t806;
t768 = -pkin(2) * t836 + qJDD(1) * pkin(7) + t779;
t822 = -g(3) + qJDD(2);
t761 = t833 * t768 + t829 * t822;
t804 = (-pkin(3) * t833 - qJ(4) * t829) * qJD(1);
t835 = qJD(3) ^ 2;
t853 = qJD(1) * t833;
t759 = -pkin(3) * t835 + qJDD(3) * qJ(4) + t804 * t853 + t761;
t823 = sin(pkin(11));
t825 = cos(pkin(11));
t854 = qJD(1) * t829;
t800 = qJD(3) * t823 + t825 * t854;
t734 = -0.2e1 * qJD(4) * t800 + t825 * t753 - t823 * t759;
t784 = qJDD(3) * t823 + t807 * t825;
t799 = qJD(3) * t825 - t823 * t854;
t728 = (-t799 * t853 - t784) * pkin(8) + (t799 * t800 - t808) * pkin(4) + t734;
t735 = 0.2e1 * qJD(4) * t799 + t823 * t753 + t825 * t759;
t783 = qJDD(3) * t825 - t807 * t823;
t785 = -pkin(4) * t853 - pkin(8) * t800;
t798 = t799 ^ 2;
t733 = -pkin(4) * t798 + pkin(8) * t783 + t785 * t853 + t735;
t828 = sin(qJ(5));
t832 = cos(qJ(5));
t720 = t832 * t728 - t828 * t733;
t775 = t799 * t832 - t800 * t828;
t746 = qJD(5) * t775 + t783 * t828 + t784 * t832;
t776 = t799 * t828 + t800 * t832;
t802 = qJDD(5) - t808;
t815 = qJD(5) - t853;
t718 = (t775 * t815 - t746) * pkin(9) + (t775 * t776 + t802) * pkin(5) + t720;
t721 = t828 * t728 + t832 * t733;
t745 = -qJD(5) * t776 + t783 * t832 - t784 * t828;
t764 = pkin(5) * t815 - pkin(9) * t776;
t774 = t775 ^ 2;
t719 = -pkin(5) * t774 + pkin(9) * t745 - t764 * t815 + t721;
t827 = sin(qJ(6));
t831 = cos(qJ(6));
t717 = t718 * t827 + t719 * t831;
t760 = -t829 * t768 + t822 * t833;
t757 = -qJDD(3) * pkin(3) - qJ(4) * t835 + t804 * t854 + qJDD(4) - t760;
t741 = -pkin(4) * t783 - pkin(8) * t798 + t800 * t785 + t757;
t722 = -pkin(5) * t745 - pkin(9) * t774 + t764 * t776 + t741;
t756 = t775 * t827 + t776 * t831;
t729 = -qJD(6) * t756 + t745 * t831 - t746 * t827;
t755 = t775 * t831 - t776 * t827;
t730 = qJD(6) * t755 + t745 * t827 + t746 * t831;
t814 = qJD(6) + t815;
t736 = Ifges(7,5) * t756 + Ifges(7,6) * t755 + Ifges(7,3) * t814;
t738 = Ifges(7,1) * t756 + Ifges(7,4) * t755 + Ifges(7,5) * t814;
t796 = qJDD(6) + t802;
t705 = -mrSges(7,1) * t722 + mrSges(7,3) * t717 + Ifges(7,4) * t730 + Ifges(7,2) * t729 + Ifges(7,6) * t796 - t736 * t756 + t738 * t814;
t716 = t718 * t831 - t719 * t827;
t737 = Ifges(7,4) * t756 + Ifges(7,2) * t755 + Ifges(7,6) * t814;
t706 = mrSges(7,2) * t722 - mrSges(7,3) * t716 + Ifges(7,1) * t730 + Ifges(7,4) * t729 + Ifges(7,5) * t796 + t736 * t755 - t737 * t814;
t749 = Ifges(6,5) * t776 + Ifges(6,6) * t775 + Ifges(6,3) * t815;
t751 = Ifges(6,1) * t776 + Ifges(6,4) * t775 + Ifges(6,5) * t815;
t742 = -mrSges(7,2) * t814 + mrSges(7,3) * t755;
t743 = mrSges(7,1) * t814 - mrSges(7,3) * t756;
t843 = m(7) * t722 - t729 * mrSges(7,1) + t730 * mrSges(7,2) - t755 * t742 + t756 * t743;
t740 = -mrSges(7,1) * t755 + mrSges(7,2) * t756;
t710 = m(7) * t716 + mrSges(7,1) * t796 - mrSges(7,3) * t730 - t740 * t756 + t742 * t814;
t711 = m(7) * t717 - mrSges(7,2) * t796 + mrSges(7,3) * t729 + t740 * t755 - t743 * t814;
t846 = -t710 * t827 + t831 * t711;
t693 = -mrSges(6,1) * t741 + mrSges(6,3) * t721 + Ifges(6,4) * t746 + Ifges(6,2) * t745 + Ifges(6,6) * t802 - pkin(5) * t843 + pkin(9) * t846 + t831 * t705 + t827 * t706 - t776 * t749 + t815 * t751;
t704 = t831 * t710 + t827 * t711;
t750 = Ifges(6,4) * t776 + Ifges(6,2) * t775 + Ifges(6,6) * t815;
t694 = mrSges(6,2) * t741 - mrSges(6,3) * t720 + Ifges(6,1) * t746 + Ifges(6,4) * t745 + Ifges(6,5) * t802 - pkin(9) * t704 - t705 * t827 + t706 * t831 + t749 * t775 - t750 * t815;
t769 = Ifges(5,5) * t800 + Ifges(5,6) * t799 - Ifges(5,3) * t853;
t771 = Ifges(5,1) * t800 + Ifges(5,4) * t799 - Ifges(5,5) * t853;
t762 = -mrSges(6,2) * t815 + mrSges(6,3) * t775;
t763 = mrSges(6,1) * t815 - mrSges(6,3) * t776;
t839 = m(6) * t741 - t745 * mrSges(6,1) + t746 * mrSges(6,2) - t775 * t762 + t776 * t763 + t843;
t758 = -mrSges(6,1) * t775 + mrSges(6,2) * t776;
t702 = m(6) * t720 + mrSges(6,1) * t802 - mrSges(6,3) * t746 - t758 * t776 + t762 * t815 + t704;
t703 = m(6) * t721 - mrSges(6,2) * t802 + mrSges(6,3) * t745 + t758 * t775 - t763 * t815 + t846;
t847 = -t702 * t828 + t832 * t703;
t677 = -mrSges(5,1) * t757 + mrSges(5,3) * t735 + Ifges(5,4) * t784 + Ifges(5,2) * t783 - Ifges(5,6) * t808 - pkin(4) * t839 + pkin(8) * t847 + t832 * t693 + t828 * t694 - t800 * t769 - t771 * t853;
t698 = t832 * t702 + t828 * t703;
t770 = Ifges(5,4) * t800 + Ifges(5,2) * t799 - Ifges(5,6) * t853;
t678 = mrSges(5,2) * t757 - mrSges(5,3) * t734 + Ifges(5,1) * t784 + Ifges(5,4) * t783 - Ifges(5,5) * t808 - pkin(8) * t698 - t693 * t828 + t694 * t832 + t769 * t799 + t770 * t853;
t780 = -mrSges(5,1) * t799 + mrSges(5,2) * t800;
t844 = mrSges(5,2) * t853 + mrSges(5,3) * t799;
t696 = m(5) * t734 - t808 * mrSges(5,1) - t784 * mrSges(5,3) - t800 * t780 - t844 * t853 + t698;
t782 = -mrSges(5,1) * t853 - mrSges(5,3) * t800;
t697 = m(5) * t735 + mrSges(5,2) * t808 + mrSges(5,3) * t783 + t780 * t799 + t782 * t853 + t847;
t692 = -t696 * t823 + t825 * t697;
t714 = m(5) * t757 - t783 * mrSges(5,1) + t784 * mrSges(5,2) + t800 * t782 - t799 * t844 + t839;
t794 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t829 + Ifges(4,2) * t833) * qJD(1);
t795 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t829 + Ifges(4,4) * t833) * qJD(1);
t856 = mrSges(4,1) * t760 - mrSges(4,2) * t761 + Ifges(4,5) * t807 + Ifges(4,6) * t808 + Ifges(4,3) * qJDD(3) - pkin(3) * t714 + qJ(4) * t692 + t825 * t677 + t823 * t678 + (t794 * t829 - t795 * t833) * qJD(1);
t805 = (-mrSges(4,1) * t833 + mrSges(4,2) * t829) * qJD(1);
t810 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t854;
t690 = m(4) * t761 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t808 - qJD(3) * t810 + t805 * t853 + t692;
t811 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t853;
t713 = m(4) * t760 + qJDD(3) * mrSges(4,1) - t807 * mrSges(4,3) + qJD(3) * t811 - t805 * t854 - t714;
t848 = t833 * t690 - t713 * t829;
t681 = m(3) * t779 - mrSges(3,1) * t836 - qJDD(1) * mrSges(3,2) + t848;
t691 = t696 * t825 + t697 * t823;
t840 = -m(4) * t767 + t808 * mrSges(4,1) - mrSges(4,2) * t807 - t810 * t854 + t811 * t853 - t691;
t686 = m(3) * t778 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t836 + t840;
t674 = t824 * t681 + t826 * t686;
t671 = m(2) * t812 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t836 + t674;
t849 = t826 * t681 - t686 * t824;
t672 = m(2) * t813 - mrSges(2,1) * t836 - qJDD(1) * mrSges(2,2) + t849;
t855 = t834 * t671 + t830 * t672;
t684 = t829 * t690 + t833 * t713;
t682 = m(3) * t822 + t684;
t850 = -t671 * t830 + t834 * t672;
t842 = -mrSges(7,1) * t716 + mrSges(7,2) * t717 - Ifges(7,5) * t730 - Ifges(7,6) * t729 - Ifges(7,3) * t796 - t756 * t737 + t755 * t738;
t793 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t829 + Ifges(4,6) * t833) * qJD(1);
t667 = mrSges(4,2) * t767 - mrSges(4,3) * t760 + Ifges(4,1) * t807 + Ifges(4,4) * t808 + Ifges(4,5) * qJDD(3) - qJ(4) * t691 - qJD(3) * t794 - t677 * t823 + t678 * t825 + t793 * t853;
t837 = mrSges(6,1) * t720 - mrSges(6,2) * t721 + Ifges(6,5) * t746 + Ifges(6,6) * t745 + Ifges(6,3) * t802 + pkin(5) * t704 + t776 * t750 - t775 * t751 - t842;
t676 = Ifges(4,4) * t807 + t799 * t771 - t800 * t770 - Ifges(5,6) * t783 - Ifges(5,5) * t784 + qJD(3) * t795 + mrSges(4,3) * t761 - mrSges(4,1) * t767 - mrSges(5,1) * t734 + mrSges(5,2) * t735 + Ifges(4,6) * qJDD(3) - pkin(4) * t698 - pkin(3) * t691 + (Ifges(4,2) + Ifges(5,3)) * t808 - t793 * t854 - t837;
t841 = mrSges(2,1) * t812 + mrSges(3,1) * t778 - mrSges(2,2) * t813 - mrSges(3,2) * t779 + pkin(1) * t674 + pkin(2) * t840 + pkin(7) * t848 + t829 * t667 + t833 * t676 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t665 = -mrSges(3,1) * t822 + mrSges(3,3) * t779 + t836 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t684 - t856;
t664 = mrSges(3,2) * t822 - mrSges(3,3) * t778 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t836 - pkin(7) * t684 + t667 * t833 - t676 * t829;
t663 = -mrSges(2,2) * g(3) - mrSges(2,3) * t812 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t836 - qJ(2) * t674 + t664 * t826 - t665 * t824;
t662 = mrSges(2,1) * g(3) + mrSges(2,3) * t813 + t836 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t682 + qJ(2) * t849 + t824 * t664 + t826 * t665;
t1 = [-m(1) * g(1) + t850; -m(1) * g(2) + t855; (-m(1) - m(2)) * g(3) + t682; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t855 - t830 * t662 + t834 * t663; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t850 + t834 * t662 + t830 * t663; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t841; t841; t682; t856; t714; t837; -t842;];
tauJB  = t1;
