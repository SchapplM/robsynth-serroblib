% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 19:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRPRR9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:32:59
% EndTime: 2019-05-05 19:34:40
% DurationCPUTime: 104.27s
% Computational Cost: add. (1534512->402), mult. (5115038->550), div. (0->0), fcn. (4400475->16), ass. (0->182)
t795 = sin(pkin(7));
t798 = cos(pkin(12));
t800 = cos(pkin(6));
t796 = sin(pkin(6));
t799 = cos(pkin(7));
t841 = t796 * t799;
t780 = (t795 * t800 + t798 * t841) * qJD(1) * pkin(9);
t804 = sin(qJ(1));
t808 = cos(qJ(1));
t791 = -g(1) * t808 - g(2) * t804;
t809 = qJD(1) ^ 2;
t847 = qJ(2) * t796;
t784 = -pkin(1) * t809 + qJDD(1) * t847 + t791;
t794 = sin(pkin(12));
t850 = pkin(9) * t794;
t823 = -pkin(2) * t798 - t795 * t850;
t836 = qJD(1) * t796;
t848 = pkin(9) * qJDD(1);
t819 = qJD(1) * t823 * t836 + t799 * t848;
t790 = t804 * g(1) - g(2) * t808;
t783 = qJDD(1) * pkin(1) + t809 * t847 + t790;
t833 = qJD(2) * t836;
t840 = t798 * t800;
t842 = t796 * t798;
t824 = -g(3) * t842 + t783 * t840 - 0.2e1 * t794 * t833;
t738 = (pkin(2) * qJDD(1) + qJD(1) * t780) * t800 + (-t796 * t819 - t784) * t794 + t824;
t785 = (pkin(2) * t800 - t841 * t850) * qJD(1);
t845 = t794 * t800;
t834 = t783 * t845 + (t784 + 0.2e1 * t833) * t798;
t739 = (-qJD(1) * t785 + t795 * t848) * t800 + (-g(3) * t794 + t798 * t819) * t796 + t834;
t832 = -t800 * g(3) + qJDD(2);
t746 = (-t783 + t823 * qJDD(1) + (-t780 * t798 + t785 * t794) * qJD(1)) * t796 + t832;
t803 = sin(qJ(3));
t807 = cos(qJ(3));
t838 = t799 * t807;
t843 = t795 * t807;
t704 = t738 * t838 - t803 * t739 + t746 * t843;
t811 = t800 * t843 + (-t794 * t803 + t798 * t838) * t796;
t768 = t811 * qJD(1);
t839 = t799 * t803;
t844 = t795 * t803;
t812 = t800 * t844 + (t794 * t807 + t798 * t839) * t796;
t761 = t768 * qJD(3) + qJDD(1) * t812;
t769 = t812 * qJD(1);
t820 = -t795 * t842 + t799 * t800;
t778 = qJDD(1) * t820 + qJDD(3);
t781 = qJD(1) * t820 + qJD(3);
t694 = (t768 * t781 - t761) * qJ(4) + (t768 * t769 + t778) * pkin(3) + t704;
t705 = t738 * t839 + t807 * t739 + t746 * t844;
t760 = -t769 * qJD(3) + qJDD(1) * t811;
t765 = pkin(3) * t781 - qJ(4) * t769;
t767 = t768 ^ 2;
t697 = -pkin(3) * t767 + qJ(4) * t760 - t765 * t781 + t705;
t793 = sin(pkin(13));
t797 = cos(pkin(13));
t758 = t768 * t793 + t769 * t797;
t686 = -0.2e1 * qJD(4) * t758 + t797 * t694 - t793 * t697;
t849 = Ifges(3,3) * t800;
t846 = t794 * t796;
t757 = t768 * t797 - t769 * t793;
t687 = 0.2e1 * qJD(4) * t757 + t793 * t694 + t797 * t697;
t729 = -mrSges(5,1) * t757 + mrSges(5,2) * t758;
t733 = t760 * t797 - t761 * t793;
t748 = mrSges(5,1) * t781 - mrSges(5,3) * t758;
t730 = -pkin(4) * t757 - pkin(10) * t758;
t777 = t781 ^ 2;
t685 = -pkin(4) * t777 + pkin(10) * t778 + t730 * t757 + t687;
t717 = -t795 * t738 + t799 * t746;
t703 = -t760 * pkin(3) - t767 * qJ(4) + t769 * t765 + qJDD(4) + t717;
t734 = t760 * t793 + t761 * t797;
t689 = (-t757 * t781 - t734) * pkin(10) + (t758 * t781 - t733) * pkin(4) + t703;
t802 = sin(qJ(5));
t806 = cos(qJ(5));
t681 = t806 * t685 + t802 * t689;
t745 = t758 * t806 + t781 * t802;
t715 = -qJD(5) * t745 - t734 * t802 + t778 * t806;
t744 = -t758 * t802 + t781 * t806;
t719 = -mrSges(6,1) * t744 + mrSges(6,2) * t745;
t756 = qJD(5) - t757;
t724 = mrSges(6,1) * t756 - mrSges(6,3) * t745;
t732 = qJDD(5) - t733;
t720 = -pkin(5) * t744 - pkin(11) * t745;
t755 = t756 ^ 2;
t679 = -pkin(5) * t755 + pkin(11) * t732 + t720 * t744 + t681;
t684 = -t778 * pkin(4) - t777 * pkin(10) + t758 * t730 - t686;
t716 = qJD(5) * t744 + t734 * t806 + t778 * t802;
t682 = (-t744 * t756 - t716) * pkin(11) + (t745 * t756 - t715) * pkin(5) + t684;
t801 = sin(qJ(6));
t805 = cos(qJ(6));
t676 = -t679 * t801 + t682 * t805;
t721 = -t745 * t801 + t756 * t805;
t692 = qJD(6) * t721 + t716 * t805 + t732 * t801;
t722 = t745 * t805 + t756 * t801;
t706 = -mrSges(7,1) * t721 + mrSges(7,2) * t722;
t743 = qJD(6) - t744;
t707 = -mrSges(7,2) * t743 + mrSges(7,3) * t721;
t714 = qJDD(6) - t715;
t674 = m(7) * t676 + mrSges(7,1) * t714 - mrSges(7,3) * t692 - t706 * t722 + t707 * t743;
t677 = t679 * t805 + t682 * t801;
t691 = -qJD(6) * t722 - t716 * t801 + t732 * t805;
t708 = mrSges(7,1) * t743 - mrSges(7,3) * t722;
t675 = m(7) * t677 - mrSges(7,2) * t714 + mrSges(7,3) * t691 + t706 * t721 - t708 * t743;
t828 = -t674 * t801 + t805 * t675;
t667 = m(6) * t681 - mrSges(6,2) * t732 + mrSges(6,3) * t715 + t719 * t744 - t724 * t756 + t828;
t680 = -t685 * t802 + t689 * t806;
t723 = -mrSges(6,2) * t756 + mrSges(6,3) * t744;
t678 = -pkin(5) * t732 - pkin(11) * t755 + t720 * t745 - t680;
t816 = -m(7) * t678 + t691 * mrSges(7,1) - mrSges(7,2) * t692 + t721 * t707 - t708 * t722;
t672 = m(6) * t680 + mrSges(6,1) * t732 - mrSges(6,3) * t716 - t719 * t745 + t723 * t756 + t816;
t829 = t806 * t667 - t672 * t802;
t659 = m(5) * t687 - mrSges(5,2) * t778 + mrSges(5,3) * t733 + t729 * t757 - t748 * t781 + t829;
t747 = -mrSges(5,2) * t781 + mrSges(5,3) * t757;
t668 = t674 * t805 + t675 * t801;
t810 = -m(6) * t684 + t715 * mrSges(6,1) - mrSges(6,2) * t716 + t744 * t723 - t724 * t745 - t668;
t664 = m(5) * t686 + mrSges(5,1) * t778 - mrSges(5,3) * t734 - t729 * t758 + t747 * t781 + t810;
t654 = t793 * t659 + t797 * t664;
t759 = -mrSges(4,1) * t768 + mrSges(4,2) * t769;
t764 = -mrSges(4,2) * t781 + mrSges(4,3) * t768;
t652 = m(4) * t704 + mrSges(4,1) * t778 - mrSges(4,3) * t761 - t759 * t769 + t764 * t781 + t654;
t766 = mrSges(4,1) * t781 - mrSges(4,3) * t769;
t830 = t797 * t659 - t664 * t793;
t653 = m(4) * t705 - mrSges(4,2) * t778 + mrSges(4,3) * t760 + t759 * t768 - t766 * t781 + t830;
t662 = t802 * t667 + t806 * t672;
t813 = m(5) * t703 - mrSges(5,1) * t733 + t734 * mrSges(5,2) - t747 * t757 + t758 * t748 + t662;
t661 = m(4) * t717 - mrSges(4,1) * t760 + mrSges(4,2) * t761 - t764 * t768 + t766 * t769 + t813;
t639 = t652 * t838 + t653 * t839 - t795 * t661;
t762 = -t794 * t784 + t824;
t826 = -mrSges(3,1) * t798 + mrSges(3,2) * t794;
t782 = t826 * t836;
t821 = -mrSges(3,2) * t800 + mrSges(3,3) * t842;
t787 = t821 * qJD(1);
t822 = mrSges(3,1) * t800 - mrSges(3,3) * t846;
t635 = m(3) * t762 + t822 * qJDD(1) + (-t782 * t846 + t787 * t800) * qJD(1) + t639;
t638 = t652 * t843 + t653 * t844 + t799 * t661;
t770 = -t796 * t783 + t832;
t786 = t822 * qJD(1);
t637 = m(3) * t770 + (t826 * qJDD(1) + (t786 * t794 - t787 * t798) * qJD(1)) * t796 + t638;
t645 = -t803 * t652 + t807 * t653;
t763 = -g(3) * t846 + t834;
t644 = m(3) * t763 + t821 * qJDD(1) + (t782 * t842 - t786 * t800) * qJD(1) + t645;
t625 = t635 * t840 - t637 * t796 + t644 * t845;
t623 = m(2) * t790 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t809 + t625;
t631 = -t635 * t794 + t798 * t644;
t630 = m(2) * t791 - mrSges(2,1) * t809 - qJDD(1) * mrSges(2,2) + t631;
t837 = t808 * t623 + t804 * t630;
t624 = t635 * t842 + t800 * t637 + t644 * t846;
t831 = -t623 * t804 + t808 * t630;
t825 = Ifges(3,5) * t794 + Ifges(3,6) * t798;
t698 = Ifges(7,5) * t722 + Ifges(7,6) * t721 + Ifges(7,3) * t743;
t700 = Ifges(7,1) * t722 + Ifges(7,4) * t721 + Ifges(7,5) * t743;
t669 = -mrSges(7,1) * t678 + mrSges(7,3) * t677 + Ifges(7,4) * t692 + Ifges(7,2) * t691 + Ifges(7,6) * t714 - t698 * t722 + t700 * t743;
t699 = Ifges(7,4) * t722 + Ifges(7,2) * t721 + Ifges(7,6) * t743;
t670 = mrSges(7,2) * t678 - mrSges(7,3) * t676 + Ifges(7,1) * t692 + Ifges(7,4) * t691 + Ifges(7,5) * t714 + t698 * t721 - t699 * t743;
t709 = Ifges(6,5) * t745 + Ifges(6,6) * t744 + Ifges(6,3) * t756;
t710 = Ifges(6,4) * t745 + Ifges(6,2) * t744 + Ifges(6,6) * t756;
t655 = mrSges(6,2) * t684 - mrSges(6,3) * t680 + Ifges(6,1) * t716 + Ifges(6,4) * t715 + Ifges(6,5) * t732 - pkin(11) * t668 - t669 * t801 + t670 * t805 + t709 * t744 - t710 * t756;
t711 = Ifges(6,1) * t745 + Ifges(6,4) * t744 + Ifges(6,5) * t756;
t656 = -mrSges(6,1) * t684 - mrSges(7,1) * t676 + mrSges(7,2) * t677 + mrSges(6,3) * t681 + Ifges(6,4) * t716 - Ifges(7,5) * t692 + Ifges(6,2) * t715 + Ifges(6,6) * t732 - Ifges(7,6) * t691 - Ifges(7,3) * t714 - pkin(5) * t668 - t699 * t722 + t700 * t721 - t709 * t745 + t711 * t756;
t725 = Ifges(5,5) * t758 + Ifges(5,6) * t757 + Ifges(5,3) * t781;
t726 = Ifges(5,4) * t758 + Ifges(5,2) * t757 + Ifges(5,6) * t781;
t640 = mrSges(5,2) * t703 - mrSges(5,3) * t686 + Ifges(5,1) * t734 + Ifges(5,4) * t733 + Ifges(5,5) * t778 - pkin(10) * t662 + t655 * t806 - t656 * t802 + t725 * t757 - t726 * t781;
t727 = Ifges(5,1) * t758 + Ifges(5,4) * t757 + Ifges(5,5) * t781;
t646 = Ifges(5,4) * t734 + Ifges(5,2) * t733 + Ifges(5,6) * t778 - t758 * t725 + t781 * t727 - mrSges(5,1) * t703 + mrSges(5,3) * t687 - Ifges(6,5) * t716 - Ifges(6,6) * t715 - Ifges(6,3) * t732 - t745 * t710 + t744 * t711 - mrSges(6,1) * t680 + mrSges(6,2) * t681 - t801 * t670 - t805 * t669 - pkin(5) * t816 - pkin(11) * t828 - pkin(4) * t662;
t749 = Ifges(4,5) * t769 + Ifges(4,6) * t768 + Ifges(4,3) * t781;
t751 = Ifges(4,1) * t769 + Ifges(4,4) * t768 + Ifges(4,5) * t781;
t626 = -mrSges(4,1) * t717 + mrSges(4,3) * t705 + Ifges(4,4) * t761 + Ifges(4,2) * t760 + Ifges(4,6) * t778 - pkin(3) * t813 + qJ(4) * t830 + t793 * t640 + t797 * t646 - t769 * t749 + t781 * t751;
t750 = Ifges(4,4) * t769 + Ifges(4,2) * t768 + Ifges(4,6) * t781;
t627 = mrSges(4,2) * t717 - mrSges(4,3) * t704 + Ifges(4,1) * t761 + Ifges(4,4) * t760 + Ifges(4,5) * t778 - qJ(4) * t654 + t640 * t797 - t646 * t793 + t749 * t768 - t750 * t781;
t818 = pkin(9) * t645 + t626 * t807 + t627 * t803;
t632 = Ifges(4,5) * t761 + Ifges(4,6) * t760 + t769 * t750 - t768 * t751 + mrSges(4,1) * t704 - mrSges(4,2) * t705 + Ifges(5,5) * t734 + Ifges(5,6) * t733 + t758 * t726 - t757 * t727 + mrSges(5,1) * t686 - mrSges(5,2) * t687 + t802 * t655 + t806 * t656 + pkin(4) * t810 + pkin(10) * t829 + pkin(3) * t654 + (Ifges(4,3) + Ifges(5,3)) * t778;
t773 = (t796 * t825 + t849) * qJD(1);
t815 = Ifges(3,5) * t800 + (Ifges(3,1) * t794 + Ifges(3,4) * t798) * t796;
t775 = t815 * qJD(1);
t814 = Ifges(3,6) * t800 + (Ifges(3,4) * t794 + Ifges(3,2) * t798) * t796;
t620 = -mrSges(3,1) * t770 + mrSges(3,3) * t763 - pkin(2) * t638 - t795 * t632 + (-t773 * t846 + t775 * t800) * qJD(1) + t818 * t799 + t814 * qJDD(1);
t774 = t814 * qJD(1);
t621 = mrSges(3,2) * t770 - mrSges(3,3) * t762 - t803 * t626 + t807 * t627 + (t773 * t842 - t774 * t800) * qJD(1) + (-t638 * t795 - t639 * t799) * pkin(9) + t815 * qJDD(1);
t817 = qJ(2) * t631 + t620 * t798 + t621 * t794;
t619 = qJDD(1) * t849 + mrSges(3,1) * t762 - mrSges(3,2) * t763 + pkin(2) * t639 + t799 * t632 + t818 * t795 + (t825 * qJDD(1) + (t774 * t794 - t775 * t798) * qJD(1)) * t796;
t618 = -mrSges(2,2) * g(3) - mrSges(2,3) * t790 + Ifges(2,5) * qJDD(1) - t809 * Ifges(2,6) - t794 * t620 + t798 * t621 + (-t624 * t796 - t625 * t800) * qJ(2);
t617 = mrSges(2,1) * g(3) + mrSges(2,3) * t791 + t809 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t624 - t796 * t619 + t800 * t817;
t1 = [-m(1) * g(1) + t831; -m(1) * g(2) + t837; (-m(1) - m(2)) * g(3) + t624; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t837 - t804 * t617 + t808 * t618; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t831 + t808 * t617 + t804 * t618; -mrSges(1,1) * g(2) + mrSges(2,1) * t790 + mrSges(1,2) * g(1) - mrSges(2,2) * t791 + Ifges(2,3) * qJDD(1) + pkin(1) * t625 + t800 * t619 + t796 * t817;];
tauB  = t1;
