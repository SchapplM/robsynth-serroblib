% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 14:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:23:22
% EndTime: 2019-05-07 14:23:51
% DurationCPUTime: 18.18s
% Computational Cost: add. (287653->377), mult. (616482->472), div. (0->0), fcn. (481741->12), ass. (0->158)
t849 = Ifges(4,1) + Ifges(5,1);
t840 = Ifges(4,4) - Ifges(5,5);
t839 = Ifges(4,5) + Ifges(5,4);
t848 = Ifges(4,2) + Ifges(5,3);
t847 = Ifges(5,6) - Ifges(4,6);
t846 = -Ifges(4,3) - Ifges(5,2);
t845 = 2 * qJD(4);
t844 = cos(qJ(3));
t799 = sin(pkin(6));
t808 = cos(qJ(2));
t825 = qJD(1) * t808;
t822 = t799 * t825;
t788 = qJD(3) - t822;
t843 = pkin(3) * t788;
t800 = cos(pkin(6));
t842 = t800 * g(3);
t841 = -mrSges(4,3) - mrSges(5,2);
t796 = qJD(1) * t800 + qJD(2);
t803 = sin(qJ(3));
t804 = sin(qJ(2));
t826 = qJD(1) * t799;
t823 = t804 * t826;
t765 = -t796 * t844 + t803 * t823;
t837 = t765 * t788;
t836 = t799 * t804;
t835 = t799 * t808;
t834 = t800 * t804;
t833 = t800 * t808;
t805 = sin(qJ(1));
t809 = cos(qJ(1));
t789 = t805 * g(1) - g(2) * t809;
t810 = qJD(1) ^ 2;
t776 = pkin(8) * t799 * t810 + qJDD(1) * pkin(1) + t789;
t790 = -g(1) * t809 - g(2) * t805;
t824 = qJDD(1) * t799;
t777 = -pkin(1) * t810 + pkin(8) * t824 + t790;
t827 = t776 * t834 + t808 * t777;
t741 = -g(3) * t836 + t827;
t774 = mrSges(3,1) * t796 - mrSges(3,3) * t823;
t778 = (-mrSges(3,1) * t808 + mrSges(3,2) * t804) * t826;
t781 = -qJD(2) * t823 + t808 * t824;
t795 = qJDD(1) * t800 + qJDD(2);
t779 = (-pkin(2) * t808 - pkin(9) * t804) * t826;
t794 = t796 ^ 2;
t717 = -t794 * pkin(2) + t795 * pkin(9) + (-g(3) * t804 + t779 * t825) * t799 + t827;
t780 = (qJD(2) * t825 + qJDD(1) * t804) * t799;
t718 = -t781 * pkin(2) - t780 * pkin(9) - t842 + (-t776 + (pkin(2) * t804 - pkin(9) * t808) * t796 * qJD(1)) * t799;
t693 = t844 * t717 + t803 * t718;
t766 = t803 * t796 + t823 * t844;
t738 = t766 * qJD(3) + t803 * t780 - t795 * t844;
t751 = mrSges(4,1) * t788 - mrSges(4,3) * t766;
t773 = qJDD(3) - t781;
t745 = pkin(3) * t765 - qJ(4) * t766;
t787 = t788 ^ 2;
t686 = -pkin(3) * t787 + t773 * qJ(4) - t765 * t745 + t788 * t845 + t693;
t752 = -mrSges(5,1) * t788 + mrSges(5,2) * t766;
t692 = -t803 * t717 + t718 * t844;
t687 = -t773 * pkin(3) - t787 * qJ(4) + t766 * t745 + qJDD(4) - t692;
t739 = -t765 * qJD(3) + t780 * t844 + t803 * t795;
t680 = (-t739 - t837) * pkin(10) + (t765 * t766 - t773) * pkin(4) + t687;
t754 = -pkin(4) * t788 - pkin(10) * t766;
t764 = t765 ^ 2;
t683 = -pkin(4) * t764 + pkin(10) * t738 + t754 * t788 + t686;
t802 = sin(qJ(5));
t807 = cos(qJ(5));
t678 = t802 * t680 + t807 * t683;
t744 = t765 * t802 + t766 * t807;
t701 = -qJD(5) * t744 + t738 * t807 - t739 * t802;
t743 = t765 * t807 - t766 * t802;
t711 = -mrSges(6,1) * t743 + mrSges(6,2) * t744;
t784 = qJD(5) - t788;
t724 = mrSges(6,1) * t784 - mrSges(6,3) * t744;
t772 = qJDD(5) - t773;
t712 = -pkin(5) * t743 - pkin(11) * t744;
t783 = t784 ^ 2;
t675 = -pkin(5) * t783 + pkin(11) * t772 + t712 * t743 + t678;
t740 = -g(3) * t835 + t776 * t833 - t804 * t777;
t716 = -t795 * pkin(2) - t794 * pkin(9) + t779 * t823 - t740;
t814 = t738 * pkin(3) + t716 + (-t739 + t837) * qJ(4);
t684 = -t738 * pkin(4) - t764 * pkin(10) - t814 + (t754 - t843 + t845) * t766;
t702 = qJD(5) * t743 + t738 * t802 + t739 * t807;
t676 = t684 + (-t743 * t784 - t702) * pkin(11) + (t744 * t784 - t701) * pkin(5);
t801 = sin(qJ(6));
t806 = cos(qJ(6));
t672 = -t675 * t801 + t676 * t806;
t721 = -t744 * t801 + t784 * t806;
t691 = qJD(6) * t721 + t702 * t806 + t772 * t801;
t700 = qJDD(6) - t701;
t722 = t744 * t806 + t784 * t801;
t703 = -mrSges(7,1) * t721 + mrSges(7,2) * t722;
t742 = qJD(6) - t743;
t704 = -mrSges(7,2) * t742 + mrSges(7,3) * t721;
t670 = m(7) * t672 + mrSges(7,1) * t700 - mrSges(7,3) * t691 - t703 * t722 + t704 * t742;
t673 = t675 * t806 + t676 * t801;
t690 = -qJD(6) * t722 - t702 * t801 + t772 * t806;
t705 = mrSges(7,1) * t742 - mrSges(7,3) * t722;
t671 = m(7) * t673 - mrSges(7,2) * t700 + mrSges(7,3) * t690 + t703 * t721 - t705 * t742;
t818 = -t670 * t801 + t806 * t671;
t662 = m(6) * t678 - mrSges(6,2) * t772 + mrSges(6,3) * t701 + t711 * t743 - t724 * t784 + t818;
t677 = t680 * t807 - t683 * t802;
t723 = -mrSges(6,2) * t784 + mrSges(6,3) * t743;
t674 = -pkin(5) * t772 - pkin(11) * t783 + t712 * t744 - t677;
t813 = -m(7) * t674 + t690 * mrSges(7,1) - mrSges(7,2) * t691 + t721 * t704 - t705 * t722;
t666 = m(6) * t677 + mrSges(6,1) * t772 - mrSges(6,3) * t702 - t711 * t744 + t723 * t784 + t813;
t819 = t807 * t662 - t802 * t666;
t817 = m(5) * t686 + t773 * mrSges(5,3) + t788 * t752 + t819;
t746 = mrSges(5,1) * t765 - mrSges(5,3) * t766;
t828 = -mrSges(4,1) * t765 - mrSges(4,2) * t766 - t746;
t653 = m(4) * t693 - t773 * mrSges(4,2) + t738 * t841 - t788 * t751 + t765 * t828 + t817;
t750 = -mrSges(4,2) * t788 - mrSges(4,3) * t765;
t656 = t802 * t662 + t807 * t666;
t753 = -mrSges(5,2) * t765 + mrSges(5,3) * t788;
t812 = -m(5) * t687 + t773 * mrSges(5,1) + t788 * t753 - t656;
t655 = m(4) * t692 + t773 * mrSges(4,1) + t739 * t841 + t788 * t750 + t766 * t828 + t812;
t820 = t844 * t653 - t655 * t803;
t645 = m(3) * t741 - mrSges(3,2) * t795 + mrSges(3,3) * t781 - t774 * t796 + t778 * t822 + t820;
t648 = t803 * t653 + t844 * t655;
t759 = -t799 * t776 - t842;
t775 = -mrSges(3,2) * t796 + mrSges(3,3) * t822;
t647 = m(3) * t759 - t781 * mrSges(3,1) + t780 * mrSges(3,2) + (t774 * t804 - t775 * t808) * t826 + t648;
t688 = (-(2 * qJD(4)) + t843) * t766 + t814;
t663 = t806 * t670 + t801 * t671;
t816 = -m(6) * t684 + t701 * mrSges(6,1) - t702 * mrSges(6,2) + t743 * t723 - t744 * t724 - t663;
t660 = m(5) * t688 + t738 * mrSges(5,1) - t739 * mrSges(5,3) - t766 * t752 + t765 * t753 + t816;
t811 = -m(4) * t716 - t738 * mrSges(4,1) - t739 * mrSges(4,2) - t765 * t750 - t766 * t751 - t660;
t659 = m(3) * t740 + t795 * mrSges(3,1) - t780 * mrSges(3,3) + t796 * t775 - t778 * t823 + t811;
t636 = t645 * t834 - t647 * t799 + t659 * t833;
t634 = m(2) * t789 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t810 + t636;
t641 = t808 * t645 - t659 * t804;
t640 = m(2) * t790 - mrSges(2,1) * t810 - qJDD(1) * mrSges(2,2) + t641;
t832 = t809 * t634 + t805 * t640;
t831 = t765 * t848 - t766 * t840 + t788 * t847;
t830 = -t765 * t847 - t766 * t839 + t788 * t846;
t829 = -t840 * t765 + t766 * t849 + t839 * t788;
t635 = t645 * t836 + t800 * t647 + t659 * t835;
t821 = -t634 * t805 + t809 * t640;
t694 = Ifges(7,5) * t722 + Ifges(7,6) * t721 + Ifges(7,3) * t742;
t696 = Ifges(7,1) * t722 + Ifges(7,4) * t721 + Ifges(7,5) * t742;
t664 = -mrSges(7,1) * t674 + mrSges(7,3) * t673 + Ifges(7,4) * t691 + Ifges(7,2) * t690 + Ifges(7,6) * t700 - t694 * t722 + t696 * t742;
t695 = Ifges(7,4) * t722 + Ifges(7,2) * t721 + Ifges(7,6) * t742;
t665 = mrSges(7,2) * t674 - mrSges(7,3) * t672 + Ifges(7,1) * t691 + Ifges(7,4) * t690 + Ifges(7,5) * t700 + t694 * t721 - t695 * t742;
t706 = Ifges(6,5) * t744 + Ifges(6,6) * t743 + Ifges(6,3) * t784;
t707 = Ifges(6,4) * t744 + Ifges(6,2) * t743 + Ifges(6,6) * t784;
t649 = mrSges(6,2) * t684 - mrSges(6,3) * t677 + Ifges(6,1) * t702 + Ifges(6,4) * t701 + Ifges(6,5) * t772 - pkin(11) * t663 - t664 * t801 + t665 * t806 + t706 * t743 - t707 * t784;
t708 = Ifges(6,1) * t744 + Ifges(6,4) * t743 + Ifges(6,5) * t784;
t650 = -mrSges(6,1) * t684 - mrSges(7,1) * t672 + mrSges(7,2) * t673 + mrSges(6,3) * t678 + Ifges(6,4) * t702 - Ifges(7,5) * t691 + Ifges(6,2) * t701 + Ifges(6,6) * t772 - Ifges(7,6) * t690 - Ifges(7,3) * t700 - pkin(5) * t663 - t695 * t722 + t696 * t721 - t706 * t744 + t708 * t784;
t632 = -mrSges(4,1) * t716 - mrSges(5,1) * t688 + mrSges(5,2) * t686 + mrSges(4,3) * t693 - pkin(3) * t660 - pkin(4) * t816 - pkin(10) * t819 - t802 * t649 - t807 * t650 - t738 * t848 + t840 * t739 + t830 * t766 - t773 * t847 + t829 * t788;
t637 = mrSges(4,2) * t716 + mrSges(5,2) * t687 - mrSges(4,3) * t692 - mrSges(5,3) * t688 - pkin(10) * t656 - qJ(4) * t660 + t807 * t649 - t802 * t650 - t840 * t738 + t739 * t849 + t830 * t765 + t839 * t773 + t831 * t788;
t756 = Ifges(3,3) * t796 + (Ifges(3,5) * t804 + Ifges(3,6) * t808) * t826;
t757 = Ifges(3,6) * t796 + (Ifges(3,4) * t804 + Ifges(3,2) * t808) * t826;
t630 = mrSges(3,2) * t759 - mrSges(3,3) * t740 + Ifges(3,1) * t780 + Ifges(3,4) * t781 + Ifges(3,5) * t795 - pkin(9) * t648 - t803 * t632 + t637 * t844 + t756 * t822 - t796 * t757;
t758 = Ifges(3,5) * t796 + (Ifges(3,1) * t804 + Ifges(3,4) * t808) * t826;
t631 = (mrSges(5,2) * pkin(3) - t839) * t739 + (mrSges(5,2) * qJ(4) - t847) * t738 - t756 * t823 + pkin(11) * t818 + (qJ(4) * t746 - t829) * t765 + (pkin(3) * t746 + t831) * t766 - qJ(4) * t817 + t846 * t773 + t806 * t664 + t801 * t665 + Ifges(3,6) * t795 + t796 * t758 + pkin(4) * t656 + Ifges(6,3) * t772 + Ifges(3,4) * t780 + Ifges(3,2) * t781 - mrSges(3,1) * t759 + t744 * t707 + mrSges(3,3) * t741 - t743 * t708 + Ifges(6,6) * t701 + Ifges(6,5) * t702 - mrSges(4,1) * t692 + mrSges(4,2) * t693 - mrSges(5,3) * t686 + mrSges(5,1) * t687 + mrSges(6,1) * t677 - mrSges(6,2) * t678 - pkin(3) * t812 + pkin(5) * t813 - pkin(2) * t648;
t815 = pkin(8) * t641 + t630 * t804 + t631 * t808;
t629 = Ifges(3,5) * t780 + Ifges(3,6) * t781 + Ifges(3,3) * t795 + mrSges(3,1) * t740 - mrSges(3,2) * t741 + t803 * t637 + t844 * t632 + pkin(2) * t811 + pkin(9) * t820 + (t757 * t804 - t758 * t808) * t826;
t628 = -mrSges(2,2) * g(3) - mrSges(2,3) * t789 + Ifges(2,5) * qJDD(1) - t810 * Ifges(2,6) + t808 * t630 - t804 * t631 + (-t635 * t799 - t636 * t800) * pkin(8);
t627 = mrSges(2,1) * g(3) + mrSges(2,3) * t790 + t810 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t635 - t799 * t629 + t800 * t815;
t1 = [-m(1) * g(1) + t821; -m(1) * g(2) + t832; (-m(1) - m(2)) * g(3) + t635; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t832 - t805 * t627 + t809 * t628; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t821 + t809 * t627 + t805 * t628; -mrSges(1,1) * g(2) + mrSges(2,1) * t789 + mrSges(1,2) * g(1) - mrSges(2,2) * t790 + Ifges(2,3) * qJDD(1) + pkin(1) * t636 + t800 * t629 + t799 * t815;];
tauB  = t1;
