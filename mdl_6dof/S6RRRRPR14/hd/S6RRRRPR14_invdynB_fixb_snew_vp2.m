% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 02:41
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR14_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 02:14:53
% EndTime: 2019-05-08 02:18:15
% DurationCPUTime: 148.47s
% Computational Cost: add. (2368159->416), mult. (5872133->556), div. (0->0), fcn. (4979041->16), ass. (0->179)
t793 = cos(pkin(6));
t786 = qJD(1) * t793 + qJD(2);
t789 = sin(pkin(7));
t792 = cos(pkin(7));
t790 = sin(pkin(6));
t801 = cos(qJ(2));
t822 = qJD(1) * t801;
t818 = t790 * t822;
t771 = (t786 * t789 + t792 * t818) * pkin(10);
t797 = sin(qJ(2));
t824 = qJD(1) * t790;
t837 = pkin(10) * t789;
t775 = (-pkin(2) * t801 - t797 * t837) * t824;
t821 = qJD(1) * qJD(2);
t781 = (qJDD(1) * t797 + t801 * t821) * t790;
t785 = qJDD(1) * t793 + qJDD(2);
t798 = sin(qJ(1));
t802 = cos(qJ(1));
t783 = t798 * g(1) - g(2) * t802;
t803 = qJD(1) ^ 2;
t838 = pkin(9) * t790;
t778 = qJDD(1) * pkin(1) + t803 * t838 + t783;
t784 = -g(1) * t802 - g(2) * t798;
t779 = -pkin(1) * t803 + qJDD(1) * t838 + t784;
t827 = t793 * t801;
t813 = t778 * t827 - t797 * t779;
t823 = qJD(1) * t797;
t836 = pkin(10) * t792;
t727 = -t781 * t836 + t785 * pkin(2) + t786 * t771 + (-g(3) * t801 - t775 * t823) * t790 + t813;
t819 = t790 * t823;
t774 = pkin(2) * t786 - t819 * t836;
t782 = (qJDD(1) * t801 - t797 * t821) * t790;
t811 = t782 * t792 + t785 * t789;
t828 = t793 * t797;
t825 = t778 * t828 + t801 * t779;
t728 = -t786 * t774 + (-g(3) * t797 + t775 * t822) * t790 + t811 * pkin(10) + t825;
t835 = t793 * g(3);
t733 = -t781 * t837 - t782 * pkin(2) - t835 + (-t778 + (-t771 * t801 + t774 * t797) * qJD(1)) * t790;
t796 = sin(qJ(3));
t800 = cos(qJ(3));
t697 = -t796 * t728 + (t727 * t792 + t733 * t789) * t800;
t829 = t792 * t801;
t834 = t789 * t796;
t761 = t786 * t834 + (t796 * t829 + t797 * t800) * t824;
t746 = -t761 * qJD(3) - t796 * t781 + t800 * t811;
t833 = t789 * t800;
t760 = (-t796 * t797 + t800 * t829) * t824 + t786 * t833;
t839 = cos(qJ(4));
t832 = t790 * t797;
t831 = t790 * t801;
t830 = t792 * t796;
t698 = t727 * t830 + t800 * t728 + t733 * t834;
t748 = -mrSges(4,1) * t760 + mrSges(4,2) * t761;
t772 = t786 * t792 - t789 * t818 + qJD(3);
t754 = mrSges(4,1) * t772 - mrSges(4,3) * t761;
t762 = -t782 * t789 + t785 * t792 + qJDD(3);
t749 = -pkin(3) * t760 - pkin(11) * t761;
t770 = t772 ^ 2;
t689 = -pkin(3) * t770 + pkin(11) * t762 + t749 * t760 + t698;
t708 = -t789 * t727 + t792 * t733;
t747 = t760 * qJD(3) + t800 * t781 + t796 * t811;
t691 = (-t760 * t772 - t747) * pkin(11) + (t761 * t772 - t746) * pkin(3) + t708;
t795 = sin(qJ(4));
t682 = t689 * t839 + t795 * t691;
t752 = t761 * t839 + t795 * t772;
t714 = qJD(4) * t752 + t747 * t795 - t762 * t839;
t751 = t761 * t795 - t772 * t839;
t730 = mrSges(5,1) * t751 + mrSges(5,2) * t752;
t758 = qJD(4) - t760;
t741 = mrSges(5,1) * t758 - mrSges(5,3) * t752;
t745 = qJDD(4) - t746;
t729 = pkin(4) * t751 - qJ(5) * t752;
t757 = t758 ^ 2;
t677 = -pkin(4) * t757 + qJ(5) * t745 - t729 * t751 + t682;
t688 = -t762 * pkin(3) - t770 * pkin(11) + t761 * t749 - t697;
t715 = -t751 * qJD(4) + t747 * t839 + t795 * t762;
t680 = (t751 * t758 - t715) * qJ(5) + (t752 * t758 + t714) * pkin(4) + t688;
t788 = sin(pkin(13));
t791 = cos(pkin(13));
t739 = t752 * t791 + t758 * t788;
t672 = -0.2e1 * qJD(5) * t739 - t788 * t677 + t791 * t680;
t701 = t715 * t791 + t745 * t788;
t738 = -t752 * t788 + t758 * t791;
t670 = (t738 * t751 - t701) * pkin(12) + (t738 * t739 + t714) * pkin(5) + t672;
t673 = 0.2e1 * qJD(5) * t738 + t791 * t677 + t788 * t680;
t700 = -t715 * t788 + t745 * t791;
t719 = pkin(5) * t751 - pkin(12) * t739;
t737 = t738 ^ 2;
t671 = -pkin(5) * t737 + pkin(12) * t700 - t719 * t751 + t673;
t794 = sin(qJ(6));
t799 = cos(qJ(6));
t668 = t670 * t799 - t671 * t794;
t709 = t738 * t799 - t739 * t794;
t685 = qJD(6) * t709 + t700 * t794 + t701 * t799;
t710 = t738 * t794 + t739 * t799;
t696 = -mrSges(7,1) * t709 + mrSges(7,2) * t710;
t750 = qJD(6) + t751;
t702 = -mrSges(7,2) * t750 + mrSges(7,3) * t709;
t713 = qJDD(6) + t714;
t666 = m(7) * t668 + mrSges(7,1) * t713 - mrSges(7,3) * t685 - t696 * t710 + t702 * t750;
t669 = t670 * t794 + t671 * t799;
t684 = -qJD(6) * t710 + t700 * t799 - t701 * t794;
t703 = mrSges(7,1) * t750 - mrSges(7,3) * t710;
t667 = m(7) * t669 - mrSges(7,2) * t713 + mrSges(7,3) * t684 + t696 * t709 - t703 * t750;
t658 = t799 * t666 + t794 * t667;
t711 = -mrSges(6,1) * t738 + mrSges(6,2) * t739;
t717 = -mrSges(6,2) * t751 + mrSges(6,3) * t738;
t656 = m(6) * t672 + mrSges(6,1) * t714 - mrSges(6,3) * t701 - t711 * t739 + t717 * t751 + t658;
t718 = mrSges(6,1) * t751 - mrSges(6,3) * t739;
t814 = -t666 * t794 + t799 * t667;
t657 = m(6) * t673 - mrSges(6,2) * t714 + mrSges(6,3) * t700 + t711 * t738 - t718 * t751 + t814;
t815 = -t656 * t788 + t791 * t657;
t653 = m(5) * t682 - mrSges(5,2) * t745 - mrSges(5,3) * t714 - t730 * t751 - t741 * t758 + t815;
t681 = -t795 * t689 + t691 * t839;
t740 = -mrSges(5,2) * t758 - mrSges(5,3) * t751;
t676 = -t745 * pkin(4) - t757 * qJ(5) + t752 * t729 + qJDD(5) - t681;
t674 = -t700 * pkin(5) - t737 * pkin(12) + t739 * t719 + t676;
t806 = m(7) * t674 - t684 * mrSges(7,1) + mrSges(7,2) * t685 - t709 * t702 + t703 * t710;
t804 = -m(6) * t676 + t700 * mrSges(6,1) - mrSges(6,2) * t701 + t738 * t717 - t718 * t739 - t806;
t662 = m(5) * t681 + mrSges(5,1) * t745 - mrSges(5,3) * t715 - t730 * t752 + t740 * t758 + t804;
t816 = t653 * t839 - t662 * t795;
t642 = m(4) * t698 - mrSges(4,2) * t762 + mrSges(4,3) * t746 + t748 * t760 - t754 * t772 + t816;
t645 = t795 * t653 + t662 * t839;
t753 = -mrSges(4,2) * t772 + mrSges(4,3) * t760;
t644 = m(4) * t708 - mrSges(4,1) * t746 + mrSges(4,2) * t747 - t753 * t760 + t754 * t761 + t645;
t654 = t656 * t791 + t657 * t788;
t805 = -m(5) * t688 - t714 * mrSges(5,1) - mrSges(5,2) * t715 - t751 * t740 - t741 * t752 - t654;
t650 = m(4) * t697 + mrSges(4,1) * t762 - mrSges(4,3) * t747 - t748 * t761 + t753 * t772 + t805;
t631 = t792 * t800 * t650 + t642 * t830 - t644 * t789;
t755 = -g(3) * t831 + t813;
t777 = -mrSges(3,2) * t786 + mrSges(3,3) * t818;
t780 = (-mrSges(3,1) * t801 + mrSges(3,2) * t797) * t824;
t627 = m(3) * t755 + mrSges(3,1) * t785 - mrSges(3,3) * t781 + t777 * t786 - t780 * t819 + t631;
t630 = t642 * t834 + t792 * t644 + t650 * t833;
t766 = -t790 * t778 - t835;
t776 = mrSges(3,1) * t786 - mrSges(3,3) * t819;
t629 = m(3) * t766 - t782 * mrSges(3,1) + t781 * mrSges(3,2) + (t776 * t797 - t777 * t801) * t824 + t630;
t637 = t800 * t642 - t650 * t796;
t756 = -g(3) * t832 + t825;
t636 = m(3) * t756 - mrSges(3,2) * t785 + mrSges(3,3) * t782 - t776 * t786 + t780 * t818 + t637;
t617 = t627 * t827 - t629 * t790 + t636 * t828;
t615 = m(2) * t783 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t803 + t617;
t623 = -t627 * t797 + t801 * t636;
t622 = m(2) * t784 - mrSges(2,1) * t803 - qJDD(1) * mrSges(2,2) + t623;
t826 = t802 * t615 + t798 * t622;
t616 = t627 * t831 + t793 * t629 + t636 * t832;
t817 = -t615 * t798 + t802 * t622;
t692 = Ifges(7,5) * t710 + Ifges(7,6) * t709 + Ifges(7,3) * t750;
t694 = Ifges(7,1) * t710 + Ifges(7,4) * t709 + Ifges(7,5) * t750;
t659 = -mrSges(7,1) * t674 + mrSges(7,3) * t669 + Ifges(7,4) * t685 + Ifges(7,2) * t684 + Ifges(7,6) * t713 - t692 * t710 + t694 * t750;
t693 = Ifges(7,4) * t710 + Ifges(7,2) * t709 + Ifges(7,6) * t750;
t660 = mrSges(7,2) * t674 - mrSges(7,3) * t668 + Ifges(7,1) * t685 + Ifges(7,4) * t684 + Ifges(7,5) * t713 + t692 * t709 - t693 * t750;
t704 = Ifges(6,5) * t739 + Ifges(6,6) * t738 + Ifges(6,3) * t751;
t706 = Ifges(6,1) * t739 + Ifges(6,4) * t738 + Ifges(6,5) * t751;
t646 = -mrSges(6,1) * t676 + mrSges(6,3) * t673 + Ifges(6,4) * t701 + Ifges(6,2) * t700 + Ifges(6,6) * t714 - pkin(5) * t806 + pkin(12) * t814 + t799 * t659 + t794 * t660 - t739 * t704 + t751 * t706;
t705 = Ifges(6,4) * t739 + Ifges(6,2) * t738 + Ifges(6,6) * t751;
t647 = mrSges(6,2) * t676 - mrSges(6,3) * t672 + Ifges(6,1) * t701 + Ifges(6,4) * t700 + Ifges(6,5) * t714 - pkin(12) * t658 - t659 * t794 + t660 * t799 + t704 * t738 - t705 * t751;
t720 = Ifges(5,5) * t752 - Ifges(5,6) * t751 + Ifges(5,3) * t758;
t721 = Ifges(5,4) * t752 - Ifges(5,2) * t751 + Ifges(5,6) * t758;
t632 = mrSges(5,2) * t688 - mrSges(5,3) * t681 + Ifges(5,1) * t715 - Ifges(5,4) * t714 + Ifges(5,5) * t745 - qJ(5) * t654 - t646 * t788 + t647 * t791 - t720 * t751 - t721 * t758;
t722 = Ifges(5,1) * t752 - Ifges(5,4) * t751 + Ifges(5,5) * t758;
t638 = Ifges(5,4) * t715 + Ifges(5,6) * t745 - t752 * t720 + t758 * t722 - mrSges(5,1) * t688 + mrSges(5,3) * t682 - Ifges(6,5) * t701 - Ifges(6,6) * t700 - t739 * t705 + t738 * t706 - mrSges(6,1) * t672 + mrSges(6,2) * t673 - Ifges(7,5) * t685 - Ifges(7,6) * t684 - Ifges(7,3) * t713 - t710 * t693 + t709 * t694 - mrSges(7,1) * t668 + mrSges(7,2) * t669 - pkin(5) * t658 - pkin(4) * t654 + (-Ifges(5,2) - Ifges(6,3)) * t714;
t743 = Ifges(4,4) * t761 + Ifges(4,2) * t760 + Ifges(4,6) * t772;
t744 = Ifges(4,1) * t761 + Ifges(4,4) * t760 + Ifges(4,5) * t772;
t618 = mrSges(4,1) * t697 - mrSges(4,2) * t698 + Ifges(4,5) * t747 + Ifges(4,6) * t746 + Ifges(4,3) * t762 + pkin(3) * t805 + pkin(11) * t816 + t795 * t632 + t638 * t839 + t761 * t743 - t760 * t744;
t763 = Ifges(3,3) * t786 + (Ifges(3,5) * t797 + Ifges(3,6) * t801) * t824;
t765 = Ifges(3,5) * t786 + (Ifges(3,1) * t797 + Ifges(3,4) * t801) * t824;
t742 = Ifges(4,5) * t761 + Ifges(4,6) * t760 + Ifges(4,3) * t772;
t619 = mrSges(4,2) * t708 - mrSges(4,3) * t697 + Ifges(4,1) * t747 + Ifges(4,4) * t746 + Ifges(4,5) * t762 - pkin(11) * t645 + t632 * t839 - t795 * t638 + t760 * t742 - t772 * t743;
t624 = Ifges(4,4) * t747 + Ifges(4,2) * t746 + Ifges(4,6) * t762 - t761 * t742 + t772 * t744 - mrSges(4,1) * t708 + mrSges(4,3) * t698 - Ifges(5,5) * t715 + Ifges(5,6) * t714 - Ifges(5,3) * t745 - t752 * t721 - t751 * t722 - mrSges(5,1) * t681 + mrSges(5,2) * t682 - t788 * t647 - t791 * t646 - pkin(4) * t804 - qJ(5) * t815 - pkin(3) * t645;
t807 = pkin(10) * t637 + t619 * t796 + t624 * t800;
t612 = -mrSges(3,1) * t766 + mrSges(3,3) * t756 + Ifges(3,4) * t781 + Ifges(3,2) * t782 + Ifges(3,6) * t785 - pkin(2) * t630 - t789 * t618 - t763 * t819 + t786 * t765 + t792 * t807;
t764 = Ifges(3,6) * t786 + (Ifges(3,4) * t797 + Ifges(3,2) * t801) * t824;
t613 = t763 * t818 + mrSges(3,2) * t766 - mrSges(3,3) * t755 + Ifges(3,1) * t781 + Ifges(3,4) * t782 + Ifges(3,5) * t785 + t800 * t619 - t796 * t624 - t786 * t764 + (-t630 * t789 - t631 * t792) * pkin(10);
t808 = pkin(9) * t623 + t612 * t801 + t613 * t797;
t611 = mrSges(3,1) * t755 - mrSges(3,2) * t756 + Ifges(3,5) * t781 + Ifges(3,6) * t782 + Ifges(3,3) * t785 + pkin(2) * t631 + t792 * t618 + (t764 * t797 - t765 * t801) * t824 + t807 * t789;
t610 = -mrSges(2,2) * g(3) - mrSges(2,3) * t783 + Ifges(2,5) * qJDD(1) - t803 * Ifges(2,6) - t797 * t612 + t801 * t613 + (-t616 * t790 - t617 * t793) * pkin(9);
t609 = mrSges(2,1) * g(3) + mrSges(2,3) * t784 + t803 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t616 - t790 * t611 + t793 * t808;
t1 = [-m(1) * g(1) + t817; -m(1) * g(2) + t826; (-m(1) - m(2)) * g(3) + t616; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(8) * t826 - t798 * t609 + t802 * t610; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(8) * t817 + t802 * t609 + t798 * t610; -mrSges(1,1) * g(2) + mrSges(2,1) * t783 + mrSges(1,2) * g(1) - mrSges(2,2) * t784 + Ifges(2,3) * qJDD(1) + pkin(1) * t617 + t793 * t611 + t790 * t808;];
tauB  = t1;
