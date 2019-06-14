% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 12:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPPRR11_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR11_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:00:00
% EndTime: 2019-05-06 12:00:17
% DurationCPUTime: 16.76s
% Computational Cost: add. (256591->382), mult. (605823->479), div. (0->0), fcn. (449112->12), ass. (0->159)
t831 = -2 * qJD(3);
t830 = Ifges(3,1) + Ifges(4,2);
t824 = Ifges(3,4) + Ifges(4,6);
t823 = Ifges(3,5) - Ifges(4,4);
t829 = Ifges(3,2) + Ifges(4,3);
t822 = Ifges(3,6) - Ifges(4,5);
t828 = Ifges(3,3) + Ifges(4,1);
t780 = cos(pkin(6));
t772 = t780 * qJD(1) + qJD(2);
t783 = sin(qJ(2));
t778 = sin(pkin(6));
t810 = qJD(1) * t778;
t803 = t783 * t810;
t827 = (pkin(2) * t772 + t831) * t803;
t784 = sin(qJ(1));
t788 = cos(qJ(1));
t767 = t784 * g(1) - t788 * g(2);
t789 = qJD(1) ^ 2;
t753 = t789 * t778 * pkin(8) + qJDD(1) * pkin(1) + t767;
t768 = -t788 * g(1) - t784 * g(2);
t807 = qJDD(1) * t778;
t754 = -t789 * pkin(1) + pkin(8) * t807 + t768;
t787 = cos(qJ(2));
t817 = t780 * t783;
t819 = t778 * t783;
t713 = -g(3) * t819 + t753 * t817 + t787 * t754;
t755 = (-pkin(2) * t787 - qJ(3) * t783) * t810;
t770 = t772 ^ 2;
t771 = t780 * qJDD(1) + qJDD(2);
t809 = qJD(1) * t787;
t804 = t778 * t809;
t700 = t770 * pkin(2) - t771 * qJ(3) - t755 * t804 + t772 * t831 - t713;
t826 = t780 * g(3);
t825 = mrSges(3,1) - mrSges(4,2);
t821 = -pkin(2) - qJ(4);
t820 = t778 ^ 2 * t789;
t818 = t778 * t787;
t816 = t780 * t787;
t730 = -t778 * t753 - t826;
t748 = t772 * mrSges(3,1) - mrSges(3,3) * t803;
t749 = -t772 * mrSges(3,2) + mrSges(3,3) * t804;
t752 = mrSges(4,1) * t803 + t772 * mrSges(4,2);
t758 = (qJD(2) * t809 + qJDD(1) * t783) * t778;
t759 = -qJD(2) * t803 + t787 * t807;
t701 = -t759 * pkin(2) + (-t772 * t804 - t758) * qJ(3) + t730 + t827;
t751 = -mrSges(4,1) * t804 - t772 * mrSges(4,3);
t750 = pkin(3) * t803 - t772 * qJ(4);
t806 = t787 ^ 2 * t820;
t682 = -pkin(3) * t806 - t826 - t758 * qJ(3) + t821 * t759 + (-t753 + (-qJ(3) * t772 * t787 - t750 * t783) * qJD(1)) * t778 + t827;
t811 = g(3) * t818 + t783 * t754;
t797 = -t770 * qJ(3) + t755 * t803 + qJDD(3) + t811;
t685 = t758 * pkin(3) + t821 * t771 + (-pkin(3) * t772 * t810 - qJ(4) * t783 * t820 - t753 * t780) * t787 + t797;
t777 = sin(pkin(11));
t779 = cos(pkin(11));
t742 = t779 * t772 - t777 * t804;
t668 = -0.2e1 * qJD(4) * t742 - t777 * t682 + t779 * t685;
t722 = -t777 * t759 + t779 * t771;
t741 = -t777 * t772 - t779 * t804;
t665 = (t741 * t803 - t722) * pkin(9) + (t741 * t742 + t758) * pkin(4) + t668;
t669 = 0.2e1 * qJD(4) * t741 + t779 * t682 + t777 * t685;
t721 = -t779 * t759 - t777 * t771;
t723 = pkin(4) * t803 - t742 * pkin(9);
t740 = t741 ^ 2;
t667 = -t740 * pkin(4) + t721 * pkin(9) - t723 * t803 + t669;
t782 = sin(qJ(5));
t786 = cos(qJ(5));
t662 = t782 * t665 + t786 * t667;
t716 = t782 * t741 + t786 * t742;
t689 = -t716 * qJD(5) + t786 * t721 - t782 * t722;
t715 = t786 * t741 - t782 * t742;
t698 = -t715 * mrSges(6,1) + t716 * mrSges(6,2);
t763 = qJD(5) + t803;
t706 = t763 * mrSges(6,1) - t716 * mrSges(6,3);
t747 = qJDD(5) + t758;
t699 = -t715 * pkin(5) - t716 * pkin(10);
t761 = t763 ^ 2;
t660 = -t761 * pkin(5) + t747 * pkin(10) + t715 * t699 + t662;
t681 = t759 * pkin(3) - qJ(4) * t806 + t772 * t750 + qJDD(4) - t700;
t671 = -t721 * pkin(4) - t740 * pkin(9) + t742 * t723 + t681;
t690 = t715 * qJD(5) + t782 * t721 + t786 * t722;
t663 = t671 + (-t715 * t763 - t690) * pkin(10) + (t716 * t763 - t689) * pkin(5);
t781 = sin(qJ(6));
t785 = cos(qJ(6));
t657 = -t781 * t660 + t785 * t663;
t703 = -t781 * t716 + t785 * t763;
t674 = t703 * qJD(6) + t785 * t690 + t781 * t747;
t704 = t785 * t716 + t781 * t763;
t686 = -t703 * mrSges(7,1) + t704 * mrSges(7,2);
t688 = qJDD(6) - t689;
t714 = qJD(6) - t715;
t691 = -t714 * mrSges(7,2) + t703 * mrSges(7,3);
t655 = m(7) * t657 + t688 * mrSges(7,1) - t674 * mrSges(7,3) - t704 * t686 + t714 * t691;
t658 = t785 * t660 + t781 * t663;
t673 = -t704 * qJD(6) - t781 * t690 + t785 * t747;
t692 = t714 * mrSges(7,1) - t704 * mrSges(7,3);
t656 = m(7) * t658 - t688 * mrSges(7,2) + t673 * mrSges(7,3) + t703 * t686 - t714 * t692;
t799 = -t781 * t655 + t785 * t656;
t646 = m(6) * t662 - t747 * mrSges(6,2) + t689 * mrSges(6,3) + t715 * t698 - t763 * t706 + t799;
t661 = t786 * t665 - t782 * t667;
t705 = -t763 * mrSges(6,2) + t715 * mrSges(6,3);
t659 = -t747 * pkin(5) - t761 * pkin(10) + t716 * t699 - t661;
t794 = -m(7) * t659 + t673 * mrSges(7,1) - t674 * mrSges(7,2) + t703 * t691 - t704 * t692;
t651 = m(6) * t661 + t747 * mrSges(6,1) - t690 * mrSges(6,3) - t716 * t698 + t763 * t705 + t794;
t639 = t782 * t646 + t786 * t651;
t717 = -t741 * mrSges(5,1) + t742 * mrSges(5,2);
t719 = -mrSges(5,2) * t803 + t741 * mrSges(5,3);
t637 = m(5) * t668 + t758 * mrSges(5,1) - t722 * mrSges(5,3) - t742 * t717 + t719 * t803 + t639;
t720 = mrSges(5,1) * t803 - t742 * mrSges(5,3);
t800 = t786 * t646 - t782 * t651;
t638 = m(5) * t669 - t758 * mrSges(5,2) + t721 * mrSges(5,3) + t741 * t717 - t720 * t803 + t800;
t801 = -t777 * t637 + t779 * t638;
t798 = m(4) * t701 - t758 * mrSges(4,3) + t751 * t804 + t801;
t630 = m(3) * t730 + t758 * mrSges(3,2) - t825 * t759 + (-t749 * t787 + (t748 - t752) * t783) * t810 + t798;
t805 = t753 * t816;
t712 = t805 - t811;
t756 = (mrSges(4,2) * t787 - mrSges(4,3) * t783) * t810;
t757 = (-mrSges(3,1) * t787 + mrSges(3,2) * t783) * t810;
t633 = t779 * t637 + t777 * t638;
t702 = -t771 * pkin(2) + t797 - t805;
t795 = -m(4) * t702 - t758 * mrSges(4,1) - t633;
t631 = m(3) * t712 - t758 * mrSges(3,3) + (t749 - t751) * t772 + t825 * t771 + (-t756 - t757) * t803 + t795;
t647 = t785 * t655 + t781 * t656;
t792 = m(6) * t671 - t689 * mrSges(6,1) + t690 * mrSges(6,2) - t715 * t705 + t716 * t706 + t647;
t791 = -m(5) * t681 + t721 * mrSges(5,1) - t722 * mrSges(5,2) + t741 * t719 - t742 * t720 - t792;
t790 = -m(4) * t700 + t771 * mrSges(4,3) + t772 * t752 + t756 * t804 - t791;
t643 = m(3) * t713 + t757 * t804 + t790 + (mrSges(3,3) + mrSges(4,1)) * t759 - t771 * mrSges(3,2) - t772 * t748;
t620 = -t778 * t630 + t631 * t816 + t643 * t817;
t618 = m(2) * t767 + qJDD(1) * mrSges(2,1) - t789 * mrSges(2,2) + t620;
t626 = -t783 * t631 + t787 * t643;
t625 = m(2) * t768 - t789 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t626;
t815 = t788 * t618 + t784 * t625;
t814 = (t823 * t783 + t822 * t787) * t810 + t828 * t772;
t813 = (-t824 * t783 - t829 * t787) * t810 - t822 * t772;
t812 = (t830 * t783 + t824 * t787) * t810 + t823 * t772;
t619 = t780 * t630 + t631 * t818 + t643 * t819;
t802 = -t784 * t618 + t788 * t625;
t675 = Ifges(7,5) * t704 + Ifges(7,6) * t703 + Ifges(7,3) * t714;
t677 = Ifges(7,1) * t704 + Ifges(7,4) * t703 + Ifges(7,5) * t714;
t648 = -mrSges(7,1) * t659 + mrSges(7,3) * t658 + Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t688 - t704 * t675 + t714 * t677;
t676 = Ifges(7,4) * t704 + Ifges(7,2) * t703 + Ifges(7,6) * t714;
t649 = mrSges(7,2) * t659 - mrSges(7,3) * t657 + Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t688 + t703 * t675 - t714 * t676;
t693 = Ifges(6,5) * t716 + Ifges(6,6) * t715 + Ifges(6,3) * t763;
t694 = Ifges(6,4) * t716 + Ifges(6,2) * t715 + Ifges(6,6) * t763;
t634 = mrSges(6,2) * t671 - mrSges(6,3) * t661 + Ifges(6,1) * t690 + Ifges(6,4) * t689 + Ifges(6,5) * t747 - pkin(10) * t647 - t781 * t648 + t785 * t649 + t715 * t693 - t763 * t694;
t695 = Ifges(6,1) * t716 + Ifges(6,4) * t715 + Ifges(6,5) * t763;
t635 = -mrSges(6,1) * t671 - mrSges(7,1) * t657 + mrSges(7,2) * t658 + mrSges(6,3) * t662 + Ifges(6,4) * t690 - Ifges(7,5) * t674 + Ifges(6,2) * t689 + Ifges(6,6) * t747 - Ifges(7,6) * t673 - Ifges(7,3) * t688 - pkin(5) * t647 - t704 * t676 + t703 * t677 - t716 * t693 + t763 * t695;
t707 = Ifges(5,5) * t742 + Ifges(5,6) * t741 + Ifges(5,3) * t803;
t709 = Ifges(5,1) * t742 + Ifges(5,4) * t741 + Ifges(5,5) * t803;
t621 = -mrSges(5,1) * t681 + mrSges(5,3) * t669 + Ifges(5,4) * t722 + Ifges(5,2) * t721 + Ifges(5,6) * t758 - pkin(4) * t792 + pkin(9) * t800 + t782 * t634 + t786 * t635 - t742 * t707 + t709 * t803;
t708 = Ifges(5,4) * t742 + Ifges(5,2) * t741 + Ifges(5,6) * t803;
t622 = mrSges(5,2) * t681 - mrSges(5,3) * t668 + Ifges(5,1) * t722 + Ifges(5,4) * t721 + Ifges(5,5) * t758 - pkin(9) * t639 + t786 * t634 - t782 * t635 + t741 * t707 - t708 * t803;
t632 = t759 * mrSges(4,2) - t752 * t803 + t798;
t615 = -mrSges(3,1) * t730 - mrSges(4,1) * t700 + mrSges(4,2) * t701 + mrSges(3,3) * t713 - pkin(2) * t632 - pkin(3) * t791 - qJ(4) * t801 - t779 * t621 - t777 * t622 + t824 * t758 + t829 * t759 + t822 * t771 + t812 * t772 - t814 * t803;
t616 = -t715 * t695 + t716 * t694 - mrSges(3,3) * t712 + Ifges(6,5) * t690 - mrSges(4,3) * t701 + mrSges(4,1) * t702 + Ifges(6,6) * t689 + mrSges(5,1) * t668 - mrSges(5,2) * t669 + t823 * t771 + t824 * t759 + mrSges(6,1) * t661 - mrSges(6,2) * t662 + pkin(3) * t633 - qJ(3) * t632 + t785 * t648 + t813 * t772 + pkin(5) * t794 + (Ifges(5,3) + t830) * t758 + t781 * t649 + Ifges(5,6) * t721 + Ifges(5,5) * t722 + mrSges(3,2) * t730 - t741 * t709 + t742 * t708 + Ifges(6,3) * t747 + pkin(4) * t639 + t814 * t804 + pkin(10) * t799;
t796 = pkin(8) * t626 + t615 * t787 + t616 * t783;
t614 = mrSges(3,1) * t712 - mrSges(3,2) * t713 + mrSges(4,2) * t702 - mrSges(4,3) * t700 + t779 * t622 - t777 * t621 - qJ(4) * t633 + pkin(2) * (-t772 * t751 + t795) + qJ(3) * t790 + (-pkin(2) * mrSges(4,2) + t828) * t771 + (qJ(3) * mrSges(4,1) + t822) * t759 + t823 * t758 + (-t812 * t787 + (-pkin(2) * t756 - t813) * t783) * t810;
t613 = -mrSges(2,2) * g(3) - mrSges(2,3) * t767 + Ifges(2,5) * qJDD(1) - t789 * Ifges(2,6) - t783 * t615 + t787 * t616 + (-t619 * t778 - t620 * t780) * pkin(8);
t612 = mrSges(2,1) * g(3) + mrSges(2,3) * t768 + t789 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t619 - t778 * t614 + t780 * t796;
t1 = [-m(1) * g(1) + t802; -m(1) * g(2) + t815; (-m(1) - m(2)) * g(3) + t619; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t815 - t784 * t612 + t788 * t613; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t802 + t788 * t612 + t784 * t613; -mrSges(1,1) * g(2) + mrSges(2,1) * t767 + mrSges(1,2) * g(1) - mrSges(2,2) * t768 + Ifges(2,3) * qJDD(1) + pkin(1) * t620 + t780 * t614 + t778 * t796;];
tauB  = t1;
