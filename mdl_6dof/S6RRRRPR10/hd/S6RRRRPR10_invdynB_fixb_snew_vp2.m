% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-07 23:04
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 22:51:22
% EndTime: 2019-05-07 22:51:55
% DurationCPUTime: 21.47s
% Computational Cost: add. (350364->379), mult. (751349->472), div. (0->0), fcn. (601093->12), ass. (0->157)
t826 = Ifges(5,1) + Ifges(6,2);
t820 = Ifges(5,4) + Ifges(6,6);
t819 = Ifges(5,5) - Ifges(6,4);
t825 = -Ifges(5,2) - Ifges(6,3);
t818 = Ifges(5,6) - Ifges(6,5);
t824 = Ifges(5,3) + Ifges(6,1);
t823 = -2 * qJD(5);
t822 = cos(qJ(4));
t780 = cos(pkin(6));
t821 = g(3) * t780;
t776 = qJD(1) * t780 + qJD(2);
t783 = sin(qJ(3));
t787 = cos(qJ(3));
t784 = sin(qJ(2));
t779 = sin(pkin(6));
t806 = qJD(1) * t779;
t803 = t784 * t806;
t753 = t776 * t787 - t783 * t803;
t754 = t776 * t783 + t787 * t803;
t782 = sin(qJ(4));
t738 = -t822 * t753 + t754 * t782;
t788 = cos(qJ(2));
t805 = qJD(1) * t788;
t802 = t779 * t805;
t770 = qJD(3) - t802;
t768 = -qJD(4) - t770;
t817 = t738 * t768;
t816 = t779 * t784;
t815 = t779 * t788;
t814 = t780 * t784;
t813 = t780 * t788;
t785 = sin(qJ(1));
t789 = cos(qJ(1));
t771 = t785 * g(1) - g(2) * t789;
t790 = qJD(1) ^ 2;
t761 = pkin(8) * t779 * t790 + qJDD(1) * pkin(1) + t771;
t772 = -g(1) * t789 - g(2) * t785;
t804 = qJDD(1) * t779;
t762 = -pkin(1) * t790 + pkin(8) * t804 + t772;
t807 = t761 * t814 + t788 * t762;
t734 = -g(3) * t816 + t807;
t759 = mrSges(3,1) * t776 - mrSges(3,3) * t803;
t763 = (-mrSges(3,1) * t788 + mrSges(3,2) * t784) * t806;
t766 = -qJD(2) * t803 + t788 * t804;
t775 = qJDD(1) * t780 + qJDD(2);
t764 = (-pkin(2) * t788 - pkin(9) * t784) * t806;
t774 = t776 ^ 2;
t713 = -pkin(2) * t774 + pkin(9) * t775 + (-g(3) * t784 + t764 * t805) * t779 + t807;
t765 = (qJD(2) * t805 + qJDD(1) * t784) * t779;
t714 = -pkin(2) * t766 - pkin(9) * t765 - t821 + (-t761 + (pkin(2) * t784 - pkin(9) * t788) * t776 * qJD(1)) * t779;
t679 = -t713 * t783 + t787 * t714;
t732 = qJD(3) * t753 + t765 * t787 + t775 * t783;
t758 = qJDD(3) - t766;
t671 = (t753 * t770 - t732) * pkin(10) + (t753 * t754 + t758) * pkin(3) + t679;
t680 = t787 * t713 + t783 * t714;
t731 = -qJD(3) * t754 - t765 * t783 + t775 * t787;
t743 = pkin(3) * t770 - pkin(10) * t754;
t752 = t753 ^ 2;
t674 = -pkin(3) * t752 + pkin(10) * t731 - t743 * t770 + t680;
t668 = t822 * t671 - t782 * t674;
t691 = -t738 * qJD(4) + t782 * t731 + t822 * t732;
t739 = t782 * t753 + t822 * t754;
t707 = mrSges(5,1) * t738 + mrSges(5,2) * t739;
t718 = mrSges(6,1) * t738 + mrSges(6,3) * t768;
t720 = mrSges(5,2) * t768 - mrSges(5,3) * t738;
t757 = qJDD(4) + t758;
t706 = pkin(4) * t738 - qJ(5) * t739;
t767 = t768 ^ 2;
t665 = -t757 * pkin(4) - t767 * qJ(5) + t739 * t706 + qJDD(5) - t668;
t660 = (t738 * t739 - t757) * pkin(11) + (t691 - t817) * pkin(5) + t665;
t690 = qJD(4) * t739 - t822 * t731 + t732 * t782;
t722 = pkin(5) * t739 + pkin(11) * t768;
t737 = t738 ^ 2;
t733 = -g(3) * t815 + t761 * t813 - t784 * t762;
t712 = -pkin(2) * t775 - pkin(9) * t774 + t764 * t803 - t733;
t678 = -pkin(3) * t731 - pkin(10) * t752 + t754 * t743 + t712;
t792 = (-t691 - t817) * qJ(5) + t678 + (-pkin(4) * t768 + t823) * t739;
t663 = t792 + (pkin(4) + pkin(11)) * t690 - pkin(5) * t737 - t722 * t739;
t781 = sin(qJ(6));
t786 = cos(qJ(6));
t658 = t660 * t786 - t663 * t781;
t716 = t738 * t786 + t768 * t781;
t677 = qJD(6) * t716 + t690 * t781 + t757 * t786;
t689 = qJDD(6) + t691;
t717 = t738 * t781 - t768 * t786;
t694 = -mrSges(7,1) * t716 + mrSges(7,2) * t717;
t736 = qJD(6) + t739;
t695 = -mrSges(7,2) * t736 + mrSges(7,3) * t716;
t656 = m(7) * t658 + mrSges(7,1) * t689 - mrSges(7,3) * t677 - t694 * t717 + t695 * t736;
t659 = t660 * t781 + t663 * t786;
t676 = -qJD(6) * t717 + t690 * t786 - t757 * t781;
t696 = mrSges(7,1) * t736 - mrSges(7,3) * t717;
t657 = m(7) * t659 - mrSges(7,2) * t689 + mrSges(7,3) * t676 + t694 * t716 - t696 * t736;
t648 = t656 * t786 + t657 * t781;
t708 = -mrSges(6,2) * t738 - mrSges(6,3) * t739;
t796 = -m(6) * t665 - t691 * mrSges(6,1) - t739 * t708 - t648;
t646 = m(5) * t668 - mrSges(5,3) * t691 - t707 * t739 + (t718 - t720) * t768 + (mrSges(5,1) - mrSges(6,2)) * t757 + t796;
t669 = t782 * t671 + t822 * t674;
t721 = -mrSges(5,1) * t768 - mrSges(5,3) * t739;
t795 = -pkin(4) * t767 + qJ(5) * t757 - t706 * t738 + t669;
t664 = 0.2e1 * qJD(5) * t768 - t795;
t719 = mrSges(6,1) * t739 - mrSges(6,2) * t768;
t662 = -pkin(5) * t690 - pkin(11) * t737 + (t823 - t722) * t768 + t795;
t797 = -m(7) * t662 + mrSges(7,1) * t676 - t677 * mrSges(7,2) + t695 * t716 - t717 * t696;
t794 = -m(6) * t664 + t757 * mrSges(6,3) - t768 * t719 - t797;
t653 = m(5) * t669 - mrSges(5,2) * t757 + t721 * t768 + (-t707 - t708) * t738 + (-mrSges(5,3) - mrSges(6,1)) * t690 + t794;
t641 = t822 * t646 + t782 * t653;
t740 = -mrSges(4,1) * t753 + mrSges(4,2) * t754;
t741 = -mrSges(4,2) * t770 + mrSges(4,3) * t753;
t639 = m(4) * t679 + mrSges(4,1) * t758 - mrSges(4,3) * t732 - t740 * t754 + t741 * t770 + t641;
t742 = mrSges(4,1) * t770 - mrSges(4,3) * t754;
t799 = -t646 * t782 + t822 * t653;
t640 = m(4) * t680 - mrSges(4,2) * t758 + mrSges(4,3) * t731 + t740 * t753 - t742 * t770 + t799;
t800 = -t639 * t783 + t787 * t640;
t631 = m(3) * t734 - mrSges(3,2) * t775 + mrSges(3,3) * t766 - t759 * t776 + t763 * t802 + t800;
t634 = t787 * t639 + t783 * t640;
t747 = -t761 * t779 - t821;
t760 = -mrSges(3,2) * t776 + mrSges(3,3) * t802;
t633 = m(3) * t747 - mrSges(3,1) * t766 + mrSges(3,2) * t765 + (t759 * t784 - t760 * t788) * t806 + t634;
t667 = pkin(4) * t690 + t792;
t811 = -t781 * t656 + t786 * t657;
t647 = m(6) * t667 - t690 * mrSges(6,2) - t691 * mrSges(6,3) - t738 * t718 - t739 * t719 + t811;
t793 = m(5) * t678 + t690 * mrSges(5,1) + t691 * mrSges(5,2) + t738 * t720 + t739 * t721 + t647;
t791 = -m(4) * t712 + t731 * mrSges(4,1) - t732 * mrSges(4,2) + t753 * t741 - t754 * t742 - t793;
t644 = m(3) * t733 + t775 * mrSges(3,1) - t765 * mrSges(3,3) + t776 * t760 - t763 * t803 + t791;
t622 = t631 * t814 - t633 * t779 + t644 * t813;
t620 = m(2) * t771 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t790 + t622;
t626 = t788 * t631 - t644 * t784;
t625 = m(2) * t772 - mrSges(2,1) * t790 - qJDD(1) * mrSges(2,2) + t626;
t812 = t789 * t620 + t785 * t625;
t810 = t738 * t818 - t739 * t819 + t768 * t824;
t809 = t738 * t825 + t739 * t820 - t768 * t818;
t808 = t820 * t738 - t739 * t826 + t819 * t768;
t621 = t631 * t816 + t780 * t633 + t644 * t815;
t801 = -t620 * t785 + t789 * t625;
t681 = Ifges(7,5) * t717 + Ifges(7,6) * t716 + Ifges(7,3) * t736;
t683 = Ifges(7,1) * t717 + Ifges(7,4) * t716 + Ifges(7,5) * t736;
t649 = -mrSges(7,1) * t662 + mrSges(7,3) * t659 + Ifges(7,4) * t677 + Ifges(7,2) * t676 + Ifges(7,6) * t689 - t681 * t717 + t683 * t736;
t682 = Ifges(7,4) * t717 + Ifges(7,2) * t716 + Ifges(7,6) * t736;
t650 = mrSges(7,2) * t662 - mrSges(7,3) * t658 + Ifges(7,1) * t677 + Ifges(7,4) * t676 + Ifges(7,5) * t689 + t681 * t716 - t682 * t736;
t627 = -mrSges(5,1) * t678 - mrSges(6,1) * t664 + mrSges(6,2) * t667 + mrSges(5,3) * t669 - pkin(4) * t647 - pkin(5) * t797 - pkin(11) * t811 - t786 * t649 - t781 * t650 + t690 * t825 + t820 * t691 + t810 * t739 + t818 * t757 + t808 * t768;
t635 = mrSges(6,1) * t665 + mrSges(7,1) * t658 + mrSges(5,2) * t678 - mrSges(7,2) * t659 - mrSges(5,3) * t668 - mrSges(6,3) * t667 + Ifges(7,5) * t677 + Ifges(7,6) * t676 + Ifges(7,3) * t689 + pkin(5) * t648 - qJ(5) * t647 + t717 * t682 - t716 * t683 + t809 * t768 + t819 * t757 + t810 * t738 + t826 * t691 - t820 * t690;
t725 = Ifges(4,5) * t754 + Ifges(4,6) * t753 + Ifges(4,3) * t770;
t727 = Ifges(4,1) * t754 + Ifges(4,4) * t753 + Ifges(4,5) * t770;
t616 = -mrSges(4,1) * t712 + mrSges(4,3) * t680 + Ifges(4,4) * t732 + Ifges(4,2) * t731 + Ifges(4,6) * t758 - pkin(3) * t793 + pkin(10) * t799 + t822 * t627 + t782 * t635 - t754 * t725 + t770 * t727;
t726 = Ifges(4,4) * t754 + Ifges(4,2) * t753 + Ifges(4,6) * t770;
t618 = mrSges(4,2) * t712 - mrSges(4,3) * t679 + Ifges(4,1) * t732 + Ifges(4,4) * t731 + Ifges(4,5) * t758 - pkin(10) * t641 - t782 * t627 + t822 * t635 + t753 * t725 - t770 * t726;
t744 = Ifges(3,3) * t776 + (Ifges(3,5) * t784 + Ifges(3,6) * t788) * t806;
t745 = Ifges(3,6) * t776 + (Ifges(3,4) * t784 + Ifges(3,2) * t788) * t806;
t615 = mrSges(3,2) * t747 - mrSges(3,3) * t733 + Ifges(3,1) * t765 + Ifges(3,4) * t766 + Ifges(3,5) * t775 - pkin(9) * t634 - t616 * t783 + t618 * t787 + t744 * t802 - t745 * t776;
t746 = Ifges(3,5) * t776 + (Ifges(3,1) * t784 + Ifges(3,4) * t788) * t806;
t617 = -t819 * t691 + (mrSges(6,1) * qJ(5) + t818) * t690 - Ifges(4,6) * t731 - qJ(5) * t794 - pkin(4) * (t718 * t768 + t796) + (qJ(5) * t708 + t808) * t738 - t809 * t739 + (pkin(4) * mrSges(6,2) - t824) * t757 - t744 * t803 + mrSges(6,3) * t664 - mrSges(6,2) * t665 - mrSges(4,1) * t679 + mrSges(4,2) * t680 + pkin(11) * t648 - mrSges(5,1) * t668 - pkin(3) * t641 + mrSges(5,2) * t669 - pkin(2) * t634 - Ifges(4,5) * t732 + mrSges(3,3) * t734 - mrSges(3,1) * t747 + t753 * t727 - t754 * t726 - Ifges(4,3) * t758 + Ifges(3,4) * t765 + Ifges(3,2) * t766 + Ifges(3,6) * t775 + t776 * t746 + t781 * t649 - t786 * t650;
t798 = pkin(8) * t626 + t615 * t784 + t617 * t788;
t614 = Ifges(3,5) * t765 + Ifges(3,6) * t766 + Ifges(3,3) * t775 + mrSges(3,1) * t733 - mrSges(3,2) * t734 + t783 * t618 + t787 * t616 + pkin(2) * t791 + pkin(9) * t800 + (t745 * t784 - t746 * t788) * t806;
t613 = -mrSges(2,2) * g(3) - mrSges(2,3) * t771 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t790 + t615 * t788 - t617 * t784 + (-t621 * t779 - t622 * t780) * pkin(8);
t612 = mrSges(2,1) * g(3) + mrSges(2,3) * t772 + Ifges(2,5) * t790 + Ifges(2,6) * qJDD(1) - pkin(1) * t621 - t614 * t779 + t780 * t798;
t1 = [-m(1) * g(1) + t801; -m(1) * g(2) + t812; (-m(1) - m(2)) * g(3) + t621; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t812 - t785 * t612 + t789 * t613; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t801 + t789 * t612 + t785 * t613; -mrSges(1,1) * g(2) + mrSges(2,1) * t771 + mrSges(1,2) * g(1) - mrSges(2,2) * t772 + Ifges(2,3) * qJDD(1) + pkin(1) * t622 + t614 * t780 + t779 * t798;];
tauB  = t1;
