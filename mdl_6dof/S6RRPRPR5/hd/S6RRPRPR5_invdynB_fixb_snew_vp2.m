% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 14:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:57:51
% EndTime: 2019-05-06 13:58:41
% DurationCPUTime: 48.63s
% Computational Cost: add. (730919->396), mult. (1936726->517), div. (0->0), fcn. (1543432->14), ass. (0->161)
t795 = -2 * qJD(3);
t755 = sin(pkin(11));
t758 = cos(pkin(11));
t762 = sin(qJ(2));
t765 = cos(qJ(2));
t756 = sin(pkin(6));
t786 = qJD(1) * t756;
t733 = (t755 * t762 - t758 * t765) * t786;
t784 = qJD(1) * qJD(2);
t742 = (qJDD(1) * t762 + t765 * t784) * t756;
t759 = cos(pkin(6));
t749 = qJDD(1) * t759 + qJDD(2);
t750 = qJD(1) * t759 + qJD(2);
t763 = sin(qJ(1));
t766 = cos(qJ(1));
t746 = t763 * g(1) - g(2) * t766;
t767 = qJD(1) ^ 2;
t793 = pkin(8) * t756;
t739 = qJDD(1) * pkin(1) + t767 * t793 + t746;
t747 = -g(1) * t766 - g(2) * t763;
t740 = -pkin(1) * t767 + qJDD(1) * t793 + t747;
t788 = t759 * t765;
t774 = t739 * t788 - t762 * t740;
t792 = t756 ^ 2 * t767;
t680 = t749 * pkin(2) - t742 * qJ(3) + (pkin(2) * t762 * t792 + (qJ(3) * qJD(1) * t750 - g(3)) * t756) * t765 + t774;
t789 = t759 * t762;
t791 = t756 * t762;
t710 = -g(3) * t791 + t739 * t789 + t765 * t740;
t782 = t762 * t786;
t736 = pkin(2) * t750 - qJ(3) * t782;
t743 = (qJDD(1) * t765 - t762 * t784) * t756;
t783 = t765 ^ 2 * t792;
t684 = -pkin(2) * t783 + qJ(3) * t743 - t736 * t750 + t710;
t734 = (t755 * t765 + t758 * t762) * t786;
t660 = t758 * t680 - t755 * t684 + t734 * t795;
t794 = cos(qJ(4));
t790 = t756 * t765;
t661 = t755 * t680 + t758 * t684 + t733 * t795;
t711 = mrSges(4,1) * t733 + mrSges(4,2) * t734;
t715 = -t742 * t755 + t743 * t758;
t721 = mrSges(4,1) * t750 - mrSges(4,3) * t734;
t712 = pkin(3) * t733 - pkin(9) * t734;
t748 = t750 ^ 2;
t654 = -pkin(3) * t748 + pkin(9) * t749 - t712 * t733 + t661;
t725 = -t759 * g(3) - t756 * t739;
t694 = -t743 * pkin(2) - qJ(3) * t783 + t736 * t782 + qJDD(3) + t725;
t716 = t742 * t758 + t743 * t755;
t663 = (t733 * t750 - t716) * pkin(9) + (t734 * t750 - t715) * pkin(3) + t694;
t761 = sin(qJ(4));
t647 = t654 * t794 + t761 * t663;
t719 = t734 * t794 + t761 * t750;
t691 = qJD(4) * t719 + t716 * t761 - t749 * t794;
t718 = t734 * t761 - t750 * t794;
t696 = mrSges(5,1) * t718 + mrSges(5,2) * t719;
t732 = qJD(4) + t733;
t704 = mrSges(5,1) * t732 - mrSges(5,3) * t719;
t714 = qJDD(4) - t715;
t695 = pkin(4) * t718 - qJ(5) * t719;
t731 = t732 ^ 2;
t642 = -pkin(4) * t731 + qJ(5) * t714 - t695 * t718 + t647;
t653 = -t749 * pkin(3) - t748 * pkin(9) + t734 * t712 - t660;
t692 = -t718 * qJD(4) + t716 * t794 + t761 * t749;
t645 = (t718 * t732 - t692) * qJ(5) + (t719 * t732 + t691) * pkin(4) + t653;
t754 = sin(pkin(12));
t757 = cos(pkin(12));
t702 = t719 * t757 + t732 * t754;
t637 = -0.2e1 * qJD(5) * t702 - t754 * t642 + t757 * t645;
t672 = t692 * t757 + t714 * t754;
t701 = -t719 * t754 + t732 * t757;
t635 = (t701 * t718 - t672) * pkin(10) + (t701 * t702 + t691) * pkin(5) + t637;
t638 = 0.2e1 * qJD(5) * t701 + t757 * t642 + t754 * t645;
t671 = -t692 * t754 + t714 * t757;
t683 = pkin(5) * t718 - pkin(10) * t702;
t700 = t701 ^ 2;
t636 = -pkin(5) * t700 + pkin(10) * t671 - t683 * t718 + t638;
t760 = sin(qJ(6));
t764 = cos(qJ(6));
t633 = t635 * t764 - t636 * t760;
t673 = t701 * t764 - t702 * t760;
t650 = qJD(6) * t673 + t671 * t760 + t672 * t764;
t674 = t701 * t760 + t702 * t764;
t659 = -mrSges(7,1) * t673 + mrSges(7,2) * t674;
t717 = qJD(6) + t718;
t664 = -mrSges(7,2) * t717 + mrSges(7,3) * t673;
t690 = qJDD(6) + t691;
t631 = m(7) * t633 + mrSges(7,1) * t690 - mrSges(7,3) * t650 - t659 * t674 + t664 * t717;
t634 = t635 * t760 + t636 * t764;
t649 = -qJD(6) * t674 + t671 * t764 - t672 * t760;
t665 = mrSges(7,1) * t717 - mrSges(7,3) * t674;
t632 = m(7) * t634 - mrSges(7,2) * t690 + mrSges(7,3) * t649 + t659 * t673 - t665 * t717;
t623 = t764 * t631 + t760 * t632;
t675 = -mrSges(6,1) * t701 + mrSges(6,2) * t702;
t681 = -mrSges(6,2) * t718 + mrSges(6,3) * t701;
t621 = m(6) * t637 + mrSges(6,1) * t691 - mrSges(6,3) * t672 - t675 * t702 + t681 * t718 + t623;
t682 = mrSges(6,1) * t718 - mrSges(6,3) * t702;
t776 = -t631 * t760 + t764 * t632;
t622 = m(6) * t638 - mrSges(6,2) * t691 + mrSges(6,3) * t671 + t675 * t701 - t682 * t718 + t776;
t777 = -t621 * t754 + t757 * t622;
t618 = m(5) * t647 - mrSges(5,2) * t714 - mrSges(5,3) * t691 - t696 * t718 - t704 * t732 + t777;
t646 = -t761 * t654 + t663 * t794;
t703 = -mrSges(5,2) * t732 - mrSges(5,3) * t718;
t641 = -t714 * pkin(4) - t731 * qJ(5) + t719 * t695 + qJDD(5) - t646;
t639 = -t671 * pkin(5) - t700 * pkin(10) + t702 * t683 + t641;
t771 = m(7) * t639 - t649 * mrSges(7,1) + mrSges(7,2) * t650 - t673 * t664 + t665 * t674;
t768 = -m(6) * t641 + t671 * mrSges(6,1) - mrSges(6,2) * t672 + t701 * t681 - t682 * t702 - t771;
t630 = m(5) * t646 + mrSges(5,1) * t714 - mrSges(5,3) * t692 - t696 * t719 + t703 * t732 + t768;
t778 = t618 * t794 - t630 * t761;
t608 = m(4) * t661 - mrSges(4,2) * t749 + mrSges(4,3) * t715 - t711 * t733 - t721 * t750 + t778;
t720 = -mrSges(4,2) * t750 - mrSges(4,3) * t733;
t619 = t621 * t757 + t622 * t754;
t769 = -m(5) * t653 - t691 * mrSges(5,1) - mrSges(5,2) * t692 - t718 * t703 - t704 * t719 - t619;
t615 = m(4) * t660 + mrSges(4,1) * t749 - mrSges(4,3) * t716 - t711 * t734 + t720 * t750 + t769;
t604 = t755 * t608 + t758 * t615;
t709 = -g(3) * t790 + t774;
t781 = t765 * t786;
t738 = -mrSges(3,2) * t750 + mrSges(3,3) * t781;
t741 = (-mrSges(3,1) * t765 + mrSges(3,2) * t762) * t786;
t602 = m(3) * t709 + mrSges(3,1) * t749 - mrSges(3,3) * t742 + t738 * t750 - t741 * t782 + t604;
t737 = mrSges(3,1) * t750 - mrSges(3,3) * t782;
t779 = t758 * t608 - t615 * t755;
t603 = m(3) * t710 - mrSges(3,2) * t749 + mrSges(3,3) * t743 - t737 * t750 + t741 * t781 + t779;
t611 = t761 * t618 + t630 * t794;
t770 = m(4) * t694 - t715 * mrSges(4,1) + t716 * mrSges(4,2) + t733 * t720 + t734 * t721 + t611;
t610 = m(3) * t725 - t743 * mrSges(3,1) + t742 * mrSges(3,2) + (t737 * t762 - t738 * t765) * t786 + t770;
t590 = t602 * t788 + t603 * t789 - t610 * t756;
t588 = m(2) * t746 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t767 + t590;
t595 = -t602 * t762 + t765 * t603;
t594 = m(2) * t747 - mrSges(2,1) * t767 - qJDD(1) * mrSges(2,2) + t595;
t787 = t766 * t588 + t763 * t594;
t589 = t602 * t790 + t603 * t791 + t759 * t610;
t780 = -t588 * t763 + t766 * t594;
t655 = Ifges(7,5) * t674 + Ifges(7,6) * t673 + Ifges(7,3) * t717;
t657 = Ifges(7,1) * t674 + Ifges(7,4) * t673 + Ifges(7,5) * t717;
t624 = -mrSges(7,1) * t639 + mrSges(7,3) * t634 + Ifges(7,4) * t650 + Ifges(7,2) * t649 + Ifges(7,6) * t690 - t655 * t674 + t657 * t717;
t656 = Ifges(7,4) * t674 + Ifges(7,2) * t673 + Ifges(7,6) * t717;
t625 = mrSges(7,2) * t639 - mrSges(7,3) * t633 + Ifges(7,1) * t650 + Ifges(7,4) * t649 + Ifges(7,5) * t690 + t655 * t673 - t656 * t717;
t666 = Ifges(6,5) * t702 + Ifges(6,6) * t701 + Ifges(6,3) * t718;
t668 = Ifges(6,1) * t702 + Ifges(6,4) * t701 + Ifges(6,5) * t718;
t612 = -mrSges(6,1) * t641 + mrSges(6,3) * t638 + Ifges(6,4) * t672 + Ifges(6,2) * t671 + Ifges(6,6) * t691 - pkin(5) * t771 + pkin(10) * t776 + t764 * t624 + t760 * t625 - t702 * t666 + t718 * t668;
t667 = Ifges(6,4) * t702 + Ifges(6,2) * t701 + Ifges(6,6) * t718;
t613 = mrSges(6,2) * t641 - mrSges(6,3) * t637 + Ifges(6,1) * t672 + Ifges(6,4) * t671 + Ifges(6,5) * t691 - pkin(10) * t623 - t624 * t760 + t625 * t764 + t666 * t701 - t667 * t718;
t685 = Ifges(5,5) * t719 - Ifges(5,6) * t718 + Ifges(5,3) * t732;
t686 = Ifges(5,4) * t719 - Ifges(5,2) * t718 + Ifges(5,6) * t732;
t596 = mrSges(5,2) * t653 - mrSges(5,3) * t646 + Ifges(5,1) * t692 - Ifges(5,4) * t691 + Ifges(5,5) * t714 - qJ(5) * t619 - t612 * t754 + t613 * t757 - t685 * t718 - t686 * t732;
t687 = Ifges(5,1) * t719 - Ifges(5,4) * t718 + Ifges(5,5) * t732;
t605 = Ifges(5,4) * t692 + Ifges(5,6) * t714 - t719 * t685 + t732 * t687 - mrSges(5,1) * t653 + mrSges(5,3) * t647 - Ifges(6,5) * t672 - Ifges(6,6) * t671 - t702 * t667 + t701 * t668 - mrSges(6,1) * t637 + mrSges(6,2) * t638 - Ifges(7,5) * t650 - Ifges(7,6) * t649 - Ifges(7,3) * t690 - t674 * t656 + t673 * t657 - mrSges(7,1) * t633 + mrSges(7,2) * t634 - pkin(5) * t623 - pkin(4) * t619 + (-Ifges(5,2) - Ifges(6,3)) * t691;
t705 = Ifges(4,5) * t734 - Ifges(4,6) * t733 + Ifges(4,3) * t750;
t706 = Ifges(4,4) * t734 - Ifges(4,2) * t733 + Ifges(4,6) * t750;
t586 = mrSges(4,2) * t694 - mrSges(4,3) * t660 + Ifges(4,1) * t716 + Ifges(4,4) * t715 + Ifges(4,5) * t749 - pkin(9) * t611 + t596 * t794 - t761 * t605 - t733 * t705 - t750 * t706;
t707 = Ifges(4,1) * t734 - Ifges(4,4) * t733 + Ifges(4,5) * t750;
t591 = Ifges(4,4) * t716 + Ifges(4,2) * t715 + Ifges(4,6) * t749 - t734 * t705 + t750 * t707 - mrSges(4,1) * t694 + mrSges(4,3) * t661 - Ifges(5,5) * t692 + Ifges(5,6) * t691 - Ifges(5,3) * t714 - t719 * t686 - t718 * t687 - mrSges(5,1) * t646 + mrSges(5,2) * t647 - t754 * t613 - t757 * t612 - pkin(4) * t768 - qJ(5) * t777 - pkin(3) * t611;
t722 = Ifges(3,3) * t750 + (Ifges(3,5) * t762 + Ifges(3,6) * t765) * t786;
t724 = Ifges(3,5) * t750 + (Ifges(3,1) * t762 + Ifges(3,4) * t765) * t786;
t583 = -mrSges(3,1) * t725 + mrSges(3,3) * t710 + Ifges(3,4) * t742 + Ifges(3,2) * t743 + Ifges(3,6) * t749 - pkin(2) * t770 + qJ(3) * t779 + t755 * t586 + t758 * t591 - t722 * t782 + t750 * t724;
t723 = Ifges(3,6) * t750 + (Ifges(3,4) * t762 + Ifges(3,2) * t765) * t786;
t584 = mrSges(3,2) * t725 - mrSges(3,3) * t709 + Ifges(3,1) * t742 + Ifges(3,4) * t743 + Ifges(3,5) * t749 - qJ(3) * t604 + t586 * t758 - t591 * t755 + t722 * t781 - t723 * t750;
t772 = pkin(8) * t595 + t583 * t765 + t584 * t762;
t585 = Ifges(3,5) * t742 + Ifges(3,6) * t743 + mrSges(3,1) * t709 - mrSges(3,2) * t710 + Ifges(4,5) * t716 + Ifges(4,6) * t715 + t734 * t706 + t733 * t707 + mrSges(4,1) * t660 - mrSges(4,2) * t661 + t761 * t596 + t794 * t605 + pkin(3) * t769 + pkin(9) * t778 + pkin(2) * t604 + (Ifges(3,3) + Ifges(4,3)) * t749 + (t723 * t762 - t724 * t765) * t786;
t582 = -mrSges(2,2) * g(3) - mrSges(2,3) * t746 + Ifges(2,5) * qJDD(1) - t767 * Ifges(2,6) - t762 * t583 + t765 * t584 + (-t589 * t756 - t590 * t759) * pkin(8);
t581 = mrSges(2,1) * g(3) + mrSges(2,3) * t747 + t767 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t589 - t756 * t585 + t759 * t772;
t1 = [-m(1) * g(1) + t780; -m(1) * g(2) + t787; (-m(1) - m(2)) * g(3) + t589; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(7) * t787 - t763 * t581 + t766 * t582; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(7) * t780 + t766 * t581 + t763 * t582; -mrSges(1,1) * g(2) + mrSges(2,1) * t746 + mrSges(1,2) * g(1) - mrSges(2,2) * t747 + Ifges(2,3) * qJDD(1) + pkin(1) * t590 + t759 * t585 + t756 * t772;];
tauB  = t1;
