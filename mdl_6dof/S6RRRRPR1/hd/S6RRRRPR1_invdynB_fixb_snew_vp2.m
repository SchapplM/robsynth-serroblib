% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 19:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:27:50
% EndTime: 2019-05-07 19:28:23
% DurationCPUTime: 32.20s
% Computational Cost: add. (512345->388), mult. (1185958->490), div. (0->0), fcn. (906230->12), ass. (0->151)
t762 = 2 * qJD(5);
t742 = qJD(1) ^ 2;
t761 = pkin(2) * t742;
t736 = sin(qJ(1));
t741 = cos(qJ(1));
t722 = -g(1) * t741 - g(2) * t736;
t711 = -pkin(1) * t742 + qJDD(1) * pkin(7) + t722;
t735 = sin(qJ(2));
t760 = t735 * t711;
t740 = cos(qJ(2));
t756 = qJD(1) * qJD(2);
t716 = qJDD(1) * t735 + t740 * t756;
t680 = qJDD(2) * pkin(2) - t716 * pkin(8) - t760 + (pkin(8) * t756 + t735 * t761 - g(3)) * t740;
t699 = -g(3) * t735 + t740 * t711;
t717 = qJDD(1) * t740 - t735 * t756;
t758 = qJD(1) * t735;
t720 = qJD(2) * pkin(2) - pkin(8) * t758;
t729 = t740 ^ 2;
t681 = pkin(8) * t717 - qJD(2) * t720 - t729 * t761 + t699;
t734 = sin(qJ(3));
t739 = cos(qJ(3));
t660 = t739 * t680 - t734 * t681;
t708 = (-t734 * t735 + t739 * t740) * qJD(1);
t684 = qJD(3) * t708 + t716 * t739 + t717 * t734;
t709 = (t734 * t740 + t735 * t739) * qJD(1);
t727 = qJDD(2) + qJDD(3);
t728 = qJD(2) + qJD(3);
t641 = (t708 * t728 - t684) * pkin(9) + (t708 * t709 + t727) * pkin(3) + t660;
t661 = t734 * t680 + t739 * t681;
t683 = -qJD(3) * t709 - t716 * t734 + t717 * t739;
t702 = pkin(3) * t728 - pkin(9) * t709;
t704 = t708 ^ 2;
t643 = -pkin(3) * t704 + pkin(9) * t683 - t702 * t728 + t661;
t733 = sin(qJ(4));
t738 = cos(qJ(4));
t623 = t738 * t641 - t733 * t643;
t695 = t708 * t738 - t709 * t733;
t659 = qJD(4) * t695 + t683 * t733 + t684 * t738;
t696 = t708 * t733 + t709 * t738;
t724 = qJDD(4) + t727;
t725 = qJD(4) + t728;
t620 = (t695 * t725 - t659) * qJ(5) + (t695 * t696 + t724) * pkin(4) + t623;
t624 = t733 * t641 + t738 * t643;
t658 = -qJD(4) * t696 + t683 * t738 - t684 * t733;
t687 = pkin(4) * t725 - qJ(5) * t696;
t694 = t695 ^ 2;
t622 = -pkin(4) * t694 + qJ(5) * t658 - t687 * t725 + t624;
t730 = sin(pkin(11));
t731 = cos(pkin(11));
t674 = t695 * t731 - t696 * t730;
t617 = t730 * t620 + t731 * t622 + t674 * t762;
t633 = t658 * t731 - t659 * t730;
t675 = t695 * t730 + t696 * t731;
t652 = -mrSges(6,1) * t674 + mrSges(6,2) * t675;
t665 = mrSges(6,1) * t725 - mrSges(6,3) * t675;
t653 = -pkin(5) * t674 - pkin(10) * t675;
t723 = t725 ^ 2;
t615 = -pkin(5) * t723 + pkin(10) * t724 + t653 * t674 + t617;
t721 = t736 * g(1) - t741 * g(2);
t748 = -qJDD(1) * pkin(1) - t721;
t685 = -t717 * pkin(2) + t720 * t758 + (-pkin(8) * t729 - pkin(7)) * t742 + t748;
t655 = -t683 * pkin(3) - t704 * pkin(9) + t709 * t702 + t685;
t626 = -t658 * pkin(4) - t694 * qJ(5) + t696 * t687 + qJDD(5) + t655;
t634 = t658 * t730 + t659 * t731;
t618 = t626 + (-t674 * t725 - t634) * pkin(10) + (t675 * t725 - t633) * pkin(5);
t732 = sin(qJ(6));
t737 = cos(qJ(6));
t612 = -t615 * t732 + t618 * t737;
t662 = -t675 * t732 + t725 * t737;
t629 = qJD(6) * t662 + t634 * t737 + t724 * t732;
t632 = qJDD(6) - t633;
t663 = t675 * t737 + t725 * t732;
t644 = -mrSges(7,1) * t662 + mrSges(7,2) * t663;
t673 = qJD(6) - t674;
t645 = -mrSges(7,2) * t673 + mrSges(7,3) * t662;
t610 = m(7) * t612 + mrSges(7,1) * t632 - mrSges(7,3) * t629 - t644 * t663 + t645 * t673;
t613 = t615 * t737 + t618 * t732;
t628 = -qJD(6) * t663 - t634 * t732 + t724 * t737;
t646 = mrSges(7,1) * t673 - mrSges(7,3) * t663;
t611 = m(7) * t613 - mrSges(7,2) * t632 + mrSges(7,3) * t628 + t644 * t662 - t646 * t673;
t750 = -t610 * t732 + t737 * t611;
t601 = m(6) * t617 - mrSges(6,2) * t724 + mrSges(6,3) * t633 + t652 * t674 - t665 * t725 + t750;
t749 = -t731 * t620 + t730 * t622;
t616 = -0.2e1 * qJD(5) * t675 - t749;
t664 = -mrSges(6,2) * t725 + mrSges(6,3) * t674;
t614 = -t724 * pkin(5) - t723 * pkin(10) + (t762 + t653) * t675 + t749;
t746 = -m(7) * t614 + t628 * mrSges(7,1) - mrSges(7,2) * t629 + t662 * t645 - t646 * t663;
t606 = m(6) * t616 + mrSges(6,1) * t724 - mrSges(6,3) * t634 - t652 * t675 + t664 * t725 + t746;
t596 = t730 * t601 + t731 * t606;
t676 = -mrSges(5,1) * t695 + mrSges(5,2) * t696;
t686 = -mrSges(5,2) * t725 + mrSges(5,3) * t695;
t594 = m(5) * t623 + mrSges(5,1) * t724 - mrSges(5,3) * t659 - t676 * t696 + t686 * t725 + t596;
t688 = mrSges(5,1) * t725 - mrSges(5,3) * t696;
t751 = t731 * t601 - t606 * t730;
t595 = m(5) * t624 - mrSges(5,2) * t724 + mrSges(5,3) * t658 + t676 * t695 - t688 * t725 + t751;
t588 = t738 * t594 + t733 * t595;
t697 = -mrSges(4,1) * t708 + mrSges(4,2) * t709;
t700 = -mrSges(4,2) * t728 + mrSges(4,3) * t708;
t586 = m(4) * t660 + mrSges(4,1) * t727 - mrSges(4,3) * t684 - t697 * t709 + t700 * t728 + t588;
t701 = mrSges(4,1) * t728 - mrSges(4,3) * t709;
t752 = -t594 * t733 + t738 * t595;
t587 = m(4) * t661 - mrSges(4,2) * t727 + mrSges(4,3) * t683 + t697 * t708 - t701 * t728 + t752;
t581 = t739 * t586 + t734 * t587;
t698 = -t740 * g(3) - t760;
t715 = (-mrSges(3,1) * t740 + mrSges(3,2) * t735) * qJD(1);
t757 = qJD(1) * t740;
t719 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t757;
t579 = m(3) * t698 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t716 + qJD(2) * t719 - t715 * t758 + t581;
t718 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t758;
t753 = -t586 * t734 + t739 * t587;
t580 = m(3) * t699 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t717 - qJD(2) * t718 + t715 * t757 + t753;
t754 = -t579 * t735 + t740 * t580;
t573 = m(2) * t722 - mrSges(2,1) * t742 - qJDD(1) * mrSges(2,2) + t754;
t710 = -t742 * pkin(7) + t748;
t602 = t737 * t610 + t732 * t611;
t747 = m(6) * t626 - t633 * mrSges(6,1) + t634 * mrSges(6,2) - t674 * t664 + t675 * t665 + t602;
t745 = m(5) * t655 - t658 * mrSges(5,1) + t659 * mrSges(5,2) - t695 * t686 + t696 * t688 + t747;
t744 = m(4) * t685 - t683 * mrSges(4,1) + t684 * mrSges(4,2) - t708 * t700 + t709 * t701 + t745;
t743 = -m(3) * t710 + t717 * mrSges(3,1) - t716 * mrSges(3,2) - t718 * t758 + t719 * t757 - t744;
t598 = m(2) * t721 + qJDD(1) * mrSges(2,1) - t742 * mrSges(2,2) + t743;
t759 = t736 * t573 + t741 * t598;
t574 = t740 * t579 + t735 * t580;
t755 = t741 * t573 - t598 * t736;
t707 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t735 + Ifges(3,4) * t740) * qJD(1);
t706 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t735 + Ifges(3,2) * t740) * qJD(1);
t705 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t735 + Ifges(3,6) * t740) * qJD(1);
t691 = Ifges(4,1) * t709 + Ifges(4,4) * t708 + Ifges(4,5) * t728;
t690 = Ifges(4,4) * t709 + Ifges(4,2) * t708 + Ifges(4,6) * t728;
t689 = Ifges(4,5) * t709 + Ifges(4,6) * t708 + Ifges(4,3) * t728;
t668 = Ifges(5,1) * t696 + Ifges(5,4) * t695 + Ifges(5,5) * t725;
t667 = Ifges(5,4) * t696 + Ifges(5,2) * t695 + Ifges(5,6) * t725;
t666 = Ifges(5,5) * t696 + Ifges(5,6) * t695 + Ifges(5,3) * t725;
t649 = Ifges(6,1) * t675 + Ifges(6,4) * t674 + Ifges(6,5) * t725;
t648 = Ifges(6,4) * t675 + Ifges(6,2) * t674 + Ifges(6,6) * t725;
t647 = Ifges(6,5) * t675 + Ifges(6,6) * t674 + Ifges(6,3) * t725;
t637 = Ifges(7,1) * t663 + Ifges(7,4) * t662 + Ifges(7,5) * t673;
t636 = Ifges(7,4) * t663 + Ifges(7,2) * t662 + Ifges(7,6) * t673;
t635 = Ifges(7,5) * t663 + Ifges(7,6) * t662 + Ifges(7,3) * t673;
t604 = mrSges(7,2) * t614 - mrSges(7,3) * t612 + Ifges(7,1) * t629 + Ifges(7,4) * t628 + Ifges(7,5) * t632 + t635 * t662 - t636 * t673;
t603 = -mrSges(7,1) * t614 + mrSges(7,3) * t613 + Ifges(7,4) * t629 + Ifges(7,2) * t628 + Ifges(7,6) * t632 - t635 * t663 + t637 * t673;
t590 = -mrSges(6,1) * t626 - mrSges(7,1) * t612 + mrSges(7,2) * t613 + mrSges(6,3) * t617 + Ifges(6,4) * t634 - Ifges(7,5) * t629 + Ifges(6,2) * t633 + Ifges(6,6) * t724 - Ifges(7,6) * t628 - Ifges(7,3) * t632 - pkin(5) * t602 - t636 * t663 + t637 * t662 - t647 * t675 + t649 * t725;
t589 = mrSges(6,2) * t626 - mrSges(6,3) * t616 + Ifges(6,1) * t634 + Ifges(6,4) * t633 + Ifges(6,5) * t724 - pkin(10) * t602 - t603 * t732 + t604 * t737 + t647 * t674 - t648 * t725;
t582 = mrSges(5,2) * t655 - mrSges(5,3) * t623 + Ifges(5,1) * t659 + Ifges(5,4) * t658 + Ifges(5,5) * t724 - qJ(5) * t596 + t589 * t731 - t590 * t730 + t666 * t695 - t667 * t725;
t575 = -mrSges(5,1) * t655 + mrSges(5,3) * t624 + Ifges(5,4) * t659 + Ifges(5,2) * t658 + Ifges(5,6) * t724 - pkin(4) * t747 + qJ(5) * t751 + t730 * t589 + t731 * t590 - t696 * t666 + t725 * t668;
t570 = mrSges(4,2) * t685 - mrSges(4,3) * t660 + Ifges(4,1) * t684 + Ifges(4,4) * t683 + Ifges(4,5) * t727 - pkin(9) * t588 - t575 * t733 + t582 * t738 + t689 * t708 - t690 * t728;
t569 = -mrSges(4,1) * t685 + mrSges(4,3) * t661 + Ifges(4,4) * t684 + Ifges(4,2) * t683 + Ifges(4,6) * t727 - pkin(3) * t745 + pkin(9) * t752 + t738 * t575 + t733 * t582 - t709 * t689 + t728 * t691;
t568 = -Ifges(3,3) * qJDD(2) - pkin(4) * t596 - pkin(10) * t750 + mrSges(2,1) * g(3) + Ifges(2,6) * qJDD(1) - pkin(3) * t588 - pkin(1) * t574 - pkin(2) * t581 + (-Ifges(5,3) - Ifges(6,3)) * t724 + (-t706 * t735 + t707 * t740) * qJD(1) - mrSges(6,1) * t616 + mrSges(6,2) * t617 - pkin(5) * t746 - mrSges(5,1) * t623 + mrSges(5,2) * t624 - Ifges(6,6) * t633 - Ifges(6,5) * t634 - Ifges(5,6) * t658 - Ifges(5,5) * t659 - mrSges(4,1) * t660 + mrSges(4,2) * t661 + t674 * t649 - t675 * t648 - Ifges(4,6) * t683 - Ifges(4,5) * t684 + t695 * t668 - t696 * t667 - mrSges(3,1) * t698 + mrSges(3,2) * t699 + t708 * t691 - t709 * t690 - Ifges(3,5) * t716 - Ifges(3,6) * t717 + mrSges(2,3) * t722 - Ifges(4,3) * t727 - t732 * t604 - t737 * t603 + t742 * Ifges(2,5);
t567 = mrSges(3,2) * t710 - mrSges(3,3) * t698 + Ifges(3,1) * t716 + Ifges(3,4) * t717 + Ifges(3,5) * qJDD(2) - pkin(8) * t581 - qJD(2) * t706 - t569 * t734 + t570 * t739 + t705 * t757;
t566 = -mrSges(3,1) * t710 + mrSges(3,3) * t699 + Ifges(3,4) * t716 + Ifges(3,2) * t717 + Ifges(3,6) * qJDD(2) - pkin(2) * t744 + pkin(8) * t753 + qJD(2) * t707 + t739 * t569 + t734 * t570 - t705 * t758;
t565 = -mrSges(2,2) * g(3) - mrSges(2,3) * t721 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t742 - pkin(7) * t574 - t566 * t735 + t567 * t740;
t1 = [-m(1) * g(1) + t755; -m(1) * g(2) + t759; (-m(1) - m(2)) * g(3) + t574; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t759 + t741 * t565 - t736 * t568; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t755 + t736 * t565 + t741 * t568; -mrSges(1,1) * g(2) + mrSges(2,1) * t721 + mrSges(1,2) * g(1) - mrSges(2,2) * t722 + Ifges(2,3) * qJDD(1) + pkin(1) * t743 + pkin(7) * t754 + t740 * t566 + t735 * t567;];
tauB  = t1;
