% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:30:24
% EndTime: 2019-05-06 17:30:39
% DurationCPUTime: 12.14s
% Computational Cost: add. (176465->365), mult. (402438->447), div. (0->0), fcn. (289140->10), ass. (0->139)
t770 = Ifges(6,1) + Ifges(7,1);
t766 = Ifges(6,4) + Ifges(7,4);
t765 = Ifges(6,5) + Ifges(7,5);
t769 = Ifges(6,2) + Ifges(7,2);
t764 = -Ifges(6,6) - Ifges(7,6);
t768 = -Ifges(6,3) - Ifges(7,3);
t733 = sin(qJ(2));
t737 = cos(qJ(2));
t755 = qJD(1) * qJD(2);
t719 = qJDD(1) * t733 + t737 * t755;
t734 = sin(qJ(1));
t738 = cos(qJ(1));
t725 = -g(1) * t738 - g(2) * t734;
t740 = qJD(1) ^ 2;
t714 = -pkin(1) * t740 + qJDD(1) * pkin(7) + t725;
t763 = t714 * t733;
t767 = pkin(2) * t740;
t675 = qJDD(2) * pkin(2) - qJ(3) * t719 - t763 + (qJ(3) * t755 + t733 * t767 - g(3)) * t737;
t699 = -g(3) * t733 + t737 * t714;
t720 = qJDD(1) * t737 - t733 * t755;
t758 = qJD(1) * t733;
t721 = qJD(2) * pkin(2) - qJ(3) * t758;
t728 = t737 ^ 2;
t676 = qJ(3) * t720 - qJD(2) * t721 - t728 * t767 + t699;
t729 = sin(pkin(10));
t730 = cos(pkin(10));
t708 = (t729 * t737 + t730 * t733) * qJD(1);
t651 = -0.2e1 * qJD(3) * t708 + t675 * t730 - t729 * t676;
t757 = qJD(1) * t737;
t707 = -t729 * t758 + t730 * t757;
t652 = 0.2e1 * qJD(3) * t707 + t729 * t675 + t730 * t676;
t687 = -mrSges(4,1) * t707 + mrSges(4,2) * t708;
t693 = -t729 * t719 + t720 * t730;
t701 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t708;
t689 = -pkin(3) * t707 - pkin(8) * t708;
t739 = qJD(2) ^ 2;
t637 = -pkin(3) * t739 + qJDD(2) * pkin(8) + t689 * t707 + t652;
t724 = g(1) * t734 - t738 * g(2);
t745 = -qJDD(1) * pkin(1) - t724;
t679 = -pkin(2) * t720 + qJDD(3) + t721 * t758 + (-qJ(3) * t728 - pkin(7)) * t740 + t745;
t694 = t719 * t730 + t720 * t729;
t641 = (-qJD(2) * t707 - t694) * pkin(8) + (qJD(2) * t708 - t693) * pkin(3) + t679;
t732 = sin(qJ(4));
t736 = cos(qJ(4));
t627 = -t637 * t732 + t736 * t641;
t696 = qJD(2) * t736 - t708 * t732;
t668 = qJD(4) * t696 + qJDD(2) * t732 + t694 * t736;
t692 = qJDD(4) - t693;
t697 = qJD(2) * t732 + t708 * t736;
t706 = qJD(4) - t707;
t623 = (t696 * t706 - t668) * pkin(9) + (t696 * t697 + t692) * pkin(4) + t627;
t628 = t736 * t637 + t732 * t641;
t667 = -qJD(4) * t697 + qJDD(2) * t736 - t694 * t732;
t682 = pkin(4) * t706 - pkin(9) * t697;
t695 = t696 ^ 2;
t625 = -pkin(4) * t695 + pkin(9) * t667 - t682 * t706 + t628;
t731 = sin(qJ(5));
t735 = cos(qJ(5));
t617 = t735 * t623 - t625 * t731;
t673 = t696 * t735 - t697 * t731;
t634 = qJD(5) * t673 + t667 * t731 + t668 * t735;
t674 = t696 * t731 + t697 * t735;
t653 = -mrSges(7,1) * t673 + mrSges(7,2) * t674;
t654 = -mrSges(6,1) * t673 + mrSges(6,2) * t674;
t702 = qJD(5) + t706;
t657 = -mrSges(6,2) * t702 + mrSges(6,3) * t673;
t690 = qJDD(5) + t692;
t614 = -0.2e1 * qJD(6) * t674 + (t673 * t702 - t634) * qJ(6) + (t673 * t674 + t690) * pkin(5) + t617;
t656 = -mrSges(7,2) * t702 + mrSges(7,3) * t673;
t754 = m(7) * t614 + t690 * mrSges(7,1) + t702 * t656;
t606 = m(6) * t617 + mrSges(6,1) * t690 + t657 * t702 + (-t653 - t654) * t674 + (-mrSges(6,3) - mrSges(7,3)) * t634 + t754;
t618 = t731 * t623 + t735 * t625;
t633 = -qJD(5) * t674 + t667 * t735 - t668 * t731;
t659 = mrSges(7,1) * t702 - mrSges(7,3) * t674;
t660 = mrSges(6,1) * t702 - mrSges(6,3) * t674;
t658 = pkin(5) * t702 - qJ(6) * t674;
t672 = t673 ^ 2;
t616 = -pkin(5) * t672 + qJ(6) * t633 + 0.2e1 * qJD(6) * t673 - t658 * t702 + t618;
t753 = m(7) * t616 + t633 * mrSges(7,3) + t673 * t653;
t609 = m(6) * t618 + mrSges(6,3) * t633 + t654 * t673 + (-t659 - t660) * t702 + (-mrSges(6,2) - mrSges(7,2)) * t690 + t753;
t604 = t735 * t606 + t731 * t609;
t677 = -mrSges(5,1) * t696 + mrSges(5,2) * t697;
t680 = -mrSges(5,2) * t706 + mrSges(5,3) * t696;
t601 = m(5) * t627 + mrSges(5,1) * t692 - mrSges(5,3) * t668 - t677 * t697 + t680 * t706 + t604;
t681 = mrSges(5,1) * t706 - mrSges(5,3) * t697;
t748 = -t606 * t731 + t735 * t609;
t602 = m(5) * t628 - mrSges(5,2) * t692 + mrSges(5,3) * t667 + t677 * t696 - t681 * t706 + t748;
t749 = -t601 * t732 + t736 * t602;
t595 = m(4) * t652 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t693 - qJD(2) * t701 + t687 * t707 + t749;
t700 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t707;
t636 = -qJDD(2) * pkin(3) - pkin(8) * t739 + t708 * t689 - t651;
t626 = -pkin(4) * t667 - pkin(9) * t695 + t697 * t682 + t636;
t620 = -pkin(5) * t633 - qJ(6) * t672 + t658 * t674 + qJDD(6) + t626;
t746 = m(7) * t620 - t633 * mrSges(7,1) + t634 * mrSges(7,2) - t673 * t656 + t674 * t659;
t743 = m(6) * t626 - t633 * mrSges(6,1) + mrSges(6,2) * t634 - t673 * t657 + t660 * t674 + t746;
t741 = -m(5) * t636 + t667 * mrSges(5,1) - mrSges(5,2) * t668 + t696 * t680 - t681 * t697 - t743;
t611 = m(4) * t651 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t694 + qJD(2) * t700 - t687 * t708 + t741;
t590 = t729 * t595 + t730 * t611;
t698 = -g(3) * t737 - t763;
t718 = (-mrSges(3,1) * t737 + mrSges(3,2) * t733) * qJD(1);
t723 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t757;
t588 = m(3) * t698 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t719 + qJD(2) * t723 - t718 * t758 + t590;
t722 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t758;
t750 = t730 * t595 - t611 * t729;
t589 = m(3) * t699 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t720 - qJD(2) * t722 + t718 * t757 + t750;
t751 = -t588 * t733 + t737 * t589;
t580 = m(2) * t725 - mrSges(2,1) * t740 - qJDD(1) * mrSges(2,2) + t751;
t713 = -pkin(7) * t740 + t745;
t596 = t736 * t601 + t732 * t602;
t744 = m(4) * t679 - t693 * mrSges(4,1) + mrSges(4,2) * t694 - t707 * t700 + t701 * t708 + t596;
t742 = -m(3) * t713 + t720 * mrSges(3,1) - mrSges(3,2) * t719 - t722 * t758 + t723 * t757 - t744;
t592 = m(2) * t724 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t740 + t742;
t762 = t734 * t580 + t738 * t592;
t581 = t737 * t588 + t733 * t589;
t761 = t764 * t673 - t765 * t674 + t768 * t702;
t760 = -t769 * t673 - t766 * t674 + t764 * t702;
t759 = t766 * t673 + t770 * t674 + t765 * t702;
t752 = t738 * t580 - t592 * t734;
t711 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t733 + Ifges(3,4) * t737) * qJD(1);
t710 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t733 + Ifges(3,2) * t737) * qJD(1);
t709 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t733 + Ifges(3,6) * t737) * qJD(1);
t685 = Ifges(4,1) * t708 + Ifges(4,4) * t707 + Ifges(4,5) * qJD(2);
t684 = Ifges(4,4) * t708 + Ifges(4,2) * t707 + Ifges(4,6) * qJD(2);
t683 = Ifges(4,5) * t708 + Ifges(4,6) * t707 + Ifges(4,3) * qJD(2);
t663 = Ifges(5,1) * t697 + Ifges(5,4) * t696 + Ifges(5,5) * t706;
t662 = Ifges(5,4) * t697 + Ifges(5,2) * t696 + Ifges(5,6) * t706;
t661 = Ifges(5,5) * t697 + Ifges(5,6) * t696 + Ifges(5,3) * t706;
t612 = -mrSges(7,3) * t634 - t653 * t674 + t754;
t603 = mrSges(6,2) * t626 + mrSges(7,2) * t620 - mrSges(6,3) * t617 - mrSges(7,3) * t614 - qJ(6) * t612 + t766 * t633 + t770 * t634 - t761 * t673 + t765 * t690 + t760 * t702;
t597 = -mrSges(6,1) * t626 + mrSges(6,3) * t618 - mrSges(7,1) * t620 + mrSges(7,3) * t616 - pkin(5) * t746 + qJ(6) * t753 + (-qJ(6) * t659 + t759) * t702 + (-mrSges(7,2) * qJ(6) - t764) * t690 + t761 * t674 + t766 * t634 + t769 * t633;
t584 = mrSges(5,2) * t636 - mrSges(5,3) * t627 + Ifges(5,1) * t668 + Ifges(5,4) * t667 + Ifges(5,5) * t692 - pkin(9) * t604 - t597 * t731 + t603 * t735 + t661 * t696 - t662 * t706;
t583 = -mrSges(5,1) * t636 + mrSges(5,3) * t628 + Ifges(5,4) * t668 + Ifges(5,2) * t667 + Ifges(5,6) * t692 - pkin(4) * t743 + pkin(9) * t748 + t735 * t597 + t731 * t603 - t697 * t661 + t706 * t663;
t582 = t759 * t673 + t760 * t674 - t708 * t683 + t696 * t663 - t697 * t662 - Ifges(5,3) * t692 + Ifges(4,2) * t693 + Ifges(4,4) * t694 + qJD(2) * t685 - mrSges(4,1) * t679 - Ifges(5,6) * t667 - Ifges(5,5) * t668 + mrSges(4,3) * t652 - mrSges(5,1) * t627 + mrSges(5,2) * t628 + mrSges(6,2) * t618 - mrSges(6,1) * t617 - mrSges(7,1) * t614 + mrSges(7,2) * t616 - pkin(5) * t612 - pkin(4) * t604 + Ifges(4,6) * qJDD(2) - pkin(3) * t596 - t765 * t634 + t768 * t690 + t764 * t633;
t577 = mrSges(4,2) * t679 - mrSges(4,3) * t651 + Ifges(4,1) * t694 + Ifges(4,4) * t693 + Ifges(4,5) * qJDD(2) - pkin(8) * t596 - qJD(2) * t684 - t583 * t732 + t584 * t736 + t683 * t707;
t576 = mrSges(3,2) * t713 - mrSges(3,3) * t698 + Ifges(3,1) * t719 + Ifges(3,4) * t720 + Ifges(3,5) * qJDD(2) - qJ(3) * t590 - qJD(2) * t710 + t577 * t730 - t582 * t729 + t709 * t757;
t575 = -pkin(1) * t581 + mrSges(2,3) * t725 - pkin(2) * t590 - Ifges(3,5) * t719 - Ifges(3,6) * t720 - mrSges(3,1) * t698 + mrSges(3,2) * t699 - t732 * t584 - t736 * t583 - pkin(3) * t741 - pkin(8) * t749 - Ifges(4,5) * t694 - Ifges(4,6) * t693 - mrSges(4,1) * t651 + mrSges(4,2) * t652 + mrSges(2,1) * g(3) - t708 * t684 + t707 * t685 + t740 * Ifges(2,5) + Ifges(2,6) * qJDD(1) + (-Ifges(3,3) - Ifges(4,3)) * qJDD(2) + (-t710 * t733 + t711 * t737) * qJD(1);
t574 = -mrSges(3,1) * t713 + mrSges(3,3) * t699 + Ifges(3,4) * t719 + Ifges(3,2) * t720 + Ifges(3,6) * qJDD(2) - pkin(2) * t744 + qJ(3) * t750 + qJD(2) * t711 + t729 * t577 + t730 * t582 - t709 * t758;
t573 = -mrSges(2,2) * g(3) - mrSges(2,3) * t724 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t740 - pkin(7) * t581 - t574 * t733 + t576 * t737;
t1 = [-m(1) * g(1) + t752; -m(1) * g(2) + t762; (-m(1) - m(2)) * g(3) + t581; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t762 + t738 * t573 - t734 * t575; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t752 + t734 * t573 + t738 * t575; -mrSges(1,1) * g(2) + mrSges(2,1) * t724 + mrSges(1,2) * g(1) - mrSges(2,2) * t725 + Ifges(2,3) * qJDD(1) + pkin(1) * t742 + pkin(7) * t751 + t737 * t574 + t733 * t576;];
tauB  = t1;
