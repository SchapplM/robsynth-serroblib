% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 03:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:31:03
% EndTime: 2019-05-06 03:31:26
% DurationCPUTime: 23.89s
% Computational Cost: add. (386130->365), mult. (923400->456), div. (0->0), fcn. (737726->12), ass. (0->150)
t738 = qJD(1) ^ 2;
t727 = cos(pkin(11));
t765 = pkin(2) * t727;
t726 = sin(pkin(11));
t764 = mrSges(3,2) * t726;
t724 = t727 ^ 2;
t763 = t724 * t738;
t732 = sin(qJ(1));
t737 = cos(qJ(1));
t714 = -t737 * g(1) - t732 * g(2);
t710 = -t738 * pkin(1) + qJDD(1) * qJ(2) + t714;
t759 = qJD(1) * qJD(2);
t757 = -t727 * g(3) - 0.2e1 * t726 * t759;
t681 = (-pkin(7) * qJDD(1) + t738 * t765 - t710) * t726 + t757;
t700 = -t726 * g(3) + (t710 + 0.2e1 * t759) * t727;
t758 = qJDD(1) * t727;
t682 = -pkin(2) * t763 + pkin(7) * t758 + t700;
t731 = sin(qJ(3));
t736 = cos(qJ(3));
t661 = t736 * t681 - t731 * t682;
t746 = t726 * t736 + t727 * t731;
t745 = -t726 * t731 + t727 * t736;
t708 = t745 * qJD(1);
t760 = t708 * qJD(3);
t698 = t746 * qJDD(1) + t760;
t709 = t746 * qJD(1);
t636 = (-t698 + t760) * pkin(8) + (t708 * t709 + qJDD(3)) * pkin(3) + t661;
t662 = t731 * t681 + t736 * t682;
t697 = -t709 * qJD(3) + t745 * qJDD(1);
t703 = qJD(3) * pkin(3) - t709 * pkin(8);
t707 = t708 ^ 2;
t642 = -t707 * pkin(3) + t697 * pkin(8) - qJD(3) * t703 + t662;
t730 = sin(qJ(4));
t735 = cos(qJ(4));
t629 = t730 * t636 + t735 * t642;
t690 = t730 * t708 + t735 * t709;
t656 = -t690 * qJD(4) + t735 * t697 - t730 * t698;
t689 = t735 * t708 - t730 * t709;
t671 = -t689 * mrSges(5,1) + t690 * mrSges(5,2);
t725 = qJD(3) + qJD(4);
t679 = t725 * mrSges(5,1) - t690 * mrSges(5,3);
t722 = qJDD(3) + qJDD(4);
t672 = -t689 * pkin(4) - t690 * pkin(9);
t721 = t725 ^ 2;
t621 = -t721 * pkin(4) + t722 * pkin(9) + t689 * t672 + t629;
t723 = t726 ^ 2;
t713 = t732 * g(1) - t737 * g(2);
t750 = qJDD(2) - t713;
t696 = (-pkin(1) - t765) * qJDD(1) + (-qJ(2) + (-t723 - t724) * pkin(7)) * t738 + t750;
t650 = -t697 * pkin(3) - t707 * pkin(8) + t709 * t703 + t696;
t657 = t689 * qJD(4) + t730 * t697 + t735 * t698;
t627 = (-t689 * t725 - t657) * pkin(9) + (t690 * t725 - t656) * pkin(4) + t650;
t729 = sin(qJ(5));
t734 = cos(qJ(5));
t616 = -t729 * t621 + t734 * t627;
t674 = -t729 * t690 + t734 * t725;
t639 = t674 * qJD(5) + t734 * t657 + t729 * t722;
t655 = qJDD(5) - t656;
t675 = t734 * t690 + t729 * t725;
t685 = qJD(5) - t689;
t614 = (t674 * t685 - t639) * pkin(10) + (t674 * t675 + t655) * pkin(5) + t616;
t617 = t734 * t621 + t729 * t627;
t638 = -t675 * qJD(5) - t729 * t657 + t734 * t722;
t665 = t685 * pkin(5) - t675 * pkin(10);
t673 = t674 ^ 2;
t615 = -t673 * pkin(5) + t638 * pkin(10) - t685 * t665 + t617;
t728 = sin(qJ(6));
t733 = cos(qJ(6));
t612 = t733 * t614 - t728 * t615;
t658 = t733 * t674 - t728 * t675;
t624 = t658 * qJD(6) + t728 * t638 + t733 * t639;
t659 = t728 * t674 + t733 * t675;
t634 = -t658 * mrSges(7,1) + t659 * mrSges(7,2);
t683 = qJD(6) + t685;
t643 = -t683 * mrSges(7,2) + t658 * mrSges(7,3);
t652 = qJDD(6) + t655;
t610 = m(7) * t612 + t652 * mrSges(7,1) - t624 * mrSges(7,3) - t659 * t634 + t683 * t643;
t613 = t728 * t614 + t733 * t615;
t623 = -t659 * qJD(6) + t733 * t638 - t728 * t639;
t644 = t683 * mrSges(7,1) - t659 * mrSges(7,3);
t611 = m(7) * t613 - t652 * mrSges(7,2) + t623 * mrSges(7,3) + t658 * t634 - t683 * t644;
t602 = t733 * t610 + t728 * t611;
t660 = -t674 * mrSges(6,1) + t675 * mrSges(6,2);
t663 = -t685 * mrSges(6,2) + t674 * mrSges(6,3);
t600 = m(6) * t616 + t655 * mrSges(6,1) - t639 * mrSges(6,3) - t675 * t660 + t685 * t663 + t602;
t664 = t685 * mrSges(6,1) - t675 * mrSges(6,3);
t751 = -t728 * t610 + t733 * t611;
t601 = m(6) * t617 - t655 * mrSges(6,2) + t638 * mrSges(6,3) + t674 * t660 - t685 * t664 + t751;
t752 = -t729 * t600 + t734 * t601;
t595 = m(5) * t629 - t722 * mrSges(5,2) + t656 * mrSges(5,3) + t689 * t671 - t725 * t679 + t752;
t628 = t735 * t636 - t730 * t642;
t678 = -t725 * mrSges(5,2) + t689 * mrSges(5,3);
t620 = -t722 * pkin(4) - t721 * pkin(9) + t690 * t672 - t628;
t618 = -t638 * pkin(5) - t673 * pkin(10) + t675 * t665 + t620;
t742 = m(7) * t618 - t623 * mrSges(7,1) + t624 * mrSges(7,2) - t658 * t643 + t659 * t644;
t740 = -m(6) * t620 + t638 * mrSges(6,1) - t639 * mrSges(6,2) + t674 * t663 - t675 * t664 - t742;
t606 = m(5) * t628 + t722 * mrSges(5,1) - t657 * mrSges(5,3) - t690 * t671 + t725 * t678 + t740;
t588 = t730 * t595 + t735 * t606;
t694 = -t708 * mrSges(4,1) + t709 * mrSges(4,2);
t701 = -qJD(3) * mrSges(4,2) + t708 * mrSges(4,3);
t586 = m(4) * t661 + qJDD(3) * mrSges(4,1) - t698 * mrSges(4,3) + qJD(3) * t701 - t709 * t694 + t588;
t702 = qJD(3) * mrSges(4,1) - t709 * mrSges(4,3);
t753 = t735 * t595 - t730 * t606;
t587 = m(4) * t662 - qJDD(3) * mrSges(4,2) + t697 * mrSges(4,3) - qJD(3) * t702 + t708 * t694 + t753;
t581 = t736 * t586 + t731 * t587;
t699 = -t726 * t710 + t757;
t744 = mrSges(3,3) * qJDD(1) + t738 * (-mrSges(3,1) * t727 + t764);
t579 = m(3) * t699 - t744 * t726 + t581;
t754 = -t731 * t586 + t736 * t587;
t580 = m(3) * t700 + t744 * t727 + t754;
t755 = -t726 * t579 + t727 * t580;
t573 = m(2) * t714 - t738 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t755;
t706 = -qJDD(1) * pkin(1) - t738 * qJ(2) + t750;
t596 = t734 * t600 + t729 * t601;
t743 = m(5) * t650 - t656 * mrSges(5,1) + t657 * mrSges(5,2) - t689 * t678 + t690 * t679 + t596;
t741 = m(4) * t696 - t697 * mrSges(4,1) + t698 * mrSges(4,2) - t708 * t701 + t709 * t702 + t743;
t739 = -m(3) * t706 + mrSges(3,1) * t758 - t741 + (t723 * t738 + t763) * mrSges(3,3);
t592 = -t738 * mrSges(2,2) + m(2) * t713 + t739 + (mrSges(2,1) - t764) * qJDD(1);
t762 = t732 * t573 + t737 * t592;
t574 = t727 * t579 + t726 * t580;
t747 = Ifges(3,5) * t726 + Ifges(3,6) * t727;
t761 = t738 * t747;
t756 = t737 * t573 - t732 * t592;
t749 = Ifges(3,1) * t726 + Ifges(3,4) * t727;
t748 = Ifges(3,4) * t726 + Ifges(3,2) * t727;
t688 = Ifges(4,1) * t709 + Ifges(4,4) * t708 + Ifges(4,5) * qJD(3);
t687 = Ifges(4,4) * t709 + Ifges(4,2) * t708 + Ifges(4,6) * qJD(3);
t686 = Ifges(4,5) * t709 + Ifges(4,6) * t708 + Ifges(4,3) * qJD(3);
t668 = Ifges(5,1) * t690 + Ifges(5,4) * t689 + Ifges(5,5) * t725;
t667 = Ifges(5,4) * t690 + Ifges(5,2) * t689 + Ifges(5,6) * t725;
t666 = Ifges(5,5) * t690 + Ifges(5,6) * t689 + Ifges(5,3) * t725;
t647 = Ifges(6,1) * t675 + Ifges(6,4) * t674 + Ifges(6,5) * t685;
t646 = Ifges(6,4) * t675 + Ifges(6,2) * t674 + Ifges(6,6) * t685;
t645 = Ifges(6,5) * t675 + Ifges(6,6) * t674 + Ifges(6,3) * t685;
t632 = Ifges(7,1) * t659 + Ifges(7,4) * t658 + Ifges(7,5) * t683;
t631 = Ifges(7,4) * t659 + Ifges(7,2) * t658 + Ifges(7,6) * t683;
t630 = Ifges(7,5) * t659 + Ifges(7,6) * t658 + Ifges(7,3) * t683;
t604 = mrSges(7,2) * t618 - mrSges(7,3) * t612 + Ifges(7,1) * t624 + Ifges(7,4) * t623 + Ifges(7,5) * t652 + t658 * t630 - t683 * t631;
t603 = -mrSges(7,1) * t618 + mrSges(7,3) * t613 + Ifges(7,4) * t624 + Ifges(7,2) * t623 + Ifges(7,6) * t652 - t659 * t630 + t683 * t632;
t590 = mrSges(6,2) * t620 - mrSges(6,3) * t616 + Ifges(6,1) * t639 + Ifges(6,4) * t638 + Ifges(6,5) * t655 - pkin(10) * t602 - t728 * t603 + t733 * t604 + t674 * t645 - t685 * t646;
t589 = -mrSges(6,1) * t620 + mrSges(6,3) * t617 + Ifges(6,4) * t639 + Ifges(6,2) * t638 + Ifges(6,6) * t655 - pkin(5) * t742 + pkin(10) * t751 + t733 * t603 + t728 * t604 - t675 * t645 + t685 * t647;
t582 = Ifges(5,4) * t657 + Ifges(5,2) * t656 + Ifges(5,6) * t722 - t690 * t666 + t725 * t668 - mrSges(5,1) * t650 + mrSges(5,3) * t629 - Ifges(6,5) * t639 - Ifges(6,6) * t638 - Ifges(6,3) * t655 - t675 * t646 + t674 * t647 - mrSges(6,1) * t616 + mrSges(6,2) * t617 - Ifges(7,5) * t624 - Ifges(7,6) * t623 - Ifges(7,3) * t652 - t659 * t631 + t658 * t632 - mrSges(7,1) * t612 + mrSges(7,2) * t613 - pkin(5) * t602 - pkin(4) * t596;
t575 = mrSges(5,2) * t650 - mrSges(5,3) * t628 + Ifges(5,1) * t657 + Ifges(5,4) * t656 + Ifges(5,5) * t722 - pkin(9) * t596 - t729 * t589 + t734 * t590 + t689 * t666 - t725 * t667;
t570 = mrSges(4,2) * t696 - mrSges(4,3) * t661 + Ifges(4,1) * t698 + Ifges(4,4) * t697 + Ifges(4,5) * qJDD(3) - pkin(8) * t588 - qJD(3) * t687 + t735 * t575 - t730 * t582 + t708 * t686;
t569 = -mrSges(4,1) * t696 + mrSges(4,3) * t662 + Ifges(4,4) * t698 + Ifges(4,2) * t697 + Ifges(4,6) * qJDD(3) - pkin(3) * t743 + pkin(8) * t753 + qJD(3) * t688 + t730 * t575 + t735 * t582 - t709 * t686;
t568 = (-t726 * t748 + t727 * t749 + Ifges(2,5)) * t738 - pkin(9) * t752 + (Ifges(2,6) - t747) * qJDD(1) + mrSges(2,1) * g(3) - t734 * t589 - t729 * t590 - Ifges(5,3) * t722 + t708 * t688 - t709 * t687 + mrSges(2,3) * t714 - mrSges(3,1) * t699 + mrSges(3,2) * t700 - pkin(1) * t574 + t689 * t668 - t690 * t667 - Ifges(4,6) * t697 - Ifges(4,5) * t698 - mrSges(4,1) * t661 + mrSges(4,2) * t662 - Ifges(5,6) * t656 - Ifges(5,5) * t657 + mrSges(5,2) * t629 - mrSges(5,1) * t628 - pkin(3) * t588 - pkin(2) * t581 - pkin(4) * t740 - Ifges(4,3) * qJDD(3);
t567 = mrSges(3,2) * t706 - mrSges(3,3) * t699 - pkin(7) * t581 + t749 * qJDD(1) - t731 * t569 + t736 * t570 + t727 * t761;
t566 = -mrSges(3,1) * t706 + mrSges(3,3) * t700 - pkin(2) * t741 + pkin(7) * t754 + t748 * qJDD(1) + t736 * t569 + t731 * t570 - t726 * t761;
t565 = -mrSges(2,2) * g(3) - mrSges(2,3) * t713 + Ifges(2,5) * qJDD(1) - t738 * Ifges(2,6) - qJ(2) * t574 - t726 * t566 + t727 * t567;
t1 = [-m(1) * g(1) + t756; -m(1) * g(2) + t762; (-m(1) - m(2)) * g(3) + t574; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t762 + t737 * t565 - t732 * t568; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t756 + t732 * t565 + t737 * t568; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t713 - mrSges(2,2) * t714 + t726 * t567 + t727 * t566 + pkin(1) * (-qJDD(1) * t764 + t739) + qJ(2) * t755;];
tauB  = t1;
