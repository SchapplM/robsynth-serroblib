% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 09:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRRPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:42:29
% EndTime: 2019-05-07 09:43:11
% DurationCPUTime: 33.47s
% Computational Cost: add. (494079->387), mult. (1155364->490), div. (0->0), fcn. (879362->12), ass. (0->149)
t741 = qJD(1) ^ 2;
t759 = pkin(2) * t741;
t735 = sin(qJ(1));
t740 = cos(qJ(1));
t721 = -g(1) * t740 - g(2) * t735;
t710 = -pkin(1) * t741 + qJDD(1) * pkin(7) + t721;
t734 = sin(qJ(2));
t758 = t734 * t710;
t739 = cos(qJ(2));
t754 = qJD(1) * qJD(2);
t715 = qJDD(1) * t734 + t739 * t754;
t677 = qJDD(2) * pkin(2) - t715 * pkin(8) - t758 + (pkin(8) * t754 + t734 * t759 - g(3)) * t739;
t698 = -g(3) * t734 + t739 * t710;
t716 = qJDD(1) * t739 - t734 * t754;
t756 = qJD(1) * t734;
t719 = qJD(2) * pkin(2) - pkin(8) * t756;
t728 = t739 ^ 2;
t678 = pkin(8) * t716 - qJD(2) * t719 - t728 * t759 + t698;
t733 = sin(qJ(3));
t738 = cos(qJ(3));
t654 = t738 * t677 - t733 * t678;
t707 = (-t733 * t734 + t738 * t739) * qJD(1);
t681 = qJD(3) * t707 + t715 * t738 + t716 * t733;
t708 = (t733 * t739 + t734 * t738) * qJD(1);
t726 = qJDD(2) + qJDD(3);
t727 = qJD(2) + qJD(3);
t639 = (t707 * t727 - t681) * qJ(4) + (t707 * t708 + t726) * pkin(3) + t654;
t655 = t733 * t677 + t738 * t678;
t680 = -qJD(3) * t708 - t715 * t733 + t716 * t738;
t700 = pkin(3) * t727 - qJ(4) * t708;
t703 = t707 ^ 2;
t641 = -pkin(3) * t703 + qJ(4) * t680 - t700 * t727 + t655;
t729 = sin(pkin(11));
t730 = cos(pkin(11));
t695 = t707 * t729 + t708 * t730;
t621 = -0.2e1 * qJD(4) * t695 + t730 * t639 - t729 * t641;
t659 = t680 * t729 + t681 * t730;
t694 = t707 * t730 - t708 * t729;
t618 = (t694 * t727 - t659) * pkin(9) + (t694 * t695 + t726) * pkin(4) + t621;
t622 = 0.2e1 * qJD(4) * t694 + t729 * t639 + t730 * t641;
t658 = t680 * t730 - t681 * t729;
t685 = pkin(4) * t727 - pkin(9) * t695;
t691 = t694 ^ 2;
t620 = -pkin(4) * t691 + pkin(9) * t658 - t685 * t727 + t622;
t732 = sin(qJ(5));
t737 = cos(qJ(5));
t615 = t732 * t618 + t737 * t620;
t672 = t694 * t732 + t695 * t737;
t631 = -qJD(5) * t672 + t658 * t737 - t659 * t732;
t671 = t694 * t737 - t695 * t732;
t650 = -mrSges(6,1) * t671 + mrSges(6,2) * t672;
t724 = qJD(5) + t727;
t663 = mrSges(6,1) * t724 - mrSges(6,3) * t672;
t723 = qJDD(5) + t726;
t651 = -pkin(5) * t671 - pkin(10) * t672;
t722 = t724 ^ 2;
t613 = -pkin(5) * t722 + pkin(10) * t723 + t651 * t671 + t615;
t720 = t735 * g(1) - t740 * g(2);
t747 = -qJDD(1) * pkin(1) - t720;
t682 = -t716 * pkin(2) + t719 * t756 + (-pkin(8) * t728 - pkin(7)) * t741 + t747;
t653 = -t680 * pkin(3) - t703 * qJ(4) + t708 * t700 + qJDD(4) + t682;
t627 = -t658 * pkin(4) - t691 * pkin(9) + t695 * t685 + t653;
t632 = qJD(5) * t671 + t658 * t732 + t659 * t737;
t616 = (-t671 * t724 - t632) * pkin(10) + (t672 * t724 - t631) * pkin(5) + t627;
t731 = sin(qJ(6));
t736 = cos(qJ(6));
t610 = -t613 * t731 + t616 * t736;
t660 = -t672 * t731 + t724 * t736;
t625 = qJD(6) * t660 + t632 * t736 + t723 * t731;
t630 = qJDD(6) - t631;
t661 = t672 * t736 + t724 * t731;
t642 = -mrSges(7,1) * t660 + mrSges(7,2) * t661;
t667 = qJD(6) - t671;
t643 = -mrSges(7,2) * t667 + mrSges(7,3) * t660;
t608 = m(7) * t610 + mrSges(7,1) * t630 - mrSges(7,3) * t625 - t642 * t661 + t643 * t667;
t611 = t613 * t736 + t616 * t731;
t624 = -qJD(6) * t661 - t632 * t731 + t723 * t736;
t644 = mrSges(7,1) * t667 - mrSges(7,3) * t661;
t609 = m(7) * t611 - mrSges(7,2) * t630 + mrSges(7,3) * t624 + t642 * t660 - t644 * t667;
t748 = -t608 * t731 + t736 * t609;
t599 = m(6) * t615 - mrSges(6,2) * t723 + mrSges(6,3) * t631 + t650 * t671 - t663 * t724 + t748;
t614 = t618 * t737 - t620 * t732;
t662 = -mrSges(6,2) * t724 + mrSges(6,3) * t671;
t612 = -pkin(5) * t723 - pkin(10) * t722 + t651 * t672 - t614;
t745 = -m(7) * t612 + t624 * mrSges(7,1) - mrSges(7,2) * t625 + t660 * t643 - t644 * t661;
t604 = m(6) * t614 + mrSges(6,1) * t723 - mrSges(6,3) * t632 - t650 * t672 + t662 * t724 + t745;
t594 = t732 * t599 + t737 * t604;
t673 = -mrSges(5,1) * t694 + mrSges(5,2) * t695;
t683 = -mrSges(5,2) * t727 + mrSges(5,3) * t694;
t592 = m(5) * t621 + mrSges(5,1) * t726 - mrSges(5,3) * t659 - t673 * t695 + t683 * t727 + t594;
t684 = mrSges(5,1) * t727 - mrSges(5,3) * t695;
t749 = t737 * t599 - t604 * t732;
t593 = m(5) * t622 - mrSges(5,2) * t726 + mrSges(5,3) * t658 + t673 * t694 - t684 * t727 + t749;
t586 = t730 * t592 + t729 * t593;
t696 = -mrSges(4,1) * t707 + mrSges(4,2) * t708;
t699 = -mrSges(4,2) * t727 + mrSges(4,3) * t707;
t584 = m(4) * t654 + mrSges(4,1) * t726 - mrSges(4,3) * t681 - t696 * t708 + t699 * t727 + t586;
t701 = mrSges(4,1) * t727 - mrSges(4,3) * t708;
t750 = -t592 * t729 + t730 * t593;
t585 = m(4) * t655 - mrSges(4,2) * t726 + mrSges(4,3) * t680 + t696 * t707 - t701 * t727 + t750;
t579 = t738 * t584 + t733 * t585;
t697 = -t739 * g(3) - t758;
t714 = (-mrSges(3,1) * t739 + mrSges(3,2) * t734) * qJD(1);
t755 = qJD(1) * t739;
t718 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t755;
t577 = m(3) * t697 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t715 + qJD(2) * t718 - t714 * t756 + t579;
t717 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t756;
t751 = -t584 * t733 + t738 * t585;
t578 = m(3) * t698 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t716 - qJD(2) * t717 + t714 * t755 + t751;
t752 = -t577 * t734 + t739 * t578;
t571 = m(2) * t721 - mrSges(2,1) * t741 - qJDD(1) * mrSges(2,2) + t752;
t709 = -t741 * pkin(7) + t747;
t600 = t736 * t608 + t731 * t609;
t746 = m(6) * t627 - t631 * mrSges(6,1) + t632 * mrSges(6,2) - t671 * t662 + t672 * t663 + t600;
t744 = m(5) * t653 - t658 * mrSges(5,1) + t659 * mrSges(5,2) - t694 * t683 + t695 * t684 + t746;
t743 = m(4) * t682 - t680 * mrSges(4,1) + t681 * mrSges(4,2) - t707 * t699 + t708 * t701 + t744;
t742 = -m(3) * t709 + t716 * mrSges(3,1) - t715 * mrSges(3,2) - t717 * t756 + t718 * t755 - t743;
t596 = m(2) * t720 + qJDD(1) * mrSges(2,1) - t741 * mrSges(2,2) + t742;
t757 = t735 * t571 + t740 * t596;
t572 = t739 * t577 + t734 * t578;
t753 = t740 * t571 - t596 * t735;
t706 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t734 + Ifges(3,4) * t739) * qJD(1);
t705 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t734 + Ifges(3,2) * t739) * qJD(1);
t704 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t734 + Ifges(3,6) * t739) * qJD(1);
t690 = Ifges(4,1) * t708 + Ifges(4,4) * t707 + Ifges(4,5) * t727;
t689 = Ifges(4,4) * t708 + Ifges(4,2) * t707 + Ifges(4,6) * t727;
t688 = Ifges(4,5) * t708 + Ifges(4,6) * t707 + Ifges(4,3) * t727;
t666 = Ifges(5,1) * t695 + Ifges(5,4) * t694 + Ifges(5,5) * t727;
t665 = Ifges(5,4) * t695 + Ifges(5,2) * t694 + Ifges(5,6) * t727;
t664 = Ifges(5,5) * t695 + Ifges(5,6) * t694 + Ifges(5,3) * t727;
t647 = Ifges(6,1) * t672 + Ifges(6,4) * t671 + Ifges(6,5) * t724;
t646 = Ifges(6,4) * t672 + Ifges(6,2) * t671 + Ifges(6,6) * t724;
t645 = Ifges(6,5) * t672 + Ifges(6,6) * t671 + Ifges(6,3) * t724;
t635 = Ifges(7,1) * t661 + Ifges(7,4) * t660 + Ifges(7,5) * t667;
t634 = Ifges(7,4) * t661 + Ifges(7,2) * t660 + Ifges(7,6) * t667;
t633 = Ifges(7,5) * t661 + Ifges(7,6) * t660 + Ifges(7,3) * t667;
t602 = mrSges(7,2) * t612 - mrSges(7,3) * t610 + Ifges(7,1) * t625 + Ifges(7,4) * t624 + Ifges(7,5) * t630 + t633 * t660 - t634 * t667;
t601 = -mrSges(7,1) * t612 + mrSges(7,3) * t611 + Ifges(7,4) * t625 + Ifges(7,2) * t624 + Ifges(7,6) * t630 - t633 * t661 + t635 * t667;
t588 = -mrSges(6,1) * t627 - mrSges(7,1) * t610 + mrSges(7,2) * t611 + mrSges(6,3) * t615 + Ifges(6,4) * t632 - Ifges(7,5) * t625 + Ifges(6,2) * t631 + Ifges(6,6) * t723 - Ifges(7,6) * t624 - Ifges(7,3) * t630 - pkin(5) * t600 - t634 * t661 + t635 * t660 - t645 * t672 + t647 * t724;
t587 = mrSges(6,2) * t627 - mrSges(6,3) * t614 + Ifges(6,1) * t632 + Ifges(6,4) * t631 + Ifges(6,5) * t723 - pkin(10) * t600 - t601 * t731 + t602 * t736 + t645 * t671 - t646 * t724;
t580 = mrSges(5,2) * t653 - mrSges(5,3) * t621 + Ifges(5,1) * t659 + Ifges(5,4) * t658 + Ifges(5,5) * t726 - pkin(9) * t594 + t587 * t737 - t588 * t732 + t664 * t694 - t665 * t727;
t573 = -mrSges(5,1) * t653 + mrSges(5,3) * t622 + Ifges(5,4) * t659 + Ifges(5,2) * t658 + Ifges(5,6) * t726 - pkin(4) * t746 + pkin(9) * t749 + t732 * t587 + t737 * t588 - t695 * t664 + t727 * t666;
t568 = mrSges(4,2) * t682 - mrSges(4,3) * t654 + Ifges(4,1) * t681 + Ifges(4,4) * t680 + Ifges(4,5) * t726 - qJ(4) * t586 - t573 * t729 + t580 * t730 + t688 * t707 - t689 * t727;
t567 = -mrSges(4,1) * t682 + mrSges(4,3) * t655 + Ifges(4,4) * t681 + Ifges(4,2) * t680 + Ifges(4,6) * t726 - pkin(3) * t744 + qJ(4) * t750 + t730 * t573 + t729 * t580 - t708 * t688 + t727 * t690;
t566 = -pkin(1) * t572 - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) - pkin(3) * t586 - mrSges(5,1) * t621 + mrSges(5,2) * t622 - pkin(10) * t748 - pkin(5) * t745 + Ifges(2,6) * qJDD(1) - pkin(4) * t594 - pkin(2) * t579 - mrSges(6,1) * t614 + mrSges(6,2) * t615 + (-Ifges(4,3) - Ifges(5,3)) * t726 + (-t705 * t734 + t706 * t739) * qJD(1) - Ifges(6,6) * t631 - Ifges(6,5) * t632 - mrSges(4,1) * t654 + mrSges(4,2) * t655 - Ifges(5,6) * t658 - Ifges(5,5) * t659 + t671 * t647 - t672 * t646 - Ifges(4,6) * t680 - Ifges(4,5) * t681 + t694 * t666 - t695 * t665 - mrSges(3,1) * t697 + mrSges(3,2) * t698 + t707 * t690 - t708 * t689 - Ifges(3,5) * t715 - Ifges(3,6) * t716 + mrSges(2,3) * t721 - Ifges(6,3) * t723 - t731 * t602 - t736 * t601 + t741 * Ifges(2,5);
t565 = mrSges(3,2) * t709 - mrSges(3,3) * t697 + Ifges(3,1) * t715 + Ifges(3,4) * t716 + Ifges(3,5) * qJDD(2) - pkin(8) * t579 - qJD(2) * t705 - t567 * t733 + t568 * t738 + t704 * t755;
t564 = -mrSges(3,1) * t709 + mrSges(3,3) * t698 + Ifges(3,4) * t715 + Ifges(3,2) * t716 + Ifges(3,6) * qJDD(2) - pkin(2) * t743 + pkin(8) * t751 + qJD(2) * t706 + t738 * t567 + t733 * t568 - t704 * t756;
t563 = -mrSges(2,2) * g(3) - mrSges(2,3) * t720 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t741 - pkin(7) * t572 - t564 * t734 + t565 * t739;
t1 = [-m(1) * g(1) + t753; -m(1) * g(2) + t757; (-m(1) - m(2)) * g(3) + t572; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t757 + t740 * t563 - t735 * t566; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t753 + t735 * t563 + t740 * t566; -mrSges(1,1) * g(2) + mrSges(2,1) * t720 + mrSges(1,2) * g(1) - mrSges(2,2) * t721 + Ifges(2,3) * qJDD(1) + pkin(1) * t742 + pkin(7) * t752 + t739 * t564 + t734 * t565;];
tauB  = t1;
