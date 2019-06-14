% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RRPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 22:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RRPRRR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:39:02
% EndTime: 2019-05-06 22:39:38
% DurationCPUTime: 34.36s
% Computational Cost: add. (561374->385), mult. (1236895->484), div. (0->0), fcn. (921128->12), ass. (0->148)
t725 = sin(qJ(1));
t730 = cos(qJ(1));
t713 = -g(1) * t730 - g(2) * t725;
t732 = qJD(1) ^ 2;
t697 = -pkin(1) * t732 + qJDD(1) * pkin(7) + t713;
t724 = sin(qJ(2));
t729 = cos(qJ(2));
t680 = -g(3) * t724 + t697 * t729;
t706 = (-mrSges(3,1) * t729 + mrSges(3,2) * t724) * qJD(1);
t745 = qJD(1) * qJD(2);
t716 = t724 * t745;
t708 = qJDD(1) * t729 - t716;
t747 = qJD(1) * t724;
t710 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t747;
t712 = t725 * g(1) - g(2) * t730;
t696 = -qJDD(1) * pkin(1) - t732 * pkin(7) - t712;
t744 = t729 * t745;
t707 = qJDD(1) * t724 + t744;
t663 = (-t707 - t744) * qJ(3) + (-t708 + t716) * pkin(2) + t696;
t705 = (-pkin(2) * t729 - qJ(3) * t724) * qJD(1);
t731 = qJD(2) ^ 2;
t746 = qJD(1) * t729;
t666 = -pkin(2) * t731 + qJDD(2) * qJ(3) + t705 * t746 + t680;
t719 = sin(pkin(11));
t720 = cos(pkin(11));
t702 = qJD(2) * t719 + t720 * t747;
t643 = -0.2e1 * qJD(3) * t702 + t663 * t720 - t719 * t666;
t685 = qJDD(2) * t719 + t707 * t720;
t701 = qJD(2) * t720 - t719 * t747;
t630 = (-t701 * t746 - t685) * pkin(8) + (t701 * t702 - t708) * pkin(3) + t643;
t644 = 0.2e1 * qJD(3) * t701 + t663 * t719 + t666 * t720;
t684 = qJDD(2) * t720 - t707 * t719;
t686 = -pkin(3) * t746 - pkin(8) * t702;
t700 = t701 ^ 2;
t632 = -pkin(3) * t700 + pkin(8) * t684 + t686 * t746 + t644;
t723 = sin(qJ(4));
t728 = cos(qJ(4));
t614 = t630 * t728 - t723 * t632;
t676 = t701 * t728 - t702 * t723;
t651 = qJD(4) * t676 + t684 * t723 + t685 * t728;
t677 = t701 * t723 + t702 * t728;
t704 = qJDD(4) - t708;
t715 = qJD(4) - t746;
t611 = (t676 * t715 - t651) * pkin(9) + (t676 * t677 + t704) * pkin(4) + t614;
t615 = t630 * t723 + t632 * t728;
t650 = -qJD(4) * t677 + t684 * t728 - t685 * t723;
t669 = pkin(4) * t715 - pkin(9) * t677;
t675 = t676 ^ 2;
t613 = -pkin(4) * t675 + pkin(9) * t650 - t669 * t715 + t615;
t722 = sin(qJ(5));
t727 = cos(qJ(5));
t601 = t611 * t727 - t722 * t613;
t658 = t676 * t727 - t677 * t722;
t627 = qJD(5) * t658 + t650 * t722 + t651 * t727;
t659 = t676 * t722 + t677 * t727;
t698 = qJDD(5) + t704;
t714 = qJD(5) + t715;
t599 = (t658 * t714 - t627) * pkin(10) + (t658 * t659 + t698) * pkin(5) + t601;
t602 = t611 * t722 + t613 * t727;
t626 = -qJD(5) * t659 + t650 * t727 - t651 * t722;
t649 = pkin(5) * t714 - pkin(10) * t659;
t657 = t658 ^ 2;
t600 = -pkin(5) * t657 + pkin(10) * t626 - t649 * t714 + t602;
t721 = sin(qJ(6));
t726 = cos(qJ(6));
t597 = t599 * t726 - t600 * t721;
t640 = t658 * t726 - t659 * t721;
t608 = qJD(6) * t640 + t626 * t721 + t627 * t726;
t641 = t658 * t721 + t659 * t726;
t623 = -mrSges(7,1) * t640 + mrSges(7,2) * t641;
t709 = qJD(6) + t714;
t633 = -mrSges(7,2) * t709 + mrSges(7,3) * t640;
t690 = qJDD(6) + t698;
t593 = m(7) * t597 + mrSges(7,1) * t690 - mrSges(7,3) * t608 - t623 * t641 + t633 * t709;
t598 = t599 * t721 + t600 * t726;
t607 = -qJD(6) * t641 + t626 * t726 - t627 * t721;
t634 = mrSges(7,1) * t709 - mrSges(7,3) * t641;
t594 = m(7) * t598 - mrSges(7,2) * t690 + mrSges(7,3) * t607 + t623 * t640 - t634 * t709;
t587 = t593 * t726 + t594 * t721;
t642 = -mrSges(6,1) * t658 + mrSges(6,2) * t659;
t646 = -mrSges(6,2) * t714 + mrSges(6,3) * t658;
t585 = m(6) * t601 + mrSges(6,1) * t698 - mrSges(6,3) * t627 - t642 * t659 + t646 * t714 + t587;
t647 = mrSges(6,1) * t714 - mrSges(6,3) * t659;
t738 = -t593 * t721 + t594 * t726;
t586 = m(6) * t602 - mrSges(6,2) * t698 + mrSges(6,3) * t626 + t642 * t658 - t647 * t714 + t738;
t581 = t585 * t727 + t586 * t722;
t660 = -mrSges(5,1) * t676 + mrSges(5,2) * t677;
t667 = -mrSges(5,2) * t715 + mrSges(5,3) * t676;
t579 = m(5) * t614 + mrSges(5,1) * t704 - mrSges(5,3) * t651 - t660 * t677 + t667 * t715 + t581;
t668 = mrSges(5,1) * t715 - mrSges(5,3) * t677;
t739 = -t585 * t722 + t586 * t727;
t580 = m(5) * t615 - mrSges(5,2) * t704 + mrSges(5,3) * t650 + t660 * t676 - t668 * t715 + t739;
t573 = t579 * t728 + t580 * t723;
t678 = -mrSges(4,1) * t701 + mrSges(4,2) * t702;
t682 = mrSges(4,2) * t746 + mrSges(4,3) * t701;
t571 = m(4) * t643 - mrSges(4,1) * t708 - mrSges(4,3) * t685 - t678 * t702 - t682 * t746 + t573;
t683 = -mrSges(4,1) * t746 - mrSges(4,3) * t702;
t740 = -t579 * t723 + t580 * t728;
t572 = m(4) * t644 + mrSges(4,2) * t708 + mrSges(4,3) * t684 + t678 * t701 + t683 * t746 + t740;
t741 = -t571 * t719 + t572 * t720;
t566 = m(3) * t680 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t708 - qJD(2) * t710 + t706 * t746 + t741;
t679 = -g(3) * t729 - t697 * t724;
t711 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t746;
t665 = -qJDD(2) * pkin(2) - qJ(3) * t731 + t705 * t747 + qJDD(3) - t679;
t645 = -pkin(3) * t684 - pkin(8) * t700 + t686 * t702 + t665;
t622 = -pkin(4) * t650 - pkin(9) * t675 + t669 * t677 + t645;
t604 = -pkin(5) * t626 - pkin(10) * t657 + t649 * t659 + t622;
t737 = m(7) * t604 - mrSges(7,1) * t607 + mrSges(7,2) * t608 - t633 * t640 + t634 * t641;
t736 = m(6) * t622 - mrSges(6,1) * t626 + mrSges(6,2) * t627 - t646 * t658 + t647 * t659 + t737;
t734 = m(5) * t645 - mrSges(5,1) * t650 + mrSges(5,2) * t651 - t667 * t676 + t668 * t677 + t736;
t733 = -m(4) * t665 + mrSges(4,1) * t684 - mrSges(4,2) * t685 + t682 * t701 - t683 * t702 - t734;
t596 = m(3) * t679 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t707 + qJD(2) * t711 - t706 * t747 + t733;
t742 = t566 * t729 - t596 * t724;
t560 = m(2) * t713 - mrSges(2,1) * t732 - qJDD(1) * mrSges(2,2) + t742;
t567 = t571 * t720 + t572 * t719;
t735 = -m(3) * t696 + mrSges(3,1) * t708 - mrSges(3,2) * t707 - t710 * t747 + t711 * t746 - t567;
t563 = m(2) * t712 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t732 + t735;
t748 = t560 * t725 + t563 * t730;
t561 = t566 * t724 + t596 * t729;
t743 = t560 * t730 - t563 * t725;
t695 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t724 + Ifges(3,4) * t729) * qJD(1);
t694 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t724 + Ifges(3,2) * t729) * qJD(1);
t693 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t724 + Ifges(3,6) * t729) * qJD(1);
t672 = Ifges(4,1) * t702 + Ifges(4,4) * t701 - Ifges(4,5) * t746;
t671 = Ifges(4,4) * t702 + Ifges(4,2) * t701 - Ifges(4,6) * t746;
t670 = Ifges(4,5) * t702 + Ifges(4,6) * t701 - Ifges(4,3) * t746;
t654 = Ifges(5,1) * t677 + Ifges(5,4) * t676 + Ifges(5,5) * t715;
t653 = Ifges(5,4) * t677 + Ifges(5,2) * t676 + Ifges(5,6) * t715;
t652 = Ifges(5,5) * t677 + Ifges(5,6) * t676 + Ifges(5,3) * t715;
t637 = Ifges(6,1) * t659 + Ifges(6,4) * t658 + Ifges(6,5) * t714;
t636 = Ifges(6,4) * t659 + Ifges(6,2) * t658 + Ifges(6,6) * t714;
t635 = Ifges(6,5) * t659 + Ifges(6,6) * t658 + Ifges(6,3) * t714;
t618 = Ifges(7,1) * t641 + Ifges(7,4) * t640 + Ifges(7,5) * t709;
t617 = Ifges(7,4) * t641 + Ifges(7,2) * t640 + Ifges(7,6) * t709;
t616 = Ifges(7,5) * t641 + Ifges(7,6) * t640 + Ifges(7,3) * t709;
t589 = mrSges(7,2) * t604 - mrSges(7,3) * t597 + Ifges(7,1) * t608 + Ifges(7,4) * t607 + Ifges(7,5) * t690 + t616 * t640 - t617 * t709;
t588 = -mrSges(7,1) * t604 + mrSges(7,3) * t598 + Ifges(7,4) * t608 + Ifges(7,2) * t607 + Ifges(7,6) * t690 - t616 * t641 + t618 * t709;
t575 = mrSges(6,2) * t622 - mrSges(6,3) * t601 + Ifges(6,1) * t627 + Ifges(6,4) * t626 + Ifges(6,5) * t698 - pkin(10) * t587 - t588 * t721 + t589 * t726 + t635 * t658 - t636 * t714;
t574 = -mrSges(6,1) * t622 + mrSges(6,3) * t602 + Ifges(6,4) * t627 + Ifges(6,2) * t626 + Ifges(6,6) * t698 - pkin(5) * t737 + pkin(10) * t738 + t726 * t588 + t721 * t589 - t659 * t635 + t714 * t637;
t569 = mrSges(5,2) * t645 - mrSges(5,3) * t614 + Ifges(5,1) * t651 + Ifges(5,4) * t650 + Ifges(5,5) * t704 - pkin(9) * t581 - t574 * t722 + t575 * t727 + t652 * t676 - t653 * t715;
t568 = -mrSges(5,1) * t645 + mrSges(5,3) * t615 + Ifges(5,4) * t651 + Ifges(5,2) * t650 + Ifges(5,6) * t704 - pkin(4) * t736 + pkin(9) * t739 + t727 * t574 + t722 * t575 - t677 * t652 + t715 * t654;
t557 = -pkin(4) * t581 - pkin(3) * t573 - t693 * t747 - pkin(2) * t567 + (Ifges(3,2) + Ifges(4,3)) * t708 + Ifges(3,6) * qJDD(2) - Ifges(5,3) * t704 + Ifges(3,4) * t707 - Ifges(6,3) * t698 + t701 * t672 - t702 * t671 - Ifges(4,5) * t685 - Ifges(7,3) * t690 + qJD(2) * t695 - mrSges(3,1) * t696 - t677 * t653 + mrSges(3,3) * t680 - Ifges(4,6) * t684 + t676 * t654 + t658 * t637 - t659 * t636 - Ifges(5,6) * t650 - Ifges(5,5) * t651 - mrSges(4,1) * t643 + mrSges(4,2) * t644 + t640 * t618 - t641 * t617 - Ifges(6,6) * t626 - Ifges(6,5) * t627 - mrSges(5,1) * t614 + mrSges(5,2) * t615 - Ifges(7,6) * t607 - Ifges(7,5) * t608 + mrSges(6,2) * t602 - mrSges(6,1) * t601 + mrSges(7,2) * t598 - mrSges(7,1) * t597 - pkin(5) * t587;
t556 = mrSges(4,2) * t665 - mrSges(4,3) * t643 + Ifges(4,1) * t685 + Ifges(4,4) * t684 - Ifges(4,5) * t708 - pkin(8) * t573 - t568 * t723 + t569 * t728 + t670 * t701 + t671 * t746;
t555 = -mrSges(4,1) * t665 + mrSges(4,3) * t644 + Ifges(4,4) * t685 + Ifges(4,2) * t684 - Ifges(4,6) * t708 - pkin(3) * t734 + pkin(8) * t740 + t728 * t568 + t723 * t569 - t702 * t670 - t672 * t746;
t554 = mrSges(3,2) * t696 - mrSges(3,3) * t679 + Ifges(3,1) * t707 + Ifges(3,4) * t708 + Ifges(3,5) * qJDD(2) - qJ(3) * t567 - qJD(2) * t694 - t555 * t719 + t556 * t720 + t693 * t746;
t553 = Ifges(2,6) * qJDD(1) + t732 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t713 - Ifges(3,5) * t707 - Ifges(3,6) * t708 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t679 + mrSges(3,2) * t680 - t719 * t556 - t720 * t555 - pkin(2) * t733 - qJ(3) * t741 - pkin(1) * t561 + (-t694 * t724 + t695 * t729) * qJD(1);
t552 = -mrSges(2,2) * g(3) - mrSges(2,3) * t712 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t732 - pkin(7) * t561 + t554 * t729 - t557 * t724;
t1 = [-m(1) * g(1) + t743; -m(1) * g(2) + t748; (-m(1) - m(2)) * g(3) + t561; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t748 + t552 * t730 - t553 * t725; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t743 + t725 * t552 + t730 * t553; -mrSges(1,1) * g(2) + mrSges(2,1) * t712 + mrSges(1,2) * g(1) - mrSges(2,2) * t713 + Ifges(2,3) * qJDD(1) + pkin(1) * t735 + pkin(7) * t742 + t724 * t554 + t729 * t557;];
tauB  = t1;
