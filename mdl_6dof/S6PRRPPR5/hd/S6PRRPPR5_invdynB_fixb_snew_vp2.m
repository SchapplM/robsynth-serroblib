% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-05-05 03:29
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:25:01
% EndTime: 2019-05-05 03:25:11
% DurationCPUTime: 7.27s
% Computational Cost: add. (110015->325), mult. (229527->402), div. (0->0), fcn. (143666->12), ass. (0->140)
t742 = -2 * qJD(4);
t741 = Ifges(4,1) + Ifges(5,2);
t737 = Ifges(4,4) + Ifges(5,6);
t736 = Ifges(4,5) - Ifges(5,4);
t740 = Ifges(4,2) + Ifges(5,3);
t735 = Ifges(4,6) - Ifges(5,5);
t739 = (Ifges(4,3) + Ifges(5,1));
t690 = sin(pkin(10));
t693 = cos(pkin(10));
t672 = t690 * g(1) - t693 * g(2);
t673 = -t693 * g(1) - t690 * g(2);
t688 = -g(3) + qJDD(1);
t700 = cos(qJ(2));
t694 = cos(pkin(6));
t697 = sin(qJ(2));
t730 = t694 * t697;
t691 = sin(pkin(6));
t731 = t691 * t697;
t626 = t672 * t730 + t700 * t673 + t688 * t731;
t702 = qJD(2) ^ 2;
t622 = -t702 * pkin(2) + qJDD(2) * pkin(8) + t626;
t641 = -t691 * t672 + t694 * t688;
t696 = sin(qJ(3));
t699 = cos(qJ(3));
t617 = t699 * t622 + t696 * t641;
t665 = (-pkin(3) * t699 - qJ(4) * t696) * qJD(2);
t701 = qJD(3) ^ 2;
t723 = qJD(2) * t699;
t606 = t701 * pkin(3) - qJDD(3) * qJ(4) + (qJD(3) * t742) - t665 * t723 - t617;
t625 = -t697 * t673 + (t672 * t694 + t688 * t691) * t700;
t738 = t702 * pkin(8);
t734 = -pkin(3) - qJ(5);
t733 = qJ(5) * t702;
t707 = -qJDD(2) * pkin(2) - t625;
t621 = t707 - t738;
t721 = qJD(2) * qJD(3);
t718 = t699 * t721;
t668 = t696 * qJDD(2) + t718;
t719 = t696 * t721;
t669 = t699 * qJDD(2) - t719;
t722 = t696 * qJD(2);
t674 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t722;
t675 = -(qJD(3) * mrSges(4,2)) + mrSges(4,3) * t723;
t677 = -mrSges(5,1) * t723 - (qJD(3) * mrSges(5,3));
t705 = pkin(3) * t719 + t722 * t742 + (-t668 - t718) * qJ(4) + t707;
t608 = -t669 * pkin(3) + t705 - t738;
t678 = mrSges(5,1) * t722 + (qJD(3) * mrSges(5,2));
t619 = t696 * t622;
t713 = -t701 * qJ(4) + t665 * t722 + qJDD(4) + t619;
t603 = t668 * pkin(4) + t734 * qJDD(3) + (-pkin(4) * t721 - t696 * t733 - t641) * t699 + t713;
t676 = pkin(4) * t722 - qJD(3) * qJ(5);
t687 = t699 ^ 2;
t605 = -t676 * t722 + (-pkin(4) * t687 - pkin(8)) * t702 + t734 * t669 + t705;
t689 = sin(pkin(11));
t692 = cos(pkin(11));
t658 = t692 * qJD(3) - t689 * t723;
t595 = -0.2e1 * qJD(5) * t658 + t692 * t603 - t689 * t605;
t639 = t692 * qJDD(3) - t689 * t669;
t657 = -t689 * qJD(3) - t692 * t723;
t593 = (t657 * t722 - t639) * pkin(9) + (t657 * t658 + t668) * pkin(5) + t595;
t596 = 0.2e1 * qJD(5) * t657 + t689 * t603 + t692 * t605;
t638 = -t689 * qJDD(3) - t692 * t669;
t640 = pkin(5) * t722 - t658 * pkin(9);
t656 = t657 ^ 2;
t594 = -t656 * pkin(5) + t638 * pkin(9) - t640 * t722 + t596;
t695 = sin(qJ(6));
t698 = cos(qJ(6));
t591 = t698 * t593 - t695 * t594;
t631 = t698 * t657 - t695 * t658;
t611 = t631 * qJD(6) + t695 * t638 + t698 * t639;
t632 = t695 * t657 + t698 * t658;
t618 = -t631 * mrSges(7,1) + t632 * mrSges(7,2);
t681 = qJD(6) + t722;
t623 = -t681 * mrSges(7,2) + t631 * mrSges(7,3);
t662 = qJDD(6) + t668;
t589 = m(7) * t591 + t662 * mrSges(7,1) - t611 * mrSges(7,3) - t632 * t618 + t681 * t623;
t592 = t695 * t593 + t698 * t594;
t610 = -t632 * qJD(6) + t698 * t638 - t695 * t639;
t624 = t681 * mrSges(7,1) - t632 * mrSges(7,3);
t590 = m(7) * t592 - t662 * mrSges(7,2) + t610 * mrSges(7,3) + t631 * t618 - t681 * t624;
t580 = t698 * t589 + t695 * t590;
t633 = -t657 * mrSges(6,1) + t658 * mrSges(6,2);
t636 = -mrSges(6,2) * t722 + t657 * mrSges(6,3);
t578 = m(6) * t595 + t668 * mrSges(6,1) - t639 * mrSges(6,3) - t658 * t633 + t636 * t722 + t580;
t637 = mrSges(6,1) * t722 - t658 * mrSges(6,3);
t715 = -t695 * t589 + t698 * t590;
t579 = m(6) * t596 - t668 * mrSges(6,2) + t638 * mrSges(6,3) + t657 * t633 - t637 * t722 + t715;
t727 = -t689 * t578 + t692 * t579;
t712 = -m(5) * t608 - t669 * mrSges(5,2) + t678 * t722 - t727;
t703 = -m(4) * t621 + t675 * t723 + t669 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t668 + (-t674 * t696 - t677 * t699) * qJD(2) + t712;
t571 = m(3) * t625 + qJDD(2) * mrSges(3,1) - t702 * mrSges(3,2) + t703;
t732 = t571 * t700;
t729 = t699 * t641;
t616 = -t619 + t729;
t666 = (mrSges(5,2) * t699 - mrSges(5,3) * t696) * qJD(2);
t667 = (-mrSges(4,1) * t699 + mrSges(4,2) * t696) * qJD(2);
t575 = t692 * t578 + t689 * t579;
t607 = -qJDD(3) * pkin(3) + t713 - t729;
t708 = -m(5) * t607 - t668 * mrSges(5,1) - t575;
t573 = m(4) * t616 - t668 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t675 - t677) * qJD(3) + (-t666 - t667) * t722 + t708;
t602 = t669 * pkin(4) + qJD(3) * t676 - t687 * t733 + qJDD(5) - t606;
t598 = -t638 * pkin(5) - t656 * pkin(9) + t658 * t640 + t602;
t709 = m(7) * t598 - t610 * mrSges(7,1) + t611 * mrSges(7,2) - t631 * t623 + t632 * t624;
t706 = -m(6) * t602 + t638 * mrSges(6,1) - t639 * mrSges(6,2) + t657 * t636 - t658 * t637 - t709;
t704 = -m(5) * t606 + qJDD(3) * mrSges(5,3) + qJD(3) * t678 + t666 * t723 - t706;
t585 = t704 + t667 * t723 - qJD(3) * t674 + m(4) * t617 - qJDD(3) * mrSges(4,2) + (mrSges(4,3) + mrSges(5,1)) * t669;
t716 = -t696 * t573 + t699 * t585;
t563 = m(3) * t626 - t702 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t716;
t566 = t699 * t573 + t696 * t585;
t565 = m(3) * t641 + t566;
t554 = t563 * t730 - t691 * t565 + t694 * t732;
t552 = m(2) * t672 + t554;
t559 = t700 * t563 - t697 * t571;
t558 = m(2) * t673 + t559;
t728 = t693 * t552 + t690 * t558;
t726 = (t739 * qJD(3)) + (t736 * t696 + t735 * t699) * qJD(2);
t725 = -t735 * qJD(3) + (-t737 * t696 - t740 * t699) * qJD(2);
t724 = t736 * qJD(3) + (t741 * t696 + t737 * t699) * qJD(2);
t553 = t563 * t731 + t694 * t565 + t691 * t732;
t717 = -t690 * t552 + t693 * t558;
t612 = Ifges(7,5) * t632 + Ifges(7,6) * t631 + Ifges(7,3) * t681;
t614 = Ifges(7,1) * t632 + Ifges(7,4) * t631 + Ifges(7,5) * t681;
t581 = -mrSges(7,1) * t598 + mrSges(7,3) * t592 + Ifges(7,4) * t611 + Ifges(7,2) * t610 + Ifges(7,6) * t662 - t632 * t612 + t681 * t614;
t613 = Ifges(7,4) * t632 + Ifges(7,2) * t631 + Ifges(7,6) * t681;
t582 = mrSges(7,2) * t598 - mrSges(7,3) * t591 + Ifges(7,1) * t611 + Ifges(7,4) * t610 + Ifges(7,5) * t662 + t631 * t612 - t681 * t613;
t627 = Ifges(6,5) * t658 + Ifges(6,6) * t657 + Ifges(6,3) * t722;
t629 = Ifges(6,1) * t658 + Ifges(6,4) * t657 + Ifges(6,5) * t722;
t567 = -mrSges(6,1) * t602 + mrSges(6,3) * t596 + Ifges(6,4) * t639 + Ifges(6,2) * t638 + Ifges(6,6) * t668 - pkin(5) * t709 + pkin(9) * t715 + t698 * t581 + t695 * t582 - t658 * t627 + t629 * t722;
t628 = Ifges(6,4) * t658 + Ifges(6,2) * t657 + Ifges(6,6) * t722;
t568 = mrSges(6,2) * t602 - mrSges(6,3) * t595 + Ifges(6,1) * t639 + Ifges(6,4) * t638 + Ifges(6,5) * t668 - pkin(9) * t580 - t695 * t581 + t698 * t582 + t657 * t627 - t628 * t722;
t574 = -t668 * mrSges(5,3) + t677 * t723 - t712;
t550 = -mrSges(4,1) * t621 - mrSges(5,1) * t606 + mrSges(5,2) * t608 + mrSges(4,3) * t617 - pkin(3) * t574 - pkin(4) * t706 - qJ(5) * t727 + t724 * qJD(3) + t735 * qJDD(3) - t692 * t567 - t689 * t568 + t737 * t668 + t740 * t669 - t726 * t722;
t555 = t736 * qJDD(3) + t737 * t669 + t725 * qJD(3) + t726 * t723 + t658 * t628 + Ifges(7,3) * t662 - t657 * t629 + t632 * t613 + Ifges(6,6) * t638 + Ifges(6,5) * t639 + mrSges(4,2) * t621 - t631 * t614 + Ifges(7,6) * t610 + Ifges(7,5) * t611 - mrSges(4,3) * t616 + mrSges(5,1) * t607 - mrSges(5,3) * t608 + mrSges(6,1) * t595 - mrSges(6,2) * t596 + mrSges(7,1) * t591 - mrSges(7,2) * t592 + pkin(5) * t580 + pkin(4) * t575 - qJ(4) * t574 + (Ifges(6,3) + t741) * t668;
t548 = mrSges(3,2) * t641 - mrSges(3,3) * t625 + Ifges(3,5) * qJDD(2) - t702 * Ifges(3,6) - pkin(8) * t566 - t696 * t550 + t699 * t555;
t549 = -pkin(2) * t566 + mrSges(3,3) * t626 - mrSges(3,1) * t641 + Ifges(3,6) * qJDD(2) - pkin(3) * (-qJD(3) * t677 + t708) - qJ(4) * t704 + t689 * t567 + qJ(5) * t575 - mrSges(4,1) * t616 + mrSges(4,2) * t617 - mrSges(5,2) * t607 + mrSges(5,3) * t606 - t692 * t568 + t702 * Ifges(3,5) + (-qJ(4) * mrSges(5,1) - t735) * t669 - t736 * t668 + (pkin(3) * mrSges(5,2) - t739) * qJDD(3) + (t724 * t699 + (pkin(3) * t666 + t725) * t696) * qJD(2);
t710 = pkin(7) * t559 + t548 * t697 + t549 * t700;
t547 = mrSges(3,1) * t625 - mrSges(3,2) * t626 + Ifges(3,3) * qJDD(2) + pkin(2) * t703 + pkin(8) * t716 + t699 * t550 + t696 * t555;
t546 = mrSges(2,2) * t688 - mrSges(2,3) * t672 + t700 * t548 - t697 * t549 + (-t553 * t691 - t554 * t694) * pkin(7);
t545 = -mrSges(2,1) * t688 + mrSges(2,3) * t673 - pkin(1) * t553 - t691 * t547 + t710 * t694;
t1 = [-m(1) * g(1) + t717; -m(1) * g(2) + t728; -m(1) * g(3) + m(2) * t688 + t553; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t728 - t690 * t545 + t693 * t546; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t717 + t693 * t545 + t690 * t546; -mrSges(1,1) * g(2) + mrSges(2,1) * t672 + mrSges(1,2) * g(1) - mrSges(2,2) * t673 + pkin(1) * t554 + t694 * t547 + t710 * t691;];
tauB  = t1;
