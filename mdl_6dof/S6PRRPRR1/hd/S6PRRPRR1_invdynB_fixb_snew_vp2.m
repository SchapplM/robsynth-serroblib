% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:15:15
% EndTime: 2019-05-05 04:15:37
% DurationCPUTime: 21.41s
% Computational Cost: add. (351454->341), mult. (771190->443), div. (0->0), fcn. (574298->14), ass. (0->141)
t672 = sin(pkin(11));
t675 = cos(pkin(11));
t660 = g(1) * t672 - g(2) * t675;
t661 = -g(1) * t675 - g(2) * t672;
t670 = -g(3) + qJDD(1);
t673 = sin(pkin(6));
t676 = cos(pkin(6));
t680 = sin(qJ(2));
t684 = cos(qJ(2));
t629 = -t680 * t661 + (t660 * t676 + t670 * t673) * t684;
t685 = qJD(2) ^ 2;
t688 = -qJDD(2) * pkin(2) - t629;
t621 = -t685 * pkin(8) + t688;
t679 = sin(qJ(3));
t683 = cos(qJ(3));
t699 = qJD(2) * qJD(3);
t698 = t683 * t699;
t658 = qJDD(2) * t679 + t698;
t659 = qJDD(2) * t683 - t679 * t699;
t701 = qJD(2) * t679;
t663 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t701;
t700 = qJD(2) * t683;
t664 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t700;
t662 = qJD(3) * pkin(3) - qJ(4) * t701;
t669 = t683 ^ 2;
t613 = -t659 * pkin(3) + qJDD(4) + t662 * t701 + (-qJ(4) * t669 - pkin(8)) * t685 + t688;
t671 = sin(pkin(12));
t674 = cos(pkin(12));
t635 = -t658 * t671 + t659 * t674;
t636 = t658 * t674 + t659 * t671;
t646 = (-t671 * t679 + t674 * t683) * qJD(2);
t639 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t646;
t647 = (t671 * t683 + t674 * t679) * qJD(2);
t640 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t647;
t703 = t676 * t680;
t704 = t673 * t680;
t630 = t660 * t703 + t684 * t661 + t670 * t704;
t622 = -pkin(2) * t685 + qJDD(2) * pkin(8) + t630;
t642 = -t660 * t673 + t670 * t676;
t614 = -t679 * t622 + t683 * t642;
t601 = (-t658 + t698) * qJ(4) + (t679 * t683 * t685 + qJDD(3)) * pkin(3) + t614;
t615 = t683 * t622 + t679 * t642;
t603 = -pkin(3) * t669 * t685 + qJ(4) * t659 - qJD(3) * t662 + t615;
t582 = -0.2e1 * qJD(4) * t647 + t674 * t601 - t671 * t603;
t579 = (qJD(3) * t646 - t636) * pkin(9) + (t646 * t647 + qJDD(3)) * pkin(4) + t582;
t583 = 0.2e1 * qJD(4) * t646 + t671 * t601 + t674 * t603;
t641 = qJD(3) * pkin(4) - pkin(9) * t647;
t645 = t646 ^ 2;
t581 = -pkin(4) * t645 + pkin(9) * t635 - qJD(3) * t641 + t583;
t678 = sin(qJ(5));
t682 = cos(qJ(5));
t576 = t678 * t579 + t682 * t581;
t627 = t646 * t682 - t647 * t678;
t628 = t646 * t678 + t647 * t682;
t612 = -pkin(5) * t627 - pkin(10) * t628;
t668 = qJD(3) + qJD(5);
t666 = t668 ^ 2;
t667 = qJDD(3) + qJDD(5);
t574 = -pkin(5) * t666 + pkin(10) * t667 + t612 * t627 + t576;
t588 = -t635 * pkin(4) - t645 * pkin(9) + t647 * t641 + t613;
t596 = -qJD(5) * t628 + t635 * t682 - t636 * t678;
t597 = qJD(5) * t627 + t635 * t678 + t636 * t682;
t577 = (-t627 * t668 - t597) * pkin(10) + (t628 * t668 - t596) * pkin(5) + t588;
t677 = sin(qJ(6));
t681 = cos(qJ(6));
t571 = -t574 * t677 + t577 * t681;
t616 = -t628 * t677 + t668 * t681;
t586 = qJD(6) * t616 + t597 * t681 + t667 * t677;
t595 = qJDD(6) - t596;
t617 = t628 * t681 + t668 * t677;
t602 = -mrSges(7,1) * t616 + mrSges(7,2) * t617;
t623 = qJD(6) - t627;
t604 = -mrSges(7,2) * t623 + mrSges(7,3) * t616;
t569 = m(7) * t571 + mrSges(7,1) * t595 - mrSges(7,3) * t586 - t602 * t617 + t604 * t623;
t572 = t574 * t681 + t577 * t677;
t585 = -qJD(6) * t617 - t597 * t677 + t667 * t681;
t605 = mrSges(7,1) * t623 - mrSges(7,3) * t617;
t570 = m(7) * t572 - mrSges(7,2) * t595 + mrSges(7,3) * t585 + t602 * t616 - t605 * t623;
t561 = t681 * t569 + t677 * t570;
t619 = -mrSges(6,2) * t668 + mrSges(6,3) * t627;
t620 = mrSges(6,1) * t668 - mrSges(6,3) * t628;
t691 = m(6) * t588 - t596 * mrSges(6,1) + t597 * mrSges(6,2) - t627 * t619 + t628 * t620 + t561;
t687 = m(5) * t613 - t635 * mrSges(5,1) + mrSges(5,2) * t636 - t646 * t639 + t640 * t647 + t691;
t686 = -m(4) * t621 + t659 * mrSges(4,1) - mrSges(4,2) * t658 - t663 * t701 + t664 * t700 - t687;
t557 = m(3) * t629 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t685 + t686;
t705 = t557 * t684;
t611 = -mrSges(6,1) * t627 + mrSges(6,2) * t628;
t693 = -t569 * t677 + t681 * t570;
t560 = m(6) * t576 - mrSges(6,2) * t667 + mrSges(6,3) * t596 + t611 * t627 - t620 * t668 + t693;
t575 = t579 * t682 - t581 * t678;
t573 = -pkin(5) * t667 - pkin(10) * t666 + t612 * t628 - t575;
t689 = -m(7) * t573 + t585 * mrSges(7,1) - mrSges(7,2) * t586 + t616 * t604 - t605 * t617;
t565 = m(6) * t575 + mrSges(6,1) * t667 - mrSges(6,3) * t597 - t611 * t628 + t619 * t668 + t689;
t554 = t678 * t560 + t682 * t565;
t633 = -mrSges(5,1) * t646 + mrSges(5,2) * t647;
t552 = m(5) * t582 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t636 + qJD(3) * t639 - t633 * t647 + t554;
t694 = t682 * t560 - t565 * t678;
t553 = m(5) * t583 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t635 - qJD(3) * t640 + t633 * t646 + t694;
t546 = t674 * t552 + t671 * t553;
t657 = (-mrSges(4,1) * t683 + mrSges(4,2) * t679) * qJD(2);
t544 = m(4) * t614 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t658 + qJD(3) * t664 - t657 * t701 + t546;
t695 = -t552 * t671 + t674 * t553;
t545 = m(4) * t615 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t659 - qJD(3) * t663 + t657 * t700 + t695;
t696 = -t544 * t679 + t683 * t545;
t536 = m(3) * t630 - mrSges(3,1) * t685 - qJDD(2) * mrSges(3,2) + t696;
t539 = t683 * t544 + t679 * t545;
t538 = m(3) * t642 + t539;
t527 = t536 * t703 - t538 * t673 + t676 * t705;
t525 = m(2) * t660 + t527;
t531 = t684 * t536 - t557 * t680;
t530 = m(2) * t661 + t531;
t702 = t675 * t525 + t672 * t530;
t526 = t536 * t704 + t676 * t538 + t673 * t705;
t697 = -t525 * t672 + t675 * t530;
t589 = Ifges(7,5) * t617 + Ifges(7,6) * t616 + Ifges(7,3) * t623;
t591 = Ifges(7,1) * t617 + Ifges(7,4) * t616 + Ifges(7,5) * t623;
t562 = -mrSges(7,1) * t573 + mrSges(7,3) * t572 + Ifges(7,4) * t586 + Ifges(7,2) * t585 + Ifges(7,6) * t595 - t589 * t617 + t591 * t623;
t590 = Ifges(7,4) * t617 + Ifges(7,2) * t616 + Ifges(7,6) * t623;
t563 = mrSges(7,2) * t573 - mrSges(7,3) * t571 + Ifges(7,1) * t586 + Ifges(7,4) * t585 + Ifges(7,5) * t595 + t589 * t616 - t590 * t623;
t606 = Ifges(6,5) * t628 + Ifges(6,6) * t627 + Ifges(6,3) * t668;
t607 = Ifges(6,4) * t628 + Ifges(6,2) * t627 + Ifges(6,6) * t668;
t547 = mrSges(6,2) * t588 - mrSges(6,3) * t575 + Ifges(6,1) * t597 + Ifges(6,4) * t596 + Ifges(6,5) * t667 - pkin(10) * t561 - t562 * t677 + t563 * t681 + t606 * t627 - t607 * t668;
t608 = Ifges(6,1) * t628 + Ifges(6,4) * t627 + Ifges(6,5) * t668;
t548 = -mrSges(6,1) * t588 - mrSges(7,1) * t571 + mrSges(7,2) * t572 + mrSges(6,3) * t576 + Ifges(6,4) * t597 - Ifges(7,5) * t586 + Ifges(6,2) * t596 + Ifges(6,6) * t667 - Ifges(7,6) * t585 - Ifges(7,3) * t595 - pkin(5) * t561 - t590 * t617 + t591 * t616 - t606 * t628 + t608 * t668;
t624 = Ifges(5,5) * t647 + Ifges(5,6) * t646 + Ifges(5,3) * qJD(3);
t626 = Ifges(5,1) * t647 + Ifges(5,4) * t646 + Ifges(5,5) * qJD(3);
t532 = -mrSges(5,1) * t613 + mrSges(5,3) * t583 + Ifges(5,4) * t636 + Ifges(5,2) * t635 + Ifges(5,6) * qJDD(3) - pkin(4) * t691 + pkin(9) * t694 + qJD(3) * t626 + t678 * t547 + t682 * t548 - t647 * t624;
t625 = Ifges(5,4) * t647 + Ifges(5,2) * t646 + Ifges(5,6) * qJD(3);
t540 = mrSges(5,2) * t613 - mrSges(5,3) * t582 + Ifges(5,1) * t636 + Ifges(5,4) * t635 + Ifges(5,5) * qJDD(3) - pkin(9) * t554 - qJD(3) * t625 + t547 * t682 - t548 * t678 + t624 * t646;
t649 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t679 + Ifges(4,6) * t683) * qJD(2);
t651 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t679 + Ifges(4,4) * t683) * qJD(2);
t521 = -mrSges(4,1) * t621 + mrSges(4,3) * t615 + Ifges(4,4) * t658 + Ifges(4,2) * t659 + Ifges(4,6) * qJDD(3) - pkin(3) * t687 + qJ(4) * t695 + qJD(3) * t651 + t674 * t532 + t671 * t540 - t649 * t701;
t650 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t679 + Ifges(4,2) * t683) * qJD(2);
t522 = mrSges(4,2) * t621 - mrSges(4,3) * t614 + Ifges(4,1) * t658 + Ifges(4,4) * t659 + Ifges(4,5) * qJDD(3) - qJ(4) * t546 - qJD(3) * t650 - t532 * t671 + t540 * t674 + t649 * t700;
t520 = mrSges(3,2) * t642 - mrSges(3,3) * t629 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t685 - pkin(8) * t539 - t521 * t679 + t522 * t683;
t523 = Ifges(3,6) * qJDD(2) - mrSges(6,1) * t575 + mrSges(6,2) * t576 - mrSges(4,1) * t614 + mrSges(4,2) * t615 - pkin(4) * t554 - pkin(5) * t689 + t627 * t608 - Ifges(6,6) * t596 - Ifges(6,5) * t597 - pkin(3) * t546 - Ifges(6,3) * t667 - mrSges(3,1) * t642 + t646 * t626 - t647 * t625 - t628 * t607 + mrSges(3,3) * t630 - Ifges(5,6) * t635 - Ifges(5,5) * t636 + t685 * Ifges(3,5) - mrSges(5,1) * t582 + mrSges(5,2) * t583 - pkin(2) * t539 - pkin(10) * t693 - t677 * t563 - t681 * t562 + (-t650 * t679 + t651 * t683) * qJD(2) + (-Ifges(4,3) - Ifges(5,3)) * qJDD(3) - Ifges(4,5) * t658 - Ifges(4,6) * t659;
t690 = pkin(7) * t531 + t520 * t680 + t523 * t684;
t519 = mrSges(3,1) * t629 - mrSges(3,2) * t630 + Ifges(3,3) * qJDD(2) + pkin(2) * t686 + pkin(8) * t696 + t683 * t521 + t679 * t522;
t518 = mrSges(2,2) * t670 - mrSges(2,3) * t660 + t684 * t520 - t680 * t523 + (-t526 * t673 - t527 * t676) * pkin(7);
t517 = -mrSges(2,1) * t670 + mrSges(2,3) * t661 - pkin(1) * t526 - t673 * t519 + t676 * t690;
t1 = [-m(1) * g(1) + t697; -m(1) * g(2) + t702; -m(1) * g(3) + m(2) * t670 + t526; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t702 - t672 * t517 + t675 * t518; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t697 + t675 * t517 + t672 * t518; -mrSges(1,1) * g(2) + mrSges(2,1) * t660 + mrSges(1,2) * g(1) - mrSges(2,2) * t661 + pkin(1) * t527 + t676 * t519 + t673 * t690;];
tauB  = t1;
