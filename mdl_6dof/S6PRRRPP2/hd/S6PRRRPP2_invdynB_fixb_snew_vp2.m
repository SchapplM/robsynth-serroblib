% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-05-05 06:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:42:26
% EndTime: 2019-05-05 06:42:33
% DurationCPUTime: 5.39s
% Computational Cost: add. (61958->296), mult. (117114->356), div. (0->0), fcn. (76024->10), ass. (0->121)
t707 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t689 = Ifges(5,4) - Ifges(6,5) - Ifges(7,4);
t706 = Ifges(5,5) + Ifges(6,4) - Ifges(7,5);
t705 = -Ifges(5,2) - Ifges(6,3) - Ifges(7,2);
t687 = -Ifges(6,6) + Ifges(7,6) + Ifges(5,6);
t704 = Ifges(7,3) + Ifges(5,3) + Ifges(6,2);
t662 = sin(qJ(4));
t663 = sin(qJ(3));
t692 = qJD(2) * t663;
t701 = cos(qJ(4));
t639 = -qJD(3) * t701 + t662 * t692;
t665 = cos(qJ(3));
t690 = qJD(2) * qJD(3);
t682 = t665 * t690;
t643 = t663 * qJDD(2) + t682;
t608 = -t639 * qJD(4) + t662 * qJDD(3) + t643 * t701;
t658 = sin(pkin(10));
t660 = cos(pkin(10));
t645 = t658 * g(1) - t660 * g(2);
t646 = -t660 * g(1) - t658 * g(2);
t657 = -g(3) + qJDD(1);
t666 = cos(qJ(2));
t661 = cos(pkin(6));
t664 = sin(qJ(2));
t695 = t661 * t664;
t659 = sin(pkin(6));
t696 = t659 * t664;
t587 = t645 * t695 + t666 * t646 + t657 * t696;
t668 = qJD(2) ^ 2;
t581 = -t668 * pkin(2) + qJDD(2) * pkin(8) + t587;
t624 = -t659 * t645 + t661 * t657;
t576 = -t663 * t581 + t665 * t624;
t642 = (-pkin(3) * t665 - pkin(9) * t663) * qJD(2);
t667 = qJD(3) ^ 2;
t672 = qJDD(3) * pkin(3) + t667 * pkin(9) - t642 * t692 + t576;
t691 = t665 * qJD(2);
t652 = -qJD(4) + t691;
t698 = t639 * t652;
t703 = -(t608 + t698) * qJ(5) - t672;
t586 = -t664 * t646 + (t645 * t661 + t657 * t659) * t666;
t702 = -2 * qJD(5);
t700 = -mrSges(5,3) - mrSges(6,2);
t577 = t665 * t581 + t663 * t624;
t573 = -t667 * pkin(3) + qJDD(3) * pkin(9) + t642 * t691 + t577;
t580 = -qJDD(2) * pkin(2) - t668 * pkin(8) - t586;
t654 = t663 * t690;
t644 = t665 * qJDD(2) - t654;
t575 = (-t643 - t682) * pkin(9) + (-t644 + t654) * pkin(3) + t580;
t568 = -t662 * t573 + t575 * t701;
t617 = -t652 * mrSges(7,2) + t639 * mrSges(7,3);
t618 = t652 * mrSges(5,2) - t639 * mrSges(5,3);
t636 = -qJDD(4) + t644;
t640 = t662 * qJD(3) + t692 * t701;
t611 = t639 * pkin(4) - t640 * qJ(5);
t651 = t652 ^ 2;
t566 = t636 * pkin(4) - t651 * qJ(5) + t640 * t611 + qJDD(5) - t568;
t623 = -t639 * mrSges(6,2) - t652 * mrSges(6,3);
t559 = -0.2e1 * qJD(6) * t640 + (-t608 + t698) * qJ(6) + (t639 * t640 + t636) * pkin(5) + t566;
t613 = -t639 * mrSges(7,1) + t640 * mrSges(7,2);
t677 = -m(7) * t559 + t608 * mrSges(7,3) + t640 * t613;
t671 = -m(6) * t566 - t636 * mrSges(6,1) - t652 * t623 + t677;
t612 = t639 * mrSges(6,1) - t640 * mrSges(6,3);
t693 = -t639 * mrSges(5,1) - t640 * mrSges(5,2) - t612;
t555 = m(5) * t568 + (-t617 - t618) * t652 + t693 * t640 + (-mrSges(5,1) - mrSges(7,1)) * t636 + t700 * t608 + t671;
t569 = t701 * t573 + t662 * t575;
t607 = t640 * qJD(4) - qJDD(3) * t701 + t662 * t643;
t620 = t652 * mrSges(7,1) - t640 * mrSges(7,3);
t621 = -t652 * mrSges(5,1) - t640 * mrSges(5,3);
t565 = -t651 * pkin(4) - t636 * qJ(5) - t639 * t611 + t652 * t702 + t569;
t622 = t652 * mrSges(6,1) + t640 * mrSges(6,2);
t619 = t652 * pkin(5) - t640 * qJ(6);
t635 = t639 ^ 2;
t561 = -t635 * pkin(5) + t607 * qJ(6) + 0.2e1 * qJD(6) * t639 - t652 * t619 + t565;
t686 = m(7) * t561 + t607 * mrSges(7,3) + t639 * t613;
t674 = m(6) * t565 - t636 * mrSges(6,3) - t652 * t622 + t686;
t556 = m(5) * t569 + (-t620 + t621) * t652 + t693 * t639 + (mrSges(5,2) - mrSges(7,2)) * t636 + t700 * t607 + t674;
t551 = t555 * t701 + t662 * t556;
t647 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t692;
t648 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t691;
t670 = -m(4) * t580 + t644 * mrSges(4,1) - t643 * mrSges(4,2) - t647 * t692 + t648 * t691 - t551;
t547 = m(3) * t586 + qJDD(2) * mrSges(3,1) - t668 * mrSges(3,2) + t670;
t699 = t547 * t666;
t697 = t652 * t617;
t641 = (-mrSges(4,1) * t665 + mrSges(4,2) * t663) * qJD(2);
t679 = -t662 * t555 + t701 * t556;
t550 = m(4) * t577 - qJDD(3) * mrSges(4,2) + t644 * mrSges(4,3) - qJD(3) * t647 + t641 * t691 + t679;
t567 = t640 * t702 + (-t640 * t652 + t607) * pkin(4) + t703;
t563 = -t635 * qJ(6) + qJDD(6) + (-pkin(4) - pkin(5)) * t607 + (pkin(4) * t652 + (2 * qJD(5)) + t619) * t640 - t703;
t676 = -m(7) * t563 + t607 * mrSges(7,1) - t608 * mrSges(7,2) + t639 * t617 - t640 * t620;
t557 = m(6) * t567 + t607 * mrSges(6,1) - t608 * mrSges(6,3) - t640 * t622 + t639 * t623 + t676;
t669 = m(5) * t672 - t607 * mrSges(5,1) - t608 * mrSges(5,2) - t639 * t618 - t640 * t621 - t557;
t553 = m(4) * t576 + qJDD(3) * mrSges(4,1) - t643 * mrSges(4,3) + qJD(3) * t648 - t641 * t692 + t669;
t680 = t665 * t550 - t663 * t553;
t540 = m(3) * t587 - t668 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t680;
t543 = t663 * t550 + t665 * t553;
t542 = m(3) * t624 + t543;
t529 = t540 * t695 - t659 * t542 + t661 * t699;
t527 = m(2) * t645 + t529;
t535 = t666 * t540 - t664 * t547;
t534 = m(2) * t646 + t535;
t694 = t660 * t527 + t658 * t534;
t528 = t540 * t696 + t661 * t542 + t659 * t699;
t685 = t687 * t639 - t706 * t640 + t704 * t652;
t684 = t705 * t639 + t689 * t640 - t687 * t652;
t683 = t689 * t639 - t707 * t640 + t706 * t652;
t681 = -t658 * t527 + t660 * t534;
t536 = mrSges(5,1) * t672 + mrSges(5,3) * t569 - mrSges(6,1) * t567 + mrSges(6,2) * t565 + mrSges(7,1) * t563 - mrSges(7,3) * t561 - pkin(5) * t676 - qJ(6) * t686 - pkin(4) * t557 + (qJ(6) * t620 + t683) * t652 + t685 * t640 + (qJ(6) * mrSges(7,2) - t687) * t636 + t689 * t608 + t705 * t607;
t558 = t636 * mrSges(7,1) - t677 + t697;
t544 = -mrSges(5,2) * t672 + mrSges(6,2) * t566 + mrSges(7,2) * t563 - mrSges(5,3) * t568 - mrSges(6,3) * t567 - mrSges(7,3) * t559 - qJ(5) * t557 - qJ(6) * t558 - t689 * t607 + t707 * t608 - t636 * t706 + t685 * t639 + t684 * t652;
t627 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t663 + Ifges(4,6) * t665) * qJD(2);
t628 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t663 + Ifges(4,2) * t665) * qJD(2);
t530 = mrSges(4,2) * t580 - mrSges(4,3) * t576 + Ifges(4,1) * t643 + Ifges(4,4) * t644 + Ifges(4,5) * qJDD(3) - pkin(9) * t551 - qJD(3) * t628 - t662 * t536 + t544 * t701 + t627 * t691;
t629 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t663 + Ifges(4,4) * t665) * qJD(2);
t531 = Ifges(4,4) * t643 + Ifges(4,2) * t644 + mrSges(7,1) * t559 - mrSges(7,2) * t561 - mrSges(6,3) * t565 - pkin(3) * t551 + pkin(5) * t558 + mrSges(4,3) * t577 - mrSges(4,1) * t580 + qJD(3) * t629 + mrSges(6,1) * t566 - mrSges(5,1) * t568 + mrSges(5,2) * t569 - t627 * t692 - pkin(4) * (t671 - t697) - qJ(5) * (-t652 * t620 + t674) + Ifges(4,6) * qJDD(3) + (pkin(4) * t612 - t684) * t640 + (qJ(5) * t612 + t683) * t639 + (pkin(4) * mrSges(7,1) + qJ(5) * mrSges(7,2) + t704) * t636 + (pkin(4) * mrSges(6,2) - t706) * t608 + (qJ(5) * mrSges(6,2) + t687) * t607;
t524 = mrSges(3,2) * t624 - mrSges(3,3) * t586 + Ifges(3,5) * qJDD(2) - t668 * Ifges(3,6) - pkin(8) * t543 + t665 * t530 - t663 * t531;
t525 = Ifges(3,6) * qJDD(2) + t668 * Ifges(3,5) - mrSges(3,1) * t624 + mrSges(3,3) * t587 - Ifges(4,5) * t643 - Ifges(4,6) * t644 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t576 + mrSges(4,2) * t577 - t662 * t544 - t701 * t536 - pkin(3) * t669 - pkin(9) * t679 - pkin(2) * t543 + (-t663 * t628 + t665 * t629) * qJD(2);
t673 = pkin(7) * t535 + t524 * t664 + t525 * t666;
t523 = mrSges(3,1) * t586 - mrSges(3,2) * t587 + Ifges(3,3) * qJDD(2) + pkin(2) * t670 + pkin(8) * t680 + t663 * t530 + t665 * t531;
t522 = mrSges(2,2) * t657 - mrSges(2,3) * t645 + t666 * t524 - t664 * t525 + (-t528 * t659 - t529 * t661) * pkin(7);
t521 = -mrSges(2,1) * t657 + mrSges(2,3) * t646 - pkin(1) * t528 - t659 * t523 + t661 * t673;
t1 = [-m(1) * g(1) + t681; -m(1) * g(2) + t694; -m(1) * g(3) + m(2) * t657 + t528; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t694 - t658 * t521 + t660 * t522; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t681 + t660 * t521 + t658 * t522; -mrSges(1,1) * g(2) + mrSges(2,1) * t645 + mrSges(1,2) * g(1) - mrSges(2,2) * t646 + pkin(1) * t529 + t661 * t523 + t659 * t673;];
tauB  = t1;
