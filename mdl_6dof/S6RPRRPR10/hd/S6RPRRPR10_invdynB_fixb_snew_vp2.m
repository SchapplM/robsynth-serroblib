% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-05-05 23:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:53:27
% EndTime: 2019-05-05 23:53:33
% DurationCPUTime: 3.69s
% Computational Cost: add. (39200->317), mult. (74862->370), div. (0->0), fcn. (45254->8), ass. (0->126)
t703 = Ifges(5,1) + Ifges(6,1);
t692 = Ifges(5,4) - Ifges(6,5);
t690 = Ifges(5,5) + Ifges(6,4);
t702 = Ifges(5,2) + Ifges(6,3);
t689 = Ifges(5,6) - Ifges(6,6);
t701 = -Ifges(5,3) - Ifges(6,2);
t656 = sin(qJ(1));
t659 = cos(qJ(1));
t640 = -t659 * g(1) - t656 * g(2);
t700 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t640;
t654 = sin(qJ(4));
t658 = cos(qJ(3));
t682 = qJD(1) * t658;
t696 = cos(qJ(4));
t630 = -t696 * qJD(3) + t654 * t682;
t655 = sin(qJ(3));
t680 = qJD(1) * qJD(3);
t678 = t655 * t680;
t635 = t658 * qJDD(1) - t678;
t595 = -t630 * qJD(4) + t654 * qJDD(3) + t696 * t635;
t639 = t656 * g(1) - t659 * g(2);
t661 = qJD(1) ^ 2;
t668 = -t661 * qJ(2) + qJDD(2) - t639;
t697 = -pkin(1) - pkin(7);
t614 = t697 * qJDD(1) + t668;
t603 = t655 * g(3) + t658 * t614;
t633 = (pkin(3) * t655 - pkin(8) * t658) * qJD(1);
t660 = qJD(3) ^ 2;
t667 = qJDD(3) * pkin(3) + t660 * pkin(8) - t633 * t682 + t603;
t681 = t655 * qJD(1);
t643 = qJD(4) + t681;
t688 = t630 * t643;
t699 = (-t595 + t688) * qJ(5) - t667;
t698 = 2 * qJD(5);
t695 = mrSges(2,1) - mrSges(3,2);
t694 = -mrSges(5,3) - mrSges(6,2);
t693 = -Ifges(3,4) + Ifges(2,5);
t691 = (Ifges(3,5) - Ifges(2,6));
t604 = -t658 * g(3) + t655 * t614;
t632 = (mrSges(4,1) * t655 + mrSges(4,2) * t658) * qJD(1);
t677 = t658 * t680;
t634 = -t655 * qJDD(1) - t677;
t637 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t682;
t613 = t697 * t661 - t700;
t575 = (-t635 + t678) * pkin(8) + (-t634 + t677) * pkin(3) + t613;
t579 = -t660 * pkin(3) + qJDD(3) * pkin(8) - t633 * t681 + t604;
t563 = t654 * t575 + t696 * t579;
t631 = t654 * qJD(3) + t696 * t682;
t594 = t631 * qJD(4) - t696 * qJDD(3) + t654 * t635;
t606 = t643 * mrSges(5,1) - t631 * mrSges(5,3);
t629 = qJDD(4) - t634;
t600 = t630 * pkin(4) - t631 * qJ(5);
t642 = t643 ^ 2;
t559 = -t642 * pkin(4) + t629 * qJ(5) - t630 * t600 + t643 * t698 + t563;
t607 = -t643 * mrSges(6,1) + t631 * mrSges(6,2);
t562 = t696 * t575 - t654 * t579;
t560 = -t629 * pkin(4) - t642 * qJ(5) + t631 * t600 + qJDD(5) - t562;
t554 = (-t595 - t688) * pkin(9) + (t630 * t631 - t629) * pkin(5) + t560;
t612 = -t643 * pkin(5) - t631 * pkin(9);
t628 = t630 ^ 2;
t555 = -t628 * pkin(5) + t594 * pkin(9) + t643 * t612 + t559;
t653 = sin(qJ(6));
t657 = cos(qJ(6));
t552 = t657 * t554 - t653 * t555;
t596 = t657 * t630 - t653 * t631;
t567 = t596 * qJD(6) + t653 * t594 + t657 * t595;
t597 = t653 * t630 + t657 * t631;
t573 = -t596 * mrSges(7,1) + t597 * mrSges(7,2);
t641 = qJD(6) - t643;
t580 = -t641 * mrSges(7,2) + t596 * mrSges(7,3);
t624 = qJDD(6) - t629;
t550 = m(7) * t552 + t624 * mrSges(7,1) - t567 * mrSges(7,3) - t597 * t573 + t641 * t580;
t553 = t653 * t554 + t657 * t555;
t566 = -t597 * qJD(6) + t657 * t594 - t653 * t595;
t581 = t641 * mrSges(7,1) - t597 * mrSges(7,3);
t551 = m(7) * t553 - t624 * mrSges(7,2) + t566 * mrSges(7,3) + t596 * t573 - t641 * t581;
t673 = -t653 * t550 + t657 * t551;
t669 = m(6) * t559 + t629 * mrSges(6,3) + t643 * t607 + t673;
t601 = t630 * mrSges(6,1) - t631 * mrSges(6,3);
t683 = -t630 * mrSges(5,1) - t631 * mrSges(5,2) - t601;
t541 = m(5) * t563 - t629 * mrSges(5,2) + t694 * t594 - t643 * t606 + t683 * t630 + t669;
t605 = -t643 * mrSges(5,2) - t630 * mrSges(5,3);
t544 = t657 * t550 + t653 * t551;
t608 = -t630 * mrSges(6,2) + t643 * mrSges(6,3);
t664 = -m(6) * t560 + t629 * mrSges(6,1) + t643 * t608 - t544;
t543 = m(5) * t562 + t629 * mrSges(5,1) + t694 * t595 + t643 * t605 + t683 * t631 + t664;
t674 = t696 * t541 - t654 * t543;
t537 = m(4) * t604 - qJDD(3) * mrSges(4,2) + t634 * mrSges(4,3) - qJD(3) * t637 - t632 * t681 + t674;
t636 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t681;
t561 = -0.2e1 * qJD(5) * t631 + (t631 * t643 + t594) * pkin(4) + t699;
t557 = -t628 * pkin(9) + (-pkin(4) - pkin(5)) * t594 + (-pkin(4) * t643 + t612 + t698) * t631 - t699;
t671 = -m(7) * t557 + t566 * mrSges(7,1) - t567 * mrSges(7,2) + t596 * t580 - t597 * t581;
t548 = m(6) * t561 + t594 * mrSges(6,1) - t595 * mrSges(6,3) - t631 * t607 + t630 * t608 + t671;
t662 = m(5) * t667 - t594 * mrSges(5,1) - t595 * mrSges(5,2) - t630 * t605 - t631 * t606 - t548;
t547 = m(4) * t603 + qJDD(3) * mrSges(4,1) - t635 * mrSges(4,3) + qJD(3) * t636 - t632 * t682 + t662;
t532 = t655 * t537 + t658 * t547;
t616 = -qJDD(1) * pkin(1) + t668;
t666 = -m(3) * t616 + (t661 * mrSges(3,3)) - t532;
t530 = m(2) * t639 - (t661 * mrSges(2,2)) + t695 * qJDD(1) + t666;
t615 = t661 * pkin(1) + t700;
t538 = t654 * t541 + t696 * t543;
t665 = -m(4) * t613 + t634 * mrSges(4,1) - t635 * mrSges(4,2) - t636 * t681 - t637 * t682 - t538;
t663 = -m(3) * t615 + (t661 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t665;
t535 = m(2) * t640 - (t661 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t663;
t687 = t659 * t530 + t656 * t535;
t686 = t702 * t630 - t692 * t631 - t689 * t643;
t685 = t689 * t630 - t690 * t631 + t701 * t643;
t684 = -t692 * t630 + t703 * t631 + t690 * t643;
t676 = -t656 * t530 + t659 * t535;
t675 = t658 * t537 - t655 * t547;
t620 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t658 - Ifges(4,4) * t655) * qJD(1);
t619 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t658 - Ifges(4,2) * t655) * qJD(1);
t618 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t658 - Ifges(4,6) * t655) * qJD(1);
t570 = Ifges(7,1) * t597 + Ifges(7,4) * t596 + Ifges(7,5) * t641;
t569 = Ifges(7,4) * t597 + Ifges(7,2) * t596 + Ifges(7,6) * t641;
t568 = Ifges(7,5) * t597 + Ifges(7,6) * t596 + Ifges(7,3) * t641;
t546 = mrSges(7,2) * t557 - mrSges(7,3) * t552 + Ifges(7,1) * t567 + Ifges(7,4) * t566 + Ifges(7,5) * t624 + t596 * t568 - t641 * t569;
t545 = -mrSges(7,1) * t557 + mrSges(7,3) * t553 + Ifges(7,4) * t567 + Ifges(7,2) * t566 + Ifges(7,6) * t624 - t597 * t568 + t641 * t570;
t531 = -m(3) * g(3) + t675;
t528 = -mrSges(5,2) * t667 + mrSges(6,2) * t560 - mrSges(5,3) * t562 - mrSges(6,3) * t561 - pkin(9) * t544 - qJ(5) * t548 - t653 * t545 + t657 * t546 - t692 * t594 + t703 * t595 + t690 * t629 + t685 * t630 + t686 * t643;
t527 = mrSges(5,1) * t667 - mrSges(6,1) * t561 + mrSges(6,2) * t559 + mrSges(5,3) * t563 - pkin(4) * t548 - pkin(5) * t671 - pkin(9) * t673 - t657 * t545 - t653 * t546 - t702 * t594 + t692 * t595 + t689 * t629 + t685 * t631 + t684 * t643;
t526 = (qJ(5) * mrSges(6,2) + t689) * t594 + (pkin(4) * mrSges(6,2) - t690) * t595 - t618 * t682 + (qJ(5) * t601 - t684) * t630 + (pkin(4) * t601 + t686) * t631 + Ifges(4,2) * t634 + Ifges(4,4) * t635 - qJ(5) * t669 + qJD(3) * t620 + Ifges(7,3) * t624 - mrSges(4,1) * t613 - t596 * t570 - pkin(4) * t664 + t597 * t569 + mrSges(4,3) * t604 + Ifges(7,5) * t567 + mrSges(6,1) * t560 - mrSges(5,1) * t562 + mrSges(5,2) * t563 + Ifges(7,6) * t566 - mrSges(7,2) * t553 - mrSges(6,3) * t559 + Ifges(4,6) * qJDD(3) + mrSges(7,1) * t552 + pkin(5) * t544 - pkin(3) * t538 + t701 * t629;
t525 = mrSges(4,2) * t613 - mrSges(4,3) * t603 + Ifges(4,1) * t635 + Ifges(4,4) * t634 + Ifges(4,5) * qJDD(3) - pkin(8) * t538 - qJD(3) * t619 - t654 * t527 + t696 * t528 - t618 * t681;
t524 = -qJ(2) * t531 - mrSges(2,3) * t639 + pkin(2) * t532 + mrSges(3,1) * t616 + pkin(8) * t674 + t654 * t528 + t696 * t527 + pkin(3) * t662 + Ifges(4,5) * t635 + Ifges(4,6) * t634 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t603 - mrSges(4,2) * t604 + (t691 * t661) + t693 * qJDD(1) + (t658 * t619 + t655 * t620) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t523 = -mrSges(3,1) * t615 + mrSges(2,3) * t640 - pkin(1) * t531 - pkin(2) * t665 - pkin(7) * t675 + t695 * g(3) - t691 * qJDD(1) - t655 * t525 - t658 * t526 + t693 * t661;
t1 = [-m(1) * g(1) + t676; -m(1) * g(2) + t687; (-m(1) - m(2) - m(3)) * g(3) + t675; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t687 - t656 * t523 + t659 * t524; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t676 + t659 * t523 + t656 * t524; qJ(2) * t663 + pkin(1) * t666 - pkin(7) * t532 + mrSges(2,1) * t639 - mrSges(2,2) * t640 + t658 * t525 - t655 * t526 + mrSges(3,2) * t616 - mrSges(3,3) * t615 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
