% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRPRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:52:23
% EndTime: 2019-05-05 03:52:35
% DurationCPUTime: 10.03s
% Computational Cost: add. (156046->316), mult. (320431->396), div. (0->0), fcn. (224770->12), ass. (0->131)
t691 = Ifges(6,1) + Ifges(7,1);
t685 = Ifges(6,4) - Ifges(7,5);
t690 = -Ifges(6,5) - Ifges(7,4);
t689 = Ifges(6,2) + Ifges(7,3);
t683 = Ifges(6,6) - Ifges(7,6);
t688 = -Ifges(6,3) - Ifges(7,2);
t648 = sin(pkin(10));
t651 = cos(pkin(10));
t638 = t648 * g(1) - t651 * g(2);
t639 = -t651 * g(1) - t648 * g(2);
t646 = -g(3) + qJDD(1);
t649 = sin(pkin(6));
t652 = cos(pkin(6));
t655 = sin(qJ(2));
t657 = cos(qJ(2));
t597 = -t655 * t639 + (t638 * t652 + t646 * t649) * t657;
t687 = cos(qJ(5));
t686 = -mrSges(6,3) - mrSges(7,2);
t659 = qJD(2) ^ 2;
t680 = t652 * t655;
t681 = t649 * t655;
t598 = t638 * t680 + t657 * t639 + t646 * t681;
t592 = -t659 * pkin(2) + qJDD(2) * pkin(8) + t598;
t615 = -t649 * t638 + t652 * t646;
t654 = sin(qJ(3));
t656 = cos(qJ(3));
t583 = t656 * t592 + t654 * t615;
t634 = (-pkin(3) * t656 - qJ(4) * t654) * qJD(2);
t658 = qJD(3) ^ 2;
t673 = t656 * qJD(2);
t564 = -t658 * pkin(3) + qJDD(3) * qJ(4) + t634 * t673 + t583;
t591 = -qJDD(2) * pkin(2) - t659 * pkin(8) - t597;
t672 = qJD(2) * qJD(3);
t670 = t656 * t672;
t636 = t654 * qJDD(2) + t670;
t645 = t654 * t672;
t637 = t656 * qJDD(2) - t645;
t570 = (-t636 - t670) * qJ(4) + (-t637 + t645) * pkin(3) + t591;
t647 = sin(pkin(11));
t650 = cos(pkin(11));
t674 = qJD(2) * t654;
t629 = t647 * qJD(3) + t650 * t674;
t559 = -0.2e1 * qJD(4) * t629 - t647 * t564 + t650 * t570;
t613 = t647 * qJDD(3) + t650 * t636;
t628 = t650 * qJD(3) - t647 * t674;
t556 = (-t628 * t673 - t613) * pkin(9) + (t628 * t629 - t637) * pkin(4) + t559;
t560 = 0.2e1 * qJD(4) * t628 + t650 * t564 + t647 * t570;
t612 = t650 * qJDD(3) - t647 * t636;
t614 = -pkin(4) * t673 - t629 * pkin(9);
t627 = t628 ^ 2;
t558 = -t627 * pkin(4) + t612 * pkin(9) + t614 * t673 + t560;
t653 = sin(qJ(5));
t552 = t653 * t556 + t687 * t558;
t605 = t653 * t628 + t687 * t629;
t571 = t605 * qJD(5) - t687 * t612 + t653 * t613;
t644 = qJD(5) - t673;
t594 = t644 * mrSges(6,1) - t605 * mrSges(6,3);
t604 = -t687 * t628 + t653 * t629;
t631 = qJDD(5) - t637;
t584 = t604 * pkin(5) - t605 * qJ(6);
t643 = t644 ^ 2;
t549 = -t643 * pkin(5) + t631 * qJ(6) + 0.2e1 * qJD(6) * t644 - t604 * t584 + t552;
t595 = -t644 * mrSges(7,1) + t605 * mrSges(7,2);
t671 = m(7) * t549 + t631 * mrSges(7,3) + t644 * t595;
t585 = t604 * mrSges(7,1) - t605 * mrSges(7,3);
t675 = -t604 * mrSges(6,1) - t605 * mrSges(6,2) - t585;
t542 = m(6) * t552 - t631 * mrSges(6,2) + t686 * t571 - t644 * t594 + t675 * t604 + t671;
t551 = t687 * t556 - t653 * t558;
t572 = -t604 * qJD(5) + t653 * t612 + t687 * t613;
t593 = -t644 * mrSges(6,2) - t604 * mrSges(6,3);
t550 = -t631 * pkin(5) - t643 * qJ(6) + t605 * t584 + qJDD(6) - t551;
t596 = -t604 * mrSges(7,2) + t644 * mrSges(7,3);
t665 = -m(7) * t550 + t631 * mrSges(7,1) + t644 * t596;
t544 = m(6) * t551 + t631 * mrSges(6,1) + t686 * t572 + t644 * t593 + t675 * t605 + t665;
t537 = t653 * t542 + t687 * t544;
t606 = -t628 * mrSges(5,1) + t629 * mrSges(5,2);
t610 = mrSges(5,2) * t673 + t628 * mrSges(5,3);
t535 = m(5) * t559 - t637 * mrSges(5,1) - t613 * mrSges(5,3) - t629 * t606 - t610 * t673 + t537;
t611 = -mrSges(5,1) * t673 - t629 * mrSges(5,3);
t666 = t687 * t542 - t653 * t544;
t536 = m(5) * t560 + t637 * mrSges(5,2) + t612 * mrSges(5,3) + t628 * t606 + t611 * t673 + t666;
t533 = t650 * t535 + t647 * t536;
t640 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t674;
t641 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t673;
t661 = -m(4) * t591 + t637 * mrSges(4,1) - t636 * mrSges(4,2) - t640 * t674 + t641 * t673 - t533;
t529 = m(3) * t597 + qJDD(2) * mrSges(3,1) - t659 * mrSges(3,2) + t661;
t682 = t529 * t657;
t635 = (-mrSges(4,1) * t656 + mrSges(4,2) * t654) * qJD(2);
t667 = -t647 * t535 + t650 * t536;
t532 = m(4) * t583 - qJDD(3) * mrSges(4,2) + t637 * mrSges(4,3) - qJD(3) * t640 + t635 * t673 + t667;
t582 = -t654 * t592 + t656 * t615;
t563 = -qJDD(3) * pkin(3) - t658 * qJ(4) + t634 * t674 + qJDD(4) - t582;
t561 = -t612 * pkin(4) - t627 * pkin(9) + t629 * t614 + t563;
t554 = -0.2e1 * qJD(6) * t605 + (t604 * t644 - t572) * qJ(6) + (t605 * t644 + t571) * pkin(5) + t561;
t547 = m(7) * t554 + t571 * mrSges(7,1) - t572 * mrSges(7,3) - t605 * t595 + t604 * t596;
t662 = m(6) * t561 + t571 * mrSges(6,1) + t572 * mrSges(6,2) + t604 * t593 + t605 * t594 + t547;
t660 = -m(5) * t563 + t612 * mrSges(5,1) - t613 * mrSges(5,2) + t628 * t610 - t629 * t611 - t662;
t546 = m(4) * t582 + qJDD(3) * mrSges(4,1) - t636 * mrSges(4,3) + qJD(3) * t641 - t635 * t674 + t660;
t668 = t656 * t532 - t654 * t546;
t523 = m(3) * t598 - t659 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t668;
t526 = t654 * t532 + t656 * t546;
t525 = m(3) * t615 + t526;
t512 = t523 * t680 - t649 * t525 + t652 * t682;
t510 = m(2) * t638 + t512;
t516 = t657 * t523 - t655 * t529;
t515 = m(2) * t639 + t516;
t679 = t651 * t510 + t648 * t515;
t678 = t689 * t604 - t685 * t605 - t683 * t644;
t677 = t683 * t604 + t690 * t605 + t688 * t644;
t676 = -t685 * t604 + t691 * t605 - t690 * t644;
t511 = t523 * t681 + t652 * t525 + t649 * t682;
t669 = -t648 * t510 + t651 * t515;
t538 = -mrSges(6,1) * t561 - mrSges(7,1) * t554 + mrSges(7,2) * t549 + mrSges(6,3) * t552 - pkin(5) * t547 - t689 * t571 + t685 * t572 + t677 * t605 + t683 * t631 + t676 * t644;
t539 = mrSges(6,2) * t561 + mrSges(7,2) * t550 - mrSges(6,3) * t551 - mrSges(7,3) * t554 - qJ(6) * t547 - t685 * t571 + t691 * t572 + t677 * t604 - t631 * t690 + t678 * t644;
t599 = Ifges(5,5) * t629 + Ifges(5,6) * t628 - Ifges(5,3) * t673;
t601 = Ifges(5,1) * t629 + Ifges(5,4) * t628 - Ifges(5,5) * t673;
t518 = -mrSges(5,1) * t563 + mrSges(5,3) * t560 + Ifges(5,4) * t613 + Ifges(5,2) * t612 - Ifges(5,6) * t637 - pkin(4) * t662 + pkin(9) * t666 + t687 * t538 + t653 * t539 - t629 * t599 - t601 * t673;
t600 = Ifges(5,4) * t629 + Ifges(5,2) * t628 - Ifges(5,6) * t673;
t519 = mrSges(5,2) * t563 - mrSges(5,3) * t559 + Ifges(5,1) * t613 + Ifges(5,4) * t612 - Ifges(5,5) * t637 - pkin(9) * t537 - t653 * t538 + t687 * t539 + t628 * t599 + t600 * t673;
t621 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t654 + Ifges(4,6) * t656) * qJD(2);
t622 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t654 + Ifges(4,2) * t656) * qJD(2);
t508 = mrSges(4,2) * t591 - mrSges(4,3) * t582 + Ifges(4,1) * t636 + Ifges(4,4) * t637 + Ifges(4,5) * qJDD(3) - qJ(4) * t533 - qJD(3) * t622 - t647 * t518 + t650 * t519 + t621 * t673;
t623 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t654 + Ifges(4,4) * t656) * qJD(2);
t517 = t628 * t601 - t629 * t600 + Ifges(4,4) * t636 - Ifges(5,6) * t612 - Ifges(5,5) * t613 + qJD(3) * t623 + mrSges(4,3) * t583 - mrSges(4,1) * t591 - mrSges(5,1) * t559 + mrSges(5,2) * t560 - mrSges(6,1) * t551 + mrSges(6,2) * t552 - mrSges(7,3) * t549 + mrSges(7,1) * t550 - pkin(4) * t537 - pkin(3) * t533 - t621 * t674 + (qJ(6) * t585 - t676) * t604 + Ifges(4,6) * qJDD(3) - pkin(5) * t665 - qJ(6) * t671 + (pkin(5) * t585 + t678) * t605 + (pkin(5) * mrSges(7,2) + t690) * t572 + (qJ(6) * mrSges(7,2) + t683) * t571 + t688 * t631 + (Ifges(5,3) + Ifges(4,2)) * t637;
t506 = mrSges(3,2) * t615 - mrSges(3,3) * t597 + Ifges(3,5) * qJDD(2) - t659 * Ifges(3,6) - pkin(8) * t526 + t656 * t508 - t654 * t517;
t507 = Ifges(3,6) * qJDD(2) + t659 * Ifges(3,5) - mrSges(3,1) * t615 + mrSges(3,3) * t598 - Ifges(4,5) * t636 - Ifges(4,6) * t637 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t582 + mrSges(4,2) * t583 - t647 * t519 - t650 * t518 - pkin(3) * t660 - qJ(4) * t667 - pkin(2) * t526 + (-t654 * t622 + t656 * t623) * qJD(2);
t663 = pkin(7) * t516 + t506 * t655 + t507 * t657;
t505 = mrSges(3,1) * t597 - mrSges(3,2) * t598 + Ifges(3,3) * qJDD(2) + pkin(2) * t661 + pkin(8) * t668 + t654 * t508 + t656 * t517;
t504 = mrSges(2,2) * t646 - mrSges(2,3) * t638 + t657 * t506 - t655 * t507 + (-t511 * t649 - t512 * t652) * pkin(7);
t503 = -mrSges(2,1) * t646 + mrSges(2,3) * t639 - pkin(1) * t511 - t649 * t505 + t663 * t652;
t1 = [-m(1) * g(1) + t669; -m(1) * g(2) + t679; -m(1) * g(3) + m(2) * t646 + t511; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t679 - t648 * t503 + t651 * t504; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t669 + t651 * t503 + t648 * t504; -mrSges(1,1) * g(2) + mrSges(2,1) * t638 + mrSges(1,2) * g(1) - mrSges(2,2) * t639 + pkin(1) * t512 + t652 * t505 + t663 * t649;];
tauB  = t1;
