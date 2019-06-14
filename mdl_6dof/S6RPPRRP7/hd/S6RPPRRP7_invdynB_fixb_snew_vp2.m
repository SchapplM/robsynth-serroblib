% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-05-05 15:05
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:03:34
% EndTime: 2019-05-05 15:03:40
% DurationCPUTime: 4.22s
% Computational Cost: add. (41454->292), mult. (93001->340), div. (0->0), fcn. (63150->8), ass. (0->123)
t695 = Ifges(6,1) + Ifges(7,1);
t685 = Ifges(6,4) + Ifges(7,4);
t684 = Ifges(6,5) + Ifges(7,5);
t694 = Ifges(6,2) + Ifges(7,2);
t693 = Ifges(6,6) + Ifges(7,6);
t692 = Ifges(6,3) + Ifges(7,3);
t640 = sin(qJ(1));
t643 = cos(qJ(1));
t620 = t640 * g(1) - t643 * g(2);
t645 = qJD(1) ^ 2;
t652 = -t645 * qJ(2) + qJDD(2) - t620;
t681 = -pkin(1) - qJ(3);
t691 = -(2 * qJD(1) * qJD(3)) + t681 * qJDD(1) + t652;
t636 = sin(pkin(9));
t630 = t636 ^ 2;
t637 = cos(pkin(9));
t674 = t637 ^ 2 + t630;
t665 = t674 * mrSges(4,3);
t621 = -t643 * g(1) - t640 * g(2);
t690 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t621;
t639 = sin(qJ(4));
t642 = cos(qJ(4));
t656 = t636 * t642 + t637 * t639;
t616 = t656 * qJD(1);
t655 = -t636 * t639 + t637 * t642;
t617 = t655 * qJD(1);
t671 = t617 * qJD(4);
t601 = -t656 * qJDD(1) - t671;
t689 = pkin(3) * t645;
t688 = mrSges(2,1) - mrSges(3,2);
t687 = -mrSges(6,2) - mrSges(7,2);
t686 = -Ifges(3,4) + Ifges(2,5);
t683 = -Ifges(2,6) + Ifges(3,5);
t680 = mrSges(4,2) * t637;
t598 = t636 * g(3) + t691 * t637;
t581 = (-pkin(7) * qJDD(1) - t636 * t689) * t637 + t598;
t599 = -g(3) * t637 + t691 * t636;
t669 = qJDD(1) * t636;
t582 = -pkin(7) * t669 - t630 * t689 + t599;
t559 = t639 * t581 + t642 * t582;
t595 = mrSges(5,1) * t616 + mrSges(5,2) * t617;
t611 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t617;
t600 = pkin(4) * t616 - pkin(8) * t617;
t644 = qJD(4) ^ 2;
t554 = -pkin(4) * t644 + qJDD(4) * pkin(8) - t600 * t616 + t559;
t651 = qJDD(3) + t690;
t589 = pkin(3) * t669 + (-t674 * pkin(7) + t681) * t645 + t651;
t672 = t616 * qJD(4);
t602 = t655 * qJDD(1) - t672;
t557 = (-t602 + t672) * pkin(8) + (-t601 + t671) * pkin(4) + t589;
t638 = sin(qJ(5));
t641 = cos(qJ(5));
t550 = -t638 * t554 + t641 * t557;
t605 = qJD(4) * t641 - t617 * t638;
t573 = qJD(5) * t605 + qJDD(4) * t638 + t602 * t641;
t606 = qJD(4) * t638 + t617 * t641;
t576 = -mrSges(7,1) * t605 + mrSges(7,2) * t606;
t577 = -mrSges(6,1) * t605 + mrSges(6,2) * t606;
t614 = qJD(5) + t616;
t584 = -mrSges(6,2) * t614 + mrSges(6,3) * t605;
t597 = qJDD(5) - t601;
t546 = -0.2e1 * qJD(6) * t606 + (t605 * t614 - t573) * qJ(6) + (t605 * t606 + t597) * pkin(5) + t550;
t583 = -mrSges(7,2) * t614 + mrSges(7,3) * t605;
t667 = m(7) * t546 + t597 * mrSges(7,1) + t614 * t583;
t538 = m(6) * t550 + t597 * mrSges(6,1) + t614 * t584 + (-t576 - t577) * t606 + (-mrSges(6,3) - mrSges(7,3)) * t573 + t667;
t551 = t641 * t554 + t638 * t557;
t572 = -qJD(5) * t606 + qJDD(4) * t641 - t602 * t638;
t585 = pkin(5) * t614 - qJ(6) * t606;
t604 = t605 ^ 2;
t548 = -pkin(5) * t604 + qJ(6) * t572 + 0.2e1 * qJD(6) * t605 - t585 * t614 + t551;
t666 = m(7) * t548 + t572 * mrSges(7,3) + t605 * t576;
t586 = mrSges(7,1) * t614 - mrSges(7,3) * t606;
t675 = -mrSges(6,1) * t614 + mrSges(6,3) * t606 - t586;
t541 = m(6) * t551 + t572 * mrSges(6,3) + t605 * t577 + t687 * t597 + t675 * t614 + t666;
t661 = -t538 * t638 + t641 * t541;
t534 = m(5) * t559 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t601 - qJD(4) * t611 - t595 * t616 + t661;
t558 = t581 * t642 - t639 * t582;
t610 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t616;
t553 = -qJDD(4) * pkin(4) - pkin(8) * t644 + t617 * t600 - t558;
t549 = -pkin(5) * t572 - qJ(6) * t604 + t585 * t606 + qJDD(6) + t553;
t660 = m(7) * t549 - t572 * mrSges(7,1) - t605 * t583;
t647 = -m(6) * t553 + t572 * mrSges(6,1) + t687 * t573 + t605 * t584 + t675 * t606 - t660;
t543 = m(5) * t558 + qJDD(4) * mrSges(5,1) - t602 * mrSges(5,3) + qJD(4) * t610 - t617 * t595 + t647;
t527 = t639 * t534 + t642 * t543;
t654 = -qJDD(1) * mrSges(4,3) - t645 * (mrSges(4,1) * t636 + t680);
t525 = m(4) * t598 + t654 * t637 + t527;
t662 = t642 * t534 - t639 * t543;
t526 = m(4) * t599 + t654 * t636 + t662;
t522 = t637 * t525 + t636 * t526;
t615 = -qJDD(1) * pkin(1) + t652;
t650 = -m(3) * t615 + t645 * mrSges(3,3) - t522;
t520 = m(2) * t620 - t645 * mrSges(2,2) + t688 * qJDD(1) + t650;
t613 = t645 * pkin(1) - t690;
t609 = t681 * t645 + t651;
t536 = t641 * t538 + t638 * t541;
t649 = m(5) * t589 - t601 * mrSges(5,1) + t602 * mrSges(5,2) + t616 * t610 + t617 * t611 + t536;
t648 = -m(4) * t609 - mrSges(4,1) * t669 - qJDD(1) * t680 - t649;
t646 = -m(3) * t613 + t645 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t648;
t531 = t646 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t665) * t645 + m(2) * t621;
t679 = t643 * t520 + t640 * t531;
t678 = t693 * t605 + t684 * t606 + t692 * t614;
t677 = -t694 * t605 - t685 * t606 - t693 * t614;
t676 = t685 * t605 + t695 * t606 + t684 * t614;
t657 = Ifges(4,5) * t637 - Ifges(4,6) * t636;
t673 = t645 * t657;
t664 = -t520 * t640 + t643 * t531;
t663 = -t636 * t525 + t637 * t526;
t659 = Ifges(4,1) * t637 - Ifges(4,4) * t636;
t658 = Ifges(4,4) * t637 - Ifges(4,2) * t636;
t592 = Ifges(5,1) * t617 - Ifges(5,4) * t616 + Ifges(5,5) * qJD(4);
t591 = Ifges(5,4) * t617 - Ifges(5,2) * t616 + Ifges(5,6) * qJD(4);
t590 = Ifges(5,5) * t617 - Ifges(5,6) * t616 + Ifges(5,3) * qJD(4);
t544 = -t573 * mrSges(7,3) - t606 * t576 + t667;
t535 = mrSges(6,2) * t553 + mrSges(7,2) * t549 - mrSges(6,3) * t550 - mrSges(7,3) * t546 - qJ(6) * t544 + t685 * t572 + t695 * t573 + t684 * t597 + t678 * t605 + t677 * t614;
t528 = -mrSges(6,1) * t553 + mrSges(6,3) * t551 - mrSges(7,1) * t549 + mrSges(7,3) * t548 - pkin(5) * t660 + qJ(6) * t666 + (-qJ(6) * t586 + t676) * t614 + (-pkin(5) * t586 - t678) * t606 + (-mrSges(7,2) * qJ(6) + t693) * t597 + (-mrSges(7,2) * pkin(5) + t685) * t573 + t694 * t572;
t523 = -mrSges(5,1) * t589 - mrSges(6,1) * t550 - mrSges(7,1) * t546 + mrSges(6,2) * t551 + mrSges(7,2) * t548 + mrSges(5,3) * t559 + Ifges(5,4) * t602 + Ifges(5,2) * t601 + Ifges(5,6) * qJDD(4) - pkin(4) * t536 - pkin(5) * t544 + qJD(4) * t592 - t617 * t590 + t677 * t606 + t676 * t605 - t692 * t597 - t684 * t573 - t693 * t572;
t521 = -m(3) * g(3) + t663;
t518 = mrSges(5,2) * t589 - mrSges(5,3) * t558 + Ifges(5,1) * t602 + Ifges(5,4) * t601 + Ifges(5,5) * qJDD(4) - pkin(8) * t536 - qJD(4) * t591 - t528 * t638 + t535 * t641 - t590 * t616;
t517 = mrSges(4,2) * t609 - mrSges(4,3) * t598 - pkin(7) * t527 + t659 * qJDD(1) + t642 * t518 - t639 * t523 - t636 * t673;
t516 = -mrSges(4,1) * t609 + mrSges(4,3) * t599 - pkin(3) * t649 + pkin(7) * t662 + t658 * qJDD(1) + t639 * t518 + t642 * t523 - t637 * t673;
t515 = Ifges(5,3) * qJDD(4) + pkin(8) * t661 + t638 * t535 + t641 * t528 - mrSges(2,3) * t620 + mrSges(3,1) * t615 + t616 * t592 + t617 * t591 + pkin(4) * t647 + Ifges(5,5) * t602 + mrSges(4,1) * t598 - mrSges(4,2) * t599 + Ifges(5,6) * t601 + mrSges(5,1) * t558 - mrSges(5,2) * t559 + pkin(3) * t527 - qJ(2) * t521 + pkin(2) * t522 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t657 + t686) * qJDD(1) + (t636 * t659 + t637 * t658 + t683) * t645;
t514 = mrSges(2,3) * t621 - mrSges(3,1) * t613 - t636 * t517 - t637 * t516 - pkin(2) * t648 - qJ(3) * t663 - pkin(1) * t521 - t683 * qJDD(1) + t688 * g(3) + (-pkin(2) * t665 + t686) * t645;
t1 = [-m(1) * g(1) + t664; -m(1) * g(2) + t679; (-m(1) - m(2) - m(3)) * g(3) + t663; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t679 - t640 * t514 + t643 * t515; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t664 + t643 * t514 + t640 * t515; pkin(1) * t650 + qJ(2) * (-t645 * t665 + t646) + t637 * t517 - t636 * t516 - qJ(3) * t522 + mrSges(2,1) * t620 - mrSges(2,2) * t621 + mrSges(3,2) * t615 - mrSges(3,3) * t613 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
