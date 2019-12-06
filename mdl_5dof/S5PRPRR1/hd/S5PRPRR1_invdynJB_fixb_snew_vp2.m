% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRPRR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:42:58
% EndTime: 2019-12-05 15:43:03
% DurationCPUTime: 4.63s
% Computational Cost: add. (50159->238), mult. (114294->302), div. (0->0), fcn. (81450->10), ass. (0->111)
t651 = qJD(2) ^ 2;
t643 = cos(pkin(9));
t680 = pkin(3) * t643;
t641 = sin(pkin(9));
t679 = mrSges(4,2) * t641;
t637 = t643 ^ 2;
t678 = t637 * t651;
t642 = sin(pkin(8));
t644 = cos(pkin(8));
t623 = t642 * g(1) - t644 * g(2);
t624 = -t644 * g(1) - t642 * g(2);
t647 = sin(qJ(2));
t650 = cos(qJ(2));
t611 = t647 * t623 + t650 * t624;
t608 = -t651 * pkin(2) + qJDD(2) * qJ(3) + t611;
t640 = -g(3) + qJDD(1);
t673 = qJD(2) * qJD(3);
t676 = t643 * t640 - 0.2e1 * t641 * t673;
t589 = (-pkin(6) * qJDD(2) + t651 * t680 - t608) * t641 + t676;
t593 = t641 * t640 + (t608 + 0.2e1 * t673) * t643;
t672 = qJDD(2) * t643;
t590 = -pkin(3) * t678 + pkin(6) * t672 + t593;
t646 = sin(qJ(4));
t649 = cos(qJ(4));
t571 = t649 * t589 - t646 * t590;
t659 = t641 * t649 + t643 * t646;
t658 = -t641 * t646 + t643 * t649;
t616 = t658 * qJD(2);
t674 = t616 * qJD(4);
t607 = t659 * qJDD(2) + t674;
t617 = t659 * qJD(2);
t567 = (-t607 + t674) * pkin(7) + (t616 * t617 + qJDD(4)) * pkin(4) + t571;
t572 = t646 * t589 + t649 * t590;
t606 = -t617 * qJD(4) + t658 * qJDD(2);
t614 = qJD(4) * pkin(4) - t617 * pkin(7);
t615 = t616 ^ 2;
t568 = -t615 * pkin(4) + t606 * pkin(7) - qJD(4) * t614 + t572;
t645 = sin(qJ(5));
t648 = cos(qJ(5));
t565 = t648 * t567 - t645 * t568;
t599 = t648 * t616 - t645 * t617;
t579 = t599 * qJD(5) + t645 * t606 + t648 * t607;
t600 = t645 * t616 + t648 * t617;
t585 = -t599 * mrSges(6,1) + t600 * mrSges(6,2);
t638 = qJD(4) + qJD(5);
t594 = -t638 * mrSges(6,2) + t599 * mrSges(6,3);
t635 = qJDD(4) + qJDD(5);
t562 = m(6) * t565 + t635 * mrSges(6,1) - t579 * mrSges(6,3) - t600 * t585 + t638 * t594;
t566 = t645 * t567 + t648 * t568;
t578 = -t600 * qJD(5) + t648 * t606 - t645 * t607;
t595 = t638 * mrSges(6,1) - t600 * mrSges(6,3);
t563 = m(6) * t566 - t635 * mrSges(6,2) + t578 * mrSges(6,3) + t599 * t585 - t638 * t595;
t552 = t648 * t562 + t645 * t563;
t603 = -t616 * mrSges(5,1) + t617 * mrSges(5,2);
t612 = -qJD(4) * mrSges(5,2) + t616 * mrSges(5,3);
t550 = m(5) * t571 + qJDD(4) * mrSges(5,1) - t607 * mrSges(5,3) + qJD(4) * t612 - t617 * t603 + t552;
t613 = qJD(4) * mrSges(5,1) - t617 * mrSges(5,3);
t666 = -t645 * t562 + t648 * t563;
t551 = m(5) * t572 - qJDD(4) * mrSges(5,2) + t606 * mrSges(5,3) - qJD(4) * t613 + t616 * t603 + t666;
t546 = t649 * t550 + t646 * t551;
t592 = -t641 * t608 + t676;
t657 = mrSges(4,3) * qJDD(2) + t651 * (-mrSges(4,1) * t643 + t679);
t544 = m(4) * t592 - t657 * t641 + t546;
t667 = -t646 * t550 + t649 * t551;
t545 = m(4) * t593 + t657 * t643 + t667;
t668 = -t641 * t544 + t643 * t545;
t536 = m(3) * t611 - t651 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t668;
t610 = t650 * t623 - t647 * t624;
t661 = qJDD(3) - t610;
t605 = -qJDD(2) * pkin(2) - t651 * qJ(3) + t661;
t636 = t641 ^ 2;
t591 = (-pkin(2) - t680) * qJDD(2) + (-qJ(3) + (-t636 - t637) * pkin(6)) * t651 + t661;
t570 = -t606 * pkin(4) - t615 * pkin(7) + t617 * t614 + t591;
t660 = m(6) * t570 - t578 * mrSges(6,1) + t579 * mrSges(6,2) - t599 * t594 + t600 * t595;
t654 = m(5) * t591 - t606 * mrSges(5,1) + t607 * mrSges(5,2) - t616 * t612 + t617 * t613 + t660;
t653 = -m(4) * t605 + mrSges(4,1) * t672 - t654 + (t636 * t651 + t678) * mrSges(4,3);
t556 = t653 + (mrSges(3,1) - t679) * qJDD(2) - t651 * mrSges(3,2) + m(3) * t610;
t533 = t647 * t536 + t650 * t556;
t531 = m(2) * t623 + t533;
t669 = t650 * t536 - t647 * t556;
t532 = m(2) * t624 + t669;
t677 = t644 * t531 + t642 * t532;
t538 = t643 * t544 + t641 * t545;
t662 = Ifges(4,5) * t641 + Ifges(4,6) * t643;
t675 = t651 * t662;
t671 = m(3) * t640 + t538;
t670 = -t642 * t531 + t644 * t532;
t665 = m(2) * t640 + t671;
t664 = Ifges(4,1) * t641 + Ifges(4,4) * t643;
t663 = Ifges(4,4) * t641 + Ifges(4,2) * t643;
t580 = Ifges(6,5) * t600 + Ifges(6,6) * t599 + Ifges(6,3) * t638;
t582 = Ifges(6,1) * t600 + Ifges(6,4) * t599 + Ifges(6,5) * t638;
t553 = -mrSges(6,1) * t570 + mrSges(6,3) * t566 + Ifges(6,4) * t579 + Ifges(6,2) * t578 + Ifges(6,6) * t635 - t600 * t580 + t638 * t582;
t581 = Ifges(6,4) * t600 + Ifges(6,2) * t599 + Ifges(6,6) * t638;
t554 = mrSges(6,2) * t570 - mrSges(6,3) * t565 + Ifges(6,1) * t579 + Ifges(6,4) * t578 + Ifges(6,5) * t635 + t599 * t580 - t638 * t581;
t596 = Ifges(5,5) * t617 + Ifges(5,6) * t616 + Ifges(5,3) * qJD(4);
t598 = Ifges(5,1) * t617 + Ifges(5,4) * t616 + Ifges(5,5) * qJD(4);
t539 = -mrSges(5,1) * t591 + mrSges(5,3) * t572 + Ifges(5,4) * t607 + Ifges(5,2) * t606 + Ifges(5,6) * qJDD(4) - pkin(4) * t660 + pkin(7) * t666 + qJD(4) * t598 + t648 * t553 + t645 * t554 - t617 * t596;
t597 = Ifges(5,4) * t617 + Ifges(5,2) * t616 + Ifges(5,6) * qJD(4);
t540 = mrSges(5,2) * t591 - mrSges(5,3) * t571 + Ifges(5,1) * t607 + Ifges(5,4) * t606 + Ifges(5,5) * qJDD(4) - pkin(7) * t552 - qJD(4) * t597 - t645 * t553 + t648 * t554 + t616 * t596;
t525 = -mrSges(4,1) * t605 + mrSges(4,3) * t593 - pkin(3) * t654 + pkin(6) * t667 + t663 * qJDD(2) + t649 * t539 + t646 * t540 - t641 * t675;
t527 = mrSges(4,2) * t605 - mrSges(4,3) * t592 - pkin(6) * t546 + t664 * qJDD(2) - t646 * t539 + t649 * t540 + t643 * t675;
t558 = qJDD(2) * t679 - t653;
t656 = mrSges(3,1) * t610 - mrSges(3,2) * t611 + Ifges(3,3) * qJDD(2) - pkin(2) * t558 + qJ(3) * t668 + t643 * t525 + t641 * t527;
t655 = -mrSges(6,1) * t565 + mrSges(6,2) * t566 - Ifges(6,5) * t579 - Ifges(6,6) * t578 - Ifges(6,3) * t635 - t600 * t581 + t599 * t582;
t652 = mrSges(5,1) * t571 - mrSges(5,2) * t572 + Ifges(5,5) * t607 + Ifges(5,6) * t606 + Ifges(5,3) * qJDD(4) + pkin(4) * t552 + t617 * t597 - t616 * t598 - t655;
t523 = -t652 + (Ifges(3,6) - t662) * qJDD(2) - mrSges(3,1) * t640 + mrSges(3,3) * t611 + mrSges(4,2) * t593 - mrSges(4,1) * t592 - pkin(3) * t546 - pkin(2) * t538 + (-t641 * t663 + t643 * t664 + Ifges(3,5)) * t651;
t522 = mrSges(3,2) * t640 - mrSges(3,3) * t610 + Ifges(3,5) * qJDD(2) - t651 * Ifges(3,6) - qJ(3) * t538 - t641 * t525 + t643 * t527;
t521 = mrSges(2,2) * t640 - mrSges(2,3) * t623 - pkin(5) * t533 + t650 * t522 - t647 * t523;
t520 = -mrSges(2,1) * t640 + mrSges(2,3) * t624 - pkin(1) * t671 + pkin(5) * t669 + t647 * t522 + t650 * t523;
t1 = [-m(1) * g(1) + t670; -m(1) * g(2) + t677; -m(1) * g(3) + t665; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t677 - t642 * t520 + t644 * t521; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t670 + t644 * t520 + t642 * t521; -mrSges(1,1) * g(2) + mrSges(2,1) * t623 + mrSges(1,2) * g(1) - mrSges(2,2) * t624 + pkin(1) * t533 + t656; t665; t656; t558; t652; -t655;];
tauJB = t1;
