% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-05-05 06:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRRRPP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:34:04
% EndTime: 2019-05-05 06:34:16
% DurationCPUTime: 10.24s
% Computational Cost: add. (165679->317), mult. (327931->396), div. (0->0), fcn. (230310->12), ass. (0->133)
t693 = Ifges(6,1) + Ifges(7,1);
t688 = Ifges(6,4) - Ifges(7,5);
t687 = Ifges(6,5) + Ifges(7,4);
t692 = Ifges(6,2) + Ifges(7,3);
t691 = -Ifges(7,2) - Ifges(6,3);
t686 = Ifges(6,6) - Ifges(7,6);
t648 = sin(pkin(10));
t650 = cos(pkin(10));
t639 = g(1) * t648 - g(2) * t650;
t640 = -g(1) * t650 - g(2) * t648;
t646 = -g(3) + qJDD(1);
t649 = sin(pkin(6));
t651 = cos(pkin(6));
t654 = sin(qJ(2));
t657 = cos(qJ(2));
t599 = -t654 * t640 + (t639 * t651 + t646 * t649) * t657;
t690 = -2 * qJD(5);
t689 = -mrSges(6,3) - mrSges(7,2);
t685 = cos(pkin(11));
t659 = qJD(2) ^ 2;
t682 = t651 * t654;
t683 = t649 * t654;
t600 = t639 * t682 + t657 * t640 + t646 * t683;
t594 = -pkin(2) * t659 + qJDD(2) * pkin(8) + t600;
t618 = -t639 * t649 + t646 * t651;
t653 = sin(qJ(3));
t656 = cos(qJ(3));
t585 = t656 * t594 + t653 * t618;
t636 = (-pkin(3) * t656 - pkin(9) * t653) * qJD(2);
t658 = qJD(3) ^ 2;
t675 = qJD(2) * t656;
t566 = -pkin(3) * t658 + qJDD(3) * pkin(9) + t636 * t675 + t585;
t593 = -qJDD(2) * pkin(2) - t659 * pkin(8) - t599;
t674 = qJD(2) * qJD(3);
t671 = t656 * t674;
t637 = qJDD(2) * t653 + t671;
t672 = t653 * t674;
t638 = qJDD(2) * t656 - t672;
t569 = (-t637 - t671) * pkin(9) + (-t638 + t672) * pkin(3) + t593;
t652 = sin(qJ(4));
t655 = cos(qJ(4));
t561 = -t652 * t566 + t655 * t569;
t676 = qJD(2) * t653;
t633 = qJD(3) * t655 - t652 * t676;
t610 = qJD(4) * t633 + qJDD(3) * t652 + t637 * t655;
t630 = qJDD(4) - t638;
t634 = qJD(3) * t652 + t655 * t676;
t645 = qJD(4) - t675;
t558 = (t633 * t645 - t610) * qJ(5) + (t633 * t634 + t630) * pkin(4) + t561;
t562 = t655 * t566 + t652 * t569;
t609 = -qJD(4) * t634 + qJDD(3) * t655 - t637 * t652;
t616 = pkin(4) * t645 - qJ(5) * t634;
t629 = t633 ^ 2;
t560 = -pkin(4) * t629 + qJ(5) * t609 - t616 * t645 + t562;
t647 = sin(pkin(11));
t611 = -t685 * t633 + t647 * t634;
t554 = t647 * t558 + t685 * t560 + t611 * t690;
t579 = -t685 * t609 + t647 * t610;
t612 = t647 * t633 + t685 * t634;
t596 = mrSges(6,1) * t645 - mrSges(6,3) * t612;
t586 = pkin(5) * t611 - qJ(6) * t612;
t644 = t645 ^ 2;
t551 = -pkin(5) * t644 + qJ(6) * t630 + 0.2e1 * qJD(6) * t645 - t586 * t611 + t554;
t597 = -mrSges(7,1) * t645 + mrSges(7,2) * t612;
t673 = m(7) * t551 + t630 * mrSges(7,3) + t645 * t597;
t587 = mrSges(7,1) * t611 - mrSges(7,3) * t612;
t677 = -mrSges(6,1) * t611 - mrSges(6,2) * t612 - t587;
t544 = m(6) * t554 - t630 * mrSges(6,2) + t689 * t579 - t645 * t596 + t677 * t611 + t673;
t664 = t685 * t558 - t647 * t560;
t553 = t612 * t690 + t664;
t580 = t647 * t609 + t685 * t610;
t595 = -mrSges(6,2) * t645 - mrSges(6,3) * t611;
t552 = -t630 * pkin(5) - t644 * qJ(6) + qJDD(6) + ((2 * qJD(5)) + t586) * t612 - t664;
t598 = -mrSges(7,2) * t611 + mrSges(7,3) * t645;
t666 = -m(7) * t552 + t630 * mrSges(7,1) + t645 * t598;
t546 = m(6) * t553 + t630 * mrSges(6,1) + t689 * t580 + t645 * t595 + t677 * t612 + t666;
t539 = t647 * t544 + t685 * t546;
t613 = -mrSges(5,1) * t633 + mrSges(5,2) * t634;
t615 = -mrSges(5,2) * t645 + mrSges(5,3) * t633;
t537 = m(5) * t561 + mrSges(5,1) * t630 - mrSges(5,3) * t610 - t613 * t634 + t615 * t645 + t539;
t617 = mrSges(5,1) * t645 - mrSges(5,3) * t634;
t667 = t685 * t544 - t546 * t647;
t538 = m(5) * t562 - mrSges(5,2) * t630 + mrSges(5,3) * t609 + t613 * t633 - t617 * t645 + t667;
t535 = t537 * t655 + t538 * t652;
t641 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t676;
t642 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t675;
t661 = -m(4) * t593 + t638 * mrSges(4,1) - mrSges(4,2) * t637 - t641 * t676 + t642 * t675 - t535;
t531 = m(3) * t599 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t659 + t661;
t684 = t531 * t657;
t635 = (-mrSges(4,1) * t656 + mrSges(4,2) * t653) * qJD(2);
t668 = -t537 * t652 + t655 * t538;
t534 = m(4) * t585 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t638 - qJD(3) * t641 + t635 * t675 + t668;
t584 = -t653 * t594 + t656 * t618;
t565 = -qJDD(3) * pkin(3) - t658 * pkin(9) + t636 * t676 - t584;
t563 = -t609 * pkin(4) - t629 * qJ(5) + t634 * t616 + qJDD(5) + t565;
t556 = -0.2e1 * qJD(6) * t612 + (t611 * t645 - t580) * qJ(6) + (t612 * t645 + t579) * pkin(5) + t563;
t549 = m(7) * t556 + t579 * mrSges(7,1) - t580 * mrSges(7,3) - t612 * t597 + t611 * t598;
t662 = m(6) * t563 + t579 * mrSges(6,1) + t580 * mrSges(6,2) + t611 * t595 + t612 * t596 + t549;
t660 = -m(5) * t565 + t609 * mrSges(5,1) - t610 * mrSges(5,2) + t633 * t615 - t634 * t617 - t662;
t548 = m(4) * t584 + qJDD(3) * mrSges(4,1) - t637 * mrSges(4,3) + qJD(3) * t642 - t635 * t676 + t660;
t669 = t656 * t534 - t548 * t653;
t525 = m(3) * t600 - mrSges(3,1) * t659 - qJDD(2) * mrSges(3,2) + t669;
t528 = t653 * t534 + t656 * t548;
t527 = m(3) * t618 + t528;
t514 = t525 * t682 - t527 * t649 + t651 * t684;
t512 = m(2) * t639 + t514;
t518 = t657 * t525 - t531 * t654;
t517 = m(2) * t640 + t518;
t681 = t650 * t512 + t648 * t517;
t680 = t692 * t611 - t688 * t612 - t686 * t645;
t679 = t686 * t611 - t687 * t612 + t691 * t645;
t678 = -t688 * t611 + t693 * t612 + t687 * t645;
t513 = t525 * t683 + t651 * t527 + t649 * t684;
t670 = -t512 * t648 + t650 * t517;
t540 = -mrSges(6,1) * t563 - mrSges(7,1) * t556 + mrSges(7,2) * t551 + mrSges(6,3) * t554 - pkin(5) * t549 - t692 * t579 + t688 * t580 + t679 * t612 + t686 * t630 + t678 * t645;
t541 = mrSges(6,2) * t563 + mrSges(7,2) * t552 - mrSges(6,3) * t553 - mrSges(7,3) * t556 - qJ(6) * t549 - t688 * t579 + t693 * t580 + t679 * t611 + t687 * t630 + t680 * t645;
t602 = Ifges(5,5) * t634 + Ifges(5,6) * t633 + Ifges(5,3) * t645;
t604 = Ifges(5,1) * t634 + Ifges(5,4) * t633 + Ifges(5,5) * t645;
t520 = -mrSges(5,1) * t565 + mrSges(5,3) * t562 + Ifges(5,4) * t610 + Ifges(5,2) * t609 + Ifges(5,6) * t630 - pkin(4) * t662 + qJ(5) * t667 + t685 * t540 + t647 * t541 - t634 * t602 + t645 * t604;
t603 = Ifges(5,4) * t634 + Ifges(5,2) * t633 + Ifges(5,6) * t645;
t521 = mrSges(5,2) * t565 - mrSges(5,3) * t561 + Ifges(5,1) * t610 + Ifges(5,4) * t609 + Ifges(5,5) * t630 - qJ(5) * t539 - t647 * t540 + t685 * t541 + t633 * t602 - t645 * t603;
t622 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t653 + Ifges(4,6) * t656) * qJD(2);
t623 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t653 + Ifges(4,2) * t656) * qJD(2);
t510 = mrSges(4,2) * t593 - mrSges(4,3) * t584 + Ifges(4,1) * t637 + Ifges(4,4) * t638 + Ifges(4,5) * qJDD(3) - pkin(9) * t535 - qJD(3) * t623 - t520 * t652 + t521 * t655 + t622 * t675;
t624 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t653 + Ifges(4,4) * t656) * qJD(2);
t519 = (mrSges(7,2) * qJ(6) + t686) * t579 + (mrSges(7,2) * pkin(5) - t687) * t580 - mrSges(5,1) * t561 + mrSges(5,2) * t562 - qJ(6) * t673 - t622 * t676 + mrSges(7,1) * t552 - mrSges(6,1) * t553 + mrSges(6,2) * t554 - pkin(5) * t666 + Ifges(4,6) * qJDD(3) - pkin(3) * t535 - pkin(4) * t539 + (-Ifges(5,3) + t691) * t630 + (qJ(6) * t587 - t678) * t611 + (pkin(5) * t587 + t680) * t612 - mrSges(7,3) * t551 + mrSges(4,3) * t585 - mrSges(4,1) * t593 - Ifges(5,6) * t609 - Ifges(5,5) * t610 + qJD(3) * t624 + t633 * t604 - t634 * t603 + Ifges(4,4) * t637 + Ifges(4,2) * t638;
t508 = mrSges(3,2) * t618 - mrSges(3,3) * t599 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t659 - pkin(8) * t528 + t510 * t656 - t519 * t653;
t509 = Ifges(3,6) * qJDD(2) + t659 * Ifges(3,5) - mrSges(3,1) * t618 + mrSges(3,3) * t600 - Ifges(4,5) * t637 - Ifges(4,6) * t638 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t584 + mrSges(4,2) * t585 - t652 * t521 - t655 * t520 - pkin(3) * t660 - pkin(9) * t668 - pkin(2) * t528 + (-t623 * t653 + t624 * t656) * qJD(2);
t663 = pkin(7) * t518 + t508 * t654 + t509 * t657;
t507 = mrSges(3,1) * t599 - mrSges(3,2) * t600 + Ifges(3,3) * qJDD(2) + pkin(2) * t661 + pkin(8) * t669 + t653 * t510 + t656 * t519;
t506 = mrSges(2,2) * t646 - mrSges(2,3) * t639 + t657 * t508 - t654 * t509 + (-t513 * t649 - t514 * t651) * pkin(7);
t505 = -mrSges(2,1) * t646 + mrSges(2,3) * t640 - pkin(1) * t513 - t649 * t507 + t663 * t651;
t1 = [-m(1) * g(1) + t670; -m(1) * g(2) + t681; -m(1) * g(3) + m(2) * t646 + t513; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t681 - t648 * t505 + t650 * t506; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t670 + t650 * t505 + t648 * t506; -mrSges(1,1) * g(2) + mrSges(2,1) * t639 + mrSges(1,2) * g(1) - mrSges(2,2) * t640 + pkin(1) * t514 + t651 * t507 + t663 * t649;];
tauB  = t1;
