% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 21:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:55:59
% EndTime: 2019-05-05 21:56:13
% DurationCPUTime: 13.19s
% Computational Cost: add. (215986->344), mult. (462462->434), div. (0->0), fcn. (322514->12), ass. (0->136)
t685 = 2 * qJD(5);
t661 = sin(qJ(1));
t665 = cos(qJ(1));
t643 = g(1) * t661 - g(2) * t665;
t635 = qJDD(1) * pkin(1) + t643;
t644 = -g(1) * t665 - g(2) * t661;
t666 = qJD(1) ^ 2;
t637 = -pkin(1) * t666 + t644;
t655 = sin(pkin(10));
t657 = cos(pkin(10));
t619 = t635 * t655 + t637 * t657;
t616 = -pkin(2) * t666 + qJDD(1) * pkin(7) + t619;
t653 = -g(3) + qJDD(2);
t660 = sin(qJ(3));
t664 = cos(qJ(3));
t599 = -t660 * t616 + t653 * t664;
t681 = qJD(1) * qJD(3);
t679 = t664 * t681;
t638 = qJDD(1) * t660 + t679;
t593 = (-t638 + t679) * pkin(8) + (t660 * t664 * t666 + qJDD(3)) * pkin(3) + t599;
t600 = t616 * t664 + t653 * t660;
t639 = qJDD(1) * t664 - t660 * t681;
t683 = qJD(1) * t660;
t642 = qJD(3) * pkin(3) - pkin(8) * t683;
t652 = t664 ^ 2;
t594 = -pkin(3) * t652 * t666 + pkin(8) * t639 - qJD(3) * t642 + t600;
t659 = sin(qJ(4));
t663 = cos(qJ(4));
t569 = t593 * t663 - t659 * t594;
t630 = (-t659 * t660 + t663 * t664) * qJD(1);
t602 = qJD(4) * t630 + t638 * t663 + t639 * t659;
t631 = (t659 * t664 + t660 * t663) * qJD(1);
t650 = qJDD(3) + qJDD(4);
t651 = qJD(3) + qJD(4);
t561 = (t630 * t651 - t602) * qJ(5) + (t630 * t631 + t650) * pkin(4) + t569;
t570 = t593 * t659 + t594 * t663;
t601 = -qJD(4) * t631 - t638 * t659 + t639 * t663;
t621 = pkin(4) * t651 - qJ(5) * t631;
t623 = t630 ^ 2;
t563 = -pkin(4) * t623 + qJ(5) * t601 - t621 * t651 + t570;
t654 = sin(pkin(11));
t656 = cos(pkin(11));
t613 = t630 * t656 - t631 * t654;
t558 = t561 * t654 + t563 * t656 + t613 * t685;
t578 = t601 * t656 - t602 * t654;
t614 = t630 * t654 + t631 * t656;
t591 = -mrSges(6,1) * t613 + mrSges(6,2) * t614;
t604 = mrSges(6,1) * t651 - mrSges(6,3) * t614;
t592 = -pkin(5) * t613 - pkin(9) * t614;
t649 = t651 ^ 2;
t556 = -pkin(5) * t649 + pkin(9) * t650 + t592 * t613 + t558;
t618 = t657 * t635 - t637 * t655;
t671 = -qJDD(1) * pkin(2) - t618;
t595 = -t639 * pkin(3) + t642 * t683 + (-pkin(8) * t652 - pkin(7)) * t666 + t671;
t565 = -t601 * pkin(4) - t623 * qJ(5) + t621 * t631 + qJDD(5) + t595;
t579 = t601 * t654 + t602 * t656;
t559 = (-t613 * t651 - t579) * pkin(9) + (t614 * t651 - t578) * pkin(5) + t565;
t658 = sin(qJ(6));
t662 = cos(qJ(6));
t553 = -t556 * t658 + t559 * t662;
t597 = -t614 * t658 + t651 * t662;
t568 = qJD(6) * t597 + t579 * t662 + t650 * t658;
t577 = qJDD(6) - t578;
t598 = t614 * t662 + t651 * t658;
t580 = -mrSges(7,1) * t597 + mrSges(7,2) * t598;
t607 = qJD(6) - t613;
t581 = -mrSges(7,2) * t607 + mrSges(7,3) * t597;
t551 = m(7) * t553 + mrSges(7,1) * t577 - mrSges(7,3) * t568 - t580 * t598 + t581 * t607;
t554 = t556 * t662 + t559 * t658;
t567 = -qJD(6) * t598 - t579 * t658 + t650 * t662;
t582 = mrSges(7,1) * t607 - mrSges(7,3) * t598;
t552 = m(7) * t554 - mrSges(7,2) * t577 + mrSges(7,3) * t567 + t580 * t597 - t582 * t607;
t673 = -t551 * t658 + t552 * t662;
t542 = m(6) * t558 - mrSges(6,2) * t650 + mrSges(6,3) * t578 + t591 * t613 - t604 * t651 + t673;
t672 = -t656 * t561 + t654 * t563;
t557 = -0.2e1 * qJD(5) * t614 - t672;
t603 = -mrSges(6,2) * t651 + mrSges(6,3) * t613;
t555 = -t650 * pkin(5) - t649 * pkin(9) + (t685 + t592) * t614 + t672;
t669 = -m(7) * t555 + mrSges(7,1) * t567 - mrSges(7,2) * t568 + t581 * t597 - t582 * t598;
t547 = m(6) * t557 + mrSges(6,1) * t650 - mrSges(6,3) * t579 - t591 * t614 + t603 * t651 + t669;
t537 = t542 * t654 + t547 * t656;
t617 = -mrSges(5,1) * t630 + mrSges(5,2) * t631;
t620 = -mrSges(5,2) * t651 + mrSges(5,3) * t630;
t535 = m(5) * t569 + mrSges(5,1) * t650 - mrSges(5,3) * t602 - t617 * t631 + t620 * t651 + t537;
t622 = mrSges(5,1) * t651 - mrSges(5,3) * t631;
t674 = t542 * t656 - t547 * t654;
t536 = m(5) * t570 - mrSges(5,2) * t650 + mrSges(5,3) * t601 + t617 * t630 - t622 * t651 + t674;
t529 = t535 * t663 + t536 * t659;
t636 = (-mrSges(4,1) * t664 + mrSges(4,2) * t660) * qJD(1);
t682 = qJD(1) * t664;
t641 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t682;
t527 = m(4) * t599 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t638 + qJD(3) * t641 - t636 * t683 + t529;
t640 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t683;
t675 = -t535 * t659 + t536 * t663;
t528 = m(4) * t600 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t639 - qJD(3) * t640 + t636 * t682 + t675;
t676 = -t527 * t660 + t528 * t664;
t521 = m(3) * t619 - mrSges(3,1) * t666 - qJDD(1) * mrSges(3,2) + t676;
t615 = -t666 * pkin(7) + t671;
t543 = t551 * t662 + t552 * t658;
t670 = m(6) * t565 - mrSges(6,1) * t578 + mrSges(6,2) * t579 - t603 * t613 + t604 * t614 + t543;
t668 = m(5) * t595 - mrSges(5,1) * t601 + mrSges(5,2) * t602 - t620 * t630 + t622 * t631 + t670;
t667 = -m(4) * t615 + mrSges(4,1) * t639 - mrSges(4,2) * t638 - t640 * t683 + t641 * t682 - t668;
t539 = m(3) * t618 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t666 + t667;
t517 = t521 * t655 + t539 * t657;
t515 = m(2) * t643 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t666 + t517;
t677 = t521 * t657 - t539 * t655;
t516 = m(2) * t644 - mrSges(2,1) * t666 - qJDD(1) * mrSges(2,2) + t677;
t684 = t515 * t665 + t516 * t661;
t522 = t527 * t664 + t528 * t660;
t680 = m(3) * t653 + t522;
t678 = -t515 * t661 + t516 * t665;
t629 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t660 + Ifges(4,4) * t664) * qJD(1);
t628 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t660 + Ifges(4,2) * t664) * qJD(1);
t627 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t660 + Ifges(4,6) * t664) * qJD(1);
t610 = Ifges(5,1) * t631 + Ifges(5,4) * t630 + Ifges(5,5) * t651;
t609 = Ifges(5,4) * t631 + Ifges(5,2) * t630 + Ifges(5,6) * t651;
t608 = Ifges(5,5) * t631 + Ifges(5,6) * t630 + Ifges(5,3) * t651;
t585 = Ifges(6,1) * t614 + Ifges(6,4) * t613 + Ifges(6,5) * t651;
t584 = Ifges(6,4) * t614 + Ifges(6,2) * t613 + Ifges(6,6) * t651;
t583 = Ifges(6,5) * t614 + Ifges(6,6) * t613 + Ifges(6,3) * t651;
t573 = Ifges(7,1) * t598 + Ifges(7,4) * t597 + Ifges(7,5) * t607;
t572 = Ifges(7,4) * t598 + Ifges(7,2) * t597 + Ifges(7,6) * t607;
t571 = Ifges(7,5) * t598 + Ifges(7,6) * t597 + Ifges(7,3) * t607;
t545 = mrSges(7,2) * t555 - mrSges(7,3) * t553 + Ifges(7,1) * t568 + Ifges(7,4) * t567 + Ifges(7,5) * t577 + t571 * t597 - t572 * t607;
t544 = -mrSges(7,1) * t555 + mrSges(7,3) * t554 + Ifges(7,4) * t568 + Ifges(7,2) * t567 + Ifges(7,6) * t577 - t571 * t598 + t573 * t607;
t531 = -mrSges(6,1) * t565 - mrSges(7,1) * t553 + mrSges(7,2) * t554 + mrSges(6,3) * t558 + Ifges(6,4) * t579 - Ifges(7,5) * t568 + Ifges(6,2) * t578 + Ifges(6,6) * t650 - Ifges(7,6) * t567 - Ifges(7,3) * t577 - pkin(5) * t543 - t572 * t598 + t573 * t597 - t583 * t614 + t585 * t651;
t530 = mrSges(6,2) * t565 - mrSges(6,3) * t557 + Ifges(6,1) * t579 + Ifges(6,4) * t578 + Ifges(6,5) * t650 - pkin(9) * t543 - t544 * t658 + t545 * t662 + t583 * t613 - t584 * t651;
t523 = mrSges(5,2) * t595 - mrSges(5,3) * t569 + Ifges(5,1) * t602 + Ifges(5,4) * t601 + Ifges(5,5) * t650 - qJ(5) * t537 + t530 * t656 - t531 * t654 + t608 * t630 - t609 * t651;
t518 = -mrSges(5,1) * t595 + mrSges(5,3) * t570 + Ifges(5,4) * t602 + Ifges(5,2) * t601 + Ifges(5,6) * t650 - pkin(4) * t670 + qJ(5) * t674 + t654 * t530 + t656 * t531 - t631 * t608 + t651 * t610;
t511 = Ifges(3,6) * qJDD(1) - pkin(2) * t522 - Ifges(4,3) * qJDD(3) - pkin(3) * t529 + (-Ifges(5,3) - Ifges(6,3)) * t650 + (-t628 * t660 + t629 * t664) * qJD(1) + t666 * Ifges(3,5) - t658 * t545 - t662 * t544 - mrSges(3,1) * t653 - Ifges(4,5) * t638 - Ifges(4,6) * t639 + t630 * t610 - t631 * t609 + t613 * t585 - t614 * t584 + mrSges(3,3) * t619 - Ifges(5,5) * t602 - mrSges(4,1) * t599 + mrSges(4,2) * t600 - Ifges(5,6) * t601 - Ifges(6,6) * t578 - Ifges(6,5) * t579 - mrSges(5,1) * t569 + mrSges(5,2) * t570 + mrSges(6,2) * t558 - mrSges(6,1) * t557 - pkin(4) * t537 - pkin(5) * t669 - pkin(9) * t673;
t510 = mrSges(4,2) * t615 - mrSges(4,3) * t599 + Ifges(4,1) * t638 + Ifges(4,4) * t639 + Ifges(4,5) * qJDD(3) - pkin(8) * t529 - qJD(3) * t628 - t518 * t659 + t523 * t663 + t627 * t682;
t509 = -mrSges(4,1) * t615 + mrSges(4,3) * t600 + Ifges(4,4) * t638 + Ifges(4,2) * t639 + Ifges(4,6) * qJDD(3) - pkin(3) * t668 + pkin(8) * t675 + qJD(3) * t629 + t663 * t518 + t659 * t523 - t627 * t683;
t508 = mrSges(3,2) * t653 - mrSges(3,3) * t618 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t666 - pkin(7) * t522 - t509 * t660 + t510 * t664;
t507 = -mrSges(2,2) * g(3) - mrSges(2,3) * t643 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t666 - qJ(2) * t517 + t508 * t657 - t511 * t655;
t506 = mrSges(2,1) * g(3) + mrSges(2,3) * t644 + t666 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t680 + qJ(2) * t677 + t655 * t508 + t657 * t511;
t1 = [-m(1) * g(1) + t678; -m(1) * g(2) + t684; (-m(1) - m(2)) * g(3) + t680; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t684 - t506 * t661 + t507 * t665; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t678 + t665 * t506 + t661 * t507; pkin(1) * t517 + mrSges(2,1) * t643 - mrSges(2,2) * t644 + t660 * t510 + t664 * t509 + pkin(2) * t667 + pkin(7) * t676 + mrSges(3,1) * t618 - mrSges(3,2) * t619 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
