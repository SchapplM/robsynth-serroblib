% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRPR2
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
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:04:52
% EndTime: 2019-05-05 22:05:08
% DurationCPUTime: 15.10s
% Computational Cost: add. (257616->342), mult. (518873->429), div. (0->0), fcn. (352653->12), ass. (0->134)
t641 = sin(qJ(1));
t645 = cos(qJ(1));
t625 = t641 * g(1) - g(2) * t645;
t617 = qJDD(1) * pkin(1) + t625;
t626 = -g(1) * t645 - g(2) * t641;
t647 = qJD(1) ^ 2;
t619 = -pkin(1) * t647 + t626;
t635 = sin(pkin(10));
t637 = cos(pkin(10));
t598 = t635 * t617 + t637 * t619;
t583 = -pkin(2) * t647 + qJDD(1) * pkin(7) + t598;
t633 = -g(3) + qJDD(2);
t640 = sin(qJ(3));
t644 = cos(qJ(3));
t576 = t644 * t583 + t640 * t633;
t618 = (-mrSges(4,1) * t644 + mrSges(4,2) * t640) * qJD(1);
t660 = qJD(1) * qJD(3);
t630 = t640 * t660;
t622 = qJDD(1) * t644 - t630;
t662 = qJD(1) * t640;
t623 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t662;
t597 = t637 * t617 - t635 * t619;
t582 = -qJDD(1) * pkin(2) - t647 * pkin(7) - t597;
t658 = t644 * t660;
t621 = qJDD(1) * t640 + t658;
t565 = (-t621 - t658) * pkin(8) + (-t622 + t630) * pkin(3) + t582;
t620 = (-pkin(3) * t644 - pkin(8) * t640) * qJD(1);
t646 = qJD(3) ^ 2;
t661 = qJD(1) * t644;
t573 = -pkin(3) * t646 + qJDD(3) * pkin(8) + t620 * t661 + t576;
t639 = sin(qJ(4));
t643 = cos(qJ(4));
t554 = t643 * t565 - t639 * t573;
t615 = qJD(3) * t643 - t639 * t662;
t593 = qJD(4) * t615 + qJDD(3) * t639 + t621 * t643;
t614 = qJDD(4) - t622;
t616 = qJD(3) * t639 + t643 * t662;
t628 = qJD(4) - t661;
t541 = (t615 * t628 - t593) * qJ(5) + (t615 * t616 + t614) * pkin(4) + t554;
t555 = t639 * t565 + t643 * t573;
t592 = -qJD(4) * t616 + qJDD(3) * t643 - t621 * t639;
t601 = pkin(4) * t628 - qJ(5) * t616;
t613 = t615 ^ 2;
t547 = -pkin(4) * t613 + qJ(5) * t592 - t601 * t628 + t555;
t634 = sin(pkin(11));
t636 = cos(pkin(11));
t596 = t615 * t634 + t616 * t636;
t535 = -0.2e1 * qJD(5) * t596 + t636 * t541 - t634 * t547;
t567 = t592 * t634 + t593 * t636;
t595 = t615 * t636 - t616 * t634;
t533 = (t595 * t628 - t567) * pkin(9) + (t595 * t596 + t614) * pkin(5) + t535;
t536 = 0.2e1 * qJD(5) * t595 + t634 * t541 + t636 * t547;
t566 = t592 * t636 - t593 * t634;
t579 = pkin(5) * t628 - pkin(9) * t596;
t594 = t595 ^ 2;
t534 = -pkin(5) * t594 + pkin(9) * t566 - t579 * t628 + t536;
t638 = sin(qJ(6));
t642 = cos(qJ(6));
t531 = t533 * t642 - t534 * t638;
t570 = t595 * t642 - t596 * t638;
t546 = qJD(6) * t570 + t566 * t638 + t567 * t642;
t571 = t595 * t638 + t596 * t642;
t556 = -mrSges(7,1) * t570 + mrSges(7,2) * t571;
t627 = qJD(6) + t628;
t557 = -mrSges(7,2) * t627 + mrSges(7,3) * t570;
t610 = qJDD(6) + t614;
t527 = m(7) * t531 + mrSges(7,1) * t610 - mrSges(7,3) * t546 - t556 * t571 + t557 * t627;
t532 = t533 * t638 + t534 * t642;
t545 = -qJD(6) * t571 + t566 * t642 - t567 * t638;
t558 = mrSges(7,1) * t627 - mrSges(7,3) * t571;
t528 = m(7) * t532 - mrSges(7,2) * t610 + mrSges(7,3) * t545 + t556 * t570 - t558 * t627;
t521 = t642 * t527 + t638 * t528;
t574 = -mrSges(6,1) * t595 + mrSges(6,2) * t596;
t577 = -mrSges(6,2) * t628 + mrSges(6,3) * t595;
t519 = m(6) * t535 + mrSges(6,1) * t614 - mrSges(6,3) * t567 - t574 * t596 + t577 * t628 + t521;
t578 = mrSges(6,1) * t628 - mrSges(6,3) * t596;
t652 = -t527 * t638 + t642 * t528;
t520 = m(6) * t536 - mrSges(6,2) * t614 + mrSges(6,3) * t566 + t574 * t595 - t578 * t628 + t652;
t515 = t636 * t519 + t634 * t520;
t599 = -mrSges(5,1) * t615 + mrSges(5,2) * t616;
t600 = -mrSges(5,2) * t628 + mrSges(5,3) * t615;
t513 = m(5) * t554 + mrSges(5,1) * t614 - mrSges(5,3) * t593 - t599 * t616 + t600 * t628 + t515;
t602 = mrSges(5,1) * t628 - mrSges(5,3) * t616;
t653 = -t519 * t634 + t636 * t520;
t514 = m(5) * t555 - mrSges(5,2) * t614 + mrSges(5,3) * t592 + t599 * t615 - t602 * t628 + t653;
t654 = -t513 * t639 + t643 * t514;
t508 = m(4) * t576 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t622 - qJD(3) * t623 + t618 * t661 + t654;
t575 = -t640 * t583 + t633 * t644;
t624 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t661;
t572 = -qJDD(3) * pkin(3) - pkin(8) * t646 + t620 * t662 - t575;
t551 = -pkin(4) * t592 - qJ(5) * t613 + t616 * t601 + qJDD(5) + t572;
t538 = -pkin(5) * t566 - pkin(9) * t594 + t579 * t596 + t551;
t651 = m(7) * t538 - t545 * mrSges(7,1) + t546 * mrSges(7,2) - t570 * t557 + t571 * t558;
t650 = m(6) * t551 - t566 * mrSges(6,1) + t567 * mrSges(6,2) - t595 * t577 + t596 * t578 + t651;
t648 = -m(5) * t572 + t592 * mrSges(5,1) - t593 * mrSges(5,2) + t615 * t600 - t616 * t602 - t650;
t530 = m(4) * t575 + qJDD(3) * mrSges(4,1) - t621 * mrSges(4,3) + qJD(3) * t624 - t618 * t662 + t648;
t655 = t644 * t508 - t530 * t640;
t502 = m(3) * t598 - mrSges(3,1) * t647 - qJDD(1) * mrSges(3,2) + t655;
t509 = t513 * t643 + t514 * t639;
t649 = -m(4) * t582 + t622 * mrSges(4,1) - mrSges(4,2) * t621 - t623 * t662 + t624 * t661 - t509;
t505 = m(3) * t597 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t647 + t649;
t496 = t635 * t502 + t637 * t505;
t494 = m(2) * t625 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t647 + t496;
t656 = t637 * t502 - t505 * t635;
t495 = m(2) * t626 - mrSges(2,1) * t647 - qJDD(1) * mrSges(2,2) + t656;
t663 = t645 * t494 + t641 * t495;
t503 = t640 * t508 + t644 * t530;
t659 = m(3) * t633 + t503;
t657 = -t494 * t641 + t645 * t495;
t609 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t640 + Ifges(4,4) * t644) * qJD(1);
t608 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t640 + Ifges(4,2) * t644) * qJD(1);
t607 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t640 + Ifges(4,6) * t644) * qJD(1);
t586 = Ifges(5,1) * t616 + Ifges(5,4) * t615 + Ifges(5,5) * t628;
t585 = Ifges(5,4) * t616 + Ifges(5,2) * t615 + Ifges(5,6) * t628;
t584 = Ifges(5,5) * t616 + Ifges(5,6) * t615 + Ifges(5,3) * t628;
t564 = Ifges(6,1) * t596 + Ifges(6,4) * t595 + Ifges(6,5) * t628;
t563 = Ifges(6,4) * t596 + Ifges(6,2) * t595 + Ifges(6,6) * t628;
t562 = Ifges(6,5) * t596 + Ifges(6,6) * t595 + Ifges(6,3) * t628;
t550 = Ifges(7,1) * t571 + Ifges(7,4) * t570 + Ifges(7,5) * t627;
t549 = Ifges(7,4) * t571 + Ifges(7,2) * t570 + Ifges(7,6) * t627;
t548 = Ifges(7,5) * t571 + Ifges(7,6) * t570 + Ifges(7,3) * t627;
t523 = mrSges(7,2) * t538 - mrSges(7,3) * t531 + Ifges(7,1) * t546 + Ifges(7,4) * t545 + Ifges(7,5) * t610 + t548 * t570 - t549 * t627;
t522 = -mrSges(7,1) * t538 + mrSges(7,3) * t532 + Ifges(7,4) * t546 + Ifges(7,2) * t545 + Ifges(7,6) * t610 - t548 * t571 + t550 * t627;
t511 = mrSges(6,2) * t551 - mrSges(6,3) * t535 + Ifges(6,1) * t567 + Ifges(6,4) * t566 + Ifges(6,5) * t614 - pkin(9) * t521 - t522 * t638 + t523 * t642 + t562 * t595 - t563 * t628;
t510 = -mrSges(6,1) * t551 + mrSges(6,3) * t536 + Ifges(6,4) * t567 + Ifges(6,2) * t566 + Ifges(6,6) * t614 - pkin(5) * t651 + pkin(9) * t652 + t642 * t522 + t638 * t523 - t596 * t562 + t628 * t564;
t499 = mrSges(5,2) * t572 - mrSges(5,3) * t554 + Ifges(5,1) * t593 + Ifges(5,4) * t592 + Ifges(5,5) * t614 - qJ(5) * t515 - t510 * t634 + t511 * t636 + t584 * t615 - t585 * t628;
t498 = -mrSges(5,1) * t572 + mrSges(5,3) * t555 + Ifges(5,4) * t593 + Ifges(5,2) * t592 + Ifges(5,6) * t614 - pkin(4) * t650 + qJ(5) * t653 + t636 * t510 + t634 * t511 - t616 * t584 + t628 * t586;
t497 = Ifges(4,6) * qJDD(3) - t607 * t662 + (-Ifges(5,3) - Ifges(6,3)) * t614 - t616 * t585 + Ifges(4,4) * t621 + Ifges(4,2) * t622 + qJD(3) * t609 - Ifges(7,3) * t610 + t615 * t586 + t595 * t564 - t596 * t563 - mrSges(4,1) * t582 - Ifges(5,6) * t592 - Ifges(5,5) * t593 - t571 * t549 + mrSges(4,3) * t576 - Ifges(6,6) * t566 - Ifges(6,5) * t567 + t570 * t550 - mrSges(5,1) * t554 + mrSges(5,2) * t555 - Ifges(7,6) * t545 - Ifges(7,5) * t546 + mrSges(6,2) * t536 + mrSges(7,2) * t532 - mrSges(6,1) * t535 - mrSges(7,1) * t531 - pkin(5) * t521 - pkin(4) * t515 - pkin(3) * t509;
t490 = mrSges(4,2) * t582 - mrSges(4,3) * t575 + Ifges(4,1) * t621 + Ifges(4,4) * t622 + Ifges(4,5) * qJDD(3) - pkin(8) * t509 - qJD(3) * t608 - t498 * t639 + t499 * t643 + t607 * t661;
t489 = Ifges(3,6) * qJDD(1) + t647 * Ifges(3,5) - mrSges(3,1) * t633 + mrSges(3,3) * t598 - Ifges(4,5) * t621 - Ifges(4,6) * t622 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t575 + mrSges(4,2) * t576 - t639 * t499 - t643 * t498 - pkin(3) * t648 - pkin(8) * t654 - pkin(2) * t503 + (-t608 * t640 + t609 * t644) * qJD(1);
t488 = mrSges(3,2) * t633 - mrSges(3,3) * t597 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t647 - pkin(7) * t503 + t490 * t644 - t497 * t640;
t487 = -mrSges(2,2) * g(3) - mrSges(2,3) * t625 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t647 - qJ(2) * t496 + t488 * t637 - t489 * t635;
t486 = mrSges(2,1) * g(3) + mrSges(2,3) * t626 + t647 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t659 + qJ(2) * t656 + t635 * t488 + t637 * t489;
t1 = [-m(1) * g(1) + t657; -m(1) * g(2) + t663; (-m(1) - m(2)) * g(3) + t659; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t663 - t641 * t486 + t645 * t487; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t657 + t645 * t486 + t641 * t487; pkin(1) * t496 + mrSges(2,1) * t625 - mrSges(2,2) * t626 + pkin(7) * t655 + t640 * t490 + t644 * t497 + pkin(2) * t649 + mrSges(3,1) * t597 - mrSges(3,2) * t598 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
