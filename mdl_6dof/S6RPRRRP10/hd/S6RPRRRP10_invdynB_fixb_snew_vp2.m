% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-05-06 01:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPRRRP10_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP10_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:55:57
% EndTime: 2019-05-06 01:56:07
% DurationCPUTime: 4.85s
% Computational Cost: add. (55703->314), mult. (107563->370), div. (0->0), fcn. (68552->8), ass. (0->123)
t678 = Ifges(6,1) + Ifges(7,1);
t669 = Ifges(6,4) - Ifges(7,5);
t668 = Ifges(7,4) + Ifges(6,5);
t677 = Ifges(6,2) + Ifges(7,3);
t666 = Ifges(6,6) - Ifges(7,6);
t676 = -Ifges(6,3) - Ifges(7,2);
t638 = sin(qJ(1));
t641 = cos(qJ(1));
t622 = -t641 * g(1) - t638 * g(2);
t675 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t622;
t674 = -pkin(1) - pkin(7);
t673 = cos(qJ(5));
t672 = mrSges(2,1) - mrSges(3,2);
t671 = -mrSges(6,3) - mrSges(7,2);
t670 = -Ifges(3,4) + Ifges(2,5);
t667 = (Ifges(3,5) - Ifges(2,6));
t621 = t638 * g(1) - t641 * g(2);
t643 = qJD(1) ^ 2;
t649 = -t643 * qJ(2) + qJDD(2) - t621;
t597 = t674 * qJDD(1) + t649;
t637 = sin(qJ(3));
t640 = cos(qJ(3));
t590 = -g(3) * t640 + t637 * t597;
t615 = (mrSges(4,1) * t637 + mrSges(4,2) * t640) * qJD(1);
t659 = qJD(1) * qJD(3);
t626 = t640 * t659;
t617 = -t637 * qJDD(1) - t626;
t660 = qJD(1) * t640;
t620 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t660;
t628 = t637 * qJD(1);
t596 = t674 * t643 - t675;
t656 = t637 * t659;
t618 = qJDD(1) * t640 - t656;
t567 = (-t618 + t656) * pkin(8) + (-t617 + t626) * pkin(3) + t596;
t616 = (pkin(3) * t637 - pkin(8) * t640) * qJD(1);
t642 = qJD(3) ^ 2;
t572 = -pkin(3) * t642 + qJDD(3) * pkin(8) - t616 * t628 + t590;
t636 = sin(qJ(4));
t639 = cos(qJ(4));
t546 = t639 * t567 - t636 * t572;
t613 = qJD(3) * t639 - t636 * t660;
t585 = qJD(4) * t613 + qJDD(3) * t636 + t618 * t639;
t612 = qJDD(4) - t617;
t614 = qJD(3) * t636 + t639 * t660;
t625 = t628 + qJD(4);
t542 = (t613 * t625 - t585) * pkin(9) + (t613 * t614 + t612) * pkin(4) + t546;
t547 = t636 * t567 + t639 * t572;
t584 = -qJD(4) * t614 + qJDD(3) * t639 - t618 * t636;
t595 = pkin(4) * t625 - pkin(9) * t614;
t611 = t613 ^ 2;
t544 = -pkin(4) * t611 + pkin(9) * t584 - t595 * t625 + t547;
t635 = sin(qJ(5));
t540 = t635 * t542 + t673 * t544;
t587 = t635 * t613 + t673 * t614;
t551 = qJD(5) * t587 - t673 * t584 + t585 * t635;
t624 = qJD(5) + t625;
t575 = mrSges(6,1) * t624 - mrSges(6,3) * t587;
t586 = -t673 * t613 + t614 * t635;
t607 = qJDD(5) + t612;
t562 = pkin(5) * t586 - qJ(6) * t587;
t623 = t624 ^ 2;
t535 = -pkin(5) * t623 + qJ(6) * t607 + 0.2e1 * qJD(6) * t624 - t562 * t586 + t540;
t576 = -mrSges(7,1) * t624 + mrSges(7,2) * t587;
t657 = m(7) * t535 + t607 * mrSges(7,3) + t624 * t576;
t563 = mrSges(7,1) * t586 - mrSges(7,3) * t587;
t661 = -mrSges(6,1) * t586 - mrSges(6,2) * t587 - t563;
t530 = m(6) * t540 - t607 * mrSges(6,2) + t671 * t551 - t624 * t575 + t661 * t586 + t657;
t539 = t673 * t542 - t635 * t544;
t552 = -t586 * qJD(5) + t635 * t584 + t673 * t585;
t574 = -mrSges(6,2) * t624 - mrSges(6,3) * t586;
t536 = -t607 * pkin(5) - t623 * qJ(6) + t587 * t562 + qJDD(6) - t539;
t573 = -mrSges(7,2) * t586 + mrSges(7,3) * t624;
t651 = -m(7) * t536 + t607 * mrSges(7,1) + t624 * t573;
t532 = m(6) * t539 + t607 * mrSges(6,1) + t671 * t552 + t624 * t574 + t661 * t587 + t651;
t526 = t635 * t530 + t673 * t532;
t588 = -mrSges(5,1) * t613 + mrSges(5,2) * t614;
t591 = -mrSges(5,2) * t625 + mrSges(5,3) * t613;
t522 = m(5) * t546 + mrSges(5,1) * t612 - mrSges(5,3) * t585 - t588 * t614 + t591 * t625 + t526;
t592 = mrSges(5,1) * t625 - mrSges(5,3) * t614;
t652 = t673 * t530 - t532 * t635;
t523 = m(5) * t547 - mrSges(5,2) * t612 + mrSges(5,3) * t584 + t588 * t613 - t592 * t625 + t652;
t653 = -t522 * t636 + t639 * t523;
t517 = m(4) * t590 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t617 - qJD(3) * t620 - t615 * t628 + t653;
t589 = g(3) * t637 + t597 * t640;
t619 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t628;
t571 = -qJDD(3) * pkin(3) - pkin(8) * t642 + t616 * t660 - t589;
t545 = -pkin(4) * t584 - pkin(9) * t611 + t614 * t595 + t571;
t538 = -0.2e1 * qJD(6) * t587 + (t586 * t624 - t552) * qJ(6) + (t587 * t624 + t551) * pkin(5) + t545;
t533 = m(7) * t538 + t551 * mrSges(7,1) - t552 * mrSges(7,3) + t586 * t573 - t587 * t576;
t645 = m(6) * t545 + t551 * mrSges(6,1) + t552 * mrSges(6,2) + t586 * t574 + t587 * t575 + t533;
t644 = -m(5) * t571 + t584 * mrSges(5,1) - t585 * mrSges(5,2) + t613 * t591 - t614 * t592 - t645;
t527 = m(4) * t589 + qJDD(3) * mrSges(4,1) - t618 * mrSges(4,3) + qJD(3) * t619 - t615 * t660 + t644;
t512 = t637 * t517 + t640 * t527;
t599 = -qJDD(1) * pkin(1) + t649;
t648 = -m(3) * t599 + (t643 * mrSges(3,3)) - t512;
t510 = m(2) * t621 - (t643 * mrSges(2,2)) + t672 * qJDD(1) + t648;
t598 = t643 * pkin(1) + t675;
t518 = t639 * t522 + t636 * t523;
t647 = -m(4) * t596 + mrSges(4,1) * t617 - t618 * mrSges(4,2) - t619 * t628 - t620 * t660 - t518;
t646 = -m(3) * t598 + (t643 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t647;
t515 = m(2) * t622 - (mrSges(2,1) * t643) - qJDD(1) * mrSges(2,2) + t646;
t665 = t641 * t510 + t638 * t515;
t664 = t677 * t586 - t669 * t587 - t666 * t624;
t663 = t666 * t586 - t668 * t587 + t676 * t624;
t662 = -t669 * t586 + t678 * t587 + t668 * t624;
t655 = -t510 * t638 + t641 * t515;
t654 = t640 * t517 - t637 * t527;
t606 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t640 - Ifges(4,4) * t637) * qJD(1);
t605 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t640 - Ifges(4,2) * t637) * qJD(1);
t604 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t640 - Ifges(4,6) * t637) * qJD(1);
t580 = Ifges(5,1) * t614 + Ifges(5,4) * t613 + Ifges(5,5) * t625;
t579 = Ifges(5,4) * t614 + Ifges(5,2) * t613 + Ifges(5,6) * t625;
t578 = Ifges(5,5) * t614 + Ifges(5,6) * t613 + Ifges(5,3) * t625;
t525 = mrSges(6,2) * t545 + mrSges(7,2) * t536 - mrSges(6,3) * t539 - mrSges(7,3) * t538 - qJ(6) * t533 - t669 * t551 + t678 * t552 + t663 * t586 + t668 * t607 + t664 * t624;
t524 = -mrSges(6,1) * t545 - mrSges(7,1) * t538 + mrSges(7,2) * t535 + mrSges(6,3) * t540 - pkin(5) * t533 - t677 * t551 + t669 * t552 + t663 * t587 + t666 * t607 + t662 * t624;
t511 = -m(3) * g(3) + t654;
t508 = mrSges(5,2) * t571 - mrSges(5,3) * t546 + Ifges(5,1) * t585 + Ifges(5,4) * t584 + Ifges(5,5) * t612 - pkin(9) * t526 - t635 * t524 + t673 * t525 + t613 * t578 - t625 * t579;
t507 = -mrSges(5,1) * t571 + mrSges(5,3) * t547 + Ifges(5,4) * t585 + Ifges(5,2) * t584 + Ifges(5,6) * t612 - pkin(4) * t645 + pkin(9) * t652 + t673 * t524 + t635 * t525 - t614 * t578 + t625 * t580;
t506 = (mrSges(7,2) * qJ(6) + t666) * t551 + (mrSges(7,2) * pkin(5) - t668) * t552 - t604 * t660 + (qJ(6) * t563 - t662) * t586 + (pkin(5) * t563 + t664) * t587 - qJ(6) * t657 - pkin(5) * t651 + t676 * t607 + Ifges(4,6) * qJDD(3) + Ifges(4,4) * t618 - Ifges(5,3) * t612 + t613 * t580 - t614 * t579 + Ifges(4,2) * t617 - mrSges(4,1) * t596 + qJD(3) * t606 + mrSges(4,3) * t590 - Ifges(5,6) * t584 - Ifges(5,5) * t585 + mrSges(5,2) * t547 + mrSges(6,2) * t540 - mrSges(5,1) * t546 - mrSges(6,1) * t539 + mrSges(7,1) * t536 - mrSges(7,3) * t535 - pkin(4) * t526 - pkin(3) * t518;
t505 = mrSges(4,2) * t596 - mrSges(4,3) * t589 + Ifges(4,1) * t618 + Ifges(4,4) * t617 + Ifges(4,5) * qJDD(3) - pkin(8) * t518 - qJD(3) * t605 - t507 * t636 + t508 * t639 - t604 * t628;
t504 = -qJ(2) * t511 - mrSges(2,3) * t621 + pkin(2) * t512 + mrSges(3,1) * t599 + t636 * t508 + t639 * t507 + pkin(3) * t644 + pkin(8) * t653 + mrSges(4,1) * t589 - mrSges(4,2) * t590 + Ifges(4,5) * t618 + Ifges(4,6) * t617 + Ifges(4,3) * qJDD(3) + (t667 * t643) + t670 * qJDD(1) + (t605 * t640 + t606 * t637) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t503 = -mrSges(3,1) * t598 + mrSges(2,3) * t622 - pkin(1) * t511 - pkin(2) * t647 - pkin(7) * t654 + t672 * g(3) - t667 * qJDD(1) - t637 * t505 - t640 * t506 + t670 * t643;
t1 = [-m(1) * g(1) + t655; -m(1) * g(2) + t665; (-m(1) - m(2) - m(3)) * g(3) + t654; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t665 - t638 * t503 + t641 * t504; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t655 + t641 * t503 + t638 * t504; pkin(1) * t648 + qJ(2) * t646 + t640 * t505 - t637 * t506 - pkin(7) * t512 + mrSges(2,1) * t621 - mrSges(2,2) * t622 - mrSges(3,3) * t598 + mrSges(3,2) * t599 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
