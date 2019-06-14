% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% tauJB [(6+6)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 22:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S6PRPPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:00:18
% EndTime: 2019-05-04 22:00:23
% DurationCPUTime: 4.71s
% Computational Cost: add. (66042->252), mult. (114096->314), div. (0->0), fcn. (68201->12), ass. (0->116)
t632 = sin(pkin(10));
t635 = cos(pkin(10));
t617 = t632 * g(1) - t635 * g(2);
t618 = -t635 * g(1) - t632 * g(2);
t628 = -g(3) + qJDD(1);
t633 = sin(pkin(6));
t636 = cos(pkin(6));
t639 = sin(qJ(2));
t642 = cos(qJ(2));
t586 = -t639 * t618 + (t617 * t636 + t628 * t633) * t642;
t644 = qJD(2) ^ 2;
t665 = t636 * t639;
t666 = t633 * t639;
t587 = t617 * t665 + t642 * t618 + t628 * t666;
t653 = qJDD(2) * qJ(3) + (2 * qJD(3) * qJD(2)) + t587;
t671 = -pkin(2) - pkin(3);
t581 = t644 * t671 + t653;
t648 = -t644 * qJ(3) + qJDD(3) - t586;
t583 = qJDD(2) * t671 + t648;
t631 = sin(pkin(11));
t634 = cos(pkin(11));
t577 = t634 * t581 + t631 * t583;
t575 = -t644 * pkin(4) - qJDD(2) * pkin(8) + t577;
t601 = -t633 * t617 + t636 * t628;
t599 = qJDD(4) - t601;
t638 = sin(qJ(5));
t641 = cos(qJ(5));
t572 = t641 * t575 + t638 * t599;
t614 = (pkin(5) * t641 + pkin(9) * t638) * qJD(2);
t643 = qJD(5) ^ 2;
t661 = t641 * qJD(2);
t569 = -t643 * pkin(5) + qJDD(5) * pkin(9) - t614 * t661 + t572;
t576 = -t631 * t581 + t634 * t583;
t574 = qJDD(2) * pkin(4) - t644 * pkin(8) - t576;
t660 = qJD(2) * qJD(5);
t658 = t641 * t660;
t615 = -t638 * qJDD(2) - t658;
t659 = t638 * t660;
t616 = -t641 * qJDD(2) + t659;
t570 = (-t615 + t658) * pkin(9) + (-t616 - t659) * pkin(5) + t574;
t637 = sin(qJ(6));
t640 = cos(qJ(6));
t565 = -t637 * t569 + t640 * t570;
t662 = qJD(2) * t638;
t611 = t640 * qJD(5) + t637 * t662;
t594 = t611 * qJD(6) + t637 * qJDD(5) + t640 * t615;
t612 = t637 * qJD(5) - t640 * t662;
t595 = -t611 * mrSges(7,1) + t612 * mrSges(7,2);
t622 = qJD(6) + t661;
t597 = -t622 * mrSges(7,2) + t611 * mrSges(7,3);
t608 = qJDD(6) - t616;
t563 = m(7) * t565 + t608 * mrSges(7,1) - t594 * mrSges(7,3) - t612 * t595 + t622 * t597;
t566 = t640 * t569 + t637 * t570;
t593 = -t612 * qJD(6) + t640 * qJDD(5) - t637 * t615;
t598 = t622 * mrSges(7,1) - t612 * mrSges(7,3);
t564 = m(7) * t566 - t608 * mrSges(7,2) + t593 * mrSges(7,3) + t611 * t595 - t622 * t598;
t558 = -t637 * t563 + t640 * t564;
t664 = t641 * t599;
t568 = -qJDD(5) * pkin(5) - t643 * pkin(9) - t664 + (-qJD(2) * t614 + t575) * t638;
t588 = Ifges(7,5) * t612 + Ifges(7,6) * t611 + Ifges(7,3) * t622;
t590 = Ifges(7,1) * t612 + Ifges(7,4) * t611 + Ifges(7,5) * t622;
t559 = -mrSges(7,1) * t568 + mrSges(7,3) * t566 + Ifges(7,4) * t594 + Ifges(7,2) * t593 + Ifges(7,6) * t608 - t612 * t588 + t622 * t590;
t589 = Ifges(7,4) * t612 + Ifges(7,2) * t611 + Ifges(7,6) * t622;
t560 = mrSges(7,2) * t568 - mrSges(7,3) * t565 + Ifges(7,1) * t594 + Ifges(7,4) * t593 + Ifges(7,5) * t608 + t611 * t588 - t622 * t589;
t567 = -m(7) * t568 + t593 * mrSges(7,1) - t594 * mrSges(7,2) + t611 * t597 - t612 * t598;
t571 = -t638 * t575 + t664;
t604 = (Ifges(6,6) * qJD(5)) + (-Ifges(6,4) * t638 - Ifges(6,2) * t641) * qJD(2);
t605 = (Ifges(6,5) * qJD(5)) + (-Ifges(6,1) * t638 - Ifges(6,4) * t641) * qJD(2);
t672 = mrSges(6,1) * t571 - mrSges(6,2) * t572 + Ifges(6,5) * t615 + Ifges(6,6) * t616 + Ifges(6,3) * qJDD(5) + pkin(5) * t567 + pkin(9) * t558 - (t638 * t604 - t641 * t605) * qJD(2) + t640 * t559 + t637 * t560;
t670 = -mrSges(3,1) - mrSges(4,1);
t669 = Ifges(4,4) + Ifges(3,5);
t668 = Ifges(3,6) - Ifges(4,6);
t613 = (mrSges(6,1) * t641 - mrSges(6,2) * t638) * qJD(2);
t619 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t662;
t556 = m(6) * t572 - qJDD(5) * mrSges(6,2) + t616 * mrSges(6,3) - qJD(5) * t619 - t613 * t661 + t558;
t620 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t661;
t561 = m(6) * t571 + qJDD(5) * mrSges(6,1) - t615 * mrSges(6,3) + qJD(5) * t620 + t613 * t662 + t567;
t655 = t641 * t556 - t638 * t561;
t548 = m(5) * t577 - t644 * mrSges(5,1) + qJDD(2) * mrSges(5,2) + t655;
t557 = t640 * t563 + t637 * t564;
t647 = -m(6) * t574 + t616 * mrSges(6,1) - t615 * mrSges(6,2) + t619 * t662 - t620 * t661 - t557;
t554 = m(5) * t576 - qJDD(2) * mrSges(5,1) - t644 * mrSges(5,2) + t647;
t544 = t631 * t548 + t634 * t554;
t585 = -qJDD(2) * pkin(2) + t648;
t543 = m(4) * t585 - qJDD(2) * mrSges(4,1) - t644 * mrSges(4,3) + t544;
t542 = m(3) * t586 + qJDD(2) * mrSges(3,1) - t644 * mrSges(3,2) - t543;
t667 = t542 * t642;
t584 = -t644 * pkin(2) + t653;
t656 = t634 * t548 - t631 * t554;
t650 = m(4) * t584 + qJDD(2) * mrSges(4,3) + t656;
t541 = m(3) * t587 - qJDD(2) * mrSges(3,2) + t644 * t670 + t650;
t553 = t638 * t556 + t641 * t561;
t552 = m(5) * t599 + t553;
t551 = m(4) * t601 - t552;
t550 = m(3) * t601 + t551;
t529 = t541 * t665 - t633 * t550 + t636 * t667;
t527 = m(2) * t617 + t529;
t534 = t642 * t541 - t639 * t542;
t533 = m(2) * t618 + t534;
t663 = t635 * t527 + t632 * t533;
t528 = t541 * t666 + t636 * t550 + t633 * t667;
t657 = -t632 * t527 + t635 * t533;
t654 = m(2) * t628 + t528;
t603 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t638 - Ifges(6,6) * t641) * qJD(2);
t545 = mrSges(6,2) * t574 - mrSges(6,3) * t571 + Ifges(6,1) * t615 + Ifges(6,4) * t616 + Ifges(6,5) * qJDD(5) - pkin(9) * t557 - qJD(5) * t604 - t637 * t559 + t640 * t560 - t603 * t661;
t646 = mrSges(7,1) * t565 - mrSges(7,2) * t566 + Ifges(7,5) * t594 + Ifges(7,6) * t593 + Ifges(7,3) * t608 + t612 * t589 - t611 * t590;
t546 = -mrSges(6,1) * t574 + mrSges(6,3) * t572 + Ifges(6,4) * t615 + Ifges(6,2) * t616 + Ifges(6,6) * qJDD(5) - pkin(5) * t557 + qJD(5) * t605 + t603 * t662 - t646;
t530 = mrSges(5,2) * t599 - mrSges(5,3) * t576 - Ifges(5,5) * qJDD(2) - t644 * Ifges(5,6) - pkin(8) * t553 + t641 * t545 - t638 * t546;
t535 = -mrSges(5,1) * t599 + mrSges(5,3) * t577 + t644 * Ifges(5,5) - Ifges(5,6) * qJDD(2) - pkin(4) * t553 - t672;
t523 = mrSges(4,2) * t584 + mrSges(3,3) * t587 - pkin(2) * t551 + pkin(3) * t552 - qJ(4) * t656 + qJDD(2) * t668 - t631 * t530 - t634 * t535 + t601 * t670 + t644 * t669;
t525 = mrSges(4,2) * t585 - mrSges(3,3) * t586 - qJ(3) * t551 - qJ(4) * t544 + t634 * t530 - t631 * t535 - t668 * t644 + (mrSges(3,2) - mrSges(4,3)) * t601 + t669 * qJDD(2);
t649 = pkin(7) * t534 + t523 * t642 + t525 * t639;
t524 = -pkin(2) * t543 + qJ(3) * (-t644 * mrSges(4,1) + t650) + mrSges(3,1) * t586 - mrSges(3,2) * t587 - pkin(3) * t544 - mrSges(4,1) * t585 + mrSges(4,3) * t584 - pkin(8) * t655 - mrSges(5,1) * t576 + mrSges(5,2) * t577 - t638 * t545 - t641 * t546 - pkin(4) * t647 + (Ifges(3,3) + Ifges(4,2) + Ifges(5,3)) * qJDD(2);
t522 = mrSges(2,2) * t628 - mrSges(2,3) * t617 - t639 * t523 + t642 * t525 + (-t528 * t633 - t529 * t636) * pkin(7);
t521 = -mrSges(2,1) * t628 + mrSges(2,3) * t618 - pkin(1) * t528 - t633 * t524 + t636 * t649;
t1 = [-m(1) * g(1) + t657; -m(1) * g(2) + t663; -m(1) * g(3) + t654; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t663 - t632 * t521 + t635 * t522; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t657 + t635 * t521 + t632 * t522; -mrSges(1,1) * g(2) + mrSges(2,1) * t617 + mrSges(1,2) * g(1) - mrSges(2,2) * t618 + pkin(1) * t529 + t636 * t524 + t633 * t649; t654; t524; t543; t552; t672; t646;];
tauJB  = t1;
