% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RPRRR7_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:37
% EndTime: 2019-12-31 19:03:41
% DurationCPUTime: 4.82s
% Computational Cost: add. (59579->270), mult. (114389->338), div. (0->0), fcn. (71793->10), ass. (0->112)
t640 = sin(qJ(1));
t644 = cos(qJ(1));
t624 = t640 * g(1) - t644 * g(2);
t615 = qJDD(1) * pkin(1) + t624;
t625 = -t644 * g(1) - t640 * g(2);
t646 = qJD(1) ^ 2;
t617 = -t646 * pkin(1) + t625;
t635 = sin(pkin(9));
t636 = cos(pkin(9));
t594 = t636 * t615 - t635 * t617;
t582 = -qJDD(1) * pkin(2) - t646 * pkin(6) - t594;
t639 = sin(qJ(3));
t643 = cos(qJ(3));
t659 = qJD(1) * qJD(3);
t658 = t643 * t659;
t619 = t639 * qJDD(1) + t658;
t629 = t639 * t659;
t620 = t643 * qJDD(1) - t629;
t570 = (-t619 - t658) * pkin(7) + (-t620 + t629) * pkin(3) + t582;
t595 = t635 * t615 + t636 * t617;
t583 = -t646 * pkin(2) + qJDD(1) * pkin(6) + t595;
t634 = -g(3) + qJDD(2);
t577 = t643 * t583 + t639 * t634;
t618 = (-pkin(3) * t643 - pkin(7) * t639) * qJD(1);
t645 = qJD(3) ^ 2;
t660 = t643 * qJD(1);
t574 = -t645 * pkin(3) + qJDD(3) * pkin(7) + t618 * t660 + t577;
t638 = sin(qJ(4));
t642 = cos(qJ(4));
t557 = t642 * t570 - t638 * t574;
t661 = qJD(1) * t639;
t613 = t642 * qJD(3) - t638 * t661;
t590 = t613 * qJD(4) + t638 * qJDD(3) + t642 * t619;
t612 = qJDD(4) - t620;
t614 = t638 * qJD(3) + t642 * t661;
t627 = qJD(4) - t660;
t554 = (t613 * t627 - t590) * pkin(8) + (t613 * t614 + t612) * pkin(4) + t557;
t558 = t638 * t570 + t642 * t574;
t589 = -t614 * qJD(4) + t642 * qJDD(3) - t638 * t619;
t599 = t627 * pkin(4) - t614 * pkin(8);
t611 = t613 ^ 2;
t555 = -t611 * pkin(4) + t589 * pkin(8) - t627 * t599 + t558;
t637 = sin(qJ(5));
t641 = cos(qJ(5));
t553 = t637 * t554 + t641 * t555;
t576 = -t639 * t583 + t643 * t634;
t573 = -qJDD(3) * pkin(3) - t645 * pkin(7) + t618 * t661 - t576;
t556 = -t589 * pkin(4) - t611 * pkin(8) + t614 * t599 + t573;
t593 = t637 * t613 + t641 * t614;
t563 = -t593 * qJD(5) + t641 * t589 - t637 * t590;
t592 = t641 * t613 - t637 * t614;
t564 = t592 * qJD(5) + t637 * t589 + t641 * t590;
t626 = qJD(5) + t627;
t567 = Ifges(6,5) * t593 + Ifges(6,6) * t592 + Ifges(6,3) * t626;
t569 = Ifges(6,1) * t593 + Ifges(6,4) * t592 + Ifges(6,5) * t626;
t608 = qJDD(5) + t612;
t541 = -mrSges(6,1) * t556 + mrSges(6,3) * t553 + Ifges(6,4) * t564 + Ifges(6,2) * t563 + Ifges(6,6) * t608 - t593 * t567 + t626 * t569;
t552 = t641 * t554 - t637 * t555;
t568 = Ifges(6,4) * t593 + Ifges(6,2) * t592 + Ifges(6,6) * t626;
t542 = mrSges(6,2) * t556 - mrSges(6,3) * t552 + Ifges(6,1) * t564 + Ifges(6,4) * t563 + Ifges(6,5) * t608 + t592 * t567 - t626 * t568;
t584 = Ifges(5,5) * t614 + Ifges(5,6) * t613 + Ifges(5,3) * t627;
t586 = Ifges(5,1) * t614 + Ifges(5,4) * t613 + Ifges(5,5) * t627;
t578 = -t626 * mrSges(6,2) + t592 * mrSges(6,3);
t579 = t626 * mrSges(6,1) - t593 * mrSges(6,3);
t652 = m(6) * t556 - t563 * mrSges(6,1) + t564 * mrSges(6,2) - t592 * t578 + t593 * t579;
t575 = -t592 * mrSges(6,1) + t593 * mrSges(6,2);
t548 = m(6) * t552 + t608 * mrSges(6,1) - t564 * mrSges(6,3) - t593 * t575 + t626 * t578;
t549 = m(6) * t553 - t608 * mrSges(6,2) + t563 * mrSges(6,3) + t592 * t575 - t626 * t579;
t654 = -t637 * t548 + t641 * t549;
t521 = -mrSges(5,1) * t573 + mrSges(5,3) * t558 + Ifges(5,4) * t590 + Ifges(5,2) * t589 + Ifges(5,6) * t612 - pkin(4) * t652 + pkin(8) * t654 + t641 * t541 + t637 * t542 - t614 * t584 + t627 * t586;
t540 = t641 * t548 + t637 * t549;
t585 = Ifges(5,4) * t614 + Ifges(5,2) * t613 + Ifges(5,6) * t627;
t528 = mrSges(5,2) * t573 - mrSges(5,3) * t557 + Ifges(5,1) * t590 + Ifges(5,4) * t589 + Ifges(5,5) * t612 - pkin(8) * t540 - t637 * t541 + t641 * t542 + t613 * t584 - t627 * t585;
t596 = -t613 * mrSges(5,1) + t614 * mrSges(5,2);
t597 = -t627 * mrSges(5,2) + t613 * mrSges(5,3);
t538 = m(5) * t557 + t612 * mrSges(5,1) - t590 * mrSges(5,3) - t614 * t596 + t627 * t597 + t540;
t598 = t627 * mrSges(5,1) - t614 * mrSges(5,3);
t539 = m(5) * t558 - t612 * mrSges(5,2) + t589 * mrSges(5,3) + t613 * t596 - t627 * t598 + t654;
t536 = -t638 * t538 + t642 * t539;
t550 = -m(5) * t573 + t589 * mrSges(5,1) - t590 * mrSges(5,2) + t613 * t597 - t614 * t598 - t652;
t606 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t639 + Ifges(4,2) * t643) * qJD(1);
t607 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t639 + Ifges(4,4) * t643) * qJD(1);
t663 = mrSges(4,1) * t576 - mrSges(4,2) * t577 + Ifges(4,5) * t619 + Ifges(4,6) * t620 + Ifges(4,3) * qJDD(3) + pkin(3) * t550 + pkin(7) * t536 + t642 * t521 + t638 * t528 + (t639 * t606 - t643 * t607) * qJD(1);
t616 = (-mrSges(4,1) * t643 + mrSges(4,2) * t639) * qJD(1);
t622 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t661;
t534 = m(4) * t577 - qJDD(3) * mrSges(4,2) + t620 * mrSges(4,3) - qJD(3) * t622 + t616 * t660 + t536;
t623 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t660;
t544 = m(4) * t576 + qJDD(3) * mrSges(4,1) - t619 * mrSges(4,3) + qJD(3) * t623 - t616 * t661 + t550;
t655 = t643 * t534 - t639 * t544;
t524 = m(3) * t595 - t646 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t655;
t535 = t642 * t538 + t638 * t539;
t649 = -m(4) * t582 + t620 * mrSges(4,1) - t619 * mrSges(4,2) - t622 * t661 + t623 * t660 - t535;
t530 = m(3) * t594 + qJDD(1) * mrSges(3,1) - t646 * mrSges(3,2) + t649;
t518 = t635 * t524 + t636 * t530;
t515 = m(2) * t624 + qJDD(1) * mrSges(2,1) - t646 * mrSges(2,2) + t518;
t656 = t636 * t524 - t635 * t530;
t516 = m(2) * t625 - t646 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t656;
t662 = t644 * t515 + t640 * t516;
t527 = t639 * t534 + t643 * t544;
t525 = m(3) * t634 + t527;
t657 = -t640 * t515 + t644 * t516;
t651 = -mrSges(6,1) * t552 + mrSges(6,2) * t553 - Ifges(6,5) * t564 - Ifges(6,6) * t563 - Ifges(6,3) * t608 - t593 * t568 + t592 * t569;
t605 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t639 + Ifges(4,6) * t643) * qJD(1);
t511 = mrSges(4,2) * t582 - mrSges(4,3) * t576 + Ifges(4,1) * t619 + Ifges(4,4) * t620 + Ifges(4,5) * qJDD(3) - pkin(7) * t535 - qJD(3) * t606 - t638 * t521 + t642 * t528 + t605 * t660;
t647 = mrSges(5,1) * t557 - mrSges(5,2) * t558 + Ifges(5,5) * t590 + Ifges(5,6) * t589 + Ifges(5,3) * t612 + pkin(4) * t540 + t614 * t585 - t613 * t586 - t651;
t520 = -mrSges(4,1) * t582 + mrSges(4,3) * t577 + Ifges(4,4) * t619 + Ifges(4,2) * t620 + Ifges(4,6) * qJDD(3) - pkin(3) * t535 + qJD(3) * t607 - t605 * t661 - t647;
t650 = mrSges(2,1) * t624 + mrSges(3,1) * t594 - mrSges(2,2) * t625 - mrSges(3,2) * t595 + pkin(1) * t518 + pkin(2) * t649 + pkin(6) * t655 + t639 * t511 + t643 * t520 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t509 = -mrSges(3,1) * t634 + mrSges(3,3) * t595 + t646 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t527 - t663;
t508 = mrSges(3,2) * t634 - mrSges(3,3) * t594 + Ifges(3,5) * qJDD(1) - t646 * Ifges(3,6) - pkin(6) * t527 + t643 * t511 - t639 * t520;
t507 = -mrSges(2,2) * g(3) - mrSges(2,3) * t624 + Ifges(2,5) * qJDD(1) - t646 * Ifges(2,6) - qJ(2) * t518 + t636 * t508 - t635 * t509;
t506 = mrSges(2,1) * g(3) + mrSges(2,3) * t625 + t646 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t525 + qJ(2) * t656 + t635 * t508 + t636 * t509;
t1 = [-m(1) * g(1) + t657; -m(1) * g(2) + t662; (-m(1) - m(2)) * g(3) + t525; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t662 - t640 * t506 + t644 * t507; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t657 + t644 * t506 + t640 * t507; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t650; t650; t525; t663; t647; -t651;];
tauJB = t1;
