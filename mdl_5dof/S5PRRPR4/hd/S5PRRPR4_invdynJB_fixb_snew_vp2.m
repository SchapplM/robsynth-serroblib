% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PRRPR4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:22:04
% EndTime: 2019-12-05 16:22:12
% DurationCPUTime: 4.90s
% Computational Cost: add. (53775->260), mult. (118226->335), div. (0->0), fcn. (79032->10), ass. (0->107)
t631 = sin(pkin(8));
t656 = cos(pkin(8));
t618 = -t656 * g(1) - t631 * g(2);
t629 = -g(3) + qJDD(1);
t635 = sin(qJ(2));
t638 = cos(qJ(2));
t598 = t638 * t618 + t635 * t629;
t639 = qJD(2) ^ 2;
t593 = -pkin(2) * t639 + qJDD(2) * pkin(6) + t598;
t617 = g(1) * t631 - t656 * g(2);
t634 = sin(qJ(3));
t637 = cos(qJ(3));
t579 = -t634 * t593 - t637 * t617;
t652 = qJD(2) * qJD(3);
t650 = t637 * t652;
t615 = qJDD(2) * t634 + t650;
t574 = (-t615 + t650) * qJ(4) + (t634 * t637 * t639 + qJDD(3)) * pkin(3) + t579;
t580 = t637 * t593 - t634 * t617;
t616 = qJDD(2) * t637 - t634 * t652;
t654 = qJD(2) * t634;
t619 = qJD(3) * pkin(3) - qJ(4) * t654;
t628 = t637 ^ 2;
t575 = -pkin(3) * t628 * t639 + qJ(4) * t616 - qJD(3) * t619 + t580;
t630 = sin(pkin(9));
t632 = cos(pkin(9));
t603 = (t630 * t637 + t632 * t634) * qJD(2);
t554 = -0.2e1 * qJD(4) * t603 + t632 * t574 - t630 * t575;
t590 = t615 * t632 + t616 * t630;
t602 = (-t630 * t634 + t632 * t637) * qJD(2);
t552 = (qJD(3) * t602 - t590) * pkin(7) + (t602 * t603 + qJDD(3)) * pkin(4) + t554;
t555 = 0.2e1 * qJD(4) * t602 + t630 * t574 + t632 * t575;
t589 = -t615 * t630 + t616 * t632;
t596 = qJD(3) * pkin(4) - pkin(7) * t603;
t601 = t602 ^ 2;
t553 = -pkin(4) * t601 + pkin(7) * t589 - qJD(3) * t596 + t555;
t633 = sin(qJ(5));
t636 = cos(qJ(5));
t550 = t552 * t636 - t553 * t633;
t584 = t602 * t636 - t603 * t633;
t564 = t584 * qJD(5) + t589 * t633 + t590 * t636;
t585 = t602 * t633 + t603 * t636;
t570 = -t584 * mrSges(6,1) + mrSges(6,2) * t585;
t627 = qJD(3) + qJD(5);
t577 = -mrSges(6,2) * t627 + t584 * mrSges(6,3);
t626 = qJDD(3) + qJDD(5);
t546 = m(6) * t550 + mrSges(6,1) * t626 - t564 * mrSges(6,3) - t570 * t585 + t577 * t627;
t551 = t552 * t633 + t553 * t636;
t563 = -qJD(5) * t585 + t589 * t636 - t590 * t633;
t578 = mrSges(6,1) * t627 - mrSges(6,3) * t585;
t547 = m(6) * t551 - mrSges(6,2) * t626 + t563 * mrSges(6,3) + t584 * t570 - t578 * t627;
t537 = t636 * t546 + t633 * t547;
t587 = -mrSges(5,1) * t602 + mrSges(5,2) * t603;
t594 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t602;
t535 = m(5) * t554 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t590 + qJD(3) * t594 - t587 * t603 + t537;
t595 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t603;
t646 = -t546 * t633 + t636 * t547;
t536 = m(5) * t555 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t589 - qJD(3) * t595 + t587 * t602 + t646;
t531 = t632 * t535 + t630 * t536;
t582 = Ifges(5,4) * t603 + Ifges(5,2) * t602 + Ifges(5,6) * qJD(3);
t583 = Ifges(5,1) * t603 + Ifges(5,4) * t602 + Ifges(5,5) * qJD(3);
t605 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t634 + Ifges(4,2) * t637) * qJD(2);
t606 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t634 + Ifges(4,4) * t637) * qJD(2);
t566 = Ifges(6,4) * t585 + Ifges(6,2) * t584 + Ifges(6,6) * t627;
t567 = Ifges(6,1) * t585 + Ifges(6,4) * t584 + Ifges(6,5) * t627;
t642 = -mrSges(6,1) * t550 + mrSges(6,2) * t551 - Ifges(6,5) * t564 - Ifges(6,6) * t563 - Ifges(6,3) * t626 - t585 * t566 + t584 * t567;
t658 = mrSges(4,1) * t579 + mrSges(5,1) * t554 - mrSges(4,2) * t580 - mrSges(5,2) * t555 + Ifges(4,5) * t615 + Ifges(5,5) * t590 + Ifges(4,6) * t616 + Ifges(5,6) * t589 + pkin(3) * t531 + pkin(4) * t537 + (t605 * t634 - t606 * t637) * qJD(2) + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t603 * t582 - t602 * t583 - t642;
t614 = (-mrSges(4,1) * t637 + mrSges(4,2) * t634) * qJD(2);
t653 = qJD(2) * t637;
t621 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t653;
t529 = m(4) * t579 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t615 + qJD(3) * t621 - t614 * t654 + t531;
t620 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t654;
t647 = -t630 * t535 + t632 * t536;
t530 = m(4) * t580 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t616 - qJD(3) * t620 + t614 * t653 + t647;
t525 = -t529 * t634 + t637 * t530;
t521 = m(3) * t598 - mrSges(3,1) * t639 - qJDD(2) * mrSges(3,2) + t525;
t597 = -t635 * t618 + t629 * t638;
t643 = -qJDD(2) * pkin(2) - t597;
t576 = -pkin(3) * t616 + qJDD(4) + t619 * t654 + (-qJ(4) * t628 - pkin(6)) * t639 + t643;
t557 = -pkin(4) * t589 - pkin(7) * t601 + t596 * t603 + t576;
t645 = m(6) * t557 - t563 * mrSges(6,1) + t564 * mrSges(6,2) - t584 * t577 + t585 * t578;
t548 = m(5) * t576 - t589 * mrSges(5,1) + mrSges(5,2) * t590 - t602 * t594 + t595 * t603 + t645;
t592 = -pkin(6) * t639 + t643;
t542 = -m(4) * t592 + t616 * mrSges(4,1) - mrSges(4,2) * t615 - t620 * t654 + t621 * t653 - t548;
t541 = m(3) * t597 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t639 + t542;
t648 = t638 * t521 - t541 * t635;
t517 = m(2) * t618 + t648;
t524 = t637 * t529 + t634 * t530;
t523 = (m(2) + m(3)) * t617 - t524;
t655 = t631 * t517 + t656 * t523;
t518 = t635 * t521 + t638 * t541;
t651 = m(2) * t629 + t518;
t649 = t656 * t517 - t523 * t631;
t565 = Ifges(6,5) * t585 + Ifges(6,6) * t584 + Ifges(6,3) * t627;
t538 = -mrSges(6,1) * t557 + mrSges(6,3) * t551 + Ifges(6,4) * t564 + Ifges(6,2) * t563 + Ifges(6,6) * t626 - t565 * t585 + t627 * t567;
t539 = mrSges(6,2) * t557 - mrSges(6,3) * t550 + Ifges(6,1) * t564 + Ifges(6,4) * t563 + Ifges(6,5) * t626 + t584 * t565 - t566 * t627;
t581 = Ifges(5,5) * t603 + Ifges(5,6) * t602 + Ifges(5,3) * qJD(3);
t526 = -mrSges(5,1) * t576 + mrSges(5,3) * t555 + Ifges(5,4) * t590 + Ifges(5,2) * t589 + Ifges(5,6) * qJDD(3) - pkin(4) * t645 + pkin(7) * t646 + qJD(3) * t583 + t636 * t538 + t633 * t539 - t603 * t581;
t527 = mrSges(5,2) * t576 - mrSges(5,3) * t554 + Ifges(5,1) * t590 + Ifges(5,4) * t589 + Ifges(5,5) * qJDD(3) - pkin(7) * t537 - qJD(3) * t582 - t538 * t633 + t539 * t636 + t581 * t602;
t604 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t634 + Ifges(4,6) * t637) * qJD(2);
t513 = -mrSges(4,1) * t592 + mrSges(4,3) * t580 + Ifges(4,4) * t615 + Ifges(4,2) * t616 + Ifges(4,6) * qJDD(3) - pkin(3) * t548 + qJ(4) * t647 + qJD(3) * t606 + t632 * t526 + t630 * t527 - t604 * t654;
t514 = mrSges(4,2) * t592 - mrSges(4,3) * t579 + Ifges(4,1) * t615 + Ifges(4,4) * t616 + Ifges(4,5) * qJDD(3) - qJ(4) * t531 - qJD(3) * t605 - t526 * t630 + t527 * t632 + t604 * t653;
t641 = mrSges(3,1) * t597 - mrSges(3,2) * t598 + Ifges(3,3) * qJDD(2) + pkin(2) * t542 + pkin(6) * t525 + t513 * t637 + t514 * t634;
t512 = mrSges(3,1) * t617 + mrSges(3,3) * t598 + t639 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t524 - t658;
t511 = -mrSges(3,2) * t617 - mrSges(3,3) * t597 + Ifges(3,5) * qJDD(2) - Ifges(3,6) * t639 - pkin(6) * t524 - t513 * t634 + t514 * t637;
t510 = -mrSges(2,1) * t629 + mrSges(2,3) * t618 - pkin(1) * t518 - t641;
t509 = mrSges(2,2) * t629 - mrSges(2,3) * t617 - pkin(5) * t518 + t511 * t638 - t512 * t635;
t1 = [-m(1) * g(1) + t649; -m(1) * g(2) + t655; -m(1) * g(3) + t651; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t655 + t656 * t509 - t631 * t510; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t649 + t631 * t509 + t656 * t510; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t617 - mrSges(2,2) * t618 + t635 * t511 + t638 * t512 + pkin(1) * (m(3) * t617 - t524) + pkin(5) * t648; t651; t641; t658; t548; -t642;];
tauJB = t1;
