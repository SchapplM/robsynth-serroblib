% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPPR1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:23
% EndTime: 2019-12-05 18:18:26
% DurationCPUTime: 2.88s
% Computational Cost: add. (45658->207), mult. (62159->259), div. (0->0), fcn. (35962->10), ass. (0->98)
t579 = qJD(1) + qJD(2);
t575 = t579 ^ 2;
t584 = cos(pkin(9));
t617 = pkin(4) * t584;
t582 = sin(pkin(9));
t616 = mrSges(5,2) * t582;
t578 = t584 ^ 2;
t615 = t575 * t578;
t576 = qJDD(1) + qJDD(2);
t614 = t576 * t584;
t603 = Ifges(5,5) * t582 + Ifges(5,6) * t584;
t613 = t575 * t603;
t588 = sin(qJ(1));
t591 = cos(qJ(1));
t563 = t591 * g(2) + t588 * g(3);
t557 = qJDD(1) * pkin(1) + t563;
t562 = t588 * g(2) - t591 * g(3);
t592 = qJD(1) ^ 2;
t558 = -t592 * pkin(1) + t562;
t587 = sin(qJ(2));
t590 = cos(qJ(2));
t546 = t590 * t557 - t587 * t558;
t543 = t576 * pkin(2) + t546;
t547 = t587 * t557 + t590 * t558;
t544 = -t575 * pkin(2) + t547;
t583 = sin(pkin(8));
t585 = cos(pkin(8));
t531 = t583 * t543 + t585 * t544;
t528 = -t575 * pkin(3) + t576 * qJ(4) + t531;
t581 = -g(1) + qJDD(3);
t611 = qJD(4) * t579;
t612 = t584 * t581 - 0.2e1 * t582 * t611;
t521 = (-pkin(7) * t576 + t575 * t617 - t528) * t582 + t612;
t525 = t582 * t581 + (t528 + 0.2e1 * t611) * t584;
t522 = -pkin(4) * t615 + pkin(7) * t614 + t525;
t586 = sin(qJ(5));
t589 = cos(qJ(5));
t519 = t589 * t521 - t586 * t522;
t598 = -t582 * t586 + t584 * t589;
t550 = t598 * t579;
t599 = t582 * t589 + t584 * t586;
t551 = t599 * t579;
t537 = -t550 * mrSges(6,1) + t551 * mrSges(6,2);
t539 = t550 * qJD(5) + t599 * t576;
t548 = -qJD(5) * mrSges(6,2) + t550 * mrSges(6,3);
t517 = m(6) * t519 + qJDD(5) * mrSges(6,1) - t539 * mrSges(6,3) + qJD(5) * t548 - t551 * t537;
t520 = t586 * t521 + t589 * t522;
t538 = -t551 * qJD(5) + t598 * t576;
t549 = qJD(5) * mrSges(6,1) - t551 * mrSges(6,3);
t518 = m(6) * t520 - qJDD(5) * mrSges(6,2) + t538 * mrSges(6,3) - qJD(5) * t549 + t550 * t537;
t507 = t589 * t517 + t586 * t518;
t524 = -t582 * t528 + t612;
t601 = mrSges(5,3) * t576 + (-mrSges(5,1) * t584 + t616) * t575;
t505 = m(5) * t524 - t601 * t582 + t507;
t606 = -t586 * t517 + t589 * t518;
t506 = m(5) * t525 + t601 * t584 + t606;
t607 = -t582 * t505 + t584 * t506;
t498 = m(4) * t531 - t575 * mrSges(4,1) - t576 * mrSges(4,2) + t607;
t530 = t585 * t543 - t583 * t544;
t602 = qJDD(4) - t530;
t527 = -t576 * pkin(3) - t575 * qJ(4) + t602;
t577 = t582 ^ 2;
t523 = (-pkin(3) - t617) * t576 + (-qJ(4) + (-t577 - t578) * pkin(7)) * t575 + t602;
t597 = m(6) * t523 - t538 * mrSges(6,1) + t539 * mrSges(6,2) - t550 * t548 + t551 * t549;
t595 = -m(5) * t527 + mrSges(5,1) * t614 - t597 + (t575 * t577 + t615) * mrSges(5,3);
t511 = m(4) * t530 - t575 * mrSges(4,2) + (mrSges(4,1) - t616) * t576 + t595;
t493 = t583 * t498 + t585 * t511;
t488 = m(3) * t546 + t576 * mrSges(3,1) - t575 * mrSges(3,2) + t493;
t608 = t585 * t498 - t583 * t511;
t489 = m(3) * t547 - t575 * mrSges(3,1) - t576 * mrSges(3,2) + t608;
t483 = t590 * t488 + t587 * t489;
t501 = t584 * t505 + t582 * t506;
t499 = m(4) * t581 + t501;
t609 = -t587 * t488 + t590 * t489;
t480 = m(2) * t562 - t592 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t609;
t481 = m(2) * t563 + qJDD(1) * mrSges(2,1) - t592 * mrSges(2,2) + t483;
t610 = t591 * t480 - t588 * t481;
t605 = Ifges(5,1) * t582 + Ifges(5,4) * t584;
t604 = Ifges(5,4) * t582 + Ifges(5,2) * t584;
t600 = -t588 * t480 - t591 * t481;
t532 = Ifges(6,5) * t551 + Ifges(6,6) * t550 + Ifges(6,3) * qJD(5);
t534 = Ifges(6,1) * t551 + Ifges(6,4) * t550 + Ifges(6,5) * qJD(5);
t508 = -mrSges(6,1) * t523 + mrSges(6,3) * t520 + Ifges(6,4) * t539 + Ifges(6,2) * t538 + Ifges(6,6) * qJDD(5) + qJD(5) * t534 - t551 * t532;
t533 = Ifges(6,4) * t551 + Ifges(6,2) * t550 + Ifges(6,6) * qJD(5);
t509 = mrSges(6,2) * t523 - mrSges(6,3) * t519 + Ifges(6,1) * t539 + Ifges(6,4) * t538 + Ifges(6,5) * qJDD(5) - qJD(5) * t533 + t550 * t532;
t491 = -mrSges(5,1) * t527 + mrSges(5,3) * t525 - pkin(4) * t597 + pkin(7) * t606 + t589 * t508 + t586 * t509 + t604 * t576 - t582 * t613;
t495 = mrSges(5,2) * t527 - mrSges(5,3) * t524 - pkin(7) * t507 - t586 * t508 + t589 * t509 + t605 * t576 + t584 * t613;
t513 = t576 * t616 - t595;
t596 = mrSges(3,1) * t546 + mrSges(4,1) * t530 - mrSges(3,2) * t547 - mrSges(4,2) * t531 + pkin(2) * t493 - pkin(3) * t513 + qJ(4) * t607 + t584 * t491 + t582 * t495 + (Ifges(3,3) + Ifges(4,3)) * t576;
t594 = mrSges(6,1) * t519 - mrSges(6,2) * t520 + Ifges(6,5) * t539 + Ifges(6,6) * t538 + Ifges(6,3) * qJDD(5) + t551 * t533 - t550 * t534;
t593 = mrSges(2,1) * t563 - mrSges(2,2) * t562 + Ifges(2,3) * qJDD(1) + pkin(1) * t483 + t596;
t484 = -mrSges(4,1) * t581 - mrSges(5,1) * t524 + mrSges(5,2) * t525 + mrSges(4,3) * t531 - pkin(3) * t501 - pkin(4) * t507 + (Ifges(4,6) - t603) * t576 - t594 + (-t582 * t604 + t584 * t605 + Ifges(4,5)) * t575;
t478 = mrSges(4,2) * t581 - mrSges(4,3) * t530 + Ifges(4,5) * t576 - t575 * Ifges(4,6) - qJ(4) * t501 - t582 * t491 + t584 * t495;
t477 = -mrSges(3,2) * g(1) - mrSges(3,3) * t546 + Ifges(3,5) * t576 - t575 * Ifges(3,6) - qJ(3) * t493 + t585 * t478 - t583 * t484;
t476 = mrSges(3,1) * g(1) + mrSges(3,3) * t547 + t575 * Ifges(3,5) + Ifges(3,6) * t576 - pkin(2) * t499 + qJ(3) * t608 + t583 * t478 + t585 * t484;
t475 = -mrSges(2,2) * g(1) - mrSges(2,3) * t563 + Ifges(2,5) * qJDD(1) - t592 * Ifges(2,6) - pkin(6) * t483 - t587 * t476 + t590 * t477;
t474 = Ifges(2,6) * qJDD(1) + t592 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t562 + t587 * t477 + t590 * t476 - pkin(1) * (-m(3) * g(1) + t499) + pkin(6) * t609;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t499; -m(1) * g(2) + t600; -m(1) * g(3) + t610; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t593; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t610 - t591 * t474 - t588 * t475; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t600 - t588 * t474 + t591 * t475; t593; t596; t499; t513; t594;];
tauJB = t1;
