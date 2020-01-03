% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRP9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRP9_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP9_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:04:12
% EndTime: 2019-12-31 22:04:19
% DurationCPUTime: 4.18s
% Computational Cost: add. (44124->288), mult. (87480->351), div. (0->0), fcn. (58133->8), ass. (0->112)
t593 = Ifges(5,1) + Ifges(6,1);
t588 = Ifges(5,4) - Ifges(6,5);
t587 = Ifges(6,4) + Ifges(5,5);
t592 = Ifges(5,2) + Ifges(6,3);
t586 = Ifges(5,6) - Ifges(6,6);
t591 = -Ifges(5,3) - Ifges(6,2);
t590 = cos(qJ(4));
t589 = -mrSges(5,3) - mrSges(6,2);
t562 = sin(qJ(1));
t565 = cos(qJ(1));
t552 = -t565 * g(1) - t562 * g(2);
t567 = qJD(1) ^ 2;
t537 = -t567 * pkin(1) + qJDD(1) * pkin(6) + t552;
t561 = sin(qJ(2));
t564 = cos(qJ(2));
t525 = -t561 * g(3) + t564 * t537;
t545 = (-mrSges(3,1) * t564 + mrSges(3,2) * t561) * qJD(1);
t578 = qJD(1) * qJD(2);
t556 = t561 * t578;
t548 = t564 * qJDD(1) - t556;
t580 = qJD(1) * t561;
t549 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t580;
t551 = t562 * g(1) - t565 * g(2);
t536 = -qJDD(1) * pkin(1) - t567 * pkin(6) - t551;
t576 = t564 * t578;
t547 = t561 * qJDD(1) + t576;
t500 = (-t547 - t576) * pkin(7) + (-t548 + t556) * pkin(2) + t536;
t546 = (-pkin(2) * t564 - pkin(7) * t561) * qJD(1);
t566 = qJD(2) ^ 2;
t579 = t564 * qJD(1);
t505 = -t566 * pkin(2) + qJDD(2) * pkin(7) + t546 * t579 + t525;
t560 = sin(qJ(3));
t563 = cos(qJ(3));
t484 = t563 * t500 - t560 * t505;
t543 = t563 * qJD(2) - t560 * t580;
t518 = t543 * qJD(3) + t560 * qJDD(2) + t563 * t547;
t542 = qJDD(3) - t548;
t544 = t560 * qJD(2) + t563 * t580;
t555 = qJD(3) - t579;
t475 = (t543 * t555 - t518) * pkin(8) + (t543 * t544 + t542) * pkin(3) + t484;
t485 = t560 * t500 + t563 * t505;
t517 = -t544 * qJD(3) + t563 * qJDD(2) - t560 * t547;
t526 = t555 * pkin(3) - t544 * pkin(8);
t541 = t543 ^ 2;
t477 = -t541 * pkin(3) + t517 * pkin(8) - t555 * t526 + t485;
t559 = sin(qJ(4));
t473 = t559 * t475 + t590 * t477;
t520 = t559 * t543 + t590 * t544;
t482 = t520 * qJD(4) - t590 * t517 + t559 * t518;
t554 = qJD(4) + t555;
t508 = t554 * mrSges(5,1) - t520 * mrSges(5,3);
t519 = -t590 * t543 + t559 * t544;
t538 = qJDD(4) + t542;
t495 = t519 * pkin(4) - t520 * qJ(5);
t553 = t554 ^ 2;
t468 = -t553 * pkin(4) + t538 * qJ(5) + 0.2e1 * qJD(5) * t554 - t519 * t495 + t473;
t509 = -t554 * mrSges(6,1) + t520 * mrSges(6,2);
t577 = m(6) * t468 + t538 * mrSges(6,3) + t554 * t509;
t496 = t519 * mrSges(6,1) - t520 * mrSges(6,3);
t581 = -t519 * mrSges(5,1) - t520 * mrSges(5,2) - t496;
t463 = m(5) * t473 - t538 * mrSges(5,2) + t589 * t482 - t554 * t508 + t581 * t519 + t577;
t472 = t590 * t475 - t559 * t477;
t483 = -t519 * qJD(4) + t559 * t517 + t590 * t518;
t507 = -t554 * mrSges(5,2) - t519 * mrSges(5,3);
t469 = -t538 * pkin(4) - t553 * qJ(5) + t520 * t495 + qJDD(5) - t472;
t506 = -t519 * mrSges(6,2) + t554 * mrSges(6,3);
t571 = -m(6) * t469 + t538 * mrSges(6,1) + t554 * t506;
t465 = m(5) * t472 + t538 * mrSges(5,1) + t589 * t483 + t554 * t507 + t581 * t520 + t571;
t458 = t559 * t463 + t590 * t465;
t521 = -t543 * mrSges(4,1) + t544 * mrSges(4,2);
t522 = -t555 * mrSges(4,2) + t543 * mrSges(4,3);
t456 = m(4) * t484 + t542 * mrSges(4,1) - t518 * mrSges(4,3) - t544 * t521 + t555 * t522 + t458;
t523 = t555 * mrSges(4,1) - t544 * mrSges(4,3);
t572 = t590 * t463 - t559 * t465;
t457 = m(4) * t485 - t542 * mrSges(4,2) + t517 * mrSges(4,3) + t543 * t521 - t555 * t523 + t572;
t573 = -t560 * t456 + t563 * t457;
t451 = m(3) * t525 - qJDD(2) * mrSges(3,2) + t548 * mrSges(3,3) - qJD(2) * t549 + t545 * t579 + t573;
t524 = -t564 * g(3) - t561 * t537;
t550 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t579;
t504 = -qJDD(2) * pkin(2) - t566 * pkin(7) + t546 * t580 - t524;
t478 = -t517 * pkin(3) - t541 * pkin(8) + t544 * t526 + t504;
t471 = -0.2e1 * qJD(5) * t520 + (t519 * t554 - t483) * qJ(5) + (t520 * t554 + t482) * pkin(4) + t478;
t466 = m(6) * t471 + t482 * mrSges(6,1) - t483 * mrSges(6,3) + t519 * t506 - t520 * t509;
t570 = m(5) * t478 + t482 * mrSges(5,1) + t483 * mrSges(5,2) + t519 * t507 + t520 * t508 + t466;
t568 = -m(4) * t504 + t517 * mrSges(4,1) - t518 * mrSges(4,2) + t543 * t522 - t544 * t523 - t570;
t460 = m(3) * t524 + qJDD(2) * mrSges(3,1) - t547 * mrSges(3,3) + qJD(2) * t550 - t545 * t580 + t568;
t574 = t564 * t451 - t561 * t460;
t445 = m(2) * t552 - t567 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t574;
t452 = t563 * t456 + t560 * t457;
t569 = -m(3) * t536 + t548 * mrSges(3,1) - t547 * mrSges(3,2) - t549 * t580 + t550 * t579 - t452;
t448 = m(2) * t551 + qJDD(1) * mrSges(2,1) - t567 * mrSges(2,2) + t569;
t585 = t562 * t445 + t565 * t448;
t446 = t561 * t451 + t564 * t460;
t584 = t592 * t519 - t588 * t520 - t586 * t554;
t583 = t586 * t519 - t587 * t520 + t591 * t554;
t582 = -t588 * t519 + t593 * t520 + t587 * t554;
t575 = t565 * t445 - t562 * t448;
t535 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t561 + Ifges(3,4) * t564) * qJD(1);
t534 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t561 + Ifges(3,2) * t564) * qJD(1);
t533 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t561 + Ifges(3,6) * t564) * qJD(1);
t513 = Ifges(4,1) * t544 + Ifges(4,4) * t543 + Ifges(4,5) * t555;
t512 = Ifges(4,4) * t544 + Ifges(4,2) * t543 + Ifges(4,6) * t555;
t511 = Ifges(4,5) * t544 + Ifges(4,6) * t543 + Ifges(4,3) * t555;
t454 = mrSges(5,2) * t478 + mrSges(6,2) * t469 - mrSges(5,3) * t472 - mrSges(6,3) * t471 - qJ(5) * t466 - t588 * t482 + t593 * t483 + t583 * t519 + t587 * t538 + t584 * t554;
t453 = -mrSges(5,1) * t478 - mrSges(6,1) * t471 + mrSges(6,2) * t468 + mrSges(5,3) * t473 - pkin(4) * t466 - t592 * t482 + t588 * t483 + t583 * t520 + t586 * t538 + t582 * t554;
t442 = mrSges(4,2) * t504 - mrSges(4,3) * t484 + Ifges(4,1) * t518 + Ifges(4,4) * t517 + Ifges(4,5) * t542 - pkin(8) * t458 - t559 * t453 + t590 * t454 + t543 * t511 - t555 * t512;
t441 = -mrSges(4,1) * t504 + mrSges(4,3) * t485 + Ifges(4,4) * t518 + Ifges(4,2) * t517 + Ifges(4,6) * t542 - pkin(3) * t570 + pkin(8) * t572 + t590 * t453 + t559 * t454 - t544 * t511 + t555 * t513;
t440 = t591 * t538 - pkin(2) * t452 + qJD(2) * t535 - mrSges(3,1) * t536 - Ifges(4,3) * t542 + t543 * t513 - t544 * t512 + Ifges(3,4) * t547 + Ifges(3,2) * t548 + mrSges(3,3) * t525 - mrSges(4,1) * t484 + mrSges(4,2) * t485 + Ifges(3,6) * qJDD(2) - mrSges(6,3) * t468 + mrSges(6,1) * t469 - Ifges(4,6) * t517 - Ifges(4,5) * t518 - pkin(4) * t571 - pkin(3) * t458 - mrSges(5,1) * t472 + mrSges(5,2) * t473 - qJ(5) * t577 - t533 * t580 + (qJ(5) * t496 - t582) * t519 + (pkin(4) * t496 + t584) * t520 + (qJ(5) * mrSges(6,2) + t586) * t482 + (pkin(4) * mrSges(6,2) - t587) * t483;
t439 = mrSges(3,2) * t536 - mrSges(3,3) * t524 + Ifges(3,1) * t547 + Ifges(3,4) * t548 + Ifges(3,5) * qJDD(2) - pkin(7) * t452 - qJD(2) * t534 - t560 * t441 + t563 * t442 + t533 * t579;
t438 = Ifges(2,6) * qJDD(1) + t567 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t552 - Ifges(3,5) * t547 - Ifges(3,6) * t548 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t524 + mrSges(3,2) * t525 - t560 * t442 - t563 * t441 - pkin(2) * t568 - pkin(7) * t573 - pkin(1) * t446 + (-t561 * t534 + t564 * t535) * qJD(1);
t437 = -mrSges(2,2) * g(3) - mrSges(2,3) * t551 + Ifges(2,5) * qJDD(1) - t567 * Ifges(2,6) - pkin(6) * t446 + t564 * t439 - t561 * t440;
t1 = [-m(1) * g(1) + t575; -m(1) * g(2) + t585; (-m(1) - m(2)) * g(3) + t446; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t585 + t565 * t437 - t562 * t438; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t575 + t562 * t437 + t565 * t438; -mrSges(1,1) * g(2) + mrSges(2,1) * t551 + mrSges(1,2) * g(1) - mrSges(2,2) * t552 + Ifges(2,3) * qJDD(1) + pkin(1) * t569 + pkin(6) * t574 + t561 * t439 + t564 * t440;];
tauB = t1;
