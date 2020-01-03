% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:13:30
% EndTime: 2019-12-31 21:13:37
% DurationCPUTime: 6.30s
% Computational Cost: add. (93725->314), mult. (210036->400), div. (0->0), fcn. (148948->10), ass. (0->123)
t601 = 2 * qJD(4);
t583 = qJD(1) ^ 2;
t600 = pkin(2) * t583;
t578 = sin(qJ(1));
t582 = cos(qJ(1));
t567 = -t582 * g(1) - t578 * g(2);
t556 = -t583 * pkin(1) + qJDD(1) * pkin(6) + t567;
t577 = sin(qJ(2));
t599 = t577 * t556;
t581 = cos(qJ(2));
t595 = qJD(1) * qJD(2);
t561 = t577 * qJDD(1) + t581 * t595;
t523 = qJDD(2) * pkin(2) - t561 * pkin(7) - t599 + (pkin(7) * t595 + t577 * t600 - g(3)) * t581;
t544 = -t577 * g(3) + t581 * t556;
t562 = t581 * qJDD(1) - t577 * t595;
t597 = qJD(1) * t577;
t565 = qJD(2) * pkin(2) - pkin(7) * t597;
t572 = t581 ^ 2;
t524 = t562 * pkin(7) - qJD(2) * t565 - t572 * t600 + t544;
t576 = sin(qJ(3));
t580 = cos(qJ(3));
t503 = t580 * t523 - t576 * t524;
t553 = (-t576 * t577 + t580 * t581) * qJD(1);
t529 = t553 * qJD(3) + t580 * t561 + t576 * t562;
t554 = (t576 * t581 + t577 * t580) * qJD(1);
t570 = qJDD(2) + qJDD(3);
t571 = qJD(2) + qJD(3);
t491 = (t553 * t571 - t529) * qJ(4) + (t553 * t554 + t570) * pkin(3) + t503;
t504 = t576 * t523 + t580 * t524;
t528 = -t554 * qJD(3) - t576 * t561 + t580 * t562;
t546 = t571 * pkin(3) - t554 * qJ(4);
t549 = t553 ^ 2;
t493 = -t549 * pkin(3) + t528 * qJ(4) - t571 * t546 + t504;
t573 = sin(pkin(9));
t574 = cos(pkin(9));
t540 = t574 * t553 - t573 * t554;
t488 = t573 * t491 + t574 * t493 + t540 * t601;
t508 = t574 * t528 - t573 * t529;
t541 = t573 * t553 + t574 * t554;
t518 = -t540 * mrSges(5,1) + t541 * mrSges(5,2);
t532 = t571 * mrSges(5,1) - t541 * mrSges(5,3);
t519 = -t540 * pkin(4) - t541 * pkin(8);
t569 = t571 ^ 2;
t486 = -t569 * pkin(4) + t570 * pkin(8) + t540 * t519 + t488;
t566 = t578 * g(1) - t582 * g(2);
t588 = -qJDD(1) * pkin(1) - t566;
t530 = -t562 * pkin(2) + t565 * t597 + (-pkin(7) * t572 - pkin(6)) * t583 + t588;
t498 = -t528 * pkin(3) - t549 * qJ(4) + t554 * t546 + qJDD(4) + t530;
t509 = t573 * t528 + t574 * t529;
t489 = (-t540 * t571 - t509) * pkin(8) + (t541 * t571 - t508) * pkin(4) + t498;
t575 = sin(qJ(5));
t579 = cos(qJ(5));
t483 = -t575 * t486 + t579 * t489;
t526 = -t575 * t541 + t579 * t571;
t496 = t526 * qJD(5) + t579 * t509 + t575 * t570;
t507 = qJDD(5) - t508;
t527 = t579 * t541 + t575 * t571;
t510 = -t526 * mrSges(6,1) + t527 * mrSges(6,2);
t534 = qJD(5) - t540;
t511 = -t534 * mrSges(6,2) + t526 * mrSges(6,3);
t481 = m(6) * t483 + t507 * mrSges(6,1) - t496 * mrSges(6,3) - t527 * t510 + t534 * t511;
t484 = t579 * t486 + t575 * t489;
t495 = -t527 * qJD(5) - t575 * t509 + t579 * t570;
t512 = t534 * mrSges(6,1) - t527 * mrSges(6,3);
t482 = m(6) * t484 - t507 * mrSges(6,2) + t495 * mrSges(6,3) + t526 * t510 - t534 * t512;
t590 = -t575 * t481 + t579 * t482;
t472 = m(5) * t488 - t570 * mrSges(5,2) + t508 * mrSges(5,3) + t540 * t518 - t571 * t532 + t590;
t589 = -t574 * t491 + t573 * t493;
t487 = -0.2e1 * qJD(4) * t541 - t589;
t531 = -t571 * mrSges(5,2) + t540 * mrSges(5,3);
t485 = -t570 * pkin(4) - t569 * pkin(8) + (t601 + t519) * t541 + t589;
t586 = -m(6) * t485 + t495 * mrSges(6,1) - t496 * mrSges(6,2) + t526 * t511 - t527 * t512;
t477 = m(5) * t487 + t570 * mrSges(5,1) - t509 * mrSges(5,3) - t541 * t518 + t571 * t531 + t586;
t467 = t573 * t472 + t574 * t477;
t542 = -t553 * mrSges(4,1) + t554 * mrSges(4,2);
t545 = -t571 * mrSges(4,2) + t553 * mrSges(4,3);
t465 = m(4) * t503 + t570 * mrSges(4,1) - t529 * mrSges(4,3) - t554 * t542 + t571 * t545 + t467;
t547 = t571 * mrSges(4,1) - t554 * mrSges(4,3);
t591 = t574 * t472 - t573 * t477;
t466 = m(4) * t504 - t570 * mrSges(4,2) + t528 * mrSges(4,3) + t553 * t542 - t571 * t547 + t591;
t459 = t580 * t465 + t576 * t466;
t543 = -t581 * g(3) - t599;
t560 = (-mrSges(3,1) * t581 + mrSges(3,2) * t577) * qJD(1);
t596 = qJD(1) * t581;
t564 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t596;
t457 = m(3) * t543 + qJDD(2) * mrSges(3,1) - t561 * mrSges(3,3) + qJD(2) * t564 - t560 * t597 + t459;
t563 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t597;
t592 = -t576 * t465 + t580 * t466;
t458 = m(3) * t544 - qJDD(2) * mrSges(3,2) + t562 * mrSges(3,3) - qJD(2) * t563 + t560 * t596 + t592;
t593 = -t577 * t457 + t581 * t458;
t451 = m(2) * t567 - t583 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t593;
t555 = -t583 * pkin(6) + t588;
t473 = t579 * t481 + t575 * t482;
t587 = m(5) * t498 - t508 * mrSges(5,1) + t509 * mrSges(5,2) - t540 * t531 + t541 * t532 + t473;
t585 = m(4) * t530 - t528 * mrSges(4,1) + t529 * mrSges(4,2) - t553 * t545 + t554 * t547 + t587;
t584 = -m(3) * t555 + t562 * mrSges(3,1) - t561 * mrSges(3,2) - t563 * t597 + t564 * t596 - t585;
t469 = m(2) * t566 + qJDD(1) * mrSges(2,1) - t583 * mrSges(2,2) + t584;
t598 = t578 * t451 + t582 * t469;
t452 = t581 * t457 + t577 * t458;
t594 = t582 * t451 - t578 * t469;
t552 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t577 + Ifges(3,4) * t581) * qJD(1);
t551 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t577 + Ifges(3,2) * t581) * qJD(1);
t550 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t577 + Ifges(3,6) * t581) * qJD(1);
t537 = Ifges(4,1) * t554 + Ifges(4,4) * t553 + Ifges(4,5) * t571;
t536 = Ifges(4,4) * t554 + Ifges(4,2) * t553 + Ifges(4,6) * t571;
t535 = Ifges(4,5) * t554 + Ifges(4,6) * t553 + Ifges(4,3) * t571;
t515 = Ifges(5,1) * t541 + Ifges(5,4) * t540 + Ifges(5,5) * t571;
t514 = Ifges(5,4) * t541 + Ifges(5,2) * t540 + Ifges(5,6) * t571;
t513 = Ifges(5,5) * t541 + Ifges(5,6) * t540 + Ifges(5,3) * t571;
t501 = Ifges(6,1) * t527 + Ifges(6,4) * t526 + Ifges(6,5) * t534;
t500 = Ifges(6,4) * t527 + Ifges(6,2) * t526 + Ifges(6,6) * t534;
t499 = Ifges(6,5) * t527 + Ifges(6,6) * t526 + Ifges(6,3) * t534;
t475 = mrSges(6,2) * t485 - mrSges(6,3) * t483 + Ifges(6,1) * t496 + Ifges(6,4) * t495 + Ifges(6,5) * t507 + t526 * t499 - t534 * t500;
t474 = -mrSges(6,1) * t485 + mrSges(6,3) * t484 + Ifges(6,4) * t496 + Ifges(6,2) * t495 + Ifges(6,6) * t507 - t527 * t499 + t534 * t501;
t461 = -mrSges(5,1) * t498 - mrSges(6,1) * t483 + mrSges(6,2) * t484 + mrSges(5,3) * t488 + Ifges(5,4) * t509 - Ifges(6,5) * t496 + Ifges(5,2) * t508 + Ifges(5,6) * t570 - Ifges(6,6) * t495 - Ifges(6,3) * t507 - pkin(4) * t473 - t527 * t500 + t526 * t501 - t541 * t513 + t571 * t515;
t460 = mrSges(5,2) * t498 - mrSges(5,3) * t487 + Ifges(5,1) * t509 + Ifges(5,4) * t508 + Ifges(5,5) * t570 - pkin(8) * t473 - t575 * t474 + t579 * t475 + t540 * t513 - t571 * t514;
t453 = mrSges(4,2) * t530 - mrSges(4,3) * t503 + Ifges(4,1) * t529 + Ifges(4,4) * t528 + Ifges(4,5) * t570 - qJ(4) * t467 + t574 * t460 - t573 * t461 + t553 * t535 - t571 * t536;
t448 = -mrSges(4,1) * t530 + mrSges(4,3) * t504 + Ifges(4,4) * t529 + Ifges(4,2) * t528 + Ifges(4,6) * t570 - pkin(3) * t587 + qJ(4) * t591 + t573 * t460 + t574 * t461 - t554 * t535 + t571 * t537;
t447 = -pkin(8) * t590 - pkin(4) * t586 + Ifges(2,6) * qJDD(1) - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) - pkin(1) * t452 + (-Ifges(4,3) - Ifges(5,3)) * t570 + (-t577 * t551 + t581 * t552) * qJD(1) - t579 * t474 + t583 * Ifges(2,5) - t575 * t475 - Ifges(3,5) * t561 - Ifges(3,6) * t562 + mrSges(2,3) * t567 + t553 * t537 - t554 * t536 + t540 * t515 - t541 * t514 - mrSges(3,1) * t543 + mrSges(3,2) * t544 - Ifges(4,6) * t528 - Ifges(4,5) * t529 - Ifges(5,6) * t508 - Ifges(5,5) * t509 - mrSges(4,1) * t503 + mrSges(4,2) * t504 - mrSges(5,1) * t487 + mrSges(5,2) * t488 - pkin(3) * t467 - pkin(2) * t459;
t446 = mrSges(3,2) * t555 - mrSges(3,3) * t543 + Ifges(3,1) * t561 + Ifges(3,4) * t562 + Ifges(3,5) * qJDD(2) - pkin(7) * t459 - qJD(2) * t551 - t576 * t448 + t580 * t453 + t550 * t596;
t445 = -mrSges(3,1) * t555 + mrSges(3,3) * t544 + Ifges(3,4) * t561 + Ifges(3,2) * t562 + Ifges(3,6) * qJDD(2) - pkin(2) * t585 + pkin(7) * t592 + qJD(2) * t552 + t580 * t448 + t576 * t453 - t550 * t597;
t444 = -mrSges(2,2) * g(3) - mrSges(2,3) * t566 + Ifges(2,5) * qJDD(1) - t583 * Ifges(2,6) - pkin(6) * t452 - t577 * t445 + t581 * t446;
t1 = [-m(1) * g(1) + t594; -m(1) * g(2) + t598; (-m(1) - m(2)) * g(3) + t452; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t598 + t582 * t444 - t578 * t447; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t594 + t578 * t444 + t582 * t447; -mrSges(1,1) * g(2) + mrSges(2,1) * t566 + mrSges(1,2) * g(1) - mrSges(2,2) * t567 + Ifges(2,3) * qJDD(1) + pkin(1) * t584 + pkin(6) * t593 + t581 * t445 + t577 * t446;];
tauB = t1;
