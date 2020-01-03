% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:04:09
% EndTime: 2019-12-31 21:04:12
% DurationCPUTime: 2.15s
% Computational Cost: add. (15211->267), mult. (29634->312), div. (0->0), fcn. (17364->6), ass. (0->100)
t598 = Ifges(4,1) + Ifges(5,1) + Ifges(6,1);
t583 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t597 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t596 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t581 = -Ifges(5,6) + Ifges(6,6) + Ifges(4,6);
t595 = Ifges(6,3) + Ifges(4,3) + Ifges(5,2);
t558 = sin(qJ(3));
t559 = sin(qJ(2));
t586 = qJD(1) * t559;
t592 = cos(qJ(3));
t539 = -t592 * qJD(2) + t558 * t586;
t561 = cos(qJ(2));
t584 = qJD(1) * qJD(2);
t576 = t561 * t584;
t543 = t559 * qJDD(1) + t576;
t508 = -t539 * qJD(3) + t558 * qJDD(2) + t592 * t543;
t560 = sin(qJ(1));
t562 = cos(qJ(1));
t549 = -t562 * g(1) - t560 * g(2);
t564 = qJD(1) ^ 2;
t531 = -t564 * pkin(1) + qJDD(1) * pkin(6) + t549;
t522 = -t561 * g(3) - t559 * t531;
t542 = (-pkin(2) * t561 - pkin(7) * t559) * qJD(1);
t563 = qJD(2) ^ 2;
t568 = qJDD(2) * pkin(2) + t563 * pkin(7) - t542 * t586 + t522;
t585 = t561 * qJD(1);
t551 = -qJD(3) + t585;
t590 = t539 * t551;
t594 = -(t508 + t590) * qJ(4) - t568;
t593 = -2 * qJD(4);
t591 = -mrSges(4,3) - mrSges(5,2);
t515 = -t551 * mrSges(6,2) + t539 * mrSges(6,3);
t589 = t551 * t515;
t523 = -t559 * g(3) + t561 * t531;
t541 = (-mrSges(3,1) * t561 + mrSges(3,2) * t559) * qJD(1);
t553 = t559 * t584;
t544 = t561 * qJDD(1) - t553;
t545 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t586;
t548 = t560 * g(1) - t562 * g(2);
t530 = -qJDD(1) * pkin(1) - t564 * pkin(6) - t548;
t479 = (-t543 - t576) * pkin(7) + (-t544 + t553) * pkin(2) + t530;
t483 = -t563 * pkin(2) + qJDD(2) * pkin(7) + t542 * t585 + t523;
t476 = t592 * t479 - t558 * t483;
t516 = t551 * mrSges(4,2) - t539 * mrSges(4,3);
t538 = -qJDD(3) + t544;
t540 = t558 * qJD(2) + t592 * t586;
t511 = t539 * pkin(3) - t540 * qJ(4);
t550 = t551 ^ 2;
t475 = t538 * pkin(3) - t550 * qJ(4) + t540 * t511 + qJDD(4) - t476;
t521 = -t539 * mrSges(5,2) - t551 * mrSges(5,3);
t467 = -0.2e1 * qJD(5) * t540 + (-t508 + t590) * qJ(5) + (t539 * t540 + t538) * pkin(4) + t475;
t513 = -t539 * mrSges(6,1) + t540 * mrSges(6,2);
t571 = -m(6) * t467 + t508 * mrSges(6,3) + t540 * t513;
t567 = -m(5) * t475 - t538 * mrSges(5,1) - t551 * t521 + t571;
t512 = t539 * mrSges(5,1) - t540 * mrSges(5,3);
t587 = -t539 * mrSges(4,1) - t540 * mrSges(4,2) - t512;
t463 = m(4) * t476 + (-t515 - t516) * t551 + t587 * t540 + (-mrSges(4,1) - mrSges(6,1)) * t538 + t591 * t508 + t567;
t477 = t558 * t479 + t592 * t483;
t507 = t540 * qJD(3) - t592 * qJDD(2) + t558 * t543;
t518 = t551 * mrSges(6,1) - t540 * mrSges(6,3);
t519 = -t551 * mrSges(4,1) - t540 * mrSges(4,3);
t473 = -t550 * pkin(3) - t538 * qJ(4) - t539 * t511 + t551 * t593 + t477;
t520 = t551 * mrSges(5,1) + t540 * mrSges(5,2);
t517 = t551 * pkin(4) - t540 * qJ(5);
t537 = t539 ^ 2;
t469 = -t537 * pkin(4) + t507 * qJ(5) + 0.2e1 * qJD(5) * t539 - t551 * t517 + t473;
t580 = m(6) * t469 + t507 * mrSges(6,3) + t539 * t513;
t569 = m(5) * t473 - t538 * mrSges(5,3) - t551 * t520 + t580;
t464 = m(4) * t477 + (-t518 + t519) * t551 + t587 * t539 + (mrSges(4,2) - mrSges(6,2)) * t538 + t591 * t507 + t569;
t573 = -t558 * t463 + t592 * t464;
t458 = m(3) * t523 - qJDD(2) * mrSges(3,2) + t544 * mrSges(3,3) - qJD(2) * t545 + t541 * t585 + t573;
t546 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t585;
t474 = t540 * t593 + (-t540 * t551 + t507) * pkin(3) + t594;
t471 = -t537 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t507 + (pkin(3) * t551 + (2 * qJD(4)) + t517) * t540 - t594;
t570 = -m(6) * t471 + t507 * mrSges(6,1) - t508 * mrSges(6,2) + t539 * t515 - t540 * t518;
t465 = m(5) * t474 + t507 * mrSges(5,1) - t508 * mrSges(5,3) - t540 * t520 + t539 * t521 + t570;
t565 = m(4) * t568 - t507 * mrSges(4,1) - t508 * mrSges(4,2) - t539 * t516 - t540 * t519 - t465;
t461 = m(3) * t522 + qJDD(2) * mrSges(3,1) - t543 * mrSges(3,3) + qJD(2) * t546 - t541 * t586 + t565;
t574 = t561 * t458 - t559 * t461;
t451 = m(2) * t549 - t564 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t574;
t459 = t592 * t463 + t558 * t464;
t566 = -m(3) * t530 + t544 * mrSges(3,1) - t543 * mrSges(3,2) - t545 * t586 + t546 * t585 - t459;
t455 = m(2) * t548 + qJDD(1) * mrSges(2,1) - t564 * mrSges(2,2) + t566;
t588 = t560 * t451 + t562 * t455;
t452 = t559 * t458 + t561 * t461;
t579 = t581 * t539 - t597 * t540 + t595 * t551;
t578 = t596 * t539 + t583 * t540 - t581 * t551;
t577 = t583 * t539 - t598 * t540 + t597 * t551;
t575 = t562 * t451 - t560 * t455;
t529 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t559 + Ifges(3,4) * t561) * qJD(1);
t528 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t559 + Ifges(3,2) * t561) * qJD(1);
t527 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t559 + Ifges(3,6) * t561) * qJD(1);
t466 = t538 * mrSges(6,1) - t571 + t589;
t453 = -mrSges(4,2) * t568 + mrSges(5,2) * t475 + mrSges(6,2) * t471 - mrSges(4,3) * t476 - mrSges(5,3) * t474 - mrSges(6,3) * t467 - qJ(4) * t465 - qJ(5) * t466 - t583 * t507 + t598 * t508 - t538 * t597 + t579 * t539 + t578 * t551;
t448 = mrSges(4,1) * t568 + mrSges(4,3) * t477 - mrSges(5,1) * t474 + mrSges(5,2) * t473 + mrSges(6,1) * t471 - mrSges(6,3) * t469 - pkin(4) * t570 - qJ(5) * t580 - pkin(3) * t465 + (qJ(5) * t518 + t577) * t551 + t579 * t540 + (qJ(5) * mrSges(6,2) - t581) * t538 + t583 * t508 + t596 * t507;
t447 = -qJ(4) * (-t551 * t518 + t569) - pkin(3) * (t567 - t589) + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t543 + Ifges(3,2) * t544 + mrSges(3,3) * t523 + qJD(2) * t529 - mrSges(3,1) * t530 - mrSges(5,3) * t473 + mrSges(5,1) * t475 - mrSges(4,1) * t476 + mrSges(4,2) * t477 + pkin(4) * t466 + mrSges(6,1) * t467 - mrSges(6,2) * t469 - pkin(2) * t459 - t527 * t586 + (pkin(3) * t512 - t578) * t540 + (qJ(4) * t512 + t577) * t539 + (pkin(3) * mrSges(6,1) + qJ(4) * mrSges(6,2) + t595) * t538 + (pkin(3) * mrSges(5,2) - t597) * t508 + (qJ(4) * mrSges(5,2) + t581) * t507;
t446 = mrSges(3,2) * t530 - mrSges(3,3) * t522 + Ifges(3,1) * t543 + Ifges(3,4) * t544 + Ifges(3,5) * qJDD(2) - pkin(7) * t459 - qJD(2) * t528 - t558 * t448 + t592 * t453 + t527 * t585;
t445 = Ifges(2,6) * qJDD(1) + t564 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t549 - Ifges(3,5) * t543 - Ifges(3,6) * t544 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t522 + mrSges(3,2) * t523 - t558 * t453 - t592 * t448 - pkin(2) * t565 - pkin(7) * t573 - pkin(1) * t452 + (-t559 * t528 + t561 * t529) * qJD(1);
t444 = -mrSges(2,2) * g(3) - mrSges(2,3) * t548 + Ifges(2,5) * qJDD(1) - t564 * Ifges(2,6) - pkin(6) * t452 + t561 * t446 - t559 * t447;
t1 = [-m(1) * g(1) + t575; -m(1) * g(2) + t588; (-m(1) - m(2)) * g(3) + t452; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t588 + t562 * t444 - t560 * t445; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t575 + t560 * t444 + t562 * t445; -mrSges(1,1) * g(2) + mrSges(2,1) * t548 + mrSges(1,2) * g(1) - mrSges(2,2) * t549 + Ifges(2,3) * qJDD(1) + pkin(1) * t566 + pkin(6) * t574 + t559 * t446 + t561 * t447;];
tauB = t1;
