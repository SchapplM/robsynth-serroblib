% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP8
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:43
% EndTime: 2019-12-31 21:07:46
% DurationCPUTime: 2.10s
% Computational Cost: add. (15209->268), mult. (29540->309), div. (0->0), fcn. (17274->6), ass. (0->103)
t589 = Ifges(4,1) + Ifges(5,2) + Ifges(6,3);
t575 = Ifges(4,4) + Ifges(5,6) - Ifges(6,6);
t574 = Ifges(4,5) - Ifges(5,4) + Ifges(6,5);
t588 = -Ifges(4,2) - Ifges(5,3) - Ifges(6,2);
t587 = -Ifges(4,6) + Ifges(5,5) + Ifges(6,4);
t586 = Ifges(4,3) + Ifges(5,1) + Ifges(6,1);
t585 = -2 * qJD(4);
t584 = 2 * qJD(5);
t583 = cos(qJ(3));
t582 = -mrSges(6,1) - mrSges(4,3);
t548 = sin(qJ(3));
t549 = sin(qJ(2));
t578 = qJD(1) * t549;
t531 = -t583 * qJD(2) + t548 * t578;
t551 = cos(qJ(2));
t577 = t551 * qJD(1);
t542 = -qJD(3) + t577;
t581 = t531 * t542;
t550 = sin(qJ(1));
t552 = cos(qJ(1));
t540 = -t552 * g(1) - t550 * g(2);
t554 = qJD(1) ^ 2;
t523 = -t554 * pkin(1) + qJDD(1) * pkin(6) + t540;
t513 = -t549 * g(3) + t551 * t523;
t533 = (-mrSges(3,1) * t551 + mrSges(3,2) * t549) * qJD(1);
t576 = qJD(1) * qJD(2);
t568 = t549 * t576;
t536 = t551 * qJDD(1) - t568;
t537 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t578;
t539 = t550 * g(1) - t552 * g(2);
t522 = -qJDD(1) * pkin(1) - t554 * pkin(6) - t539;
t567 = t551 * t576;
t535 = t549 * qJDD(1) + t567;
t471 = (-t535 - t567) * pkin(7) + (-t536 + t568) * pkin(2) + t522;
t534 = (-pkin(2) * t551 - pkin(7) * t549) * qJD(1);
t553 = qJD(2) ^ 2;
t475 = -t553 * pkin(2) + qJDD(2) * pkin(7) + t534 * t577 + t513;
t468 = t583 * t471 - t548 * t475;
t497 = -t531 * qJD(3) + t548 * qJDD(2) + t583 * t535;
t532 = t548 * qJD(2) + t583 * t578;
t501 = -t532 * mrSges(6,2) + t531 * mrSges(6,3);
t503 = t531 * mrSges(4,1) + t532 * mrSges(4,2);
t505 = t542 * mrSges(4,2) - t531 * mrSges(4,3);
t509 = t531 * mrSges(5,1) + t542 * mrSges(5,3);
t530 = qJDD(3) - t536;
t502 = t531 * pkin(3) - t532 * qJ(4);
t541 = t542 ^ 2;
t467 = -t530 * pkin(3) - t541 * qJ(4) + t532 * t502 + qJDD(4) - t468;
t504 = -t531 * mrSges(5,2) - t532 * mrSges(5,3);
t460 = t542 * t584 + (t531 * t532 - t530) * qJ(5) + (t497 - t581) * pkin(4) + t467;
t510 = -t531 * mrSges(6,1) - t542 * mrSges(6,2);
t563 = m(6) * t460 - t530 * mrSges(6,3) + t542 * t510;
t560 = -m(5) * t467 - t497 * mrSges(5,1) - t532 * t504 - t563;
t455 = m(4) * t468 + (-t505 + t509) * t542 + (-t501 - t503) * t532 + (mrSges(4,1) - mrSges(5,2)) * t530 + t582 * t497 + t560;
t469 = t548 * t471 + t583 * t475;
t496 = t532 * qJD(3) - t583 * qJDD(2) + t548 * t535;
t506 = -t542 * mrSges(4,1) - t532 * mrSges(4,3);
t558 = -t541 * pkin(3) + t530 * qJ(4) - t531 * t502 + t469;
t465 = 0.2e1 * qJD(4) * t542 - t558;
t511 = t532 * mrSges(5,1) - t542 * mrSges(5,2);
t507 = t532 * pkin(4) + t542 * qJ(5);
t529 = t531 ^ 2;
t464 = -t496 * pkin(4) - t529 * qJ(5) + qJDD(5) + (t585 - t507) * t542 + t558;
t508 = t532 * mrSges(6,1) + t542 * mrSges(6,3);
t569 = -m(6) * t464 - t530 * mrSges(6,2) + t542 * t508;
t561 = -m(5) * t465 + t530 * mrSges(5,3) - t542 * t511 - t569;
t579 = -t501 - t504;
t457 = m(4) * t469 - t530 * mrSges(4,2) + t542 * t506 + (-t503 + t579) * t531 + (-mrSges(5,1) + t582) * t496 + t561;
t564 = -t548 * t455 + t583 * t457;
t451 = m(3) * t513 - qJDD(2) * mrSges(3,2) + t536 * mrSges(3,3) - qJD(2) * t537 + t533 * t577 + t564;
t512 = -t551 * g(3) - t549 * t523;
t538 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t577;
t474 = -qJDD(2) * pkin(2) - t553 * pkin(7) + t534 * t578 - t512;
t557 = (-t497 - t581) * qJ(4) + t474 + (-t542 * pkin(3) + t585) * t532;
t466 = t496 * pkin(3) + t557;
t462 = -t529 * pkin(4) + t531 * t584 - t532 * t507 + (pkin(3) + qJ(5)) * t496 + t557;
t562 = m(6) * t462 - t497 * mrSges(6,2) + t496 * mrSges(6,3) - t532 * t508 + t531 * t510;
t559 = -m(5) * t466 + t496 * mrSges(5,2) + t531 * t509 - t562;
t555 = -m(4) * t474 - t496 * mrSges(4,1) - t531 * t505 + (-t506 + t511) * t532 + (-mrSges(4,2) + mrSges(5,3)) * t497 + t559;
t454 = m(3) * t512 + qJDD(2) * mrSges(3,1) - t535 * mrSges(3,3) + qJD(2) * t538 - t533 * t578 + t555;
t565 = t551 * t451 - t549 * t454;
t444 = m(2) * t540 - t554 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t565;
t452 = t583 * t455 + t548 * t457;
t556 = -m(3) * t522 + t536 * mrSges(3,1) - t535 * mrSges(3,2) - t537 * t578 + t538 * t577 - t452;
t448 = m(2) * t539 + qJDD(1) * mrSges(2,1) - t554 * mrSges(2,2) + t556;
t580 = t550 * t444 + t552 * t448;
t445 = t549 * t451 + t551 * t454;
t572 = -t587 * t531 - t574 * t532 + t586 * t542;
t571 = t588 * t531 + t575 * t532 + t587 * t542;
t570 = t575 * t531 - t589 * t532 + t574 * t542;
t566 = t552 * t444 - t550 * t448;
t521 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t549 + Ifges(3,4) * t551) * qJD(1);
t520 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t549 + Ifges(3,2) * t551) * qJD(1);
t519 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t549 + Ifges(3,6) * t551) * qJD(1);
t459 = t497 * mrSges(6,1) + t532 * t501 + t563;
t458 = -t497 * mrSges(5,3) - t532 * t511 - t559;
t446 = mrSges(5,1) * t467 + mrSges(6,1) * t460 + mrSges(4,2) * t474 - mrSges(6,2) * t462 - mrSges(4,3) * t468 - mrSges(5,3) * t466 + pkin(4) * t459 - qJ(4) * t458 - t575 * t496 + t589 * t497 + t574 * t530 + t572 * t531 + t571 * t542;
t441 = -mrSges(4,1) * t474 + mrSges(4,3) * t469 - mrSges(5,1) * t465 + mrSges(5,2) * t466 + mrSges(6,1) * t464 - mrSges(6,3) * t462 - pkin(4) * (t531 * t501 + t569) - qJ(5) * t562 - pkin(3) * t458 + t570 * t542 + t572 * t532 - t587 * t530 + t575 * t497 + (-pkin(4) * mrSges(6,1) + t588) * t496;
t440 = -qJ(4) * t561 + Ifges(3,6) * qJDD(2) - pkin(3) * (t542 * t509 + t560) + Ifges(3,4) * t535 + Ifges(3,2) * t536 + mrSges(3,3) * t513 + qJD(2) * t521 - mrSges(3,1) * t522 - t519 * t578 - mrSges(4,1) * t468 + mrSges(4,2) * t469 - mrSges(6,2) * t464 + mrSges(5,3) * t465 - mrSges(5,2) * t467 + qJ(5) * t459 + mrSges(6,3) * t460 - pkin(2) * t452 + (pkin(3) * t501 - t571) * t532 + (pkin(3) * mrSges(5,2) - t586) * t530 + (pkin(3) * mrSges(6,1) - t574) * t497 + (-qJ(4) * t579 + t570) * t531 + (-qJ(4) * (-mrSges(5,1) - mrSges(6,1)) - t587) * t496;
t439 = mrSges(3,2) * t522 - mrSges(3,3) * t512 + Ifges(3,1) * t535 + Ifges(3,4) * t536 + Ifges(3,5) * qJDD(2) - pkin(7) * t452 - qJD(2) * t520 - t548 * t441 + t583 * t446 + t519 * t577;
t438 = Ifges(2,6) * qJDD(1) + t554 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t540 - Ifges(3,5) * t535 - Ifges(3,6) * t536 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t512 + mrSges(3,2) * t513 - t548 * t446 - t583 * t441 - pkin(2) * t555 - pkin(7) * t564 - pkin(1) * t445 + (-t549 * t520 + t551 * t521) * qJD(1);
t437 = -mrSges(2,2) * g(3) - mrSges(2,3) * t539 + Ifges(2,5) * qJDD(1) - t554 * Ifges(2,6) - pkin(6) * t445 + t551 * t439 - t549 * t440;
t1 = [-m(1) * g(1) + t566; -m(1) * g(2) + t580; (-m(1) - m(2)) * g(3) + t445; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t580 + t552 * t437 - t550 * t438; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t566 + t550 * t437 + t552 * t438; -mrSges(1,1) * g(2) + mrSges(2,1) * t539 + mrSges(1,2) * g(1) - mrSges(2,2) * t540 + Ifges(2,3) * qJDD(1) + pkin(1) * t556 + pkin(6) * t565 + t549 * t439 + t551 * t440;];
tauB = t1;
