% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:41:55
% EndTime: 2019-05-06 08:41:59
% DurationCPUTime: 2.51s
% Computational Cost: add. (13250->316), mult. (29305->367), div. (0->0), fcn. (16664->8), ass. (0->130)
t591 = -2 * qJD(4);
t590 = Ifges(5,1) + Ifges(6,1);
t579 = Ifges(3,4) + Ifges(4,6);
t578 = Ifges(5,4) - Ifges(6,5);
t577 = Ifges(3,5) - Ifges(4,4);
t576 = Ifges(5,5) + Ifges(6,4);
t589 = Ifges(3,2) + Ifges(4,3);
t588 = Ifges(4,2) + Ifges(3,1);
t587 = Ifges(5,2) + Ifges(6,3);
t586 = -Ifges(6,2) - Ifges(5,3);
t575 = Ifges(3,6) - Ifges(4,5);
t574 = Ifges(5,6) - Ifges(6,6);
t585 = Ifges(3,3) + Ifges(4,1);
t536 = cos(qJ(2));
t533 = sin(qJ(2));
t561 = qJD(1) * qJD(2);
t557 = t533 * t561;
t506 = qJDD(1) * t536 - t557;
t529 = sin(pkin(9));
t530 = cos(pkin(9));
t476 = qJDD(2) * t530 - t506 * t529;
t563 = qJD(1) * t536;
t496 = qJD(2) * t529 + t530 * t563;
t564 = qJD(1) * t533;
t559 = t496 * t564;
t584 = (-t476 + t559) * qJ(5);
t511 = pkin(3) * t564 - qJD(2) * qJ(4);
t528 = t536 ^ 2;
t539 = qJD(1) ^ 2;
t558 = t536 * t561;
t505 = qJDD(1) * t533 + t558;
t534 = sin(qJ(1));
t537 = cos(qJ(1));
t556 = g(1) * t534 - t537 * g(2);
t548 = -qJDD(1) * pkin(1) - t556;
t542 = pkin(2) * t557 - 0.2e1 * qJD(3) * t564 + (-t505 - t558) * qJ(3) + t548;
t432 = -t511 * t564 + (-pkin(3) * t528 - pkin(7)) * t539 + (-pkin(2) - qJ(4)) * t506 + t542;
t552 = -g(1) * t537 - g(2) * t534;
t493 = -pkin(1) * t539 + qJDD(1) * pkin(7) + t552;
t467 = -t536 * g(3) - t533 * t493;
t502 = (-t536 * pkin(2) - t533 * qJ(3)) * qJD(1);
t538 = qJD(2) ^ 2;
t447 = -qJDD(2) * pkin(2) - qJ(3) * t538 + t502 * t564 + qJDD(3) - t467;
t442 = (-t533 * t536 * t539 - qJDD(2)) * qJ(4) + (t505 - t558) * pkin(3) + t447;
t497 = qJD(2) * t530 - t529 * t563;
t426 = -t529 * t432 + t442 * t530 + t497 * t591;
t583 = 2 * qJD(5);
t582 = pkin(7) * t539;
t581 = mrSges(3,1) - mrSges(4,2);
t580 = -mrSges(5,3) - mrSges(6,2);
t573 = t533 ^ 2 * t539;
t572 = t587 * t496 - t578 * t497 - t574 * t564;
t571 = t574 * t496 - t576 * t497 + t586 * t564;
t570 = -t578 * t496 + t590 * t497 + t576 * t564;
t463 = mrSges(6,1) * t496 - mrSges(6,3) * t497;
t569 = -mrSges(5,1) * t496 - mrSges(5,2) * t497 - t463;
t468 = -t533 * g(3) + t536 * t493;
t568 = t585 * qJD(2) + (t533 * t577 + t536 * t575) * qJD(1);
t567 = t575 * qJD(2) + (t533 * t579 + t536 * t589) * qJD(1);
t566 = t577 * qJD(2) + (t588 * t533 + t536 * t579) * qJD(1);
t512 = -mrSges(4,1) * t563 - qJD(2) * mrSges(4,3);
t565 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t563 - t512;
t560 = qJD(3) * qJD(2);
t427 = t530 * t432 + t529 * t442 + t496 * t591;
t473 = mrSges(5,1) * t564 - mrSges(5,3) * t497;
t475 = qJDD(2) * t529 + t530 * t506;
t462 = pkin(4) * t496 - qJ(5) * t497;
t423 = -pkin(4) * t573 + t505 * qJ(5) - t496 * t462 + t564 * t583 + t427;
t474 = -mrSges(6,1) * t564 + mrSges(6,2) * t497;
t424 = -pkin(4) * t505 - qJ(5) * t573 + t497 * t462 + qJDD(5) - t426;
t420 = (-t476 - t559) * pkin(8) + (t496 * t497 - t505) * pkin(5) + t424;
t477 = -pkin(5) * t564 - pkin(8) * t497;
t494 = t496 ^ 2;
t421 = -pkin(5) * t494 + pkin(8) * t475 + t477 * t564 + t423;
t532 = sin(qJ(6));
t535 = cos(qJ(6));
t418 = t420 * t535 - t421 * t532;
t460 = t496 * t535 - t497 * t532;
t435 = qJD(6) * t460 + t475 * t532 + t476 * t535;
t461 = t496 * t532 + t497 * t535;
t444 = -mrSges(7,1) * t460 + mrSges(7,2) * t461;
t516 = qJD(6) - t564;
t448 = -mrSges(7,2) * t516 + mrSges(7,3) * t460;
t501 = qJDD(6) - t505;
t415 = m(7) * t418 + mrSges(7,1) * t501 - mrSges(7,3) * t435 - t444 * t461 + t448 * t516;
t419 = t420 * t532 + t421 * t535;
t434 = -qJD(6) * t461 + t475 * t535 - t476 * t532;
t449 = mrSges(7,1) * t516 - mrSges(7,3) * t461;
t416 = m(7) * t419 - mrSges(7,2) * t501 + mrSges(7,3) * t434 + t444 * t460 - t449 * t516;
t554 = -t415 * t532 + t535 * t416;
t547 = m(6) * t423 + t505 * mrSges(6,3) + t474 * t564 + t554;
t407 = m(5) * t427 - mrSges(5,2) * t505 - t473 * t564 + t580 * t475 + t569 * t496 + t547;
t472 = -mrSges(5,2) * t564 - mrSges(5,3) * t496;
t410 = t415 * t535 + t416 * t532;
t471 = -mrSges(6,2) * t496 + mrSges(6,3) * t564;
t544 = -m(6) * t424 + t505 * mrSges(6,1) + t471 * t564 - t410;
t408 = m(5) * t426 + mrSges(5,1) * t505 + t472 * t564 + t580 * t476 + t569 * t497 + t544;
t555 = t530 * t407 - t408 * t529;
t550 = t538 * pkin(2) - qJDD(2) * qJ(3) - t502 * t563 - t468;
t405 = t529 * t407 + t530 * t408;
t445 = -pkin(2) * t506 + t542 - t582;
t549 = m(4) * t445 + t555;
t546 = m(4) * t447 + t505 * mrSges(4,1) + t405;
t521 = -0.2e1 * t560;
t543 = qJ(4) * t528 * t539 - pkin(3) * t506 - qJD(2) * t511 - qJDD(4) + t550;
t425 = -pkin(8) * t494 + t521 + (-pkin(4) - pkin(5)) * t475 - t584 + (-pkin(4) * t564 + t477 + t583) * t497 + t543;
t545 = -m(7) * t425 + t434 * mrSges(7,1) - t435 * mrSges(7,2) + t460 * t448 - t461 * t449;
t438 = 0.2e1 * t560 - t543;
t440 = Ifges(7,4) * t461 + Ifges(7,2) * t460 + Ifges(7,6) * t516;
t441 = Ifges(7,1) * t461 + Ifges(7,4) * t460 + Ifges(7,5) * t516;
t541 = mrSges(7,1) * t418 - mrSges(7,2) * t419 + Ifges(7,5) * t435 + Ifges(7,6) * t434 + Ifges(7,3) * t501 + t461 * t440 - t460 * t441;
t429 = -0.2e1 * qJD(5) * t497 + t584 + (t497 * t564 + t475) * pkin(4) + t438;
t417 = m(6) * t429 + t475 * mrSges(6,1) - t476 * mrSges(6,3) + t496 * t471 - t497 * t474 + t545;
t413 = m(5) * t438 + t475 * mrSges(5,1) + t476 * mrSges(5,2) + t496 * t472 + t497 * t473 + t417;
t446 = t521 + t550;
t503 = (t536 * mrSges(4,2) - t533 * mrSges(4,3)) * qJD(1);
t513 = mrSges(4,1) * t564 + qJD(2) * mrSges(4,2);
t540 = -m(4) * t446 + qJDD(2) * mrSges(4,3) + qJD(2) * t513 + t503 * t563 + t413;
t509 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t564;
t504 = (-t536 * mrSges(3,1) + t533 * mrSges(3,2)) * qJD(1);
t492 = t548 - t582;
t439 = Ifges(7,5) * t461 + Ifges(7,6) * t460 + Ifges(7,3) * t516;
t412 = mrSges(7,2) * t425 - mrSges(7,3) * t418 + Ifges(7,1) * t435 + Ifges(7,4) * t434 + Ifges(7,5) * t501 + t439 * t460 - t440 * t516;
t411 = -mrSges(7,1) * t425 + mrSges(7,3) * t419 + Ifges(7,4) * t435 + Ifges(7,2) * t434 + Ifges(7,6) * t501 - t439 * t461 + t441 * t516;
t409 = mrSges(6,2) * t476 + t463 * t497 - t544;
t404 = qJDD(2) * mrSges(4,2) + qJD(2) * t512 + t503 * t564 + t546;
t403 = mrSges(4,2) * t506 - mrSges(4,3) * t505 + (t512 * t536 - t513 * t533) * qJD(1) + t549;
t402 = mrSges(5,2) * t438 + mrSges(6,2) * t424 - mrSges(5,3) * t426 - mrSges(6,3) * t429 - pkin(8) * t410 - qJ(5) * t417 - t532 * t411 + t535 * t412 - t578 * t475 + t590 * t476 + t571 * t496 + t576 * t505 + t572 * t564;
t401 = -mrSges(5,1) * t438 - mrSges(6,1) * t429 + mrSges(6,2) * t423 + mrSges(5,3) * t427 - pkin(4) * t417 - pkin(5) * t545 - pkin(8) * t554 - t535 * t411 - t532 * t412 - t587 * t475 + t578 * t476 + t571 * t497 + t574 * t505 + t570 * t564;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t556 - mrSges(2,2) * t552 + t533 * (-qJ(3) * t403 - t567 * qJD(2) - t541 - pkin(4) * t409 - pkin(5) * t410 + (-t586 + t588) * t505 - t572 * t497 + (-qJ(5) * t463 + t570) * t496 + t576 * t476 + t579 * t506 + qJ(5) * t547 + t568 * t563 + (-qJ(5) * mrSges(6,2) - t574) * t475 + t577 * qJDD(2) + pkin(3) * t405 + mrSges(6,3) * t423 - mrSges(6,1) * t424 + mrSges(5,1) * t426 - mrSges(5,2) * t427 - mrSges(4,3) * t445 + mrSges(4,1) * t447 - mrSges(3,3) * t467 + mrSges(3,2) * t492) + t536 * (-mrSges(3,1) * t492 - mrSges(4,1) * t446 + mrSges(4,2) * t445 + mrSges(3,3) * t468 - pkin(2) * t403 + pkin(3) * t413 - qJ(4) * t555 + t566 * qJD(2) + t575 * qJDD(2) - t530 * t401 - t529 * t402 + t579 * t505 + t589 * t506 - t568 * t564) + pkin(1) * (-m(3) * t492 + t581 * t506 + (-mrSges(3,2) + mrSges(4,3)) * t505 + (t565 * t536 + (-t509 + t513) * t533) * qJD(1) - t549) + pkin(7) * (t536 * (t540 + t504 * t563 - qJDD(2) * mrSges(3,2) + (mrSges(3,3) + mrSges(4,1)) * t506 + m(3) * t468 - qJD(2) * t509) + (-m(3) * t467 + t505 * mrSges(3,3) - t581 * qJDD(2) - t565 * qJD(2) + (t503 + t504) * t564 + t546) * t533); mrSges(3,1) * t467 - mrSges(3,2) * t468 + mrSges(4,2) * t447 - mrSges(4,3) * t446 + t530 * t402 - t529 * t401 - qJ(4) * t405 - pkin(2) * t404 + qJ(3) * t540 + (qJ(3) * mrSges(4,1) + t575) * t506 + t577 * t505 + t585 * qJDD(2) + (t567 * t533 - t566 * t536) * qJD(1); t404; t413; t409; t541;];
tauJ  = t1;
