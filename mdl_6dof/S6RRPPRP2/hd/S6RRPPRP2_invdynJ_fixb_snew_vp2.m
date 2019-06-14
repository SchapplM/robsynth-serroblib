% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-05-06 09:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:09:20
% EndTime: 2019-05-06 09:09:25
% DurationCPUTime: 2.94s
% Computational Cost: add. (16988->314), mult. (39199->362), div. (0->0), fcn. (25841->8), ass. (0->128)
t586 = Ifges(6,4) + Ifges(7,4);
t602 = Ifges(6,2) + Ifges(7,2);
t595 = Ifges(6,6) + Ifges(7,6);
t601 = -2 * qJD(4);
t600 = Ifges(4,1) + Ifges(5,2);
t599 = -Ifges(5,1) - Ifges(4,3);
t598 = Ifges(6,1) + Ifges(7,1);
t585 = Ifges(4,5) - Ifges(5,4);
t597 = Ifges(6,5) + Ifges(7,5);
t596 = -Ifges(4,2) - Ifges(5,3);
t584 = Ifges(4,6) - Ifges(5,5);
t583 = -Ifges(5,6) - Ifges(4,4);
t594 = Ifges(6,3) + Ifges(7,3);
t538 = sin(pkin(9));
t543 = cos(qJ(2));
t569 = qJD(1) * t543;
t540 = sin(qJ(2));
t570 = qJD(1) * t540;
t581 = cos(pkin(9));
t517 = t538 * t570 - t581 * t569;
t539 = sin(qJ(5));
t542 = cos(qJ(5));
t499 = -qJD(2) * t539 + t517 * t542;
t500 = qJD(2) * t542 + t517 * t539;
t518 = (t538 * t543 + t581 * t540) * qJD(1);
t515 = qJD(5) + t518;
t593 = t602 * t499 + t500 * t586 + t595 * t515;
t487 = pkin(3) * t517 - qJ(4) * t518;
t545 = qJD(2) ^ 2;
t565 = qJD(1) * qJD(2);
t527 = qJDD(1) * t540 + t543 * t565;
t546 = qJD(1) ^ 2;
t541 = sin(qJ(1));
t544 = cos(qJ(1));
t557 = -g(1) * t544 - g(2) * t541;
t524 = -pkin(1) * t546 + qJDD(1) * pkin(7) + t557;
t579 = t540 * t524;
t589 = pkin(2) * t546;
t467 = qJDD(2) * pkin(2) - t527 * qJ(3) - t579 + (qJ(3) * t565 + t540 * t589 - g(3)) * t543;
t503 = -g(3) * t540 + t543 * t524;
t528 = qJDD(1) * t543 - t540 * t565;
t529 = qJD(2) * pkin(2) - qJ(3) * t570;
t537 = t543 ^ 2;
t468 = qJ(3) * t528 - qJD(2) * t529 - t537 * t589 + t503;
t576 = t538 * t467 + t581 * t468;
t592 = pkin(3) * t545 - qJDD(2) * qJ(4) + qJD(2) * t601 + t517 * t487 - t576;
t441 = -0.2e1 * qJD(3) * t518 + t581 * t467 - t538 * t468;
t496 = t538 * t527 - t581 * t528;
t461 = -qJD(5) * t500 - qJDD(2) * t539 + t496 * t542;
t475 = -mrSges(7,2) * t515 + mrSges(7,3) * t499;
t591 = -t461 * mrSges(7,1) - t499 * t475;
t476 = -mrSges(6,2) * t515 + mrSges(6,3) * t499;
t590 = -(mrSges(6,1) + mrSges(7,1)) * t461 - (t475 + t476) * t499;
t560 = t541 * g(1) - t544 * g(2);
t555 = -qJDD(1) * pkin(1) - t560;
t474 = -t528 * pkin(2) + qJDD(3) + t529 * t570 + (-qJ(3) * t537 - pkin(7)) * t546 + t555;
t497 = t581 * t527 + t538 * t528;
t506 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t518;
t568 = qJD(2) * t517;
t547 = (-t497 + t568) * qJ(4) + t474 + (qJD(2) * pkin(3) + t601) * t518;
t440 = t496 * pkin(3) + t547;
t508 = mrSges(5,1) * t518 + qJD(2) * mrSges(5,2);
t438 = -qJDD(2) * pkin(3) - t545 * qJ(4) + t518 * t487 + qJDD(4) - t441;
t432 = (t517 * t518 - qJDD(2)) * pkin(8) + (t497 + t568) * pkin(4) + t438;
t509 = pkin(4) * t518 - qJD(2) * pkin(8);
t516 = t517 ^ 2;
t436 = -t516 * pkin(4) - t518 * t509 + (pkin(3) + pkin(8)) * t496 + t547;
t426 = t542 * t432 - t539 * t436;
t462 = qJD(5) * t499 + qJDD(2) * t542 + t496 * t539;
t469 = -mrSges(7,1) * t499 + mrSges(7,2) * t500;
t470 = -mrSges(6,1) * t499 + mrSges(6,2) * t500;
t495 = qJDD(5) + t497;
t421 = -0.2e1 * qJD(6) * t500 + (t499 * t515 - t462) * qJ(6) + (t499 * t500 + t495) * pkin(5) + t426;
t563 = m(7) * t421 + t495 * mrSges(7,1) + t515 * t475;
t414 = m(6) * t426 + t495 * mrSges(6,1) + t515 * t476 + (-t469 - t470) * t500 + (-mrSges(6,3) - mrSges(7,3)) * t462 + t563;
t427 = t539 * t432 + t542 * t436;
t478 = mrSges(7,1) * t515 - mrSges(7,3) * t500;
t479 = mrSges(6,1) * t515 - mrSges(6,3) * t500;
t477 = pkin(5) * t515 - qJ(6) * t500;
t498 = t499 ^ 2;
t423 = -pkin(5) * t498 + qJ(6) * t461 + 0.2e1 * qJD(6) * t499 - t477 * t515 + t427;
t562 = m(7) * t423 + t461 * mrSges(7,3) + t499 * t469;
t416 = m(6) * t427 + t461 * mrSges(6,3) + t499 * t470 + (-t478 - t479) * t515 + (-mrSges(6,2) - mrSges(7,2)) * t495 + t562;
t558 = -t539 * t414 + t542 * t416;
t552 = -m(5) * t440 + t497 * mrSges(5,3) + t518 * t508 - t558;
t507 = mrSges(5,1) * t517 - qJD(2) * mrSges(5,3);
t571 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t517 - t507;
t588 = mrSges(4,1) - mrSges(5,2);
t406 = m(4) * t474 + t497 * mrSges(4,2) + t588 * t496 + t518 * t506 + t571 * t517 - t552;
t488 = mrSges(4,1) * t517 + mrSges(4,2) * t518;
t410 = t542 * t414 + t539 * t416;
t489 = -mrSges(5,2) * t517 - mrSges(5,3) * t518;
t551 = -m(5) * t438 - t497 * mrSges(5,1) - t518 * t489 - t410;
t405 = m(4) * t441 - t497 * mrSges(4,3) + t571 * qJD(2) + t588 * qJDD(2) - t518 * t488 + t551;
t567 = qJD(3) * t517;
t512 = -0.2e1 * t567;
t442 = t512 + t576;
t437 = 0.2e1 * t567 + t592;
t434 = -pkin(4) * t496 - pkin(8) * t516 + qJD(2) * t509 + t512 - t592;
t429 = -pkin(5) * t461 - qJ(6) * t498 + t477 * t500 + qJDD(6) + t434;
t561 = m(7) * t429 + t462 * mrSges(7,2) + t500 * t478;
t554 = -m(6) * t434 - t462 * mrSges(6,2) - t500 * t479 - t561;
t550 = -m(5) * t437 + qJDD(2) * mrSges(5,3) + qJD(2) * t508 - t554;
t413 = t550 + (-t488 - t489) * t517 - qJDD(2) * mrSges(4,2) + (-mrSges(4,3) - mrSges(5,1)) * t496 - qJD(2) * t506 + m(4) * t442 + t590;
t402 = t581 * t405 + t538 * t413;
t578 = -t595 * t499 - t597 * t500 - t594 * t515;
t577 = -t586 * t499 - t598 * t500 - t597 * t515;
t574 = t599 * qJD(2) + t584 * t517 - t585 * t518;
t573 = t584 * qJD(2) + t596 * t517 - t583 * t518;
t572 = t585 * qJD(2) + t583 * t517 + t600 * t518;
t559 = -t405 * t538 + t581 * t413;
t418 = -t462 * mrSges(7,3) - t500 * t469 + t563;
t548 = mrSges(6,1) * t426 + mrSges(7,1) * t421 - mrSges(6,2) * t427 - mrSges(7,2) * t423 + pkin(5) * t418 + t595 * t461 + t597 * t462 + t594 * t495 + t577 * t499 + t593 * t500;
t531 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t569;
t530 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t570;
t526 = (-mrSges(3,1) * t543 + mrSges(3,2) * t540) * qJD(1);
t523 = -t546 * pkin(7) + t555;
t521 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t540 + Ifges(3,4) * t543) * qJD(1);
t520 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t540 + Ifges(3,2) * t543) * qJD(1);
t502 = -t543 * g(3) - t579;
t424 = t561 + t591;
t409 = mrSges(6,2) * t434 + mrSges(7,2) * t429 - mrSges(6,3) * t426 - mrSges(7,3) * t421 - qJ(6) * t418 + t586 * t461 + t598 * t462 + t597 * t495 - t578 * t499 - t593 * t515;
t408 = qJDD(2) * mrSges(5,2) + qJD(2) * t507 - t551;
t407 = -t496 * mrSges(5,2) - t517 * t507 - t552;
t403 = -mrSges(6,1) * t434 + mrSges(6,3) * t427 - mrSges(7,1) * t429 + mrSges(7,3) * t423 - pkin(5) * t424 + qJ(6) * t562 + (-qJ(6) * t478 - t577) * t515 + t578 * t500 + (-mrSges(7,2) * qJ(6) + t595) * t495 + t586 * t462 + t602 * t461;
t401 = mrSges(5,1) * t438 + mrSges(4,2) * t474 - mrSges(4,3) * t441 - mrSges(5,3) * t440 + pkin(4) * t410 - qJ(4) * t407 - t573 * qJD(2) + t585 * qJDD(2) + t583 * t496 + t600 * t497 + t574 * t517 + t548;
t400 = -mrSges(4,1) * t474 + mrSges(4,3) * t442 - mrSges(5,1) * t437 + mrSges(5,2) * t440 - t539 * t409 - t542 * t403 - pkin(4) * (t554 - t590) - pkin(8) * t558 - pkin(3) * t407 + t574 * t518 - t583 * t497 + t596 * t496 + t584 * qJDD(2) + t572 * qJD(2);
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t560 - mrSges(2,2) * t557 + t540 * (mrSges(3,2) * t523 - mrSges(3,3) * t502 + Ifges(3,1) * t527 + Ifges(3,4) * t528 + Ifges(3,5) * qJDD(2) - qJ(3) * t402 - qJD(2) * t520 - t538 * t400 + t581 * t401) + t543 * (-mrSges(3,1) * t523 + mrSges(3,3) * t503 + Ifges(3,4) * t527 + Ifges(3,2) * t528 + Ifges(3,6) * qJDD(2) - pkin(2) * t406 + qJ(3) * t559 + qJD(2) * t521 + t581 * t400 + t538 * t401) + pkin(1) * (-m(3) * t523 + t528 * mrSges(3,1) - t527 * mrSges(3,2) + (-t530 * t540 + t531 * t543) * qJD(1) - t406) + pkin(7) * (t543 * (m(3) * t503 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t528 - qJD(2) * t530 + t526 * t569 + t559) - t540 * (m(3) * t502 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t527 + qJD(2) * t531 - t526 * t570 + t402)); t542 * t409 - t539 * t403 + Ifges(3,5) * t527 + Ifges(3,6) * t528 + qJ(4) * (-t461 * mrSges(6,1) - t499 * t476 + t550 + t591) + mrSges(3,1) * t502 - mrSges(3,2) * t503 - mrSges(5,3) * t437 + mrSges(5,2) * t438 + mrSges(4,1) * t441 - mrSges(4,2) * t442 - pkin(8) * t410 - pkin(3) * t408 + pkin(2) * t402 + t573 * t518 + (-qJ(4) * t489 + t572) * t517 + t585 * t497 + (-mrSges(5,1) * qJ(4) - t584) * t496 + (t540 * t520 - t543 * t521) * qJD(1) + (Ifges(3,3) - t599) * qJDD(2); t406; t408; t548; t424;];
tauJ  = t1;
