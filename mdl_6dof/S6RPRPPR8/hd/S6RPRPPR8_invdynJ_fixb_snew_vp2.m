% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-05-05 17:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:22:22
% EndTime: 2019-05-05 17:22:25
% DurationCPUTime: 1.62s
% Computational Cost: add. (3949->264), mult. (7889->293), div. (0->0), fcn. (3515->6), ass. (0->116)
t508 = sin(qJ(3));
t511 = cos(qJ(3));
t537 = Ifges(5,5) - Ifges(6,4) - Ifges(4,4);
t572 = t511 * (Ifges(4,1) + Ifges(5,1) + Ifges(6,2)) + t508 * t537;
t571 = (Ifges(5,3) + Ifges(6,1) + Ifges(4,2)) * t508 + t511 * t537;
t538 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t536 = Ifges(4,6) - Ifges(5,6) + Ifges(6,5);
t566 = t511 * (t571 * qJD(1) - t536 * qJD(3));
t514 = qJD(1) ^ 2;
t509 = sin(qJ(1));
t512 = cos(qJ(1));
t533 = g(1) * t509 - t512 * g(2);
t525 = -qJ(2) * t514 + qJDD(2) - t533;
t445 = (-pkin(1) - pkin(7)) * qJDD(1) + t525;
t440 = -g(3) * t511 + t508 * t445;
t565 = -qJDD(3) * qJ(4) - t440;
t472 = (t508 * pkin(3) - t511 * qJ(4)) * qJD(1);
t543 = qJD(1) * t511;
t564 = t472 * t543 + qJDD(4);
t551 = t445 * t511;
t557 = g(3) * t508;
t439 = t551 + t557;
t513 = qJD(3) ^ 2;
t552 = qJ(4) * t513;
t427 = -qJDD(3) * pkin(3) - t439 - t552 + t564;
t544 = qJD(1) * t508;
t487 = -mrSges(5,2) * t544 + qJD(3) * mrSges(5,3);
t563 = m(5) * t427 - qJDD(3) * mrSges(5,1) - qJD(3) * t487;
t562 = -2 * qJD(5);
t561 = -pkin(3) - pkin(8);
t560 = -pkin(4) - pkin(8);
t558 = pkin(3) * t513;
t556 = mrSges(4,1) - mrSges(6,2);
t555 = -mrSges(4,3) - mrSges(5,2);
t554 = -pkin(5) - qJ(4);
t540 = qJD(1) * qJD(3);
t476 = qJDD(1) * t508 + t511 * t540;
t553 = t476 * mrSges(6,2);
t481 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t544;
t550 = t481 * t508;
t549 = t508 ^ 2 * t514;
t492 = t508 * t540;
t477 = qJDD(1) * t511 - t492;
t483 = -qJD(3) * pkin(4) - qJ(5) * t543;
t539 = qJD(2) * qJD(1);
t495 = -0.2e1 * t539;
t545 = -t512 * g(1) - t509 * g(2);
t531 = t514 * pkin(1) - qJDD(1) * qJ(2) - t545;
t529 = -t514 * pkin(7) - t531;
t526 = t476 * pkin(3) - t477 * qJ(4) + t529;
t542 = qJD(4) * t511;
t516 = -qJ(5) * t549 + 0.2e1 * qJD(1) * t542 + t483 * t543 + qJDD(5) + t495 - t526;
t416 = t516 + t560 * t476 + pkin(5) * t477 + (t508 * t554 + t511 * t561) * t540;
t475 = (t511 * pkin(5) - t508 * pkin(8)) * qJD(1);
t518 = t511 * t514 * t508 * pkin(4) + t543 * t562 + (-t477 - t492) * qJ(5) - t557 + t564;
t419 = t554 * t513 + (-qJD(1) * t475 - t445) * t511 + (-pkin(3) + t560) * qJDD(3) + t518;
t507 = sin(qJ(6));
t510 = cos(qJ(6));
t414 = t416 * t510 - t419 * t507;
t469 = -qJD(3) * t510 - t507 * t544;
t437 = qJD(6) * t469 - qJDD(3) * t507 + t476 * t510;
t470 = -qJD(3) * t507 + t510 * t544;
t438 = -mrSges(7,1) * t469 + mrSges(7,2) * t470;
t489 = qJD(6) + t543;
t441 = -mrSges(7,2) * t489 + mrSges(7,3) * t469;
t467 = qJDD(6) + t477;
t411 = m(7) * t414 + mrSges(7,1) * t467 - mrSges(7,3) * t437 - t438 * t470 + t441 * t489;
t415 = t416 * t507 + t419 * t510;
t436 = -qJD(6) * t470 - qJDD(3) * t510 - t476 * t507;
t442 = mrSges(7,1) * t489 - mrSges(7,3) * t470;
t412 = m(7) * t415 - mrSges(7,2) * t467 + mrSges(7,3) * t436 + t438 * t469 - t442 * t489;
t548 = -t507 * t411 + t510 * t412;
t547 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t544 - t481;
t484 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t543;
t486 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t543;
t546 = -t484 - t486;
t541 = t562 + t472;
t534 = t572 * qJD(1) + t538 * qJD(3);
t473 = (t508 * mrSges(5,1) - t511 * mrSges(5,3)) * qJD(1);
t532 = qJD(1) * (-t473 - (mrSges(4,1) * t508 + mrSges(4,2) * t511) * qJD(1));
t530 = pkin(3) * t511 + qJ(4) * t508;
t485 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t543;
t493 = 0.2e1 * qJD(4) * qJD(3);
t426 = -t472 * t544 + t493 - t558 - t565;
t521 = pkin(4) * t549 - qJ(5) * t476 + t565;
t421 = t558 + (-0.2e1 * qJD(4) - t483) * qJD(3) + t541 * t544 + t521;
t471 = (t511 * mrSges(6,1) + t508 * mrSges(6,2)) * qJD(1);
t418 = qJDD(3) * pkin(5) + qJD(3) * t483 + t493 + t561 * t513 + (t475 - t541) * t544 - t521;
t524 = -m(7) * t418 + mrSges(7,1) * t436 - t437 * mrSges(7,2) + t441 * t469 - t470 * t442;
t517 = -m(6) * t421 + qJDD(3) * mrSges(6,1) + t476 * mrSges(6,3) + qJD(3) * t484 + t471 * t544 - t524;
t515 = m(5) * t426 + qJDD(3) * mrSges(5,3) + qJD(3) * t486 + t517;
t422 = -t552 - t551 + (-pkin(3) - pkin(4)) * qJDD(3) + t518;
t527 = -m(6) * t422 + t471 * t543 - t548;
t528 = (m(4) * t439 + t556 * qJDD(3) + t547 * qJD(3) + t511 * t532 + (mrSges(6,3) + t555) * t477 + t527 - t563) * t511 + (m(4) * t440 - qJDD(3) * mrSges(4,2) - qJD(3) * t485 + t476 * t555 + t508 * t532 + t515) * t508;
t405 = t510 * t411 + t507 * t412;
t420 = -pkin(4) * t476 - t530 * t540 + t516;
t523 = m(6) * t420 + t477 * mrSges(6,1) + t405;
t430 = Ifges(7,4) * t470 + Ifges(7,2) * t469 + Ifges(7,6) * t489;
t431 = Ifges(7,1) * t470 + Ifges(7,4) * t469 + Ifges(7,5) * t489;
t522 = mrSges(7,1) * t414 - mrSges(7,2) * t415 + Ifges(7,5) * t437 + Ifges(7,6) * t436 + Ifges(7,3) * t467 + t470 * t430 - t469 * t431;
t520 = qJDD(3) * mrSges(6,2) + qJD(3) * t481 - t527;
t494 = 0.2e1 * t539;
t424 = t494 + (qJD(3) * t530 - 0.2e1 * t542) * qJD(1) + t526;
t519 = m(5) * t424 + t476 * mrSges(5,1) + t487 * t544 - t523;
t447 = -qJDD(1) * pkin(1) + t525;
t446 = t495 + t531;
t444 = t494 + t529;
t429 = Ifges(7,5) * t470 + Ifges(7,6) * t469 + Ifges(7,3) * t489;
t408 = mrSges(7,2) * t418 - mrSges(7,3) * t414 + Ifges(7,1) * t437 + Ifges(7,4) * t436 + Ifges(7,5) * t467 + t429 * t469 - t430 * t489;
t407 = -mrSges(7,1) * t418 + mrSges(7,3) * t415 + Ifges(7,4) * t437 + Ifges(7,2) * t436 + Ifges(7,6) * t467 - t429 * t470 + t431 * t489;
t404 = -mrSges(6,3) * t477 + t520;
t403 = t553 + (t484 * t511 + t550) * qJD(1) + t523;
t402 = t473 * t543 + (mrSges(5,2) - mrSges(6,3)) * t477 + t520 + t563;
t401 = -t553 - t477 * mrSges(5,3) + (t511 * t546 - t550) * qJD(1) + t519;
t399 = m(3) * t447 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t514 + t528;
t1 = [mrSges(2,1) * t533 - mrSges(2,2) * t545 + mrSges(3,2) * t447 - mrSges(3,3) * t446 + t511 * (mrSges(6,1) * t420 + mrSges(4,2) * t444 + mrSges(5,2) * t427 - mrSges(4,3) * t439 - mrSges(5,3) * t424 - mrSges(6,3) * t422 + pkin(5) * t405 - qJ(4) * t401 - qJ(5) * t404 + t522) - t508 * (-mrSges(4,1) * t444 - mrSges(5,1) * t424 + mrSges(5,2) * t426 - mrSges(6,2) * t420 + mrSges(4,3) * t440 + mrSges(6,3) * t421 - pkin(3) * t401 + pkin(4) * t403 + pkin(8) * t405 - qJ(5) * t517 + t507 * t407 - t510 * t408) - pkin(7) * t528 - pkin(1) * t399 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t508 * t536 + t511 * t538) * qJDD(3) + (-t508 * t534 + t566) * qJD(3) + t572 * t477 + t571 * t476 + (-m(3) * t446 + m(4) * t444 + t514 * mrSges(3,2) + t519 + qJDD(1) * mrSges(3,3) + (mrSges(4,2) - mrSges(5,3)) * t477 + t556 * t476 + (t547 * t508 + (t485 + t546) * t511) * qJD(1)) * qJ(2); t399; -mrSges(6,1) * t421 + mrSges(6,2) * t422 + mrSges(5,3) * t426 - mrSges(5,1) * t427 - pkin(3) * t402 - pkin(5) * t524 - pkin(8) * t548 - t507 * t408 + qJ(4) * t515 - t510 * t407 - pkin(4) * t404 + mrSges(4,1) * t439 - mrSges(4,2) * t440 + t538 * t477 + (-qJ(4) * mrSges(5,2) - t536) * t476 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * qJDD(3) + (-t566 + (-qJ(4) * t473 + t534) * t508) * qJD(1); t402; t403; t522;];
tauJ  = t1;
