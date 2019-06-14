% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-05-07 05:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:52:33
% EndTime: 2019-05-07 05:52:41
% DurationCPUTime: 4.88s
% Computational Cost: add. (50980->330), mult. (103911->403), div. (0->0), fcn. (70766->10), ass. (0->132)
t590 = Ifges(4,1) + Ifges(5,1);
t582 = Ifges(4,4) - Ifges(5,5);
t581 = Ifges(4,5) + Ifges(5,4);
t589 = Ifges(4,2) + Ifges(5,3);
t588 = Ifges(5,2) + Ifges(4,3);
t580 = -Ifges(5,6) + Ifges(4,6);
t550 = sin(qJ(3));
t551 = sin(qJ(2));
t574 = qJD(1) * t551;
t585 = cos(qJ(3));
t527 = -qJD(2) * t585 + t550 * t574;
t554 = cos(qJ(2));
t572 = qJD(1) * qJD(2);
t570 = t554 * t572;
t531 = qJDD(1) * t551 + t570;
t494 = -t527 * qJD(3) + t550 * qJDD(2) + t531 * t585;
t528 = t550 * qJD(2) + t574 * t585;
t501 = mrSges(5,1) * t527 - mrSges(5,3) * t528;
t557 = qJD(1) ^ 2;
t552 = sin(qJ(1));
t555 = cos(qJ(1));
t569 = g(1) * t552 - t555 * g(2);
t519 = -qJDD(1) * pkin(1) - pkin(7) * t557 - t569;
t540 = t551 * t572;
t532 = qJDD(1) * t554 - t540;
t468 = (-t531 - t570) * pkin(8) + (-t532 + t540) * pkin(2) + t519;
t565 = -g(1) * t555 - g(2) * t552;
t520 = -pkin(1) * t557 + qJDD(1) * pkin(7) + t565;
t509 = -g(3) * t551 + t554 * t520;
t530 = (-pkin(2) * t554 - pkin(8) * t551) * qJD(1);
t556 = qJD(2) ^ 2;
t573 = qJD(1) * t554;
t472 = -pkin(2) * t556 + qJDD(2) * pkin(8) + t530 * t573 + t509;
t451 = t468 * t585 - t550 * t472;
t500 = pkin(3) * t527 - qJ(4) * t528;
t526 = -qJDD(3) + t532;
t538 = -qJD(3) + t573;
t537 = t538 ^ 2;
t444 = t526 * pkin(3) - t537 * qJ(4) + t528 * t500 + qJDD(4) - t451;
t579 = t527 * t538;
t429 = (-t494 + t579) * qJ(5) + (t527 * t528 + t526) * pkin(4) + t444;
t452 = t550 * t468 + t585 * t472;
t586 = -2 * qJD(4);
t442 = -pkin(3) * t537 - t526 * qJ(4) - t527 * t500 + t538 * t586 + t452;
t493 = qJD(3) * t528 - qJDD(2) * t585 + t531 * t550;
t504 = pkin(4) * t538 - qJ(5) * t528;
t525 = t527 ^ 2;
t433 = -pkin(4) * t525 + qJ(5) * t493 - t504 * t538 + t442;
t546 = sin(pkin(10));
t547 = cos(pkin(10));
t497 = t527 * t546 + t528 * t547;
t423 = -0.2e1 * qJD(5) * t497 + t547 * t429 - t433 * t546;
t461 = t493 * t546 + t494 * t547;
t496 = t527 * t547 - t528 * t546;
t421 = (t496 * t538 - t461) * pkin(9) + (t496 * t497 + t526) * pkin(5) + t423;
t424 = 0.2e1 * qJD(5) * t496 + t546 * t429 + t547 * t433;
t460 = t493 * t547 - t494 * t546;
t475 = pkin(5) * t538 - pkin(9) * t497;
t495 = t496 ^ 2;
t422 = -pkin(5) * t495 + pkin(9) * t460 - t475 * t538 + t424;
t549 = sin(qJ(6));
t553 = cos(qJ(6));
t419 = t421 * t553 - t422 * t549;
t464 = t496 * t553 - t497 * t549;
t439 = qJD(6) * t464 + t460 * t549 + t461 * t553;
t465 = t496 * t549 + t497 * t553;
t450 = -mrSges(7,1) * t464 + mrSges(7,2) * t465;
t536 = qJD(6) + t538;
t453 = -mrSges(7,2) * t536 + mrSges(7,3) * t464;
t524 = qJDD(6) + t526;
t415 = m(7) * t419 + mrSges(7,1) * t524 - mrSges(7,3) * t439 - t450 * t465 + t453 * t536;
t420 = t421 * t549 + t422 * t553;
t438 = -qJD(6) * t465 + t460 * t553 - t461 * t549;
t454 = mrSges(7,1) * t536 - mrSges(7,3) * t465;
t416 = m(7) * t420 - mrSges(7,2) * t524 + mrSges(7,3) * t438 + t450 * t464 - t454 * t536;
t408 = t553 * t415 + t549 * t416;
t466 = -mrSges(6,1) * t496 + mrSges(6,2) * t497;
t473 = -mrSges(6,2) * t538 + mrSges(6,3) * t496;
t406 = m(6) * t423 + mrSges(6,1) * t526 - mrSges(6,3) * t461 - t466 * t497 + t473 * t538 + t408;
t474 = mrSges(6,1) * t538 - mrSges(6,3) * t497;
t566 = -t415 * t549 + t553 * t416;
t407 = m(6) * t424 - mrSges(6,2) * t526 + mrSges(6,3) * t460 + t466 * t496 - t474 * t538 + t566;
t404 = t406 * t547 + t407 * t546;
t507 = -mrSges(5,2) * t527 - mrSges(5,3) * t538;
t561 = -m(5) * t444 - t526 * mrSges(5,1) - t538 * t507 - t404;
t403 = mrSges(5,2) * t494 + t501 * t528 - t561;
t458 = Ifges(6,4) * t497 + Ifges(6,2) * t496 + Ifges(6,6) * t538;
t459 = Ifges(6,1) * t497 + Ifges(6,4) * t496 + Ifges(6,5) * t538;
t446 = Ifges(7,4) * t465 + Ifges(7,2) * t464 + Ifges(7,6) * t536;
t447 = Ifges(7,1) * t465 + Ifges(7,4) * t464 + Ifges(7,5) * t536;
t560 = -mrSges(7,1) * t419 + mrSges(7,2) * t420 - Ifges(7,5) * t439 - Ifges(7,6) * t438 - Ifges(7,3) * t524 - t465 * t446 + t464 * t447;
t506 = mrSges(5,1) * t538 + mrSges(5,2) * t528;
t567 = -t406 * t546 + t547 * t407;
t563 = m(5) * t442 - t526 * mrSges(5,3) - t538 * t506 + t567;
t576 = -t582 * t527 + t528 * t590 - t581 * t538;
t578 = t527 * t589 - t528 * t582 + t538 * t580;
t587 = -t580 * t493 + t581 * t494 - (Ifges(6,3) + t588) * t526 + t576 * t527 - t578 * t528 + mrSges(4,1) * t451 - mrSges(5,1) * t444 - mrSges(6,1) * t423 - mrSges(4,2) * t452 + mrSges(6,2) * t424 + mrSges(5,3) * t442 - Ifges(6,5) * t461 - Ifges(6,6) * t460 - pkin(3) * t403 - pkin(4) * t404 - pkin(5) * t408 + qJ(4) * (-mrSges(5,2) * t493 - t501 * t527 + t563) - t497 * t458 + t496 * t459 + t560;
t584 = pkin(3) * t538;
t583 = -mrSges(4,3) - mrSges(5,2);
t577 = t527 * t580 - t528 * t581 + t538 * t588;
t575 = -mrSges(4,1) * t527 - mrSges(4,2) * t528 - t501;
t508 = -t554 * g(3) - t551 * t520;
t505 = -mrSges(4,1) * t538 - mrSges(4,3) * t528;
t400 = m(4) * t452 + mrSges(4,2) * t526 + t493 * t583 + t505 * t538 + t527 * t575 + t563;
t503 = mrSges(4,2) * t538 - mrSges(4,3) * t527;
t401 = m(4) * t451 - mrSges(4,1) * t526 + t494 * t583 - t503 * t538 + t528 * t575 + t561;
t568 = t585 * t400 - t401 * t550;
t471 = -qJDD(2) * pkin(2) - t556 * pkin(8) + t530 * t574 - t508;
t562 = t493 * pkin(3) + t471 + (-t494 - t579) * qJ(4);
t432 = -pkin(4) * t493 - qJ(5) * t525 + qJDD(5) - t562 + ((2 * qJD(4)) + t504 + t584) * t528;
t426 = -pkin(5) * t460 - pkin(9) * t495 + t475 * t497 + t432;
t564 = m(7) * t426 - t438 * mrSges(7,1) + t439 * mrSges(7,2) - t464 * t453 + t465 * t454;
t396 = t550 * t400 + t401 * t585;
t417 = m(6) * t432 - t460 * mrSges(6,1) + t461 * mrSges(6,2) - t496 * t473 + t497 * t474 + t564;
t443 = (t586 - t584) * t528 + t562;
t411 = m(5) * t443 + t493 * mrSges(5,1) - t494 * mrSges(5,3) - t528 * t506 + t527 * t507 - t417;
t559 = -m(4) * t471 - t493 * mrSges(4,1) - t494 * mrSges(4,2) - t527 * t503 - t528 * t505 - t411;
t534 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t573;
t533 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t574;
t529 = (-mrSges(3,1) * t554 + mrSges(3,2) * t551) * qJD(1);
t518 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t551 + Ifges(3,4) * t554) * qJD(1);
t517 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t551 + Ifges(3,2) * t554) * qJD(1);
t516 = Ifges(3,3) * qJD(2) + (t551 * Ifges(3,5) + t554 * Ifges(3,6)) * qJD(1);
t457 = Ifges(6,5) * t497 + Ifges(6,6) * t496 + Ifges(6,3) * t538;
t445 = Ifges(7,5) * t465 + Ifges(7,6) * t464 + Ifges(7,3) * t536;
t410 = mrSges(7,2) * t426 - mrSges(7,3) * t419 + Ifges(7,1) * t439 + Ifges(7,4) * t438 + Ifges(7,5) * t524 + t445 * t464 - t446 * t536;
t409 = -mrSges(7,1) * t426 + mrSges(7,3) * t420 + Ifges(7,4) * t439 + Ifges(7,2) * t438 + Ifges(7,6) * t524 - t445 * t465 + t447 * t536;
t398 = mrSges(6,2) * t432 - mrSges(6,3) * t423 + Ifges(6,1) * t461 + Ifges(6,4) * t460 + Ifges(6,5) * t526 - pkin(9) * t408 - t409 * t549 + t410 * t553 + t457 * t496 - t458 * t538;
t397 = -mrSges(6,1) * t432 + mrSges(6,3) * t424 + Ifges(6,4) * t461 + Ifges(6,2) * t460 + Ifges(6,6) * t526 - pkin(5) * t564 + pkin(9) * t566 + t553 * t409 + t549 * t410 - t497 * t457 + t538 * t459;
t395 = mrSges(4,2) * t471 + mrSges(5,2) * t444 - mrSges(4,3) * t451 - mrSges(5,3) * t443 - qJ(4) * t411 - qJ(5) * t404 - t397 * t546 + t398 * t547 - t582 * t493 + t494 * t590 - t581 * t526 + t577 * t527 - t578 * t538;
t394 = -mrSges(4,1) * t471 - mrSges(5,1) * t443 + mrSges(5,2) * t442 + mrSges(4,3) * t452 - pkin(3) * t411 + pkin(4) * t417 - qJ(5) * t567 - t547 * t397 - t546 * t398 - t493 * t589 + t582 * t494 - t580 * t526 + t577 * t528 - t576 * t538;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t569 - mrSges(2,2) * t565 + t551 * (mrSges(3,2) * t519 - mrSges(3,3) * t508 + Ifges(3,1) * t531 + Ifges(3,4) * t532 + Ifges(3,5) * qJDD(2) - pkin(8) * t396 - qJD(2) * t517 - t550 * t394 + t585 * t395 + t516 * t573) + t554 * (-mrSges(3,1) * t519 + mrSges(3,3) * t509 + Ifges(3,4) * t531 + Ifges(3,2) * t532 + Ifges(3,6) * qJDD(2) - pkin(2) * t396 + qJD(2) * t518 - t516 * t574 - t587) + pkin(1) * (-m(3) * t519 + t532 * mrSges(3,1) - t531 * mrSges(3,2) + (-t533 * t551 + t534 * t554) * qJD(1) - t396) + pkin(7) * (t554 * (m(3) * t509 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t532 - qJD(2) * t533 + t529 * t573 + t568) - t551 * (m(3) * t508 + qJDD(2) * mrSges(3,1) - t531 * mrSges(3,3) + qJD(2) * t534 - t529 * t574 + t559)); Ifges(3,5) * t531 + Ifges(3,6) * t532 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t508 - mrSges(3,2) * t509 + t550 * t395 + t585 * t394 + pkin(2) * t559 + pkin(8) * t568 + (t551 * t517 - t554 * t518) * qJD(1); t587; t403; t417; -t560;];
tauJ  = t1;
