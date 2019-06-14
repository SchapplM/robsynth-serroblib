% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-05-06 10:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:52:59
% EndTime: 2019-05-06 10:53:05
% DurationCPUTime: 4.20s
% Computational Cost: add. (35853->333), mult. (81943->410), div. (0->0), fcn. (54260->10), ass. (0->131)
t582 = Ifges(3,1) + Ifges(4,1);
t577 = Ifges(3,4) - Ifges(4,5);
t576 = Ifges(3,5) + Ifges(4,4);
t581 = Ifges(3,2) + Ifges(4,3);
t575 = Ifges(3,6) - Ifges(4,6);
t580 = (Ifges(3,3) + Ifges(4,2));
t579 = 2 * qJD(3);
t578 = mrSges(3,3) + mrSges(4,2);
t551 = cos(qJ(2));
t554 = qJD(1) ^ 2;
t574 = t551 ^ 2 * t554;
t548 = sin(qJ(1));
t552 = cos(qJ(1));
t562 = -g(1) * t552 - g(2) * t548;
t511 = -pkin(1) * t554 + qJDD(1) * pkin(7) + t562;
t547 = sin(qJ(2));
t490 = -g(3) * t547 + t551 * t511;
t515 = (-pkin(2) * t551 - qJ(3) * t547) * qJD(1);
t553 = qJD(2) ^ 2;
t568 = qJD(1) * t551;
t471 = -pkin(2) * t553 + qJDD(2) * qJ(3) + (qJD(2) * t579) + t515 * t568 + t490;
t567 = qJD(1) * qJD(2);
t565 = t547 * t567;
t519 = qJDD(1) * t551 - t565;
t569 = qJD(1) * t547;
t522 = -(qJD(2) * pkin(3)) - qJ(4) * t569;
t467 = -pkin(3) * t574 - qJ(4) * t519 + qJD(2) * t522 + t471;
t489 = -t551 * g(3) - t547 * t511;
t474 = -qJDD(2) * pkin(2) - t553 * qJ(3) + t515 * t569 + qJDD(3) - t489;
t566 = t551 * t567;
t518 = qJDD(1) * t547 + t566;
t468 = (-t518 + t566) * qJ(4) + (-t547 * t551 * t554 - qJDD(2)) * pkin(3) + t474;
t542 = sin(pkin(10));
t543 = cos(pkin(10));
t503 = (-t542 * t551 + t543 * t547) * qJD(1);
t437 = -0.2e1 * qJD(4) * t503 - t542 * t467 + t543 * t468;
t488 = t518 * t543 - t519 * t542;
t502 = (-t542 * t547 - t543 * t551) * qJD(1);
t434 = (-qJD(2) * t502 - t488) * pkin(8) + (t502 * t503 - qJDD(2)) * pkin(4) + t437;
t438 = 0.2e1 * qJD(4) * t502 + t543 * t467 + t542 * t468;
t487 = -t518 * t542 - t519 * t543;
t493 = -qJD(2) * pkin(4) - pkin(8) * t503;
t501 = t502 ^ 2;
t436 = -pkin(4) * t501 + pkin(8) * t487 + qJD(2) * t493 + t438;
t546 = sin(qJ(5));
t550 = cos(qJ(5));
t431 = t546 * t434 + t550 * t436;
t481 = t502 * t550 - t503 * t546;
t482 = t502 * t546 + t503 * t550;
t466 = -pkin(5) * t481 - pkin(9) * t482;
t536 = -qJD(2) + qJD(5);
t534 = t536 ^ 2;
t535 = -qJDD(2) + qJDD(5);
t429 = -pkin(5) * t534 + pkin(9) * t535 + t466 * t481 + t431;
t570 = t548 * g(1) - t552 * g(2);
t510 = -qJDD(1) * pkin(1) - t554 * pkin(7) - t570;
t560 = -t519 * pkin(2) + t510 + (-t518 - t566) * qJ(3);
t454 = -pkin(2) * t565 + t519 * pkin(3) - qJ(4) * t574 + qJDD(4) - t560 + (t522 + t579) * t569;
t440 = -t487 * pkin(4) - t501 * pkin(8) + t503 * t493 + t454;
t451 = -qJD(5) * t482 + t487 * t550 - t488 * t546;
t452 = qJD(5) * t481 + t487 * t546 + t488 * t550;
t432 = (-t481 * t536 - t452) * pkin(9) + (t482 * t536 - t451) * pkin(5) + t440;
t545 = sin(qJ(6));
t549 = cos(qJ(6));
t426 = -t429 * t545 + t432 * t549;
t472 = -t482 * t545 + t536 * t549;
t443 = qJD(6) * t472 + t452 * t549 + t535 * t545;
t450 = qJDD(6) - t451;
t473 = t482 * t549 + t536 * t545;
t453 = -mrSges(7,1) * t472 + mrSges(7,2) * t473;
t477 = qJD(6) - t481;
t455 = -mrSges(7,2) * t477 + mrSges(7,3) * t472;
t423 = m(7) * t426 + mrSges(7,1) * t450 - mrSges(7,3) * t443 - t453 * t473 + t455 * t477;
t427 = t429 * t549 + t432 * t545;
t442 = -qJD(6) * t473 - t452 * t545 + t535 * t549;
t456 = mrSges(7,1) * t477 - mrSges(7,3) * t473;
t424 = m(7) * t427 - mrSges(7,2) * t450 + mrSges(7,3) * t442 + t453 * t472 - t456 * t477;
t415 = -t423 * t545 + t549 * t424;
t465 = -mrSges(6,1) * t481 + mrSges(6,2) * t482;
t476 = mrSges(6,1) * t536 - mrSges(6,3) * t482;
t412 = m(6) * t431 - mrSges(6,2) * t535 + mrSges(6,3) * t451 + t465 * t481 - t476 * t536 + t415;
t430 = t434 * t550 - t436 * t546;
t428 = -pkin(5) * t535 - pkin(9) * t534 + t466 * t482 - t430;
t425 = -m(7) * t428 + t442 * mrSges(7,1) - mrSges(7,2) * t443 + t472 * t455 - t456 * t473;
t475 = -mrSges(6,2) * t536 + mrSges(6,3) * t481;
t419 = m(6) * t430 + mrSges(6,1) * t535 - mrSges(6,3) * t452 - t465 * t482 + t475 * t536 + t425;
t408 = t546 * t412 + t550 * t419;
t414 = t549 * t423 + t545 * t424;
t573 = (t580 * qJD(2)) + (t547 * t576 + t551 * t575) * qJD(1);
t572 = -t575 * qJD(2) + (-t547 * t577 - t581 * t551) * qJD(1);
t571 = t576 * qJD(2) + (t547 * t582 + t551 * t577) * qJD(1);
t485 = -mrSges(5,1) * t502 + mrSges(5,2) * t503;
t491 = qJD(2) * mrSges(5,2) + mrSges(5,3) * t502;
t406 = m(5) * t437 - qJDD(2) * mrSges(5,1) - mrSges(5,3) * t488 - qJD(2) * t491 - t485 * t503 + t408;
t492 = -qJD(2) * mrSges(5,1) - mrSges(5,3) * t503;
t563 = t550 * t412 - t419 * t546;
t407 = m(5) * t438 + qJDD(2) * mrSges(5,2) + mrSges(5,3) * t487 + qJD(2) * t492 + t485 * t502 + t563;
t564 = -t542 * t406 + t543 * t407;
t402 = t543 * t406 + t542 * t407;
t561 = m(6) * t440 - t451 * mrSges(6,1) + t452 * mrSges(6,2) - t481 * t475 + t482 * t476 + t414;
t516 = (-mrSges(4,1) * t551 - mrSges(4,3) * t547) * qJD(1);
t524 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t569;
t559 = m(4) * t471 + qJDD(2) * mrSges(4,3) + qJD(2) * t524 + t516 * t568 + t564;
t526 = mrSges(4,2) * t568 + qJD(2) * mrSges(4,3);
t558 = m(4) * t474 - qJDD(2) * mrSges(4,1) - qJD(2) * t526 + t402;
t413 = m(5) * t454 - t487 * mrSges(5,1) + t488 * mrSges(5,2) - t502 * t491 + t503 * t492 + t561;
t445 = Ifges(7,4) * t473 + Ifges(7,2) * t472 + Ifges(7,6) * t477;
t446 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t477;
t557 = mrSges(7,1) * t426 - mrSges(7,2) * t427 + Ifges(7,5) * t443 + Ifges(7,6) * t442 + Ifges(7,3) * t450 + t445 * t473 - t446 * t472;
t469 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t569 + t560;
t556 = m(4) * t469 - t413;
t444 = Ifges(7,5) * t473 + Ifges(7,6) * t472 + Ifges(7,3) * t477;
t416 = -mrSges(7,1) * t428 + mrSges(7,3) * t427 + Ifges(7,4) * t443 + Ifges(7,2) * t442 + Ifges(7,6) * t450 - t444 * t473 + t446 * t477;
t417 = mrSges(7,2) * t428 - mrSges(7,3) * t426 + Ifges(7,1) * t443 + Ifges(7,4) * t442 + Ifges(7,5) * t450 + t444 * t472 - t445 * t477;
t458 = Ifges(6,4) * t482 + Ifges(6,2) * t481 + (Ifges(6,6) * t536);
t459 = Ifges(6,1) * t482 + Ifges(6,4) * t481 + (Ifges(6,5) * t536);
t555 = mrSges(6,1) * t430 - mrSges(6,2) * t431 + Ifges(6,5) * t452 + Ifges(6,6) * t451 + Ifges(6,3) * t535 + pkin(5) * t425 + pkin(9) * t415 + t549 * t416 + t545 * t417 + t482 * t458 - t481 * t459;
t525 = -(qJD(2) * mrSges(3,2)) + mrSges(3,3) * t568;
t523 = (qJD(2) * mrSges(3,1)) - mrSges(3,3) * t569;
t517 = (-mrSges(3,1) * t551 + mrSges(3,2) * t547) * qJD(1);
t480 = Ifges(5,1) * t503 + Ifges(5,4) * t502 - (Ifges(5,5) * qJD(2));
t479 = Ifges(5,4) * t503 + Ifges(5,2) * t502 - (Ifges(5,6) * qJD(2));
t478 = Ifges(5,5) * t503 + Ifges(5,6) * t502 - (Ifges(5,3) * qJD(2));
t457 = Ifges(6,5) * t482 + Ifges(6,6) * t481 + Ifges(6,3) * t536;
t409 = t556 + (-t547 * t524 - t551 * t526) * qJD(1) - t519 * mrSges(4,1) - t518 * mrSges(4,3);
t404 = -mrSges(6,1) * t440 + mrSges(6,3) * t431 + Ifges(6,4) * t452 + Ifges(6,2) * t451 + Ifges(6,6) * t535 - pkin(5) * t414 - t457 * t482 + t459 * t536 - t557;
t403 = mrSges(6,2) * t440 - mrSges(6,3) * t430 + Ifges(6,1) * t452 + Ifges(6,4) * t451 + Ifges(6,5) * t535 - pkin(9) * t414 - t416 * t545 + t417 * t549 + t457 * t481 - t458 * t536;
t401 = t518 * mrSges(4,2) + t516 * t569 + t558;
t400 = mrSges(5,2) * t454 - mrSges(5,3) * t437 + Ifges(5,1) * t488 + Ifges(5,4) * t487 - Ifges(5,5) * qJDD(2) - pkin(8) * t408 + qJD(2) * t479 + t403 * t550 - t404 * t546 + t478 * t502;
t399 = -mrSges(5,1) * t454 + mrSges(5,3) * t438 + Ifges(5,4) * t488 + Ifges(5,2) * t487 - Ifges(5,6) * qJDD(2) - pkin(4) * t561 + pkin(8) * t563 - qJD(2) * t480 + t546 * t403 + t550 * t404 - t503 * t478;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t570 - mrSges(2,2) * t562 + t547 * (mrSges(3,2) * t510 + mrSges(4,2) * t474 - mrSges(3,3) * t489 - mrSges(4,3) * t469 - qJ(3) * t409 - qJ(4) * t402 + t572 * qJD(2) + t576 * qJDD(2) - t542 * t399 + t543 * t400 + t582 * t518 + t577 * t519 + t573 * t568) + t551 * (-mrSges(3,1) * t510 - mrSges(4,1) * t469 + mrSges(4,2) * t471 + mrSges(3,3) * t490 - pkin(2) * t409 + pkin(3) * t413 - qJ(4) * t564 + t571 * qJD(2) + t575 * qJDD(2) - t543 * t399 - t542 * t400 + t577 * t518 + t581 * t519 - t573 * t569) + pkin(1) * (-t556 + ((t525 + t526) * t551 + (-t523 + t524) * t547) * qJD(1) + (mrSges(3,1) + mrSges(4,1)) * t519 + (-mrSges(3,2) + mrSges(4,3)) * t518 - m(3) * t510) + pkin(7) * (t551 * (m(3) * t490 - qJDD(2) * mrSges(3,2) - qJD(2) * t523 + t517 * t568 + t578 * t519 + t559) + (-m(3) * t489 - qJDD(2) * mrSges(3,1) - qJD(2) * t525 + t578 * t518 + (t516 + t517) * t569 + t558) * t547); -t555 + (-t572 * t547 - t571 * t551) * qJD(1) + qJ(3) * t559 + (mrSges(4,2) * qJ(3) + t575) * t519 - pkin(2) * t401 + t576 * t518 + (Ifges(5,3) + t580) * qJDD(2) + t502 * t480 - t503 * t479 + mrSges(3,1) * t489 - mrSges(3,2) * t490 - Ifges(5,6) * t487 - Ifges(5,5) * t488 + mrSges(4,3) * t471 - mrSges(4,1) * t474 - mrSges(5,1) * t437 + mrSges(5,2) * t438 - pkin(4) * t408 - pkin(3) * t402; t401; t413; t555; t557;];
tauJ  = t1;
