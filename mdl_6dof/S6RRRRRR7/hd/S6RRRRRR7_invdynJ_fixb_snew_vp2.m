% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRRRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 12:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 12:20:33
% EndTime: 2019-05-08 12:21:01
% DurationCPUTime: 20.18s
% Computational Cost: add. (330366->364), mult. (703920->471), div. (0->0), fcn. (573709->14), ass. (0->155)
t547 = sin(pkin(6));
t553 = sin(qJ(2));
t559 = cos(qJ(2));
t577 = qJD(1) * qJD(2);
t538 = (-qJDD(1) * t559 + t553 * t577) * t547;
t585 = pkin(8) * t547;
t548 = cos(pkin(6));
t584 = t548 * g(3);
t561 = qJD(1) ^ 2;
t554 = sin(qJ(1));
t560 = cos(qJ(1));
t574 = t554 * g(1) - g(2) * t560;
t533 = qJDD(1) * pkin(1) + t561 * t585 + t574;
t583 = t533 * t548;
t582 = t547 * t553;
t581 = t547 * t559;
t579 = qJD(1) * t547;
t536 = (-pkin(2) * t559 - pkin(9) * t553) * t579;
t544 = qJD(1) * t548 + qJD(2);
t542 = t544 ^ 2;
t543 = qJDD(1) * t548 + qJDD(2);
t578 = qJD(1) * t559;
t570 = -g(1) * t560 - g(2) * t554;
t534 = -pkin(1) * t561 + qJDD(1) * t585 + t570;
t580 = t559 * t534 + t553 * t583;
t484 = -t542 * pkin(2) + t543 * pkin(9) + (-g(3) * t553 + t536 * t578) * t547 + t580;
t537 = (qJDD(1) * t553 + t559 * t577) * t547;
t485 = t538 * pkin(2) - t537 * pkin(9) - t584 + (-t533 + (pkin(2) * t553 - pkin(9) * t559) * t544 * qJD(1)) * t547;
t552 = sin(qJ(3));
t558 = cos(qJ(3));
t461 = t558 * t484 + t552 * t485;
t576 = t553 * t579;
t525 = t544 * t558 - t552 * t576;
t526 = t544 * t552 + t558 * t576;
t509 = -pkin(3) * t525 - pkin(10) * t526;
t530 = qJDD(3) + t538;
t575 = t547 * t578;
t541 = qJD(3) - t575;
t539 = t541 ^ 2;
t453 = -pkin(3) * t539 + pkin(10) * t530 + t509 * t525 + t461;
t506 = -g(3) * t581 - t553 * t534 + t559 * t583;
t483 = -t543 * pkin(2) - t542 * pkin(9) + t536 * t576 - t506;
t504 = -t526 * qJD(3) - t552 * t537 + t543 * t558;
t505 = qJD(3) * t525 + t537 * t558 + t543 * t552;
t457 = (-t525 * t541 - t505) * pkin(10) + (t526 * t541 - t504) * pkin(3) + t483;
t551 = sin(qJ(4));
t557 = cos(qJ(4));
t436 = -t551 * t453 + t557 * t457;
t511 = -t526 * t551 + t541 * t557;
t471 = qJD(4) * t511 + t505 * t557 + t530 * t551;
t502 = qJDD(4) - t504;
t512 = t526 * t557 + t541 * t551;
t524 = qJD(4) - t525;
t427 = (t511 * t524 - t471) * pkin(11) + (t511 * t512 + t502) * pkin(4) + t436;
t437 = t557 * t453 + t551 * t457;
t470 = -qJD(4) * t512 - t505 * t551 + t530 * t557;
t493 = pkin(4) * t524 - pkin(11) * t512;
t510 = t511 ^ 2;
t435 = -pkin(4) * t510 + pkin(11) * t470 - t493 * t524 + t437;
t550 = sin(qJ(5));
t556 = cos(qJ(5));
t421 = t556 * t427 - t550 * t435;
t487 = t511 * t556 - t512 * t550;
t450 = qJD(5) * t487 + t470 * t550 + t471 * t556;
t488 = t511 * t550 + t512 * t556;
t497 = qJDD(5) + t502;
t522 = qJD(5) + t524;
t418 = (t487 * t522 - t450) * pkin(12) + (t487 * t488 + t497) * pkin(5) + t421;
t422 = t550 * t427 + t556 * t435;
t449 = -qJD(5) * t488 + t470 * t556 - t471 * t550;
t474 = pkin(5) * t522 - pkin(12) * t488;
t486 = t487 ^ 2;
t419 = -pkin(5) * t486 + pkin(12) * t449 - t474 * t522 + t422;
t549 = sin(qJ(6));
t555 = cos(qJ(6));
t416 = t418 * t555 - t419 * t549;
t466 = t487 * t555 - t488 * t549;
t433 = qJD(6) * t466 + t449 * t549 + t450 * t555;
t467 = t487 * t549 + t488 * t555;
t445 = -mrSges(7,1) * t466 + mrSges(7,2) * t467;
t519 = qJD(6) + t522;
t458 = -mrSges(7,2) * t519 + mrSges(7,3) * t466;
t496 = qJDD(6) + t497;
t412 = m(7) * t416 + mrSges(7,1) * t496 - mrSges(7,3) * t433 - t445 * t467 + t458 * t519;
t417 = t418 * t549 + t419 * t555;
t432 = -qJD(6) * t467 + t449 * t555 - t450 * t549;
t459 = mrSges(7,1) * t519 - mrSges(7,3) * t467;
t413 = m(7) * t417 - mrSges(7,2) * t496 + mrSges(7,3) * t432 + t445 * t466 - t459 * t519;
t404 = t555 * t412 + t549 * t413;
t468 = -mrSges(6,1) * t487 + mrSges(6,2) * t488;
t472 = -mrSges(6,2) * t522 + mrSges(6,3) * t487;
t401 = m(6) * t421 + mrSges(6,1) * t497 - mrSges(6,3) * t450 - t468 * t488 + t472 * t522 + t404;
t473 = mrSges(6,1) * t522 - mrSges(6,3) * t488;
t571 = -t412 * t549 + t555 * t413;
t402 = m(6) * t422 - mrSges(6,2) * t497 + mrSges(6,3) * t449 + t468 * t487 - t473 * t522 + t571;
t397 = t556 * t401 + t550 * t402;
t489 = -mrSges(5,1) * t511 + mrSges(5,2) * t512;
t491 = -mrSges(5,2) * t524 + mrSges(5,3) * t511;
t395 = m(5) * t436 + mrSges(5,1) * t502 - mrSges(5,3) * t471 - t489 * t512 + t491 * t524 + t397;
t492 = mrSges(5,1) * t524 - mrSges(5,3) * t512;
t572 = -t401 * t550 + t556 * t402;
t396 = m(5) * t437 - mrSges(5,2) * t502 + mrSges(5,3) * t470 + t489 * t511 - t492 * t524 + t572;
t391 = -t395 * t551 + t557 * t396;
t508 = -mrSges(4,1) * t525 + mrSges(4,2) * t526;
t514 = mrSges(4,1) * t541 - mrSges(4,3) * t526;
t389 = m(4) * t461 - mrSges(4,2) * t530 + mrSges(4,3) * t504 + t508 * t525 - t514 * t541 + t391;
t460 = -t552 * t484 + t485 * t558;
t452 = -pkin(3) * t530 - pkin(10) * t539 + t526 * t509 - t460;
t438 = -pkin(4) * t470 - pkin(11) * t510 + t512 * t493 + t452;
t424 = -pkin(5) * t449 - pkin(12) * t486 + t474 * t488 + t438;
t569 = m(7) * t424 - t432 * mrSges(7,1) + t433 * mrSges(7,2) - t466 * t458 + t467 * t459;
t566 = m(6) * t438 - t449 * mrSges(6,1) + mrSges(6,2) * t450 - t487 * t472 + t473 * t488 + t569;
t414 = -m(5) * t452 + t470 * mrSges(5,1) - mrSges(5,2) * t471 + t511 * t491 - t492 * t512 - t566;
t513 = -mrSges(4,2) * t541 + mrSges(4,3) * t525;
t408 = m(4) * t460 + mrSges(4,1) * t530 - mrSges(4,3) * t505 - t508 * t526 + t513 * t541 + t414;
t385 = t552 * t389 + t558 * t408;
t573 = t558 * t389 - t408 * t552;
t390 = t395 * t557 + t396 * t551;
t441 = Ifges(7,4) * t467 + Ifges(7,2) * t466 + Ifges(7,6) * t519;
t442 = Ifges(7,1) * t467 + Ifges(7,4) * t466 + Ifges(7,5) * t519;
t567 = -mrSges(7,1) * t416 + mrSges(7,2) * t417 - Ifges(7,5) * t433 - Ifges(7,6) * t432 - Ifges(7,3) * t496 - t467 * t441 + t466 * t442;
t565 = -m(4) * t483 + t504 * mrSges(4,1) - mrSges(4,2) * t505 + t525 * t513 - t514 * t526 - t390;
t463 = Ifges(6,4) * t488 + Ifges(6,2) * t487 + Ifges(6,6) * t522;
t464 = Ifges(6,1) * t488 + Ifges(6,4) * t487 + Ifges(6,5) * t522;
t564 = -mrSges(6,1) * t421 + mrSges(6,2) * t422 - Ifges(6,5) * t450 - Ifges(6,6) * t449 - Ifges(6,3) * t497 - pkin(5) * t404 - t488 * t463 + t487 * t464 + t567;
t440 = Ifges(7,5) * t467 + Ifges(7,6) * t466 + Ifges(7,3) * t519;
t405 = -mrSges(7,1) * t424 + mrSges(7,3) * t417 + Ifges(7,4) * t433 + Ifges(7,2) * t432 + Ifges(7,6) * t496 - t440 * t467 + t442 * t519;
t406 = mrSges(7,2) * t424 - mrSges(7,3) * t416 + Ifges(7,1) * t433 + Ifges(7,4) * t432 + Ifges(7,5) * t496 + t440 * t466 - t441 * t519;
t462 = Ifges(6,5) * t488 + Ifges(6,6) * t487 + Ifges(6,3) * t522;
t392 = -mrSges(6,1) * t438 + mrSges(6,3) * t422 + Ifges(6,4) * t450 + Ifges(6,2) * t449 + Ifges(6,6) * t497 - pkin(5) * t569 + pkin(12) * t571 + t555 * t405 + t549 * t406 - t488 * t462 + t522 * t464;
t393 = mrSges(6,2) * t438 - mrSges(6,3) * t421 + Ifges(6,1) * t450 + Ifges(6,4) * t449 + Ifges(6,5) * t497 - pkin(12) * t404 - t405 * t549 + t406 * t555 + t462 * t487 - t463 * t522;
t475 = Ifges(5,5) * t512 + Ifges(5,6) * t511 + Ifges(5,3) * t524;
t477 = Ifges(5,1) * t512 + Ifges(5,4) * t511 + Ifges(5,5) * t524;
t382 = -mrSges(5,1) * t452 + mrSges(5,3) * t437 + Ifges(5,4) * t471 + Ifges(5,2) * t470 + Ifges(5,6) * t502 - pkin(4) * t566 + pkin(11) * t572 + t556 * t392 + t550 * t393 - t512 * t475 + t524 * t477;
t476 = Ifges(5,4) * t512 + Ifges(5,2) * t511 + Ifges(5,6) * t524;
t383 = mrSges(5,2) * t452 - mrSges(5,3) * t436 + Ifges(5,1) * t471 + Ifges(5,4) * t470 + Ifges(5,5) * t502 - pkin(11) * t397 - t392 * t550 + t393 * t556 + t475 * t511 - t476 * t524;
t499 = Ifges(4,4) * t526 + Ifges(4,2) * t525 + Ifges(4,6) * t541;
t500 = Ifges(4,1) * t526 + Ifges(4,4) * t525 + Ifges(4,5) * t541;
t563 = mrSges(4,1) * t460 - mrSges(4,2) * t461 + Ifges(4,5) * t505 + Ifges(4,6) * t504 + Ifges(4,3) * t530 + pkin(3) * t414 + pkin(10) * t391 + t557 * t382 + t551 * t383 + t526 * t499 - t525 * t500;
t562 = mrSges(5,1) * t436 - mrSges(5,2) * t437 + Ifges(5,5) * t471 + Ifges(5,6) * t470 + Ifges(5,3) * t502 + pkin(4) * t397 + t512 * t476 - t511 * t477 - t564;
t535 = (-mrSges(3,1) * t559 + mrSges(3,2) * t553) * t579;
t532 = -mrSges(3,2) * t544 + mrSges(3,3) * t575;
t531 = mrSges(3,1) * t544 - mrSges(3,3) * t576;
t518 = -t547 * t533 - t584;
t517 = Ifges(3,5) * t544 + (Ifges(3,1) * t553 + Ifges(3,4) * t559) * t579;
t516 = Ifges(3,6) * t544 + (Ifges(3,4) * t553 + Ifges(3,2) * t559) * t579;
t515 = Ifges(3,3) * t544 + (Ifges(3,5) * t553 + Ifges(3,6) * t559) * t579;
t507 = -g(3) * t582 + t580;
t498 = Ifges(4,5) * t526 + Ifges(4,6) * t525 + Ifges(4,3) * t541;
t386 = m(3) * t506 + mrSges(3,1) * t543 - mrSges(3,3) * t537 + t532 * t544 - t535 * t576 + t565;
t384 = m(3) * t507 - mrSges(3,2) * t543 - mrSges(3,3) * t538 - t531 * t544 + t535 * t575 + t573;
t381 = -mrSges(4,1) * t483 + mrSges(4,3) * t461 + Ifges(4,4) * t505 + Ifges(4,2) * t504 + Ifges(4,6) * t530 - pkin(3) * t390 - t526 * t498 + t541 * t500 - t562;
t380 = mrSges(4,2) * t483 - mrSges(4,3) * t460 + Ifges(4,1) * t505 + Ifges(4,4) * t504 + Ifges(4,5) * t530 - pkin(10) * t390 - t382 * t551 + t383 * t557 + t498 * t525 - t499 * t541;
t379 = Ifges(3,5) * t537 - Ifges(3,6) * t538 + Ifges(3,3) * t543 + mrSges(3,1) * t506 - mrSges(3,2) * t507 + t552 * t380 + t558 * t381 + pkin(2) * t565 + pkin(9) * t573 + (t516 * t553 - t517 * t559) * t579;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t574 - mrSges(2,2) * t570 + (mrSges(3,2) * t518 - mrSges(3,3) * t506 + Ifges(3,1) * t537 - Ifges(3,4) * t538 + Ifges(3,5) * t543 - pkin(9) * t385 + t380 * t558 - t381 * t552 + t515 * t575 - t516 * t544) * t582 + (-mrSges(3,1) * t518 + mrSges(3,3) * t507 + Ifges(3,4) * t537 - Ifges(3,2) * t538 + Ifges(3,6) * t543 - pkin(2) * t385 - t515 * t576 + t544 * t517 - t563) * t581 + t548 * t379 + pkin(1) * ((t384 * t553 + t386 * t559) * t548 + (-m(3) * t518 - t538 * mrSges(3,1) - t537 * mrSges(3,2) + (-t531 * t553 + t532 * t559) * t579 - t385) * t547) + (t384 * t559 - t386 * t553) * t585; t379; t563; t562; -t564; -t567;];
tauJ  = t1;
