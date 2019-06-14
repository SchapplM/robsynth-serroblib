% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:03:05
% EndTime: 2019-05-06 23:03:22
% DurationCPUTime: 16.78s
% Computational Cost: add. (259957->362), mult. (591520->471), div. (0->0), fcn. (487005->14), ass. (0->152)
t548 = sin(pkin(6));
t584 = pkin(8) * t548;
t550 = cos(pkin(6));
t583 = t550 * g(3);
t561 = qJD(1) ^ 2;
t555 = sin(qJ(1));
t560 = cos(qJ(1));
t573 = g(1) * t555 - g(2) * t560;
t532 = qJDD(1) * pkin(1) + t561 * t584 + t573;
t582 = t532 * t550;
t554 = sin(qJ(2));
t581 = t548 * t554;
t559 = cos(qJ(2));
t580 = t548 * t559;
t578 = qJD(1) * t548;
t534 = (-pkin(2) * t559 - qJ(3) * t554) * t578;
t544 = qJD(1) * t550 + qJD(2);
t542 = t544 ^ 2;
t543 = qJDD(1) * t550 + qJDD(2);
t577 = qJD(1) * t559;
t568 = -g(1) * t560 - g(2) * t555;
t576 = qJDD(1) * t548;
t533 = -pkin(1) * t561 + pkin(8) * t576 + t568;
t579 = t533 * t559 + t554 * t582;
t492 = -t542 * pkin(2) + t543 * qJ(3) + (-g(3) * t554 + t534 * t577) * t548 + t579;
t536 = (qJD(2) * t577 + qJDD(1) * t554) * t548;
t575 = t554 * t578;
t537 = -qJD(2) * t575 + t559 * t576;
t493 = -t537 * pkin(2) - t583 - t536 * qJ(3) + (-t532 + (pkin(2) * t554 - qJ(3) * t559) * t544 * qJD(1)) * t548;
t547 = sin(pkin(12));
t549 = cos(pkin(12));
t525 = t544 * t547 + t549 * t575;
t467 = -0.2e1 * qJD(3) * t525 - t547 * t492 + t493 * t549;
t512 = t536 * t549 + t543 * t547;
t524 = t544 * t549 - t547 * t575;
t574 = t548 * t577;
t455 = (-t524 * t574 - t512) * pkin(9) + (t524 * t525 - t537) * pkin(3) + t467;
t468 = 0.2e1 * qJD(3) * t524 + t492 * t549 + t493 * t547;
t511 = -t536 * t547 + t543 * t549;
t513 = -pkin(3) * t574 - pkin(9) * t525;
t522 = t524 ^ 2;
t458 = -pkin(3) * t522 + pkin(9) * t511 + t513 * t574 + t468;
t553 = sin(qJ(4));
t558 = cos(qJ(4));
t435 = t455 * t558 - t553 * t458;
t505 = t524 * t558 - t525 * t553;
t476 = qJD(4) * t505 + t511 * t553 + t512 * t558;
t506 = t524 * t553 + t525 * t558;
t529 = qJDD(4) - t537;
t540 = qJD(4) - t574;
t432 = (t505 * t540 - t476) * pkin(10) + (t505 * t506 + t529) * pkin(4) + t435;
t436 = t455 * t553 + t458 * t558;
t475 = -qJD(4) * t506 + t511 * t558 - t512 * t553;
t496 = pkin(4) * t540 - pkin(10) * t506;
t504 = t505 ^ 2;
t434 = -pkin(4) * t504 + pkin(10) * t475 - t496 * t540 + t436;
t552 = sin(qJ(5));
t557 = cos(qJ(5));
t429 = t432 * t552 + t434 * t557;
t486 = t505 * t552 + t506 * t557;
t451 = -qJD(5) * t486 + t475 * t557 - t476 * t552;
t485 = t505 * t557 - t506 * t552;
t469 = -mrSges(6,1) * t485 + mrSges(6,2) * t486;
t539 = qJD(5) + t540;
t478 = mrSges(6,1) * t539 - mrSges(6,3) * t486;
t528 = qJDD(5) + t529;
t470 = -pkin(5) * t485 - pkin(11) * t486;
t538 = t539 ^ 2;
t426 = -pkin(5) * t538 + pkin(11) * t528 + t470 * t485 + t429;
t502 = -g(3) * t580 - t533 * t554 + t559 * t582;
t491 = -t543 * pkin(2) - t542 * qJ(3) + t534 * t575 + qJDD(3) - t502;
t471 = -t511 * pkin(3) - t522 * pkin(9) + t513 * t525 + t491;
t441 = -t475 * pkin(4) - t504 * pkin(10) + t496 * t506 + t471;
t452 = qJD(5) * t485 + t475 * t552 + t476 * t557;
t430 = t441 + (-t485 * t539 - t452) * pkin(11) + (t486 * t539 - t451) * pkin(5);
t551 = sin(qJ(6));
t556 = cos(qJ(6));
t423 = -t426 * t551 + t430 * t556;
t472 = -t486 * t551 + t539 * t556;
t439 = qJD(6) * t472 + t452 * t556 + t528 * t551;
t446 = qJDD(6) - t451;
t473 = t486 * t556 + t539 * t551;
t459 = -mrSges(7,1) * t472 + mrSges(7,2) * t473;
t484 = qJD(6) - t485;
t460 = -mrSges(7,2) * t484 + mrSges(7,3) * t472;
t419 = m(7) * t423 + mrSges(7,1) * t446 - mrSges(7,3) * t439 - t459 * t473 + t460 * t484;
t424 = t426 * t556 + t430 * t551;
t438 = -qJD(6) * t473 - t452 * t551 + t528 * t556;
t461 = mrSges(7,1) * t484 - mrSges(7,3) * t473;
t420 = m(7) * t424 - mrSges(7,2) * t446 + mrSges(7,3) * t438 + t459 * t472 - t461 * t484;
t569 = -t419 * t551 + t420 * t556;
t405 = m(6) * t429 - mrSges(6,2) * t528 + mrSges(6,3) * t451 + t469 * t485 - t478 * t539 + t569;
t428 = t432 * t557 - t434 * t552;
t477 = -mrSges(6,2) * t539 + mrSges(6,3) * t485;
t425 = -pkin(5) * t528 - pkin(11) * t538 + t470 * t486 - t428;
t566 = -m(7) * t425 + mrSges(7,1) * t438 - mrSges(7,2) * t439 + t460 * t472 - t461 * t473;
t415 = m(6) * t428 + mrSges(6,1) * t528 - mrSges(6,3) * t452 - t469 * t486 + t477 * t539 + t566;
t401 = t405 * t552 + t415 * t557;
t487 = -mrSges(5,1) * t505 + mrSges(5,2) * t506;
t494 = -mrSges(5,2) * t540 + mrSges(5,3) * t505;
t399 = m(5) * t435 + mrSges(5,1) * t529 - mrSges(5,3) * t476 - t487 * t506 + t494 * t540 + t401;
t495 = mrSges(5,1) * t540 - mrSges(5,3) * t506;
t570 = t405 * t557 - t415 * t552;
t400 = m(5) * t436 - mrSges(5,2) * t529 + mrSges(5,3) * t475 + t487 * t505 - t495 * t540 + t570;
t393 = t399 * t558 + t400 * t553;
t507 = -mrSges(4,1) * t524 + mrSges(4,2) * t525;
t509 = mrSges(4,2) * t574 + mrSges(4,3) * t524;
t391 = m(4) * t467 - mrSges(4,1) * t537 - mrSges(4,3) * t512 - t507 * t525 - t509 * t574 + t393;
t510 = -mrSges(4,1) * t574 - mrSges(4,3) * t525;
t571 = -t399 * t553 + t400 * t558;
t392 = m(4) * t468 + mrSges(4,2) * t537 + mrSges(4,3) * t511 + t507 * t524 + t510 * t574 + t571;
t386 = t391 * t549 + t392 * t547;
t408 = t419 * t556 + t420 * t551;
t572 = -t391 * t547 + t392 * t549;
t567 = m(6) * t441 - mrSges(6,1) * t451 + mrSges(6,2) * t452 - t477 * t485 + t478 * t486 + t408;
t447 = Ifges(7,5) * t473 + Ifges(7,6) * t472 + Ifges(7,3) * t484;
t449 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t484;
t412 = -mrSges(7,1) * t425 + mrSges(7,3) * t424 + Ifges(7,4) * t439 + Ifges(7,2) * t438 + Ifges(7,6) * t446 - t447 * t473 + t449 * t484;
t448 = Ifges(7,4) * t473 + Ifges(7,2) * t472 + Ifges(7,6) * t484;
t413 = mrSges(7,2) * t425 - mrSges(7,3) * t423 + Ifges(7,1) * t439 + Ifges(7,4) * t438 + Ifges(7,5) * t446 + t447 * t472 - t448 * t484;
t463 = Ifges(6,4) * t486 + Ifges(6,2) * t485 + Ifges(6,6) * t539;
t464 = Ifges(6,1) * t486 + Ifges(6,4) * t485 + Ifges(6,5) * t539;
t565 = -mrSges(6,1) * t428 + mrSges(6,2) * t429 - Ifges(6,5) * t452 - Ifges(6,6) * t451 - Ifges(6,3) * t528 - pkin(5) * t566 - pkin(11) * t569 - t412 * t556 - t413 * t551 - t463 * t486 + t485 * t464;
t564 = m(5) * t471 - mrSges(5,1) * t475 + mrSges(5,2) * t476 - t494 * t505 + t495 * t506 + t567;
t563 = mrSges(7,1) * t423 - mrSges(7,2) * t424 + Ifges(7,5) * t439 + Ifges(7,6) * t438 + Ifges(7,3) * t446 + t448 * t473 - t449 * t472;
t406 = m(4) * t491 - mrSges(4,1) * t511 + mrSges(4,2) * t512 - t509 * t524 + t510 * t525 + t564;
t480 = Ifges(5,4) * t506 + Ifges(5,2) * t505 + Ifges(5,6) * t540;
t481 = Ifges(5,1) * t506 + Ifges(5,4) * t505 + Ifges(5,5) * t540;
t562 = mrSges(5,1) * t435 - mrSges(5,2) * t436 + Ifges(5,5) * t476 + Ifges(5,6) * t475 + Ifges(5,3) * t529 + pkin(4) * t401 + t506 * t480 - t505 * t481 - t565;
t535 = (-mrSges(3,1) * t559 + mrSges(3,2) * t554) * t578;
t531 = -mrSges(3,2) * t544 + mrSges(3,3) * t574;
t530 = mrSges(3,1) * t544 - mrSges(3,3) * t575;
t517 = -t548 * t532 - t583;
t516 = Ifges(3,5) * t544 + (Ifges(3,1) * t554 + Ifges(3,4) * t559) * t578;
t515 = Ifges(3,6) * t544 + (t554 * Ifges(3,4) + Ifges(3,2) * t559) * t578;
t514 = Ifges(3,3) * t544 + (Ifges(3,5) * t554 + Ifges(3,6) * t559) * t578;
t503 = -g(3) * t581 + t579;
t499 = Ifges(4,1) * t525 + Ifges(4,4) * t524 - Ifges(4,5) * t574;
t498 = Ifges(4,4) * t525 + Ifges(4,2) * t524 - Ifges(4,6) * t574;
t497 = Ifges(4,5) * t525 + Ifges(4,6) * t524 - Ifges(4,3) * t574;
t479 = Ifges(5,5) * t506 + Ifges(5,6) * t505 + Ifges(5,3) * t540;
t462 = Ifges(6,5) * t486 + Ifges(6,6) * t485 + Ifges(6,3) * t539;
t402 = m(3) * t502 + mrSges(3,1) * t543 - mrSges(3,3) * t536 + t531 * t544 - t535 * t575 - t406;
t395 = -mrSges(6,1) * t441 + mrSges(6,3) * t429 + Ifges(6,4) * t452 + Ifges(6,2) * t451 + Ifges(6,6) * t528 - pkin(5) * t408 - t462 * t486 + t464 * t539 - t563;
t394 = mrSges(6,2) * t441 - mrSges(6,3) * t428 + Ifges(6,1) * t452 + Ifges(6,4) * t451 + Ifges(6,5) * t528 - pkin(11) * t408 - t412 * t551 + t413 * t556 + t462 * t485 - t463 * t539;
t387 = mrSges(5,2) * t471 - mrSges(5,3) * t435 + Ifges(5,1) * t476 + Ifges(5,4) * t475 + Ifges(5,5) * t529 - pkin(10) * t401 + t394 * t557 - t395 * t552 + t479 * t505 - t480 * t540;
t385 = m(3) * t503 - mrSges(3,2) * t543 + mrSges(3,3) * t537 - t530 * t544 + t535 * t574 + t572;
t384 = -mrSges(5,1) * t471 + mrSges(5,3) * t436 + Ifges(5,4) * t476 + Ifges(5,2) * t475 + Ifges(5,6) * t529 - pkin(4) * t567 + pkin(10) * t570 + t552 * t394 + t557 * t395 - t506 * t479 + t540 * t481;
t383 = mrSges(4,2) * t491 - mrSges(4,3) * t467 + Ifges(4,1) * t512 + Ifges(4,4) * t511 - Ifges(4,5) * t537 - pkin(9) * t393 - t384 * t553 + t387 * t558 + t497 * t524 + t498 * t574;
t382 = -mrSges(4,1) * t491 + mrSges(4,3) * t468 + Ifges(4,4) * t512 + Ifges(4,2) * t511 - Ifges(4,6) * t537 - pkin(3) * t564 + pkin(9) * t571 + t558 * t384 + t553 * t387 - t525 * t497 - t499 * t574;
t381 = Ifges(3,5) * t536 + Ifges(3,6) * t537 + Ifges(3,3) * t543 + mrSges(3,1) * t502 - mrSges(3,2) * t503 + t547 * t383 + t549 * t382 - pkin(2) * t406 + qJ(3) * t572 + (t515 * t554 - t516 * t559) * t578;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t573 - mrSges(2,2) * t568 + (mrSges(3,2) * t517 - mrSges(3,3) * t502 + Ifges(3,1) * t536 + Ifges(3,4) * t537 + Ifges(3,5) * t543 - qJ(3) * t386 - t382 * t547 + t383 * t549 + t514 * t574 - t515 * t544) * t581 + (-t514 * t575 - pkin(2) * t386 - t562 + Ifges(3,6) * t543 + t544 * t516 + Ifges(3,4) * t536 + t524 * t499 - t525 * t498 - Ifges(4,5) * t512 - mrSges(3,1) * t517 - Ifges(4,6) * t511 + mrSges(3,3) * t503 - mrSges(4,1) * t467 + mrSges(4,2) * t468 + (Ifges(3,2) + Ifges(4,3)) * t537 - pkin(3) * t393) * t580 + t550 * t381 + pkin(1) * ((t385 * t554 + t402 * t559) * t550 + (-m(3) * t517 + t537 * mrSges(3,1) - t536 * mrSges(3,2) + (-t530 * t554 + t531 * t559) * t578 - t386) * t548) + (t385 * t559 - t402 * t554) * t584; t381; t406; t562; -t565; t563;];
tauJ  = t1;
