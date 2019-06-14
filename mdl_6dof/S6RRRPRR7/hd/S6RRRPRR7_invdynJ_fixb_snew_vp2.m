% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 11:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:31:23
% EndTime: 2019-05-07 11:31:50
% DurationCPUTime: 18.06s
% Computational Cost: add. (276887->362), mult. (613635->471), div. (0->0), fcn. (502604->14), ass. (0->152)
t554 = sin(qJ(2));
t559 = cos(qJ(2));
t548 = sin(pkin(6));
t578 = qJD(1) * t548;
t535 = (-pkin(2) * t559 - pkin(9) * t554) * t578;
t550 = cos(pkin(6));
t544 = qJD(1) * t550 + qJD(2);
t542 = t544 ^ 2;
t543 = qJDD(1) * t550 + qJDD(2);
t577 = qJD(1) * t559;
t561 = qJD(1) ^ 2;
t555 = sin(qJ(1));
t560 = cos(qJ(1));
t568 = -g(1) * t560 - g(2) * t555;
t576 = qJDD(1) * t548;
t533 = -pkin(1) * t561 + pkin(8) * t576 + t568;
t573 = t555 * g(1) - g(2) * t560;
t585 = pkin(8) * t548;
t532 = qJDD(1) * pkin(1) + t561 * t585 + t573;
t582 = t532 * t550;
t579 = t559 * t533 + t554 * t582;
t492 = -t542 * pkin(2) + t543 * pkin(9) + (-g(3) * t554 + t535 * t577) * t548 + t579;
t536 = (qJD(2) * t577 + qJDD(1) * t554) * t548;
t575 = t554 * t578;
t537 = -qJD(2) * t575 + t559 * t576;
t584 = t550 * g(3);
t493 = -t537 * pkin(2) - t536 * pkin(9) - t584 + (-t532 + (pkin(2) * t554 - pkin(9) * t559) * t544 * qJD(1)) * t548;
t553 = sin(qJ(3));
t558 = cos(qJ(3));
t470 = -t553 * t492 + t558 * t493;
t524 = t544 * t558 - t553 * t575;
t504 = qJD(3) * t524 + t536 * t558 + t543 * t553;
t525 = t544 * t553 + t558 * t575;
t529 = qJDD(3) - t537;
t574 = t548 * t577;
t540 = qJD(3) - t574;
t455 = (t524 * t540 - t504) * qJ(4) + (t524 * t525 + t529) * pkin(3) + t470;
t471 = t558 * t492 + t553 * t493;
t503 = -qJD(3) * t525 - t536 * t553 + t543 * t558;
t514 = pkin(3) * t540 - qJ(4) * t525;
t523 = t524 ^ 2;
t458 = -pkin(3) * t523 + qJ(4) * t503 - t514 * t540 + t471;
t547 = sin(pkin(12));
t549 = cos(pkin(12));
t511 = t524 * t547 + t525 * t549;
t435 = -0.2e1 * qJD(4) * t511 + t549 * t455 - t547 * t458;
t478 = t503 * t547 + t504 * t549;
t510 = t524 * t549 - t525 * t547;
t432 = (t510 * t540 - t478) * pkin(10) + (t510 * t511 + t529) * pkin(4) + t435;
t436 = 0.2e1 * qJD(4) * t510 + t547 * t455 + t549 * t458;
t477 = t503 * t549 - t504 * t547;
t496 = pkin(4) * t540 - pkin(10) * t511;
t509 = t510 ^ 2;
t434 = -pkin(4) * t509 + pkin(10) * t477 - t496 * t540 + t436;
t552 = sin(qJ(5));
t557 = cos(qJ(5));
t429 = t552 * t432 + t557 * t434;
t486 = t510 * t552 + t511 * t557;
t453 = -qJD(5) * t486 + t477 * t557 - t478 * t552;
t485 = t510 * t557 - t511 * t552;
t468 = -mrSges(6,1) * t485 + mrSges(6,2) * t486;
t539 = qJD(5) + t540;
t475 = mrSges(6,1) * t539 - mrSges(6,3) * t486;
t528 = qJDD(5) + t529;
t469 = -pkin(5) * t485 - pkin(11) * t486;
t538 = t539 ^ 2;
t426 = -pkin(5) * t538 + pkin(11) * t528 + t469 * t485 + t429;
t580 = t548 * t559;
t505 = -g(3) * t580 - t554 * t533 + t559 * t582;
t491 = -t543 * pkin(2) - t542 * pkin(9) + t535 * t575 - t505;
t462 = -t503 * pkin(3) - t523 * qJ(4) + t525 * t514 + qJDD(4) + t491;
t438 = -t477 * pkin(4) - t509 * pkin(10) + t511 * t496 + t462;
t454 = qJD(5) * t485 + t477 * t552 + t478 * t557;
t430 = t438 + (t486 * t539 - t453) * pkin(5) + (-t485 * t539 - t454) * pkin(11);
t551 = sin(qJ(6));
t556 = cos(qJ(6));
t423 = -t426 * t551 + t430 * t556;
t472 = -t486 * t551 + t539 * t556;
t441 = qJD(6) * t472 + t454 * t556 + t528 * t551;
t452 = qJDD(6) - t453;
t473 = t486 * t556 + t539 * t551;
t459 = -mrSges(7,1) * t472 + mrSges(7,2) * t473;
t484 = qJD(6) - t485;
t460 = -mrSges(7,2) * t484 + mrSges(7,3) * t472;
t419 = m(7) * t423 + mrSges(7,1) * t452 - mrSges(7,3) * t441 - t459 * t473 + t460 * t484;
t424 = t426 * t556 + t430 * t551;
t440 = -qJD(6) * t473 - t454 * t551 + t528 * t556;
t461 = mrSges(7,1) * t484 - mrSges(7,3) * t473;
t420 = m(7) * t424 - mrSges(7,2) * t452 + mrSges(7,3) * t440 + t459 * t472 - t461 * t484;
t569 = -t419 * t551 + t556 * t420;
t404 = m(6) * t429 - mrSges(6,2) * t528 + mrSges(6,3) * t453 + t468 * t485 - t475 * t539 + t569;
t428 = t432 * t557 - t434 * t552;
t474 = -mrSges(6,2) * t539 + mrSges(6,3) * t485;
t425 = -pkin(5) * t528 - pkin(11) * t538 + t469 * t486 - t428;
t566 = -m(7) * t425 + t440 * mrSges(7,1) - mrSges(7,2) * t441 + t472 * t460 - t461 * t473;
t415 = m(6) * t428 + mrSges(6,1) * t528 - mrSges(6,3) * t454 - t468 * t486 + t474 * t539 + t566;
t401 = t552 * t404 + t557 * t415;
t487 = -mrSges(5,1) * t510 + mrSges(5,2) * t511;
t494 = -mrSges(5,2) * t540 + mrSges(5,3) * t510;
t399 = m(5) * t435 + mrSges(5,1) * t529 - mrSges(5,3) * t478 - t487 * t511 + t494 * t540 + t401;
t495 = mrSges(5,1) * t540 - mrSges(5,3) * t511;
t570 = t557 * t404 - t415 * t552;
t400 = m(5) * t436 - mrSges(5,2) * t529 + mrSges(5,3) * t477 + t487 * t510 - t495 * t540 + t570;
t393 = t549 * t399 + t547 * t400;
t480 = Ifges(5,4) * t511 + Ifges(5,2) * t510 + Ifges(5,6) * t540;
t481 = Ifges(5,1) * t511 + Ifges(5,4) * t510 + Ifges(5,5) * t540;
t498 = Ifges(4,4) * t525 + Ifges(4,2) * t524 + Ifges(4,6) * t540;
t499 = Ifges(4,1) * t525 + Ifges(4,4) * t524 + Ifges(4,5) * t540;
t442 = Ifges(7,5) * t473 + Ifges(7,6) * t472 + Ifges(7,3) * t484;
t444 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t484;
t412 = -mrSges(7,1) * t425 + mrSges(7,3) * t424 + Ifges(7,4) * t441 + Ifges(7,2) * t440 + Ifges(7,6) * t452 - t442 * t473 + t444 * t484;
t443 = Ifges(7,4) * t473 + Ifges(7,2) * t472 + Ifges(7,6) * t484;
t413 = mrSges(7,2) * t425 - mrSges(7,3) * t423 + Ifges(7,1) * t441 + Ifges(7,4) * t440 + Ifges(7,5) * t452 + t442 * t472 - t443 * t484;
t464 = Ifges(6,4) * t486 + Ifges(6,2) * t485 + Ifges(6,6) * t539;
t465 = Ifges(6,1) * t486 + Ifges(6,4) * t485 + Ifges(6,5) * t539;
t565 = -mrSges(6,1) * t428 + mrSges(6,2) * t429 - Ifges(6,5) * t454 - Ifges(6,6) * t453 - Ifges(6,3) * t528 - pkin(5) * t566 - pkin(11) * t569 - t556 * t412 - t551 * t413 - t486 * t464 + t485 * t465;
t586 = mrSges(4,1) * t470 + mrSges(5,1) * t435 - mrSges(4,2) * t471 - mrSges(5,2) * t436 + Ifges(4,5) * t504 + Ifges(5,5) * t478 + Ifges(4,6) * t503 + Ifges(5,6) * t477 + pkin(3) * t393 + pkin(4) * t401 + t511 * t480 - t510 * t481 + t525 * t498 - t524 * t499 + (Ifges(4,3) + Ifges(5,3)) * t529 - t565;
t581 = t548 * t554;
t512 = -mrSges(4,1) * t524 + mrSges(4,2) * t525;
t513 = -mrSges(4,2) * t540 + mrSges(4,3) * t524;
t391 = m(4) * t470 + mrSges(4,1) * t529 - mrSges(4,3) * t504 - t512 * t525 + t513 * t540 + t393;
t515 = mrSges(4,1) * t540 - mrSges(4,3) * t525;
t571 = -t399 * t547 + t549 * t400;
t392 = m(4) * t471 - mrSges(4,2) * t529 + mrSges(4,3) * t503 + t512 * t524 - t515 * t540 + t571;
t386 = t558 * t391 + t553 * t392;
t408 = t556 * t419 + t551 * t420;
t572 = -t391 * t553 + t558 * t392;
t567 = m(6) * t438 - t453 * mrSges(6,1) + t454 * mrSges(6,2) - t485 * t474 + t486 * t475 + t408;
t406 = m(5) * t462 - t477 * mrSges(5,1) + t478 * mrSges(5,2) - t510 * t494 + t511 * t495 + t567;
t564 = mrSges(7,1) * t423 - mrSges(7,2) * t424 + Ifges(7,5) * t441 + Ifges(7,6) * t440 + Ifges(7,3) * t452 + t443 * t473 - t444 * t472;
t563 = -m(4) * t491 + t503 * mrSges(4,1) - t504 * mrSges(4,2) + t524 * t513 - t525 * t515 - t406;
t534 = (-mrSges(3,1) * t559 + mrSges(3,2) * t554) * t578;
t531 = -mrSges(3,2) * t544 + mrSges(3,3) * t574;
t530 = mrSges(3,1) * t544 - mrSges(3,3) * t575;
t519 = -t548 * t532 - t584;
t518 = Ifges(3,5) * t544 + (Ifges(3,1) * t554 + Ifges(3,4) * t559) * t578;
t517 = Ifges(3,6) * t544 + (Ifges(3,4) * t554 + Ifges(3,2) * t559) * t578;
t516 = Ifges(3,3) * t544 + (Ifges(3,5) * t554 + Ifges(3,6) * t559) * t578;
t506 = -g(3) * t581 + t579;
t497 = Ifges(4,5) * t525 + Ifges(4,6) * t524 + Ifges(4,3) * t540;
t479 = Ifges(5,5) * t511 + Ifges(5,6) * t510 + Ifges(5,3) * t540;
t463 = Ifges(6,5) * t486 + Ifges(6,6) * t485 + Ifges(6,3) * t539;
t405 = m(3) * t505 + t543 * mrSges(3,1) - t536 * mrSges(3,3) + t544 * t531 - t534 * t575 + t563;
t395 = -mrSges(6,1) * t438 + mrSges(6,3) * t429 + Ifges(6,4) * t454 + Ifges(6,2) * t453 + Ifges(6,6) * t528 - pkin(5) * t408 - t463 * t486 + t465 * t539 - t564;
t394 = mrSges(6,2) * t438 - mrSges(6,3) * t428 + Ifges(6,1) * t454 + Ifges(6,4) * t453 + Ifges(6,5) * t528 - pkin(11) * t408 - t412 * t551 + t413 * t556 + t463 * t485 - t464 * t539;
t387 = mrSges(5,2) * t462 - mrSges(5,3) * t435 + Ifges(5,1) * t478 + Ifges(5,4) * t477 + Ifges(5,5) * t529 - pkin(10) * t401 + t394 * t557 - t395 * t552 + t479 * t510 - t480 * t540;
t385 = m(3) * t506 - mrSges(3,2) * t543 + mrSges(3,3) * t537 - t530 * t544 + t534 * t574 + t572;
t384 = -mrSges(5,1) * t462 + mrSges(5,3) * t436 + Ifges(5,4) * t478 + Ifges(5,2) * t477 + Ifges(5,6) * t529 - pkin(4) * t567 + pkin(10) * t570 + t552 * t394 + t557 * t395 - t511 * t479 + t540 * t481;
t383 = mrSges(4,2) * t491 - mrSges(4,3) * t470 + Ifges(4,1) * t504 + Ifges(4,4) * t503 + Ifges(4,5) * t529 - qJ(4) * t393 - t384 * t547 + t387 * t549 + t497 * t524 - t498 * t540;
t382 = -mrSges(4,1) * t491 + mrSges(4,3) * t471 + Ifges(4,4) * t504 + Ifges(4,2) * t503 + Ifges(4,6) * t529 - pkin(3) * t406 + qJ(4) * t571 + t549 * t384 + t547 * t387 - t525 * t497 + t540 * t499;
t381 = Ifges(3,5) * t536 + Ifges(3,6) * t537 + Ifges(3,3) * t543 + mrSges(3,1) * t505 - mrSges(3,2) * t506 + t553 * t383 + t558 * t382 + pkin(2) * t563 + pkin(9) * t572 + (t517 * t554 - t518 * t559) * t578;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t573 - mrSges(2,2) * t568 + (mrSges(3,2) * t519 - mrSges(3,3) * t505 + Ifges(3,1) * t536 + Ifges(3,4) * t537 + Ifges(3,5) * t543 - pkin(9) * t386 - t382 * t553 + t383 * t558 + t516 * t574 - t517 * t544) * t581 + (-mrSges(3,1) * t519 + mrSges(3,3) * t506 + Ifges(3,4) * t536 + Ifges(3,2) * t537 + Ifges(3,6) * t543 - pkin(2) * t386 - t516 * t575 + t544 * t518 - t586) * t580 + t550 * t381 + pkin(1) * ((t385 * t554 + t405 * t559) * t550 + (-m(3) * t519 + t537 * mrSges(3,1) - t536 * mrSges(3,2) + (-t530 * t554 + t531 * t559) * t578 - t386) * t548) + (t385 * t559 - t405 * t554) * t585; t381; t586; t406; -t565; t564;];
tauJ  = t1;
