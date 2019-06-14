% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 10:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:40:13
% EndTime: 2019-05-07 10:40:24
% DurationCPUTime: 5.55s
% Computational Cost: add. (53074->329), mult. (109133->400), div. (0->0), fcn. (76419->10), ass. (0->135)
t579 = Ifges(4,4) + Ifges(5,6);
t590 = Ifges(4,2) + Ifges(5,3);
t586 = Ifges(4,6) - Ifges(5,5);
t589 = -2 * qJD(4);
t588 = Ifges(4,1) + Ifges(5,2);
t587 = Ifges(4,5) - Ifges(5,4);
t585 = Ifges(4,3) + Ifges(5,1);
t545 = sin(qJ(3));
t550 = cos(qJ(2));
t572 = qJD(1) * t550;
t546 = sin(qJ(2));
t573 = qJD(1) * t546;
t582 = cos(qJ(3));
t520 = t545 * t573 - t572 * t582;
t521 = (t545 * t550 + t546 * t582) * qJD(1);
t541 = qJD(2) + qJD(3);
t584 = t590 * t520 - t579 * t521 - t586 * t541;
t570 = qJD(1) * qJD(2);
t527 = qJDD(1) * t546 + t550 * t570;
t552 = qJD(1) ^ 2;
t547 = sin(qJ(1));
t551 = cos(qJ(1));
t565 = -g(1) * t551 - g(2) * t547;
t523 = -pkin(1) * t552 + qJDD(1) * pkin(7) + t565;
t577 = t546 * t523;
t581 = pkin(2) * t552;
t468 = qJDD(2) * pkin(2) - t527 * pkin(8) - t577 + (pkin(8) * t570 + t546 * t581 - g(3)) * t550;
t505 = -g(3) * t546 + t550 * t523;
t528 = qJDD(1) * t550 - t546 * t570;
t531 = qJD(2) * pkin(2) - pkin(8) * t573;
t542 = t550 ^ 2;
t469 = pkin(8) * t528 - qJD(2) * t531 - t542 * t581 + t505;
t447 = t545 * t468 + t469 * t582;
t496 = pkin(3) * t520 - qJ(4) * t521;
t539 = t541 ^ 2;
t540 = qJDD(2) + qJDD(3);
t437 = pkin(3) * t539 - t540 * qJ(4) + t520 * t496 + t541 * t589 - t447;
t483 = qJD(3) * t521 + t527 * t545 - t528 * t582;
t484 = -t520 * qJD(3) + t527 * t582 + t545 * t528;
t569 = t547 * g(1) - t551 * g(2);
t564 = -qJDD(1) * pkin(1) - t569;
t485 = -t528 * pkin(2) + t531 * t573 + (-pkin(8) * t542 - pkin(7)) * t552 + t564;
t507 = mrSges(4,1) * t541 - mrSges(4,3) * t521;
t578 = t520 * t541;
t554 = (-t484 + t578) * qJ(4) + t485 + (pkin(3) * t541 + t589) * t521;
t435 = t483 * pkin(3) + t554;
t509 = mrSges(5,1) * t521 + mrSges(5,2) * t541;
t510 = pkin(4) * t521 - pkin(9) * t541;
t516 = t520 ^ 2;
t423 = -t516 * pkin(4) - t521 * t510 + (pkin(3) + pkin(9)) * t483 + t554;
t446 = t468 * t582 - t545 * t469;
t439 = -t540 * pkin(3) - t539 * qJ(4) + t521 * t496 + qJDD(4) - t446;
t426 = (t520 * t521 - t540) * pkin(9) + (t484 + t578) * pkin(4) + t439;
t544 = sin(qJ(5));
t549 = cos(qJ(5));
t418 = -t544 * t423 + t549 * t426;
t502 = t520 * t549 - t541 * t544;
t453 = qJD(5) * t502 + t483 * t544 + t540 * t549;
t482 = qJDD(5) + t484;
t503 = t520 * t544 + t541 * t549;
t515 = qJD(5) + t521;
t415 = (t502 * t515 - t453) * pkin(10) + (t502 * t503 + t482) * pkin(5) + t418;
t419 = t549 * t423 + t544 * t426;
t452 = -qJD(5) * t503 + t483 * t549 - t540 * t544;
t488 = pkin(5) * t515 - pkin(10) * t503;
t501 = t502 ^ 2;
t416 = -pkin(5) * t501 + pkin(10) * t452 - t488 * t515 + t419;
t543 = sin(qJ(6));
t548 = cos(qJ(6));
t413 = t415 * t548 - t416 * t543;
t460 = t502 * t548 - t503 * t543;
t433 = qJD(6) * t460 + t452 * t543 + t453 * t548;
t461 = t502 * t543 + t503 * t548;
t444 = -mrSges(7,1) * t460 + mrSges(7,2) * t461;
t513 = qJD(6) + t515;
t454 = -mrSges(7,2) * t513 + mrSges(7,3) * t460;
t472 = qJDD(6) + t482;
t409 = m(7) * t413 + mrSges(7,1) * t472 - mrSges(7,3) * t433 - t444 * t461 + t454 * t513;
t414 = t415 * t543 + t416 * t548;
t432 = -qJD(6) * t461 + t452 * t548 - t453 * t543;
t455 = mrSges(7,1) * t513 - mrSges(7,3) * t461;
t410 = m(7) * t414 - mrSges(7,2) * t472 + mrSges(7,3) * t432 + t444 * t460 - t455 * t513;
t400 = t548 * t409 + t543 * t410;
t465 = -mrSges(6,1) * t502 + mrSges(6,2) * t503;
t486 = -mrSges(6,2) * t515 + mrSges(6,3) * t502;
t397 = m(6) * t418 + mrSges(6,1) * t482 - mrSges(6,3) * t453 - t465 * t503 + t486 * t515 + t400;
t487 = mrSges(6,1) * t515 - mrSges(6,3) * t503;
t566 = -t409 * t543 + t548 * t410;
t398 = m(6) * t419 - mrSges(6,2) * t482 + mrSges(6,3) * t452 + t465 * t502 - t487 * t515 + t566;
t567 = -t544 * t397 + t549 * t398;
t562 = -m(5) * t435 + t484 * mrSges(5,3) + t521 * t509 - t567;
t508 = mrSges(5,1) * t520 - mrSges(5,3) * t541;
t574 = -mrSges(4,2) * t541 - mrSges(4,3) * t520 - t508;
t580 = mrSges(4,1) - mrSges(5,2);
t583 = m(4) * t485 + t484 * mrSges(4,2) + t580 * t483 + t521 * t507 + t574 * t520 - t562;
t497 = mrSges(4,1) * t520 + mrSges(4,2) * t521;
t395 = t549 * t397 + t544 * t398;
t498 = -mrSges(5,2) * t520 - mrSges(5,3) * t521;
t560 = -m(5) * t439 - t484 * mrSges(5,1) - t521 * t498 - t395;
t391 = m(4) * t446 - t484 * mrSges(4,3) - t521 * t497 + t540 * t580 + t541 * t574 + t560;
t428 = -pkin(4) * t483 - pkin(9) * t516 + t541 * t510 - t437;
t421 = -pkin(5) * t452 - pkin(10) * t501 + t488 * t503 + t428;
t561 = m(7) * t421 - t432 * mrSges(7,1) + t433 * mrSges(7,2) - t460 * t454 + t461 * t455;
t557 = -m(6) * t428 + t452 * mrSges(6,1) - t453 * mrSges(6,2) + t502 * t486 - t503 * t487 - t561;
t556 = -m(5) * t437 + t540 * mrSges(5,3) + t541 * t509 - t557;
t405 = (-mrSges(4,3) - mrSges(5,1)) * t483 - t541 * t507 - t540 * mrSges(4,2) + m(4) * t447 + t556 + (-t497 - t498) * t520;
t386 = t391 * t582 + t545 * t405;
t576 = t586 * t520 - t587 * t521 - t585 * t541;
t575 = -t579 * t520 + t588 * t521 + t587 * t541;
t568 = -t391 * t545 + t405 * t582;
t441 = Ifges(7,4) * t461 + Ifges(7,2) * t460 + Ifges(7,6) * t513;
t442 = Ifges(7,1) * t461 + Ifges(7,4) * t460 + Ifges(7,5) * t513;
t559 = mrSges(7,1) * t413 - mrSges(7,2) * t414 + Ifges(7,5) * t433 + Ifges(7,6) * t432 + Ifges(7,3) * t472 + t461 * t441 - t460 * t442;
t457 = Ifges(6,4) * t503 + Ifges(6,2) * t502 + Ifges(6,6) * t515;
t458 = Ifges(6,1) * t503 + Ifges(6,4) * t502 + Ifges(6,5) * t515;
t555 = mrSges(6,1) * t418 - mrSges(6,2) * t419 + Ifges(6,5) * t453 + Ifges(6,6) * t452 + Ifges(6,3) * t482 + pkin(5) * t400 + t503 * t457 - t502 * t458 + t559;
t440 = Ifges(7,5) * t461 + Ifges(7,6) * t460 + Ifges(7,3) * t513;
t401 = -mrSges(7,1) * t421 + mrSges(7,3) * t414 + Ifges(7,4) * t433 + Ifges(7,2) * t432 + Ifges(7,6) * t472 - t440 * t461 + t442 * t513;
t402 = mrSges(7,2) * t421 - mrSges(7,3) * t413 + Ifges(7,1) * t433 + Ifges(7,4) * t432 + Ifges(7,5) * t472 + t440 * t460 - t441 * t513;
t456 = Ifges(6,5) * t503 + Ifges(6,6) * t502 + Ifges(6,3) * t515;
t387 = -mrSges(6,1) * t428 + mrSges(6,3) * t419 + Ifges(6,4) * t453 + Ifges(6,2) * t452 + Ifges(6,6) * t482 - pkin(5) * t561 + pkin(10) * t566 + t548 * t401 + t543 * t402 - t503 * t456 + t515 * t458;
t389 = mrSges(6,2) * t428 - mrSges(6,3) * t418 + Ifges(6,1) * t453 + Ifges(6,4) * t452 + Ifges(6,5) * t482 - pkin(10) * t400 - t401 * t543 + t402 * t548 + t456 * t502 - t457 * t515;
t394 = t540 * mrSges(5,2) + t541 * t508 - t560;
t553 = -mrSges(4,2) * t447 - mrSges(5,3) * t437 - pkin(3) * t394 - pkin(9) * t395 - t544 * t387 + t549 * t389 + t575 * t520 + qJ(4) * (-t520 * t498 + t556) + mrSges(5,2) * t439 + mrSges(4,1) * t446 + t585 * t540 - t584 * t521 + t587 * t484 + (-mrSges(5,1) * qJ(4) - t586) * t483;
t530 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t572;
t529 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t573;
t526 = (-mrSges(3,1) * t550 + mrSges(3,2) * t546) * qJD(1);
t522 = -t552 * pkin(7) + t564;
t519 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t546 + Ifges(3,4) * t550) * qJD(1);
t518 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t546 + Ifges(3,2) * t550) * qJD(1);
t504 = -t550 * g(3) - t577;
t392 = -t483 * mrSges(5,2) - t520 * t508 - t562;
t385 = mrSges(5,1) * t439 + mrSges(4,2) * t485 - mrSges(4,3) * t446 - mrSges(5,3) * t435 + pkin(4) * t395 - qJ(4) * t392 - t579 * t483 + t588 * t484 + t576 * t520 + t587 * t540 + t584 * t541 + t555;
t384 = -mrSges(4,1) * t485 - mrSges(5,1) * t437 + mrSges(5,2) * t435 + mrSges(4,3) * t447 - pkin(3) * t392 - pkin(4) * t557 - pkin(9) * t567 - t549 * t387 - t544 * t389 - t590 * t483 + t579 * t484 + t576 * t521 + t586 * t540 + t575 * t541;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t569 - mrSges(2,2) * t565 + t546 * (mrSges(3,2) * t522 - mrSges(3,3) * t504 + Ifges(3,1) * t527 + Ifges(3,4) * t528 + Ifges(3,5) * qJDD(2) - pkin(8) * t386 - qJD(2) * t518 - t545 * t384 + t582 * t385) + t550 * (-mrSges(3,1) * t522 + mrSges(3,3) * t505 + Ifges(3,4) * t527 + Ifges(3,2) * t528 + Ifges(3,6) * qJDD(2) - pkin(2) * t583 + pkin(8) * t568 + qJD(2) * t519 + t582 * t384 + t545 * t385) + pkin(1) * (-m(3) * t522 + t528 * mrSges(3,1) - t527 * mrSges(3,2) + (-t529 * t546 + t530 * t550) * qJD(1) - t583) + pkin(7) * (t550 * (m(3) * t505 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t528 - qJD(2) * t529 + t526 * t572 + t568) - t546 * (m(3) * t504 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t527 + qJD(2) * t530 - t526 * t573 + t386)); Ifges(3,5) * t527 + Ifges(3,6) * t528 + mrSges(3,1) * t504 - mrSges(3,2) * t505 + t553 + (t546 * t518 - t550 * t519) * qJD(1) + Ifges(3,3) * qJDD(2) + pkin(2) * t386; t553; t394; t555; t559;];
tauJ  = t1;
