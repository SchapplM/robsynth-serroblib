% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR8
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
% Datum: 2019-05-07 12:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:16:33
% EndTime: 2019-05-07 12:16:55
% DurationCPUTime: 16.22s
% Computational Cost: add. (253063->362), mult. (557292->472), div. (0->0), fcn. (453628->14), ass. (0->152)
t549 = sin(pkin(6));
t555 = sin(qJ(2));
t560 = cos(qJ(2));
t577 = qJD(1) * qJD(2);
t540 = (-qJDD(1) * t560 + t555 * t577) * t549;
t580 = qJD(1) * t549;
t538 = (-pkin(2) * t560 - pkin(9) * t555) * t580;
t551 = cos(pkin(6));
t545 = qJD(1) * t551 + qJD(2);
t543 = t545 ^ 2;
t544 = qJDD(1) * t551 + qJDD(2);
t579 = qJD(1) * t560;
t562 = qJD(1) ^ 2;
t556 = sin(qJ(1));
t561 = cos(qJ(1));
t569 = -g(1) * t561 - g(2) * t556;
t587 = pkin(8) * t549;
t536 = -pkin(1) * t562 + qJDD(1) * t587 + t569;
t574 = t556 * g(1) - g(2) * t561;
t535 = qJDD(1) * pkin(1) + t562 * t587 + t574;
t584 = t535 * t551;
t581 = t560 * t536 + t555 * t584;
t492 = -t543 * pkin(2) + t544 * pkin(9) + (-g(3) * t555 + t538 * t579) * t549 + t581;
t539 = (qJDD(1) * t555 + t560 * t577) * t549;
t586 = t551 * g(3);
t493 = t540 * pkin(2) - t539 * pkin(9) - t586 + (-t535 + (pkin(2) * t555 - pkin(9) * t560) * t545 * qJD(1)) * t549;
t554 = sin(qJ(3));
t559 = cos(qJ(3));
t463 = -t554 * t492 + t559 * t493;
t576 = t555 * t580;
t528 = t545 * t559 - t554 * t576;
t507 = qJD(3) * t528 + t539 * t559 + t544 * t554;
t529 = t545 * t554 + t559 * t576;
t532 = qJDD(3) + t540;
t575 = t549 * t579;
t542 = qJD(3) - t575;
t448 = (t528 * t542 - t507) * qJ(4) + (t528 * t529 + t532) * pkin(3) + t463;
t464 = t559 * t492 + t554 * t493;
t506 = -qJD(3) * t529 - t539 * t554 + t544 * t559;
t518 = pkin(3) * t542 - qJ(4) * t529;
t527 = t528 ^ 2;
t455 = -pkin(3) * t527 + qJ(4) * t506 - t518 * t542 + t464;
t548 = sin(pkin(12));
t550 = cos(pkin(12));
t515 = t528 * t548 + t529 * t550;
t436 = -0.2e1 * qJD(4) * t515 + t448 * t550 - t548 * t455;
t514 = t528 * t550 - t548 * t529;
t437 = 0.2e1 * qJD(4) * t514 + t548 * t448 + t550 * t455;
t487 = -pkin(4) * t514 - pkin(10) * t515;
t541 = t542 ^ 2;
t435 = -pkin(4) * t541 + pkin(10) * t532 + t487 * t514 + t437;
t582 = t549 * t560;
t509 = -g(3) * t582 - t555 * t536 + t560 * t584;
t491 = -t544 * pkin(2) - t543 * pkin(9) + t538 * t576 - t509;
t457 = -t506 * pkin(3) - t527 * qJ(4) + t529 * t518 + qJDD(4) + t491;
t480 = t506 * t550 - t548 * t507;
t481 = t506 * t548 + t507 * t550;
t440 = (-t514 * t542 - t481) * pkin(10) + (t515 * t542 - t480) * pkin(4) + t457;
t553 = sin(qJ(5));
t558 = cos(qJ(5));
t430 = -t553 * t435 + t558 * t440;
t495 = -t515 * t553 + t542 * t558;
t460 = qJD(5) * t495 + t481 * t558 + t532 * t553;
t479 = qJDD(5) - t480;
t496 = t515 * t558 + t542 * t553;
t513 = qJD(5) - t514;
t428 = (t495 * t513 - t460) * pkin(11) + (t495 * t496 + t479) * pkin(5) + t430;
t431 = t558 * t435 + t553 * t440;
t459 = -qJD(5) * t496 - t481 * t553 + t532 * t558;
t476 = pkin(5) * t513 - pkin(11) * t496;
t494 = t495 ^ 2;
t429 = -pkin(5) * t494 + pkin(11) * t459 - t476 * t513 + t431;
t552 = sin(qJ(6));
t557 = cos(qJ(6));
t426 = t428 * t557 - t429 * t552;
t470 = t495 * t557 - t496 * t552;
t445 = qJD(6) * t470 + t459 * t552 + t460 * t557;
t471 = t495 * t552 + t496 * t557;
t456 = -mrSges(7,1) * t470 + mrSges(7,2) * t471;
t508 = qJD(6) + t513;
t461 = -mrSges(7,2) * t508 + mrSges(7,3) * t470;
t477 = qJDD(6) + t479;
t422 = m(7) * t426 + mrSges(7,1) * t477 - mrSges(7,3) * t445 - t456 * t471 + t461 * t508;
t427 = t428 * t552 + t429 * t557;
t444 = -qJD(6) * t471 + t459 * t557 - t460 * t552;
t462 = mrSges(7,1) * t508 - mrSges(7,3) * t471;
t423 = m(7) * t427 - mrSges(7,2) * t477 + mrSges(7,3) * t444 + t456 * t470 - t462 * t508;
t414 = t557 * t422 + t552 * t423;
t472 = -mrSges(6,1) * t495 + mrSges(6,2) * t496;
t474 = -mrSges(6,2) * t513 + mrSges(6,3) * t495;
t412 = m(6) * t430 + mrSges(6,1) * t479 - mrSges(6,3) * t460 - t472 * t496 + t474 * t513 + t414;
t475 = mrSges(6,1) * t513 - mrSges(6,3) * t496;
t571 = -t422 * t552 + t557 * t423;
t413 = m(6) * t431 - mrSges(6,2) * t479 + mrSges(6,3) * t459 + t472 * t495 - t475 * t513 + t571;
t408 = -t412 * t553 + t558 * t413;
t486 = -mrSges(5,1) * t514 + mrSges(5,2) * t515;
t498 = mrSges(5,1) * t542 - mrSges(5,3) * t515;
t405 = m(5) * t437 - mrSges(5,2) * t532 + mrSges(5,3) * t480 + t486 * t514 - t498 * t542 + t408;
t434 = -pkin(4) * t532 - pkin(10) * t541 + t515 * t487 - t436;
t432 = -pkin(5) * t459 - pkin(11) * t494 + t476 * t496 + t434;
t567 = m(7) * t432 - t444 * mrSges(7,1) + mrSges(7,2) * t445 - t470 * t461 + t462 * t471;
t424 = -m(6) * t434 + t459 * mrSges(6,1) - mrSges(6,2) * t460 + t495 * t474 - t475 * t496 - t567;
t497 = -mrSges(5,2) * t542 + mrSges(5,3) * t514;
t418 = m(5) * t436 + mrSges(5,1) * t532 - mrSges(5,3) * t481 - t486 * t515 + t497 * t542 + t424;
t399 = t548 * t405 + t550 * t418;
t449 = Ifges(7,5) * t471 + Ifges(7,6) * t470 + Ifges(7,3) * t508;
t451 = Ifges(7,1) * t471 + Ifges(7,4) * t470 + Ifges(7,5) * t508;
t415 = -mrSges(7,1) * t432 + mrSges(7,3) * t427 + Ifges(7,4) * t445 + Ifges(7,2) * t444 + Ifges(7,6) * t477 - t449 * t471 + t451 * t508;
t450 = Ifges(7,4) * t471 + Ifges(7,2) * t470 + Ifges(7,6) * t508;
t416 = mrSges(7,2) * t432 - mrSges(7,3) * t426 + Ifges(7,1) * t445 + Ifges(7,4) * t444 + Ifges(7,5) * t477 + t449 * t470 - t450 * t508;
t465 = Ifges(6,5) * t496 + Ifges(6,6) * t495 + Ifges(6,3) * t513;
t467 = Ifges(6,1) * t496 + Ifges(6,4) * t495 + Ifges(6,5) * t513;
t400 = -mrSges(6,1) * t434 + mrSges(6,3) * t431 + Ifges(6,4) * t460 + Ifges(6,2) * t459 + Ifges(6,6) * t479 - pkin(5) * t567 + pkin(11) * t571 + t557 * t415 + t552 * t416 - t496 * t465 + t513 * t467;
t466 = Ifges(6,4) * t496 + Ifges(6,2) * t495 + Ifges(6,6) * t513;
t401 = mrSges(6,2) * t434 - mrSges(6,3) * t430 + Ifges(6,1) * t460 + Ifges(6,4) * t459 + Ifges(6,5) * t479 - pkin(11) * t414 - t415 * t552 + t416 * t557 + t465 * t495 - t466 * t513;
t483 = Ifges(5,4) * t515 + Ifges(5,2) * t514 + Ifges(5,6) * t542;
t484 = Ifges(5,1) * t515 + Ifges(5,4) * t514 + Ifges(5,5) * t542;
t501 = Ifges(4,4) * t529 + Ifges(4,2) * t528 + Ifges(4,6) * t542;
t502 = Ifges(4,1) * t529 + Ifges(4,4) * t528 + Ifges(4,5) * t542;
t588 = Ifges(4,5) * t507 + Ifges(4,6) * t506 + t529 * t501 - t528 * t502 + mrSges(4,1) * t463 - mrSges(4,2) * t464 + Ifges(5,5) * t481 + Ifges(5,6) * t480 + t515 * t483 - t514 * t484 + mrSges(5,1) * t436 - mrSges(5,2) * t437 + t553 * t401 + t558 * t400 + pkin(4) * t424 + pkin(10) * t408 + pkin(3) * t399 + (Ifges(4,3) + Ifges(5,3)) * t532;
t583 = t549 * t555;
t516 = -mrSges(4,1) * t528 + mrSges(4,2) * t529;
t517 = -mrSges(4,2) * t542 + mrSges(4,3) * t528;
t397 = m(4) * t463 + mrSges(4,1) * t532 - mrSges(4,3) * t507 - t516 * t529 + t517 * t542 + t399;
t519 = mrSges(4,1) * t542 - mrSges(4,3) * t529;
t572 = t550 * t405 - t418 * t548;
t398 = m(4) * t464 - mrSges(4,2) * t532 + mrSges(4,3) * t506 + t516 * t528 - t519 * t542 + t572;
t392 = t559 * t397 + t554 * t398;
t407 = t558 * t412 + t553 * t413;
t573 = -t397 * t554 + t559 * t398;
t566 = -mrSges(7,1) * t426 + mrSges(7,2) * t427 - Ifges(7,5) * t445 - Ifges(7,6) * t444 - Ifges(7,3) * t477 - t471 * t450 + t470 * t451;
t406 = m(5) * t457 - t480 * mrSges(5,1) + mrSges(5,2) * t481 - t514 * t497 + t498 * t515 + t407;
t565 = -m(4) * t491 + t506 * mrSges(4,1) - mrSges(4,2) * t507 + t528 * t517 - t519 * t529 - t406;
t564 = mrSges(6,1) * t430 - mrSges(6,2) * t431 + Ifges(6,5) * t460 + Ifges(6,6) * t459 + Ifges(6,3) * t479 + pkin(5) * t414 + t496 * t466 - t495 * t467 - t566;
t537 = (-mrSges(3,1) * t560 + mrSges(3,2) * t555) * t580;
t534 = -mrSges(3,2) * t545 + mrSges(3,3) * t575;
t533 = mrSges(3,1) * t545 - mrSges(3,3) * t576;
t523 = -t549 * t535 - t586;
t522 = Ifges(3,5) * t545 + (Ifges(3,1) * t555 + Ifges(3,4) * t560) * t580;
t521 = Ifges(3,6) * t545 + (Ifges(3,4) * t555 + Ifges(3,2) * t560) * t580;
t520 = Ifges(3,3) * t545 + (Ifges(3,5) * t555 + Ifges(3,6) * t560) * t580;
t510 = -g(3) * t583 + t581;
t500 = Ifges(4,5) * t529 + Ifges(4,6) * t528 + Ifges(4,3) * t542;
t482 = Ifges(5,5) * t515 + Ifges(5,6) * t514 + Ifges(5,3) * t542;
t402 = m(3) * t509 + mrSges(3,1) * t544 - mrSges(3,3) * t539 + t534 * t545 - t537 * t576 + t565;
t393 = -mrSges(5,1) * t457 + mrSges(5,3) * t437 + Ifges(5,4) * t481 + Ifges(5,2) * t480 + Ifges(5,6) * t532 - pkin(4) * t407 - t515 * t482 + t542 * t484 - t564;
t391 = m(3) * t510 - mrSges(3,2) * t544 - mrSges(3,3) * t540 - t533 * t545 + t537 * t575 + t573;
t390 = mrSges(5,2) * t457 - mrSges(5,3) * t436 + Ifges(5,1) * t481 + Ifges(5,4) * t480 + Ifges(5,5) * t532 - pkin(10) * t407 - t400 * t553 + t401 * t558 + t482 * t514 - t483 * t542;
t389 = mrSges(4,2) * t491 - mrSges(4,3) * t463 + Ifges(4,1) * t507 + Ifges(4,4) * t506 + Ifges(4,5) * t532 - qJ(4) * t399 + t390 * t550 - t393 * t548 + t500 * t528 - t501 * t542;
t388 = -mrSges(4,1) * t491 + mrSges(4,3) * t464 + Ifges(4,4) * t507 + Ifges(4,2) * t506 + Ifges(4,6) * t532 - pkin(3) * t406 + qJ(4) * t572 + t548 * t390 + t550 * t393 - t529 * t500 + t542 * t502;
t387 = Ifges(3,5) * t539 - Ifges(3,6) * t540 + Ifges(3,3) * t544 + mrSges(3,1) * t509 - mrSges(3,2) * t510 + t554 * t389 + t559 * t388 + pkin(2) * t565 + pkin(9) * t573 + (t521 * t555 - t522 * t560) * t580;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t574 - mrSges(2,2) * t569 + (mrSges(3,2) * t523 - mrSges(3,3) * t509 + Ifges(3,1) * t539 - Ifges(3,4) * t540 + Ifges(3,5) * t544 - pkin(9) * t392 - t388 * t554 + t389 * t559 + t520 * t575 - t521 * t545) * t583 + (-mrSges(3,1) * t523 + mrSges(3,3) * t510 + Ifges(3,4) * t539 - Ifges(3,2) * t540 + Ifges(3,6) * t544 - pkin(2) * t392 - t520 * t576 + t545 * t522 - t588) * t582 + t551 * t387 + pkin(1) * ((t391 * t555 + t402 * t560) * t551 + (-m(3) * t523 - t540 * mrSges(3,1) - t539 * mrSges(3,2) + (-t533 * t555 + t534 * t560) * t580 - t392) * t549) + (t391 * t560 - t402 * t555) * t587; t387; t588; t406; t564; -t566;];
tauJ  = t1;
