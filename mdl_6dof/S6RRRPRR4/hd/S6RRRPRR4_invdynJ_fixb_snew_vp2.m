% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 10:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 10:26:02
% EndTime: 2019-05-07 10:26:13
% DurationCPUTime: 10.40s
% Computational Cost: add. (168875->350), mult. (348161->445), div. (0->0), fcn. (257346->12), ass. (0->138)
t559 = cos(qJ(3));
t540 = qJD(1) ^ 2;
t558 = pkin(2) * t540;
t535 = sin(qJ(1));
t539 = cos(qJ(1));
t548 = -g(1) * t539 - g(2) * t535;
t513 = -pkin(1) * t540 + qJDD(1) * pkin(7) + t548;
t534 = sin(qJ(2));
t557 = t534 * t513;
t538 = cos(qJ(2));
t554 = qJD(1) * qJD(2);
t517 = qJDD(1) * t534 + t538 * t554;
t476 = qJDD(2) * pkin(2) - t517 * pkin(8) - t557 + (pkin(8) * t554 + t534 * t558 - g(3)) * t538;
t501 = -g(3) * t534 + t538 * t513;
t518 = qJDD(1) * t538 - t534 * t554;
t556 = qJD(1) * t534;
t521 = qJD(2) * pkin(2) - pkin(8) * t556;
t528 = t538 ^ 2;
t477 = pkin(8) * t518 - qJD(2) * t521 - t528 * t558 + t501;
t533 = sin(qJ(3));
t456 = t533 * t476 + t477 * t559;
t511 = (t533 * t538 + t534 * t559) * qJD(1);
t484 = qJD(3) * t511 + t517 * t533 - t518 * t559;
t555 = qJD(1) * t538;
t510 = t533 * t556 - t555 * t559;
t494 = mrSges(4,1) * t510 + mrSges(4,2) * t511;
t527 = qJD(2) + qJD(3);
t503 = mrSges(4,1) * t527 - mrSges(4,3) * t511;
t526 = qJDD(2) + qJDD(3);
t485 = -t510 * qJD(3) + t517 * t559 + t533 * t518;
t553 = t535 * g(1) - t539 * g(2);
t546 = -qJDD(1) * pkin(1) - t553;
t486 = -t518 * pkin(2) + t521 * t556 + (-pkin(8) * t528 - pkin(7)) * t540 + t546;
t441 = (t510 * t527 - t485) * qJ(4) + (t511 * t527 + t484) * pkin(3) + t486;
t493 = pkin(3) * t510 - qJ(4) * t511;
t525 = t527 ^ 2;
t444 = -pkin(3) * t525 + qJ(4) * t526 - t493 * t510 + t456;
t529 = sin(pkin(11));
t530 = cos(pkin(11));
t499 = t511 * t530 + t527 * t529;
t427 = -0.2e1 * qJD(4) * t499 + t530 * t441 - t529 * t444;
t470 = t485 * t530 + t526 * t529;
t498 = -t511 * t529 + t527 * t530;
t417 = (t498 * t510 - t470) * pkin(9) + (t498 * t499 + t484) * pkin(4) + t427;
t428 = 0.2e1 * qJD(4) * t498 + t529 * t441 + t530 * t444;
t469 = -t485 * t529 + t526 * t530;
t489 = pkin(4) * t510 - pkin(9) * t499;
t497 = t498 ^ 2;
t419 = -pkin(4) * t497 + pkin(9) * t469 - t489 * t510 + t428;
t532 = sin(qJ(5));
t537 = cos(qJ(5));
t411 = t537 * t417 - t532 * t419;
t467 = t498 * t537 - t499 * t532;
t440 = qJD(5) * t467 + t469 * t532 + t470 * t537;
t468 = t498 * t532 + t499 * t537;
t483 = qJDD(5) + t484;
t506 = qJD(5) + t510;
t409 = (t467 * t506 - t440) * pkin(10) + (t467 * t468 + t483) * pkin(5) + t411;
t412 = t532 * t417 + t537 * t419;
t439 = -qJD(5) * t468 + t469 * t537 - t470 * t532;
t459 = pkin(5) * t506 - pkin(10) * t468;
t466 = t467 ^ 2;
t410 = -pkin(5) * t466 + pkin(10) * t439 - t459 * t506 + t412;
t531 = sin(qJ(6));
t536 = cos(qJ(6));
t407 = t409 * t536 - t410 * t531;
t451 = t467 * t536 - t468 * t531;
t425 = qJD(6) * t451 + t439 * t531 + t440 * t536;
t452 = t467 * t531 + t468 * t536;
t435 = -mrSges(7,1) * t451 + mrSges(7,2) * t452;
t505 = qJD(6) + t506;
t445 = -mrSges(7,2) * t505 + mrSges(7,3) * t451;
t479 = qJDD(6) + t483;
t400 = m(7) * t407 + mrSges(7,1) * t479 - mrSges(7,3) * t425 - t435 * t452 + t445 * t505;
t408 = t409 * t531 + t410 * t536;
t424 = -qJD(6) * t452 + t439 * t536 - t440 * t531;
t446 = mrSges(7,1) * t505 - mrSges(7,3) * t452;
t401 = m(7) * t408 - mrSges(7,2) * t479 + mrSges(7,3) * t424 + t435 * t451 - t446 * t505;
t394 = t536 * t400 + t531 * t401;
t453 = -mrSges(6,1) * t467 + mrSges(6,2) * t468;
t457 = -mrSges(6,2) * t506 + mrSges(6,3) * t467;
t392 = m(6) * t411 + mrSges(6,1) * t483 - mrSges(6,3) * t440 - t453 * t468 + t457 * t506 + t394;
t458 = mrSges(6,1) * t506 - mrSges(6,3) * t468;
t549 = -t400 * t531 + t536 * t401;
t393 = m(6) * t412 - mrSges(6,2) * t483 + mrSges(6,3) * t439 + t453 * t467 - t458 * t506 + t549;
t388 = t537 * t392 + t532 * t393;
t472 = -mrSges(5,1) * t498 + mrSges(5,2) * t499;
t487 = -mrSges(5,2) * t510 + mrSges(5,3) * t498;
t386 = m(5) * t427 + mrSges(5,1) * t484 - mrSges(5,3) * t470 - t472 * t499 + t487 * t510 + t388;
t488 = mrSges(5,1) * t510 - mrSges(5,3) * t499;
t550 = -t392 * t532 + t537 * t393;
t387 = m(5) * t428 - mrSges(5,2) * t484 + mrSges(5,3) * t469 + t472 * t498 - t488 * t510 + t550;
t551 = -t386 * t529 + t530 * t387;
t378 = m(4) * t456 - mrSges(4,2) * t526 - mrSges(4,3) * t484 - t494 * t510 - t503 * t527 + t551;
t455 = t476 * t559 - t533 * t477;
t443 = -t526 * pkin(3) - t525 * qJ(4) + t511 * t493 + qJDD(4) - t455;
t429 = -t469 * pkin(4) - t497 * pkin(9) + t499 * t489 + t443;
t414 = -t439 * pkin(5) - t466 * pkin(10) + t468 * t459 + t429;
t547 = m(7) * t414 - t424 * mrSges(7,1) + t425 * mrSges(7,2) - t451 * t445 + t452 * t446;
t542 = m(6) * t429 - t439 * mrSges(6,1) + mrSges(6,2) * t440 - t467 * t457 + t458 * t468 + t547;
t405 = m(5) * t443 - t469 * mrSges(5,1) + mrSges(5,2) * t470 - t498 * t487 + t488 * t499 + t542;
t502 = -mrSges(4,2) * t527 - mrSges(4,3) * t510;
t403 = m(4) * t455 + mrSges(4,1) * t526 - mrSges(4,3) * t485 - t494 * t511 + t502 * t527 - t405;
t375 = t533 * t378 + t403 * t559;
t380 = t530 * t386 + t529 * t387;
t552 = t378 * t559 - t403 * t533;
t431 = Ifges(7,4) * t452 + Ifges(7,2) * t451 + Ifges(7,6) * t505;
t432 = Ifges(7,1) * t452 + Ifges(7,4) * t451 + Ifges(7,5) * t505;
t545 = -mrSges(7,1) * t407 + mrSges(7,2) * t408 - Ifges(7,5) * t425 - Ifges(7,6) * t424 - Ifges(7,3) * t479 - t452 * t431 + t451 * t432;
t430 = Ifges(7,5) * t452 + Ifges(7,6) * t451 + Ifges(7,3) * t505;
t395 = -mrSges(7,1) * t414 + mrSges(7,3) * t408 + Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t479 - t430 * t452 + t432 * t505;
t396 = mrSges(7,2) * t414 - mrSges(7,3) * t407 + Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t479 + t430 * t451 - t431 * t505;
t447 = Ifges(6,5) * t468 + Ifges(6,6) * t467 + Ifges(6,3) * t506;
t449 = Ifges(6,1) * t468 + Ifges(6,4) * t467 + Ifges(6,5) * t506;
t381 = -mrSges(6,1) * t429 + mrSges(6,3) * t412 + Ifges(6,4) * t440 + Ifges(6,2) * t439 + Ifges(6,6) * t483 - pkin(5) * t547 + pkin(10) * t549 + t536 * t395 + t531 * t396 - t468 * t447 + t506 * t449;
t448 = Ifges(6,4) * t468 + Ifges(6,2) * t467 + Ifges(6,6) * t506;
t382 = mrSges(6,2) * t429 - mrSges(6,3) * t411 + Ifges(6,1) * t440 + Ifges(6,4) * t439 + Ifges(6,5) * t483 - pkin(10) * t394 - t395 * t531 + t396 * t536 + t447 * t467 - t448 * t506;
t460 = Ifges(5,5) * t499 + Ifges(5,6) * t498 + Ifges(5,3) * t510;
t462 = Ifges(5,1) * t499 + Ifges(5,4) * t498 + Ifges(5,5) * t510;
t372 = -mrSges(5,1) * t443 + mrSges(5,3) * t428 + Ifges(5,4) * t470 + Ifges(5,2) * t469 + Ifges(5,6) * t484 - pkin(4) * t542 + pkin(9) * t550 + t537 * t381 + t532 * t382 - t499 * t460 + t510 * t462;
t461 = Ifges(5,4) * t499 + Ifges(5,2) * t498 + Ifges(5,6) * t510;
t374 = mrSges(5,2) * t443 - mrSges(5,3) * t427 + Ifges(5,1) * t470 + Ifges(5,4) * t469 + Ifges(5,5) * t484 - pkin(9) * t388 - t381 * t532 + t382 * t537 + t460 * t498 - t461 * t510;
t491 = Ifges(4,4) * t511 - Ifges(4,2) * t510 + Ifges(4,6) * t527;
t492 = Ifges(4,1) * t511 - Ifges(4,4) * t510 + Ifges(4,5) * t527;
t544 = mrSges(4,1) * t455 - mrSges(4,2) * t456 + Ifges(4,5) * t485 - Ifges(4,6) * t484 + Ifges(4,3) * t526 - pkin(3) * t405 + qJ(4) * t551 + t530 * t372 + t529 * t374 + t511 * t491 + t510 * t492;
t543 = m(4) * t486 + t484 * mrSges(4,1) + t485 * mrSges(4,2) + t510 * t502 + t511 * t503 + t380;
t541 = mrSges(6,1) * t411 - mrSges(6,2) * t412 + Ifges(6,5) * t440 + Ifges(6,6) * t439 + Ifges(6,3) * t483 + pkin(5) * t394 + t468 * t448 - t467 * t449 - t545;
t520 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t555;
t519 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t556;
t516 = (-mrSges(3,1) * t538 + mrSges(3,2) * t534) * qJD(1);
t512 = -t540 * pkin(7) + t546;
t509 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t534 + Ifges(3,4) * t538) * qJD(1);
t508 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t534 + Ifges(3,2) * t538) * qJD(1);
t500 = -t538 * g(3) - t557;
t490 = Ifges(4,5) * t511 - Ifges(4,6) * t510 + Ifges(4,3) * t527;
t370 = (-Ifges(5,3) - Ifges(4,2)) * t484 + t527 * t492 + Ifges(4,6) * t526 - t511 * t490 + t498 * t462 - t499 * t461 + Ifges(4,4) * t485 - mrSges(4,1) * t486 - Ifges(5,6) * t469 - Ifges(5,5) * t470 + mrSges(4,3) * t456 + mrSges(5,2) * t428 - mrSges(5,1) * t427 - pkin(4) * t388 - pkin(3) * t380 - t541;
t369 = mrSges(4,2) * t486 - mrSges(4,3) * t455 + Ifges(4,1) * t485 - Ifges(4,4) * t484 + Ifges(4,5) * t526 - qJ(4) * t380 - t372 * t529 + t374 * t530 - t490 * t510 - t491 * t527;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t553 - mrSges(2,2) * t548 + t534 * (mrSges(3,2) * t512 - mrSges(3,3) * t500 + Ifges(3,1) * t517 + Ifges(3,4) * t518 + Ifges(3,5) * qJDD(2) - pkin(8) * t375 - qJD(2) * t508 + t559 * t369 - t533 * t370) + t538 * (-mrSges(3,1) * t512 + mrSges(3,3) * t501 + Ifges(3,4) * t517 + Ifges(3,2) * t518 + Ifges(3,6) * qJDD(2) - pkin(2) * t543 + pkin(8) * t552 + qJD(2) * t509 + t533 * t369 + t559 * t370) + pkin(1) * (-m(3) * t512 + t518 * mrSges(3,1) - t517 * mrSges(3,2) + (-t519 * t534 + t520 * t538) * qJD(1) - t543) + pkin(7) * (t538 * (m(3) * t501 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t518 - qJD(2) * t519 + t516 * t555 + t552) - t534 * (m(3) * t500 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t517 + qJD(2) * t520 - t516 * t556 + t375)); pkin(2) * t375 + (t534 * t508 - t538 * t509) * qJD(1) + t544 + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t518 + Ifges(3,5) * t517 + mrSges(3,1) * t500 - mrSges(3,2) * t501; t544; t405; t541; -t545;];
tauJ  = t1;
