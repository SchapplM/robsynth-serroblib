% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-06 13:20
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:15:38
% EndTime: 2019-05-06 13:15:47
% DurationCPUTime: 8.76s
% Computational Cost: add. (127165->349), mult. (295817->446), div. (0->0), fcn. (216391->12), ass. (0->135)
t530 = sin(qJ(2));
t534 = cos(qJ(2));
t550 = qJD(1) * qJD(2);
t516 = qJDD(1) * t530 + t534 * t550;
t537 = qJD(1) ^ 2;
t531 = sin(qJ(1));
t535 = cos(qJ(1));
t543 = -g(1) * t535 - g(2) * t531;
t513 = -pkin(1) * t537 + qJDD(1) * pkin(7) + t543;
t554 = t530 * t513;
t556 = pkin(2) * t537;
t476 = qJDD(2) * pkin(2) - t516 * qJ(3) - t554 + (qJ(3) * t550 + t530 * t556 - g(3)) * t534;
t498 = -g(3) * t530 + t534 * t513;
t517 = qJDD(1) * t534 - t530 * t550;
t553 = qJD(1) * t530;
t518 = qJD(2) * pkin(2) - qJ(3) * t553;
t523 = t534 ^ 2;
t477 = qJ(3) * t517 - qJD(2) * t518 - t523 * t556 + t498;
t525 = sin(pkin(10));
t527 = cos(pkin(10));
t507 = (t525 * t534 + t527 * t530) * qJD(1);
t452 = -0.2e1 * qJD(3) * t507 + t476 * t527 - t525 * t477;
t552 = qJD(1) * t534;
t506 = -t525 * t553 + t527 * t552;
t453 = 0.2e1 * qJD(3) * t506 + t525 * t476 + t527 * t477;
t489 = -pkin(3) * t506 - pkin(8) * t507;
t536 = qJD(2) ^ 2;
t439 = -pkin(3) * t536 + qJDD(2) * pkin(8) + t489 * t506 + t453;
t549 = t531 * g(1) - t535 * g(2);
t541 = -qJDD(1) * pkin(1) - t549;
t480 = -t517 * pkin(2) + qJDD(3) + t518 * t553 + (-qJ(3) * t523 - pkin(7)) * t537 + t541;
t492 = -t525 * t516 + t517 * t527;
t493 = t516 * t527 + t517 * t525;
t442 = (-qJD(2) * t506 - t493) * pkin(8) + (qJD(2) * t507 - t492) * pkin(3) + t480;
t529 = sin(qJ(4));
t533 = cos(qJ(4));
t429 = -t529 * t439 + t533 * t442;
t495 = qJD(2) * t533 - t507 * t529;
t467 = qJD(4) * t495 + qJDD(2) * t529 + t493 * t533;
t491 = qJDD(4) - t492;
t496 = qJD(2) * t529 + t507 * t533;
t505 = qJD(4) - t506;
t418 = (t495 * t505 - t467) * qJ(5) + (t495 * t496 + t491) * pkin(4) + t429;
t430 = t533 * t439 + t529 * t442;
t466 = -qJD(4) * t496 + qJDD(2) * t533 - t493 * t529;
t482 = pkin(4) * t505 - qJ(5) * t496;
t494 = t495 ^ 2;
t420 = -pkin(4) * t494 + qJ(5) * t466 - t482 * t505 + t430;
t524 = sin(pkin(11));
t526 = cos(pkin(11));
t475 = t495 * t524 + t496 * t526;
t412 = -0.2e1 * qJD(5) * t475 + t526 * t418 - t524 * t420;
t447 = t466 * t524 + t467 * t526;
t474 = t495 * t526 - t496 * t524;
t410 = (t474 * t505 - t447) * pkin(9) + (t474 * t475 + t491) * pkin(5) + t412;
t413 = 0.2e1 * qJD(5) * t474 + t524 * t418 + t526 * t420;
t446 = t466 * t526 - t467 * t524;
t459 = pkin(5) * t505 - pkin(9) * t475;
t471 = t474 ^ 2;
t411 = -pkin(5) * t471 + pkin(9) * t446 - t459 * t505 + t413;
t528 = sin(qJ(6));
t532 = cos(qJ(6));
t408 = t410 * t532 - t411 * t528;
t454 = t474 * t532 - t475 * t528;
t426 = qJD(6) * t454 + t446 * t528 + t447 * t532;
t455 = t474 * t528 + t475 * t532;
t436 = -mrSges(7,1) * t454 + mrSges(7,2) * t455;
t501 = qJD(6) + t505;
t443 = -mrSges(7,2) * t501 + mrSges(7,3) * t454;
t490 = qJDD(6) + t491;
t402 = m(7) * t408 + mrSges(7,1) * t490 - mrSges(7,3) * t426 - t436 * t455 + t443 * t501;
t409 = t410 * t528 + t411 * t532;
t425 = -qJD(6) * t455 + t446 * t532 - t447 * t528;
t444 = mrSges(7,1) * t501 - mrSges(7,3) * t455;
t403 = m(7) * t409 - mrSges(7,2) * t490 + mrSges(7,3) * t425 + t436 * t454 - t444 * t501;
t396 = t532 * t402 + t528 * t403;
t456 = -mrSges(6,1) * t474 + mrSges(6,2) * t475;
t457 = -mrSges(6,2) * t505 + mrSges(6,3) * t474;
t394 = m(6) * t412 + mrSges(6,1) * t491 - mrSges(6,3) * t447 - t456 * t475 + t457 * t505 + t396;
t458 = mrSges(6,1) * t505 - mrSges(6,3) * t475;
t545 = -t402 * t528 + t532 * t403;
t395 = m(6) * t413 - mrSges(6,2) * t491 + mrSges(6,3) * t446 + t456 * t474 - t458 * t505 + t545;
t390 = t526 * t394 + t524 * t395;
t449 = Ifges(6,4) * t475 + Ifges(6,2) * t474 + Ifges(6,6) * t505;
t450 = Ifges(6,1) * t475 + Ifges(6,4) * t474 + Ifges(6,5) * t505;
t461 = Ifges(5,4) * t496 + Ifges(5,2) * t495 + Ifges(5,6) * t505;
t462 = Ifges(5,1) * t496 + Ifges(5,4) * t495 + Ifges(5,5) * t505;
t432 = Ifges(7,4) * t455 + Ifges(7,2) * t454 + Ifges(7,6) * t501;
t433 = Ifges(7,1) * t455 + Ifges(7,4) * t454 + Ifges(7,5) * t501;
t540 = -mrSges(7,1) * t408 + mrSges(7,2) * t409 - Ifges(7,5) * t426 - Ifges(7,6) * t425 - Ifges(7,3) * t490 - t455 * t432 + t454 * t433;
t557 = mrSges(5,1) * t429 + mrSges(6,1) * t412 - mrSges(5,2) * t430 - mrSges(6,2) * t413 + Ifges(5,5) * t467 + Ifges(6,5) * t447 + Ifges(5,6) * t466 + Ifges(6,6) * t446 + pkin(4) * t390 + pkin(5) * t396 + t475 * t449 - t474 * t450 + t496 * t461 - t495 * t462 + (Ifges(5,3) + Ifges(6,3)) * t491 - t540;
t487 = -mrSges(4,1) * t506 + mrSges(4,2) * t507;
t500 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t507;
t478 = -mrSges(5,1) * t495 + mrSges(5,2) * t496;
t481 = -mrSges(5,2) * t505 + mrSges(5,3) * t495;
t388 = m(5) * t429 + mrSges(5,1) * t491 - mrSges(5,3) * t467 - t478 * t496 + t481 * t505 + t390;
t483 = mrSges(5,1) * t505 - mrSges(5,3) * t496;
t546 = -t394 * t524 + t526 * t395;
t389 = m(5) * t430 - mrSges(5,2) * t491 + mrSges(5,3) * t466 + t478 * t495 - t483 * t505 + t546;
t547 = -t388 * t529 + t533 * t389;
t380 = m(4) * t453 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t492 - qJD(2) * t500 + t487 * t506 + t547;
t499 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t506;
t438 = -qJDD(2) * pkin(3) - pkin(8) * t536 + t507 * t489 - t452;
t428 = -pkin(4) * t466 - qJ(5) * t494 + t496 * t482 + qJDD(5) + t438;
t415 = -pkin(5) * t446 - pkin(9) * t471 + t459 * t475 + t428;
t542 = m(7) * t415 - t425 * mrSges(7,1) + t426 * mrSges(7,2) - t454 * t443 + t455 * t444;
t406 = m(6) * t428 - t446 * mrSges(6,1) + mrSges(6,2) * t447 - t474 * t457 + t458 * t475 + t542;
t539 = -m(5) * t438 + t466 * mrSges(5,1) - mrSges(5,2) * t467 + t495 * t481 - t483 * t496 - t406;
t405 = m(4) * t452 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t493 + qJD(2) * t499 - t487 * t507 + t539;
t377 = t525 * t380 + t527 * t405;
t382 = t533 * t388 + t529 * t389;
t548 = t527 * t380 - t405 * t525;
t381 = m(4) * t480 - t492 * mrSges(4,1) + t493 * mrSges(4,2) - t506 * t499 + t507 * t500 + t382;
t520 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t552;
t519 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t553;
t515 = (-mrSges(3,1) * t534 + mrSges(3,2) * t530) * qJD(1);
t512 = -t537 * pkin(7) + t541;
t510 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t530 + Ifges(3,4) * t534) * qJD(1);
t509 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t530 + Ifges(3,2) * t534) * qJD(1);
t497 = -t534 * g(3) - t554;
t486 = Ifges(4,1) * t507 + Ifges(4,4) * t506 + Ifges(4,5) * qJD(2);
t485 = Ifges(4,4) * t507 + Ifges(4,2) * t506 + Ifges(4,6) * qJD(2);
t484 = Ifges(4,5) * t507 + Ifges(4,6) * t506 + Ifges(4,3) * qJD(2);
t460 = Ifges(5,5) * t496 + Ifges(5,6) * t495 + Ifges(5,3) * t505;
t448 = Ifges(6,5) * t475 + Ifges(6,6) * t474 + Ifges(6,3) * t505;
t431 = Ifges(7,5) * t455 + Ifges(7,6) * t454 + Ifges(7,3) * t501;
t398 = mrSges(7,2) * t415 - mrSges(7,3) * t408 + Ifges(7,1) * t426 + Ifges(7,4) * t425 + Ifges(7,5) * t490 + t431 * t454 - t432 * t501;
t397 = -mrSges(7,1) * t415 + mrSges(7,3) * t409 + Ifges(7,4) * t426 + Ifges(7,2) * t425 + Ifges(7,6) * t490 - t431 * t455 + t433 * t501;
t384 = mrSges(6,2) * t428 - mrSges(6,3) * t412 + Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t491 - pkin(9) * t396 - t397 * t528 + t398 * t532 + t448 * t474 - t449 * t505;
t383 = -mrSges(6,1) * t428 + mrSges(6,3) * t413 + Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t491 - pkin(5) * t542 + pkin(9) * t545 + t532 * t397 + t528 * t398 - t475 * t448 + t505 * t450;
t376 = mrSges(5,2) * t438 - mrSges(5,3) * t429 + Ifges(5,1) * t467 + Ifges(5,4) * t466 + Ifges(5,5) * t491 - qJ(5) * t390 - t383 * t524 + t384 * t526 + t460 * t495 - t461 * t505;
t375 = -mrSges(5,1) * t438 + mrSges(5,3) * t430 + Ifges(5,4) * t467 + Ifges(5,2) * t466 + Ifges(5,6) * t491 - pkin(4) * t406 + qJ(5) * t546 + t526 * t383 + t524 * t384 - t496 * t460 + t505 * t462;
t374 = -mrSges(4,1) * t480 + mrSges(4,3) * t453 + Ifges(4,4) * t493 + Ifges(4,2) * t492 + Ifges(4,6) * qJDD(2) - pkin(3) * t382 + qJD(2) * t486 - t507 * t484 - t557;
t373 = mrSges(4,2) * t480 - mrSges(4,3) * t452 + Ifges(4,1) * t493 + Ifges(4,4) * t492 + Ifges(4,5) * qJDD(2) - pkin(8) * t382 - qJD(2) * t485 - t375 * t529 + t376 * t533 + t484 * t506;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t549 - mrSges(2,2) * t543 + t530 * (mrSges(3,2) * t512 - mrSges(3,3) * t497 + Ifges(3,1) * t516 + Ifges(3,4) * t517 + Ifges(3,5) * qJDD(2) - qJ(3) * t377 - qJD(2) * t509 + t527 * t373 - t525 * t374) + t534 * (-mrSges(3,1) * t512 + mrSges(3,3) * t498 + Ifges(3,4) * t516 + Ifges(3,2) * t517 + Ifges(3,6) * qJDD(2) - pkin(2) * t381 + qJ(3) * t548 + qJD(2) * t510 + t525 * t373 + t527 * t374) + pkin(1) * (-m(3) * t512 + t517 * mrSges(3,1) - t516 * mrSges(3,2) + (-t519 * t530 + t520 * t534) * qJD(1) - t381) + pkin(7) * (t534 * (m(3) * t498 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t517 - qJD(2) * t519 + t515 * t552 + t548) - t530 * (m(3) * t497 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t516 + qJD(2) * t520 - t515 * t553 + t377)); Ifges(3,5) * t516 + Ifges(3,6) * t517 + mrSges(3,1) * t497 - mrSges(3,2) * t498 + Ifges(4,5) * t493 + Ifges(4,6) * t492 + t507 * t485 - t506 * t486 + mrSges(4,1) * t452 - mrSges(4,2) * t453 + t529 * t376 + t533 * t375 + pkin(3) * t539 + pkin(8) * t547 + pkin(2) * t377 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t530 * t509 - t534 * t510) * qJD(1); t381; t557; t406; -t540;];
tauJ  = t1;
