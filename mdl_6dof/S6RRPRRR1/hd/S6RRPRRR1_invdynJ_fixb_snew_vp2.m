% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 19:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:30:36
% EndTime: 2019-05-06 19:30:46
% DurationCPUTime: 9.59s
% Computational Cost: add. (136794->349), mult. (329922->446), div. (0->0), fcn. (253207->12), ass. (0->138)
t537 = qJD(1) ^ 2;
t555 = pkin(2) * t537;
t531 = sin(qJ(1));
t536 = cos(qJ(1));
t545 = -g(1) * t536 - g(2) * t531;
t508 = -pkin(1) * t537 + qJDD(1) * pkin(7) + t545;
t530 = sin(qJ(2));
t554 = t508 * t530;
t535 = cos(qJ(2));
t551 = qJD(1) * qJD(2);
t511 = qJDD(1) * t530 + t535 * t551;
t477 = qJDD(2) * pkin(2) - qJ(3) * t511 - t554 + (qJ(3) * t551 + t530 * t555 - g(3)) * t535;
t494 = -g(3) * t530 + t535 * t508;
t512 = qJDD(1) * t535 - t530 * t551;
t553 = qJD(1) * t530;
t513 = qJD(2) * pkin(2) - qJ(3) * t553;
t524 = t535 ^ 2;
t478 = qJ(3) * t512 - qJD(2) * t513 - t524 * t555 + t494;
t525 = sin(pkin(11));
t526 = cos(pkin(11));
t503 = (t535 * t525 + t530 * t526) * qJD(1);
t453 = -0.2e1 * qJD(3) * t503 + t526 * t477 - t478 * t525;
t492 = t511 * t526 + t512 * t525;
t502 = (-t530 * t525 + t535 * t526) * qJD(1);
t439 = (qJD(2) * t502 - t492) * pkin(8) + (t502 * t503 + qJDD(2)) * pkin(3) + t453;
t454 = 0.2e1 * qJD(3) * t502 + t525 * t477 + t526 * t478;
t491 = -t511 * t525 + t512 * t526;
t497 = qJD(2) * pkin(3) - pkin(8) * t503;
t501 = t502 ^ 2;
t441 = -pkin(3) * t501 + pkin(8) * t491 - qJD(2) * t497 + t454;
t529 = sin(qJ(4));
t534 = cos(qJ(4));
t419 = t534 * t439 - t441 * t529;
t487 = t502 * t534 - t503 * t529;
t460 = qJD(4) * t487 + t491 * t529 + t492 * t534;
t488 = t502 * t529 + t503 * t534;
t522 = qJDD(2) + qJDD(4);
t523 = qJD(2) + qJD(4);
t415 = (t487 * t523 - t460) * pkin(9) + (t487 * t488 + t522) * pkin(4) + t419;
t420 = t529 * t439 + t534 * t441;
t459 = -qJD(4) * t488 + t491 * t534 - t492 * t529;
t482 = pkin(4) * t523 - pkin(9) * t488;
t483 = t487 ^ 2;
t417 = -pkin(4) * t483 + pkin(9) * t459 - t482 * t523 + t420;
t528 = sin(qJ(5));
t533 = cos(qJ(5));
t412 = t528 * t415 + t533 * t417;
t472 = t487 * t528 + t488 * t533;
t430 = -qJD(5) * t472 + t459 * t533 - t460 * t528;
t471 = t487 * t533 - t488 * t528;
t449 = -mrSges(6,1) * t471 + mrSges(6,2) * t472;
t520 = qJD(5) + t523;
t464 = mrSges(6,1) * t520 - mrSges(6,3) * t472;
t519 = qJDD(5) + t522;
t450 = -pkin(5) * t471 - pkin(10) * t472;
t518 = t520 ^ 2;
t409 = -pkin(5) * t518 + pkin(10) * t519 + t450 * t471 + t412;
t550 = g(1) * t531 - t536 * g(2);
t544 = -qJDD(1) * pkin(1) - t550;
t479 = -pkin(2) * t512 + qJDD(3) + t513 * t553 + (-qJ(3) * t524 - pkin(7)) * t537 + t544;
t452 = -pkin(3) * t491 - pkin(8) * t501 + t503 * t497 + t479;
t425 = -pkin(4) * t459 - pkin(9) * t483 + t488 * t482 + t452;
t431 = qJD(5) * t471 + t459 * t528 + t460 * t533;
t413 = t425 + (-t471 * t520 - t431) * pkin(10) + (t472 * t520 - t430) * pkin(5);
t527 = sin(qJ(6));
t532 = cos(qJ(6));
t406 = -t409 * t527 + t413 * t532;
t461 = -t472 * t527 + t520 * t532;
t423 = qJD(6) * t461 + t431 * t532 + t519 * t527;
t429 = qJDD(6) - t430;
t462 = t472 * t532 + t520 * t527;
t442 = -mrSges(7,1) * t461 + mrSges(7,2) * t462;
t468 = qJD(6) - t471;
t443 = -mrSges(7,2) * t468 + mrSges(7,3) * t461;
t403 = m(7) * t406 + mrSges(7,1) * t429 - mrSges(7,3) * t423 - t442 * t462 + t443 * t468;
t407 = t409 * t532 + t413 * t527;
t422 = -qJD(6) * t462 - t431 * t527 + t519 * t532;
t444 = mrSges(7,1) * t468 - mrSges(7,3) * t462;
t404 = m(7) * t407 - mrSges(7,2) * t429 + mrSges(7,3) * t422 + t442 * t461 - t444 * t468;
t546 = -t403 * t527 + t532 * t404;
t391 = m(6) * t412 - mrSges(6,2) * t519 + mrSges(6,3) * t430 + t449 * t471 - t464 * t520 + t546;
t411 = t415 * t533 - t417 * t528;
t463 = -mrSges(6,2) * t520 + mrSges(6,3) * t471;
t408 = -pkin(5) * t519 - pkin(10) * t518 + t450 * t472 - t411;
t543 = -m(7) * t408 + t422 * mrSges(7,1) - mrSges(7,2) * t423 + t461 * t443 - t444 * t462;
t399 = m(6) * t411 + mrSges(6,1) * t519 - mrSges(6,3) * t431 - t449 * t472 + t463 * t520 + t543;
t387 = t528 * t391 + t533 * t399;
t473 = -mrSges(5,1) * t487 + mrSges(5,2) * t488;
t480 = -mrSges(5,2) * t523 + mrSges(5,3) * t487;
t384 = m(5) * t419 + mrSges(5,1) * t522 - mrSges(5,3) * t460 - t473 * t488 + t480 * t523 + t387;
t481 = mrSges(5,1) * t523 - mrSges(5,3) * t488;
t547 = t533 * t391 - t399 * t528;
t385 = m(5) * t420 - mrSges(5,2) * t522 + mrSges(5,3) * t459 + t473 * t487 - t481 * t523 + t547;
t378 = t534 * t384 + t529 * t385;
t490 = -mrSges(4,1) * t502 + mrSges(4,2) * t503;
t495 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t502;
t376 = m(4) * t453 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t492 + qJD(2) * t495 - t490 * t503 + t378;
t496 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t503;
t548 = -t529 * t384 + t534 * t385;
t377 = m(4) * t454 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t491 - qJD(2) * t496 + t490 * t502 + t548;
t371 = t526 * t376 + t525 * t377;
t393 = t532 * t403 + t527 * t404;
t552 = qJD(1) * t535;
t549 = -t376 * t525 + t526 * t377;
t542 = -m(6) * t425 + t430 * mrSges(6,1) - t431 * mrSges(6,2) + t471 * t463 - t472 * t464 - t393;
t432 = Ifges(7,5) * t462 + Ifges(7,6) * t461 + Ifges(7,3) * t468;
t434 = Ifges(7,1) * t462 + Ifges(7,4) * t461 + Ifges(7,5) * t468;
t396 = -mrSges(7,1) * t408 + mrSges(7,3) * t407 + Ifges(7,4) * t423 + Ifges(7,2) * t422 + Ifges(7,6) * t429 - t432 * t462 + t434 * t468;
t433 = Ifges(7,4) * t462 + Ifges(7,2) * t461 + Ifges(7,6) * t468;
t397 = mrSges(7,2) * t408 - mrSges(7,3) * t406 + Ifges(7,1) * t423 + Ifges(7,4) * t422 + Ifges(7,5) * t429 + t432 * t461 - t433 * t468;
t446 = Ifges(6,4) * t472 + Ifges(6,2) * t471 + Ifges(6,6) * t520;
t447 = Ifges(6,1) * t472 + Ifges(6,4) * t471 + Ifges(6,5) * t520;
t541 = mrSges(6,1) * t411 - mrSges(6,2) * t412 + Ifges(6,5) * t431 + Ifges(6,6) * t430 + Ifges(6,3) * t519 + pkin(5) * t543 + pkin(10) * t546 + t532 * t396 + t527 * t397 + t472 * t446 - t471 * t447;
t540 = mrSges(7,1) * t406 - mrSges(7,2) * t407 + Ifges(7,5) * t423 + Ifges(7,6) * t422 + Ifges(7,3) * t429 + t433 * t462 - t434 * t461;
t539 = -m(5) * t452 + t459 * mrSges(5,1) - t460 * mrSges(5,2) + t487 * t480 - t488 * t481 + t542;
t466 = Ifges(5,4) * t488 + Ifges(5,2) * t487 + Ifges(5,6) * t523;
t467 = Ifges(5,1) * t488 + Ifges(5,4) * t487 + Ifges(5,5) * t523;
t538 = mrSges(5,1) * t419 - mrSges(5,2) * t420 + Ifges(5,5) * t460 + Ifges(5,6) * t459 + Ifges(5,3) * t522 + pkin(4) * t387 + t488 * t466 - t487 * t467 + t541;
t388 = m(4) * t479 - t491 * mrSges(4,1) + t492 * mrSges(4,2) - t502 * t495 + t503 * t496 - t539;
t515 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t552;
t514 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t553;
t510 = (-t535 * mrSges(3,1) + t530 * mrSges(3,2)) * qJD(1);
t507 = -pkin(7) * t537 + t544;
t506 = Ifges(3,5) * qJD(2) + (t530 * Ifges(3,1) + t535 * Ifges(3,4)) * qJD(1);
t505 = Ifges(3,6) * qJD(2) + (t530 * Ifges(3,4) + t535 * Ifges(3,2)) * qJD(1);
t493 = -g(3) * t535 - t554;
t486 = Ifges(4,1) * t503 + Ifges(4,4) * t502 + Ifges(4,5) * qJD(2);
t485 = Ifges(4,4) * t503 + Ifges(4,2) * t502 + Ifges(4,6) * qJD(2);
t484 = Ifges(4,5) * t503 + Ifges(4,6) * t502 + Ifges(4,3) * qJD(2);
t465 = Ifges(5,5) * t488 + Ifges(5,6) * t487 + Ifges(5,3) * t523;
t445 = Ifges(6,5) * t472 + Ifges(6,6) * t471 + Ifges(6,3) * t520;
t380 = -mrSges(6,1) * t425 + mrSges(6,3) * t412 + Ifges(6,4) * t431 + Ifges(6,2) * t430 + Ifges(6,6) * t519 - pkin(5) * t393 - t445 * t472 + t447 * t520 - t540;
t379 = mrSges(6,2) * t425 - mrSges(6,3) * t411 + Ifges(6,1) * t431 + Ifges(6,4) * t430 + Ifges(6,5) * t519 - pkin(10) * t393 - t396 * t527 + t397 * t532 + t445 * t471 - t446 * t520;
t372 = mrSges(5,2) * t452 - mrSges(5,3) * t419 + Ifges(5,1) * t460 + Ifges(5,4) * t459 + Ifges(5,5) * t522 - pkin(9) * t387 + t379 * t533 - t380 * t528 + t465 * t487 - t466 * t523;
t370 = -mrSges(5,1) * t452 + mrSges(5,3) * t420 + Ifges(5,4) * t460 + Ifges(5,2) * t459 + Ifges(5,6) * t522 + pkin(4) * t542 + pkin(9) * t547 + t528 * t379 + t533 * t380 - t488 * t465 + t523 * t467;
t369 = mrSges(4,2) * t479 - mrSges(4,3) * t453 + Ifges(4,1) * t492 + Ifges(4,4) * t491 + Ifges(4,5) * qJDD(2) - pkin(8) * t378 - qJD(2) * t485 - t370 * t529 + t372 * t534 + t484 * t502;
t368 = -mrSges(4,1) * t479 + mrSges(4,3) * t454 + Ifges(4,4) * t492 + Ifges(4,2) * t491 + Ifges(4,6) * qJDD(2) + pkin(3) * t539 + pkin(8) * t548 + qJD(2) * t486 + t534 * t370 + t529 * t372 - t503 * t484;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t550 - mrSges(2,2) * t545 + t530 * (mrSges(3,2) * t507 - mrSges(3,3) * t493 + Ifges(3,1) * t511 + Ifges(3,4) * t512 + Ifges(3,5) * qJDD(2) - qJ(3) * t371 - qJD(2) * t505 - t368 * t525 + t369 * t526) + t535 * (-mrSges(3,1) * t507 + mrSges(3,3) * t494 + Ifges(3,4) * t511 + Ifges(3,2) * t512 + Ifges(3,6) * qJDD(2) - pkin(2) * t388 + qJ(3) * t549 + qJD(2) * t506 + t526 * t368 + t525 * t369) + pkin(1) * (-m(3) * t507 - t511 * mrSges(3,2) + t512 * mrSges(3,1) + (-t530 * t514 + t535 * t515) * qJD(1) - t388) + pkin(7) * (t535 * (m(3) * t494 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t512 - qJD(2) * t514 + t510 * t552 + t549) - t530 * (m(3) * t493 + qJDD(2) * mrSges(3,1) - t511 * mrSges(3,3) + qJD(2) * t515 - t510 * t553 + t371)); Ifges(3,5) * t511 + Ifges(3,6) * t512 - t502 * t486 + t503 * t485 + Ifges(4,5) * t492 + mrSges(3,1) * t493 - mrSges(3,2) * t494 + Ifges(4,6) * t491 + mrSges(4,1) * t453 - mrSges(4,2) * t454 + pkin(2) * t371 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t530 * t505 - t535 * t506) * qJD(1) + t538 + pkin(3) * t378; t388; t538; t541; t540;];
tauJ  = t1;
