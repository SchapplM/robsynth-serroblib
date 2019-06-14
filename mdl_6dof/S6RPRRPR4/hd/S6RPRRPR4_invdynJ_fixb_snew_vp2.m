% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 22:31
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:27:47
% EndTime: 2019-05-05 22:27:54
% DurationCPUTime: 7.27s
% Computational Cost: add. (107403->318), mult. (261715->406), div. (0->0), fcn. (208062->12), ass. (0->133)
t530 = qJD(1) ^ 2;
t555 = cos(qJ(4));
t554 = pkin(2) * t530;
t553 = pkin(7) * qJDD(1);
t526 = sin(qJ(1));
t529 = cos(qJ(1));
t541 = -g(1) * t529 - g(2) * t526;
t507 = -pkin(1) * t530 + qJDD(1) * qJ(2) + t541;
t520 = sin(pkin(10));
t522 = cos(pkin(10));
t548 = qJD(1) * qJD(2);
t546 = -t522 * g(3) - 0.2e1 * t520 * t548;
t482 = (t522 * t554 - t507 - t553) * t520 + t546;
t497 = -g(3) * t520 + (t507 + 0.2e1 * t548) * t522;
t517 = t522 ^ 2;
t483 = -t517 * t554 + t522 * t553 + t497;
t525 = sin(qJ(3));
t528 = cos(qJ(3));
t461 = t528 * t482 - t525 * t483;
t538 = t520 * t528 + t522 * t525;
t537 = -t520 * t525 + t522 * t528;
t505 = t537 * qJD(1);
t549 = t505 * qJD(3);
t495 = t538 * qJDD(1) + t549;
t506 = t538 * qJD(1);
t435 = (-t495 + t549) * pkin(8) + (t505 * t506 + qJDD(3)) * pkin(3) + t461;
t462 = t525 * t482 + t528 * t483;
t494 = -t506 * qJD(3) + t537 * qJDD(1);
t500 = qJD(3) * pkin(3) - pkin(8) * t506;
t504 = t505 ^ 2;
t438 = -pkin(3) * t504 + pkin(8) * t494 - qJD(3) * t500 + t462;
t524 = sin(qJ(4));
t428 = t524 * t435 + t555 * t438;
t489 = t524 * t505 + t555 * t506;
t456 = qJD(4) * t489 - t555 * t494 + t495 * t524;
t488 = -t555 * t505 + t506 * t524;
t471 = mrSges(5,1) * t488 + mrSges(5,2) * t489;
t518 = qJD(3) + qJD(4);
t480 = mrSges(5,1) * t518 - mrSges(5,3) * t489;
t515 = qJDD(3) + qJDD(4);
t470 = pkin(4) * t488 - qJ(5) * t489;
t514 = t518 ^ 2;
t419 = -pkin(4) * t514 + qJ(5) * t515 - t470 * t488 + t428;
t547 = t526 * g(1) - t529 * g(2);
t540 = qJDD(2) - t547;
t551 = -t520 ^ 2 - t517;
t493 = (-pkin(2) * t522 - pkin(1)) * qJDD(1) + (t551 * pkin(7) - qJ(2)) * t530 + t540;
t450 = -t494 * pkin(3) - t504 * pkin(8) + t506 * t500 + t493;
t457 = -t488 * qJD(4) + t524 * t494 + t555 * t495;
t422 = (t488 * t518 - t457) * qJ(5) + (t489 * t518 + t456) * pkin(4) + t450;
t519 = sin(pkin(11));
t521 = cos(pkin(11));
t476 = t489 * t521 + t518 * t519;
t414 = -0.2e1 * qJD(5) * t476 - t519 * t419 + t521 * t422;
t447 = t457 * t521 + t515 * t519;
t475 = -t489 * t519 + t518 * t521;
t412 = (t475 * t488 - t447) * pkin(9) + (t475 * t476 + t456) * pkin(5) + t414;
t415 = 0.2e1 * qJD(5) * t475 + t521 * t419 + t519 * t422;
t446 = -t457 * t519 + t515 * t521;
t465 = pkin(5) * t488 - pkin(9) * t476;
t474 = t475 ^ 2;
t413 = -pkin(5) * t474 + pkin(9) * t446 - t465 * t488 + t415;
t523 = sin(qJ(6));
t527 = cos(qJ(6));
t410 = t412 * t527 - t413 * t523;
t458 = t475 * t527 - t476 * t523;
t425 = qJD(6) * t458 + t446 * t523 + t447 * t527;
t459 = t475 * t523 + t476 * t527;
t433 = -mrSges(7,1) * t458 + mrSges(7,2) * t459;
t484 = qJD(6) + t488;
t439 = -mrSges(7,2) * t484 + mrSges(7,3) * t458;
t455 = qJDD(6) + t456;
t406 = m(7) * t410 + mrSges(7,1) * t455 - mrSges(7,3) * t425 - t433 * t459 + t439 * t484;
t411 = t412 * t523 + t413 * t527;
t424 = -qJD(6) * t459 + t446 * t527 - t447 * t523;
t440 = mrSges(7,1) * t484 - mrSges(7,3) * t459;
t407 = m(7) * t411 - mrSges(7,2) * t455 + mrSges(7,3) * t424 + t433 * t458 - t440 * t484;
t398 = t527 * t406 + t523 * t407;
t460 = -mrSges(6,1) * t475 + mrSges(6,2) * t476;
t463 = -mrSges(6,2) * t488 + mrSges(6,3) * t475;
t396 = m(6) * t414 + mrSges(6,1) * t456 - mrSges(6,3) * t447 - t460 * t476 + t463 * t488 + t398;
t464 = mrSges(6,1) * t488 - mrSges(6,3) * t476;
t542 = -t406 * t523 + t527 * t407;
t397 = m(6) * t415 - mrSges(6,2) * t456 + mrSges(6,3) * t446 + t460 * t475 - t464 * t488 + t542;
t543 = -t396 * t519 + t521 * t397;
t390 = m(5) * t428 - mrSges(5,2) * t515 - mrSges(5,3) * t456 - t471 * t488 - t480 * t518 + t543;
t427 = t555 * t435 - t524 * t438;
t418 = -t515 * pkin(4) - t514 * qJ(5) + t489 * t470 + qJDD(5) - t427;
t416 = -t446 * pkin(5) - t474 * pkin(9) + t476 * t465 + t418;
t535 = m(7) * t416 - t424 * mrSges(7,1) + mrSges(7,2) * t425 - t458 * t439 + t440 * t459;
t409 = m(6) * t418 - t446 * mrSges(6,1) + mrSges(6,2) * t447 - t475 * t463 + t464 * t476 + t535;
t479 = -mrSges(5,2) * t518 - mrSges(5,3) * t488;
t402 = m(5) * t427 + mrSges(5,1) * t515 - mrSges(5,3) * t457 - t471 * t489 + t479 * t518 - t409;
t382 = t524 * t390 + t555 * t402;
t492 = -mrSges(4,1) * t505 + mrSges(4,2) * t506;
t498 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t505;
t380 = m(4) * t461 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t495 + qJD(3) * t498 - t492 * t506 + t382;
t499 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t506;
t544 = t555 * t390 - t402 * t524;
t381 = m(4) * t462 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t494 - qJD(3) * t499 + t492 * t505 + t544;
t552 = t528 * t380 + t525 * t381;
t392 = t521 * t396 + t519 * t397;
t545 = -t525 * t380 + t528 * t381;
t539 = -mrSges(3,1) * t522 + mrSges(3,2) * t520;
t536 = mrSges(3,3) * qJDD(1) + t530 * t539;
t534 = m(5) * t450 + t456 * mrSges(5,1) + t457 * mrSges(5,2) + t488 * t479 + t489 * t480 + t392;
t429 = Ifges(7,5) * t459 + Ifges(7,6) * t458 + Ifges(7,3) * t484;
t431 = Ifges(7,1) * t459 + Ifges(7,4) * t458 + Ifges(7,5) * t484;
t399 = -mrSges(7,1) * t416 + mrSges(7,3) * t411 + Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t455 - t429 * t459 + t431 * t484;
t430 = Ifges(7,4) * t459 + Ifges(7,2) * t458 + Ifges(7,6) * t484;
t400 = mrSges(7,2) * t416 - mrSges(7,3) * t410 + Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t455 + t429 * t458 - t430 * t484;
t441 = Ifges(6,5) * t476 + Ifges(6,6) * t475 + Ifges(6,3) * t488;
t443 = Ifges(6,1) * t476 + Ifges(6,4) * t475 + Ifges(6,5) * t488;
t384 = -mrSges(6,1) * t418 + mrSges(6,3) * t415 + Ifges(6,4) * t447 + Ifges(6,2) * t446 + Ifges(6,6) * t456 - pkin(5) * t535 + pkin(9) * t542 + t527 * t399 + t523 * t400 - t476 * t441 + t488 * t443;
t442 = Ifges(6,4) * t476 + Ifges(6,2) * t475 + Ifges(6,6) * t488;
t386 = mrSges(6,2) * t418 - mrSges(6,3) * t414 + Ifges(6,1) * t447 + Ifges(6,4) * t446 + Ifges(6,5) * t456 - pkin(9) * t398 - t399 * t523 + t400 * t527 + t441 * t475 - t442 * t488;
t467 = Ifges(5,4) * t489 - Ifges(5,2) * t488 + Ifges(5,6) * t518;
t468 = Ifges(5,1) * t489 - Ifges(5,4) * t488 + Ifges(5,5) * t518;
t533 = mrSges(5,1) * t427 - mrSges(5,2) * t428 + Ifges(5,5) * t457 - Ifges(5,6) * t456 + Ifges(5,3) * t515 - pkin(4) * t409 + qJ(5) * t543 + t521 * t384 + t519 * t386 + t489 * t467 + t468 * t488;
t532 = mrSges(7,1) * t410 - mrSges(7,2) * t411 + Ifges(7,5) * t425 + Ifges(7,6) * t424 + Ifges(7,3) * t455 + t459 * t430 - t458 * t431;
t531 = m(4) * t493 - t494 * mrSges(4,1) + t495 * mrSges(4,2) - t505 * t498 + t506 * t499 + t534;
t503 = -qJDD(1) * pkin(1) - t530 * qJ(2) + t540;
t496 = -t520 * t507 + t546;
t487 = Ifges(4,1) * t506 + Ifges(4,4) * t505 + Ifges(4,5) * qJD(3);
t486 = Ifges(4,4) * t506 + Ifges(4,2) * t505 + Ifges(4,6) * qJD(3);
t485 = Ifges(4,5) * t506 + Ifges(4,6) * t505 + Ifges(4,3) * qJD(3);
t466 = Ifges(5,5) * t489 - Ifges(5,6) * t488 + Ifges(5,3) * t518;
t387 = t551 * t530 * mrSges(3,3) + m(3) * t503 + t539 * qJDD(1) + t531;
t376 = -t532 + Ifges(5,6) * t515 + t518 * t468 - t489 * t466 + t475 * t443 - t476 * t442 + Ifges(5,4) * t457 - Ifges(6,6) * t446 - Ifges(6,5) * t447 - mrSges(5,1) * t450 + mrSges(5,3) * t428 + mrSges(6,2) * t415 - mrSges(6,1) * t414 - pkin(4) * t392 + (-Ifges(5,2) - Ifges(6,3)) * t456 - pkin(5) * t398;
t375 = mrSges(5,2) * t450 - mrSges(5,3) * t427 + Ifges(5,1) * t457 - Ifges(5,4) * t456 + Ifges(5,5) * t515 - qJ(5) * t392 - t384 * t519 + t386 * t521 - t466 * t488 - t467 * t518;
t374 = mrSges(4,2) * t493 - mrSges(4,3) * t461 + Ifges(4,1) * t495 + Ifges(4,4) * t494 + Ifges(4,5) * qJDD(3) - pkin(8) * t382 - qJD(3) * t486 + t555 * t375 - t524 * t376 + t505 * t485;
t373 = -mrSges(4,1) * t493 + mrSges(4,3) * t462 + Ifges(4,4) * t495 + Ifges(4,2) * t494 + Ifges(4,6) * qJDD(3) - pkin(3) * t534 + pkin(8) * t544 + qJD(3) * t487 + t524 * t375 + t555 * t376 - t506 * t485;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t547 - mrSges(2,2) * t541 + t520 * (mrSges(3,2) * t503 - mrSges(3,3) * t496 + t528 * t374 - t525 * t373 - pkin(7) * t552 + (Ifges(3,1) * t520 + Ifges(3,4) * t522) * qJDD(1)) + t522 * (-mrSges(3,1) * t503 + mrSges(3,3) * t497 + t525 * t374 + t528 * t373 - pkin(2) * t531 + pkin(7) * t545 + (Ifges(3,4) * t520 + Ifges(3,2) * t522) * qJDD(1)) - pkin(1) * t387 + qJ(2) * ((m(3) * t497 + t536 * t522 + t545) * t522 + (-m(3) * t496 + t536 * t520 - t552) * t520); t387; mrSges(4,1) * t461 - mrSges(4,2) * t462 + Ifges(4,5) * t495 + Ifges(4,6) * t494 + Ifges(4,3) * qJDD(3) + pkin(3) * t382 + t486 * t506 - t487 * t505 + t533; t533; t409; t532;];
tauJ  = t1;
