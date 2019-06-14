% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPRR6
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
% Datum: 2019-05-07 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:00:07
% EndTime: 2019-05-07 11:00:35
% DurationCPUTime: 13.03s
% Computational Cost: add. (197856->349), mult. (419879->440), div. (0->0), fcn. (310907->12), ass. (0->137)
t535 = qJD(1) ^ 2;
t528 = sin(qJ(1));
t533 = cos(qJ(1));
t547 = g(1) * t528 - t533 * g(2);
t503 = -qJDD(1) * pkin(1) - pkin(7) * t535 - t547;
t527 = sin(qJ(2));
t532 = cos(qJ(2));
t549 = qJD(1) * qJD(2);
t548 = t532 * t549;
t512 = qJDD(1) * t527 + t548;
t519 = t527 * t549;
t513 = qJDD(1) * t532 - t519;
t468 = (-t512 - t548) * pkin(8) + (-t513 + t519) * pkin(2) + t503;
t542 = -g(1) * t533 - g(2) * t528;
t504 = -pkin(1) * t535 + qJDD(1) * pkin(7) + t542;
t493 = -g(3) * t527 + t532 * t504;
t511 = (-pkin(2) * t532 - pkin(8) * t527) * qJD(1);
t534 = qJD(2) ^ 2;
t551 = qJD(1) * t532;
t471 = -pkin(2) * t534 + qJDD(2) * pkin(8) + t511 * t551 + t493;
t526 = sin(qJ(3));
t531 = cos(qJ(3));
t449 = t531 * t468 - t471 * t526;
t550 = t527 * qJD(1);
t508 = qJD(2) * t531 - t526 * t550;
t484 = qJD(3) * t508 + qJDD(2) * t526 + t512 * t531;
t507 = qJDD(3) - t513;
t509 = qJD(2) * t526 + t531 * t550;
t518 = qJD(3) - t551;
t435 = (t508 * t518 - t484) * qJ(4) + (t508 * t509 + t507) * pkin(3) + t449;
t450 = t526 * t468 + t531 * t471;
t483 = -qJD(3) * t509 + qJDD(2) * t531 - t512 * t526;
t490 = pkin(3) * t518 - qJ(4) * t509;
t506 = t508 ^ 2;
t437 = -pkin(3) * t506 + qJ(4) * t483 - t490 * t518 + t450;
t522 = sin(pkin(11));
t523 = cos(pkin(11));
t487 = t508 * t522 + t509 * t523;
t416 = -0.2e1 * qJD(4) * t487 + t523 * t435 - t437 * t522;
t459 = t483 * t522 + t484 * t523;
t486 = t508 * t523 - t509 * t522;
t406 = (t486 * t518 - t459) * pkin(9) + (t486 * t487 + t507) * pkin(4) + t416;
t417 = 0.2e1 * qJD(4) * t486 + t522 * t435 + t523 * t437;
t458 = t483 * t523 - t484 * t522;
t474 = pkin(4) * t518 - pkin(9) * t487;
t485 = t486 ^ 2;
t414 = -pkin(4) * t485 + pkin(9) * t458 - t474 * t518 + t417;
t525 = sin(qJ(5));
t530 = cos(qJ(5));
t400 = t530 * t406 - t414 * t525;
t463 = t486 * t530 - t487 * t525;
t431 = qJD(5) * t463 + t458 * t525 + t459 * t530;
t464 = t486 * t525 + t487 * t530;
t505 = qJDD(5) + t507;
t517 = qJD(5) + t518;
t397 = (t463 * t517 - t431) * pkin(10) + (t463 * t464 + t505) * pkin(5) + t400;
t401 = t525 * t406 + t530 * t414;
t430 = -qJD(5) * t464 + t458 * t530 - t459 * t525;
t453 = pkin(5) * t517 - pkin(10) * t464;
t462 = t463 ^ 2;
t398 = -pkin(5) * t462 + pkin(10) * t430 - t453 * t517 + t401;
t524 = sin(qJ(6));
t529 = cos(qJ(6));
t395 = t397 * t529 - t398 * t524;
t445 = t463 * t529 - t464 * t524;
t412 = qJD(6) * t445 + t430 * t524 + t431 * t529;
t446 = t463 * t524 + t464 * t529;
t425 = -mrSges(7,1) * t445 + mrSges(7,2) * t446;
t514 = qJD(6) + t517;
t438 = -mrSges(7,2) * t514 + mrSges(7,3) * t445;
t499 = qJDD(6) + t505;
t391 = m(7) * t395 + mrSges(7,1) * t499 - mrSges(7,3) * t412 - t425 * t446 + t438 * t514;
t396 = t397 * t524 + t398 * t529;
t411 = -qJD(6) * t446 + t430 * t529 - t431 * t524;
t439 = mrSges(7,1) * t514 - mrSges(7,3) * t446;
t392 = m(7) * t396 - mrSges(7,2) * t499 + mrSges(7,3) * t411 + t425 * t445 - t439 * t514;
t385 = t529 * t391 + t524 * t392;
t447 = -mrSges(6,1) * t463 + mrSges(6,2) * t464;
t451 = -mrSges(6,2) * t517 + mrSges(6,3) * t463;
t382 = m(6) * t400 + mrSges(6,1) * t505 - mrSges(6,3) * t431 - t447 * t464 + t451 * t517 + t385;
t452 = mrSges(6,1) * t517 - mrSges(6,3) * t464;
t543 = -t391 * t524 + t529 * t392;
t383 = m(6) * t401 - mrSges(6,2) * t505 + mrSges(6,3) * t430 + t447 * t463 - t452 * t517 + t543;
t378 = t530 * t382 + t525 * t383;
t465 = -mrSges(5,1) * t486 + mrSges(5,2) * t487;
t472 = -mrSges(5,2) * t518 + mrSges(5,3) * t486;
t376 = m(5) * t416 + mrSges(5,1) * t507 - mrSges(5,3) * t459 - t465 * t487 + t472 * t518 + t378;
t473 = mrSges(5,1) * t518 - mrSges(5,3) * t487;
t544 = -t382 * t525 + t530 * t383;
t377 = m(5) * t417 - mrSges(5,2) * t507 + mrSges(5,3) * t458 + t465 * t486 - t473 * t518 + t544;
t370 = t523 * t376 + t522 * t377;
t456 = Ifges(5,4) * t487 + Ifges(5,2) * t486 + Ifges(5,6) * t518;
t457 = Ifges(5,1) * t487 + Ifges(5,4) * t486 + Ifges(5,5) * t518;
t476 = Ifges(4,4) * t509 + Ifges(4,2) * t508 + Ifges(4,6) * t518;
t477 = Ifges(4,1) * t509 + Ifges(4,4) * t508 + Ifges(4,5) * t518;
t441 = Ifges(6,4) * t464 + Ifges(6,2) * t463 + Ifges(6,6) * t517;
t442 = Ifges(6,1) * t464 + Ifges(6,4) * t463 + Ifges(6,5) * t517;
t421 = Ifges(7,4) * t446 + Ifges(7,2) * t445 + Ifges(7,6) * t514;
t422 = Ifges(7,1) * t446 + Ifges(7,4) * t445 + Ifges(7,5) * t514;
t539 = -mrSges(7,1) * t395 + mrSges(7,2) * t396 - Ifges(7,5) * t412 - Ifges(7,6) * t411 - Ifges(7,3) * t499 - t446 * t421 + t445 * t422;
t538 = -mrSges(6,1) * t400 + mrSges(6,2) * t401 - Ifges(6,5) * t431 - Ifges(6,6) * t430 - Ifges(6,3) * t505 - pkin(5) * t385 - t464 * t441 + t463 * t442 + t539;
t553 = mrSges(4,1) * t449 + mrSges(5,1) * t416 - mrSges(4,2) * t450 - mrSges(5,2) * t417 + Ifges(4,5) * t484 + Ifges(5,5) * t459 + Ifges(4,6) * t483 + Ifges(5,6) * t458 + pkin(3) * t370 + pkin(4) * t378 + t487 * t456 - t486 * t457 + t509 * t476 - t508 * t477 + (Ifges(4,3) + Ifges(5,3)) * t507 - t538;
t492 = -t532 * g(3) - t527 * t504;
t488 = -mrSges(4,1) * t508 + mrSges(4,2) * t509;
t489 = -mrSges(4,2) * t518 + mrSges(4,3) * t508;
t368 = m(4) * t449 + mrSges(4,1) * t507 - mrSges(4,3) * t484 - t488 * t509 + t489 * t518 + t370;
t491 = mrSges(4,1) * t518 - mrSges(4,3) * t509;
t545 = -t376 * t522 + t523 * t377;
t369 = m(4) * t450 - mrSges(4,2) * t507 + mrSges(4,3) * t483 + t488 * t508 - t491 * t518 + t545;
t546 = -t368 * t526 + t531 * t369;
t470 = -qJDD(2) * pkin(2) - pkin(8) * t534 + t511 * t550 - t492;
t448 = -pkin(3) * t483 - qJ(4) * t506 + t509 * t490 + qJDD(4) + t470;
t419 = -pkin(4) * t458 - pkin(9) * t485 + t487 * t474 + t448;
t403 = -pkin(5) * t430 - pkin(10) * t462 + t453 * t464 + t419;
t541 = m(7) * t403 - t411 * mrSges(7,1) + t412 * mrSges(7,2) - t445 * t438 + t446 * t439;
t364 = t368 * t531 + t369 * t526;
t540 = m(6) * t419 - t430 * mrSges(6,1) + t431 * mrSges(6,2) - t463 * t451 + t464 * t452 + t541;
t393 = m(5) * t448 - t458 * mrSges(5,1) + t459 * mrSges(5,2) - t486 * t472 + t487 * t473 + t540;
t537 = -m(4) * t470 + t483 * mrSges(4,1) - t484 * mrSges(4,2) + t508 * t489 - t509 * t491 - t393;
t516 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t551;
t515 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t550;
t510 = (-t532 * mrSges(3,1) + t527 * mrSges(3,2)) * qJD(1);
t502 = Ifges(3,5) * qJD(2) + (t527 * Ifges(3,1) + t532 * Ifges(3,4)) * qJD(1);
t501 = Ifges(3,6) * qJD(2) + (t527 * Ifges(3,4) + t532 * Ifges(3,2)) * qJD(1);
t475 = Ifges(4,5) * t509 + Ifges(4,6) * t508 + Ifges(4,3) * t518;
t455 = Ifges(5,5) * t487 + Ifges(5,6) * t486 + Ifges(5,3) * t518;
t440 = Ifges(6,5) * t464 + Ifges(6,6) * t463 + Ifges(6,3) * t517;
t420 = Ifges(7,5) * t446 + Ifges(7,6) * t445 + Ifges(7,3) * t514;
t387 = mrSges(7,2) * t403 - mrSges(7,3) * t395 + Ifges(7,1) * t412 + Ifges(7,4) * t411 + Ifges(7,5) * t499 + t420 * t445 - t421 * t514;
t386 = -mrSges(7,1) * t403 + mrSges(7,3) * t396 + Ifges(7,4) * t412 + Ifges(7,2) * t411 + Ifges(7,6) * t499 - t420 * t446 + t422 * t514;
t372 = mrSges(6,2) * t419 - mrSges(6,3) * t400 + Ifges(6,1) * t431 + Ifges(6,4) * t430 + Ifges(6,5) * t505 - pkin(10) * t385 - t386 * t524 + t387 * t529 + t440 * t463 - t441 * t517;
t371 = -mrSges(6,1) * t419 + mrSges(6,3) * t401 + Ifges(6,4) * t431 + Ifges(6,2) * t430 + Ifges(6,6) * t505 - pkin(5) * t541 + pkin(10) * t543 + t529 * t386 + t524 * t387 - t464 * t440 + t517 * t442;
t366 = mrSges(5,2) * t448 - mrSges(5,3) * t416 + Ifges(5,1) * t459 + Ifges(5,4) * t458 + Ifges(5,5) * t507 - pkin(9) * t378 - t371 * t525 + t372 * t530 + t455 * t486 - t456 * t518;
t365 = -mrSges(5,1) * t448 + mrSges(5,3) * t417 + Ifges(5,4) * t459 + Ifges(5,2) * t458 + Ifges(5,6) * t507 - pkin(4) * t540 + pkin(9) * t544 + t530 * t371 + t525 * t372 - t487 * t455 + t518 * t457;
t363 = mrSges(4,2) * t470 - mrSges(4,3) * t449 + Ifges(4,1) * t484 + Ifges(4,4) * t483 + Ifges(4,5) * t507 - qJ(4) * t370 - t365 * t522 + t366 * t523 + t475 * t508 - t476 * t518;
t362 = -mrSges(4,1) * t470 + mrSges(4,3) * t450 + Ifges(4,4) * t484 + Ifges(4,2) * t483 + Ifges(4,6) * t507 - pkin(3) * t393 + qJ(4) * t545 + t523 * t365 + t522 * t366 - t509 * t475 + t518 * t477;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t547 - mrSges(2,2) * t542 + t527 * (mrSges(3,2) * t503 - mrSges(3,3) * t492 + Ifges(3,1) * t512 + Ifges(3,4) * t513 + Ifges(3,5) * qJDD(2) - pkin(8) * t364 - qJD(2) * t501 - t526 * t362 + t531 * t363) + t532 * (-mrSges(3,1) * t503 + mrSges(3,3) * t493 + Ifges(3,4) * t512 + Ifges(3,2) * t513 + Ifges(3,6) * qJDD(2) - pkin(2) * t364 + qJD(2) * t502 - t553) + pkin(1) * (-m(3) * t503 + mrSges(3,1) * t513 - mrSges(3,2) * t512 + (-t515 * t527 + t516 * t532) * qJD(1) - t364) + pkin(7) * (t532 * (m(3) * t493 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t513 - qJD(2) * t515 + t510 * t551 + t546) - t527 * (m(3) * t492 + qJDD(2) * mrSges(3,1) - t512 * mrSges(3,3) + qJD(2) * t516 - t510 * t550 + t537)); Ifges(3,5) * t512 + Ifges(3,6) * t513 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t492 - mrSges(3,2) * t493 + t526 * t363 + t531 * t362 + pkin(2) * t537 + pkin(8) * t546 + (t527 * t501 - t532 * t502) * qJD(1); t553; t393; -t538; -t539;];
tauJ  = t1;
