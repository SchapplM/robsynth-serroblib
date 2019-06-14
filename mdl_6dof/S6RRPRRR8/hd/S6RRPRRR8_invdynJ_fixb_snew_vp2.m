% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRRR8
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
% Datum: 2019-05-06 22:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:38:50
% EndTime: 2019-05-06 22:39:02
% DurationCPUTime: 11.77s
% Computational Cost: add. (181773->349), mult. (400703->440), div. (0->0), fcn. (299159->12), ass. (0->137)
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
t468 = (-t512 - t548) * qJ(3) + (-t513 + t519) * pkin(2) + t503;
t542 = -g(1) * t533 - g(2) * t528;
t504 = -pkin(1) * t535 + qJDD(1) * pkin(7) + t542;
t485 = -g(3) * t527 + t532 * t504;
t510 = (-pkin(2) * t532 - qJ(3) * t527) * qJD(1);
t534 = qJD(2) ^ 2;
t550 = qJD(1) * t532;
t471 = -pkin(2) * t534 + qJDD(2) * qJ(3) + t510 * t550 + t485;
t522 = sin(pkin(11));
t523 = cos(pkin(11));
t551 = qJD(1) * t527;
t508 = qJD(2) * t522 + t523 * t551;
t448 = -0.2e1 * qJD(3) * t508 + t523 * t468 - t471 * t522;
t490 = qJDD(2) * t522 + t512 * t523;
t507 = qJD(2) * t523 - t522 * t551;
t435 = (-t507 * t550 - t490) * pkin(8) + (t507 * t508 - t513) * pkin(3) + t448;
t449 = 0.2e1 * qJD(3) * t507 + t522 * t468 + t523 * t471;
t489 = qJDD(2) * t523 - t512 * t522;
t491 = -pkin(3) * t550 - pkin(8) * t508;
t506 = t507 ^ 2;
t437 = -pkin(3) * t506 + pkin(8) * t489 + t491 * t550 + t449;
t526 = sin(qJ(4));
t531 = cos(qJ(4));
t416 = t531 * t435 - t437 * t526;
t481 = t507 * t531 - t508 * t526;
t456 = qJD(4) * t481 + t489 * t526 + t490 * t531;
t482 = t507 * t526 + t508 * t531;
t509 = qJDD(4) - t513;
t518 = qJD(4) - t550;
t412 = (t481 * t518 - t456) * pkin(9) + (t481 * t482 + t509) * pkin(4) + t416;
t417 = t526 * t435 + t531 * t437;
t455 = -qJD(4) * t482 + t489 * t531 - t490 * t526;
t474 = pkin(4) * t518 - pkin(9) * t482;
t480 = t481 ^ 2;
t414 = -pkin(4) * t480 + pkin(9) * t455 - t474 * t518 + t417;
t525 = sin(qJ(5));
t530 = cos(qJ(5));
t400 = t530 * t412 - t414 * t525;
t463 = t481 * t530 - t482 * t525;
t431 = qJD(5) * t463 + t455 * t525 + t456 * t530;
t464 = t481 * t525 + t482 * t530;
t505 = qJDD(5) + t509;
t517 = qJD(5) + t518;
t397 = (t463 * t517 - t431) * pkin(10) + (t463 * t464 + t505) * pkin(5) + t400;
t401 = t525 * t412 + t530 * t414;
t430 = -qJD(5) * t464 + t455 * t530 - t456 * t525;
t454 = pkin(5) * t517 - pkin(10) * t464;
t462 = t463 ^ 2;
t398 = -pkin(5) * t462 + pkin(10) * t430 - t454 * t517 + t401;
t524 = sin(qJ(6));
t529 = cos(qJ(6));
t395 = t397 * t529 - t398 * t524;
t445 = t463 * t529 - t464 * t524;
t409 = qJD(6) * t445 + t430 * t524 + t431 * t529;
t446 = t463 * t524 + t464 * t529;
t425 = -mrSges(7,1) * t445 + mrSges(7,2) * t446;
t514 = qJD(6) + t517;
t438 = -mrSges(7,2) * t514 + mrSges(7,3) * t445;
t497 = qJDD(6) + t505;
t391 = m(7) * t395 + mrSges(7,1) * t497 - mrSges(7,3) * t409 - t425 * t446 + t438 * t514;
t396 = t397 * t524 + t398 * t529;
t408 = -qJD(6) * t446 + t430 * t529 - t431 * t524;
t439 = mrSges(7,1) * t514 - mrSges(7,3) * t446;
t392 = m(7) * t396 - mrSges(7,2) * t497 + mrSges(7,3) * t408 + t425 * t445 - t439 * t514;
t385 = t529 * t391 + t524 * t392;
t447 = -mrSges(6,1) * t463 + mrSges(6,2) * t464;
t451 = -mrSges(6,2) * t517 + mrSges(6,3) * t463;
t382 = m(6) * t400 + mrSges(6,1) * t505 - mrSges(6,3) * t431 - t447 * t464 + t451 * t517 + t385;
t452 = mrSges(6,1) * t517 - mrSges(6,3) * t464;
t543 = -t391 * t524 + t529 * t392;
t383 = m(6) * t401 - mrSges(6,2) * t505 + mrSges(6,3) * t430 + t447 * t463 - t452 * t517 + t543;
t378 = t530 * t382 + t525 * t383;
t465 = -mrSges(5,1) * t481 + mrSges(5,2) * t482;
t472 = -mrSges(5,2) * t518 + mrSges(5,3) * t481;
t376 = m(5) * t416 + mrSges(5,1) * t509 - mrSges(5,3) * t456 - t465 * t482 + t472 * t518 + t378;
t473 = mrSges(5,1) * t518 - mrSges(5,3) * t482;
t544 = -t382 * t525 + t530 * t383;
t377 = m(5) * t417 - mrSges(5,2) * t509 + mrSges(5,3) * t455 + t465 * t481 - t473 * t518 + t544;
t370 = t531 * t376 + t526 * t377;
t484 = -t532 * g(3) - t527 * t504;
t483 = -mrSges(4,1) * t507 + mrSges(4,2) * t508;
t487 = mrSges(4,2) * t550 + mrSges(4,3) * t507;
t368 = m(4) * t448 - mrSges(4,1) * t513 - mrSges(4,3) * t490 - t483 * t508 - t487 * t550 + t370;
t488 = -mrSges(4,1) * t550 - mrSges(4,3) * t508;
t545 = -t376 * t526 + t531 * t377;
t369 = m(4) * t449 + mrSges(4,2) * t513 + mrSges(4,3) * t489 + t483 * t507 + t488 * t550 + t545;
t546 = -t368 * t522 + t523 * t369;
t470 = -qJDD(2) * pkin(2) - qJ(3) * t534 + t510 * t551 + qJDD(3) - t484;
t450 = -pkin(3) * t489 - pkin(8) * t506 + t508 * t491 + t470;
t424 = -pkin(4) * t455 - pkin(9) * t480 + t482 * t474 + t450;
t403 = -pkin(5) * t430 - pkin(10) * t462 + t454 * t464 + t424;
t541 = m(7) * t403 - t408 * mrSges(7,1) + t409 * mrSges(7,2) - t445 * t438 + t446 * t439;
t364 = t368 * t523 + t369 * t522;
t540 = m(6) * t424 - t430 * mrSges(6,1) + t431 * mrSges(6,2) - t463 * t451 + t464 * t452 + t541;
t419 = Ifges(7,4) * t446 + Ifges(7,2) * t445 + Ifges(7,6) * t514;
t420 = Ifges(7,1) * t446 + Ifges(7,4) * t445 + Ifges(7,5) * t514;
t539 = -mrSges(7,1) * t395 + mrSges(7,2) * t396 - Ifges(7,5) * t409 - Ifges(7,6) * t408 - Ifges(7,3) * t497 - t446 * t419 + t445 * t420;
t538 = m(5) * t450 - t455 * mrSges(5,1) + t456 * mrSges(5,2) - t481 * t472 + t482 * t473 + t540;
t441 = Ifges(6,4) * t464 + Ifges(6,2) * t463 + Ifges(6,6) * t517;
t442 = Ifges(6,1) * t464 + Ifges(6,4) * t463 + Ifges(6,5) * t517;
t537 = -mrSges(6,1) * t400 + mrSges(6,2) * t401 - Ifges(6,5) * t431 - Ifges(6,6) * t430 - Ifges(6,3) * t505 - pkin(5) * t385 - t464 * t441 + t463 * t442 + t539;
t393 = m(4) * t470 - t489 * mrSges(4,1) + t490 * mrSges(4,2) - t507 * t487 + t508 * t488 + t538;
t458 = Ifges(5,4) * t482 + Ifges(5,2) * t481 + Ifges(5,6) * t518;
t459 = Ifges(5,1) * t482 + Ifges(5,4) * t481 + Ifges(5,5) * t518;
t536 = mrSges(5,1) * t416 - mrSges(5,2) * t417 + Ifges(5,5) * t456 + Ifges(5,6) * t455 + Ifges(5,3) * t509 + pkin(4) * t378 + t482 * t458 - t481 * t459 - t537;
t516 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t550;
t515 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t551;
t511 = (-mrSges(3,1) * t532 + mrSges(3,2) * t527) * qJD(1);
t502 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t527 + Ifges(3,4) * t532) * qJD(1);
t501 = Ifges(3,6) * qJD(2) + (t527 * Ifges(3,4) + Ifges(3,2) * t532) * qJD(1);
t477 = Ifges(4,1) * t508 + Ifges(4,4) * t507 - Ifges(4,5) * t550;
t476 = Ifges(4,4) * t508 + Ifges(4,2) * t507 - Ifges(4,6) * t550;
t475 = Ifges(4,5) * t508 + Ifges(4,6) * t507 - Ifges(4,3) * t550;
t457 = Ifges(5,5) * t482 + Ifges(5,6) * t481 + Ifges(5,3) * t518;
t440 = Ifges(6,5) * t464 + Ifges(6,6) * t463 + Ifges(6,3) * t517;
t418 = Ifges(7,5) * t446 + Ifges(7,6) * t445 + Ifges(7,3) * t514;
t387 = mrSges(7,2) * t403 - mrSges(7,3) * t395 + Ifges(7,1) * t409 + Ifges(7,4) * t408 + Ifges(7,5) * t497 + t418 * t445 - t419 * t514;
t386 = -mrSges(7,1) * t403 + mrSges(7,3) * t396 + Ifges(7,4) * t409 + Ifges(7,2) * t408 + Ifges(7,6) * t497 - t418 * t446 + t420 * t514;
t372 = mrSges(6,2) * t424 - mrSges(6,3) * t400 + Ifges(6,1) * t431 + Ifges(6,4) * t430 + Ifges(6,5) * t505 - pkin(10) * t385 - t386 * t524 + t387 * t529 + t440 * t463 - t441 * t517;
t371 = -mrSges(6,1) * t424 + mrSges(6,3) * t401 + Ifges(6,4) * t431 + Ifges(6,2) * t430 + Ifges(6,6) * t505 - pkin(5) * t541 + pkin(10) * t543 + t529 * t386 + t524 * t387 - t464 * t440 + t517 * t442;
t366 = mrSges(5,2) * t450 - mrSges(5,3) * t416 + Ifges(5,1) * t456 + Ifges(5,4) * t455 + Ifges(5,5) * t509 - pkin(9) * t378 - t371 * t525 + t372 * t530 + t457 * t481 - t458 * t518;
t365 = -mrSges(5,1) * t450 + mrSges(5,3) * t417 + Ifges(5,4) * t456 + Ifges(5,2) * t455 + Ifges(5,6) * t509 - pkin(4) * t540 + pkin(9) * t544 + t530 * t371 + t525 * t372 - t482 * t457 + t518 * t459;
t363 = mrSges(4,2) * t470 - mrSges(4,3) * t448 + Ifges(4,1) * t490 + Ifges(4,4) * t489 - Ifges(4,5) * t513 - pkin(8) * t370 - t365 * t526 + t366 * t531 + t475 * t507 + t476 * t550;
t362 = -mrSges(4,1) * t470 + mrSges(4,3) * t449 + Ifges(4,4) * t490 + Ifges(4,2) * t489 - Ifges(4,6) * t513 - pkin(3) * t538 + pkin(8) * t545 + t531 * t365 + t526 * t366 - t508 * t475 - t477 * t550;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t547 - mrSges(2,2) * t542 + t527 * (mrSges(3,2) * t503 - mrSges(3,3) * t484 + Ifges(3,1) * t512 + Ifges(3,4) * t513 + Ifges(3,5) * qJDD(2) - qJ(3) * t364 - qJD(2) * t501 - t522 * t362 + t523 * t363) + t532 * (-mrSges(3,1) * t503 - mrSges(4,1) * t448 + mrSges(4,2) * t449 + mrSges(3,3) * t485 + Ifges(3,4) * t512 - Ifges(4,5) * t490 + Ifges(3,6) * qJDD(2) - Ifges(4,6) * t489 - pkin(2) * t364 - pkin(3) * t370 + qJD(2) * t502 - t508 * t476 + t507 * t477 - t536 + (Ifges(3,2) + Ifges(4,3)) * t513) + pkin(1) * (-m(3) * t503 + mrSges(3,1) * t513 - mrSges(3,2) * t512 + (-t515 * t527 + t516 * t532) * qJD(1) - t364) + pkin(7) * (t532 * (m(3) * t485 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t513 - qJD(2) * t515 + t511 * t550 + t546) - t527 * (m(3) * t484 + qJDD(2) * mrSges(3,1) - t512 * mrSges(3,3) + qJD(2) * t516 - t511 * t551 - t393)); Ifges(3,5) * t512 + Ifges(3,6) * t513 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t484 - mrSges(3,2) * t485 + t522 * t363 + t523 * t362 - pkin(2) * t393 + qJ(3) * t546 + (t527 * t501 - t532 * t502) * qJD(1); t393; t536; -t537; -t539;];
tauJ  = t1;
