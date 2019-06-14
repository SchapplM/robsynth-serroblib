% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR1
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
% Datum: 2019-05-06 12:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:48:30
% EndTime: 2019-05-06 12:48:39
% DurationCPUTime: 8.46s
% Computational Cost: add. (120662->348), mult. (283492->447), div. (0->0), fcn. (213116->12), ass. (0->135)
t551 = cos(qJ(4));
t534 = qJD(1) ^ 2;
t550 = pkin(2) * t534;
t530 = sin(qJ(1));
t533 = cos(qJ(1));
t540 = -t533 * g(1) - t530 * g(2);
t509 = -t534 * pkin(1) + qJDD(1) * pkin(7) + t540;
t529 = sin(qJ(2));
t549 = t529 * t509;
t532 = cos(qJ(2));
t546 = qJD(1) * qJD(2);
t512 = t529 * qJDD(1) + t532 * t546;
t472 = qJDD(2) * pkin(2) - t512 * qJ(3) - t549 + (qJ(3) * t546 + t529 * t550 - g(3)) * t532;
t494 = -t529 * g(3) + t532 * t509;
t513 = t532 * qJDD(1) - t529 * t546;
t548 = qJD(1) * t529;
t514 = qJD(2) * pkin(2) - qJ(3) * t548;
t522 = t532 ^ 2;
t473 = t513 * qJ(3) - qJD(2) * t514 - t522 * t550 + t494;
t524 = sin(pkin(10));
t526 = cos(pkin(10));
t504 = (t532 * t524 + t529 * t526) * qJD(1);
t444 = -0.2e1 * qJD(3) * t504 + t526 * t472 - t524 * t473;
t492 = t526 * t512 + t524 * t513;
t503 = (-t529 * t524 + t532 * t526) * qJD(1);
t431 = (qJD(2) * t503 - t492) * pkin(8) + (t503 * t504 + qJDD(2)) * pkin(3) + t444;
t445 = 0.2e1 * qJD(3) * t503 + t524 * t472 + t526 * t473;
t491 = -t524 * t512 + t526 * t513;
t497 = qJD(2) * pkin(3) - t504 * pkin(8);
t502 = t503 ^ 2;
t435 = -t502 * pkin(3) + t491 * pkin(8) - qJD(2) * t497 + t445;
t528 = sin(qJ(4));
t422 = t528 * t431 + t551 * t435;
t487 = t528 * t503 + t504 * t551;
t455 = t487 * qJD(4) - t551 * t491 + t528 * t492;
t486 = -t551 * t503 + t528 * t504;
t468 = t486 * mrSges(5,1) + t487 * mrSges(5,2);
t521 = qJD(2) + qJD(4);
t481 = t521 * mrSges(5,1) - t487 * mrSges(5,3);
t520 = qJDD(2) + qJDD(4);
t467 = t486 * pkin(4) - t487 * qJ(5);
t519 = t521 ^ 2;
t416 = -t519 * pkin(4) + t520 * qJ(5) - t486 * t467 + t422;
t545 = t530 * g(1) - t533 * g(2);
t539 = -qJDD(1) * pkin(1) - t545;
t479 = -t513 * pkin(2) + qJDD(3) + t514 * t548 + (-qJ(3) * t522 - pkin(7)) * t534 + t539;
t443 = -t491 * pkin(3) - t502 * pkin(8) + t504 * t497 + t479;
t456 = -t486 * qJD(4) + t528 * t491 + t492 * t551;
t419 = (t486 * t521 - t456) * qJ(5) + (t487 * t521 + t455) * pkin(4) + t443;
t523 = sin(pkin(11));
t525 = cos(pkin(11));
t478 = t525 * t487 + t523 * t521;
t411 = -0.2e1 * qJD(5) * t478 - t523 * t416 + t525 * t419;
t449 = t525 * t456 + t523 * t520;
t477 = -t523 * t487 + t525 * t521;
t409 = (t477 * t486 - t449) * pkin(9) + (t477 * t478 + t455) * pkin(5) + t411;
t412 = 0.2e1 * qJD(5) * t477 + t525 * t416 + t523 * t419;
t448 = -t523 * t456 + t525 * t520;
t462 = t486 * pkin(5) - t478 * pkin(9);
t476 = t477 ^ 2;
t410 = -t476 * pkin(5) + t448 * pkin(9) - t486 * t462 + t412;
t527 = sin(qJ(6));
t531 = cos(qJ(6));
t407 = t531 * t409 - t527 * t410;
t457 = t531 * t477 - t527 * t478;
t425 = t457 * qJD(6) + t527 * t448 + t531 * t449;
t458 = t527 * t477 + t531 * t478;
t432 = -t457 * mrSges(7,1) + t458 * mrSges(7,2);
t482 = qJD(6) + t486;
t436 = -t482 * mrSges(7,2) + t457 * mrSges(7,3);
t454 = qJDD(6) + t455;
t403 = m(7) * t407 + t454 * mrSges(7,1) - t425 * mrSges(7,3) - t458 * t432 + t482 * t436;
t408 = t527 * t409 + t531 * t410;
t424 = -t458 * qJD(6) + t531 * t448 - t527 * t449;
t437 = t482 * mrSges(7,1) - t458 * mrSges(7,3);
t404 = m(7) * t408 - t454 * mrSges(7,2) + t424 * mrSges(7,3) + t457 * t432 - t482 * t437;
t395 = t531 * t403 + t527 * t404;
t459 = -t477 * mrSges(6,1) + t478 * mrSges(6,2);
t460 = -t486 * mrSges(6,2) + t477 * mrSges(6,3);
t393 = m(6) * t411 + t455 * mrSges(6,1) - t449 * mrSges(6,3) - t478 * t459 + t486 * t460 + t395;
t461 = t486 * mrSges(6,1) - t478 * mrSges(6,3);
t541 = -t527 * t403 + t531 * t404;
t394 = m(6) * t412 - t455 * mrSges(6,2) + t448 * mrSges(6,3) + t477 * t459 - t486 * t461 + t541;
t542 = -t523 * t393 + t525 * t394;
t386 = m(5) * t422 - t520 * mrSges(5,2) - t455 * mrSges(5,3) - t486 * t468 - t521 * t481 + t542;
t421 = t551 * t431 - t528 * t435;
t415 = -t520 * pkin(4) - t519 * qJ(5) + t487 * t467 + qJDD(5) - t421;
t413 = -t448 * pkin(5) - t476 * pkin(9) + t478 * t462 + t415;
t538 = m(7) * t413 - t424 * mrSges(7,1) + t425 * mrSges(7,2) - t457 * t436 + t458 * t437;
t406 = m(6) * t415 - t448 * mrSges(6,1) + t449 * mrSges(6,2) - t477 * t460 + t478 * t461 + t538;
t480 = -t521 * mrSges(5,2) - t486 * mrSges(5,3);
t399 = m(5) * t421 + t520 * mrSges(5,1) - t456 * mrSges(5,3) - t487 * t468 + t521 * t480 - t406;
t379 = t528 * t386 + t551 * t399;
t489 = -t503 * mrSges(4,1) + t504 * mrSges(4,2);
t495 = -qJD(2) * mrSges(4,2) + t503 * mrSges(4,3);
t377 = m(4) * t444 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t492 + qJD(2) * t495 - t489 * t504 + t379;
t496 = qJD(2) * mrSges(4,1) - t504 * mrSges(4,3);
t543 = t551 * t386 - t528 * t399;
t378 = m(4) * t445 - qJDD(2) * mrSges(4,2) + t491 * mrSges(4,3) - qJD(2) * t496 + t503 * t489 + t543;
t372 = t526 * t377 + t524 * t378;
t389 = t525 * t393 + t523 * t394;
t547 = qJD(1) * t532;
t544 = -t377 * t524 + t526 * t378;
t537 = m(5) * t443 + t455 * mrSges(5,1) + t456 * mrSges(5,2) + t486 * t480 + t487 * t481 + t389;
t426 = Ifges(7,5) * t458 + Ifges(7,6) * t457 + Ifges(7,3) * t482;
t428 = Ifges(7,1) * t458 + Ifges(7,4) * t457 + Ifges(7,5) * t482;
t396 = -mrSges(7,1) * t413 + mrSges(7,3) * t408 + Ifges(7,4) * t425 + Ifges(7,2) * t424 + Ifges(7,6) * t454 - t458 * t426 + t482 * t428;
t427 = Ifges(7,4) * t458 + Ifges(7,2) * t457 + Ifges(7,6) * t482;
t397 = mrSges(7,2) * t413 - mrSges(7,3) * t407 + Ifges(7,1) * t425 + Ifges(7,4) * t424 + Ifges(7,5) * t454 + t457 * t426 - t482 * t427;
t438 = Ifges(6,5) * t478 + Ifges(6,6) * t477 + Ifges(6,3) * t486;
t440 = Ifges(6,1) * t478 + Ifges(6,4) * t477 + Ifges(6,5) * t486;
t381 = -mrSges(6,1) * t415 + mrSges(6,3) * t412 + Ifges(6,4) * t449 + Ifges(6,2) * t448 + Ifges(6,6) * t455 - pkin(5) * t538 + pkin(9) * t541 + t531 * t396 + t527 * t397 - t478 * t438 + t486 * t440;
t439 = Ifges(6,4) * t478 + Ifges(6,2) * t477 + Ifges(6,6) * t486;
t383 = mrSges(6,2) * t415 - mrSges(6,3) * t411 + Ifges(6,1) * t449 + Ifges(6,4) * t448 + Ifges(6,5) * t455 - pkin(9) * t395 - t527 * t396 + t531 * t397 + t477 * t438 - t486 * t439;
t464 = Ifges(5,4) * t487 - Ifges(5,2) * t486 + Ifges(5,6) * t521;
t465 = Ifges(5,1) * t487 - Ifges(5,4) * t486 + Ifges(5,5) * t521;
t536 = mrSges(5,1) * t421 - mrSges(5,2) * t422 + Ifges(5,5) * t456 - Ifges(5,6) * t455 + Ifges(5,3) * t520 - pkin(4) * t406 + qJ(5) * t542 + t525 * t381 + t523 * t383 + t487 * t464 + t486 * t465;
t535 = mrSges(7,1) * t407 - mrSges(7,2) * t408 + Ifges(7,5) * t425 + Ifges(7,6) * t424 + Ifges(7,3) * t454 + t458 * t427 - t457 * t428;
t387 = m(4) * t479 - t491 * mrSges(4,1) + t492 * mrSges(4,2) - t503 * t495 + t504 * t496 + t537;
t516 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t547;
t515 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t548;
t511 = (-t532 * mrSges(3,1) + t529 * mrSges(3,2)) * qJD(1);
t508 = -t534 * pkin(7) + t539;
t507 = Ifges(3,5) * qJD(2) + (t529 * Ifges(3,1) + t532 * Ifges(3,4)) * qJD(1);
t506 = Ifges(3,6) * qJD(2) + (t529 * Ifges(3,4) + t532 * Ifges(3,2)) * qJD(1);
t493 = -t532 * g(3) - t549;
t485 = Ifges(4,1) * t504 + Ifges(4,4) * t503 + Ifges(4,5) * qJD(2);
t484 = Ifges(4,4) * t504 + Ifges(4,2) * t503 + Ifges(4,6) * qJD(2);
t483 = Ifges(4,5) * t504 + Ifges(4,6) * t503 + Ifges(4,3) * qJD(2);
t463 = Ifges(5,5) * t487 - Ifges(5,6) * t486 + Ifges(5,3) * t521;
t373 = -t535 + Ifges(5,6) * t520 + t521 * t465 - t487 * t463 + t477 * t440 - t478 * t439 + Ifges(5,4) * t456 - Ifges(6,5) * t449 - mrSges(5,1) * t443 - Ifges(6,6) * t448 + mrSges(5,3) * t422 + mrSges(6,2) * t412 - mrSges(6,1) * t411 + (-Ifges(5,2) - Ifges(6,3)) * t455 - pkin(5) * t395 - pkin(4) * t389;
t371 = mrSges(5,2) * t443 - mrSges(5,3) * t421 + Ifges(5,1) * t456 - Ifges(5,4) * t455 + Ifges(5,5) * t520 - qJ(5) * t389 - t381 * t523 + t383 * t525 - t463 * t486 - t464 * t521;
t370 = mrSges(4,2) * t479 - mrSges(4,3) * t444 + Ifges(4,1) * t492 + Ifges(4,4) * t491 + Ifges(4,5) * qJDD(2) - pkin(8) * t379 - qJD(2) * t484 + t371 * t551 - t528 * t373 + t503 * t483;
t369 = -mrSges(4,1) * t479 + mrSges(4,3) * t445 + Ifges(4,4) * t492 + Ifges(4,2) * t491 + Ifges(4,6) * qJDD(2) - pkin(3) * t537 + pkin(8) * t543 + qJD(2) * t485 + t528 * t371 + t373 * t551 - t504 * t483;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t545 - mrSges(2,2) * t540 + t529 * (mrSges(3,2) * t508 - mrSges(3,3) * t493 + Ifges(3,1) * t512 + Ifges(3,4) * t513 + Ifges(3,5) * qJDD(2) - qJ(3) * t372 - qJD(2) * t506 - t524 * t369 + t526 * t370) + t532 * (-mrSges(3,1) * t508 + mrSges(3,3) * t494 + Ifges(3,4) * t512 + Ifges(3,2) * t513 + Ifges(3,6) * qJDD(2) - pkin(2) * t387 + qJ(3) * t544 + qJD(2) * t507 + t526 * t369 + t524 * t370) + pkin(1) * ((-t515 * t529 + t516 * t532) * qJD(1) - t387 - t512 * mrSges(3,2) + t513 * mrSges(3,1) - m(3) * t508) + pkin(7) * (t532 * (m(3) * t494 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t513 - qJD(2) * t515 + t511 * t547 + t544) - t529 * (m(3) * t493 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t512 + qJD(2) * t516 - t511 * t548 + t372)); pkin(3) * t379 + t536 + pkin(2) * t372 + Ifges(3,5) * t512 + Ifges(3,6) * t513 - t503 * t485 + t504 * t484 + Ifges(4,6) * t491 + Ifges(4,5) * t492 + mrSges(3,1) * t493 - mrSges(3,2) * t494 + mrSges(4,1) * t444 - mrSges(4,2) * t445 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t529 * t506 - t532 * t507) * qJD(1); t387; t536; t406; t535;];
tauJ  = t1;
