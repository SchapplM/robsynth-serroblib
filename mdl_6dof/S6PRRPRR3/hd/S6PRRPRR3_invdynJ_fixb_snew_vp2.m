% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR3
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-05 04:54
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:46:22
% EndTime: 2019-05-05 04:46:28
% DurationCPUTime: 5.37s
% Computational Cost: add. (57952->299), mult. (140799->397), div. (0->0), fcn. (112265->16), ass. (0->133)
t541 = -2 * qJD(4);
t502 = sin(pkin(12));
t506 = cos(pkin(12));
t491 = g(1) * t502 - g(2) * t506;
t500 = -g(3) + qJDD(1);
t504 = sin(pkin(6));
t508 = cos(pkin(6));
t540 = t491 * t508 + t500 * t504;
t501 = sin(pkin(13));
t505 = cos(pkin(13));
t511 = sin(qJ(3));
t515 = cos(qJ(3));
t503 = sin(pkin(7));
t531 = qJD(2) * t503;
t479 = (t501 * t511 - t505 * t515) * t531;
t492 = -g(1) * t506 - g(2) * t502;
t512 = sin(qJ(2));
t516 = cos(qJ(2));
t461 = -t492 * t512 + t540 * t516;
t517 = qJD(2) ^ 2;
t539 = pkin(9) * t503;
t453 = qJDD(2) * pkin(2) + t517 * t539 + t461;
t462 = t516 * t492 + t540 * t512;
t454 = -pkin(2) * t517 + qJDD(2) * t539 + t462;
t507 = cos(pkin(7));
t533 = t507 * t515;
t481 = -t491 * t504 + t500 * t508;
t538 = t481 * t503;
t425 = t453 * t533 - t511 * t454 + t515 * t538;
t529 = qJD(2) * qJD(3);
t486 = (qJDD(2) * t511 + t515 * t529) * t503;
t496 = qJDD(2) * t507 + qJDD(3);
t497 = qJD(2) * t507 + qJD(3);
t526 = t515 * t531;
t536 = t503 ^ 2 * t517;
t417 = (t497 * t526 - t486) * qJ(4) + (t511 * t515 * t536 + t496) * pkin(3) + t425;
t534 = t507 * t511;
t426 = t453 * t534 + t515 * t454 + t511 * t538;
t527 = t511 * t531;
t482 = pkin(3) * t497 - qJ(4) * t527;
t487 = (qJDD(2) * t515 - t511 * t529) * t503;
t528 = t515 ^ 2 * t536;
t418 = -pkin(3) * t528 + qJ(4) * t487 - t482 * t497 + t426;
t480 = (t501 * t515 + t505 * t511) * t531;
t407 = t505 * t417 - t501 * t418 + t480 * t541;
t408 = t501 * t417 + t505 * t418 + t479 * t541;
t455 = mrSges(5,1) * t479 + mrSges(5,2) * t480;
t459 = -t486 * t501 + t487 * t505;
t467 = mrSges(5,1) * t497 - mrSges(5,3) * t480;
t456 = pkin(4) * t479 - pkin(10) * t480;
t495 = t497 ^ 2;
t406 = -pkin(4) * t495 + pkin(10) * t496 - t456 * t479 + t408;
t438 = -t503 * t453 + t507 * t481;
t427 = -t487 * pkin(3) - qJ(4) * t528 + t482 * t527 + qJDD(4) + t438;
t460 = t486 * t505 + t487 * t501;
t410 = (t479 * t497 - t460) * pkin(10) + (t480 * t497 - t459) * pkin(4) + t427;
t510 = sin(qJ(5));
t514 = cos(qJ(5));
t403 = t514 * t406 + t510 * t410;
t464 = -t480 * t510 + t497 * t514;
t465 = t480 * t514 + t497 * t510;
t440 = -pkin(5) * t464 - pkin(11) * t465;
t458 = qJDD(5) - t459;
t478 = qJD(5) + t479;
t477 = t478 ^ 2;
t400 = -pkin(5) * t477 + pkin(11) * t458 + t440 * t464 + t403;
t405 = -t496 * pkin(4) - t495 * pkin(10) + t480 * t456 - t407;
t436 = -qJD(5) * t465 - t460 * t510 + t496 * t514;
t437 = qJD(5) * t464 + t460 * t514 + t496 * t510;
t401 = (-t464 * t478 - t437) * pkin(11) + (t465 * t478 - t436) * pkin(5) + t405;
t509 = sin(qJ(6));
t513 = cos(qJ(6));
t396 = -t400 * t509 + t401 * t513;
t442 = -t465 * t509 + t478 * t513;
t413 = qJD(6) * t442 + t437 * t513 + t458 * t509;
t443 = t465 * t513 + t478 * t509;
t423 = -mrSges(7,1) * t442 + mrSges(7,2) * t443;
t463 = qJD(6) - t464;
t428 = -mrSges(7,2) * t463 + mrSges(7,3) * t442;
t435 = qJDD(6) - t436;
t394 = m(7) * t396 + mrSges(7,1) * t435 - mrSges(7,3) * t413 - t423 * t443 + t428 * t463;
t397 = t400 * t513 + t401 * t509;
t412 = -qJD(6) * t443 - t437 * t509 + t458 * t513;
t429 = mrSges(7,1) * t463 - mrSges(7,3) * t443;
t395 = m(7) * t397 - mrSges(7,2) * t435 + mrSges(7,3) * t412 + t423 * t442 - t429 * t463;
t388 = -t394 * t509 + t513 * t395;
t439 = -mrSges(6,1) * t464 + mrSges(6,2) * t465;
t445 = mrSges(6,1) * t478 - mrSges(6,3) * t465;
t386 = m(6) * t403 - mrSges(6,2) * t458 + mrSges(6,3) * t436 + t439 * t464 - t445 * t478 + t388;
t402 = -t406 * t510 + t410 * t514;
t399 = -pkin(5) * t458 - pkin(11) * t477 + t440 * t465 - t402;
t398 = -m(7) * t399 + t412 * mrSges(7,1) - mrSges(7,2) * t413 + t442 * t428 - t429 * t443;
t444 = -mrSges(6,2) * t478 + mrSges(6,3) * t464;
t392 = m(6) * t402 + mrSges(6,1) * t458 - mrSges(6,3) * t437 - t439 * t465 + t444 * t478 + t398;
t523 = t514 * t386 - t392 * t510;
t378 = m(5) * t408 - mrSges(5,2) * t496 + mrSges(5,3) * t459 - t455 * t479 - t467 * t497 + t523;
t466 = -mrSges(5,2) * t497 - mrSges(5,3) * t479;
t387 = t394 * t513 + t395 * t509;
t520 = -m(6) * t405 + t436 * mrSges(6,1) - mrSges(6,2) * t437 + t464 * t444 - t445 * t465 - t387;
t383 = m(5) * t407 + mrSges(5,1) * t496 - mrSges(5,3) * t460 - t455 * t480 + t466 * t497 + t520;
t373 = t501 * t378 + t505 * t383;
t484 = -mrSges(4,2) * t497 + mrSges(4,3) * t526;
t485 = (-mrSges(4,1) * t515 + mrSges(4,2) * t511) * t531;
t371 = m(4) * t425 + mrSges(4,1) * t496 - mrSges(4,3) * t486 + t484 * t497 - t485 * t527 + t373;
t483 = mrSges(4,1) * t497 - mrSges(4,3) * t527;
t524 = t505 * t378 - t383 * t501;
t372 = m(4) * t426 - mrSges(4,2) * t496 + mrSges(4,3) * t487 - t483 * t497 + t485 * t526 + t524;
t532 = t371 * t533 + t372 * t534;
t381 = t510 * t386 + t514 * t392;
t525 = -t371 * t511 + t515 * t372;
t380 = m(5) * t427 - t459 * mrSges(5,1) + t460 * mrSges(5,2) + t479 * t466 + t480 * t467 + t381;
t420 = Ifges(7,4) * t443 + Ifges(7,2) * t442 + Ifges(7,6) * t463;
t421 = Ifges(7,1) * t443 + Ifges(7,4) * t442 + Ifges(7,5) * t463;
t519 = mrSges(7,1) * t396 - mrSges(7,2) * t397 + Ifges(7,5) * t413 + Ifges(7,6) * t412 + Ifges(7,3) * t435 + t420 * t443 - t421 * t442;
t419 = Ifges(7,5) * t443 + Ifges(7,6) * t442 + Ifges(7,3) * t463;
t389 = -mrSges(7,1) * t399 + mrSges(7,3) * t397 + Ifges(7,4) * t413 + Ifges(7,2) * t412 + Ifges(7,6) * t435 - t419 * t443 + t421 * t463;
t390 = mrSges(7,2) * t399 - mrSges(7,3) * t396 + Ifges(7,1) * t413 + Ifges(7,4) * t412 + Ifges(7,5) * t435 + t419 * t442 - t420 * t463;
t431 = Ifges(6,4) * t465 + Ifges(6,2) * t464 + Ifges(6,6) * t478;
t432 = Ifges(6,1) * t465 + Ifges(6,4) * t464 + Ifges(6,5) * t478;
t518 = mrSges(6,1) * t402 - mrSges(6,2) * t403 + Ifges(6,5) * t437 + Ifges(6,6) * t436 + Ifges(6,3) * t458 + pkin(5) * t398 + pkin(11) * t388 + t513 * t389 + t509 * t390 + t465 * t431 - t464 * t432;
t472 = Ifges(4,5) * t497 + (Ifges(4,1) * t511 + Ifges(4,4) * t515) * t531;
t471 = Ifges(4,6) * t497 + (Ifges(4,4) * t511 + Ifges(4,2) * t515) * t531;
t448 = Ifges(5,1) * t480 - Ifges(5,4) * t479 + Ifges(5,5) * t497;
t447 = Ifges(5,4) * t480 - Ifges(5,2) * t479 + Ifges(5,6) * t497;
t446 = Ifges(5,5) * t480 - Ifges(5,6) * t479 + Ifges(5,3) * t497;
t430 = Ifges(6,5) * t465 + Ifges(6,6) * t464 + Ifges(6,3) * t478;
t379 = m(4) * t438 - t487 * mrSges(4,1) + t486 * mrSges(4,2) + (t483 * t511 - t484 * t515) * t531 + t380;
t375 = -mrSges(6,1) * t405 + mrSges(6,3) * t403 + Ifges(6,4) * t437 + Ifges(6,2) * t436 + Ifges(6,6) * t458 - pkin(5) * t387 - t430 * t465 + t432 * t478 - t519;
t374 = mrSges(6,2) * t405 - mrSges(6,3) * t402 + Ifges(6,1) * t437 + Ifges(6,4) * t436 + Ifges(6,5) * t458 - pkin(11) * t387 - t389 * t509 + t390 * t513 + t430 * t464 - t431 * t478;
t367 = -mrSges(5,1) * t427 + mrSges(5,3) * t408 + Ifges(5,4) * t460 + Ifges(5,2) * t459 + Ifges(5,6) * t496 - pkin(4) * t381 - t480 * t446 + t497 * t448 - t518;
t366 = mrSges(5,2) * t427 - mrSges(5,3) * t407 + Ifges(5,1) * t460 + Ifges(5,4) * t459 + Ifges(5,5) * t496 - pkin(10) * t381 + t374 * t514 - t375 * t510 - t446 * t479 - t447 * t497;
t365 = Ifges(4,5) * t486 + Ifges(4,6) * t487 + mrSges(4,1) * t425 - mrSges(4,2) * t426 + Ifges(5,5) * t460 + Ifges(5,6) * t459 + t480 * t447 + t479 * t448 + mrSges(5,1) * t407 - mrSges(5,2) * t408 + t510 * t374 + t514 * t375 + pkin(4) * t520 + pkin(10) * t523 + pkin(3) * t373 + (Ifges(4,3) + Ifges(5,3)) * t496 + (t471 * t511 - t472 * t515) * t531;
t1 = [m(2) * t500 + t508 * (m(3) * t481 + t507 * t379 + (t371 * t515 + t372 * t511) * t503) + (t512 * (m(3) * t462 - mrSges(3,1) * t517 - qJDD(2) * mrSges(3,2) + t525) + t516 * (m(3) * t461 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t517 - t379 * t503 + t532)) * t504; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t461 - mrSges(3,2) * t462 + t507 * t365 + pkin(2) * t532 + (t511 * (mrSges(4,2) * t438 - mrSges(4,3) * t425 + Ifges(4,1) * t486 + Ifges(4,4) * t487 + Ifges(4,5) * t496 - qJ(4) * t373 + t366 * t505 - t367 * t501 - t471 * t497) + t515 * (-mrSges(4,1) * t438 + mrSges(4,3) * t426 + Ifges(4,4) * t486 + Ifges(4,2) * t487 + Ifges(4,6) * t496 - pkin(3) * t380 + qJ(4) * t524 + t501 * t366 + t505 * t367 + t497 * t472) - pkin(2) * t379 + pkin(9) * t525) * t503; t365; t380; t518; t519;];
tauJ  = t1;
