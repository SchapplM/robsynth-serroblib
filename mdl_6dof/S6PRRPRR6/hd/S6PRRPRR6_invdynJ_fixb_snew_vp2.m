% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR6
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
% Datum: 2019-05-05 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:41:13
% EndTime: 2019-05-05 05:41:20
% DurationCPUTime: 6.48s
% Computational Cost: add. (72719->300), mult. (157786->397), div. (0->0), fcn. (127962->16), ass. (0->132)
t493 = sin(pkin(12));
t497 = cos(pkin(12));
t483 = t493 * g(1) - t497 * g(2);
t491 = -g(3) + qJDD(1);
t495 = sin(pkin(6));
t499 = cos(pkin(6));
t527 = t483 * t499 + t491 * t495;
t484 = -t497 * g(1) - t493 * g(2);
t503 = sin(qJ(2));
t507 = cos(qJ(2));
t451 = -t503 * t484 + t527 * t507;
t494 = sin(pkin(7));
t508 = qJD(2) ^ 2;
t449 = t508 * t494 * pkin(9) + qJDD(2) * pkin(2) + t451;
t452 = t507 * t484 + t527 * t503;
t518 = qJDD(2) * t494;
t450 = -t508 * pkin(2) + pkin(9) * t518 + t452;
t498 = cos(pkin(7));
t502 = sin(qJ(3));
t506 = cos(qJ(3));
t470 = -t495 * t483 + t499 * t491;
t525 = t470 * t494;
t418 = -t502 * t450 + (t449 * t498 + t525) * t506;
t520 = qJD(2) * t494;
t474 = (-t506 * pkin(3) - t502 * qJ(4)) * t520;
t490 = t498 * qJD(2) + qJD(3);
t488 = t490 ^ 2;
t489 = t498 * qJDD(2) + qJDD(3);
t517 = t502 * t520;
t412 = -t489 * pkin(3) - t488 * qJ(4) + t474 * t517 + qJDD(4) - t418;
t492 = sin(pkin(13));
t496 = cos(pkin(13));
t468 = t496 * t490 - t492 * t517;
t519 = qJD(2) * t506;
t516 = t494 * t519;
t454 = mrSges(5,2) * t516 + t468 * mrSges(5,3);
t469 = t492 * t490 + t496 * t517;
t455 = -mrSges(5,1) * t516 - t469 * mrSges(5,3);
t476 = (qJD(3) * t519 + qJDD(2) * t502) * t494;
t456 = -t492 * t476 + t496 * t489;
t457 = t496 * t476 + t492 * t489;
t522 = t498 * t502;
t419 = t449 * t522 + t506 * t450 + t502 * t525;
t413 = -t488 * pkin(3) + t489 * qJ(4) + t474 * t516 + t419;
t466 = t498 * t470;
t477 = -qJD(3) * t517 + t506 * t518;
t416 = -t477 * pkin(3) - t476 * qJ(4) + t466 + (-t449 + (pkin(3) * t502 - qJ(4) * t506) * t490 * qJD(2)) * t494;
t401 = -0.2e1 * qJD(4) * t469 - t492 * t413 + t496 * t416;
t398 = (-t468 * t516 - t457) * pkin(10) + (t468 * t469 - t477) * pkin(4) + t401;
t402 = 0.2e1 * qJD(4) * t468 + t496 * t413 + t492 * t416;
t458 = -pkin(4) * t516 - t469 * pkin(10);
t467 = t468 ^ 2;
t400 = -t467 * pkin(4) + t456 * pkin(10) + t458 * t516 + t402;
t501 = sin(qJ(5));
t505 = cos(qJ(5));
t395 = t501 * t398 + t505 * t400;
t444 = t505 * t468 - t501 * t469;
t445 = t501 * t468 + t505 * t469;
t431 = -t444 * pkin(5) - t445 * pkin(11);
t471 = qJDD(5) - t477;
t482 = qJD(5) - t516;
t481 = t482 ^ 2;
t393 = -t481 * pkin(5) + t471 * pkin(11) + t444 * t431 + t395;
t403 = -t456 * pkin(4) - t467 * pkin(10) + t469 * t458 + t412;
t422 = -t445 * qJD(5) + t505 * t456 - t501 * t457;
t423 = t444 * qJD(5) + t501 * t456 + t505 * t457;
t396 = (-t444 * t482 - t423) * pkin(11) + (t445 * t482 - t422) * pkin(5) + t403;
t500 = sin(qJ(6));
t504 = cos(qJ(6));
t390 = -t500 * t393 + t504 * t396;
t433 = -t500 * t445 + t504 * t482;
t406 = t433 * qJD(6) + t504 * t423 + t500 * t471;
t434 = t504 * t445 + t500 * t482;
t417 = -t433 * mrSges(7,1) + t434 * mrSges(7,2);
t421 = qJDD(6) - t422;
t443 = qJD(6) - t444;
t424 = -t443 * mrSges(7,2) + t433 * mrSges(7,3);
t387 = m(7) * t390 + t421 * mrSges(7,1) - t406 * mrSges(7,3) - t434 * t417 + t443 * t424;
t391 = t504 * t393 + t500 * t396;
t405 = -t434 * qJD(6) - t500 * t423 + t504 * t471;
t425 = t443 * mrSges(7,1) - t434 * mrSges(7,3);
t388 = m(7) * t391 - t421 * mrSges(7,2) + t405 * mrSges(7,3) + t433 * t417 - t443 * t425;
t378 = t504 * t387 + t500 * t388;
t435 = -t482 * mrSges(6,2) + t444 * mrSges(6,3);
t436 = t482 * mrSges(6,1) - t445 * mrSges(6,3);
t511 = m(6) * t403 - t422 * mrSges(6,1) + t423 * mrSges(6,2) - t444 * t435 + t445 * t436 + t378;
t377 = m(5) * t412 - t456 * mrSges(5,1) + t457 * mrSges(5,2) - t468 * t454 + t469 * t455 + t511;
t473 = -t490 * mrSges(4,2) + mrSges(4,3) * t516;
t475 = (-t506 * mrSges(4,1) + t502 * mrSges(4,2)) * t520;
t373 = m(4) * t418 + t489 * mrSges(4,1) - t476 * mrSges(4,3) + t490 * t473 - t475 * t517 - t377;
t526 = t373 * t506;
t472 = t490 * mrSges(4,1) - mrSges(4,3) * t517;
t379 = -t500 * t387 + t504 * t388;
t430 = -t444 * mrSges(6,1) + t445 * mrSges(6,2);
t376 = m(6) * t395 - t471 * mrSges(6,2) + t422 * mrSges(6,3) + t444 * t430 - t482 * t436 + t379;
t394 = t505 * t398 - t501 * t400;
t392 = -t471 * pkin(5) - t481 * pkin(11) + t445 * t431 - t394;
t389 = -m(7) * t392 + t405 * mrSges(7,1) - t406 * mrSges(7,2) + t433 * t424 - t434 * t425;
t383 = m(6) * t394 + t471 * mrSges(6,1) - t423 * mrSges(6,3) - t445 * t430 + t482 * t435 + t389;
t371 = t501 * t376 + t505 * t383;
t448 = -t468 * mrSges(5,1) + t469 * mrSges(5,2);
t369 = m(5) * t401 - t477 * mrSges(5,1) - t457 * mrSges(5,3) - t469 * t448 - t454 * t516 + t371;
t513 = t505 * t376 - t501 * t383;
t370 = m(5) * t402 + t477 * mrSges(5,2) + t456 * mrSges(5,3) + t468 * t448 + t455 * t516 + t513;
t514 = -t492 * t369 + t496 * t370;
t361 = m(4) * t419 - t489 * mrSges(4,2) + t477 * mrSges(4,3) - t490 * t472 + t475 * t516 + t514;
t521 = t361 * t522 + t498 * t526;
t363 = t496 * t369 + t492 * t370;
t515 = t506 * t361 - t502 * t373;
t408 = Ifges(7,4) * t434 + Ifges(7,2) * t433 + Ifges(7,6) * t443;
t409 = Ifges(7,1) * t434 + Ifges(7,4) * t433 + Ifges(7,5) * t443;
t510 = mrSges(7,1) * t390 - mrSges(7,2) * t391 + Ifges(7,5) * t406 + Ifges(7,6) * t405 + Ifges(7,3) * t421 + t434 * t408 - t433 * t409;
t407 = Ifges(7,5) * t434 + Ifges(7,6) * t433 + Ifges(7,3) * t443;
t380 = -mrSges(7,1) * t392 + mrSges(7,3) * t391 + Ifges(7,4) * t406 + Ifges(7,2) * t405 + Ifges(7,6) * t421 - t434 * t407 + t443 * t409;
t381 = mrSges(7,2) * t392 - mrSges(7,3) * t390 + Ifges(7,1) * t406 + Ifges(7,4) * t405 + Ifges(7,5) * t421 + t433 * t407 - t443 * t408;
t427 = Ifges(6,4) * t445 + Ifges(6,2) * t444 + Ifges(6,6) * t482;
t428 = Ifges(6,1) * t445 + Ifges(6,4) * t444 + Ifges(6,5) * t482;
t509 = mrSges(6,1) * t394 - mrSges(6,2) * t395 + Ifges(6,5) * t423 + Ifges(6,6) * t422 + Ifges(6,3) * t471 + pkin(5) * t389 + pkin(11) * t379 + t504 * t380 + t500 * t381 + t445 * t427 - t444 * t428;
t462 = Ifges(4,5) * t490 + (t502 * Ifges(4,1) + t506 * Ifges(4,4)) * t520;
t461 = Ifges(4,6) * t490 + (t502 * Ifges(4,4) + Ifges(4,2) * t506) * t520;
t439 = Ifges(5,1) * t469 + Ifges(5,4) * t468 - Ifges(5,5) * t516;
t438 = Ifges(5,4) * t469 + Ifges(5,2) * t468 - Ifges(5,6) * t516;
t437 = Ifges(5,5) * t469 + Ifges(5,6) * t468 - Ifges(5,3) * t516;
t432 = -t494 * t449 + t466;
t426 = Ifges(6,5) * t445 + Ifges(6,6) * t444 + Ifges(6,3) * t482;
t365 = -mrSges(6,1) * t403 + mrSges(6,3) * t395 + Ifges(6,4) * t423 + Ifges(6,2) * t422 + Ifges(6,6) * t471 - pkin(5) * t378 - t445 * t426 + t482 * t428 - t510;
t364 = mrSges(6,2) * t403 - mrSges(6,3) * t394 + Ifges(6,1) * t423 + Ifges(6,4) * t422 + Ifges(6,5) * t471 - pkin(11) * t378 - t500 * t380 + t504 * t381 + t444 * t426 - t482 * t427;
t362 = m(4) * t432 - t477 * mrSges(4,1) + t476 * mrSges(4,2) + (t472 * t502 - t473 * t506) * t520 + t363;
t358 = mrSges(5,2) * t412 - mrSges(5,3) * t401 + Ifges(5,1) * t457 + Ifges(5,4) * t456 - Ifges(5,5) * t477 - pkin(10) * t371 + t505 * t364 - t501 * t365 + t468 * t437 + t438 * t516;
t357 = -mrSges(5,1) * t412 + mrSges(5,3) * t402 + Ifges(5,4) * t457 + Ifges(5,2) * t456 - Ifges(5,6) * t477 - pkin(4) * t511 + pkin(10) * t513 + t501 * t364 + t505 * t365 - t469 * t437 - t439 * t516;
t356 = Ifges(4,5) * t476 + Ifges(4,6) * t477 + Ifges(4,3) * t489 + mrSges(4,1) * t418 - mrSges(4,2) * t419 + t492 * t358 + t496 * t357 - pkin(3) * t377 + qJ(4) * t514 + (t502 * t461 - t506 * t462) * t520;
t1 = [m(2) * t491 + t499 * (m(3) * t470 + t498 * t362 + (t361 * t502 + t526) * t494) + (t503 * (m(3) * t452 - t508 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t515) + t507 * (m(3) * t451 + qJDD(2) * mrSges(3,1) - t508 * mrSges(3,2) - t494 * t362 + t521)) * t495; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t451 - mrSges(3,2) * t452 + t498 * t356 + pkin(2) * t521 + (t502 * (mrSges(4,2) * t432 - mrSges(4,3) * t418 + Ifges(4,1) * t476 + Ifges(4,4) * t477 + Ifges(4,5) * t489 - qJ(4) * t363 - t492 * t357 + t496 * t358 - t490 * t461) + t506 * (-mrSges(4,1) * t432 - mrSges(5,1) * t401 + mrSges(5,2) * t402 + mrSges(4,3) * t419 + Ifges(4,4) * t476 - Ifges(5,5) * t457 + Ifges(4,6) * t489 - Ifges(5,6) * t456 - pkin(3) * t363 - pkin(4) * t371 - t469 * t438 + t468 * t439 + t490 * t462 - t509 + (Ifges(4,2) + Ifges(5,3)) * t477) - pkin(2) * t362 + pkin(9) * t515) * t494; t356; t377; t509; t510;];
tauJ  = t1;
