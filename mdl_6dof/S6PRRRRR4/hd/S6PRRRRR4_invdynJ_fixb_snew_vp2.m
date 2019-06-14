% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 11:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:25:33
% EndTime: 2019-05-05 11:25:40
% DurationCPUTime: 6.99s
% Computational Cost: add. (84063->302), mult. (173582->396), div. (0->0), fcn. (141466->16), ass. (0->135)
t500 = sin(pkin(13));
t503 = cos(pkin(13));
t491 = g(1) * t500 - g(2) * t503;
t499 = -g(3) + qJDD(1);
t502 = sin(pkin(6));
t505 = cos(pkin(6));
t539 = t491 * t505 + t499 * t502;
t492 = -g(1) * t503 - g(2) * t500;
t510 = sin(qJ(2));
t515 = cos(qJ(2));
t461 = -t492 * t510 + t515 * t539;
t501 = sin(pkin(7));
t516 = qJD(2) ^ 2;
t458 = pkin(9) * t501 * t516 + qJDD(2) * pkin(2) + t461;
t462 = t515 * t492 + t510 * t539;
t530 = qJDD(2) * t501;
t459 = -pkin(2) * t516 + pkin(9) * t530 + t462;
t504 = cos(pkin(7));
t509 = sin(qJ(3));
t514 = cos(qJ(3));
t476 = -t491 * t502 + t499 * t505;
t537 = t476 * t501;
t429 = -t509 * t459 + t514 * (t458 * t504 + t537);
t498 = qJD(2) * t504 + qJD(3);
t531 = qJD(2) * t514;
t528 = t501 * t531;
t480 = -mrSges(4,2) * t498 + mrSges(4,3) * t528;
t532 = qJD(2) * t501;
t481 = (-mrSges(4,1) * t514 + mrSges(4,2) * t509) * t532;
t483 = (qJD(3) * t531 + qJDD(2) * t509) * t501;
t497 = qJDD(2) * t504 + qJDD(3);
t482 = (-pkin(3) * t514 - pkin(10) * t509) * t532;
t496 = t498 ^ 2;
t529 = t509 * t532;
t416 = -pkin(3) * t497 - pkin(10) * t496 + t482 * t529 - t429;
t508 = sin(qJ(4));
t513 = cos(qJ(4));
t475 = t498 * t508 + t513 * t529;
t450 = -qJD(4) * t475 - t483 * t508 + t497 * t513;
t474 = t498 * t513 - t508 * t529;
t451 = qJD(4) * t474 + t483 * t513 + t497 * t508;
t490 = qJD(4) - t528;
t463 = -mrSges(5,2) * t490 + mrSges(5,3) * t474;
t464 = mrSges(5,1) * t490 - mrSges(5,3) * t475;
t534 = t504 * t509;
t430 = t458 * t534 + t514 * t459 + t509 * t537;
t417 = -pkin(3) * t496 + pkin(10) * t497 + t482 * t528 + t430;
t471 = t504 * t476;
t484 = -qJD(3) * t529 + t514 * t530;
t427 = -pkin(3) * t484 - pkin(10) * t483 + t471 + (-t458 + (pkin(3) * t509 - pkin(10) * t514) * t498 * qJD(2)) * t501;
t405 = -t417 * t508 + t513 * t427;
t478 = qJDD(4) - t484;
t402 = (t474 * t490 - t451) * pkin(11) + (t474 * t475 + t478) * pkin(4) + t405;
t406 = t513 * t417 + t508 * t427;
t465 = pkin(4) * t490 - pkin(11) * t475;
t473 = t474 ^ 2;
t404 = -pkin(4) * t473 + pkin(11) * t450 - t465 * t490 + t406;
t507 = sin(qJ(5));
t512 = cos(qJ(5));
t399 = t507 * t402 + t512 * t404;
t456 = t474 * t512 - t475 * t507;
t457 = t474 * t507 + t475 * t512;
t438 = -pkin(5) * t456 - pkin(12) * t457;
t477 = qJDD(5) + t478;
t489 = qJD(5) + t490;
t488 = t489 ^ 2;
t396 = -pkin(5) * t488 + pkin(12) * t477 + t438 * t456 + t399;
t407 = -pkin(4) * t450 - pkin(11) * t473 + t475 * t465 + t416;
t422 = -qJD(5) * t457 + t450 * t512 - t451 * t507;
t423 = qJD(5) * t456 + t450 * t507 + t451 * t512;
t400 = (-t456 * t489 - t423) * pkin(12) + (t457 * t489 - t422) * pkin(5) + t407;
t506 = sin(qJ(6));
t511 = cos(qJ(6));
t393 = -t396 * t506 + t400 * t511;
t440 = -t457 * t506 + t489 * t511;
t410 = qJD(6) * t440 + t423 * t511 + t477 * t506;
t421 = qJDD(6) - t422;
t441 = t457 * t511 + t489 * t506;
t428 = -mrSges(7,1) * t440 + mrSges(7,2) * t441;
t453 = qJD(6) - t456;
t431 = -mrSges(7,2) * t453 + mrSges(7,3) * t440;
t389 = m(7) * t393 + mrSges(7,1) * t421 - mrSges(7,3) * t410 - t428 * t441 + t431 * t453;
t394 = t396 * t511 + t400 * t506;
t409 = -qJD(6) * t441 - t423 * t506 + t477 * t511;
t432 = mrSges(7,1) * t453 - mrSges(7,3) * t441;
t390 = m(7) * t394 - mrSges(7,2) * t421 + mrSges(7,3) * t409 + t428 * t440 - t432 * t453;
t378 = t511 * t389 + t506 * t390;
t442 = -mrSges(6,2) * t489 + mrSges(6,3) * t456;
t443 = mrSges(6,1) * t489 - mrSges(6,3) * t457;
t521 = m(6) * t407 - t422 * mrSges(6,1) + mrSges(6,2) * t423 - t456 * t442 + t443 * t457 + t378;
t518 = -m(5) * t416 + t450 * mrSges(5,1) - mrSges(5,2) * t451 + t474 * t463 - t464 * t475 - t521;
t373 = m(4) * t429 + mrSges(4,1) * t497 - mrSges(4,3) * t483 + t480 * t498 - t481 * t529 + t518;
t538 = t373 * t514;
t479 = mrSges(4,1) * t498 - mrSges(4,3) * t529;
t437 = -mrSges(6,1) * t456 + mrSges(6,2) * t457;
t524 = -t389 * t506 + t511 * t390;
t376 = m(6) * t399 - mrSges(6,2) * t477 + mrSges(6,3) * t422 + t437 * t456 - t443 * t489 + t524;
t398 = t402 * t512 - t404 * t507;
t395 = -pkin(5) * t477 - pkin(12) * t488 + t438 * t457 - t398;
t522 = -m(7) * t395 + t409 * mrSges(7,1) - mrSges(7,2) * t410 + t440 * t431 - t432 * t441;
t385 = m(6) * t398 + mrSges(6,1) * t477 - mrSges(6,3) * t423 - t437 * t457 + t442 * t489 + t522;
t371 = t507 * t376 + t512 * t385;
t460 = -mrSges(5,1) * t474 + mrSges(5,2) * t475;
t369 = m(5) * t405 + mrSges(5,1) * t478 - mrSges(5,3) * t451 - t460 * t475 + t463 * t490 + t371;
t525 = t512 * t376 - t385 * t507;
t370 = m(5) * t406 - mrSges(5,2) * t478 + mrSges(5,3) * t450 + t460 * t474 - t464 * t490 + t525;
t526 = -t369 * t508 + t513 * t370;
t361 = m(4) * t430 - mrSges(4,2) * t497 + mrSges(4,3) * t484 - t479 * t498 + t481 * t528 + t526;
t533 = t361 * t534 + t504 * t538;
t363 = t513 * t369 + t508 * t370;
t527 = t514 * t361 - t373 * t509;
t411 = Ifges(7,5) * t441 + Ifges(7,6) * t440 + Ifges(7,3) * t453;
t413 = Ifges(7,1) * t441 + Ifges(7,4) * t440 + Ifges(7,5) * t453;
t382 = -mrSges(7,1) * t395 + mrSges(7,3) * t394 + Ifges(7,4) * t410 + Ifges(7,2) * t409 + Ifges(7,6) * t421 - t411 * t441 + t413 * t453;
t412 = Ifges(7,4) * t441 + Ifges(7,2) * t440 + Ifges(7,6) * t453;
t383 = mrSges(7,2) * t395 - mrSges(7,3) * t393 + Ifges(7,1) * t410 + Ifges(7,4) * t409 + Ifges(7,5) * t421 + t411 * t440 - t412 * t453;
t434 = Ifges(6,4) * t457 + Ifges(6,2) * t456 + Ifges(6,6) * t489;
t435 = Ifges(6,1) * t457 + Ifges(6,4) * t456 + Ifges(6,5) * t489;
t520 = -mrSges(6,1) * t398 + mrSges(6,2) * t399 - Ifges(6,5) * t423 - Ifges(6,6) * t422 - Ifges(6,3) * t477 - pkin(5) * t522 - pkin(12) * t524 - t511 * t382 - t506 * t383 - t457 * t434 + t456 * t435;
t519 = mrSges(7,1) * t393 - mrSges(7,2) * t394 + Ifges(7,5) * t410 + Ifges(7,6) * t409 + Ifges(7,3) * t421 + t412 * t441 - t413 * t440;
t445 = Ifges(5,4) * t475 + Ifges(5,2) * t474 + Ifges(5,6) * t490;
t446 = Ifges(5,1) * t475 + Ifges(5,4) * t474 + Ifges(5,5) * t490;
t517 = mrSges(5,1) * t405 - mrSges(5,2) * t406 + Ifges(5,5) * t451 + Ifges(5,6) * t450 + Ifges(5,3) * t478 + pkin(4) * t371 + t475 * t445 - t474 * t446 - t520;
t469 = Ifges(4,5) * t498 + (Ifges(4,1) * t509 + Ifges(4,4) * t514) * t532;
t468 = Ifges(4,6) * t498 + (Ifges(4,4) * t509 + Ifges(4,2) * t514) * t532;
t444 = Ifges(5,5) * t475 + Ifges(5,6) * t474 + Ifges(5,3) * t490;
t439 = -t458 * t501 + t471;
t433 = Ifges(6,5) * t457 + Ifges(6,6) * t456 + Ifges(6,3) * t489;
t365 = -mrSges(6,1) * t407 + mrSges(6,3) * t399 + Ifges(6,4) * t423 + Ifges(6,2) * t422 + Ifges(6,6) * t477 - pkin(5) * t378 - t433 * t457 + t435 * t489 - t519;
t364 = mrSges(6,2) * t407 - mrSges(6,3) * t398 + Ifges(6,1) * t423 + Ifges(6,4) * t422 + Ifges(6,5) * t477 - pkin(12) * t378 - t382 * t506 + t383 * t511 + t433 * t456 - t434 * t489;
t362 = m(4) * t439 - mrSges(4,1) * t484 + mrSges(4,2) * t483 + (t479 * t509 - t480 * t514) * t532 + t363;
t358 = mrSges(5,2) * t416 - mrSges(5,3) * t405 + Ifges(5,1) * t451 + Ifges(5,4) * t450 + Ifges(5,5) * t478 - pkin(11) * t371 + t364 * t512 - t365 * t507 + t444 * t474 - t445 * t490;
t357 = -mrSges(5,1) * t416 + mrSges(5,3) * t406 + Ifges(5,4) * t451 + Ifges(5,2) * t450 + Ifges(5,6) * t478 - pkin(4) * t521 + pkin(11) * t525 + t507 * t364 + t512 * t365 - t475 * t444 + t490 * t446;
t356 = Ifges(4,5) * t483 + Ifges(4,6) * t484 + Ifges(4,3) * t497 + mrSges(4,1) * t429 - mrSges(4,2) * t430 + t508 * t358 + t513 * t357 + pkin(3) * t518 + pkin(10) * t526 + (t468 * t509 - t469 * t514) * t532;
t1 = [m(2) * t499 + t505 * (m(3) * t476 + t362 * t504 + (t361 * t509 + t538) * t501) + (t510 * (m(3) * t462 - mrSges(3,1) * t516 - qJDD(2) * mrSges(3,2) + t527) + t515 * (m(3) * t461 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t516 - t362 * t501 + t533)) * t502; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t461 - mrSges(3,2) * t462 + t504 * t356 + pkin(2) * t533 + (t509 * (mrSges(4,2) * t439 - mrSges(4,3) * t429 + Ifges(4,1) * t483 + Ifges(4,4) * t484 + Ifges(4,5) * t497 - pkin(10) * t363 - t357 * t508 + t358 * t513 - t468 * t498) + t514 * (-mrSges(4,1) * t439 + mrSges(4,3) * t430 + Ifges(4,4) * t483 + Ifges(4,2) * t484 + Ifges(4,6) * t497 - pkin(3) * t363 + t498 * t469 - t517) - pkin(2) * t362 + pkin(9) * t527) * t501; t356; t517; -t520; t519;];
tauJ  = t1;
