% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 02:47
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:44:28
% EndTime: 2019-05-06 02:44:32
% DurationCPUTime: 3.62s
% Computational Cost: add. (45031->294), mult. (87974->373), div. (0->0), fcn. (60897->12), ass. (0->122)
t483 = sin(qJ(1));
t488 = cos(qJ(1));
t503 = t483 * g(1) - g(2) * t488;
t459 = qJDD(1) * pkin(1) + t503;
t489 = qJD(1) ^ 2;
t498 = -g(1) * t488 - g(2) * t483;
t461 = -pkin(1) * t489 + t498;
t477 = sin(pkin(11));
t478 = cos(pkin(11));
t438 = t477 * t459 + t478 * t461;
t434 = -pkin(2) * t489 + qJDD(1) * pkin(7) + t438;
t476 = -g(3) + qJDD(2);
t482 = sin(qJ(3));
t487 = cos(qJ(3));
t421 = -t434 * t482 + t487 * t476;
t505 = qJD(1) * qJD(3);
t504 = t487 * t505;
t462 = qJDD(1) * t482 + t504;
t403 = (-t462 + t504) * pkin(8) + (t482 * t487 * t489 + qJDD(3)) * pkin(3) + t421;
t422 = t487 * t434 + t482 * t476;
t463 = qJDD(1) * t487 - t482 * t505;
t507 = qJD(1) * t482;
t466 = qJD(3) * pkin(3) - pkin(8) * t507;
t475 = t487 ^ 2;
t404 = -pkin(3) * t475 * t489 + pkin(8) * t463 - qJD(3) * t466 + t422;
t481 = sin(qJ(4));
t486 = cos(qJ(4));
t389 = t481 * t403 + t486 * t404;
t454 = (t487 * t481 + t482 * t486) * qJD(1);
t423 = -t454 * qJD(4) - t481 * t462 + t463 * t486;
t506 = qJD(1) * t487;
t453 = -t481 * t507 + t486 * t506;
t435 = -mrSges(5,1) * t453 + mrSges(5,2) * t454;
t474 = qJD(3) + qJD(4);
t443 = mrSges(5,1) * t474 - mrSges(5,3) * t454;
t473 = qJDD(3) + qJDD(4);
t436 = -pkin(4) * t453 - pkin(9) * t454;
t472 = t474 ^ 2;
t382 = -pkin(4) * t472 + pkin(9) * t473 + t436 * t453 + t389;
t437 = t459 * t478 - t477 * t461;
t497 = -qJDD(1) * pkin(2) - t437;
t409 = -pkin(3) * t463 + t466 * t507 + (-pkin(8) * t475 - pkin(7)) * t489 + t497;
t424 = qJD(4) * t453 + t462 * t486 + t463 * t481;
t385 = (-t453 * t474 - t424) * pkin(9) + (t454 * t474 - t423) * pkin(4) + t409;
t480 = sin(qJ(5));
t485 = cos(qJ(5));
t372 = -t382 * t480 + t485 * t385;
t440 = -t454 * t480 + t474 * t485;
t397 = qJD(5) * t440 + t424 * t485 + t473 * t480;
t420 = qJDD(5) - t423;
t441 = t454 * t485 + t474 * t480;
t446 = qJD(5) - t453;
t370 = (t440 * t446 - t397) * pkin(10) + (t440 * t441 + t420) * pkin(5) + t372;
t373 = t485 * t382 + t480 * t385;
t396 = -qJD(5) * t441 - t424 * t480 + t473 * t485;
t427 = pkin(5) * t446 - pkin(10) * t441;
t439 = t440 ^ 2;
t371 = -pkin(5) * t439 + pkin(10) * t396 - t427 * t446 + t373;
t479 = sin(qJ(6));
t484 = cos(qJ(6));
t368 = t370 * t484 - t371 * t479;
t410 = t440 * t484 - t441 * t479;
t379 = qJD(6) * t410 + t396 * t479 + t397 * t484;
t411 = t440 * t479 + t441 * t484;
t394 = -mrSges(7,1) * t410 + mrSges(7,2) * t411;
t444 = qJD(6) + t446;
t401 = -mrSges(7,2) * t444 + mrSges(7,3) * t410;
t415 = qJDD(6) + t420;
t364 = m(7) * t368 + mrSges(7,1) * t415 - mrSges(7,3) * t379 - t394 * t411 + t401 * t444;
t369 = t370 * t479 + t371 * t484;
t378 = -qJD(6) * t411 + t396 * t484 - t397 * t479;
t402 = mrSges(7,1) * t444 - mrSges(7,3) * t411;
t365 = m(7) * t369 - mrSges(7,2) * t415 + mrSges(7,3) * t378 + t394 * t410 - t402 * t444;
t356 = t484 * t364 + t479 * t365;
t412 = -mrSges(6,1) * t440 + mrSges(6,2) * t441;
t425 = -mrSges(6,2) * t446 + mrSges(6,3) * t440;
t354 = m(6) * t372 + mrSges(6,1) * t420 - mrSges(6,3) * t397 - t412 * t441 + t425 * t446 + t356;
t426 = mrSges(6,1) * t446 - mrSges(6,3) * t441;
t499 = -t364 * t479 + t484 * t365;
t355 = m(6) * t373 - mrSges(6,2) * t420 + mrSges(6,3) * t396 + t412 * t440 - t426 * t446 + t499;
t500 = -t354 * t480 + t485 * t355;
t348 = m(5) * t389 - mrSges(5,2) * t473 + mrSges(5,3) * t423 + t435 * t453 - t443 * t474 + t500;
t388 = t403 * t486 - t481 * t404;
t442 = -mrSges(5,2) * t474 + mrSges(5,3) * t453;
t381 = -pkin(4) * t473 - pkin(9) * t472 + t454 * t436 - t388;
t374 = -pkin(5) * t396 - pkin(10) * t439 + t427 * t441 + t381;
t496 = m(7) * t374 - t378 * mrSges(7,1) + mrSges(7,2) * t379 - t410 * t401 + t402 * t411;
t492 = -m(6) * t381 + t396 * mrSges(6,1) - mrSges(6,2) * t397 + t440 * t425 - t426 * t441 - t496;
t360 = m(5) * t388 + mrSges(5,1) * t473 - mrSges(5,3) * t424 - t435 * t454 + t442 * t474 + t492;
t341 = t481 * t348 + t486 * t360;
t350 = t485 * t354 + t480 * t355;
t460 = (-t487 * mrSges(4,1) + t482 * mrSges(4,2)) * qJD(1);
t465 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t506;
t339 = m(4) * t421 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t462 + qJD(3) * t465 - t460 * t507 + t341;
t464 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t507;
t501 = t486 * t348 - t360 * t481;
t340 = m(4) * t422 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t463 - qJD(3) * t464 + t460 * t506 + t501;
t502 = -t339 * t482 + t487 * t340;
t391 = Ifges(7,4) * t411 + Ifges(7,2) * t410 + Ifges(7,6) * t444;
t392 = Ifges(7,1) * t411 + Ifges(7,4) * t410 + Ifges(7,5) * t444;
t495 = -mrSges(7,1) * t368 + mrSges(7,2) * t369 - Ifges(7,5) * t379 - Ifges(7,6) * t378 - Ifges(7,3) * t415 - t411 * t391 + t410 * t392;
t494 = m(5) * t409 - t423 * mrSges(5,1) + t424 * mrSges(5,2) - t453 * t442 + t454 * t443 + t350;
t390 = Ifges(7,5) * t411 + Ifges(7,6) * t410 + Ifges(7,3) * t444;
t357 = -mrSges(7,1) * t374 + mrSges(7,3) * t369 + Ifges(7,4) * t379 + Ifges(7,2) * t378 + Ifges(7,6) * t415 - t390 * t411 + t392 * t444;
t358 = mrSges(7,2) * t374 - mrSges(7,3) * t368 + Ifges(7,1) * t379 + Ifges(7,4) * t378 + Ifges(7,5) * t415 + t390 * t410 - t391 * t444;
t405 = Ifges(6,5) * t441 + Ifges(6,6) * t440 + Ifges(6,3) * t446;
t407 = Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * t446;
t343 = -mrSges(6,1) * t381 + mrSges(6,3) * t373 + Ifges(6,4) * t397 + Ifges(6,2) * t396 + Ifges(6,6) * t420 - pkin(5) * t496 + pkin(10) * t499 + t484 * t357 + t479 * t358 - t441 * t405 + t446 * t407;
t406 = Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * t446;
t345 = mrSges(6,2) * t381 - mrSges(6,3) * t372 + Ifges(6,1) * t397 + Ifges(6,4) * t396 + Ifges(6,5) * t420 - pkin(10) * t356 - t357 * t479 + t358 * t484 + t405 * t440 - t406 * t446;
t430 = Ifges(5,4) * t454 + Ifges(5,2) * t453 + Ifges(5,6) * t474;
t431 = Ifges(5,1) * t454 + Ifges(5,4) * t453 + Ifges(5,5) * t474;
t493 = mrSges(5,1) * t388 - mrSges(5,2) * t389 + Ifges(5,5) * t424 + Ifges(5,6) * t423 + Ifges(5,3) * t473 + pkin(4) * t492 + pkin(9) * t500 + t485 * t343 + t480 * t345 + t454 * t430 - t431 * t453;
t433 = -pkin(7) * t489 + t497;
t491 = -m(4) * t433 + t463 * mrSges(4,1) - mrSges(4,2) * t462 - t464 * t507 + t465 * t506 - t494;
t490 = mrSges(6,1) * t372 - mrSges(6,2) * t373 + Ifges(6,5) * t397 + Ifges(6,6) * t396 + Ifges(6,3) * t420 + pkin(5) * t356 + t441 * t406 - t440 * t407 - t495;
t452 = Ifges(4,5) * qJD(3) + (t482 * Ifges(4,1) + t487 * Ifges(4,4)) * qJD(1);
t451 = Ifges(4,6) * qJD(3) + (t482 * Ifges(4,4) + t487 * Ifges(4,2)) * qJD(1);
t429 = Ifges(5,5) * t454 + Ifges(5,6) * t453 + Ifges(5,3) * t474;
t337 = -mrSges(5,1) * t409 + mrSges(5,3) * t389 + Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * t473 - pkin(4) * t350 - t454 * t429 + t474 * t431 - t490;
t336 = mrSges(5,2) * t409 - mrSges(5,3) * t388 + Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * t473 - pkin(9) * t350 - t343 * t480 + t345 * t485 + t429 * t453 - t430 * t474;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t503 - mrSges(2,2) * t498 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t437 - mrSges(3,2) * t438 + t482 * (mrSges(4,2) * t433 - mrSges(4,3) * t421 + Ifges(4,1) * t462 + Ifges(4,4) * t463 + Ifges(4,5) * qJDD(3) - pkin(8) * t341 - qJD(3) * t451 + t486 * t336 - t481 * t337) + t487 * (-mrSges(4,1) * t433 + mrSges(4,3) * t422 + Ifges(4,4) * t462 + Ifges(4,2) * t463 + Ifges(4,6) * qJDD(3) - pkin(3) * t494 + pkin(8) * t501 + qJD(3) * t452 + t481 * t336 + t486 * t337) + pkin(2) * t491 + pkin(7) * t502 + pkin(1) * (t477 * (m(3) * t438 - mrSges(3,1) * t489 - qJDD(1) * mrSges(3,2) + t502) + t478 * (m(3) * t437 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t489 + t491)); m(3) * t476 + t339 * t487 + t340 * t482; Ifges(4,3) * qJDD(3) + t493 + mrSges(4,1) * t421 - mrSges(4,2) * t422 + Ifges(4,5) * t462 + Ifges(4,6) * t463 + pkin(3) * t341 + (t451 * t482 - t452 * t487) * qJD(1); t493; t490; -t495;];
tauJ  = t1;
