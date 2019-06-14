% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-05-05 19:23
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:20:14
% EndTime: 2019-05-05 19:20:18
% DurationCPUTime: 3.25s
% Computational Cost: add. (31370->290), mult. (68810->366), div. (0->0), fcn. (46354->10), ass. (0->114)
t486 = qJD(1) ^ 2;
t480 = sin(qJ(1));
t484 = cos(qJ(1));
t500 = g(1) * t480 - t484 * g(2);
t491 = -qJ(2) * t486 + qJDD(2) - t500;
t507 = -pkin(1) - pkin(7);
t447 = t507 * qJDD(1) + t491;
t479 = sin(qJ(3));
t483 = cos(qJ(3));
t437 = t479 * g(3) + t483 * t447;
t503 = qJD(1) * qJD(3);
t501 = t479 * t503;
t463 = qJDD(1) * t483 - t501;
t414 = (-t463 - t501) * qJ(4) + (-t479 * t483 * t486 + qJDD(3)) * pkin(3) + t437;
t438 = -g(3) * t483 + t479 * t447;
t462 = -qJDD(1) * t479 - t483 * t503;
t505 = qJD(1) * t483;
t465 = qJD(3) * pkin(3) - qJ(4) * t505;
t474 = t479 ^ 2;
t415 = -pkin(3) * t474 * t486 + qJ(4) * t462 - qJD(3) * t465 + t438;
t475 = sin(pkin(10));
t476 = cos(pkin(10));
t506 = qJD(1) * t479;
t455 = -t475 * t506 + t476 * t505;
t395 = -0.2e1 * qJD(4) * t455 + t414 * t476 - t475 * t415;
t495 = -g(1) * t484 - g(2) * t480;
t492 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t495;
t454 = -t475 * t505 - t476 * t506;
t396 = 0.2e1 * qJD(4) * t454 + t475 * t414 + t476 * t415;
t429 = -mrSges(5,1) * t454 + mrSges(5,2) * t455;
t435 = t462 * t476 - t475 * t463;
t446 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t455;
t431 = -pkin(4) * t454 - pkin(8) * t455;
t485 = qJD(3) ^ 2;
t385 = -pkin(4) * t485 + qJDD(3) * pkin(8) + t431 * t454 + t396;
t417 = -pkin(3) * t462 + qJDD(4) + t465 * t505 + (-qJ(4) * t474 + t507) * t486 + t492;
t436 = t462 * t475 + t463 * t476;
t394 = (-qJD(3) * t454 - t436) * pkin(8) + (qJD(3) * t455 - t435) * pkin(4) + t417;
t478 = sin(qJ(5));
t482 = cos(qJ(5));
t380 = -t385 * t478 + t482 * t394;
t440 = qJD(3) * t482 - t455 * t478;
t410 = qJD(5) * t440 + qJDD(3) * t478 + t436 * t482;
t434 = qJDD(5) - t435;
t441 = qJD(3) * t478 + t455 * t482;
t452 = qJD(5) - t454;
t378 = (t440 * t452 - t410) * pkin(9) + (t440 * t441 + t434) * pkin(5) + t380;
t381 = t482 * t385 + t478 * t394;
t409 = -qJD(5) * t441 + qJDD(3) * t482 - t436 * t478;
t424 = pkin(5) * t452 - pkin(9) * t441;
t439 = t440 ^ 2;
t379 = -pkin(5) * t439 + pkin(9) * t409 - t424 * t452 + t381;
t477 = sin(qJ(6));
t481 = cos(qJ(6));
t376 = t378 * t481 - t379 * t477;
t418 = t440 * t481 - t441 * t477;
t390 = qJD(6) * t418 + t409 * t477 + t410 * t481;
t419 = t440 * t477 + t441 * t481;
t401 = -mrSges(7,1) * t418 + mrSges(7,2) * t419;
t449 = qJD(6) + t452;
t402 = -mrSges(7,2) * t449 + mrSges(7,3) * t418;
t432 = qJDD(6) + t434;
t373 = m(7) * t376 + mrSges(7,1) * t432 - mrSges(7,3) * t390 - t401 * t419 + t402 * t449;
t377 = t378 * t477 + t379 * t481;
t389 = -qJD(6) * t419 + t409 * t481 - t410 * t477;
t403 = mrSges(7,1) * t449 - mrSges(7,3) * t419;
t374 = m(7) * t377 - mrSges(7,2) * t432 + mrSges(7,3) * t389 + t401 * t418 - t403 * t449;
t365 = t481 * t373 + t477 * t374;
t420 = -mrSges(6,1) * t440 + mrSges(6,2) * t441;
t422 = -mrSges(6,2) * t452 + mrSges(6,3) * t440;
t363 = m(6) * t380 + mrSges(6,1) * t434 - mrSges(6,3) * t410 - t420 * t441 + t422 * t452 + t365;
t423 = mrSges(6,1) * t452 - mrSges(6,3) * t441;
t497 = -t373 * t477 + t481 * t374;
t364 = m(6) * t381 - mrSges(6,2) * t434 + mrSges(6,3) * t409 + t420 * t440 - t423 * t452 + t497;
t498 = -t478 * t363 + t482 * t364;
t357 = m(5) * t396 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t435 - qJD(3) * t446 + t429 * t454 + t498;
t445 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t454;
t384 = -qJDD(3) * pkin(4) - pkin(8) * t485 + t455 * t431 - t395;
t382 = -pkin(5) * t409 - pkin(9) * t439 + t424 * t441 + t384;
t490 = m(7) * t382 - t389 * mrSges(7,1) + t390 * mrSges(7,2) - t418 * t402 + t419 * t403;
t488 = -m(6) * t384 + t409 * mrSges(6,1) - t410 * mrSges(6,2) + t440 * t422 - t441 * t423 - t490;
t369 = m(5) * t395 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t436 + qJD(3) * t445 - t429 * t455 + t488;
t352 = t475 * t357 + t476 * t369;
t359 = t482 * t363 + t478 * t364;
t499 = t476 * t357 - t369 * t475;
t461 = (mrSges(4,1) * t479 + mrSges(4,2) * t483) * qJD(1);
t464 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t506;
t466 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t505;
t494 = (m(4) * t437 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t463 + qJD(3) * t464 - t461 * t505 + t352) * t483 + (m(4) * t438 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t462 - qJD(3) * t466 - t461 * t506 + t499) * t479;
t358 = m(5) * t417 - mrSges(5,1) * t435 + t436 * mrSges(5,2) - t445 * t454 + t455 * t446 + t359;
t398 = Ifges(7,4) * t419 + Ifges(7,2) * t418 + Ifges(7,6) * t449;
t399 = Ifges(7,1) * t419 + Ifges(7,4) * t418 + Ifges(7,5) * t449;
t489 = -mrSges(7,1) * t376 + mrSges(7,2) * t377 - Ifges(7,5) * t390 - Ifges(7,6) * t389 - Ifges(7,3) * t432 - t419 * t398 + t418 * t399;
t405 = Ifges(6,4) * t441 + Ifges(6,2) * t440 + Ifges(6,6) * t452;
t406 = Ifges(6,1) * t441 + Ifges(6,4) * t440 + Ifges(6,5) * t452;
t487 = mrSges(6,1) * t380 - mrSges(6,2) * t381 + Ifges(6,5) * t410 + Ifges(6,6) * t409 + Ifges(6,3) * t434 + pkin(5) * t365 + t441 * t405 - t440 * t406 - t489;
t458 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t483 - Ifges(4,4) * t479) * qJD(1);
t457 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t483 - Ifges(4,2) * t479) * qJD(1);
t453 = -qJDD(1) * pkin(1) + t491;
t448 = pkin(1) * t486 - t492;
t444 = t507 * t486 + t492;
t427 = Ifges(5,1) * t455 + Ifges(5,4) * t454 + Ifges(5,5) * qJD(3);
t426 = Ifges(5,4) * t455 + Ifges(5,2) * t454 + Ifges(5,6) * qJD(3);
t425 = Ifges(5,5) * t455 + Ifges(5,6) * t454 + Ifges(5,3) * qJD(3);
t404 = Ifges(6,5) * t441 + Ifges(6,6) * t440 + Ifges(6,3) * t452;
t397 = Ifges(7,5) * t419 + Ifges(7,6) * t418 + Ifges(7,3) * t449;
t367 = mrSges(7,2) * t382 - mrSges(7,3) * t376 + Ifges(7,1) * t390 + Ifges(7,4) * t389 + Ifges(7,5) * t432 + t397 * t418 - t398 * t449;
t366 = -mrSges(7,1) * t382 + mrSges(7,3) * t377 + Ifges(7,4) * t390 + Ifges(7,2) * t389 + Ifges(7,6) * t432 - t397 * t419 + t399 * t449;
t354 = mrSges(6,2) * t384 - mrSges(6,3) * t380 + Ifges(6,1) * t410 + Ifges(6,4) * t409 + Ifges(6,5) * t434 - pkin(9) * t365 - t366 * t477 + t367 * t481 + t404 * t440 - t405 * t452;
t353 = -mrSges(6,1) * t384 + mrSges(6,3) * t381 + Ifges(6,4) * t410 + Ifges(6,2) * t409 + Ifges(6,6) * t434 - pkin(5) * t490 + pkin(9) * t497 + t481 * t366 + t477 * t367 - t441 * t404 + t452 * t406;
t349 = -mrSges(5,1) * t417 + mrSges(5,3) * t396 + Ifges(5,4) * t436 + Ifges(5,2) * t435 + Ifges(5,6) * qJDD(3) - pkin(4) * t359 + qJD(3) * t427 - t455 * t425 - t487;
t348 = m(3) * t453 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t486 + t494;
t347 = mrSges(5,2) * t417 - mrSges(5,3) * t395 + Ifges(5,1) * t436 + Ifges(5,4) * t435 + Ifges(5,5) * qJDD(3) - pkin(8) * t359 - qJD(3) * t426 - t353 * t478 + t354 * t482 + t425 * t454;
t1 = [mrSges(2,1) * t500 - mrSges(2,2) * t495 + mrSges(3,2) * t453 - mrSges(3,3) * t448 + t483 * (mrSges(4,2) * t444 - mrSges(4,3) * t437 + Ifges(4,1) * t463 + Ifges(4,4) * t462 + Ifges(4,5) * qJDD(3) - qJ(4) * t352 - qJD(3) * t457 + t347 * t476 - t349 * t475) - t479 * (-mrSges(4,1) * t444 + mrSges(4,3) * t438 + Ifges(4,4) * t463 + Ifges(4,2) * t462 + Ifges(4,6) * qJDD(3) - pkin(3) * t358 + qJ(4) * t499 + qJD(3) * t458 + t475 * t347 + t476 * t349) - pkin(7) * t494 - pkin(1) * t348 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t448 + m(4) * t444 - mrSges(4,1) * t462 + mrSges(3,2) * t486 + mrSges(4,2) * t463 + t358 + qJDD(1) * mrSges(3,3) + (t464 * t479 + t466 * t483) * qJD(1)) * qJ(2); t348; Ifges(4,5) * t463 + Ifges(4,6) * t462 + mrSges(4,1) * t437 - mrSges(4,2) * t438 + Ifges(5,5) * t436 + Ifges(5,6) * t435 + t455 * t426 - t454 * t427 + mrSges(5,1) * t395 - mrSges(5,2) * t396 + t478 * t354 + t482 * t353 + pkin(4) * t488 + pkin(8) * t498 + pkin(3) * t352 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t457 * t483 + t458 * t479) * qJD(1); t358; t487; -t489;];
tauJ  = t1;
