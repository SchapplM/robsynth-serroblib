% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR7
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
% Datum: 2019-05-05 19:14
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:12:46
% EndTime: 2019-05-05 19:12:49
% DurationCPUTime: 3.36s
% Computational Cost: add. (33358->289), mult. (74238->368), div. (0->0), fcn. (51778->10), ass. (0->114)
t470 = sin(qJ(1));
t474 = cos(qJ(1));
t484 = -g(1) * t474 - g(2) * t470;
t481 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t484;
t494 = -pkin(1) - pkin(7);
t475 = qJD(1) ^ 2;
t488 = g(1) * t470 - t474 * g(2);
t480 = -qJ(2) * t475 + qJDD(2) - t488;
t438 = t494 * qJDD(1) + t480;
t469 = sin(qJ(3));
t473 = cos(qJ(3));
t430 = t469 * g(3) + t473 * t438;
t491 = qJD(1) * qJD(3);
t489 = t469 * t491;
t452 = qJDD(1) * t473 - t489;
t410 = (-t452 - t489) * qJ(4) + (-t469 * t473 * t475 + qJDD(3)) * pkin(3) + t430;
t431 = -g(3) * t473 + t469 * t438;
t451 = -qJDD(1) * t469 - t473 * t491;
t492 = qJD(1) * t473;
t454 = qJD(3) * pkin(3) - qJ(4) * t492;
t464 = t469 ^ 2;
t411 = -pkin(3) * t464 * t475 + qJ(4) * t451 - qJD(3) * t454 + t431;
t465 = sin(pkin(10));
t466 = cos(pkin(10));
t445 = (-t465 * t469 + t466 * t473) * qJD(1);
t387 = -0.2e1 * qJD(4) * t445 + t466 * t410 - t411 * t465;
t429 = t451 * t465 + t452 * t466;
t444 = (-t465 * t473 - t466 * t469) * qJD(1);
t377 = (qJD(3) * t444 - t429) * pkin(8) + (t444 * t445 + qJDD(3)) * pkin(4) + t387;
t388 = 0.2e1 * qJD(4) * t444 + t465 * t410 + t466 * t411;
t428 = t451 * t466 - t452 * t465;
t437 = qJD(3) * pkin(4) - pkin(8) * t445;
t443 = t444 ^ 2;
t379 = -pkin(4) * t443 + pkin(8) * t428 - qJD(3) * t437 + t388;
t468 = sin(qJ(5));
t472 = cos(qJ(5));
t374 = t468 * t377 + t472 * t379;
t423 = t444 * t468 + t445 * t472;
t396 = -qJD(5) * t423 + t428 * t472 - t429 * t468;
t422 = t444 * t472 - t445 * t468;
t405 = -mrSges(6,1) * t422 + mrSges(6,2) * t423;
t462 = qJD(3) + qJD(5);
t417 = mrSges(6,1) * t462 - mrSges(6,3) * t423;
t461 = qJDD(3) + qJDD(5);
t406 = -pkin(5) * t422 - pkin(9) * t423;
t460 = t462 ^ 2;
t371 = -pkin(5) * t460 + pkin(9) * t461 + t406 * t422 + t374;
t413 = -pkin(3) * t451 + qJDD(4) + t454 * t492 + (-qJ(4) * t464 + t494) * t475 + t481;
t390 = -pkin(4) * t428 - pkin(8) * t443 + t445 * t437 + t413;
t397 = qJD(5) * t422 + t428 * t468 + t429 * t472;
t375 = t390 + (t423 * t462 - t396) * pkin(5) + (-t422 * t462 - t397) * pkin(9);
t467 = sin(qJ(6));
t471 = cos(qJ(6));
t368 = -t371 * t467 + t375 * t471;
t414 = -t423 * t467 + t462 * t471;
t382 = qJD(6) * t414 + t397 * t471 + t461 * t467;
t395 = qJDD(6) - t396;
t415 = t423 * t471 + t462 * t467;
t398 = -mrSges(7,1) * t414 + mrSges(7,2) * t415;
t418 = qJD(6) - t422;
t399 = -mrSges(7,2) * t418 + mrSges(7,3) * t414;
t365 = m(7) * t368 + mrSges(7,1) * t395 - mrSges(7,3) * t382 - t398 * t415 + t399 * t418;
t369 = t371 * t471 + t375 * t467;
t381 = -qJD(6) * t415 - t397 * t467 + t461 * t471;
t400 = mrSges(7,1) * t418 - mrSges(7,3) * t415;
t366 = m(7) * t369 - mrSges(7,2) * t395 + mrSges(7,3) * t381 + t398 * t414 - t400 * t418;
t485 = -t365 * t467 + t471 * t366;
t352 = m(6) * t374 - mrSges(6,2) * t461 + mrSges(6,3) * t396 + t405 * t422 - t417 * t462 + t485;
t373 = t377 * t472 - t379 * t468;
t416 = -mrSges(6,2) * t462 + mrSges(6,3) * t422;
t370 = -pkin(5) * t461 - pkin(9) * t460 + t406 * t423 - t373;
t479 = -m(7) * t370 + t381 * mrSges(7,1) - mrSges(7,2) * t382 + t414 * t399 - t400 * t415;
t361 = m(6) * t373 + mrSges(6,1) * t461 - mrSges(6,3) * t397 - t405 * t423 + t416 * t462 + t479;
t349 = t468 * t352 + t472 * t361;
t426 = -mrSges(5,1) * t444 + mrSges(5,2) * t445;
t435 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t444;
t347 = m(5) * t387 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t429 + qJD(3) * t435 - t426 * t445 + t349;
t436 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t445;
t486 = t472 * t352 - t361 * t468;
t348 = m(5) * t388 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t428 - qJD(3) * t436 + t426 * t444 + t486;
t341 = t466 * t347 + t465 * t348;
t355 = t471 * t365 + t467 * t366;
t493 = qJD(1) * t469;
t487 = -t347 * t465 + t466 * t348;
t450 = (mrSges(4,1) * t469 + mrSges(4,2) * t473) * qJD(1);
t453 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t493;
t455 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t492;
t483 = (m(4) * t430 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t452 + qJD(3) * t453 - t450 * t492 + t341) * t473 + (m(4) * t431 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t451 - qJD(3) * t455 - t450 * t493 + t487) * t469;
t478 = m(6) * t390 - t396 * mrSges(6,1) + t397 * mrSges(6,2) - t422 * t416 + t423 * t417 + t355;
t383 = Ifges(7,5) * t415 + Ifges(7,6) * t414 + Ifges(7,3) * t418;
t385 = Ifges(7,1) * t415 + Ifges(7,4) * t414 + Ifges(7,5) * t418;
t358 = -mrSges(7,1) * t370 + mrSges(7,3) * t369 + Ifges(7,4) * t382 + Ifges(7,2) * t381 + Ifges(7,6) * t395 - t383 * t415 + t385 * t418;
t384 = Ifges(7,4) * t415 + Ifges(7,2) * t414 + Ifges(7,6) * t418;
t359 = mrSges(7,2) * t370 - mrSges(7,3) * t368 + Ifges(7,1) * t382 + Ifges(7,4) * t381 + Ifges(7,5) * t395 + t383 * t414 - t384 * t418;
t402 = Ifges(6,4) * t423 + Ifges(6,2) * t422 + Ifges(6,6) * t462;
t403 = Ifges(6,1) * t423 + Ifges(6,4) * t422 + Ifges(6,5) * t462;
t477 = mrSges(6,1) * t373 - mrSges(6,2) * t374 + Ifges(6,5) * t397 + Ifges(6,6) * t396 + Ifges(6,3) * t461 + pkin(5) * t479 + pkin(9) * t485 + t471 * t358 + t467 * t359 + t423 * t402 - t422 * t403;
t476 = mrSges(7,1) * t368 - mrSges(7,2) * t369 + Ifges(7,5) * t382 + Ifges(7,6) * t381 + Ifges(7,3) * t395 + t384 * t415 - t385 * t414;
t353 = m(5) * t413 - t428 * mrSges(5,1) + t429 * mrSges(5,2) - t444 * t435 + t445 * t436 + t478;
t448 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t473 - Ifges(4,4) * t469) * qJD(1);
t447 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t473 - Ifges(4,2) * t469) * qJD(1);
t442 = -qJDD(1) * pkin(1) + t480;
t439 = pkin(1) * t475 - t481;
t434 = t494 * t475 + t481;
t421 = Ifges(5,1) * t445 + Ifges(5,4) * t444 + (Ifges(5,5) * qJD(3));
t420 = Ifges(5,4) * t445 + Ifges(5,2) * t444 + (Ifges(5,6) * qJD(3));
t419 = Ifges(5,5) * t445 + Ifges(5,6) * t444 + (Ifges(5,3) * qJD(3));
t401 = Ifges(6,5) * t423 + Ifges(6,6) * t422 + Ifges(6,3) * t462;
t343 = -mrSges(6,1) * t390 + mrSges(6,3) * t374 + Ifges(6,4) * t397 + Ifges(6,2) * t396 + Ifges(6,6) * t461 - pkin(5) * t355 - t401 * t423 + t403 * t462 - t476;
t342 = mrSges(6,2) * t390 - mrSges(6,3) * t373 + Ifges(6,1) * t397 + Ifges(6,4) * t396 + Ifges(6,5) * t461 - pkin(9) * t355 - t358 * t467 + t359 * t471 + t401 * t422 - t402 * t462;
t338 = mrSges(5,2) * t413 - mrSges(5,3) * t387 + Ifges(5,1) * t429 + Ifges(5,4) * t428 + Ifges(5,5) * qJDD(3) - pkin(8) * t349 - qJD(3) * t420 + t342 * t472 - t343 * t468 + t419 * t444;
t337 = m(3) * t442 + qJDD(1) * mrSges(3,2) - (mrSges(3,3) * t475) + t483;
t336 = -mrSges(5,1) * t413 + mrSges(5,3) * t388 + Ifges(5,4) * t429 + Ifges(5,2) * t428 + Ifges(5,6) * qJDD(3) - pkin(4) * t478 + pkin(8) * t486 + qJD(3) * t421 + t468 * t342 + t472 * t343 - t445 * t419;
t1 = [mrSges(2,1) * t488 - mrSges(2,2) * t484 + mrSges(3,2) * t442 - mrSges(3,3) * t439 + t473 * (mrSges(4,2) * t434 - mrSges(4,3) * t430 + Ifges(4,1) * t452 + Ifges(4,4) * t451 + Ifges(4,5) * qJDD(3) - qJ(4) * t341 - qJD(3) * t447 - t336 * t465 + t338 * t466) - t469 * (-mrSges(4,1) * t434 + mrSges(4,3) * t431 + Ifges(4,4) * t452 + Ifges(4,2) * t451 + Ifges(4,6) * qJDD(3) - pkin(3) * t353 + qJ(4) * t487 + qJD(3) * t448 + t466 * t336 + t465 * t338) - pkin(7) * t483 - pkin(1) * t337 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t439 + m(4) * t434 - t451 * mrSges(4,1) + t475 * mrSges(3,2) + t452 * mrSges(4,2) + t353 + qJDD(1) * mrSges(3,3) + (t453 * t469 + t455 * t473) * qJD(1)) * qJ(2); t337; t477 + pkin(4) * t349 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t473 * t447 + t469 * t448) * qJD(1) + pkin(3) * t341 + Ifges(4,6) * t451 + Ifges(4,5) * t452 + mrSges(5,1) * t387 - mrSges(5,2) * t388 + Ifges(5,6) * t428 + Ifges(5,5) * t429 + mrSges(4,1) * t430 - mrSges(4,2) * t431 - t444 * t421 + t445 * t420; t353; t477; t476;];
tauJ  = t1;
