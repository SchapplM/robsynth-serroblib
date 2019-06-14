% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRRR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR7_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:05:41
% EndTime: 2019-05-05 16:05:45
% DurationCPUTime: 3.17s
% Computational Cost: add. (29123->258), mult. (67053->326), div. (0->0), fcn. (49684->10), ass. (0->114)
t480 = qJD(1) ^ 2;
t475 = sin(qJ(1));
t479 = cos(qJ(1));
t499 = g(1) * t475 - t479 * g(2);
t488 = -qJ(2) * t480 + qJDD(2) - t499;
t508 = -pkin(1) - qJ(3);
t512 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t508 + t488;
t471 = cos(pkin(10));
t511 = t471 ^ 2;
t494 = -g(1) * t479 - g(2) * t475;
t510 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t494;
t509 = pkin(3) * t480;
t507 = t471 * mrSges(4,2);
t470 = sin(pkin(10));
t436 = t470 * g(3) + t512 * t471;
t418 = (-pkin(7) * qJDD(1) - t470 * t509) * t471 + t436;
t437 = -g(3) * t471 + t512 * t470;
t466 = t470 ^ 2;
t501 = qJDD(1) * t470;
t421 = -pkin(7) * t501 - t466 * t509 + t437;
t474 = sin(qJ(4));
t478 = cos(qJ(4));
t404 = t478 * t418 - t421 * t474;
t491 = -t470 * t474 + t471 * t478;
t492 = -t470 * t478 - t471 * t474;
t451 = t492 * qJD(1);
t503 = qJD(4) * t451;
t439 = qJDD(1) * t491 + t503;
t452 = t491 * qJD(1);
t385 = (-t439 + t503) * pkin(8) + (t451 * t452 + qJDD(4)) * pkin(4) + t404;
t405 = t474 * t418 + t478 * t421;
t438 = -qJD(4) * t452 + qJDD(1) * t492;
t446 = qJD(4) * pkin(4) - pkin(8) * t452;
t450 = t451 ^ 2;
t387 = -pkin(4) * t450 + pkin(8) * t438 - qJD(4) * t446 + t405;
t473 = sin(qJ(5));
t477 = cos(qJ(5));
t383 = t473 * t385 + t477 * t387;
t431 = t451 * t473 + t452 * t477;
t402 = -qJD(5) * t431 + t438 * t477 - t439 * t473;
t430 = t451 * t477 - t452 * t473;
t413 = -mrSges(6,1) * t430 + mrSges(6,2) * t431;
t468 = qJD(4) + qJD(5);
t423 = mrSges(6,1) * t468 - mrSges(6,3) * t431;
t465 = qJDD(4) + qJDD(5);
t414 = -pkin(5) * t430 - pkin(9) * t431;
t464 = t468 ^ 2;
t379 = -pkin(5) * t464 + pkin(9) * t465 + t414 * t430 + t383;
t487 = qJDD(3) + t510;
t505 = -t466 - t511;
t425 = pkin(3) * t501 + (pkin(7) * t505 + t508) * t480 + t487;
t396 = -pkin(4) * t438 - pkin(8) * t450 + t452 * t446 + t425;
t403 = qJD(5) * t430 + t438 * t473 + t439 * t477;
t380 = (-t430 * t468 - t403) * pkin(9) + (t431 * t468 - t402) * pkin(5) + t396;
t472 = sin(qJ(6));
t476 = cos(qJ(6));
t376 = -t379 * t472 + t380 * t476;
t419 = -t431 * t472 + t468 * t476;
t390 = qJD(6) * t419 + t403 * t476 + t465 * t472;
t401 = qJDD(6) - t402;
t420 = t431 * t476 + t468 * t472;
t406 = -mrSges(7,1) * t419 + mrSges(7,2) * t420;
t426 = qJD(6) - t430;
t407 = -mrSges(7,2) * t426 + mrSges(7,3) * t419;
t373 = m(7) * t376 + mrSges(7,1) * t401 - mrSges(7,3) * t390 - t406 * t420 + t407 * t426;
t377 = t379 * t476 + t380 * t472;
t389 = -qJD(6) * t420 - t403 * t472 + t465 * t476;
t408 = mrSges(7,1) * t426 - mrSges(7,3) * t420;
t374 = m(7) * t377 - mrSges(7,2) * t401 + mrSges(7,3) * t389 + t406 * t419 - t408 * t426;
t495 = -t373 * t472 + t476 * t374;
t361 = m(6) * t383 - mrSges(6,2) * t465 + mrSges(6,3) * t402 + t413 * t430 - t423 * t468 + t495;
t382 = t385 * t477 - t387 * t473;
t422 = -mrSges(6,2) * t468 + mrSges(6,3) * t430;
t378 = -pkin(5) * t465 - pkin(9) * t464 + t414 * t431 - t382;
t486 = -m(7) * t378 + t389 * mrSges(7,1) - mrSges(7,2) * t390 + t419 * t407 - t408 * t420;
t369 = m(6) * t382 + mrSges(6,1) * t465 - mrSges(6,3) * t403 - t413 * t431 + t422 * t468 + t486;
t358 = t473 * t361 + t477 * t369;
t434 = -mrSges(5,1) * t451 + mrSges(5,2) * t452;
t444 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t451;
t356 = m(5) * t404 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t439 + qJD(4) * t444 - t434 * t452 + t358;
t445 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t452;
t496 = t477 * t361 - t369 * t473;
t357 = m(5) * t405 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t438 - qJD(4) * t445 + t434 * t451 + t496;
t506 = t478 * t356 + t474 * t357;
t363 = t476 * t373 + t472 * t374;
t498 = t505 * mrSges(4,3);
t497 = -t356 * t474 + t478 * t357;
t490 = -mrSges(4,3) * qJDD(1) - t480 * (t470 * mrSges(4,1) + t507);
t493 = (m(4) * t436 + t471 * t490 + t506) * t471 + (m(4) * t437 + t470 * t490 + t497) * t470;
t485 = m(6) * t396 - mrSges(6,1) * t402 + t403 * mrSges(6,2) - t422 * t430 + t431 * t423 + t363;
t391 = Ifges(7,5) * t420 + Ifges(7,6) * t419 + Ifges(7,3) * t426;
t393 = Ifges(7,1) * t420 + Ifges(7,4) * t419 + Ifges(7,5) * t426;
t366 = -mrSges(7,1) * t378 + mrSges(7,3) * t377 + Ifges(7,4) * t390 + Ifges(7,2) * t389 + Ifges(7,6) * t401 - t391 * t420 + t393 * t426;
t392 = Ifges(7,4) * t420 + Ifges(7,2) * t419 + Ifges(7,6) * t426;
t367 = mrSges(7,2) * t378 - mrSges(7,3) * t376 + Ifges(7,1) * t390 + Ifges(7,4) * t389 + Ifges(7,5) * t401 + t391 * t419 - t392 * t426;
t410 = Ifges(6,4) * t431 + Ifges(6,2) * t430 + Ifges(6,6) * t468;
t411 = Ifges(6,1) * t431 + Ifges(6,4) * t430 + Ifges(6,5) * t468;
t484 = mrSges(6,1) * t382 - mrSges(6,2) * t383 + Ifges(6,5) * t403 + Ifges(6,6) * t402 + Ifges(6,3) * t465 + pkin(5) * t486 + pkin(9) * t495 + t476 * t366 + t472 * t367 + t431 * t410 - t411 * t430;
t483 = mrSges(7,1) * t376 - mrSges(7,2) * t377 + Ifges(7,5) * t390 + Ifges(7,6) * t389 + Ifges(7,3) * t401 + t392 * t420 - t393 * t419;
t482 = m(5) * t425 - mrSges(5,1) * t438 + t439 * mrSges(5,2) - t444 * t451 + t452 * t445 + t485;
t443 = t480 * t508 + t487;
t481 = m(4) * t443 + mrSges(4,1) * t501 + qJDD(1) * t507 + t482;
t449 = -qJDD(1) * pkin(1) + t488;
t448 = pkin(1) * t480 - t510;
t429 = Ifges(5,1) * t452 + Ifges(5,4) * t451 + Ifges(5,5) * qJD(4);
t428 = Ifges(5,4) * t452 + Ifges(5,2) * t451 + Ifges(5,6) * qJD(4);
t427 = Ifges(5,5) * t452 + Ifges(5,6) * t451 + Ifges(5,3) * qJD(4);
t409 = Ifges(6,5) * t431 + Ifges(6,6) * t430 + Ifges(6,3) * t468;
t352 = -mrSges(6,1) * t396 + mrSges(6,3) * t383 + Ifges(6,4) * t403 + Ifges(6,2) * t402 + Ifges(6,6) * t465 - pkin(5) * t363 - t409 * t431 + t411 * t468 - t483;
t351 = mrSges(6,2) * t396 - mrSges(6,3) * t382 + Ifges(6,1) * t403 + Ifges(6,4) * t402 + Ifges(6,5) * t465 - pkin(9) * t363 - t366 * t472 + t367 * t476 + t409 * t430 - t410 * t468;
t348 = mrSges(5,2) * t425 - mrSges(5,3) * t404 + Ifges(5,1) * t439 + Ifges(5,4) * t438 + Ifges(5,5) * qJDD(4) - pkin(8) * t358 - qJD(4) * t428 + t351 * t477 - t352 * t473 + t427 * t451;
t347 = m(3) * t449 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t480 + t493;
t346 = -mrSges(5,1) * t425 + mrSges(5,3) * t405 + Ifges(5,4) * t439 + Ifges(5,2) * t438 + Ifges(5,6) * qJDD(4) - pkin(4) * t485 + pkin(8) * t496 + qJD(4) * t429 + t473 * t351 + t477 * t352 - t452 * t427;
t1 = [mrSges(2,1) * t499 - mrSges(2,2) * t494 + mrSges(3,2) * t449 - mrSges(3,3) * t448 + t471 * (mrSges(4,2) * t443 - mrSges(4,3) * t436 - pkin(7) * t506 - t474 * t346 + t478 * t348) - t470 * (-mrSges(4,1) * t443 + mrSges(4,3) * t437 - pkin(3) * t482 + pkin(7) * t497 + t478 * t346 + t474 * t348) - qJ(3) * t493 - pkin(1) * t347 + qJ(2) * (t481 + (mrSges(3,2) + t498) * t480 - m(3) * t448) + (Ifges(4,1) * t511 + qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (-0.2e1 * Ifges(4,4) * t471 + Ifges(4,2) * t470) * t470) * qJDD(1); t347; t480 * t498 + t481; mrSges(5,1) * t404 - mrSges(5,2) * t405 + Ifges(5,5) * t439 + Ifges(5,6) * t438 + Ifges(5,3) * qJDD(4) + pkin(4) * t358 + t428 * t452 - t429 * t451 + t484; t484; t483;];
tauJ  = t1;
