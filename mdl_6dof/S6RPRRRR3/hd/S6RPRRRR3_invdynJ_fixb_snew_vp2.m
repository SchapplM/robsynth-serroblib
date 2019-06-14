% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRR3
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
% Datum: 2019-05-06 03:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:57:41
% EndTime: 2019-05-06 02:57:45
% DurationCPUTime: 4.08s
% Computational Cost: add. (52797->294), mult. (104140->369), div. (0->0), fcn. (71760->12), ass. (0->122)
t466 = sin(qJ(1));
t471 = cos(qJ(1));
t486 = t466 * g(1) - g(2) * t471;
t445 = qJDD(1) * pkin(1) + t486;
t473 = qJD(1) ^ 2;
t481 = -g(1) * t471 - g(2) * t466;
t447 = -pkin(1) * t473 + t481;
t460 = sin(pkin(11));
t461 = cos(pkin(11));
t422 = t461 * t445 - t460 * t447;
t409 = -qJDD(1) * pkin(2) - t473 * pkin(7) - t422;
t465 = sin(qJ(3));
t470 = cos(qJ(3));
t488 = qJD(1) * qJD(3);
t487 = t470 * t488;
t449 = qJDD(1) * t465 + t487;
t457 = t465 * t488;
t450 = qJDD(1) * t470 - t457;
t394 = (-t449 - t487) * pkin(8) + (-t450 + t457) * pkin(3) + t409;
t423 = t460 * t445 + t461 * t447;
t410 = -pkin(2) * t473 + qJDD(1) * pkin(7) + t423;
t459 = -g(3) + qJDD(2);
t403 = t470 * t410 + t465 * t459;
t448 = (-pkin(3) * t470 - pkin(8) * t465) * qJD(1);
t472 = qJD(3) ^ 2;
t489 = qJD(1) * t470;
t400 = -pkin(3) * t472 + qJDD(3) * pkin(8) + t448 * t489 + t403;
t464 = sin(qJ(4));
t469 = cos(qJ(4));
t378 = t469 * t394 - t464 * t400;
t490 = qJD(1) * t465;
t443 = qJD(3) * t469 - t464 * t490;
t418 = qJD(4) * t443 + qJDD(3) * t464 + t449 * t469;
t442 = qJDD(4) - t450;
t444 = qJD(3) * t464 + t469 * t490;
t455 = qJD(4) - t489;
t368 = (t443 * t455 - t418) * pkin(9) + (t443 * t444 + t442) * pkin(4) + t378;
t379 = t464 * t394 + t469 * t400;
t417 = -qJD(4) * t444 + qJDD(3) * t469 - t449 * t464;
t427 = pkin(4) * t455 - pkin(9) * t444;
t441 = t443 ^ 2;
t370 = -pkin(4) * t441 + pkin(9) * t417 - t427 * t455 + t379;
t463 = sin(qJ(5));
t468 = cos(qJ(5));
t356 = t468 * t368 - t463 * t370;
t420 = t443 * t468 - t444 * t463;
t386 = qJD(5) * t420 + t417 * t463 + t418 * t468;
t421 = t443 * t463 + t444 * t468;
t438 = qJDD(5) + t442;
t454 = qJD(5) + t455;
t353 = (t420 * t454 - t386) * pkin(10) + (t420 * t421 + t438) * pkin(5) + t356;
t357 = t463 * t368 + t468 * t370;
t385 = -qJD(5) * t421 + t417 * t468 - t418 * t463;
t406 = pkin(5) * t454 - pkin(10) * t421;
t419 = t420 ^ 2;
t354 = -pkin(5) * t419 + pkin(10) * t385 - t406 * t454 + t357;
t462 = sin(qJ(6));
t467 = cos(qJ(6));
t351 = t353 * t467 - t354 * t462;
t397 = t420 * t467 - t421 * t462;
t365 = qJD(6) * t397 + t385 * t462 + t386 * t467;
t398 = t420 * t462 + t421 * t467;
t380 = -mrSges(7,1) * t397 + mrSges(7,2) * t398;
t451 = qJD(6) + t454;
t387 = -mrSges(7,2) * t451 + mrSges(7,3) * t397;
t431 = qJDD(6) + t438;
t348 = m(7) * t351 + mrSges(7,1) * t431 - mrSges(7,3) * t365 - t380 * t398 + t387 * t451;
t352 = t353 * t462 + t354 * t467;
t364 = -qJD(6) * t398 + t385 * t467 - t386 * t462;
t388 = mrSges(7,1) * t451 - mrSges(7,3) * t398;
t349 = m(7) * t352 - mrSges(7,2) * t431 + mrSges(7,3) * t364 + t380 * t397 - t388 * t451;
t341 = t467 * t348 + t462 * t349;
t401 = -mrSges(6,1) * t420 + mrSges(6,2) * t421;
t404 = -mrSges(6,2) * t454 + mrSges(6,3) * t420;
t338 = m(6) * t356 + mrSges(6,1) * t438 - mrSges(6,3) * t386 - t401 * t421 + t404 * t454 + t341;
t405 = mrSges(6,1) * t454 - mrSges(6,3) * t421;
t482 = -t348 * t462 + t467 * t349;
t339 = m(6) * t357 - mrSges(6,2) * t438 + mrSges(6,3) * t385 + t401 * t420 - t405 * t454 + t482;
t334 = t468 * t338 + t463 * t339;
t446 = (-mrSges(4,1) * t470 + mrSges(4,2) * t465) * qJD(1);
t452 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t490;
t424 = -mrSges(5,1) * t443 + mrSges(5,2) * t444;
t425 = -mrSges(5,2) * t455 + mrSges(5,3) * t443;
t332 = m(5) * t378 + mrSges(5,1) * t442 - mrSges(5,3) * t418 - t424 * t444 + t425 * t455 + t334;
t426 = mrSges(5,1) * t455 - mrSges(5,3) * t444;
t483 = -t338 * t463 + t468 * t339;
t333 = m(5) * t379 - mrSges(5,2) * t442 + mrSges(5,3) * t417 + t424 * t443 - t426 * t455 + t483;
t484 = -t332 * t464 + t469 * t333;
t327 = m(4) * t403 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t450 - qJD(3) * t452 + t446 * t489 + t484;
t402 = -t465 * t410 + t459 * t470;
t453 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t489;
t399 = -qJDD(3) * pkin(3) - pkin(8) * t472 + t448 * t490 - t402;
t377 = -pkin(4) * t417 - pkin(9) * t441 + t444 * t427 + t399;
t359 = -pkin(5) * t385 - pkin(10) * t419 + t406 * t421 + t377;
t480 = m(7) * t359 - t364 * mrSges(7,1) + t365 * mrSges(7,2) - t397 * t387 + t398 * t388;
t478 = m(6) * t377 - t385 * mrSges(6,1) + t386 * mrSges(6,2) - t420 * t404 + t421 * t405 + t480;
t475 = -m(5) * t399 + t417 * mrSges(5,1) - t418 * mrSges(5,2) + t443 * t425 - t444 * t426 - t478;
t344 = m(4) * t402 + qJDD(3) * mrSges(4,1) - t449 * mrSges(4,3) + qJD(3) * t453 - t446 * t490 + t475;
t485 = t470 * t327 - t344 * t465;
t328 = t332 * t469 + t333 * t464;
t373 = Ifges(7,4) * t398 + Ifges(7,2) * t397 + Ifges(7,6) * t451;
t374 = Ifges(7,1) * t398 + Ifges(7,4) * t397 + Ifges(7,5) * t451;
t479 = -mrSges(7,1) * t351 + mrSges(7,2) * t352 - Ifges(7,5) * t365 - Ifges(7,6) * t364 - Ifges(7,3) * t431 - t398 * t373 + t397 * t374;
t477 = -m(4) * t409 + t450 * mrSges(4,1) - mrSges(4,2) * t449 - t452 * t490 + t453 * t489 - t328;
t392 = Ifges(6,4) * t421 + Ifges(6,2) * t420 + Ifges(6,6) * t454;
t393 = Ifges(6,1) * t421 + Ifges(6,4) * t420 + Ifges(6,5) * t454;
t476 = -mrSges(6,1) * t356 + mrSges(6,2) * t357 - Ifges(6,5) * t386 - Ifges(6,6) * t385 - Ifges(6,3) * t438 - pkin(5) * t341 - t421 * t392 + t420 * t393 + t479;
t412 = Ifges(5,4) * t444 + Ifges(5,2) * t443 + Ifges(5,6) * t455;
t413 = Ifges(5,1) * t444 + Ifges(5,4) * t443 + Ifges(5,5) * t455;
t474 = mrSges(5,1) * t378 - mrSges(5,2) * t379 + Ifges(5,5) * t418 + Ifges(5,6) * t417 + Ifges(5,3) * t442 + pkin(4) * t334 + t444 * t412 - t443 * t413 - t476;
t437 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t465 + Ifges(4,4) * t470) * qJD(1);
t436 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t465 + Ifges(4,2) * t470) * qJD(1);
t411 = Ifges(5,5) * t444 + Ifges(5,6) * t443 + Ifges(5,3) * t455;
t391 = Ifges(6,5) * t421 + Ifges(6,6) * t420 + Ifges(6,3) * t454;
t372 = Ifges(7,5) * t398 + Ifges(7,6) * t397 + Ifges(7,3) * t451;
t343 = mrSges(7,2) * t359 - mrSges(7,3) * t351 + Ifges(7,1) * t365 + Ifges(7,4) * t364 + Ifges(7,5) * t431 + t372 * t397 - t373 * t451;
t342 = -mrSges(7,1) * t359 + mrSges(7,3) * t352 + Ifges(7,4) * t365 + Ifges(7,2) * t364 + Ifges(7,6) * t431 - t372 * t398 + t374 * t451;
t330 = mrSges(6,2) * t377 - mrSges(6,3) * t356 + Ifges(6,1) * t386 + Ifges(6,4) * t385 + Ifges(6,5) * t438 - pkin(10) * t341 - t342 * t462 + t343 * t467 + t391 * t420 - t392 * t454;
t329 = -mrSges(6,1) * t377 + mrSges(6,3) * t357 + Ifges(6,4) * t386 + Ifges(6,2) * t385 + Ifges(6,6) * t438 - pkin(5) * t480 + pkin(10) * t482 + t467 * t342 + t462 * t343 - t421 * t391 + t454 * t393;
t325 = mrSges(5,2) * t399 - mrSges(5,3) * t378 + Ifges(5,1) * t418 + Ifges(5,4) * t417 + Ifges(5,5) * t442 - pkin(9) * t334 - t329 * t463 + t330 * t468 + t411 * t443 - t412 * t455;
t324 = -mrSges(5,1) * t399 + mrSges(5,3) * t379 + Ifges(5,4) * t418 + Ifges(5,2) * t417 + Ifges(5,6) * t442 - pkin(4) * t478 + pkin(9) * t483 + t468 * t329 + t463 * t330 - t444 * t411 + t455 * t413;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t486 - mrSges(2,2) * t481 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t422 - mrSges(3,2) * t423 + t465 * (mrSges(4,2) * t409 - mrSges(4,3) * t402 + Ifges(4,1) * t449 + Ifges(4,4) * t450 + Ifges(4,5) * qJDD(3) - pkin(8) * t328 - qJD(3) * t436 - t464 * t324 + t469 * t325) + t470 * (-mrSges(4,1) * t409 + mrSges(4,3) * t403 + Ifges(4,4) * t449 + Ifges(4,2) * t450 + Ifges(4,6) * qJDD(3) - pkin(3) * t328 + qJD(3) * t437 - t474) + pkin(2) * t477 + pkin(7) * t485 + pkin(1) * (t460 * (m(3) * t423 - mrSges(3,1) * t473 - qJDD(1) * mrSges(3,2) + t485) + t461 * (m(3) * t422 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t473 + t477)); m(3) * t459 + t327 * t465 + t344 * t470; Ifges(4,5) * t449 + Ifges(4,6) * t450 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t402 - mrSges(4,2) * t403 + t464 * t325 + t469 * t324 + pkin(3) * t475 + pkin(8) * t484 + (t436 * t465 - t437 * t470) * qJD(1); t474; -t476; -t479;];
tauJ  = t1;
