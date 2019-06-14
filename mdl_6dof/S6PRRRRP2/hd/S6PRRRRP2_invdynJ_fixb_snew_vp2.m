% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:35:12
% EndTime: 2019-05-05 09:35:16
% DurationCPUTime: 2.30s
% Computational Cost: add. (19774->267), mult. (38631->333), div. (0->0), fcn. (27528->12), ass. (0->114)
t521 = Ifges(6,1) + Ifges(7,1);
t514 = Ifges(6,4) - Ifges(7,5);
t513 = -Ifges(6,5) - Ifges(7,4);
t520 = Ifges(6,2) + Ifges(7,3);
t512 = Ifges(6,6) - Ifges(7,6);
t519 = -Ifges(6,3) - Ifges(7,2);
t476 = sin(pkin(11));
t478 = cos(pkin(11));
t462 = g(1) * t476 - g(2) * t478;
t475 = -g(3) + qJDD(1);
t477 = sin(pkin(6));
t479 = cos(pkin(6));
t518 = t462 * t479 + t475 * t477;
t481 = sin(qJ(4));
t482 = sin(qJ(3));
t484 = cos(qJ(4));
t485 = cos(qJ(3));
t452 = (t481 * t482 - t484 * t485) * qJD(2);
t463 = -g(1) * t478 - g(2) * t476;
t483 = sin(qJ(2));
t486 = cos(qJ(2));
t434 = -t483 * t463 + t486 * t518;
t502 = qJD(2) * qJD(3);
t500 = t485 * t502;
t460 = qJDD(2) * t482 + t500;
t461 = qJDD(2) * t485 - t482 * t502;
t424 = -qJD(4) * t452 + t460 * t484 + t461 * t481;
t453 = (t481 * t485 + t482 * t484) * qJD(2);
t473 = qJD(3) + qJD(4);
t480 = sin(qJ(5));
t516 = cos(qJ(5));
t439 = t453 * t480 - t516 * t473;
t472 = qJDD(3) + qJDD(4);
t397 = -t439 * qJD(5) + t516 * t424 + t480 * t472;
t440 = t516 * t453 + t480 * t473;
t412 = mrSges(7,1) * t439 - mrSges(7,3) * t440;
t435 = t486 * t463 + t518 * t483;
t487 = qJD(2) ^ 2;
t430 = -pkin(2) * t487 + qJDD(2) * pkin(8) + t435;
t445 = -t462 * t477 + t475 * t479;
t407 = -t430 * t482 + t485 * t445;
t392 = (-t460 + t500) * pkin(9) + (t482 * t485 * t487 + qJDD(3)) * pkin(3) + t407;
t408 = t485 * t430 + t482 * t445;
t504 = qJD(2) * t482;
t467 = qJD(3) * pkin(3) - pkin(9) * t504;
t474 = t485 ^ 2;
t393 = -pkin(3) * t474 * t487 + pkin(9) * t461 - qJD(3) * t467 + t408;
t388 = t481 * t392 + t484 * t393;
t438 = pkin(4) * t452 - pkin(10) * t453;
t471 = t473 ^ 2;
t383 = -pkin(4) * t471 + pkin(10) * t472 - t438 * t452 + t388;
t493 = -qJDD(2) * pkin(2) - t434;
t398 = -pkin(3) * t461 + t467 * t504 + (-pkin(9) * t474 - pkin(8)) * t487 + t493;
t423 = -qJD(4) * t453 - t460 * t481 + t461 * t484;
t385 = (t452 * t473 - t424) * pkin(10) + (t453 * t473 - t423) * pkin(4) + t398;
t379 = -t480 * t383 + t516 * t385;
t411 = pkin(5) * t439 - qJ(6) * t440;
t421 = qJDD(5) - t423;
t447 = qJD(5) + t452;
t446 = t447 ^ 2;
t377 = -t421 * pkin(5) - t446 * qJ(6) + t440 * t411 + qJDD(6) - t379;
t425 = -mrSges(7,2) * t439 + mrSges(7,3) * t447;
t496 = -m(7) * t377 + t421 * mrSges(7,1) + t447 * t425;
t373 = mrSges(7,2) * t397 + t412 * t440 - t496;
t380 = t516 * t383 + t480 * t385;
t376 = -pkin(5) * t446 + qJ(6) * t421 + 0.2e1 * qJD(6) * t447 - t411 * t439 + t380;
t396 = qJD(5) * t440 + t424 * t480 - t516 * t472;
t428 = -mrSges(7,1) * t447 + mrSges(7,2) * t440;
t501 = m(7) * t376 + t421 * mrSges(7,3) + t447 * t428;
t506 = t514 * t439 - t521 * t440 + t513 * t447;
t507 = t520 * t439 - t514 * t440 - t512 * t447;
t517 = -t512 * t396 - t513 * t397 - t519 * t421 - t506 * t439 - t507 * t440 + mrSges(6,1) * t379 - mrSges(7,1) * t377 - mrSges(6,2) * t380 + mrSges(7,3) * t376 - pkin(5) * t373 + qJ(6) * (-mrSges(7,2) * t396 - t412 * t439 + t501);
t515 = -mrSges(6,3) - mrSges(7,2);
t437 = mrSges(5,1) * t452 + mrSges(5,2) * t453;
t444 = mrSges(5,1) * t473 - mrSges(5,3) * t453;
t427 = mrSges(6,1) * t447 - mrSges(6,3) * t440;
t505 = -mrSges(6,1) * t439 - mrSges(6,2) * t440 - t412;
t368 = m(6) * t380 - mrSges(6,2) * t421 + t515 * t396 - t427 * t447 + t505 * t439 + t501;
t426 = -mrSges(6,2) * t447 - mrSges(6,3) * t439;
t370 = m(6) * t379 + mrSges(6,1) * t421 + t515 * t397 + t426 * t447 + t505 * t440 + t496;
t497 = t516 * t368 - t370 * t480;
t357 = m(5) * t388 - mrSges(5,2) * t472 + mrSges(5,3) * t423 - t437 * t452 - t444 * t473 + t497;
t387 = t392 * t484 - t481 * t393;
t443 = -mrSges(5,2) * t473 - mrSges(5,3) * t452;
t382 = -pkin(4) * t472 - pkin(10) * t471 + t453 * t438 - t387;
t378 = -0.2e1 * qJD(6) * t440 + (t439 * t447 - t397) * qJ(6) + (t440 * t447 + t396) * pkin(5) + t382;
t374 = m(7) * t378 + mrSges(7,1) * t396 - t397 * mrSges(7,3) + t425 * t439 - t440 * t428;
t489 = -m(6) * t382 - t396 * mrSges(6,1) - mrSges(6,2) * t397 - t439 * t426 - t427 * t440 - t374;
t365 = m(5) * t387 + mrSges(5,1) * t472 - mrSges(5,3) * t424 - t437 * t453 + t443 * t473 + t489;
t354 = t481 * t357 + t484 * t365;
t363 = t480 * t368 + t516 * t370;
t508 = t512 * t439 + t513 * t440 + t519 * t447;
t503 = qJD(2) * t485;
t459 = (-mrSges(4,1) * t485 + mrSges(4,2) * t482) * qJD(2);
t465 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t503;
t352 = m(4) * t407 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t460 + qJD(3) * t465 - t459 * t504 + t354;
t464 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t504;
t498 = t484 * t357 - t365 * t481;
t353 = m(4) * t408 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t461 - qJD(3) * t464 + t459 * t503 + t498;
t499 = -t352 * t482 + t485 * t353;
t492 = m(5) * t398 - t423 * mrSges(5,1) + t424 * mrSges(5,2) + t452 * t443 + t453 * t444 + t363;
t359 = -mrSges(6,1) * t382 - mrSges(7,1) * t378 + mrSges(7,2) * t376 + mrSges(6,3) * t380 - pkin(5) * t374 - t520 * t396 + t514 * t397 + t512 * t421 + t508 * t440 - t506 * t447;
t361 = mrSges(6,2) * t382 + mrSges(7,2) * t377 - mrSges(6,3) * t379 - mrSges(7,3) * t378 - qJ(6) * t374 - t514 * t396 + t521 * t397 - t513 * t421 + t508 * t439 + t507 * t447;
t432 = Ifges(5,4) * t453 - Ifges(5,2) * t452 + Ifges(5,6) * t473;
t433 = Ifges(5,1) * t453 - Ifges(5,4) * t452 + Ifges(5,5) * t473;
t491 = mrSges(5,1) * t387 - mrSges(5,2) * t388 + Ifges(5,5) * t424 + Ifges(5,6) * t423 + Ifges(5,3) * t472 + pkin(4) * t489 + pkin(10) * t497 + t516 * t359 + t480 * t361 + t453 * t432 + t433 * t452;
t429 = -pkin(8) * t487 + t493;
t488 = -m(4) * t429 + t461 * mrSges(4,1) - mrSges(4,2) * t460 - t464 * t504 + t465 * t503 - t492;
t451 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t482 + Ifges(4,4) * t485) * qJD(2);
t450 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t482 + Ifges(4,2) * t485) * qJD(2);
t431 = Ifges(5,5) * t453 - Ifges(5,6) * t452 + Ifges(5,3) * t473;
t350 = -mrSges(5,1) * t398 + mrSges(5,3) * t388 + Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * t472 - pkin(4) * t363 - t453 * t431 + t473 * t433 - t517;
t349 = mrSges(5,2) * t398 - mrSges(5,3) * t387 + Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * t472 - pkin(10) * t363 - t480 * t359 + t516 * t361 - t452 * t431 - t473 * t432;
t1 = [m(2) * t475 + t479 * (m(3) * t445 + t352 * t485 + t353 * t482) + (t483 * (m(3) * t435 - mrSges(3,1) * t487 - qJDD(2) * mrSges(3,2) + t499) + t486 * (m(3) * t434 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t487 + t488)) * t477; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t434 - mrSges(3,2) * t435 + t482 * (mrSges(4,2) * t429 - mrSges(4,3) * t407 + Ifges(4,1) * t460 + Ifges(4,4) * t461 + Ifges(4,5) * qJDD(3) - pkin(9) * t354 - qJD(3) * t450 + t349 * t484 - t350 * t481) + t485 * (-mrSges(4,1) * t429 + mrSges(4,3) * t408 + Ifges(4,4) * t460 + Ifges(4,2) * t461 + Ifges(4,6) * qJDD(3) - pkin(3) * t492 + pkin(9) * t498 + qJD(3) * t451 + t481 * t349 + t484 * t350) + pkin(2) * t488 + pkin(8) * t499; Ifges(4,3) * qJDD(3) + t491 + mrSges(4,1) * t407 - mrSges(4,2) * t408 + Ifges(4,5) * t460 + Ifges(4,6) * t461 + pkin(3) * t354 + (t450 * t482 - t451 * t485) * qJD(2); t491; t517; t373;];
tauJ  = t1;
