% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Datum: 2019-05-05 13:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:54:56
% EndTime: 2019-05-05 13:54:59
% DurationCPUTime: 2.86s
% Computational Cost: add. (28458->269), mult. (65842->347), div. (0->0), fcn. (46749->12), ass. (0->120)
t477 = qJD(1) ^ 2;
t503 = cos(qJ(4));
t469 = cos(pkin(10));
t502 = pkin(3) * t469;
t466 = sin(pkin(10));
t501 = mrSges(4,2) * t466;
t463 = t469 ^ 2;
t500 = t463 * t477;
t473 = sin(qJ(1));
t475 = cos(qJ(1));
t490 = t473 * g(1) - g(2) * t475;
t450 = qJDD(1) * pkin(1) + t490;
t485 = -g(1) * t475 - g(2) * t473;
t451 = -pkin(1) * t477 + t485;
t467 = sin(pkin(9));
t470 = cos(pkin(9));
t434 = t467 * t450 + t470 * t451;
t425 = -pkin(2) * t477 + qJDD(1) * qJ(3) + t434;
t464 = -g(3) + qJDD(2);
t494 = qJD(1) * qJD(3);
t498 = t469 * t464 - 0.2e1 * t466 * t494;
t405 = (-pkin(7) * qJDD(1) + t477 * t502 - t425) * t466 + t498;
t411 = t466 * t464 + (t425 + 0.2e1 * t494) * t469;
t492 = qJDD(1) * t469;
t406 = -pkin(3) * t500 + pkin(7) * t492 + t411;
t472 = sin(qJ(4));
t390 = t472 * t405 + t503 * t406;
t491 = t469 * t503;
t497 = qJD(1) * t466;
t443 = -qJD(1) * t491 + t472 * t497;
t482 = t503 * t466 + t469 * t472;
t444 = t482 * qJD(1);
t428 = mrSges(5,1) * t443 + mrSges(5,2) * t444;
t493 = qJDD(1) * t466;
t496 = qJD(4) * t444;
t431 = -qJDD(1) * t491 + t472 * t493 + t496;
t441 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t444;
t427 = pkin(4) * t443 - qJ(5) * t444;
t476 = qJD(4) ^ 2;
t382 = -pkin(4) * t476 + qJDD(4) * qJ(5) - t427 * t443 + t390;
t462 = t466 ^ 2;
t433 = t450 * t470 - t467 * t451;
t484 = qJDD(3) - t433;
t407 = (-pkin(2) - t502) * qJDD(1) + (-qJ(3) + (-t462 - t463) * pkin(7)) * t477 + t484;
t495 = t443 * qJD(4);
t432 = t482 * qJDD(1) - t495;
t385 = (-t432 + t495) * qJ(5) + (t431 + t496) * pkin(4) + t407;
t465 = sin(pkin(11));
t468 = cos(pkin(11));
t439 = qJD(4) * t465 + t444 * t468;
t377 = -0.2e1 * qJD(5) * t439 - t382 * t465 + t468 * t385;
t419 = qJDD(4) * t465 + t432 * t468;
t438 = qJD(4) * t468 - t444 * t465;
t375 = (t438 * t443 - t419) * pkin(8) + (t438 * t439 + t431) * pkin(5) + t377;
t378 = 0.2e1 * qJD(5) * t438 + t468 * t382 + t465 * t385;
t417 = pkin(5) * t443 - pkin(8) * t439;
t418 = qJDD(4) * t468 - t432 * t465;
t437 = t438 ^ 2;
t376 = -pkin(5) * t437 + pkin(8) * t418 - t417 * t443 + t378;
t471 = sin(qJ(6));
t474 = cos(qJ(6));
t373 = t375 * t474 - t376 * t471;
t408 = t438 * t474 - t439 * t471;
t388 = qJD(6) * t408 + t418 * t471 + t419 * t474;
t409 = t438 * t471 + t439 * t474;
t395 = -mrSges(7,1) * t408 + mrSges(7,2) * t409;
t442 = qJD(6) + t443;
t396 = -mrSges(7,2) * t442 + mrSges(7,3) * t408;
t430 = qJDD(6) + t431;
t370 = m(7) * t373 + mrSges(7,1) * t430 - mrSges(7,3) * t388 - t395 * t409 + t396 * t442;
t374 = t375 * t471 + t376 * t474;
t387 = -qJD(6) * t409 + t418 * t474 - t419 * t471;
t397 = mrSges(7,1) * t442 - mrSges(7,3) * t409;
t371 = m(7) * t374 - mrSges(7,2) * t430 + mrSges(7,3) * t387 + t395 * t408 - t397 * t442;
t362 = t474 * t370 + t471 * t371;
t412 = -mrSges(6,1) * t438 + mrSges(6,2) * t439;
t415 = -mrSges(6,2) * t443 + mrSges(6,3) * t438;
t360 = m(6) * t377 + mrSges(6,1) * t431 - mrSges(6,3) * t419 - t412 * t439 + t415 * t443 + t362;
t416 = mrSges(6,1) * t443 - mrSges(6,3) * t439;
t486 = -t370 * t471 + t474 * t371;
t361 = m(6) * t378 - mrSges(6,2) * t431 + mrSges(6,3) * t418 + t412 * t438 - t416 * t443 + t486;
t487 = -t360 * t465 + t468 * t361;
t355 = m(5) * t390 - qJDD(4) * mrSges(5,2) - mrSges(5,3) * t431 - qJD(4) * t441 - t428 * t443 + t487;
t389 = t503 * t405 - t472 * t406;
t381 = -qJDD(4) * pkin(4) - t476 * qJ(5) + t444 * t427 + qJDD(5) - t389;
t379 = -t418 * pkin(5) - t437 * pkin(8) + t439 * t417 + t381;
t481 = m(7) * t379 - t387 * mrSges(7,1) + mrSges(7,2) * t388 - t408 * t396 + t397 * t409;
t372 = m(6) * t381 - t418 * mrSges(6,1) + mrSges(6,2) * t419 - t438 * t415 + t416 * t439 + t481;
t440 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t443;
t366 = m(5) * t389 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t432 + qJD(4) * t440 - t428 * t444 - t372;
t499 = t472 * t355 + t503 * t366;
t356 = t468 * t360 + t465 * t361;
t410 = -t425 * t466 + t498;
t483 = mrSges(4,3) * qJDD(1) + t477 * (-mrSges(4,1) * t469 + t501);
t348 = m(4) * t410 - t483 * t466 + t499;
t488 = t503 * t355 - t366 * t472;
t349 = m(4) * t411 + t483 * t469 + t488;
t489 = -t348 * t466 + t469 * t349;
t480 = m(5) * t407 + t431 * mrSges(5,1) + mrSges(5,2) * t432 + t443 * t440 + t441 * t444 + t356;
t421 = -qJDD(1) * pkin(2) - qJ(3) * t477 + t484;
t479 = -m(4) * t421 + mrSges(4,1) * t492 - t480 + (t462 * t477 + t500) * mrSges(4,3);
t392 = Ifges(7,4) * t409 + Ifges(7,2) * t408 + Ifges(7,6) * t442;
t393 = Ifges(7,1) * t409 + Ifges(7,4) * t408 + Ifges(7,5) * t442;
t478 = mrSges(7,1) * t373 - mrSges(7,2) * t374 + Ifges(7,5) * t388 + Ifges(7,6) * t387 + Ifges(7,3) * t430 + t409 * t392 - t408 * t393;
t449 = (Ifges(4,5) * t466 + Ifges(4,6) * t469) * qJD(1);
t424 = Ifges(5,1) * t444 - Ifges(5,4) * t443 + Ifges(5,5) * qJD(4);
t423 = Ifges(5,4) * t444 - Ifges(5,2) * t443 + Ifges(5,6) * qJD(4);
t422 = Ifges(5,5) * t444 - Ifges(5,6) * t443 + Ifges(5,3) * qJD(4);
t401 = Ifges(6,1) * t439 + Ifges(6,4) * t438 + Ifges(6,5) * t443;
t400 = Ifges(6,4) * t439 + Ifges(6,2) * t438 + Ifges(6,6) * t443;
t399 = Ifges(6,5) * t439 + Ifges(6,6) * t438 + Ifges(6,3) * t443;
t391 = Ifges(7,5) * t409 + Ifges(7,6) * t408 + Ifges(7,3) * t442;
t364 = mrSges(7,2) * t379 - mrSges(7,3) * t373 + Ifges(7,1) * t388 + Ifges(7,4) * t387 + Ifges(7,5) * t430 + t391 * t408 - t392 * t442;
t363 = -mrSges(7,1) * t379 + mrSges(7,3) * t374 + Ifges(7,4) * t388 + Ifges(7,2) * t387 + Ifges(7,6) * t430 - t391 * t409 + t393 * t442;
t352 = mrSges(4,2) * t493 - t479;
t351 = mrSges(6,2) * t381 - mrSges(6,3) * t377 + Ifges(6,1) * t419 + Ifges(6,4) * t418 + Ifges(6,5) * t431 - pkin(8) * t362 - t363 * t471 + t364 * t474 + t399 * t438 - t400 * t443;
t350 = -mrSges(6,1) * t381 + mrSges(6,3) * t378 + Ifges(6,4) * t419 + Ifges(6,2) * t418 + Ifges(6,6) * t431 - pkin(5) * t481 + pkin(8) * t486 + t474 * t363 + t471 * t364 - t439 * t399 + t443 * t401;
t346 = Ifges(5,6) * qJDD(4) + (-Ifges(5,2) - Ifges(6,3)) * t431 - t444 * t422 + t438 * t401 - t439 * t400 + Ifges(5,4) * t432 + qJD(4) * t424 - Ifges(6,6) * t418 - Ifges(6,5) * t419 - mrSges(5,1) * t407 + mrSges(5,3) * t390 - mrSges(6,1) * t377 + mrSges(6,2) * t378 - pkin(5) * t362 - pkin(4) * t356 - t478;
t345 = mrSges(5,2) * t407 - mrSges(5,3) * t389 + Ifges(5,1) * t432 - Ifges(5,4) * t431 + Ifges(5,5) * qJDD(4) - qJ(5) * t356 - qJD(4) * t423 - t350 * t465 + t351 * t468 - t422 * t443;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t490 - mrSges(2,2) * t485 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t433 - mrSges(3,2) * t434 + t466 * (t469 * qJD(1) * t449 + mrSges(4,2) * t421 - mrSges(4,3) * t410 + t503 * t345 - t472 * t346 - pkin(7) * t499 + (Ifges(4,1) * t466 + Ifges(4,4) * t469) * qJDD(1)) + t469 * (-t449 * t497 - mrSges(4,1) * t421 + mrSges(4,3) * t411 + t472 * t345 + t503 * t346 - pkin(3) * t480 + pkin(7) * t488 + (Ifges(4,4) * t466 + Ifges(4,2) * t469) * qJDD(1)) - pkin(2) * t352 + qJ(3) * t489 + pkin(1) * (t467 * (m(3) * t434 - mrSges(3,1) * t477 - qJDD(1) * mrSges(3,2) + t489) + t470 * (t479 + m(3) * t433 - mrSges(3,2) * t477 + (mrSges(3,1) - t501) * qJDD(1))); m(3) * t464 + t348 * t469 + t349 * t466; t352; mrSges(5,1) * t389 - mrSges(5,2) * t390 + Ifges(5,5) * t432 - Ifges(5,6) * t431 + Ifges(5,3) * qJDD(4) - pkin(4) * t372 + qJ(5) * t487 + t468 * t350 + t465 * t351 + t444 * t423 + t443 * t424; t372; t478;];
tauJ  = t1;
