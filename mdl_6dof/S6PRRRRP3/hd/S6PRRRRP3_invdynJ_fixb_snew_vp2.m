% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRRRP3
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
% Datum: 2019-05-05 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:43:15
% EndTime: 2019-05-05 09:43:20
% DurationCPUTime: 2.86s
% Computational Cost: add. (22212->268), mult. (43200->328), div. (0->0), fcn. (30807->12), ass. (0->113)
t481 = sin(qJ(4));
t485 = cos(qJ(4));
t482 = sin(qJ(3));
t505 = qJD(2) * t482;
t462 = qJD(3) * t481 + t485 * t505;
t486 = cos(qJ(3));
t503 = qJD(2) * qJD(3);
t500 = t486 * t503;
t465 = qJDD(2) * t482 + t500;
t435 = -qJD(4) * t462 + qJDD(3) * t485 - t465 * t481;
t461 = qJD(3) * t485 - t481 * t505;
t436 = qJD(4) * t461 + qJDD(3) * t481 + t465 * t485;
t480 = sin(qJ(5));
t484 = cos(qJ(5));
t438 = t461 * t484 - t462 * t480;
t396 = qJD(5) * t438 + t435 * t480 + t436 * t484;
t439 = t461 * t480 + t462 * t484;
t415 = -mrSges(7,1) * t438 + mrSges(7,2) * t439;
t476 = sin(pkin(11));
t478 = cos(pkin(11));
t468 = -g(1) * t478 - g(2) * t476;
t483 = sin(qJ(2));
t487 = cos(qJ(2));
t467 = g(1) * t476 - g(2) * t478;
t475 = -g(3) + qJDD(1);
t477 = sin(pkin(6));
t479 = cos(pkin(6));
t512 = t467 * t479 + t475 * t477;
t428 = t487 * t468 + t483 * t512;
t489 = qJD(2) ^ 2;
t421 = -pkin(2) * t489 + qJDD(2) * pkin(8) + t428;
t445 = -t467 * t477 + t475 * t479;
t414 = t486 * t421 + t482 * t445;
t464 = (-pkin(3) * t486 - pkin(9) * t482) * qJD(2);
t488 = qJD(3) ^ 2;
t504 = qJD(2) * t486;
t400 = -pkin(3) * t488 + qJDD(3) * pkin(9) + t464 * t504 + t414;
t427 = -t483 * t468 + t487 * t512;
t420 = -qJDD(2) * pkin(2) - t489 * pkin(8) - t427;
t474 = t482 * t503;
t466 = qJDD(2) * t486 - t474;
t403 = (-t465 - t500) * pkin(9) + (-t466 + t474) * pkin(3) + t420;
t382 = -t481 * t400 + t485 * t403;
t458 = qJDD(4) - t466;
t473 = qJD(4) - t504;
t379 = (t461 * t473 - t436) * pkin(10) + (t461 * t462 + t458) * pkin(4) + t382;
t383 = t485 * t400 + t481 * t403;
t444 = pkin(4) * t473 - pkin(10) * t462;
t457 = t461 ^ 2;
t381 = -pkin(4) * t457 + pkin(10) * t435 - t444 * t473 + t383;
t373 = t484 * t379 - t480 * t381;
t454 = qJDD(5) + t458;
t472 = qJD(5) + t473;
t368 = -0.2e1 * qJD(6) * t439 + (t438 * t472 - t396) * qJ(6) + (t438 * t439 + t454) * pkin(5) + t373;
t422 = -mrSges(7,2) * t472 + mrSges(7,3) * t438;
t502 = m(7) * t368 + t454 * mrSges(7,1) + t472 * t422;
t365 = -t396 * mrSges(7,3) - t439 * t415 + t502;
t374 = t480 * t379 + t484 * t381;
t395 = -qJD(5) * t439 + t435 * t484 - t436 * t480;
t424 = pkin(5) * t472 - qJ(6) * t439;
t437 = t438 ^ 2;
t370 = -pkin(5) * t437 + qJ(6) * t395 + 0.2e1 * qJD(6) * t438 - t424 * t472 + t374;
t510 = Ifges(6,4) + Ifges(7,4);
t516 = Ifges(6,5) + Ifges(7,5);
t517 = Ifges(6,1) + Ifges(7,1);
t511 = t510 * t438 + t517 * t439 + t516 * t472;
t515 = Ifges(6,6) + Ifges(7,6);
t519 = Ifges(6,2) + Ifges(7,2);
t513 = t519 * t438 + t510 * t439 + t515 * t472;
t514 = Ifges(6,3) + Ifges(7,3);
t520 = mrSges(6,1) * t373 + mrSges(7,1) * t368 - mrSges(6,2) * t374 - mrSges(7,2) * t370 + pkin(5) * t365 + t515 * t395 + t396 * t516 - t438 * t511 + t513 * t439 + t514 * t454;
t416 = -mrSges(6,1) * t438 + mrSges(6,2) * t439;
t423 = -mrSges(6,2) * t472 + mrSges(6,3) * t438;
t359 = m(6) * t373 + t454 * mrSges(6,1) + t472 * t423 + (-t415 - t416) * t439 + (-mrSges(6,3) - mrSges(7,3)) * t396 + t502;
t425 = mrSges(7,1) * t472 - mrSges(7,3) * t439;
t426 = mrSges(6,1) * t472 - mrSges(6,3) * t439;
t501 = m(7) * t370 + t395 * mrSges(7,3) + t438 * t415;
t362 = m(6) * t374 + t395 * mrSges(6,3) + t438 * t416 + (-t425 - t426) * t472 + (-mrSges(6,2) - mrSges(7,2)) * t454 + t501;
t357 = t484 * t359 + t480 * t362;
t430 = Ifges(5,4) * t462 + Ifges(5,2) * t461 + Ifges(5,6) * t473;
t431 = Ifges(5,1) * t462 + Ifges(5,4) * t461 + Ifges(5,5) * t473;
t518 = mrSges(5,1) * t382 - mrSges(5,2) * t383 + Ifges(5,5) * t436 + Ifges(5,6) * t435 + Ifges(5,3) * t458 + pkin(4) * t357 + t462 * t430 - t461 * t431 + t520;
t507 = -t515 * t438 - t439 * t516 - t514 * t472;
t463 = (-mrSges(4,1) * t486 + mrSges(4,2) * t482) * qJD(2);
t469 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t505;
t440 = -mrSges(5,1) * t461 + mrSges(5,2) * t462;
t442 = -mrSges(5,2) * t473 + mrSges(5,3) * t461;
t354 = m(5) * t382 + mrSges(5,1) * t458 - mrSges(5,3) * t436 - t440 * t462 + t442 * t473 + t357;
t443 = mrSges(5,1) * t473 - mrSges(5,3) * t462;
t496 = -t359 * t480 + t484 * t362;
t355 = m(5) * t383 - mrSges(5,2) * t458 + mrSges(5,3) * t435 + t440 * t461 - t443 * t473 + t496;
t497 = -t354 * t481 + t485 * t355;
t350 = m(4) * t414 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t466 - qJD(3) * t469 + t463 * t504 + t497;
t413 = -t482 * t421 + t445 * t486;
t470 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t504;
t399 = -qJDD(3) * pkin(3) - pkin(9) * t488 + t464 * t505 - t413;
t384 = -pkin(4) * t435 - pkin(10) * t457 + t462 * t444 + t399;
t376 = -pkin(5) * t395 - qJ(6) * t437 + t424 * t439 + qJDD(6) + t384;
t371 = m(7) * t376 - t395 * mrSges(7,1) + t396 * mrSges(7,2) - t438 * t422 + t439 * t425;
t494 = m(6) * t384 - t395 * mrSges(6,1) + t396 * mrSges(6,2) - t438 * t423 + t439 * t426 + t371;
t491 = -m(5) * t399 + t435 * mrSges(5,1) - t436 * mrSges(5,2) + t461 * t442 - t462 * t443 - t494;
t363 = m(4) * t413 + qJDD(3) * mrSges(4,1) - t465 * mrSges(4,3) + qJD(3) * t470 - t463 * t505 + t491;
t498 = t486 * t350 - t363 * t482;
t351 = t354 * t485 + t355 * t481;
t493 = -m(4) * t420 + t466 * mrSges(4,1) - mrSges(4,2) * t465 - t469 * t505 + t470 * t504 - t351;
t453 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t482 + Ifges(4,4) * t486) * qJD(2);
t452 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t482 + Ifges(4,2) * t486) * qJD(2);
t429 = Ifges(5,5) * t462 + Ifges(5,6) * t461 + Ifges(5,3) * t473;
t356 = mrSges(6,2) * t384 + mrSges(7,2) * t376 - mrSges(6,3) * t373 - mrSges(7,3) * t368 - qJ(6) * t365 + t510 * t395 + t517 * t396 - t507 * t438 + t516 * t454 - t513 * t472;
t352 = -mrSges(6,1) * t384 + mrSges(6,3) * t374 - mrSges(7,1) * t376 + mrSges(7,3) * t370 - pkin(5) * t371 + qJ(6) * t501 + (-qJ(6) * t425 + t511) * t472 + (-mrSges(7,2) * qJ(6) + t515) * t454 + t507 * t439 + t510 * t396 + t519 * t395;
t348 = mrSges(5,2) * t399 - mrSges(5,3) * t382 + Ifges(5,1) * t436 + Ifges(5,4) * t435 + Ifges(5,5) * t458 - pkin(10) * t357 - t352 * t480 + t356 * t484 + t429 * t461 - t430 * t473;
t347 = -mrSges(5,1) * t399 + mrSges(5,3) * t383 + Ifges(5,4) * t436 + Ifges(5,2) * t435 + Ifges(5,6) * t458 - pkin(4) * t494 + pkin(10) * t496 + t484 * t352 + t480 * t356 - t462 * t429 + t473 * t431;
t1 = [m(2) * t475 + t479 * (m(3) * t445 + t350 * t482 + t363 * t486) + (t483 * (m(3) * t428 - mrSges(3,1) * t489 - qJDD(2) * mrSges(3,2) + t498) + t487 * (m(3) * t427 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t489 + t493)) * t477; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t427 - mrSges(3,2) * t428 + t482 * (mrSges(4,2) * t420 - mrSges(4,3) * t413 + Ifges(4,1) * t465 + Ifges(4,4) * t466 + Ifges(4,5) * qJDD(3) - pkin(9) * t351 - qJD(3) * t452 - t347 * t481 + t348 * t485) + t486 * (-mrSges(4,1) * t420 + mrSges(4,3) * t414 + Ifges(4,4) * t465 + Ifges(4,2) * t466 + Ifges(4,6) * qJDD(3) - pkin(3) * t351 + qJD(3) * t453 - t518) + pkin(2) * t493 + pkin(8) * t498; Ifges(4,5) * t465 + Ifges(4,6) * t466 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t413 - mrSges(4,2) * t414 + t481 * t348 + t485 * t347 + pkin(3) * t491 + pkin(9) * t497 + (t452 * t482 - t453 * t486) * qJD(2); t518; t520; t371;];
tauJ  = t1;
