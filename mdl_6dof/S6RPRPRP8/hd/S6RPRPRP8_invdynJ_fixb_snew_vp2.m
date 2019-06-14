% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
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
% Datum: 2019-05-05 18:03
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:00:44
% EndTime: 2019-05-05 18:00:48
% DurationCPUTime: 1.89s
% Computational Cost: add. (13562->265), mult. (29319->324), div. (0->0), fcn. (18729->8), ass. (0->105)
t505 = -2 * qJD(4);
t504 = Ifges(6,1) + Ifges(7,1);
t497 = Ifges(6,4) - Ifges(7,5);
t496 = -Ifges(6,5) - Ifges(7,4);
t503 = Ifges(6,2) + Ifges(7,3);
t495 = Ifges(6,6) - Ifges(7,6);
t502 = -Ifges(6,3) - Ifges(7,2);
t469 = qJD(1) ^ 2;
t465 = sin(qJ(1));
t467 = cos(qJ(1));
t482 = g(1) * t465 - t467 * g(2);
t472 = -qJ(2) * t469 + qJDD(2) - t482;
t500 = -pkin(1) - pkin(7);
t434 = t500 * qJDD(1) + t472;
t464 = sin(qJ(3));
t466 = cos(qJ(3));
t425 = t464 * g(3) + t466 * t434;
t486 = qJD(1) * qJD(3);
t483 = t464 * t486;
t449 = qJDD(1) * t466 - t483;
t401 = (-t449 - t483) * qJ(4) + (-t464 * t466 * t469 + qJDD(3)) * pkin(3) + t425;
t426 = -g(3) * t466 + t464 * t434;
t448 = -qJDD(1) * t464 - t466 * t486;
t488 = qJD(1) * t466;
t451 = qJD(3) * pkin(3) - qJ(4) * t488;
t460 = t464 ^ 2;
t402 = -pkin(3) * t460 * t469 + qJ(4) * t448 - qJD(3) * t451 + t426;
t461 = sin(pkin(9));
t462 = cos(pkin(9));
t489 = qJD(1) * t464;
t442 = -t461 * t489 + t462 * t488;
t382 = t401 * t462 - t461 * t402 + t442 * t505;
t477 = -g(1) * t467 - g(2) * t465;
t473 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t477;
t441 = (t461 * t466 + t462 * t464) * qJD(1);
t424 = t448 * t461 + t449 * t462;
t463 = sin(qJ(5));
t499 = cos(qJ(5));
t427 = -t499 * qJD(3) + t442 * t463;
t395 = -t427 * qJD(5) + t463 * qJDD(3) + t499 * t424;
t428 = t463 * qJD(3) + t499 * t442;
t406 = mrSges(7,1) * t427 - mrSges(7,3) * t428;
t383 = t461 * t401 + t462 * t402 + t441 * t505;
t418 = pkin(4) * t441 - pkin(8) * t442;
t468 = qJD(3) ^ 2;
t379 = -pkin(4) * t468 + qJDD(3) * pkin(8) - t418 * t441 + t383;
t404 = -pkin(3) * t448 + qJDD(4) + t451 * t488 + (-qJ(4) * t460 + t500) * t469 + t473;
t423 = t448 * t462 - t449 * t461;
t381 = (qJD(3) * t441 - t424) * pkin(8) + (qJD(3) * t442 - t423) * pkin(4) + t404;
t375 = -t463 * t379 + t499 * t381;
t405 = pkin(5) * t427 - qJ(6) * t428;
t422 = qJDD(5) - t423;
t439 = qJD(5) + t441;
t438 = t439 ^ 2;
t373 = -t422 * pkin(5) - t438 * qJ(6) + t428 * t405 + qJDD(6) - t375;
t409 = -mrSges(7,2) * t427 + mrSges(7,3) * t439;
t478 = -m(7) * t373 + t422 * mrSges(7,1) + t439 * t409;
t369 = mrSges(7,2) * t395 + t406 * t428 - t478;
t376 = t499 * t379 + t463 * t381;
t372 = -pkin(5) * t438 + qJ(6) * t422 + 0.2e1 * qJD(6) * t439 - t405 * t427 + t376;
t394 = qJD(5) * t428 - t499 * qJDD(3) + t424 * t463;
t412 = -mrSges(7,1) * t439 + mrSges(7,2) * t428;
t484 = m(7) * t372 + t422 * mrSges(7,3) + t439 * t412;
t491 = t497 * t427 - t504 * t428 + t496 * t439;
t492 = t503 * t427 - t497 * t428 - t495 * t439;
t501 = -t495 * t394 - t496 * t395 - t502 * t422 - t491 * t427 - t492 * t428 + mrSges(6,1) * t375 - mrSges(7,1) * t373 - mrSges(6,2) * t376 + mrSges(7,3) * t372 - pkin(5) * t369 + qJ(6) * (-mrSges(7,2) * t394 - t406 * t427 + t484);
t498 = -mrSges(6,3) - mrSges(7,2);
t417 = mrSges(5,1) * t441 + mrSges(5,2) * t442;
t433 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t442;
t411 = mrSges(6,1) * t439 - mrSges(6,3) * t428;
t490 = -mrSges(6,1) * t427 - mrSges(6,2) * t428 - t406;
t365 = m(6) * t376 - mrSges(6,2) * t422 + t498 * t394 - t411 * t439 + t490 * t427 + t484;
t410 = -mrSges(6,2) * t439 - mrSges(6,3) * t427;
t367 = m(6) * t375 + mrSges(6,1) * t422 + t498 * t395 + t410 * t439 + t490 * t428 + t478;
t480 = t499 * t365 - t367 * t463;
t356 = m(5) * t383 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t423 - qJD(3) * t433 - t417 * t441 + t480;
t432 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t441;
t378 = -qJDD(3) * pkin(4) - pkin(8) * t468 + t442 * t418 - t382;
t374 = -0.2e1 * qJD(6) * t428 + (t427 * t439 - t395) * qJ(6) + (t428 * t439 + t394) * pkin(5) + t378;
t370 = m(7) * t374 + mrSges(7,1) * t394 - t395 * mrSges(7,3) + t409 * t427 - t428 * t412;
t470 = -m(6) * t378 - t394 * mrSges(6,1) - mrSges(6,2) * t395 - t427 * t410 - t411 * t428 - t370;
t362 = m(5) * t382 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t424 + qJD(3) * t432 - t417 * t442 + t470;
t353 = t461 * t356 + t462 * t362;
t360 = t463 * t365 + t499 * t367;
t493 = t495 * t427 + t496 * t428 + t502 * t439;
t481 = t462 * t356 - t362 * t461;
t447 = (t464 * mrSges(4,1) + t466 * mrSges(4,2)) * qJD(1);
t450 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t489;
t452 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t488;
t476 = (m(4) * t425 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t449 + qJD(3) * t450 - t447 * t488 + t353) * t466 + (m(4) * t426 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t448 - qJD(3) * t452 - t447 * t489 + t481) * t464;
t359 = m(5) * t404 - mrSges(5,1) * t423 + t424 * mrSges(5,2) + t432 * t441 + t442 * t433 + t360;
t445 = Ifges(4,5) * qJD(3) + (t466 * Ifges(4,1) - t464 * Ifges(4,4)) * qJD(1);
t444 = Ifges(4,6) * qJD(3) + (t466 * Ifges(4,4) - t464 * Ifges(4,2)) * qJD(1);
t440 = -qJDD(1) * pkin(1) + t472;
t435 = pkin(1) * t469 - t473;
t431 = t500 * t469 + t473;
t415 = Ifges(5,1) * t442 - Ifges(5,4) * t441 + Ifges(5,5) * qJD(3);
t414 = Ifges(5,4) * t442 - Ifges(5,2) * t441 + Ifges(5,6) * qJD(3);
t413 = Ifges(5,5) * t442 - Ifges(5,6) * t441 + Ifges(5,3) * qJD(3);
t358 = mrSges(6,2) * t378 + mrSges(7,2) * t373 - mrSges(6,3) * t375 - mrSges(7,3) * t374 - qJ(6) * t370 - t497 * t394 + t504 * t395 - t496 * t422 + t493 * t427 + t492 * t439;
t357 = -mrSges(6,1) * t378 - mrSges(7,1) * t374 + mrSges(7,2) * t372 + mrSges(6,3) * t376 - pkin(5) * t370 - t503 * t394 + t497 * t395 + t495 * t422 + t493 * t428 - t491 * t439;
t350 = -mrSges(5,1) * t404 + mrSges(5,3) * t383 + Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * qJDD(3) - pkin(4) * t360 + qJD(3) * t415 - t442 * t413 - t501;
t349 = mrSges(5,2) * t404 - mrSges(5,3) * t382 + Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * qJDD(3) - pkin(8) * t360 - qJD(3) * t414 - t463 * t357 + t499 * t358 - t441 * t413;
t348 = m(3) * t440 + qJDD(1) * mrSges(3,2) - mrSges(3,3) * t469 + t476;
t1 = [mrSges(2,1) * t482 - mrSges(2,2) * t477 + mrSges(3,2) * t440 - mrSges(3,3) * t435 + t466 * (mrSges(4,2) * t431 - mrSges(4,3) * t425 + Ifges(4,1) * t449 + Ifges(4,4) * t448 + Ifges(4,5) * qJDD(3) - qJ(4) * t353 - qJD(3) * t444 + t462 * t349 - t350 * t461) - t464 * (-mrSges(4,1) * t431 + mrSges(4,3) * t426 + Ifges(4,4) * t449 + Ifges(4,2) * t448 + Ifges(4,6) * qJDD(3) - pkin(3) * t359 + qJ(4) * t481 + qJD(3) * t445 + t461 * t349 + t462 * t350) - pkin(7) * t476 - pkin(1) * t348 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t435 + m(4) * t431 - mrSges(4,1) * t448 + mrSges(3,2) * t469 + mrSges(4,2) * t449 + t359 + qJDD(1) * mrSges(3,3) + (t450 * t464 + t452 * t466) * qJD(1)) * qJ(2); t348; Ifges(4,5) * t449 + Ifges(4,6) * t448 + mrSges(4,1) * t425 - mrSges(4,2) * t426 + Ifges(5,5) * t424 + Ifges(5,6) * t423 + t442 * t414 + t441 * t415 + mrSges(5,1) * t382 - mrSges(5,2) * t383 + t463 * t358 + t499 * t357 + pkin(4) * t470 + pkin(8) * t480 + pkin(3) * t353 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t466 * t444 + t464 * t445) * qJD(1); t359; t501; t369;];
tauJ  = t1;
