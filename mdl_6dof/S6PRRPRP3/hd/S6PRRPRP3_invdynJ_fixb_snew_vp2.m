% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-05-05 03:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:52:19
% EndTime: 2019-05-05 03:52:23
% DurationCPUTime: 2.41s
% Computational Cost: add. (18114->266), mult. (37233->329), div. (0->0), fcn. (26195->12), ass. (0->111)
t502 = Ifges(6,1) + Ifges(7,1);
t491 = Ifges(6,4) - Ifges(7,5);
t497 = -Ifges(6,5) - Ifges(7,4);
t501 = Ifges(6,2) + Ifges(7,3);
t489 = Ifges(6,6) - Ifges(7,6);
t457 = sin(pkin(11));
t460 = cos(pkin(11));
t464 = sin(qJ(3));
t482 = qJD(2) * t464;
t438 = qJD(3) * t460 - t457 * t482;
t439 = qJD(3) * t457 + t460 * t482;
t463 = sin(qJ(5));
t493 = cos(qJ(5));
t414 = -t438 * t493 + t463 * t439;
t466 = cos(qJ(3));
t480 = qJD(2) * qJD(3);
t478 = t466 * t480;
t446 = qJDD(2) * t464 + t478;
t422 = qJDD(3) * t460 - t446 * t457;
t423 = qJDD(3) * t457 + t446 * t460;
t382 = -t414 * qJD(5) + t463 * t422 + t423 * t493;
t415 = t463 * t438 + t439 * t493;
t395 = mrSges(7,1) * t414 - mrSges(7,3) * t415;
t458 = sin(pkin(10));
t461 = cos(pkin(10));
t449 = -g(1) * t461 - g(2) * t458;
t465 = sin(qJ(2));
t467 = cos(qJ(2));
t448 = g(1) * t458 - g(2) * t461;
t456 = -g(3) + qJDD(1);
t459 = sin(pkin(6));
t462 = cos(pkin(6));
t499 = t448 * t462 + t456 * t459;
t408 = t467 * t449 + t465 * t499;
t469 = qJD(2) ^ 2;
t402 = -pkin(2) * t469 + qJDD(2) * pkin(8) + t408;
t425 = -t448 * t459 + t456 * t462;
t393 = t466 * t402 + t464 * t425;
t444 = (-pkin(3) * t466 - qJ(4) * t464) * qJD(2);
t468 = qJD(3) ^ 2;
t481 = qJD(2) * t466;
t374 = -pkin(3) * t468 + qJDD(3) * qJ(4) + t444 * t481 + t393;
t407 = -t465 * t449 + t467 * t499;
t401 = -qJDD(2) * pkin(2) - t469 * pkin(8) - t407;
t455 = t464 * t480;
t447 = qJDD(2) * t466 - t455;
t380 = (-t446 - t478) * qJ(4) + (-t447 + t455) * pkin(3) + t401;
t369 = -0.2e1 * qJD(4) * t439 - t457 * t374 + t460 * t380;
t366 = (-t438 * t481 - t423) * pkin(9) + (t438 * t439 - t447) * pkin(4) + t369;
t370 = 0.2e1 * qJD(4) * t438 + t460 * t374 + t457 * t380;
t424 = -pkin(4) * t481 - pkin(9) * t439;
t437 = t438 ^ 2;
t368 = -pkin(4) * t437 + pkin(9) * t422 + t424 * t481 + t370;
t361 = t493 * t366 - t463 * t368;
t394 = pkin(5) * t414 - qJ(6) * t415;
t441 = qJDD(5) - t447;
t454 = qJD(5) - t481;
t453 = t454 ^ 2;
t360 = -t441 * pkin(5) - t453 * qJ(6) + t415 * t394 + qJDD(6) - t361;
t406 = -mrSges(7,2) * t414 + mrSges(7,3) * t454;
t474 = -m(7) * t360 + t441 * mrSges(7,1) + t454 * t406;
t356 = t382 * mrSges(7,2) + t415 * t395 - t474;
t362 = t463 * t366 + t493 * t368;
t359 = -pkin(5) * t453 + qJ(6) * t441 + 0.2e1 * qJD(6) * t454 - t394 * t414 + t362;
t381 = t415 * qJD(5) - t422 * t493 + t463 * t423;
t405 = -mrSges(7,1) * t454 + mrSges(7,2) * t415;
t479 = m(7) * t359 + t441 * mrSges(7,3) + t454 * t405;
t485 = t501 * t414 - t491 * t415 - t489 * t454;
t494 = t491 * t414 - t502 * t415 + t497 * t454;
t496 = -Ifges(6,3) - Ifges(7,2);
t500 = t497 * t382 + t494 * t414 + t489 * t381 + t496 * t441 - mrSges(6,1) * t361 + mrSges(7,1) * t360 + mrSges(6,2) * t362 - mrSges(7,3) * t359 + pkin(5) * t356 - qJ(6) * (-t381 * mrSges(7,2) - t414 * t395 + t479) + t485 * t415;
t492 = -mrSges(6,3) - mrSges(7,2);
t404 = mrSges(6,1) * t454 - mrSges(6,3) * t415;
t483 = -mrSges(6,1) * t414 - mrSges(6,2) * t415 - t395;
t350 = m(6) * t362 - t441 * mrSges(6,2) + t492 * t381 - t454 * t404 + t483 * t414 + t479;
t403 = -mrSges(6,2) * t454 - mrSges(6,3) * t414;
t352 = m(6) * t361 + t441 * mrSges(6,1) + t492 * t382 + t454 * t403 + t483 * t415 + t474;
t345 = t463 * t350 + t493 * t352;
t486 = t489 * t414 + t415 * t497 + t496 * t454;
t445 = (-mrSges(4,1) * t466 + mrSges(4,2) * t464) * qJD(2);
t450 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t482;
t416 = -mrSges(5,1) * t438 + mrSges(5,2) * t439;
t420 = mrSges(5,2) * t481 + mrSges(5,3) * t438;
t343 = m(5) * t369 - mrSges(5,1) * t447 - mrSges(5,3) * t423 - t416 * t439 - t420 * t481 + t345;
t421 = -mrSges(5,1) * t481 - mrSges(5,3) * t439;
t475 = t493 * t350 - t352 * t463;
t344 = m(5) * t370 + mrSges(5,2) * t447 + mrSges(5,3) * t422 + t416 * t438 + t421 * t481 + t475;
t476 = -t343 * t457 + t460 * t344;
t340 = m(4) * t393 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t447 - qJD(3) * t450 + t445 * t481 + t476;
t392 = -t464 * t402 + t466 * t425;
t373 = -qJDD(3) * pkin(3) - t468 * qJ(4) + t444 * t482 + qJDD(4) - t392;
t371 = -t422 * pkin(4) - t437 * pkin(9) + t439 * t424 + t373;
t364 = -0.2e1 * qJD(6) * t415 + (t414 * t454 - t382) * qJ(6) + (t415 * t454 + t381) * pkin(5) + t371;
t357 = m(7) * t364 + t381 * mrSges(7,1) - t382 * mrSges(7,3) - t415 * t405 + t414 * t406;
t472 = m(6) * t371 + t381 * mrSges(6,1) + t382 * mrSges(6,2) + t414 * t403 + t415 * t404 + t357;
t354 = m(5) * t373 - t422 * mrSges(5,1) + t423 * mrSges(5,2) - t438 * t420 + t439 * t421 + t472;
t451 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t481;
t353 = m(4) * t392 + qJDD(3) * mrSges(4,1) - t446 * mrSges(4,3) + qJD(3) * t451 - t445 * t482 - t354;
t477 = t466 * t340 - t353 * t464;
t341 = t343 * t460 + t344 * t457;
t470 = -m(4) * t401 + t447 * mrSges(4,1) - mrSges(4,2) * t446 - t450 * t482 + t451 * t481 - t341;
t433 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t464 + Ifges(4,4) * t466) * qJD(2);
t432 = Ifges(4,6) * qJD(3) + (t464 * Ifges(4,4) + Ifges(4,2) * t466) * qJD(2);
t411 = Ifges(5,1) * t439 + Ifges(5,4) * t438 - Ifges(5,5) * t481;
t410 = Ifges(5,4) * t439 + Ifges(5,2) * t438 - Ifges(5,6) * t481;
t409 = Ifges(5,5) * t439 + Ifges(5,6) * t438 - Ifges(5,3) * t481;
t347 = mrSges(6,2) * t371 + mrSges(7,2) * t360 - mrSges(6,3) * t361 - mrSges(7,3) * t364 - qJ(6) * t357 - t491 * t381 + t502 * t382 + t486 * t414 - t497 * t441 + t485 * t454;
t346 = -mrSges(6,1) * t371 - mrSges(7,1) * t364 + mrSges(7,2) * t359 + mrSges(6,3) * t362 - pkin(5) * t357 - t501 * t381 + t491 * t382 + t486 * t415 + t489 * t441 - t494 * t454;
t338 = mrSges(5,2) * t373 - mrSges(5,3) * t369 + Ifges(5,1) * t423 + Ifges(5,4) * t422 - Ifges(5,5) * t447 - pkin(9) * t345 - t463 * t346 + t347 * t493 + t438 * t409 + t410 * t481;
t337 = -mrSges(5,1) * t373 + mrSges(5,3) * t370 + Ifges(5,4) * t423 + Ifges(5,2) * t422 - Ifges(5,6) * t447 - pkin(4) * t472 + pkin(9) * t475 + t346 * t493 + t463 * t347 - t439 * t409 - t411 * t481;
t1 = [m(2) * t456 + t462 * (m(3) * t425 + t340 * t464 + t353 * t466) + (t465 * (m(3) * t408 - mrSges(3,1) * t469 - qJDD(2) * mrSges(3,2) + t477) + t467 * (m(3) * t407 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t469 + t470)) * t459; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t407 - mrSges(3,2) * t408 + t464 * (mrSges(4,2) * t401 - mrSges(4,3) * t392 + Ifges(4,1) * t446 + Ifges(4,4) * t447 + Ifges(4,5) * qJDD(3) - qJ(4) * t341 - qJD(3) * t432 - t337 * t457 + t338 * t460) + t466 * (Ifges(4,6) * qJDD(3) + Ifges(4,4) * t446 + t438 * t411 - t439 * t410 - Ifges(5,6) * t422 - Ifges(5,5) * t423 + qJD(3) * t433 + mrSges(4,3) * t393 - mrSges(4,1) * t401 - mrSges(5,1) * t369 + mrSges(5,2) * t370 - pkin(4) * t345 - pkin(3) * t341 + (Ifges(5,3) + Ifges(4,2)) * t447 + t500) + pkin(2) * t470 + pkin(8) * t477; Ifges(4,5) * t446 + Ifges(4,6) * t447 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t392 - mrSges(4,2) * t393 + t457 * t338 + t460 * t337 - pkin(3) * t354 + qJ(4) * t476 + (t432 * t464 - t433 * t466) * qJD(2); t354; -t500; t356;];
tauJ  = t1;
