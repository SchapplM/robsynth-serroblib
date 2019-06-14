% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRRRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:13:02
% EndTime: 2019-05-06 01:13:07
% DurationCPUTime: 2.84s
% Computational Cost: add. (21555->270), mult. (41884->327), div. (0->0), fcn. (27643->10), ass. (0->111)
t464 = sin(qJ(4));
t468 = cos(qJ(4));
t465 = sin(qJ(3));
t489 = qJD(1) * t465;
t446 = qJD(3) * t464 + t468 * t489;
t469 = cos(qJ(3));
t487 = qJD(1) * qJD(3);
t484 = t469 * t487;
t451 = qJDD(1) * t465 + t484;
t419 = -qJD(4) * t446 + qJDD(3) * t468 - t451 * t464;
t445 = qJD(3) * t468 - t464 * t489;
t420 = qJD(4) * t445 + qJDD(3) * t464 + t451 * t468;
t463 = sin(qJ(5));
t467 = cos(qJ(5));
t422 = t445 * t467 - t446 * t463;
t382 = qJD(5) * t422 + t419 * t463 + t420 * t467;
t423 = t445 * t463 + t446 * t467;
t399 = -mrSges(7,1) * t422 + mrSges(7,2) * t423;
t466 = sin(qJ(1));
t470 = cos(qJ(1));
t483 = t466 * g(1) - g(2) * t470;
t447 = qJDD(1) * pkin(1) + t483;
t472 = qJD(1) ^ 2;
t478 = -g(1) * t470 - g(2) * t466;
t449 = -pkin(1) * t472 + t478;
t461 = sin(pkin(10));
t462 = cos(pkin(10));
t424 = t462 * t447 - t461 * t449;
t411 = -qJDD(1) * pkin(2) - t472 * pkin(7) - t424;
t458 = t465 * t487;
t452 = qJDD(1) * t469 - t458;
t392 = (-t451 - t484) * pkin(8) + (-t452 + t458) * pkin(3) + t411;
t425 = t461 * t447 + t462 * t449;
t412 = -pkin(2) * t472 + qJDD(1) * pkin(7) + t425;
t460 = -g(3) + qJDD(2);
t403 = t469 * t412 + t465 * t460;
t450 = (-pkin(3) * t469 - pkin(8) * t465) * qJD(1);
t471 = qJD(3) ^ 2;
t488 = qJD(1) * t469;
t398 = -pkin(3) * t471 + qJDD(3) * pkin(8) + t450 * t488 + t403;
t369 = t468 * t392 - t464 * t398;
t444 = qJDD(4) - t452;
t456 = qJD(4) - t488;
t365 = (t445 * t456 - t420) * pkin(9) + (t445 * t446 + t444) * pkin(4) + t369;
t370 = t464 * t392 + t468 * t398;
t429 = pkin(4) * t456 - pkin(9) * t446;
t443 = t445 ^ 2;
t367 = -pkin(4) * t443 + pkin(9) * t419 - t429 * t456 + t370;
t359 = t467 * t365 - t463 * t367;
t440 = qJDD(5) + t444;
t455 = qJD(5) + t456;
t354 = -0.2e1 * qJD(6) * t423 + (t422 * t455 - t382) * qJ(6) + (t422 * t423 + t440) * pkin(5) + t359;
t404 = -mrSges(7,2) * t455 + mrSges(7,3) * t422;
t486 = m(7) * t354 + t440 * mrSges(7,1) + t455 * t404;
t351 = -t382 * mrSges(7,3) - t423 * t399 + t486;
t360 = t463 * t365 + t467 * t367;
t381 = -qJD(5) * t423 + t419 * t467 - t420 * t463;
t406 = pkin(5) * t455 - qJ(6) * t423;
t421 = t422 ^ 2;
t357 = -pkin(5) * t421 + qJ(6) * t381 + 0.2e1 * qJD(6) * t422 - t406 * t455 + t360;
t492 = Ifges(6,4) + Ifges(7,4);
t497 = Ifges(6,5) + Ifges(7,5);
t498 = Ifges(6,1) + Ifges(7,1);
t493 = t492 * t422 + t498 * t423 + t497 * t455;
t496 = Ifges(6,6) + Ifges(7,6);
t500 = Ifges(6,2) + Ifges(7,2);
t494 = t500 * t422 + t423 * t492 + t496 * t455;
t495 = Ifges(6,3) + Ifges(7,3);
t501 = mrSges(6,1) * t359 + mrSges(7,1) * t354 - mrSges(6,2) * t360 - mrSges(7,2) * t357 + pkin(5) * t351 + t496 * t381 + t497 * t382 - t493 * t422 + t494 * t423 + t495 * t440;
t400 = -mrSges(6,1) * t422 + mrSges(6,2) * t423;
t405 = -mrSges(6,2) * t455 + mrSges(6,3) * t422;
t345 = m(6) * t359 + t440 * mrSges(6,1) + t455 * t405 + (-t399 - t400) * t423 + (-mrSges(6,3) - mrSges(7,3)) * t382 + t486;
t407 = mrSges(7,1) * t455 - mrSges(7,3) * t423;
t408 = mrSges(6,1) * t455 - mrSges(6,3) * t423;
t485 = m(7) * t357 + t381 * mrSges(7,3) + t422 * t399;
t348 = m(6) * t360 + t381 * mrSges(6,3) + t422 * t400 + (-t407 - t408) * t455 + (-mrSges(6,2) - mrSges(7,2)) * t440 + t485;
t343 = t467 * t345 + t463 * t348;
t414 = Ifges(5,4) * t446 + Ifges(5,2) * t445 + Ifges(5,6) * t456;
t415 = Ifges(5,1) * t446 + Ifges(5,4) * t445 + Ifges(5,5) * t456;
t499 = mrSges(5,1) * t369 - mrSges(5,2) * t370 + Ifges(5,5) * t420 + Ifges(5,6) * t419 + Ifges(5,3) * t444 + pkin(4) * t343 + t446 * t414 - t445 * t415 + t501;
t491 = -t496 * t422 - t497 * t423 - t495 * t455;
t448 = (-mrSges(4,1) * t469 + mrSges(4,2) * t465) * qJD(1);
t453 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t489;
t426 = -mrSges(5,1) * t445 + mrSges(5,2) * t446;
t427 = -mrSges(5,2) * t456 + mrSges(5,3) * t445;
t340 = m(5) * t369 + mrSges(5,1) * t444 - mrSges(5,3) * t420 - t426 * t446 + t427 * t456 + t343;
t428 = mrSges(5,1) * t456 - mrSges(5,3) * t446;
t479 = -t345 * t463 + t467 * t348;
t341 = m(5) * t370 - mrSges(5,2) * t444 + mrSges(5,3) * t419 + t426 * t445 - t428 * t456 + t479;
t480 = -t340 * t464 + t468 * t341;
t336 = m(4) * t403 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t452 - qJD(3) * t453 + t448 * t488 + t480;
t402 = -t465 * t412 + t460 * t469;
t454 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t488;
t397 = -qJDD(3) * pkin(3) - pkin(8) * t471 + t450 * t489 - t402;
t368 = -pkin(4) * t419 - pkin(9) * t443 + t446 * t429 + t397;
t362 = -pkin(5) * t381 - qJ(6) * t421 + t406 * t423 + qJDD(6) + t368;
t355 = m(7) * t362 - t381 * mrSges(7,1) + t382 * mrSges(7,2) - t422 * t404 + t423 * t407;
t477 = m(6) * t368 - t381 * mrSges(6,1) + t382 * mrSges(6,2) - t422 * t405 + t423 * t408 + t355;
t474 = -m(5) * t397 + t419 * mrSges(5,1) - t420 * mrSges(5,2) + t445 * t427 - t446 * t428 - t477;
t349 = m(4) * t402 + qJDD(3) * mrSges(4,1) - t451 * mrSges(4,3) + qJD(3) * t454 - t448 * t489 + t474;
t481 = t469 * t336 - t349 * t465;
t337 = t340 * t468 + t341 * t464;
t476 = -m(4) * t411 + t452 * mrSges(4,1) - mrSges(4,2) * t451 - t453 * t489 + t454 * t488 - t337;
t439 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t465 + Ifges(4,4) * t469) * qJD(1);
t438 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t465 + Ifges(4,2) * t469) * qJD(1);
t413 = Ifges(5,5) * t446 + Ifges(5,6) * t445 + Ifges(5,3) * t456;
t342 = mrSges(6,2) * t368 + mrSges(7,2) * t362 - mrSges(6,3) * t359 - mrSges(7,3) * t354 - qJ(6) * t351 + t492 * t381 + t498 * t382 - t491 * t422 + t497 * t440 - t494 * t455;
t338 = -mrSges(6,1) * t368 + mrSges(6,3) * t360 - mrSges(7,1) * t362 + mrSges(7,3) * t357 - pkin(5) * t355 + qJ(6) * t485 + (-qJ(6) * t407 + t493) * t455 + (-mrSges(7,2) * qJ(6) + t496) * t440 + t491 * t423 + t492 * t382 + t500 * t381;
t334 = mrSges(5,2) * t397 - mrSges(5,3) * t369 + Ifges(5,1) * t420 + Ifges(5,4) * t419 + Ifges(5,5) * t444 - pkin(9) * t343 - t338 * t463 + t342 * t467 + t413 * t445 - t414 * t456;
t333 = -mrSges(5,1) * t397 + mrSges(5,3) * t370 + Ifges(5,4) * t420 + Ifges(5,2) * t419 + Ifges(5,6) * t444 - pkin(4) * t477 + pkin(9) * t479 + t467 * t338 + t463 * t342 - t446 * t413 + t456 * t415;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t483 - mrSges(2,2) * t478 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t424 - mrSges(3,2) * t425 + t465 * (mrSges(4,2) * t411 - mrSges(4,3) * t402 + Ifges(4,1) * t451 + Ifges(4,4) * t452 + Ifges(4,5) * qJDD(3) - pkin(8) * t337 - qJD(3) * t438 - t464 * t333 + t468 * t334) + t469 * (-mrSges(4,1) * t411 + mrSges(4,3) * t403 + Ifges(4,4) * t451 + Ifges(4,2) * t452 + Ifges(4,6) * qJDD(3) - pkin(3) * t337 + qJD(3) * t439 - t499) + pkin(2) * t476 + pkin(7) * t481 + pkin(1) * (t461 * (m(3) * t425 - mrSges(3,1) * t472 - qJDD(1) * mrSges(3,2) + t481) + t462 * (m(3) * t424 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t472 + t476)); m(3) * t460 + t336 * t465 + t349 * t469; Ifges(4,5) * t451 + Ifges(4,6) * t452 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t402 - mrSges(4,2) * t403 + t464 * t334 + t468 * t333 + pkin(3) * t474 + pkin(8) * t480 + (t438 * t465 - t439 * t469) * qJD(1); t499; t501; t355;];
tauJ  = t1;
