% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:38:00
% EndTime: 2019-12-31 21:38:07
% DurationCPUTime: 5.08s
% Computational Cost: add. (54655->287), mult. (119825->382), div. (0->0), fcn. (93217->12), ass. (0->123)
t421 = sin(pkin(5));
t426 = sin(qJ(2));
t429 = cos(qJ(2));
t444 = qJD(1) * qJD(2);
t410 = (-qJDD(1) * t429 + t426 * t444) * t421;
t453 = cos(qJ(3));
t452 = t421 * pkin(7);
t423 = cos(pkin(5));
t451 = t423 * g(3);
t431 = qJD(1) ^ 2;
t427 = sin(qJ(1));
t430 = cos(qJ(1));
t441 = t427 * g(1) - t430 * g(2);
t405 = qJDD(1) * pkin(1) + t431 * t452 + t441;
t450 = t405 * t423;
t449 = t421 * t426;
t448 = t421 * t429;
t446 = qJD(1) * t421;
t408 = (-t429 * pkin(2) - t426 * pkin(8)) * t446;
t417 = t423 * qJD(1) + qJD(2);
t415 = t417 ^ 2;
t416 = t423 * qJDD(1) + qJDD(2);
t445 = qJD(1) * t429;
t438 = -t430 * g(1) - t427 * g(2);
t406 = -t431 * pkin(1) + qJDD(1) * t452 + t438;
t447 = t429 * t406 + t426 * t450;
t361 = -t415 * pkin(2) + t416 * pkin(8) + (-g(3) * t426 + t408 * t445) * t421 + t447;
t409 = (qJDD(1) * t426 + t429 * t444) * t421;
t362 = t410 * pkin(2) - t409 * pkin(8) - t451 + (-t405 + (pkin(2) * t426 - pkin(8) * t429) * t417 * qJD(1)) * t421;
t425 = sin(qJ(3));
t345 = t453 * t361 + t425 * t362;
t443 = t426 * t446;
t398 = -t453 * t417 + t425 * t443;
t399 = t425 * t417 + t453 * t443;
t382 = t398 * pkin(3) - t399 * qJ(4);
t402 = qJDD(3) + t410;
t442 = t421 * t445;
t414 = qJD(3) - t442;
t413 = t414 ^ 2;
t337 = -t413 * pkin(3) + t402 * qJ(4) - t398 * t382 + t345;
t380 = -g(3) * t448 - t426 * t406 + t429 * t450;
t360 = -t416 * pkin(2) - t415 * pkin(8) + t408 * t443 - t380;
t378 = t399 * qJD(3) + t425 * t409 - t453 * t416;
t379 = -t398 * qJD(3) + t453 * t409 + t425 * t416;
t340 = (t398 * t414 - t379) * qJ(4) + (t399 * t414 + t378) * pkin(3) + t360;
t420 = sin(pkin(10));
t422 = cos(pkin(10));
t388 = t422 * t399 + t420 * t414;
t332 = -0.2e1 * qJD(4) * t388 - t420 * t337 + t422 * t340;
t367 = t422 * t379 + t420 * t402;
t387 = -t420 * t399 + t422 * t414;
t330 = (t398 * t387 - t367) * pkin(9) + (t387 * t388 + t378) * pkin(4) + t332;
t333 = 0.2e1 * qJD(4) * t387 + t422 * t337 + t420 * t340;
t366 = -t420 * t379 + t422 * t402;
t371 = t398 * pkin(4) - t388 * pkin(9);
t386 = t387 ^ 2;
t331 = -t386 * pkin(4) + t366 * pkin(9) - t398 * t371 + t333;
t424 = sin(qJ(5));
t428 = cos(qJ(5));
t328 = t428 * t330 - t424 * t331;
t363 = t428 * t387 - t424 * t388;
t343 = t363 * qJD(5) + t424 * t366 + t428 * t367;
t364 = t424 * t387 + t428 * t388;
t350 = -t363 * mrSges(6,1) + t364 * mrSges(6,2);
t397 = qJD(5) + t398;
t351 = -t397 * mrSges(6,2) + t363 * mrSges(6,3);
t376 = qJDD(5) + t378;
t325 = m(6) * t328 + t376 * mrSges(6,1) - t343 * mrSges(6,3) - t364 * t350 + t397 * t351;
t329 = t424 * t330 + t428 * t331;
t342 = -t364 * qJD(5) + t428 * t366 - t424 * t367;
t352 = t397 * mrSges(6,1) - t364 * mrSges(6,3);
t326 = m(6) * t329 - t376 * mrSges(6,2) + t342 * mrSges(6,3) + t363 * t350 - t397 * t352;
t317 = t428 * t325 + t424 * t326;
t368 = -t387 * mrSges(5,1) + t388 * mrSges(5,2);
t437 = -t398 * mrSges(5,2) + t387 * mrSges(5,3);
t315 = m(5) * t332 + t378 * mrSges(5,1) - t367 * mrSges(5,3) - t388 * t368 + t398 * t437 + t317;
t370 = t398 * mrSges(5,1) - t388 * mrSges(5,3);
t439 = -t424 * t325 + t428 * t326;
t316 = m(5) * t333 - t378 * mrSges(5,2) + t366 * mrSges(5,3) + t387 * t368 - t398 * t370 + t439;
t313 = -t420 * t315 + t422 * t316;
t383 = t398 * mrSges(4,1) + t399 * mrSges(4,2);
t390 = t414 * mrSges(4,1) - t399 * mrSges(4,3);
t311 = m(4) * t345 - t402 * mrSges(4,2) - t378 * mrSges(4,3) - t398 * t383 - t414 * t390 + t313;
t344 = -t425 * t361 + t453 * t362;
t336 = -t402 * pkin(3) - t413 * qJ(4) + t399 * t382 + qJDD(4) - t344;
t334 = -t366 * pkin(4) - t386 * pkin(9) + t388 * t371 + t336;
t435 = m(6) * t334 - t342 * mrSges(6,1) + t343 * mrSges(6,2) - t363 * t351 + t364 * t352;
t327 = m(5) * t336 - t366 * mrSges(5,1) + t367 * mrSges(5,2) + t388 * t370 - t387 * t437 + t435;
t389 = -t414 * mrSges(4,2) - t398 * mrSges(4,3);
t321 = m(4) * t344 + t402 * mrSges(4,1) - t379 * mrSges(4,3) - t399 * t383 + t414 * t389 - t327;
t305 = t425 * t311 + t453 * t321;
t440 = t453 * t311 - t425 * t321;
t312 = t422 * t315 + t420 * t316;
t434 = -m(4) * t360 - t378 * mrSges(4,1) - t379 * mrSges(4,2) - t398 * t389 - t399 * t390 - t312;
t347 = Ifges(6,4) * t364 + Ifges(6,2) * t363 + Ifges(6,6) * t397;
t348 = Ifges(6,1) * t364 + Ifges(6,4) * t363 + Ifges(6,5) * t397;
t433 = mrSges(6,1) * t328 - mrSges(6,2) * t329 + Ifges(6,5) * t343 + Ifges(6,6) * t342 + Ifges(6,3) * t376 + t364 * t347 - t363 * t348;
t346 = Ifges(6,5) * t364 + Ifges(6,6) * t363 + Ifges(6,3) * t397;
t318 = -mrSges(6,1) * t334 + mrSges(6,3) * t329 + Ifges(6,4) * t343 + Ifges(6,2) * t342 + Ifges(6,6) * t376 - t364 * t346 + t397 * t348;
t319 = mrSges(6,2) * t334 - mrSges(6,3) * t328 + Ifges(6,1) * t343 + Ifges(6,4) * t342 + Ifges(6,5) * t376 + t363 * t346 - t397 * t347;
t353 = Ifges(5,5) * t388 + Ifges(5,6) * t387 + Ifges(5,3) * t398;
t355 = Ifges(5,1) * t388 + Ifges(5,4) * t387 + Ifges(5,5) * t398;
t306 = -mrSges(5,1) * t336 + mrSges(5,3) * t333 + Ifges(5,4) * t367 + Ifges(5,2) * t366 + Ifges(5,6) * t378 - pkin(4) * t435 + pkin(9) * t439 + t428 * t318 + t424 * t319 - t388 * t353 + t398 * t355;
t354 = Ifges(5,4) * t388 + Ifges(5,2) * t387 + Ifges(5,6) * t398;
t307 = mrSges(5,2) * t336 - mrSges(5,3) * t332 + Ifges(5,1) * t367 + Ifges(5,4) * t366 + Ifges(5,5) * t378 - pkin(9) * t317 - t424 * t318 + t428 * t319 + t387 * t353 - t398 * t354;
t373 = Ifges(4,4) * t399 - Ifges(4,2) * t398 + Ifges(4,6) * t414;
t374 = Ifges(4,1) * t399 - Ifges(4,4) * t398 + Ifges(4,5) * t414;
t432 = mrSges(4,1) * t344 - mrSges(4,2) * t345 + Ifges(4,5) * t379 - Ifges(4,6) * t378 + Ifges(4,3) * t402 - pkin(3) * t327 + qJ(4) * t313 + t422 * t306 + t420 * t307 + t399 * t373 + t398 * t374;
t407 = (-t429 * mrSges(3,1) + t426 * mrSges(3,2)) * t446;
t404 = -t417 * mrSges(3,2) + mrSges(3,3) * t442;
t403 = t417 * mrSges(3,1) - mrSges(3,3) * t443;
t394 = -t421 * t405 - t451;
t393 = Ifges(3,5) * t417 + (t426 * Ifges(3,1) + t429 * Ifges(3,4)) * t446;
t392 = Ifges(3,6) * t417 + (t426 * Ifges(3,4) + t429 * Ifges(3,2)) * t446;
t391 = Ifges(3,3) * t417 + (t426 * Ifges(3,5) + t429 * Ifges(3,6)) * t446;
t381 = -g(3) * t449 + t447;
t372 = Ifges(4,5) * t399 - Ifges(4,6) * t398 + Ifges(4,3) * t414;
t308 = m(3) * t380 + t416 * mrSges(3,1) - t409 * mrSges(3,3) + t417 * t404 - t407 * t443 + t434;
t304 = m(3) * t381 - t416 * mrSges(3,2) - t410 * mrSges(3,3) - t417 * t403 + t407 * t442 + t440;
t303 = (-Ifges(4,2) - Ifges(5,3)) * t378 - t433 + t414 * t374 + Ifges(4,6) * t402 - t399 * t372 + t387 * t355 - t388 * t354 + Ifges(4,4) * t379 - Ifges(5,6) * t366 - Ifges(5,5) * t367 - mrSges(4,1) * t360 + mrSges(4,3) * t345 - mrSges(5,1) * t332 + mrSges(5,2) * t333 - pkin(4) * t317 - pkin(3) * t312;
t302 = mrSges(4,2) * t360 - mrSges(4,3) * t344 + Ifges(4,1) * t379 - Ifges(4,4) * t378 + Ifges(4,5) * t402 - qJ(4) * t312 - t420 * t306 + t422 * t307 - t398 * t372 - t414 * t373;
t301 = Ifges(3,5) * t409 - Ifges(3,6) * t410 + Ifges(3,3) * t416 + mrSges(3,1) * t380 - mrSges(3,2) * t381 + t425 * t302 + t453 * t303 + pkin(2) * t434 + pkin(8) * t440 + (t426 * t392 - t429 * t393) * t446;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t441 - mrSges(2,2) * t438 + (mrSges(3,2) * t394 - mrSges(3,3) * t380 + Ifges(3,1) * t409 - Ifges(3,4) * t410 + Ifges(3,5) * t416 - pkin(8) * t305 + t453 * t302 - t425 * t303 + t391 * t442 - t417 * t392) * t449 + (-mrSges(3,1) * t394 + mrSges(3,3) * t381 + Ifges(3,4) * t409 - Ifges(3,2) * t410 + Ifges(3,6) * t416 - pkin(2) * t305 - t391 * t443 + t417 * t393 - t432) * t448 + t423 * t301 + pkin(1) * ((t426 * t304 + t429 * t308) * t423 + (-m(3) * t394 - t410 * mrSges(3,1) - t409 * mrSges(3,2) + (-t403 * t426 + t404 * t429) * t446 - t305) * t421) + (t429 * t304 - t426 * t308) * t452; t301; t432; t327; t433;];
tauJ = t1;
