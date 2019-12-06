% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:19
% EndTime: 2019-12-05 18:38:23
% DurationCPUTime: 3.73s
% Computational Cost: add. (36555->276), mult. (84103->356), div. (0->0), fcn. (60046->10), ass. (0->109)
t435 = qJD(1) ^ 2;
t451 = pkin(2) * t435;
t430 = sin(qJ(1));
t434 = cos(qJ(1));
t442 = -t434 * g(1) - t430 * g(2);
t408 = -t435 * pkin(1) + qJDD(1) * pkin(6) + t442;
t429 = sin(qJ(2));
t450 = t429 * t408;
t433 = cos(qJ(2));
t447 = qJD(1) * qJD(2);
t411 = t429 * qJDD(1) + t433 * t447;
t375 = qJDD(2) * pkin(2) - t411 * pkin(7) - t450 + (pkin(7) * t447 + t429 * t451 - g(3)) * t433;
t396 = -t429 * g(3) + t433 * t408;
t412 = t433 * qJDD(1) - t429 * t447;
t449 = qJD(1) * t429;
t415 = qJD(2) * pkin(2) - pkin(7) * t449;
t424 = t433 ^ 2;
t376 = t412 * pkin(7) - qJD(2) * t415 - t424 * t451 + t396;
t428 = sin(qJ(3));
t432 = cos(qJ(3));
t355 = t432 * t375 - t428 * t376;
t405 = (-t429 * t428 + t433 * t432) * qJD(1);
t380 = t405 * qJD(3) + t432 * t411 + t428 * t412;
t406 = (t433 * t428 + t429 * t432) * qJD(1);
t422 = qJDD(2) + qJDD(3);
t423 = qJD(2) + qJD(3);
t343 = (t405 * t423 - t380) * qJ(4) + (t405 * t406 + t422) * pkin(3) + t355;
t356 = t428 * t375 + t432 * t376;
t379 = -t406 * qJD(3) - t428 * t411 + t432 * t412;
t398 = t423 * pkin(3) - t406 * qJ(4);
t401 = t405 ^ 2;
t345 = -t401 * pkin(3) + t379 * qJ(4) - t423 * t398 + t356;
t425 = sin(pkin(9));
t426 = cos(pkin(9));
t393 = t425 * t405 + t426 * t406;
t331 = -0.2e1 * qJD(4) * t393 + t426 * t343 - t425 * t345;
t361 = t425 * t379 + t426 * t380;
t392 = t426 * t405 - t425 * t406;
t328 = (t392 * t423 - t361) * pkin(8) + (t392 * t393 + t422) * pkin(4) + t331;
t332 = 0.2e1 * qJD(4) * t392 + t425 * t343 + t426 * t345;
t360 = t426 * t379 - t425 * t380;
t384 = t423 * pkin(4) - t393 * pkin(8);
t390 = t392 ^ 2;
t329 = -t390 * pkin(4) + t360 * pkin(8) - t423 * t384 + t332;
t427 = sin(qJ(5));
t431 = cos(qJ(5));
t326 = t431 * t328 - t427 * t329;
t368 = t431 * t392 - t427 * t393;
t339 = t368 * qJD(5) + t427 * t360 + t431 * t361;
t369 = t427 * t392 + t431 * t393;
t350 = -t368 * mrSges(6,1) + t369 * mrSges(6,2);
t420 = qJD(5) + t423;
t362 = -t420 * mrSges(6,2) + t368 * mrSges(6,3);
t419 = qJDD(5) + t422;
t322 = m(6) * t326 + t419 * mrSges(6,1) - t339 * mrSges(6,3) - t369 * t350 + t420 * t362;
t327 = t427 * t328 + t431 * t329;
t338 = -t369 * qJD(5) + t431 * t360 - t427 * t361;
t363 = t420 * mrSges(6,1) - t369 * mrSges(6,3);
t323 = m(6) * t327 - t419 * mrSges(6,2) + t338 * mrSges(6,3) + t368 * t350 - t420 * t363;
t316 = t431 * t322 + t427 * t323;
t370 = -t392 * mrSges(5,1) + t393 * mrSges(5,2);
t382 = -t423 * mrSges(5,2) + t392 * mrSges(5,3);
t313 = m(5) * t331 + t422 * mrSges(5,1) - t361 * mrSges(5,3) - t393 * t370 + t423 * t382 + t316;
t383 = t423 * mrSges(5,1) - t393 * mrSges(5,3);
t443 = -t427 * t322 + t431 * t323;
t314 = m(5) * t332 - t422 * mrSges(5,2) + t360 * mrSges(5,3) + t392 * t370 - t423 * t383 + t443;
t309 = t426 * t313 + t425 * t314;
t394 = -t405 * mrSges(4,1) + t406 * mrSges(4,2);
t397 = -t423 * mrSges(4,2) + t405 * mrSges(4,3);
t306 = m(4) * t355 + t422 * mrSges(4,1) - t380 * mrSges(4,3) - t406 * t394 + t423 * t397 + t309;
t399 = t423 * mrSges(4,1) - t406 * mrSges(4,3);
t444 = -t425 * t313 + t426 * t314;
t307 = m(4) * t356 - t422 * mrSges(4,2) + t379 * mrSges(4,3) + t405 * t394 - t423 * t399 + t444;
t300 = t432 * t306 + t428 * t307;
t448 = qJD(1) * t433;
t446 = t430 * g(1) - t434 * g(2);
t445 = -t428 * t306 + t432 * t307;
t441 = -qJDD(1) * pkin(1) - t446;
t381 = -t412 * pkin(2) + t415 * t449 + (-pkin(7) * t424 - pkin(6)) * t435 + t441;
t352 = -t379 * pkin(3) - t401 * qJ(4) + t406 * t398 + qJDD(4) + t381;
t334 = -t360 * pkin(4) - t390 * pkin(8) + t393 * t384 + t352;
t440 = -m(6) * t334 + t338 * mrSges(6,1) - t339 * mrSges(6,2) + t368 * t362 - t369 * t363;
t347 = Ifges(6,4) * t369 + Ifges(6,2) * t368 + Ifges(6,6) * t420;
t348 = Ifges(6,1) * t369 + Ifges(6,4) * t368 + Ifges(6,5) * t420;
t439 = mrSges(6,1) * t326 - mrSges(6,2) * t327 + Ifges(6,5) * t339 + Ifges(6,6) * t338 + Ifges(6,3) * t419 + t369 * t347 - t368 * t348;
t438 = -m(5) * t352 + t360 * mrSges(5,1) - t361 * mrSges(5,2) + t392 * t382 - t393 * t383 + t440;
t437 = m(4) * t381 - t379 * mrSges(4,1) + t380 * mrSges(4,2) - t405 * t397 + t406 * t399 - t438;
t365 = Ifges(5,4) * t393 + Ifges(5,2) * t392 + Ifges(5,6) * t423;
t366 = Ifges(5,1) * t393 + Ifges(5,4) * t392 + Ifges(5,5) * t423;
t388 = Ifges(4,4) * t406 + Ifges(4,2) * t405 + Ifges(4,6) * t423;
t389 = Ifges(4,1) * t406 + Ifges(4,4) * t405 + Ifges(4,5) * t423;
t436 = mrSges(4,1) * t355 + mrSges(5,1) * t331 - mrSges(4,2) * t356 - mrSges(5,2) * t332 + Ifges(4,5) * t380 + Ifges(5,5) * t361 + Ifges(4,6) * t379 + Ifges(5,6) * t360 + pkin(3) * t309 + pkin(4) * t316 + t393 * t365 - t392 * t366 + t406 * t388 - t405 * t389 + t439 + (Ifges(5,3) + Ifges(4,3)) * t422;
t414 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t448;
t413 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t449;
t410 = (-t433 * mrSges(3,1) + t429 * mrSges(3,2)) * qJD(1);
t407 = -t435 * pkin(6) + t441;
t404 = Ifges(3,5) * qJD(2) + (t429 * Ifges(3,1) + t433 * Ifges(3,4)) * qJD(1);
t403 = Ifges(3,6) * qJD(2) + (t429 * Ifges(3,4) + t433 * Ifges(3,2)) * qJD(1);
t395 = -t433 * g(3) - t450;
t387 = Ifges(4,5) * t406 + Ifges(4,6) * t405 + Ifges(4,3) * t423;
t364 = Ifges(5,5) * t393 + Ifges(5,6) * t392 + Ifges(5,3) * t423;
t346 = Ifges(6,5) * t369 + Ifges(6,6) * t368 + Ifges(6,3) * t420;
t318 = mrSges(6,2) * t334 - mrSges(6,3) * t326 + Ifges(6,1) * t339 + Ifges(6,4) * t338 + Ifges(6,5) * t419 + t368 * t346 - t420 * t347;
t317 = -mrSges(6,1) * t334 + mrSges(6,3) * t327 + Ifges(6,4) * t339 + Ifges(6,2) * t338 + Ifges(6,6) * t419 - t369 * t346 + t420 * t348;
t302 = mrSges(5,2) * t352 - mrSges(5,3) * t331 + Ifges(5,1) * t361 + Ifges(5,4) * t360 + Ifges(5,5) * t422 - pkin(8) * t316 - t427 * t317 + t431 * t318 + t392 * t364 - t423 * t365;
t301 = -mrSges(5,1) * t352 + mrSges(5,3) * t332 + Ifges(5,4) * t361 + Ifges(5,2) * t360 + Ifges(5,6) * t422 + pkin(4) * t440 + pkin(8) * t443 + t431 * t317 + t427 * t318 - t393 * t364 + t423 * t366;
t299 = mrSges(4,2) * t381 - mrSges(4,3) * t355 + Ifges(4,1) * t380 + Ifges(4,4) * t379 + Ifges(4,5) * t422 - qJ(4) * t309 - t425 * t301 + t426 * t302 + t405 * t387 - t423 * t388;
t298 = -mrSges(4,1) * t381 + mrSges(4,3) * t356 + Ifges(4,4) * t380 + Ifges(4,2) * t379 + Ifges(4,6) * t422 + pkin(3) * t438 + qJ(4) * t444 + t426 * t301 + t425 * t302 - t406 * t387 + t423 * t389;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t446 - mrSges(2,2) * t442 + t429 * (mrSges(3,2) * t407 - mrSges(3,3) * t395 + Ifges(3,1) * t411 + Ifges(3,4) * t412 + Ifges(3,5) * qJDD(2) - pkin(7) * t300 - qJD(2) * t403 - t428 * t298 + t432 * t299) + t433 * (-mrSges(3,1) * t407 + mrSges(3,3) * t396 + Ifges(3,4) * t411 + Ifges(3,2) * t412 + Ifges(3,6) * qJDD(2) - pkin(2) * t437 + pkin(7) * t445 + qJD(2) * t404 + t432 * t298 + t428 * t299) + pkin(1) * (-t437 - m(3) * t407 - t411 * mrSges(3,2) + t412 * mrSges(3,1) + (-t429 * t413 + t433 * t414) * qJD(1)) + pkin(6) * (t433 * (m(3) * t396 - qJDD(2) * mrSges(3,2) + t412 * mrSges(3,3) - qJD(2) * t413 + t410 * t448 + t445) - t429 * (m(3) * t395 + qJDD(2) * mrSges(3,1) - t411 * mrSges(3,3) + qJD(2) * t414 - t410 * t449 + t300)); t436 + pkin(2) * t300 + Ifges(3,3) * qJDD(2) + (t429 * t403 - t433 * t404) * qJD(1) + mrSges(3,1) * t395 - mrSges(3,2) * t396 + Ifges(3,5) * t411 + Ifges(3,6) * t412; t436; -t438; t439;];
tauJ = t1;
