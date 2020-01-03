% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR12_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR12_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:29:47
% EndTime: 2019-12-31 18:29:49
% DurationCPUTime: 2.21s
% Computational Cost: add. (19195->248), mult. (47168->323), div. (0->0), fcn. (34198->10), ass. (0->107)
t415 = qJD(1) ^ 2;
t438 = cos(qJ(3));
t437 = pkin(2) * t415;
t436 = pkin(6) * qJDD(1);
t411 = sin(qJ(1));
t413 = cos(qJ(1));
t423 = -t413 * g(1) - t411 * g(2);
t395 = -t415 * pkin(1) + qJDD(1) * qJ(2) + t423;
t406 = sin(pkin(8));
t408 = cos(pkin(8));
t430 = qJD(1) * qJD(2);
t427 = -t408 * g(3) - 0.2e1 * t406 * t430;
t364 = (t408 * t437 - t395 - t436) * t406 + t427;
t382 = -t406 * g(3) + (t395 + 0.2e1 * t430) * t408;
t404 = t408 ^ 2;
t371 = -t404 * t437 + t408 * t436 + t382;
t410 = sin(qJ(3));
t350 = t410 * t364 + t438 * t371;
t429 = t408 * t438;
t433 = qJD(1) * t406;
t393 = -qJD(1) * t429 + t410 * t433;
t419 = t438 * t406 + t408 * t410;
t394 = t419 * qJD(1);
t376 = t393 * mrSges(4,1) + t394 * mrSges(4,2);
t431 = t394 * qJD(3);
t379 = t431 + (t406 * t410 - t429) * qJDD(1);
t389 = qJD(3) * mrSges(4,1) - t394 * mrSges(4,3);
t375 = t393 * pkin(3) - t394 * qJ(4);
t414 = qJD(3) ^ 2;
t340 = -t414 * pkin(3) + qJDD(3) * qJ(4) - t393 * t375 + t350;
t428 = t411 * g(1) - t413 * g(2);
t422 = qJDD(2) - t428;
t434 = -t406 ^ 2 - t404;
t378 = (-pkin(2) * t408 - pkin(1)) * qJDD(1) + (t434 * pkin(6) - qJ(2)) * t415 + t422;
t432 = t393 * qJD(3);
t380 = t419 * qJDD(1) - t432;
t343 = (-t380 + t432) * qJ(4) + (t379 + t431) * pkin(3) + t378;
t405 = sin(pkin(9));
t407 = cos(pkin(9));
t387 = t405 * qJD(3) + t407 * t394;
t332 = -0.2e1 * qJD(4) * t387 - t405 * t340 + t407 * t343;
t370 = t405 * qJDD(3) + t407 * t380;
t386 = t407 * qJD(3) - t405 * t394;
t330 = (t386 * t393 - t370) * pkin(7) + (t386 * t387 + t379) * pkin(4) + t332;
t333 = 0.2e1 * qJD(4) * t386 + t407 * t340 + t405 * t343;
t368 = t393 * pkin(4) - t387 * pkin(7);
t369 = t407 * qJDD(3) - t405 * t380;
t385 = t386 ^ 2;
t331 = -t385 * pkin(4) + t369 * pkin(7) - t393 * t368 + t333;
t409 = sin(qJ(5));
t412 = cos(qJ(5));
t328 = t412 * t330 - t409 * t331;
t357 = t412 * t386 - t409 * t387;
t339 = t357 * qJD(5) + t409 * t369 + t412 * t370;
t358 = t409 * t386 + t412 * t387;
t348 = -t357 * mrSges(6,1) + t358 * mrSges(6,2);
t391 = qJD(5) + t393;
t351 = -t391 * mrSges(6,2) + t357 * mrSges(6,3);
t377 = qJDD(5) + t379;
t325 = m(6) * t328 + t377 * mrSges(6,1) - t339 * mrSges(6,3) - t358 * t348 + t391 * t351;
t329 = t409 * t330 + t412 * t331;
t338 = -t358 * qJD(5) + t412 * t369 - t409 * t370;
t352 = t391 * mrSges(6,1) - t358 * mrSges(6,3);
t326 = m(6) * t329 - t377 * mrSges(6,2) + t338 * mrSges(6,3) + t357 * t348 - t391 * t352;
t317 = t412 * t325 + t409 * t326;
t359 = -t386 * mrSges(5,1) + t387 * mrSges(5,2);
t366 = -t393 * mrSges(5,2) + t386 * mrSges(5,3);
t315 = m(5) * t332 + t379 * mrSges(5,1) - t370 * mrSges(5,3) - t387 * t359 + t393 * t366 + t317;
t367 = t393 * mrSges(5,1) - t387 * mrSges(5,3);
t424 = -t409 * t325 + t412 * t326;
t316 = m(5) * t333 - t379 * mrSges(5,2) + t369 * mrSges(5,3) + t386 * t359 - t393 * t367 + t424;
t425 = -t405 * t315 + t407 * t316;
t310 = m(4) * t350 - qJDD(3) * mrSges(4,2) - t379 * mrSges(4,3) - qJD(3) * t389 - t393 * t376 + t425;
t349 = t438 * t364 - t410 * t371;
t337 = -qJDD(3) * pkin(3) - t414 * qJ(4) + t394 * t375 + qJDD(4) - t349;
t334 = -t369 * pkin(4) - t385 * pkin(7) + t387 * t368 + t337;
t418 = m(6) * t334 - t338 * mrSges(6,1) + t339 * mrSges(6,2) - t357 * t351 + t358 * t352;
t327 = m(5) * t337 - t369 * mrSges(5,1) + t370 * mrSges(5,2) - t386 * t366 + t387 * t367 + t418;
t388 = -qJD(3) * mrSges(4,2) - t393 * mrSges(4,3);
t321 = m(4) * t349 + qJDD(3) * mrSges(4,1) - t380 * mrSges(4,3) + qJD(3) * t388 - t394 * t376 - t327;
t435 = t410 * t310 + t438 * t321;
t311 = t407 * t315 + t405 * t316;
t426 = t438 * t310 - t410 * t321;
t421 = -t408 * mrSges(3,1) + t406 * mrSges(3,2);
t420 = mrSges(3,3) * qJDD(1) + t415 * t421;
t417 = m(4) * t378 + t379 * mrSges(4,1) + t380 * mrSges(4,2) + t393 * t388 + t394 * t389 + t311;
t345 = Ifges(6,4) * t358 + Ifges(6,2) * t357 + Ifges(6,6) * t391;
t346 = Ifges(6,1) * t358 + Ifges(6,4) * t357 + Ifges(6,5) * t391;
t416 = mrSges(6,1) * t328 - mrSges(6,2) * t329 + Ifges(6,5) * t339 + Ifges(6,6) * t338 + Ifges(6,3) * t377 + t358 * t345 - t357 * t346;
t397 = (Ifges(3,5) * t406 + Ifges(3,6) * t408) * qJD(1);
t392 = -qJDD(1) * pkin(1) - t415 * qJ(2) + t422;
t381 = -t406 * t395 + t427;
t374 = Ifges(4,1) * t394 - Ifges(4,4) * t393 + Ifges(4,5) * qJD(3);
t373 = Ifges(4,4) * t394 - Ifges(4,2) * t393 + Ifges(4,6) * qJD(3);
t372 = Ifges(4,5) * t394 - Ifges(4,6) * t393 + Ifges(4,3) * qJD(3);
t355 = Ifges(5,1) * t387 + Ifges(5,4) * t386 + Ifges(5,5) * t393;
t354 = Ifges(5,4) * t387 + Ifges(5,2) * t386 + Ifges(5,6) * t393;
t353 = Ifges(5,5) * t387 + Ifges(5,6) * t386 + Ifges(5,3) * t393;
t344 = Ifges(6,5) * t358 + Ifges(6,6) * t357 + Ifges(6,3) * t391;
t319 = mrSges(6,2) * t334 - mrSges(6,3) * t328 + Ifges(6,1) * t339 + Ifges(6,4) * t338 + Ifges(6,5) * t377 + t357 * t344 - t391 * t345;
t318 = -mrSges(6,1) * t334 + mrSges(6,3) * t329 + Ifges(6,4) * t339 + Ifges(6,2) * t338 + Ifges(6,6) * t377 - t358 * t344 + t391 * t346;
t307 = t434 * t415 * mrSges(3,3) + m(3) * t392 + t421 * qJDD(1) + t417;
t306 = mrSges(5,2) * t337 - mrSges(5,3) * t332 + Ifges(5,1) * t370 + Ifges(5,4) * t369 + Ifges(5,5) * t379 - pkin(7) * t317 - t409 * t318 + t412 * t319 + t386 * t353 - t393 * t354;
t305 = -mrSges(5,1) * t337 + mrSges(5,3) * t333 + Ifges(5,4) * t370 + Ifges(5,2) * t369 + Ifges(5,6) * t379 - pkin(4) * t418 + pkin(7) * t424 + t412 * t318 + t409 * t319 - t387 * t353 + t393 * t355;
t304 = -t416 - t394 * t372 + t386 * t355 - t387 * t354 + Ifges(4,4) * t380 - mrSges(4,1) * t378 - Ifges(5,5) * t370 + qJD(3) * t374 - Ifges(5,6) * t369 + mrSges(4,3) * t350 - mrSges(5,1) * t332 + mrSges(5,2) * t333 - pkin(4) * t317 - pkin(3) * t311 + Ifges(4,6) * qJDD(3) + (-Ifges(4,2) - Ifges(5,3)) * t379;
t303 = mrSges(4,2) * t378 - mrSges(4,3) * t349 + Ifges(4,1) * t380 - Ifges(4,4) * t379 + Ifges(4,5) * qJDD(3) - qJ(4) * t311 - qJD(3) * t373 - t405 * t305 + t407 * t306 - t393 * t372;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t428 - mrSges(2,2) * t423 + t406 * (t408 * qJD(1) * t397 + mrSges(3,2) * t392 - mrSges(3,3) * t381 + t438 * t303 - t410 * t304 - pkin(6) * t435 + (Ifges(3,1) * t406 + Ifges(3,4) * t408) * qJDD(1)) + t408 * (-t397 * t433 - mrSges(3,1) * t392 + mrSges(3,3) * t382 + t410 * t303 + t438 * t304 - pkin(2) * t417 + pkin(6) * t426 + (Ifges(3,4) * t406 + Ifges(3,2) * t408) * qJDD(1)) - pkin(1) * t307 + qJ(2) * ((m(3) * t382 + t420 * t408 + t426) * t408 + (-m(3) * t381 + t420 * t406 - t435) * t406); t307; mrSges(4,1) * t349 - mrSges(4,2) * t350 + Ifges(4,5) * t380 - Ifges(4,6) * t379 + Ifges(4,3) * qJDD(3) - pkin(3) * t327 + qJ(4) * t425 + t407 * t305 + t405 * t306 + t394 * t373 + t393 * t374; t327; t416;];
tauJ = t1;
