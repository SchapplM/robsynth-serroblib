% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR9
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR9_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:47
% EndTime: 2019-12-31 21:22:53
% DurationCPUTime: 3.53s
% Computational Cost: add. (34814->274), mult. (72553->350), div. (0->0), fcn. (50325->10), ass. (0->108)
t415 = qJD(1) ^ 2;
t409 = sin(qJ(1));
t413 = cos(qJ(1));
t424 = g(1) * t409 - t413 * g(2);
t386 = -qJDD(1) * pkin(1) - pkin(6) * t415 - t424;
t408 = sin(qJ(2));
t412 = cos(qJ(2));
t426 = qJD(1) * qJD(2);
t425 = t412 * t426;
t395 = qJDD(1) * t408 + t425;
t401 = t408 * t426;
t396 = qJDD(1) * t412 - t401;
t353 = (-t395 - t425) * pkin(7) + (-t396 + t401) * pkin(2) + t386;
t420 = -g(1) * t413 - g(2) * t409;
t387 = -pkin(1) * t415 + qJDD(1) * pkin(6) + t420;
t378 = -g(3) * t408 + t412 * t387;
t394 = (-t412 * pkin(2) - t408 * pkin(7)) * qJD(1);
t414 = qJD(2) ^ 2;
t428 = qJD(1) * t412;
t356 = -pkin(2) * t414 + qJDD(2) * pkin(7) + t394 * t428 + t378;
t407 = sin(qJ(3));
t411 = cos(qJ(3));
t337 = t411 * t353 - t356 * t407;
t427 = t408 * qJD(1);
t391 = qJD(2) * t411 - t407 * t427;
t369 = qJD(3) * t391 + qJDD(2) * t407 + t395 * t411;
t390 = qJDD(3) - t396;
t392 = qJD(2) * t407 + t411 * t427;
t400 = qJD(3) - t428;
t327 = (t391 * t400 - t369) * qJ(4) + (t391 * t392 + t390) * pkin(3) + t337;
t338 = t407 * t353 + t411 * t356;
t368 = -qJD(3) * t392 + qJDD(2) * t411 - t395 * t407;
t375 = pkin(3) * t400 - qJ(4) * t392;
t389 = t391 ^ 2;
t329 = -pkin(3) * t389 + qJ(4) * t368 - t375 * t400 + t338;
t404 = sin(pkin(9));
t405 = cos(pkin(9));
t372 = t391 * t404 + t392 * t405;
t314 = -0.2e1 * qJD(4) * t372 + t405 * t327 - t329 * t404;
t346 = t368 * t404 + t369 * t405;
t371 = t391 * t405 - t392 * t404;
t312 = (t371 * t400 - t346) * pkin(8) + (t371 * t372 + t390) * pkin(4) + t314;
t315 = 0.2e1 * qJD(4) * t371 + t404 * t327 + t405 * t329;
t345 = t368 * t405 - t369 * t404;
t359 = pkin(4) * t400 - pkin(8) * t372;
t370 = t371 ^ 2;
t313 = -pkin(4) * t370 + pkin(8) * t345 - t359 * t400 + t315;
t406 = sin(qJ(5));
t410 = cos(qJ(5));
t310 = t312 * t410 - t313 * t406;
t348 = t371 * t410 - t372 * t406;
t323 = qJD(5) * t348 + t345 * t406 + t346 * t410;
t349 = t371 * t406 + t372 * t410;
t335 = -mrSges(6,1) * t348 + mrSges(6,2) * t349;
t399 = qJD(5) + t400;
t339 = -mrSges(6,2) * t399 + mrSges(6,3) * t348;
t388 = qJDD(5) + t390;
t306 = m(6) * t310 + mrSges(6,1) * t388 - mrSges(6,3) * t323 - t335 * t349 + t339 * t399;
t311 = t312 * t406 + t313 * t410;
t322 = -qJD(5) * t349 + t345 * t410 - t346 * t406;
t340 = mrSges(6,1) * t399 - mrSges(6,3) * t349;
t307 = m(6) * t311 - mrSges(6,2) * t388 + mrSges(6,3) * t322 + t335 * t348 - t340 * t399;
t300 = t410 * t306 + t406 * t307;
t350 = -mrSges(5,1) * t371 + mrSges(5,2) * t372;
t357 = -mrSges(5,2) * t400 + mrSges(5,3) * t371;
t298 = m(5) * t314 + mrSges(5,1) * t390 - mrSges(5,3) * t346 - t350 * t372 + t357 * t400 + t300;
t358 = mrSges(5,1) * t400 - mrSges(5,3) * t372;
t421 = -t306 * t406 + t410 * t307;
t299 = m(5) * t315 - mrSges(5,2) * t390 + mrSges(5,3) * t345 + t350 * t371 - t358 * t400 + t421;
t294 = t405 * t298 + t404 * t299;
t343 = Ifges(5,4) * t372 + Ifges(5,2) * t371 + Ifges(5,6) * t400;
t344 = Ifges(5,1) * t372 + Ifges(5,4) * t371 + Ifges(5,5) * t400;
t361 = Ifges(4,4) * t392 + Ifges(4,2) * t391 + Ifges(4,6) * t400;
t362 = Ifges(4,1) * t392 + Ifges(4,4) * t391 + Ifges(4,5) * t400;
t331 = Ifges(6,4) * t349 + Ifges(6,2) * t348 + Ifges(6,6) * t399;
t332 = Ifges(6,1) * t349 + Ifges(6,4) * t348 + Ifges(6,5) * t399;
t418 = -mrSges(6,1) * t310 + mrSges(6,2) * t311 - Ifges(6,5) * t323 - Ifges(6,6) * t322 - Ifges(6,3) * t388 - t349 * t331 + t348 * t332;
t430 = mrSges(4,1) * t337 + mrSges(5,1) * t314 - mrSges(4,2) * t338 - mrSges(5,2) * t315 + Ifges(4,5) * t369 + Ifges(5,5) * t346 + Ifges(4,6) * t368 + Ifges(5,6) * t345 + pkin(3) * t294 + pkin(4) * t300 + t372 * t343 - t371 * t344 + t392 * t361 - t391 * t362 + (Ifges(4,3) + Ifges(5,3)) * t390 - t418;
t377 = -t412 * g(3) - t408 * t387;
t373 = -mrSges(4,1) * t391 + mrSges(4,2) * t392;
t374 = -mrSges(4,2) * t400 + mrSges(4,3) * t391;
t292 = m(4) * t337 + mrSges(4,1) * t390 - mrSges(4,3) * t369 - t373 * t392 + t374 * t400 + t294;
t376 = mrSges(4,1) * t400 - mrSges(4,3) * t392;
t422 = -t298 * t404 + t405 * t299;
t293 = m(4) * t338 - mrSges(4,2) * t390 + mrSges(4,3) * t368 + t373 * t391 - t376 * t400 + t422;
t423 = -t292 * t407 + t411 * t293;
t355 = -qJDD(2) * pkin(2) - pkin(7) * t414 + t394 * t427 - t377;
t336 = -pkin(3) * t368 - qJ(4) * t389 + t392 * t375 + qJDD(4) + t355;
t317 = -pkin(4) * t345 - pkin(8) * t370 + t359 * t372 + t336;
t419 = m(6) * t317 - t322 * mrSges(6,1) + t323 * mrSges(6,2) - t348 * t339 + t349 * t340;
t288 = t292 * t411 + t293 * t407;
t308 = m(5) * t336 - t345 * mrSges(5,1) + t346 * mrSges(5,2) - t371 * t357 + t372 * t358 + t419;
t417 = -m(4) * t355 + t368 * mrSges(4,1) - t369 * mrSges(4,2) + t391 * t374 - t392 * t376 - t308;
t398 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t428;
t397 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t427;
t393 = (-t412 * mrSges(3,1) + t408 * mrSges(3,2)) * qJD(1);
t385 = Ifges(3,5) * qJD(2) + (t408 * Ifges(3,1) + t412 * Ifges(3,4)) * qJD(1);
t384 = Ifges(3,6) * qJD(2) + (t408 * Ifges(3,4) + t412 * Ifges(3,2)) * qJD(1);
t360 = Ifges(4,5) * t392 + Ifges(4,6) * t391 + Ifges(4,3) * t400;
t342 = Ifges(5,5) * t372 + Ifges(5,6) * t371 + Ifges(5,3) * t400;
t330 = Ifges(6,5) * t349 + Ifges(6,6) * t348 + Ifges(6,3) * t399;
t302 = mrSges(6,2) * t317 - mrSges(6,3) * t310 + Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * t388 + t330 * t348 - t331 * t399;
t301 = -mrSges(6,1) * t317 + mrSges(6,3) * t311 + Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * t388 - t330 * t349 + t332 * t399;
t290 = mrSges(5,2) * t336 - mrSges(5,3) * t314 + Ifges(5,1) * t346 + Ifges(5,4) * t345 + Ifges(5,5) * t390 - pkin(8) * t300 - t301 * t406 + t302 * t410 + t342 * t371 - t343 * t400;
t289 = -mrSges(5,1) * t336 + mrSges(5,3) * t315 + Ifges(5,4) * t346 + Ifges(5,2) * t345 + Ifges(5,6) * t390 - pkin(4) * t419 + pkin(8) * t421 + t410 * t301 + t406 * t302 - t372 * t342 + t400 * t344;
t287 = mrSges(4,2) * t355 - mrSges(4,3) * t337 + Ifges(4,1) * t369 + Ifges(4,4) * t368 + Ifges(4,5) * t390 - qJ(4) * t294 - t289 * t404 + t290 * t405 + t360 * t391 - t361 * t400;
t286 = -mrSges(4,1) * t355 + mrSges(4,3) * t338 + Ifges(4,4) * t369 + Ifges(4,2) * t368 + Ifges(4,6) * t390 - pkin(3) * t308 + qJ(4) * t422 + t405 * t289 + t404 * t290 - t392 * t360 + t400 * t362;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t424 - mrSges(2,2) * t420 + t408 * (mrSges(3,2) * t386 - mrSges(3,3) * t377 + Ifges(3,1) * t395 + Ifges(3,4) * t396 + Ifges(3,5) * qJDD(2) - pkin(7) * t288 - qJD(2) * t384 - t286 * t407 + t287 * t411) + t412 * (-mrSges(3,1) * t386 + mrSges(3,3) * t378 + Ifges(3,4) * t395 + Ifges(3,2) * t396 + Ifges(3,6) * qJDD(2) - pkin(2) * t288 + qJD(2) * t385 - t430) + pkin(1) * (-m(3) * t386 + mrSges(3,1) * t396 - mrSges(3,2) * t395 + (-t397 * t408 + t398 * t412) * qJD(1) - t288) + pkin(6) * (t412 * (m(3) * t378 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t396 - qJD(2) * t397 + t393 * t428 + t423) - t408 * (m(3) * t377 + qJDD(2) * mrSges(3,1) - t395 * mrSges(3,3) + qJD(2) * t398 - t393 * t427 + t417)); Ifges(3,5) * t395 + Ifges(3,6) * t396 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t377 - mrSges(3,2) * t378 + t407 * t287 + t411 * t286 + pkin(2) * t417 + pkin(7) * t423 + (t408 * t384 - t412 * t385) * qJD(1); t430; t308; -t418;];
tauJ = t1;
