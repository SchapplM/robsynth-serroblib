% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR13_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:32:38
% EndTime: 2019-12-31 20:32:43
% DurationCPUTime: 3.07s
% Computational Cost: add. (31430->274), mult. (68402->350), div. (0->0), fcn. (47582->10), ass. (0->108)
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
t353 = (-t395 - t425) * qJ(3) + (-t396 + t401) * pkin(2) + t386;
t420 = -g(1) * t413 - g(2) * t409;
t387 = -pkin(1) * t415 + qJDD(1) * pkin(6) + t420;
t370 = -g(3) * t408 + t412 * t387;
t393 = (-t412 * pkin(2) - t408 * qJ(3)) * qJD(1);
t414 = qJD(2) ^ 2;
t428 = qJD(1) * t412;
t356 = -pkin(2) * t414 + qJDD(2) * qJ(3) + t393 * t428 + t370;
t404 = sin(pkin(9));
t405 = cos(pkin(9));
t427 = t408 * qJD(1);
t391 = qJD(2) * t404 + t405 * t427;
t336 = -0.2e1 * qJD(3) * t391 + t405 * t353 - t356 * t404;
t375 = qJDD(2) * t404 + t395 * t405;
t390 = qJD(2) * t405 - t404 * t427;
t327 = (-t390 * t428 - t375) * pkin(7) + (t390 * t391 - t396) * pkin(3) + t336;
t337 = 0.2e1 * qJD(3) * t390 + t404 * t353 + t405 * t356;
t374 = qJDD(2) * t405 - t395 * t404;
t376 = -pkin(3) * t428 - pkin(7) * t391;
t389 = t390 ^ 2;
t329 = -pkin(3) * t389 + pkin(7) * t374 + t376 * t428 + t337;
t407 = sin(qJ(4));
t411 = cos(qJ(4));
t314 = t411 * t327 - t329 * t407;
t366 = t390 * t411 - t391 * t407;
t343 = qJD(4) * t366 + t374 * t407 + t375 * t411;
t367 = t390 * t407 + t391 * t411;
t392 = qJDD(4) - t396;
t400 = qJD(4) - t428;
t312 = (t366 * t400 - t343) * pkin(8) + (t366 * t367 + t392) * pkin(4) + t314;
t315 = t407 * t327 + t411 * t329;
t342 = -qJD(4) * t367 + t374 * t411 - t375 * t407;
t359 = pkin(4) * t400 - pkin(8) * t367;
t365 = t366 ^ 2;
t313 = -pkin(4) * t365 + pkin(8) * t342 - t359 * t400 + t315;
t406 = sin(qJ(5));
t410 = cos(qJ(5));
t310 = t312 * t410 - t313 * t406;
t348 = t366 * t410 - t367 * t406;
t323 = qJD(5) * t348 + t342 * t406 + t343 * t410;
t349 = t366 * t406 + t367 * t410;
t335 = -mrSges(6,1) * t348 + mrSges(6,2) * t349;
t399 = qJD(5) + t400;
t339 = -mrSges(6,2) * t399 + mrSges(6,3) * t348;
t388 = qJDD(5) + t392;
t306 = m(6) * t310 + mrSges(6,1) * t388 - mrSges(6,3) * t323 - t335 * t349 + t339 * t399;
t311 = t312 * t406 + t313 * t410;
t322 = -qJD(5) * t349 + t342 * t410 - t343 * t406;
t340 = mrSges(6,1) * t399 - mrSges(6,3) * t349;
t307 = m(6) * t311 - mrSges(6,2) * t388 + mrSges(6,3) * t322 + t335 * t348 - t340 * t399;
t300 = t410 * t306 + t406 * t307;
t350 = -mrSges(5,1) * t366 + mrSges(5,2) * t367;
t357 = -mrSges(5,2) * t400 + mrSges(5,3) * t366;
t298 = m(5) * t314 + mrSges(5,1) * t392 - mrSges(5,3) * t343 - t350 * t367 + t357 * t400 + t300;
t358 = mrSges(5,1) * t400 - mrSges(5,3) * t367;
t421 = -t306 * t406 + t410 * t307;
t299 = m(5) * t315 - mrSges(5,2) * t392 + mrSges(5,3) * t342 + t350 * t366 - t358 * t400 + t421;
t294 = t411 * t298 + t407 * t299;
t369 = -t412 * g(3) - t408 * t387;
t368 = -mrSges(4,1) * t390 + mrSges(4,2) * t391;
t372 = mrSges(4,2) * t428 + mrSges(4,3) * t390;
t292 = m(4) * t336 - mrSges(4,1) * t396 - mrSges(4,3) * t375 - t368 * t391 - t372 * t428 + t294;
t373 = -mrSges(4,1) * t428 - mrSges(4,3) * t391;
t422 = -t298 * t407 + t411 * t299;
t293 = m(4) * t337 + mrSges(4,2) * t396 + mrSges(4,3) * t374 + t368 * t390 + t373 * t428 + t422;
t423 = -t292 * t404 + t405 * t293;
t355 = -qJDD(2) * pkin(2) - qJ(3) * t414 + t393 * t427 + qJDD(3) - t369;
t338 = -pkin(3) * t374 - pkin(7) * t389 + t391 * t376 + t355;
t317 = -pkin(4) * t342 - pkin(8) * t365 + t359 * t367 + t338;
t419 = m(6) * t317 - t322 * mrSges(6,1) + t323 * mrSges(6,2) - t348 * t339 + t349 * t340;
t288 = t292 * t405 + t293 * t404;
t331 = Ifges(6,4) * t349 + Ifges(6,2) * t348 + Ifges(6,6) * t399;
t332 = Ifges(6,1) * t349 + Ifges(6,4) * t348 + Ifges(6,5) * t399;
t418 = -mrSges(6,1) * t310 + mrSges(6,2) * t311 - Ifges(6,5) * t323 - Ifges(6,6) * t322 - Ifges(6,3) * t388 - t349 * t331 + t348 * t332;
t417 = m(5) * t338 - t342 * mrSges(5,1) + t343 * mrSges(5,2) - t366 * t357 + t367 * t358 + t419;
t308 = m(4) * t355 - t374 * mrSges(4,1) + t375 * mrSges(4,2) - t390 * t372 + t391 * t373 + t417;
t345 = Ifges(5,4) * t367 + Ifges(5,2) * t366 + Ifges(5,6) * t400;
t346 = Ifges(5,1) * t367 + Ifges(5,4) * t366 + Ifges(5,5) * t400;
t416 = mrSges(5,1) * t314 - mrSges(5,2) * t315 + Ifges(5,5) * t343 + Ifges(5,6) * t342 + Ifges(5,3) * t392 + pkin(4) * t300 + t367 * t345 - t366 * t346 - t418;
t398 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t428;
t397 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t427;
t394 = (-t412 * mrSges(3,1) + t408 * mrSges(3,2)) * qJD(1);
t385 = Ifges(3,5) * qJD(2) + (t408 * Ifges(3,1) + t412 * Ifges(3,4)) * qJD(1);
t384 = Ifges(3,6) * qJD(2) + (t408 * Ifges(3,4) + Ifges(3,2) * t412) * qJD(1);
t362 = Ifges(4,1) * t391 + Ifges(4,4) * t390 - Ifges(4,5) * t428;
t361 = Ifges(4,4) * t391 + Ifges(4,2) * t390 - Ifges(4,6) * t428;
t360 = Ifges(4,5) * t391 + Ifges(4,6) * t390 - Ifges(4,3) * t428;
t344 = Ifges(5,5) * t367 + Ifges(5,6) * t366 + Ifges(5,3) * t400;
t330 = Ifges(6,5) * t349 + Ifges(6,6) * t348 + Ifges(6,3) * t399;
t302 = mrSges(6,2) * t317 - mrSges(6,3) * t310 + Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * t388 + t330 * t348 - t331 * t399;
t301 = -mrSges(6,1) * t317 + mrSges(6,3) * t311 + Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * t388 - t330 * t349 + t332 * t399;
t290 = mrSges(5,2) * t338 - mrSges(5,3) * t314 + Ifges(5,1) * t343 + Ifges(5,4) * t342 + Ifges(5,5) * t392 - pkin(8) * t300 - t301 * t406 + t302 * t410 + t344 * t366 - t345 * t400;
t289 = -mrSges(5,1) * t338 + mrSges(5,3) * t315 + Ifges(5,4) * t343 + Ifges(5,2) * t342 + Ifges(5,6) * t392 - pkin(4) * t419 + pkin(8) * t421 + t410 * t301 + t406 * t302 - t367 * t344 + t400 * t346;
t287 = mrSges(4,2) * t355 - mrSges(4,3) * t336 + Ifges(4,1) * t375 + Ifges(4,4) * t374 - Ifges(4,5) * t396 - pkin(7) * t294 - t289 * t407 + t290 * t411 + t360 * t390 + t361 * t428;
t286 = -mrSges(4,1) * t355 + mrSges(4,3) * t337 + Ifges(4,4) * t375 + Ifges(4,2) * t374 - Ifges(4,6) * t396 - pkin(3) * t417 + pkin(7) * t422 + t411 * t289 + t407 * t290 - t391 * t360 - t362 * t428;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t424 - mrSges(2,2) * t420 + t408 * (mrSges(3,2) * t386 - mrSges(3,3) * t369 + Ifges(3,1) * t395 + Ifges(3,4) * t396 + Ifges(3,5) * qJDD(2) - qJ(3) * t288 - qJD(2) * t384 - t404 * t286 + t405 * t287) + t412 * (-mrSges(3,1) * t386 - mrSges(4,1) * t336 + mrSges(4,2) * t337 + mrSges(3,3) * t370 + Ifges(3,4) * t395 - Ifges(4,5) * t375 + Ifges(3,6) * qJDD(2) - Ifges(4,6) * t374 - pkin(2) * t288 - pkin(3) * t294 + qJD(2) * t385 - t391 * t361 + t390 * t362 - t416 + (Ifges(3,2) + Ifges(4,3)) * t396) + pkin(1) * (-m(3) * t386 + mrSges(3,1) * t396 - mrSges(3,2) * t395 + (-t397 * t408 + t398 * t412) * qJD(1) - t288) + pkin(6) * (t412 * (m(3) * t370 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t396 - qJD(2) * t397 + t394 * t428 + t423) - t408 * (m(3) * t369 + qJDD(2) * mrSges(3,1) - t395 * mrSges(3,3) + qJD(2) * t398 - t394 * t427 - t308)); Ifges(3,5) * t395 + Ifges(3,6) * t396 + Ifges(3,3) * qJDD(2) + mrSges(3,1) * t369 - mrSges(3,2) * t370 + t404 * t287 + t405 * t286 - pkin(2) * t308 + qJ(3) * t423 + (t408 * t384 - t412 * t385) * qJD(1); t308; t416; -t418;];
tauJ = t1;
