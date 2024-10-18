% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR11
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
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:34
% EndTime: 2024-09-27 21:45:35
% DurationCPUTime: 0.34s
% Computational Cost: add. (14161->213), mult. (35247->286), div. (0->0), fcn. (26732->12), ass. (0->99)
t410 = qJD(2) ^ 2;
t398 = sin(pkin(10));
t400 = cos(pkin(10));
t384 = g(1) * t398 - g(2) * t400;
t385 = -g(1) * t400 - g(2) * t398;
t405 = sin(qJ(2));
t409 = cos(qJ(2));
t416 = t409 * t384 - t385 * t405;
t399 = sin(pkin(5));
t428 = pkin(7) * t399;
t359 = qJDD(2) * pkin(2) + t410 * t428 + t416;
t397 = -g(3) + qJDD(1);
t401 = cos(pkin(5));
t429 = t359 * t401 + t397 * t399;
t426 = t399 ^ 2 * t410;
t424 = t405 * t384 + t409 * t385;
t360 = -pkin(2) * t410 + qJDD(2) * t428 + t424;
t404 = sin(qJ(3));
t408 = cos(qJ(3));
t338 = -t360 * t404 + t429 * t408;
t422 = qJD(2) * qJD(3);
t375 = (qJDD(2) * t404 + t408 * t422) * t399;
t390 = t401 * qJDD(2) + qJDD(3);
t391 = t401 * qJD(2) + qJD(3);
t423 = qJD(2) * t399;
t419 = t408 * t423;
t329 = (t391 * t419 - t375) * pkin(8) + (t404 * t408 * t426 + t390) * pkin(3) + t338;
t339 = t408 * t360 + t429 * t404;
t420 = t404 * t423;
t374 = pkin(3) * t391 - pkin(8) * t420;
t376 = (qJDD(2) * t408 - t404 * t422) * t399;
t421 = t408 ^ 2 * t426;
t330 = -pkin(3) * t421 + pkin(8) * t376 - t374 * t391 + t339;
t403 = sin(qJ(4));
t407 = cos(qJ(4));
t316 = t407 * t329 - t330 * t403;
t369 = (-t403 * t404 + t407 * t408) * t423;
t345 = qJD(4) * t369 + t375 * t407 + t376 * t403;
t370 = (t408 * t403 + t404 * t407) * t423;
t388 = qJDD(4) + t390;
t389 = qJD(4) + t391;
t313 = (t369 * t389 - t345) * pkin(9) + (t369 * t370 + t388) * pkin(4) + t316;
t317 = t403 * t329 + t407 * t330;
t344 = -qJD(4) * t370 - t375 * t403 + t376 * t407;
t363 = pkin(4) * t389 - pkin(9) * t370;
t368 = t369 ^ 2;
t314 = -pkin(4) * t368 + pkin(9) * t344 - t363 * t389 + t317;
t402 = sin(qJ(5));
t406 = cos(qJ(5));
t311 = t313 * t406 - t314 * t402;
t352 = t369 * t406 - t370 * t402;
t324 = qJD(5) * t352 + t344 * t402 + t345 * t406;
t353 = t369 * t402 + t370 * t406;
t335 = -mrSges(6,1) * t352 + mrSges(6,2) * t353;
t383 = qJD(5) + t389;
t346 = -mrSges(6,2) * t383 + mrSges(6,3) * t352;
t382 = qJDD(5) + t388;
t308 = m(6) * t311 + mrSges(6,1) * t382 - mrSges(6,3) * t324 - t335 * t353 + t346 * t383;
t312 = t313 * t402 + t314 * t406;
t323 = -qJD(5) * t353 + t344 * t406 - t345 * t402;
t347 = mrSges(6,1) * t383 - mrSges(6,3) * t353;
t309 = m(6) * t312 - mrSges(6,2) * t382 + mrSges(6,3) * t323 + t335 * t352 - t347 * t383;
t301 = t406 * t308 + t402 * t309;
t355 = -mrSges(5,1) * t369 + mrSges(5,2) * t370;
t361 = -mrSges(5,2) * t389 + mrSges(5,3) * t369;
t298 = m(5) * t316 + mrSges(5,1) * t388 - mrSges(5,3) * t345 - t355 * t370 + t361 * t389 + t301;
t362 = mrSges(5,1) * t389 - mrSges(5,3) * t370;
t417 = -t308 * t402 + t406 * t309;
t299 = m(5) * t317 - mrSges(5,2) * t388 + mrSges(5,3) * t344 + t355 * t369 - t362 * t389 + t417;
t294 = t407 * t298 + t403 * t299;
t418 = -t298 * t403 + t407 * t299;
t354 = -t359 * t399 + t401 * t397;
t372 = -mrSges(4,2) * t391 + mrSges(4,3) * t419;
t373 = (-t408 * mrSges(4,1) + t404 * mrSges(4,2)) * t423;
t292 = m(4) * t338 + mrSges(4,1) * t390 - mrSges(4,3) * t375 + t372 * t391 - t373 * t420 + t294;
t371 = mrSges(4,1) * t391 - mrSges(4,3) * t420;
t293 = m(4) * t339 - mrSges(4,2) * t390 + mrSges(4,3) * t376 - t371 * t391 + t373 * t419 + t418;
t415 = t292 * t408 + t293 * t404;
t337 = -pkin(3) * t376 - pkin(8) * t421 + t374 * t420 + t354;
t319 = -pkin(4) * t344 - pkin(9) * t368 + t363 * t370 + t337;
t414 = m(6) * t319 - mrSges(6,1) * t323 + t324 * mrSges(6,2) - t346 * t352 + t353 * t347;
t332 = Ifges(6,4) * t353 + Ifges(6,2) * t352 + Ifges(6,6) * t383;
t333 = Ifges(6,1) * t353 + Ifges(6,4) * t352 + Ifges(6,5) * t383;
t413 = mrSges(6,1) * t311 - mrSges(6,2) * t312 + Ifges(6,5) * t324 + Ifges(6,6) * t323 + Ifges(6,3) * t382 + t353 * t332 - t352 * t333;
t412 = m(5) * t337 - mrSges(5,1) * t344 + t345 * mrSges(5,2) - t369 * t361 + t370 * t362 + t414;
t349 = Ifges(5,4) * t370 + Ifges(5,2) * t369 + Ifges(5,6) * t389;
t350 = Ifges(5,1) * t370 + Ifges(5,4) * t369 + Ifges(5,5) * t389;
t411 = mrSges(5,1) * t316 - mrSges(5,2) * t317 + Ifges(5,5) * t345 + Ifges(5,6) * t344 + Ifges(5,3) * t388 + pkin(4) * t301 + t370 * t349 - t369 * t350 + t413;
t366 = Ifges(4,5) * t391 + (t404 * Ifges(4,1) + t408 * Ifges(4,4)) * t423;
t365 = Ifges(4,6) * t391 + (t404 * Ifges(4,4) + t408 * Ifges(4,2)) * t423;
t348 = Ifges(5,5) * t370 + Ifges(5,6) * t369 + Ifges(5,3) * t389;
t331 = Ifges(6,5) * t353 + Ifges(6,6) * t352 + Ifges(6,3) * t383;
t304 = t412 + m(4) * t354 - mrSges(4,1) * t376 + mrSges(4,2) * t375 + (t371 * t404 - t372 * t408) * t423;
t303 = mrSges(6,2) * t319 - mrSges(6,3) * t311 + Ifges(6,1) * t324 + Ifges(6,4) * t323 + Ifges(6,5) * t382 + t331 * t352 - t332 * t383;
t302 = -mrSges(6,1) * t319 + mrSges(6,3) * t312 + Ifges(6,4) * t324 + Ifges(6,2) * t323 + Ifges(6,6) * t382 - t331 * t353 + t333 * t383;
t291 = mrSges(5,2) * t337 - mrSges(5,3) * t316 + Ifges(5,1) * t345 + Ifges(5,4) * t344 + Ifges(5,5) * t388 - pkin(9) * t301 - t302 * t402 + t303 * t406 + t348 * t369 - t349 * t389;
t290 = -mrSges(5,1) * t337 + mrSges(5,3) * t317 + Ifges(5,4) * t345 + Ifges(5,2) * t344 + Ifges(5,6) * t388 - pkin(4) * t414 + pkin(9) * t417 + t406 * t302 + t402 * t303 - t370 * t348 + t389 * t350;
t289 = (t404 * t365 - t408 * t366) * t423 + t411 + Ifges(4,3) * t390 + Ifges(4,5) * t375 + Ifges(4,6) * t376 + mrSges(4,1) * t338 - mrSges(4,2) * t339 + pkin(3) * t294;
t1 = [t304 * t401 + t415 * t399 + (m(2) + m(3)) * t397; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t416 - mrSges(3,2) * t424 + (pkin(2) * t415 + t289) * t401 + (t404 * (mrSges(4,2) * t354 - mrSges(4,3) * t338 + Ifges(4,1) * t375 + Ifges(4,4) * t376 + Ifges(4,5) * t390 - pkin(8) * t294 - t290 * t403 + t291 * t407 - t365 * t391) + t408 * (-mrSges(4,1) * t354 + mrSges(4,3) * t339 + Ifges(4,4) * t375 + Ifges(4,2) * t376 + Ifges(4,6) * t390 - pkin(3) * t412 + pkin(8) * t418 + t407 * t290 + t403 * t291 + t391 * t366) - pkin(2) * t304 + pkin(7) * (-t292 * t404 + t293 * t408)) * t399; t289; t411; t413;];
tauJ = t1;
