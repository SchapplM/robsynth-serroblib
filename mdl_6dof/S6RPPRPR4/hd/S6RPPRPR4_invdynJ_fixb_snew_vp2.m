% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-05 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:15:10
% EndTime: 2019-05-05 14:15:12
% DurationCPUTime: 1.44s
% Computational Cost: add. (12153->238), mult. (24070->302), div. (0->0), fcn. (13177->10), ass. (0->101)
t385 = sin(pkin(10));
t387 = cos(pkin(10));
t390 = sin(qJ(4));
t393 = cos(qJ(4));
t364 = (t390 * t385 - t393 * t387) * qJD(1);
t416 = 2 * qJD(5);
t415 = -pkin(1) - pkin(2);
t396 = qJD(1) ^ 2;
t391 = sin(qJ(1));
t394 = cos(qJ(1));
t406 = -t394 * g(1) - t391 * g(2);
t401 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t406;
t356 = t396 * t415 + t401;
t410 = t391 * g(1) - t394 * g(2);
t400 = -t396 * qJ(2) + qJDD(2) - t410;
t359 = qJDD(1) * t415 + t400;
t386 = sin(pkin(9));
t388 = cos(pkin(9));
t338 = t388 * t356 + t386 * t359;
t335 = -t396 * pkin(3) - qJDD(1) * pkin(7) + t338;
t383 = g(3) + qJDD(3);
t331 = -t390 * t335 + t393 * t383;
t412 = qJD(1) * qJD(4);
t411 = t393 * t412;
t372 = -t390 * qJDD(1) - t411;
t321 = (-t372 - t411) * qJ(5) + (t390 * t393 * t396 + qJDD(4)) * pkin(4) + t331;
t332 = t393 * t335 + t390 * t383;
t373 = -t393 * qJDD(1) + t390 * t412;
t413 = t390 * qJD(1);
t374 = qJD(4) * pkin(4) + qJ(5) * t413;
t382 = t393 ^ 2;
t322 = -t382 * t396 * pkin(4) + t373 * qJ(5) - qJD(4) * t374 + t332;
t317 = t385 * t321 + t387 * t322 + t364 * t416;
t365 = (t393 * t385 + t390 * t387) * qJD(1);
t345 = -t364 * mrSges(6,1) - t365 * mrSges(6,2);
t349 = -t385 * t372 + t387 * t373;
t358 = qJD(4) * mrSges(6,1) + t365 * mrSges(6,3);
t346 = -t364 * pkin(5) + t365 * pkin(8);
t395 = qJD(4) ^ 2;
t315 = -t395 * pkin(5) + qJDD(4) * pkin(8) + t364 * t346 + t317;
t337 = -t386 * t356 + t388 * t359;
t405 = qJDD(1) * pkin(3) - t337;
t323 = -t374 * t413 - t373 * pkin(4) + qJDD(5) + (-qJ(5) * t382 - pkin(7)) * t396 + t405;
t350 = t387 * t372 + t385 * t373;
t318 = (-qJD(4) * t364 - t350) * pkin(8) + (-qJD(4) * t365 - t349) * pkin(5) + t323;
t389 = sin(qJ(6));
t392 = cos(qJ(6));
t312 = -t389 * t315 + t392 * t318;
t353 = t392 * qJD(4) + t389 * t365;
t330 = t353 * qJD(6) + t389 * qJDD(4) + t392 * t350;
t354 = qJD(4) * t389 - t365 * t392;
t336 = -mrSges(7,1) * t353 + mrSges(7,2) * t354;
t362 = qJD(6) - t364;
t339 = -mrSges(7,2) * t362 + mrSges(7,3) * t353;
t348 = qJDD(6) - t349;
t310 = m(7) * t312 + mrSges(7,1) * t348 - mrSges(7,3) * t330 - t336 * t354 + t339 * t362;
t313 = t315 * t392 + t318 * t389;
t329 = -qJD(6) * t354 + qJDD(4) * t392 - t350 * t389;
t340 = mrSges(7,1) * t362 - mrSges(7,3) * t354;
t311 = m(7) * t313 - mrSges(7,2) * t348 + mrSges(7,3) * t329 + t336 * t353 - t340 * t362;
t407 = -t310 * t389 + t311 * t392;
t300 = m(6) * t317 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t349 - qJD(4) * t358 + t345 * t364 + t407;
t403 = -t387 * t321 + t385 * t322;
t316 = t365 * t416 - t403;
t357 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t364;
t314 = -qJDD(4) * pkin(5) - t395 * pkin(8) + (-(2 * qJD(5)) - t346) * t365 + t403;
t399 = -m(7) * t314 + mrSges(7,1) * t329 - mrSges(7,2) * t330 + t339 * t353 - t340 * t354;
t306 = m(6) * t316 + qJDD(4) * mrSges(6,1) - mrSges(6,3) * t350 + qJD(4) * t357 + t345 * t365 + t399;
t296 = t300 * t385 + t306 * t387;
t302 = t310 * t392 + t311 * t389;
t414 = qJD(1) * t393;
t371 = (mrSges(5,1) * t393 - mrSges(5,2) * t390) * qJD(1);
t376 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t414;
t294 = m(5) * t331 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t372 + qJD(4) * t376 + t371 * t413 + t296;
t375 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t413;
t408 = t300 * t387 - t306 * t385;
t295 = m(5) * t332 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t373 - qJD(4) * t375 - t371 * t414 + t408;
t409 = -t294 * t390 + t295 * t393;
t290 = m(4) * t338 - mrSges(4,1) * t396 + qJDD(1) * mrSges(4,2) + t409;
t301 = m(6) * t323 - mrSges(6,1) * t349 + mrSges(6,2) * t350 - t357 * t364 - t358 * t365 + t302;
t334 = -t396 * pkin(7) + t405;
t397 = -m(5) * t334 + mrSges(5,1) * t373 - mrSges(5,2) * t372 + t375 * t413 - t376 * t414 - t301;
t297 = m(4) * t337 - qJDD(1) * mrSges(4,1) - mrSges(4,2) * t396 + t397;
t404 = t290 * t386 + t297 * t388;
t325 = Ifges(7,4) * t354 + Ifges(7,2) * t353 + Ifges(7,6) * t362;
t326 = Ifges(7,1) * t354 + Ifges(7,4) * t353 + Ifges(7,5) * t362;
t398 = mrSges(7,1) * t312 - mrSges(7,2) * t313 + Ifges(7,5) * t330 + Ifges(7,6) * t329 + Ifges(7,3) * t348 + t325 * t354 - t326 * t353;
t368 = Ifges(5,5) * qJD(4) + (-Ifges(5,1) * t390 - Ifges(5,4) * t393) * qJD(1);
t367 = Ifges(5,6) * qJD(4) + (-Ifges(5,4) * t390 - Ifges(5,2) * t393) * qJD(1);
t363 = -qJDD(1) * pkin(1) + t400;
t360 = -pkin(1) * t396 + t401;
t343 = -Ifges(6,1) * t365 + Ifges(6,4) * t364 + Ifges(6,5) * qJD(4);
t342 = -Ifges(6,4) * t365 + Ifges(6,2) * t364 + Ifges(6,6) * qJD(4);
t341 = -Ifges(6,5) * t365 + Ifges(6,6) * t364 + Ifges(6,3) * qJD(4);
t324 = Ifges(7,5) * t354 + Ifges(7,6) * t353 + Ifges(7,3) * t362;
t304 = mrSges(7,2) * t314 - mrSges(7,3) * t312 + Ifges(7,1) * t330 + Ifges(7,4) * t329 + Ifges(7,5) * t348 + t324 * t353 - t325 * t362;
t303 = -mrSges(7,1) * t314 + mrSges(7,3) * t313 + Ifges(7,4) * t330 + Ifges(7,2) * t329 + Ifges(7,6) * t348 - t324 * t354 + t326 * t362;
t292 = -mrSges(6,1) * t323 + mrSges(6,3) * t317 + Ifges(6,4) * t350 + Ifges(6,2) * t349 + Ifges(6,6) * qJDD(4) - pkin(5) * t302 + qJD(4) * t343 + t341 * t365 - t398;
t291 = mrSges(6,2) * t323 - mrSges(6,3) * t316 + Ifges(6,1) * t350 + Ifges(6,4) * t349 + Ifges(6,5) * qJDD(4) - pkin(8) * t302 - qJD(4) * t342 - t303 * t389 + t304 * t392 + t341 * t364;
t289 = m(3) * t363 - qJDD(1) * mrSges(3,1) - mrSges(3,3) * t396 + t404;
t1 = [-pkin(1) * t289 + qJ(2) * (m(3) * t360 - mrSges(3,1) * t396 + t290 * t388 - t297 * t386) - mrSges(2,2) * t406 + mrSges(2,1) * t410 - pkin(2) * t404 - mrSges(3,1) * t363 + mrSges(3,3) * t360 - t390 * (mrSges(5,2) * t334 - mrSges(5,3) * t331 + Ifges(5,1) * t372 + Ifges(5,4) * t373 + Ifges(5,5) * qJDD(4) - qJ(5) * t296 - qJD(4) * t367 + t291 * t387 - t292 * t385) - t393 * (-mrSges(5,1) * t334 + mrSges(5,3) * t332 + Ifges(5,4) * t372 + Ifges(5,2) * t373 + Ifges(5,6) * qJDD(4) - pkin(4) * t301 + qJ(5) * t408 + qJD(4) * t368 + t385 * t291 + t387 * t292) - pkin(3) * t397 - pkin(7) * t409 - mrSges(4,1) * t337 + mrSges(4,2) * t338 + (mrSges(3,3) * qJ(2) + Ifges(3,2) + Ifges(2,3) + Ifges(4,3)) * qJDD(1); t289; m(4) * t383 + t294 * t393 + t295 * t390; Ifges(5,5) * t372 + Ifges(5,6) * t373 + mrSges(5,1) * t331 - mrSges(5,2) * t332 + Ifges(6,5) * t350 + Ifges(6,6) * t349 - t365 * t342 - t364 * t343 + mrSges(6,1) * t316 - mrSges(6,2) * t317 + t389 * t304 + t392 * t303 + pkin(5) * t399 + pkin(8) * t407 + pkin(4) * t296 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (-t367 * t390 + t368 * t393) * qJD(1); t301; t398;];
tauJ  = t1;
