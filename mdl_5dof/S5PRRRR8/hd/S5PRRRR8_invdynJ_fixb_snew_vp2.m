% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRRR8
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRRR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:15:31
% EndTime: 2019-12-05 17:15:33
% DurationCPUTime: 1.24s
% Computational Cost: add. (8852->217), mult. (17346->285), div. (0->0), fcn. (12112->12), ass. (0->95)
t373 = sin(pkin(10));
t375 = cos(pkin(10));
t361 = t373 * g(1) - t375 * g(2);
t372 = -g(3) + qJDD(1);
t374 = sin(pkin(5));
t376 = cos(pkin(5));
t403 = t361 * t376 + t372 * t374;
t378 = sin(qJ(4));
t379 = sin(qJ(3));
t382 = cos(qJ(4));
t383 = cos(qJ(3));
t351 = (t379 * t378 - t383 * t382) * qJD(2);
t362 = -t375 * g(1) - t373 * g(2);
t380 = sin(qJ(2));
t384 = cos(qJ(2));
t334 = -t380 * t362 + t384 * t403;
t335 = t384 * t362 + t403 * t380;
t385 = qJD(2) ^ 2;
t330 = -t385 * pkin(2) + qJDD(2) * pkin(7) + t335;
t345 = -t374 * t361 + t376 * t372;
t316 = -t379 * t330 + t383 * t345;
t398 = qJD(2) * qJD(3);
t397 = t383 * t398;
t359 = t379 * qJDD(2) + t397;
t306 = (-t359 + t397) * pkin(8) + (t379 * t383 * t385 + qJDD(3)) * pkin(3) + t316;
t317 = t383 * t330 + t379 * t345;
t360 = t383 * qJDD(2) - t379 * t398;
t399 = t379 * qJD(2);
t366 = qJD(3) * pkin(3) - pkin(8) * t399;
t371 = t383 ^ 2;
t307 = -t371 * t385 * pkin(3) + t360 * pkin(8) - qJD(3) * t366 + t317;
t303 = t378 * t306 + t382 * t307;
t352 = (t383 * t378 + t379 * t382) * qJD(2);
t325 = -t352 * qJD(4) - t378 * t359 + t382 * t360;
t337 = t351 * mrSges(5,1) + t352 * mrSges(5,2);
t370 = qJD(3) + qJD(4);
t344 = t370 * mrSges(5,1) - t352 * mrSges(5,3);
t369 = qJDD(3) + qJDD(4);
t338 = t351 * pkin(4) - t352 * pkin(9);
t368 = t370 ^ 2;
t299 = -t368 * pkin(4) + t369 * pkin(9) - t351 * t338 + t303;
t390 = -qJDD(2) * pkin(2) - t334;
t311 = -t360 * pkin(3) + t366 * t399 + (-pkin(8) * t371 - pkin(7)) * t385 + t390;
t326 = -t351 * qJD(4) + t382 * t359 + t378 * t360;
t300 = (t351 * t370 - t326) * pkin(9) + (t352 * t370 - t325) * pkin(4) + t311;
t377 = sin(qJ(5));
t381 = cos(qJ(5));
t296 = -t377 * t299 + t381 * t300;
t339 = -t377 * t352 + t381 * t370;
t310 = t339 * qJD(5) + t381 * t326 + t377 * t369;
t340 = t381 * t352 + t377 * t370;
t318 = -t339 * mrSges(6,1) + t340 * mrSges(6,2);
t323 = qJDD(5) - t325;
t346 = qJD(5) + t351;
t327 = -t346 * mrSges(6,2) + t339 * mrSges(6,3);
t293 = m(6) * t296 + t323 * mrSges(6,1) - t310 * mrSges(6,3) - t340 * t318 + t346 * t327;
t297 = t381 * t299 + t377 * t300;
t309 = -t340 * qJD(5) - t377 * t326 + t381 * t369;
t328 = t346 * mrSges(6,1) - t340 * mrSges(6,3);
t294 = m(6) * t297 - t323 * mrSges(6,2) + t309 * mrSges(6,3) + t339 * t318 - t346 * t328;
t394 = -t377 * t293 + t381 * t294;
t281 = m(5) * t303 - t369 * mrSges(5,2) + t325 * mrSges(5,3) - t351 * t337 - t370 * t344 + t394;
t302 = t382 * t306 - t378 * t307;
t343 = -t370 * mrSges(5,2) - t351 * mrSges(5,3);
t298 = -t369 * pkin(4) - t368 * pkin(9) + t352 * t338 - t302;
t391 = -m(6) * t298 + t309 * mrSges(6,1) - t310 * mrSges(6,2) + t339 * t327 - t340 * t328;
t289 = m(5) * t302 + t369 * mrSges(5,1) - t326 * mrSges(5,3) - t352 * t337 + t370 * t343 + t391;
t278 = t378 * t281 + t382 * t289;
t283 = t381 * t293 + t377 * t294;
t400 = qJD(2) * t383;
t358 = (-t383 * mrSges(4,1) + t379 * mrSges(4,2)) * qJD(2);
t364 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t400;
t276 = m(4) * t316 + qJDD(3) * mrSges(4,1) - t359 * mrSges(4,3) + qJD(3) * t364 - t358 * t399 + t278;
t363 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t399;
t395 = t382 * t281 - t378 * t289;
t277 = m(4) * t317 - qJDD(3) * mrSges(4,2) + t360 * mrSges(4,3) - qJD(3) * t363 + t358 * t400 + t395;
t396 = -t379 * t276 + t383 * t277;
t389 = m(5) * t311 - t325 * mrSges(5,1) + t326 * mrSges(5,2) + t351 * t343 + t352 * t344 + t283;
t312 = Ifges(6,5) * t340 + Ifges(6,6) * t339 + Ifges(6,3) * t346;
t314 = Ifges(6,1) * t340 + Ifges(6,4) * t339 + Ifges(6,5) * t346;
t286 = -mrSges(6,1) * t298 + mrSges(6,3) * t297 + Ifges(6,4) * t310 + Ifges(6,2) * t309 + Ifges(6,6) * t323 - t340 * t312 + t346 * t314;
t313 = Ifges(6,4) * t340 + Ifges(6,2) * t339 + Ifges(6,6) * t346;
t287 = mrSges(6,2) * t298 - mrSges(6,3) * t296 + Ifges(6,1) * t310 + Ifges(6,4) * t309 + Ifges(6,5) * t323 + t339 * t312 - t346 * t313;
t332 = Ifges(5,4) * t352 - Ifges(5,2) * t351 + Ifges(5,6) * t370;
t333 = Ifges(5,1) * t352 - Ifges(5,4) * t351 + Ifges(5,5) * t370;
t388 = mrSges(5,1) * t302 - mrSges(5,2) * t303 + Ifges(5,5) * t326 + Ifges(5,6) * t325 + Ifges(5,3) * t369 + pkin(4) * t391 + pkin(9) * t394 + t381 * t286 + t377 * t287 + t352 * t332 + t351 * t333;
t387 = mrSges(6,1) * t296 - mrSges(6,2) * t297 + Ifges(6,5) * t310 + Ifges(6,6) * t309 + Ifges(6,3) * t323 + t340 * t313 - t339 * t314;
t329 = -t385 * pkin(7) + t390;
t386 = -m(4) * t329 + t360 * mrSges(4,1) - t359 * mrSges(4,2) - t363 * t399 + t364 * t400 - t389;
t350 = Ifges(4,5) * qJD(3) + (t379 * Ifges(4,1) + t383 * Ifges(4,4)) * qJD(2);
t349 = Ifges(4,6) * qJD(3) + (t379 * Ifges(4,4) + t383 * Ifges(4,2)) * qJD(2);
t331 = Ifges(5,5) * t352 - Ifges(5,6) * t351 + Ifges(5,3) * t370;
t274 = -mrSges(5,1) * t311 + mrSges(5,3) * t303 + Ifges(5,4) * t326 + Ifges(5,2) * t325 + Ifges(5,6) * t369 - pkin(4) * t283 - t352 * t331 + t370 * t333 - t387;
t273 = mrSges(5,2) * t311 - mrSges(5,3) * t302 + Ifges(5,1) * t326 + Ifges(5,4) * t325 + Ifges(5,5) * t369 - pkin(9) * t283 - t377 * t286 + t381 * t287 - t351 * t331 - t370 * t332;
t1 = [m(2) * t372 + t376 * (m(3) * t345 + t383 * t276 + t379 * t277) + (t380 * (m(3) * t335 - t385 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t396) + t384 * (m(3) * t334 + qJDD(2) * mrSges(3,1) - t385 * mrSges(3,2) + t386)) * t374; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t334 - mrSges(3,2) * t335 + t379 * (mrSges(4,2) * t329 - mrSges(4,3) * t316 + Ifges(4,1) * t359 + Ifges(4,4) * t360 + Ifges(4,5) * qJDD(3) - pkin(8) * t278 - qJD(3) * t349 + t382 * t273 - t378 * t274) + t383 * (-mrSges(4,1) * t329 + mrSges(4,3) * t317 + Ifges(4,4) * t359 + Ifges(4,2) * t360 + Ifges(4,6) * qJDD(3) - pkin(3) * t389 + pkin(8) * t395 + qJD(3) * t350 + t378 * t273 + t382 * t274) + pkin(2) * t386 + pkin(7) * t396; Ifges(4,5) * t359 + Ifges(4,6) * t360 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t316 - mrSges(4,2) * t317 + (t379 * t349 - t383 * t350) * qJD(2) + t388 + pkin(3) * t278; t388; t387;];
tauJ = t1;
