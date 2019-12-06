% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPPR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:56
% EndTime: 2019-12-05 17:28:57
% DurationCPUTime: 1.13s
% Computational Cost: add. (4605->176), mult. (10702->256), div. (0->0), fcn. (6478->10), ass. (0->94)
t370 = sin(qJ(1));
t372 = cos(qJ(1));
t400 = t372 * g(2) + t370 * g(3);
t345 = qJDD(1) * pkin(1) + t400;
t373 = qJD(1) ^ 2;
t392 = t370 * g(2) - g(3) * t372;
t346 = -pkin(1) * t373 + t392;
t365 = sin(pkin(7));
t368 = cos(pkin(7));
t328 = t365 * t345 + t368 * t346;
t409 = -pkin(2) * t373 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t328;
t327 = t345 * t368 - t365 * t346;
t378 = -qJ(3) * t373 + qJDD(3) - t327;
t364 = sin(pkin(8));
t367 = cos(pkin(8));
t386 = -pkin(3) * t367 - qJ(4) * t364;
t399 = qJD(1) * t364;
t408 = (-pkin(2) + t386) * qJDD(1) + t378 - 0.2e1 * qJD(4) * t399;
t362 = -g(1) + qJDD(2);
t313 = t362 * t367 - t409 * t364;
t407 = mrSges(4,2) * t364;
t366 = cos(pkin(9));
t406 = Ifges(5,4) * t366;
t363 = sin(pkin(9));
t405 = Ifges(5,2) * t363;
t361 = t364 ^ 2;
t404 = t361 * t373;
t403 = t363 * t364;
t402 = t364 * t366;
t314 = t364 * t362 + t409 * t367;
t342 = t386 * qJD(1);
t398 = qJD(1) * t367;
t307 = t342 * t398 + t314;
t384 = -pkin(4) * t367 - pkin(6) * t402;
t401 = t408 * t366;
t299 = t384 * qJDD(1) + (-t307 + (-pkin(4) * t361 * t366 + pkin(6) * t364 * t367) * t373) * t363 + t401;
t302 = t366 * t307 + t408 * t363;
t341 = t384 * qJD(1);
t394 = t363 ^ 2 * t404;
t396 = qJDD(1) * t364;
t300 = -pkin(6) * t363 * t396 - pkin(4) * t394 + t341 * t398 + t302;
t369 = sin(qJ(5));
t371 = cos(qJ(5));
t297 = t299 * t371 - t300 * t369;
t380 = (-t363 * t371 - t366 * t369) * t364;
t332 = qJD(1) * t380;
t379 = (-t363 * t369 + t366 * t371) * t364;
t333 = qJD(1) * t379;
t317 = -mrSges(6,1) * t332 + mrSges(6,2) * t333;
t320 = qJD(5) * t332 + qJDD(1) * t379;
t350 = qJD(5) - t398;
t325 = -mrSges(6,2) * t350 + mrSges(6,3) * t332;
t395 = qJDD(1) * t367;
t349 = qJDD(5) - t395;
t295 = m(6) * t297 + mrSges(6,1) * t349 - mrSges(6,3) * t320 - t317 * t333 + t325 * t350;
t298 = t299 * t369 + t300 * t371;
t319 = -qJD(5) * t333 + qJDD(1) * t380;
t326 = mrSges(6,1) * t350 - mrSges(6,3) * t333;
t296 = m(6) * t298 - mrSges(6,2) * t349 + mrSges(6,3) * t319 + t317 * t332 - t326 * t350;
t288 = t371 * t295 + t369 * t296;
t301 = -t307 * t363 + t401;
t387 = mrSges(5,1) * t363 + mrSges(5,2) * t366;
t334 = t387 * t399;
t381 = mrSges(5,2) * t367 - mrSges(5,3) * t403;
t336 = t381 * qJD(1);
t382 = -mrSges(5,1) * t367 - mrSges(5,3) * t402;
t286 = m(5) * t301 + t382 * qJDD(1) + (-t334 * t402 - t336 * t367) * qJD(1) + t288;
t337 = t382 * qJD(1);
t390 = -t369 * t295 + t371 * t296;
t287 = m(5) * t302 + t381 * qJDD(1) + (-t334 * t403 + t337 * t367) * qJD(1) + t390;
t343 = (-mrSges(4,1) * t367 + t407) * qJD(1);
t283 = m(4) * t314 - t286 * t363 + t287 * t366 + (qJDD(1) * mrSges(4,3) + qJD(1) * t343) * t367;
t306 = t342 * t399 + qJDD(4) - t313;
t304 = -pkin(6) * t394 + (pkin(4) * qJDD(1) * t363 + qJD(1) * t341 * t366) * t364 + t306;
t377 = -m(6) * t304 + mrSges(6,1) * t319 - t320 * mrSges(6,2) + t325 * t332 - t333 * t326;
t376 = m(5) * t306 - t377;
t385 = t336 * t363 + t337 * t366;
t291 = m(4) * t313 + ((-mrSges(4,3) - t387) * qJDD(1) + (-t343 - t385) * qJD(1)) * t364 - t376;
t391 = t367 * t283 - t291 * t364;
t285 = t286 * t366 + t287 * t363;
t383 = -Ifges(5,5) * t366 + Ifges(5,6) * t363 + Ifges(4,4);
t323 = -qJDD(1) * pkin(2) + t378;
t375 = -m(4) * t323 + mrSges(4,1) * t395 - t285 + (t367 ^ 2 * t373 + t404) * mrSges(4,3);
t309 = Ifges(6,4) * t333 + Ifges(6,2) * t332 + Ifges(6,6) * t350;
t310 = Ifges(6,1) * t333 + Ifges(6,4) * t332 + Ifges(6,5) * t350;
t374 = mrSges(6,1) * t297 - mrSges(6,2) * t298 + Ifges(6,5) * t320 + Ifges(6,6) * t319 + Ifges(6,3) * t349 + t333 * t309 - t332 * t310;
t344 = (Ifges(4,5) * t364 + Ifges(4,6) * t367) * qJD(1);
t331 = (-Ifges(5,5) * t367 + (Ifges(5,1) * t366 - Ifges(5,4) * t363) * t364) * qJD(1);
t330 = (-Ifges(5,6) * t367 + (-t405 + t406) * t364) * qJD(1);
t308 = Ifges(6,5) * t333 + Ifges(6,6) * t332 + Ifges(6,3) * t350;
t290 = mrSges(6,2) * t304 - mrSges(6,3) * t297 + Ifges(6,1) * t320 + Ifges(6,4) * t319 + Ifges(6,5) * t349 + t308 * t332 - t309 * t350;
t289 = -mrSges(6,1) * t304 + mrSges(6,3) * t298 + Ifges(6,4) * t320 + Ifges(6,2) * t319 + Ifges(6,6) * t349 - t308 * t333 + t310 * t350;
t284 = mrSges(4,2) * t396 - t375;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t400 - mrSges(2,2) * t392 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t327 - mrSges(3,2) * t328 + t364 * (t344 * t398 + mrSges(4,2) * t323 - mrSges(4,3) * t313 + t366 * (mrSges(5,2) * t306 - mrSges(5,3) * t301 - pkin(6) * t288 - t369 * t289 + t371 * t290 + t330 * t398) - t363 * (-mrSges(5,1) * t306 + mrSges(5,3) * t302 + pkin(4) * t377 + pkin(6) * t390 + t371 * t289 + t369 * t290 - t331 * t398) - qJ(4) * t285 + (t383 * t367 + (Ifges(5,1) * t366 ^ 2 + Ifges(4,1) + (t405 - 0.2e1 * t406) * t363) * t364) * qJDD(1)) + t367 * (-mrSges(4,1) * t323 - mrSges(5,1) * t301 + mrSges(5,2) * t302 + mrSges(4,3) * t314 - pkin(3) * t285 - pkin(4) * t288 + (Ifges(4,2) + Ifges(5,3)) * t395 + (t383 * qJDD(1) + (-t330 * t366 - t331 * t363 - t344) * qJD(1)) * t364 - t374) - pkin(2) * t284 + qJ(3) * t391 + pkin(1) * (t365 * (m(3) * t328 - mrSges(3,1) * t373 - qJDD(1) * mrSges(3,2) + t391) + t368 * (m(3) * t327 - mrSges(3,2) * t373 + (mrSges(3,1) - t407) * qJDD(1) + t375)); m(3) * t362 + t283 * t364 + t291 * t367; t284; (t385 * qJD(1) + t387 * qJDD(1)) * t364 + t376; t374;];
tauJ = t1;
