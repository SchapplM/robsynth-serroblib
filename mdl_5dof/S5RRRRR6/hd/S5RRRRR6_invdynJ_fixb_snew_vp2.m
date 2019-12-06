% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:21
% EndTime: 2019-12-05 19:00:22
% DurationCPUTime: 1.63s
% Computational Cost: add. (23359->223), mult. (29931->289), div. (0->0), fcn. (19448->10), ass. (0->95)
t378 = qJD(1) + qJD(2);
t374 = t378 ^ 2;
t406 = pkin(3) * t374;
t382 = sin(qJ(3));
t405 = t378 * t382;
t387 = cos(qJ(3));
t404 = t378 * t387;
t384 = sin(qJ(1));
t389 = cos(qJ(1));
t402 = t389 * g(2) + t384 * g(3);
t361 = qJDD(1) * pkin(1) + t402;
t400 = t384 * g(2) - t389 * g(3);
t362 = -qJD(1) ^ 2 * pkin(1) + t400;
t383 = sin(qJ(2));
t388 = cos(qJ(2));
t342 = t383 * t361 + t388 * t362;
t376 = qJDD(1) + qJDD(2);
t339 = -t374 * pkin(2) + t376 * pkin(7) + t342;
t403 = t382 * t339;
t401 = qJD(3) * t378;
t356 = t382 * t376 + t387 * t401;
t317 = qJDD(3) * pkin(3) - t356 * pkin(8) - t403 + (pkin(8) * t401 + t382 * t406 - g(1)) * t387;
t329 = -t382 * g(1) + t387 * t339;
t357 = t387 * t376 - t382 * t401;
t365 = qJD(3) * pkin(3) - pkin(8) * t405;
t379 = t387 ^ 2;
t318 = t357 * pkin(8) - qJD(3) * t365 - t379 * t406 + t329;
t381 = sin(qJ(4));
t386 = cos(qJ(4));
t299 = t386 * t317 - t381 * t318;
t350 = (-t381 * t382 + t386 * t387) * t378;
t325 = t350 * qJD(4) + t386 * t356 + t381 * t357;
t351 = (t381 * t387 + t382 * t386) * t378;
t375 = qJDD(3) + qJDD(4);
t377 = qJD(3) + qJD(4);
t294 = (t350 * t377 - t325) * pkin(9) + (t350 * t351 + t375) * pkin(4) + t299;
t300 = t381 * t317 + t386 * t318;
t324 = -t351 * qJD(4) - t381 * t356 + t386 * t357;
t345 = t377 * pkin(4) - t351 * pkin(9);
t346 = t350 ^ 2;
t295 = -t346 * pkin(4) + t324 * pkin(9) - t377 * t345 + t300;
t380 = sin(qJ(5));
t385 = cos(qJ(5));
t292 = t385 * t294 - t380 * t295;
t334 = t385 * t350 - t380 * t351;
t306 = t334 * qJD(5) + t380 * t324 + t385 * t325;
t335 = t380 * t350 + t385 * t351;
t313 = -t334 * mrSges(6,1) + t335 * mrSges(6,2);
t370 = qJD(5) + t377;
t326 = -t370 * mrSges(6,2) + t334 * mrSges(6,3);
t369 = qJDD(5) + t375;
t289 = m(6) * t292 + t369 * mrSges(6,1) - t306 * mrSges(6,3) - t335 * t313 + t370 * t326;
t293 = t380 * t294 + t385 * t295;
t305 = -t335 * qJD(5) + t385 * t324 - t380 * t325;
t327 = t370 * mrSges(6,1) - t335 * mrSges(6,3);
t290 = m(6) * t293 - t369 * mrSges(6,2) + t305 * mrSges(6,3) + t334 * t313 - t370 * t327;
t282 = t385 * t289 + t380 * t290;
t337 = -t350 * mrSges(5,1) + t351 * mrSges(5,2);
t343 = -t377 * mrSges(5,2) + t350 * mrSges(5,3);
t279 = m(5) * t299 + t375 * mrSges(5,1) - t325 * mrSges(5,3) - t351 * t337 + t377 * t343 + t282;
t344 = t377 * mrSges(5,1) - t351 * mrSges(5,3);
t397 = -t380 * t289 + t385 * t290;
t280 = m(5) * t300 - t375 * mrSges(5,2) + t324 * mrSges(5,3) + t350 * t337 - t377 * t344 + t397;
t275 = t386 * t279 + t381 * t280;
t328 = -t387 * g(1) - t403;
t355 = (-mrSges(4,1) * t387 + mrSges(4,2) * t382) * t378;
t363 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t405;
t364 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t404;
t398 = -t381 * t279 + t386 * t280;
t399 = -t382 * (m(4) * t328 + qJDD(3) * mrSges(4,1) - t356 * mrSges(4,3) + qJD(3) * t364 - t355 * t405 + t275) + t387 * (m(4) * t329 - qJDD(3) * mrSges(4,2) + t357 * mrSges(4,3) - qJD(3) * t363 + t355 * t404 + t398);
t341 = t388 * t361 - t383 * t362;
t395 = -t376 * pkin(2) - t341;
t319 = -t357 * pkin(3) + t365 * t405 + (-pkin(8) * t379 - pkin(7)) * t374 + t395;
t297 = -t324 * pkin(4) - t346 * pkin(9) + t351 * t345 + t319;
t396 = m(6) * t297 - t305 * mrSges(6,1) + t306 * mrSges(6,2) - t334 * t326 + t335 * t327;
t308 = Ifges(6,5) * t335 + Ifges(6,6) * t334 + Ifges(6,3) * t370;
t310 = Ifges(6,1) * t335 + Ifges(6,4) * t334 + Ifges(6,5) * t370;
t283 = -mrSges(6,1) * t297 + mrSges(6,3) * t293 + Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t369 - t335 * t308 + t370 * t310;
t309 = Ifges(6,4) * t335 + Ifges(6,2) * t334 + Ifges(6,6) * t370;
t284 = mrSges(6,2) * t297 - mrSges(6,3) * t292 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t369 + t334 * t308 - t370 * t309;
t330 = Ifges(5,5) * t351 + Ifges(5,6) * t350 + Ifges(5,3) * t377;
t332 = Ifges(5,1) * t351 + Ifges(5,4) * t350 + Ifges(5,5) * t377;
t271 = -mrSges(5,1) * t319 + mrSges(5,3) * t300 + Ifges(5,4) * t325 + Ifges(5,2) * t324 + Ifges(5,6) * t375 - pkin(4) * t396 + pkin(9) * t397 + t385 * t283 + t380 * t284 - t351 * t330 + t377 * t332;
t331 = Ifges(5,4) * t351 + Ifges(5,2) * t350 + Ifges(5,6) * t377;
t272 = mrSges(5,2) * t319 - mrSges(5,3) * t299 + Ifges(5,1) * t325 + Ifges(5,4) * t324 + Ifges(5,5) * t375 - pkin(9) * t282 - t380 * t283 + t385 * t284 + t350 * t330 - t377 * t331;
t338 = -t374 * pkin(7) + t395;
t347 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t382 + Ifges(4,6) * t387) * t378;
t348 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t382 + Ifges(4,2) * t387) * t378;
t349 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t382 + Ifges(4,4) * t387) * t378;
t392 = m(5) * t319 - t324 * mrSges(5,1) + t325 * mrSges(5,2) - t350 * t343 + t351 * t344 + t396;
t390 = -m(4) * t338 + t357 * mrSges(4,1) - t356 * mrSges(4,2) - t363 * t405 + t364 * t404 - t392;
t394 = -mrSges(3,2) * t342 + t387 * (-mrSges(4,1) * t338 + mrSges(4,3) * t329 + Ifges(4,4) * t356 + Ifges(4,2) * t357 + Ifges(4,6) * qJDD(3) - pkin(3) * t392 + pkin(8) * t398 + qJD(3) * t349 + t386 * t271 + t381 * t272 - t347 * t405) + t382 * (mrSges(4,2) * t338 - mrSges(4,3) * t328 + Ifges(4,1) * t356 + Ifges(4,4) * t357 + Ifges(4,5) * qJDD(3) - pkin(8) * t275 - qJD(3) * t348 - t381 * t271 + t386 * t272 + t347 * t404) + pkin(7) * t399 + pkin(2) * t390 + mrSges(3,1) * t341 + Ifges(3,3) * t376;
t393 = mrSges(6,1) * t292 - mrSges(6,2) * t293 + Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t369 + t335 * t309 - t334 * t310;
t391 = mrSges(5,1) * t299 - mrSges(5,2) * t300 + Ifges(5,5) * t325 + Ifges(5,6) * t324 + Ifges(5,3) * t375 + pkin(4) * t282 + t351 * t331 - t350 * t332 + t393;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t402 - mrSges(2,2) * t400 + pkin(1) * (t383 * (m(3) * t342 - t374 * mrSges(3,1) - t376 * mrSges(3,2) + t399) + t388 * (m(3) * t341 + t376 * mrSges(3,1) - t374 * mrSges(3,2) + t390)) + t394; t394; t391 + Ifges(4,3) * qJDD(3) + (t382 * t348 - t387 * t349) * t378 + Ifges(4,5) * t356 + Ifges(4,6) * t357 + mrSges(4,1) * t328 - mrSges(4,2) * t329 + pkin(3) * t275; t391; t393;];
tauJ = t1;
