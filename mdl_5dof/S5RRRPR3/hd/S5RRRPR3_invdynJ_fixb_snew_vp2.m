% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR3
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:36
% EndTime: 2019-12-05 18:42:38
% DurationCPUTime: 1.47s
% Computational Cost: add. (20826->222), mult. (27803->290), div. (0->0), fcn. (17752->10), ass. (0->92)
t373 = qJD(1) + qJD(2);
t369 = t373 ^ 2;
t399 = pkin(3) * t369;
t378 = sin(qJ(3));
t398 = t373 * t378;
t382 = cos(qJ(3));
t397 = t373 * t382;
t380 = sin(qJ(1));
t384 = cos(qJ(1));
t395 = t384 * g(2) + t380 * g(3);
t359 = qJDD(1) * pkin(1) + t395;
t393 = t380 * g(2) - t384 * g(3);
t360 = -qJD(1) ^ 2 * pkin(1) + t393;
t379 = sin(qJ(2));
t383 = cos(qJ(2));
t338 = t379 * t359 + t383 * t360;
t371 = qJDD(1) + qJDD(2);
t332 = -t369 * pkin(2) + t371 * pkin(7) + t338;
t396 = t378 * t332;
t394 = qJD(3) * t373;
t354 = t378 * t371 + t382 * t394;
t316 = qJDD(3) * pkin(3) - t354 * qJ(4) - t396 + (qJ(4) * t394 + t378 * t399 - g(1)) * t382;
t322 = -t378 * g(1) + t382 * t332;
t355 = t382 * t371 - t378 * t394;
t361 = qJD(3) * pkin(3) - qJ(4) * t398;
t374 = t382 ^ 2;
t317 = t355 * qJ(4) - qJD(3) * t361 - t374 * t399 + t322;
t375 = sin(pkin(9));
t376 = cos(pkin(9));
t346 = (t375 * t382 + t376 * t378) * t373;
t296 = -0.2e1 * qJD(4) * t346 + t376 * t316 - t375 * t317;
t335 = t376 * t354 + t375 * t355;
t345 = (-t375 * t378 + t376 * t382) * t373;
t294 = (qJD(3) * t345 - t335) * pkin(8) + (t345 * t346 + qJDD(3)) * pkin(4) + t296;
t297 = 0.2e1 * qJD(4) * t345 + t375 * t316 + t376 * t317;
t334 = -t375 * t354 + t376 * t355;
t341 = qJD(3) * pkin(4) - t346 * pkin(8);
t344 = t345 ^ 2;
t295 = -t344 * pkin(4) + t334 * pkin(8) - qJD(3) * t341 + t297;
t377 = sin(qJ(5));
t381 = cos(qJ(5));
t292 = t381 * t294 - t377 * t295;
t326 = t381 * t345 - t377 * t346;
t306 = t326 * qJD(5) + t377 * t334 + t381 * t335;
t327 = t377 * t345 + t381 * t346;
t312 = -t326 * mrSges(6,1) + t327 * mrSges(6,2);
t372 = qJD(3) + qJD(5);
t319 = -t372 * mrSges(6,2) + t326 * mrSges(6,3);
t370 = qJDD(3) + qJDD(5);
t288 = m(6) * t292 + t370 * mrSges(6,1) - t306 * mrSges(6,3) - t327 * t312 + t372 * t319;
t293 = t377 * t294 + t381 * t295;
t305 = -t327 * qJD(5) + t381 * t334 - t377 * t335;
t320 = t372 * mrSges(6,1) - t327 * mrSges(6,3);
t289 = m(6) * t293 - t370 * mrSges(6,2) + t305 * mrSges(6,3) + t326 * t312 - t372 * t320;
t281 = t381 * t288 + t377 * t289;
t330 = -t345 * mrSges(5,1) + t346 * mrSges(5,2);
t339 = -qJD(3) * mrSges(5,2) + t345 * mrSges(5,3);
t279 = m(5) * t296 + qJDD(3) * mrSges(5,1) - t335 * mrSges(5,3) + qJD(3) * t339 - t346 * t330 + t281;
t340 = qJD(3) * mrSges(5,1) - t346 * mrSges(5,3);
t390 = -t377 * t288 + t381 * t289;
t280 = m(5) * t297 - qJDD(3) * mrSges(5,2) + t334 * mrSges(5,3) - qJD(3) * t340 + t345 * t330 + t390;
t275 = t376 * t279 + t375 * t280;
t321 = -t382 * g(1) - t396;
t353 = (-mrSges(4,1) * t382 + mrSges(4,2) * t378) * t373;
t362 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t398;
t363 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t397;
t391 = -t375 * t279 + t376 * t280;
t392 = -t378 * (m(4) * t321 + qJDD(3) * mrSges(4,1) - t354 * mrSges(4,3) + qJD(3) * t363 - t353 * t398 + t275) + t382 * (m(4) * t322 - qJDD(3) * mrSges(4,2) + t355 * mrSges(4,3) - qJD(3) * t362 + t353 * t397 + t391);
t337 = t383 * t359 - t379 * t360;
t388 = -t371 * pkin(2) - t337;
t318 = -t355 * pkin(3) + qJDD(4) + t361 * t398 + (-qJ(4) * t374 - pkin(7)) * t369 + t388;
t299 = -t334 * pkin(4) - t344 * pkin(8) + t346 * t341 + t318;
t389 = m(6) * t299 - t305 * mrSges(6,1) + t306 * mrSges(6,2) - t326 * t319 + t327 * t320;
t307 = Ifges(6,5) * t327 + Ifges(6,6) * t326 + Ifges(6,3) * t372;
t309 = Ifges(6,1) * t327 + Ifges(6,4) * t326 + Ifges(6,5) * t372;
t282 = -mrSges(6,1) * t299 + mrSges(6,3) * t293 + Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t370 - t327 * t307 + t372 * t309;
t308 = Ifges(6,4) * t327 + Ifges(6,2) * t326 + Ifges(6,6) * t372;
t283 = mrSges(6,2) * t299 - mrSges(6,3) * t292 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t370 + t326 * t307 - t372 * t308;
t323 = Ifges(5,5) * t346 + Ifges(5,6) * t345 + Ifges(5,3) * qJD(3);
t325 = Ifges(5,1) * t346 + Ifges(5,4) * t345 + Ifges(5,5) * qJD(3);
t271 = -mrSges(5,1) * t318 + mrSges(5,3) * t297 + Ifges(5,4) * t335 + Ifges(5,2) * t334 + Ifges(5,6) * qJDD(3) - pkin(4) * t389 + pkin(8) * t390 + qJD(3) * t325 + t381 * t282 + t377 * t283 - t346 * t323;
t324 = Ifges(5,4) * t346 + Ifges(5,2) * t345 + Ifges(5,6) * qJD(3);
t272 = mrSges(5,2) * t318 - mrSges(5,3) * t296 + Ifges(5,1) * t335 + Ifges(5,4) * t334 + Ifges(5,5) * qJDD(3) - pkin(8) * t281 - qJD(3) * t324 - t377 * t282 + t381 * t283 + t345 * t323;
t290 = m(5) * t318 - t334 * mrSges(5,1) + t335 * mrSges(5,2) - t345 * t339 + t346 * t340 + t389;
t331 = -t369 * pkin(7) + t388;
t347 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t378 + Ifges(4,6) * t382) * t373;
t348 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t378 + Ifges(4,2) * t382) * t373;
t349 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t378 + Ifges(4,4) * t382) * t373;
t385 = -m(4) * t331 + t355 * mrSges(4,1) - t354 * mrSges(4,2) - t362 * t398 + t363 * t397 - t290;
t387 = -mrSges(3,2) * t338 + t382 * (-mrSges(4,1) * t331 + mrSges(4,3) * t322 + Ifges(4,4) * t354 + Ifges(4,2) * t355 + Ifges(4,6) * qJDD(3) - pkin(3) * t290 + qJ(4) * t391 + qJD(3) * t349 + t376 * t271 + t375 * t272 - t347 * t398) + t378 * (mrSges(4,2) * t331 - mrSges(4,3) * t321 + Ifges(4,1) * t354 + Ifges(4,4) * t355 + Ifges(4,5) * qJDD(3) - qJ(4) * t275 - qJD(3) * t348 - t375 * t271 + t376 * t272 + t347 * t397) + pkin(7) * t392 + pkin(2) * t385 + mrSges(3,1) * t337 + Ifges(3,3) * t371;
t386 = mrSges(6,1) * t292 - mrSges(6,2) * t293 + Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t370 + t327 * t308 - t326 * t309;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t395 - mrSges(2,2) * t393 + pkin(1) * (t379 * (m(3) * t338 - t369 * mrSges(3,1) - t371 * mrSges(3,2) + t392) + t383 * (m(3) * t337 + t371 * mrSges(3,1) - t369 * mrSges(3,2) + t385)) + t387; t387; t386 + mrSges(5,1) * t296 - mrSges(5,2) * t297 + mrSges(4,1) * t321 - mrSges(4,2) * t322 + Ifges(5,6) * t334 + Ifges(5,5) * t335 - t345 * t325 + t346 * t324 + Ifges(4,5) * t354 + Ifges(4,6) * t355 + pkin(4) * t281 + pkin(3) * t275 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t378 * t348 - t382 * t349) * t373; t290; t386;];
tauJ = t1;
