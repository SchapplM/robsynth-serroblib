% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:35:34
% EndTime: 2019-12-05 16:35:37
% DurationCPUTime: 1.24s
% Computational Cost: add. (6994->218), mult. (14286->283), div. (0->0), fcn. (9727->12), ass. (0->97)
t346 = sin(pkin(9));
t349 = cos(pkin(9));
t337 = t346 * g(1) - t349 * g(2);
t344 = -g(3) + qJDD(1);
t347 = sin(pkin(5));
t350 = cos(pkin(5));
t378 = t337 * t350 + t344 * t347;
t338 = -t349 * g(1) - t346 * g(2);
t353 = sin(qJ(2));
t356 = cos(qJ(2));
t302 = -t353 * t338 + t356 * t378;
t377 = 2 * qJD(4);
t355 = cos(qJ(3));
t358 = qJD(2) ^ 2;
t375 = t355 ^ 2 * t358;
t319 = -t347 * t337 + t350 * t344;
t373 = t355 * t319;
t303 = t356 * t338 + t378 * t353;
t299 = -t358 * pkin(2) + qJDD(2) * pkin(7) + t303;
t352 = sin(qJ(3));
t287 = t355 * t299 + t352 * t319;
t372 = qJD(2) * t355;
t371 = t352 * qJD(2);
t370 = qJD(2) * qJD(3);
t333 = (-t355 * pkin(3) - t352 * qJ(4)) * qJD(2);
t357 = qJD(3) ^ 2;
t283 = -t357 * pkin(3) + qJDD(3) * qJ(4) + t333 * t372 + t287;
t298 = -qJDD(2) * pkin(2) - t358 * pkin(7) - t302;
t368 = t355 * t370;
t335 = t352 * qJDD(2) + t368;
t369 = t352 * t370;
t336 = t355 * qJDD(2) - t369;
t285 = (-t335 - t368) * qJ(4) + (-t336 + t369) * pkin(3) + t298;
t345 = sin(pkin(10));
t348 = cos(pkin(10));
t328 = t348 * qJD(3) - t345 * t371;
t279 = t348 * t283 + t345 * t285 + t328 * t377;
t334 = (-t355 * mrSges(4,1) + t352 * mrSges(4,2)) * qJD(2);
t339 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t371;
t329 = t345 * qJD(3) + t348 * t371;
t308 = -t328 * mrSges(5,1) + t329 * mrSges(5,2);
t316 = -mrSges(5,1) * t372 - t329 * mrSges(5,3);
t317 = t348 * qJDD(3) - t345 * t335;
t309 = -t328 * pkin(4) - t329 * pkin(8);
t277 = -pkin(4) * t375 - t336 * pkin(8) + t328 * t309 + t279;
t318 = t345 * qJDD(3) + t348 * t335;
t296 = t352 * t299;
t362 = -qJDD(3) * pkin(3) - t357 * qJ(4) + t333 * t371 + qJDD(4) + t296;
t280 = -t317 * pkin(4) - t318 * pkin(8) + (-t319 + (-pkin(4) * t329 + pkin(8) * t328) * qJD(2)) * t355 + t362;
t351 = sin(qJ(5));
t354 = cos(qJ(5));
t274 = -t351 * t277 + t354 * t280;
t310 = -t351 * t329 - t354 * t372;
t294 = t310 * qJD(5) + t354 * t318 - t351 * t336;
t311 = t354 * t329 - t351 * t372;
t295 = -t310 * mrSges(6,1) + t311 * mrSges(6,2);
t326 = qJD(5) - t328;
t300 = -t326 * mrSges(6,2) + t310 * mrSges(6,3);
t314 = qJDD(5) - t317;
t272 = m(6) * t274 + t314 * mrSges(6,1) - t294 * mrSges(6,3) - t311 * t295 + t326 * t300;
t275 = t354 * t277 + t351 * t280;
t293 = -t311 * qJD(5) - t351 * t318 - t354 * t336;
t301 = t326 * mrSges(6,1) - t311 * mrSges(6,3);
t273 = m(6) * t275 - t314 * mrSges(6,2) + t293 * mrSges(6,3) + t310 * t295 - t326 * t301;
t365 = -t351 * t272 + t354 * t273;
t265 = m(5) * t279 + t336 * mrSges(5,2) + t317 * mrSges(5,3) + t328 * t308 + t316 * t372 + t365;
t364 = t345 * t283 - t348 * t285;
t278 = -0.2e1 * qJD(4) * t329 - t364;
t315 = mrSges(5,2) * t372 + t328 * mrSges(5,3);
t276 = -pkin(8) * t375 + t336 * pkin(4) + (t377 + t309) * t329 + t364;
t361 = -m(6) * t276 + t293 * mrSges(6,1) - t294 * mrSges(6,2) + t310 * t300 - t311 * t301;
t270 = m(5) * t278 - t336 * mrSges(5,1) - t318 * mrSges(5,3) - t329 * t308 - t315 * t372 + t361;
t366 = t348 * t265 - t345 * t270;
t261 = m(4) * t287 - qJDD(3) * mrSges(4,2) + t336 * mrSges(4,3) - qJD(3) * t339 + t334 * t372 + t366;
t267 = t354 * t272 + t351 * t273;
t282 = t362 - t373;
t266 = m(5) * t282 - t317 * mrSges(5,1) + t318 * mrSges(5,2) - t328 * t315 + t329 * t316 + t267;
t286 = -t296 + t373;
t340 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t372;
t263 = m(4) * t286 + qJDD(3) * mrSges(4,1) - t335 * mrSges(4,3) + qJD(3) * t340 - t334 * t371 - t266;
t367 = t355 * t261 - t352 * t263;
t262 = t345 * t265 + t348 * t270;
t360 = -m(4) * t298 + t336 * mrSges(4,1) - t335 * mrSges(4,2) - t339 * t371 + t340 * t372 - t262;
t289 = Ifges(6,4) * t311 + Ifges(6,2) * t310 + Ifges(6,6) * t326;
t290 = Ifges(6,1) * t311 + Ifges(6,4) * t310 + Ifges(6,5) * t326;
t359 = mrSges(6,1) * t274 - mrSges(6,2) * t275 + Ifges(6,5) * t294 + Ifges(6,6) * t293 + Ifges(6,3) * t314 + t311 * t289 - t310 * t290;
t325 = Ifges(4,5) * qJD(3) + (t352 * Ifges(4,1) + t355 * Ifges(4,4)) * qJD(2);
t324 = Ifges(4,6) * qJD(3) + (t352 * Ifges(4,4) + Ifges(4,2) * t355) * qJD(2);
t306 = Ifges(5,1) * t329 + Ifges(5,4) * t328 - Ifges(5,5) * t372;
t305 = Ifges(5,4) * t329 + Ifges(5,2) * t328 - Ifges(5,6) * t372;
t304 = Ifges(5,5) * t329 + Ifges(5,6) * t328 - Ifges(5,3) * t372;
t288 = Ifges(6,5) * t311 + Ifges(6,6) * t310 + Ifges(6,3) * t326;
t269 = mrSges(6,2) * t276 - mrSges(6,3) * t274 + Ifges(6,1) * t294 + Ifges(6,4) * t293 + Ifges(6,5) * t314 + t310 * t288 - t326 * t289;
t268 = -mrSges(6,1) * t276 + mrSges(6,3) * t275 + Ifges(6,4) * t294 + Ifges(6,2) * t293 + Ifges(6,6) * t314 - t311 * t288 + t326 * t290;
t259 = -mrSges(5,1) * t282 + mrSges(5,3) * t279 + Ifges(5,4) * t318 + Ifges(5,2) * t317 - Ifges(5,6) * t336 - pkin(4) * t267 - t329 * t304 - t306 * t372 - t359;
t258 = mrSges(5,2) * t282 - mrSges(5,3) * t278 + Ifges(5,1) * t318 + Ifges(5,4) * t317 - Ifges(5,5) * t336 - pkin(8) * t267 - t351 * t268 + t354 * t269 + t328 * t304 + t305 * t372;
t1 = [m(2) * t344 + t350 * (m(3) * t319 + t352 * t261 + t355 * t263) + (t353 * (m(3) * t303 - t358 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t367) + t356 * (m(3) * t302 + qJDD(2) * mrSges(3,1) - t358 * mrSges(3,2) + t360)) * t347; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t302 - mrSges(3,2) * t303 + t352 * (mrSges(4,2) * t298 - mrSges(4,3) * t286 + Ifges(4,1) * t335 + Ifges(4,4) * t336 + Ifges(4,5) * qJDD(3) - qJ(4) * t262 - qJD(3) * t324 + t348 * t258 - t345 * t259) + t355 * (Ifges(4,4) * t335 + Ifges(4,6) * qJDD(3) + qJD(3) * t325 - mrSges(4,1) * t298 + mrSges(4,3) * t287 - Ifges(5,5) * t318 - Ifges(5,6) * t317 - t329 * t305 + t328 * t306 - mrSges(5,1) * t278 + mrSges(5,2) * t279 - t351 * t269 - t354 * t268 - pkin(4) * t361 - pkin(8) * t365 - pkin(3) * t262 + (Ifges(4,2) + Ifges(5,3)) * t336) + pkin(2) * t360 + pkin(7) * t367; Ifges(4,5) * t335 + Ifges(4,6) * t336 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t286 - mrSges(4,2) * t287 + t345 * t258 + t348 * t259 - pkin(3) * t266 + qJ(4) * t366 + (t352 * t324 - t355 * t325) * qJD(2); t266; t359;];
tauJ = t1;
