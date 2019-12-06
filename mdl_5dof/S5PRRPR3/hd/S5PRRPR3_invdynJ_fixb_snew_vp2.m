% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:22
% EndTime: 2019-12-05 16:19:23
% DurationCPUTime: 0.96s
% Computational Cost: add. (7100->206), mult. (15692->272), div. (0->0), fcn. (10526->10), ass. (0->86)
t358 = qJD(2) ^ 2;
t349 = sin(pkin(8));
t351 = cos(pkin(8));
t336 = t349 * g(1) - t351 * g(2);
t337 = -t351 * g(1) - t349 * g(2);
t354 = sin(qJ(2));
t357 = cos(qJ(2));
t369 = t354 * t336 + t357 * t337;
t315 = -t358 * pkin(2) + qJDD(2) * pkin(6) + t369;
t347 = -g(3) + qJDD(1);
t353 = sin(qJ(3));
t356 = cos(qJ(3));
t305 = -t353 * t315 + t356 * t347;
t366 = qJD(2) * qJD(3);
t365 = t356 * t366;
t334 = t353 * qJDD(2) + t365;
t300 = (-t334 + t365) * qJ(4) + (t353 * t356 * t358 + qJDD(3)) * pkin(3) + t305;
t306 = t356 * t315 + t353 * t347;
t335 = t356 * qJDD(2) - t353 * t366;
t367 = t353 * qJD(2);
t338 = qJD(3) * pkin(3) - qJ(4) * t367;
t346 = t356 ^ 2;
t301 = -t346 * t358 * pkin(3) + t335 * qJ(4) - qJD(3) * t338 + t306;
t348 = sin(pkin(9));
t350 = cos(pkin(9));
t325 = (t356 * t348 + t353 * t350) * qJD(2);
t282 = -0.2e1 * qJD(4) * t325 + t350 * t300 - t348 * t301;
t317 = t350 * t334 + t348 * t335;
t324 = (-t353 * t348 + t356 * t350) * qJD(2);
t280 = (qJD(3) * t324 - t317) * pkin(7) + (t324 * t325 + qJDD(3)) * pkin(4) + t282;
t283 = 0.2e1 * qJD(4) * t324 + t348 * t300 + t350 * t301;
t316 = -t348 * t334 + t350 * t335;
t320 = qJD(3) * pkin(4) - t325 * pkin(7);
t323 = t324 ^ 2;
t281 = -t323 * pkin(4) + t316 * pkin(7) - qJD(3) * t320 + t283;
t352 = sin(qJ(5));
t355 = cos(qJ(5));
t278 = t355 * t280 - t352 * t281;
t310 = t355 * t324 - t352 * t325;
t291 = t310 * qJD(5) + t352 * t316 + t355 * t317;
t311 = t352 * t324 + t355 * t325;
t296 = -t310 * mrSges(6,1) + t311 * mrSges(6,2);
t345 = qJD(3) + qJD(5);
t303 = -t345 * mrSges(6,2) + t310 * mrSges(6,3);
t344 = qJDD(3) + qJDD(5);
t274 = m(6) * t278 + t344 * mrSges(6,1) - t291 * mrSges(6,3) - t311 * t296 + t345 * t303;
t279 = t352 * t280 + t355 * t281;
t290 = -t311 * qJD(5) + t355 * t316 - t352 * t317;
t304 = t345 * mrSges(6,1) - t311 * mrSges(6,3);
t275 = m(6) * t279 - t344 * mrSges(6,2) + t290 * mrSges(6,3) + t310 * t296 - t345 * t304;
t268 = t355 * t274 + t352 * t275;
t312 = -t324 * mrSges(5,1) + t325 * mrSges(5,2);
t318 = -qJD(3) * mrSges(5,2) + t324 * mrSges(5,3);
t266 = m(5) * t282 + qJDD(3) * mrSges(5,1) - t317 * mrSges(5,3) + qJD(3) * t318 - t325 * t312 + t268;
t319 = qJD(3) * mrSges(5,1) - t325 * mrSges(5,3);
t363 = -t352 * t274 + t355 * t275;
t267 = m(5) * t283 - qJDD(3) * mrSges(5,2) + t316 * mrSges(5,3) - qJD(3) * t319 + t324 * t312 + t363;
t262 = t350 * t266 + t348 * t267;
t368 = qJD(2) * t356;
t364 = -t348 * t266 + t350 * t267;
t362 = t357 * t336 - t354 * t337;
t361 = -qJDD(2) * pkin(2) - t362;
t302 = -t335 * pkin(3) + qJDD(4) + t338 * t367 + (-qJ(4) * t346 - pkin(6)) * t358 + t361;
t285 = -t316 * pkin(4) - t323 * pkin(7) + t325 * t320 + t302;
t360 = m(6) * t285 - t290 * mrSges(6,1) + t291 * mrSges(6,2) - t310 * t303 + t311 * t304;
t293 = Ifges(6,4) * t311 + Ifges(6,2) * t310 + Ifges(6,6) * t345;
t294 = Ifges(6,1) * t311 + Ifges(6,4) * t310 + Ifges(6,5) * t345;
t359 = mrSges(6,1) * t278 - mrSges(6,2) * t279 + Ifges(6,5) * t291 + Ifges(6,6) * t290 + Ifges(6,3) * t344 + t311 * t293 - t310 * t294;
t276 = m(5) * t302 - t316 * mrSges(5,1) + t317 * mrSges(5,2) - t324 * t318 + t325 * t319 + t360;
t340 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t368;
t339 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t367;
t333 = (-t356 * mrSges(4,1) + t353 * mrSges(4,2)) * qJD(2);
t328 = Ifges(4,5) * qJD(3) + (t353 * Ifges(4,1) + t356 * Ifges(4,4)) * qJD(2);
t327 = Ifges(4,6) * qJD(3) + (t353 * Ifges(4,4) + t356 * Ifges(4,2)) * qJD(2);
t314 = -t358 * pkin(6) + t361;
t309 = Ifges(5,1) * t325 + Ifges(5,4) * t324 + Ifges(5,5) * qJD(3);
t308 = Ifges(5,4) * t325 + Ifges(5,2) * t324 + Ifges(5,6) * qJD(3);
t307 = Ifges(5,5) * t325 + Ifges(5,6) * t324 + Ifges(5,3) * qJD(3);
t292 = Ifges(6,5) * t311 + Ifges(6,6) * t310 + Ifges(6,3) * t345;
t270 = mrSges(6,2) * t285 - mrSges(6,3) * t278 + Ifges(6,1) * t291 + Ifges(6,4) * t290 + Ifges(6,5) * t344 + t310 * t292 - t345 * t293;
t269 = -mrSges(6,1) * t285 + mrSges(6,3) * t279 + Ifges(6,4) * t291 + Ifges(6,2) * t290 + Ifges(6,6) * t344 - t311 * t292 + t345 * t294;
t261 = m(4) * t306 - qJDD(3) * mrSges(4,2) + t335 * mrSges(4,3) - qJD(3) * t339 + t333 * t368 + t364;
t260 = m(4) * t305 + qJDD(3) * mrSges(4,1) - t334 * mrSges(4,3) + qJD(3) * t340 - t333 * t367 + t262;
t259 = mrSges(5,2) * t302 - mrSges(5,3) * t282 + Ifges(5,1) * t317 + Ifges(5,4) * t316 + Ifges(5,5) * qJDD(3) - pkin(7) * t268 - qJD(3) * t308 - t352 * t269 + t355 * t270 + t324 * t307;
t258 = -mrSges(5,1) * t302 + mrSges(5,3) * t283 + Ifges(5,4) * t317 + Ifges(5,2) * t316 + Ifges(5,6) * qJDD(3) - pkin(4) * t360 + pkin(7) * t363 + qJD(3) * t309 + t355 * t269 + t352 * t270 - t325 * t307;
t1 = [t356 * t260 + t353 * t261 + (m(2) + m(3)) * t347; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t362 - mrSges(3,2) * t369 + t353 * (mrSges(4,2) * t314 - mrSges(4,3) * t305 + Ifges(4,1) * t334 + Ifges(4,4) * t335 + Ifges(4,5) * qJDD(3) - qJ(4) * t262 - qJD(3) * t327 - t348 * t258 + t350 * t259) + t356 * (-mrSges(4,1) * t314 + mrSges(4,3) * t306 + Ifges(4,4) * t334 + Ifges(4,2) * t335 + Ifges(4,6) * qJDD(3) - pkin(3) * t276 + qJ(4) * t364 + qJD(3) * t328 + t350 * t258 + t348 * t259) + pkin(6) * (-t353 * t260 + t356 * t261) + (-m(4) * t314 + t335 * mrSges(4,1) - t334 * mrSges(4,2) - t276 + (-t339 * t353 + t340 * t356) * qJD(2)) * pkin(2); (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t353 * t327 - t356 * t328) * qJD(2) + t359 + Ifges(4,5) * t334 + Ifges(4,6) * t335 - t324 * t309 + t325 * t308 + Ifges(5,6) * t316 + Ifges(5,5) * t317 + mrSges(4,1) * t305 - mrSges(4,2) * t306 + mrSges(5,1) * t282 - mrSges(5,2) * t283 + pkin(4) * t268 + pkin(3) * t262; t276; t359;];
tauJ = t1;
