% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR15
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR15_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:49
% EndTime: 2019-12-31 18:36:50
% DurationCPUTime: 1.02s
% Computational Cost: add. (6792->214), mult. (14031->273), div. (0->0), fcn. (8380->8), ass. (0->86)
t347 = sin(qJ(1));
t350 = cos(qJ(1));
t358 = -t350 * g(1) - t347 * g(2);
t369 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t358;
t368 = (-pkin(1) - pkin(6));
t352 = qJD(1) ^ 2;
t314 = (t368 * t352) - t369;
t346 = sin(qJ(3));
t349 = cos(qJ(3));
t365 = qJD(1) * qJD(3);
t362 = t349 * t365;
t335 = t346 * qJDD(1) + t362;
t363 = t346 * t365;
t336 = t349 * qJDD(1) - t363;
t297 = (-t336 + t363) * qJ(4) + (t335 + t362) * pkin(3) + t314;
t361 = t347 * g(1) - t350 * g(2);
t355 = -t352 * qJ(2) + qJDD(2) - t361;
t320 = t368 * qJDD(1) + t355;
t311 = -t349 * g(3) + t346 * t320;
t333 = (t346 * pkin(3) - t349 * qJ(4)) * qJD(1);
t351 = qJD(3) ^ 2;
t367 = t346 * qJD(1);
t300 = -t351 * pkin(3) + qJDD(3) * qJ(4) - t333 * t367 + t311;
t343 = sin(pkin(8));
t344 = cos(pkin(8));
t366 = t349 * qJD(1);
t331 = t343 * qJD(3) + t344 * t366;
t284 = -0.2e1 * qJD(4) * t331 + t344 * t297 - t343 * t300;
t318 = t343 * qJDD(3) + t344 * t336;
t330 = t344 * qJD(3) - t343 * t366;
t282 = (t330 * t367 - t318) * pkin(7) + (t330 * t331 + t335) * pkin(4) + t284;
t285 = 0.2e1 * qJD(4) * t330 + t343 * t297 + t344 * t300;
t317 = t344 * qJDD(3) - t343 * t336;
t319 = pkin(4) * t367 - t331 * pkin(7);
t329 = t330 ^ 2;
t283 = -t329 * pkin(4) + t317 * pkin(7) - t319 * t367 + t285;
t345 = sin(qJ(5));
t348 = cos(qJ(5));
t280 = t348 * t282 - t345 * t283;
t307 = t348 * t330 - t345 * t331;
t289 = t307 * qJD(5) + t345 * t317 + t348 * t318;
t308 = t345 * t330 + t348 * t331;
t294 = -t307 * mrSges(6,1) + t308 * mrSges(6,2);
t339 = qJD(5) + t367;
t301 = -t339 * mrSges(6,2) + t307 * mrSges(6,3);
t332 = qJDD(5) + t335;
t277 = m(6) * t280 + t332 * mrSges(6,1) - t289 * mrSges(6,3) - t308 * t294 + t339 * t301;
t281 = t345 * t282 + t348 * t283;
t288 = -t308 * qJD(5) + t348 * t317 - t345 * t318;
t302 = t339 * mrSges(6,1) - t308 * mrSges(6,3);
t278 = m(6) * t281 - t332 * mrSges(6,2) + t288 * mrSges(6,3) + t307 * t294 - t339 * t302;
t270 = t348 * t277 + t345 * t278;
t309 = -t330 * mrSges(5,1) + t331 * mrSges(5,2);
t315 = -mrSges(5,2) * t367 + t330 * mrSges(5,3);
t268 = m(5) * t284 + t335 * mrSges(5,1) - t318 * mrSges(5,3) - t331 * t309 + t315 * t367 + t270;
t316 = mrSges(5,1) * t367 - t331 * mrSges(5,3);
t359 = -t345 * t277 + t348 * t278;
t269 = m(5) * t285 - t335 * mrSges(5,2) + t317 * mrSges(5,3) + t330 * t309 - t316 * t367 + t359;
t264 = t344 * t268 + t343 * t269;
t360 = -t343 * t268 + t344 * t269;
t310 = t346 * g(3) + t349 * t320;
t299 = -qJDD(3) * pkin(3) - t351 * qJ(4) + t333 * t366 + qJDD(4) - t310;
t286 = -t317 * pkin(4) - t329 * pkin(7) + t331 * t319 + t299;
t354 = m(6) * t286 - t288 * mrSges(6,1) + t289 * mrSges(6,2) - t307 * t301 + t308 * t302;
t279 = m(5) * t299 - t317 * mrSges(5,1) + t318 * mrSges(5,2) - t330 * t315 + t331 * t316 + t354;
t334 = (t346 * mrSges(4,1) + t349 * mrSges(4,2)) * qJD(1);
t337 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t367;
t338 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t366;
t357 = t346 * (m(4) * t311 - qJDD(3) * mrSges(4,2) - t335 * mrSges(4,3) - qJD(3) * t338 - t334 * t367 + t360) + t349 * (m(4) * t310 + qJDD(3) * mrSges(4,1) - t336 * mrSges(4,3) + qJD(3) * t337 - t334 * t366 - t279);
t291 = Ifges(6,4) * t308 + Ifges(6,2) * t307 + Ifges(6,6) * t339;
t292 = Ifges(6,1) * t308 + Ifges(6,4) * t307 + Ifges(6,5) * t339;
t353 = mrSges(6,1) * t280 - mrSges(6,2) * t281 + Ifges(6,5) * t289 + Ifges(6,6) * t288 + Ifges(6,3) * t332 + t308 * t291 - t307 * t292;
t328 = (Ifges(4,5) * qJD(3)) + (t349 * Ifges(4,1) - t346 * Ifges(4,4)) * qJD(1);
t327 = (Ifges(4,6) * qJD(3)) + (t349 * Ifges(4,4) - Ifges(4,2) * t346) * qJD(1);
t322 = -qJDD(1) * pkin(1) + t355;
t321 = t352 * pkin(1) + t369;
t305 = Ifges(5,1) * t331 + Ifges(5,4) * t330 + Ifges(5,5) * t367;
t304 = Ifges(5,4) * t331 + Ifges(5,2) * t330 + Ifges(5,6) * t367;
t303 = Ifges(5,5) * t331 + Ifges(5,6) * t330 + Ifges(5,3) * t367;
t290 = Ifges(6,5) * t308 + Ifges(6,6) * t307 + Ifges(6,3) * t339;
t272 = mrSges(6,2) * t286 - mrSges(6,3) * t280 + Ifges(6,1) * t289 + Ifges(6,4) * t288 + Ifges(6,5) * t332 + t307 * t290 - t339 * t291;
t271 = -mrSges(6,1) * t286 + mrSges(6,3) * t281 + Ifges(6,4) * t289 + Ifges(6,2) * t288 + Ifges(6,6) * t332 - t308 * t290 + t339 * t292;
t262 = mrSges(5,2) * t299 - mrSges(5,3) * t284 + Ifges(5,1) * t318 + Ifges(5,4) * t317 + Ifges(5,5) * t335 - pkin(7) * t270 - t345 * t271 + t348 * t272 + t330 * t303 - t304 * t367;
t261 = m(3) * t322 + qJDD(1) * mrSges(3,2) - (t352 * mrSges(3,3)) + t357;
t260 = -mrSges(5,1) * t299 + mrSges(5,3) * t285 + Ifges(5,4) * t318 + Ifges(5,2) * t317 + Ifges(5,6) * t335 - pkin(4) * t354 + pkin(7) * t359 + t348 * t271 + t345 * t272 - t331 * t303 + t305 * t367;
t1 = [mrSges(2,1) * t361 - mrSges(2,2) * t358 + mrSges(3,2) * t322 - mrSges(3,3) * t321 + t349 * (mrSges(4,2) * t314 - mrSges(4,3) * t310 + Ifges(4,1) * t336 - Ifges(4,4) * t335 + Ifges(4,5) * qJDD(3) - qJ(4) * t264 - qJD(3) * t327 - t343 * t260 + t344 * t262) - t346 * (-mrSges(4,1) * t314 - mrSges(5,1) * t284 + mrSges(5,2) * t285 + mrSges(4,3) * t311 + Ifges(4,4) * t336 - Ifges(5,5) * t318 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t317 - pkin(3) * t264 - pkin(4) * t270 + qJD(3) * t328 - t331 * t304 + t330 * t305 - t353 + (-Ifges(4,2) - Ifges(5,3)) * t335) - pkin(6) * t357 - pkin(1) * t261 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-m(3) * t321 + m(4) * t314 + t335 * mrSges(4,1) + t352 * mrSges(3,2) + t336 * mrSges(4,2) + t264 + qJDD(1) * mrSges(3,3) + (t337 * t346 + t338 * t349) * qJD(1)) * qJ(2); t261; Ifges(4,5) * t336 - Ifges(4,6) * t335 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t310 - mrSges(4,2) * t311 + t343 * t262 + t344 * t260 - pkin(3) * t279 + qJ(4) * t360 + (t349 * t327 + t346 * t328) * qJD(1); t279; t353;];
tauJ = t1;
