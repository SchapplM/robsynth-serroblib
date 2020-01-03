% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR7_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:09
% EndTime: 2019-12-31 18:19:10
% DurationCPUTime: 1.06s
% Computational Cost: add. (7028->219), mult. (14917->285), div. (0->0), fcn. (9288->10), ass. (0->92)
t349 = sin(pkin(9));
t351 = cos(pkin(9));
t354 = sin(qJ(3));
t357 = cos(qJ(3));
t325 = (t354 * t349 - t357 * t351) * qJD(1);
t376 = 2 * qJD(4);
t355 = sin(qJ(1));
t358 = cos(qJ(1));
t371 = t355 * g(1) - t358 * g(2);
t336 = qJDD(1) * pkin(1) + t371;
t360 = qJD(1) ^ 2;
t367 = -t358 * g(1) - t355 * g(2);
t338 = -t360 * pkin(1) + t367;
t350 = sin(pkin(8));
t352 = cos(pkin(8));
t316 = t350 * t336 + t352 * t338;
t309 = -t360 * pkin(2) + qJDD(1) * pkin(6) + t316;
t348 = -g(3) + qJDD(2);
t300 = -t354 * t309 + t357 * t348;
t373 = qJD(1) * qJD(3);
t372 = t357 * t373;
t339 = t354 * qJDD(1) + t372;
t289 = (-t339 + t372) * qJ(4) + (t354 * t357 * t360 + qJDD(3)) * pkin(3) + t300;
t301 = t357 * t309 + t354 * t348;
t340 = t357 * qJDD(1) - t354 * t373;
t374 = t354 * qJD(1);
t341 = qJD(3) * pkin(3) - qJ(4) * t374;
t347 = t357 ^ 2;
t290 = -t347 * t360 * pkin(3) + t340 * qJ(4) - qJD(3) * t341 + t301;
t285 = t349 * t289 + t351 * t290 - t325 * t376;
t326 = (t357 * t349 + t354 * t351) * qJD(1);
t311 = t325 * mrSges(5,1) + t326 * mrSges(5,2);
t317 = -t349 * t339 + t351 * t340;
t322 = qJD(3) * mrSges(5,1) - t326 * mrSges(5,3);
t312 = t325 * pkin(4) - t326 * pkin(7);
t359 = qJD(3) ^ 2;
t283 = -t359 * pkin(4) + qJDD(3) * pkin(7) - t325 * t312 + t285;
t315 = t352 * t336 - t350 * t338;
t364 = -qJDD(1) * pkin(2) - t315;
t291 = -t340 * pkin(3) + qJDD(4) + t341 * t374 + (-qJ(4) * t347 - pkin(6)) * t360 + t364;
t318 = t351 * t339 + t349 * t340;
t286 = (qJD(3) * t325 - t318) * pkin(7) + (qJD(3) * t326 - t317) * pkin(4) + t291;
t353 = sin(qJ(5));
t356 = cos(qJ(5));
t280 = -t353 * t283 + t356 * t286;
t319 = t356 * qJD(3) - t353 * t326;
t298 = t319 * qJD(5) + t353 * qJDD(3) + t356 * t318;
t320 = t353 * qJD(3) + t356 * t326;
t299 = -t319 * mrSges(6,1) + t320 * mrSges(6,2);
t324 = qJD(5) + t325;
t302 = -t324 * mrSges(6,2) + t319 * mrSges(6,3);
t314 = qJDD(5) - t317;
t278 = m(6) * t280 + t314 * mrSges(6,1) - t298 * mrSges(6,3) - t320 * t299 + t324 * t302;
t281 = t356 * t283 + t353 * t286;
t297 = -t320 * qJD(5) + t356 * qJDD(3) - t353 * t318;
t303 = t324 * mrSges(6,1) - t320 * mrSges(6,3);
t279 = m(6) * t281 - t314 * mrSges(6,2) + t297 * mrSges(6,3) + t319 * t299 - t324 * t303;
t368 = -t353 * t278 + t356 * t279;
t268 = m(5) * t285 - qJDD(3) * mrSges(5,2) + t317 * mrSges(5,3) - qJD(3) * t322 - t325 * t311 + t368;
t366 = -t351 * t289 + t349 * t290;
t284 = -0.2e1 * qJD(4) * t326 - t366;
t321 = -qJD(3) * mrSges(5,2) - t325 * mrSges(5,3);
t282 = -qJDD(3) * pkin(4) - t359 * pkin(7) + (t376 + t312) * t326 + t366;
t363 = -m(6) * t282 + t297 * mrSges(6,1) - t298 * mrSges(6,2) + t319 * t302 - t320 * t303;
t274 = m(5) * t284 + qJDD(3) * mrSges(5,1) - t318 * mrSges(5,3) + qJD(3) * t321 - t326 * t311 + t363;
t265 = t349 * t268 + t351 * t274;
t270 = t356 * t278 + t353 * t279;
t375 = qJD(1) * t357;
t337 = (-t357 * mrSges(4,1) + t354 * mrSges(4,2)) * qJD(1);
t343 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t375;
t263 = m(4) * t300 + qJDD(3) * mrSges(4,1) - t339 * mrSges(4,3) + qJD(3) * t343 - t337 * t374 + t265;
t342 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t374;
t369 = t351 * t268 - t349 * t274;
t264 = m(4) * t301 - qJDD(3) * mrSges(4,2) + t340 * mrSges(4,3) - qJD(3) * t342 + t337 * t375 + t369;
t370 = -t354 * t263 + t357 * t264;
t269 = m(5) * t291 - t317 * mrSges(5,1) + t318 * mrSges(5,2) + t325 * t321 + t326 * t322 + t270;
t293 = Ifges(6,4) * t320 + Ifges(6,2) * t319 + Ifges(6,6) * t324;
t294 = Ifges(6,1) * t320 + Ifges(6,4) * t319 + Ifges(6,5) * t324;
t362 = mrSges(6,1) * t280 - mrSges(6,2) * t281 + Ifges(6,5) * t298 + Ifges(6,6) * t297 + Ifges(6,3) * t314 + t320 * t293 - t319 * t294;
t308 = -t360 * pkin(6) + t364;
t361 = -m(4) * t308 + t340 * mrSges(4,1) - t339 * mrSges(4,2) - t342 * t374 + t343 * t375 - t269;
t332 = Ifges(4,5) * qJD(3) + (t354 * Ifges(4,1) + t357 * Ifges(4,4)) * qJD(1);
t331 = Ifges(4,6) * qJD(3) + (t354 * Ifges(4,4) + t357 * Ifges(4,2)) * qJD(1);
t307 = Ifges(5,1) * t326 - Ifges(5,4) * t325 + Ifges(5,5) * qJD(3);
t306 = Ifges(5,4) * t326 - Ifges(5,2) * t325 + Ifges(5,6) * qJD(3);
t305 = Ifges(5,5) * t326 - Ifges(5,6) * t325 + Ifges(5,3) * qJD(3);
t292 = Ifges(6,5) * t320 + Ifges(6,6) * t319 + Ifges(6,3) * t324;
t272 = mrSges(6,2) * t282 - mrSges(6,3) * t280 + Ifges(6,1) * t298 + Ifges(6,4) * t297 + Ifges(6,5) * t314 + t319 * t292 - t324 * t293;
t271 = -mrSges(6,1) * t282 + mrSges(6,3) * t281 + Ifges(6,4) * t298 + Ifges(6,2) * t297 + Ifges(6,6) * t314 - t320 * t292 + t324 * t294;
t261 = -mrSges(5,1) * t291 + mrSges(5,3) * t285 + Ifges(5,4) * t318 + Ifges(5,2) * t317 + Ifges(5,6) * qJDD(3) - pkin(4) * t270 + qJD(3) * t307 - t326 * t305 - t362;
t260 = mrSges(5,2) * t291 - mrSges(5,3) * t284 + Ifges(5,1) * t318 + Ifges(5,4) * t317 + Ifges(5,5) * qJDD(3) - pkin(7) * t270 - qJD(3) * t306 - t353 * t271 + t356 * t272 - t325 * t305;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t371 - mrSges(2,2) * t367 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t315 - mrSges(3,2) * t316 + t354 * (mrSges(4,2) * t308 - mrSges(4,3) * t300 + Ifges(4,1) * t339 + Ifges(4,4) * t340 + Ifges(4,5) * qJDD(3) - qJ(4) * t265 - qJD(3) * t331 + t351 * t260 - t349 * t261) + t357 * (-mrSges(4,1) * t308 + mrSges(4,3) * t301 + Ifges(4,4) * t339 + Ifges(4,2) * t340 + Ifges(4,6) * qJDD(3) - pkin(3) * t269 + qJ(4) * t369 + qJD(3) * t332 + t349 * t260 + t351 * t261) + pkin(2) * t361 + pkin(6) * t370 + pkin(1) * (t350 * (m(3) * t316 - t360 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t370) + t352 * (m(3) * t315 + qJDD(1) * mrSges(3,1) - t360 * mrSges(3,2) + t361)); m(3) * t348 + t357 * t263 + t354 * t264; Ifges(4,5) * t339 + Ifges(4,6) * t340 + mrSges(4,1) * t300 - mrSges(4,2) * t301 + Ifges(5,5) * t318 + Ifges(5,6) * t317 + t326 * t306 + t325 * t307 + mrSges(5,1) * t284 - mrSges(5,2) * t285 + t353 * t272 + t356 * t271 + pkin(4) * t363 + pkin(7) * t368 + pkin(3) * t265 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + (t354 * t331 - t357 * t332) * qJD(1); t269; t362;];
tauJ = t1;
