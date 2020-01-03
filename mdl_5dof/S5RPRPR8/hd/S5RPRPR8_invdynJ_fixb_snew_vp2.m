% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPR8
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:34
% EndTime: 2019-12-31 18:21:35
% DurationCPUTime: 1.03s
% Computational Cost: add. (8037->217), mult. (16437->280), div. (0->0), fcn. (10180->10), ass. (0->90)
t344 = sin(qJ(1));
t347 = cos(qJ(1));
t357 = t344 * g(1) - t347 * g(2);
t325 = qJDD(1) * pkin(1) + t357;
t349 = qJD(1) ^ 2;
t353 = -t347 * g(1) - t344 * g(2);
t328 = -t349 * pkin(1) + t353;
t339 = sin(pkin(8));
t341 = cos(pkin(8));
t301 = t341 * t325 - t339 * t328;
t293 = -qJDD(1) * pkin(2) - t349 * pkin(6) - t301;
t343 = sin(qJ(3));
t346 = cos(qJ(3));
t359 = qJD(1) * qJD(3);
t358 = t346 * t359;
t329 = t343 * qJDD(1) + t358;
t335 = t343 * t359;
t330 = t346 * qJDD(1) - t335;
t282 = (-t329 - t358) * qJ(4) + (-t330 + t335) * pkin(3) + t293;
t302 = t339 * t325 + t341 * t328;
t294 = -t349 * pkin(2) + qJDD(1) * pkin(6) + t302;
t337 = -g(3) + qJDD(2);
t288 = t346 * t294 + t343 * t337;
t326 = (-t346 * pkin(3) - t343 * qJ(4)) * qJD(1);
t348 = qJD(3) ^ 2;
t360 = t346 * qJD(1);
t286 = -t348 * pkin(3) + qJDD(3) * qJ(4) + t326 * t360 + t288;
t338 = sin(pkin(9));
t340 = cos(pkin(9));
t361 = t343 * qJD(1);
t322 = t338 * qJD(3) + t340 * t361;
t270 = -0.2e1 * qJD(4) * t322 + t340 * t282 - t338 * t286;
t308 = t338 * qJDD(3) + t340 * t329;
t321 = t340 * qJD(3) - t338 * t361;
t268 = (-t321 * t360 - t308) * pkin(7) + (t321 * t322 - t330) * pkin(4) + t270;
t271 = 0.2e1 * qJD(4) * t321 + t338 * t282 + t340 * t286;
t307 = t340 * qJDD(3) - t338 * t329;
t309 = -pkin(4) * t360 - t322 * pkin(7);
t320 = t321 ^ 2;
t269 = -t320 * pkin(4) + t307 * pkin(7) + t309 * t360 + t271;
t342 = sin(qJ(5));
t345 = cos(qJ(5));
t266 = t345 * t268 - t342 * t269;
t299 = t345 * t321 - t342 * t322;
t275 = t299 * qJD(5) + t342 * t307 + t345 * t308;
t300 = t342 * t321 + t345 * t322;
t285 = -t299 * mrSges(6,1) + t300 * mrSges(6,2);
t333 = qJD(5) - t360;
t289 = -t333 * mrSges(6,2) + t299 * mrSges(6,3);
t324 = qJDD(5) - t330;
t263 = m(6) * t266 + t324 * mrSges(6,1) - t275 * mrSges(6,3) - t300 * t285 + t333 * t289;
t267 = t342 * t268 + t345 * t269;
t274 = -t300 * qJD(5) + t345 * t307 - t342 * t308;
t290 = t333 * mrSges(6,1) - t300 * mrSges(6,3);
t264 = m(6) * t267 - t324 * mrSges(6,2) + t274 * mrSges(6,3) + t299 * t285 - t333 * t290;
t256 = t345 * t263 + t342 * t264;
t327 = (-t346 * mrSges(4,1) + t343 * mrSges(4,2)) * qJD(1);
t331 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t361;
t303 = -t321 * mrSges(5,1) + t322 * mrSges(5,2);
t305 = mrSges(5,2) * t360 + t321 * mrSges(5,3);
t254 = m(5) * t270 - t330 * mrSges(5,1) - t308 * mrSges(5,3) - t322 * t303 - t305 * t360 + t256;
t306 = -mrSges(5,1) * t360 - t322 * mrSges(5,3);
t354 = -t342 * t263 + t345 * t264;
t255 = m(5) * t271 + t330 * mrSges(5,2) + t307 * mrSges(5,3) + t321 * t303 + t306 * t360 + t354;
t355 = -t338 * t254 + t340 * t255;
t251 = m(4) * t288 - qJDD(3) * mrSges(4,2) + t330 * mrSges(4,3) - qJD(3) * t331 + t327 * t360 + t355;
t287 = -t343 * t294 + t346 * t337;
t284 = -qJDD(3) * pkin(3) - t348 * qJ(4) + t326 * t361 + qJDD(4) - t287;
t272 = -t307 * pkin(4) - t320 * pkin(7) + t322 * t309 + t284;
t352 = m(6) * t272 - t274 * mrSges(6,1) + t275 * mrSges(6,2) - t299 * t289 + t300 * t290;
t265 = m(5) * t284 - t307 * mrSges(5,1) + t308 * mrSges(5,2) - t321 * t305 + t322 * t306 + t352;
t332 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t360;
t259 = m(4) * t287 + qJDD(3) * mrSges(4,1) - t329 * mrSges(4,3) + qJD(3) * t332 - t327 * t361 - t265;
t356 = t346 * t251 - t343 * t259;
t252 = t340 * t254 + t338 * t255;
t351 = -m(4) * t293 + t330 * mrSges(4,1) - t329 * mrSges(4,2) - t331 * t361 + t332 * t360 - t252;
t279 = Ifges(6,4) * t300 + Ifges(6,2) * t299 + Ifges(6,6) * t333;
t280 = Ifges(6,1) * t300 + Ifges(6,4) * t299 + Ifges(6,5) * t333;
t350 = mrSges(6,1) * t266 - mrSges(6,2) * t267 + Ifges(6,5) * t275 + Ifges(6,6) * t274 + Ifges(6,3) * t324 + t300 * t279 - t299 * t280;
t318 = Ifges(4,5) * qJD(3) + (t343 * Ifges(4,1) + t346 * Ifges(4,4)) * qJD(1);
t317 = Ifges(4,6) * qJD(3) + (t343 * Ifges(4,4) + Ifges(4,2) * t346) * qJD(1);
t297 = Ifges(5,1) * t322 + Ifges(5,4) * t321 - Ifges(5,5) * t360;
t296 = Ifges(5,4) * t322 + Ifges(5,2) * t321 - Ifges(5,6) * t360;
t295 = Ifges(5,5) * t322 + Ifges(5,6) * t321 - Ifges(5,3) * t360;
t278 = Ifges(6,5) * t300 + Ifges(6,6) * t299 + Ifges(6,3) * t333;
t258 = mrSges(6,2) * t272 - mrSges(6,3) * t266 + Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t324 + t299 * t278 - t333 * t279;
t257 = -mrSges(6,1) * t272 + mrSges(6,3) * t267 + Ifges(6,4) * t275 + Ifges(6,2) * t274 + Ifges(6,6) * t324 - t300 * t278 + t333 * t280;
t249 = mrSges(5,2) * t284 - mrSges(5,3) * t270 + Ifges(5,1) * t308 + Ifges(5,4) * t307 - Ifges(5,5) * t330 - pkin(7) * t256 - t342 * t257 + t345 * t258 + t321 * t295 + t296 * t360;
t248 = -mrSges(5,1) * t284 + mrSges(5,3) * t271 + Ifges(5,4) * t308 + Ifges(5,2) * t307 - Ifges(5,6) * t330 - pkin(4) * t352 + pkin(7) * t354 + t345 * t257 + t342 * t258 - t322 * t295 - t297 * t360;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t357 - mrSges(2,2) * t353 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t301 - mrSges(3,2) * t302 + t343 * (mrSges(4,2) * t293 - mrSges(4,3) * t287 + Ifges(4,1) * t329 + Ifges(4,4) * t330 + Ifges(4,5) * qJDD(3) - qJ(4) * t252 - qJD(3) * t317 - t338 * t248 + t340 * t249) + t346 * (-mrSges(4,1) * t293 - mrSges(5,1) * t270 + mrSges(5,2) * t271 + mrSges(4,3) * t288 + Ifges(4,4) * t329 - Ifges(5,5) * t308 + Ifges(4,6) * qJDD(3) - Ifges(5,6) * t307 - pkin(3) * t252 - pkin(4) * t256 + qJD(3) * t318 - t322 * t296 + t321 * t297 - t350 + (Ifges(4,2) + Ifges(5,3)) * t330) + pkin(2) * t351 + pkin(6) * t356 + pkin(1) * (t339 * (m(3) * t302 - t349 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t356) + t341 * (m(3) * t301 + qJDD(1) * mrSges(3,1) - t349 * mrSges(3,2) + t351)); m(3) * t337 + t343 * t251 + t346 * t259; Ifges(4,5) * t329 + Ifges(4,6) * t330 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t287 - mrSges(4,2) * t288 + t338 * t249 + t340 * t248 - pkin(3) * t265 + qJ(4) * t355 + (t343 * t317 - t346 * t318) * qJD(1); t265; t350;];
tauJ = t1;
