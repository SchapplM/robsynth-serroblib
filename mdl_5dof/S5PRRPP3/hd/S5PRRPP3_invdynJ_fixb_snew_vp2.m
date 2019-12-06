% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:12:12
% EndTime: 2019-12-05 16:12:14
% DurationCPUTime: 1.01s
% Computational Cost: add. (2675->191), mult. (5463->231), div. (0->0), fcn. (3211->8), ass. (0->83)
t362 = Ifges(5,1) + Ifges(6,1);
t350 = Ifges(5,4) - Ifges(6,5);
t361 = Ifges(5,5) + Ifges(6,4);
t360 = Ifges(5,2) + Ifges(6,3);
t358 = Ifges(5,6) - Ifges(6,6);
t359 = Ifges(6,2) + Ifges(5,3);
t323 = sin(pkin(8));
t326 = sin(qJ(3));
t342 = qJD(2) * t326;
t348 = cos(pkin(8));
t301 = -t348 * qJD(3) + t323 * t342;
t302 = t323 * qJD(3) + t348 * t342;
t328 = cos(qJ(3));
t341 = qJD(2) * t328;
t355 = t360 * t301 - t350 * t302 + t358 * t341;
t354 = t350 * t301 - t362 * t302 + t361 * t341;
t353 = -2 * qJD(4);
t352 = -2 * qJD(5);
t351 = -mrSges(5,3) - mrSges(6,2);
t349 = t326 * Ifges(4,4);
t324 = sin(pkin(7));
t325 = cos(pkin(7));
t313 = -g(1) * t324 + g(2) * t325;
t347 = t313 * t328;
t331 = qJD(2) ^ 2;
t346 = t328 ^ 2 * t331;
t314 = -g(1) * t325 - g(2) * t324;
t322 = -g(3) + qJDD(1);
t327 = sin(qJ(2));
t329 = cos(qJ(2));
t293 = -t327 * t314 + t322 * t329;
t283 = -qJDD(2) * pkin(2) - pkin(6) * t331 - t293;
t340 = qJD(2) * qJD(3);
t339 = t328 * t340;
t311 = qJDD(2) * t326 + t339;
t338 = t326 * t340;
t312 = qJDD(2) * t328 - t338;
t266 = (-t311 - t339) * qJ(4) + (-t312 + t338) * pkin(3) + t283;
t294 = t329 * t314 + t327 * t322;
t284 = -pkin(2) * t331 + qJDD(2) * pkin(6) + t294;
t269 = t328 * t284 + t326 * t313;
t309 = (-pkin(3) * t328 - qJ(4) * t326) * qJD(2);
t330 = qJD(3) ^ 2;
t267 = -pkin(3) * t330 + qJDD(3) * qJ(4) + t309 * t341 + t269;
t262 = t323 * t266 + t348 * t267 + t301 * t353;
t278 = pkin(4) * t301 - qJ(5) * t302;
t258 = -pkin(4) * t346 - qJ(5) * t312 - t278 * t301 + t341 * t352 + t262;
t345 = m(6) * t258 - t312 * mrSges(6,3);
t344 = t358 * t301 - t302 * t361 + t359 * t341;
t279 = mrSges(6,1) * t301 - mrSges(6,3) * t302;
t343 = -mrSges(5,1) * t301 - mrSges(5,2) * t302 - t279;
t334 = t348 * t266 - t323 * t267;
t259 = -qJ(5) * t346 + t312 * pkin(4) + qJDD(5) + ((2 * qJD(4)) + t278) * t302 - t334;
t337 = -m(6) * t259 - t312 * mrSges(6,1);
t291 = -t348 * qJDD(3) + t311 * t323;
t292 = t323 * qJDD(3) + t348 * t311;
t281 = t326 * t284;
t333 = -qJDD(3) * pkin(3) - qJ(4) * t330 + t309 * t342 + qJDD(4) + t281;
t260 = pkin(4) * t291 - qJ(5) * t292 + t302 * t352 + (-t313 + (-pkin(4) * t302 - qJ(5) * t301) * qJD(2)) * t328 + t333;
t287 = -mrSges(6,2) * t301 - mrSges(6,3) * t341;
t290 = mrSges(6,1) * t341 + mrSges(6,2) * t302;
t256 = m(6) * t260 + mrSges(6,1) * t291 - t292 * mrSges(6,3) + t287 * t301 - t302 * t290;
t265 = t333 - t347;
t288 = mrSges(5,2) * t341 - mrSges(5,3) * t301;
t289 = -mrSges(5,1) * t341 - mrSges(5,3) * t302;
t254 = m(5) * t265 + t291 * mrSges(5,1) + mrSges(5,2) * t292 + t301 * t288 + t289 * t302 + t256;
t268 = -t281 + t347;
t310 = (-mrSges(4,1) * t328 + mrSges(4,2) * t326) * qJD(2);
t315 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t342;
t316 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t341;
t251 = m(5) * t262 + mrSges(5,2) * t312 + t343 * t301 + t351 * t291 + (t289 - t290) * t341 + t345;
t261 = t302 * t353 + t334;
t252 = m(5) * t261 - mrSges(5,1) * t312 + t343 * t302 + t351 * t292 + (-t287 - t288) * t341 + t337;
t335 = t348 * t251 - t252 * t323;
t336 = -(m(4) * t268 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t311 + qJD(3) * t316 - t310 * t342 - t254) * t326 + t328 * (m(4) * t269 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t312 - qJD(3) * t315 + t310 * t341 + t335);
t249 = t323 * t251 + t348 * t252;
t332 = -m(4) * t283 + t312 * mrSges(4,1) - t311 * mrSges(4,2) - t315 * t342 + t316 * t341 - t249;
t299 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t326 + Ifges(4,4) * t328) * qJD(2);
t298 = Ifges(4,6) * qJD(3) + (Ifges(4,2) * t328 + t349) * qJD(2);
t255 = mrSges(6,2) * t292 + t279 * t302 + t287 * t341 - t337;
t248 = mrSges(5,2) * t265 + mrSges(6,2) * t259 - mrSges(5,3) * t261 - mrSges(6,3) * t260 - qJ(5) * t256 - t350 * t291 + t362 * t292 + t344 * t301 - t312 * t361 - t355 * t341;
t247 = -mrSges(5,1) * t265 - mrSges(6,1) * t260 + mrSges(6,2) * t258 + mrSges(5,3) * t262 - pkin(4) * t256 - t360 * t291 + t350 * t292 + t344 * t302 - t312 * t358 + t354 * t341;
t1 = [m(2) * t322 + t327 * (m(3) * t294 - mrSges(3,1) * t331 - qJDD(2) * mrSges(3,2) + t336) + t329 * (m(3) * t293 + qJDD(2) * mrSges(3,1) - t331 * mrSges(3,2) + t332); Ifges(3,3) * qJDD(2) + mrSges(3,1) * t293 - mrSges(3,2) * t294 + t326 * (mrSges(4,2) * t283 - mrSges(4,3) * t268 + Ifges(4,1) * t311 + Ifges(4,5) * qJDD(3) - qJ(4) * t249 - qJD(3) * t298 - t323 * t247 + t348 * t248) + t328 * (Ifges(4,4) * t311 + Ifges(4,6) * qJDD(3) + qJD(3) * t299 - mrSges(4,1) * t283 + mrSges(4,3) * t269 - mrSges(5,1) * t261 + mrSges(5,2) * t262 + mrSges(6,1) * t259 - mrSges(6,3) * t258 + pkin(4) * t255 - qJ(5) * (-t290 * t341 + t345) - pkin(3) * t249 + t355 * t302 + (qJ(5) * t279 + t354) * t301 - t361 * t292 + (qJ(5) * mrSges(6,2) + t358) * t291) + pkin(2) * t332 + pkin(6) * t336 + (t349 + t328 * (Ifges(4,2) + t359)) * t312; Ifges(4,5) * t311 + Ifges(4,6) * t312 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t268 - mrSges(4,2) * t269 + t323 * t248 + t348 * t247 - pkin(3) * t254 + qJ(4) * t335 + (t298 * t326 - t299 * t328) * qJD(2); t254; t255;];
tauJ = t1;
