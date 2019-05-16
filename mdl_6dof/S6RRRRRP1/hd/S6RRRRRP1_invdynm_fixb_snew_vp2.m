% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 04:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:22:04
% EndTime: 2019-05-08 04:23:02
% DurationCPUTime: 22.76s
% Computational Cost: add. (396488->382), mult. (857893->468), div. (0->0), fcn. (639663->10), ass. (0->145)
t330 = sin(qJ(3));
t331 = sin(qJ(2));
t335 = cos(qJ(3));
t336 = cos(qJ(2));
t303 = (t330 * t336 + t331 * t335) * qJD(1);
t360 = qJD(1) * qJD(2);
t310 = qJDD(1) * t331 + t336 * t360;
t311 = qJDD(1) * t336 - t331 * t360;
t278 = -qJD(3) * t303 - t310 * t330 + t311 * t335;
t302 = (-t330 * t331 + t335 * t336) * qJD(1);
t279 = qJD(3) * t302 + t310 * t335 + t311 * t330;
t329 = sin(qJ(4));
t334 = cos(qJ(4));
t289 = t302 * t334 - t303 * t329;
t245 = qJD(4) * t289 + t278 * t329 + t279 * t334;
t290 = t302 * t329 + t303 * t334;
t325 = qJD(2) + qJD(3);
t322 = qJD(4) + t325;
t328 = sin(qJ(5));
t333 = cos(qJ(5));
t273 = -t290 * t328 + t322 * t333;
t324 = qJDD(2) + qJDD(3);
t321 = qJDD(4) + t324;
t220 = qJD(5) * t273 + t245 * t333 + t321 * t328;
t274 = t290 * t333 + t322 * t328;
t252 = -mrSges(7,1) * t273 + mrSges(7,2) * t274;
t332 = sin(qJ(1));
t337 = cos(qJ(1));
t317 = -g(1) * t337 - g(2) * t332;
t338 = qJD(1) ^ 2;
t305 = -pkin(1) * t338 + qJDD(1) * pkin(7) + t317;
t365 = t331 * t305;
t367 = pkin(2) * t338;
t270 = qJDD(2) * pkin(2) - t310 * pkin(8) - t365 + (pkin(8) * t360 + t331 * t367 - g(3)) * t336;
t293 = -g(3) * t331 + t336 * t305;
t362 = qJD(1) * t331;
t315 = qJD(2) * pkin(2) - pkin(8) * t362;
t327 = t336 ^ 2;
t271 = pkin(8) * t311 - qJD(2) * t315 - t327 * t367 + t293;
t250 = t335 * t270 - t330 * t271;
t206 = (t302 * t325 - t279) * pkin(9) + (t302 * t303 + t324) * pkin(3) + t250;
t251 = t330 * t270 + t335 * t271;
t296 = pkin(3) * t325 - pkin(9) * t303;
t298 = t302 ^ 2;
t218 = -pkin(3) * t298 + pkin(9) * t278 - t296 * t325 + t251;
t204 = t329 * t206 + t334 * t218;
t265 = -pkin(4) * t289 - pkin(10) * t290;
t320 = t322 ^ 2;
t198 = -pkin(4) * t320 + pkin(10) * t321 + t265 * t289 + t204;
t316 = t332 * g(1) - t337 * g(2);
t351 = -qJDD(1) * pkin(1) - t316;
t280 = -t311 * pkin(2) + t315 * t362 + (-pkin(8) * t327 - pkin(7)) * t338 + t351;
t226 = -t278 * pkin(3) - t298 * pkin(9) + t303 * t296 + t280;
t244 = -qJD(4) * t290 + t278 * t334 - t279 * t329;
t201 = (-t289 * t322 - t245) * pkin(10) + (t290 * t322 - t244) * pkin(4) + t226;
t192 = -t328 * t198 + t333 * t201;
t241 = qJDD(5) - t244;
t286 = qJD(5) - t289;
t188 = -0.2e1 * qJD(6) * t274 + (t273 * t286 - t220) * qJ(6) + (t273 * t274 + t241) * pkin(5) + t192;
t254 = -mrSges(7,2) * t286 + mrSges(7,3) * t273;
t359 = m(7) * t188 + t241 * mrSges(7,1) + t286 * t254;
t185 = -t220 * mrSges(7,3) - t274 * t252 + t359;
t193 = t333 * t198 + t328 * t201;
t219 = -qJD(5) * t274 - t245 * t328 + t321 * t333;
t234 = Ifges(6,4) * t274 + Ifges(6,2) * t273 + Ifges(6,6) * t286;
t235 = Ifges(7,1) * t274 + Ifges(7,4) * t273 + Ifges(7,5) * t286;
t236 = Ifges(6,1) * t274 + Ifges(6,4) * t273 + Ifges(6,5) * t286;
t256 = pkin(5) * t286 - qJ(6) * t274;
t272 = t273 ^ 2;
t191 = -pkin(5) * t272 + qJ(6) * t219 + 0.2e1 * qJD(6) * t273 - t256 * t286 + t193;
t233 = Ifges(7,4) * t274 + Ifges(7,2) * t273 + Ifges(7,6) * t286;
t348 = -mrSges(7,1) * t188 + mrSges(7,2) * t191 - Ifges(7,5) * t220 - Ifges(7,6) * t219 - Ifges(7,3) * t241 - t274 * t233;
t369 = mrSges(6,1) * t192 - mrSges(6,2) * t193 + Ifges(6,5) * t220 + Ifges(6,6) * t219 + Ifges(6,3) * t241 + pkin(5) * t185 + t274 * t234 - (t236 + t235) * t273 - t348;
t264 = -mrSges(5,1) * t289 + mrSges(5,2) * t290;
t282 = mrSges(5,1) * t322 - mrSges(5,3) * t290;
t253 = -mrSges(6,1) * t273 + mrSges(6,2) * t274;
t255 = -mrSges(6,2) * t286 + mrSges(6,3) * t273;
t177 = m(6) * t192 + t241 * mrSges(6,1) + t286 * t255 + (-t252 - t253) * t274 + (-mrSges(6,3) - mrSges(7,3)) * t220 + t359;
t358 = m(7) * t191 + t219 * mrSges(7,3) + t273 * t252;
t257 = mrSges(7,1) * t286 - mrSges(7,3) * t274;
t363 = -mrSges(6,1) * t286 + mrSges(6,3) * t274 - t257;
t366 = -mrSges(6,2) - mrSges(7,2);
t180 = m(6) * t193 + t219 * mrSges(6,3) + t366 * t241 + t273 * t253 + t363 * t286 + t358;
t354 = -t177 * t328 + t333 * t180;
t170 = m(5) * t204 - mrSges(5,2) * t321 + mrSges(5,3) * t244 + t264 * t289 - t282 * t322 + t354;
t203 = t206 * t334 - t329 * t218;
t281 = -mrSges(5,2) * t322 + mrSges(5,3) * t289;
t197 = -pkin(4) * t321 - pkin(10) * t320 + t290 * t265 - t203;
t195 = -pkin(5) * t219 - qJ(6) * t272 + t256 * t274 + qJDD(6) + t197;
t353 = -m(7) * t195 + t219 * mrSges(7,1) + t273 * t254;
t343 = -m(6) * t197 + t219 * mrSges(6,1) + t366 * t220 + t273 * t255 + t363 * t274 + t353;
t182 = m(5) * t203 + t321 * mrSges(5,1) - t245 * mrSges(5,3) - t290 * t264 + t322 * t281 + t343;
t163 = t329 * t170 + t334 * t182;
t291 = -mrSges(4,1) * t302 + mrSges(4,2) * t303;
t294 = -mrSges(4,2) * t325 + mrSges(4,3) * t302;
t160 = m(4) * t250 + mrSges(4,1) * t324 - mrSges(4,3) * t279 - t291 * t303 + t294 * t325 + t163;
t295 = mrSges(4,1) * t325 - mrSges(4,3) * t303;
t355 = t334 * t170 - t182 * t329;
t161 = m(4) * t251 - mrSges(4,2) * t324 + mrSges(4,3) * t278 + t291 * t302 - t295 * t325 + t355;
t155 = t335 * t160 + t330 * t161;
t292 = -t336 * g(3) - t365;
t300 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t331 + Ifges(3,2) * t336) * qJD(1);
t301 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t331 + Ifges(3,4) * t336) * qJD(1);
t284 = Ifges(4,4) * t303 + Ifges(4,2) * t302 + Ifges(4,6) * t325;
t285 = Ifges(4,1) * t303 + Ifges(4,4) * t302 + Ifges(4,5) * t325;
t231 = Ifges(7,5) * t274 + Ifges(7,6) * t273 + Ifges(7,3) * t286;
t232 = Ifges(6,5) * t274 + Ifges(6,6) * t273 + Ifges(6,3) * t286;
t349 = -mrSges(7,1) * t195 + mrSges(7,3) * t191 + Ifges(7,4) * t220 + Ifges(7,2) * t219 + Ifges(7,6) * t241 + t286 * t235;
t165 = Ifges(6,4) * t220 + Ifges(6,2) * t219 + Ifges(6,6) * t241 + t286 * t236 - mrSges(6,1) * t197 + mrSges(6,3) * t193 - pkin(5) * (t220 * mrSges(7,2) - t353) + qJ(6) * (-t241 * mrSges(7,2) - t286 * t257 + t358) + (-pkin(5) * t257 - t231 - t232) * t274 + t349;
t347 = mrSges(7,2) * t195 - mrSges(7,3) * t188 + Ifges(7,1) * t220 + Ifges(7,4) * t219 + Ifges(7,5) * t241 + t273 * t231;
t172 = mrSges(6,2) * t197 - mrSges(6,3) * t192 + Ifges(6,1) * t220 + Ifges(6,4) * t219 + Ifges(6,5) * t241 - qJ(6) * t185 + t273 * t232 + (-t233 - t234) * t286 + t347;
t260 = Ifges(5,4) * t290 + Ifges(5,2) * t289 + Ifges(5,6) * t322;
t261 = Ifges(5,1) * t290 + Ifges(5,4) * t289 + Ifges(5,5) * t322;
t345 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t245 - Ifges(5,6) * t244 - Ifges(5,3) * t321 - pkin(4) * t343 - pkin(10) * t354 - t333 * t165 - t328 * t172 - t290 * t260 + t289 * t261;
t341 = -mrSges(4,1) * t250 + mrSges(4,2) * t251 - Ifges(4,5) * t279 - Ifges(4,6) * t278 - Ifges(4,3) * t324 - pkin(3) * t163 - t303 * t284 + t302 * t285 + t345;
t368 = mrSges(3,1) * t292 - mrSges(3,2) * t293 + Ifges(3,5) * t310 + Ifges(3,6) * t311 + Ifges(3,3) * qJDD(2) + pkin(2) * t155 + (t300 * t331 - t301 * t336) * qJD(1) - t341;
t174 = t333 * t177 + t328 * t180;
t361 = qJD(1) * t336;
t309 = (-mrSges(3,1) * t336 + mrSges(3,2) * t331) * qJD(1);
t314 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t361;
t153 = m(3) * t292 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t310 + qJD(2) * t314 - t309 * t362 + t155;
t313 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t362;
t356 = -t160 * t330 + t335 * t161;
t154 = m(3) * t293 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t311 - qJD(2) * t313 + t309 * t361 + t356;
t357 = -t153 * t331 + t336 * t154;
t350 = m(5) * t226 - t244 * mrSges(5,1) + t245 * mrSges(5,2) - t289 * t281 + t290 * t282 + t174;
t259 = Ifges(5,5) * t290 + Ifges(5,6) * t289 + Ifges(5,3) * t322;
t151 = mrSges(5,2) * t226 - mrSges(5,3) * t203 + Ifges(5,1) * t245 + Ifges(5,4) * t244 + Ifges(5,5) * t321 - pkin(10) * t174 - t165 * t328 + t172 * t333 + t259 * t289 - t260 * t322;
t156 = -mrSges(5,1) * t226 + mrSges(5,3) * t204 + Ifges(5,4) * t245 + Ifges(5,2) * t244 + Ifges(5,6) * t321 - pkin(4) * t174 - t290 * t259 + t322 * t261 - t369;
t283 = Ifges(4,5) * t303 + Ifges(4,6) * t302 + Ifges(4,3) * t325;
t146 = -mrSges(4,1) * t280 + mrSges(4,3) * t251 + Ifges(4,4) * t279 + Ifges(4,2) * t278 + Ifges(4,6) * t324 - pkin(3) * t350 + pkin(9) * t355 + t329 * t151 + t334 * t156 - t303 * t283 + t325 * t285;
t147 = mrSges(4,2) * t280 - mrSges(4,3) * t250 + Ifges(4,1) * t279 + Ifges(4,4) * t278 + Ifges(4,5) * t324 - pkin(9) * t163 + t151 * t334 - t156 * t329 + t283 * t302 - t284 * t325;
t299 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t331 + Ifges(3,6) * t336) * qJD(1);
t304 = -t338 * pkin(7) + t351;
t344 = m(4) * t280 - t278 * mrSges(4,1) + mrSges(4,2) * t279 - t302 * t294 + t295 * t303 + t350;
t142 = -mrSges(3,1) * t304 + mrSges(3,3) * t293 + Ifges(3,4) * t310 + Ifges(3,2) * t311 + Ifges(3,6) * qJDD(2) - pkin(2) * t344 + pkin(8) * t356 + qJD(2) * t301 + t335 * t146 + t330 * t147 - t299 * t362;
t144 = mrSges(3,2) * t304 - mrSges(3,3) * t292 + Ifges(3,1) * t310 + Ifges(3,4) * t311 + Ifges(3,5) * qJDD(2) - pkin(8) * t155 - qJD(2) * t300 - t146 * t330 + t147 * t335 + t299 * t361;
t340 = -m(3) * t304 + t311 * mrSges(3,1) - mrSges(3,2) * t310 - t313 * t362 + t314 * t361 - t344;
t346 = mrSges(2,1) * t316 - mrSges(2,2) * t317 + Ifges(2,3) * qJDD(1) + pkin(1) * t340 + pkin(7) * t357 + t336 * t142 + t331 * t144;
t166 = m(2) * t316 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t338 + t340;
t150 = t153 * t336 + t154 * t331;
t148 = m(2) * t317 - mrSges(2,1) * t338 - qJDD(1) * mrSges(2,2) + t357;
t145 = mrSges(2,1) * g(3) + mrSges(2,3) * t317 + t338 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t150 - t368;
t140 = -mrSges(2,2) * g(3) - mrSges(2,3) * t316 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t338 - pkin(7) * t150 - t142 * t331 + t144 * t336;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t140 - t332 * t145 - pkin(6) * (t148 * t332 + t166 * t337), t140, t144, t147, t151, t172, -t233 * t286 + t347; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t332 * t140 + t337 * t145 + pkin(6) * (t148 * t337 - t166 * t332), t145, t142, t146, t156, t165, -t274 * t231 + t349; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t346, t346, t368, -t341, -t345, t369, -t273 * t235 - t348;];
m_new  = t1;
