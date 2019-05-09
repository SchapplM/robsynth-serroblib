% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-05-06 17:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:16:38
% EndTime: 2019-05-06 17:17:26
% DurationCPUTime: 21.72s
% Computational Cost: add. (341668->381), mult. (790974->470), div. (0->0), fcn. (583235->10), ass. (0->143)
t328 = sin(pkin(10));
t329 = cos(pkin(10));
t332 = sin(qJ(2));
t336 = cos(qJ(2));
t301 = (-t328 * t332 + t329 * t336) * qJD(1);
t302 = (t328 * t336 + t329 * t332) * qJD(1);
t331 = sin(qJ(4));
t335 = cos(qJ(4));
t282 = t301 * t335 - t302 * t331;
t360 = qJD(1) * qJD(2);
t312 = qJDD(1) * t332 + t336 * t360;
t313 = qJDD(1) * t336 - t332 * t360;
t290 = -t312 * t328 + t313 * t329;
t291 = t312 * t329 + t313 * t328;
t251 = qJD(4) * t282 + t290 * t331 + t291 * t335;
t283 = t301 * t331 + t302 * t335;
t324 = qJD(2) + qJD(4);
t330 = sin(qJ(5));
t334 = cos(qJ(5));
t274 = -t283 * t330 + t324 * t334;
t323 = qJDD(2) + qJDD(4);
t223 = qJD(5) * t274 + t251 * t334 + t323 * t330;
t275 = t283 * t334 + t324 * t330;
t252 = -mrSges(7,1) * t274 + mrSges(7,2) * t275;
t333 = sin(qJ(1));
t337 = cos(qJ(1));
t319 = -g(1) * t337 - g(2) * t333;
t338 = qJD(1) ^ 2;
t307 = -pkin(1) * t338 + qJDD(1) * pkin(7) + t319;
t365 = t332 * t307;
t367 = pkin(2) * t338;
t270 = qJDD(2) * pkin(2) - t312 * qJ(3) - t365 + (qJ(3) * t360 + t332 * t367 - g(3)) * t336;
t293 = -g(3) * t332 + t307 * t336;
t362 = qJD(1) * t332;
t315 = qJD(2) * pkin(2) - qJ(3) * t362;
t327 = t336 ^ 2;
t271 = qJ(3) * t313 - qJD(2) * t315 - t327 * t367 + t293;
t237 = -0.2e1 * qJD(3) * t302 + t270 * t329 - t328 * t271;
t206 = (qJD(2) * t301 - t291) * pkin(8) + (t301 * t302 + qJDD(2)) * pkin(3) + t237;
t238 = 0.2e1 * qJD(3) * t301 + t270 * t328 + t271 * t329;
t296 = qJD(2) * pkin(3) - pkin(8) * t302;
t300 = t301 ^ 2;
t209 = -pkin(3) * t300 + pkin(8) * t290 - qJD(2) * t296 + t238;
t204 = t206 * t331 + t209 * t335;
t265 = -pkin(4) * t282 - pkin(9) * t283;
t322 = t324 ^ 2;
t198 = -pkin(4) * t322 + pkin(9) * t323 + t265 * t282 + t204;
t318 = t333 * g(1) - g(2) * t337;
t351 = -qJDD(1) * pkin(1) - t318;
t273 = -t313 * pkin(2) + qJDD(3) + t315 * t362 + (-qJ(3) * t327 - pkin(7)) * t338 + t351;
t234 = -t290 * pkin(3) - t300 * pkin(8) + t296 * t302 + t273;
t250 = -qJD(4) * t283 + t290 * t335 - t291 * t331;
t201 = (-t282 * t324 - t251) * pkin(9) + (t283 * t324 - t250) * pkin(4) + t234;
t192 = -t330 * t198 + t201 * t334;
t249 = qJDD(5) - t250;
t278 = qJD(5) - t282;
t188 = -0.2e1 * qJD(6) * t275 + (t274 * t278 - t223) * qJ(6) + (t274 * t275 + t249) * pkin(5) + t192;
t254 = -mrSges(7,2) * t278 + mrSges(7,3) * t274;
t359 = m(7) * t188 + mrSges(7,1) * t249 + t254 * t278;
t185 = -t223 * mrSges(7,3) - t275 * t252 + t359;
t193 = t198 * t334 + t201 * t330;
t222 = -qJD(5) * t275 - t251 * t330 + t323 * t334;
t228 = Ifges(6,4) * t275 + Ifges(6,2) * t274 + Ifges(6,6) * t278;
t229 = Ifges(7,1) * t275 + Ifges(7,4) * t274 + Ifges(7,5) * t278;
t230 = Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t278;
t256 = pkin(5) * t278 - qJ(6) * t275;
t272 = t274 ^ 2;
t191 = -pkin(5) * t272 + qJ(6) * t222 + 0.2e1 * qJD(6) * t274 - t256 * t278 + t193;
t227 = Ifges(7,4) * t275 + Ifges(7,2) * t274 + Ifges(7,6) * t278;
t348 = -mrSges(7,1) * t188 + mrSges(7,2) * t191 - Ifges(7,5) * t223 - Ifges(7,6) * t222 - Ifges(7,3) * t249 - t227 * t275;
t369 = mrSges(6,1) * t192 - mrSges(6,2) * t193 + Ifges(6,5) * t223 + Ifges(6,6) * t222 + Ifges(6,3) * t249 + pkin(5) * t185 + t275 * t228 - (t230 + t229) * t274 - t348;
t264 = -mrSges(5,1) * t282 + mrSges(5,2) * t283;
t277 = mrSges(5,1) * t324 - mrSges(5,3) * t283;
t253 = -mrSges(6,1) * t274 + mrSges(6,2) * t275;
t255 = -mrSges(6,2) * t278 + mrSges(6,3) * t274;
t177 = m(6) * t192 + t249 * mrSges(6,1) + t278 * t255 + (-t252 - t253) * t275 + (-mrSges(6,3) - mrSges(7,3)) * t223 + t359;
t358 = m(7) * t191 + mrSges(7,3) * t222 + t252 * t274;
t257 = mrSges(7,1) * t278 - mrSges(7,3) * t275;
t363 = -mrSges(6,1) * t278 + mrSges(6,3) * t275 - t257;
t366 = -mrSges(6,2) - mrSges(7,2);
t180 = m(6) * t193 + t222 * mrSges(6,3) + t249 * t366 + t274 * t253 + t278 * t363 + t358;
t354 = -t177 * t330 + t180 * t334;
t170 = m(5) * t204 - mrSges(5,2) * t323 + mrSges(5,3) * t250 + t264 * t282 - t277 * t324 + t354;
t203 = t206 * t335 - t209 * t331;
t276 = -mrSges(5,2) * t324 + mrSges(5,3) * t282;
t197 = -pkin(4) * t323 - pkin(9) * t322 + t265 * t283 - t203;
t195 = -pkin(5) * t222 - qJ(6) * t272 + t256 * t275 + qJDD(6) + t197;
t353 = -m(7) * t195 + mrSges(7,1) * t222 + t254 * t274;
t343 = -m(6) * t197 + mrSges(6,1) * t222 + t223 * t366 + t255 * t274 + t275 * t363 + t353;
t182 = m(5) * t203 + t323 * mrSges(5,1) - t251 * mrSges(5,3) - t283 * t264 + t324 * t276 + t343;
t163 = t170 * t331 + t182 * t335;
t286 = -mrSges(4,1) * t301 + mrSges(4,2) * t302;
t294 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t301;
t160 = m(4) * t237 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t291 + qJD(2) * t294 - t286 * t302 + t163;
t295 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t302;
t355 = t170 * t335 - t182 * t331;
t161 = m(4) * t238 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t290 - qJD(2) * t295 + t286 * t301 + t355;
t155 = t160 * t329 + t161 * t328;
t292 = -t336 * g(3) - t365;
t304 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t332 + Ifges(3,2) * t336) * qJD(1);
t305 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t332 + Ifges(3,4) * t336) * qJD(1);
t280 = Ifges(4,4) * t302 + Ifges(4,2) * t301 + Ifges(4,6) * qJD(2);
t281 = Ifges(4,1) * t302 + Ifges(4,4) * t301 + Ifges(4,5) * qJD(2);
t225 = Ifges(7,5) * t275 + Ifges(7,6) * t274 + Ifges(7,3) * t278;
t226 = Ifges(6,5) * t275 + Ifges(6,6) * t274 + Ifges(6,3) * t278;
t349 = -mrSges(7,1) * t195 + mrSges(7,3) * t191 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t249 + t229 * t278;
t165 = Ifges(6,4) * t223 + Ifges(6,2) * t222 + Ifges(6,6) * t249 + t278 * t230 - mrSges(6,1) * t197 + mrSges(6,3) * t193 - pkin(5) * (t223 * mrSges(7,2) - t353) + qJ(6) * (-t249 * mrSges(7,2) - t278 * t257 + t358) + (-pkin(5) * t257 - t225 - t226) * t275 + t349;
t347 = mrSges(7,2) * t195 - mrSges(7,3) * t188 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t249 + t225 * t274;
t172 = mrSges(6,2) * t197 - mrSges(6,3) * t192 + Ifges(6,1) * t223 + Ifges(6,4) * t222 + Ifges(6,5) * t249 - qJ(6) * t185 + t274 * t226 + (-t227 - t228) * t278 + t347;
t260 = Ifges(5,4) * t283 + Ifges(5,2) * t282 + Ifges(5,6) * t324;
t261 = Ifges(5,1) * t283 + Ifges(5,4) * t282 + Ifges(5,5) * t324;
t345 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t251 - Ifges(5,6) * t250 - Ifges(5,3) * t323 - pkin(4) * t343 - pkin(9) * t354 - t165 * t334 - t172 * t330 - t260 * t283 + t282 * t261;
t341 = -mrSges(4,1) * t237 + mrSges(4,2) * t238 - Ifges(4,5) * t291 - Ifges(4,6) * t290 - Ifges(4,3) * qJDD(2) - pkin(3) * t163 - t280 * t302 + t301 * t281 + t345;
t368 = mrSges(3,1) * t292 - mrSges(3,2) * t293 + Ifges(3,5) * t312 + Ifges(3,6) * t313 + Ifges(3,3) * qJDD(2) + pkin(2) * t155 + (t304 * t332 - t305 * t336) * qJD(1) - t341;
t174 = t177 * t334 + t180 * t330;
t361 = qJD(1) * t336;
t311 = (-mrSges(3,1) * t336 + mrSges(3,2) * t332) * qJD(1);
t317 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t361;
t153 = m(3) * t292 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t312 + qJD(2) * t317 - t311 * t362 + t155;
t316 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t362;
t356 = -t160 * t328 + t161 * t329;
t154 = m(3) * t293 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t313 - qJD(2) * t316 + t311 * t361 + t356;
t357 = -t153 * t332 + t154 * t336;
t350 = m(5) * t234 - mrSges(5,1) * t250 + mrSges(5,2) * t251 - t276 * t282 + t277 * t283 + t174;
t259 = Ifges(5,5) * t283 + Ifges(5,6) * t282 + Ifges(5,3) * t324;
t151 = mrSges(5,2) * t234 - mrSges(5,3) * t203 + Ifges(5,1) * t251 + Ifges(5,4) * t250 + Ifges(5,5) * t323 - pkin(9) * t174 - t165 * t330 + t172 * t334 + t259 * t282 - t260 * t324;
t156 = -mrSges(5,1) * t234 + mrSges(5,3) * t204 + Ifges(5,4) * t251 + Ifges(5,2) * t250 + Ifges(5,6) * t323 - pkin(4) * t174 - t283 * t259 + t324 * t261 - t369;
t279 = Ifges(4,5) * t302 + Ifges(4,6) * t301 + Ifges(4,3) * qJD(2);
t146 = -mrSges(4,1) * t273 + mrSges(4,3) * t238 + Ifges(4,4) * t291 + Ifges(4,2) * t290 + Ifges(4,6) * qJDD(2) - pkin(3) * t350 + pkin(8) * t355 + qJD(2) * t281 + t331 * t151 + t335 * t156 - t302 * t279;
t147 = mrSges(4,2) * t273 - mrSges(4,3) * t237 + Ifges(4,1) * t291 + Ifges(4,4) * t290 + Ifges(4,5) * qJDD(2) - pkin(8) * t163 - qJD(2) * t280 + t151 * t335 - t156 * t331 + t279 * t301;
t303 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t332 + Ifges(3,6) * t336) * qJD(1);
t306 = -t338 * pkin(7) + t351;
t344 = m(4) * t273 - mrSges(4,1) * t290 + mrSges(4,2) * t291 - t294 * t301 + t295 * t302 + t350;
t142 = -mrSges(3,1) * t306 + mrSges(3,3) * t293 + Ifges(3,4) * t312 + Ifges(3,2) * t313 + Ifges(3,6) * qJDD(2) - pkin(2) * t344 + qJ(3) * t356 + qJD(2) * t305 + t329 * t146 + t328 * t147 - t303 * t362;
t144 = mrSges(3,2) * t306 - mrSges(3,3) * t292 + Ifges(3,1) * t312 + Ifges(3,4) * t313 + Ifges(3,5) * qJDD(2) - qJ(3) * t155 - qJD(2) * t304 - t146 * t328 + t147 * t329 + t303 * t361;
t340 = -m(3) * t306 + mrSges(3,1) * t313 - mrSges(3,2) * t312 - t316 * t362 + t317 * t361 - t344;
t346 = mrSges(2,1) * t318 - mrSges(2,2) * t319 + Ifges(2,3) * qJDD(1) + pkin(1) * t340 + pkin(7) * t357 + t142 * t336 + t144 * t332;
t166 = m(2) * t318 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t338 + t340;
t150 = t153 * t336 + t154 * t332;
t148 = m(2) * t319 - mrSges(2,1) * t338 - qJDD(1) * mrSges(2,2) + t357;
t145 = mrSges(2,1) * g(3) + mrSges(2,3) * t319 + t338 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t150 - t368;
t140 = -mrSges(2,2) * g(3) - mrSges(2,3) * t318 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t338 - pkin(7) * t150 - t142 * t332 + t144 * t336;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t337 * t140 - t333 * t145 - pkin(6) * (t148 * t333 + t166 * t337), t140, t144, t147, t151, t172, -t227 * t278 + t347; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t333 * t140 + t337 * t145 + pkin(6) * (t148 * t337 - t166 * t333), t145, t142, t146, t156, t165, -t275 * t225 + t349; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t346, t346, t368, -t341, -t345, t369, -t274 * t229 - t348;];
m_new  = t1;
