% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 19:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:26:30
% EndTime: 2019-05-07 19:27:39
% DurationCPUTime: 53.30s
% Computational Cost: add. (935424->386), mult. (2167132->488), div. (0->0), fcn. (1657951->12), ass. (0->155)
t335 = sin(qJ(2));
t340 = cos(qJ(2));
t362 = qJD(1) * qJD(2);
t311 = qJDD(1) * t335 + t340 * t362;
t336 = sin(qJ(1));
t341 = cos(qJ(1));
t318 = -g(1) * t341 - g(2) * t336;
t342 = qJD(1) ^ 2;
t306 = -pkin(1) * t342 + qJDD(1) * pkin(7) + t318;
t365 = t335 * t306;
t366 = pkin(2) * t342;
t273 = qJDD(2) * pkin(2) - t311 * pkin(8) - t365 + (pkin(8) * t362 + t335 * t366 - g(3)) * t340;
t294 = -g(3) * t335 + t340 * t306;
t312 = qJDD(1) * t340 - t335 * t362;
t364 = qJD(1) * t335;
t316 = qJD(2) * pkin(2) - pkin(8) * t364;
t329 = t340 ^ 2;
t274 = pkin(8) * t312 - qJD(2) * t316 - t329 * t366 + t294;
t334 = sin(qJ(3));
t339 = cos(qJ(3));
t252 = t339 * t273 - t334 * t274;
t303 = (-t334 * t335 + t339 * t340) * qJD(1);
t279 = qJD(3) * t303 + t311 * t339 + t312 * t334;
t304 = (t334 * t340 + t335 * t339) * qJD(1);
t326 = qJDD(2) + qJDD(3);
t327 = qJD(2) + qJD(3);
t229 = (t303 * t327 - t279) * pkin(9) + (t303 * t304 + t326) * pkin(3) + t252;
t253 = t334 * t273 + t339 * t274;
t278 = -qJD(3) * t304 - t311 * t334 + t312 * t339;
t297 = pkin(3) * t327 - pkin(9) * t304;
t299 = t303 ^ 2;
t231 = -pkin(3) * t299 + pkin(9) * t278 - t297 * t327 + t253;
t333 = sin(qJ(4));
t338 = cos(qJ(4));
t208 = t338 * t229 - t333 * t231;
t290 = t303 * t338 - t304 * t333;
t249 = qJD(4) * t290 + t278 * t333 + t279 * t338;
t291 = t303 * t333 + t304 * t338;
t323 = qJDD(4) + t326;
t324 = qJD(4) + t327;
t204 = (t290 * t324 - t249) * qJ(5) + (t290 * t291 + t323) * pkin(4) + t208;
t209 = t333 * t229 + t338 * t231;
t248 = -qJD(4) * t291 + t278 * t338 - t279 * t333;
t282 = pkin(4) * t324 - qJ(5) * t291;
t289 = t290 ^ 2;
t206 = -pkin(4) * t289 + qJ(5) * t248 - t282 * t324 + t209;
t330 = sin(pkin(11));
t331 = cos(pkin(11));
t266 = t290 * t331 - t291 * t330;
t367 = 2 * qJD(5);
t201 = t330 * t204 + t331 * t206 + t266 * t367;
t220 = t248 * t331 - t249 * t330;
t267 = t290 * t330 + t291 * t331;
t240 = -mrSges(6,1) * t266 + mrSges(6,2) * t267;
t257 = mrSges(6,1) * t324 - mrSges(6,3) * t267;
t241 = -pkin(5) * t266 - pkin(10) * t267;
t322 = t324 ^ 2;
t198 = -pkin(5) * t322 + pkin(10) * t323 + t241 * t266 + t201;
t317 = t336 * g(1) - t341 * g(2);
t354 = -qJDD(1) * pkin(1) - t317;
t280 = -t312 * pkin(2) + t316 * t364 + (-pkin(8) * t329 - pkin(7)) * t342 + t354;
t243 = -pkin(3) * t278 - pkin(9) * t299 + t304 * t297 + t280;
t211 = -pkin(4) * t248 - qJ(5) * t289 + t291 * t282 + qJDD(5) + t243;
t221 = t248 * t330 + t249 * t331;
t202 = (-t266 * t324 - t221) * pkin(10) + (t267 * t324 - t220) * pkin(5) + t211;
t332 = sin(qJ(6));
t337 = cos(qJ(6));
t195 = -t198 * t332 + t202 * t337;
t254 = -t267 * t332 + t324 * t337;
t214 = qJD(6) * t254 + t221 * t337 + t323 * t332;
t219 = qJDD(6) - t220;
t255 = t267 * t337 + t324 * t332;
t232 = -mrSges(7,1) * t254 + mrSges(7,2) * t255;
t265 = qJD(6) - t266;
t233 = -mrSges(7,2) * t265 + mrSges(7,3) * t254;
t191 = m(7) * t195 + mrSges(7,1) * t219 - mrSges(7,3) * t214 - t232 * t255 + t233 * t265;
t196 = t198 * t337 + t202 * t332;
t213 = -qJD(6) * t255 - t221 * t332 + t323 * t337;
t234 = mrSges(7,1) * t265 - mrSges(7,3) * t255;
t192 = m(7) * t196 - mrSges(7,2) * t219 + mrSges(7,3) * t213 + t232 * t254 - t234 * t265;
t357 = -t191 * t332 + t337 * t192;
t178 = m(6) * t201 - mrSges(6,2) * t323 + mrSges(6,3) * t220 + t240 * t266 - t257 * t324 + t357;
t356 = -t331 * t204 + t330 * t206;
t200 = -0.2e1 * qJD(5) * t267 - t356;
t256 = -mrSges(6,2) * t324 + mrSges(6,3) * t266;
t197 = -t323 * pkin(5) - t322 * pkin(10) + (t367 + t241) * t267 + t356;
t351 = -m(7) * t197 + t213 * mrSges(7,1) - mrSges(7,2) * t214 + t254 * t233 - t234 * t255;
t187 = m(6) * t200 + mrSges(6,1) * t323 - mrSges(6,3) * t221 - t240 * t267 + t256 * t324 + t351;
t173 = t330 * t178 + t331 * t187;
t268 = -mrSges(5,1) * t290 + mrSges(5,2) * t291;
t281 = -mrSges(5,2) * t324 + mrSges(5,3) * t290;
t170 = m(5) * t208 + mrSges(5,1) * t323 - mrSges(5,3) * t249 - t268 * t291 + t281 * t324 + t173;
t283 = mrSges(5,1) * t324 - mrSges(5,3) * t291;
t358 = t331 * t178 - t187 * t330;
t171 = m(5) * t209 - mrSges(5,2) * t323 + mrSges(5,3) * t248 + t268 * t290 - t283 * t324 + t358;
t164 = t338 * t170 + t333 * t171;
t292 = -mrSges(4,1) * t303 + mrSges(4,2) * t304;
t295 = -mrSges(4,2) * t327 + mrSges(4,3) * t303;
t161 = m(4) * t252 + mrSges(4,1) * t326 - mrSges(4,3) * t279 - t292 * t304 + t295 * t327 + t164;
t296 = mrSges(4,1) * t327 - mrSges(4,3) * t304;
t359 = -t170 * t333 + t338 * t171;
t162 = m(4) * t253 - mrSges(4,2) * t326 + mrSges(4,3) * t278 + t292 * t303 - t296 * t327 + t359;
t156 = t339 * t161 + t334 * t162;
t293 = -t340 * g(3) - t365;
t301 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t335 + Ifges(3,2) * t340) * qJD(1);
t302 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t335 + Ifges(3,4) * t340) * qJD(1);
t285 = Ifges(4,4) * t304 + Ifges(4,2) * t303 + Ifges(4,6) * t327;
t286 = Ifges(4,1) * t304 + Ifges(4,4) * t303 + Ifges(4,5) * t327;
t259 = Ifges(5,4) * t291 + Ifges(5,2) * t290 + Ifges(5,6) * t324;
t260 = Ifges(5,1) * t291 + Ifges(5,4) * t290 + Ifges(5,5) * t324;
t222 = Ifges(7,5) * t255 + Ifges(7,6) * t254 + Ifges(7,3) * t265;
t224 = Ifges(7,1) * t255 + Ifges(7,4) * t254 + Ifges(7,5) * t265;
t184 = -mrSges(7,1) * t197 + mrSges(7,3) * t196 + Ifges(7,4) * t214 + Ifges(7,2) * t213 + Ifges(7,6) * t219 - t222 * t255 + t224 * t265;
t223 = Ifges(7,4) * t255 + Ifges(7,2) * t254 + Ifges(7,6) * t265;
t185 = mrSges(7,2) * t197 - mrSges(7,3) * t195 + Ifges(7,1) * t214 + Ifges(7,4) * t213 + Ifges(7,5) * t219 + t222 * t254 - t223 * t265;
t236 = Ifges(6,4) * t267 + Ifges(6,2) * t266 + Ifges(6,6) * t324;
t237 = Ifges(6,1) * t267 + Ifges(6,4) * t266 + Ifges(6,5) * t324;
t349 = -mrSges(6,1) * t200 + mrSges(6,2) * t201 - Ifges(6,5) * t221 - Ifges(6,6) * t220 - Ifges(6,3) * t323 - pkin(5) * t351 - pkin(10) * t357 - t337 * t184 - t332 * t185 - t267 * t236 + t266 * t237;
t346 = -mrSges(5,1) * t208 + mrSges(5,2) * t209 - Ifges(5,5) * t249 - Ifges(5,6) * t248 - Ifges(5,3) * t323 - pkin(4) * t173 - t291 * t259 + t290 * t260 + t349;
t344 = mrSges(4,1) * t252 - mrSges(4,2) * t253 + Ifges(4,5) * t279 + Ifges(4,6) * t278 + Ifges(4,3) * t326 + pkin(3) * t164 + t304 * t285 - t303 * t286 - t346;
t368 = mrSges(3,1) * t293 - mrSges(3,2) * t294 + Ifges(3,5) * t311 + Ifges(3,6) * t312 + Ifges(3,3) * qJDD(2) + pkin(2) * t156 + (t301 * t335 - t302 * t340) * qJD(1) + t344;
t180 = t337 * t191 + t332 * t192;
t363 = qJD(1) * t340;
t310 = (-mrSges(3,1) * t340 + mrSges(3,2) * t335) * qJD(1);
t315 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t363;
t154 = m(3) * t293 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t311 + qJD(2) * t315 - t310 * t364 + t156;
t314 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t364;
t360 = -t161 * t334 + t339 * t162;
t155 = m(3) * t294 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t312 - qJD(2) * t314 + t310 * t363 + t360;
t361 = -t154 * t335 + t340 * t155;
t353 = m(6) * t211 - t220 * mrSges(6,1) + t221 * mrSges(6,2) - t266 * t256 + t267 * t257 + t180;
t235 = Ifges(6,5) * t267 + Ifges(6,6) * t266 + Ifges(6,3) * t324;
t165 = mrSges(6,2) * t211 - mrSges(6,3) * t200 + Ifges(6,1) * t221 + Ifges(6,4) * t220 + Ifges(6,5) * t323 - pkin(10) * t180 - t184 * t332 + t185 * t337 + t235 * t266 - t236 * t324;
t348 = mrSges(7,1) * t195 - mrSges(7,2) * t196 + Ifges(7,5) * t214 + Ifges(7,6) * t213 + Ifges(7,3) * t219 + t223 * t255 - t224 * t254;
t166 = -mrSges(6,1) * t211 + mrSges(6,3) * t201 + Ifges(6,4) * t221 + Ifges(6,2) * t220 + Ifges(6,6) * t323 - pkin(5) * t180 - t235 * t267 + t237 * t324 - t348;
t258 = Ifges(5,5) * t291 + Ifges(5,6) * t290 + Ifges(5,3) * t324;
t152 = -mrSges(5,1) * t243 + mrSges(5,3) * t209 + Ifges(5,4) * t249 + Ifges(5,2) * t248 + Ifges(5,6) * t323 - pkin(4) * t353 + qJ(5) * t358 + t330 * t165 + t331 * t166 - t291 * t258 + t324 * t260;
t157 = mrSges(5,2) * t243 - mrSges(5,3) * t208 + Ifges(5,1) * t249 + Ifges(5,4) * t248 + Ifges(5,5) * t323 - qJ(5) * t173 + t165 * t331 - t166 * t330 + t258 * t290 - t259 * t324;
t284 = Ifges(4,5) * t304 + Ifges(4,6) * t303 + Ifges(4,3) * t327;
t350 = m(5) * t243 - t248 * mrSges(5,1) + t249 * mrSges(5,2) - t290 * t281 + t291 * t283 + t353;
t147 = -mrSges(4,1) * t280 + mrSges(4,3) * t253 + Ifges(4,4) * t279 + Ifges(4,2) * t278 + Ifges(4,6) * t326 - pkin(3) * t350 + pkin(9) * t359 + t338 * t152 + t333 * t157 - t304 * t284 + t327 * t286;
t148 = mrSges(4,2) * t280 - mrSges(4,3) * t252 + Ifges(4,1) * t279 + Ifges(4,4) * t278 + Ifges(4,5) * t326 - pkin(9) * t164 - t152 * t333 + t157 * t338 + t284 * t303 - t285 * t327;
t300 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t335 + Ifges(3,6) * t340) * qJD(1);
t305 = -pkin(7) * t342 + t354;
t347 = m(4) * t280 - t278 * mrSges(4,1) + t279 * mrSges(4,2) - t303 * t295 + t304 * t296 + t350;
t143 = -mrSges(3,1) * t305 + mrSges(3,3) * t294 + Ifges(3,4) * t311 + Ifges(3,2) * t312 + Ifges(3,6) * qJDD(2) - pkin(2) * t347 + pkin(8) * t360 + qJD(2) * t302 + t339 * t147 + t334 * t148 - t300 * t364;
t145 = mrSges(3,2) * t305 - mrSges(3,3) * t293 + Ifges(3,1) * t311 + Ifges(3,4) * t312 + Ifges(3,5) * qJDD(2) - pkin(8) * t156 - qJD(2) * t301 - t147 * t334 + t148 * t339 + t300 * t363;
t345 = -m(3) * t305 + t312 * mrSges(3,1) - t311 * mrSges(3,2) - t314 * t364 + t315 * t363 - t347;
t352 = mrSges(2,1) * t317 - mrSges(2,2) * t318 + Ifges(2,3) * qJDD(1) + pkin(1) * t345 + pkin(7) * t361 + t340 * t143 + t335 * t145;
t174 = m(2) * t317 + qJDD(1) * mrSges(2,1) - t342 * mrSges(2,2) + t345;
t151 = t154 * t340 + t155 * t335;
t149 = m(2) * t318 - mrSges(2,1) * t342 - qJDD(1) * mrSges(2,2) + t361;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t318 + t342 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t368;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t317 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t342 - pkin(7) * t151 - t143 * t335 + t145 * t340;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t341 * t141 - t336 * t146 - pkin(6) * (t149 * t336 + t174 * t341), t141, t145, t148, t157, t165, t185; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t336 * t141 + t341 * t146 + pkin(6) * (t149 * t341 - t174 * t336), t146, t143, t147, t152, t166, t184; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t352, t352, t368, t344, -t346, -t349, t348;];
m_new  = t1;
