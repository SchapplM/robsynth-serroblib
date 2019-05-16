% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRP2
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
% Datum: 2019-05-08 04:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP2_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:29:33
% EndTime: 2019-05-08 04:30:29
% DurationCPUTime: 22.12s
% Computational Cost: add. (388135->383), mult. (840839->468), div. (0->0), fcn. (625711->10), ass. (0->144)
t329 = sin(qJ(2));
t333 = cos(qJ(2));
t356 = qJD(1) * qJD(2);
t306 = qJDD(1) * t329 + t333 * t356;
t330 = sin(qJ(1));
t334 = cos(qJ(1));
t313 = -g(1) * t334 - g(2) * t330;
t335 = qJD(1) ^ 2;
t301 = -pkin(1) * t335 + qJDD(1) * pkin(7) + t313;
t362 = t329 * t301;
t364 = pkin(2) * t335;
t266 = qJDD(2) * pkin(2) - t306 * pkin(8) - t362 + (pkin(8) * t356 + t329 * t364 - g(3)) * t333;
t289 = -g(3) * t329 + t333 * t301;
t307 = qJDD(1) * t333 - t329 * t356;
t358 = qJD(1) * t329;
t311 = qJD(2) * pkin(2) - pkin(8) * t358;
t325 = t333 ^ 2;
t267 = pkin(8) * t307 - qJD(2) * t311 - t325 * t364 + t289;
t328 = sin(qJ(3));
t332 = cos(qJ(3));
t246 = t332 * t266 - t328 * t267;
t298 = (-t328 * t329 + t332 * t333) * qJD(1);
t274 = qJD(3) * t298 + t306 * t332 + t307 * t328;
t299 = (t328 * t333 + t329 * t332) * qJD(1);
t322 = qJDD(2) + qJDD(3);
t323 = qJD(2) + qJD(3);
t206 = (t298 * t323 - t274) * pkin(9) + (t298 * t299 + t322) * pkin(3) + t246;
t247 = t328 * t266 + t332 * t267;
t273 = -qJD(3) * t299 - t306 * t328 + t307 * t332;
t292 = pkin(3) * t323 - pkin(9) * t299;
t294 = t298 ^ 2;
t215 = -pkin(3) * t294 + pkin(9) * t273 - t292 * t323 + t247;
t327 = sin(qJ(4));
t331 = cos(qJ(4));
t204 = t327 * t206 + t331 * t215;
t285 = t298 * t331 - t299 * t327;
t286 = t298 * t327 + t299 * t331;
t261 = -pkin(4) * t285 - pkin(10) * t286;
t320 = qJD(4) + t323;
t318 = t320 ^ 2;
t319 = qJDD(4) + t322;
t199 = -pkin(4) * t318 + pkin(10) * t319 + t261 * t285 + t204;
t312 = t330 * g(1) - t334 * g(2);
t347 = -qJDD(1) * pkin(1) - t312;
t275 = -t307 * pkin(2) + t311 * t358 + (-pkin(8) * t325 - pkin(7)) * t335 + t347;
t221 = -t273 * pkin(3) - t294 * pkin(9) + t299 * t292 + t275;
t239 = -qJD(4) * t286 + t273 * t331 - t274 * t327;
t240 = qJD(4) * t285 + t273 * t327 + t274 * t331;
t201 = (-t285 * t320 - t240) * pkin(10) + (t286 * t320 - t239) * pkin(4) + t221;
t326 = sin(qJ(5));
t365 = cos(qJ(5));
t195 = -t326 * t199 + t365 * t201;
t196 = t365 * t199 + t326 * t201;
t269 = t365 * t286 + t326 * t320;
t216 = t269 * qJD(5) + t326 * t240 - t365 * t319;
t268 = t326 * t286 - t365 * t320;
t217 = -t268 * qJD(5) + t365 * t240 + t326 * t319;
t282 = qJD(5) - t285;
t226 = Ifges(7,5) * t269 + Ifges(7,6) * t282 + Ifges(7,3) * t268;
t229 = Ifges(6,4) * t269 - Ifges(6,2) * t268 + Ifges(6,6) * t282;
t231 = Ifges(6,1) * t269 - Ifges(6,4) * t268 + Ifges(6,5) * t282;
t236 = qJDD(5) - t239;
t249 = mrSges(7,1) * t268 - mrSges(7,3) * t269;
t248 = pkin(5) * t268 - qJ(6) * t269;
t281 = t282 ^ 2;
t191 = -pkin(5) * t281 + qJ(6) * t236 + 0.2e1 * qJD(6) * t282 - t248 * t268 + t196;
t193 = -t236 * pkin(5) - t281 * qJ(6) + t269 * t248 + qJDD(6) - t195;
t230 = Ifges(7,1) * t269 + Ifges(7,4) * t282 + Ifges(7,5) * t268;
t345 = mrSges(7,1) * t193 - mrSges(7,3) * t191 - Ifges(7,4) * t217 - Ifges(7,2) * t236 - Ifges(7,6) * t216 - t268 * t230;
t251 = -mrSges(7,2) * t268 + mrSges(7,3) * t282;
t350 = -m(7) * t193 + t236 * mrSges(7,1) + t282 * t251;
t254 = -mrSges(7,1) * t282 + mrSges(7,2) * t269;
t355 = m(7) * t191 + t236 * mrSges(7,3) + t282 * t254;
t367 = -(-t229 + t226) * t269 + mrSges(6,1) * t195 - mrSges(6,2) * t196 + Ifges(6,5) * t217 - Ifges(6,6) * t216 + Ifges(6,3) * t236 + pkin(5) * (-t217 * mrSges(7,2) - t269 * t249 + t350) + qJ(6) * (-t216 * mrSges(7,2) - t268 * t249 + t355) + t268 * t231 - t345;
t260 = -mrSges(5,1) * t285 + mrSges(5,2) * t286;
t277 = mrSges(5,1) * t320 - mrSges(5,3) * t286;
t253 = mrSges(6,1) * t282 - mrSges(6,3) * t269;
t359 = -mrSges(6,1) * t268 - mrSges(6,2) * t269 - t249;
t363 = -mrSges(6,3) - mrSges(7,2);
t181 = m(6) * t196 - t236 * mrSges(6,2) + t363 * t216 - t282 * t253 + t359 * t268 + t355;
t252 = -mrSges(6,2) * t282 - mrSges(6,3) * t268;
t183 = m(6) * t195 + t236 * mrSges(6,1) + t363 * t217 + t282 * t252 + t359 * t269 + t350;
t351 = t365 * t181 - t183 * t326;
t169 = m(5) * t204 - mrSges(5,2) * t319 + mrSges(5,3) * t239 + t260 * t285 - t277 * t320 + t351;
t203 = t331 * t206 - t327 * t215;
t276 = -mrSges(5,2) * t320 + mrSges(5,3) * t285;
t198 = -t319 * pkin(4) - t318 * pkin(10) + t286 * t261 - t203;
t194 = -0.2e1 * qJD(6) * t269 + (t268 * t282 - t217) * qJ(6) + (t269 * t282 + t216) * pkin(5) + t198;
t188 = m(7) * t194 + mrSges(7,1) * t216 - t217 * mrSges(7,3) + t251 * t268 - t269 * t254;
t340 = -m(6) * t198 - t216 * mrSges(6,1) - mrSges(6,2) * t217 - t268 * t252 - t253 * t269 - t188;
t178 = m(5) * t203 + mrSges(5,1) * t319 - mrSges(5,3) * t240 - t260 * t286 + t276 * t320 + t340;
t164 = t327 * t169 + t331 * t178;
t287 = -mrSges(4,1) * t298 + mrSges(4,2) * t299;
t290 = -mrSges(4,2) * t323 + mrSges(4,3) * t298;
t161 = m(4) * t246 + mrSges(4,1) * t322 - mrSges(4,3) * t274 - t287 * t299 + t290 * t323 + t164;
t291 = mrSges(4,1) * t323 - mrSges(4,3) * t299;
t352 = t331 * t169 - t178 * t327;
t162 = m(4) * t247 - mrSges(4,2) * t322 + mrSges(4,3) * t273 + t287 * t298 - t291 * t323 + t352;
t155 = t332 * t161 + t328 * t162;
t288 = -t333 * g(3) - t362;
t296 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t329 + Ifges(3,2) * t333) * qJD(1);
t297 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t329 + Ifges(3,4) * t333) * qJD(1);
t279 = Ifges(4,4) * t299 + Ifges(4,2) * t298 + Ifges(4,6) * t323;
t280 = Ifges(4,1) * t299 + Ifges(4,4) * t298 + Ifges(4,5) * t323;
t349 = -mrSges(7,1) * t194 + mrSges(7,2) * t191;
t228 = Ifges(7,4) * t269 + Ifges(7,2) * t282 + Ifges(7,6) * t268;
t361 = -Ifges(6,5) * t269 + Ifges(6,6) * t268 - Ifges(6,3) * t282 - t228;
t171 = -mrSges(6,1) * t198 + mrSges(6,3) * t196 - pkin(5) * t188 + (t230 + t231) * t282 + t361 * t269 + (Ifges(6,6) - Ifges(7,6)) * t236 + (Ifges(6,4) - Ifges(7,5)) * t217 + (-Ifges(6,2) - Ifges(7,3)) * t216 + t349;
t344 = mrSges(7,2) * t193 - mrSges(7,3) * t194 + Ifges(7,1) * t217 + Ifges(7,4) * t236 + Ifges(7,5) * t216 + t282 * t226;
t173 = mrSges(6,2) * t198 - mrSges(6,3) * t195 + Ifges(6,1) * t217 - Ifges(6,4) * t216 + Ifges(6,5) * t236 - qJ(6) * t188 - t282 * t229 + t361 * t268 + t344;
t256 = Ifges(5,4) * t286 + Ifges(5,2) * t285 + Ifges(5,6) * t320;
t257 = Ifges(5,1) * t286 + Ifges(5,4) * t285 + Ifges(5,5) * t320;
t342 = -mrSges(5,1) * t203 + mrSges(5,2) * t204 - Ifges(5,5) * t240 - Ifges(5,6) * t239 - Ifges(5,3) * t319 - pkin(4) * t340 - pkin(10) * t351 - t365 * t171 - t326 * t173 - t286 * t256 + t285 * t257;
t339 = -mrSges(4,1) * t246 + mrSges(4,2) * t247 - Ifges(4,5) * t274 - Ifges(4,6) * t273 - Ifges(4,3) * t322 - pkin(3) * t164 - t299 * t279 + t298 * t280 + t342;
t366 = mrSges(3,1) * t288 - mrSges(3,2) * t289 + Ifges(3,5) * t306 + Ifges(3,6) * t307 + Ifges(3,3) * qJDD(2) + pkin(2) * t155 + (t296 * t329 - t297 * t333) * qJD(1) - t339;
t175 = t326 * t181 + t365 * t183;
t357 = qJD(1) * t333;
t305 = (-mrSges(3,1) * t333 + mrSges(3,2) * t329) * qJD(1);
t310 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t357;
t153 = m(3) * t288 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t306 + qJD(2) * t310 - t305 * t358 + t155;
t309 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t358;
t353 = -t161 * t328 + t332 * t162;
t154 = m(3) * t289 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t307 - qJD(2) * t309 + t305 * t357 + t353;
t354 = -t153 * t329 + t333 * t154;
t346 = m(5) * t221 - t239 * mrSges(5,1) + t240 * mrSges(5,2) - t285 * t276 + t286 * t277 + t175;
t255 = Ifges(5,5) * t286 + Ifges(5,6) * t285 + Ifges(5,3) * t320;
t156 = mrSges(5,2) * t221 - mrSges(5,3) * t203 + Ifges(5,1) * t240 + Ifges(5,4) * t239 + Ifges(5,5) * t319 - pkin(10) * t175 - t326 * t171 + t365 * t173 + t285 * t255 - t320 * t256;
t157 = -mrSges(5,1) * t221 + mrSges(5,3) * t204 + Ifges(5,4) * t240 + Ifges(5,2) * t239 + Ifges(5,6) * t319 - pkin(4) * t175 - t286 * t255 + t320 * t257 - t367;
t278 = Ifges(4,5) * t299 + Ifges(4,6) * t298 + Ifges(4,3) * t323;
t147 = -mrSges(4,1) * t275 + mrSges(4,3) * t247 + Ifges(4,4) * t274 + Ifges(4,2) * t273 + Ifges(4,6) * t322 - pkin(3) * t346 + pkin(9) * t352 + t327 * t156 + t331 * t157 - t299 * t278 + t323 * t280;
t148 = mrSges(4,2) * t275 - mrSges(4,3) * t246 + Ifges(4,1) * t274 + Ifges(4,4) * t273 + Ifges(4,5) * t322 - pkin(9) * t164 + t156 * t331 - t157 * t327 + t278 * t298 - t279 * t323;
t295 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t329 + Ifges(3,6) * t333) * qJD(1);
t300 = -t335 * pkin(7) + t347;
t341 = m(4) * t275 - t273 * mrSges(4,1) + mrSges(4,2) * t274 - t298 * t290 + t291 * t299 + t346;
t143 = -mrSges(3,1) * t300 + mrSges(3,3) * t289 + Ifges(3,4) * t306 + Ifges(3,2) * t307 + Ifges(3,6) * qJDD(2) - pkin(2) * t341 + pkin(8) * t353 + qJD(2) * t297 + t332 * t147 + t328 * t148 - t295 * t358;
t145 = mrSges(3,2) * t300 - mrSges(3,3) * t288 + Ifges(3,1) * t306 + Ifges(3,4) * t307 + Ifges(3,5) * qJDD(2) - pkin(8) * t155 - qJD(2) * t296 - t147 * t328 + t148 * t332 + t295 * t357;
t338 = -m(3) * t300 + t307 * mrSges(3,1) - mrSges(3,2) * t306 - t309 * t358 + t310 * t357 - t341;
t343 = mrSges(2,1) * t312 - mrSges(2,2) * t313 + Ifges(2,3) * qJDD(1) + pkin(1) * t338 + pkin(7) * t354 + t333 * t143 + t329 * t145;
t165 = m(2) * t312 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t335 + t338;
t151 = t153 * t333 + t154 * t329;
t149 = m(2) * t313 - mrSges(2,1) * t335 - qJDD(1) * mrSges(2,2) + t354;
t146 = mrSges(2,1) * g(3) + mrSges(2,3) * t313 + t335 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t151 - t366;
t141 = -mrSges(2,2) * g(3) - mrSges(2,3) * t312 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t335 - pkin(7) * t151 - t143 * t329 + t145 * t333;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t334 * t141 - t330 * t146 - pkin(6) * (t149 * t330 + t165 * t334), t141, t145, t148, t156, t173, -t228 * t268 + t344; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t330 * t141 + t334 * t146 + pkin(6) * (t149 * t334 - t165 * t330), t146, t143, t147, t157, t171, -t269 * t226 - t345; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t343, t343, t366, -t339, -t342, t367, Ifges(7,5) * t217 + Ifges(7,6) * t236 + Ifges(7,3) * t216 + t269 * t228 - t282 * t230 - t349;];
m_new  = t1;
