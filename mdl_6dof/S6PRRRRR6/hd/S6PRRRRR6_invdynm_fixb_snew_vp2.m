% Calculate vector of cutting torques with Newton-Euler for
% S6PRRRRR6
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 12:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6PRRRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 12:32:22
% EndTime: 2019-05-05 12:35:18
% DurationCPUTime: 172.44s
% Computational Cost: add. (3251181->369), mult. (7649050->503), div. (0->0), fcn. (6460399->18), ass. (0->173)
t321 = sin(pkin(14));
t325 = cos(pkin(14));
t314 = g(1) * t321 - g(2) * t325;
t315 = -g(1) * t325 - g(2) * t321;
t320 = -g(3) + qJDD(1);
t338 = cos(qJ(2));
t328 = cos(pkin(6));
t333 = sin(qJ(2));
t358 = t328 * t333;
t324 = sin(pkin(6));
t364 = t324 * t333;
t288 = t314 * t358 + t315 * t338 + t320 * t364;
t339 = qJD(2) ^ 2;
t323 = sin(pkin(7));
t371 = pkin(10) * t323;
t285 = -pkin(2) * t339 + qJDD(2) * t371 + t288;
t327 = cos(pkin(7));
t319 = qJD(2) * t327 + qJD(3);
t322 = sin(pkin(8));
t326 = cos(pkin(8));
t337 = cos(qJ(3));
t355 = qJD(2) * t323;
t351 = t337 * t355;
t300 = (t319 * t322 + t326 * t351) * pkin(11);
t332 = sin(qJ(3));
t370 = pkin(11) * t322;
t304 = (-pkin(3) * t337 - t332 * t370) * t355;
t354 = qJD(2) * qJD(3);
t308 = (qJDD(2) * t332 + t337 * t354) * t323;
t318 = qJDD(2) * t327 + qJDD(3);
t357 = t328 * t338;
t363 = t324 * t338;
t287 = t314 * t357 - t315 * t333 + t320 * t363;
t284 = qJDD(2) * pkin(2) + t339 * t371 + t287;
t303 = -t314 * t324 + t320 * t328;
t359 = t327 * t337;
t365 = t323 * t337;
t356 = t284 * t359 + t303 * t365;
t369 = pkin(11) * t326;
t243 = -t308 * t369 + t318 * pkin(3) + t319 * t300 + (-t304 * t355 - t285) * t332 + t356;
t360 = t327 * t332;
t366 = t323 * t332;
t261 = t284 * t360 + t285 * t337 + t303 * t366;
t352 = t332 * t355;
t302 = pkin(3) * t319 - t352 * t369;
t309 = (qJDD(2) * t337 - t332 * t354) * t323;
t348 = t309 * t326 + t318 * t322;
t244 = pkin(11) * t348 - t319 * t302 + t304 * t351 + t261;
t298 = t327 * t303;
t257 = -t308 * t370 - t309 * pkin(3) + t298 + (-t284 + (-t300 * t337 + t302 * t332) * qJD(2)) * t323;
t331 = sin(qJ(4));
t336 = cos(qJ(4));
t229 = -t331 * t244 + (t243 * t326 + t257 * t322) * t336;
t361 = t326 * t337;
t368 = t322 * t331;
t291 = t319 * t368 + (t331 * t361 + t332 * t336) * t355;
t271 = -t291 * qJD(4) - t331 * t308 + t336 * t348;
t367 = t322 * t336;
t290 = (-t331 * t332 + t336 * t361) * t355 + t319 * t367;
t362 = t326 * t331;
t230 = t243 * t362 + t244 * t336 + t257 * t368;
t273 = -mrSges(5,1) * t290 + mrSges(5,2) * t291;
t301 = t319 * t326 - t322 * t351 + qJD(4);
t280 = mrSges(5,1) * t301 - mrSges(5,3) * t291;
t292 = -t309 * t322 + t318 * t326 + qJDD(4);
t274 = -pkin(4) * t290 - pkin(12) * t291;
t299 = t301 ^ 2;
t226 = -pkin(4) * t299 + pkin(12) * t292 + t274 * t290 + t230;
t234 = -t322 * t243 + t257 * t326;
t272 = t290 * qJD(4) + t336 * t308 + t331 * t348;
t228 = (-t290 * t301 - t272) * pkin(12) + (t291 * t301 - t271) * pkin(4) + t234;
t330 = sin(qJ(5));
t335 = cos(qJ(5));
t222 = t226 * t335 + t228 * t330;
t277 = -t291 * t330 + t301 * t335;
t278 = t291 * t335 + t301 * t330;
t259 = -pkin(5) * t277 - pkin(13) * t278;
t270 = qJDD(5) - t271;
t289 = qJD(5) - t290;
t286 = t289 ^ 2;
t220 = -pkin(5) * t286 + pkin(13) * t270 + t259 * t277 + t222;
t225 = -t292 * pkin(4) - t299 * pkin(12) + t274 * t291 - t229;
t247 = -qJD(5) * t278 - t272 * t330 + t292 * t335;
t248 = qJD(5) * t277 + t272 * t335 + t292 * t330;
t223 = (-t277 * t289 - t248) * pkin(13) + (t278 * t289 - t247) * pkin(5) + t225;
t329 = sin(qJ(6));
t334 = cos(qJ(6));
t217 = -t220 * t329 + t223 * t334;
t263 = -t278 * t329 + t289 * t334;
t233 = qJD(6) * t263 + t248 * t334 + t270 * t329;
t264 = t278 * t334 + t289 * t329;
t239 = -mrSges(7,1) * t263 + mrSges(7,2) * t264;
t246 = qJDD(6) - t247;
t276 = qJD(6) - t277;
t249 = -mrSges(7,2) * t276 + mrSges(7,3) * t263;
t214 = m(7) * t217 + mrSges(7,1) * t246 - mrSges(7,3) * t233 - t239 * t264 + t249 * t276;
t218 = t220 * t334 + t223 * t329;
t232 = -qJD(6) * t264 - t248 * t329 + t270 * t334;
t250 = mrSges(7,1) * t276 - mrSges(7,3) * t264;
t215 = m(7) * t218 - mrSges(7,2) * t246 + mrSges(7,3) * t232 + t239 * t263 - t250 * t276;
t208 = -t214 * t329 + t215 * t334;
t258 = -mrSges(6,1) * t277 + mrSges(6,2) * t278;
t266 = mrSges(6,1) * t289 - mrSges(6,3) * t278;
t206 = m(6) * t222 - mrSges(6,2) * t270 + mrSges(6,3) * t247 + t258 * t277 - t266 * t289 + t208;
t221 = -t226 * t330 + t228 * t335;
t219 = -pkin(5) * t270 - pkin(13) * t286 + t259 * t278 - t221;
t216 = -m(7) * t219 + mrSges(7,1) * t232 - mrSges(7,2) * t233 + t249 * t263 - t250 * t264;
t265 = -mrSges(6,2) * t289 + mrSges(6,3) * t277;
t212 = m(6) * t221 + mrSges(6,1) * t270 - mrSges(6,3) * t248 - t258 * t278 + t265 * t289 + t216;
t350 = t206 * t335 - t212 * t330;
t197 = m(5) * t230 - mrSges(5,2) * t292 + mrSges(5,3) * t271 + t273 * t290 - t280 * t301 + t350;
t200 = t206 * t330 + t212 * t335;
t279 = -mrSges(5,2) * t301 + mrSges(5,3) * t290;
t199 = m(5) * t234 - mrSges(5,1) * t271 + mrSges(5,2) * t272 - t279 * t290 + t280 * t291 + t200;
t207 = t214 * t334 + t215 * t329;
t342 = -m(6) * t225 + mrSges(6,1) * t247 - mrSges(6,2) * t248 + t265 * t277 - t266 * t278 - t207;
t203 = m(5) * t229 + mrSges(5,1) * t292 - mrSges(5,3) * t272 - t273 * t291 + t279 * t301 + t342;
t186 = t203 * t326 * t336 + t197 * t362 - t199 * t322;
t260 = -t285 * t332 + t356;
t306 = -mrSges(4,2) * t319 + mrSges(4,3) * t351;
t307 = (-mrSges(4,1) * t337 + mrSges(4,2) * t332) * t355;
t182 = m(4) * t260 + mrSges(4,1) * t318 - mrSges(4,3) * t308 + t306 * t319 - t307 * t352 + t186;
t185 = t197 * t368 + t199 * t326 + t203 * t367;
t275 = -t323 * t284 + t298;
t305 = mrSges(4,1) * t319 - mrSges(4,3) * t352;
t184 = m(4) * t275 - t309 * mrSges(4,1) + t308 * mrSges(4,2) + (t305 * t332 - t306 * t337) * t355 + t185;
t191 = t197 * t336 - t203 * t331;
t190 = m(4) * t261 - mrSges(4,2) * t318 + mrSges(4,3) * t309 - t305 * t319 + t307 * t351 + t191;
t172 = t182 * t359 - t184 * t323 + t190 * t360;
t169 = m(3) * t287 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t339 + t172;
t176 = -t182 * t332 + t190 * t337;
t175 = m(3) * t288 - mrSges(3,1) * t339 - qJDD(2) * mrSges(3,2) + t176;
t167 = -t169 * t333 + t175 * t338;
t372 = pkin(9) * t167;
t171 = t182 * t365 + t184 * t327 + t190 * t366;
t170 = m(3) * t303 + t171;
t161 = t169 * t357 - t170 * t324 + t175 * t358;
t235 = Ifges(7,5) * t264 + Ifges(7,6) * t263 + Ifges(7,3) * t276;
t237 = Ifges(7,1) * t264 + Ifges(7,4) * t263 + Ifges(7,5) * t276;
t209 = -mrSges(7,1) * t219 + mrSges(7,3) * t218 + Ifges(7,4) * t233 + Ifges(7,2) * t232 + Ifges(7,6) * t246 - t235 * t264 + t237 * t276;
t236 = Ifges(7,4) * t264 + Ifges(7,2) * t263 + Ifges(7,6) * t276;
t210 = mrSges(7,2) * t219 - mrSges(7,3) * t217 + Ifges(7,1) * t233 + Ifges(7,4) * t232 + Ifges(7,5) * t246 + t235 * t263 - t236 * t276;
t253 = Ifges(6,5) * t278 + Ifges(6,6) * t277 + Ifges(6,3) * t289;
t254 = Ifges(6,4) * t278 + Ifges(6,2) * t277 + Ifges(6,6) * t289;
t192 = mrSges(6,2) * t225 - mrSges(6,3) * t221 + Ifges(6,1) * t248 + Ifges(6,4) * t247 + Ifges(6,5) * t270 - pkin(13) * t207 - t209 * t329 + t210 * t334 + t253 * t277 - t254 * t289;
t255 = Ifges(6,1) * t278 + Ifges(6,4) * t277 + Ifges(6,5) * t289;
t341 = mrSges(7,1) * t217 - mrSges(7,2) * t218 + Ifges(7,5) * t233 + Ifges(7,6) * t232 + Ifges(7,3) * t246 + t236 * t264 - t237 * t263;
t193 = -mrSges(6,1) * t225 + mrSges(6,3) * t222 + Ifges(6,4) * t248 + Ifges(6,2) * t247 + Ifges(6,6) * t270 - pkin(5) * t207 - t253 * t278 + t255 * t289 - t341;
t268 = Ifges(5,4) * t291 + Ifges(5,2) * t290 + Ifges(5,6) * t301;
t269 = Ifges(5,1) * t291 + Ifges(5,4) * t290 + Ifges(5,5) * t301;
t177 = mrSges(5,1) * t229 - mrSges(5,2) * t230 + Ifges(5,5) * t272 + Ifges(5,6) * t271 + Ifges(5,3) * t292 + pkin(4) * t342 + pkin(12) * t350 + t330 * t192 + t335 * t193 + t291 * t268 - t290 * t269;
t295 = Ifges(4,3) * t319 + (Ifges(4,5) * t332 + Ifges(4,6) * t337) * t355;
t297 = Ifges(4,5) * t319 + (Ifges(4,1) * t332 + Ifges(4,4) * t337) * t355;
t267 = Ifges(5,5) * t291 + Ifges(5,6) * t290 + Ifges(5,3) * t301;
t178 = mrSges(5,2) * t234 - mrSges(5,3) * t229 + Ifges(5,1) * t272 + Ifges(5,4) * t271 + Ifges(5,5) * t292 - pkin(12) * t200 + t192 * t335 - t193 * t330 + t267 * t290 - t268 * t301;
t340 = mrSges(6,1) * t221 - mrSges(6,2) * t222 + Ifges(6,5) * t248 + Ifges(6,6) * t247 + Ifges(6,3) * t270 + pkin(5) * t216 + pkin(13) * t208 + t209 * t334 + t210 * t329 + t254 * t278 - t255 * t277;
t179 = -mrSges(5,1) * t234 + mrSges(5,3) * t230 + Ifges(5,4) * t272 + Ifges(5,2) * t271 + Ifges(5,6) * t292 - pkin(4) * t200 - t267 * t291 + t269 * t301 - t340;
t344 = pkin(11) * t191 + t178 * t331 + t179 * t336;
t163 = -mrSges(4,1) * t275 + mrSges(4,3) * t261 + Ifges(4,4) * t308 + Ifges(4,2) * t309 + Ifges(4,6) * t318 - pkin(3) * t185 - t322 * t177 - t295 * t352 + t319 * t297 + t326 * t344;
t296 = Ifges(4,6) * t319 + (Ifges(4,4) * t332 + Ifges(4,2) * t337) * t355;
t164 = t295 * t351 + mrSges(4,2) * t275 - mrSges(4,3) * t260 + Ifges(4,1) * t308 + Ifges(4,4) * t309 + Ifges(4,5) * t318 + t336 * t178 - t331 * t179 - t319 * t296 + (-t185 * t322 - t186 * t326) * pkin(11);
t345 = pkin(10) * t176 + t163 * t337 + t164 * t332;
t162 = mrSges(4,1) * t260 - mrSges(4,2) * t261 + Ifges(4,5) * t308 + Ifges(4,6) * t309 + Ifges(4,3) * t318 + pkin(3) * t186 + t326 * t177 + (t296 * t332 - t297 * t337) * t355 + t344 * t322;
t153 = mrSges(3,1) * t287 - mrSges(3,2) * t288 + Ifges(3,3) * qJDD(2) + pkin(2) * t172 + t327 * t162 + t323 * t345;
t155 = -mrSges(3,1) * t303 + mrSges(3,3) * t288 + t339 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t171 - t323 * t162 + t327 * t345;
t157 = mrSges(3,2) * t303 - mrSges(3,3) * t287 + Ifges(3,5) * qJDD(2) - t339 * Ifges(3,6) - t332 * t163 + t337 * t164 + (-t171 * t323 - t172 * t327) * pkin(10);
t343 = mrSges(2,1) * t314 - mrSges(2,2) * t315 + pkin(1) * t161 + t153 * t328 + t155 * t363 + t157 * t364 + t324 * t372;
t165 = m(2) * t315 + t167;
t160 = t328 * t170 + (t169 * t338 + t175 * t333) * t324;
t158 = m(2) * t314 + t161;
t151 = mrSges(2,2) * t320 - mrSges(2,3) * t314 - t333 * t155 + t338 * t157 + (-t160 * t324 - t161 * t328) * pkin(9);
t150 = -mrSges(2,1) * t320 + mrSges(2,3) * t315 - pkin(1) * t160 - t324 * t153 + (t155 * t338 + t157 * t333 + t372) * t328;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t325 * t151 - t321 * t150 - qJ(1) * (t158 * t325 + t165 * t321), t151, t157, t164, t178, t192, t210; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t321 * t151 + t325 * t150 + qJ(1) * (-t158 * t321 + t165 * t325), t150, t155, t163, t179, t193, t209; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t343, t343, t153, t162, t177, t340, t341;];
m_new  = t1;
