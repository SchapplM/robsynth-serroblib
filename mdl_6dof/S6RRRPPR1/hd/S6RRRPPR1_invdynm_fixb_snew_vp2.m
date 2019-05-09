% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-05-07 04:11
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR1_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:06:55
% EndTime: 2019-05-07 04:07:52
% DurationCPUTime: 43.75s
% Computational Cost: add. (771810->386), mult. (1760928->489), div. (0->0), fcn. (1313839->12), ass. (0->152)
t374 = -2 * qJD(4);
t340 = sin(qJ(2));
t344 = cos(qJ(2));
t366 = qJD(1) * qJD(2);
t319 = qJDD(1) * t340 + t344 * t366;
t341 = sin(qJ(1));
t345 = cos(qJ(1));
t326 = -g(1) * t345 - g(2) * t341;
t346 = qJD(1) ^ 2;
t314 = -pkin(1) * t346 + qJDD(1) * pkin(7) + t326;
t370 = t340 * t314;
t372 = pkin(2) * t346;
t273 = qJDD(2) * pkin(2) - t319 * pkin(8) - t370 + (pkin(8) * t366 + t340 * t372 - g(3)) * t344;
t301 = -g(3) * t340 + t344 * t314;
t320 = qJDD(1) * t344 - t340 * t366;
t369 = qJD(1) * t340;
t324 = qJD(2) * pkin(2) - pkin(8) * t369;
t334 = t344 ^ 2;
t274 = pkin(8) * t320 - qJD(2) * t324 - t334 * t372 + t301;
t339 = sin(qJ(3));
t343 = cos(qJ(3));
t246 = t343 * t273 - t339 * t274;
t311 = (-t339 * t340 + t343 * t344) * qJD(1);
t285 = qJD(3) * t311 + t319 * t343 + t320 * t339;
t312 = (t339 * t344 + t340 * t343) * qJD(1);
t331 = qJDD(2) + qJDD(3);
t332 = qJD(2) + qJD(3);
t227 = (t311 * t332 - t285) * qJ(4) + (t311 * t312 + t331) * pkin(3) + t246;
t247 = t339 * t273 + t343 * t274;
t284 = -qJD(3) * t312 - t319 * t339 + t320 * t343;
t303 = pkin(3) * t332 - qJ(4) * t312;
t307 = t311 ^ 2;
t231 = -pkin(3) * t307 + qJ(4) * t284 - t303 * t332 + t247;
t336 = sin(pkin(10));
t371 = cos(pkin(10));
t298 = t336 * t311 + t371 * t312;
t211 = t371 * t227 - t336 * t231 + t298 * t374;
t297 = -t371 * t311 + t312 * t336;
t212 = t336 * t227 + t371 * t231 + t297 * t374;
t257 = -t371 * t284 + t285 * t336;
t268 = mrSges(5,1) * t297 + mrSges(5,2) * t298;
t288 = mrSges(5,1) * t332 - mrSges(5,3) * t298;
t267 = pkin(4) * t297 - qJ(5) * t298;
t330 = t332 ^ 2;
t209 = -pkin(4) * t330 + qJ(5) * t331 - t267 * t297 + t212;
t325 = t341 * g(1) - t345 * g(2);
t358 = -qJDD(1) * pkin(1) - t325;
t286 = -t320 * pkin(2) + t324 * t369 + (-pkin(8) * t334 - pkin(7)) * t346 + t358;
t233 = -t284 * pkin(3) - t307 * qJ(4) + t312 * t303 + qJDD(4) + t286;
t258 = t336 * t284 + t371 * t285;
t215 = (t297 * t332 - t258) * qJ(5) + (t298 * t332 + t257) * pkin(4) + t233;
t335 = sin(pkin(11));
t337 = cos(pkin(11));
t280 = t298 * t337 + t332 * t335;
t204 = -0.2e1 * qJD(5) * t280 - t335 * t209 + t337 * t215;
t245 = t258 * t337 + t331 * t335;
t279 = -t298 * t335 + t332 * t337;
t202 = (t279 * t297 - t245) * pkin(9) + (t279 * t280 + t257) * pkin(5) + t204;
t205 = 0.2e1 * qJD(5) * t279 + t337 * t209 + t335 * t215;
t244 = -t258 * t335 + t331 * t337;
t261 = pkin(5) * t297 - pkin(9) * t280;
t278 = t279 ^ 2;
t203 = -pkin(5) * t278 + pkin(9) * t244 - t261 * t297 + t205;
t338 = sin(qJ(6));
t342 = cos(qJ(6));
t200 = t202 * t342 - t203 * t338;
t248 = t279 * t342 - t280 * t338;
t221 = qJD(6) * t248 + t244 * t338 + t245 * t342;
t249 = t279 * t338 + t280 * t342;
t228 = -mrSges(7,1) * t248 + mrSges(7,2) * t249;
t291 = qJD(6) + t297;
t234 = -mrSges(7,2) * t291 + mrSges(7,3) * t248;
t255 = qJDD(6) + t257;
t195 = m(7) * t200 + mrSges(7,1) * t255 - mrSges(7,3) * t221 - t228 * t249 + t234 * t291;
t201 = t202 * t338 + t203 * t342;
t220 = -qJD(6) * t249 + t244 * t342 - t245 * t338;
t235 = mrSges(7,1) * t291 - mrSges(7,3) * t249;
t196 = m(7) * t201 - mrSges(7,2) * t255 + mrSges(7,3) * t220 + t228 * t248 - t235 * t291;
t187 = t342 * t195 + t338 * t196;
t256 = -mrSges(6,1) * t279 + mrSges(6,2) * t280;
t259 = -mrSges(6,2) * t297 + mrSges(6,3) * t279;
t185 = m(6) * t204 + mrSges(6,1) * t257 - mrSges(6,3) * t245 - t256 * t280 + t259 * t297 + t187;
t260 = mrSges(6,1) * t297 - mrSges(6,3) * t280;
t361 = -t195 * t338 + t342 * t196;
t186 = m(6) * t205 - mrSges(6,2) * t257 + mrSges(6,3) * t244 + t256 * t279 - t260 * t297 + t361;
t362 = -t185 * t335 + t337 * t186;
t178 = m(5) * t212 - mrSges(5,2) * t331 - mrSges(5,3) * t257 - t268 * t297 - t288 * t332 + t362;
t287 = -mrSges(5,2) * t332 - mrSges(5,3) * t297;
t208 = -t331 * pkin(4) - t330 * qJ(5) + t298 * t267 + qJDD(5) - t211;
t206 = -t244 * pkin(5) - t278 * pkin(9) + t280 * t261 + t208;
t355 = m(7) * t206 - t220 * mrSges(7,1) + mrSges(7,2) * t221 - t248 * t234 + t235 * t249;
t351 = -m(6) * t208 + t244 * mrSges(6,1) - mrSges(6,2) * t245 + t279 * t259 - t260 * t280 - t355;
t191 = m(5) * t211 + mrSges(5,1) * t331 - mrSges(5,3) * t258 - t268 * t298 + t287 * t332 + t351;
t169 = t336 * t178 + t371 * t191;
t299 = -mrSges(4,1) * t311 + mrSges(4,2) * t312;
t302 = -mrSges(4,2) * t332 + mrSges(4,3) * t311;
t166 = m(4) * t246 + mrSges(4,1) * t331 - mrSges(4,3) * t285 - t299 * t312 + t302 * t332 + t169;
t304 = mrSges(4,1) * t332 - mrSges(4,3) * t312;
t363 = t371 * t178 - t191 * t336;
t167 = m(4) * t247 - mrSges(4,2) * t331 + mrSges(4,3) * t284 + t299 * t311 - t304 * t332 + t363;
t161 = t343 * t166 + t339 * t167;
t300 = -t344 * g(3) - t370;
t309 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t340 + Ifges(3,2) * t344) * qJD(1);
t310 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t340 + Ifges(3,4) * t344) * qJD(1);
t293 = Ifges(4,4) * t312 + Ifges(4,2) * t311 + Ifges(4,6) * t332;
t294 = Ifges(4,1) * t312 + Ifges(4,4) * t311 + Ifges(4,5) * t332;
t222 = Ifges(7,5) * t249 + Ifges(7,6) * t248 + Ifges(7,3) * t291;
t224 = Ifges(7,1) * t249 + Ifges(7,4) * t248 + Ifges(7,5) * t291;
t188 = -mrSges(7,1) * t206 + mrSges(7,3) * t201 + Ifges(7,4) * t221 + Ifges(7,2) * t220 + Ifges(7,6) * t255 - t222 * t249 + t224 * t291;
t223 = Ifges(7,4) * t249 + Ifges(7,2) * t248 + Ifges(7,6) * t291;
t189 = mrSges(7,2) * t206 - mrSges(7,3) * t200 + Ifges(7,1) * t221 + Ifges(7,4) * t220 + Ifges(7,5) * t255 + t222 * t248 - t223 * t291;
t236 = Ifges(6,5) * t280 + Ifges(6,6) * t279 + Ifges(6,3) * t297;
t238 = Ifges(6,1) * t280 + Ifges(6,4) * t279 + Ifges(6,5) * t297;
t171 = -mrSges(6,1) * t208 + mrSges(6,3) * t205 + Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t257 - pkin(5) * t355 + pkin(9) * t361 + t342 * t188 + t338 * t189 - t280 * t236 + t297 * t238;
t237 = Ifges(6,4) * t280 + Ifges(6,2) * t279 + Ifges(6,6) * t297;
t173 = mrSges(6,2) * t208 - mrSges(6,3) * t204 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t257 - pkin(9) * t187 - t188 * t338 + t189 * t342 + t236 * t279 - t237 * t297;
t263 = Ifges(5,4) * t298 - Ifges(5,2) * t297 + Ifges(5,6) * t332;
t264 = Ifges(5,1) * t298 - Ifges(5,4) * t297 + Ifges(5,5) * t332;
t353 = -mrSges(5,1) * t211 + mrSges(5,2) * t212 - Ifges(5,5) * t258 + Ifges(5,6) * t257 - Ifges(5,3) * t331 - pkin(4) * t351 - qJ(5) * t362 - t337 * t171 - t335 * t173 - t298 * t263 - t297 * t264;
t350 = -mrSges(4,1) * t246 + mrSges(4,2) * t247 - Ifges(4,5) * t285 - Ifges(4,6) * t284 - Ifges(4,3) * t331 - pkin(3) * t169 - t312 * t293 + t311 * t294 + t353;
t373 = mrSges(3,1) * t300 - mrSges(3,2) * t301 + Ifges(3,5) * t319 + Ifges(3,6) * t320 + Ifges(3,3) * qJDD(2) + pkin(2) * t161 + (t309 * t340 - t310 * t344) * qJD(1) - t350;
t180 = t337 * t185 + t335 * t186;
t368 = qJD(1) * t344;
t318 = (-mrSges(3,1) * t344 + mrSges(3,2) * t340) * qJD(1);
t323 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t368;
t159 = m(3) * t300 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t319 + qJD(2) * t323 - t318 * t369 + t161;
t322 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t369;
t364 = -t166 * t339 + t343 * t167;
t160 = m(3) * t301 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t320 - qJD(2) * t322 + t318 * t368 + t364;
t365 = -t159 * t340 + t344 * t160;
t357 = m(5) * t233 + t257 * mrSges(5,1) + t258 * mrSges(5,2) + t297 * t287 + t298 * t288 + t180;
t262 = Ifges(5,5) * t298 - Ifges(5,6) * t297 + Ifges(5,3) * t332;
t157 = mrSges(5,2) * t233 - mrSges(5,3) * t211 + Ifges(5,1) * t258 - Ifges(5,4) * t257 + Ifges(5,5) * t331 - qJ(5) * t180 - t171 * t335 + t173 * t337 - t262 * t297 - t263 * t332;
t354 = -mrSges(7,1) * t200 + mrSges(7,2) * t201 - Ifges(7,5) * t221 - Ifges(7,6) * t220 - Ifges(7,3) * t255 - t249 * t223 + t248 * t224;
t348 = -mrSges(6,1) * t204 + mrSges(6,2) * t205 - Ifges(6,5) * t245 - Ifges(6,6) * t244 - pkin(5) * t187 - t280 * t237 + t279 * t238 + t354;
t162 = mrSges(5,3) * t212 - pkin(4) * t180 + t348 + Ifges(5,4) * t258 - mrSges(5,1) * t233 - t298 * t262 + Ifges(5,6) * t331 + t332 * t264 + (-Ifges(5,2) - Ifges(6,3)) * t257;
t292 = Ifges(4,5) * t312 + Ifges(4,6) * t311 + Ifges(4,3) * t332;
t152 = -mrSges(4,1) * t286 + mrSges(4,3) * t247 + Ifges(4,4) * t285 + Ifges(4,2) * t284 + Ifges(4,6) * t331 - pkin(3) * t357 + qJ(4) * t363 + t336 * t157 + t371 * t162 - t312 * t292 + t332 * t294;
t153 = mrSges(4,2) * t286 - mrSges(4,3) * t246 + Ifges(4,1) * t285 + Ifges(4,4) * t284 + Ifges(4,5) * t331 - qJ(4) * t169 + t371 * t157 - t336 * t162 + t311 * t292 - t332 * t293;
t308 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t340 + Ifges(3,6) * t344) * qJD(1);
t313 = -t346 * pkin(7) + t358;
t352 = m(4) * t286 - t284 * mrSges(4,1) + mrSges(4,2) * t285 - t311 * t302 + t304 * t312 + t357;
t148 = -mrSges(3,1) * t313 + mrSges(3,3) * t301 + Ifges(3,4) * t319 + Ifges(3,2) * t320 + Ifges(3,6) * qJDD(2) - pkin(2) * t352 + pkin(8) * t364 + qJD(2) * t310 + t343 * t152 + t339 * t153 - t308 * t369;
t150 = mrSges(3,2) * t313 - mrSges(3,3) * t300 + Ifges(3,1) * t319 + Ifges(3,4) * t320 + Ifges(3,5) * qJDD(2) - pkin(8) * t161 - qJD(2) * t309 - t152 * t339 + t153 * t343 + t308 * t368;
t349 = -m(3) * t313 + t320 * mrSges(3,1) - mrSges(3,2) * t319 - t322 * t369 + t323 * t368 - t352;
t356 = mrSges(2,1) * t325 - mrSges(2,2) * t326 + Ifges(2,3) * qJDD(1) + pkin(1) * t349 + pkin(7) * t365 + t344 * t148 + t340 * t150;
t174 = m(2) * t325 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t346 + t349;
t156 = t159 * t344 + t160 * t340;
t154 = m(2) * t326 - mrSges(2,1) * t346 - qJDD(1) * mrSges(2,2) + t365;
t151 = mrSges(2,1) * g(3) + mrSges(2,3) * t326 + t346 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t156 - t373;
t146 = -mrSges(2,2) * g(3) - mrSges(2,3) * t325 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t346 - pkin(7) * t156 - t148 * t340 + t150 * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t345 * t146 - t341 * t151 - pkin(6) * (t154 * t341 + t174 * t345), t146, t150, t153, t157, t173, t189; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t341 * t146 + t345 * t151 + pkin(6) * (t154 * t345 - t174 * t341), t151, t148, t152, t162, t171, t188; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t356, t356, t373, -t350, -t353, Ifges(6,3) * t257 - t348, -t354;];
m_new  = t1;
