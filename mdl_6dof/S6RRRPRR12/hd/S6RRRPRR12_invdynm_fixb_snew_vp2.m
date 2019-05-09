% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 15:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:56:24
% EndTime: 2019-05-07 14:59:41
% DurationCPUTime: 99.14s
% Computational Cost: add. (1795494->399), mult. (3921666->518), div. (0->0), fcn. (3171519->14), ass. (0->165)
t339 = sin(pkin(6));
t345 = sin(qJ(2));
t349 = cos(qJ(2));
t366 = qJD(1) * qJD(2);
t324 = (-qJDD(1) * t349 + t345 * t366) * t339;
t376 = cos(qJ(3));
t375 = pkin(8) * t339;
t341 = cos(pkin(6));
t374 = t341 * g(3);
t373 = t339 * t345;
t372 = t339 * t349;
t371 = t341 * t345;
t370 = t341 * t349;
t368 = qJD(1) * t339;
t322 = (-pkin(2) * t349 - pkin(9) * t345) * t368;
t334 = qJD(1) * t341 + qJD(2);
t332 = t334 ^ 2;
t333 = qJDD(1) * t341 + qJDD(2);
t367 = qJD(1) * t349;
t346 = sin(qJ(1));
t350 = cos(qJ(1));
t330 = t346 * g(1) - g(2) * t350;
t351 = qJD(1) ^ 2;
t319 = qJDD(1) * pkin(1) + t351 * t375 + t330;
t331 = -g(1) * t350 - g(2) * t346;
t320 = -pkin(1) * t351 + qJDD(1) * t375 + t331;
t369 = t319 * t371 + t349 * t320;
t269 = -t332 * pkin(2) + t333 * pkin(9) + (-g(3) * t345 + t322 * t367) * t339 + t369;
t323 = (qJDD(1) * t345 + t349 * t366) * t339;
t270 = t324 * pkin(2) - t323 * pkin(9) - t374 + (-t319 + (pkin(2) * t345 - pkin(9) * t349) * t334 * qJD(1)) * t339;
t344 = sin(qJ(3));
t249 = t269 * t376 + t344 * t270;
t365 = t345 * t368;
t312 = -t334 * t376 + t344 * t365;
t313 = t344 * t334 + t365 * t376;
t295 = pkin(3) * t312 - qJ(4) * t313;
t316 = qJDD(3) + t324;
t364 = t339 * t367;
t329 = qJD(3) - t364;
t328 = t329 ^ 2;
t236 = -pkin(3) * t328 + qJ(4) * t316 - t295 * t312 + t249;
t293 = -g(3) * t372 + t319 * t370 - t345 * t320;
t268 = -t333 * pkin(2) - t332 * pkin(9) + t322 * t365 - t293;
t291 = qJD(3) * t313 + t323 * t344 - t333 * t376;
t292 = -t312 * qJD(3) + t323 * t376 + t344 * t333;
t240 = (t312 * t329 - t292) * qJ(4) + (t313 * t329 + t291) * pkin(3) + t268;
t338 = sin(pkin(12));
t340 = cos(pkin(12));
t301 = t313 * t340 + t329 * t338;
t224 = -0.2e1 * qJD(4) * t301 - t338 * t236 + t340 * t240;
t276 = t292 * t340 + t316 * t338;
t300 = -t313 * t338 + t329 * t340;
t215 = (t300 * t312 - t276) * pkin(10) + (t300 * t301 + t291) * pkin(4) + t224;
t225 = 0.2e1 * qJD(4) * t300 + t340 * t236 + t338 * t240;
t275 = -t292 * t338 + t316 * t340;
t281 = pkin(4) * t312 - pkin(10) * t301;
t299 = t300 ^ 2;
t217 = -pkin(4) * t299 + pkin(10) * t275 - t281 * t312 + t225;
t343 = sin(qJ(5));
t348 = cos(qJ(5));
t209 = t348 * t215 - t343 * t217;
t272 = t300 * t348 - t301 * t343;
t245 = qJD(5) * t272 + t275 * t343 + t276 * t348;
t273 = t300 * t343 + t301 * t348;
t289 = qJDD(5) + t291;
t311 = qJD(5) + t312;
t206 = (t272 * t311 - t245) * pkin(11) + (t272 * t273 + t289) * pkin(5) + t209;
t210 = t343 * t215 + t348 * t217;
t244 = -qJD(5) * t273 + t275 * t348 - t276 * t343;
t259 = pkin(5) * t311 - pkin(11) * t273;
t271 = t272 ^ 2;
t207 = -pkin(5) * t271 + pkin(11) * t244 - t259 * t311 + t210;
t342 = sin(qJ(6));
t347 = cos(qJ(6));
t204 = t206 * t347 - t207 * t342;
t254 = t272 * t347 - t273 * t342;
t223 = qJD(6) * t254 + t244 * t342 + t245 * t347;
t255 = t272 * t342 + t273 * t347;
t233 = -mrSges(7,1) * t254 + mrSges(7,2) * t255;
t310 = qJD(6) + t311;
t246 = -mrSges(7,2) * t310 + mrSges(7,3) * t254;
t284 = qJDD(6) + t289;
t198 = m(7) * t204 + mrSges(7,1) * t284 - mrSges(7,3) * t223 - t233 * t255 + t246 * t310;
t205 = t206 * t342 + t207 * t347;
t222 = -qJD(6) * t255 + t244 * t347 - t245 * t342;
t247 = mrSges(7,1) * t310 - mrSges(7,3) * t255;
t199 = m(7) * t205 - mrSges(7,2) * t284 + mrSges(7,3) * t222 + t233 * t254 - t247 * t310;
t192 = t347 * t198 + t342 * t199;
t256 = -mrSges(6,1) * t272 + mrSges(6,2) * t273;
t257 = -mrSges(6,2) * t311 + mrSges(6,3) * t272;
t189 = m(6) * t209 + mrSges(6,1) * t289 - mrSges(6,3) * t245 - t256 * t273 + t257 * t311 + t192;
t258 = mrSges(6,1) * t311 - mrSges(6,3) * t273;
t361 = -t198 * t342 + t347 * t199;
t190 = m(6) * t210 - mrSges(6,2) * t289 + mrSges(6,3) * t244 + t256 * t272 - t258 * t311 + t361;
t185 = t348 * t189 + t343 * t190;
t277 = -mrSges(5,1) * t300 + mrSges(5,2) * t301;
t279 = -mrSges(5,2) * t312 + mrSges(5,3) * t300;
t183 = m(5) * t224 + mrSges(5,1) * t291 - mrSges(5,3) * t276 - t277 * t301 + t279 * t312 + t185;
t280 = mrSges(5,1) * t312 - mrSges(5,3) * t301;
t362 = -t189 * t343 + t348 * t190;
t184 = m(5) * t225 - mrSges(5,2) * t291 + mrSges(5,3) * t275 + t277 * t300 - t280 * t312 + t362;
t179 = -t183 * t338 + t340 * t184;
t296 = mrSges(4,1) * t312 + mrSges(4,2) * t313;
t303 = mrSges(4,1) * t329 - mrSges(4,3) * t313;
t177 = m(4) * t249 - mrSges(4,2) * t316 - mrSges(4,3) * t291 - t296 * t312 - t303 * t329 + t179;
t248 = -t344 * t269 + t270 * t376;
t235 = -t316 * pkin(3) - t328 * qJ(4) + t313 * t295 + qJDD(4) - t248;
t227 = -t275 * pkin(4) - t299 * pkin(10) + t301 * t281 + t235;
t212 = -t244 * pkin(5) - t271 * pkin(11) + t273 * t259 + t227;
t360 = m(7) * t212 - t222 * mrSges(7,1) + t223 * mrSges(7,2) - t254 * t246 + t255 * t247;
t356 = m(6) * t227 - t244 * mrSges(6,1) + mrSges(6,2) * t245 - t272 * t257 + t258 * t273 + t360;
t202 = -m(5) * t235 + t275 * mrSges(5,1) - mrSges(5,2) * t276 + t300 * t279 - t280 * t301 - t356;
t302 = -mrSges(4,2) * t329 - mrSges(4,3) * t312;
t201 = m(4) * t248 + mrSges(4,1) * t316 - mrSges(4,3) * t292 - t296 * t313 + t302 * t329 + t202;
t172 = t344 * t177 + t201 * t376;
t294 = -g(3) * t373 + t369;
t317 = mrSges(3,1) * t334 - mrSges(3,3) * t365;
t321 = (-mrSges(3,1) * t349 + mrSges(3,2) * t345) * t368;
t363 = t177 * t376 - t201 * t344;
t170 = m(3) * t294 - mrSges(3,2) * t333 - mrSges(3,3) * t324 - t317 * t334 + t321 * t364 + t363;
t318 = -mrSges(3,2) * t334 + mrSges(3,3) * t364;
t178 = t183 * t340 + t184 * t338;
t355 = -m(4) * t268 - t291 * mrSges(4,1) - mrSges(4,2) * t292 - t312 * t302 - t303 * t313 - t178;
t174 = m(3) * t293 + mrSges(3,1) * t333 - mrSges(3,3) * t323 + t318 * t334 - t321 * t365 + t355;
t164 = t349 * t170 - t174 * t345;
t307 = -t339 * t319 - t374;
t171 = m(3) * t307 + t324 * mrSges(3,1) + t323 * mrSges(3,2) + (t317 * t345 - t318 * t349) * t368 + t172;
t161 = t170 * t371 - t171 * t339 + t174 * t370;
t228 = Ifges(7,5) * t255 + Ifges(7,6) * t254 + Ifges(7,3) * t310;
t230 = Ifges(7,1) * t255 + Ifges(7,4) * t254 + Ifges(7,5) * t310;
t193 = -mrSges(7,1) * t212 + mrSges(7,3) * t205 + Ifges(7,4) * t223 + Ifges(7,2) * t222 + Ifges(7,6) * t284 - t228 * t255 + t230 * t310;
t229 = Ifges(7,4) * t255 + Ifges(7,2) * t254 + Ifges(7,6) * t310;
t194 = mrSges(7,2) * t212 - mrSges(7,3) * t204 + Ifges(7,1) * t223 + Ifges(7,4) * t222 + Ifges(7,5) * t284 + t228 * t254 - t229 * t310;
t250 = Ifges(6,5) * t273 + Ifges(6,6) * t272 + Ifges(6,3) * t311;
t252 = Ifges(6,1) * t273 + Ifges(6,4) * t272 + Ifges(6,5) * t311;
t180 = -mrSges(6,1) * t227 + mrSges(6,3) * t210 + Ifges(6,4) * t245 + Ifges(6,2) * t244 + Ifges(6,6) * t289 - pkin(5) * t360 + pkin(11) * t361 + t347 * t193 + t342 * t194 - t273 * t250 + t311 * t252;
t251 = Ifges(6,4) * t273 + Ifges(6,2) * t272 + Ifges(6,6) * t311;
t181 = mrSges(6,2) * t227 - mrSges(6,3) * t209 + Ifges(6,1) * t245 + Ifges(6,4) * t244 + Ifges(6,5) * t289 - pkin(11) * t192 - t193 * t342 + t194 * t347 + t250 * t272 - t251 * t311;
t260 = Ifges(5,5) * t301 + Ifges(5,6) * t300 + Ifges(5,3) * t312;
t262 = Ifges(5,1) * t301 + Ifges(5,4) * t300 + Ifges(5,5) * t312;
t166 = -mrSges(5,1) * t235 + mrSges(5,3) * t225 + Ifges(5,4) * t276 + Ifges(5,2) * t275 + Ifges(5,6) * t291 - pkin(4) * t356 + pkin(10) * t362 + t348 * t180 + t343 * t181 - t301 * t260 + t312 * t262;
t261 = Ifges(5,4) * t301 + Ifges(5,2) * t300 + Ifges(5,6) * t312;
t167 = mrSges(5,2) * t235 - mrSges(5,3) * t224 + Ifges(5,1) * t276 + Ifges(5,4) * t275 + Ifges(5,5) * t291 - pkin(10) * t185 - t180 * t343 + t181 * t348 + t260 * t300 - t261 * t312;
t285 = Ifges(4,5) * t313 - Ifges(4,6) * t312 + Ifges(4,3) * t329;
t286 = Ifges(4,4) * t313 - Ifges(4,2) * t312 + Ifges(4,6) * t329;
t157 = mrSges(4,2) * t268 - mrSges(4,3) * t248 + Ifges(4,1) * t292 - Ifges(4,4) * t291 + Ifges(4,5) * t316 - qJ(4) * t178 - t166 * t338 + t167 * t340 - t285 * t312 - t286 * t329;
t287 = Ifges(4,1) * t313 - Ifges(4,4) * t312 + Ifges(4,5) * t329;
t357 = -mrSges(7,1) * t204 + mrSges(7,2) * t205 - Ifges(7,5) * t223 - Ifges(7,6) * t222 - Ifges(7,3) * t284 - t255 * t229 + t254 * t230;
t354 = -mrSges(6,1) * t209 + mrSges(6,2) * t210 - Ifges(6,5) * t245 - Ifges(6,6) * t244 - Ifges(6,3) * t289 - pkin(5) * t192 - t273 * t251 + t272 * t252 + t357;
t352 = mrSges(5,1) * t224 - mrSges(5,2) * t225 + Ifges(5,5) * t276 + Ifges(5,6) * t275 + pkin(4) * t185 + t301 * t261 - t300 * t262 - t354;
t165 = -pkin(3) * t178 + (-Ifges(5,3) - Ifges(4,2)) * t291 - t352 + t329 * t287 + Ifges(4,6) * t316 - t313 * t285 + Ifges(4,4) * t292 - mrSges(4,1) * t268 + mrSges(4,3) * t249;
t305 = Ifges(3,6) * t334 + (Ifges(3,4) * t345 + Ifges(3,2) * t349) * t368;
t306 = Ifges(3,5) * t334 + (Ifges(3,1) * t345 + Ifges(3,4) * t349) * t368;
t152 = Ifges(3,5) * t323 - Ifges(3,6) * t324 + Ifges(3,3) * t333 + mrSges(3,1) * t293 - mrSges(3,2) * t294 + t344 * t157 + t376 * t165 + pkin(2) * t355 + pkin(9) * t363 + (t305 * t345 - t306 * t349) * t368;
t304 = Ifges(3,3) * t334 + (Ifges(3,5) * t345 + Ifges(3,6) * t349) * t368;
t154 = mrSges(3,2) * t307 - mrSges(3,3) * t293 + Ifges(3,1) * t323 - Ifges(3,4) * t324 + Ifges(3,5) * t333 - pkin(9) * t172 + t157 * t376 - t344 * t165 + t304 * t364 - t334 * t305;
t353 = mrSges(4,1) * t248 - mrSges(4,2) * t249 + Ifges(4,5) * t292 - Ifges(4,6) * t291 + Ifges(4,3) * t316 + pkin(3) * t202 + qJ(4) * t179 + t340 * t166 + t338 * t167 + t313 * t286 + t312 * t287;
t156 = -mrSges(3,1) * t307 + mrSges(3,3) * t294 + Ifges(3,4) * t323 - Ifges(3,2) * t324 + Ifges(3,6) * t333 - pkin(2) * t172 - t304 * t365 + t334 * t306 - t353;
t358 = mrSges(2,1) * t330 - mrSges(2,2) * t331 + Ifges(2,3) * qJDD(1) + pkin(1) * t161 + t341 * t152 + t154 * t373 + t156 * t372 + t164 * t375;
t162 = m(2) * t331 - mrSges(2,1) * t351 - qJDD(1) * mrSges(2,2) + t164;
t160 = t341 * t171 + (t170 * t345 + t174 * t349) * t339;
t158 = m(2) * t330 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t351 + t161;
t150 = -mrSges(2,2) * g(3) - mrSges(2,3) * t330 + Ifges(2,5) * qJDD(1) - t351 * Ifges(2,6) + t349 * t154 - t345 * t156 + (-t160 * t339 - t161 * t341) * pkin(8);
t149 = mrSges(2,1) * g(3) + mrSges(2,3) * t331 + t351 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t160 - t339 * t152 + (pkin(8) * t164 + t154 * t345 + t156 * t349) * t341;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t350 * t150 - t346 * t149 - pkin(7) * (t158 * t350 + t162 * t346), t150, t154, t157, t167, t181, t194; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t346 * t150 + t350 * t149 + pkin(7) * (-t158 * t346 + t162 * t350), t149, t156, t165, t166, t180, t193; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t358, t358, t152, t353, Ifges(5,3) * t291 + t352, -t354, -t357;];
m_new  = t1;
