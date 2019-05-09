% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 10:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR5_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 10:04:24
% EndTime: 2019-05-08 10:07:45
% DurationCPUTime: 99.53s
% Computational Cost: add. (1837162->399), mult. (3980117->515), div. (0->0), fcn. (3287594->14), ass. (0->167)
t341 = sin(pkin(6));
t380 = pkin(8) * t341;
t342 = cos(pkin(6));
t379 = g(3) * t342;
t347 = sin(qJ(2));
t378 = t341 * t347;
t353 = cos(qJ(2));
t377 = t341 * t353;
t376 = t342 * t347;
t375 = t342 * t353;
t373 = qJD(1) * t341;
t324 = (-pkin(2) * t353 - pkin(9) * t347) * t373;
t337 = qJD(1) * t342 + qJD(2);
t335 = t337 ^ 2;
t336 = qJDD(1) * t342 + qJDD(2);
t372 = qJD(1) * t353;
t348 = sin(qJ(1));
t354 = cos(qJ(1));
t332 = t348 * g(1) - g(2) * t354;
t355 = qJD(1) ^ 2;
t321 = qJDD(1) * pkin(1) + t355 * t380 + t332;
t333 = -g(1) * t354 - g(2) * t348;
t371 = qJDD(1) * t341;
t322 = -pkin(1) * t355 + pkin(8) * t371 + t333;
t374 = t321 * t376 + t353 * t322;
t281 = -pkin(2) * t335 + pkin(9) * t336 + (-g(3) * t347 + t324 * t372) * t341 + t374;
t325 = (qJD(2) * t372 + qJDD(1) * t347) * t341;
t370 = t347 * t373;
t326 = -qJD(2) * t370 + t353 * t371;
t282 = -pkin(2) * t326 - pkin(9) * t325 - t379 + (-t321 + (pkin(2) * t347 - pkin(9) * t353) * t337 * qJD(1)) * t341;
t346 = sin(qJ(3));
t352 = cos(qJ(3));
t256 = -t346 * t281 + t352 * t282;
t313 = t337 * t352 - t346 * t370;
t293 = qJD(3) * t313 + t325 * t352 + t336 * t346;
t314 = t337 * t346 + t352 * t370;
t318 = qJDD(3) - t326;
t369 = t341 * t372;
t331 = qJD(3) - t369;
t241 = (t313 * t331 - t293) * pkin(10) + (t313 * t314 + t318) * pkin(3) + t256;
t257 = t352 * t281 + t346 * t282;
t292 = -qJD(3) * t314 - t325 * t346 + t336 * t352;
t302 = pkin(3) * t331 - pkin(10) * t314;
t312 = t313 ^ 2;
t244 = -pkin(3) * t312 + pkin(10) * t292 - t302 * t331 + t257;
t345 = sin(qJ(4));
t351 = cos(qJ(4));
t221 = t351 * t241 - t244 * t345;
t297 = t313 * t351 - t314 * t345;
t262 = qJD(4) * t297 + t292 * t345 + t293 * t351;
t298 = t313 * t345 + t314 * t351;
t317 = qJDD(4) + t318;
t329 = qJD(4) + t331;
t217 = (t297 * t329 - t262) * pkin(11) + (t297 * t298 + t317) * pkin(4) + t221;
t222 = t345 * t241 + t351 * t244;
t261 = -qJD(4) * t298 + t292 * t351 - t293 * t345;
t285 = pkin(4) * t329 - pkin(11) * t298;
t296 = t297 ^ 2;
t219 = -pkin(4) * t296 + pkin(11) * t261 - t285 * t329 + t222;
t344 = sin(qJ(5));
t350 = cos(qJ(5));
t214 = t344 * t217 + t350 * t219;
t275 = t297 * t344 + t298 * t350;
t233 = -qJD(5) * t275 + t261 * t350 - t262 * t344;
t274 = t297 * t350 - t298 * t344;
t254 = -mrSges(6,1) * t274 + mrSges(6,2) * t275;
t328 = qJD(5) + t329;
t267 = mrSges(6,1) * t328 - mrSges(6,3) * t275;
t311 = qJDD(5) + t317;
t255 = -pkin(5) * t274 - pkin(12) * t275;
t327 = t328 ^ 2;
t211 = -pkin(5) * t327 + pkin(12) * t311 + t255 * t274 + t214;
t294 = -g(3) * t377 + t321 * t375 - t347 * t322;
t280 = -pkin(2) * t336 - pkin(9) * t335 + t324 * t370 - t294;
t248 = -pkin(3) * t292 - pkin(10) * t312 + t314 * t302 + t280;
t227 = -pkin(4) * t261 - pkin(11) * t296 + t298 * t285 + t248;
t234 = qJD(5) * t274 + t261 * t344 + t262 * t350;
t215 = (t275 * t328 - t233) * pkin(5) + (-t274 * t328 - t234) * pkin(12) + t227;
t343 = sin(qJ(6));
t349 = cos(qJ(6));
t208 = -t211 * t343 + t215 * t349;
t264 = -t275 * t343 + t328 * t349;
t225 = qJD(6) * t264 + t234 * t349 + t311 * t343;
t232 = qJDD(6) - t233;
t265 = t275 * t349 + t328 * t343;
t245 = -mrSges(7,1) * t264 + mrSges(7,2) * t265;
t273 = qJD(6) - t274;
t246 = -mrSges(7,2) * t273 + mrSges(7,3) * t264;
t204 = m(7) * t208 + mrSges(7,1) * t232 - mrSges(7,3) * t225 - t245 * t265 + t246 * t273;
t209 = t211 * t349 + t215 * t343;
t224 = -qJD(6) * t265 - t234 * t343 + t311 * t349;
t247 = mrSges(7,1) * t273 - mrSges(7,3) * t265;
t205 = m(7) * t209 - mrSges(7,2) * t232 + mrSges(7,3) * t224 + t245 * t264 - t247 * t273;
t365 = -t204 * t343 + t349 * t205;
t191 = m(6) * t214 - mrSges(6,2) * t311 + mrSges(6,3) * t233 + t254 * t274 - t267 * t328 + t365;
t213 = t217 * t350 - t219 * t344;
t266 = -mrSges(6,2) * t328 + mrSges(6,3) * t274;
t210 = -pkin(5) * t311 - pkin(12) * t327 + t255 * t275 - t213;
t363 = -m(7) * t210 + t224 * mrSges(7,1) - mrSges(7,2) * t225 + t264 * t246 - t247 * t265;
t200 = m(6) * t213 + mrSges(6,1) * t311 - mrSges(6,3) * t234 - t254 * t275 + t266 * t328 + t363;
t186 = t344 * t191 + t350 * t200;
t276 = -mrSges(5,1) * t297 + mrSges(5,2) * t298;
t283 = -mrSges(5,2) * t329 + mrSges(5,3) * t297;
t183 = m(5) * t221 + mrSges(5,1) * t317 - mrSges(5,3) * t262 - t276 * t298 + t283 * t329 + t186;
t284 = mrSges(5,1) * t329 - mrSges(5,3) * t298;
t366 = t350 * t191 - t200 * t344;
t184 = m(5) * t222 - mrSges(5,2) * t317 + mrSges(5,3) * t261 + t276 * t297 - t284 * t329 + t366;
t177 = t351 * t183 + t345 * t184;
t299 = -mrSges(4,1) * t313 + mrSges(4,2) * t314;
t300 = -mrSges(4,2) * t331 + mrSges(4,3) * t313;
t175 = m(4) * t256 + mrSges(4,1) * t318 - mrSges(4,3) * t293 - t299 * t314 + t300 * t331 + t177;
t301 = mrSges(4,1) * t331 - mrSges(4,3) * t314;
t367 = -t183 * t345 + t351 * t184;
t176 = m(4) * t257 - mrSges(4,2) * t318 + mrSges(4,3) * t292 + t299 * t313 - t301 * t331 + t367;
t170 = t352 * t175 + t346 * t176;
t193 = t349 * t204 + t343 * t205;
t295 = -g(3) * t378 + t374;
t319 = mrSges(3,1) * t337 - mrSges(3,3) * t370;
t323 = (-mrSges(3,1) * t353 + mrSges(3,2) * t347) * t373;
t368 = -t175 * t346 + t352 * t176;
t168 = m(3) * t295 - mrSges(3,2) * t336 + mrSges(3,3) * t326 - t319 * t337 + t323 * t369 + t368;
t320 = -mrSges(3,2) * t337 + mrSges(3,3) * t369;
t364 = m(6) * t227 - t233 * mrSges(6,1) + t234 * mrSges(6,2) - t274 * t266 + t275 * t267 + t193;
t360 = m(5) * t248 - t261 * mrSges(5,1) + t262 * mrSges(5,2) - t297 * t283 + t298 * t284 + t364;
t357 = -m(4) * t280 + t292 * mrSges(4,1) - t293 * mrSges(4,2) + t313 * t300 - t314 * t301 - t360;
t188 = m(3) * t294 + t336 * mrSges(3,1) - t325 * mrSges(3,3) + t337 * t320 - t323 * t370 + t357;
t164 = t353 * t168 - t188 * t347;
t306 = -t321 * t341 - t379;
t169 = m(3) * t306 - mrSges(3,1) * t326 + mrSges(3,2) * t325 + (t319 * t347 - t320 * t353) * t373 + t170;
t161 = t168 * t376 - t169 * t341 + t188 * t375;
t235 = Ifges(7,5) * t265 + Ifges(7,6) * t264 + Ifges(7,3) * t273;
t237 = Ifges(7,1) * t265 + Ifges(7,4) * t264 + Ifges(7,5) * t273;
t197 = -mrSges(7,1) * t210 + mrSges(7,3) * t209 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t232 - t235 * t265 + t237 * t273;
t236 = Ifges(7,4) * t265 + Ifges(7,2) * t264 + Ifges(7,6) * t273;
t198 = mrSges(7,2) * t210 - mrSges(7,3) * t208 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t232 + t235 * t264 - t236 * t273;
t249 = Ifges(6,5) * t275 + Ifges(6,6) * t274 + Ifges(6,3) * t328;
t250 = Ifges(6,4) * t275 + Ifges(6,2) * t274 + Ifges(6,6) * t328;
t178 = mrSges(6,2) * t227 - mrSges(6,3) * t213 + Ifges(6,1) * t234 + Ifges(6,4) * t233 + Ifges(6,5) * t311 - pkin(12) * t193 - t197 * t343 + t198 * t349 + t249 * t274 - t250 * t328;
t251 = Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t328;
t359 = mrSges(7,1) * t208 - mrSges(7,2) * t209 + Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t232 + t236 * t265 - t237 * t264;
t179 = -mrSges(6,1) * t227 + mrSges(6,3) * t214 + Ifges(6,4) * t234 + Ifges(6,2) * t233 + Ifges(6,6) * t311 - pkin(5) * t193 - t249 * t275 + t251 * t328 - t359;
t268 = Ifges(5,5) * t298 + Ifges(5,6) * t297 + Ifges(5,3) * t329;
t270 = Ifges(5,1) * t298 + Ifges(5,4) * t297 + Ifges(5,5) * t329;
t165 = -mrSges(5,1) * t248 + mrSges(5,3) * t222 + Ifges(5,4) * t262 + Ifges(5,2) * t261 + Ifges(5,6) * t317 - pkin(4) * t364 + pkin(11) * t366 + t344 * t178 + t350 * t179 - t298 * t268 + t329 * t270;
t269 = Ifges(5,4) * t298 + Ifges(5,2) * t297 + Ifges(5,6) * t329;
t171 = mrSges(5,2) * t248 - mrSges(5,3) * t221 + Ifges(5,1) * t262 + Ifges(5,4) * t261 + Ifges(5,5) * t317 - pkin(11) * t186 + t178 * t350 - t179 * t344 + t268 * t297 - t269 * t329;
t286 = Ifges(4,5) * t314 + Ifges(4,6) * t313 + Ifges(4,3) * t331;
t288 = Ifges(4,1) * t314 + Ifges(4,4) * t313 + Ifges(4,5) * t331;
t154 = -mrSges(4,1) * t280 + mrSges(4,3) * t257 + Ifges(4,4) * t293 + Ifges(4,2) * t292 + Ifges(4,6) * t318 - pkin(3) * t360 + pkin(10) * t367 + t351 * t165 + t345 * t171 - t314 * t286 + t331 * t288;
t287 = Ifges(4,4) * t314 + Ifges(4,2) * t313 + Ifges(4,6) * t331;
t155 = mrSges(4,2) * t280 - mrSges(4,3) * t256 + Ifges(4,1) * t293 + Ifges(4,4) * t292 + Ifges(4,5) * t318 - pkin(10) * t177 - t165 * t345 + t171 * t351 + t286 * t313 - t287 * t331;
t304 = Ifges(3,6) * t337 + (Ifges(3,4) * t347 + Ifges(3,2) * t353) * t373;
t305 = Ifges(3,5) * t337 + (Ifges(3,1) * t347 + Ifges(3,4) * t353) * t373;
t151 = Ifges(3,5) * t325 + Ifges(3,6) * t326 + Ifges(3,3) * t336 + mrSges(3,1) * t294 - mrSges(3,2) * t295 + t346 * t155 + t352 * t154 + pkin(2) * t357 + pkin(9) * t368 + (t304 * t347 - t305 * t353) * t373;
t303 = Ifges(3,3) * t337 + (Ifges(3,5) * t347 + Ifges(3,6) * t353) * t373;
t153 = mrSges(3,2) * t306 - mrSges(3,3) * t294 + Ifges(3,1) * t325 + Ifges(3,4) * t326 + Ifges(3,5) * t336 - pkin(9) * t170 - t154 * t346 + t155 * t352 + t303 * t369 - t304 * t337;
t361 = -mrSges(6,1) * t213 + mrSges(6,2) * t214 - Ifges(6,5) * t234 - Ifges(6,6) * t233 - Ifges(6,3) * t311 - pkin(5) * t363 - pkin(12) * t365 - t349 * t197 - t343 * t198 - t275 * t250 + t274 * t251;
t358 = -mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t262 - Ifges(5,6) * t261 - Ifges(5,3) * t317 - pkin(4) * t186 - t298 * t269 + t297 * t270 + t361;
t356 = mrSges(4,1) * t256 - mrSges(4,2) * t257 + Ifges(4,5) * t293 + Ifges(4,6) * t292 + Ifges(4,3) * t318 + pkin(3) * t177 + t314 * t287 - t313 * t288 - t358;
t157 = -mrSges(3,1) * t306 + mrSges(3,3) * t295 + Ifges(3,4) * t325 + Ifges(3,2) * t326 + Ifges(3,6) * t336 - pkin(2) * t170 - t303 * t370 + t337 * t305 - t356;
t362 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t161 + t342 * t151 + t153 * t378 + t157 * t377 + t164 * t380;
t162 = m(2) * t333 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t164;
t160 = t342 * t169 + (t168 * t347 + t188 * t353) * t341;
t158 = m(2) * t332 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t161;
t149 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t355 + t153 * t353 - t157 * t347 + (-t160 * t341 - t161 * t342) * pkin(8);
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + Ifges(2,5) * t355 + Ifges(2,6) * qJDD(1) - pkin(1) * t160 - t151 * t341 + (pkin(8) * t164 + t153 * t347 + t157 * t353) * t342;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t149 - t348 * t148 - pkin(7) * (t158 * t354 + t162 * t348), t149, t153, t155, t171, t178, t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t348 * t149 + t354 * t148 + pkin(7) * (-t158 * t348 + t162 * t354), t148, t157, t154, t165, t179, t197; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t362, t362, t151, t356, -t358, -t361, t359;];
m_new  = t1;
