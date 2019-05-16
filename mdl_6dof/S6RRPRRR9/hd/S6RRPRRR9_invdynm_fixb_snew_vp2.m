% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-06 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 23:01:09
% EndTime: 2019-05-06 23:03:05
% DurationCPUTime: 89.31s
% Computational Cost: add. (1619287->399), mult. (3685749->517), div. (0->0), fcn. (3031439->14), ass. (0->165)
t342 = sin(pkin(6));
t380 = pkin(8) * t342;
t344 = cos(pkin(6));
t379 = t344 * g(3);
t348 = sin(qJ(2));
t378 = t342 * t348;
t353 = cos(qJ(2));
t377 = t342 * t353;
t376 = t344 * t348;
t375 = t344 * t353;
t373 = qJD(1) * t342;
t324 = (-pkin(2) * t353 - qJ(3) * t348) * t373;
t337 = qJD(1) * t344 + qJD(2);
t335 = t337 ^ 2;
t336 = qJDD(1) * t344 + qJDD(2);
t372 = qJD(1) * t353;
t349 = sin(qJ(1));
t354 = cos(qJ(1));
t332 = t349 * g(1) - g(2) * t354;
t355 = qJD(1) ^ 2;
t322 = qJDD(1) * pkin(1) + t355 * t380 + t332;
t333 = -g(1) * t354 - g(2) * t349;
t371 = qJDD(1) * t342;
t323 = -pkin(1) * t355 + pkin(8) * t371 + t333;
t374 = t322 * t376 + t353 * t323;
t281 = -t335 * pkin(2) + t336 * qJ(3) + (-g(3) * t348 + t324 * t372) * t342 + t374;
t326 = (qJD(2) * t372 + qJDD(1) * t348) * t342;
t370 = t348 * t373;
t327 = -qJD(2) * t370 + t353 * t371;
t282 = -t327 * pkin(2) - t379 - t326 * qJ(3) + (-t322 + (pkin(2) * t348 - qJ(3) * t353) * t337 * qJD(1)) * t342;
t341 = sin(pkin(12));
t343 = cos(pkin(12));
t314 = t337 * t341 + t343 * t370;
t253 = -0.2e1 * qJD(3) * t314 - t341 * t281 + t343 * t282;
t301 = t326 * t343 + t336 * t341;
t313 = t337 * t343 - t341 * t370;
t369 = t342 * t372;
t241 = (-t313 * t369 - t301) * pkin(9) + (t313 * t314 - t327) * pkin(3) + t253;
t254 = 0.2e1 * qJD(3) * t313 + t343 * t281 + t341 * t282;
t300 = -t326 * t341 + t336 * t343;
t302 = -pkin(3) * t369 - pkin(9) * t314;
t311 = t313 ^ 2;
t244 = -pkin(3) * t311 + pkin(9) * t300 + t302 * t369 + t254;
t347 = sin(qJ(4));
t352 = cos(qJ(4));
t221 = t352 * t241 - t347 * t244;
t294 = t313 * t352 - t314 * t347;
t265 = qJD(4) * t294 + t300 * t347 + t301 * t352;
t295 = t313 * t347 + t314 * t352;
t319 = qJDD(4) - t327;
t331 = qJD(4) - t369;
t217 = (t294 * t331 - t265) * pkin(10) + (t294 * t295 + t319) * pkin(4) + t221;
t222 = t347 * t241 + t352 * t244;
t264 = -qJD(4) * t295 + t300 * t352 - t301 * t347;
t285 = pkin(4) * t331 - pkin(10) * t295;
t293 = t294 ^ 2;
t219 = -pkin(4) * t293 + pkin(10) * t264 - t285 * t331 + t222;
t346 = sin(qJ(5));
t351 = cos(qJ(5));
t214 = t346 * t217 + t351 * t219;
t275 = t294 * t346 + t295 * t351;
t237 = -qJD(5) * t275 + t264 * t351 - t265 * t346;
t274 = t294 * t351 - t295 * t346;
t255 = -mrSges(6,1) * t274 + mrSges(6,2) * t275;
t329 = qJD(5) + t331;
t267 = mrSges(6,1) * t329 - mrSges(6,3) * t275;
t318 = qJDD(5) + t319;
t256 = -pkin(5) * t274 - pkin(11) * t275;
t328 = t329 ^ 2;
t211 = -pkin(5) * t328 + pkin(11) * t318 + t256 * t274 + t214;
t291 = -g(3) * t377 + t322 * t375 - t348 * t323;
t280 = -t336 * pkin(2) - t335 * qJ(3) + t324 * t370 + qJDD(3) - t291;
t257 = -t300 * pkin(3) - t311 * pkin(9) + t314 * t302 + t280;
t227 = -t264 * pkin(4) - t293 * pkin(10) + t295 * t285 + t257;
t238 = qJD(5) * t274 + t264 * t346 + t265 * t351;
t215 = (-t274 * t329 - t238) * pkin(11) + (t275 * t329 - t237) * pkin(5) + t227;
t345 = sin(qJ(6));
t350 = cos(qJ(6));
t208 = -t211 * t345 + t215 * t350;
t259 = -t275 * t345 + t329 * t350;
t225 = qJD(6) * t259 + t238 * t350 + t318 * t345;
t232 = qJDD(6) - t237;
t260 = t275 * t350 + t329 * t345;
t245 = -mrSges(7,1) * t259 + mrSges(7,2) * t260;
t273 = qJD(6) - t274;
t246 = -mrSges(7,2) * t273 + mrSges(7,3) * t259;
t204 = m(7) * t208 + mrSges(7,1) * t232 - mrSges(7,3) * t225 - t245 * t260 + t246 * t273;
t209 = t211 * t350 + t215 * t345;
t224 = -qJD(6) * t260 - t238 * t345 + t318 * t350;
t247 = mrSges(7,1) * t273 - mrSges(7,3) * t260;
t205 = m(7) * t209 - mrSges(7,2) * t232 + mrSges(7,3) * t224 + t245 * t259 - t247 * t273;
t365 = -t204 * t345 + t350 * t205;
t191 = m(6) * t214 - mrSges(6,2) * t318 + mrSges(6,3) * t237 + t255 * t274 - t267 * t329 + t365;
t213 = t217 * t351 - t219 * t346;
t266 = -mrSges(6,2) * t329 + mrSges(6,3) * t274;
t210 = -pkin(5) * t318 - pkin(11) * t328 + t256 * t275 - t213;
t363 = -m(7) * t210 + t224 * mrSges(7,1) - mrSges(7,2) * t225 + t259 * t246 - t247 * t260;
t200 = m(6) * t213 + mrSges(6,1) * t318 - mrSges(6,3) * t238 - t255 * t275 + t266 * t329 + t363;
t186 = t346 * t191 + t351 * t200;
t276 = -mrSges(5,1) * t294 + mrSges(5,2) * t295;
t283 = -mrSges(5,2) * t331 + mrSges(5,3) * t294;
t183 = m(5) * t221 + mrSges(5,1) * t319 - mrSges(5,3) * t265 - t276 * t295 + t283 * t331 + t186;
t284 = mrSges(5,1) * t331 - mrSges(5,3) * t295;
t366 = t351 * t191 - t200 * t346;
t184 = m(5) * t222 - mrSges(5,2) * t319 + mrSges(5,3) * t264 + t276 * t294 - t284 * t331 + t366;
t177 = t352 * t183 + t347 * t184;
t296 = -mrSges(4,1) * t313 + mrSges(4,2) * t314;
t298 = mrSges(4,2) * t369 + mrSges(4,3) * t313;
t175 = m(4) * t253 - mrSges(4,1) * t327 - mrSges(4,3) * t301 - t296 * t314 - t298 * t369 + t177;
t299 = -mrSges(4,1) * t369 - mrSges(4,3) * t314;
t367 = -t183 * t347 + t352 * t184;
t176 = m(4) * t254 + mrSges(4,2) * t327 + mrSges(4,3) * t300 + t296 * t313 + t299 * t369 + t367;
t170 = t343 * t175 + t341 * t176;
t193 = t350 * t204 + t345 * t205;
t292 = -g(3) * t378 + t374;
t320 = mrSges(3,1) * t337 - mrSges(3,3) * t370;
t325 = (-mrSges(3,1) * t353 + mrSges(3,2) * t348) * t373;
t368 = -t175 * t341 + t343 * t176;
t168 = m(3) * t292 - mrSges(3,2) * t336 + mrSges(3,3) * t327 - t320 * t337 + t325 * t369 + t368;
t321 = -mrSges(3,2) * t337 + mrSges(3,3) * t369;
t364 = m(6) * t227 - t237 * mrSges(6,1) + t238 * mrSges(6,2) - t274 * t266 + t275 * t267 + t193;
t360 = m(5) * t257 - t264 * mrSges(5,1) + t265 * mrSges(5,2) - t294 * t283 + t295 * t284 + t364;
t357 = -m(4) * t280 + t300 * mrSges(4,1) - t301 * mrSges(4,2) + t313 * t298 - t314 * t299 - t360;
t188 = m(3) * t291 + t336 * mrSges(3,1) - t326 * mrSges(3,3) + t337 * t321 - t325 * t370 + t357;
t164 = t353 * t168 - t188 * t348;
t306 = -t342 * t322 - t379;
t169 = m(3) * t306 - t327 * mrSges(3,1) + t326 * mrSges(3,2) + (t320 * t348 - t321 * t353) * t373 + t170;
t161 = t168 * t376 - t169 * t342 + t188 * t375;
t233 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t273;
t235 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t273;
t197 = -mrSges(7,1) * t210 + mrSges(7,3) * t209 + Ifges(7,4) * t225 + Ifges(7,2) * t224 + Ifges(7,6) * t232 - t233 * t260 + t235 * t273;
t234 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t273;
t198 = mrSges(7,2) * t210 - mrSges(7,3) * t208 + Ifges(7,1) * t225 + Ifges(7,4) * t224 + Ifges(7,5) * t232 + t233 * t259 - t234 * t273;
t248 = Ifges(6,5) * t275 + Ifges(6,6) * t274 + Ifges(6,3) * t329;
t249 = Ifges(6,4) * t275 + Ifges(6,2) * t274 + Ifges(6,6) * t329;
t178 = mrSges(6,2) * t227 - mrSges(6,3) * t213 + Ifges(6,1) * t238 + Ifges(6,4) * t237 + Ifges(6,5) * t318 - pkin(11) * t193 - t197 * t345 + t198 * t350 + t248 * t274 - t249 * t329;
t250 = Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t329;
t359 = mrSges(7,1) * t208 - mrSges(7,2) * t209 + Ifges(7,5) * t225 + Ifges(7,6) * t224 + Ifges(7,3) * t232 + t234 * t260 - t235 * t259;
t179 = -mrSges(6,1) * t227 + mrSges(6,3) * t214 + Ifges(6,4) * t238 + Ifges(6,2) * t237 + Ifges(6,6) * t318 - pkin(5) * t193 - t248 * t275 + t250 * t329 - t359;
t268 = Ifges(5,5) * t295 + Ifges(5,6) * t294 + Ifges(5,3) * t331;
t270 = Ifges(5,1) * t295 + Ifges(5,4) * t294 + Ifges(5,5) * t331;
t165 = -mrSges(5,1) * t257 + mrSges(5,3) * t222 + Ifges(5,4) * t265 + Ifges(5,2) * t264 + Ifges(5,6) * t319 - pkin(4) * t364 + pkin(10) * t366 + t346 * t178 + t351 * t179 - t295 * t268 + t331 * t270;
t269 = Ifges(5,4) * t295 + Ifges(5,2) * t294 + Ifges(5,6) * t331;
t171 = mrSges(5,2) * t257 - mrSges(5,3) * t221 + Ifges(5,1) * t265 + Ifges(5,4) * t264 + Ifges(5,5) * t319 - pkin(10) * t186 + t178 * t351 - t179 * t346 + t268 * t294 - t269 * t331;
t286 = Ifges(4,5) * t314 + Ifges(4,6) * t313 - Ifges(4,3) * t369;
t288 = Ifges(4,1) * t314 + Ifges(4,4) * t313 - Ifges(4,5) * t369;
t154 = -mrSges(4,1) * t280 + mrSges(4,3) * t254 + Ifges(4,4) * t301 + Ifges(4,2) * t300 - Ifges(4,6) * t327 - pkin(3) * t360 + pkin(9) * t367 + t352 * t165 + t347 * t171 - t314 * t286 - t288 * t369;
t287 = Ifges(4,4) * t314 + Ifges(4,2) * t313 - Ifges(4,6) * t369;
t155 = mrSges(4,2) * t280 - mrSges(4,3) * t253 + Ifges(4,1) * t301 + Ifges(4,4) * t300 - Ifges(4,5) * t327 - pkin(9) * t177 - t165 * t347 + t171 * t352 + t286 * t313 + t287 * t369;
t304 = Ifges(3,6) * t337 + (Ifges(3,4) * t348 + Ifges(3,2) * t353) * t373;
t305 = Ifges(3,5) * t337 + (Ifges(3,1) * t348 + Ifges(3,4) * t353) * t373;
t151 = Ifges(3,5) * t326 + Ifges(3,6) * t327 + Ifges(3,3) * t336 + mrSges(3,1) * t291 - mrSges(3,2) * t292 + t341 * t155 + t343 * t154 + pkin(2) * t357 + qJ(3) * t368 + (t304 * t348 - t305 * t353) * t373;
t303 = Ifges(3,3) * t337 + (Ifges(3,5) * t348 + Ifges(3,6) * t353) * t373;
t153 = mrSges(3,2) * t306 - mrSges(3,3) * t291 + Ifges(3,1) * t326 + Ifges(3,4) * t327 + Ifges(3,5) * t336 - qJ(3) * t170 - t154 * t341 + t155 * t343 + t303 * t369 - t304 * t337;
t361 = -mrSges(6,1) * t213 + mrSges(6,2) * t214 - Ifges(6,5) * t238 - Ifges(6,6) * t237 - Ifges(6,3) * t318 - pkin(5) * t363 - pkin(11) * t365 - t350 * t197 - t345 * t198 - t275 * t249 + t274 * t250;
t358 = -mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t265 - Ifges(5,6) * t264 - Ifges(5,3) * t319 - pkin(4) * t186 - t295 * t269 + t294 * t270 + t361;
t356 = mrSges(4,1) * t253 - mrSges(4,2) * t254 + Ifges(4,5) * t301 + Ifges(4,6) * t300 + pkin(3) * t177 + t314 * t287 - t313 * t288 - t358;
t157 = (Ifges(3,2) + Ifges(4,3)) * t327 - pkin(2) * t170 - t356 + Ifges(3,6) * t336 + t337 * t305 + Ifges(3,4) * t326 - mrSges(3,1) * t306 + mrSges(3,3) * t292 - t303 * t370;
t362 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t161 + t344 * t151 + t153 * t378 + t157 * t377 + t164 * t380;
t162 = m(2) * t333 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t164;
t160 = t344 * t169 + (t168 * t348 + t188 * t353) * t342;
t158 = m(2) * t332 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t161;
t149 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) + t353 * t153 - t348 * t157 + (-t160 * t342 - t161 * t344) * pkin(8);
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t160 - t342 * t151 + (pkin(8) * t164 + t153 * t348 + t157 * t353) * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t149 - t349 * t148 - pkin(7) * (t158 * t354 + t162 * t349), t149, t153, t155, t171, t178, t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t349 * t149 + t354 * t148 + pkin(7) * (-t158 * t349 + t162 * t354), t148, t157, t154, t165, t179, t197; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t362, t362, t151, -Ifges(4,3) * t327 + t356, -t358, -t361, t359;];
m_new  = t1;
