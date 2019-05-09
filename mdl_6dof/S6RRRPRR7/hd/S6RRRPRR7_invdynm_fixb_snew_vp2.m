% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR7
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
% Datum: 2019-05-07 11:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 11:28:07
% EndTime: 2019-05-07 11:31:23
% DurationCPUTime: 95.15s
% Computational Cost: add. (1691625->398), mult. (3749540->517), div. (0->0), fcn. (3068822->14), ass. (0->165)
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
t325 = (-pkin(2) * t353 - pkin(9) * t348) * t373;
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
t281 = -t335 * pkin(2) + t336 * pkin(9) + (-g(3) * t348 + t325 * t372) * t342 + t374;
t326 = (qJD(2) * t372 + qJDD(1) * t348) * t342;
t370 = t348 * t373;
t327 = -qJD(2) * t370 + t353 * t371;
t282 = -t327 * pkin(2) - t326 * pkin(9) - t379 + (-t322 + (pkin(2) * t348 - pkin(9) * t353) * t337 * qJD(1)) * t342;
t347 = sin(qJ(3));
t352 = cos(qJ(3));
t256 = -t347 * t281 + t352 * t282;
t313 = t337 * t352 - t347 * t370;
t293 = qJD(3) * t313 + t326 * t352 + t336 * t347;
t314 = t337 * t347 + t352 * t370;
t319 = qJDD(3) - t327;
t369 = t342 * t372;
t331 = qJD(3) - t369;
t241 = (t313 * t331 - t293) * qJ(4) + (t313 * t314 + t319) * pkin(3) + t256;
t257 = t352 * t281 + t347 * t282;
t292 = -qJD(3) * t314 - t326 * t347 + t336 * t352;
t303 = pkin(3) * t331 - qJ(4) * t314;
t312 = t313 ^ 2;
t244 = -pkin(3) * t312 + qJ(4) * t292 - t303 * t331 + t257;
t341 = sin(pkin(12));
t343 = cos(pkin(12));
t300 = t313 * t341 + t314 * t343;
t221 = -0.2e1 * qJD(4) * t300 + t343 * t241 - t341 * t244;
t267 = t292 * t341 + t293 * t343;
t299 = t313 * t343 - t314 * t341;
t217 = (t299 * t331 - t267) * pkin(10) + (t299 * t300 + t319) * pkin(4) + t221;
t222 = 0.2e1 * qJD(4) * t299 + t341 * t241 + t343 * t244;
t266 = t292 * t343 - t293 * t341;
t285 = pkin(4) * t331 - pkin(10) * t300;
t298 = t299 ^ 2;
t219 = -pkin(4) * t298 + pkin(10) * t266 - t285 * t331 + t222;
t346 = sin(qJ(5));
t351 = cos(qJ(5));
t214 = t346 * t217 + t351 * t219;
t275 = t299 * t346 + t300 * t351;
t239 = -qJD(5) * t275 + t266 * t351 - t267 * t346;
t274 = t299 * t351 - t300 * t346;
t254 = -mrSges(6,1) * t274 + mrSges(6,2) * t275;
t329 = qJD(5) + t331;
t262 = mrSges(6,1) * t329 - mrSges(6,3) * t275;
t318 = qJDD(5) + t319;
t255 = -pkin(5) * t274 - pkin(11) * t275;
t328 = t329 ^ 2;
t211 = -pkin(5) * t328 + pkin(11) * t318 + t255 * t274 + t214;
t294 = -g(3) * t377 + t322 * t375 - t348 * t323;
t280 = -t336 * pkin(2) - t335 * pkin(9) + t325 * t370 - t294;
t248 = -t292 * pkin(3) - t312 * qJ(4) + t314 * t303 + qJDD(4) + t280;
t224 = -t266 * pkin(4) - t298 * pkin(10) + t300 * t285 + t248;
t240 = qJD(5) * t274 + t266 * t346 + t267 * t351;
t215 = (t275 * t329 - t239) * pkin(5) + (-t274 * t329 - t240) * pkin(11) + t224;
t345 = sin(qJ(6));
t350 = cos(qJ(6));
t208 = -t211 * t345 + t215 * t350;
t259 = -t275 * t345 + t329 * t350;
t227 = qJD(6) * t259 + t240 * t350 + t318 * t345;
t238 = qJDD(6) - t239;
t260 = t275 * t350 + t329 * t345;
t245 = -mrSges(7,1) * t259 + mrSges(7,2) * t260;
t273 = qJD(6) - t274;
t246 = -mrSges(7,2) * t273 + mrSges(7,3) * t259;
t204 = m(7) * t208 + mrSges(7,1) * t238 - mrSges(7,3) * t227 - t245 * t260 + t246 * t273;
t209 = t211 * t350 + t215 * t345;
t226 = -qJD(6) * t260 - t240 * t345 + t318 * t350;
t247 = mrSges(7,1) * t273 - mrSges(7,3) * t260;
t205 = m(7) * t209 - mrSges(7,2) * t238 + mrSges(7,3) * t226 + t245 * t259 - t247 * t273;
t365 = -t204 * t345 + t350 * t205;
t189 = m(6) * t214 - mrSges(6,2) * t318 + mrSges(6,3) * t239 + t254 * t274 - t262 * t329 + t365;
t213 = t217 * t351 - t219 * t346;
t261 = -mrSges(6,2) * t329 + mrSges(6,3) * t274;
t210 = -pkin(5) * t318 - pkin(11) * t328 + t255 * t275 - t213;
t363 = -m(7) * t210 + t226 * mrSges(7,1) - mrSges(7,2) * t227 + t259 * t246 - t247 * t260;
t200 = m(6) * t213 + mrSges(6,1) * t318 - mrSges(6,3) * t240 - t254 * t275 + t261 * t329 + t363;
t186 = t346 * t189 + t351 * t200;
t276 = -mrSges(5,1) * t299 + mrSges(5,2) * t300;
t283 = -mrSges(5,2) * t331 + mrSges(5,3) * t299;
t183 = m(5) * t221 + mrSges(5,1) * t319 - mrSges(5,3) * t267 - t276 * t300 + t283 * t331 + t186;
t284 = mrSges(5,1) * t331 - mrSges(5,3) * t300;
t366 = t351 * t189 - t200 * t346;
t184 = m(5) * t222 - mrSges(5,2) * t319 + mrSges(5,3) * t266 + t276 * t299 - t284 * t331 + t366;
t177 = t343 * t183 + t341 * t184;
t301 = -mrSges(4,1) * t313 + mrSges(4,2) * t314;
t302 = -mrSges(4,2) * t331 + mrSges(4,3) * t313;
t175 = m(4) * t256 + mrSges(4,1) * t319 - mrSges(4,3) * t293 - t301 * t314 + t302 * t331 + t177;
t304 = mrSges(4,1) * t331 - mrSges(4,3) * t314;
t367 = -t183 * t341 + t343 * t184;
t176 = m(4) * t257 - mrSges(4,2) * t319 + mrSges(4,3) * t292 + t301 * t313 - t304 * t331 + t367;
t170 = t352 * t175 + t347 * t176;
t193 = t350 * t204 + t345 * t205;
t295 = -g(3) * t378 + t374;
t320 = mrSges(3,1) * t337 - mrSges(3,3) * t370;
t324 = (-mrSges(3,1) * t353 + mrSges(3,2) * t348) * t373;
t368 = -t175 * t347 + t352 * t176;
t168 = m(3) * t295 - mrSges(3,2) * t336 + mrSges(3,3) * t327 - t320 * t337 + t324 * t369 + t368;
t321 = -mrSges(3,2) * t337 + mrSges(3,3) * t369;
t364 = m(6) * t224 - t239 * mrSges(6,1) + t240 * mrSges(6,2) - t274 * t261 + t275 * t262 + t193;
t360 = m(5) * t248 - t266 * mrSges(5,1) + t267 * mrSges(5,2) - t299 * t283 + t300 * t284 + t364;
t357 = -m(4) * t280 + t292 * mrSges(4,1) - t293 * mrSges(4,2) + t313 * t302 - t314 * t304 - t360;
t191 = m(3) * t294 + t336 * mrSges(3,1) - t326 * mrSges(3,3) + t337 * t321 - t324 * t370 + t357;
t164 = t353 * t168 - t191 * t348;
t308 = -t342 * t322 - t379;
t169 = m(3) * t308 - t327 * mrSges(3,1) + t326 * mrSges(3,2) + (t320 * t348 - t321 * t353) * t373 + t170;
t161 = t168 * t376 - t169 * t342 + t191 * t375;
t228 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t273;
t230 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t273;
t197 = -mrSges(7,1) * t210 + mrSges(7,3) * t209 + Ifges(7,4) * t227 + Ifges(7,2) * t226 + Ifges(7,6) * t238 - t228 * t260 + t230 * t273;
t229 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t273;
t198 = mrSges(7,2) * t210 - mrSges(7,3) * t208 + Ifges(7,1) * t227 + Ifges(7,4) * t226 + Ifges(7,5) * t238 + t228 * t259 - t229 * t273;
t249 = Ifges(6,5) * t275 + Ifges(6,6) * t274 + Ifges(6,3) * t329;
t250 = Ifges(6,4) * t275 + Ifges(6,2) * t274 + Ifges(6,6) * t329;
t178 = mrSges(6,2) * t224 - mrSges(6,3) * t213 + Ifges(6,1) * t240 + Ifges(6,4) * t239 + Ifges(6,5) * t318 - pkin(11) * t193 - t197 * t345 + t198 * t350 + t249 * t274 - t250 * t329;
t251 = Ifges(6,1) * t275 + Ifges(6,4) * t274 + Ifges(6,5) * t329;
t359 = mrSges(7,1) * t208 - mrSges(7,2) * t209 + Ifges(7,5) * t227 + Ifges(7,6) * t226 + Ifges(7,3) * t238 + t229 * t260 - t230 * t259;
t179 = -mrSges(6,1) * t224 + mrSges(6,3) * t214 + Ifges(6,4) * t240 + Ifges(6,2) * t239 + Ifges(6,6) * t318 - pkin(5) * t193 - t249 * t275 + t251 * t329 - t359;
t268 = Ifges(5,5) * t300 + Ifges(5,6) * t299 + Ifges(5,3) * t331;
t270 = Ifges(5,1) * t300 + Ifges(5,4) * t299 + Ifges(5,5) * t331;
t165 = -mrSges(5,1) * t248 + mrSges(5,3) * t222 + Ifges(5,4) * t267 + Ifges(5,2) * t266 + Ifges(5,6) * t319 - pkin(4) * t364 + pkin(10) * t366 + t346 * t178 + t351 * t179 - t300 * t268 + t331 * t270;
t269 = Ifges(5,4) * t300 + Ifges(5,2) * t299 + Ifges(5,6) * t331;
t171 = mrSges(5,2) * t248 - mrSges(5,3) * t221 + Ifges(5,1) * t267 + Ifges(5,4) * t266 + Ifges(5,5) * t319 - pkin(10) * t186 + t178 * t351 - t179 * t346 + t268 * t299 - t269 * t331;
t286 = Ifges(4,5) * t314 + Ifges(4,6) * t313 + Ifges(4,3) * t331;
t288 = Ifges(4,1) * t314 + Ifges(4,4) * t313 + Ifges(4,5) * t331;
t154 = -mrSges(4,1) * t280 + mrSges(4,3) * t257 + Ifges(4,4) * t293 + Ifges(4,2) * t292 + Ifges(4,6) * t319 - pkin(3) * t360 + qJ(4) * t367 + t343 * t165 + t341 * t171 - t314 * t286 + t331 * t288;
t287 = Ifges(4,4) * t314 + Ifges(4,2) * t313 + Ifges(4,6) * t331;
t155 = mrSges(4,2) * t280 - mrSges(4,3) * t256 + Ifges(4,1) * t293 + Ifges(4,4) * t292 + Ifges(4,5) * t319 - qJ(4) * t177 - t165 * t341 + t171 * t343 + t286 * t313 - t287 * t331;
t306 = Ifges(3,6) * t337 + (Ifges(3,4) * t348 + Ifges(3,2) * t353) * t373;
t307 = Ifges(3,5) * t337 + (Ifges(3,1) * t348 + Ifges(3,4) * t353) * t373;
t151 = Ifges(3,5) * t326 + Ifges(3,6) * t327 + Ifges(3,3) * t336 + mrSges(3,1) * t294 - mrSges(3,2) * t295 + t347 * t155 + t352 * t154 + pkin(2) * t357 + pkin(9) * t368 + (t306 * t348 - t307 * t353) * t373;
t305 = Ifges(3,3) * t337 + (Ifges(3,5) * t348 + Ifges(3,6) * t353) * t373;
t153 = mrSges(3,2) * t308 - mrSges(3,3) * t294 + Ifges(3,1) * t326 + Ifges(3,4) * t327 + Ifges(3,5) * t336 - pkin(9) * t170 - t154 * t347 + t155 * t352 + t305 * t369 - t306 * t337;
t361 = -mrSges(6,1) * t213 + mrSges(6,2) * t214 - Ifges(6,5) * t240 - Ifges(6,6) * t239 - Ifges(6,3) * t318 - pkin(5) * t363 - pkin(11) * t365 - t350 * t197 - t345 * t198 - t275 * t250 + t274 * t251;
t358 = -mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t267 - Ifges(5,6) * t266 - Ifges(5,3) * t319 - pkin(4) * t186 - t300 * t269 + t299 * t270 + t361;
t356 = mrSges(4,1) * t256 - mrSges(4,2) * t257 + Ifges(4,5) * t293 + Ifges(4,6) * t292 + Ifges(4,3) * t319 + pkin(3) * t177 + t314 * t287 - t313 * t288 - t358;
t157 = -mrSges(3,1) * t308 + mrSges(3,3) * t295 + Ifges(3,4) * t326 + Ifges(3,2) * t327 + Ifges(3,6) * t336 - pkin(2) * t170 - t305 * t370 + t337 * t307 - t356;
t362 = mrSges(2,1) * t332 - mrSges(2,2) * t333 + Ifges(2,3) * qJDD(1) + pkin(1) * t161 + t344 * t151 + t153 * t378 + t157 * t377 + t164 * t380;
t162 = m(2) * t333 - mrSges(2,1) * t355 - qJDD(1) * mrSges(2,2) + t164;
t160 = t344 * t169 + (t168 * t348 + t191 * t353) * t342;
t158 = m(2) * t332 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t355 + t161;
t149 = -mrSges(2,2) * g(3) - mrSges(2,3) * t332 + Ifges(2,5) * qJDD(1) - t355 * Ifges(2,6) + t353 * t153 - t348 * t157 + (-t160 * t342 - t161 * t344) * pkin(8);
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t333 + t355 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t160 - t342 * t151 + (pkin(8) * t164 + t153 * t348 + t157 * t353) * t344;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t354 * t149 - t349 * t148 - pkin(7) * (t158 * t354 + t162 * t349), t149, t153, t155, t171, t178, t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t349 * t149 + t354 * t148 + pkin(7) * (-t158 * t349 + t162 * t354), t148, t157, t154, t165, t179, t197; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t362, t362, t151, t356, -t358, -t361, t359;];
m_new  = t1;
