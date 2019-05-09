% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-07 21:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:10:21
% EndTime: 2019-05-07 21:13:16
% DurationCPUTime: 94.69s
% Computational Cost: add. (1754017->399), mult. (3843675->517), div. (0->0), fcn. (3153304->14), ass. (0->167)
t381 = 2 * qJD(5);
t341 = sin(pkin(6));
t380 = pkin(8) * t341;
t343 = cos(pkin(6));
t379 = t343 * g(3);
t347 = sin(qJ(2));
t378 = t341 * t347;
t352 = cos(qJ(2));
t377 = t341 * t352;
t376 = t343 * t347;
t375 = t343 * t352;
t373 = qJD(1) * t341;
t324 = (-pkin(2) * t352 - pkin(9) * t347) * t373;
t336 = qJD(1) * t343 + qJD(2);
t334 = t336 ^ 2;
t335 = qJDD(1) * t343 + qJDD(2);
t372 = qJD(1) * t352;
t348 = sin(qJ(1));
t353 = cos(qJ(1));
t331 = t348 * g(1) - g(2) * t353;
t354 = qJD(1) ^ 2;
t321 = qJDD(1) * pkin(1) + t354 * t380 + t331;
t332 = -g(1) * t353 - g(2) * t348;
t371 = qJDD(1) * t341;
t322 = -pkin(1) * t354 + pkin(8) * t371 + t332;
t374 = t321 * t376 + t352 * t322;
t282 = -t334 * pkin(2) + t335 * pkin(9) + (-g(3) * t347 + t324 * t372) * t341 + t374;
t325 = (qJD(2) * t372 + qJDD(1) * t347) * t341;
t370 = t347 * t373;
t326 = -qJD(2) * t370 + t352 * t371;
t283 = -t326 * pkin(2) - t325 * pkin(9) - t379 + (-t321 + (pkin(2) * t347 - pkin(9) * t352) * t336 * qJD(1)) * t341;
t346 = sin(qJ(3));
t351 = cos(qJ(3));
t256 = -t346 * t282 + t351 * t283;
t313 = t336 * t351 - t346 * t370;
t294 = qJD(3) * t313 + t325 * t351 + t335 * t346;
t314 = t336 * t346 + t351 * t370;
t318 = qJDD(3) - t326;
t369 = t341 * t372;
t330 = qJD(3) - t369;
t234 = (t313 * t330 - t294) * pkin(10) + (t313 * t314 + t318) * pkin(3) + t256;
t257 = t351 * t282 + t346 * t283;
t293 = -qJD(3) * t314 - t325 * t346 + t335 * t351;
t303 = pkin(3) * t330 - pkin(10) * t314;
t312 = t313 ^ 2;
t244 = -pkin(3) * t312 + pkin(10) * t293 - t303 * t330 + t257;
t345 = sin(qJ(4));
t350 = cos(qJ(4));
t221 = t350 * t234 - t345 * t244;
t298 = t313 * t350 - t314 * t345;
t262 = qJD(4) * t298 + t293 * t345 + t294 * t350;
t299 = t313 * t345 + t314 * t350;
t317 = qJDD(4) + t318;
t328 = qJD(4) + t330;
t217 = (t298 * t328 - t262) * qJ(5) + (t298 * t299 + t317) * pkin(4) + t221;
t222 = t345 * t234 + t350 * t244;
t261 = -qJD(4) * t299 + t293 * t350 - t294 * t345;
t285 = pkin(4) * t328 - qJ(5) * t299;
t297 = t298 ^ 2;
t219 = -pkin(4) * t297 + qJ(5) * t261 - t285 * t328 + t222;
t340 = sin(pkin(12));
t342 = cos(pkin(12));
t275 = t298 * t342 - t299 * t340;
t214 = t340 * t217 + t342 * t219 + t275 * t381;
t240 = t261 * t342 - t262 * t340;
t276 = t298 * t340 + t299 * t342;
t254 = -mrSges(6,1) * t275 + mrSges(6,2) * t276;
t267 = mrSges(6,1) * t328 - mrSges(6,3) * t276;
t255 = -pkin(5) * t275 - pkin(11) * t276;
t327 = t328 ^ 2;
t211 = -pkin(5) * t327 + pkin(11) * t317 + t255 * t275 + t214;
t295 = -g(3) * t377 + t321 * t375 - t347 * t322;
t281 = -t335 * pkin(2) - t334 * pkin(9) + t324 * t370 - t295;
t248 = -t293 * pkin(3) - t312 * pkin(10) + t314 * t303 + t281;
t224 = -t261 * pkin(4) - t297 * qJ(5) + t299 * t285 + qJDD(5) + t248;
t241 = t261 * t340 + t262 * t342;
t215 = t224 + (-t275 * t328 - t241) * pkin(11) + (t276 * t328 - t240) * pkin(5);
t344 = sin(qJ(6));
t349 = cos(qJ(6));
t208 = -t211 * t344 + t215 * t349;
t264 = -t276 * t344 + t328 * t349;
t227 = qJD(6) * t264 + t241 * t349 + t317 * t344;
t239 = qJDD(6) - t240;
t265 = t276 * t349 + t328 * t344;
t245 = -mrSges(7,1) * t264 + mrSges(7,2) * t265;
t274 = qJD(6) - t275;
t246 = -mrSges(7,2) * t274 + mrSges(7,3) * t264;
t204 = m(7) * t208 + mrSges(7,1) * t239 - mrSges(7,3) * t227 - t245 * t265 + t246 * t274;
t209 = t211 * t349 + t215 * t344;
t226 = -qJD(6) * t265 - t241 * t344 + t317 * t349;
t247 = mrSges(7,1) * t274 - mrSges(7,3) * t265;
t205 = m(7) * t209 - mrSges(7,2) * t239 + mrSges(7,3) * t226 + t245 * t264 - t247 * t274;
t365 = -t204 * t344 + t349 * t205;
t191 = m(6) * t214 - mrSges(6,2) * t317 + mrSges(6,3) * t240 + t254 * t275 - t267 * t328 + t365;
t364 = -t342 * t217 + t340 * t219;
t213 = -0.2e1 * qJD(5) * t276 - t364;
t266 = -mrSges(6,2) * t328 + mrSges(6,3) * t275;
t210 = -t317 * pkin(5) - t327 * pkin(11) + (t381 + t255) * t276 + t364;
t362 = -m(7) * t210 + t226 * mrSges(7,1) - mrSges(7,2) * t227 + t264 * t246 - t247 * t265;
t200 = m(6) * t213 + mrSges(6,1) * t317 - mrSges(6,3) * t241 - t254 * t276 + t266 * t328 + t362;
t186 = t340 * t191 + t342 * t200;
t277 = -mrSges(5,1) * t298 + mrSges(5,2) * t299;
t284 = -mrSges(5,2) * t328 + mrSges(5,3) * t298;
t183 = m(5) * t221 + mrSges(5,1) * t317 - mrSges(5,3) * t262 - t277 * t299 + t284 * t328 + t186;
t286 = mrSges(5,1) * t328 - mrSges(5,3) * t299;
t366 = t342 * t191 - t200 * t340;
t184 = m(5) * t222 - mrSges(5,2) * t317 + mrSges(5,3) * t261 + t277 * t298 - t286 * t328 + t366;
t177 = t350 * t183 + t345 * t184;
t300 = -mrSges(4,1) * t313 + mrSges(4,2) * t314;
t301 = -mrSges(4,2) * t330 + mrSges(4,3) * t313;
t175 = m(4) * t256 + mrSges(4,1) * t318 - mrSges(4,3) * t294 - t300 * t314 + t301 * t330 + t177;
t302 = mrSges(4,1) * t330 - mrSges(4,3) * t314;
t367 = -t183 * t345 + t350 * t184;
t176 = m(4) * t257 - mrSges(4,2) * t318 + mrSges(4,3) * t293 + t300 * t313 - t302 * t330 + t367;
t170 = t351 * t175 + t346 * t176;
t193 = t349 * t204 + t344 * t205;
t296 = -g(3) * t378 + t374;
t319 = mrSges(3,1) * t336 - mrSges(3,3) * t370;
t323 = (-mrSges(3,1) * t352 + mrSges(3,2) * t347) * t373;
t368 = -t175 * t346 + t351 * t176;
t168 = m(3) * t296 - mrSges(3,2) * t335 + mrSges(3,3) * t326 - t319 * t336 + t323 * t369 + t368;
t320 = -mrSges(3,2) * t336 + mrSges(3,3) * t369;
t363 = m(6) * t224 - t240 * mrSges(6,1) + t241 * mrSges(6,2) - t275 * t266 + t276 * t267 + t193;
t359 = m(5) * t248 - t261 * mrSges(5,1) + t262 * mrSges(5,2) - t298 * t284 + t299 * t286 + t363;
t356 = -m(4) * t281 + t293 * mrSges(4,1) - t294 * mrSges(4,2) + t313 * t301 - t314 * t302 - t359;
t188 = m(3) * t295 + t335 * mrSges(3,1) - t325 * mrSges(3,3) + t336 * t320 - t323 * t370 + t356;
t164 = t352 * t168 - t188 * t347;
t307 = -t341 * t321 - t379;
t169 = m(3) * t307 - t326 * mrSges(3,1) + t325 * mrSges(3,2) + (t319 * t347 - t320 * t352) * t373 + t170;
t161 = t168 * t376 - t169 * t341 + t188 * t375;
t228 = Ifges(7,5) * t265 + Ifges(7,6) * t264 + Ifges(7,3) * t274;
t230 = Ifges(7,1) * t265 + Ifges(7,4) * t264 + Ifges(7,5) * t274;
t197 = -mrSges(7,1) * t210 + mrSges(7,3) * t209 + Ifges(7,4) * t227 + Ifges(7,2) * t226 + Ifges(7,6) * t239 - t228 * t265 + t230 * t274;
t229 = Ifges(7,4) * t265 + Ifges(7,2) * t264 + Ifges(7,6) * t274;
t198 = mrSges(7,2) * t210 - mrSges(7,3) * t208 + Ifges(7,1) * t227 + Ifges(7,4) * t226 + Ifges(7,5) * t239 + t228 * t264 - t229 * t274;
t249 = Ifges(6,5) * t276 + Ifges(6,6) * t275 + Ifges(6,3) * t328;
t250 = Ifges(6,4) * t276 + Ifges(6,2) * t275 + Ifges(6,6) * t328;
t178 = mrSges(6,2) * t224 - mrSges(6,3) * t213 + Ifges(6,1) * t241 + Ifges(6,4) * t240 + Ifges(6,5) * t317 - pkin(11) * t193 - t197 * t344 + t198 * t349 + t249 * t275 - t250 * t328;
t251 = Ifges(6,1) * t276 + Ifges(6,4) * t275 + Ifges(6,5) * t328;
t358 = mrSges(7,1) * t208 - mrSges(7,2) * t209 + Ifges(7,5) * t227 + Ifges(7,6) * t226 + Ifges(7,3) * t239 + t229 * t265 - t230 * t264;
t179 = -mrSges(6,1) * t224 + mrSges(6,3) * t214 + Ifges(6,4) * t241 + Ifges(6,2) * t240 + Ifges(6,6) * t317 - pkin(5) * t193 - t249 * t276 + t251 * t328 - t358;
t268 = Ifges(5,5) * t299 + Ifges(5,6) * t298 + Ifges(5,3) * t328;
t270 = Ifges(5,1) * t299 + Ifges(5,4) * t298 + Ifges(5,5) * t328;
t165 = -mrSges(5,1) * t248 + mrSges(5,3) * t222 + Ifges(5,4) * t262 + Ifges(5,2) * t261 + Ifges(5,6) * t317 - pkin(4) * t363 + qJ(5) * t366 + t340 * t178 + t342 * t179 - t299 * t268 + t328 * t270;
t269 = Ifges(5,4) * t299 + Ifges(5,2) * t298 + Ifges(5,6) * t328;
t171 = mrSges(5,2) * t248 - mrSges(5,3) * t221 + Ifges(5,1) * t262 + Ifges(5,4) * t261 + Ifges(5,5) * t317 - qJ(5) * t186 + t178 * t342 - t179 * t340 + t268 * t298 - t269 * t328;
t287 = Ifges(4,5) * t314 + Ifges(4,6) * t313 + Ifges(4,3) * t330;
t289 = Ifges(4,1) * t314 + Ifges(4,4) * t313 + Ifges(4,5) * t330;
t154 = -mrSges(4,1) * t281 + mrSges(4,3) * t257 + Ifges(4,4) * t294 + Ifges(4,2) * t293 + Ifges(4,6) * t318 - pkin(3) * t359 + pkin(10) * t367 + t350 * t165 + t345 * t171 - t314 * t287 + t330 * t289;
t288 = Ifges(4,4) * t314 + Ifges(4,2) * t313 + Ifges(4,6) * t330;
t155 = mrSges(4,2) * t281 - mrSges(4,3) * t256 + Ifges(4,1) * t294 + Ifges(4,4) * t293 + Ifges(4,5) * t318 - pkin(10) * t177 - t165 * t345 + t171 * t350 + t287 * t313 - t288 * t330;
t305 = Ifges(3,6) * t336 + (Ifges(3,4) * t347 + Ifges(3,2) * t352) * t373;
t306 = Ifges(3,5) * t336 + (Ifges(3,1) * t347 + Ifges(3,4) * t352) * t373;
t151 = Ifges(3,5) * t325 + Ifges(3,6) * t326 + Ifges(3,3) * t335 + mrSges(3,1) * t295 - mrSges(3,2) * t296 + t346 * t155 + t351 * t154 + pkin(2) * t356 + pkin(9) * t368 + (t305 * t347 - t306 * t352) * t373;
t304 = Ifges(3,3) * t336 + (Ifges(3,5) * t347 + Ifges(3,6) * t352) * t373;
t153 = mrSges(3,2) * t307 - mrSges(3,3) * t295 + Ifges(3,1) * t325 + Ifges(3,4) * t326 + Ifges(3,5) * t335 - pkin(9) * t170 - t154 * t346 + t155 * t351 + t304 * t369 - t305 * t336;
t360 = -mrSges(6,1) * t213 + mrSges(6,2) * t214 - Ifges(6,5) * t241 - Ifges(6,6) * t240 - Ifges(6,3) * t317 - pkin(5) * t362 - pkin(11) * t365 - t349 * t197 - t344 * t198 - t276 * t250 + t275 * t251;
t357 = -mrSges(5,1) * t221 + mrSges(5,2) * t222 - Ifges(5,5) * t262 - Ifges(5,6) * t261 - Ifges(5,3) * t317 - pkin(4) * t186 - t299 * t269 + t298 * t270 + t360;
t355 = mrSges(4,1) * t256 - mrSges(4,2) * t257 + Ifges(4,5) * t294 + Ifges(4,6) * t293 + Ifges(4,3) * t318 + pkin(3) * t177 + t314 * t288 - t313 * t289 - t357;
t157 = -mrSges(3,1) * t307 + mrSges(3,3) * t296 + Ifges(3,4) * t325 + Ifges(3,2) * t326 + Ifges(3,6) * t335 - pkin(2) * t170 - t304 * t370 + t336 * t306 - t355;
t361 = mrSges(2,1) * t331 - mrSges(2,2) * t332 + Ifges(2,3) * qJDD(1) + pkin(1) * t161 + t343 * t151 + t153 * t378 + t157 * t377 + t164 * t380;
t162 = m(2) * t332 - mrSges(2,1) * t354 - qJDD(1) * mrSges(2,2) + t164;
t160 = t343 * t169 + (t168 * t347 + t188 * t352) * t341;
t158 = m(2) * t331 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t354 + t161;
t149 = -mrSges(2,2) * g(3) - mrSges(2,3) * t331 + Ifges(2,5) * qJDD(1) - t354 * Ifges(2,6) + t352 * t153 - t347 * t157 + (-t160 * t341 - t161 * t343) * pkin(8);
t148 = mrSges(2,1) * g(3) + mrSges(2,3) * t332 + t354 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t160 - t341 * t151 + (pkin(8) * t164 + t153 * t347 + t157 * t352) * t343;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t353 * t149 - t348 * t148 - pkin(7) * (t158 * t353 + t162 * t348), t149, t153, t155, t171, t178, t198; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t348 * t149 + t353 * t148 + pkin(7) * (-t158 * t348 + t162 * t353), t148, t157, t154, t165, t179, t197; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t361, t361, t151, t355, -t357, -t360, t358;];
m_new  = t1;
