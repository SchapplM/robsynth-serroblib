% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 19:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 18:06:19
% EndTime: 2019-05-08 18:20:55
% DurationCPUTime: 570.88s
% Computational Cost: add. (9953773->430), mult. (25520788->585), div. (0->0), fcn. (22118114->18), ass. (0->196)
t361 = cos(pkin(6));
t353 = t361 * qJD(1) + qJD(2);
t357 = sin(pkin(7));
t360 = cos(pkin(7));
t358 = sin(pkin(6));
t372 = cos(qJ(2));
t391 = qJD(1) * t372;
t388 = t358 * t391;
t337 = (t353 * t357 + t360 * t388) * pkin(11);
t366 = sin(qJ(2));
t393 = qJD(1) * t358;
t412 = pkin(11) * t357;
t341 = (-pkin(2) * t372 - t366 * t412) * t393;
t390 = qJD(1) * qJD(2);
t347 = (qJDD(1) * t366 + t372 * t390) * t358;
t352 = t361 * qJDD(1) + qJDD(2);
t367 = sin(qJ(1));
t373 = cos(qJ(1));
t350 = t367 * g(1) - t373 * g(2);
t374 = qJD(1) ^ 2;
t413 = pkin(10) * t358;
t344 = qJDD(1) * pkin(1) + t374 * t413 + t350;
t351 = -t373 * g(1) - t367 * g(2);
t345 = -t374 * pkin(1) + qJDD(1) * t413 + t351;
t396 = t361 * t372;
t386 = t344 * t396 - t366 * t345;
t392 = qJD(1) * t366;
t411 = pkin(11) * t360;
t302 = -t347 * t411 + t352 * pkin(2) + t353 * t337 + (-g(3) * t372 - t341 * t392) * t358 + t386;
t389 = t358 * t392;
t340 = t353 * pkin(2) - t389 * t411;
t348 = (qJDD(1) * t372 - t366 * t390) * t358;
t382 = t348 * t360 + t352 * t357;
t397 = t361 * t366;
t394 = t344 * t397 + t372 * t345;
t303 = -t353 * t340 + (-g(3) * t366 + t341 * t391) * t358 + t382 * pkin(11) + t394;
t408 = t361 * g(3);
t308 = -t347 * t412 - t348 * pkin(2) - t408 + (-t344 + (-t337 * t372 + t340 * t366) * qJD(1)) * t358;
t365 = sin(qJ(3));
t371 = cos(qJ(3));
t399 = t360 * t371;
t404 = t357 * t371;
t276 = t302 * t399 - t365 * t303 + t308 * t404;
t398 = t360 * t372;
t328 = t353 * t404 + (-t365 * t366 + t371 * t398) * t393;
t317 = t328 * qJD(3) + t371 * t347 + t382 * t365;
t405 = t357 * t365;
t329 = t353 * t405 + (t365 * t398 + t366 * t371) * t393;
t356 = sin(pkin(8));
t410 = pkin(12) * t356;
t318 = -t328 * pkin(3) - t329 * t410;
t338 = t360 * t353 - t357 * t388 + qJD(3);
t359 = cos(pkin(8));
t383 = t328 * t359 + t338 * t356;
t321 = t383 * pkin(12);
t330 = -t357 * t348 + t360 * t352 + qJDD(3);
t409 = pkin(12) * t359;
t260 = t330 * pkin(3) - t317 * t409 - t329 * t318 + t338 * t321 + t276;
t400 = t360 * t365;
t277 = t302 * t400 + t371 * t303 + t308 * t405;
t323 = t338 * pkin(3) - t329 * t409;
t316 = -t329 * qJD(3) - t365 * t347 + t382 * t371;
t384 = t316 * t359 + t330 * t356;
t261 = t384 * pkin(12) + t328 * t318 - t338 * t323 + t277;
t291 = -t357 * t302 + t360 * t308;
t268 = -t316 * pkin(3) - t317 * t410 - t328 * t321 + t329 * t323 + t291;
t364 = sin(qJ(4));
t370 = cos(qJ(4));
t247 = -t364 * t261 + (t260 * t359 + t268 * t356) * t370;
t312 = t370 * t329 + t383 * t364;
t281 = -t312 * qJD(4) - t364 * t317 + t384 * t370;
t311 = -t364 * t329 + t383 * t370;
t282 = t311 * qJD(4) + t370 * t317 + t384 * t364;
t292 = -t311 * mrSges(5,1) + t312 * mrSges(5,2);
t322 = -t356 * t328 + t359 * t338 + qJD(4);
t297 = -t322 * mrSges(5,2) + t311 * mrSges(5,3);
t304 = -t356 * t316 + t359 * t330 + qJDD(4);
t401 = t359 * t364;
t406 = t356 * t364;
t248 = t260 * t401 + t370 * t261 + t268 * t406;
t293 = -t311 * pkin(4) - t312 * pkin(13);
t320 = t322 ^ 2;
t244 = -t320 * pkin(4) + t304 * pkin(13) + t311 * t293 + t248;
t249 = -t356 * t260 + t359 * t268;
t246 = (-t311 * t322 - t282) * pkin(13) + (t312 * t322 - t281) * pkin(4) + t249;
t363 = sin(qJ(5));
t369 = cos(qJ(5));
t240 = t369 * t244 + t363 * t246;
t295 = -t363 * t312 + t369 * t322;
t296 = t369 * t312 + t363 * t322;
t279 = -t295 * pkin(5) - t296 * pkin(14);
t280 = qJDD(5) - t281;
t310 = qJD(5) - t311;
t309 = t310 ^ 2;
t238 = -t309 * pkin(5) + t280 * pkin(14) + t295 * t279 + t240;
t243 = -t304 * pkin(4) - t320 * pkin(13) + t312 * t293 - t247;
t264 = -t296 * qJD(5) - t363 * t282 + t369 * t304;
t265 = t295 * qJD(5) + t369 * t282 + t363 * t304;
t241 = (-t295 * t310 - t265) * pkin(14) + (t296 * t310 - t264) * pkin(5) + t243;
t362 = sin(qJ(6));
t368 = cos(qJ(6));
t234 = -t362 * t238 + t368 * t241;
t284 = -t362 * t296 + t368 * t310;
t252 = t284 * qJD(6) + t368 * t265 + t362 * t280;
t263 = qJDD(6) - t264;
t285 = t368 * t296 + t362 * t310;
t269 = -t284 * mrSges(7,1) + t285 * mrSges(7,2);
t294 = qJD(6) - t295;
t270 = -t294 * mrSges(7,2) + t284 * mrSges(7,3);
t232 = m(7) * t234 + t263 * mrSges(7,1) - t252 * mrSges(7,3) - t285 * t269 + t294 * t270;
t235 = t368 * t238 + t362 * t241;
t251 = -t285 * qJD(6) - t362 * t265 + t368 * t280;
t271 = t294 * mrSges(7,1) - t285 * mrSges(7,3);
t233 = m(7) * t235 - t263 * mrSges(7,2) + t251 * mrSges(7,3) + t284 * t269 - t294 * t271;
t225 = t368 * t232 + t362 * t233;
t286 = -t310 * mrSges(6,2) + t295 * mrSges(6,3);
t287 = t310 * mrSges(6,1) - t296 * mrSges(6,3);
t377 = -m(6) * t243 + t264 * mrSges(6,1) - t265 * mrSges(6,2) + t295 * t286 - t296 * t287 - t225;
t221 = m(5) * t247 + t304 * mrSges(5,1) - t282 * mrSges(5,3) - t312 * t292 + t322 * t297 + t377;
t407 = t221 * t370;
t403 = t358 * t366;
t402 = t358 * t372;
t226 = -t362 * t232 + t368 * t233;
t278 = -t295 * mrSges(6,1) + t296 * mrSges(6,2);
t224 = m(6) * t240 - t280 * mrSges(6,2) + t264 * mrSges(6,3) + t295 * t278 - t310 * t287 + t226;
t239 = -t363 * t244 + t369 * t246;
t237 = -t280 * pkin(5) - t309 * pkin(14) + t296 * t279 - t239;
t236 = -m(7) * t237 + t251 * mrSges(7,1) - t252 * mrSges(7,2) + t284 * t270 - t285 * t271;
t230 = m(6) * t239 + t280 * mrSges(6,1) - t265 * mrSges(6,3) - t296 * t278 + t310 * t286 + t236;
t218 = t363 * t224 + t369 * t230;
t298 = t322 * mrSges(5,1) - t312 * mrSges(5,3);
t387 = t369 * t224 - t363 * t230;
t215 = m(5) * t248 - t304 * mrSges(5,2) + t281 * mrSges(5,3) + t311 * t292 - t322 * t298 + t387;
t217 = m(5) * t249 - t281 * mrSges(5,1) + t282 * mrSges(5,2) - t311 * t297 + t312 * t298 + t218;
t204 = t215 * t401 - t356 * t217 + t359 * t407;
t319 = -t328 * mrSges(4,1) + t329 * mrSges(4,2);
t324 = -t338 * mrSges(4,2) + t328 * mrSges(4,3);
t200 = m(4) * t276 + t330 * mrSges(4,1) - t317 * mrSges(4,3) - t329 * t319 + t338 * t324 + t204;
t203 = t215 * t406 + t359 * t217 + t356 * t407;
t325 = t338 * mrSges(4,1) - t329 * mrSges(4,3);
t202 = m(4) * t291 - t316 * mrSges(4,1) + t317 * mrSges(4,2) - t328 * t324 + t329 * t325 + t203;
t209 = t370 * t215 - t364 * t221;
t208 = m(4) * t277 - t330 * mrSges(4,2) + t316 * mrSges(4,3) + t328 * t319 - t338 * t325 + t209;
t189 = t200 * t404 + t360 * t202 + t208 * t405;
t190 = t200 * t399 - t357 * t202 + t208 * t400;
t326 = -g(3) * t402 + t386;
t343 = -t353 * mrSges(3,2) + mrSges(3,3) * t388;
t346 = (-mrSges(3,1) * t372 + mrSges(3,2) * t366) * t393;
t187 = m(3) * t326 + t352 * mrSges(3,1) - t347 * mrSges(3,3) + t353 * t343 - t346 * t389 + t190;
t194 = -t365 * t200 + t371 * t208;
t327 = -g(3) * t403 + t394;
t342 = t353 * mrSges(3,1) - mrSges(3,3) * t389;
t193 = m(3) * t327 - t352 * mrSges(3,2) + t348 * mrSges(3,3) - t353 * t342 + t346 * t388 + t194;
t185 = -t366 * t187 + t372 * t193;
t334 = -t358 * t344 - t408;
t188 = m(3) * t334 - t348 * mrSges(3,1) + t347 * mrSges(3,2) + (t342 * t366 - t343 * t372) * t393 + t189;
t179 = t187 * t396 - t358 * t188 + t193 * t397;
t253 = Ifges(7,5) * t285 + Ifges(7,6) * t284 + Ifges(7,3) * t294;
t255 = Ifges(7,1) * t285 + Ifges(7,4) * t284 + Ifges(7,5) * t294;
t227 = -mrSges(7,1) * t237 + mrSges(7,3) * t235 + Ifges(7,4) * t252 + Ifges(7,2) * t251 + Ifges(7,6) * t263 - t285 * t253 + t294 * t255;
t254 = Ifges(7,4) * t285 + Ifges(7,2) * t284 + Ifges(7,6) * t294;
t228 = mrSges(7,2) * t237 - mrSges(7,3) * t234 + Ifges(7,1) * t252 + Ifges(7,4) * t251 + Ifges(7,5) * t263 + t284 * t253 - t294 * t254;
t272 = Ifges(6,5) * t296 + Ifges(6,6) * t295 + Ifges(6,3) * t310;
t273 = Ifges(6,4) * t296 + Ifges(6,2) * t295 + Ifges(6,6) * t310;
t210 = mrSges(6,2) * t243 - mrSges(6,3) * t239 + Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t280 - pkin(14) * t225 - t362 * t227 + t368 * t228 + t295 * t272 - t310 * t273;
t274 = Ifges(6,1) * t296 + Ifges(6,4) * t295 + Ifges(6,5) * t310;
t376 = mrSges(7,1) * t234 - mrSges(7,2) * t235 + Ifges(7,5) * t252 + Ifges(7,6) * t251 + Ifges(7,3) * t263 + t285 * t254 - t284 * t255;
t211 = -mrSges(6,1) * t243 + mrSges(6,3) * t240 + Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t280 - pkin(5) * t225 - t296 * t272 + t310 * t274 - t376;
t289 = Ifges(5,4) * t312 + Ifges(5,2) * t311 + Ifges(5,6) * t322;
t290 = Ifges(5,1) * t312 + Ifges(5,4) * t311 + Ifges(5,5) * t322;
t195 = mrSges(5,1) * t247 - mrSges(5,2) * t248 + Ifges(5,5) * t282 + Ifges(5,6) * t281 + Ifges(5,3) * t304 + pkin(4) * t377 + pkin(13) * t387 + t363 * t210 + t369 * t211 + t312 * t289 - t311 * t290;
t313 = Ifges(4,5) * t329 + Ifges(4,6) * t328 + Ifges(4,3) * t338;
t315 = Ifges(4,1) * t329 + Ifges(4,4) * t328 + Ifges(4,5) * t338;
t288 = Ifges(5,5) * t312 + Ifges(5,6) * t311 + Ifges(5,3) * t322;
t196 = mrSges(5,2) * t249 - mrSges(5,3) * t247 + Ifges(5,1) * t282 + Ifges(5,4) * t281 + Ifges(5,5) * t304 - pkin(13) * t218 + t369 * t210 - t363 * t211 + t311 * t288 - t322 * t289;
t375 = mrSges(6,1) * t239 - mrSges(6,2) * t240 + Ifges(6,5) * t265 + Ifges(6,6) * t264 + Ifges(6,3) * t280 + pkin(5) * t236 + pkin(14) * t226 + t368 * t227 + t362 * t228 + t296 * t273 - t295 * t274;
t197 = -mrSges(5,1) * t249 + mrSges(5,3) * t248 + Ifges(5,4) * t282 + Ifges(5,2) * t281 + Ifges(5,6) * t304 - pkin(4) * t218 - t312 * t288 + t322 * t290 - t375;
t379 = pkin(12) * t209 + t196 * t364 + t197 * t370;
t181 = -mrSges(4,1) * t291 + mrSges(4,3) * t277 + Ifges(4,4) * t317 + Ifges(4,2) * t316 + Ifges(4,6) * t330 - pkin(3) * t203 - t356 * t195 - t329 * t313 + t338 * t315 + t379 * t359;
t314 = Ifges(4,4) * t329 + Ifges(4,2) * t328 + Ifges(4,6) * t338;
t182 = mrSges(4,2) * t291 - mrSges(4,3) * t276 + Ifges(4,1) * t317 + Ifges(4,4) * t316 + Ifges(4,5) * t330 + t370 * t196 - t364 * t197 + t328 * t313 - t338 * t314 + (-t203 * t356 - t204 * t359) * pkin(12);
t380 = pkin(11) * t194 + t181 * t371 + t182 * t365;
t180 = mrSges(4,1) * t276 - mrSges(4,2) * t277 + Ifges(4,5) * t317 + Ifges(4,6) * t316 + Ifges(4,3) * t330 + pkin(3) * t204 + t359 * t195 + t329 * t314 - t328 * t315 + t379 * t356;
t332 = Ifges(3,6) * t353 + (Ifges(3,4) * t366 + Ifges(3,2) * t372) * t393;
t333 = Ifges(3,5) * t353 + (Ifges(3,1) * t366 + Ifges(3,4) * t372) * t393;
t171 = mrSges(3,1) * t326 - mrSges(3,2) * t327 + Ifges(3,5) * t347 + Ifges(3,6) * t348 + Ifges(3,3) * t352 + pkin(2) * t190 + t360 * t180 + (t332 * t366 - t333 * t372) * t393 + t380 * t357;
t331 = Ifges(3,3) * t353 + (Ifges(3,5) * t366 + Ifges(3,6) * t372) * t393;
t173 = -mrSges(3,1) * t334 + mrSges(3,3) * t327 + Ifges(3,4) * t347 + Ifges(3,2) * t348 + Ifges(3,6) * t352 - pkin(2) * t189 - t357 * t180 - t331 * t389 + t353 * t333 + t380 * t360;
t175 = t331 * t388 + mrSges(3,2) * t334 - mrSges(3,3) * t326 + Ifges(3,1) * t347 + Ifges(3,4) * t348 + Ifges(3,5) * t352 - t365 * t181 + t371 * t182 - t353 * t332 + (-t189 * t357 - t190 * t360) * pkin(11);
t378 = mrSges(2,1) * t350 - mrSges(2,2) * t351 + Ifges(2,3) * qJDD(1) + pkin(1) * t179 + t361 * t171 + t173 * t402 + t175 * t403 + t185 * t413;
t183 = m(2) * t351 - t374 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t185;
t178 = t361 * t188 + (t187 * t372 + t193 * t366) * t358;
t176 = m(2) * t350 + qJDD(1) * mrSges(2,1) - t374 * mrSges(2,2) + t179;
t169 = -mrSges(2,2) * g(3) - mrSges(2,3) * t350 + Ifges(2,5) * qJDD(1) - t374 * Ifges(2,6) - t366 * t173 + t372 * t175 + (-t178 * t358 - t179 * t361) * pkin(10);
t168 = mrSges(2,1) * g(3) + mrSges(2,3) * t351 + t374 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t178 - t358 * t171 + (pkin(10) * t185 + t173 * t372 + t175 * t366) * t361;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t373 * t169 - t367 * t168 - pkin(9) * (t373 * t176 + t367 * t183), t169, t175, t182, t196, t210, t228; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t367 * t169 + t373 * t168 + pkin(9) * (-t367 * t176 + t373 * t183), t168, t173, t181, t197, t211, t227; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t378, t378, t171, t180, t195, t375, t376;];
m_new  = t1;
