% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-05-09 12:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRRR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-09 11:07:26
% EndTime: 2019-05-09 11:27:42
% DurationCPUTime: 533.51s
% Computational Cost: add. (8736372->429), mult. (23694513->587), div. (0->0), fcn. (20428384->18), ass. (0->196)
t368 = cos(pkin(6));
t358 = t368 * qJD(1) + qJD(2);
t363 = sin(pkin(7));
t367 = cos(pkin(7));
t364 = sin(pkin(6));
t377 = cos(qJ(2));
t396 = qJD(1) * t377;
t393 = t364 * t396;
t340 = (t358 * t363 + t367 * t393) * qJ(3);
t372 = sin(qJ(2));
t398 = qJD(1) * t364;
t414 = qJ(3) * t363;
t346 = (-pkin(2) * t377 - t372 * t414) * t398;
t395 = qJD(1) * qJD(2);
t352 = (qJDD(1) * t372 + t377 * t395) * t364;
t357 = t368 * qJDD(1) + qJDD(2);
t373 = sin(qJ(1));
t378 = cos(qJ(1));
t355 = t373 * g(1) - t378 * g(2);
t379 = qJD(1) ^ 2;
t418 = pkin(10) * t364;
t349 = qJDD(1) * pkin(1) + t379 * t418 + t355;
t356 = -t378 * g(1) - t373 * g(2);
t350 = -t379 * pkin(1) + qJDD(1) * t418 + t356;
t401 = t368 * t377;
t391 = t349 * t401 - t372 * t350;
t397 = qJD(1) * t372;
t413 = qJ(3) * t367;
t305 = -t352 * t413 + t357 * pkin(2) + t358 * t340 + (-g(3) * t377 - t346 * t397) * t364 + t391;
t394 = t364 * t397;
t345 = t358 * pkin(2) - t394 * t413;
t353 = (qJDD(1) * t377 - t372 * t395) * t364;
t387 = t353 * t367 + t357 * t363;
t402 = t368 * t372;
t399 = t349 * t402 + t377 * t350;
t306 = -t358 * t345 + (-g(3) * t372 + t346 * t396) * t364 + t387 * qJ(3) + t399;
t415 = t368 * g(3);
t310 = -t352 * t414 - t353 * pkin(2) - t415 + (-t349 + (-t340 * t377 + t345 * t372) * qJD(1)) * t364;
t361 = sin(pkin(14));
t365 = cos(pkin(14));
t403 = t367 * t377;
t411 = t361 * t363;
t334 = t358 * t411 + (t361 * t403 + t365 * t372) * t398;
t405 = t365 * t367;
t408 = t363 * t365;
t279 = -0.2e1 * qJD(3) * t334 + t305 * t405 - t306 * t361 + t310 * t408;
t333 = t358 * t408 + (-t361 * t372 + t365 * t403) * t398;
t362 = sin(pkin(8));
t417 = pkin(11) * t362;
t319 = -pkin(3) * t333 - t334 * t417;
t343 = t367 * t358 - t363 * t393;
t366 = cos(pkin(8));
t388 = t333 * t366 + t343 * t362;
t322 = t388 * pkin(11);
t328 = t365 * t352 + t387 * t361;
t335 = -t363 * t353 + t357 * t367;
t416 = pkin(11) * t366;
t263 = pkin(3) * t335 - t319 * t334 + t322 * t343 - t328 * t416 + t279;
t410 = t361 * t367;
t280 = 0.2e1 * qJD(3) * t333 + t305 * t410 + t365 * t306 + t310 * t411;
t324 = pkin(3) * t343 - t334 * t416;
t327 = -t361 * t352 + t387 * t365;
t389 = t327 * t366 + t335 * t362;
t264 = t389 * pkin(11) + t319 * t333 - t324 * t343 + t280;
t294 = -t305 * t363 + t367 * t310 + qJDD(3);
t268 = -pkin(3) * t327 - t322 * t333 + t324 * t334 - t328 * t417 + t294;
t371 = sin(qJ(4));
t376 = cos(qJ(4));
t250 = -t371 * t264 + (t263 * t366 + t268 * t362) * t376;
t314 = t376 * t334 + t388 * t371;
t292 = -t314 * qJD(4) - t371 * t328 + t389 * t376;
t313 = -t371 * t334 + t388 * t376;
t293 = t313 * qJD(4) + t376 * t328 + t389 * t371;
t295 = -mrSges(5,1) * t313 + mrSges(5,2) * t314;
t323 = -t333 * t362 + t343 * t366 + qJD(4);
t300 = -mrSges(5,2) * t323 + mrSges(5,3) * t313;
t318 = -t327 * t362 + t335 * t366 + qJDD(4);
t404 = t366 * t371;
t409 = t362 * t371;
t251 = t263 * t404 + t376 * t264 + t268 * t409;
t296 = -pkin(4) * t313 - pkin(12) * t314;
t321 = t323 ^ 2;
t247 = -pkin(4) * t321 + pkin(12) * t318 + t296 * t313 + t251;
t252 = -t263 * t362 + t366 * t268;
t249 = (-t313 * t323 - t293) * pkin(12) + (t314 * t323 - t292) * pkin(4) + t252;
t370 = sin(qJ(5));
t375 = cos(qJ(5));
t243 = t375 * t247 + t370 * t249;
t298 = -t370 * t314 + t375 * t323;
t299 = t375 * t314 + t370 * t323;
t282 = -pkin(5) * t298 - pkin(13) * t299;
t291 = qJDD(5) - t292;
t312 = qJD(5) - t313;
t311 = t312 ^ 2;
t241 = -pkin(5) * t311 + pkin(13) * t291 + t282 * t298 + t243;
t246 = -t318 * pkin(4) - t321 * pkin(12) + t314 * t296 - t250;
t271 = -t299 * qJD(5) - t370 * t293 + t375 * t318;
t272 = t298 * qJD(5) + t375 * t293 + t370 * t318;
t244 = (-t298 * t312 - t272) * pkin(13) + (t299 * t312 - t271) * pkin(5) + t246;
t369 = sin(qJ(6));
t374 = cos(qJ(6));
t237 = -t369 * t241 + t374 * t244;
t284 = -t369 * t299 + t374 * t312;
t255 = t284 * qJD(6) + t374 * t272 + t369 * t291;
t285 = t374 * t299 + t369 * t312;
t265 = -mrSges(7,1) * t284 + mrSges(7,2) * t285;
t270 = qJDD(6) - t271;
t297 = qJD(6) - t298;
t273 = -mrSges(7,2) * t297 + mrSges(7,3) * t284;
t235 = m(7) * t237 + mrSges(7,1) * t270 - mrSges(7,3) * t255 - t265 * t285 + t273 * t297;
t238 = t374 * t241 + t369 * t244;
t254 = -t285 * qJD(6) - t369 * t272 + t374 * t291;
t274 = mrSges(7,1) * t297 - mrSges(7,3) * t285;
t236 = m(7) * t238 - mrSges(7,2) * t270 + mrSges(7,3) * t254 + t265 * t284 - t274 * t297;
t228 = t374 * t235 + t369 * t236;
t286 = -mrSges(6,2) * t312 + mrSges(6,3) * t298;
t287 = mrSges(6,1) * t312 - mrSges(6,3) * t299;
t382 = -m(6) * t246 + t271 * mrSges(6,1) - t272 * mrSges(6,2) + t298 * t286 - t299 * t287 - t228;
t224 = m(5) * t250 + t318 * mrSges(5,1) - t293 * mrSges(5,3) - t314 * t295 + t323 * t300 + t382;
t412 = t224 * t376;
t407 = t364 * t372;
t406 = t364 * t377;
t229 = -t369 * t235 + t374 * t236;
t281 = -mrSges(6,1) * t298 + mrSges(6,2) * t299;
t227 = m(6) * t243 - t291 * mrSges(6,2) + t271 * mrSges(6,3) + t298 * t281 - t312 * t287 + t229;
t242 = -t370 * t247 + t375 * t249;
t240 = -t291 * pkin(5) - t311 * pkin(13) + t299 * t282 - t242;
t239 = -m(7) * t240 + t254 * mrSges(7,1) - mrSges(7,2) * t255 + t284 * t273 - t274 * t285;
t233 = m(6) * t242 + mrSges(6,1) * t291 - mrSges(6,3) * t272 - t281 * t299 + t286 * t312 + t239;
t221 = t370 * t227 + t375 * t233;
t301 = mrSges(5,1) * t323 - mrSges(5,3) * t314;
t392 = t375 * t227 - t370 * t233;
t218 = m(5) * t251 - t318 * mrSges(5,2) + t292 * mrSges(5,3) + t313 * t295 - t323 * t301 + t392;
t220 = m(5) * t252 - mrSges(5,1) * t292 + mrSges(5,2) * t293 - t300 * t313 + t301 * t314 + t221;
t207 = t218 * t404 - t220 * t362 + t366 * t412;
t320 = -mrSges(4,1) * t333 + mrSges(4,2) * t334;
t325 = -mrSges(4,2) * t343 + mrSges(4,3) * t333;
t203 = m(4) * t279 + mrSges(4,1) * t335 - mrSges(4,3) * t328 - t320 * t334 + t325 * t343 + t207;
t206 = t218 * t409 + t366 * t220 + t362 * t412;
t326 = mrSges(4,1) * t343 - mrSges(4,3) * t334;
t205 = m(4) * t294 - mrSges(4,1) * t327 + mrSges(4,2) * t328 - t325 * t333 + t326 * t334 + t206;
t212 = t376 * t218 - t371 * t224;
t211 = m(4) * t280 - t335 * mrSges(4,2) + t327 * mrSges(4,3) + t333 * t320 - t343 * t326 + t212;
t192 = t203 * t408 + t367 * t205 + t211 * t411;
t193 = t203 * t405 - t363 * t205 + t211 * t410;
t329 = -g(3) * t406 + t391;
t348 = -t358 * mrSges(3,2) + mrSges(3,3) * t393;
t351 = (-mrSges(3,1) * t377 + mrSges(3,2) * t372) * t398;
t190 = m(3) * t329 + t357 * mrSges(3,1) - t352 * mrSges(3,3) + t358 * t348 - t351 * t394 + t193;
t197 = -t361 * t203 + t365 * t211;
t330 = -g(3) * t407 + t399;
t347 = t358 * mrSges(3,1) - mrSges(3,3) * t394;
t196 = m(3) * t330 - t357 * mrSges(3,2) + t353 * mrSges(3,3) - t358 * t347 + t351 * t393 + t197;
t188 = -t372 * t190 + t377 * t196;
t339 = -t364 * t349 - t415;
t191 = m(3) * t339 - t353 * mrSges(3,1) + t352 * mrSges(3,2) + (t347 * t372 - t348 * t377) * t398 + t192;
t182 = t190 * t401 - t191 * t364 + t196 * t402;
t256 = Ifges(7,5) * t285 + Ifges(7,6) * t284 + Ifges(7,3) * t297;
t258 = Ifges(7,1) * t285 + Ifges(7,4) * t284 + Ifges(7,5) * t297;
t230 = -mrSges(7,1) * t240 + mrSges(7,3) * t238 + Ifges(7,4) * t255 + Ifges(7,2) * t254 + Ifges(7,6) * t270 - t256 * t285 + t258 * t297;
t257 = Ifges(7,4) * t285 + Ifges(7,2) * t284 + Ifges(7,6) * t297;
t231 = mrSges(7,2) * t240 - mrSges(7,3) * t237 + Ifges(7,1) * t255 + Ifges(7,4) * t254 + Ifges(7,5) * t270 + t256 * t284 - t257 * t297;
t275 = Ifges(6,5) * t299 + Ifges(6,6) * t298 + Ifges(6,3) * t312;
t276 = Ifges(6,4) * t299 + Ifges(6,2) * t298 + Ifges(6,6) * t312;
t213 = mrSges(6,2) * t246 - mrSges(6,3) * t242 + Ifges(6,1) * t272 + Ifges(6,4) * t271 + Ifges(6,5) * t291 - pkin(13) * t228 - t369 * t230 + t374 * t231 + t298 * t275 - t312 * t276;
t277 = Ifges(6,1) * t299 + Ifges(6,4) * t298 + Ifges(6,5) * t312;
t381 = mrSges(7,1) * t237 - mrSges(7,2) * t238 + Ifges(7,5) * t255 + Ifges(7,6) * t254 + Ifges(7,3) * t270 + t257 * t285 - t258 * t284;
t214 = -mrSges(6,1) * t246 + mrSges(6,3) * t243 + Ifges(6,4) * t272 + Ifges(6,2) * t271 + Ifges(6,6) * t291 - pkin(5) * t228 - t275 * t299 + t277 * t312 - t381;
t288 = Ifges(5,5) * t314 + Ifges(5,6) * t313 + Ifges(5,3) * t323;
t289 = Ifges(5,4) * t314 + Ifges(5,2) * t313 + Ifges(5,6) * t323;
t199 = mrSges(5,2) * t252 - mrSges(5,3) * t250 + Ifges(5,1) * t293 + Ifges(5,4) * t292 + Ifges(5,5) * t318 - pkin(12) * t221 + t375 * t213 - t370 * t214 + t313 * t288 - t323 * t289;
t290 = Ifges(5,1) * t314 + Ifges(5,4) * t313 + Ifges(5,5) * t323;
t380 = mrSges(6,1) * t242 - mrSges(6,2) * t243 + Ifges(6,5) * t272 + Ifges(6,6) * t271 + Ifges(6,3) * t291 + pkin(5) * t239 + pkin(13) * t229 + t374 * t230 + t369 * t231 + t299 * t276 - t298 * t277;
t200 = -mrSges(5,1) * t252 + mrSges(5,3) * t251 + Ifges(5,4) * t293 + Ifges(5,2) * t292 + Ifges(5,6) * t318 - pkin(4) * t221 - t314 * t288 + t323 * t290 - t380;
t385 = pkin(11) * t212 + t199 * t371 + t200 * t376;
t198 = mrSges(5,1) * t250 - mrSges(5,2) * t251 + Ifges(5,5) * t293 + Ifges(5,6) * t292 + Ifges(5,3) * t318 + pkin(4) * t382 + pkin(12) * t392 + t370 * t213 + t375 * t214 + t314 * t289 - t313 * t290;
t315 = Ifges(4,5) * t334 + Ifges(4,6) * t333 + Ifges(4,3) * t343;
t317 = Ifges(4,1) * t334 + Ifges(4,4) * t333 + Ifges(4,5) * t343;
t184 = -mrSges(4,1) * t294 + mrSges(4,3) * t280 + Ifges(4,4) * t328 + Ifges(4,2) * t327 + Ifges(4,6) * t335 - pkin(3) * t206 - t362 * t198 - t334 * t315 + t343 * t317 + t385 * t366;
t316 = Ifges(4,4) * t334 + Ifges(4,2) * t333 + Ifges(4,6) * t343;
t185 = mrSges(4,2) * t294 - mrSges(4,3) * t279 + Ifges(4,1) * t328 + Ifges(4,4) * t327 + Ifges(4,5) * t335 + t376 * t199 - t371 * t200 + t333 * t315 - t343 * t316 + (-t206 * t362 - t207 * t366) * pkin(11);
t384 = qJ(3) * t197 + t184 * t365 + t185 * t361;
t183 = mrSges(4,1) * t279 - mrSges(4,2) * t280 + Ifges(4,5) * t328 + Ifges(4,6) * t327 + Ifges(4,3) * t335 + pkin(3) * t207 + t366 * t198 + t334 * t316 - t333 * t317 + t385 * t362;
t337 = Ifges(3,6) * t358 + (Ifges(3,4) * t372 + Ifges(3,2) * t377) * t398;
t338 = Ifges(3,5) * t358 + (Ifges(3,1) * t372 + Ifges(3,4) * t377) * t398;
t174 = mrSges(3,1) * t329 - mrSges(3,2) * t330 + Ifges(3,5) * t352 + Ifges(3,6) * t353 + Ifges(3,3) * t357 + pkin(2) * t193 + t367 * t183 + (t337 * t372 - t338 * t377) * t398 + t384 * t363;
t336 = Ifges(3,3) * t358 + (Ifges(3,5) * t372 + Ifges(3,6) * t377) * t398;
t176 = -mrSges(3,1) * t339 + mrSges(3,3) * t330 + Ifges(3,4) * t352 + Ifges(3,2) * t353 + Ifges(3,6) * t357 - pkin(2) * t192 - t363 * t183 - t336 * t394 + t358 * t338 + t384 * t367;
t178 = t336 * t393 + mrSges(3,2) * t339 - mrSges(3,3) * t329 + Ifges(3,1) * t352 + Ifges(3,4) * t353 + Ifges(3,5) * t357 - t361 * t184 + t365 * t185 - t358 * t337 + (-t192 * t363 - t193 * t367) * qJ(3);
t383 = mrSges(2,1) * t355 - mrSges(2,2) * t356 + Ifges(2,3) * qJDD(1) + pkin(1) * t182 + t368 * t174 + t176 * t406 + t178 * t407 + t188 * t418;
t186 = m(2) * t356 - t379 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t188;
t181 = t368 * t191 + (t190 * t377 + t196 * t372) * t364;
t179 = m(2) * t355 + qJDD(1) * mrSges(2,1) - t379 * mrSges(2,2) + t182;
t172 = -mrSges(2,2) * g(3) - mrSges(2,3) * t355 + Ifges(2,5) * qJDD(1) - t379 * Ifges(2,6) - t372 * t176 + t377 * t178 + (-t181 * t364 - t182 * t368) * pkin(10);
t171 = mrSges(2,1) * g(3) + mrSges(2,3) * t356 + t379 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t181 - t364 * t174 + (pkin(10) * t188 + t176 * t377 + t178 * t372) * t368;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t378 * t172 - t373 * t171 - pkin(9) * (t378 * t179 + t373 * t186), t172, t178, t185, t199, t213, t231; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t373 * t172 + t378 * t171 + pkin(9) * (-t373 * t179 + t378 * t186), t171, t176, t184, t200, t214, t230; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t383, t383, t174, t183, t198, t380, t381;];
m_new  = t1;
