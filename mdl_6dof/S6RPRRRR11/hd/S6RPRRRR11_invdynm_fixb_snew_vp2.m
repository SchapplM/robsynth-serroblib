% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 05:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 05:31:05
% EndTime: 2019-05-06 05:34:22
% DurationCPUTime: 185.09s
% Computational Cost: add. (2984283->402), mult. (9253685->550), div. (0->0), fcn. (7981344->16), ass. (0->185)
t349 = sin(pkin(13));
t351 = sin(pkin(6));
t352 = cos(pkin(13));
t354 = cos(pkin(6));
t358 = sin(qJ(3));
t353 = cos(pkin(7));
t363 = cos(qJ(3));
t396 = t353 * t363;
t350 = sin(pkin(7));
t401 = t350 * t363;
t369 = t351 * (-t349 * t358 + t352 * t396) + t354 * t401;
t322 = t369 * qJD(1);
t397 = t353 * t358;
t402 = t350 * t358;
t371 = t354 * t402 + (t349 * t363 + t352 * t397) * t351;
t323 = t371 * qJD(1);
t308 = -t323 * qJD(3) + qJDD(1) * t369;
t399 = t351 * t353;
t334 = (t350 * t354 + t352 * t399) * qJD(1) * pkin(9);
t359 = sin(qJ(1));
t364 = cos(qJ(1));
t346 = -g(1) * t364 - g(2) * t359;
t365 = qJD(1) ^ 2;
t405 = qJ(2) * t351;
t338 = -pkin(1) * t365 + qJDD(1) * t405 + t346;
t408 = pkin(9) * t349;
t383 = -pkin(2) * t352 - t350 * t408;
t395 = qJD(1) * t351;
t406 = pkin(9) * qJDD(1);
t378 = qJD(1) * t383 * t395 + t353 * t406;
t345 = t359 * g(1) - g(2) * t364;
t337 = qJDD(1) * pkin(1) + t365 * t405 + t345;
t391 = qJD(2) * t395;
t398 = t352 * t354;
t400 = t351 * t352;
t384 = -g(3) * t400 + t337 * t398 - 0.2e1 * t349 * t391;
t288 = (pkin(2) * qJDD(1) + qJD(1) * t334) * t354 + (-t351 * t378 - t338) * t349 + t384;
t339 = (pkin(2) * t354 - t399 * t408) * qJD(1);
t403 = t349 * t354;
t392 = t337 * t403 + (t338 + 0.2e1 * t391) * t352;
t289 = (-qJD(1) * t339 + t350 * t406) * t354 + (-g(3) * t349 + t352 * t378) * t351 + t392;
t390 = -t354 * g(3) + qJDD(2);
t298 = (-t337 + t383 * qJDD(1) + (-t334 * t352 + t339 * t349) * qJD(1)) * t351 + t390;
t256 = -t358 * t289 + (t288 * t353 + t298 * t350) * t363;
t407 = Ifges(3,3) * t354;
t404 = t349 * t351;
t257 = t288 * t397 + t363 * t289 + t298 * t402;
t307 = -pkin(3) * t322 - pkin(10) * t323;
t379 = -t350 * t400 + t353 * t354;
t335 = qJD(1) * t379 + qJD(3);
t331 = t335 ^ 2;
t332 = qJDD(1) * t379 + qJDD(3);
t248 = -pkin(3) * t331 + pkin(10) * t332 + t307 * t322 + t257;
t267 = -t350 * t288 + t353 * t298;
t309 = t322 * qJD(3) + qJDD(1) * t371;
t253 = (-t322 * t335 - t309) * pkin(10) + (t323 * t335 - t308) * pkin(3) + t267;
t357 = sin(qJ(4));
t362 = cos(qJ(4));
t238 = t362 * t248 + t357 * t253;
t315 = -t357 * t323 + t335 * t362;
t316 = t323 * t362 + t335 * t357;
t291 = -pkin(4) * t315 - pkin(11) * t316;
t305 = qJDD(4) - t308;
t321 = qJD(4) - t322;
t320 = t321 ^ 2;
t233 = -pkin(4) * t320 + pkin(11) * t305 + t291 * t315 + t238;
t247 = -t332 * pkin(3) - t331 * pkin(10) + t323 * t307 - t256;
t283 = -t316 * qJD(4) - t357 * t309 + t332 * t362;
t284 = qJD(4) * t315 + t309 * t362 + t332 * t357;
t236 = (-t315 * t321 - t284) * pkin(11) + (t316 * t321 - t283) * pkin(4) + t247;
t356 = sin(qJ(5));
t361 = cos(qJ(5));
t228 = -t356 * t233 + t361 * t236;
t296 = -t316 * t356 + t321 * t361;
t260 = qJD(5) * t296 + t284 * t361 + t305 * t356;
t281 = qJDD(5) - t283;
t297 = t316 * t361 + t321 * t356;
t312 = qJD(5) - t315;
t226 = (t296 * t312 - t260) * pkin(12) + (t296 * t297 + t281) * pkin(5) + t228;
t229 = t361 * t233 + t356 * t236;
t259 = -qJD(5) * t297 - t284 * t356 + t305 * t361;
t274 = pkin(5) * t312 - pkin(12) * t297;
t295 = t296 ^ 2;
t227 = -pkin(5) * t295 + pkin(12) * t259 - t274 * t312 + t229;
t355 = sin(qJ(6));
t360 = cos(qJ(6));
t224 = t226 * t360 - t227 * t355;
t268 = t296 * t360 - t297 * t355;
t243 = qJD(6) * t268 + t259 * t355 + t260 * t360;
t269 = t296 * t355 + t297 * t360;
t255 = -mrSges(7,1) * t268 + mrSges(7,2) * t269;
t310 = qJD(6) + t312;
t261 = -mrSges(7,2) * t310 + mrSges(7,3) * t268;
t276 = qJDD(6) + t281;
t220 = m(7) * t224 + mrSges(7,1) * t276 - mrSges(7,3) * t243 - t255 * t269 + t261 * t310;
t225 = t226 * t355 + t227 * t360;
t242 = -qJD(6) * t269 + t259 * t360 - t260 * t355;
t262 = mrSges(7,1) * t310 - mrSges(7,3) * t269;
t221 = m(7) * t225 - mrSges(7,2) * t276 + mrSges(7,3) * t242 + t255 * t268 - t262 * t310;
t212 = t360 * t220 + t355 * t221;
t270 = -mrSges(6,1) * t296 + mrSges(6,2) * t297;
t272 = -mrSges(6,2) * t312 + mrSges(6,3) * t296;
t210 = m(6) * t228 + mrSges(6,1) * t281 - mrSges(6,3) * t260 - t270 * t297 + t272 * t312 + t212;
t273 = mrSges(6,1) * t312 - mrSges(6,3) * t297;
t388 = -t220 * t355 + t360 * t221;
t211 = m(6) * t229 - mrSges(6,2) * t281 + mrSges(6,3) * t259 + t270 * t296 - t273 * t312 + t388;
t208 = -t210 * t356 + t361 * t211;
t290 = -mrSges(5,1) * t315 + mrSges(5,2) * t316;
t300 = mrSges(5,1) * t321 - mrSges(5,3) * t316;
t206 = m(5) * t238 - mrSges(5,2) * t305 + mrSges(5,3) * t283 + t290 * t315 - t300 * t321 + t208;
t237 = -t357 * t248 + t253 * t362;
t232 = -pkin(4) * t305 - pkin(11) * t320 + t316 * t291 - t237;
t230 = -pkin(5) * t259 - pkin(12) * t295 + t274 * t297 + t232;
t376 = m(7) * t230 - t242 * mrSges(7,1) + mrSges(7,2) * t243 - t268 * t261 + t262 * t269;
t222 = -m(6) * t232 + t259 * mrSges(6,1) - mrSges(6,2) * t260 + t296 * t272 - t273 * t297 - t376;
t299 = -mrSges(5,2) * t321 + mrSges(5,3) * t315;
t216 = m(5) * t237 + mrSges(5,1) * t305 - mrSges(5,3) * t284 - t290 * t316 + t299 * t321 + t222;
t198 = t357 * t206 + t362 * t216;
t306 = -mrSges(4,1) * t322 + mrSges(4,2) * t323;
t318 = mrSges(4,1) * t335 - mrSges(4,3) * t323;
t389 = t362 * t206 - t216 * t357;
t195 = m(4) * t257 - mrSges(4,2) * t332 + mrSges(4,3) * t308 + t306 * t322 - t318 * t335 + t389;
t317 = -mrSges(4,2) * t335 + mrSges(4,3) * t322;
t197 = m(4) * t267 - mrSges(4,1) * t308 + mrSges(4,2) * t309 - t317 * t322 + t318 * t323 + t198;
t207 = t210 * t361 + t211 * t356;
t368 = -m(5) * t247 + t283 * mrSges(5,1) - mrSges(5,2) * t284 + t315 * t299 - t300 * t316 - t207;
t203 = m(4) * t256 + mrSges(4,1) * t332 - mrSges(4,3) * t309 - t306 * t323 + t317 * t335 + t368;
t184 = t195 * t402 + t353 * t197 + t203 * t401;
t185 = t195 * t397 - t350 * t197 + t203 * t396;
t313 = -t349 * t338 + t384;
t387 = -mrSges(3,1) * t352 + mrSges(3,2) * t349;
t336 = t387 * t395;
t381 = -mrSges(3,2) * t354 + mrSges(3,3) * t400;
t341 = t381 * qJD(1);
t382 = mrSges(3,1) * t354 - mrSges(3,3) * t404;
t182 = m(3) * t313 + t382 * qJDD(1) + (-t336 * t404 + t341 * t354) * qJD(1) + t185;
t190 = t363 * t195 - t358 * t203;
t314 = -g(3) * t404 + t392;
t340 = t382 * qJD(1);
t189 = m(3) * t314 + t381 * qJDD(1) + (t336 * t400 - t340 * t354) * qJD(1) + t190;
t179 = -t182 * t349 + t352 * t189;
t324 = -t351 * t337 + t390;
t183 = m(3) * t324 + (t387 * qJDD(1) + (t340 * t349 - t341 * t352) * qJD(1)) * t351 + t184;
t174 = t182 * t398 - t183 * t351 + t189 * t403;
t386 = Ifges(3,5) * t349 + Ifges(3,6) * t352;
t250 = Ifges(7,5) * t269 + Ifges(7,6) * t268 + Ifges(7,3) * t310;
t252 = Ifges(7,1) * t269 + Ifges(7,4) * t268 + Ifges(7,5) * t310;
t213 = -mrSges(7,1) * t230 + mrSges(7,3) * t225 + Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t276 - t250 * t269 + t252 * t310;
t251 = Ifges(7,4) * t269 + Ifges(7,2) * t268 + Ifges(7,6) * t310;
t214 = mrSges(7,2) * t230 - mrSges(7,3) * t224 + Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t276 + t250 * t268 - t251 * t310;
t263 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t312;
t265 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t312;
t199 = -mrSges(6,1) * t232 + mrSges(6,3) * t229 + Ifges(6,4) * t260 + Ifges(6,2) * t259 + Ifges(6,6) * t281 - pkin(5) * t376 + pkin(12) * t388 + t360 * t213 + t355 * t214 - t297 * t263 + t312 * t265;
t264 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t312;
t200 = mrSges(6,2) * t232 - mrSges(6,3) * t228 + Ifges(6,1) * t260 + Ifges(6,4) * t259 + Ifges(6,5) * t281 - pkin(12) * t212 - t213 * t355 + t214 * t360 + t263 * t296 - t264 * t312;
t277 = Ifges(5,5) * t316 + Ifges(5,6) * t315 + Ifges(5,3) * t321;
t278 = Ifges(5,4) * t316 + Ifges(5,2) * t315 + Ifges(5,6) * t321;
t186 = mrSges(5,2) * t247 - mrSges(5,3) * t237 + Ifges(5,1) * t284 + Ifges(5,4) * t283 + Ifges(5,5) * t305 - pkin(11) * t207 - t199 * t356 + t200 * t361 + t277 * t315 - t278 * t321;
t279 = Ifges(5,1) * t316 + Ifges(5,4) * t315 + Ifges(5,5) * t321;
t372 = -mrSges(7,1) * t224 + mrSges(7,2) * t225 - Ifges(7,5) * t243 - Ifges(7,6) * t242 - Ifges(7,3) * t276 - t269 * t251 + t268 * t252;
t366 = mrSges(6,1) * t228 - mrSges(6,2) * t229 + Ifges(6,5) * t260 + Ifges(6,6) * t259 + Ifges(6,3) * t281 + pkin(5) * t212 + t297 * t264 - t296 * t265 - t372;
t191 = -mrSges(5,1) * t247 + mrSges(5,3) * t238 + Ifges(5,4) * t284 + Ifges(5,2) * t283 + Ifges(5,6) * t305 - pkin(4) * t207 - t316 * t277 + t321 * t279 - t366;
t301 = Ifges(4,5) * t323 + Ifges(4,6) * t322 + Ifges(4,3) * t335;
t302 = Ifges(4,4) * t323 + Ifges(4,2) * t322 + Ifges(4,6) * t335;
t176 = mrSges(4,2) * t267 - mrSges(4,3) * t256 + Ifges(4,1) * t309 + Ifges(4,4) * t308 + Ifges(4,5) * t332 - pkin(10) * t198 + t186 * t362 - t191 * t357 + t301 * t322 - t302 * t335;
t303 = Ifges(4,1) * t323 + Ifges(4,4) * t322 + Ifges(4,5) * t335;
t367 = mrSges(5,1) * t237 - mrSges(5,2) * t238 + Ifges(5,5) * t284 + Ifges(5,6) * t283 + Ifges(5,3) * t305 + pkin(4) * t222 + pkin(11) * t208 + t361 * t199 + t356 * t200 + t316 * t278 - t315 * t279;
t180 = -mrSges(4,1) * t267 + mrSges(4,3) * t257 + Ifges(4,4) * t309 + Ifges(4,2) * t308 + Ifges(4,6) * t332 - pkin(3) * t198 - t323 * t301 + t335 * t303 - t367;
t377 = pkin(9) * t190 + t176 * t358 + t180 * t363;
t375 = Ifges(3,5) * t354 + (Ifges(3,1) * t349 + Ifges(3,4) * t352) * t351;
t374 = Ifges(3,6) * t354 + (Ifges(3,4) * t349 + Ifges(3,2) * t352) * t351;
t175 = mrSges(4,1) * t256 - mrSges(4,2) * t257 + Ifges(4,5) * t309 + Ifges(4,6) * t308 + Ifges(4,3) * t332 + pkin(3) * t368 + pkin(10) * t389 + t357 * t186 + t362 * t191 + t323 * t302 - t322 * t303;
t328 = t374 * qJD(1);
t329 = t375 * qJD(1);
t166 = qJDD(1) * t407 + mrSges(3,1) * t313 - mrSges(3,2) * t314 + pkin(2) * t185 + t353 * t175 + t377 * t350 + (t386 * qJDD(1) + (t328 * t349 - t329 * t352) * qJD(1)) * t351;
t327 = (t351 * t386 + t407) * qJD(1);
t168 = -mrSges(3,1) * t324 + mrSges(3,3) * t314 - pkin(2) * t184 - t350 * t175 + (-t327 * t404 + t329 * t354) * qJD(1) + t377 * t353 + t374 * qJDD(1);
t170 = mrSges(3,2) * t324 - mrSges(3,3) * t313 + t363 * t176 - t358 * t180 + (t327 * t400 - t328 * t354) * qJD(1) + (-t184 * t350 - t185 * t353) * pkin(9) + t375 * qJDD(1);
t373 = mrSges(2,1) * t345 - mrSges(2,2) * t346 + Ifges(2,3) * qJDD(1) + pkin(1) * t174 + t354 * t166 + t168 * t400 + t170 * t404 + t179 * t405;
t177 = m(2) * t346 - mrSges(2,1) * t365 - qJDD(1) * mrSges(2,2) + t179;
t173 = t354 * t183 + (t182 * t352 + t189 * t349) * t351;
t171 = m(2) * t345 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t365 + t174;
t164 = -mrSges(2,2) * g(3) - mrSges(2,3) * t345 + Ifges(2,5) * qJDD(1) - t365 * Ifges(2,6) - t349 * t168 + t352 * t170 + (-t173 * t351 - t174 * t354) * qJ(2);
t163 = mrSges(2,1) * g(3) + mrSges(2,3) * t346 + t365 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t173 - t351 * t166 + (qJ(2) * t179 + t168 * t352 + t170 * t349) * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t364 * t164 - t359 * t163 - pkin(8) * (t171 * t364 + t177 * t359), t164, t170, t176, t186, t200, t214; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t359 * t164 + t364 * t163 + pkin(8) * (-t171 * t359 + t177 * t364), t163, t168, t180, t191, t199, t213; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t373, t373, t166, t175, t367, t366, -t372;];
m_new  = t1;
