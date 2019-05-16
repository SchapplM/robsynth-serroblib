% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-05-08 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR8_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 13:33:53
% EndTime: 2019-05-08 13:39:50
% DurationCPUTime: 218.59s
% Computational Cost: add. (3920891->415), mult. (9612242->551), div. (0->0), fcn. (8203156->16), ass. (0->182)
t358 = cos(pkin(6));
t352 = qJD(1) * t358 + qJD(2);
t355 = sin(pkin(7));
t357 = cos(pkin(7));
t356 = sin(pkin(6));
t369 = cos(qJ(2));
t390 = qJD(1) * t369;
t387 = t356 * t390;
t380 = t352 * t355 + t357 * t387;
t334 = t380 * pkin(10);
t363 = sin(qJ(2));
t392 = qJD(1) * t356;
t404 = pkin(10) * t355;
t339 = (-pkin(2) * t369 - t363 * t404) * t392;
t389 = qJD(1) * qJD(2);
t345 = (qJDD(1) * t363 + t369 * t389) * t356;
t351 = qJDD(1) * t358 + qJDD(2);
t364 = sin(qJ(1));
t370 = cos(qJ(1));
t349 = t364 * g(1) - g(2) * t370;
t371 = qJD(1) ^ 2;
t405 = pkin(9) * t356;
t342 = qJDD(1) * pkin(1) + t371 * t405 + t349;
t350 = -g(1) * t370 - g(2) * t364;
t343 = -pkin(1) * t371 + qJDD(1) * t405 + t350;
t395 = t358 * t369;
t383 = t342 * t395 - t343 * t363;
t391 = qJD(1) * t363;
t403 = pkin(10) * t357;
t291 = -t345 * t403 + pkin(2) * t351 + t334 * t352 + (-g(3) * t369 - t339 * t391) * t356 + t383;
t388 = t356 * t391;
t337 = pkin(2) * t352 - t388 * t403;
t346 = (qJDD(1) * t369 - t363 * t389) * t356;
t381 = t346 * t357 + t351 * t355;
t396 = t358 * t363;
t393 = t342 * t396 + t369 * t343;
t292 = -t337 * t352 + (-g(3) * t363 + t339 * t390) * t356 + t381 * pkin(10) + t393;
t402 = g(3) * t358;
t298 = -t345 * t404 - pkin(2) * t346 - t402 + (-t342 + (-t334 * t369 + t337 * t363) * qJD(1)) * t356;
t362 = sin(qJ(3));
t368 = cos(qJ(3));
t263 = -t362 * t292 + (t291 * t357 + t298 * t355) * t368;
t324 = -t362 * t388 + t368 * t380;
t397 = t357 * t362;
t400 = t355 * t362;
t325 = t352 * t400 + (t363 * t368 + t369 * t397) * t392;
t309 = -t325 * qJD(3) - t362 * t345 + t368 * t381;
t310 = qJD(3) * t324 + t345 * t368 + t362 * t381;
t311 = -mrSges(4,1) * t324 + mrSges(4,2) * t325;
t335 = t352 * t357 - t355 * t387 + qJD(3);
t316 = -mrSges(4,2) * t335 + mrSges(4,3) * t324;
t326 = -t346 * t355 + t351 * t357 + qJDD(3);
t312 = -pkin(3) * t324 - pkin(11) * t325;
t333 = t335 ^ 2;
t247 = -pkin(3) * t326 - pkin(11) * t333 + t325 * t312 - t263;
t361 = sin(qJ(4));
t367 = cos(qJ(4));
t315 = t325 * t367 + t335 * t361;
t276 = -qJD(4) * t315 - t310 * t361 + t326 * t367;
t314 = -t325 * t361 + t335 * t367;
t277 = qJD(4) * t314 + t310 * t367 + t326 * t361;
t323 = qJD(4) - t324;
t300 = -mrSges(5,2) * t323 + mrSges(5,3) * t314;
t301 = mrSges(5,1) * t323 - mrSges(5,3) * t315;
t264 = t291 * t397 + t368 * t292 + t298 * t400;
t248 = -pkin(3) * t333 + pkin(11) * t326 + t312 * t324 + t264;
t274 = -t291 * t355 + t357 * t298;
t251 = (-t324 * t335 - t310) * pkin(11) + (t325 * t335 - t309) * pkin(3) + t274;
t240 = -t248 * t361 + t367 * t251;
t308 = qJDD(4) - t309;
t237 = (t314 * t323 - t277) * pkin(12) + (t314 * t315 + t308) * pkin(4) + t240;
t241 = t367 * t248 + t361 * t251;
t302 = pkin(4) * t323 - pkin(12) * t315;
t313 = t314 ^ 2;
t239 = -pkin(4) * t313 + pkin(12) * t276 - t302 * t323 + t241;
t360 = sin(qJ(5));
t366 = cos(qJ(5));
t234 = t360 * t237 + t366 * t239;
t293 = t314 * t366 - t315 * t360;
t294 = t314 * t360 + t315 * t366;
t273 = -pkin(5) * t293 - pkin(13) * t294;
t307 = qJDD(5) + t308;
t321 = qJD(5) + t323;
t320 = t321 ^ 2;
t231 = -pkin(5) * t320 + pkin(13) * t307 + t273 * t293 + t234;
t242 = -pkin(4) * t276 - pkin(12) * t313 + t315 * t302 + t247;
t256 = -qJD(5) * t294 + t276 * t366 - t277 * t360;
t257 = qJD(5) * t293 + t276 * t360 + t277 * t366;
t235 = (-t293 * t321 - t257) * pkin(13) + (t294 * t321 - t256) * pkin(5) + t242;
t359 = sin(qJ(6));
t365 = cos(qJ(6));
t228 = -t231 * t359 + t235 * t365;
t278 = -t294 * t359 + t321 * t365;
t245 = qJD(6) * t278 + t257 * t365 + t307 * t359;
t255 = qJDD(6) - t256;
t279 = t294 * t365 + t321 * t359;
t265 = -mrSges(7,1) * t278 + mrSges(7,2) * t279;
t290 = qJD(6) - t293;
t266 = -mrSges(7,2) * t290 + mrSges(7,3) * t278;
t224 = m(7) * t228 + mrSges(7,1) * t255 - mrSges(7,3) * t245 - t265 * t279 + t266 * t290;
t229 = t231 * t365 + t235 * t359;
t244 = -qJD(6) * t279 - t257 * t359 + t307 * t365;
t267 = mrSges(7,1) * t290 - mrSges(7,3) * t279;
t225 = m(7) * t229 - mrSges(7,2) * t255 + mrSges(7,3) * t244 + t265 * t278 - t267 * t290;
t213 = t365 * t224 + t359 * t225;
t280 = -mrSges(6,2) * t321 + mrSges(6,3) * t293;
t281 = mrSges(6,1) * t321 - mrSges(6,3) * t294;
t376 = m(6) * t242 - t256 * mrSges(6,1) + mrSges(6,2) * t257 - t293 * t280 + t281 * t294 + t213;
t373 = -m(5) * t247 + t276 * mrSges(5,1) - mrSges(5,2) * t277 + t314 * t300 - t301 * t315 - t376;
t208 = m(4) * t263 + mrSges(4,1) * t326 - mrSges(4,3) * t310 - t311 * t325 + t316 * t335 + t373;
t401 = t208 * t368;
t399 = t356 * t363;
t398 = t356 * t369;
t272 = -mrSges(6,1) * t293 + mrSges(6,2) * t294;
t384 = -t224 * t359 + t365 * t225;
t211 = m(6) * t234 - mrSges(6,2) * t307 + mrSges(6,3) * t256 + t272 * t293 - t281 * t321 + t384;
t233 = t237 * t366 - t239 * t360;
t230 = -pkin(5) * t307 - pkin(13) * t320 + t273 * t294 - t233;
t378 = -m(7) * t230 + t244 * mrSges(7,1) - mrSges(7,2) * t245 + t278 * t266 - t267 * t279;
t220 = m(6) * t233 + mrSges(6,1) * t307 - mrSges(6,3) * t257 - t272 * t294 + t280 * t321 + t378;
t205 = t360 * t211 + t366 * t220;
t295 = -mrSges(5,1) * t314 + mrSges(5,2) * t315;
t203 = m(5) * t240 + mrSges(5,1) * t308 - mrSges(5,3) * t277 - t295 * t315 + t300 * t323 + t205;
t385 = t366 * t211 - t220 * t360;
t204 = m(5) * t241 - mrSges(5,2) * t308 + mrSges(5,3) * t276 + t295 * t314 - t301 * t323 + t385;
t197 = t367 * t203 + t361 * t204;
t317 = mrSges(4,1) * t335 - mrSges(4,3) * t325;
t386 = -t361 * t203 + t367 * t204;
t194 = m(4) * t264 - mrSges(4,2) * t326 + mrSges(4,3) * t309 + t311 * t324 - t317 * t335 + t386;
t196 = m(4) * t274 - mrSges(4,1) * t309 + mrSges(4,2) * t310 - t316 * t324 + t317 * t325 + t197;
t183 = t194 * t400 + t357 * t196 + t355 * t401;
t184 = t194 * t397 - t196 * t355 + t357 * t401;
t318 = -g(3) * t398 + t383;
t341 = -mrSges(3,2) * t352 + mrSges(3,3) * t387;
t344 = (-mrSges(3,1) * t369 + mrSges(3,2) * t363) * t392;
t181 = m(3) * t318 + mrSges(3,1) * t351 - mrSges(3,3) * t345 + t341 * t352 - t344 * t388 + t184;
t190 = t368 * t194 - t208 * t362;
t319 = -g(3) * t399 + t393;
t340 = mrSges(3,1) * t352 - mrSges(3,3) * t388;
t189 = m(3) * t319 - mrSges(3,2) * t351 + mrSges(3,3) * t346 - t340 * t352 + t344 * t387 + t190;
t178 = -t181 * t363 + t369 * t189;
t330 = -t342 * t356 - t402;
t182 = m(3) * t330 - mrSges(3,1) * t346 + mrSges(3,2) * t345 + (t340 * t363 - t341 * t369) * t392 + t183;
t173 = t181 * t395 - t182 * t356 + t189 * t396;
t258 = Ifges(7,5) * t279 + Ifges(7,6) * t278 + Ifges(7,3) * t290;
t260 = Ifges(7,1) * t279 + Ifges(7,4) * t278 + Ifges(7,5) * t290;
t217 = -mrSges(7,1) * t230 + mrSges(7,3) * t229 + Ifges(7,4) * t245 + Ifges(7,2) * t244 + Ifges(7,6) * t255 - t258 * t279 + t260 * t290;
t259 = Ifges(7,4) * t279 + Ifges(7,2) * t278 + Ifges(7,6) * t290;
t218 = mrSges(7,2) * t230 - mrSges(7,3) * t228 + Ifges(7,1) * t245 + Ifges(7,4) * t244 + Ifges(7,5) * t255 + t258 * t278 - t259 * t290;
t268 = Ifges(6,5) * t294 + Ifges(6,6) * t293 + Ifges(6,3) * t321;
t269 = Ifges(6,4) * t294 + Ifges(6,2) * t293 + Ifges(6,6) * t321;
t198 = mrSges(6,2) * t242 - mrSges(6,3) * t233 + Ifges(6,1) * t257 + Ifges(6,4) * t256 + Ifges(6,5) * t307 - pkin(13) * t213 - t217 * t359 + t218 * t365 + t268 * t293 - t269 * t321;
t270 = Ifges(6,1) * t294 + Ifges(6,4) * t293 + Ifges(6,5) * t321;
t374 = mrSges(7,1) * t228 - mrSges(7,2) * t229 + Ifges(7,5) * t245 + Ifges(7,6) * t244 + Ifges(7,3) * t255 + t259 * t279 - t260 * t278;
t199 = -mrSges(6,1) * t242 + mrSges(6,3) * t234 + Ifges(6,4) * t257 + Ifges(6,2) * t256 + Ifges(6,6) * t307 - pkin(5) * t213 - t268 * t294 + t270 * t321 - t374;
t282 = Ifges(5,5) * t315 + Ifges(5,6) * t314 + Ifges(5,3) * t323;
t284 = Ifges(5,1) * t315 + Ifges(5,4) * t314 + Ifges(5,5) * t323;
t185 = -mrSges(5,1) * t247 + mrSges(5,3) * t241 + Ifges(5,4) * t277 + Ifges(5,2) * t276 + Ifges(5,6) * t308 - pkin(4) * t376 + pkin(12) * t385 + t360 * t198 + t366 * t199 - t315 * t282 + t323 * t284;
t283 = Ifges(5,4) * t315 + Ifges(5,2) * t314 + Ifges(5,6) * t323;
t186 = mrSges(5,2) * t247 - mrSges(5,3) * t240 + Ifges(5,1) * t277 + Ifges(5,4) * t276 + Ifges(5,5) * t308 - pkin(12) * t205 + t198 * t366 - t199 * t360 + t282 * t314 - t283 * t323;
t303 = Ifges(4,5) * t325 + Ifges(4,6) * t324 + Ifges(4,3) * t335;
t304 = Ifges(4,4) * t325 + Ifges(4,2) * t324 + Ifges(4,6) * t335;
t175 = mrSges(4,2) * t274 - mrSges(4,3) * t263 + Ifges(4,1) * t310 + Ifges(4,4) * t309 + Ifges(4,5) * t326 - pkin(11) * t197 - t185 * t361 + t186 * t367 + t303 * t324 - t304 * t335;
t305 = Ifges(4,1) * t325 + Ifges(4,4) * t324 + Ifges(4,5) * t335;
t375 = -mrSges(6,1) * t233 + mrSges(6,2) * t234 - Ifges(6,5) * t257 - Ifges(6,6) * t256 - Ifges(6,3) * t307 - pkin(5) * t378 - pkin(13) * t384 - t365 * t217 - t359 * t218 - t294 * t269 + t293 * t270;
t372 = mrSges(5,1) * t240 - mrSges(5,2) * t241 + Ifges(5,5) * t277 + Ifges(5,6) * t276 + Ifges(5,3) * t308 + pkin(4) * t205 + t315 * t283 - t314 * t284 - t375;
t179 = -mrSges(4,1) * t274 + mrSges(4,3) * t264 + Ifges(4,4) * t310 + Ifges(4,2) * t309 + Ifges(4,6) * t326 - pkin(3) * t197 - t325 * t303 + t335 * t305 - t372;
t379 = pkin(10) * t190 + t175 * t362 + t179 * t368;
t174 = mrSges(4,1) * t263 - mrSges(4,2) * t264 + Ifges(4,5) * t310 + Ifges(4,6) * t309 + Ifges(4,3) * t326 + pkin(3) * t373 + pkin(11) * t386 + t367 * t185 + t361 * t186 + t325 * t304 - t324 * t305;
t328 = Ifges(3,6) * t352 + (Ifges(3,4) * t363 + Ifges(3,2) * t369) * t392;
t329 = Ifges(3,5) * t352 + (Ifges(3,1) * t363 + Ifges(3,4) * t369) * t392;
t165 = mrSges(3,1) * t318 - mrSges(3,2) * t319 + Ifges(3,5) * t345 + Ifges(3,6) * t346 + Ifges(3,3) * t351 + pkin(2) * t184 + t174 * t357 + (t328 * t363 - t329 * t369) * t392 + t379 * t355;
t327 = Ifges(3,3) * t352 + (Ifges(3,5) * t363 + Ifges(3,6) * t369) * t392;
t167 = -mrSges(3,1) * t330 + mrSges(3,3) * t319 + Ifges(3,4) * t345 + Ifges(3,2) * t346 + Ifges(3,6) * t351 - pkin(2) * t183 - t174 * t355 - t327 * t388 + t329 * t352 + t357 * t379;
t169 = t327 * t387 + mrSges(3,2) * t330 - mrSges(3,3) * t318 + Ifges(3,1) * t345 + Ifges(3,4) * t346 + Ifges(3,5) * t351 + t175 * t368 - t179 * t362 - t328 * t352 + (-t183 * t355 - t184 * t357) * pkin(10);
t377 = mrSges(2,1) * t349 - mrSges(2,2) * t350 + Ifges(2,3) * qJDD(1) + pkin(1) * t173 + t358 * t165 + t167 * t398 + t169 * t399 + t178 * t405;
t176 = m(2) * t350 - mrSges(2,1) * t371 - qJDD(1) * mrSges(2,2) + t178;
t172 = t182 * t358 + (t181 * t369 + t189 * t363) * t356;
t170 = m(2) * t349 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t371 + t173;
t163 = -mrSges(2,2) * g(3) - mrSges(2,3) * t349 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t371 - t167 * t363 + t169 * t369 + (-t172 * t356 - t173 * t358) * pkin(9);
t162 = mrSges(2,1) * g(3) + mrSges(2,3) * t350 + Ifges(2,5) * t371 + Ifges(2,6) * qJDD(1) - pkin(1) * t172 - t356 * t165 + (pkin(9) * t178 + t167 * t369 + t169 * t363) * t358;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t370 * t163 - t364 * t162 - pkin(8) * (t170 * t370 + t176 * t364), t163, t169, t175, t186, t198, t218; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t364 * t163 + t370 * t162 + pkin(8) * (-t364 * t170 + t370 * t176), t162, t167, t179, t185, t199, t217; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t377, t377, t165, t174, t372, -t375, t374;];
m_new  = t1;
