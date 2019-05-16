% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 23:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:19:30
% EndTime: 2019-05-05 23:22:27
% DurationCPUTime: 173.15s
% Computational Cost: add. (2760131->402), mult. (8662945->552), div. (0->0), fcn. (7455718->16), ass. (0->185)
t350 = sin(pkin(12));
t352 = sin(pkin(6));
t354 = cos(pkin(12));
t356 = cos(pkin(6));
t359 = sin(qJ(3));
t355 = cos(pkin(7));
t363 = cos(qJ(3));
t399 = t355 * t363;
t351 = sin(pkin(7));
t404 = t351 * t363;
t370 = (-t350 * t359 + t354 * t399) * t352 + t356 * t404;
t322 = t370 * qJD(1);
t400 = t355 * t359;
t405 = t351 * t359;
t372 = t356 * t405 + (t350 * t363 + t354 * t400) * t352;
t323 = t372 * qJD(1);
t311 = -t323 * qJD(3) + t370 * qJDD(1);
t402 = t352 * t355;
t334 = (t351 * t356 + t354 * t402) * qJD(1) * pkin(9);
t360 = sin(qJ(1));
t364 = cos(qJ(1));
t346 = -t364 * g(1) - t360 * g(2);
t365 = qJD(1) ^ 2;
t408 = qJ(2) * t352;
t338 = -t365 * pkin(1) + qJDD(1) * t408 + t346;
t411 = pkin(9) * t350;
t384 = -pkin(2) * t354 - t351 * t411;
t398 = qJD(1) * t352;
t409 = pkin(9) * qJDD(1);
t379 = qJD(1) * t384 * t398 + t355 * t409;
t345 = t360 * g(1) - t364 * g(2);
t337 = qJDD(1) * pkin(1) + t365 * t408 + t345;
t394 = qJD(2) * t398;
t401 = t354 * t356;
t403 = t352 * t354;
t385 = -g(3) * t403 + t337 * t401 - 0.2e1 * t350 * t394;
t290 = (pkin(2) * qJDD(1) + qJD(1) * t334) * t356 + (-t379 * t352 - t338) * t350 + t385;
t339 = (pkin(2) * t356 - t402 * t411) * qJD(1);
t406 = t350 * t356;
t395 = t337 * t406 + (t338 + 0.2e1 * t394) * t354;
t291 = (-qJD(1) * t339 + t351 * t409) * t356 + (-g(3) * t350 + t379 * t354) * t352 + t395;
t393 = -t356 * g(3) + qJDD(2);
t300 = (-t337 + t384 * qJDD(1) + (-t334 * t354 + t339 * t350) * qJD(1)) * t352 + t393;
t257 = -t359 * t291 + (t290 * t355 + t300 * t351) * t363;
t412 = 2 * qJD(5);
t410 = Ifges(3,3) * t356;
t407 = t350 * t352;
t258 = t290 * t400 + t363 * t291 + t300 * t405;
t310 = -t322 * pkin(3) - t323 * pkin(10);
t380 = -t351 * t403 + t355 * t356;
t335 = t380 * qJD(1) + qJD(3);
t331 = t335 ^ 2;
t332 = t380 * qJDD(1) + qJDD(3);
t248 = -t331 * pkin(3) + t332 * pkin(10) + t322 * t310 + t258;
t274 = -t351 * t290 + t355 * t300;
t312 = t322 * qJD(3) + t372 * qJDD(1);
t251 = (-t322 * t335 - t312) * pkin(10) + (t323 * t335 - t311) * pkin(3) + t274;
t358 = sin(qJ(4));
t362 = cos(qJ(4));
t240 = -t358 * t248 + t362 * t251;
t316 = -t358 * t323 + t362 * t335;
t286 = t316 * qJD(4) + t362 * t312 + t358 * t332;
t308 = qJDD(4) - t311;
t317 = t362 * t323 + t358 * t335;
t321 = qJD(4) - t322;
t237 = (t316 * t321 - t286) * qJ(5) + (t316 * t317 + t308) * pkin(4) + t240;
t241 = t362 * t248 + t358 * t251;
t285 = -t317 * qJD(4) - t358 * t312 + t362 * t332;
t302 = t321 * pkin(4) - t317 * qJ(5);
t315 = t316 ^ 2;
t239 = -t315 * pkin(4) + t285 * qJ(5) - t321 * t302 + t241;
t349 = sin(pkin(13));
t353 = cos(pkin(13));
t294 = t353 * t316 - t349 * t317;
t234 = t349 * t237 + t353 * t239 + t294 * t412;
t264 = t353 * t285 - t349 * t286;
t295 = t349 * t316 + t353 * t317;
t272 = -t294 * mrSges(6,1) + t295 * mrSges(6,2);
t278 = t321 * mrSges(6,1) - t295 * mrSges(6,3);
t273 = -t294 * pkin(5) - t295 * pkin(11);
t320 = t321 ^ 2;
t231 = -t320 * pkin(5) + t308 * pkin(11) + t294 * t273 + t234;
t247 = -t332 * pkin(3) - t331 * pkin(10) + t323 * t310 - t257;
t242 = -t285 * pkin(4) - t315 * qJ(5) + t317 * t302 + qJDD(5) + t247;
t265 = t349 * t285 + t353 * t286;
t235 = (-t294 * t321 - t265) * pkin(11) + (t295 * t321 - t264) * pkin(5) + t242;
t357 = sin(qJ(6));
t361 = cos(qJ(6));
t228 = -t357 * t231 + t361 * t235;
t275 = -t357 * t295 + t361 * t321;
t245 = t275 * qJD(6) + t361 * t265 + t357 * t308;
t276 = t361 * t295 + t357 * t321;
t259 = -t275 * mrSges(7,1) + t276 * mrSges(7,2);
t263 = qJDD(6) - t264;
t293 = qJD(6) - t294;
t266 = -t293 * mrSges(7,2) + t275 * mrSges(7,3);
t224 = m(7) * t228 + t263 * mrSges(7,1) - t245 * mrSges(7,3) - t276 * t259 + t293 * t266;
t229 = t361 * t231 + t357 * t235;
t244 = -t276 * qJD(6) - t357 * t265 + t361 * t308;
t267 = t293 * mrSges(7,1) - t276 * mrSges(7,3);
t225 = m(7) * t229 - t263 * mrSges(7,2) + t244 * mrSges(7,3) + t275 * t259 - t293 * t267;
t390 = -t357 * t224 + t361 * t225;
t210 = m(6) * t234 - t308 * mrSges(6,2) + t264 * mrSges(6,3) + t294 * t272 - t321 * t278 + t390;
t387 = -t353 * t237 + t349 * t239;
t233 = -0.2e1 * qJD(5) * t295 - t387;
t277 = -t321 * mrSges(6,2) + t294 * mrSges(6,3);
t230 = -t308 * pkin(5) - t320 * pkin(11) + (t412 + t273) * t295 + t387;
t377 = -m(7) * t230 + t244 * mrSges(7,1) - t245 * mrSges(7,2) + t275 * t266 - t276 * t267;
t220 = m(6) * t233 + t308 * mrSges(6,1) - t265 * mrSges(6,3) - t295 * t272 + t321 * t277 + t377;
t205 = t349 * t210 + t353 * t220;
t296 = -t316 * mrSges(5,1) + t317 * mrSges(5,2);
t301 = -t321 * mrSges(5,2) + t316 * mrSges(5,3);
t203 = m(5) * t240 + t308 * mrSges(5,1) - t286 * mrSges(5,3) - t317 * t296 + t321 * t301 + t205;
t303 = t321 * mrSges(5,1) - t317 * mrSges(5,3);
t391 = t353 * t210 - t349 * t220;
t204 = m(5) * t241 - t308 * mrSges(5,2) + t285 * mrSges(5,3) + t316 * t296 - t321 * t303 + t391;
t197 = t362 * t203 + t358 * t204;
t213 = t361 * t224 + t357 * t225;
t309 = -t322 * mrSges(4,1) + t323 * mrSges(4,2);
t319 = t335 * mrSges(4,1) - t323 * mrSges(4,3);
t392 = -t358 * t203 + t362 * t204;
t194 = m(4) * t258 - t332 * mrSges(4,2) + t311 * mrSges(4,3) + t322 * t309 - t335 * t319 + t392;
t318 = -t335 * mrSges(4,2) + t322 * mrSges(4,3);
t196 = m(4) * t274 - t311 * mrSges(4,1) + t312 * mrSges(4,2) - t322 * t318 + t323 * t319 + t197;
t373 = m(6) * t242 - t264 * mrSges(6,1) + t265 * mrSges(6,2) - t294 * t277 + t295 * t278 + t213;
t367 = -m(5) * t247 + t285 * mrSges(5,1) - t286 * mrSges(5,2) + t316 * t301 - t317 * t303 - t373;
t211 = m(4) * t257 + t332 * mrSges(4,1) - t312 * mrSges(4,3) - t323 * t309 + t335 * t318 + t367;
t183 = t194 * t405 + t355 * t196 + t211 * t404;
t184 = t194 * t400 - t351 * t196 + t211 * t399;
t313 = -t350 * t338 + t385;
t389 = -mrSges(3,1) * t354 + mrSges(3,2) * t350;
t336 = t389 * t398;
t382 = -mrSges(3,2) * t356 + mrSges(3,3) * t403;
t341 = t382 * qJD(1);
t383 = mrSges(3,1) * t356 - mrSges(3,3) * t407;
t181 = m(3) * t313 + t383 * qJDD(1) + (-t336 * t407 + t341 * t356) * qJD(1) + t184;
t190 = t363 * t194 - t359 * t211;
t314 = -g(3) * t407 + t395;
t340 = t383 * qJD(1);
t189 = m(3) * t314 + t382 * qJDD(1) + (t336 * t403 - t340 * t356) * qJD(1) + t190;
t178 = -t350 * t181 + t354 * t189;
t324 = -t352 * t337 + t393;
t182 = m(3) * t324 + (t389 * qJDD(1) + (t340 * t350 - t341 * t354) * qJD(1)) * t352 + t183;
t173 = t181 * t401 - t352 * t182 + t189 * t406;
t388 = Ifges(3,5) * t350 + Ifges(3,6) * t354;
t252 = Ifges(7,5) * t276 + Ifges(7,6) * t275 + Ifges(7,3) * t293;
t254 = Ifges(7,1) * t276 + Ifges(7,4) * t275 + Ifges(7,5) * t293;
t217 = -mrSges(7,1) * t230 + mrSges(7,3) * t229 + Ifges(7,4) * t245 + Ifges(7,2) * t244 + Ifges(7,6) * t263 - t276 * t252 + t293 * t254;
t253 = Ifges(7,4) * t276 + Ifges(7,2) * t275 + Ifges(7,6) * t293;
t218 = mrSges(7,2) * t230 - mrSges(7,3) * t228 + Ifges(7,1) * t245 + Ifges(7,4) * t244 + Ifges(7,5) * t263 + t275 * t252 - t293 * t253;
t268 = Ifges(6,5) * t295 + Ifges(6,6) * t294 + Ifges(6,3) * t321;
t269 = Ifges(6,4) * t295 + Ifges(6,2) * t294 + Ifges(6,6) * t321;
t198 = mrSges(6,2) * t242 - mrSges(6,3) * t233 + Ifges(6,1) * t265 + Ifges(6,4) * t264 + Ifges(6,5) * t308 - pkin(11) * t213 - t357 * t217 + t361 * t218 + t294 * t268 - t321 * t269;
t270 = Ifges(6,1) * t295 + Ifges(6,4) * t294 + Ifges(6,5) * t321;
t368 = mrSges(7,1) * t228 - mrSges(7,2) * t229 + Ifges(7,5) * t245 + Ifges(7,6) * t244 + Ifges(7,3) * t263 + t276 * t253 - t275 * t254;
t199 = -mrSges(6,1) * t242 + mrSges(6,3) * t234 + Ifges(6,4) * t265 + Ifges(6,2) * t264 + Ifges(6,6) * t308 - pkin(5) * t213 - t295 * t268 + t321 * t270 - t368;
t279 = Ifges(5,5) * t317 + Ifges(5,6) * t316 + Ifges(5,3) * t321;
t281 = Ifges(5,1) * t317 + Ifges(5,4) * t316 + Ifges(5,5) * t321;
t185 = -mrSges(5,1) * t247 + mrSges(5,3) * t241 + Ifges(5,4) * t286 + Ifges(5,2) * t285 + Ifges(5,6) * t308 - pkin(4) * t373 + qJ(5) * t391 + t349 * t198 + t353 * t199 - t317 * t279 + t321 * t281;
t280 = Ifges(5,4) * t317 + Ifges(5,2) * t316 + Ifges(5,6) * t321;
t186 = mrSges(5,2) * t247 - mrSges(5,3) * t240 + Ifges(5,1) * t286 + Ifges(5,4) * t285 + Ifges(5,5) * t308 - qJ(5) * t205 + t353 * t198 - t349 * t199 + t316 * t279 - t321 * t280;
t304 = Ifges(4,5) * t323 + Ifges(4,6) * t322 + Ifges(4,3) * t335;
t305 = Ifges(4,4) * t323 + Ifges(4,2) * t322 + Ifges(4,6) * t335;
t175 = mrSges(4,2) * t274 - mrSges(4,3) * t257 + Ifges(4,1) * t312 + Ifges(4,4) * t311 + Ifges(4,5) * t332 - pkin(10) * t197 - t358 * t185 + t362 * t186 + t322 * t304 - t335 * t305;
t306 = Ifges(4,1) * t323 + Ifges(4,4) * t322 + Ifges(4,5) * t335;
t369 = -mrSges(6,1) * t233 + mrSges(6,2) * t234 - Ifges(6,5) * t265 - Ifges(6,6) * t264 - Ifges(6,3) * t308 - pkin(5) * t377 - pkin(11) * t390 - t361 * t217 - t357 * t218 - t295 * t269 + t294 * t270;
t366 = mrSges(5,1) * t240 - mrSges(5,2) * t241 + Ifges(5,5) * t286 + Ifges(5,6) * t285 + Ifges(5,3) * t308 + pkin(4) * t205 + t317 * t280 - t316 * t281 - t369;
t179 = -mrSges(4,1) * t274 + mrSges(4,3) * t258 + Ifges(4,4) * t312 + Ifges(4,2) * t311 + Ifges(4,6) * t332 - pkin(3) * t197 - t323 * t304 + t335 * t306 - t366;
t378 = pkin(9) * t190 + t175 * t359 + t179 * t363;
t376 = Ifges(3,5) * t356 + (Ifges(3,1) * t350 + Ifges(3,4) * t354) * t352;
t375 = Ifges(3,6) * t356 + (Ifges(3,4) * t350 + Ifges(3,2) * t354) * t352;
t174 = mrSges(4,1) * t257 - mrSges(4,2) * t258 + Ifges(4,5) * t312 + Ifges(4,6) * t311 + Ifges(4,3) * t332 + pkin(3) * t367 + pkin(10) * t392 + t362 * t185 + t358 * t186 + t323 * t305 - t322 * t306;
t328 = t375 * qJD(1);
t329 = t376 * qJD(1);
t165 = qJDD(1) * t410 + mrSges(3,1) * t313 - mrSges(3,2) * t314 + pkin(2) * t184 + t355 * t174 + t378 * t351 + (t388 * qJDD(1) + (t328 * t350 - t329 * t354) * qJD(1)) * t352;
t327 = (t388 * t352 + t410) * qJD(1);
t167 = -mrSges(3,1) * t324 + mrSges(3,3) * t314 - pkin(2) * t183 - t351 * t174 + (-t327 * t407 + t329 * t356) * qJD(1) + t378 * t355 + t375 * qJDD(1);
t169 = mrSges(3,2) * t324 - mrSges(3,3) * t313 + t363 * t175 - t359 * t179 + (t327 * t403 - t328 * t356) * qJD(1) + (-t183 * t351 - t184 * t355) * pkin(9) + t376 * qJDD(1);
t374 = mrSges(2,1) * t345 - mrSges(2,2) * t346 + Ifges(2,3) * qJDD(1) + pkin(1) * t173 + t356 * t165 + t167 * t403 + t169 * t407 + t178 * t408;
t176 = m(2) * t346 - t365 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t178;
t172 = t356 * t182 + (t181 * t354 + t189 * t350) * t352;
t170 = m(2) * t345 + qJDD(1) * mrSges(2,1) - t365 * mrSges(2,2) + t173;
t163 = -mrSges(2,2) * g(3) - mrSges(2,3) * t345 + Ifges(2,5) * qJDD(1) - t365 * Ifges(2,6) - t350 * t167 + t354 * t169 + (-t172 * t352 - t173 * t356) * qJ(2);
t162 = mrSges(2,1) * g(3) + mrSges(2,3) * t346 + t365 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t172 - t352 * t165 + (qJ(2) * t178 + t167 * t354 + t169 * t350) * t356;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t364 * t163 - t360 * t162 - pkin(8) * (t364 * t170 + t360 * t176), t163, t169, t175, t186, t198, t218; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t360 * t163 + t364 * t162 + pkin(8) * (-t360 * t170 + t364 * t176), t162, t167, t179, t185, t199, t217; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t374, t374, t165, t174, t366, -t369, t368;];
m_new  = t1;
