% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 19:44
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 19:29:52
% EndTime: 2019-05-05 19:32:38
% DurationCPUTime: 162.66s
% Computational Cost: add. (2524041->400), mult. (8413941->551), div. (0->0), fcn. (7239919->16), ass. (0->183)
t351 = sin(pkin(7));
t354 = cos(pkin(12));
t356 = cos(pkin(6));
t352 = sin(pkin(6));
t355 = cos(pkin(7));
t396 = t352 * t355;
t334 = (t351 * t356 + t354 * t396) * qJD(1) * pkin(9);
t360 = sin(qJ(1));
t364 = cos(qJ(1));
t346 = -t364 * g(1) - t360 * g(2);
t365 = qJD(1) ^ 2;
t402 = qJ(2) * t352;
t338 = -t365 * pkin(1) + qJDD(1) * t402 + t346;
t350 = sin(pkin(12));
t405 = pkin(9) * t350;
t381 = -pkin(2) * t354 - t351 * t405;
t392 = qJD(1) * t352;
t403 = pkin(9) * qJDD(1);
t377 = qJD(1) * t381 * t392 + t355 * t403;
t345 = t360 * g(1) - t364 * g(2);
t337 = qJDD(1) * pkin(1) + t365 * t402 + t345;
t389 = qJD(2) * t392;
t395 = t354 * t356;
t397 = t352 * t354;
t382 = -g(3) * t397 + t337 * t395 - 0.2e1 * t350 * t389;
t291 = (pkin(2) * qJDD(1) + qJD(1) * t334) * t356 + (-t377 * t352 - t338) * t350 + t382;
t339 = (pkin(2) * t356 - t396 * t405) * qJD(1);
t400 = t350 * t356;
t390 = t337 * t400 + (t338 + 0.2e1 * t389) * t354;
t292 = (-qJD(1) * t339 + t351 * t403) * t356 + (-g(3) * t350 + t377 * t354) * t352 + t390;
t388 = -t356 * g(3) + qJDD(2);
t299 = (-t337 + t381 * qJDD(1) + (-t334 * t354 + t339 * t350) * qJD(1)) * t352 + t388;
t359 = sin(qJ(3));
t363 = cos(qJ(3));
t393 = t355 * t363;
t398 = t351 * t363;
t254 = t291 * t393 - t359 * t292 + t299 * t398;
t370 = t356 * t398 + (-t350 * t359 + t354 * t393) * t352;
t321 = t370 * qJD(1);
t394 = t355 * t359;
t399 = t351 * t359;
t371 = t356 * t399 + (t350 * t363 + t354 * t394) * t352;
t314 = t321 * qJD(3) + t371 * qJDD(1);
t322 = t371 * qJD(1);
t378 = -t351 * t397 + t355 * t356;
t332 = t378 * qJDD(1) + qJDD(3);
t335 = t378 * qJD(1) + qJD(3);
t244 = (t321 * t335 - t314) * qJ(4) + (t321 * t322 + t332) * pkin(3) + t254;
t255 = t291 * t394 + t363 * t292 + t299 * t399;
t313 = -t322 * qJD(3) + t370 * qJDD(1);
t318 = t335 * pkin(3) - t322 * qJ(4);
t320 = t321 ^ 2;
t247 = -t320 * pkin(3) + t313 * qJ(4) - t335 * t318 + t255;
t349 = sin(pkin(13));
t353 = cos(pkin(13));
t311 = t349 * t321 + t353 * t322;
t236 = -0.2e1 * qJD(4) * t311 + t353 * t244 - t349 * t247;
t404 = Ifges(3,3) * t356;
t401 = t350 * t352;
t310 = t353 * t321 - t349 * t322;
t237 = 0.2e1 * qJD(4) * t310 + t349 * t244 + t353 * t247;
t280 = -t310 * mrSges(5,1) + t311 * mrSges(5,2);
t286 = t353 * t313 - t349 * t314;
t301 = t335 * mrSges(5,1) - t311 * mrSges(5,3);
t281 = -t310 * pkin(4) - t311 * pkin(10);
t331 = t335 ^ 2;
t234 = -t331 * pkin(4) + t332 * pkin(10) + t310 * t281 + t237;
t268 = -t351 * t291 + t355 * t299;
t253 = -t313 * pkin(3) - t320 * qJ(4) + t322 * t318 + qJDD(4) + t268;
t287 = t349 * t313 + t353 * t314;
t239 = (-t310 * t335 - t287) * pkin(10) + (t311 * t335 - t286) * pkin(4) + t253;
t358 = sin(qJ(5));
t362 = cos(qJ(5));
t230 = t362 * t234 + t358 * t239;
t297 = -t358 * t311 + t362 * t335;
t298 = t362 * t311 + t358 * t335;
t271 = -t297 * pkin(5) - t298 * pkin(11);
t285 = qJDD(5) - t286;
t309 = qJD(5) - t310;
t308 = t309 ^ 2;
t228 = -t308 * pkin(5) + t285 * pkin(11) + t297 * t271 + t230;
t233 = -t332 * pkin(4) - t331 * pkin(10) + t311 * t281 - t236;
t265 = -t298 * qJD(5) - t358 * t287 + t362 * t332;
t266 = t297 * qJD(5) + t362 * t287 + t358 * t332;
t231 = (-t297 * t309 - t266) * pkin(11) + (t298 * t309 - t265) * pkin(5) + t233;
t357 = sin(qJ(6));
t361 = cos(qJ(6));
t224 = -t357 * t228 + t361 * t231;
t272 = -t357 * t298 + t361 * t309;
t242 = t272 * qJD(6) + t361 * t266 + t357 * t285;
t273 = t361 * t298 + t357 * t309;
t256 = -t272 * mrSges(7,1) + t273 * mrSges(7,2);
t296 = qJD(6) - t297;
t257 = -t296 * mrSges(7,2) + t272 * mrSges(7,3);
t264 = qJDD(6) - t265;
t222 = m(7) * t224 + t264 * mrSges(7,1) - t242 * mrSges(7,3) - t273 * t256 + t296 * t257;
t225 = t361 * t228 + t357 * t231;
t241 = -t273 * qJD(6) - t357 * t266 + t361 * t285;
t258 = t296 * mrSges(7,1) - t273 * mrSges(7,3);
t223 = m(7) * t225 - t264 * mrSges(7,2) + t241 * mrSges(7,3) + t272 * t256 - t296 * t258;
t216 = -t357 * t222 + t361 * t223;
t270 = -t297 * mrSges(6,1) + t298 * mrSges(6,2);
t275 = t309 * mrSges(6,1) - t298 * mrSges(6,3);
t213 = m(6) * t230 - t285 * mrSges(6,2) + t265 * mrSges(6,3) + t297 * t270 - t309 * t275 + t216;
t229 = -t358 * t234 + t362 * t239;
t227 = -t285 * pkin(5) - t308 * pkin(11) + t298 * t271 - t229;
t226 = -m(7) * t227 + t241 * mrSges(7,1) - t242 * mrSges(7,2) + t272 * t257 - t273 * t258;
t274 = -t309 * mrSges(6,2) + t297 * mrSges(6,3);
t220 = m(6) * t229 + t285 * mrSges(6,1) - t266 * mrSges(6,3) - t298 * t270 + t309 * t274 + t226;
t386 = t362 * t213 - t358 * t220;
t204 = m(5) * t237 - t332 * mrSges(5,2) + t286 * mrSges(5,3) + t310 * t280 - t335 * t301 + t386;
t300 = -t335 * mrSges(5,2) + t310 * mrSges(5,3);
t215 = t361 * t222 + t357 * t223;
t368 = -m(6) * t233 + t265 * mrSges(6,1) - t266 * mrSges(6,2) + t297 * t274 - t298 * t275 - t215;
t210 = m(5) * t236 + t332 * mrSges(5,1) - t287 * mrSges(5,3) - t311 * t280 + t335 * t300 + t368;
t197 = t349 * t204 + t353 * t210;
t208 = t358 * t213 + t362 * t220;
t312 = -t321 * mrSges(4,1) + t322 * mrSges(4,2);
t317 = -t335 * mrSges(4,2) + t321 * mrSges(4,3);
t195 = m(4) * t254 + t332 * mrSges(4,1) - t314 * mrSges(4,3) - t322 * t312 + t335 * t317 + t197;
t319 = t335 * mrSges(4,1) - t322 * mrSges(4,3);
t387 = t353 * t204 - t349 * t210;
t196 = m(4) * t255 - t332 * mrSges(4,2) + t313 * mrSges(4,3) + t321 * t312 - t335 * t319 + t387;
t372 = m(5) * t253 - t286 * mrSges(5,1) + t287 * mrSges(5,2) - t310 * t300 + t311 * t301 + t208;
t206 = m(4) * t268 - t313 * mrSges(4,1) + t314 * mrSges(4,2) - t321 * t317 + t322 * t319 + t372;
t182 = t195 * t398 + t196 * t399 + t355 * t206;
t183 = t195 * t393 + t196 * t394 - t351 * t206;
t315 = -t350 * t338 + t382;
t384 = -mrSges(3,1) * t354 + mrSges(3,2) * t350;
t336 = t384 * t392;
t379 = -mrSges(3,2) * t356 + mrSges(3,3) * t397;
t341 = t379 * qJD(1);
t380 = mrSges(3,1) * t356 - mrSges(3,3) * t401;
t180 = m(3) * t315 + t380 * qJDD(1) + (-t336 * t401 + t341 * t356) * qJD(1) + t183;
t188 = -t359 * t195 + t363 * t196;
t316 = -g(3) * t401 + t390;
t340 = t380 * qJD(1);
t187 = m(3) * t316 + t379 * qJDD(1) + (t336 * t397 - t340 * t356) * qJD(1) + t188;
t177 = -t350 * t180 + t354 * t187;
t323 = -t352 * t337 + t388;
t181 = m(3) * t323 + (t384 * qJDD(1) + (t340 * t350 - t341 * t354) * qJD(1)) * t352 + t182;
t172 = t180 * t395 - t352 * t181 + t187 * t400;
t383 = Ifges(3,5) * t350 + Ifges(3,6) * t354;
t248 = Ifges(7,5) * t273 + Ifges(7,6) * t272 + Ifges(7,3) * t296;
t250 = Ifges(7,1) * t273 + Ifges(7,4) * t272 + Ifges(7,5) * t296;
t217 = -mrSges(7,1) * t227 + mrSges(7,3) * t225 + Ifges(7,4) * t242 + Ifges(7,2) * t241 + Ifges(7,6) * t264 - t273 * t248 + t296 * t250;
t249 = Ifges(7,4) * t273 + Ifges(7,2) * t272 + Ifges(7,6) * t296;
t218 = mrSges(7,2) * t227 - mrSges(7,3) * t224 + Ifges(7,1) * t242 + Ifges(7,4) * t241 + Ifges(7,5) * t264 + t272 * t248 - t296 * t249;
t259 = Ifges(6,5) * t298 + Ifges(6,6) * t297 + Ifges(6,3) * t309;
t260 = Ifges(6,4) * t298 + Ifges(6,2) * t297 + Ifges(6,6) * t309;
t199 = mrSges(6,2) * t233 - mrSges(6,3) * t229 + Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t285 - pkin(11) * t215 - t357 * t217 + t361 * t218 + t297 * t259 - t309 * t260;
t261 = Ifges(6,1) * t298 + Ifges(6,4) * t297 + Ifges(6,5) * t309;
t367 = mrSges(7,1) * t224 - mrSges(7,2) * t225 + Ifges(7,5) * t242 + Ifges(7,6) * t241 + Ifges(7,3) * t264 + t273 * t249 - t272 * t250;
t201 = -mrSges(6,1) * t233 + mrSges(6,3) * t230 + Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t285 - pkin(5) * t215 - t298 * t259 + t309 * t261 - t367;
t276 = Ifges(5,5) * t311 + Ifges(5,6) * t310 + Ifges(5,3) * t335;
t277 = Ifges(5,4) * t311 + Ifges(5,2) * t310 + Ifges(5,6) * t335;
t184 = mrSges(5,2) * t253 - mrSges(5,3) * t236 + Ifges(5,1) * t287 + Ifges(5,4) * t286 + Ifges(5,5) * t332 - pkin(10) * t208 + t362 * t199 - t358 * t201 + t310 * t276 - t335 * t277;
t278 = Ifges(5,1) * t311 + Ifges(5,4) * t310 + Ifges(5,5) * t335;
t366 = mrSges(6,1) * t229 - mrSges(6,2) * t230 + Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t285 + pkin(5) * t226 + pkin(11) * t216 + t361 * t217 + t357 * t218 + t298 * t260 - t297 * t261;
t189 = -mrSges(5,1) * t253 + mrSges(5,3) * t237 + Ifges(5,4) * t287 + Ifges(5,2) * t286 + Ifges(5,6) * t332 - pkin(4) * t208 - t311 * t276 + t335 * t278 - t366;
t302 = Ifges(4,5) * t322 + Ifges(4,6) * t321 + Ifges(4,3) * t335;
t304 = Ifges(4,1) * t322 + Ifges(4,4) * t321 + Ifges(4,5) * t335;
t173 = -mrSges(4,1) * t268 + mrSges(4,3) * t255 + Ifges(4,4) * t314 + Ifges(4,2) * t313 + Ifges(4,6) * t332 - pkin(3) * t372 + qJ(4) * t387 + t349 * t184 + t353 * t189 - t322 * t302 + t335 * t304;
t303 = Ifges(4,4) * t322 + Ifges(4,2) * t321 + Ifges(4,6) * t335;
t174 = mrSges(4,2) * t268 - mrSges(4,3) * t254 + Ifges(4,1) * t314 + Ifges(4,4) * t313 + Ifges(4,5) * t332 - qJ(4) * t197 + t353 * t184 - t349 * t189 + t321 * t302 - t335 * t303;
t376 = pkin(9) * t188 + t173 * t363 + t174 * t359;
t375 = Ifges(3,5) * t356 + (Ifges(3,1) * t350 + Ifges(3,4) * t354) * t352;
t374 = Ifges(3,6) * t356 + (Ifges(3,4) * t350 + Ifges(3,2) * t354) * t352;
t369 = mrSges(5,1) * t236 - mrSges(5,2) * t237 + Ifges(5,5) * t287 + Ifges(5,6) * t286 + Ifges(5,3) * t332 + pkin(4) * t368 + pkin(10) * t386 + t358 * t199 + t362 * t201 + t311 * t277 - t310 * t278;
t178 = mrSges(4,1) * t254 - mrSges(4,2) * t255 + Ifges(4,5) * t314 + Ifges(4,6) * t313 + Ifges(4,3) * t332 + pkin(3) * t197 + t322 * t303 - t321 * t304 + t369;
t327 = t374 * qJD(1);
t328 = t375 * qJD(1);
t164 = qJDD(1) * t404 + mrSges(3,1) * t315 - mrSges(3,2) * t316 + pkin(2) * t183 + t355 * t178 + t376 * t351 + (t383 * qJDD(1) + (t327 * t350 - t328 * t354) * qJD(1)) * t352;
t326 = (t383 * t352 + t404) * qJD(1);
t166 = -mrSges(3,1) * t323 + mrSges(3,3) * t316 - pkin(2) * t182 - t351 * t178 + (-t326 * t401 + t328 * t356) * qJD(1) + t376 * t355 + t374 * qJDD(1);
t168 = mrSges(3,2) * t323 - mrSges(3,3) * t315 - t359 * t173 + t363 * t174 + (t326 * t397 - t327 * t356) * qJD(1) + (-t182 * t351 - t183 * t355) * pkin(9) + t375 * qJDD(1);
t373 = mrSges(2,1) * t345 - mrSges(2,2) * t346 + Ifges(2,3) * qJDD(1) + pkin(1) * t172 + t356 * t164 + t166 * t397 + t168 * t401 + t177 * t402;
t175 = m(2) * t346 - t365 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t177;
t171 = t356 * t181 + (t180 * t354 + t187 * t350) * t352;
t169 = m(2) * t345 + qJDD(1) * mrSges(2,1) - t365 * mrSges(2,2) + t172;
t162 = -mrSges(2,2) * g(3) - mrSges(2,3) * t345 + Ifges(2,5) * qJDD(1) - t365 * Ifges(2,6) - t350 * t166 + t354 * t168 + (-t171 * t352 - t172 * t356) * qJ(2);
t161 = mrSges(2,1) * g(3) + mrSges(2,3) * t346 + t365 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t171 - t352 * t164 + (qJ(2) * t177 + t166 * t354 + t168 * t350) * t356;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t364 * t162 - t360 * t161 - pkin(8) * (t364 * t169 + t360 * t175), t162, t168, t174, t184, t199, t218; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t360 * t162 + t364 * t161 + pkin(8) * (-t360 * t169 + t364 * t175), t161, t166, t173, t189, t201, t217; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t373, t373, t164, t178, t369, t366, t367;];
m_new  = t1;
