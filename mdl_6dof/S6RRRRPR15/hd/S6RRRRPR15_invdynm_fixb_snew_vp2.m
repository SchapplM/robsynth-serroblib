% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6]';
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
% Datum: 2019-05-08 03:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR15_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR15_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 03:19:40
% EndTime: 2019-05-08 03:22:36
% DurationCPUTime: 78.03s
% Computational Cost: add. (1384219->420), mult. (3421391->534), div. (0->0), fcn. (2865839->14), ass. (0->176)
t354 = cos(pkin(6));
t348 = qJD(1) * t354 + qJD(2);
t351 = sin(pkin(7));
t353 = cos(pkin(7));
t352 = sin(pkin(6));
t362 = cos(qJ(2));
t387 = qJD(1) * t362;
t384 = t352 * t387;
t381 = t353 * t384;
t330 = (t348 * t351 + t381) * pkin(10);
t358 = sin(qJ(2));
t389 = qJD(1) * t352;
t404 = pkin(10) * t351;
t334 = (-pkin(2) * t362 - t358 * t404) * t389;
t386 = qJD(1) * qJD(2);
t340 = (qJDD(1) * t358 + t362 * t386) * t352;
t347 = qJDD(1) * t354 + qJDD(2);
t359 = sin(qJ(1));
t363 = cos(qJ(1));
t345 = t359 * g(1) - g(2) * t363;
t364 = qJD(1) ^ 2;
t405 = pkin(9) * t352;
t337 = qJDD(1) * pkin(1) + t364 * t405 + t345;
t346 = -g(1) * t363 - g(2) * t359;
t338 = -pkin(1) * t364 + qJDD(1) * t405 + t346;
t393 = t354 * t362;
t382 = t337 * t393 - t358 * t338;
t388 = qJD(1) * t358;
t403 = pkin(10) * t353;
t280 = -t340 * t403 + t347 * pkin(2) + t348 * t330 + (-g(3) * t362 - t334 * t388) * t352 + t382;
t385 = t352 * t388;
t333 = pkin(2) * t348 - t385 * t403;
t341 = (qJDD(1) * t362 - t358 * t386) * t352;
t379 = t341 * t353 + t347 * t351;
t394 = t354 * t358;
t390 = t337 * t394 + t362 * t338;
t281 = -t348 * t333 + (-g(3) * t358 + t334 * t387) * t352 + t379 * pkin(10) + t390;
t402 = t354 * g(3);
t287 = -t340 * t404 - t341 * pkin(2) - t402 + (-t337 + (-t330 * t362 + t333 * t358) * qJD(1)) * t352;
t357 = sin(qJ(3));
t361 = cos(qJ(3));
t245 = -t357 * t281 + (t280 * t353 + t287 * t351) * t361;
t395 = t353 * t357;
t399 = t351 * t357;
t246 = t280 * t395 + t361 * t281 + t287 * t399;
t398 = t351 * t361;
t319 = t348 * t398 - t357 * t385 + t361 * t381;
t320 = t348 * t399 + (t358 * t361 + t362 * t395) * t389;
t306 = -pkin(3) * t319 - pkin(11) * t320;
t321 = -t341 * t351 + t347 * t353 + qJDD(3);
t331 = t348 * t353 - t351 * t384 + qJD(3);
t329 = t331 ^ 2;
t237 = -pkin(3) * t329 + pkin(11) * t321 + t306 * t319 + t246;
t252 = -t351 * t280 + t353 * t287;
t303 = -t320 * qJD(3) - t357 * t340 + t361 * t379;
t304 = t319 * qJD(3) + t361 * t340 + t357 * t379;
t239 = (-t319 * t331 - t304) * pkin(11) + (t320 * t331 - t303) * pkin(3) + t252;
t356 = sin(qJ(4));
t406 = cos(qJ(4));
t232 = -t356 * t237 + t239 * t406;
t233 = t406 * t237 + t356 * t239;
t311 = t320 * t406 + t356 * t331;
t262 = t311 * qJD(4) + t356 * t304 - t321 * t406;
t310 = t320 * t356 - t331 * t406;
t263 = -t310 * qJD(4) + t304 * t406 + t356 * t321;
t317 = qJD(4) - t319;
t271 = Ifges(5,4) * t311 - Ifges(5,2) * t310 + Ifges(5,6) * t317;
t284 = -mrSges(6,2) * t310 - mrSges(6,3) * t311;
t291 = mrSges(6,1) * t310 - mrSges(6,3) * t317;
t302 = qJDD(4) - t303;
t282 = pkin(4) * t310 - qJ(5) * t311;
t316 = t317 ^ 2;
t230 = -t302 * pkin(4) - t316 * qJ(5) + t311 * t282 + qJDD(5) - t232;
t400 = t310 * t317;
t224 = (t310 * t311 - t302) * pkin(12) + (t263 + t400) * pkin(5) + t230;
t295 = pkin(5) * t311 - pkin(12) * t317;
t309 = t310 ^ 2;
t236 = -t321 * pkin(3) - t329 * pkin(11) + t320 * t306 - t245;
t407 = -2 * qJD(5);
t366 = (-t263 + t400) * qJ(5) + t236 + (t317 * pkin(4) + t407) * t311;
t227 = -t309 * pkin(5) - t311 * t295 + (pkin(4) + pkin(12)) * t262 + t366;
t355 = sin(qJ(6));
t360 = cos(qJ(6));
t221 = t224 * t360 - t227 * t355;
t289 = t310 * t360 - t317 * t355;
t244 = qJD(6) * t289 + t262 * t355 + t302 * t360;
t290 = t310 * t355 + t317 * t360;
t254 = -mrSges(7,1) * t289 + mrSges(7,2) * t290;
t261 = qJDD(6) + t263;
t308 = qJD(6) + t311;
t266 = -mrSges(7,2) * t308 + mrSges(7,3) * t289;
t218 = m(7) * t221 + mrSges(7,1) * t261 - mrSges(7,3) * t244 - t254 * t290 + t266 * t308;
t222 = t224 * t355 + t227 * t360;
t243 = -qJD(6) * t290 + t262 * t360 - t302 * t355;
t267 = mrSges(7,1) * t308 - mrSges(7,3) * t290;
t219 = m(7) * t222 - mrSges(7,2) * t261 + mrSges(7,3) * t243 + t254 * t289 - t267 * t308;
t207 = t360 * t218 + t355 * t219;
t374 = -t316 * pkin(4) + t302 * qJ(5) - t310 * t282 + t233;
t226 = -t262 * pkin(5) - t309 * pkin(12) + ((2 * qJD(5)) + t295) * t317 + t374;
t248 = Ifges(7,5) * t290 + Ifges(7,6) * t289 + Ifges(7,3) * t308;
t250 = Ifges(7,1) * t290 + Ifges(7,4) * t289 + Ifges(7,5) * t308;
t210 = -mrSges(7,1) * t226 + mrSges(7,3) * t222 + Ifges(7,4) * t244 + Ifges(7,2) * t243 + Ifges(7,6) * t261 - t248 * t290 + t250 * t308;
t249 = Ifges(7,4) * t290 + Ifges(7,2) * t289 + Ifges(7,6) * t308;
t211 = mrSges(7,2) * t226 - mrSges(7,3) * t221 + Ifges(7,1) * t244 + Ifges(7,4) * t243 + Ifges(7,5) * t261 + t248 * t289 - t249 * t308;
t228 = t317 * t407 - t374;
t268 = Ifges(6,5) * t317 - Ifges(6,6) * t311 + Ifges(6,3) * t310;
t370 = -mrSges(6,2) * t230 + mrSges(6,3) * t228 - Ifges(6,1) * t302 + Ifges(6,4) * t263 - Ifges(6,5) * t262 + pkin(12) * t207 + t355 * t210 - t360 * t211 + t311 * t268;
t223 = -m(7) * t226 + t243 * mrSges(7,1) - t244 * mrSges(7,2) + t289 * t266 - t290 * t267;
t292 = mrSges(6,1) * t311 + mrSges(6,2) * t317;
t371 = -m(6) * t228 + t302 * mrSges(6,3) + t317 * t292 - t223;
t375 = -m(6) * t230 - t263 * mrSges(6,1) - t311 * t284 - t207;
t270 = Ifges(6,4) * t317 - Ifges(6,2) * t311 + Ifges(6,6) * t310;
t391 = Ifges(5,1) * t311 - Ifges(5,4) * t310 + Ifges(5,5) * t317 - t270;
t408 = t310 * t391 + mrSges(5,1) * t232 - mrSges(5,2) * t233 + Ifges(5,5) * t263 - Ifges(5,6) * t262 + Ifges(5,3) * t302 + pkin(4) * (-t302 * mrSges(6,2) - t317 * t291 + t375) + qJ(5) * (-t262 * mrSges(6,1) - t310 * t284 + t371) + t311 * t271 - t370;
t401 = Ifges(5,4) + Ifges(6,6);
t397 = t352 * t358;
t396 = t352 * t362;
t283 = mrSges(5,1) * t310 + mrSges(5,2) * t311;
t293 = -mrSges(5,2) * t317 - mrSges(5,3) * t310;
t204 = m(5) * t232 - t263 * mrSges(5,3) - t311 * t283 + (-t291 + t293) * t317 + (mrSges(5,1) - mrSges(6,2)) * t302 + t375;
t294 = mrSges(5,1) * t317 - mrSges(5,3) * t311;
t214 = m(5) * t233 - t302 * mrSges(5,2) - t317 * t294 + (-t283 - t284) * t310 + (-mrSges(5,3) - mrSges(6,1)) * t262 + t371;
t199 = t406 * t204 + t356 * t214;
t208 = -t355 * t218 + t360 * t219;
t272 = Ifges(6,1) * t317 - Ifges(6,4) * t311 + Ifges(6,5) * t310;
t392 = -Ifges(5,5) * t311 + Ifges(5,6) * t310 - Ifges(5,3) * t317 - t272;
t305 = -mrSges(4,1) * t319 + mrSges(4,2) * t320;
t313 = mrSges(4,1) * t331 - mrSges(4,3) * t320;
t383 = -t204 * t356 + t406 * t214;
t196 = m(4) * t246 - mrSges(4,2) * t321 + mrSges(4,3) * t303 + t305 * t319 - t313 * t331 + t383;
t312 = -mrSges(4,2) * t331 + mrSges(4,3) * t319;
t198 = m(4) * t252 - mrSges(4,1) * t303 + mrSges(4,2) * t304 - t312 * t319 + t313 * t320 + t199;
t231 = t262 * pkin(4) + t366;
t377 = -m(6) * t231 + t262 * mrSges(6,2) + t310 * t291 - t208;
t367 = -m(5) * t236 - t262 * mrSges(5,1) - t310 * t293 + (t292 - t294) * t311 + (-mrSges(5,2) + mrSges(6,3)) * t263 + t377;
t202 = m(4) * t245 + t321 * mrSges(4,1) - t304 * mrSges(4,3) - t320 * t305 + t331 * t312 + t367;
t185 = t196 * t399 + t353 * t198 + t202 * t398;
t186 = t353 * t361 * t202 + t196 * t395 - t198 * t351;
t314 = -g(3) * t396 + t382;
t336 = -mrSges(3,2) * t348 + mrSges(3,3) * t384;
t339 = (-mrSges(3,1) * t362 + mrSges(3,2) * t358) * t389;
t183 = m(3) * t314 + mrSges(3,1) * t347 - mrSges(3,3) * t340 + t336 * t348 - t339 * t385 + t186;
t192 = t361 * t196 - t202 * t357;
t315 = -g(3) * t397 + t390;
t335 = mrSges(3,1) * t348 - mrSges(3,3) * t385;
t191 = m(3) * t315 - mrSges(3,2) * t347 + mrSges(3,3) * t341 - t335 * t348 + t339 * t384 + t192;
t180 = -t183 * t358 + t362 * t191;
t325 = -t352 * t337 - t402;
t184 = m(3) * t325 - t341 * mrSges(3,1) + t340 * mrSges(3,2) + (t335 * t358 - t336 * t362) * t389 + t185;
t175 = t183 * t393 - t184 * t352 + t191 * t394;
t206 = -t263 * mrSges(6,3) - t311 * t292 - t377;
t369 = -mrSges(6,1) * t228 + mrSges(6,2) * t231 - pkin(5) * t223 - pkin(12) * t208 - t360 * t210 - t355 * t211;
t187 = -mrSges(5,1) * t236 + mrSges(5,3) * t233 - pkin(4) * t206 + t391 * t317 + t392 * t311 + (Ifges(5,6) - Ifges(6,5)) * t302 + t401 * t263 + (-Ifges(5,2) - Ifges(6,3)) * t262 + t369;
t372 = mrSges(7,1) * t221 - mrSges(7,2) * t222 + Ifges(7,5) * t244 + Ifges(7,6) * t243 + Ifges(7,3) * t261 + t290 * t249 - t289 * t250;
t368 = mrSges(6,1) * t230 - mrSges(6,3) * t231 + pkin(5) * t207 + t372;
t188 = t368 + (-t271 + t268) * t317 + t392 * t310 + (Ifges(5,5) - Ifges(6,4)) * t302 + (Ifges(5,1) + Ifges(6,2)) * t263 - t401 * t262 + mrSges(5,2) * t236 - mrSges(5,3) * t232 - qJ(5) * t206;
t297 = Ifges(4,5) * t320 + Ifges(4,6) * t319 + Ifges(4,3) * t331;
t298 = Ifges(4,4) * t320 + Ifges(4,2) * t319 + Ifges(4,6) * t331;
t177 = mrSges(4,2) * t252 - mrSges(4,3) * t245 + Ifges(4,1) * t304 + Ifges(4,4) * t303 + Ifges(4,5) * t321 - pkin(11) * t199 - t356 * t187 + t188 * t406 + t319 * t297 - t331 * t298;
t299 = Ifges(4,1) * t320 + Ifges(4,4) * t319 + Ifges(4,5) * t331;
t181 = -mrSges(4,1) * t252 + mrSges(4,3) * t246 + Ifges(4,4) * t304 + Ifges(4,2) * t303 + Ifges(4,6) * t321 - pkin(3) * t199 - t320 * t297 + t331 * t299 - t408;
t376 = pkin(10) * t192 + t177 * t357 + t181 * t361;
t176 = mrSges(4,1) * t245 - mrSges(4,2) * t246 + Ifges(4,5) * t304 + Ifges(4,6) * t303 + Ifges(4,3) * t321 + pkin(3) * t367 + pkin(11) * t383 + t187 * t406 + t356 * t188 + t320 * t298 - t319 * t299;
t323 = Ifges(3,6) * t348 + (Ifges(3,4) * t358 + Ifges(3,2) * t362) * t389;
t324 = Ifges(3,5) * t348 + (Ifges(3,1) * t358 + Ifges(3,4) * t362) * t389;
t167 = mrSges(3,1) * t314 - mrSges(3,2) * t315 + Ifges(3,5) * t340 + Ifges(3,6) * t341 + Ifges(3,3) * t347 + pkin(2) * t186 + t353 * t176 + (t323 * t358 - t324 * t362) * t389 + t376 * t351;
t322 = Ifges(3,3) * t348 + (Ifges(3,5) * t358 + Ifges(3,6) * t362) * t389;
t169 = -mrSges(3,1) * t325 + mrSges(3,3) * t315 + Ifges(3,4) * t340 + Ifges(3,2) * t341 + Ifges(3,6) * t347 - pkin(2) * t185 - t351 * t176 - t322 * t385 + t348 * t324 + t353 * t376;
t171 = t322 * t384 + mrSges(3,2) * t325 - mrSges(3,3) * t314 + Ifges(3,1) * t340 + Ifges(3,4) * t341 + Ifges(3,5) * t347 + t361 * t177 - t357 * t181 - t348 * t323 + (-t185 * t351 - t186 * t353) * pkin(10);
t373 = mrSges(2,1) * t345 - mrSges(2,2) * t346 + Ifges(2,3) * qJDD(1) + pkin(1) * t175 + t354 * t167 + t169 * t396 + t171 * t397 + t180 * t405;
t178 = m(2) * t346 - mrSges(2,1) * t364 - qJDD(1) * mrSges(2,2) + t180;
t174 = t354 * t184 + (t183 * t362 + t191 * t358) * t352;
t172 = m(2) * t345 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t364 + t175;
t165 = -mrSges(2,2) * g(3) - mrSges(2,3) * t345 + Ifges(2,5) * qJDD(1) - t364 * Ifges(2,6) - t358 * t169 + t362 * t171 + (-t174 * t352 - t175 * t354) * pkin(9);
t164 = mrSges(2,1) * g(3) + mrSges(2,3) * t346 + t364 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t174 - t352 * t167 + (pkin(9) * t180 + t169 * t362 + t171 * t358) * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t363 * t165 - t359 * t164 - pkin(8) * (t172 * t363 + t178 * t359), t165, t171, t177, t188, -t310 * t270 - t370, t211; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t359 * t165 + t363 * t164 + pkin(8) * (-t172 * t359 + t178 * t363), t164, t169, t181, t187, Ifges(6,4) * t302 - Ifges(6,2) * t263 + Ifges(6,6) * t262 - t317 * t268 + t310 * t272 - t368, t210; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t373, t373, t167, t176, t408, Ifges(6,5) * t302 - Ifges(6,6) * t263 + Ifges(6,3) * t262 + t317 * t270 + t311 * t272 - t369, t372;];
m_new  = t1;
