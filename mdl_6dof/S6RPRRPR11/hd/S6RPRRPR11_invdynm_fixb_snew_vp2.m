% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR11
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
% Datum: 2019-05-06 00:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:05:58
% EndTime: 2019-05-06 00:09:03
% DurationCPUTime: 180.11s
% Computational Cost: add. (2879283->402), mult. (9011985->552), div. (0->0), fcn. (7741630->16), ass. (0->183)
t348 = sin(pkin(12));
t350 = sin(pkin(6));
t352 = cos(pkin(12));
t354 = cos(pkin(6));
t357 = sin(qJ(3));
t353 = cos(pkin(7));
t360 = cos(qJ(3));
t393 = t353 * t360;
t349 = sin(pkin(7));
t398 = t349 * t360;
t366 = t350 * (-t348 * t357 + t352 * t393) + t354 * t398;
t318 = t366 * qJD(1);
t394 = t353 * t357;
t399 = t349 * t357;
t368 = t354 * t399 + (t348 * t360 + t352 * t394) * t350;
t319 = t368 * qJD(1);
t307 = -t319 * qJD(3) + qJDD(1) * t366;
t396 = t350 * t353;
t332 = (t349 * t354 + t352 * t396) * qJD(1) * pkin(9);
t358 = sin(qJ(1));
t361 = cos(qJ(1));
t344 = -g(1) * t361 - g(2) * t358;
t362 = qJD(1) ^ 2;
t402 = qJ(2) * t350;
t336 = -pkin(1) * t362 + qJDD(1) * t402 + t344;
t405 = pkin(9) * t348;
t380 = -pkin(2) * t352 - t349 * t405;
t392 = qJD(1) * t350;
t403 = pkin(9) * qJDD(1);
t375 = qJD(1) * t380 * t392 + t353 * t403;
t343 = g(1) * t358 - g(2) * t361;
t335 = qJDD(1) * pkin(1) + t362 * t402 + t343;
t388 = qJD(2) * t392;
t395 = t352 * t354;
t397 = t350 * t352;
t381 = -g(3) * t397 + t335 * t395 - 0.2e1 * t348 * t388;
t286 = (pkin(2) * qJDD(1) + qJD(1) * t332) * t354 + (-t350 * t375 - t336) * t348 + t381;
t337 = (pkin(2) * t354 - t396 * t405) * qJD(1);
t400 = t348 * t354;
t389 = t335 * t400 + (t336 + 0.2e1 * t388) * t352;
t287 = (-qJD(1) * t337 + t349 * t403) * t354 + (-g(3) * t348 + t352 * t375) * t350 + t389;
t387 = -t354 * g(3) + qJDD(2);
t298 = (-t335 + t380 * qJDD(1) + (-t332 * t352 + t337 * t348) * qJD(1)) * t350 + t387;
t255 = -t357 * t287 + (t286 * t353 + t298 * t349) * t360;
t406 = cos(qJ(4));
t404 = Ifges(3,3) * t354;
t401 = t348 * t350;
t256 = t286 * t394 + t287 * t360 + t298 * t399;
t306 = -pkin(3) * t318 - pkin(10) * t319;
t376 = -t349 * t397 + t353 * t354;
t333 = qJD(1) * t376 + qJD(3);
t329 = t333 ^ 2;
t330 = qJDD(1) * t376 + qJDD(3);
t247 = -pkin(3) * t329 + pkin(10) * t330 + t306 * t318 + t256;
t264 = -t349 * t286 + t298 * t353;
t308 = t318 * qJD(3) + qJDD(1) * t368;
t253 = (-t318 * t333 - t308) * pkin(10) + (t319 * t333 - t307) * pkin(3) + t264;
t356 = sin(qJ(4));
t237 = t247 * t406 + t253 * t356;
t312 = t319 * t356 - t333 * t406;
t313 = t319 * t406 + t333 * t356;
t288 = pkin(4) * t312 - qJ(5) * t313;
t304 = qJDD(4) - t307;
t317 = qJD(4) - t318;
t316 = t317 ^ 2;
t232 = -pkin(4) * t316 + qJ(5) * t304 - t288 * t312 + t237;
t246 = -t330 * pkin(3) - t329 * pkin(10) + t306 * t319 - t255;
t281 = qJD(4) * t313 + t308 * t356 - t330 * t406;
t282 = -qJD(4) * t312 + t308 * t406 + t330 * t356;
t235 = (t312 * t317 - t282) * qJ(5) + (t313 * t317 + t281) * pkin(4) + t246;
t347 = sin(pkin(13));
t351 = cos(pkin(13));
t297 = t313 * t351 + t317 * t347;
t227 = -0.2e1 * qJD(5) * t297 - t347 * t232 + t235 * t351;
t266 = t282 * t351 + t304 * t347;
t296 = -t313 * t347 + t317 * t351;
t225 = (t296 * t312 - t266) * pkin(11) + (t296 * t297 + t281) * pkin(5) + t227;
t228 = 0.2e1 * qJD(5) * t296 + t232 * t351 + t235 * t347;
t265 = -t282 * t347 + t304 * t351;
t273 = pkin(5) * t312 - pkin(11) * t297;
t295 = t296 ^ 2;
t226 = -pkin(5) * t295 + pkin(11) * t265 - t273 * t312 + t228;
t355 = sin(qJ(6));
t359 = cos(qJ(6));
t223 = t225 * t359 - t226 * t355;
t267 = t296 * t359 - t297 * t355;
t243 = qJD(6) * t267 + t265 * t355 + t266 * t359;
t268 = t296 * t355 + t297 * t359;
t254 = -mrSges(7,1) * t267 + mrSges(7,2) * t268;
t309 = qJD(6) + t312;
t257 = -mrSges(7,2) * t309 + mrSges(7,3) * t267;
t279 = qJDD(6) + t281;
t219 = m(7) * t223 + mrSges(7,1) * t279 - mrSges(7,3) * t243 - t254 * t268 + t257 * t309;
t224 = t225 * t355 + t226 * t359;
t242 = -qJD(6) * t268 + t265 * t359 - t266 * t355;
t258 = mrSges(7,1) * t309 - mrSges(7,3) * t268;
t220 = m(7) * t224 - mrSges(7,2) * t279 + mrSges(7,3) * t242 + t254 * t267 - t258 * t309;
t211 = t219 * t359 + t220 * t355;
t269 = -mrSges(6,1) * t296 + mrSges(6,2) * t297;
t271 = -mrSges(6,2) * t312 + mrSges(6,3) * t296;
t209 = m(6) * t227 + mrSges(6,1) * t281 - mrSges(6,3) * t266 - t269 * t297 + t271 * t312 + t211;
t272 = mrSges(6,1) * t312 - mrSges(6,3) * t297;
t385 = -t219 * t355 + t220 * t359;
t210 = m(6) * t228 - mrSges(6,2) * t281 + mrSges(6,3) * t265 + t269 * t296 - t272 * t312 + t385;
t207 = -t209 * t347 + t210 * t351;
t289 = mrSges(5,1) * t312 + mrSges(5,2) * t313;
t300 = mrSges(5,1) * t317 - mrSges(5,3) * t313;
t205 = m(5) * t237 - mrSges(5,2) * t304 - mrSges(5,3) * t281 - t289 * t312 - t300 * t317 + t207;
t236 = -t247 * t356 + t253 * t406;
t231 = -pkin(4) * t304 - qJ(5) * t316 + t288 * t313 + qJDD(5) - t236;
t229 = -pkin(5) * t265 - pkin(11) * t295 + t273 * t297 + t231;
t373 = m(7) * t229 - mrSges(7,1) * t242 + mrSges(7,2) * t243 - t257 * t267 + t258 * t268;
t221 = -m(6) * t231 + mrSges(6,1) * t265 - mrSges(6,2) * t266 + t271 * t296 - t272 * t297 - t373;
t299 = -mrSges(5,2) * t317 - mrSges(5,3) * t312;
t215 = m(5) * t236 + mrSges(5,1) * t304 - mrSges(5,3) * t282 - t289 * t313 + t299 * t317 + t221;
t197 = t205 * t356 + t215 * t406;
t305 = -mrSges(4,1) * t318 + mrSges(4,2) * t319;
t315 = mrSges(4,1) * t333 - mrSges(4,3) * t319;
t386 = t205 * t406 - t215 * t356;
t194 = m(4) * t256 - mrSges(4,2) * t330 + mrSges(4,3) * t307 + t305 * t318 - t315 * t333 + t386;
t314 = -mrSges(4,2) * t333 + mrSges(4,3) * t318;
t196 = m(4) * t264 - mrSges(4,1) * t307 + mrSges(4,2) * t308 - t314 * t318 + t315 * t319 + t197;
t206 = t209 * t351 + t210 * t347;
t365 = -m(5) * t246 - mrSges(5,1) * t281 - mrSges(5,2) * t282 - t299 * t312 - t300 * t313 - t206;
t202 = m(4) * t255 + mrSges(4,1) * t330 - mrSges(4,3) * t308 - t305 * t319 + t314 * t333 + t365;
t183 = t194 * t399 + t196 * t353 + t202 * t398;
t184 = t194 * t394 - t349 * t196 + t202 * t393;
t310 = -t348 * t336 + t381;
t384 = -mrSges(3,1) * t352 + mrSges(3,2) * t348;
t334 = t384 * t392;
t378 = -mrSges(3,2) * t354 + mrSges(3,3) * t397;
t339 = t378 * qJD(1);
t379 = mrSges(3,1) * t354 - mrSges(3,3) * t401;
t181 = m(3) * t310 + t379 * qJDD(1) + (-t334 * t401 + t339 * t354) * qJD(1) + t184;
t189 = t194 * t360 - t357 * t202;
t311 = -g(3) * t401 + t389;
t338 = t379 * qJD(1);
t188 = m(3) * t311 + t378 * qJDD(1) + (t334 * t397 - t338 * t354) * qJD(1) + t189;
t178 = -t181 * t348 + t188 * t352;
t320 = -t350 * t335 + t387;
t182 = m(3) * t320 + (t384 * qJDD(1) + (t338 * t348 - t339 * t352) * qJD(1)) * t350 + t183;
t173 = t181 * t395 - t182 * t350 + t188 * t400;
t383 = Ifges(3,5) * t348 + Ifges(3,6) * t352;
t249 = Ifges(7,5) * t268 + Ifges(7,6) * t267 + Ifges(7,3) * t309;
t251 = Ifges(7,1) * t268 + Ifges(7,4) * t267 + Ifges(7,5) * t309;
t212 = -mrSges(7,1) * t229 + mrSges(7,3) * t224 + Ifges(7,4) * t243 + Ifges(7,2) * t242 + Ifges(7,6) * t279 - t249 * t268 + t251 * t309;
t250 = Ifges(7,4) * t268 + Ifges(7,2) * t267 + Ifges(7,6) * t309;
t213 = mrSges(7,2) * t229 - mrSges(7,3) * t223 + Ifges(7,1) * t243 + Ifges(7,4) * t242 + Ifges(7,5) * t279 + t249 * t267 - t250 * t309;
t259 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t312;
t261 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t312;
t198 = -mrSges(6,1) * t231 + mrSges(6,3) * t228 + Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t281 - pkin(5) * t373 + pkin(11) * t385 + t359 * t212 + t355 * t213 - t297 * t259 + t312 * t261;
t260 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t312;
t199 = mrSges(6,2) * t231 - mrSges(6,3) * t227 + Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t281 - pkin(11) * t211 - t212 * t355 + t213 * t359 + t259 * t296 - t260 * t312;
t275 = Ifges(5,5) * t313 - Ifges(5,6) * t312 + Ifges(5,3) * t317;
t276 = Ifges(5,4) * t313 - Ifges(5,2) * t312 + Ifges(5,6) * t317;
t185 = mrSges(5,2) * t246 - mrSges(5,3) * t236 + Ifges(5,1) * t282 - Ifges(5,4) * t281 + Ifges(5,5) * t304 - qJ(5) * t206 - t198 * t347 + t199 * t351 - t275 * t312 - t276 * t317;
t277 = Ifges(5,1) * t313 - Ifges(5,4) * t312 + Ifges(5,5) * t317;
t369 = -mrSges(7,1) * t223 + mrSges(7,2) * t224 - Ifges(7,5) * t243 - Ifges(7,6) * t242 - Ifges(7,3) * t279 - t250 * t268 + t267 * t251;
t364 = -mrSges(6,1) * t227 + mrSges(6,2) * t228 - Ifges(6,5) * t266 - Ifges(6,6) * t265 - pkin(5) * t211 - t297 * t260 + t296 * t261 + t369;
t190 = -pkin(4) * t206 + t364 + (-Ifges(5,2) - Ifges(6,3)) * t281 + t317 * t277 - t313 * t275 + Ifges(5,6) * t304 + Ifges(5,4) * t282 - mrSges(5,1) * t246 + mrSges(5,3) * t237;
t301 = Ifges(4,5) * t319 + Ifges(4,6) * t318 + Ifges(4,3) * t333;
t302 = Ifges(4,4) * t319 + Ifges(4,2) * t318 + Ifges(4,6) * t333;
t175 = mrSges(4,2) * t264 - mrSges(4,3) * t255 + Ifges(4,1) * t308 + Ifges(4,4) * t307 + Ifges(4,5) * t330 - pkin(10) * t197 + t185 * t406 - t190 * t356 + t301 * t318 - t302 * t333;
t303 = Ifges(4,1) * t319 + Ifges(4,4) * t318 + Ifges(4,5) * t333;
t363 = mrSges(5,1) * t236 - mrSges(5,2) * t237 + Ifges(5,5) * t282 - Ifges(5,6) * t281 + Ifges(5,3) * t304 + pkin(4) * t221 + qJ(5) * t207 + t198 * t351 + t199 * t347 + t276 * t313 + t277 * t312;
t179 = -mrSges(4,1) * t264 + mrSges(4,3) * t256 + Ifges(4,4) * t308 + Ifges(4,2) * t307 + Ifges(4,6) * t330 - pkin(3) * t197 - t301 * t319 + t303 * t333 - t363;
t374 = pkin(9) * t189 + t175 * t357 + t179 * t360;
t372 = Ifges(3,5) * t354 + (Ifges(3,1) * t348 + Ifges(3,4) * t352) * t350;
t371 = Ifges(3,6) * t354 + (Ifges(3,4) * t348 + Ifges(3,2) * t352) * t350;
t174 = mrSges(4,1) * t255 - mrSges(4,2) * t256 + Ifges(4,5) * t308 + Ifges(4,6) * t307 + Ifges(4,3) * t330 + pkin(3) * t365 + pkin(10) * t386 + t356 * t185 + t190 * t406 + t319 * t302 - t318 * t303;
t324 = t371 * qJD(1);
t325 = t372 * qJD(1);
t165 = qJDD(1) * t404 + mrSges(3,1) * t310 - mrSges(3,2) * t311 + pkin(2) * t184 + t353 * t174 + t374 * t349 + (t383 * qJDD(1) + (t324 * t348 - t325 * t352) * qJD(1)) * t350;
t323 = (t350 * t383 + t404) * qJD(1);
t167 = -mrSges(3,1) * t320 + mrSges(3,3) * t311 - pkin(2) * t183 - t349 * t174 + (-t323 * t401 + t325 * t354) * qJD(1) + t374 * t353 + t371 * qJDD(1);
t169 = mrSges(3,2) * t320 - mrSges(3,3) * t310 + t360 * t175 - t357 * t179 + (t323 * t397 - t324 * t354) * qJD(1) + (-t183 * t349 - t184 * t353) * pkin(9) + t372 * qJDD(1);
t370 = mrSges(2,1) * t343 - mrSges(2,2) * t344 + Ifges(2,3) * qJDD(1) + pkin(1) * t173 + t165 * t354 + t167 * t397 + t169 * t401 + t178 * t402;
t176 = m(2) * t344 - mrSges(2,1) * t362 - qJDD(1) * mrSges(2,2) + t178;
t172 = t354 * t182 + (t181 * t352 + t188 * t348) * t350;
t170 = m(2) * t343 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t362 + t173;
t163 = -mrSges(2,2) * g(3) - mrSges(2,3) * t343 + Ifges(2,5) * qJDD(1) - t362 * Ifges(2,6) - t348 * t167 + t352 * t169 + (-t172 * t350 - t173 * t354) * qJ(2);
t162 = mrSges(2,1) * g(3) + mrSges(2,3) * t344 + t362 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t172 - t350 * t165 + (qJ(2) * t178 + t167 * t352 + t169 * t348) * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t361 * t163 - t358 * t162 - pkin(8) * (t170 * t361 + t176 * t358), t163, t169, t175, t185, t199, t213; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t358 * t163 + t361 * t162 + pkin(8) * (-t170 * t358 + t176 * t361), t162, t167, t179, t190, t198, t212; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t370, t370, t165, t174, t363, Ifges(6,3) * t281 - t364, -t369;];
m_new  = t1;
