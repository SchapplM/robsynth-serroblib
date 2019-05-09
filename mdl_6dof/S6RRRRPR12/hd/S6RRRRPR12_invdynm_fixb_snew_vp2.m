% Calculate vector of cutting torques with Newton-Euler for
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-05-08 00:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRRPR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 00:10:52
% EndTime: 2019-05-08 00:16:03
% DurationCPUTime: 208.64s
% Computational Cost: add. (3764683->416), mult. (9328442->557), div. (0->0), fcn. (7922680->16), ass. (0->182)
t355 = cos(pkin(6));
t347 = t355 * qJD(1) + qJD(2);
t351 = sin(pkin(7));
t354 = cos(pkin(7));
t352 = sin(pkin(6));
t364 = cos(qJ(2));
t388 = qJD(1) * t364;
t384 = t352 * t388;
t331 = (t347 * t351 + t354 * t384) * pkin(10);
t359 = sin(qJ(2));
t390 = qJD(1) * t352;
t402 = pkin(10) * t351;
t335 = (-pkin(2) * t364 - t359 * t402) * t390;
t387 = qJD(1) * qJD(2);
t341 = (qJDD(1) * t359 + t364 * t387) * t352;
t346 = t355 * qJDD(1) + qJDD(2);
t360 = sin(qJ(1));
t365 = cos(qJ(1));
t344 = t360 * g(1) - t365 * g(2);
t366 = qJD(1) ^ 2;
t403 = pkin(9) * t352;
t338 = qJDD(1) * pkin(1) + t366 * t403 + t344;
t345 = -t365 * g(1) - t360 * g(2);
t339 = -t366 * pkin(1) + qJDD(1) * t403 + t345;
t392 = t355 * t364;
t380 = t338 * t392 - t359 * t339;
t389 = qJD(1) * t359;
t401 = pkin(10) * t354;
t291 = -t341 * t401 + t346 * pkin(2) + t347 * t331 + (-g(3) * t364 - t335 * t389) * t352 + t380;
t385 = t352 * t389;
t334 = t347 * pkin(2) - t385 * t401;
t342 = (qJDD(1) * t364 - t359 * t387) * t352;
t377 = t342 * t354 + t346 * t351;
t393 = t355 * t359;
t391 = t338 * t393 + t364 * t339;
t292 = -t347 * t334 + (-g(3) * t359 + t335 * t388) * t352 + t377 * pkin(10) + t391;
t400 = t355 * g(3);
t298 = -t341 * t402 - t342 * pkin(2) - t400 + (-t338 + (-t331 * t364 + t334 * t359) * qJD(1)) * t352;
t358 = sin(qJ(3));
t363 = cos(qJ(3));
t256 = -t358 * t292 + (t291 * t354 + t298 * t351) * t363;
t394 = t354 * t364;
t399 = t351 * t358;
t322 = t347 * t399 + (t358 * t394 + t359 * t363) * t390;
t308 = -t322 * qJD(3) - t358 * t341 + t363 * t377;
t398 = t351 * t363;
t321 = (-t358 * t359 + t363 * t394) * t390 + t347 * t398;
t404 = 2 * qJD(5);
t397 = t352 * t359;
t396 = t352 * t364;
t395 = t354 * t358;
t257 = t291 * t395 + t363 * t292 + t298 * t399;
t311 = -t321 * pkin(3) - t322 * pkin(11);
t323 = -t351 * t342 + t354 * t346 + qJDD(3);
t332 = t354 * t347 - t351 * t384 + qJD(3);
t330 = t332 ^ 2;
t247 = -t330 * pkin(3) + t323 * pkin(11) + t321 * t311 + t257;
t273 = -t351 * t291 + t354 * t298;
t309 = t321 * qJD(3) + t363 * t341 + t358 * t377;
t250 = (-t321 * t332 - t309) * pkin(11) + (t322 * t332 - t308) * pkin(3) + t273;
t357 = sin(qJ(4));
t362 = cos(qJ(4));
t239 = -t357 * t247 + t362 * t250;
t313 = -t357 * t322 + t362 * t332;
t276 = t313 * qJD(4) + t362 * t309 + t357 * t323;
t307 = qJDD(4) - t308;
t314 = t362 * t322 + t357 * t332;
t320 = qJD(4) - t321;
t236 = (t313 * t320 - t276) * qJ(5) + (t313 * t314 + t307) * pkin(4) + t239;
t240 = t362 * t247 + t357 * t250;
t275 = -t314 * qJD(4) - t357 * t309 + t362 * t323;
t301 = t320 * pkin(4) - t314 * qJ(5);
t312 = t313 ^ 2;
t238 = -t312 * pkin(4) + t275 * qJ(5) - t320 * t301 + t240;
t350 = sin(pkin(13));
t353 = cos(pkin(13));
t293 = t353 * t313 - t350 * t314;
t233 = t350 * t236 + t353 * t238 + t293 * t404;
t262 = t353 * t275 - t350 * t276;
t294 = t350 * t313 + t353 * t314;
t271 = -mrSges(6,1) * t293 + mrSges(6,2) * t294;
t280 = t320 * mrSges(6,1) - t294 * mrSges(6,3);
t272 = -pkin(5) * t293 - pkin(12) * t294;
t319 = t320 ^ 2;
t230 = -t319 * pkin(5) + t307 * pkin(12) + t293 * t272 + t233;
t246 = -t323 * pkin(3) - t330 * pkin(11) + t322 * t311 - t256;
t241 = -t275 * pkin(4) - t312 * qJ(5) + t314 * t301 + qJDD(5) + t246;
t263 = t350 * t275 + t353 * t276;
t234 = (-t293 * t320 - t263) * pkin(12) + (t294 * t320 - t262) * pkin(5) + t241;
t356 = sin(qJ(6));
t361 = cos(qJ(6));
t227 = -t356 * t230 + t361 * t234;
t277 = -t356 * t294 + t361 * t320;
t244 = t277 * qJD(6) + t361 * t263 + t356 * t307;
t261 = qJDD(6) - t262;
t278 = t361 * t294 + t356 * t320;
t264 = -mrSges(7,1) * t277 + mrSges(7,2) * t278;
t290 = qJD(6) - t293;
t265 = -mrSges(7,2) * t290 + mrSges(7,3) * t277;
t223 = m(7) * t227 + mrSges(7,1) * t261 - mrSges(7,3) * t244 - t264 * t278 + t265 * t290;
t228 = t361 * t230 + t356 * t234;
t243 = -t278 * qJD(6) - t356 * t263 + t361 * t307;
t266 = mrSges(7,1) * t290 - mrSges(7,3) * t278;
t224 = m(7) * t228 - mrSges(7,2) * t261 + mrSges(7,3) * t243 + t264 * t277 - t266 * t290;
t381 = -t356 * t223 + t361 * t224;
t207 = m(6) * t233 - t307 * mrSges(6,2) + t262 * mrSges(6,3) + t293 * t271 - t320 * t280 + t381;
t379 = -t353 * t236 + t350 * t238;
t232 = -0.2e1 * qJD(5) * t294 - t379;
t279 = -t320 * mrSges(6,2) + t293 * mrSges(6,3);
t229 = -t307 * pkin(5) - t319 * pkin(12) + (t404 + t272) * t294 + t379;
t373 = -m(7) * t229 + t243 * mrSges(7,1) - t244 * mrSges(7,2) + t277 * t265 - t278 * t266;
t219 = m(6) * t232 + t307 * mrSges(6,1) - t263 * mrSges(6,3) - t294 * t271 + t320 * t279 + t373;
t204 = t350 * t207 + t353 * t219;
t295 = -t313 * mrSges(5,1) + t314 * mrSges(5,2);
t300 = -t320 * mrSges(5,2) + t313 * mrSges(5,3);
t202 = m(5) * t239 + t307 * mrSges(5,1) - t276 * mrSges(5,3) - t314 * t295 + t320 * t300 + t204;
t302 = t320 * mrSges(5,1) - t314 * mrSges(5,3);
t382 = t353 * t207 - t350 * t219;
t203 = m(5) * t240 - t307 * mrSges(5,2) + t275 * mrSges(5,3) + t313 * t295 - t320 * t302 + t382;
t196 = t362 * t202 + t357 * t203;
t212 = t361 * t223 + t356 * t224;
t310 = -t321 * mrSges(4,1) + t322 * mrSges(4,2);
t316 = t332 * mrSges(4,1) - t322 * mrSges(4,3);
t383 = -t357 * t202 + t362 * t203;
t193 = m(4) * t257 - t323 * mrSges(4,2) + t308 * mrSges(4,3) + t321 * t310 - t332 * t316 + t383;
t315 = -t332 * mrSges(4,2) + t321 * mrSges(4,3);
t195 = m(4) * t273 - t308 * mrSges(4,1) + t309 * mrSges(4,2) - t321 * t315 + t322 * t316 + t196;
t371 = m(6) * t241 - t262 * mrSges(6,1) + t263 * mrSges(6,2) - t293 * t279 + t294 * t280 + t212;
t368 = -m(5) * t246 + t275 * mrSges(5,1) - t276 * mrSges(5,2) + t313 * t300 - t314 * t302 - t371;
t210 = m(4) * t256 + t323 * mrSges(4,1) - t309 * mrSges(4,3) - t322 * t310 + t332 * t315 + t368;
t182 = t193 * t399 + t354 * t195 + t210 * t398;
t183 = t354 * t363 * t210 + t193 * t395 - t351 * t195;
t317 = -g(3) * t396 + t380;
t337 = -t347 * mrSges(3,2) + mrSges(3,3) * t384;
t340 = (-mrSges(3,1) * t364 + mrSges(3,2) * t359) * t390;
t180 = m(3) * t317 + t346 * mrSges(3,1) - t341 * mrSges(3,3) + t347 * t337 - t340 * t385 + t183;
t189 = t363 * t193 - t358 * t210;
t318 = -g(3) * t397 + t391;
t336 = t347 * mrSges(3,1) - mrSges(3,3) * t385;
t188 = m(3) * t318 - t346 * mrSges(3,2) + t342 * mrSges(3,3) - t347 * t336 + t340 * t384 + t189;
t177 = -t359 * t180 + t364 * t188;
t327 = -t352 * t338 - t400;
t181 = m(3) * t327 - t342 * mrSges(3,1) + t341 * mrSges(3,2) + (t336 * t359 - t337 * t364) * t390 + t182;
t172 = t180 * t392 - t352 * t181 + t188 * t393;
t251 = Ifges(7,5) * t278 + Ifges(7,6) * t277 + Ifges(7,3) * t290;
t253 = Ifges(7,1) * t278 + Ifges(7,4) * t277 + Ifges(7,5) * t290;
t216 = -mrSges(7,1) * t229 + mrSges(7,3) * t228 + Ifges(7,4) * t244 + Ifges(7,2) * t243 + Ifges(7,6) * t261 - t251 * t278 + t253 * t290;
t252 = Ifges(7,4) * t278 + Ifges(7,2) * t277 + Ifges(7,6) * t290;
t217 = mrSges(7,2) * t229 - mrSges(7,3) * t227 + Ifges(7,1) * t244 + Ifges(7,4) * t243 + Ifges(7,5) * t261 + t251 * t277 - t252 * t290;
t267 = Ifges(6,5) * t294 + Ifges(6,6) * t293 + Ifges(6,3) * t320;
t268 = Ifges(6,4) * t294 + Ifges(6,2) * t293 + Ifges(6,6) * t320;
t197 = mrSges(6,2) * t241 - mrSges(6,3) * t232 + Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t307 - pkin(12) * t212 - t356 * t216 + t361 * t217 + t293 * t267 - t320 * t268;
t269 = Ifges(6,1) * t294 + Ifges(6,4) * t293 + Ifges(6,5) * t320;
t369 = mrSges(7,1) * t227 - mrSges(7,2) * t228 + Ifges(7,5) * t244 + Ifges(7,6) * t243 + Ifges(7,3) * t261 + t278 * t252 - t277 * t253;
t198 = -mrSges(6,1) * t241 + mrSges(6,3) * t233 + Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t307 - pkin(5) * t212 - t294 * t267 + t320 * t269 - t369;
t281 = Ifges(5,5) * t314 + Ifges(5,6) * t313 + Ifges(5,3) * t320;
t283 = Ifges(5,1) * t314 + Ifges(5,4) * t313 + Ifges(5,5) * t320;
t184 = -mrSges(5,1) * t246 + mrSges(5,3) * t240 + Ifges(5,4) * t276 + Ifges(5,2) * t275 + Ifges(5,6) * t307 - pkin(4) * t371 + qJ(5) * t382 + t350 * t197 + t353 * t198 - t314 * t281 + t320 * t283;
t282 = Ifges(5,4) * t314 + Ifges(5,2) * t313 + Ifges(5,6) * t320;
t185 = mrSges(5,2) * t246 - mrSges(5,3) * t239 + Ifges(5,1) * t276 + Ifges(5,4) * t275 + Ifges(5,5) * t307 - qJ(5) * t204 + t353 * t197 - t350 * t198 + t313 * t281 - t320 * t282;
t303 = Ifges(4,5) * t322 + Ifges(4,6) * t321 + Ifges(4,3) * t332;
t304 = Ifges(4,4) * t322 + Ifges(4,2) * t321 + Ifges(4,6) * t332;
t174 = mrSges(4,2) * t273 - mrSges(4,3) * t256 + Ifges(4,1) * t309 + Ifges(4,4) * t308 + Ifges(4,5) * t323 - pkin(11) * t196 - t357 * t184 + t362 * t185 + t321 * t303 - t332 * t304;
t305 = Ifges(4,1) * t322 + Ifges(4,4) * t321 + Ifges(4,5) * t332;
t370 = -mrSges(6,1) * t232 + mrSges(6,2) * t233 - Ifges(6,5) * t263 - Ifges(6,6) * t262 - Ifges(6,3) * t307 - pkin(5) * t373 - pkin(12) * t381 - t361 * t216 - t356 * t217 - t294 * t268 + t293 * t269;
t367 = mrSges(5,1) * t239 - mrSges(5,2) * t240 + Ifges(5,5) * t276 + Ifges(5,6) * t275 + Ifges(5,3) * t307 + pkin(4) * t204 + t314 * t282 - t313 * t283 - t370;
t178 = -mrSges(4,1) * t273 + mrSges(4,3) * t257 + Ifges(4,4) * t309 + Ifges(4,2) * t308 + Ifges(4,6) * t323 - pkin(3) * t196 - t322 * t303 + t332 * t305 - t367;
t374 = pkin(10) * t189 + t174 * t358 + t178 * t363;
t173 = mrSges(4,1) * t256 - mrSges(4,2) * t257 + Ifges(4,5) * t309 + Ifges(4,6) * t308 + Ifges(4,3) * t323 + pkin(3) * t368 + pkin(11) * t383 + t362 * t184 + t357 * t185 + t322 * t304 - t321 * t305;
t325 = Ifges(3,6) * t347 + (Ifges(3,4) * t359 + Ifges(3,2) * t364) * t390;
t326 = Ifges(3,5) * t347 + (Ifges(3,1) * t359 + Ifges(3,4) * t364) * t390;
t164 = mrSges(3,1) * t317 - mrSges(3,2) * t318 + Ifges(3,5) * t341 + Ifges(3,6) * t342 + Ifges(3,3) * t346 + pkin(2) * t183 + t354 * t173 + (t325 * t359 - t326 * t364) * t390 + t374 * t351;
t324 = Ifges(3,3) * t347 + (Ifges(3,5) * t359 + Ifges(3,6) * t364) * t390;
t166 = -mrSges(3,1) * t327 + mrSges(3,3) * t318 + Ifges(3,4) * t341 + Ifges(3,2) * t342 + Ifges(3,6) * t346 - pkin(2) * t182 - t351 * t173 - t324 * t385 + t347 * t326 + t354 * t374;
t168 = t324 * t384 + mrSges(3,2) * t327 - mrSges(3,3) * t317 + Ifges(3,1) * t341 + Ifges(3,4) * t342 + Ifges(3,5) * t346 + t363 * t174 - t358 * t178 - t347 * t325 + (-t182 * t351 - t183 * t354) * pkin(10);
t372 = mrSges(2,1) * t344 - mrSges(2,2) * t345 + Ifges(2,3) * qJDD(1) + pkin(1) * t172 + t355 * t164 + t166 * t396 + t168 * t397 + t177 * t403;
t175 = m(2) * t345 - t366 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t177;
t171 = t355 * t181 + (t180 * t364 + t188 * t359) * t352;
t169 = m(2) * t344 + qJDD(1) * mrSges(2,1) - t366 * mrSges(2,2) + t172;
t162 = -mrSges(2,2) * g(3) - mrSges(2,3) * t344 + Ifges(2,5) * qJDD(1) - t366 * Ifges(2,6) - t359 * t166 + t364 * t168 + (-t171 * t352 - t172 * t355) * pkin(9);
t161 = mrSges(2,1) * g(3) + mrSges(2,3) * t345 + t366 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t171 - t352 * t164 + (pkin(9) * t177 + t166 * t364 + t168 * t359) * t355;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t365 * t162 - t360 * t161 - pkin(8) * (t365 * t169 + t360 * t175), t162, t168, t174, t185, t197, t217; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t360 * t162 + t365 * t161 + pkin(8) * (-t360 * t169 + t365 * t175), t161, t166, t178, t184, t198, t216; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t372, t372, t164, t173, t367, -t370, t369;];
m_new  = t1;
