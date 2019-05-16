% Calculate vector of cutting torques with Newton-Euler for
% S6RRPRPR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-05-06 16:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPRPR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR13_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 16:32:14
% EndTime: 2019-05-06 16:33:05
% DurationCPUTime: 27.96s
% Computational Cost: add. (491306->406), mult. (1113990->507), div. (0->0), fcn. (821165->12), ass. (0->162)
t403 = -2 * qJD(3);
t354 = cos(pkin(6));
t345 = qJD(1) * t354 + qJD(2);
t357 = sin(qJ(2));
t352 = sin(pkin(6));
t388 = qJD(1) * t352;
t382 = t357 * t388;
t402 = (pkin(2) * t345 + t403) * t382;
t358 = sin(qJ(1));
t362 = cos(qJ(1));
t340 = t358 * g(1) - g(2) * t362;
t363 = qJD(1) ^ 2;
t400 = pkin(8) * t352;
t323 = qJDD(1) * pkin(1) + t363 * t400 + t340;
t341 = -g(1) * t362 - g(2) * t358;
t385 = qJDD(1) * t352;
t324 = -pkin(1) * t363 + pkin(8) * t385 + t341;
t361 = cos(qJ(2));
t393 = t354 * t357;
t395 = t352 * t357;
t282 = -g(3) * t395 + t323 * t393 + t361 * t324;
t325 = (-pkin(2) * t361 - qJ(3) * t357) * t388;
t343 = t345 ^ 2;
t344 = qJDD(1) * t354 + qJDD(2);
t387 = qJD(1) * t361;
t381 = t352 * t387;
t251 = t343 * pkin(2) - t344 * qJ(3) - t325 * t381 + t345 * t403 - t282;
t401 = -pkin(2) - pkin(9);
t399 = t354 * g(3);
t398 = mrSges(3,1) - mrSges(4,2);
t397 = Ifges(3,4) + Ifges(4,6);
t396 = t352 ^ 2 * t363;
t394 = t352 * t361;
t392 = t354 * t361;
t328 = pkin(3) * t382 - pkin(9) * t345;
t329 = (qJD(2) * t387 + qJDD(1) * t357) * t352;
t330 = -qJD(2) * t382 + t361 * t385;
t384 = t361 ^ 2 * t396;
t245 = -pkin(3) * t384 - t399 - t329 * qJ(3) + t401 * t330 + (-t323 + (-qJ(3) * t345 * t361 - t328 * t357) * qJD(1)) * t352 + t402;
t389 = g(3) * t394 + t357 * t324;
t377 = -t343 * qJ(3) + t325 * t382 + qJDD(3) + t389;
t247 = t329 * pkin(3) + t401 * t344 + (-pkin(3) * t345 * t388 - pkin(9) * t357 * t396 - t323 * t354) * t361 + t377;
t356 = sin(qJ(4));
t360 = cos(qJ(4));
t229 = t360 * t245 + t356 * t247;
t308 = t345 * t356 + t360 * t381;
t309 = t345 * t360 - t356 * t381;
t283 = pkin(4) * t308 - qJ(5) * t309;
t318 = qJDD(4) + t329;
t335 = qJD(4) + t382;
t333 = t335 ^ 2;
t223 = -pkin(4) * t333 + qJ(5) * t318 - t283 * t308 + t229;
t244 = t330 * pkin(3) - pkin(9) * t384 + t345 * t328 - t251;
t279 = qJD(4) * t309 + t360 * t330 + t344 * t356;
t280 = -qJD(4) * t308 - t330 * t356 + t344 * t360;
t226 = (t308 * t335 - t280) * qJ(5) + (t309 * t335 + t279) * pkin(4) + t244;
t351 = sin(pkin(11));
t353 = cos(pkin(11));
t289 = t309 * t353 + t335 * t351;
t218 = -0.2e1 * qJD(5) * t289 - t351 * t223 + t353 * t226;
t263 = t280 * t353 + t318 * t351;
t288 = -t309 * t351 + t335 * t353;
t216 = (t288 * t308 - t263) * pkin(10) + (t288 * t289 + t279) * pkin(5) + t218;
t219 = 0.2e1 * qJD(5) * t288 + t353 * t223 + t351 * t226;
t262 = -t280 * t351 + t318 * t353;
t269 = pkin(5) * t308 - pkin(10) * t289;
t287 = t288 ^ 2;
t217 = -pkin(5) * t287 + pkin(10) * t262 - t269 * t308 + t219;
t355 = sin(qJ(6));
t359 = cos(qJ(6));
t214 = t216 * t359 - t217 * t355;
t259 = t288 * t359 - t289 * t355;
t235 = qJD(6) * t259 + t262 * t355 + t263 * t359;
t260 = t288 * t355 + t289 * t359;
t240 = -mrSges(7,1) * t259 + mrSges(7,2) * t260;
t307 = qJD(6) + t308;
t248 = -mrSges(7,2) * t307 + mrSges(7,3) * t259;
t277 = qJDD(6) + t279;
t209 = m(7) * t214 + mrSges(7,1) * t277 - mrSges(7,3) * t235 - t240 * t260 + t248 * t307;
t215 = t216 * t355 + t217 * t359;
t234 = -qJD(6) * t260 + t262 * t359 - t263 * t355;
t249 = mrSges(7,1) * t307 - mrSges(7,3) * t260;
t210 = m(7) * t215 - mrSges(7,2) * t277 + mrSges(7,3) * t234 + t240 * t259 - t249 * t307;
t202 = t359 * t209 + t355 * t210;
t264 = -mrSges(6,1) * t288 + mrSges(6,2) * t289;
t267 = -mrSges(6,2) * t308 + mrSges(6,3) * t288;
t200 = m(6) * t218 + mrSges(6,1) * t279 - mrSges(6,3) * t263 - t264 * t289 + t267 * t308 + t202;
t268 = mrSges(6,1) * t308 - mrSges(6,3) * t289;
t379 = -t209 * t355 + t359 * t210;
t201 = m(6) * t219 - mrSges(6,2) * t279 + mrSges(6,3) * t262 + t264 * t288 - t268 * t308 + t379;
t195 = t353 * t200 + t351 * t201;
t298 = Ifges(4,1) * t345 + (-Ifges(4,4) * t357 - Ifges(4,5) * t361) * t388;
t391 = Ifges(3,3) * t345 + (Ifges(3,5) * t357 + Ifges(3,6) * t361) * t388 + t298;
t296 = Ifges(4,5) * t345 + (-Ifges(4,6) * t357 - Ifges(4,3) * t361) * t388;
t390 = t296 - Ifges(3,6) * t345 - (Ifges(3,4) * t357 + Ifges(3,2) * t361) * t388;
t383 = t323 * t392;
t281 = t383 - t389;
t320 = -mrSges(3,2) * t345 + mrSges(3,3) * t381;
t321 = -mrSges(4,1) * t381 - mrSges(4,3) * t345;
t326 = (mrSges(4,2) * t361 - mrSges(4,3) * t357) * t388;
t327 = (-mrSges(3,1) * t361 + mrSges(3,2) * t357) * t388;
t284 = mrSges(5,1) * t308 + mrSges(5,2) * t309;
t291 = mrSges(5,1) * t335 - mrSges(5,3) * t309;
t380 = -t200 * t351 + t353 * t201;
t192 = m(5) * t229 - mrSges(5,2) * t318 - mrSges(5,3) * t279 - t284 * t308 - t291 * t335 + t380;
t228 = -t356 * t245 + t247 * t360;
t290 = -mrSges(5,2) * t335 - mrSges(5,3) * t308;
t222 = -pkin(4) * t318 - qJ(5) * t333 + t309 * t283 + qJDD(5) - t228;
t220 = -pkin(5) * t262 - pkin(10) * t287 + t269 * t289 + t222;
t375 = m(7) * t220 - t234 * mrSges(7,1) + mrSges(7,2) * t235 - t259 * t248 + t249 * t260;
t366 = -m(6) * t222 + t262 * mrSges(6,1) - mrSges(6,2) * t263 + t288 * t267 - t268 * t289 - t375;
t205 = m(5) * t228 + mrSges(5,1) * t318 - mrSges(5,3) * t280 - t284 * t309 + t290 * t335 + t366;
t182 = t356 * t192 + t360 * t205;
t258 = -t344 * pkin(2) + t377 - t383;
t376 = -m(4) * t258 - t329 * mrSges(4,1) - t182;
t180 = m(3) * t281 - t329 * mrSges(3,3) + (t320 - t321) * t345 + t398 * t344 + (-t326 - t327) * t382 + t376;
t319 = mrSges(3,1) * t345 - mrSges(3,3) * t382;
t322 = mrSges(4,1) * t382 + mrSges(4,2) * t345;
t372 = m(5) * t244 + t279 * mrSges(5,1) + t280 * mrSges(5,2) + t308 * t290 + t309 * t291 + t195;
t367 = -m(4) * t251 + t344 * mrSges(4,3) + t345 * t322 + t326 * t381 + t372;
t190 = -t345 * t319 - t344 * mrSges(3,2) + m(3) * t282 + (mrSges(3,3) + mrSges(4,1)) * t330 + t367 + t327 * t381;
t176 = -t180 * t357 + t361 * t190;
t183 = t360 * t192 - t356 * t205;
t299 = -t352 * t323 - t399;
t252 = -t330 * pkin(2) + (-t345 * t381 - t329) * qJ(3) + t299 + t402;
t378 = m(4) * t252 - t329 * mrSges(4,3) + t321 * t381 + t183;
t179 = m(3) * t299 + t329 * mrSges(3,2) - t398 * t330 + (-t320 * t361 + (t319 - t322) * t357) * t388 + t378;
t171 = -t179 * t352 + t180 * t392 + t190 * t393;
t295 = Ifges(3,5) * t345 + (Ifges(3,1) * t357 + Ifges(3,4) * t361) * t388;
t236 = Ifges(7,5) * t260 + Ifges(7,6) * t259 + Ifges(7,3) * t307;
t238 = Ifges(7,1) * t260 + Ifges(7,4) * t259 + Ifges(7,5) * t307;
t203 = -mrSges(7,1) * t220 + mrSges(7,3) * t215 + Ifges(7,4) * t235 + Ifges(7,2) * t234 + Ifges(7,6) * t277 - t236 * t260 + t238 * t307;
t237 = Ifges(7,4) * t260 + Ifges(7,2) * t259 + Ifges(7,6) * t307;
t204 = mrSges(7,2) * t220 - mrSges(7,3) * t214 + Ifges(7,1) * t235 + Ifges(7,4) * t234 + Ifges(7,5) * t277 + t236 * t259 - t237 * t307;
t253 = Ifges(6,5) * t289 + Ifges(6,6) * t288 + Ifges(6,3) * t308;
t255 = Ifges(6,1) * t289 + Ifges(6,4) * t288 + Ifges(6,5) * t308;
t185 = -mrSges(6,1) * t222 + mrSges(6,3) * t219 + Ifges(6,4) * t263 + Ifges(6,2) * t262 + Ifges(6,6) * t279 - pkin(5) * t375 + pkin(10) * t379 + t359 * t203 + t355 * t204 - t289 * t253 + t308 * t255;
t254 = Ifges(6,4) * t289 + Ifges(6,2) * t288 + Ifges(6,6) * t308;
t187 = mrSges(6,2) * t222 - mrSges(6,3) * t218 + Ifges(6,1) * t263 + Ifges(6,4) * t262 + Ifges(6,5) * t279 - pkin(10) * t202 - t203 * t355 + t204 * t359 + t253 * t288 - t254 * t308;
t271 = Ifges(5,5) * t309 - Ifges(5,6) * t308 + Ifges(5,3) * t335;
t272 = Ifges(5,4) * t309 - Ifges(5,2) * t308 + Ifges(5,6) * t335;
t173 = mrSges(5,2) * t244 - mrSges(5,3) * t228 + Ifges(5,1) * t280 - Ifges(5,4) * t279 + Ifges(5,5) * t318 - qJ(5) * t195 - t185 * t351 + t187 * t353 - t271 * t308 - t272 * t335;
t273 = Ifges(5,1) * t309 - Ifges(5,4) * t308 + Ifges(5,5) * t335;
t371 = -mrSges(7,1) * t214 + mrSges(7,2) * t215 - Ifges(7,5) * t235 - Ifges(7,6) * t234 - Ifges(7,3) * t277 - t260 * t237 + t259 * t238;
t364 = -mrSges(6,1) * t218 + mrSges(6,2) * t219 - Ifges(6,5) * t263 - Ifges(6,6) * t262 - pkin(5) * t202 - t289 * t254 + t288 * t255 + t371;
t177 = t335 * t273 + Ifges(5,6) * t318 - t309 * t271 + Ifges(5,4) * t280 - mrSges(5,1) * t244 + mrSges(5,3) * t229 + t364 - pkin(4) * t195 + (-Ifges(5,2) - Ifges(6,3)) * t279;
t297 = Ifges(4,4) * t345 + (-Ifges(4,2) * t357 - Ifges(4,6) * t361) * t388;
t370 = mrSges(4,2) * t258 - mrSges(4,3) * t251 + Ifges(4,1) * t344 - Ifges(4,4) * t329 - Ifges(4,5) * t330 - pkin(9) * t182 + t360 * t173 - t356 * t177 + t297 * t381;
t163 = qJ(3) * (mrSges(4,1) * t330 + t367) + Ifges(3,3) * t344 + Ifges(3,6) * t330 + Ifges(3,5) * t329 + mrSges(3,1) * t281 - mrSges(3,2) * t282 + t370 + pkin(2) * (-t344 * mrSges(4,2) - t345 * t321 + t376) + (-t295 * t361 + (-pkin(2) * t326 - t390) * t357) * t388;
t181 = t330 * mrSges(4,2) - t322 * t382 + t378;
t368 = -mrSges(4,1) * t251 + mrSges(4,2) * t252 + pkin(3) * t372 - pkin(9) * t183 - t356 * t173 - t360 * t177;
t165 = -mrSges(3,1) * t299 + mrSges(3,3) * t282 - pkin(2) * t181 + (t295 - t297) * t345 + (Ifges(3,6) - Ifges(4,5)) * t344 + (Ifges(3,2) + Ifges(4,3)) * t330 + t397 * t329 - t391 * t382 + t368;
t369 = -mrSges(5,1) * t228 + mrSges(5,2) * t229 - Ifges(5,5) * t280 + Ifges(5,6) * t279 - Ifges(5,3) * t318 - pkin(4) * t366 - qJ(5) * t380 - t353 * t185 - t351 * t187 - t309 * t272 - t308 * t273;
t365 = -mrSges(4,1) * t258 + mrSges(4,3) * t252 - pkin(3) * t182 + t369;
t167 = mrSges(3,2) * t299 - mrSges(3,3) * t281 + t390 * t345 + (Ifges(3,5) - Ifges(4,4)) * t344 + t397 * t330 + (Ifges(3,1) + Ifges(4,2)) * t329 - qJ(3) * t181 + t391 * t381 - t365;
t374 = mrSges(2,1) * t340 - mrSges(2,2) * t341 + Ifges(2,3) * qJDD(1) + pkin(1) * t171 + t354 * t163 + t165 * t394 + t167 * t395 + t176 * t400;
t174 = m(2) * t341 - mrSges(2,1) * t363 - qJDD(1) * mrSges(2,2) + t176;
t170 = t354 * t179 + (t180 * t361 + t190 * t357) * t352;
t168 = m(2) * t340 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t363 + t171;
t161 = -mrSges(2,2) * g(3) - mrSges(2,3) * t340 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - t357 * t165 + t361 * t167 + (-t170 * t352 - t171 * t354) * pkin(8);
t160 = mrSges(2,1) * g(3) + mrSges(2,3) * t341 + t363 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t170 - t352 * t163 + (pkin(8) * t176 + t165 * t361 + t167 * t357) * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t362 * t161 - t358 * t160 - pkin(7) * (t168 * t362 + t174 * t358), t161, t167, -t296 * t382 + t370, t173, t187, t204; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t358 * t161 + t362 * t160 + pkin(7) * (-t168 * t358 + t174 * t362), t160, t165, Ifges(4,4) * t344 - Ifges(4,2) * t329 - Ifges(4,6) * t330 - t345 * t296 - t298 * t381 + t365, t177, t185, t203; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t374, t374, t163, Ifges(4,5) * t344 - Ifges(4,6) * t329 - Ifges(4,3) * t330 + t345 * t297 + t298 * t382 - t368, -t369, Ifges(6,3) * t279 - t364, -t371;];
m_new  = t1;
