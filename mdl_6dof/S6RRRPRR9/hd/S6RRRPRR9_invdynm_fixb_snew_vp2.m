% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-05-07 13:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:54:08
% EndTime: 2019-05-07 12:59:06
% DurationCPUTime: 199.70s
% Computational Cost: add. (3453321->414), mult. (9059003->556), div. (0->0), fcn. (7693825->16), ass. (0->181)
t354 = cos(pkin(6));
t346 = t354 * qJD(1) + qJD(2);
t350 = sin(pkin(7));
t353 = cos(pkin(7));
t351 = sin(pkin(6));
t363 = cos(qJ(2));
t382 = qJD(1) * t363;
t378 = t351 * t382;
t330 = (t346 * t350 + t353 * t378) * pkin(10);
t358 = sin(qJ(2));
t384 = qJD(1) * t351;
t397 = pkin(10) * t350;
t334 = (-pkin(2) * t363 - t358 * t397) * t384;
t380 = qJD(1) * qJD(2);
t340 = (qJDD(1) * t358 + t363 * t380) * t351;
t345 = t354 * qJDD(1) + qJDD(2);
t359 = sin(qJ(1));
t364 = cos(qJ(1));
t343 = t359 * g(1) - t364 * g(2);
t365 = qJD(1) ^ 2;
t398 = pkin(9) * t351;
t337 = qJDD(1) * pkin(1) + t365 * t398 + t343;
t344 = -t364 * g(1) - t359 * g(2);
t338 = -t365 * pkin(1) + qJDD(1) * t398 + t344;
t386 = t354 * t363;
t374 = t337 * t386 - t358 * t338;
t383 = qJD(1) * t358;
t396 = pkin(10) * t353;
t289 = -t340 * t396 + t345 * pkin(2) + t346 * t330 + (-g(3) * t363 - t334 * t383) * t351 + t374;
t379 = t351 * t383;
t333 = t346 * pkin(2) - t379 * t396;
t341 = (qJDD(1) * t363 - t358 * t380) * t351;
t373 = t341 * t353 + t345 * t350;
t387 = t354 * t358;
t385 = t337 * t387 + t363 * t338;
t290 = -t346 * t333 + (-g(3) * t358 + t334 * t382) * t351 + t373 * pkin(10) + t385;
t395 = t354 * g(3);
t294 = -t340 * t397 - t341 * pkin(2) - t395 + (-t337 + (-t330 * t363 + t333 * t358) * qJD(1)) * t351;
t357 = sin(qJ(3));
t362 = cos(qJ(3));
t389 = t353 * t362;
t393 = t350 * t362;
t252 = t289 * t389 - t357 * t290 + t294 * t393;
t388 = t353 * t363;
t319 = t346 * t393 + (-t357 * t358 + t362 * t388) * t384;
t305 = t319 * qJD(3) + t362 * t340 + t373 * t357;
t394 = t350 * t357;
t320 = t346 * t394 + (t357 * t388 + t358 * t362) * t384;
t322 = -t350 * t341 + t353 * t345 + qJDD(3);
t331 = t353 * t346 - t350 * t378 + qJD(3);
t242 = (t319 * t331 - t305) * qJ(4) + (t319 * t320 + t322) * pkin(3) + t252;
t390 = t353 * t357;
t253 = t289 * t390 + t362 * t290 + t294 * t394;
t304 = -t320 * qJD(3) - t357 * t340 + t373 * t362;
t314 = t331 * pkin(3) - t320 * qJ(4);
t318 = t319 ^ 2;
t245 = -t318 * pkin(3) + t304 * qJ(4) - t331 * t314 + t253;
t349 = sin(pkin(13));
t352 = cos(pkin(13));
t311 = t349 * t319 + t352 * t320;
t234 = -0.2e1 * qJD(4) * t311 + t352 * t242 - t349 * t245;
t392 = t351 * t358;
t391 = t351 * t363;
t310 = t352 * t319 - t349 * t320;
t235 = 0.2e1 * qJD(4) * t310 + t349 * t242 + t352 * t245;
t276 = t352 * t304 - t349 * t305;
t284 = -t310 * mrSges(5,1) + t311 * mrSges(5,2);
t299 = t331 * mrSges(5,1) - t311 * mrSges(5,3);
t285 = -t310 * pkin(4) - t311 * pkin(11);
t329 = t331 ^ 2;
t232 = -t329 * pkin(4) + t322 * pkin(11) + t310 * t285 + t235;
t265 = -t350 * t289 + t353 * t294;
t251 = -t304 * pkin(3) - t318 * qJ(4) + t320 * t314 + qJDD(4) + t265;
t277 = t349 * t304 + t352 * t305;
t237 = (-t310 * t331 - t277) * pkin(11) + (t311 * t331 - t276) * pkin(4) + t251;
t356 = sin(qJ(5));
t361 = cos(qJ(5));
t228 = t361 * t232 + t356 * t237;
t296 = -t356 * t311 + t361 * t331;
t297 = t361 * t311 + t356 * t331;
t269 = -pkin(5) * t296 - pkin(12) * t297;
t273 = qJDD(5) - t276;
t309 = qJD(5) - t310;
t308 = t309 ^ 2;
t226 = -pkin(5) * t308 + pkin(12) * t273 + t269 * t296 + t228;
t231 = -t322 * pkin(4) - t329 * pkin(11) + t311 * t285 - t234;
t257 = -t297 * qJD(5) - t356 * t277 + t361 * t322;
t258 = t296 * qJD(5) + t361 * t277 + t356 * t322;
t229 = (-t296 * t309 - t258) * pkin(12) + (t297 * t309 - t257) * pkin(5) + t231;
t355 = sin(qJ(6));
t360 = cos(qJ(6));
t222 = -t355 * t226 + t360 * t229;
t274 = -t355 * t297 + t360 * t309;
t240 = t274 * qJD(6) + t360 * t258 + t355 * t273;
t275 = t360 * t297 + t355 * t309;
t254 = -mrSges(7,1) * t274 + mrSges(7,2) * t275;
t256 = qJDD(6) - t257;
t295 = qJD(6) - t296;
t259 = -mrSges(7,2) * t295 + mrSges(7,3) * t274;
t220 = m(7) * t222 + mrSges(7,1) * t256 - mrSges(7,3) * t240 - t254 * t275 + t259 * t295;
t223 = t360 * t226 + t355 * t229;
t239 = -t275 * qJD(6) - t355 * t258 + t360 * t273;
t260 = mrSges(7,1) * t295 - mrSges(7,3) * t275;
t221 = m(7) * t223 - mrSges(7,2) * t256 + mrSges(7,3) * t239 + t254 * t274 - t260 * t295;
t214 = -t355 * t220 + t360 * t221;
t268 = -mrSges(6,1) * t296 + mrSges(6,2) * t297;
t279 = mrSges(6,1) * t309 - mrSges(6,3) * t297;
t211 = m(6) * t228 - t273 * mrSges(6,2) + t257 * mrSges(6,3) + t296 * t268 - t309 * t279 + t214;
t227 = -t356 * t232 + t361 * t237;
t225 = -t273 * pkin(5) - t308 * pkin(12) + t297 * t269 - t227;
t224 = -m(7) * t225 + t239 * mrSges(7,1) - mrSges(7,2) * t240 + t274 * t259 - t260 * t275;
t278 = -mrSges(6,2) * t309 + mrSges(6,3) * t296;
t218 = m(6) * t227 + mrSges(6,1) * t273 - mrSges(6,3) * t258 - t268 * t297 + t278 * t309 + t224;
t376 = t361 * t211 - t356 * t218;
t202 = m(5) * t235 - t322 * mrSges(5,2) + t276 * mrSges(5,3) + t310 * t284 - t331 * t299 + t376;
t298 = -t331 * mrSges(5,2) + t310 * mrSges(5,3);
t213 = t360 * t220 + t355 * t221;
t368 = -m(6) * t231 + t257 * mrSges(6,1) - t258 * mrSges(6,2) + t296 * t278 - t297 * t279 - t213;
t208 = m(5) * t234 + t322 * mrSges(5,1) - t277 * mrSges(5,3) - t311 * t284 + t331 * t298 + t368;
t195 = t349 * t202 + t352 * t208;
t206 = t356 * t211 + t361 * t218;
t312 = -t319 * mrSges(4,1) + t320 * mrSges(4,2);
t313 = -t331 * mrSges(4,2) + t319 * mrSges(4,3);
t193 = m(4) * t252 + t322 * mrSges(4,1) - t305 * mrSges(4,3) - t320 * t312 + t331 * t313 + t195;
t315 = t331 * mrSges(4,1) - t320 * mrSges(4,3);
t377 = t352 * t202 - t349 * t208;
t194 = m(4) * t253 - t322 * mrSges(4,2) + t304 * mrSges(4,3) + t319 * t312 - t331 * t315 + t377;
t370 = m(5) * t251 - t276 * mrSges(5,1) + t277 * mrSges(5,2) - t310 * t298 + t311 * t299 + t206;
t204 = m(4) * t265 - t304 * mrSges(4,1) + t305 * mrSges(4,2) - t319 * t313 + t320 * t315 + t370;
t180 = t193 * t393 + t194 * t394 + t353 * t204;
t181 = t193 * t389 + t194 * t390 - t350 * t204;
t316 = -g(3) * t391 + t374;
t336 = -t346 * mrSges(3,2) + mrSges(3,3) * t378;
t339 = (-mrSges(3,1) * t363 + mrSges(3,2) * t358) * t384;
t178 = m(3) * t316 + t345 * mrSges(3,1) - t340 * mrSges(3,3) + t346 * t336 - t339 * t379 + t181;
t186 = -t357 * t193 + t362 * t194;
t317 = -g(3) * t392 + t385;
t335 = t346 * mrSges(3,1) - mrSges(3,3) * t379;
t185 = m(3) * t317 - t345 * mrSges(3,2) + t341 * mrSges(3,3) - t346 * t335 + t339 * t378 + t186;
t175 = -t358 * t178 + t363 * t185;
t326 = -t351 * t337 - t395;
t179 = m(3) * t326 - t341 * mrSges(3,1) + t340 * mrSges(3,2) + (t335 * t358 - t336 * t363) * t384 + t180;
t170 = t178 * t386 - t351 * t179 + t185 * t387;
t246 = Ifges(7,5) * t275 + Ifges(7,6) * t274 + Ifges(7,3) * t295;
t248 = Ifges(7,1) * t275 + Ifges(7,4) * t274 + Ifges(7,5) * t295;
t215 = -mrSges(7,1) * t225 + mrSges(7,3) * t223 + Ifges(7,4) * t240 + Ifges(7,2) * t239 + Ifges(7,6) * t256 - t246 * t275 + t248 * t295;
t247 = Ifges(7,4) * t275 + Ifges(7,2) * t274 + Ifges(7,6) * t295;
t216 = mrSges(7,2) * t225 - mrSges(7,3) * t222 + Ifges(7,1) * t240 + Ifges(7,4) * t239 + Ifges(7,5) * t256 + t246 * t274 - t247 * t295;
t261 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t309;
t262 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t309;
t197 = mrSges(6,2) * t231 - mrSges(6,3) * t227 + Ifges(6,1) * t258 + Ifges(6,4) * t257 + Ifges(6,5) * t273 - pkin(12) * t213 - t355 * t215 + t360 * t216 + t296 * t261 - t309 * t262;
t263 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t309;
t367 = mrSges(7,1) * t222 - mrSges(7,2) * t223 + Ifges(7,5) * t240 + Ifges(7,6) * t239 + Ifges(7,3) * t256 + t247 * t275 - t274 * t248;
t199 = -mrSges(6,1) * t231 + mrSges(6,3) * t228 + Ifges(6,4) * t258 + Ifges(6,2) * t257 + Ifges(6,6) * t273 - pkin(5) * t213 - t261 * t297 + t263 * t309 - t367;
t280 = Ifges(5,5) * t311 + Ifges(5,6) * t310 + Ifges(5,3) * t331;
t281 = Ifges(5,4) * t311 + Ifges(5,2) * t310 + Ifges(5,6) * t331;
t182 = mrSges(5,2) * t251 - mrSges(5,3) * t234 + Ifges(5,1) * t277 + Ifges(5,4) * t276 + Ifges(5,5) * t322 - pkin(11) * t206 + t361 * t197 - t356 * t199 + t310 * t280 - t331 * t281;
t282 = Ifges(5,1) * t311 + Ifges(5,4) * t310 + Ifges(5,5) * t331;
t366 = mrSges(6,1) * t227 - mrSges(6,2) * t228 + Ifges(6,5) * t258 + Ifges(6,6) * t257 + Ifges(6,3) * t273 + pkin(5) * t224 + pkin(12) * t214 + t360 * t215 + t355 * t216 + t297 * t262 - t296 * t263;
t187 = -mrSges(5,1) * t251 + mrSges(5,3) * t235 + Ifges(5,4) * t277 + Ifges(5,2) * t276 + Ifges(5,6) * t322 - pkin(4) * t206 - t311 * t280 + t331 * t282 - t366;
t300 = Ifges(4,5) * t320 + Ifges(4,6) * t319 + Ifges(4,3) * t331;
t302 = Ifges(4,1) * t320 + Ifges(4,4) * t319 + Ifges(4,5) * t331;
t171 = -mrSges(4,1) * t265 + mrSges(4,3) * t253 + Ifges(4,4) * t305 + Ifges(4,2) * t304 + Ifges(4,6) * t322 - pkin(3) * t370 + qJ(4) * t377 + t349 * t182 + t352 * t187 - t320 * t300 + t331 * t302;
t301 = Ifges(4,4) * t320 + Ifges(4,2) * t319 + Ifges(4,6) * t331;
t172 = mrSges(4,2) * t265 - mrSges(4,3) * t252 + Ifges(4,1) * t305 + Ifges(4,4) * t304 + Ifges(4,5) * t322 - qJ(4) * t195 + t352 * t182 - t349 * t187 + t319 * t300 - t331 * t301;
t372 = pkin(10) * t186 + t171 * t362 + t172 * t357;
t369 = mrSges(5,1) * t234 - mrSges(5,2) * t235 + Ifges(5,5) * t277 + Ifges(5,6) * t276 + Ifges(5,3) * t322 + pkin(4) * t368 + pkin(11) * t376 + t356 * t197 + t361 * t199 + t311 * t281 - t310 * t282;
t176 = mrSges(4,1) * t252 - mrSges(4,2) * t253 + Ifges(4,5) * t305 + Ifges(4,6) * t304 + Ifges(4,3) * t322 + pkin(3) * t195 + t320 * t301 - t319 * t302 + t369;
t324 = Ifges(3,6) * t346 + (Ifges(3,4) * t358 + Ifges(3,2) * t363) * t384;
t325 = Ifges(3,5) * t346 + (Ifges(3,1) * t358 + Ifges(3,4) * t363) * t384;
t162 = mrSges(3,1) * t316 - mrSges(3,2) * t317 + Ifges(3,5) * t340 + Ifges(3,6) * t341 + Ifges(3,3) * t345 + pkin(2) * t181 + t353 * t176 + (t324 * t358 - t325 * t363) * t384 + t372 * t350;
t323 = Ifges(3,3) * t346 + (Ifges(3,5) * t358 + Ifges(3,6) * t363) * t384;
t164 = -mrSges(3,1) * t326 + mrSges(3,3) * t317 + Ifges(3,4) * t340 + Ifges(3,2) * t341 + Ifges(3,6) * t345 - pkin(2) * t180 - t350 * t176 - t323 * t379 + t346 * t325 + t372 * t353;
t166 = t323 * t378 + mrSges(3,2) * t326 - mrSges(3,3) * t316 + Ifges(3,1) * t340 + Ifges(3,4) * t341 + Ifges(3,5) * t345 - t357 * t171 + t362 * t172 - t346 * t324 + (-t180 * t350 - t181 * t353) * pkin(10);
t371 = mrSges(2,1) * t343 - mrSges(2,2) * t344 + Ifges(2,3) * qJDD(1) + pkin(1) * t170 + t354 * t162 + t164 * t391 + t166 * t392 + t175 * t398;
t173 = m(2) * t344 - t365 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t175;
t169 = t354 * t179 + (t178 * t363 + t185 * t358) * t351;
t167 = m(2) * t343 + qJDD(1) * mrSges(2,1) - t365 * mrSges(2,2) + t170;
t160 = -mrSges(2,2) * g(3) - mrSges(2,3) * t343 + Ifges(2,5) * qJDD(1) - t365 * Ifges(2,6) - t358 * t164 + t363 * t166 + (-t169 * t351 - t170 * t354) * pkin(9);
t159 = mrSges(2,1) * g(3) + mrSges(2,3) * t344 + t365 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t169 - t351 * t162 + (pkin(9) * t175 + t164 * t363 + t166 * t358) * t354;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t364 * t160 - t359 * t159 - pkin(8) * (t364 * t167 + t359 * t173), t160, t166, t172, t182, t197, t216; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t359 * t160 + t364 * t159 + pkin(8) * (-t359 * t167 + t364 * t173), t159, t164, t171, t187, t199, t215; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t371, t371, t162, t176, t369, t366, t367;];
m_new  = t1;
