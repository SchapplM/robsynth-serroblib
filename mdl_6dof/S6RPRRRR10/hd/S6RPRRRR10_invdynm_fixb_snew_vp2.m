% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR10
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
% Datum: 2019-05-06 04:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR10_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR10_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 04:41:50
% EndTime: 2019-05-06 04:44:55
% DurationCPUTime: 177.80s
% Computational Cost: add. (2885175->403), mult. (8946745->549), div. (0->0), fcn. (7736194->16), ass. (0->188)
t356 = sin(pkin(13));
t358 = sin(pkin(6));
t359 = cos(pkin(13));
t361 = cos(pkin(6));
t370 = cos(qJ(3));
t360 = cos(pkin(7));
t365 = sin(qJ(3));
t405 = t360 * t365;
t357 = sin(pkin(7));
t409 = t357 * t365;
t377 = t361 * t409 + (t356 * t370 + t359 * t405) * t358;
t328 = t377 * qJD(1);
t407 = t358 * t360;
t388 = t357 * t361 + t359 * t407;
t384 = t388 * t370;
t411 = t356 * t358;
t402 = t365 * t411;
t314 = -t328 * qJD(3) + (t384 - t402) * qJDD(1);
t383 = t388 * qJD(1);
t339 = pkin(9) * t383;
t366 = sin(qJ(1));
t371 = cos(qJ(1));
t353 = -g(1) * t371 - g(2) * t366;
t372 = qJD(1) ^ 2;
t413 = qJ(2) * t358;
t343 = -pkin(1) * t372 + qJDD(1) * t413 + t353;
t416 = pkin(9) * t356;
t391 = -pkin(2) * t359 - t357 * t416;
t403 = qJD(1) * t358;
t414 = pkin(9) * qJDD(1);
t386 = qJD(1) * t391 * t403 + t360 * t414;
t352 = t366 * g(1) - g(2) * t371;
t342 = qJDD(1) * pkin(1) + t372 * t413 + t352;
t400 = qJD(2) * t403;
t406 = t359 * t361;
t408 = t358 * t359;
t392 = -g(3) * t408 + t342 * t406 - 0.2e1 * t356 * t400;
t293 = (pkin(2) * qJDD(1) + qJD(1) * t339) * t361 + (-t358 * t386 - t343) * t356 + t392;
t344 = (pkin(2) * t361 - t407 * t416) * qJD(1);
t410 = t356 * t361;
t401 = t342 * t410 + (t343 + 0.2e1 * t400) * t359;
t294 = (-qJD(1) * t344 + t357 * t414) * t361 + (-g(3) * t356 + t359 * t386) * t358 + t401;
t399 = -t361 * g(3) + qJDD(2);
t302 = (-t342 + t391 * qJDD(1) + (-t339 * t359 + t344 * t356) * qJD(1)) * t358 + t399;
t266 = -t365 * t294 + (t293 * t360 + t302 * t357) * t370;
t415 = Ifges(3,3) * t361;
t349 = qJD(1) * t402;
t327 = t370 * t383 - t349;
t312 = -mrSges(4,1) * t327 + mrSges(4,2) * t328;
t315 = t327 * qJD(3) + qJDD(1) * t377;
t387 = -t357 * t408 + t360 * t361;
t340 = qJD(1) * t387 + qJD(3);
t321 = -mrSges(4,2) * t340 + mrSges(4,3) * t327;
t337 = qJDD(1) * t387 + qJDD(3);
t313 = -pkin(3) * t327 - pkin(10) * t328;
t336 = t340 ^ 2;
t250 = -t337 * pkin(3) - t336 * pkin(10) + t328 * t313 - t266;
t364 = sin(qJ(4));
t369 = cos(qJ(4));
t320 = t328 * t369 + t340 * t364;
t288 = -qJD(4) * t320 - t315 * t364 + t337 * t369;
t319 = -t328 * t364 + t340 * t369;
t289 = qJD(4) * t319 + t315 * t369 + t337 * t364;
t326 = -qJD(1) * t384 + qJD(4) + t349;
t303 = -mrSges(5,2) * t326 + mrSges(5,3) * t319;
t304 = mrSges(5,1) * t326 - mrSges(5,3) * t320;
t267 = t293 * t405 + t370 * t294 + t302 * t409;
t251 = -pkin(3) * t336 + pkin(10) * t337 + t313 * t327 + t267;
t277 = -t357 * t293 + t360 * t302;
t254 = (-t327 * t340 - t315) * pkin(10) + (t328 * t340 - t314) * pkin(3) + t277;
t243 = -t364 * t251 + t369 * t254;
t311 = qJDD(4) - t314;
t240 = (t319 * t326 - t289) * pkin(11) + (t319 * t320 + t311) * pkin(4) + t243;
t244 = t369 * t251 + t364 * t254;
t305 = pkin(4) * t326 - pkin(11) * t320;
t318 = t319 ^ 2;
t242 = -pkin(4) * t318 + pkin(11) * t288 - t305 * t326 + t244;
t363 = sin(qJ(5));
t368 = cos(qJ(5));
t237 = t363 * t240 + t368 * t242;
t296 = t319 * t368 - t320 * t363;
t297 = t319 * t363 + t320 * t368;
t276 = -pkin(5) * t296 - pkin(12) * t297;
t310 = qJDD(5) + t311;
t324 = qJD(5) + t326;
t323 = t324 ^ 2;
t234 = -pkin(5) * t323 + pkin(12) * t310 + t276 * t296 + t237;
t245 = -t288 * pkin(4) - t318 * pkin(11) + t320 * t305 + t250;
t263 = -qJD(5) * t297 + t288 * t368 - t289 * t363;
t264 = qJD(5) * t296 + t288 * t363 + t289 * t368;
t238 = (-t296 * t324 - t264) * pkin(12) + (t297 * t324 - t263) * pkin(5) + t245;
t362 = sin(qJ(6));
t367 = cos(qJ(6));
t231 = -t234 * t362 + t238 * t367;
t278 = -t297 * t362 + t324 * t367;
t248 = qJD(6) * t278 + t264 * t367 + t310 * t362;
t262 = qJDD(6) - t263;
t279 = t297 * t367 + t324 * t362;
t268 = -mrSges(7,1) * t278 + mrSges(7,2) * t279;
t295 = qJD(6) - t296;
t269 = -mrSges(7,2) * t295 + mrSges(7,3) * t278;
t227 = m(7) * t231 + mrSges(7,1) * t262 - mrSges(7,3) * t248 - t268 * t279 + t269 * t295;
t232 = t234 * t367 + t238 * t362;
t247 = -qJD(6) * t279 - t264 * t362 + t310 * t367;
t270 = mrSges(7,1) * t295 - mrSges(7,3) * t279;
t228 = m(7) * t232 - mrSges(7,2) * t262 + mrSges(7,3) * t247 + t268 * t278 - t270 * t295;
t216 = t367 * t227 + t362 * t228;
t280 = -mrSges(6,2) * t324 + mrSges(6,3) * t296;
t281 = mrSges(6,1) * t324 - mrSges(6,3) * t297;
t378 = m(6) * t245 - t263 * mrSges(6,1) + mrSges(6,2) * t264 - t296 * t280 + t281 * t297 + t216;
t374 = -m(5) * t250 + t288 * mrSges(5,1) - mrSges(5,2) * t289 + t319 * t303 - t304 * t320 - t378;
t211 = m(4) * t266 + mrSges(4,1) * t337 - mrSges(4,3) * t315 - t312 * t328 + t321 * t340 + t374;
t412 = t211 * t370;
t275 = -mrSges(6,1) * t296 + mrSges(6,2) * t297;
t396 = -t227 * t362 + t367 * t228;
t214 = m(6) * t237 - mrSges(6,2) * t310 + mrSges(6,3) * t263 + t275 * t296 - t281 * t324 + t396;
t236 = t240 * t368 - t242 * t363;
t233 = -pkin(5) * t310 - pkin(12) * t323 + t276 * t297 - t236;
t382 = -m(7) * t233 + t247 * mrSges(7,1) - mrSges(7,2) * t248 + t278 * t269 - t270 * t279;
t223 = m(6) * t236 + mrSges(6,1) * t310 - mrSges(6,3) * t264 - t275 * t297 + t280 * t324 + t382;
t208 = t363 * t214 + t368 * t223;
t298 = -mrSges(5,1) * t319 + mrSges(5,2) * t320;
t206 = m(5) * t243 + mrSges(5,1) * t311 - mrSges(5,3) * t289 - t298 * t320 + t303 * t326 + t208;
t397 = t368 * t214 - t223 * t363;
t207 = m(5) * t244 - mrSges(5,2) * t311 + mrSges(5,3) * t288 + t298 * t319 - t304 * t326 + t397;
t200 = t369 * t206 + t364 * t207;
t322 = mrSges(4,1) * t340 - mrSges(4,3) * t328;
t398 = -t206 * t364 + t369 * t207;
t197 = m(4) * t267 - mrSges(4,2) * t337 + mrSges(4,3) * t314 + t312 * t327 - t322 * t340 + t398;
t199 = m(4) * t277 - mrSges(4,1) * t314 + mrSges(4,2) * t315 - t321 * t327 + t322 * t328 + t200;
t186 = t197 * t409 + t360 * t199 + t357 * t412;
t187 = t197 * t405 - t357 * t199 + t360 * t412;
t316 = -t356 * t343 + t392;
t395 = -mrSges(3,1) * t359 + mrSges(3,2) * t356;
t341 = t395 * t403;
t389 = -mrSges(3,2) * t361 + mrSges(3,3) * t408;
t346 = t389 * qJD(1);
t390 = mrSges(3,1) * t361 - mrSges(3,3) * t411;
t184 = m(3) * t316 + t390 * qJDD(1) + (-t341 * t411 + t346 * t361) * qJD(1) + t187;
t193 = t370 * t197 - t365 * t211;
t317 = -g(3) * t411 + t401;
t345 = t390 * qJD(1);
t192 = m(3) * t317 + t389 * qJDD(1) + (t341 * t408 - t345 * t361) * qJD(1) + t193;
t181 = -t184 * t356 + t359 * t192;
t329 = -t358 * t342 + t399;
t185 = m(3) * t329 + (t395 * qJDD(1) + (t345 * t356 - t346 * t359) * qJD(1)) * t358 + t186;
t176 = t184 * t406 - t185 * t358 + t192 * t410;
t394 = Ifges(3,5) * t356 + Ifges(3,6) * t359;
t255 = Ifges(7,5) * t279 + Ifges(7,6) * t278 + Ifges(7,3) * t295;
t257 = Ifges(7,1) * t279 + Ifges(7,4) * t278 + Ifges(7,5) * t295;
t220 = -mrSges(7,1) * t233 + mrSges(7,3) * t232 + Ifges(7,4) * t248 + Ifges(7,2) * t247 + Ifges(7,6) * t262 - t255 * t279 + t257 * t295;
t256 = Ifges(7,4) * t279 + Ifges(7,2) * t278 + Ifges(7,6) * t295;
t221 = mrSges(7,2) * t233 - mrSges(7,3) * t231 + Ifges(7,1) * t248 + Ifges(7,4) * t247 + Ifges(7,5) * t262 + t255 * t278 - t256 * t295;
t271 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t324;
t272 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t324;
t201 = mrSges(6,2) * t245 - mrSges(6,3) * t236 + Ifges(6,1) * t264 + Ifges(6,4) * t263 + Ifges(6,5) * t310 - pkin(12) * t216 - t220 * t362 + t221 * t367 + t271 * t296 - t272 * t324;
t273 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t324;
t375 = mrSges(7,1) * t231 - mrSges(7,2) * t232 + Ifges(7,5) * t248 + Ifges(7,6) * t247 + Ifges(7,3) * t262 + t256 * t279 - t257 * t278;
t202 = -mrSges(6,1) * t245 + mrSges(6,3) * t237 + Ifges(6,4) * t264 + Ifges(6,2) * t263 + Ifges(6,6) * t310 - pkin(5) * t216 - t271 * t297 + t273 * t324 - t375;
t282 = Ifges(5,5) * t320 + Ifges(5,6) * t319 + Ifges(5,3) * t326;
t284 = Ifges(5,1) * t320 + Ifges(5,4) * t319 + Ifges(5,5) * t326;
t188 = -mrSges(5,1) * t250 + mrSges(5,3) * t244 + Ifges(5,4) * t289 + Ifges(5,2) * t288 + Ifges(5,6) * t311 - pkin(4) * t378 + pkin(11) * t397 + t363 * t201 + t368 * t202 - t320 * t282 + t326 * t284;
t283 = Ifges(5,4) * t320 + Ifges(5,2) * t319 + Ifges(5,6) * t326;
t189 = mrSges(5,2) * t250 - mrSges(5,3) * t243 + Ifges(5,1) * t289 + Ifges(5,4) * t288 + Ifges(5,5) * t311 - pkin(11) * t208 + t201 * t368 - t202 * t363 + t282 * t319 - t283 * t326;
t306 = Ifges(4,5) * t328 + Ifges(4,6) * t327 + Ifges(4,3) * t340;
t307 = Ifges(4,4) * t328 + Ifges(4,2) * t327 + Ifges(4,6) * t340;
t178 = mrSges(4,2) * t277 - mrSges(4,3) * t266 + Ifges(4,1) * t315 + Ifges(4,4) * t314 + Ifges(4,5) * t337 - pkin(10) * t200 - t188 * t364 + t189 * t369 + t306 * t327 - t307 * t340;
t308 = Ifges(4,1) * t328 + Ifges(4,4) * t327 + Ifges(4,5) * t340;
t376 = -mrSges(6,1) * t236 + mrSges(6,2) * t237 - Ifges(6,5) * t264 - Ifges(6,6) * t263 - Ifges(6,3) * t310 - pkin(5) * t382 - pkin(12) * t396 - t367 * t220 - t362 * t221 - t297 * t272 + t296 * t273;
t373 = mrSges(5,1) * t243 - mrSges(5,2) * t244 + Ifges(5,5) * t289 + Ifges(5,6) * t288 + Ifges(5,3) * t311 + pkin(4) * t208 + t320 * t283 - t319 * t284 - t376;
t182 = -mrSges(4,1) * t277 + mrSges(4,3) * t267 + Ifges(4,4) * t315 + Ifges(4,2) * t314 + Ifges(4,6) * t337 - pkin(3) * t200 - t328 * t306 + t340 * t308 - t373;
t385 = pkin(9) * t193 + t178 * t365 + t182 * t370;
t381 = Ifges(3,5) * t361 + (Ifges(3,1) * t356 + Ifges(3,4) * t359) * t358;
t380 = Ifges(3,6) * t361 + (Ifges(3,4) * t356 + Ifges(3,2) * t359) * t358;
t177 = mrSges(4,1) * t266 - mrSges(4,2) * t267 + Ifges(4,5) * t315 + Ifges(4,6) * t314 + Ifges(4,3) * t337 + pkin(3) * t374 + pkin(10) * t398 + t369 * t188 + t364 * t189 + t328 * t307 - t327 * t308;
t333 = t380 * qJD(1);
t334 = t381 * qJD(1);
t168 = qJDD(1) * t415 + mrSges(3,1) * t316 - mrSges(3,2) * t317 + pkin(2) * t187 + t360 * t177 + t385 * t357 + (t394 * qJDD(1) + (t333 * t356 - t334 * t359) * qJD(1)) * t358;
t332 = (t358 * t394 + t415) * qJD(1);
t170 = -mrSges(3,1) * t329 + mrSges(3,3) * t317 - pkin(2) * t186 - t357 * t177 + (-t332 * t411 + t334 * t361) * qJD(1) + t385 * t360 + t380 * qJDD(1);
t172 = mrSges(3,2) * t329 - mrSges(3,3) * t316 + t370 * t178 - t365 * t182 + (t332 * t408 - t333 * t361) * qJD(1) + (-t186 * t357 - t187 * t360) * pkin(9) + t381 * qJDD(1);
t379 = mrSges(2,1) * t352 - mrSges(2,2) * t353 + Ifges(2,3) * qJDD(1) + pkin(1) * t176 + t361 * t168 + t170 * t408 + t172 * t411 + t181 * t413;
t179 = m(2) * t353 - mrSges(2,1) * t372 - qJDD(1) * mrSges(2,2) + t181;
t175 = t361 * t185 + (t184 * t359 + t192 * t356) * t358;
t173 = m(2) * t352 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t372 + t176;
t166 = -mrSges(2,2) * g(3) - mrSges(2,3) * t352 + Ifges(2,5) * qJDD(1) - t372 * Ifges(2,6) - t356 * t170 + t359 * t172 + (-t175 * t358 - t176 * t361) * qJ(2);
t165 = mrSges(2,1) * g(3) + mrSges(2,3) * t353 + t372 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t175 - t358 * t168 + (qJ(2) * t181 + t170 * t359 + t172 * t356) * t361;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t371 * t166 - t366 * t165 - pkin(8) * (t173 * t371 + t179 * t366), t166, t172, t178, t189, t201, t221; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t366 * t166 + t371 * t165 + pkin(8) * (-t173 * t366 + t179 * t371), t165, t170, t182, t188, t202, t220; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t379, t379, t168, t177, t373, -t376, t375;];
m_new  = t1;
