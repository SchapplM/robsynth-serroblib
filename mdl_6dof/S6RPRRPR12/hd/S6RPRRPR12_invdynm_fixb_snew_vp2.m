% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-05-06 00:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRPR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:41:30
% EndTime: 2019-05-06 00:42:49
% DurationCPUTime: 64.29s
% Computational Cost: add. (1018707->407), mult. (3177569->532), div. (0->0), fcn. (2694143->14), ass. (0->180)
t353 = sin(pkin(12));
t356 = cos(pkin(12));
t361 = sin(qJ(3));
t357 = cos(pkin(7));
t364 = cos(qJ(3));
t403 = t357 * t364;
t423 = -t353 * t361 + t356 * t403;
t355 = sin(pkin(6));
t358 = cos(pkin(6));
t404 = t357 * t361;
t354 = sin(pkin(7));
t409 = t354 * t361;
t374 = t358 * t409 + (t353 * t364 + t356 * t404) * t355;
t323 = t374 * qJD(1);
t408 = t354 * t364;
t398 = t358 * t408;
t309 = -t323 * qJD(3) + (t355 * t423 + t398) * qJDD(1);
t406 = t355 * t357;
t336 = (t354 * t358 + t356 * t406) * qJD(1) * pkin(9);
t362 = sin(qJ(1));
t365 = cos(qJ(1));
t350 = -t365 * g(1) - t362 * g(2);
t366 = qJD(1) ^ 2;
t414 = qJ(2) * t355;
t340 = -t366 * pkin(1) + qJDD(1) * t414 + t350;
t418 = pkin(9) * t353;
t388 = -pkin(2) * t356 - t354 * t418;
t400 = qJD(1) * t355;
t415 = pkin(9) * qJDD(1);
t383 = qJD(1) * t388 * t400 + t357 * t415;
t349 = t362 * g(1) - t365 * g(2);
t339 = qJDD(1) * pkin(1) + t366 * t414 + t349;
t395 = qJD(2) * t400;
t405 = t356 * t358;
t407 = t355 * t356;
t389 = -g(3) * t407 + t339 * t405 - 0.2e1 * t353 * t395;
t284 = (pkin(2) * qJDD(1) + qJD(1) * t336) * t358 + (-t383 * t355 - t340) * t353 + t389;
t341 = (pkin(2) * t358 - t406 * t418) * qJD(1);
t411 = t353 * t358;
t396 = t339 * t411 + (t340 + 0.2e1 * t395) * t356;
t285 = (-qJD(1) * t341 + t354 * t415) * t358 + (-g(3) * t353 + t383 * t356) * t355 + t396;
t394 = -t358 * g(3) + qJDD(2);
t294 = (-t339 + t388 * qJDD(1) + (-t336 * t356 + t341 * t353) * qJD(1)) * t355 + t394;
t244 = -t361 * t285 + (t284 * t357 + t294 * t354) * t364;
t245 = t284 * t404 + t364 * t285 + t294 * t409;
t322 = qJD(1) * t398 + t400 * t423;
t308 = -t322 * pkin(3) - t323 * pkin(10);
t384 = -t354 * t407 + t357 * t358;
t337 = t384 * qJD(1) + qJD(3);
t333 = t337 ^ 2;
t334 = t384 * qJDD(1) + qJDD(3);
t241 = -t333 * pkin(3) + t334 * pkin(10) + t322 * t308 + t245;
t256 = -t354 * t284 + t357 * t294;
t310 = t322 * qJD(3) + t374 * qJDD(1);
t243 = (-t322 * t337 - t310) * pkin(10) + (t323 * t337 - t309) * pkin(3) + t256;
t360 = sin(qJ(4));
t419 = cos(qJ(4));
t236 = -t360 * t241 + t419 * t243;
t237 = t419 * t241 + t360 * t243;
t316 = t360 * t323 - t419 * t337;
t317 = t419 * t323 + t360 * t337;
t321 = qJD(4) - t322;
t267 = Ifges(5,4) * t317 - Ifges(5,2) * t316 + Ifges(5,6) * t321;
t278 = t317 * qJD(4) + t360 * t310 - t419 * t334;
t279 = -t316 * qJD(4) + t419 * t310 + t360 * t334;
t288 = -t316 * mrSges(6,2) - t317 * mrSges(6,3);
t295 = t316 * mrSges(6,1) - t321 * mrSges(6,3);
t306 = qJDD(4) - t309;
t286 = t316 * pkin(4) - t317 * qJ(5);
t320 = t321 ^ 2;
t234 = -t306 * pkin(4) - t320 * qJ(5) + t317 * t286 + qJDD(5) - t236;
t413 = t316 * t321;
t228 = (t316 * t317 - t306) * pkin(11) + (t279 + t413) * pkin(5) + t234;
t299 = t317 * pkin(5) - t321 * pkin(11);
t315 = t316 ^ 2;
t240 = -t334 * pkin(3) - t333 * pkin(10) + t323 * t308 - t244;
t420 = -2 * qJD(5);
t368 = (-t279 + t413) * qJ(5) + t240 + (t321 * pkin(4) + t420) * t317;
t231 = -t315 * pkin(5) - t317 * t299 + (pkin(4) + pkin(11)) * t278 + t368;
t359 = sin(qJ(6));
t363 = cos(qJ(6));
t225 = t363 * t228 - t359 * t231;
t292 = t363 * t316 - t359 * t321;
t251 = t292 * qJD(6) + t359 * t278 + t363 * t306;
t293 = t359 * t316 + t363 * t321;
t258 = -t292 * mrSges(7,1) + t293 * mrSges(7,2);
t312 = qJD(6) + t317;
t261 = -t312 * mrSges(7,2) + t292 * mrSges(7,3);
t275 = qJDD(6) + t279;
t222 = m(7) * t225 + t275 * mrSges(7,1) - t251 * mrSges(7,3) - t293 * t258 + t312 * t261;
t226 = t359 * t228 + t363 * t231;
t250 = -t293 * qJD(6) + t363 * t278 - t359 * t306;
t262 = t312 * mrSges(7,1) - t293 * mrSges(7,3);
t223 = m(7) * t226 - t275 * mrSges(7,2) + t250 * mrSges(7,3) + t292 * t258 - t312 * t262;
t211 = t363 * t222 + t359 * t223;
t377 = -t320 * pkin(4) + t306 * qJ(5) - t316 * t286 + t237;
t230 = -t278 * pkin(5) - t315 * pkin(11) + ((2 * qJD(5)) + t299) * t321 + t377;
t252 = Ifges(7,5) * t293 + Ifges(7,6) * t292 + Ifges(7,3) * t312;
t254 = Ifges(7,1) * t293 + Ifges(7,4) * t292 + Ifges(7,5) * t312;
t214 = -mrSges(7,1) * t230 + mrSges(7,3) * t226 + Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t275 - t293 * t252 + t312 * t254;
t253 = Ifges(7,4) * t293 + Ifges(7,2) * t292 + Ifges(7,6) * t312;
t215 = mrSges(7,2) * t230 - mrSges(7,3) * t225 + Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t275 + t292 * t252 - t312 * t253;
t232 = t321 * t420 - t377;
t264 = Ifges(6,5) * t321 - Ifges(6,6) * t317 + Ifges(6,3) * t316;
t372 = -mrSges(6,2) * t234 + mrSges(6,3) * t232 - Ifges(6,1) * t306 + Ifges(6,4) * t279 - Ifges(6,5) * t278 + pkin(11) * t211 + t359 * t214 - t363 * t215 + t317 * t264;
t227 = -m(7) * t230 + t250 * mrSges(7,1) - t251 * mrSges(7,2) + t292 * t261 - t293 * t262;
t296 = t317 * mrSges(6,1) + t321 * mrSges(6,2);
t373 = -m(6) * t232 + t306 * mrSges(6,3) + t321 * t296 - t227;
t378 = -m(6) * t234 - t279 * mrSges(6,1) - t317 * t288 - t211;
t266 = Ifges(6,4) * t321 - Ifges(6,2) * t317 + Ifges(6,6) * t316;
t401 = Ifges(5,1) * t317 - Ifges(5,4) * t316 + Ifges(5,5) * t321 - t266;
t422 = t401 * t316 + mrSges(5,1) * t236 - mrSges(5,2) * t237 + Ifges(5,5) * t279 - Ifges(5,6) * t278 + Ifges(5,3) * t306 + pkin(4) * (-t306 * mrSges(6,2) - t321 * t295 + t378) + qJ(5) * (-t278 * mrSges(6,1) - t316 * t288 + t373) + t317 * t267 - t372;
t417 = Ifges(5,4) + Ifges(6,6);
t416 = Ifges(3,3) * t358;
t412 = t353 * t355;
t287 = t316 * mrSges(5,1) + t317 * mrSges(5,2);
t297 = -t321 * mrSges(5,2) - t316 * mrSges(5,3);
t208 = m(5) * t236 - t279 * mrSges(5,3) - t317 * t287 + (-t295 + t297) * t321 + (mrSges(5,1) - mrSges(6,2)) * t306 + t378;
t298 = t321 * mrSges(5,1) - t317 * mrSges(5,3);
t218 = m(5) * t237 - t306 * mrSges(5,2) - t321 * t298 + (-t287 - t288) * t316 + (-mrSges(5,3) - mrSges(6,1)) * t278 + t373;
t203 = t419 * t208 + t360 * t218;
t212 = -t359 * t222 + t363 * t223;
t268 = Ifges(6,1) * t321 - Ifges(6,4) * t317 + Ifges(6,5) * t316;
t402 = -Ifges(5,5) * t317 + Ifges(5,6) * t316 - Ifges(5,3) * t321 - t268;
t307 = -t322 * mrSges(4,1) + t323 * mrSges(4,2);
t319 = t337 * mrSges(4,1) - t323 * mrSges(4,3);
t393 = -t360 * t208 + t419 * t218;
t200 = m(4) * t245 - t334 * mrSges(4,2) + t309 * mrSges(4,3) + t322 * t307 - t337 * t319 + t393;
t318 = -t337 * mrSges(4,2) + t322 * mrSges(4,3);
t202 = m(4) * t256 - t309 * mrSges(4,1) + t310 * mrSges(4,2) - t322 * t318 + t323 * t319 + t203;
t235 = t278 * pkin(4) + t368;
t382 = -m(6) * t235 + t278 * mrSges(6,2) + t316 * t295 - t212;
t369 = -m(5) * t240 - t278 * mrSges(5,1) - t316 * t297 + (t296 - t298) * t317 + (-mrSges(5,2) + mrSges(6,3)) * t279 + t382;
t206 = m(4) * t244 + t334 * mrSges(4,1) - t310 * mrSges(4,3) - t323 * t307 + t337 * t318 + t369;
t189 = t200 * t409 + t357 * t202 + t206 * t408;
t190 = t200 * t404 - t354 * t202 + t206 * t403;
t313 = -t353 * t340 + t389;
t392 = -mrSges(3,1) * t356 + mrSges(3,2) * t353;
t338 = t392 * t400;
t386 = -mrSges(3,2) * t358 + mrSges(3,3) * t407;
t343 = t386 * qJD(1);
t387 = mrSges(3,1) * t358 - mrSges(3,3) * t412;
t187 = m(3) * t313 + t387 * qJDD(1) + (-t338 * t412 + t343 * t358) * qJD(1) + t190;
t196 = t364 * t200 - t361 * t206;
t314 = -g(3) * t412 + t396;
t342 = t387 * qJD(1);
t195 = m(3) * t314 + t386 * qJDD(1) + (t338 * t407 - t342 * t358) * qJD(1) + t196;
t184 = -t353 * t187 + t356 * t195;
t324 = -t355 * t339 + t394;
t188 = m(3) * t324 + (t392 * qJDD(1) + (t342 * t353 - t343 * t356) * qJD(1)) * t355 + t189;
t179 = t187 * t405 - t355 * t188 + t195 * t411;
t391 = Ifges(3,5) * t353 + Ifges(3,6) * t356;
t210 = -t279 * mrSges(6,3) - t317 * t296 - t382;
t371 = -mrSges(6,1) * t232 + mrSges(6,2) * t235 - pkin(5) * t227 - pkin(11) * t212 - t363 * t214 - t359 * t215;
t191 = -mrSges(5,1) * t240 + mrSges(5,3) * t237 - pkin(4) * t210 + t401 * t321 + t402 * t317 + (Ifges(5,6) - Ifges(6,5)) * t306 + t417 * t279 + (-Ifges(5,2) - Ifges(6,3)) * t278 + t371;
t375 = mrSges(7,1) * t225 - mrSges(7,2) * t226 + Ifges(7,5) * t251 + Ifges(7,6) * t250 + Ifges(7,3) * t275 + t293 * t253 - t292 * t254;
t370 = mrSges(6,1) * t234 - mrSges(6,3) * t235 + pkin(5) * t211 + t375;
t192 = t370 + (-t267 + t264) * t321 + t402 * t316 + (Ifges(5,5) - Ifges(6,4)) * t306 + (Ifges(5,1) + Ifges(6,2)) * t279 - t417 * t278 - qJ(5) * t210 + mrSges(5,2) * t240 - mrSges(5,3) * t236;
t301 = Ifges(4,5) * t323 + Ifges(4,6) * t322 + Ifges(4,3) * t337;
t302 = Ifges(4,4) * t323 + Ifges(4,2) * t322 + Ifges(4,6) * t337;
t181 = mrSges(4,2) * t256 - mrSges(4,3) * t244 + Ifges(4,1) * t310 + Ifges(4,4) * t309 + Ifges(4,5) * t334 - pkin(10) * t203 - t360 * t191 + t419 * t192 + t322 * t301 - t337 * t302;
t303 = Ifges(4,1) * t323 + Ifges(4,4) * t322 + Ifges(4,5) * t337;
t185 = -mrSges(4,1) * t256 + mrSges(4,3) * t245 + Ifges(4,4) * t310 + Ifges(4,2) * t309 + Ifges(4,6) * t334 - pkin(3) * t203 - t323 * t301 + t337 * t303 - t422;
t381 = pkin(9) * t196 + t181 * t361 + t185 * t364;
t380 = Ifges(3,5) * t358 + (Ifges(3,1) * t353 + Ifges(3,4) * t356) * t355;
t379 = Ifges(3,6) * t358 + (Ifges(3,4) * t353 + Ifges(3,2) * t356) * t355;
t180 = mrSges(4,1) * t244 - mrSges(4,2) * t245 + Ifges(4,5) * t310 + Ifges(4,6) * t309 + Ifges(4,3) * t334 + pkin(3) * t369 + pkin(10) * t393 + t419 * t191 + t360 * t192 + t323 * t302 - t322 * t303;
t328 = t379 * qJD(1);
t329 = t380 * qJD(1);
t171 = qJDD(1) * t416 + mrSges(3,1) * t313 - mrSges(3,2) * t314 + pkin(2) * t190 + t357 * t180 + t381 * t354 + (t391 * qJDD(1) + (t328 * t353 - t329 * t356) * qJD(1)) * t355;
t327 = (t391 * t355 + t416) * qJD(1);
t173 = -mrSges(3,1) * t324 + mrSges(3,3) * t314 - pkin(2) * t189 - t354 * t180 + (-t327 * t412 + t329 * t358) * qJD(1) + t381 * t357 + t379 * qJDD(1);
t175 = mrSges(3,2) * t324 - mrSges(3,3) * t313 + t364 * t181 - t361 * t185 + (t327 * t407 - t328 * t358) * qJD(1) + (-t189 * t354 - t190 * t357) * pkin(9) + t380 * qJDD(1);
t376 = mrSges(2,1) * t349 - mrSges(2,2) * t350 + Ifges(2,3) * qJDD(1) + pkin(1) * t179 + t358 * t171 + t173 * t407 + t175 * t412 + t184 * t414;
t182 = m(2) * t350 - t366 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t184;
t178 = t358 * t188 + (t187 * t356 + t195 * t353) * t355;
t176 = m(2) * t349 + qJDD(1) * mrSges(2,1) - t366 * mrSges(2,2) + t179;
t169 = -mrSges(2,2) * g(3) - mrSges(2,3) * t349 + Ifges(2,5) * qJDD(1) - t366 * Ifges(2,6) - t353 * t173 + t356 * t175 + (-t178 * t355 - t179 * t358) * qJ(2);
t168 = mrSges(2,1) * g(3) + mrSges(2,3) * t350 + t366 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t178 - t355 * t171 + (qJ(2) * t184 + t173 * t356 + t175 * t353) * t358;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t365 * t169 - t362 * t168 - pkin(8) * (t365 * t176 + t362 * t182), t169, t175, t181, t192, -t316 * t266 - t372, t215; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t362 * t169 + t365 * t168 + pkin(8) * (-t362 * t176 + t365 * t182), t168, t173, t185, t191, Ifges(6,4) * t306 - Ifges(6,2) * t279 + Ifges(6,6) * t278 - t321 * t264 + t316 * t268 - t370, t214; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t376, t376, t171, t180, t422, Ifges(6,5) * t306 - Ifges(6,6) * t279 + Ifges(6,3) * t278 + t321 * t266 + t317 * t268 - t371, t375;];
m_new  = t1;
