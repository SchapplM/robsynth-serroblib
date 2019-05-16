% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR13
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
% Datum: 2019-05-07 15:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR13_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 15:33:23
% EndTime: 2019-05-07 15:37:33
% DurationCPUTime: 203.28s
% Computational Cost: add. (3647635->417), mult. (9128128->555), div. (0->0), fcn. (7736931->16), ass. (0->181)
t407 = cos(qJ(3));
t360 = sin(pkin(6));
t406 = pkin(9) * t360;
t359 = sin(pkin(7));
t405 = pkin(10) * t359;
t362 = cos(pkin(7));
t404 = pkin(10) * t362;
t363 = cos(pkin(6));
t403 = g(3) * t363;
t366 = sin(qJ(3));
t402 = t359 * t366;
t367 = sin(qJ(2));
t401 = t360 * t367;
t371 = cos(qJ(2));
t400 = t360 * t371;
t399 = t362 * t366;
t398 = t363 * t367;
t397 = t363 * t371;
t355 = qJD(1) * t363 + qJD(2);
t393 = qJD(1) * t371;
t388 = t360 * t393;
t383 = t362 * t388;
t335 = (t355 * t359 + t383) * pkin(10);
t395 = qJD(1) * t360;
t340 = (-pkin(2) * t371 - t367 * t405) * t395;
t392 = qJD(1) * qJD(2);
t346 = (qJDD(1) * t367 + t371 * t392) * t360;
t354 = qJDD(1) * t363 + qJDD(2);
t368 = sin(qJ(1));
t372 = cos(qJ(1));
t352 = t368 * g(1) - g(2) * t372;
t373 = qJD(1) ^ 2;
t343 = qJDD(1) * pkin(1) + t373 * t406 + t352;
t353 = -g(1) * t372 - g(2) * t368;
t344 = -pkin(1) * t373 + qJDD(1) * t406 + t353;
t384 = t343 * t397 - t344 * t367;
t394 = qJD(1) * t367;
t292 = -t346 * t404 + pkin(2) * t354 + t335 * t355 + (-g(3) * t371 - t340 * t394) * t360 + t384;
t389 = t360 * t394;
t339 = pkin(2) * t355 - t389 * t404;
t347 = (qJDD(1) * t371 - t367 * t392) * t360;
t382 = t347 * t362 + t354 * t359;
t396 = t343 * t398 + t371 * t344;
t293 = -t339 * t355 + (-g(3) * t367 + t340 * t393) * t360 + t382 * pkin(10) + t396;
t300 = -t346 * t405 - pkin(2) * t347 - t403 + (-t343 + (-t335 * t371 + t339 * t367) * qJD(1)) * t360;
t266 = t292 * t399 + t407 * t293 + t300 * t402;
t391 = t359 * t407;
t325 = -t355 * t391 + t366 * t389 - t407 * t383;
t326 = t355 * t402 + (t367 * t407 + t371 * t399) * t395;
t312 = pkin(3) * t325 - qJ(4) * t326;
t327 = -t347 * t359 + t354 * t362 + qJDD(3);
t337 = t355 * t362 - t359 * t388 + qJD(3);
t334 = t337 ^ 2;
t250 = -pkin(3) * t334 + qJ(4) * t327 - t312 * t325 + t266;
t276 = -t292 * t359 + t362 * t300;
t390 = t362 * t407;
t310 = qJD(3) * t326 + t346 * t366 - t347 * t390 - t354 * t391;
t311 = -t325 * qJD(3) + t346 * t407 + t382 * t366;
t253 = (t325 * t337 - t311) * qJ(4) + (t326 * t337 + t310) * pkin(3) + t276;
t358 = sin(pkin(13));
t361 = cos(pkin(13));
t318 = t326 * t361 + t337 * t358;
t242 = -0.2e1 * qJD(4) * t318 - t250 * t358 + t361 * t253;
t297 = t311 * t361 + t327 * t358;
t317 = -t326 * t358 + t337 * t361;
t239 = (t317 * t325 - t297) * pkin(11) + (t317 * t318 + t310) * pkin(4) + t242;
t243 = 0.2e1 * qJD(4) * t317 + t361 * t250 + t358 * t253;
t296 = -t311 * t358 + t327 * t361;
t304 = pkin(4) * t325 - pkin(11) * t318;
t316 = t317 ^ 2;
t241 = -pkin(4) * t316 + pkin(11) * t296 - t304 * t325 + t243;
t365 = sin(qJ(5));
t370 = cos(qJ(5));
t236 = t365 * t239 + t370 * t241;
t291 = t317 * t365 + t318 * t370;
t263 = -qJD(5) * t291 + t296 * t370 - t297 * t365;
t290 = t317 * t370 - t318 * t365;
t274 = -mrSges(6,1) * t290 + mrSges(6,2) * t291;
t324 = qJD(5) + t325;
t280 = mrSges(6,1) * t324 - mrSges(6,3) * t291;
t309 = qJDD(5) + t310;
t275 = -pkin(5) * t290 - pkin(12) * t291;
t323 = t324 ^ 2;
t233 = -pkin(5) * t323 + pkin(12) * t309 + t275 * t290 + t236;
t265 = t292 * t390 - t366 * t293 + t300 * t391;
t249 = -t327 * pkin(3) - t334 * qJ(4) + t326 * t312 + qJDD(4) - t265;
t244 = -t296 * pkin(4) - t316 * pkin(11) + t318 * t304 + t249;
t264 = qJD(5) * t290 + t296 * t365 + t297 * t370;
t237 = (-t290 * t324 - t264) * pkin(12) + (t291 * t324 - t263) * pkin(5) + t244;
t364 = sin(qJ(6));
t369 = cos(qJ(6));
t230 = -t233 * t364 + t237 * t369;
t277 = -t291 * t364 + t324 * t369;
t247 = qJD(6) * t277 + t264 * t369 + t309 * t364;
t262 = qJDD(6) - t263;
t278 = t291 * t369 + t324 * t364;
t267 = -mrSges(7,1) * t277 + mrSges(7,2) * t278;
t289 = qJD(6) - t290;
t268 = -mrSges(7,2) * t289 + mrSges(7,3) * t277;
t226 = m(7) * t230 + mrSges(7,1) * t262 - mrSges(7,3) * t247 - t267 * t278 + t268 * t289;
t231 = t233 * t369 + t237 * t364;
t246 = -qJD(6) * t278 - t264 * t364 + t309 * t369;
t269 = mrSges(7,1) * t289 - mrSges(7,3) * t278;
t227 = m(7) * t231 - mrSges(7,2) * t262 + mrSges(7,3) * t246 + t267 * t277 - t269 * t289;
t385 = -t226 * t364 + t369 * t227;
t210 = m(6) * t236 - mrSges(6,2) * t309 + mrSges(6,3) * t263 + t274 * t290 - t280 * t324 + t385;
t235 = t239 * t370 - t241 * t365;
t279 = -mrSges(6,2) * t324 + mrSges(6,3) * t290;
t232 = -pkin(5) * t309 - pkin(12) * t323 + t275 * t291 - t235;
t380 = -m(7) * t232 + t246 * mrSges(7,1) - mrSges(7,2) * t247 + t277 * t268 - t269 * t278;
t222 = m(6) * t235 + mrSges(6,1) * t309 - mrSges(6,3) * t264 - t274 * t291 + t279 * t324 + t380;
t207 = t365 * t210 + t370 * t222;
t295 = -mrSges(5,1) * t317 + mrSges(5,2) * t318;
t302 = -mrSges(5,2) * t325 + mrSges(5,3) * t317;
t205 = m(5) * t242 + mrSges(5,1) * t310 - mrSges(5,3) * t297 - t295 * t318 + t302 * t325 + t207;
t303 = mrSges(5,1) * t325 - mrSges(5,3) * t318;
t386 = t370 * t210 - t222 * t365;
t206 = m(5) * t243 - mrSges(5,2) * t310 + mrSges(5,3) * t296 + t295 * t317 - t303 * t325 + t386;
t199 = t361 * t205 + t358 * t206;
t215 = t369 * t226 + t364 * t227;
t313 = mrSges(4,1) * t325 + mrSges(4,2) * t326;
t320 = mrSges(4,1) * t337 - mrSges(4,3) * t326;
t387 = -t205 * t358 + t361 * t206;
t196 = m(4) * t266 - mrSges(4,2) * t327 - mrSges(4,3) * t310 - t313 * t325 - t320 * t337 + t387;
t319 = -mrSges(4,2) * t337 - mrSges(4,3) * t325;
t198 = m(4) * t276 + mrSges(4,1) * t310 + mrSges(4,2) * t311 + t319 * t325 + t320 * t326 + t199;
t378 = m(6) * t244 - t263 * mrSges(6,1) + mrSges(6,2) * t264 - t290 * t279 + t280 * t291 + t215;
t375 = -m(5) * t249 + t296 * mrSges(5,1) - mrSges(5,2) * t297 + t317 * t302 - t303 * t318 - t378;
t213 = m(4) * t265 + mrSges(4,1) * t327 - mrSges(4,3) * t311 - t313 * t326 + t319 * t337 + t375;
t185 = t196 * t402 + t362 * t198 + t213 * t391;
t186 = t196 * t399 - t198 * t359 + t213 * t390;
t321 = -g(3) * t400 + t384;
t342 = -mrSges(3,2) * t355 + mrSges(3,3) * t388;
t345 = (-mrSges(3,1) * t371 + mrSges(3,2) * t367) * t395;
t183 = m(3) * t321 + mrSges(3,1) * t354 - mrSges(3,3) * t346 + t342 * t355 - t345 * t389 + t186;
t192 = t407 * t196 - t213 * t366;
t322 = -g(3) * t401 + t396;
t341 = mrSges(3,1) * t355 - mrSges(3,3) * t389;
t191 = m(3) * t322 - mrSges(3,2) * t354 + mrSges(3,3) * t347 - t341 * t355 + t345 * t388 + t192;
t180 = -t183 * t367 + t371 * t191;
t331 = -t343 * t360 - t403;
t184 = m(3) * t331 - mrSges(3,1) * t347 + mrSges(3,2) * t346 + (t341 * t367 - t342 * t371) * t395 + t185;
t175 = t183 * t397 - t184 * t360 + t191 * t398;
t254 = Ifges(7,5) * t278 + Ifges(7,6) * t277 + Ifges(7,3) * t289;
t256 = Ifges(7,1) * t278 + Ifges(7,4) * t277 + Ifges(7,5) * t289;
t219 = -mrSges(7,1) * t232 + mrSges(7,3) * t231 + Ifges(7,4) * t247 + Ifges(7,2) * t246 + Ifges(7,6) * t262 - t254 * t278 + t256 * t289;
t255 = Ifges(7,4) * t278 + Ifges(7,2) * t277 + Ifges(7,6) * t289;
t220 = mrSges(7,2) * t232 - mrSges(7,3) * t230 + Ifges(7,1) * t247 + Ifges(7,4) * t246 + Ifges(7,5) * t262 + t254 * t277 - t255 * t289;
t270 = Ifges(6,5) * t291 + Ifges(6,6) * t290 + Ifges(6,3) * t324;
t271 = Ifges(6,4) * t291 + Ifges(6,2) * t290 + Ifges(6,6) * t324;
t200 = mrSges(6,2) * t244 - mrSges(6,3) * t235 + Ifges(6,1) * t264 + Ifges(6,4) * t263 + Ifges(6,5) * t309 - pkin(12) * t215 - t219 * t364 + t220 * t369 + t270 * t290 - t271 * t324;
t272 = Ifges(6,1) * t291 + Ifges(6,4) * t290 + Ifges(6,5) * t324;
t376 = mrSges(7,1) * t230 - mrSges(7,2) * t231 + Ifges(7,5) * t247 + Ifges(7,6) * t246 + Ifges(7,3) * t262 + t255 * t278 - t256 * t277;
t201 = -mrSges(6,1) * t244 + mrSges(6,3) * t236 + Ifges(6,4) * t264 + Ifges(6,2) * t263 + Ifges(6,6) * t309 - pkin(5) * t215 - t270 * t291 + t272 * t324 - t376;
t281 = Ifges(5,5) * t318 + Ifges(5,6) * t317 + Ifges(5,3) * t325;
t283 = Ifges(5,1) * t318 + Ifges(5,4) * t317 + Ifges(5,5) * t325;
t187 = -mrSges(5,1) * t249 + mrSges(5,3) * t243 + Ifges(5,4) * t297 + Ifges(5,2) * t296 + Ifges(5,6) * t310 - pkin(4) * t378 + pkin(11) * t386 + t365 * t200 + t370 * t201 - t318 * t281 + t325 * t283;
t282 = Ifges(5,4) * t318 + Ifges(5,2) * t317 + Ifges(5,6) * t325;
t188 = mrSges(5,2) * t249 - mrSges(5,3) * t242 + Ifges(5,1) * t297 + Ifges(5,4) * t296 + Ifges(5,5) * t310 - pkin(11) * t207 + t200 * t370 - t201 * t365 + t281 * t317 - t282 * t325;
t305 = Ifges(4,5) * t326 - Ifges(4,6) * t325 + Ifges(4,3) * t337;
t306 = Ifges(4,4) * t326 - Ifges(4,2) * t325 + Ifges(4,6) * t337;
t177 = mrSges(4,2) * t276 - mrSges(4,3) * t265 + Ifges(4,1) * t311 - Ifges(4,4) * t310 + Ifges(4,5) * t327 - qJ(4) * t199 - t187 * t358 + t188 * t361 - t305 * t325 - t306 * t337;
t307 = Ifges(4,1) * t326 - Ifges(4,4) * t325 + Ifges(4,5) * t337;
t377 = -mrSges(6,1) * t235 + mrSges(6,2) * t236 - Ifges(6,5) * t264 - Ifges(6,6) * t263 - Ifges(6,3) * t309 - pkin(5) * t380 - pkin(12) * t385 - t369 * t219 - t364 * t220 - t291 * t271 + t290 * t272;
t374 = -mrSges(5,1) * t242 + mrSges(5,2) * t243 - Ifges(5,5) * t297 - Ifges(5,6) * t296 - pkin(4) * t207 - t318 * t282 + t317 * t283 + t377;
t181 = (-Ifges(5,3) - Ifges(4,2)) * t310 - pkin(3) * t199 + t337 * t307 - t326 * t305 + Ifges(4,6) * t327 + Ifges(4,4) * t311 - mrSges(4,1) * t276 + mrSges(4,3) * t266 + t374;
t381 = pkin(10) * t192 + t177 * t366 + t181 * t407;
t176 = mrSges(4,1) * t265 - mrSges(4,2) * t266 + Ifges(4,5) * t311 - Ifges(4,6) * t310 + Ifges(4,3) * t327 + pkin(3) * t375 + qJ(4) * t387 + t361 * t187 + t358 * t188 + t326 * t306 + t325 * t307;
t329 = Ifges(3,6) * t355 + (Ifges(3,4) * t367 + Ifges(3,2) * t371) * t395;
t330 = Ifges(3,5) * t355 + (Ifges(3,1) * t367 + Ifges(3,4) * t371) * t395;
t167 = mrSges(3,1) * t321 - mrSges(3,2) * t322 + Ifges(3,5) * t346 + Ifges(3,6) * t347 + Ifges(3,3) * t354 + pkin(2) * t186 + t362 * t176 + (t329 * t367 - t330 * t371) * t395 + t381 * t359;
t328 = Ifges(3,3) * t355 + (Ifges(3,5) * t367 + Ifges(3,6) * t371) * t395;
t169 = -mrSges(3,1) * t331 + mrSges(3,3) * t322 + Ifges(3,4) * t346 + Ifges(3,2) * t347 + Ifges(3,6) * t354 - pkin(2) * t185 - t359 * t176 - t328 * t389 + t355 * t330 + t362 * t381;
t171 = t328 * t388 + mrSges(3,2) * t331 - mrSges(3,3) * t321 + t407 * t177 + Ifges(3,1) * t346 + Ifges(3,4) * t347 + Ifges(3,5) * t354 - t366 * t181 - t355 * t329 + (-t185 * t359 - t186 * t362) * pkin(10);
t379 = mrSges(2,1) * t352 - mrSges(2,2) * t353 + Ifges(2,3) * qJDD(1) + pkin(1) * t175 + t363 * t167 + t169 * t400 + t171 * t401 + t180 * t406;
t178 = m(2) * t353 - mrSges(2,1) * t373 - qJDD(1) * mrSges(2,2) + t180;
t174 = t184 * t363 + (t183 * t371 + t191 * t367) * t360;
t172 = m(2) * t352 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t373 + t175;
t165 = -mrSges(2,2) * g(3) - mrSges(2,3) * t352 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t373 - t169 * t367 + t171 * t371 + (-t174 * t360 - t175 * t363) * pkin(9);
t164 = mrSges(2,1) * g(3) + mrSges(2,3) * t353 + Ifges(2,5) * t373 + Ifges(2,6) * qJDD(1) - pkin(1) * t174 - t167 * t360 + (pkin(9) * t180 + t169 * t371 + t171 * t367) * t363;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t372 * t165 - t368 * t164 - pkin(8) * (t172 * t372 + t178 * t368), t165, t171, t177, t188, t200, t220; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t368 * t165 + t372 * t164 + pkin(8) * (-t172 * t368 + t178 * t372), t164, t169, t181, t187, t201, t219; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t379, t379, t167, t176, Ifges(5,3) * t310 - t374, -t377, t376;];
m_new  = t1;
