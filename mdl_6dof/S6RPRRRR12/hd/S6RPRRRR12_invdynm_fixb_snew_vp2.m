% Calculate vector of cutting torques with Newton-Euler for
% S6RPRRRR12
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-05-06 07:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRRRR12_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR12_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 06:24:05
% EndTime: 2019-05-06 06:31:38
% DurationCPUTime: 448.51s
% Computational Cost: add. (7279886->416), mult. (23743554->580), div. (0->0), fcn. (20876754->18), ass. (0->198)
t358 = sin(pkin(7));
t360 = cos(pkin(14));
t363 = cos(pkin(6));
t359 = sin(pkin(6));
t362 = cos(pkin(7));
t406 = t359 * t362;
t341 = (t358 * t363 + t360 * t406) * qJD(1) * pkin(10);
t368 = sin(qJ(1));
t373 = cos(qJ(1));
t353 = -g(1) * t373 - g(2) * t368;
t374 = qJD(1) ^ 2;
t415 = qJ(2) * t359;
t345 = -pkin(1) * t374 + qJDD(1) * t415 + t353;
t356 = sin(pkin(14));
t420 = pkin(10) * t356;
t389 = -pkin(2) * t360 - t358 * t420;
t401 = qJD(1) * t359;
t416 = pkin(10) * qJDD(1);
t385 = qJD(1) * t389 * t401 + t362 * t416;
t352 = t368 * g(1) - g(2) * t373;
t344 = qJDD(1) * pkin(1) + t374 * t415 + t352;
t399 = qJD(2) * t401;
t405 = t360 * t363;
t407 = t359 * t360;
t390 = -g(3) * t407 + t344 * t405 - 0.2e1 * t356 * t399;
t304 = (pkin(2) * qJDD(1) + qJD(1) * t341) * t363 + (-t385 * t359 - t345) * t356 + t390;
t346 = (pkin(2) * t363 - t406 * t420) * qJD(1);
t411 = t356 * t363;
t400 = t344 * t411 + (t345 + 0.2e1 * t399) * t360;
t305 = (-qJD(1) * t346 + t358 * t416) * t363 + (-g(3) * t356 + t385 * t360) * t359 + t400;
t398 = -g(3) * t363 + qJDD(2);
t313 = (-t344 + t389 * qJDD(1) + (-t341 * t360 + t346 * t356) * qJD(1)) * t359 + t398;
t367 = sin(qJ(3));
t372 = cos(qJ(3));
t402 = t362 * t372;
t408 = t358 * t372;
t278 = t304 * t402 - t305 * t367 + t313 * t408;
t378 = t363 * t408 + (-t356 * t367 + t360 * t402) * t359;
t330 = t378 * qJD(1);
t403 = t362 * t367;
t409 = t358 * t367;
t379 = t363 * t409 + (t356 * t372 + t360 * t403) * t359;
t331 = t379 * qJD(1);
t357 = sin(pkin(8));
t419 = pkin(11) * t357;
t318 = -pkin(3) * t330 - t331 * t419;
t321 = qJD(3) * t330 + t379 * qJDD(1);
t386 = -t358 * t407 + t362 * t363;
t342 = t386 * qJD(1) + qJD(3);
t361 = cos(pkin(8));
t392 = t330 * t361 + t342 * t357;
t323 = t392 * pkin(11);
t339 = t386 * qJDD(1) + qJDD(3);
t418 = pkin(11) * t361;
t262 = pkin(3) * t339 - t318 * t331 - t321 * t418 + t323 * t342 + t278;
t279 = t304 * t403 + t372 * t305 + t313 * t409;
t327 = pkin(3) * t342 - t331 * t418;
t320 = -qJD(3) * t331 + t378 * qJDD(1);
t393 = t320 * t361 + t339 * t357;
t263 = t393 * pkin(11) + t318 * t330 - t327 * t342 + t279;
t293 = -t304 * t358 + t362 * t313;
t266 = -pkin(3) * t320 - t321 * t419 - t323 * t330 + t327 * t331 + t293;
t366 = sin(qJ(4));
t371 = cos(qJ(4));
t249 = -t366 * t263 + (t262 * t361 + t266 * t357) * t371;
t312 = t331 * t371 + t392 * t366;
t286 = -qJD(4) * t312 - t321 * t366 + t393 * t371;
t311 = -t331 * t366 + t392 * t371;
t417 = Ifges(3,3) * t363;
t287 = qJD(4) * t311 + t321 * t371 + t393 * t366;
t294 = -mrSges(5,1) * t311 + mrSges(5,2) * t312;
t324 = -t330 * t357 + t342 * t361 + qJD(4);
t299 = -mrSges(5,2) * t324 + mrSges(5,3) * t311;
t314 = -t320 * t357 + t339 * t361 + qJDD(4);
t404 = t361 * t366;
t410 = t357 * t366;
t250 = t262 * t404 + t371 * t263 + t266 * t410;
t295 = -pkin(4) * t311 - pkin(12) * t312;
t322 = t324 ^ 2;
t246 = -pkin(4) * t322 + pkin(12) * t314 + t295 * t311 + t250;
t251 = -t262 * t357 + t361 * t266;
t248 = (-t311 * t324 - t287) * pkin(12) + (t312 * t324 - t286) * pkin(4) + t251;
t365 = sin(qJ(5));
t370 = cos(qJ(5));
t242 = t370 * t246 + t365 * t248;
t297 = -t312 * t365 + t370 * t324;
t298 = t312 * t370 + t324 * t365;
t281 = -pkin(5) * t297 - pkin(13) * t298;
t285 = qJDD(5) - t286;
t309 = qJD(5) - t311;
t308 = t309 ^ 2;
t240 = -pkin(5) * t308 + pkin(13) * t285 + t281 * t297 + t242;
t245 = -pkin(4) * t314 - pkin(12) * t322 + t312 * t295 - t249;
t270 = -qJD(5) * t298 - t287 * t365 + t314 * t370;
t271 = qJD(5) * t297 + t287 * t370 + t314 * t365;
t243 = (-t297 * t309 - t271) * pkin(13) + (t298 * t309 - t270) * pkin(5) + t245;
t364 = sin(qJ(6));
t369 = cos(qJ(6));
t236 = -t240 * t364 + t243 * t369;
t283 = -t298 * t364 + t309 * t369;
t254 = qJD(6) * t283 + t271 * t369 + t285 * t364;
t284 = t298 * t369 + t309 * t364;
t267 = -mrSges(7,1) * t283 + mrSges(7,2) * t284;
t269 = qJDD(6) - t270;
t296 = qJD(6) - t297;
t272 = -mrSges(7,2) * t296 + mrSges(7,3) * t283;
t234 = m(7) * t236 + mrSges(7,1) * t269 - mrSges(7,3) * t254 - t267 * t284 + t272 * t296;
t237 = t240 * t369 + t243 * t364;
t253 = -qJD(6) * t284 - t271 * t364 + t285 * t369;
t273 = mrSges(7,1) * t296 - mrSges(7,3) * t284;
t235 = m(7) * t237 - mrSges(7,2) * t269 + mrSges(7,3) * t253 + t267 * t283 - t273 * t296;
t227 = t234 * t369 + t235 * t364;
t288 = -mrSges(6,2) * t309 + mrSges(6,3) * t297;
t289 = mrSges(6,1) * t309 - mrSges(6,3) * t298;
t377 = -m(6) * t245 + t270 * mrSges(6,1) - mrSges(6,2) * t271 + t297 * t288 - t289 * t298 - t227;
t223 = m(5) * t249 + mrSges(5,1) * t314 - mrSges(5,3) * t287 - t294 * t312 + t299 * t324 + t377;
t414 = t223 * t371;
t412 = t356 * t359;
t228 = -t234 * t364 + t369 * t235;
t280 = -mrSges(6,1) * t297 + mrSges(6,2) * t298;
t226 = m(6) * t242 - mrSges(6,2) * t285 + mrSges(6,3) * t270 + t280 * t297 - t289 * t309 + t228;
t241 = -t246 * t365 + t248 * t370;
t239 = -pkin(5) * t285 - pkin(13) * t308 + t281 * t298 - t241;
t238 = -m(7) * t239 + t253 * mrSges(7,1) - mrSges(7,2) * t254 + t283 * t272 - t273 * t284;
t232 = m(6) * t241 + mrSges(6,1) * t285 - mrSges(6,3) * t271 - t280 * t298 + t288 * t309 + t238;
t220 = t365 * t226 + t370 * t232;
t300 = mrSges(5,1) * t324 - mrSges(5,3) * t312;
t397 = t370 * t226 - t365 * t232;
t217 = m(5) * t250 - mrSges(5,2) * t314 + mrSges(5,3) * t286 + t294 * t311 - t300 * t324 + t397;
t219 = m(5) * t251 - mrSges(5,1) * t286 + mrSges(5,2) * t287 - t299 * t311 + t300 * t312 + t220;
t206 = t217 * t404 - t219 * t357 + t361 * t414;
t319 = -mrSges(4,1) * t330 + mrSges(4,2) * t331;
t328 = -mrSges(4,2) * t342 + mrSges(4,3) * t330;
t202 = m(4) * t278 + mrSges(4,1) * t339 - mrSges(4,3) * t321 - t319 * t331 + t328 * t342 + t206;
t205 = t217 * t410 + t361 * t219 + t357 * t414;
t329 = mrSges(4,1) * t342 - mrSges(4,3) * t331;
t204 = m(4) * t293 - mrSges(4,1) * t320 + mrSges(4,2) * t321 - t328 * t330 + t329 * t331 + t205;
t211 = t371 * t217 - t223 * t366;
t210 = m(4) * t279 - mrSges(4,2) * t339 + mrSges(4,3) * t320 + t319 * t330 - t329 * t342 + t211;
t191 = t202 * t408 + t362 * t204 + t210 * t409;
t192 = t202 * t402 - t204 * t358 + t210 * t403;
t325 = -t345 * t356 + t390;
t396 = -mrSges(3,1) * t360 + mrSges(3,2) * t356;
t343 = t396 * t401;
t387 = -mrSges(3,2) * t363 + mrSges(3,3) * t407;
t348 = t387 * qJD(1);
t388 = mrSges(3,1) * t363 - mrSges(3,3) * t412;
t189 = m(3) * t325 + t388 * qJDD(1) + (-t343 * t412 + t348 * t363) * qJD(1) + t192;
t196 = -t202 * t367 + t372 * t210;
t326 = -g(3) * t412 + t400;
t347 = t388 * qJD(1);
t195 = m(3) * t326 + t387 * qJDD(1) + (t343 * t407 - t347 * t363) * qJD(1) + t196;
t187 = -t189 * t356 + t360 * t195;
t332 = -t344 * t359 + t398;
t190 = m(3) * t332 + (t396 * qJDD(1) + (t347 * t356 - t348 * t360) * qJD(1)) * t359 + t191;
t181 = t189 * t405 - t190 * t359 + t195 * t411;
t395 = Ifges(3,5) * t356 + Ifges(3,6) * t360;
t255 = Ifges(7,5) * t284 + Ifges(7,6) * t283 + Ifges(7,3) * t296;
t257 = Ifges(7,1) * t284 + Ifges(7,4) * t283 + Ifges(7,5) * t296;
t229 = -mrSges(7,1) * t239 + mrSges(7,3) * t237 + Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t269 - t255 * t284 + t257 * t296;
t256 = Ifges(7,4) * t284 + Ifges(7,2) * t283 + Ifges(7,6) * t296;
t230 = mrSges(7,2) * t239 - mrSges(7,3) * t236 + Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t269 + t255 * t283 - t256 * t296;
t274 = Ifges(6,5) * t298 + Ifges(6,6) * t297 + Ifges(6,3) * t309;
t275 = Ifges(6,4) * t298 + Ifges(6,2) * t297 + Ifges(6,6) * t309;
t212 = mrSges(6,2) * t245 - mrSges(6,3) * t241 + Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t285 - pkin(13) * t227 - t229 * t364 + t230 * t369 + t274 * t297 - t275 * t309;
t276 = Ifges(6,1) * t298 + Ifges(6,4) * t297 + Ifges(6,5) * t309;
t376 = mrSges(7,1) * t236 - mrSges(7,2) * t237 + Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t269 + t256 * t284 - t283 * t257;
t213 = -mrSges(6,1) * t245 + mrSges(6,3) * t242 + Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t285 - pkin(5) * t227 - t274 * t298 + t276 * t309 - t376;
t291 = Ifges(5,4) * t312 + Ifges(5,2) * t311 + Ifges(5,6) * t324;
t292 = Ifges(5,1) * t312 + Ifges(5,4) * t311 + Ifges(5,5) * t324;
t197 = mrSges(5,1) * t249 - mrSges(5,2) * t250 + Ifges(5,5) * t287 + Ifges(5,6) * t286 + Ifges(5,3) * t314 + pkin(4) * t377 + pkin(12) * t397 + t365 * t212 + t370 * t213 + t312 * t291 - t311 * t292;
t315 = Ifges(4,5) * t331 + Ifges(4,6) * t330 + Ifges(4,3) * t342;
t317 = Ifges(4,1) * t331 + Ifges(4,4) * t330 + Ifges(4,5) * t342;
t290 = Ifges(5,5) * t312 + Ifges(5,6) * t311 + Ifges(5,3) * t324;
t198 = mrSges(5,2) * t251 - mrSges(5,3) * t249 + Ifges(5,1) * t287 + Ifges(5,4) * t286 + Ifges(5,5) * t314 - pkin(12) * t220 + t212 * t370 - t213 * t365 + t290 * t311 - t291 * t324;
t375 = mrSges(6,1) * t241 - mrSges(6,2) * t242 + Ifges(6,5) * t271 + Ifges(6,6) * t270 + Ifges(6,3) * t285 + pkin(5) * t238 + pkin(13) * t228 + t369 * t229 + t364 * t230 + t298 * t275 - t297 * t276;
t199 = -mrSges(5,1) * t251 + mrSges(5,3) * t250 + Ifges(5,4) * t287 + Ifges(5,2) * t286 + Ifges(5,6) * t314 - pkin(4) * t220 - t312 * t290 + t324 * t292 - t375;
t383 = pkin(11) * t211 + t198 * t366 + t199 * t371;
t183 = -mrSges(4,1) * t293 + mrSges(4,3) * t279 + Ifges(4,4) * t321 + Ifges(4,2) * t320 + Ifges(4,6) * t339 - pkin(3) * t205 - t197 * t357 - t315 * t331 + t317 * t342 + t383 * t361;
t316 = Ifges(4,4) * t331 + Ifges(4,2) * t330 + Ifges(4,6) * t342;
t184 = mrSges(4,2) * t293 - mrSges(4,3) * t278 + Ifges(4,1) * t321 + Ifges(4,4) * t320 + Ifges(4,5) * t339 + t198 * t371 - t199 * t366 + t315 * t330 - t316 * t342 + (-t205 * t357 - t206 * t361) * pkin(11);
t384 = pkin(10) * t196 + t183 * t372 + t184 * t367;
t382 = Ifges(3,5) * t363 + (Ifges(3,1) * t356 + Ifges(3,4) * t360) * t359;
t381 = Ifges(3,6) * t363 + (Ifges(3,4) * t356 + Ifges(3,2) * t360) * t359;
t182 = mrSges(4,1) * t278 - mrSges(4,2) * t279 + Ifges(4,5) * t321 + Ifges(4,6) * t320 + Ifges(4,3) * t339 + pkin(3) * t206 + t197 * t361 + t316 * t331 - t317 * t330 + t383 * t357;
t336 = t381 * qJD(1);
t337 = t382 * qJD(1);
t173 = qJDD(1) * t417 + mrSges(3,1) * t325 - mrSges(3,2) * t326 + pkin(2) * t192 + t182 * t362 + t384 * t358 + (t395 * qJDD(1) + (t336 * t356 - t337 * t360) * qJD(1)) * t359;
t335 = (t395 * t359 + t417) * qJD(1);
t175 = -mrSges(3,1) * t332 + mrSges(3,3) * t326 - pkin(2) * t191 - t182 * t358 + (-t335 * t412 + t337 * t363) * qJD(1) + t384 * t362 + t381 * qJDD(1);
t177 = mrSges(3,2) * t332 - mrSges(3,3) * t325 - t183 * t367 + t184 * t372 + (t335 * t407 - t336 * t363) * qJD(1) + (-t191 * t358 - t192 * t362) * pkin(10) + t382 * qJDD(1);
t380 = mrSges(2,1) * t352 - mrSges(2,2) * t353 + Ifges(2,3) * qJDD(1) + pkin(1) * t181 + t363 * t173 + t175 * t407 + t177 * t412 + t187 * t415;
t185 = m(2) * t353 - mrSges(2,1) * t374 - qJDD(1) * mrSges(2,2) + t187;
t180 = t190 * t363 + (t189 * t360 + t195 * t356) * t359;
t178 = m(2) * t352 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t374 + t181;
t171 = -mrSges(2,2) * g(3) - mrSges(2,3) * t352 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t374 - t175 * t356 + t177 * t360 + (-t180 * t359 - t181 * t363) * qJ(2);
t170 = mrSges(2,1) * g(3) + mrSges(2,3) * t353 + Ifges(2,5) * t374 + Ifges(2,6) * qJDD(1) - pkin(1) * t180 - t359 * t173 + (qJ(2) * t187 + t175 * t360 + t177 * t356) * t363;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t373 * t171 - t368 * t170 - pkin(9) * (t178 * t373 + t185 * t368), t171, t177, t184, t198, t212, t230; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t368 * t171 + t373 * t170 + pkin(9) * (-t368 * t178 + t185 * t373), t170, t175, t183, t199, t213, t229; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t380, t380, t173, t182, t197, t375, t376;];
m_new  = t1;
