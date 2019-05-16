% Calculate vector of cutting torques with Newton-Euler for
% S6RPRPRR11
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
% Datum: 2019-05-05 20:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPRPRR11_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:11:34
% EndTime: 2019-05-05 20:14:26
% DurationCPUTime: 169.13s
% Computational Cost: add. (2666769->402), mult. (8462631->550), div. (0->0), fcn. (7269969->16), ass. (0->183)
t367 = cos(pkin(6));
t370 = sin(qJ(3));
t366 = cos(pkin(7));
t418 = cos(qJ(3));
t402 = t366 * t418;
t362 = sin(pkin(7));
t403 = t362 * t418;
t363 = sin(pkin(6));
t365 = cos(pkin(12));
t410 = t363 * t365;
t361 = sin(pkin(12));
t413 = t361 * t363;
t419 = -t367 * t403 + t370 * t413 - t402 * t410;
t417 = pkin(9) * t361;
t416 = Ifges(3,3) * t367;
t415 = pkin(9) * qJDD(1);
t414 = qJ(2) * t363;
t412 = t361 * t367;
t411 = t362 * t370;
t409 = t363 * t366;
t408 = t365 * t367;
t407 = t366 * t370;
t341 = (t362 * t367 + t365 * t409) * qJD(1) * pkin(9);
t371 = sin(qJ(1));
t374 = cos(qJ(1));
t357 = -t374 * g(1) - t371 * g(2);
t375 = qJD(1) ^ 2;
t345 = -t375 * pkin(1) + qJDD(1) * t414 + t357;
t391 = -pkin(2) * t365 - t362 * t417;
t406 = qJD(1) * t363;
t387 = qJD(1) * t391 * t406 + t366 * t415;
t356 = t371 * g(1) - t374 * g(2);
t344 = qJDD(1) * pkin(1) + t375 * t414 + t356;
t401 = qJD(2) * t406;
t392 = -g(3) * t410 + t344 * t408 - 0.2e1 * t361 * t401;
t293 = (pkin(2) * qJDD(1) + qJD(1) * t341) * t367 + (-t387 * t363 - t345) * t361 + t392;
t346 = (pkin(2) * t367 - t409 * t417) * qJD(1);
t404 = t344 * t412 + (t345 + 0.2e1 * t401) * t365;
t294 = (-qJD(1) * t346 + t362 * t415) * t367 + (-g(3) * t361 + t387 * t365) * t363 + t404;
t400 = -t367 * g(3) + qJDD(2);
t302 = (-t344 + t391 * qJDD(1) + (-t341 * t365 + t346 * t361) * qJD(1)) * t363 + t400;
t264 = t293 * t407 + t418 * t294 + t302 * t411;
t329 = t419 * qJD(1);
t380 = t367 * t411 + (t418 * t361 + t365 * t407) * t363;
t330 = t380 * qJD(1);
t314 = t329 * pkin(3) - t330 * qJ(4);
t388 = -t362 * t410 + t366 * t367;
t342 = t388 * qJD(1) + qJD(3);
t338 = t342 ^ 2;
t339 = t388 * qJDD(1) + qJDD(3);
t254 = -t338 * pkin(3) + t339 * qJ(4) - t329 * t314 + t264;
t280 = -t362 * t293 + t366 * t302;
t316 = t330 * qJD(3) + t419 * qJDD(1);
t317 = -t329 * qJD(3) + t380 * qJDD(1);
t257 = (t329 * t342 - t317) * qJ(4) + (t330 * t342 + t316) * pkin(3) + t280;
t360 = sin(pkin(13));
t364 = cos(pkin(13));
t324 = t364 * t330 + t360 * t342;
t246 = -0.2e1 * qJD(4) * t324 - t360 * t254 + t364 * t257;
t308 = t364 * t317 + t360 * t339;
t323 = -t360 * t330 + t364 * t342;
t243 = (t323 * t329 - t308) * pkin(10) + (t323 * t324 + t316) * pkin(4) + t246;
t247 = 0.2e1 * qJD(4) * t323 + t364 * t254 + t360 * t257;
t306 = t329 * pkin(4) - t324 * pkin(10);
t307 = -t360 * t317 + t364 * t339;
t320 = t323 ^ 2;
t245 = -t320 * pkin(4) + t307 * pkin(10) - t329 * t306 + t247;
t369 = sin(qJ(5));
t373 = cos(qJ(5));
t240 = t369 * t243 + t373 * t245;
t297 = t369 * t323 + t373 * t324;
t270 = -t297 * qJD(5) + t373 * t307 - t369 * t308;
t296 = t373 * t323 - t369 * t324;
t278 = -t296 * mrSges(6,1) + t297 * mrSges(6,2);
t328 = qJD(5) + t329;
t284 = t328 * mrSges(6,1) - t297 * mrSges(6,3);
t313 = qJDD(5) + t316;
t279 = -t296 * pkin(5) - t297 * pkin(11);
t327 = t328 ^ 2;
t237 = -t327 * pkin(5) + t313 * pkin(11) + t296 * t279 + t240;
t263 = t293 * t402 - t370 * t294 + t302 * t403;
t253 = -t339 * pkin(3) - t338 * qJ(4) + t330 * t314 + qJDD(4) - t263;
t248 = -t307 * pkin(4) - t320 * pkin(10) + t324 * t306 + t253;
t271 = t296 * qJD(5) + t369 * t307 + t373 * t308;
t241 = (-t296 * t328 - t271) * pkin(11) + (t297 * t328 - t270) * pkin(5) + t248;
t368 = sin(qJ(6));
t372 = cos(qJ(6));
t234 = -t368 * t237 + t372 * t241;
t281 = -t368 * t297 + t372 * t328;
t251 = t281 * qJD(6) + t372 * t271 + t368 * t313;
t282 = t372 * t297 + t368 * t328;
t265 = -t281 * mrSges(7,1) + t282 * mrSges(7,2);
t269 = qJDD(6) - t270;
t295 = qJD(6) - t296;
t272 = -t295 * mrSges(7,2) + t281 * mrSges(7,3);
t230 = m(7) * t234 + t269 * mrSges(7,1) - t251 * mrSges(7,3) - t282 * t265 + t295 * t272;
t235 = t372 * t237 + t368 * t241;
t250 = -t282 * qJD(6) - t368 * t271 + t372 * t313;
t273 = t295 * mrSges(7,1) - t282 * mrSges(7,3);
t231 = m(7) * t235 - t269 * mrSges(7,2) + t250 * mrSges(7,3) + t281 * t265 - t295 * t273;
t397 = -t368 * t230 + t372 * t231;
t214 = m(6) * t240 - t313 * mrSges(6,2) + t270 * mrSges(6,3) + t296 * t278 - t328 * t284 + t397;
t239 = t373 * t243 - t369 * t245;
t283 = -t328 * mrSges(6,2) + t296 * mrSges(6,3);
t236 = -t313 * pkin(5) - t327 * pkin(11) + t297 * t279 - t239;
t385 = -m(7) * t236 + t250 * mrSges(7,1) - t251 * mrSges(7,2) + t281 * t272 - t282 * t273;
t226 = m(6) * t239 + t313 * mrSges(6,1) - t271 * mrSges(6,3) - t297 * t278 + t328 * t283 + t385;
t211 = t369 * t214 + t373 * t226;
t298 = -t323 * mrSges(5,1) + t324 * mrSges(5,2);
t304 = -t329 * mrSges(5,2) + t323 * mrSges(5,3);
t209 = m(5) * t246 + t316 * mrSges(5,1) - t308 * mrSges(5,3) - t324 * t298 + t329 * t304 + t211;
t305 = t329 * mrSges(5,1) - t324 * mrSges(5,3);
t398 = t373 * t214 - t369 * t226;
t210 = m(5) * t247 - t316 * mrSges(5,2) + t307 * mrSges(5,3) + t323 * t298 - t329 * t305 + t398;
t203 = t364 * t209 + t360 * t210;
t219 = t372 * t230 + t368 * t231;
t315 = t329 * mrSges(4,1) + t330 * mrSges(4,2);
t326 = t342 * mrSges(4,1) - t330 * mrSges(4,3);
t399 = -t360 * t209 + t364 * t210;
t200 = m(4) * t264 - t339 * mrSges(4,2) - t316 * mrSges(4,3) - t329 * t315 - t342 * t326 + t399;
t325 = -t342 * mrSges(4,2) - t329 * mrSges(4,3);
t202 = m(4) * t280 + t316 * mrSges(4,1) + t317 * mrSges(4,2) + t329 * t325 + t330 * t326 + t203;
t381 = m(6) * t248 - t270 * mrSges(6,1) + t271 * mrSges(6,2) - t296 * t283 + t297 * t284 + t219;
t377 = -m(5) * t253 + t307 * mrSges(5,1) - t308 * mrSges(5,2) + t323 * t304 - t324 * t305 - t381;
t217 = m(4) * t263 + t339 * mrSges(4,1) - t317 * mrSges(4,3) - t330 * t315 + t342 * t325 + t377;
t189 = t200 * t411 + t366 * t202 + t217 * t403;
t190 = t200 * t407 - t362 * t202 + t217 * t402;
t321 = -t361 * t345 + t392;
t395 = -mrSges(3,1) * t365 + mrSges(3,2) * t361;
t343 = t395 * t406;
t389 = -mrSges(3,2) * t367 + mrSges(3,3) * t410;
t348 = t389 * qJD(1);
t390 = mrSges(3,1) * t367 - mrSges(3,3) * t413;
t187 = m(3) * t321 + t390 * qJDD(1) + (-t343 * t413 + t348 * t367) * qJD(1) + t190;
t196 = t418 * t200 - t370 * t217;
t322 = -g(3) * t413 + t404;
t347 = t390 * qJD(1);
t195 = m(3) * t322 + t389 * qJDD(1) + (t343 * t410 - t347 * t367) * qJD(1) + t196;
t184 = -t361 * t187 + t365 * t195;
t331 = -t363 * t344 + t400;
t188 = m(3) * t331 + (t395 * qJDD(1) + (t347 * t361 - t348 * t365) * qJD(1)) * t363 + t189;
t179 = t187 * t408 - t363 * t188 + t195 * t412;
t394 = Ifges(3,5) * t361 + Ifges(3,6) * t365;
t258 = Ifges(7,5) * t282 + Ifges(7,6) * t281 + Ifges(7,3) * t295;
t260 = Ifges(7,1) * t282 + Ifges(7,4) * t281 + Ifges(7,5) * t295;
t223 = -mrSges(7,1) * t236 + mrSges(7,3) * t235 + Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t269 - t282 * t258 + t295 * t260;
t259 = Ifges(7,4) * t282 + Ifges(7,2) * t281 + Ifges(7,6) * t295;
t224 = mrSges(7,2) * t236 - mrSges(7,3) * t234 + Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t269 + t281 * t258 - t295 * t259;
t274 = Ifges(6,5) * t297 + Ifges(6,6) * t296 + Ifges(6,3) * t328;
t275 = Ifges(6,4) * t297 + Ifges(6,2) * t296 + Ifges(6,6) * t328;
t204 = mrSges(6,2) * t248 - mrSges(6,3) * t239 + Ifges(6,1) * t271 + Ifges(6,4) * t270 + Ifges(6,5) * t313 - pkin(11) * t219 - t368 * t223 + t372 * t224 + t296 * t274 - t328 * t275;
t276 = Ifges(6,1) * t297 + Ifges(6,4) * t296 + Ifges(6,5) * t328;
t378 = mrSges(7,1) * t234 - mrSges(7,2) * t235 + Ifges(7,5) * t251 + Ifges(7,6) * t250 + Ifges(7,3) * t269 + t282 * t259 - t281 * t260;
t205 = -mrSges(6,1) * t248 + mrSges(6,3) * t240 + Ifges(6,4) * t271 + Ifges(6,2) * t270 + Ifges(6,6) * t313 - pkin(5) * t219 - t297 * t274 + t328 * t276 - t378;
t285 = Ifges(5,5) * t324 + Ifges(5,6) * t323 + Ifges(5,3) * t329;
t287 = Ifges(5,1) * t324 + Ifges(5,4) * t323 + Ifges(5,5) * t329;
t191 = -mrSges(5,1) * t253 + mrSges(5,3) * t247 + Ifges(5,4) * t308 + Ifges(5,2) * t307 + Ifges(5,6) * t316 - pkin(4) * t381 + pkin(10) * t398 + t369 * t204 + t373 * t205 - t324 * t285 + t329 * t287;
t286 = Ifges(5,4) * t324 + Ifges(5,2) * t323 + Ifges(5,6) * t329;
t192 = mrSges(5,2) * t253 - mrSges(5,3) * t246 + Ifges(5,1) * t308 + Ifges(5,4) * t307 + Ifges(5,5) * t316 - pkin(10) * t211 + t373 * t204 - t369 * t205 + t323 * t285 - t329 * t286;
t309 = Ifges(4,5) * t330 - Ifges(4,6) * t329 + Ifges(4,3) * t342;
t310 = Ifges(4,4) * t330 - Ifges(4,2) * t329 + Ifges(4,6) * t342;
t181 = mrSges(4,2) * t280 - mrSges(4,3) * t263 + Ifges(4,1) * t317 - Ifges(4,4) * t316 + Ifges(4,5) * t339 - qJ(4) * t203 - t360 * t191 + t364 * t192 - t329 * t309 - t342 * t310;
t311 = Ifges(4,1) * t330 - Ifges(4,4) * t329 + Ifges(4,5) * t342;
t379 = -mrSges(6,1) * t239 + mrSges(6,2) * t240 - Ifges(6,5) * t271 - Ifges(6,6) * t270 - Ifges(6,3) * t313 - pkin(5) * t385 - pkin(11) * t397 - t372 * t223 - t368 * t224 - t297 * t275 + t296 * t276;
t376 = -mrSges(5,1) * t246 + mrSges(5,2) * t247 - Ifges(5,5) * t308 - Ifges(5,6) * t307 - pkin(4) * t211 - t324 * t286 + t323 * t287 + t379;
t185 = Ifges(4,6) * t339 + t342 * t311 - t330 * t309 + Ifges(4,4) * t317 - mrSges(4,1) * t280 + mrSges(4,3) * t264 - pkin(3) * t203 + t376 + (-Ifges(5,3) - Ifges(4,2)) * t316;
t386 = pkin(9) * t196 + t181 * t370 + t418 * t185;
t384 = Ifges(3,5) * t367 + (Ifges(3,1) * t361 + Ifges(3,4) * t365) * t363;
t383 = Ifges(3,6) * t367 + (Ifges(3,4) * t361 + Ifges(3,2) * t365) * t363;
t180 = mrSges(4,1) * t263 - mrSges(4,2) * t264 + Ifges(4,5) * t317 - Ifges(4,6) * t316 + Ifges(4,3) * t339 + pkin(3) * t377 + qJ(4) * t399 + t364 * t191 + t360 * t192 + t330 * t310 + t329 * t311;
t335 = t383 * qJD(1);
t336 = t384 * qJD(1);
t171 = qJDD(1) * t416 + mrSges(3,1) * t321 - mrSges(3,2) * t322 + pkin(2) * t190 + t366 * t180 + t386 * t362 + (t394 * qJDD(1) + (t335 * t361 - t336 * t365) * qJD(1)) * t363;
t334 = (t394 * t363 + t416) * qJD(1);
t173 = -mrSges(3,1) * t331 + mrSges(3,3) * t322 - pkin(2) * t189 - t362 * t180 + (-t334 * t413 + t336 * t367) * qJD(1) + t386 * t366 + t383 * qJDD(1);
t175 = mrSges(3,2) * t331 - mrSges(3,3) * t321 + t418 * t181 - t370 * t185 + (t334 * t410 - t335 * t367) * qJD(1) + (-t189 * t362 - t190 * t366) * pkin(9) + t384 * qJDD(1);
t382 = mrSges(2,1) * t356 - mrSges(2,2) * t357 + Ifges(2,3) * qJDD(1) + pkin(1) * t179 + t367 * t171 + t173 * t410 + t175 * t413 + t184 * t414;
t182 = m(2) * t357 - t375 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t184;
t178 = t367 * t188 + (t187 * t365 + t195 * t361) * t363;
t176 = m(2) * t356 + qJDD(1) * mrSges(2,1) - t375 * mrSges(2,2) + t179;
t169 = -mrSges(2,2) * g(3) - mrSges(2,3) * t356 + Ifges(2,5) * qJDD(1) - t375 * Ifges(2,6) - t361 * t173 + t365 * t175 + (-t178 * t363 - t179 * t367) * qJ(2);
t168 = mrSges(2,1) * g(3) + mrSges(2,3) * t357 + t375 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t178 - t363 * t171 + (qJ(2) * t184 + t173 * t365 + t175 * t361) * t367;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t374 * t169 - t371 * t168 - pkin(8) * (t374 * t176 + t371 * t182), t169, t175, t181, t192, t204, t224; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t371 * t169 + t374 * t168 + pkin(8) * (-t371 * t176 + t374 * t182), t168, t173, t185, t191, t205, t223; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t382, t382, t171, t180, Ifges(5,3) * t316 - t376, -t379, t378;];
m_new  = t1;
