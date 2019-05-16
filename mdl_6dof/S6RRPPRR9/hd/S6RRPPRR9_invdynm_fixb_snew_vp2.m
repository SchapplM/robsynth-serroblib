% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-05-06 11:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR9_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:31:47
% EndTime: 2019-05-06 11:32:09
% DurationCPUTime: 9.75s
% Computational Cost: add. (153417->413), mult. (358351->499), div. (0->0), fcn. (239539->10), ass. (0->170)
t455 = -2 * qJD(3);
t454 = -2 * qJD(4);
t389 = sin(qJ(1));
t393 = cos(qJ(1));
t370 = t389 * g(1) - g(2) * t393;
t394 = qJD(1) ^ 2;
t384 = sin(pkin(6));
t452 = pkin(8) * t384;
t342 = qJDD(1) * pkin(1) + t394 * t452 + t370;
t371 = -g(1) * t393 - g(2) * t389;
t425 = qJDD(1) * t384;
t343 = -pkin(1) * t394 + pkin(8) * t425 + t371;
t388 = sin(qJ(2));
t385 = cos(pkin(6));
t392 = cos(qJ(2));
t437 = t385 * t392;
t439 = t384 * t392;
t281 = -g(3) * t439 + t342 * t437 - t388 * t343;
t429 = qJD(1) * t384;
t344 = (-pkin(2) * t392 - qJ(3) * t388) * t429;
t377 = qJD(1) * t385 + qJD(2);
t375 = t377 ^ 2;
t376 = qJDD(1) * t385 + qJDD(2);
t421 = t388 * t429;
t263 = -t376 * pkin(2) - t375 * qJ(3) + t344 * t421 + qJDD(3) - t281;
t428 = qJD(1) * t392;
t349 = (qJD(2) * t428 + qJDD(1) * t388) * t384;
t420 = t384 * t428;
t417 = t377 * t420;
t441 = t384 ^ 2 * t394;
t422 = t392 * t441;
t249 = t263 - (t388 * t422 + t376) * qJ(4) - (-t349 + t417) * pkin(3) + t377 * t454;
t438 = t385 * t388;
t432 = t342 * t438 + t392 * t343;
t453 = pkin(2) * t375 - t376 * qJ(3) - t344 * t420 + t377 * t455 - t432;
t451 = g(3) * t385;
t450 = mrSges(3,1) - mrSges(4,2);
t449 = -mrSges(3,3) - mrSges(5,1);
t448 = Ifges(4,6) + Ifges(3,4);
t447 = mrSges(5,2) * t349;
t446 = t349 * mrSges(5,1);
t445 = qJ(3) * t377;
t337 = pkin(3) * t421 - qJ(4) * t377;
t444 = t337 * t388;
t340 = mrSges(5,1) * t420 + mrSges(5,2) * t377;
t443 = t340 * t392;
t442 = t342 * t384;
t440 = t384 * t388;
t296 = Ifges(5,5) * t377 + (-Ifges(5,6) * t392 + Ifges(5,3) * t388) * t429;
t436 = t392 * t296;
t348 = pkin(4) * t420 - pkin(9) * t377;
t382 = t388 ^ 2;
t383 = t392 ^ 2;
t350 = -qJD(2) * t421 + t392 * t425;
t413 = t421 * t455 - t451 + (t377 * t421 - t350) * pkin(2);
t411 = -t350 * qJ(4) + t420 * t454 + t413;
t239 = (-pkin(3) * t383 - pkin(4) * t382) * t441 + (pkin(9) - qJ(3)) * t349 + (-t342 + (-t444 + (-t348 - t445) * t392) * qJD(1)) * t384 + t411;
t423 = t383 * t441;
t401 = -qJ(4) * t423 + t377 * t337 + qJDD(4) - t453;
t243 = -pkin(9) * t376 + (pkin(3) + pkin(4)) * t350 + (pkin(9) * t422 + (pkin(4) * qJD(1) * t377 - g(3)) * t384) * t388 + t401;
t387 = sin(qJ(5));
t391 = cos(qJ(5));
t236 = t391 * t239 + t387 * t243;
t317 = t377 * t391 + t387 * t421;
t279 = -qJD(5) * t317 + t349 * t391 - t376 * t387;
t316 = -t377 * t387 + t391 * t421;
t283 = -mrSges(6,1) * t316 + mrSges(6,2) * t317;
t358 = qJD(5) + t420;
t288 = mrSges(6,1) * t358 - mrSges(6,3) * t317;
t333 = qJDD(5) + t350;
t284 = -pkin(5) * t316 - pkin(10) * t317;
t355 = t358 ^ 2;
t232 = -pkin(5) * t355 + pkin(10) * t333 + t284 * t316 + t236;
t242 = -pkin(9) * t382 * t441 - pkin(4) * t349 + t377 * t348 - t249;
t280 = qJD(5) * t316 + t349 * t387 + t376 * t391;
t237 = t242 + (t317 * t358 - t279) * pkin(5) + (-t316 * t358 - t280) * pkin(10);
t386 = sin(qJ(6));
t390 = cos(qJ(6));
t229 = -t232 * t386 + t237 * t390;
t285 = -t317 * t386 + t358 * t390;
t254 = qJD(6) * t285 + t280 * t390 + t333 * t386;
t286 = t317 * t390 + t358 * t386;
t264 = -mrSges(7,1) * t285 + mrSges(7,2) * t286;
t315 = qJD(6) - t316;
t267 = -mrSges(7,2) * t315 + mrSges(7,3) * t285;
t276 = qJDD(6) - t279;
t225 = m(7) * t229 + mrSges(7,1) * t276 - mrSges(7,3) * t254 - t264 * t286 + t267 * t315;
t230 = t232 * t390 + t237 * t386;
t253 = -qJD(6) * t286 - t280 * t386 + t333 * t390;
t268 = mrSges(7,1) * t315 - mrSges(7,3) * t286;
t226 = m(7) * t230 - mrSges(7,2) * t276 + mrSges(7,3) * t253 + t264 * t285 - t268 * t315;
t418 = -t225 * t386 + t390 * t226;
t210 = m(6) * t236 - mrSges(6,2) * t333 + mrSges(6,3) * t279 + t283 * t316 - t288 * t358 + t418;
t235 = -t239 * t387 + t243 * t391;
t287 = -mrSges(6,2) * t358 + mrSges(6,3) * t316;
t231 = -pkin(5) * t333 - pkin(10) * t355 + t284 * t317 - t235;
t409 = -m(7) * t231 + t253 * mrSges(7,1) - mrSges(7,2) * t254 + t285 * t267 - t268 * t286;
t221 = m(6) * t235 + mrSges(6,1) * t333 - mrSges(6,3) * t280 - t283 * t317 + t287 * t358 + t409;
t203 = t387 * t210 + t391 * t221;
t214 = t390 * t225 + t386 * t226;
t301 = Ifges(4,1) * t377 + (-Ifges(4,4) * t388 - Ifges(4,5) * t392) * t429;
t435 = Ifges(3,3) * t377 + (Ifges(3,5) * t388 + Ifges(3,6) * t392) * t429 + t301;
t299 = Ifges(4,4) * t377 + (-Ifges(4,2) * t388 - Ifges(4,6) * t392) * t429;
t434 = -t296 + t299;
t297 = Ifges(4,5) * t377 + (-Ifges(4,6) * t388 - Ifges(4,3) * t392) * t429;
t298 = Ifges(5,4) * t377 + (-Ifges(5,2) * t392 + Ifges(5,6) * t388) * t429;
t433 = -t297 - t298;
t338 = mrSges(5,1) * t421 - mrSges(5,3) * t377;
t341 = mrSges(4,1) * t421 + mrSges(4,2) * t377;
t431 = -t338 - t341;
t345 = (mrSges(4,2) * t392 - mrSges(4,3) * t388) * t429;
t347 = (-mrSges(5,2) * t388 - mrSges(5,3) * t392) * t429;
t430 = -t345 - t347;
t424 = g(3) * t440;
t282 = -t424 + t432;
t335 = mrSges(3,1) * t377 - mrSges(3,3) * t421;
t346 = (-mrSges(3,1) * t392 + mrSges(3,2) * t388) * t429;
t256 = t424 + t453;
t251 = pkin(3) * t350 + t401 - t424;
t415 = m(5) * t251 + t376 * mrSges(5,2) + t377 * t338 + t347 * t420 + t203;
t405 = -m(4) * t256 + t376 * mrSges(4,3) + t377 * t341 + t345 * t420 + t415;
t198 = m(3) * t282 - mrSges(3,2) * t376 - t335 * t377 + (mrSges(4,1) - t449) * t350 + t346 * t420 + t405;
t336 = -mrSges(3,2) * t377 + mrSges(3,3) * t420;
t339 = -mrSges(4,1) * t420 - mrSges(4,3) * t377;
t412 = -m(6) * t242 + t279 * mrSges(6,1) - t280 * mrSges(6,2) + t316 * t287 - t317 * t288 - t214;
t404 = -m(5) * t249 + t376 * mrSges(5,3) + t377 * t340 - t412;
t400 = -m(4) * t263 - t349 * mrSges(4,1) + t404;
t206 = (-t346 + t430) * t421 + t400 + (t336 - t339) * t377 + t450 * t376 + t449 * t349 + m(3) * t281;
t189 = t392 * t198 - t206 * t388;
t419 = t391 * t210 - t221 * t387;
t302 = -t442 - t451;
t257 = -t442 + (-t349 - t417) * qJ(3) + t413;
t246 = -pkin(3) * t423 - qJ(3) * t349 + (-t342 + (-t392 * t445 - t444) * qJD(1)) * t384 + t411;
t416 = m(5) * t246 - t350 * mrSges(5,3) + t419;
t410 = m(4) * t257 - t349 * mrSges(4,3) + t339 * t420 + t416;
t195 = m(3) * t302 - t450 * t350 + (mrSges(3,2) - mrSges(5,2)) * t349 + ((-t336 - t340) * t392 + (t335 + t431) * t388) * t429 + t410;
t186 = -t195 * t384 + t198 * t438 + t206 * t437;
t294 = Ifges(3,6) * t377 + (Ifges(3,4) * t388 + Ifges(3,2) * t392) * t429;
t295 = Ifges(3,5) * t377 + (Ifges(3,1) * t388 + Ifges(3,4) * t392) * t429;
t207 = t347 * t421 - t404 + t446;
t258 = Ifges(7,5) * t286 + Ifges(7,6) * t285 + Ifges(7,3) * t315;
t260 = Ifges(7,1) * t286 + Ifges(7,4) * t285 + Ifges(7,5) * t315;
t218 = -mrSges(7,1) * t231 + mrSges(7,3) * t230 + Ifges(7,4) * t254 + Ifges(7,2) * t253 + Ifges(7,6) * t276 - t258 * t286 + t260 * t315;
t259 = Ifges(7,4) * t286 + Ifges(7,2) * t285 + Ifges(7,6) * t315;
t219 = mrSges(7,2) * t231 - mrSges(7,3) * t229 + Ifges(7,1) * t254 + Ifges(7,4) * t253 + Ifges(7,5) * t276 + t258 * t285 - t259 * t315;
t269 = Ifges(6,5) * t317 + Ifges(6,6) * t316 + Ifges(6,3) * t358;
t270 = Ifges(6,4) * t317 + Ifges(6,2) * t316 + Ifges(6,6) * t358;
t192 = mrSges(6,2) * t242 - mrSges(6,3) * t235 + Ifges(6,1) * t280 + Ifges(6,4) * t279 + Ifges(6,5) * t333 - pkin(10) * t214 - t218 * t386 + t219 * t390 + t269 * t316 - t270 * t358;
t271 = Ifges(6,1) * t317 + Ifges(6,4) * t316 + Ifges(6,5) * t358;
t399 = mrSges(7,1) * t229 - mrSges(7,2) * t230 + Ifges(7,5) * t254 + Ifges(7,6) * t253 + Ifges(7,3) * t276 + t259 * t286 - t260 * t285;
t194 = -mrSges(6,1) * t242 + mrSges(6,3) * t236 + Ifges(6,4) * t280 + Ifges(6,2) * t279 + Ifges(6,6) * t333 - pkin(5) * t214 - t269 * t317 + t271 * t358 - t399;
t403 = mrSges(5,2) * t251 - mrSges(5,3) * t249 + Ifges(5,1) * t376 - Ifges(5,4) * t350 + Ifges(5,5) * t349 - pkin(9) * t203 + t391 * t192 - t387 * t194;
t396 = mrSges(4,2) * t263 - mrSges(4,3) * t256 + Ifges(4,1) * t376 - Ifges(4,4) * t349 - Ifges(4,5) * t350 - qJ(4) * t207 + t299 * t420 + t403;
t178 = t396 + qJ(3) * t405 + (qJ(3) * (mrSges(4,1) + mrSges(5,1)) + Ifges(3,6)) * t350 + ((-t295 - t296) * t392 + (pkin(2) * t430 + t294 + t433) * t388) * t429 + pkin(2) * (-t376 * mrSges(4,2) - t377 * t339 + t400 - t446) + Ifges(3,3) * t376 + Ifges(3,5) * t349 + mrSges(3,1) * t281 - mrSges(3,2) * t282;
t199 = mrSges(4,2) * t350 - t447 + (t388 * t431 - t443) * t429 + t410;
t300 = Ifges(5,1) * t377 + (-Ifges(5,4) * t392 + Ifges(5,5) * t388) * t429;
t402 = mrSges(5,1) * t249 - mrSges(5,2) * t246 + Ifges(5,5) * t376 - Ifges(5,6) * t350 + Ifges(5,3) * t349 + pkin(4) * t412 + pkin(9) * t419 + t387 * t192 + t391 * t194 + t377 * t298 + t300 * t420;
t398 = mrSges(4,1) * t263 - mrSges(4,3) * t257 + pkin(3) * t207 + t402;
t180 = t398 + t435 * t420 + (-t294 + t297) * t377 + (Ifges(3,5) - Ifges(4,4)) * t376 + t448 * t350 + (Ifges(3,1) + Ifges(4,2)) * t349 + mrSges(3,2) * t302 - mrSges(3,3) * t281 - qJ(3) * t199;
t407 = mrSges(6,1) * t235 - mrSges(6,2) * t236 + Ifges(6,5) * t280 + Ifges(6,6) * t279 + Ifges(6,3) * t333 + pkin(5) * t409 + pkin(10) * t418 + t390 * t218 + t386 * t219 + t317 * t270 - t316 * t271;
t397 = mrSges(5,1) * t251 - mrSges(5,3) * t246 - Ifges(5,4) * t376 + Ifges(5,2) * t350 - Ifges(5,6) * t349 + pkin(4) * t203 - t300 * t421 + t407;
t395 = mrSges(4,1) * t256 - mrSges(4,2) * t257 + pkin(3) * (-mrSges(5,1) * t350 - t415) + qJ(4) * (-t447 + (-t338 * t388 - t443) * t429 + t416) - t397;
t182 = (t295 - t434) * t377 + (Ifges(3,6) - Ifges(4,5)) * t376 + (Ifges(4,3) + Ifges(3,2)) * t350 + t448 * t349 - t395 - t435 * t421 - mrSges(3,1) * t302 + mrSges(3,3) * t282 - pkin(2) * t199;
t408 = mrSges(2,1) * t370 - mrSges(2,2) * t371 + Ifges(2,3) * qJDD(1) + pkin(1) * t186 + t385 * t178 + t180 * t440 + t182 * t439 + t189 * t452;
t187 = m(2) * t371 - mrSges(2,1) * t394 - qJDD(1) * mrSges(2,2) + t189;
t185 = t195 * t385 + (t198 * t388 + t206 * t392) * t384;
t183 = m(2) * t370 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t394 + t186;
t176 = -mrSges(2,2) * g(3) - mrSges(2,3) * t370 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t394 + t180 * t392 - t182 * t388 + (-t185 * t384 - t186 * t385) * pkin(8);
t175 = mrSges(2,1) * g(3) + mrSges(2,3) * t371 + Ifges(2,5) * t394 + Ifges(2,6) * qJDD(1) - pkin(1) * t185 - t178 * t384 + (pkin(8) * t189 + t180 * t388 + t182 * t392) * t385;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t393 * t176 - t389 * t175 - pkin(7) * (t183 * t393 + t187 * t389), t176, t180, t396 + (t388 * t433 - t436) * t429, (-t388 * t298 - t436) * t429 + t403, t192, t219; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t389 * t176 + t393 * t175 + pkin(7) * (-t183 * t389 + t187 * t393), t175, t182, Ifges(4,4) * t376 - Ifges(4,2) * t349 - Ifges(4,6) * t350 - t377 * t297 - t301 * t420 - t398, -t377 * t296 - t397, t194, t218; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t408, t408, t178, Ifges(4,5) * t376 - Ifges(4,6) * t349 - Ifges(4,3) * t350 + t301 * t421 + t377 * t434 + t395, t402, t407, t399;];
m_new  = t1;
