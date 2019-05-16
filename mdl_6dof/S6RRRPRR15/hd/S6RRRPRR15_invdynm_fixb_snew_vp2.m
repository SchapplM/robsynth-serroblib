% Calculate vector of cutting torques with Newton-Euler for
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-05-07 17:28
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRRPRR15_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR15_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:05:57
% EndTime: 2019-05-07 17:08:01
% DurationCPUTime: 72.51s
% Computational Cost: add. (1248661->422), mult. (3126781->543), div. (0->0), fcn. (2597171->14), ass. (0->177)
t373 = cos(qJ(2));
t364 = sin(pkin(7));
t418 = cos(pkin(6));
t391 = qJD(1) * t418 + qJD(2);
t387 = t391 * t364;
t365 = sin(pkin(6));
t408 = qJD(1) * t365;
t417 = cos(pkin(7));
t392 = t417 * t408;
t424 = t373 * t392 + t387;
t340 = t424 * pkin(10);
t369 = sin(qJ(2));
t421 = pkin(10) * t369;
t346 = (-pkin(2) * t373 - t364 * t421) * t408;
t404 = qJD(1) * qJD(2);
t352 = (qJDD(1) * t369 + t373 * t404) * t365;
t361 = qJDD(1) * t418 + qJDD(2);
t370 = sin(qJ(1));
t374 = cos(qJ(1));
t359 = t370 * g(1) - g(2) * t374;
t375 = qJD(1) ^ 2;
t422 = pkin(9) * t365;
t349 = qJDD(1) * pkin(1) + t375 * t422 + t359;
t360 = -g(1) * t374 - g(2) * t370;
t350 = -pkin(1) * t375 + qJDD(1) * t422 + t360;
t396 = t373 * t418;
t394 = t349 * t396 - t369 * t350;
t401 = t417 * pkin(10);
t407 = qJD(1) * t369;
t276 = -t352 * t401 + t361 * pkin(2) + t391 * t340 + (-g(3) * t373 - t346 * t407) * t365 + t394;
t345 = pkin(2) * t391 - t392 * t421;
t353 = (qJDD(1) * t373 - t369 * t404) * t365;
t389 = t353 * t417 + t361 * t364;
t406 = qJD(1) * t373;
t397 = t369 * t418;
t409 = t349 * t397 + t373 * t350;
t277 = -t391 * t345 + (-g(3) * t369 + t346 * t406) * t365 + t389 * pkin(10) + t409;
t402 = t418 * g(3);
t283 = -t352 * t364 * pkin(10) - t402 - t353 * pkin(2) + (-t349 + (-t340 * t373 + t345 * t369) * qJD(1)) * t365;
t368 = sin(qJ(3));
t398 = t368 * t417;
t415 = t364 * t368;
t423 = cos(qJ(3));
t253 = t276 * t398 + t277 * t423 + t283 * t415;
t400 = t365 * t407;
t326 = t368 * t400 - t424 * t423;
t327 = t368 * t387 + (t369 * t423 + t373 * t398) * t408;
t308 = pkin(3) * t326 - qJ(4) * t327;
t331 = -t364 * t353 + t361 * t417 + qJDD(3);
t399 = t365 * t406;
t342 = t364 * t399 - t391 * t417 - qJD(3);
t339 = t342 ^ 2;
t246 = t339 * pkin(3) - t331 * qJ(4) + 0.2e1 * qJD(4) * t342 + t326 * t308 - t253;
t420 = mrSges(4,1) - mrSges(5,2);
t419 = -Ifges(5,6) - Ifges(4,4);
t416 = t326 * t342;
t414 = t365 * t369;
t413 = t365 * t373;
t393 = t417 * t423;
t403 = t364 * t423;
t252 = t276 * t393 - t368 * t277 + t283 * t403;
t248 = -t331 * pkin(3) - t339 * qJ(4) + t327 * t308 + qJDD(4) - t252;
t305 = -t326 * qJD(3) + t352 * t423 + t389 * t368;
t240 = (t326 * t327 - t331) * pkin(11) + (t305 - t416) * pkin(4) + t248;
t304 = t327 * qJD(3) + t368 * t352 - t353 * t393 - t361 * t403;
t318 = pkin(4) * t327 + pkin(11) * t342;
t325 = t326 ^ 2;
t258 = -t364 * t276 + t417 * t283;
t382 = (-t305 - t416) * qJ(4) + t258 + (-t342 * pkin(3) - 0.2e1 * qJD(4)) * t327;
t241 = -t325 * pkin(4) - t327 * t318 + (pkin(3) + pkin(11)) * t304 + t382;
t367 = sin(qJ(5));
t372 = cos(qJ(5));
t236 = t367 * t240 + t372 * t241;
t312 = t326 * t372 + t342 * t367;
t313 = t326 * t367 - t342 * t372;
t280 = -pkin(5) * t312 - pkin(12) * t313;
t303 = qJDD(5) + t305;
t324 = qJD(5) + t327;
t323 = t324 ^ 2;
t233 = -pkin(5) * t323 + pkin(12) * t303 + t280 * t312 + t236;
t243 = -t304 * pkin(4) - t325 * pkin(11) - t342 * t318 - t246;
t265 = -qJD(5) * t313 + t304 * t372 - t331 * t367;
t266 = qJD(5) * t312 + t304 * t367 + t331 * t372;
t237 = (-t312 * t324 - t266) * pkin(12) + (t313 * t324 - t265) * pkin(5) + t243;
t366 = sin(qJ(6));
t371 = cos(qJ(6));
t228 = -t233 * t366 + t237 * t371;
t287 = -t313 * t366 + t324 * t371;
t251 = qJD(6) * t287 + t266 * t371 + t303 * t366;
t288 = t313 * t371 + t324 * t366;
t260 = -mrSges(7,1) * t287 + mrSges(7,2) * t288;
t264 = qJDD(6) - t265;
t311 = qJD(6) - t312;
t267 = -mrSges(7,2) * t311 + mrSges(7,3) * t287;
t226 = m(7) * t228 + mrSges(7,1) * t264 - mrSges(7,3) * t251 - t260 * t288 + t267 * t311;
t229 = t233 * t371 + t237 * t366;
t250 = -qJD(6) * t288 - t266 * t366 + t303 * t371;
t268 = mrSges(7,1) * t311 - mrSges(7,3) * t288;
t227 = m(7) * t229 - mrSges(7,2) * t264 + mrSges(7,3) * t250 + t260 * t287 - t268 * t311;
t216 = t371 * t226 + t366 * t227;
t295 = -Ifges(5,1) * t342 - Ifges(5,4) * t327 + Ifges(5,5) * t326;
t412 = -Ifges(4,5) * t327 + Ifges(4,6) * t326 + Ifges(4,3) * t342 - t295;
t293 = -Ifges(5,4) * t342 - Ifges(5,2) * t327 + Ifges(5,6) * t326;
t411 = -Ifges(4,1) * t327 + Ifges(4,4) * t326 + Ifges(4,5) * t342 + t293;
t315 = mrSges(5,1) * t326 + mrSges(5,3) * t342;
t410 = mrSges(4,2) * t342 - mrSges(4,3) * t326 - t315;
t317 = -mrSges(4,1) * t342 - mrSges(4,3) * t327;
t279 = -mrSges(6,1) * t312 + mrSges(6,2) * t313;
t290 = mrSges(6,1) * t324 - mrSges(6,3) * t313;
t395 = -t226 * t366 + t371 * t227;
t213 = m(6) * t236 - mrSges(6,2) * t303 + mrSges(6,3) * t265 + t279 * t312 - t290 * t324 + t395;
t235 = t240 * t372 - t241 * t367;
t289 = -mrSges(6,2) * t324 + mrSges(6,3) * t312;
t232 = -pkin(5) * t303 - pkin(12) * t323 + t280 * t313 - t235;
t385 = -m(7) * t232 + t250 * mrSges(7,1) - mrSges(7,2) * t251 + t287 * t267 - t268 * t288;
t222 = m(6) * t235 + mrSges(6,1) * t303 - mrSges(6,3) * t266 - t279 * t313 + t289 * t324 + t385;
t207 = t372 * t213 - t367 * t222;
t245 = t304 * pkin(3) + t382;
t316 = mrSges(5,1) * t327 - mrSges(5,2) * t342;
t388 = m(5) * t245 - t305 * mrSges(5,3) - t327 * t316 + t207;
t201 = m(4) * t258 + t305 * mrSges(4,2) + t304 * t420 + t327 * t317 + t326 * t410 + t388;
t309 = mrSges(4,1) * t326 + mrSges(4,2) * t327;
t206 = t367 * t213 + t372 * t222;
t310 = -mrSges(5,2) * t326 - mrSges(5,3) * t327;
t384 = -m(5) * t248 - t305 * mrSges(5,1) - t327 * t310 - t206;
t204 = m(4) * t252 - t305 * mrSges(4,3) - t327 * t309 + t331 * t420 - t342 * t410 + t384;
t214 = -m(6) * t243 + t265 * mrSges(6,1) - t266 * mrSges(6,2) + t312 * t289 - t313 * t290 - t216;
t379 = -m(5) * t246 + t331 * mrSges(5,3) - t342 * t316 - t214;
t211 = m(4) * t253 + t342 * t317 - t331 * mrSges(4,2) + (-t309 - t310) * t326 + (-mrSges(4,3) - mrSges(5,1)) * t304 + t379;
t191 = t417 * t201 + t204 * t403 + t211 * t415;
t192 = -t364 * t201 + t204 * t393 + t211 * t398;
t320 = -g(3) * t413 + t394;
t348 = -mrSges(3,2) * t391 + mrSges(3,3) * t399;
t351 = (-mrSges(3,1) * t373 + mrSges(3,2) * t369) * t408;
t189 = m(3) * t320 + t361 * mrSges(3,1) - t352 * mrSges(3,3) + t348 * t391 - t351 * t400 + t192;
t196 = -t368 * t204 + t211 * t423;
t321 = -g(3) * t414 + t409;
t347 = mrSges(3,1) * t391 - mrSges(3,3) * t400;
t195 = m(3) * t321 - t361 * mrSges(3,2) + t353 * mrSges(3,3) - t347 * t391 + t351 * t399 + t196;
t186 = -t189 * t369 + t373 * t195;
t335 = -t365 * t349 - t402;
t190 = m(3) * t335 - t353 * mrSges(3,1) + t352 * mrSges(3,2) + (t347 * t369 - t348 * t373) * t408 + t191;
t181 = t189 * t396 - t190 * t365 + t195 * t397;
t294 = Ifges(4,4) * t327 - Ifges(4,2) * t326 - Ifges(4,6) * t342;
t254 = Ifges(7,5) * t288 + Ifges(7,6) * t287 + Ifges(7,3) * t311;
t256 = Ifges(7,1) * t288 + Ifges(7,4) * t287 + Ifges(7,5) * t311;
t220 = -mrSges(7,1) * t232 + mrSges(7,3) * t229 + Ifges(7,4) * t251 + Ifges(7,2) * t250 + Ifges(7,6) * t264 - t254 * t288 + t256 * t311;
t255 = Ifges(7,4) * t288 + Ifges(7,2) * t287 + Ifges(7,6) * t311;
t221 = mrSges(7,2) * t232 - mrSges(7,3) * t228 + Ifges(7,1) * t251 + Ifges(7,4) * t250 + Ifges(7,5) * t264 + t254 * t287 - t255 * t311;
t269 = Ifges(6,5) * t313 + Ifges(6,6) * t312 + Ifges(6,3) * t324;
t270 = Ifges(6,4) * t313 + Ifges(6,2) * t312 + Ifges(6,6) * t324;
t198 = mrSges(6,2) * t243 - mrSges(6,3) * t235 + Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t303 - pkin(12) * t216 - t220 * t366 + t221 * t371 + t269 * t312 - t270 * t324;
t271 = Ifges(6,1) * t313 + Ifges(6,4) * t312 + Ifges(6,5) * t324;
t377 = mrSges(7,1) * t228 - mrSges(7,2) * t229 + Ifges(7,5) * t251 + Ifges(7,6) * t250 + Ifges(7,3) * t264 + t255 * t288 - t256 * t287;
t199 = -mrSges(6,1) * t243 + mrSges(6,3) * t236 + Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t303 - pkin(5) * t216 - t269 * t313 + t271 * t324 - t377;
t291 = -Ifges(5,5) * t342 - Ifges(5,6) * t327 + Ifges(5,3) * t326;
t381 = mrSges(5,2) * t248 - mrSges(5,3) * t246 + Ifges(5,1) * t331 - Ifges(5,4) * t305 + Ifges(5,5) * t304 - pkin(11) * t206 + t372 * t198 - t367 * t199 - t327 * t291;
t182 = -mrSges(4,2) * t253 + qJ(4) * (-t304 * mrSges(5,1) + t379) + t327 * t294 + Ifges(4,3) * t331 - Ifges(4,6) * t304 + Ifges(4,5) * t305 + pkin(3) * (-t331 * mrSges(5,2) + t342 * t315 + t384) + (-qJ(4) * t310 - t411) * t326 + t381 + mrSges(4,1) * t252;
t205 = -t304 * mrSges(5,2) - t326 * t315 + t388;
t378 = -mrSges(5,1) * t246 + mrSges(5,2) * t245 - pkin(4) * t214 - pkin(11) * t207 - t367 * t198 - t372 * t199;
t183 = -mrSges(4,1) * t258 + mrSges(4,3) * t253 - pkin(3) * t205 + t411 * t342 + (Ifges(4,6) - Ifges(5,5)) * t331 + t412 * t327 - t419 * t305 + (-Ifges(4,2) - Ifges(5,3)) * t304 + t378;
t380 = mrSges(6,1) * t235 - mrSges(6,2) * t236 + Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t303 + pkin(5) * t385 + pkin(12) * t395 + t371 * t220 + t366 * t221 + t313 * t270 - t312 * t271;
t376 = mrSges(5,1) * t248 - mrSges(5,3) * t245 + pkin(4) * t206 + t380;
t187 = mrSges(4,2) * t258 - qJ(4) * t205 + t376 + (-t291 + t294) * t342 + (Ifges(4,5) - Ifges(5,4)) * t331 + t412 * t326 + (Ifges(5,2) + Ifges(4,1)) * t305 + t419 * t304 - mrSges(4,3) * t252;
t333 = Ifges(3,6) * qJD(2) + (Ifges(3,6) * t418 + (Ifges(3,4) * t369 + Ifges(3,2) * t373) * t365) * qJD(1);
t334 = Ifges(3,5) * qJD(2) + (Ifges(3,5) * t418 + (Ifges(3,1) * t369 + Ifges(3,4) * t373) * t365) * qJD(1);
t173 = mrSges(3,1) * t320 - mrSges(3,2) * t321 + t417 * t182 + Ifges(3,5) * t352 + Ifges(3,6) * t353 + Ifges(3,3) * t361 + pkin(2) * t192 + (t333 * t369 - t334 * t373) * t408 + (pkin(10) * t196 + t183 * t423 + t187 * t368) * t364;
t332 = Ifges(3,3) * qJD(2) + (Ifges(3,3) * t418 + (Ifges(3,5) * t369 + Ifges(3,6) * t373) * t365) * qJD(1);
t175 = -mrSges(3,1) * t335 + mrSges(3,3) * t321 + Ifges(3,4) * t352 + Ifges(3,2) * t353 + Ifges(3,6) * t361 - pkin(2) * t191 - t364 * t182 + t183 * t393 + t187 * t398 + t196 * t401 - t332 * t400 + t334 * t391;
t177 = Ifges(3,1) * t352 + Ifges(3,4) * t353 + Ifges(3,5) * t361 + t332 * t399 - t391 * t333 + mrSges(3,2) * t335 - mrSges(3,3) * t320 + t423 * t187 - t368 * t183 + (-t364 * t191 - t192 * t417) * pkin(10);
t383 = mrSges(2,1) * t359 - mrSges(2,2) * t360 + Ifges(2,3) * qJDD(1) + pkin(1) * t181 + t418 * t173 + t175 * t413 + t177 * t414 + t186 * t422;
t184 = m(2) * t360 - mrSges(2,1) * t375 - qJDD(1) * mrSges(2,2) + t186;
t180 = t418 * t190 + (t189 * t373 + t195 * t369) * t365;
t178 = m(2) * t359 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t375 + t181;
t171 = -mrSges(2,2) * g(3) - mrSges(2,3) * t359 + Ifges(2,5) * qJDD(1) - t375 * Ifges(2,6) - t369 * t175 + t373 * t177 + (-t180 * t365 - t181 * t418) * pkin(9);
t170 = pkin(9) * t186 * t418 + mrSges(2,1) * g(3) + mrSges(2,3) * t360 + t375 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t180 - t365 * t173 + t175 * t396 + t177 * t397;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t374 * t171 - t370 * t170 - pkin(8) * (t178 * t374 + t184 * t370), t171, t177, t187, -t326 * t293 + t381, t198, t221; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t370 * t171 + t374 * t170 + pkin(8) * (-t178 * t370 + t184 * t374), t170, t175, t183, Ifges(5,4) * t331 - Ifges(5,2) * t305 + Ifges(5,6) * t304 + t342 * t291 + t326 * t295 - t376, t199, t220; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t383, t383, t173, t182, Ifges(5,5) * t331 - Ifges(5,6) * t305 + Ifges(5,3) * t304 - t342 * t293 + t327 * t295 - t378, t380, t377;];
m_new  = t1;
