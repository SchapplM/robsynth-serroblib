% Calculate vector of cutting torques with Newton-Euler for
% S6RRPPRR7
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
% Datum: 2019-05-06 11:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RRPPRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 11:02:55
% EndTime: 2019-05-06 11:03:17
% DurationCPUTime: 9.96s
% Computational Cost: add. (153151->420), mult. (352895->499), div. (0->0), fcn. (239696->10), ass. (0->168)
t370 = sin(qJ(1));
t374 = cos(qJ(1));
t351 = t370 * g(1) - t374 * g(2);
t375 = qJD(1) ^ 2;
t365 = sin(pkin(6));
t431 = pkin(8) * t365;
t326 = qJDD(1) * pkin(1) + t375 * t431 + t351;
t366 = cos(pkin(6));
t287 = -t366 * g(3) - t365 * t326;
t369 = sin(qJ(2));
t373 = cos(qJ(2));
t409 = qJD(1) * t373;
t333 = (qJD(2) * t409 + qJDD(1) * t369) * t365;
t410 = qJD(1) * t365;
t401 = t369 * t410;
t406 = qJDD(1) * t365;
t334 = -qJD(2) * t401 + t373 * t406;
t357 = t366 * qJD(1) + qJD(2);
t402 = t365 * t409;
t398 = t357 * t402;
t397 = -t334 * pkin(2) + t287 + (-t333 - t398) * qJ(3);
t432 = pkin(2) * t357;
t246 = (-(2 * qJD(3)) + t432) * t401 + t397;
t437 = m(4) * t246 - t334 * mrSges(4,1);
t356 = t366 * qJDD(1) + qJDD(2);
t436 = Ifges(5,5) * t334 + Ifges(5,6) * t333 + Ifges(5,3) * t356;
t319 = -t357 * pkin(3) - qJ(4) * t401;
t420 = t365 ^ 2 * t375;
t404 = t373 ^ 2 * t420;
t386 = -qJ(4) * t404 + qJDD(4) - t397 + ((2 * qJD(3)) + t319) * t401;
t433 = pkin(3) + pkin(9);
t228 = t333 * pkin(4) + t433 * t334 + (pkin(4) * t373 + (-pkin(2) - pkin(9)) * t369) * t357 * t410 + t386;
t332 = (pkin(4) * t369 + pkin(9) * t373) * t410;
t352 = -t374 * g(1) - t370 * g(2);
t327 = -t375 * pkin(1) + pkin(8) * t406 + t352;
t416 = t366 * t373;
t418 = t365 * t373;
t267 = -g(3) * t418 + t326 * t416 - t369 * t327;
t328 = (-pkin(2) * t373 - qJ(3) * t369) * t410;
t434 = t357 ^ 2;
t391 = -qJ(3) * t434 + t328 * t401 + qJDD(3) - t267;
t407 = qJD(1) * qJD(4);
t399 = -0.2e1 * t365 * t407;
t382 = t369 * t399 + t391 + (-t333 + t398) * qJ(4);
t403 = t373 * t420;
t232 = -t434 * pkin(4) + (-pkin(3) * t403 - t332 * t410) * t369 + (-pkin(2) - t433) * t356 + t382;
t368 = sin(qJ(5));
t372 = cos(qJ(5));
t225 = t368 * t228 + t372 * t232;
t302 = -t368 * t357 - t372 * t402;
t265 = -t302 * qJD(5) + t368 * t334 - t372 * t356;
t301 = -t372 * t357 + t368 * t402;
t269 = -mrSges(6,1) * t301 + mrSges(6,2) * t302;
t340 = qJD(5) + t401;
t274 = mrSges(6,1) * t340 - mrSges(6,3) * t302;
t317 = qJDD(5) + t333;
t270 = -pkin(5) * t301 - pkin(10) * t302;
t338 = t340 ^ 2;
t222 = -pkin(5) * t338 + pkin(10) * t317 + t270 * t301 + t225;
t408 = qJD(3) * t357;
t339 = 0.2e1 * t408;
t417 = t366 * t369;
t412 = t326 * t417 + t373 * t327;
t396 = pkin(2) * t434 - t356 * qJ(3) - t328 * t402 - t412;
t387 = pkin(3) * t404 + t334 * qJ(4) - t357 * t319 + t396;
t430 = g(3) * t369;
t230 = t356 * pkin(4) - t434 * pkin(9) + t339 + t373 * t399 + (-t332 * t409 - t430) * t365 - t387;
t266 = t301 * qJD(5) - t372 * t334 - t368 * t356;
t226 = t230 + (-t301 * t340 - t266) * pkin(10) + (t302 * t340 - t265) * pkin(5);
t367 = sin(qJ(6));
t371 = cos(qJ(6));
t219 = -t222 * t367 + t226 * t371;
t271 = -t367 * t302 + t371 * t340;
t241 = t271 * qJD(6) + t371 * t266 + t367 * t317;
t272 = t371 * t302 + t367 * t340;
t253 = -mrSges(7,1) * t271 + mrSges(7,2) * t272;
t299 = qJD(6) - t301;
t255 = -mrSges(7,2) * t299 + mrSges(7,3) * t271;
t263 = qJDD(6) - t265;
t215 = m(7) * t219 + mrSges(7,1) * t263 - t241 * mrSges(7,3) - t253 * t272 + t255 * t299;
t220 = t222 * t371 + t226 * t367;
t240 = -t272 * qJD(6) - t367 * t266 + t371 * t317;
t256 = mrSges(7,1) * t299 - mrSges(7,3) * t272;
t216 = m(7) * t220 - mrSges(7,2) * t263 + t240 * mrSges(7,3) + t253 * t271 - t256 * t299;
t400 = -t215 * t367 + t371 * t216;
t201 = m(6) * t225 - mrSges(6,2) * t317 + mrSges(6,3) * t265 + t269 * t301 - t274 * t340 + t400;
t224 = t228 * t372 - t232 * t368;
t273 = -mrSges(6,2) * t340 + mrSges(6,3) * t301;
t221 = -pkin(5) * t317 - pkin(10) * t338 + t270 * t302 - t224;
t394 = -m(7) * t221 + t240 * mrSges(7,1) - t241 * mrSges(7,2) + t271 * t255 - t256 * t272;
t211 = m(6) * t224 + mrSges(6,1) * t317 - mrSges(6,3) * t266 - t269 * t302 + t273 * t340 + t394;
t195 = t372 * t201 - t368 * t211;
t429 = t356 * pkin(2);
t235 = -t429 + (-t369 * t403 - t356) * pkin(3) + t382;
t323 = -t357 * mrSges(5,1) + mrSges(5,3) * t402;
t331 = (mrSges(5,1) * t369 - mrSges(5,2) * t373) * t410;
t395 = m(5) * t235 + t356 * mrSges(5,2) + t357 * t323 - t331 * t401 + t195;
t192 = -t333 * mrSges(5,3) + t395;
t252 = t391 - t429;
t279 = Ifges(4,6) * t357 + (Ifges(4,5) * t369 - Ifges(4,3) * t373) * t410;
t282 = Ifges(4,2) * t357 + (Ifges(4,4) * t369 - Ifges(4,6) * t373) * t410;
t435 = mrSges(4,2) * t252 - mrSges(4,3) * t246 + Ifges(4,1) * t333 + Ifges(4,4) * t356 - Ifges(4,5) * t334 - qJ(4) * t192 + t357 * t279 + t282 * t402;
t428 = mrSges(3,3) + mrSges(4,2);
t427 = Ifges(3,5) + Ifges(5,6);
t320 = t357 * mrSges(5,2) - mrSges(5,3) * t401;
t421 = t357 * t320;
t419 = t365 * t369;
t204 = t371 * t215 + t367 * t216;
t278 = -Ifges(5,3) * t357 + (-Ifges(5,5) * t373 - Ifges(5,6) * t369) * t410;
t415 = -t278 + t282;
t284 = -Ifges(5,5) * t357 + (-Ifges(5,1) * t373 - Ifges(5,4) * t369) * t410;
t414 = -Ifges(3,6) * t357 - (Ifges(3,4) * t369 + Ifges(3,2) * t373) * t410 + t284;
t285 = Ifges(4,4) * t357 + (Ifges(4,1) * t369 - Ifges(4,5) * t373) * t410;
t413 = t285 + Ifges(3,5) * t357 + (Ifges(3,1) * t369 + Ifges(3,4) * t373) * t410;
t411 = -t357 * mrSges(3,1) + mrSges(3,3) * t401 + t320;
t405 = g(3) * t419;
t324 = -t357 * mrSges(3,2) + mrSges(3,3) * t402;
t329 = (-mrSges(4,1) * t373 - mrSges(4,3) * t369) * t410;
t330 = (-mrSges(3,1) * t373 + mrSges(3,2) * t369) * t410;
t325 = mrSges(4,2) * t402 + t357 * mrSges(4,3);
t390 = -m(4) * t252 + t356 * mrSges(4,1) + t357 * t325 - t395;
t189 = m(3) * t267 + t356 * mrSges(3,1) + t357 * t324 + (-t329 - t330) * t401 + (mrSges(5,3) - t428) * t333 + t390;
t268 = -t405 + t412;
t245 = t339 - t396 - t405;
t322 = -t357 * mrSges(4,1) + mrSges(4,2) * t401;
t202 = -m(6) * t230 + t265 * mrSges(6,1) - t266 * mrSges(6,2) + t301 * t273 - t302 * t274 - t204;
t233 = -0.2e1 * t408 + (0.2e1 * t373 * t407 + t430) * t365 + t387;
t385 = -m(5) * t233 - t334 * mrSges(5,3) - t202;
t380 = m(4) * t245 + t356 * mrSges(4,3) + t357 * t322 + t329 * t402 + t385;
t198 = t380 + (t330 - t331) * t402 + t411 * t357 + (-mrSges(3,2) + mrSges(5,1)) * t356 + t428 * t334 + m(3) * t268;
t183 = -t189 * t369 + t373 * t198;
t194 = t368 * t201 + t372 * t211;
t238 = t334 * pkin(3) - t401 * t432 + t386;
t393 = -m(5) * t238 - t333 * mrSges(5,1) + t323 * t402 - t194;
t187 = m(3) * t287 + (-mrSges(3,1) + mrSges(5,2)) * t334 + (mrSges(3,2) - mrSges(4,3)) * t333 + ((-t324 - t325) * t373 + (-t322 - t411) * t369) * t410 + t393 + t437;
t180 = -t187 * t365 + t189 * t416 + t198 * t417;
t247 = Ifges(7,5) * t272 + Ifges(7,6) * t271 + Ifges(7,3) * t299;
t249 = Ifges(7,1) * t272 + Ifges(7,4) * t271 + Ifges(7,5) * t299;
t208 = -mrSges(7,1) * t221 + mrSges(7,3) * t220 + Ifges(7,4) * t241 + Ifges(7,2) * t240 + Ifges(7,6) * t263 - t247 * t272 + t249 * t299;
t248 = Ifges(7,4) * t272 + Ifges(7,2) * t271 + Ifges(7,6) * t299;
t209 = mrSges(7,2) * t221 - mrSges(7,3) * t219 + Ifges(7,1) * t241 + Ifges(7,4) * t240 + Ifges(7,5) * t263 + t247 * t271 - t248 * t299;
t257 = Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t340;
t258 = Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t340;
t185 = mrSges(6,2) * t230 - mrSges(6,3) * t224 + Ifges(6,1) * t266 + Ifges(6,4) * t265 + Ifges(6,5) * t317 - pkin(10) * t204 - t208 * t367 + t209 * t371 + t257 * t301 - t258 * t340;
t259 = Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t340;
t381 = mrSges(7,1) * t219 - mrSges(7,2) * t220 + Ifges(7,5) * t241 + Ifges(7,6) * t240 + Ifges(7,3) * t263 + t248 * t272 - t249 * t271;
t186 = -mrSges(6,1) * t230 + mrSges(6,3) * t225 + Ifges(6,4) * t266 + Ifges(6,2) * t265 + Ifges(6,6) * t317 - pkin(5) * t204 - t257 * t302 + t259 * t340 - t381;
t281 = -Ifges(5,6) * t357 + (-Ifges(5,4) * t373 - Ifges(5,2) * t369) * t410;
t383 = -mrSges(5,1) * t233 + mrSges(5,2) * t235 - pkin(4) * t202 - pkin(9) * t195 - t368 * t185 - t372 * t186 + t281 * t402;
t377 = -mrSges(4,1) * t252 + mrSges(4,3) * t245 + Ifges(4,4) * t333 + Ifges(4,2) * t356 - Ifges(4,6) * t334 - pkin(3) * t192 + t383;
t172 = t377 + ((-qJ(3) * t331 - t413) * t373 + (-pkin(2) * t329 - t279 - t414) * t369) * t410 + (mrSges(5,1) * qJ(3) + Ifges(3,3) + Ifges(5,3)) * t356 + (mrSges(4,2) * qJ(3) + Ifges(5,5) + Ifges(3,6)) * t334 + (pkin(2) * (-mrSges(4,2) + mrSges(5,3)) + t427) * t333 + mrSges(3,1) * t267 - mrSges(3,2) * t268 + qJ(3) * (t380 + t421) + pkin(2) * t390;
t389 = t334 * mrSges(5,2) + t393;
t190 = -t333 * mrSges(4,3) + (-t325 * t373 + (-t320 - t322) * t369) * t410 + t389 + t437;
t280 = Ifges(3,3) * t357 + (Ifges(3,5) * t369 + Ifges(3,6) * t373) * t410;
t388 = -mrSges(5,2) * t238 + mrSges(5,3) * t233 + Ifges(5,1) * t334 + Ifges(5,4) * t333 + Ifges(5,5) * t356 + pkin(9) * t194 - t372 * t185 + t368 * t186;
t378 = mrSges(4,1) * t246 - mrSges(4,2) * t245 + pkin(3) * (-t320 * t401 + t389) + qJ(4) * (t356 * mrSges(5,1) - t331 * t402 + t385 + t421) - t388;
t174 = -t378 + (-t280 - t415) * t401 + (-t281 + t413) * t357 + (Ifges(3,6) - Ifges(4,6)) * t356 + (Ifges(4,3) + Ifges(3,2)) * t334 + (Ifges(3,4) - Ifges(4,5)) * t333 - mrSges(3,1) * t287 + mrSges(3,3) * t268 - pkin(2) * t190;
t384 = mrSges(6,1) * t224 - mrSges(6,2) * t225 + Ifges(6,5) * t266 + Ifges(6,6) * t265 + Ifges(6,3) * t317 + pkin(5) * t394 + pkin(10) * t400 + t371 * t208 + t367 * t209 + t302 * t258 - t301 * t259;
t379 = mrSges(5,1) * t238 - mrSges(5,3) * t235 + pkin(4) * t194 + t384;
t176 = t435 + t379 + (-t278 + t280) * t402 + (Ifges(5,4) + Ifges(3,4)) * t334 + (Ifges(3,1) + Ifges(5,2)) * t333 + t427 * t356 + t414 * t357 + mrSges(3,2) * t287 - mrSges(3,3) * t267 - qJ(3) * t190;
t392 = mrSges(2,1) * t351 - mrSges(2,2) * t352 + Ifges(2,3) * qJDD(1) + pkin(1) * t180 + t366 * t172 + t174 * t418 + t176 * t419 + t183 * t431;
t376 = -Ifges(5,4) * t334 - Ifges(5,2) * t333 - Ifges(5,6) * t356 + t278 * t402 - t357 * t284 - t379;
t181 = m(2) * t352 - mrSges(2,1) * t375 - qJDD(1) * mrSges(2,2) + t183;
t179 = t366 * t187 + (t189 * t373 + t198 * t369) * t365;
t177 = m(2) * t351 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t375 + t180;
t170 = -mrSges(2,2) * g(3) - mrSges(2,3) * t351 + Ifges(2,5) * qJDD(1) - t375 * Ifges(2,6) - t369 * t174 + t373 * t176 + (-t179 * t365 - t180 * t366) * pkin(8);
t169 = mrSges(2,1) * g(3) + mrSges(2,3) * t352 + t375 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t179 - t365 * t172 + (pkin(8) * t183 + t174 * t373 + t176 * t369) * t366;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t374 * t170 - t370 * t169 - pkin(7) * (t177 * t374 + t181 * t370), t170, t176, -t376 + t435, -t278 * t401 + t357 * t281 - t388, t185, t209; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t370 * t170 + t374 * t169 + pkin(7) * (-t177 * t370 + t181 * t374), t169, t174, (-t285 * t373 + (-t279 - t284) * t369) * t410 + t377 + t436, t376, t186, t208; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t392, t392, t172, t378 + t415 * t401 + (-t285 + t281) * t357 + Ifges(4,6) * t356 - Ifges(4,3) * t334 + Ifges(4,5) * t333, t284 * t401 - t383 - t436, t384, t381;];
m_new  = t1;
