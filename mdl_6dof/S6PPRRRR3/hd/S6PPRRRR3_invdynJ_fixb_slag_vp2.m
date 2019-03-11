% Calculate vector of inverse dynamics joint torques for
% S6PPRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:20
% EndTime: 2019-03-08 19:08:04
% DurationCPUTime: 24.04s
% Computational Cost: add. (14705->793), mult. (38955->1189), div. (0->0), fcn. (36500->18), ass. (0->372)
t432 = cos(pkin(6));
t261 = qJD(1) * t432 + qJD(2);
t273 = cos(pkin(7));
t281 = cos(qJ(3));
t271 = cos(pkin(14));
t270 = sin(pkin(6));
t400 = qJD(1) * t270;
t380 = t271 * t400;
t269 = sin(pkin(7));
t412 = t269 * t281;
t267 = sin(pkin(14));
t277 = sin(qJ(3));
t418 = t267 * t277;
t174 = t273 * t281 * t380 + t261 * t412 - t400 * t418;
t413 = t269 * t277;
t411 = t271 * t273;
t503 = (t267 * t281 + t277 * t411) * t270;
t175 = qJD(1) * t503 + t261 * t413;
t268 = sin(pkin(8));
t276 = sin(qJ(4));
t416 = t268 * t276;
t262 = pkin(10) * t416;
t272 = cos(pkin(8));
t280 = cos(qJ(4));
t409 = t272 * t280;
t238 = pkin(3) * t409 - t262;
t410 = t272 * t276;
t534 = -t238 * qJD(4) + t174 * t280 - t175 * t410;
t279 = cos(qJ(5));
t450 = pkin(11) * t279;
t414 = t268 * t280;
t239 = pkin(3) * t410 + pkin(10) * t414;
t221 = pkin(11) * t272 + t239;
t350 = -pkin(4) * t280 - pkin(11) * t276;
t222 = (-pkin(3) + t350) * t268;
t328 = t268 * (pkin(4) * t276 - pkin(11) * t280);
t230 = qJD(4) * t328;
t275 = sin(qJ(5));
t394 = qJD(5) * t279;
t395 = qJD(5) * t275;
t417 = t268 * t275;
t510 = -t175 * t417 - t221 * t395 + t222 * t394 + t275 * t230 - t279 * t534;
t496 = t239 * qJD(4) - t174 * t276 - t175 * t409;
t524 = -m(7) - m(6);
t397 = qJD(4) * t276;
t376 = t268 * t397;
t533 = -pkin(12) * t376 - t510;
t236 = -t279 * t272 + t275 * t416;
t396 = qJD(4) * t280;
t375 = t268 * t396;
t193 = -qJD(5) * t236 + t279 * t375;
t415 = t268 * t279;
t237 = t272 * t275 + t276 * t415;
t194 = qJD(5) * t237 + t275 * t375;
t532 = pkin(5) * t194 - pkin(12) * t193 + t496;
t366 = t269 * t432;
t299 = t281 * t366 + (t281 * t411 - t418) * t270;
t310 = -t270 * t271 * t269 + t273 * t432;
t306 = t310 * t268;
t531 = t299 * t272 + t306;
t399 = qJD(3) * t268;
t379 = t276 * t399;
t229 = qJD(3) * t328;
t167 = pkin(10) * t399 + t175;
t210 = t261 * t273 - t269 * t380;
t169 = qJD(3) * pkin(3) + t174;
t423 = t169 * t272;
t332 = t210 * t268 + t423;
t78 = -t276 * t167 + t280 * t332;
t65 = t275 * t229 + t279 * t78;
t530 = pkin(11) * t395 + pkin(12) * t379 + t65;
t349 = pkin(5) * t275 - pkin(12) * t279;
t398 = qJD(3) * t280;
t424 = t167 * t280;
t529 = -qJD(6) * t450 + t349 * qJD(5) - t169 * t410 - t424 - (t210 * t276 + t349 * t398) * t268;
t498 = t279 * t221 + t275 * t222;
t509 = -qJD(5) * t498 - t175 * t415 + t230 * t279 + t275 * t534;
t259 = qJDD(1) * t432 + qJDD(2);
t388 = qJDD(1) * t270;
t371 = t271 * t388;
t357 = t273 * t371;
t372 = t267 * t388;
t118 = qJD(3) * t174 + t259 * t413 + t277 * t357 + t281 * t372;
t528 = pkin(10) * qJDD(3) * t268 + qJD(4) * t423 + t118;
t274 = sin(qJ(6));
t278 = cos(qJ(6));
t344 = mrSges(7,1) * t274 + mrSges(7,2) * t278;
t480 = -t344 + mrSges(5,2) - mrSges(6,3);
t185 = t277 * t366 + t503;
t403 = t279 * t280;
t202 = (-t274 * t403 + t276 * t278) * t399;
t527 = t274 * t394 + t202;
t345 = -mrSges(7,1) * t278 + mrSges(7,2) * t274;
t319 = m(7) * pkin(5) - t345;
t346 = -mrSges(6,1) * t279 + mrSges(6,2) * t275;
t385 = m(7) * pkin(12) + mrSges(7,3);
t526 = -pkin(4) * t524 + t275 * t385 + t279 * t319 + mrSges(5,1) - t346;
t260 = qJD(3) * t272 + qJD(4);
t212 = t260 * t275 + t279 * t379;
t438 = Ifges(6,4) * t212;
t378 = t268 * t398;
t254 = qJD(5) - t378;
t512 = t254 * Ifges(6,6);
t211 = t260 * t279 - t275 * t379;
t513 = t211 * Ifges(6,2);
t138 = t438 + t512 + t513;
t200 = t272 * t210;
t112 = t200 + (qJD(3) * t350 - t169) * t268;
t79 = t276 * t332 + t424;
t72 = pkin(11) * t260 + t79;
t42 = t112 * t275 + t279 * t72;
t38 = pkin(12) * t254 + t42;
t71 = -pkin(4) * t260 - t78;
t51 = -pkin(5) * t211 - pkin(12) * t212 + t71;
t11 = t274 * t51 + t278 * t38;
t520 = t11 * mrSges(7,2);
t10 = -t274 * t38 + t278 * t51;
t521 = t10 * mrSges(7,1);
t525 = -t138 / 0.2e1 - t520 + t521;
t389 = qJD(3) * qJD(4);
t234 = (qJDD(3) * t276 + t280 * t389) * t268;
t258 = qJDD(3) * t272 + qJDD(4);
t149 = qJD(5) * t211 + t234 * t279 + t258 * t275;
t172 = -t212 * t274 + t254 * t278;
t233 = (-qJDD(3) * t280 + t276 * t389) * t268;
t225 = qJDD(5) + t233;
t82 = qJD(6) * t172 + t149 * t278 + t225 * t274;
t473 = t82 / 0.2e1;
t173 = t212 * t278 + t254 * t274;
t83 = -qJD(6) * t173 - t149 * t274 + t225 * t278;
t472 = t83 / 0.2e1;
t119 = -qJD(3) * t175 + t259 * t412 - t277 * t372 + t281 * t357;
t116 = qJDD(3) * pkin(3) + t119;
t206 = t259 * t273 - t269 * t371;
t17 = t116 * t410 - t167 * t397 + t206 * t416 + t210 * t375 + t280 * t528;
t15 = pkin(11) * t258 + t17;
t91 = -t116 * t268 + t272 * t206;
t63 = pkin(4) * t233 - pkin(11) * t234 + t91;
t7 = t112 * t394 + t279 * t15 + t275 * t63 - t395 * t72;
t5 = pkin(12) * t225 + t7;
t150 = -qJD(5) * t212 - t234 * t275 + t258 * t279;
t18 = t280 * (t116 * t272 + t206 * t268) - t167 * t396 - t210 * t376 - t528 * t276;
t16 = -pkin(4) * t258 - t18;
t9 = -pkin(5) * t150 - pkin(12) * t149 + t16;
t1 = qJD(6) * t10 + t274 * t9 + t278 * t5;
t523 = t1 * mrSges(7,2);
t2 = -qJD(6) * t11 - t274 * t5 + t278 * t9;
t522 = t2 * mrSges(7,1);
t146 = qJDD(6) - t150;
t465 = t146 / 0.2e1;
t464 = t149 / 0.2e1;
t463 = t150 / 0.2e1;
t453 = t225 / 0.2e1;
t519 = t17 * mrSges(5,2);
t518 = t18 * mrSges(5,1);
t220 = t262 + (-pkin(3) * t280 - pkin(4)) * t272;
t153 = pkin(5) * t236 - pkin(12) * t237 + t220;
t155 = -pkin(12) * t414 + t498;
t90 = t153 * t274 + t155 * t278;
t516 = -qJD(6) * t90 + t274 * t533 + t278 * t532;
t89 = t153 * t278 - t155 * t274;
t515 = qJD(6) * t89 + t274 * t532 - t278 * t533;
t514 = -pkin(5) * t376 - t509;
t511 = mrSges(6,1) + t319;
t486 = mrSges(6,2) - t385;
t120 = mrSges(6,1) * t225 - mrSges(6,3) * t149;
t36 = -mrSges(7,1) * t83 + mrSges(7,2) * t82;
t508 = -t120 + t36;
t253 = -pkin(5) * t279 - pkin(12) * t275 - pkin(4);
t391 = qJD(6) * t278;
t507 = t253 * t391 + t274 * t529 - t278 * t530;
t393 = qJD(6) * t274;
t506 = -t253 * t393 + t274 * t530 + t278 * t529;
t191 = mrSges(5,1) * t258 - mrSges(5,3) * t234;
t88 = -mrSges(6,1) * t150 + mrSges(6,2) * t149;
t505 = -t191 + t88;
t288 = -t268 * t299 + t272 * t310;
t181 = t185 * qJD(3);
t180 = t299 * qJD(3);
t488 = qJD(4) * t531 + t180;
t50 = -t181 * t410 - t185 * t397 + t280 * t488;
t504 = qJD(5) * t288 + t50;
t117 = -mrSges(7,1) * t172 + mrSges(7,2) * t173;
t442 = mrSges(6,3) * t212;
t183 = mrSges(6,1) * t254 - t442;
t501 = -t183 + t117;
t500 = t275 * t391 + t527;
t203 = (t274 * t276 + t278 * t403) * t399;
t392 = qJD(6) * t275;
t499 = t274 * t392 - t278 * t394 + t203;
t363 = mrSges(5,3) * t379;
t497 = mrSges(5,1) * t260 + mrSges(6,1) * t211 - mrSges(6,2) * t212 - t363;
t431 = cos(pkin(13));
t348 = t432 * t431;
t430 = sin(pkin(13));
t302 = -t267 * t430 + t271 * t348;
t365 = t270 * t431;
t495 = -t269 * t365 + t302 * t273;
t402 = t280 * t281;
t406 = t276 * t277;
t494 = t272 * t402 - t406;
t322 = (-mrSges(5,1) * t280 + mrSges(5,2) * t276) * t268;
t228 = qJD(3) * t322;
t493 = mrSges(4,1) * qJD(3) - t228 * t268;
t492 = t1 * t278 - t2 * t274;
t207 = qJD(6) - t211;
t128 = -t169 * t268 + t200;
t427 = t128 * t268;
t491 = -t260 * t268 * (Ifges(5,5) * t280 - Ifges(5,6) * t276) / 0.2e1 - (mrSges(5,1) * t276 + mrSges(5,2) * t280) * t427;
t487 = m(3) + m(4) + m(5) - t524;
t8 = -qJD(5) * t42 - t15 * t275 + t279 * t63;
t485 = -m(6) * t71 + t497;
t41 = t112 * t279 - t275 * t72;
t37 = -pkin(5) * t254 - t41;
t484 = -m(7) * t37 - t501;
t483 = t522 - t523;
t481 = pkin(11) * t524 + t480;
t94 = t173 * Ifges(7,5) + t172 * Ifges(7,6) + t207 * Ifges(7,3);
t478 = -mrSges(6,3) * t42 + t94 / 0.2e1 + t525;
t477 = t268 ^ 2;
t282 = qJD(3) ^ 2;
t22 = t82 * Ifges(7,4) + t83 * Ifges(7,2) + t146 * Ifges(7,6);
t476 = t22 / 0.2e1;
t475 = Ifges(7,1) * t473 + Ifges(7,4) * t472 + Ifges(7,5) * t465;
t474 = Ifges(6,1) * t464 + Ifges(6,4) * t463 + Ifges(6,5) * t453;
t435 = Ifges(7,4) * t173;
t95 = Ifges(7,2) * t172 + Ifges(7,6) * t207 + t435;
t470 = -t95 / 0.2e1;
t469 = t95 / 0.2e1;
t170 = Ifges(7,4) * t172;
t96 = Ifges(7,1) * t173 + Ifges(7,5) * t207 + t170;
t468 = -t96 / 0.2e1;
t467 = t96 / 0.2e1;
t462 = -t172 / 0.2e1;
t461 = t172 / 0.2e1;
t460 = -t173 / 0.2e1;
t459 = t173 / 0.2e1;
t458 = -t207 / 0.2e1;
t457 = t207 / 0.2e1;
t454 = t212 / 0.2e1;
t452 = pkin(3) * t268;
t6 = -pkin(5) * t225 - t8;
t447 = t275 * t6;
t446 = t279 * t7;
t443 = mrSges(6,3) * t211;
t441 = mrSges(6,3) * t275;
t440 = Ifges(5,4) * t276;
t439 = Ifges(5,4) * t280;
t437 = Ifges(6,4) * t275;
t436 = Ifges(6,4) * t279;
t434 = Ifges(7,4) * t274;
t433 = Ifges(7,4) * t278;
t428 = mrSges(4,2) * qJD(3);
t301 = t267 * t348 + t271 * t430;
t147 = t277 * t495 + t281 * t301;
t426 = t147 * t268;
t347 = t432 * t430;
t304 = -t267 * t431 - t271 * t347;
t364 = t270 * t430;
t292 = t269 * t364 + t273 * t304;
t303 = -t267 * t347 + t271 * t431;
t148 = t277 * t292 + t281 * t303;
t425 = t148 * t268;
t422 = t185 * t268;
t421 = t211 * t274;
t420 = t211 * t278;
t408 = t274 * t275;
t407 = t275 * t278;
t405 = t276 * t281;
t404 = t277 * t280;
t401 = mrSges(4,2) * qJDD(3);
t21 = Ifges(7,5) * t82 + Ifges(7,6) * t83 + Ifges(7,3) * t146;
t386 = pkin(11) * t394;
t382 = Ifges(6,5) * t149 + Ifges(6,6) * t150 + Ifges(6,3) * t225;
t381 = Ifges(5,5) * t234 - Ifges(5,6) * t233 + Ifges(5,3) * t258;
t373 = t416 / 0.2e1;
t368 = t394 / 0.2e1;
t367 = -t392 / 0.2e1;
t362 = mrSges(5,3) * t378;
t359 = t399 * t413;
t285 = t301 * t277 - t281 * t495;
t355 = -t285 * pkin(3) + pkin(10) * t426;
t286 = t277 * t303 - t281 * t292;
t354 = -t286 * pkin(3) + pkin(10) * t425;
t353 = t299 * pkin(3) + pkin(10) * t422;
t343 = Ifges(6,1) * t279 - t437;
t342 = Ifges(7,1) * t278 - t434;
t341 = Ifges(7,1) * t274 + t433;
t340 = -Ifges(6,2) * t275 + t436;
t339 = -Ifges(7,2) * t274 + t433;
t338 = Ifges(7,2) * t278 + t434;
t337 = Ifges(6,5) * t279 - Ifges(6,6) * t275;
t336 = Ifges(7,5) * t278 - Ifges(7,6) * t274;
t335 = Ifges(7,5) * t274 + Ifges(7,6) * t278;
t295 = t299 * t280;
t100 = t185 * t276 - t272 * t295 - t280 * t306;
t101 = t185 * t280 + t276 * t531;
t59 = t101 * t279 + t275 * t288;
t30 = t100 * t278 - t274 * t59;
t31 = t100 * t274 + t278 * t59;
t64 = t229 * t279 - t275 * t78;
t325 = t272 * t405 + t404;
t187 = t269 * t325 + t273 * t416;
t235 = -t268 * t412 + t272 * t273;
t152 = t187 * t279 + t235 * t275;
t186 = -t269 * t494 - t273 * t414;
t106 = t152 * t278 + t186 * t274;
t105 = -t152 * t274 + t186 * t278;
t151 = t187 * t275 - t279 * t235;
t164 = -t221 * t275 + t222 * t279;
t195 = -t237 * t274 - t278 * t414;
t327 = -t237 * t278 + t274 * t414;
t321 = (Ifges(5,2) * t280 + t440) * t268;
t313 = t276 * t477 * (Ifges(5,1) * t280 - t440);
t297 = t299 * mrSges(4,1);
t293 = -t269 * t304 + t273 * t364;
t291 = -t269 * t302 - t273 * t365;
t290 = t293 * t268;
t289 = t291 * t268;
t284 = t286 * t280;
t283 = t285 * t280;
t257 = Ifges(5,4) * t378;
t227 = -mrSges(5,2) * t260 + t362;
t218 = t253 * t274 + t278 * t450;
t217 = t253 * t278 - t274 * t450;
t205 = Ifges(6,4) * t211;
t198 = Ifges(5,1) * t379 + t260 * Ifges(5,5) + t257;
t197 = t260 * Ifges(5,6) + qJD(3) * t321;
t190 = -mrSges(5,2) * t258 - mrSges(5,3) * t233;
t182 = -mrSges(6,2) * t254 + t443;
t171 = mrSges(5,1) * t233 + mrSges(5,2) * t234;
t163 = pkin(5) * t212 - pkin(12) * t211;
t154 = pkin(5) * t414 - t164;
t141 = t273 * t376 + (t325 * qJD(4) + (t272 * t404 + t405) * qJD(3)) * t269;
t140 = t273 * t375 + (t494 * qJD(4) + (-t272 * t406 + t402) * qJD(3)) * t269;
t139 = t212 * Ifges(6,1) + t254 * Ifges(6,5) + t205;
t137 = t212 * Ifges(6,5) + t211 * Ifges(6,6) + t254 * Ifges(6,3);
t130 = mrSges(7,1) * t207 - mrSges(7,3) * t173;
t129 = -mrSges(7,2) * t207 + mrSges(7,3) * t172;
t127 = -t185 * t410 + t295;
t126 = t185 * t409 + t276 * t299;
t125 = qJD(6) * t327 - t193 * t274 + t278 * t376;
t124 = qJD(6) * t195 + t193 * t278 + t274 * t376;
t121 = -mrSges(6,2) * t225 + mrSges(6,3) * t150;
t103 = t268 * t286 + t272 * t293;
t102 = t268 * t285 + t272 * t291;
t87 = -t148 * t410 - t284;
t86 = t148 * t409 - t276 * t286;
t85 = -t147 * t410 - t283;
t84 = t147 * t409 - t276 * t285;
t74 = qJD(5) * t152 + t140 * t275 - t279 * t359;
t73 = -qJD(5) * t151 + t140 * t279 + t275 * t359;
t67 = t149 * Ifges(6,4) + t150 * Ifges(6,2) + t225 * Ifges(6,6);
t60 = -pkin(5) * t379 - t64;
t58 = t101 * t275 - t279 * t288;
t57 = t148 * t280 + (-t272 * t286 + t290) * t276;
t56 = t148 * t276 + t272 * t284 - t280 * t290;
t55 = t147 * t280 + (-t272 * t285 + t289) * t276;
t54 = t147 * t276 + t272 * t283 - t280 * t289;
t49 = t181 * t409 + t185 * t396 + t276 * t488;
t48 = -mrSges(7,2) * t146 + mrSges(7,3) * t83;
t47 = mrSges(7,1) * t146 - mrSges(7,3) * t82;
t35 = t103 * t275 + t279 * t57;
t33 = t102 * t275 + t279 * t55;
t29 = t163 * t274 + t278 * t41;
t28 = t163 * t278 - t274 * t41;
t27 = -qJD(6) * t106 + t141 * t278 - t274 * t73;
t26 = qJD(6) * t105 + t141 * t274 + t278 * t73;
t13 = -t101 * t395 + t181 * t417 + t279 * t504;
t4 = qJD(6) * t30 + t13 * t278 + t274 * t49;
t3 = -qJD(6) * t31 - t13 * t274 + t278 * t49;
t12 = [-t180 * t428 + m(3) * (t259 * t432 + (t267 ^ 2 + t271 ^ 2) * t270 ^ 2 * qJDD(1)) + m(5) * (t17 * t101 + t288 * t91 + t79 * t50) + m(6) * (t13 * t42 + t59 * t7) + m(4) * (t118 * t185 + t119 * t299 + t175 * t180 + t206 * t310) + t50 * t227 + t101 * t190 + t13 * t182 + m(7) * (t1 * t31 + t10 * t3 + t11 * t4 + t2 * t30) + t4 * t129 + t3 * t130 + t59 * t121 + t31 * t48 + t30 * t47 + t288 * t171 + qJDD(3) * t297 - t185 * t401 + m(2) * qJDD(1) + (-m(6) * t8 + m(7) * t6 + t508) * t58 + (-m(5) * t78 - t485) * t49 - (m(4) * t174 - m(5) * t427 + t493) * t181 + (-m(6) * t41 - t484) * (t101 * t394 - t181 * t415 + t275 * t504) + (-m(5) * t18 + m(6) * t16 + t505) * t100 + (-m(2) - t487) * g(3); t105 * t47 + t106 * t48 + t152 * t121 + t26 * t129 + t27 * t130 + t140 * t227 + t235 * t171 + t73 * t182 + t187 * t190 + t501 * t74 + t505 * t186 + t508 * t151 - t497 * t141 + ((mrSges(4,1) * qJDD(3) - mrSges(4,2) * t282) * t281 + (-mrSges(4,1) * t282 + t228 * t399 - t401) * t277) * t269 + m(7) * (t1 * t106 + t10 * t27 + t105 * t2 + t11 * t26 + t151 * t6 + t37 * t74) + m(6) * (t141 * t71 - t151 * t8 + t152 * t7 + t16 * t186 - t41 * t74 + t42 * t73) + m(3) * t259 + m(5) * (t128 * t359 + t140 * t79 - t141 * t78 + t17 * t187 - t18 * t186 + t235 * t91) + m(4) * (t206 * t273 + (t118 * t277 + t119 * t281 + (-t174 * t277 + t175 * t281) * qJD(3)) * t269) + (-t432 * g(3) + (-g(1) * t430 + g(2) * t431) * t270) * t487; t272 * t518 - t496 * t497 - t382 * t414 / 0.2e1 + t8 * (-mrSges(6,1) * t414 - mrSges(6,3) * t237) + t254 * (Ifges(6,5) * t193 + Ifges(6,3) * t376) / 0.2e1 + (Ifges(6,5) * t237 - Ifges(6,3) * t414) * t453 + t211 * (Ifges(6,4) * t193 + Ifges(6,6) * t376) / 0.2e1 + (Ifges(6,4) * t237 - Ifges(6,6) * t414) * t463 + (Ifges(5,4) * t234 - Ifges(5,2) * t233 + Ifges(5,6) * t258) * t414 / 0.2e1 + (Ifges(5,1) * t234 - Ifges(5,4) * t233 + Ifges(5,5) * t258) * t373 + (-m(5) * t353 - t127 * mrSges(5,1) + t185 * mrSges(4,2) - mrSges(5,3) * t422 - t297 + t524 * (t127 * pkin(4) + pkin(11) * t126 + t353) + t486 * (t127 * t275 - t185 * t415) - t511 * (t127 * t279 + t185 * t417) + t480 * t126) * g(3) + (-m(5) * t354 + mrSges(4,1) * t286 - t87 * mrSges(5,1) + t148 * mrSges(4,2) - mrSges(5,3) * t425 + t524 * (t87 * pkin(4) + pkin(11) * t86 + t354) - t511 * (t148 * t417 + t279 * t87) + t480 * t86 + t486 * (-t148 * t415 + t275 * t87)) * g(1) + (-m(5) * t355 + mrSges(4,1) * t285 - t85 * mrSges(5,1) + t147 * mrSges(4,2) - mrSges(5,3) * t426 + t524 * (t85 * pkin(4) + pkin(11) * t84 + t355) - t511 * (t147 * t417 + t279 * t85) + t480 * t84 + t486 * (-t147 * t415 + t275 * t85)) * g(2) + t41 * (mrSges(6,1) * t376 - mrSges(6,3) * t193) + (t16 * t237 + t193 * t71 - t376 * t42 + t414 * t7) * mrSges(6,2) + (Ifges(7,4) * t124 + Ifges(7,2) * t125) * t461 + t272 * t381 / 0.2e1 + t498 * t121 + (t16 * t220 + t164 * t8 + t41 * t509 + t42 * t510 + t496 * t71 + t498 * t7) * m(6) + (t1 * t195 - t10 * t124 + t11 * t125 + t2 * t327) * mrSges(7,3) + (-Ifges(7,4) * t327 + Ifges(7,2) * t195) * t472 + (-Ifges(7,5) * t327 + Ifges(7,6) * t195) * t465 + (-Ifges(7,1) * t327 + Ifges(7,4) * t195) * t473 + t6 * (-mrSges(7,1) * t195 - mrSges(7,2) * t327) - t327 * t475 + (t17 * t414 - t18 * t416 - t375 * t78 - t376 * t79) * mrSges(5,3) + t258 * (Ifges(5,3) * t272 + (Ifges(5,5) * t276 + Ifges(5,6) * t280) * t268) / 0.2e1 - t171 * t452 + (Ifges(7,5) * t124 + Ifges(7,6) * t125) * t457 + (t268 * t198 + t477 * qJD(3) * (-Ifges(5,2) * t276 + t439)) * t396 / 0.2e1 - t233 * (Ifges(5,6) * t272 + t321) / 0.2e1 - t534 * t227 + (t17 * t239 - t175 * t427 + t18 * t238 - t452 * t91 - t496 * t78 - t534 * t79) * m(5) + t313 * t389 / 0.2e1 - t197 * t376 / 0.2e1 + Ifges(4,3) * qJDD(3) + t238 * t191 + t239 * t190 + (Ifges(7,1) * t124 + Ifges(7,4) * t125) * t459 + t220 * t88 + t193 * t139 / 0.2e1 + t164 * t120 + t154 * t36 + t37 * (-mrSges(7,1) * t125 + mrSges(7,2) * t124) - t118 * mrSges(4,2) + t119 * mrSges(4,1) + t89 * t47 + t90 * t48 + t237 * t474 + t195 * t476 + t124 * t467 + t125 * t469 + t174 * t428 + t91 * t322 + (-t7 * mrSges(6,3) - t67 / 0.2e1 + t16 * mrSges(6,1) + t21 / 0.2e1 - Ifges(6,2) * t463 - Ifges(6,4) * t464 + Ifges(7,3) * t465 + Ifges(7,6) * t472 + Ifges(7,5) * t473 - Ifges(6,6) * t453 + t483) * t236 + (t137 * t373 - t491) * qJD(4) + t493 * t175 + t509 * t183 + t510 * t182 + (-t512 / 0.2e1 - t513 / 0.2e1 + t71 * mrSges(6,1) - Ifges(6,4) * t454 + Ifges(7,3) * t457 + Ifges(7,5) * t459 + Ifges(7,6) * t461 + t478) * t194 + (Ifges(6,1) * t237 - Ifges(6,5) * t414) * t464 + (Ifges(6,1) * t193 + Ifges(6,5) * t376) * t454 + t514 * t117 + t515 * t129 + t516 * t130 + (t1 * t90 + t10 * t516 + t11 * t515 + t154 * t6 + t2 * t89 + t37 * t514) * m(7) - t272 * t519 + t234 * (Ifges(5,5) * t272 + (t276 * Ifges(5,1) + t439) * t268) / 0.2e1; t279 * t523 + (-pkin(4) * t16 + (-t275 * t8 + t446 + (-t275 * t42 - t279 * t41) * qJD(5)) * pkin(11) - t41 * t64 - t42 * t65) * m(6) + t527 * t470 - t8 * t441 - t282 * t313 / 0.2e1 + (Ifges(7,4) * t203 + Ifges(7,2) * t202) * t462 + (-t42 * (-mrSges(6,2) * t276 - t280 * t441) - t41 * (mrSges(6,1) * t276 - mrSges(6,3) * t403)) * t399 + (-t386 - t64) * t183 + (-t60 + t386) * t117 + t508 * pkin(11) * t275 - t22 * t408 / 0.2e1 + (Ifges(7,5) * t203 + Ifges(7,6) * t202) * t458 + (-pkin(11) * t182 + t478) * t395 + (t362 - t227) * t78 + t518 - t519 + t254 * t71 * (mrSges(6,1) * t275 + mrSges(6,2) * t279) + t16 * t346 - ((-Ifges(5,2) * t379 + t279 * t139 + t275 * t94 + t198 + t257) * t280 + t254 * (Ifges(6,3) * t276 + t280 * t337) + t212 * (Ifges(6,5) * t276 + t280 * t343) + t211 * (Ifges(6,6) * t276 + t280 * t340) + t276 * t137) * t399 / 0.2e1 + (-t394 * t41 + t446) * mrSges(6,3) + (Ifges(7,1) * t203 + Ifges(7,4) * t202) * t460 + t274 * t96 * t367 + t381 + (t211 * t340 + t212 * t343 + t254 * t337) * qJD(5) / 0.2e1 + (t481 * t55 + t526 * t54) * g(2) + (t481 * t57 + t526 * t56) * g(1) + (t100 * t526 + t101 * t481) * g(3) + (Ifges(7,5) * t460 + Ifges(7,6) * t462 + Ifges(7,3) * t458 - t525) * t275 * t378 + t279 * t67 / 0.2e1 - t279 * t21 / 0.2e1 + t217 * t47 + t218 * t48 + (t367 * t95 + t368 * t96) * t278 - t65 * t182 - pkin(4) * t88 + (-Ifges(7,5) * t279 + t275 * t342) * t473 + t275 * t474 + t407 * t475 + (Ifges(6,2) * t279 + t437) * t463 + (Ifges(6,1) * t275 + t436) * t464 + (-Ifges(7,3) * t279 + t275 * t336) * t465 + t203 * t468 + (-Ifges(7,6) * t279 + t275 * t339) * t472 + t121 * t450 + (Ifges(6,5) * t275 + Ifges(6,6) * t279) * t453 + (-t335 * t392 + (Ifges(7,3) * t275 + t279 * t336) * qJD(5)) * t457 + (-t341 * t392 + (Ifges(7,5) * t275 + t279 * t342) * qJD(5)) * t459 + (-t338 * t392 + (Ifges(7,6) * t275 + t279 * t339) * qJD(5)) * t461 + t344 * t447 + (t363 + t485) * t79 + (t197 * t373 + t491) * qJD(3) + (t500 * mrSges(7,1) - t499 * mrSges(7,2)) * t37 + (-t1 * t408 + t10 * t499 - t11 * t500 - t2 * t407) * mrSges(7,3) + t139 * t368 + t506 * t130 + t507 * t129 + (t1 * t218 + t2 * t217 + (t37 * t394 + t447) * pkin(11) - t37 * t60 + t507 * t11 + t506 * t10) * m(7) - t279 * t522; t212 * t520 + (-pkin(5) * t6 - t10 * t28 - t11 * t29) * m(7) + t207 * t37 * t344 + t6 * t345 + t382 + (t172 * t339 + t173 * t342 + t207 * t336) * qJD(6) / 0.2e1 - (Ifges(6,1) * t211 - t438 + t94) * t212 / 0.2e1 - (-Ifges(6,2) * t212 + t139 + t205) * t211 / 0.2e1 - t254 * (Ifges(6,5) * t211 - Ifges(6,6) * t212) / 0.2e1 - t71 * (mrSges(6,1) * t212 + mrSges(6,2) * t211) - t29 * t129 - t28 * t130 - pkin(5) * t36 - t7 * mrSges(6,2) + t8 * mrSges(6,1) + t341 * t473 + t274 * t475 + t278 * t476 + t335 * t465 + t391 * t467 + t420 * t468 + t421 * t469 + t393 * t470 + t338 * t472 + t138 * t454 + (Ifges(7,3) * t212 + t211 * t336) * t458 + (Ifges(7,5) * t212 + t211 * t342) * t460 + (Ifges(7,6) * t212 + t211 * t339) * t462 + (t442 + t484) * t42 + ((-t393 + t421) * t11 + (-t391 + t420) * t10 + t492) * mrSges(7,3) + (m(7) * ((-t10 * t278 - t11 * t274) * qJD(6) + t492) - t274 * t47 + t278 * t48 - t129 * t393 - t130 * t391) * pkin(12) + (-t182 + t443) * t41 + (t486 * t59 + t511 * t58) * g(3) + (t486 * t33 - t511 * (t102 * t279 - t275 * t55)) * g(2) + (t486 * t35 - t511 * (t103 * t279 - t275 * t57)) * g(1) - t212 * t521; (Ifges(7,1) * t172 - t435) * t460 + t95 * t459 + (Ifges(7,5) * t172 - Ifges(7,6) * t173) * t458 - t10 * t129 - t37 * (mrSges(7,1) * t173 + mrSges(7,2) * t172) + t11 * t130 - g(1) * ((-t274 * t35 + t278 * t56) * mrSges(7,1) + (-t274 * t56 - t278 * t35) * mrSges(7,2)) - g(2) * ((-t274 * t33 + t278 * t54) * mrSges(7,1) + (-t274 * t54 - t278 * t33) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t30 - mrSges(7,2) * t31) + (t10 * t172 + t11 * t173) * mrSges(7,3) + t21 + (-Ifges(7,2) * t173 + t170 + t96) * t462 + t483;];
tau  = t12;
