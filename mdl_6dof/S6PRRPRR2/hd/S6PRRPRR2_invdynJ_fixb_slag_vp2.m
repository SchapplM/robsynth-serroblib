% Calculate vector of inverse dynamics joint torques for
% S6PRRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:56:00
% EndTime: 2019-03-08 21:56:45
% DurationCPUTime: 27.33s
% Computational Cost: add. (11901->787), mult. (27321->1065), div. (0->0), fcn. (21243->18), ass. (0->359)
t269 = qJ(5) + qJ(6);
t266 = sin(t269);
t267 = cos(t269);
t276 = sin(qJ(5));
t280 = cos(qJ(5));
t316 = -mrSges(6,1) * t280 + mrSges(6,2) * t276;
t417 = t280 * pkin(5);
t465 = -mrSges(5,1) - m(7) * (pkin(4) + t417) - mrSges(7,1) * t267 + mrSges(7,2) * t266 - m(6) * pkin(4) + t316;
t278 = sin(qJ(2));
t272 = sin(pkin(6));
t372 = qJD(1) * t272;
t342 = t278 * t372;
t277 = sin(qJ(3));
t366 = qJD(3) * t277;
t471 = pkin(3) * t366 - t342;
t457 = mrSges(5,2) + m(7) * (-pkin(10) - pkin(9)) - mrSges(7,3) - m(6) * pkin(9);
t270 = sin(pkin(12));
t281 = cos(qJ(3));
t396 = cos(pkin(12));
t228 = t270 * t281 + t277 * t396;
t218 = t228 * qJD(3);
t294 = -t270 * t277 + t281 * t396;
t219 = t294 * qJD(3);
t517 = pkin(4) * t218 - pkin(9) * t219 + t471;
t274 = -qJ(4) - pkin(8);
t331 = qJD(3) * t274;
t210 = qJD(4) * t281 + t277 * t331;
t211 = -qJD(4) * t277 + t281 * t331;
t282 = cos(qJ(2));
t341 = t282 * t372;
t475 = t210 * t396 + t270 * t211 - t294 * t341;
t268 = qJ(3) + pkin(12);
t264 = sin(t268);
t265 = cos(t268);
t516 = (-mrSges(6,3) + t457) * t264 + t465 * t265;
t262 = pkin(3) * t281 + pkin(2);
t153 = -pkin(4) * t294 - pkin(9) * t228 - t262;
t243 = t274 * t277;
t244 = t274 * t281;
t181 = t270 * t243 - t244 * t396;
t363 = qJD(5) * t280;
t364 = qJD(5) * t276;
t483 = t153 * t363 - t181 * t364 + t517 * t276 + t475 * t280;
t515 = -t475 * t276 + t517 * t280;
t496 = -m(6) - m(7);
t468 = m(5) - t496;
t512 = -pkin(3) * t468 - mrSges(4,1);
t169 = t280 * t181;
t387 = t219 * t280;
t511 = -pkin(10) * t387 + pkin(5) * t218 + (-t169 + (pkin(10) * t228 - t153) * t276) * qJD(5) + t515;
t298 = t219 * t276 + t228 * t363;
t510 = pkin(10) * t298 - t483;
t422 = pkin(3) * t270;
t259 = pkin(9) + t422;
t415 = pkin(10) + t259;
t330 = qJD(5) * t415;
t324 = qJD(2) * t396;
t369 = qJD(2) * t277;
t216 = -t270 * t369 + t281 * t324;
t390 = t216 * t276;
t237 = qJD(2) * pkin(8) + t342;
t321 = qJ(4) * qJD(2) + t237;
t273 = cos(pkin(6));
t371 = qJD(1) * t273;
t340 = t277 * t371;
t175 = t281 * t321 + t340;
t162 = t270 * t175;
t256 = t281 * t371;
t174 = -t277 * t321 + t256;
t100 = t174 * t396 - t162;
t367 = qJD(2) * t281;
t217 = -t270 * t367 - t277 * t324;
t356 = pkin(3) * t369;
t130 = -pkin(4) * t217 - pkin(9) * t216 + t356;
t58 = t280 * t100 + t276 * t130;
t509 = -pkin(10) * t390 + t276 * t330 + t58;
t389 = t216 * t280;
t57 = -t100 * t276 + t280 * t130;
t508 = pkin(5) * t217 + pkin(10) * t389 - t280 * t330 - t57;
t505 = t364 - t390;
t279 = cos(qJ(6));
t205 = qJD(5) - t216;
t184 = qJD(3) * t276 - t217 * t280;
t201 = -qJD(2) * t262 + qJD(4) - t341;
t115 = -pkin(4) * t216 + pkin(9) * t217 + t201;
t166 = qJD(3) * pkin(3) + t174;
t326 = t396 * t175;
t91 = t270 * t166 + t326;
t86 = qJD(3) * pkin(9) + t91;
t52 = t280 * t115 - t276 * t86;
t42 = -pkin(10) * t184 + t52;
t32 = pkin(5) * t205 + t42;
t275 = sin(qJ(6));
t183 = qJD(3) * t280 + t217 * t276;
t53 = t115 * t276 + t280 * t86;
t43 = pkin(10) * t183 + t53;
t400 = t275 * t43;
t15 = t279 * t32 - t400;
t399 = t279 * t43;
t16 = t275 * t32 + t399;
t478 = Ifges(5,6) * qJD(3);
t504 = t201 * mrSges(5,1) + t52 * mrSges(6,1) + t15 * mrSges(7,1) - t53 * mrSges(6,2) - t16 * mrSges(7,2) - t478 / 0.2e1;
t479 = Ifges(5,5) * qJD(3);
t503 = -t479 / 0.2e1 - t201 * mrSges(5,2);
t452 = m(7) * pkin(5);
t370 = qJD(2) * t272;
t335 = qJD(1) * t370;
t248 = t282 * t335;
t360 = qJDD(1) * t272;
t204 = t278 * t360 + t248;
t193 = qJDD(2) * pkin(8) + t204;
t500 = qJD(3) * t371 + t193;
t315 = mrSges(6,1) * t276 + mrSges(6,2) * t280;
t499 = -t266 * mrSges(7,1) - t267 * mrSges(7,2) - t276 * t452 - mrSges(5,3) - t315;
t381 = t272 * t278;
t220 = t273 * t281 - t277 * t381;
t322 = t279 * t183 - t184 * t275;
t362 = qJD(2) * qJD(3);
t234 = qJDD(2) * t281 - t277 * t362;
t235 = qJDD(2) * t277 + t281 * t362;
t168 = t270 * t234 + t235 * t396;
t94 = qJD(5) * t183 + qJDD(3) * t276 + t168 * t280;
t95 = -qJD(5) * t184 + qJDD(3) * t280 - t168 * t276;
t28 = qJD(6) * t322 + t275 * t95 + t279 * t94;
t451 = t28 / 0.2e1;
t112 = t183 * t275 + t184 * t279;
t29 = -qJD(6) * t112 - t275 * t94 + t279 * t95;
t450 = t29 / 0.2e1;
t443 = t94 / 0.2e1;
t442 = t95 / 0.2e1;
t167 = t234 * t396 - t270 * t235;
t161 = qJDD(5) - t167;
t155 = qJDD(6) + t161;
t437 = t155 / 0.2e1;
t436 = t161 / 0.2e1;
t495 = t234 / 0.2e1;
t385 = t228 * t280;
t87 = t280 * t153 - t181 * t276;
t65 = -pkin(5) * t294 - pkin(10) * t385 + t87;
t386 = t228 * t276;
t88 = t276 * t153 + t169;
t70 = -pkin(10) * t386 + t88;
t31 = t275 * t65 + t279 * t70;
t494 = -qJD(6) * t31 + t275 * t510 + t279 * t511;
t30 = -t275 * t70 + t279 * t65;
t493 = qJD(6) * t30 + t275 * t511 - t279 * t510;
t225 = t415 * t276;
t226 = t415 * t280;
t150 = -t225 * t275 + t226 * t279;
t492 = -qJD(6) * t150 + t275 * t509 + t279 * t508;
t149 = -t225 * t279 - t226 * t275;
t491 = qJD(6) * t149 + t275 * t508 - t279 * t509;
t200 = qJD(6) + t205;
t487 = t205 * Ifges(6,3);
t489 = t183 * Ifges(6,6);
t490 = t184 * Ifges(6,5) + t112 * Ifges(7,5) + Ifges(7,6) * t322 + t200 * Ifges(7,3) + t487 + t489;
t202 = Ifges(5,4) * t216;
t486 = t216 * Ifges(5,2);
t398 = qJDD(3) / 0.2e1;
t484 = -qJD(5) * t88 + t515;
t482 = -t452 - mrSges(6,1);
t146 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t168;
t48 = -mrSges(6,1) * t95 + mrSges(6,2) * t94;
t481 = t48 - t146;
t98 = t174 * t270 + t326;
t480 = pkin(5) * t505 - t98;
t304 = t275 * t276 - t279 * t280;
t140 = t304 * t228;
t126 = t304 * t216;
t467 = qJD(5) + qJD(6);
t170 = t467 * t304;
t477 = -t170 + t126;
t231 = t275 * t280 + t276 * t279;
t125 = t231 * t216;
t171 = t467 * t231;
t476 = -t171 + t125;
t413 = mrSges(5,3) * t217;
t474 = -qJD(3) * mrSges(5,1) - mrSges(6,1) * t183 + mrSges(6,2) * t184 - t413;
t359 = qJDD(1) * t273;
t107 = -t237 * t366 + t277 * t359 + t281 * t500;
t187 = t237 * t281 + t340;
t252 = t281 * t359;
t108 = -qJD(3) * t187 - t193 * t277 + t252;
t470 = t107 * t281 - t108 * t277;
t361 = qJD(2) * qJD(4);
t365 = qJD(3) * t281;
t76 = -t237 * t365 + qJDD(3) * pkin(3) - qJ(4) * t235 + t252 + (-t361 - t500) * t277;
t81 = qJ(4) * t234 + t281 * t361 + t107;
t39 = t270 * t76 + t396 * t81;
t37 = qJDD(3) * pkin(9) + t39;
t247 = t278 * t335;
t203 = t282 * t360 - t247;
t192 = -qJDD(2) * pkin(2) - t203;
t151 = -pkin(3) * t234 + qJDD(4) + t192;
t66 = -pkin(4) * t167 - pkin(9) * t168 + t151;
t10 = t115 * t363 + t276 * t66 + t280 * t37 - t364 * t86;
t11 = -qJD(5) * t53 - t276 * t37 + t280 * t66;
t469 = t10 * t280 - t11 * t276;
t466 = 0.2e1 * t398;
t90 = t166 * t396 - t162;
t85 = -qJD(3) * pkin(4) - t90;
t67 = -t183 * pkin(5) + t85;
t464 = -mrSges(7,1) * t67 + mrSges(7,3) * t16;
t463 = mrSges(7,2) * t67 - t15 * mrSges(7,3);
t462 = m(5) * t90 - t474;
t6 = pkin(5) * t161 - pkin(10) * t94 + t11;
t7 = pkin(10) * t95 + t10;
t2 = qJD(6) * t15 + t275 * t6 + t279 * t7;
t3 = -qJD(6) * t16 - t275 * t7 + t279 * t6;
t461 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t460 = t11 * mrSges(6,1) - t10 * mrSges(6,2);
t352 = m(4) * pkin(8) + mrSges(4,3);
t459 = mrSges(3,2) - t352 + t499;
t318 = -mrSges(4,1) * t281 + mrSges(4,2) * t277;
t293 = m(4) * pkin(2) - t318;
t458 = mrSges(3,1) + t293 - t516;
t284 = qJD(2) ^ 2;
t455 = Ifges(7,4) * t451 + Ifges(7,2) * t450 + Ifges(7,6) * t437;
t454 = Ifges(7,1) * t451 + Ifges(7,4) * t450 + Ifges(7,5) * t437;
t453 = m(5) * pkin(3);
t449 = Ifges(6,1) * t443 + Ifges(6,4) * t442 + Ifges(6,5) * t436;
t405 = Ifges(7,4) * t112;
t46 = Ifges(7,2) * t322 + Ifges(7,6) * t200 + t405;
t448 = -t46 / 0.2e1;
t447 = t46 / 0.2e1;
t105 = Ifges(7,4) * t322;
t47 = Ifges(7,1) * t112 + Ifges(7,5) * t200 + t105;
t446 = -t47 / 0.2e1;
t445 = t47 / 0.2e1;
t441 = -t322 / 0.2e1;
t440 = t322 / 0.2e1;
t439 = -t112 / 0.2e1;
t438 = t112 / 0.2e1;
t435 = -t183 / 0.2e1;
t434 = -t184 / 0.2e1;
t433 = t184 / 0.2e1;
t432 = -t200 / 0.2e1;
t431 = t200 / 0.2e1;
t430 = -t205 / 0.2e1;
t429 = -t216 / 0.2e1;
t428 = -t217 / 0.2e1;
t424 = t280 / 0.2e1;
t421 = pkin(5) * t184;
t419 = g(3) * t272;
t414 = mrSges(5,3) * t216;
t412 = mrSges(6,3) * t183;
t411 = mrSges(6,3) * t184;
t410 = Ifges(4,4) * t277;
t409 = Ifges(4,4) * t281;
t408 = Ifges(5,4) * t217;
t407 = Ifges(6,4) * t276;
t406 = Ifges(6,4) * t280;
t401 = t184 * Ifges(6,4);
t397 = cos(pkin(11));
t271 = sin(pkin(11));
t384 = t271 * t272;
t383 = t271 * t278;
t382 = t271 * t282;
t380 = t272 * t281;
t379 = t272 * t282;
t328 = t397 * t278;
t213 = t273 * t328 + t382;
t329 = t272 * t397;
t157 = t213 * t265 - t264 * t329;
t327 = t397 * t282;
t212 = -t273 * t327 + t383;
t377 = (-t157 * t266 + t212 * t267) * mrSges(7,1) + (-t157 * t267 - t212 * t266) * mrSges(7,2);
t215 = -t273 * t383 + t327;
t159 = t215 * t265 + t264 * t384;
t214 = t273 * t382 + t328;
t376 = (-t159 * t266 + t214 * t267) * mrSges(7,1) + (-t159 * t267 - t214 * t266) * mrSges(7,2);
t196 = t264 * t273 + t265 * t381;
t375 = (-t196 * t266 - t267 * t379) * mrSges(7,1) + (-t196 * t267 + t266 * t379) * mrSges(7,2);
t368 = qJD(2) * t278;
t358 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t155;
t357 = Ifges(6,5) * t94 + Ifges(6,6) * t95 + Ifges(6,3) * t161;
t56 = -mrSges(7,1) * t322 + mrSges(7,2) * t112;
t353 = t56 + t474;
t350 = mrSges(4,3) * t369;
t349 = mrSges(4,3) * t367;
t347 = t276 * t379;
t345 = t280 * t379;
t182 = Ifges(6,4) * t183;
t84 = t184 * Ifges(6,1) + t205 * Ifges(6,5) + t182;
t344 = t84 * t424;
t343 = t396 * pkin(3);
t339 = t272 * t368;
t338 = t282 * t370;
t332 = -t364 / 0.2e1;
t96 = -t167 * mrSges(5,1) + t168 * mrSges(5,2);
t128 = t210 * t270 - t396 * t211;
t180 = -t396 * t243 - t244 * t270;
t260 = -t343 - pkin(4);
t38 = -t270 * t81 + t396 * t76;
t313 = Ifges(6,1) * t280 - t407;
t312 = t281 * Ifges(4,2) + t410;
t311 = -Ifges(6,2) * t276 + t406;
t310 = Ifges(4,5) * t281 - Ifges(4,6) * t277;
t309 = Ifges(6,5) * t280 - Ifges(6,6) * t276;
t221 = t273 * t277 + t278 * t380;
t136 = t270 * t220 + t221 * t396;
t116 = -t136 * t276 - t345;
t301 = -t136 * t280 + t347;
t59 = t116 * t279 + t275 * t301;
t60 = t116 * t275 - t279 * t301;
t118 = -mrSges(6,2) * t205 + t412;
t119 = mrSges(6,1) * t205 - t411;
t307 = t118 * t280 - t119 * t276;
t302 = t358 + t461;
t299 = t85 * t315;
t297 = t228 * t364 - t387;
t238 = -qJD(2) * pkin(2) - t341;
t296 = t238 * (mrSges(4,1) * t277 + mrSges(4,2) * t281);
t295 = t277 * (Ifges(4,1) * t281 - t410);
t36 = -qJDD(3) * pkin(4) - t38;
t288 = -g(1) * t214 - g(2) * t212 + g(3) * t379;
t263 = Ifges(4,4) * t367;
t242 = -qJD(3) * mrSges(4,2) + t349;
t241 = qJD(3) * mrSges(4,1) - t350;
t239 = t260 - t417;
t232 = t318 * qJD(2);
t223 = Ifges(4,1) * t369 + Ifges(4,5) * qJD(3) + t263;
t222 = Ifges(4,6) * qJD(3) + qJD(2) * t312;
t207 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t235;
t206 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t234;
t190 = -qJD(3) * mrSges(5,2) + t414;
t186 = -t237 * t277 + t256;
t177 = -mrSges(4,1) * t234 + mrSges(4,2) * t235;
t173 = qJD(3) * t220 + t281 * t338;
t172 = -qJD(3) * t221 - t277 * t338;
t147 = -mrSges(5,1) * t216 - mrSges(5,2) * t217;
t145 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t167;
t139 = t231 * t228;
t138 = -t217 * Ifges(5,1) + t202 + t479;
t137 = -t408 + t478 + t486;
t135 = -t220 * t396 + t221 * t270;
t127 = pkin(5) * t386 + t180;
t99 = t270 * t172 + t173 * t396;
t97 = -t172 * t396 + t173 * t270;
t83 = t183 * Ifges(6,2) + t205 * Ifges(6,6) + t401;
t78 = mrSges(7,1) * t200 - mrSges(7,3) * t112;
t77 = -mrSges(7,2) * t200 + mrSges(7,3) * t322;
t75 = pkin(5) * t298 + t128;
t64 = -mrSges(6,2) * t161 + mrSges(6,3) * t95;
t63 = mrSges(6,1) * t161 - mrSges(6,3) * t94;
t55 = t140 * t467 - t231 * t219;
t54 = -t171 * t228 - t219 * t304;
t50 = qJD(5) * t301 - t276 * t99 + t280 * t339;
t49 = qJD(5) * t116 + t276 * t339 + t280 * t99;
t33 = t94 * Ifges(6,4) + t95 * Ifges(6,2) + t161 * Ifges(6,6);
t23 = -t95 * pkin(5) + t36;
t22 = -mrSges(7,2) * t155 + mrSges(7,3) * t29;
t21 = mrSges(7,1) * t155 - mrSges(7,3) * t28;
t18 = t279 * t42 - t400;
t17 = -t275 * t42 - t399;
t14 = -qJD(6) * t60 - t275 * t49 + t279 * t50;
t13 = qJD(6) * t59 + t275 * t50 + t279 * t49;
t12 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t1 = [m(2) * qJDD(1) + t116 * t63 - t301 * t64 + t49 * t118 + t50 * t119 + t13 * t77 + t136 * t145 + t14 * t78 + t172 * t241 + t173 * t242 + t99 * t190 + t221 * t206 + t220 * t207 + t59 * t21 + t60 * t22 + t353 * t97 + (t12 + t481) * t135 + (-m(2) - m(3) - m(4) - t468) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t284 - t177 - t96) * t282 + (-mrSges(3,1) * t284 - mrSges(3,2) * qJDD(2) + (t147 + t232) * qJD(2)) * t278) * t272 + m(4) * (t107 * t221 + t108 * t220 + t172 * t186 + t173 * t187 + (-t192 * t282 + t238 * t368) * t272) + m(5) * (-t135 * t38 + t136 * t39 - t90 * t97 + t91 * t99 + (-t151 * t282 + t201 * t368) * t272) + m(3) * (qJDD(1) * t273 ^ 2 + (t203 * t282 + t204 * t278) * t272) + m(7) * (t13 * t16 + t135 * t23 + t14 * t15 + t2 * t60 + t3 * t59 + t67 * t97) + m(6) * (-t10 * t301 + t11 * t116 + t135 * t36 + t49 * t53 + t50 * t52 + t85 * t97); t205 * (-Ifges(6,5) * t297 - Ifges(6,6) * t298) / 0.2e1 + t183 * (-Ifges(6,4) * t297 - Ifges(6,2) * t298) / 0.2e1 + ((t516 * t282 + (t274 * t468 + t499) * t278) * t272 - t468 * t262 * t379) * g(3) + (t296 + t310 * qJD(3) / 0.2e1) * qJD(3) + (-t139 * t2 + t140 * t3 - t15 * t54 + t16 * t55) * mrSges(7,3) + (-Ifges(7,4) * t140 - Ifges(7,2) * t139) * t450 - (t358 + t357) * t294 / 0.2e1 - (t151 * mrSges(5,1) - t39 * mrSges(5,3) - Ifges(5,4) * t168 + Ifges(6,5) * t443 + Ifges(7,5) * t451 - Ifges(5,2) * t167 - Ifges(5,6) * t466 + Ifges(6,6) * t442 + Ifges(7,6) * t450 + Ifges(6,3) * t436 + Ifges(7,3) * t437 + t460 + t461) * t294 + (-Ifges(6,1) * t297 - Ifges(6,4) * t298) * t433 + t235 * t409 / 0.2e1 + (-t10 * t386 - t11 * t385 + t297 * t52 - t298 * t53) * mrSges(6,3) + (-pkin(2) * t192 - (t238 * t278 + (-t186 * t277 + t187 * t281) * t282) * t372) * m(4) + t312 * t495 + Ifges(4,6) * t281 * t398 + (Ifges(7,4) * t54 + Ifges(7,2) * t55) * t440 + (-Ifges(7,1) * t140 - Ifges(7,4) * t139) * t451 + t23 * (mrSges(7,1) * t139 - mrSges(7,2) * t140) + (-Ifges(7,5) * t140 - Ifges(7,6) * t139) * t437 + (t295 + t281 * (-Ifges(4,2) * t277 + t409)) * t362 / 0.2e1 - t140 * t454 - t139 * t455 + (-t282 * t419 + t203 + t247) * mrSges(3,1) + (t278 * t419 - t204 + t248) * mrSges(3,2) - (t278 * t352 + t282 * t293) * t419 + (Ifges(7,1) * t54 + Ifges(7,4) * t55) * t438 + ((-m(6) * t85 - m(7) * t67 + t462 - t56) * t341 + t151 * mrSges(5,2) - t38 * mrSges(5,3) + Ifges(5,1) * t168 + Ifges(5,4) * t167 + Ifges(5,5) * t466 + t309 * t436 + t311 * t442 + t313 * t443 + t36 * t315 + t84 * t332) * t228 - t33 * t386 / 0.2e1 - t222 * t366 / 0.2e1 + t223 * t365 / 0.2e1 + Ifges(3,3) * qJDD(2) - t232 * t342 + t281 * (Ifges(4,4) * t235 + Ifges(4,2) * t234 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t262 * t96 + t181 * t145 - pkin(2) * t177 + t54 * t445 + t55 * t447 + t385 * t449 + t127 * t12 + t493 * t77 + t494 * t78 + (t127 * t23 + t15 * t494 + t16 * t493 + t2 * t31 + t3 * t30 + t67 * t75) * m(7) + (Ifges(4,1) * t235 + Ifges(4,4) * t495 + Ifges(4,5) * t466 + t241 * t341) * t277 + t87 * t63 + t88 * t64 + t483 * t118 + t484 * t119 + (t10 * t88 + t11 * t87 + t128 * t85 + t180 * t36 + t483 * t53 + t484 * t52) * m(6) + t75 * t56 + t67 * (-mrSges(7,1) * t55 + mrSges(7,2) * t54) + t481 * t180 + t474 * t128 + t475 * t190 + (-t128 * t90 - t151 * t262 - t180 * t38 + t181 * t39 + t201 * t471 + t475 * t91) * m(5) + (m(4) * ((-t186 * t281 - t187 * t277) * qJD(3) + t470) + t281 * t206 - t277 * t207 - t241 * t365 - t242 * t366) * pkin(8) + (-t186 * t365 - t187 * t366 + t470) * mrSges(4,3) + t471 * t147 - t298 * t83 / 0.2e1 + (-t468 * (-t214 * t262 - t215 * t274) + t459 * t215 + t458 * t214) * g(1) + (-t468 * (-t212 * t262 - t213 * t274) + t459 * t213 + t458 * t212) * g(2) + t30 * t21 + t31 * t22 + (Ifges(7,5) * t54 + Ifges(7,6) * t55) * t431 + (-t90 * mrSges(5,3) + t202 / 0.2e1 + t138 / 0.2e1 + Ifges(5,1) * t428 + t344 - t503) * t219 + (-t91 * mrSges(5,3) + t490 / 0.2e1 - t137 / 0.2e1 - t486 / 0.2e1 + Ifges(7,6) * t440 + Ifges(7,5) * t438 - Ifges(5,4) * t428 + Ifges(7,3) * t431 + Ifges(6,5) * t433 + t489 / 0.2e1 + t487 / 0.2e1 + t504) * t218 - t281 * t242 * t341 + t85 * (mrSges(6,1) * t298 - mrSges(6,2) * t297) + t192 * t318; (t202 + t138) * t429 + (-Ifges(7,4) * t126 - Ifges(7,2) * t125) * t441 + (-Ifges(7,1) * t126 - Ifges(7,4) * t125) * t439 + (-Ifges(7,5) * t126 - Ifges(7,6) * t125) * t432 - t67 * (mrSges(7,1) * t125 - mrSges(7,2) * t126) + (t344 + t299) * qJD(5) + (-(-t213 * t281 + t277 * t329) * mrSges(4,2) + t457 * t157 + t465 * (-t213 * t264 - t265 * t329) + (pkin(3) * t496 - mrSges(4,1) - t453) * (-t213 * t277 - t281 * t329)) * g(2) + (t125 * t16 - t126 * t15 - t2 * t304 - t231 * t3) * mrSges(7,3) + (Ifges(7,4) * t231 - Ifges(7,2) * t304) * t450 + (Ifges(7,1) * t231 - Ifges(7,4) * t304) * t451 + t23 * (mrSges(7,1) * t304 + mrSges(7,2) * t231) + (Ifges(7,5) * t231 - Ifges(7,6) * t304) * t437 - t304 * t455 + (t390 / 0.2e1 + t332) * t83 - m(5) * (t100 * t91 + t201 * t356) + (t260 * t36 - t52 * t57 - t53 * t58 - t85 * t98) * m(6) - (-Ifges(4,2) * t369 + t223 + t263) * t367 / 0.2e1 + (t183 * t311 + t184 * t313 + t205 * t309) * qJD(5) / 0.2e1 + (t270 * t39 + t38 * t396) * t453 + t231 * t454 - t284 * t295 / 0.2e1 + (t241 + t350) * t187 - t91 * t413 + (-t242 + t349) * t186 - t84 * t389 / 0.2e1 + t222 * t369 / 0.2e1 - t310 * t362 / 0.2e1 - t147 * t356 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + t260 * t48 + Ifges(4,5) * t235 + t239 * t12 + Ifges(4,6) * t234 - t100 * t190 + Ifges(5,5) * t168 + Ifges(5,6) * t167 + t149 * t21 + t150 * t22 + (Ifges(6,2) * t280 + t407) * t442 + (Ifges(6,1) * t276 + t406) * t443 - t126 * t446 - t125 * t448 + (Ifges(6,5) * t276 + Ifges(6,6) * t280) * t436 + t33 * t424 + t137 * t428 + t90 * t414 + t145 * t422 - t57 * t119 - t58 * t118 + (Ifges(5,1) * t216 + t408 + t490) * t217 / 0.2e1 + t491 * t77 + t492 * t78 + (t149 * t3 + t15 * t492 + t150 * t2 + t16 * t491 + t23 * t239 + t480 * t67) * m(7) - t107 * mrSges(4,2) + t108 * mrSges(4,1) + t480 * t56 + (m(6) * ((-t276 * t53 - t280 * t52) * qJD(5) + t469) + t280 * t64 - t276 * t63 - t119 * t363 - t118 * t364) * t259 + t38 * mrSges(5,1) - t39 * mrSges(5,2) - (Ifges(7,1) * t438 + Ifges(7,4) * t440 + Ifges(7,5) * t431 + t445 + t463) * t170 - (Ifges(7,4) * t438 + Ifges(7,2) * t440 + Ifges(7,6) * t431 + t447 + t464) * t171 + t462 * t98 + (t309 * t430 + t311 * t435 + t313 * t434 - t299 + t503) * t216 + (-Ifges(6,5) * t434 - Ifges(7,5) * t439 + Ifges(5,2) * t429 - Ifges(6,6) * t435 - Ifges(7,6) * t441 - Ifges(6,3) * t430 - Ifges(7,3) * t432 + t504) * t217 + (-g(1) * t159 - g(2) * t157 - g(3) * t196 - t505 * t53 + (-t363 + t389) * t52 + t469) * mrSges(6,3) + (-(-t215 * t281 - t277 * t384) * mrSges(4,2) + t457 * t159 + t465 * (-t215 * t264 + t265 * t384) + t512 * (-t215 * t277 + t271 * t380)) * g(1) + (mrSges(4,2) * t221 + t457 * t196 + t465 * (-t264 * t381 + t265 * t273) + t512 * t220) * g(3) + t146 * t343 + t276 * t449 - qJD(2) * t296 + t36 * t316; -t304 * t21 + t231 * t22 + t276 * t64 + t280 * t63 + t476 * t78 + t477 * t77 + t307 * qJD(5) + t353 * t217 + (-t190 - t307) * t216 + t96 + (t15 * t476 + t16 * t477 + t2 * t231 + t217 * t67 - t3 * t304 + t288) * m(7) + (t10 * t276 + t11 * t280 + t217 * t85 + t288 + t205 * (-t276 * t52 + t280 * t53)) * m(6) + (-t216 * t91 - t217 * t90 + t151 + t288) * m(5); -(Ifges(7,4) * t439 + Ifges(7,2) * t441 + Ifges(7,6) * t432 + t448 - t464) * t112 + (Ifges(7,1) * t439 + Ifges(7,4) * t441 + Ifges(7,5) * t432 + t446 - t463) * t322 + (-t375 - (-t196 * t280 + t347) * mrSges(6,2) + t482 * (-t196 * t276 - t345)) * g(3) + (-t376 - (-t159 * t280 - t214 * t276) * mrSges(6,2) + t482 * (-t159 * t276 + t214 * t280)) * g(1) + (-t377 - (-t157 * t280 - t212 * t276) * mrSges(6,2) + t482 * (-t157 * t276 + t212 * t280)) * g(2) + (t2 * t275 + t279 * t3 + (-t15 * t275 + t16 * t279) * qJD(6)) * t452 + (-Ifges(6,2) * t184 + t182 + t84) * t435 - t56 * t421 - m(7) * (t15 * t17 + t16 * t18 + t421 * t67) + t357 + t302 + t460 - t85 * (mrSges(6,1) * t184 + mrSges(6,2) * t183) + (Ifges(6,1) * t183 - t401) * t434 + (Ifges(6,5) * t183 - Ifges(6,6) * t184) * t430 + t83 * t433 - t18 * t77 - t17 * t78 + (t412 - t118) * t52 + (t411 + t119) * t53 + ((-t275 * t78 + t279 * t77) * qJD(6) + t21 * t279 + t22 * t275) * pkin(5); -t67 * (mrSges(7,1) * t112 + mrSges(7,2) * t322) + (Ifges(7,1) * t322 - t405) * t439 + t46 * t438 + (Ifges(7,5) * t322 - Ifges(7,6) * t112) * t432 - t15 * t77 + t16 * t78 - g(1) * t376 - g(2) * t377 - g(3) * t375 + (t112 * t16 + t15 * t322) * mrSges(7,3) + t302 + (-Ifges(7,2) * t112 + t105 + t47) * t441;];
tau  = t1;
