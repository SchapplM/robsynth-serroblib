% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:21
% EndTime: 2019-03-09 11:48:17
% DurationCPUTime: 32.69s
% Computational Cost: add. (15695->815), mult. (35732->1036), div. (0->0), fcn. (26556->14), ass. (0->365)
t535 = -mrSges(7,1) - mrSges(6,1);
t548 = mrSges(6,2) + mrSges(7,2);
t532 = Ifges(6,4) + Ifges(7,4);
t300 = sin(qJ(2));
t304 = cos(qJ(2));
t405 = sin(pkin(10));
t340 = qJD(1) * t405;
t406 = cos(pkin(10));
t341 = qJD(1) * t406;
t241 = -t300 * t340 + t304 * t341;
t242 = -t300 * t341 - t304 * t340;
t376 = qJD(1) * t300;
t362 = pkin(2) * t376;
t156 = -pkin(3) * t242 - pkin(8) * t241 + t362;
t297 = -qJ(3) - pkin(7);
t274 = t297 * t304;
t262 = qJD(1) * t274;
t245 = t405 * t262;
t272 = t297 * t300;
t261 = qJD(1) * t272;
t189 = t261 * t406 + t245;
t299 = sin(qJ(4));
t303 = cos(qJ(4));
t109 = t299 * t156 + t303 * t189;
t352 = t405 * pkin(2);
t280 = t352 + pkin(8);
t430 = pkin(9) + t280;
t345 = qJD(4) * t430;
t388 = t299 * t241;
t551 = -pkin(9) * t388 + t299 * t345 + t109;
t108 = t303 * t156 - t189 * t299;
t402 = t241 * t303;
t550 = pkin(4) * t242 + pkin(9) * t402 - t303 * t345 - t108;
t371 = qJD(4) * t299;
t549 = t371 - t388;
t533 = Ifges(6,1) + Ifges(7,1);
t531 = -Ifges(7,5) - Ifges(6,5);
t530 = Ifges(6,2) + Ifges(7,2);
t529 = Ifges(6,6) + Ifges(7,6);
t528 = Ifges(7,3) + Ifges(6,3);
t298 = sin(qJ(5));
t302 = cos(qJ(5));
t260 = t298 * t303 + t299 * t302;
t151 = t260 * t241;
t500 = qJD(4) + qJD(5);
t194 = t500 * t260;
t509 = t151 - t194;
t321 = t298 * t299 - t302 * t303;
t152 = t321 * t241;
t193 = t500 * t321;
t510 = -t193 + t152;
t496 = m(4) + m(5) + m(6) + m(7);
t296 = qJ(4) + qJ(5);
t291 = cos(t296);
t292 = t303 * pkin(4);
t267 = pkin(5) * t291 + t292;
t263 = pkin(3) + t267;
t284 = t292 + pkin(3);
t290 = sin(t296);
t333 = -mrSges(5,1) * t303 + mrSges(5,2) * t299;
t547 = -m(5) * pkin(3) - m(6) * t284 - m(7) * t263 + t548 * t290 + t291 * t535 + t333;
t546 = -mrSges(4,2) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t199 = qJD(2) * t303 + t242 * t299;
t200 = qJD(2) * t299 - t242 * t303;
t136 = t199 * t298 + t200 * t302;
t231 = qJD(4) - t241;
t223 = qJD(5) + t231;
t338 = t302 * t199 - t200 * t298;
t541 = t532 * t338;
t523 = t136 * t533 - t223 * t531 + t541;
t545 = t523 / 0.2e1;
t539 = t532 * t136;
t524 = t223 * t529 + t338 * t530 + t539;
t544 = t524 / 0.2e1;
t367 = qJD(1) * qJD(2);
t349 = t300 * t367;
t366 = qJDD(1) * t304;
t264 = -t349 + t366;
t265 = qJDD(1) * t300 + t304 * t367;
t190 = t264 * t406 - t405 * t265;
t187 = qJDD(4) - t190;
t184 = qJDD(5) + t187;
t191 = t264 * t405 + t265 * t406;
t128 = qJD(4) * t199 + qJDD(2) * t299 + t191 * t303;
t129 = -qJD(4) * t200 + qJDD(2) * t303 - t191 * t299;
t55 = qJD(5) * t338 + t128 * t302 + t129 * t298;
t56 = -qJD(5) * t136 - t128 * t298 + t129 * t302;
t543 = -t530 * t56 / 0.2e1 - t532 * t55 / 0.2e1 - t529 * t184 / 0.2e1;
t252 = t430 * t299;
t253 = t430 * t303;
t368 = qJD(5) * t302;
t369 = qJD(5) * t298;
t517 = -t252 * t368 - t253 * t369 + t298 * t550 - t302 * t551;
t178 = -t298 * t252 + t302 * t253;
t515 = -qJD(5) * t178 + t298 * t551 + t302 * t550;
t542 = qJ(6) * t136;
t273 = -mrSges(3,1) * t304 + mrSges(3,2) * t300;
t295 = qJ(2) + pkin(10);
t288 = sin(t295);
t289 = cos(t295);
t540 = -mrSges(4,1) * t289 - t288 * t546 + t273;
t256 = t300 * t405 - t304 * t406;
t244 = t256 * qJD(2);
t257 = t300 * t406 + t304 * t405;
t370 = qJD(4) * t303;
t319 = -t244 * t299 + t257 * t370;
t344 = t406 * t262;
t188 = t261 * t405 - t344;
t507 = pkin(4) * t549 - t188;
t301 = sin(qJ(1));
t305 = cos(qJ(1));
t538 = g(1) * t305 + g(2) * t301;
t293 = t304 * pkin(2);
t285 = t293 + pkin(1);
t268 = -qJD(1) * t285 + qJD(3);
t145 = -pkin(3) * t241 + pkin(8) * t242 + t268;
t251 = qJD(2) * pkin(2) + t261;
t182 = t405 * t251 - t344;
t167 = qJD(2) * pkin(8) + t182;
t98 = t303 * t145 - t167 * t299;
t85 = -pkin(9) * t200 + t98;
t73 = pkin(4) * t231 + t85;
t99 = t145 * t299 + t167 * t303;
t86 = pkin(9) * t199 + t99;
t83 = t302 * t86;
t33 = t298 * t73 + t83;
t511 = qJ(6) * t338;
t25 = t33 + t511;
t181 = t251 * t406 + t245;
t166 = -qJD(2) * pkin(3) - t181;
t130 = -t199 * pkin(4) + t166;
t78 = -pkin(5) * t338 + qJD(6) + t130;
t537 = -t78 * mrSges(7,1) + mrSges(6,3) * t33 + mrSges(7,3) * t25 + t544;
t484 = m(6) * pkin(4);
t476 = t128 / 0.2e1;
t475 = t129 / 0.2e1;
t463 = t187 / 0.2e1;
t536 = t264 / 0.2e1;
t473 = -t338 / 0.2e1;
t432 = qJD(2) / 0.2e1;
t526 = -t184 * t531 + t532 * t56 + t533 * t55;
t39 = -mrSges(7,2) * t184 + mrSges(7,3) * t56;
t40 = -mrSges(6,2) * t184 + mrSges(6,3) * t56;
t525 = t39 + t40;
t522 = qJ(6) * t509 - qJD(6) * t321 + t517;
t521 = pkin(5) * t242 - qJ(6) * t510 - qJD(6) * t260 + t515;
t227 = Ifges(4,4) * t241;
t520 = t199 * Ifges(5,6);
t519 = t231 * Ifges(5,3);
t518 = t241 * Ifges(4,2);
t407 = qJDD(2) / 0.2e1;
t180 = pkin(3) * t256 - pkin(8) * t257 - t285;
t197 = t272 * t405 - t274 * t406;
t192 = t303 * t197;
t125 = t299 * t180 + t192;
t397 = t257 * t299;
t103 = -pkin(9) * t397 + t125;
t124 = t303 * t180 - t197 * t299;
t396 = t257 * t303;
t95 = pkin(4) * t256 - pkin(9) * t396 + t124;
t58 = t302 * t103 + t298 * t95;
t516 = -pkin(5) * t509 + t507;
t514 = t484 + mrSges(5,1);
t513 = Ifges(4,5) * qJD(2);
t512 = Ifges(4,6) * qJD(2);
t163 = t321 * t257;
t508 = qJD(2) * mrSges(4,1) + mrSges(5,1) * t199 - mrSges(5,2) * t200 + mrSges(4,3) * t242;
t506 = t184 * t528 + t529 * t56 - t531 * t55;
t390 = t291 * t301;
t391 = t290 * t305;
t214 = -t289 * t391 + t390;
t389 = t291 * t305;
t392 = t290 * t301;
t215 = t289 * t389 + t392;
t505 = t214 * t535 + t548 * t215;
t212 = t289 * t392 + t389;
t213 = -t289 * t390 + t391;
t504 = -t212 * t535 - t548 * t213;
t286 = pkin(7) * t366;
t254 = -pkin(7) * t349 + t286;
t255 = t265 * pkin(7);
t503 = t254 * t304 + t255 * t300;
t404 = qJDD(1) * pkin(1);
t226 = -pkin(2) * t264 + qJDD(3) - t404;
t106 = -pkin(3) * t190 - pkin(8) * t191 + t226;
t373 = qJD(3) * t300;
t173 = qJDD(2) * pkin(2) - qJ(3) * t265 - qJD(1) * t373 - t255;
t374 = qJD(2) * t300;
t357 = pkin(7) * t374;
t372 = qJD(3) * t304;
t185 = qJ(3) * t264 + t286 + (-t357 + t372) * qJD(1);
t120 = t405 * t173 + t406 * t185;
t112 = qJDD(2) * pkin(8) + t120;
t30 = t299 * t106 + t303 * t112 + t145 * t370 - t167 * t371;
t31 = -qJD(4) * t99 + t303 * t106 - t112 * t299;
t502 = -t299 * t31 + t30 * t303;
t501 = mrSges(6,1) * t290 + t291 * t548;
t499 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t498 = t200 * Ifges(5,5) - t136 * t531 + t223 * t528 + t338 * t529 + t519 + t520;
t497 = 0.2e1 * t407;
t495 = -t268 * mrSges(4,2) + t181 * mrSges(4,3);
t81 = t298 * t86;
t32 = t302 * t73 - t81;
t24 = t32 - t542;
t21 = pkin(5) * t223 + t24;
t492 = t78 * mrSges(7,2) - mrSges(6,3) * t32 - mrSges(7,3) * t21;
t440 = pkin(4) * t299;
t266 = pkin(5) * t290 + t440;
t491 = -m(6) * t440 - m(7) * t266;
t490 = t31 * mrSges(5,1) - t30 * mrSges(5,2);
t488 = m(3) * pkin(1) + mrSges(2,1) - t540;
t19 = pkin(4) * t187 - pkin(9) * t128 + t31;
t23 = pkin(9) * t129 + t30;
t6 = -qJD(5) * t33 + t302 * t19 - t23 * t298;
t2 = pkin(5) * t184 - qJ(6) * t55 - qJD(6) * t136 + t6;
t5 = t298 * t19 + t302 * t23 + t73 * t368 - t369 * t86;
t3 = qJ(6) * t56 + qJD(6) * t338 + t5;
t487 = t6 * mrSges(6,1) + t2 * mrSges(7,1) - t5 * mrSges(6,2) - t3 * mrSges(7,2);
t486 = t268 * mrSges(4,1) + t98 * mrSges(5,1) + t32 * mrSges(6,1) + t21 * mrSges(7,1) - t99 * mrSges(5,2) - t33 * mrSges(6,2) - t25 * mrSges(7,2) - t182 * mrSges(4,3);
t483 = m(7) * pkin(5);
t482 = t55 / 0.2e1;
t481 = t56 / 0.2e1;
t480 = Ifges(5,1) * t476 + Ifges(5,4) * t475 + Ifges(5,5) * t463;
t306 = -pkin(9) - pkin(8);
t472 = t338 / 0.2e1;
t470 = -t136 / 0.2e1;
t469 = t136 / 0.2e1;
t464 = t184 / 0.2e1;
t460 = -t199 / 0.2e1;
t459 = -t200 / 0.2e1;
t458 = t200 / 0.2e1;
t457 = -t223 / 0.2e1;
t456 = t223 / 0.2e1;
t455 = -t231 / 0.2e1;
t453 = -t242 / 0.2e1;
t442 = pkin(4) * t200;
t439 = pkin(4) * t302;
t437 = pkin(7) * t304;
t436 = pkin(8) * t288;
t433 = g(3) * t288;
t36 = t302 * t85 - t81;
t428 = mrSges(5,3) * t199;
t427 = mrSges(5,3) * t200;
t426 = mrSges(6,3) * t338;
t425 = mrSges(6,3) * t136;
t424 = mrSges(7,3) * t338;
t423 = mrSges(7,3) * t136;
t422 = Ifges(3,4) * t300;
t421 = Ifges(3,4) * t304;
t420 = Ifges(5,4) * t299;
t419 = Ifges(5,4) * t303;
t414 = t200 * Ifges(5,4);
t413 = t242 * Ifges(4,4);
t400 = t244 * t303;
t395 = t266 * t301;
t294 = -qJ(6) + t306;
t394 = t288 * t294;
t393 = t288 * t306;
t387 = t299 * t301;
t386 = t299 * t305;
t385 = t301 * t303;
t198 = Ifges(5,4) * t199;
t123 = t200 * Ifges(5,1) + t231 * Ifges(5,5) + t198;
t384 = t303 * t123;
t383 = t303 * t305;
t375 = qJD(1) * t304;
t365 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t376) * t437;
t361 = pkin(2) * t374;
t359 = pkin(4) * t369;
t358 = pkin(4) * t368;
t354 = Ifges(5,5) * t128 + Ifges(5,6) * t129 + Ifges(5,3) * t187;
t353 = t406 * pkin(2);
t350 = t384 / 0.2e1;
t15 = -t56 * mrSges(7,1) + t55 * mrSges(7,2);
t347 = -t371 / 0.2e1;
t35 = -t298 * t85 - t83;
t346 = qJD(2) * t297;
t342 = -t190 * mrSges(4,1) + t191 * mrSges(4,2);
t57 = -t103 * t298 + t302 * t95;
t239 = t300 * t346 + t372;
t240 = t304 * t346 - t373;
t155 = t239 * t406 + t240 * t405;
t243 = t257 * qJD(2);
t157 = pkin(3) * t243 + pkin(8) * t244 + t361;
t339 = -t155 * t299 + t303 * t157;
t177 = -t302 * t252 - t253 * t298;
t281 = -t353 - pkin(3);
t336 = pkin(3) * t289 + t436;
t335 = mrSges(3,1) * t300 + mrSges(3,2) * t304;
t332 = mrSges(5,1) * t299 + mrSges(5,2) * t303;
t330 = Ifges(5,1) * t303 - t420;
t329 = t304 * Ifges(3,2) + t422;
t328 = -Ifges(5,2) * t299 + t419;
t327 = Ifges(3,5) * t304 - Ifges(3,6) * t300;
t326 = Ifges(5,5) * t303 - Ifges(5,6) * t299;
t119 = t173 * t406 - t405 * t185;
t154 = t239 * t405 - t406 * t240;
t196 = -t406 * t272 - t274 * t405;
t143 = -mrSges(5,2) * t231 + t428;
t144 = mrSges(5,1) * t231 - t427;
t324 = t143 * t303 - t144 * t299;
t323 = t263 * t289 - t394;
t322 = t289 * t284 - t393;
t269 = -t292 + t281;
t153 = pkin(4) * t397 + t196;
t320 = pkin(1) * t335;
t234 = -t289 * t386 + t385;
t232 = t289 * t387 + t383;
t318 = t257 * t371 + t400;
t46 = pkin(9) * t400 + pkin(4) * t243 + (-t192 + (pkin(9) * t257 - t180) * t299) * qJD(4) + t339;
t61 = t303 * t155 + t299 * t157 + t180 * t370 - t197 * t371;
t54 = -pkin(9) * t319 + t61;
t9 = -t103 * t369 + t298 * t46 + t302 * t54 + t95 * t368;
t317 = t166 * t332;
t316 = t300 * (Ifges(3,1) * t304 - t422);
t110 = pkin(4) * t319 + t154;
t111 = -qJDD(2) * pkin(3) - t119;
t63 = -t129 * pkin(4) + t111;
t10 = -qJD(5) * t58 - t298 * t54 + t302 * t46;
t310 = t487 + t506;
t287 = Ifges(3,4) * t375;
t283 = pkin(5) + t439;
t271 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t375;
t249 = Ifges(3,1) * t376 + Ifges(3,5) * qJD(2) + t287;
t248 = Ifges(3,6) * qJD(2) + qJD(1) * t329;
t235 = t289 * t383 + t387;
t233 = -t289 * t385 + t386;
t210 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t241;
t201 = pkin(5) * t321 + t269;
t174 = -mrSges(4,1) * t241 - mrSges(4,2) * t242;
t170 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t191;
t169 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t190;
t162 = t260 * t257;
t161 = -t242 * Ifges(4,1) + t227 + t513;
t160 = -t413 + t512 + t518;
t141 = -qJ(6) * t321 + t178;
t140 = -qJ(6) * t260 + t177;
t122 = t199 * Ifges(5,2) + t231 * Ifges(5,6) + t414;
t116 = mrSges(6,1) * t223 - t425;
t115 = mrSges(7,1) * t223 - t423;
t114 = -mrSges(6,2) * t223 + t426;
t113 = -mrSges(7,2) * t223 + t424;
t105 = t162 * pkin(5) + t153;
t104 = pkin(5) * t136 + t442;
t94 = -mrSges(5,2) * t187 + mrSges(5,3) * t129;
t93 = mrSges(5,1) * t187 - mrSges(5,3) * t128;
t88 = -mrSges(6,1) * t338 + mrSges(6,2) * t136;
t87 = -mrSges(7,1) * t338 + mrSges(7,2) * t136;
t75 = t163 * t500 + t260 * t244;
t74 = -t194 * t257 + t244 * t321;
t72 = -mrSges(5,1) * t129 + mrSges(5,2) * t128;
t62 = -qJD(4) * t125 + t339;
t59 = t128 * Ifges(5,4) + t129 * Ifges(5,2) + t187 * Ifges(5,6);
t47 = -t75 * pkin(5) + t110;
t41 = -qJ(6) * t162 + t58;
t38 = mrSges(6,1) * t184 - mrSges(6,3) * t55;
t37 = mrSges(7,1) * t184 - mrSges(7,3) * t55;
t34 = pkin(5) * t256 + qJ(6) * t163 + t57;
t27 = t36 - t542;
t26 = t35 - t511;
t20 = -t56 * pkin(5) + qJDD(6) + t63;
t16 = -mrSges(6,1) * t56 + mrSges(6,2) * t55;
t8 = qJ(6) * t75 - qJD(6) * t162 + t9;
t7 = pkin(5) * t243 - qJ(6) * t74 + qJD(6) * t163 + t10;
t1 = [(t304 * (-Ifges(3,2) * t300 + t421) + t316) * t367 / 0.2e1 - qJDD(2) * mrSges(3,2) * t437 + t174 * t361 + (-t531 * t469 + t528 * t456 + t519 / 0.2e1 + t520 / 0.2e1 - t518 / 0.2e1 - t160 / 0.2e1 + t529 * t472 + t486 - Ifges(4,6) * t432 - Ifges(4,4) * t453 + Ifges(5,5) * t458 + t498 / 0.2e1) * t243 + t231 * (-Ifges(5,5) * t318 - Ifges(5,6) * t319) / 0.2e1 + t166 * (mrSges(5,1) * t319 - mrSges(5,2) * t318) + t20 * (mrSges(7,1) * t162 - mrSges(7,2) * t163) + t63 * (mrSges(6,1) * t162 - mrSges(6,2) * t163) + (-Ifges(5,1) * t318 - Ifges(5,4) * t319) * t458 - t248 * t374 / 0.2e1 + t199 * (-Ifges(5,4) * t318 - Ifges(5,2) * t319) / 0.2e1 + t396 * t480 + m(7) * (t105 * t20 + t2 * t34 + t21 * t7 + t25 * t8 + t3 * t41 + t47 * t78) + m(6) * (t10 * t32 + t110 * t130 + t153 * t63 + t33 * t9 + t5 * t58 + t57 * t6) - t273 * t404 + t75 * t544 + t74 * t545 + t265 * t421 / 0.2e1 - t59 * t397 / 0.2e1 + (-t162 * t5 + t163 * t6 - t32 * t74 + t33 * t75) * mrSges(6,3) + (-t162 * t3 + t163 * t2 - t21 * t74 + t25 * t75) * mrSges(7,3) + (-t30 * t397 - t31 * t396 + t318 * t98 - t319 * t99) * mrSges(5,3) + (-m(4) * t119 + m(5) * t111 - t170 + t72) * t196 + (-t387 * t484 - m(7) * t395 - t235 * mrSges(5,1) - t234 * mrSges(5,2) - t496 * (t305 * t285 - t301 * t297) + t535 * t215 - t548 * t214 + t499 * t301 + (-m(5) * t336 - m(6) * t322 - m(7) * t323 - t488) * t305) * g(2) + (-t233 * mrSges(5,1) - t232 * mrSges(5,2) + t535 * t213 - t548 * t212 + (t297 * t496 + t491 + t499) * t305 + (m(4) * t285 - m(6) * (-t285 - t322) - m(7) * (-t285 - t323) - m(5) * (-t285 - t336) + t488) * t301) * g(1) - t285 * t342 - t320 * t367 - t271 * t357 + t162 * t543 + Ifges(2,3) * qJDD(1) + t304 * (Ifges(3,4) * t265 + Ifges(3,2) * t264 + Ifges(3,6) * qJDD(2)) / 0.2e1 + m(5) * (t124 * t31 + t125 * t30 + t61 * t99 + t62 * t98) + m(4) * (t120 * t197 + t155 * t182 - t226 * t285 + t268 * t361) - pkin(1) * (-mrSges(3,1) * t264 + mrSges(3,2) * t265) + t329 * t536 + t155 * t210 + t197 * t169 + t153 * t16 + t61 * t143 + t62 * t144 + t124 * t93 + t125 * t94 + t130 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) + t7 * t115 + t10 * t116 + t105 * t15 + t110 * t88 + t8 * t113 + t9 * t114 + t47 * t87 + t78 * (-mrSges(7,1) * t75 + mrSges(7,2) * t74) + (Ifges(3,1) * t265 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t265) + Ifges(3,4) * t536 + t497 * Ifges(3,5)) * t300 + (t327 * t432 - t365) * qJD(2) + (t530 * t75 + t532 * t74) * t472 + (-t162 * t530 - t163 * t532) * t481 + (t532 * t75 + t533 * t74) * t469 + (-t162 * t532 - t163 * t533) * t482 + t57 * t38 + t58 * t40 + (-t162 * t529 + t163 * t531) * t464 + (t529 * t75 - t531 * t74) * t456 + (t226 * mrSges(4,1) - t120 * mrSges(4,3) - Ifges(4,4) * t191 + Ifges(5,5) * t476 - Ifges(4,2) * t190 - Ifges(4,6) * t497 + Ifges(5,6) * t475 + Ifges(5,3) * t463 + t464 * t528 + t481 * t529 - t482 * t531 + t487 + t490) * t256 - t526 * t163 / 0.2e1 + t41 * t39 - (Ifges(4,5) * t432 + Ifges(4,1) * t453 + t227 / 0.2e1 + t161 / 0.2e1 + t350 - t495) * t244 + t34 * t37 + (-m(4) * t181 + m(5) * t166 - t508) * t154 + (t354 + t506) * t256 / 0.2e1 - t319 * t122 / 0.2e1 + (t264 * t437 + t503) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t503) + (t226 * mrSges(4,2) - t119 * mrSges(4,3) + Ifges(4,1) * t191 + Ifges(4,4) * t190 + Ifges(4,5) * t497 + t111 * t332 + t123 * t347 + t326 * t463 + t328 * t475 + t330 * t476) * t257 + Ifges(3,6) * t304 * t407 + t304 * t249 * t432; -(-Ifges(3,2) * t376 + t249 + t287) * t375 / 0.2e1 + (t199 * t328 + t200 * t330 + t231 * t326) * qJD(4) / 0.2e1 - (Ifges(4,2) * t242 + t161 + t227 + t384) * t241 / 0.2e1 - (-t531 * t456 + t533 * t469 + t532 * t472 + t492 + t545) * t193 + t538 * (t335 + t496 * pkin(2) * t300 + (-m(5) * pkin(8) + m(6) * t306 + m(7) * t294 - t546) * t289 + (mrSges(4,1) - t547) * t288) + t111 * t333 + t169 * t352 + t170 * t353 + (t248 / 0.2e1 + pkin(7) * t271) * t376 + (-t108 * t98 - t109 * t99 + t111 * t281 - t166 * t188) * m(5) + (t317 + t350) * qJD(4) + ((t119 * t406 + t120 * t405) * pkin(2) + t181 * t188 - t182 * t189 - t268 * t362) * m(4) + (-t549 * t99 + (-t370 + t402) * t98 + t502) * mrSges(5,3) + (t365 + (-t316 / 0.2e1 + t320) * qJD(1)) * qJD(1) + t299 * t480 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + t151 * t544 + t152 * t545 + (t388 / 0.2e1 + t347) * t122 + (t151 * t33 - t152 * t32 - t260 * t6 - t321 * t5) * mrSges(6,3) + (t151 * t25 - t152 * t21 - t2 * t260 - t3 * t321) * mrSges(7,3) + t63 * (mrSges(6,1) * t321 + mrSges(6,2) * t260) + t20 * (mrSges(7,1) * t321 + mrSges(7,2) * t260) + (t260 * t533 - t321 * t532) * t482 + (-t260 * t531 - t321 * t529) * t464 + (t260 * t532 - t321 * t530) * t481 + (Ifges(5,2) * t303 + t420) * t475 + (Ifges(5,1) * t299 + t419) * t476 + t160 * t453 + (Ifges(5,5) * t299 + Ifges(5,6) * t303) * t463 - t78 * (mrSges(7,1) * t151 - mrSges(7,2) * t152) - (t529 * t456 + t532 * t469 + t530 * t472 + t537) * t194 - t327 * t367 / 0.2e1 - t174 * t362 + t321 * t543 + t303 * t59 / 0.2e1 + t281 * t72 + Ifges(3,6) * t264 + Ifges(3,5) * t265 + t269 * t16 - t254 * mrSges(3,2) - t255 * mrSges(3,1) + t201 * t15 - t189 * t210 + Ifges(4,6) * t190 + Ifges(4,5) * t191 + t177 * t38 + t178 * t40 + t140 * t37 + t141 * t39 - t109 * t143 - t108 * t144 + t119 * mrSges(4,1) - t120 * mrSges(4,2) + (-t151 * t532 - t152 * t533) * t470 + (-Ifges(5,3) * t455 - Ifges(5,5) * t459 - Ifges(5,6) * t460 - t512 / 0.2e1 - t529 * t473 + t531 * t470 - t528 * t457 + t486) * t242 + (-t151 * t529 + t152 * t531) * t457 + (-t151 * t530 - t152 * t532) * t473 + t526 * t260 / 0.2e1 + t521 * t115 + t522 * t113 + (t140 * t2 + t141 * t3 + t20 * t201 + t21 * t521 + t25 * t522 + t516 * t78) * m(7) + t515 * t116 + t516 * t87 + (t130 * t507 + t177 * t6 + t178 * t5 + t269 * t63 + t32 * t515 + t33 * t517) * m(6) + t517 * t114 + (-t317 + t326 * t455 + t330 * t459 + t328 * t460 - t513 / 0.2e1 + t495) * t241 + (-mrSges(6,1) * t509 + mrSges(6,2) * t510) * t130 + t508 * t188 + t507 * t88 + (-t143 * t371 - t299 * t93 + m(5) * ((-t99 * t299 - t98 * t303) * qJD(4) + t502) - t144 * t370 + t303 * t94) * t280 + (Ifges(4,1) * t241 + t413 + t498) * t242 / 0.2e1 + (-m(5) * (t293 + t436) - m(4) * t293 - m(7) * (t293 - t394) - m(6) * (t293 - t393) + t547 * t289 + t540) * g(3); (-t210 - t324) * t241 + t324 * qJD(4) + t342 + t303 * t93 + t299 * t94 - (t37 + t38) * t321 + (t87 + t88 - t508) * t242 + t525 * t260 + (t114 + t113) * t510 + (t115 + t116) * t509 + (-g(1) * t301 + g(2) * t305) * t496 + (-t2 * t321 + t21 * t509 + t242 * t78 + t25 * t510 + t260 * t3) * m(7) + (t130 * t242 + t260 * t5 + t32 * t509 - t321 * t6 + t33 * t510) * m(6) + (t166 * t242 + t299 * t30 + t303 * t31 + t231 * (-t299 * t98 + t303 * t99)) * m(5) + (-t181 * t242 - t182 * t241 + t226) * m(4); (-Ifges(5,2) * t200 + t123 + t198) * t460 + (t358 - t36) * t114 + (t427 + t144) * t99 + t38 * t439 + (-t359 - t35) * t116 + (t358 - t27) * t113 + (t428 - t143) * t98 + t310 + (t2 * t283 + (t298 * t3 + (-t21 * t298 + t25 * t302) * qJD(5)) * pkin(4) - t104 * t78 - t21 * t26 - t25 * t27) * m(7) + (-t359 - t26) * t115 + (t298 * t5 + t302 * t6 + (-t298 * t32 + t302 * t33) * qJD(5)) * t484 - t88 * t442 + (Ifges(5,5) * t199 - Ifges(5,6) * t200) * t455 + t122 * t458 + (Ifges(5,1) * t199 - t414) * t459 + t354 + t490 - m(6) * (t130 * t442 + t32 * t35 + t33 * t36) + (-t130 * mrSges(6,1) - t529 * t457 - t532 * t470 - t530 * t473 + t537) * t136 + t283 * t37 + t523 * t473 - t166 * (mrSges(5,1) * t200 + mrSges(5,2) * t199) - t104 * t87 + (-t130 * mrSges(6,2) - t531 * t457 + t533 * t470 + t532 * t473 - t492) * t338 + (-m(7) * (-t267 * t305 - t289 * t395) - mrSges(5,2) * t233 + t514 * t232 + t504) * g(2) + (-m(7) * (-t266 * t289 * t305 + t267 * t301) + mrSges(5,2) * t235 - t514 * t234 + t505) * g(1) + (mrSges(7,1) * t290 + t332 - t491 + t501) * t433 + t525 * pkin(4) * t298; t21 * t424 + t310 + t2 * t483 - t130 * (mrSges(6,1) * t136 + mrSges(6,2) * t338) - t78 * (mrSges(7,1) * t136 + mrSges(7,2) * t338) - t24 * t113 + (t338 * t533 - t539) * t470 + t524 * t469 + (-t136 * t529 - t338 * t531) * t457 + (-(-mrSges(7,1) - t483) * t290 + t501) * t433 + (t425 + t116) * t33 + (t426 - t114) * t32 + (t423 - m(7) * (-t21 + t24) + t115) * t25 + (t212 * t483 + t504) * g(2) + (-t214 * t483 + t505) * g(1) + (-t136 * t530 + t523 + t541) * t473 + (t37 + (-m(7) * t78 - t87) * t136) * pkin(5); -t338 * t113 + t136 * t115 + (g(3) * t289 + t21 * t136 - t338 * t25 - t288 * t538 + t20) * m(7) + t15;];
tau  = t1;
