% Calculate vector of inverse dynamics joint torques for
% S5RRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:27:44
% EndTime: 2019-12-31 22:28:39
% DurationCPUTime: 28.90s
% Computational Cost: add. (12306->794), mult. (27428->1088), div. (0->0), fcn. (19180->14), ass. (0->359)
t298 = sin(qJ(2));
t303 = cos(qJ(2));
t342 = pkin(2) * t298 - pkin(7) * t303;
t246 = t342 * qJD(1);
t302 = cos(qJ(3));
t297 = sin(qJ(3));
t382 = qJD(1) * t298;
t360 = t297 * t382;
t172 = pkin(6) * t360 + t302 * t246;
t390 = t302 * t303;
t324 = pkin(3) * t298 - pkin(8) * t390;
t305 = -pkin(8) - pkin(7);
t361 = qJD(3) * t305;
t525 = -qJD(1) * t324 + t302 * t361 - t172;
t219 = t297 * t246;
t393 = t298 * t302;
t395 = t297 * t303;
t524 = t219 + (-pkin(6) * t393 - pkin(8) * t395) * qJD(1) - t297 * t361;
t296 = sin(qJ(4));
t301 = cos(qJ(4));
t241 = t296 * t302 + t297 * t301;
t487 = qJD(3) + qJD(4);
t167 = t487 * t241;
t319 = t241 * t303;
t192 = qJD(1) * t319;
t523 = t167 - t192;
t262 = t305 * t297;
t263 = t305 * t302;
t373 = qJD(4) * t301;
t374 = qJD(4) * t296;
t501 = t262 * t373 + t263 * t374 + t525 * t296 - t524 * t301;
t171 = t296 * t262 - t301 * t263;
t500 = -qJD(4) * t171 + t524 * t296 + t525 * t301;
t522 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3);
t521 = m(4) * pkin(6);
t370 = qJD(1) * qJD(2);
t251 = qJDD(1) * t298 + t303 * t370;
t520 = t251 / 0.2e1;
t325 = t296 * t297 - t301 * t302;
t166 = t487 * t325;
t318 = t325 * t303;
t193 = qJD(1) * t318;
t519 = -pkin(4) * t382 + t500 + (t166 - t193) * pkin(9);
t518 = t523 * pkin(9) - t501;
t340 = mrSges(3,1) * t298 + mrSges(3,2) * t303;
t419 = Ifges(3,4) * t298;
t517 = pkin(1) * t340 - t298 * (Ifges(3,1) * t303 - t419) / 0.2e1;
t379 = qJD(2) * t302;
t238 = -t360 + t379;
t283 = pkin(6) * t382;
t260 = -qJD(2) * pkin(2) + t283;
t176 = -pkin(3) * t238 + t260;
t358 = t302 * t382;
t239 = qJD(2) * t297 + t358;
t344 = t301 * t238 - t239 * t296;
t109 = -pkin(4) * t344 + t176;
t295 = sin(qJ(5));
t300 = cos(qJ(5));
t156 = t238 * t296 + t239 * t301;
t345 = -t156 * t295 + t300 * t344;
t146 = qJD(3) * t238 + qJDD(2) * t297 + t251 * t302;
t147 = -qJD(3) * t239 + qJDD(2) * t302 - t251 * t297;
t60 = qJD(4) * t344 + t146 * t301 + t147 * t296;
t61 = -qJD(4) * t156 - t146 * t296 + t147 * t301;
t18 = qJD(5) * t345 + t295 * t61 + t300 * t60;
t95 = t156 * t300 + t295 * t344;
t19 = -qJD(5) * t95 - t295 * t60 + t300 * t61;
t250 = qJDD(1) * t303 - t298 * t370;
t235 = qJDD(3) - t250;
t226 = qJDD(4) + t235;
t210 = qJDD(5) + t226;
t367 = Ifges(6,5) * t18 + Ifges(6,6) * t19 + Ifges(6,3) * t210;
t343 = pkin(2) * t303 + pkin(7) * t298;
t255 = -pkin(1) - t343;
t227 = t255 * qJD(1);
t381 = qJD(1) * t303;
t284 = pkin(6) * t381;
t261 = qJD(2) * pkin(7) + t284;
t161 = t302 * t227 - t261 * t297;
t123 = -pkin(8) * t239 + t161;
t271 = qJD(3) - t381;
t111 = pkin(3) * t271 + t123;
t162 = t227 * t297 + t261 * t302;
t124 = pkin(8) * t238 + t162;
t404 = qJDD(1) * pkin(1);
t165 = -pkin(2) * t250 - pkin(7) * t251 - t404;
t233 = t250 * pkin(6);
t205 = qJDD(2) * pkin(7) + t233;
t76 = -qJD(3) * t162 + t302 * t165 - t205 * t297;
t51 = pkin(3) * t235 - pkin(8) * t146 + t76;
t375 = qJD(3) * t302;
t377 = qJD(3) * t297;
t75 = t297 * t165 + t302 * t205 + t227 * t375 - t261 * t377;
t57 = pkin(8) * t147 + t75;
t12 = t111 * t373 - t124 * t374 + t296 * t51 + t301 * t57;
t10 = pkin(9) * t61 + t12;
t505 = pkin(9) * t344;
t118 = t301 * t124;
t63 = t111 * t296 + t118;
t50 = t63 + t505;
t411 = t295 * t50;
t264 = qJD(4) + t271;
t509 = pkin(9) * t156;
t116 = t296 * t124;
t62 = t301 * t111 - t116;
t49 = t62 - t509;
t45 = pkin(4) * t264 + t49;
t20 = t300 * t45 - t411;
t13 = -qJD(4) * t63 - t296 * t57 + t301 * t51;
t9 = pkin(4) * t226 - pkin(9) * t60 + t13;
t2 = qJD(5) * t20 + t10 * t300 + t295 * t9;
t406 = t300 * t50;
t21 = t295 * t45 + t406;
t3 = -qJD(5) * t21 - t10 * t295 + t300 * t9;
t513 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t323 = t367 - t513;
t366 = Ifges(5,5) * t60 + Ifges(5,6) * t61 + Ifges(5,3) * t226;
t434 = mrSges(6,3) * t21;
t435 = mrSges(6,3) * t20;
t256 = qJD(5) + t264;
t442 = -t256 / 0.2e1;
t455 = -t95 / 0.2e1;
t457 = -t345 / 0.2e1;
t89 = Ifges(6,4) * t345;
t44 = Ifges(6,1) * t95 + Ifges(6,5) * t256 + t89;
t466 = -t44 / 0.2e1;
t433 = Ifges(6,4) * t95;
t43 = Ifges(6,2) * t345 + Ifges(6,6) * t256 + t433;
t468 = -t43 / 0.2e1;
t511 = -t13 * mrSges(5,1) + t12 * mrSges(5,2);
t516 = t323 + t366 - t511 + (-mrSges(6,2) * t109 + Ifges(6,1) * t455 + Ifges(6,4) * t457 + Ifges(6,5) * t442 + t435 + t466) * t345 - (mrSges(6,1) * t109 + Ifges(6,4) * t455 + Ifges(6,2) * t457 + Ifges(6,6) * t442 - t434 + t468) * t95;
t294 = qJ(3) + qJ(4);
t287 = cos(t294);
t429 = pkin(3) * t302;
t254 = pkin(4) * t287 + t429;
t245 = pkin(2) + t254;
t288 = qJ(5) + t294;
t276 = sin(t288);
t277 = cos(t288);
t279 = pkin(2) + t429;
t286 = sin(t294);
t339 = -mrSges(4,1) * t302 + mrSges(4,2) * t297;
t515 = m(4) * pkin(2) + m(5) * t279 + m(6) * t245 + mrSges(5,1) * t287 + t277 * mrSges(6,1) - mrSges(5,2) * t286 - t276 * mrSges(6,2) - t339;
t293 = -pkin(9) + t305;
t514 = -m(4) * pkin(7) + m(5) * t305 + m(6) * t293 - t522;
t512 = -t76 * mrSges(4,1) + t75 * mrSges(4,2);
t333 = t303 * Ifges(3,2) + t419;
t510 = t21 * mrSges(6,2) + t63 * mrSges(5,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t333 / 0.2e1 - t20 * mrSges(6,1) - t62 * mrSges(5,1);
t428 = pkin(4) * t156;
t359 = t297 * t381;
t495 = -t284 + (-t359 + t377) * pkin(3);
t378 = qJD(2) * t303;
t313 = t297 * t378 + t298 * t375;
t474 = m(5) * pkin(3);
t472 = t18 / 0.2e1;
t471 = t19 / 0.2e1;
t464 = t60 / 0.2e1;
t463 = t61 / 0.2e1;
t453 = t146 / 0.2e1;
t452 = t147 / 0.2e1;
t447 = t210 / 0.2e1;
t446 = t226 / 0.2e1;
t445 = t235 / 0.2e1;
t170 = t301 * t262 + t263 * t296;
t131 = -pkin(9) * t241 + t170;
t132 = -pkin(9) * t325 + t171;
t80 = t131 * t295 + t132 * t300;
t504 = -qJD(5) * t80 + t295 * t518 + t300 * t519;
t79 = t131 * t300 - t132 * t295;
t503 = qJD(5) * t79 + t295 * t519 - t300 * t518;
t502 = mrSges(4,1) + t474;
t430 = pkin(3) * t301;
t278 = pkin(4) + t430;
t371 = qJD(5) * t300;
t372 = qJD(5) * t295;
t399 = t295 * t296;
t68 = -t123 * t296 - t118;
t52 = t68 - t505;
t69 = t301 * t123 - t116;
t53 = t69 - t509;
t499 = -t295 * t52 - t300 * t53 + t278 * t371 + (-t296 * t372 + (t300 * t301 - t399) * qJD(4)) * pkin(3);
t398 = t296 * t300;
t498 = t295 * t53 - t300 * t52 - t278 * t372 + (-t296 * t371 + (-t295 * t301 - t398) * qJD(4)) * pkin(3);
t204 = t325 * t298;
t496 = t523 * pkin(4) + t495;
t237 = t302 * t255;
t160 = -pkin(8) * t393 + t237 + (-pkin(6) * t297 - pkin(3)) * t303;
t273 = pkin(6) * t390;
t187 = t297 * t255 + t273;
t397 = t297 * t298;
t169 = -pkin(8) * t397 + t187;
t101 = t296 * t160 + t301 * t169;
t494 = qJD(2) * mrSges(3,1) + mrSges(4,1) * t238 - mrSges(4,2) * t239 - mrSges(3,3) * t382;
t336 = -mrSges(6,1) * t276 - mrSges(6,2) * t277;
t493 = mrSges(5,1) * t286 + mrSges(5,2) * t287 - t336;
t299 = sin(qJ(1));
t304 = cos(qJ(1));
t389 = t303 * t304;
t199 = -t286 * t389 + t287 * t299;
t200 = t286 * t299 + t287 * t389;
t184 = -t276 * t389 + t277 * t299;
t185 = t276 * t299 + t277 * t389;
t387 = t184 * mrSges(6,1) - t185 * mrSges(6,2);
t492 = -t199 * mrSges(5,1) + t200 * mrSges(5,2) - t387;
t392 = t299 * t303;
t197 = t286 * t392 + t287 * t304;
t198 = t286 * t304 - t287 * t392;
t182 = t276 * t392 + t277 * t304;
t183 = t276 * t304 - t277 * t392;
t388 = -t182 * mrSges(6,1) + t183 * mrSges(6,2);
t491 = t197 * mrSges(5,1) - t198 * mrSges(5,2) - t388;
t231 = Ifges(4,4) * t238;
t135 = t239 * Ifges(4,1) + t271 * Ifges(4,5) + t231;
t282 = Ifges(3,4) * t381;
t490 = Ifges(3,1) * t382 + Ifges(3,5) * qJD(2) + t302 * t135 + t282;
t234 = t251 * pkin(6);
t489 = t233 * t303 + t234 * t298;
t488 = -t297 * t76 + t302 * t75;
t486 = t239 * Ifges(4,5) + t156 * Ifges(5,5) + t95 * Ifges(6,5) + t238 * Ifges(4,6) + Ifges(5,6) * t344 + Ifges(6,6) * t345 + t271 * Ifges(4,3) + t264 * Ifges(5,3) + t256 * Ifges(6,3);
t485 = -m(3) - m(5) - m(4) - m(6);
t484 = -mrSges(5,1) * t176 + mrSges(5,3) * t63;
t483 = mrSges(5,2) * t176 - mrSges(5,3) * t62;
t427 = pkin(4) * t286;
t431 = pkin(3) * t297;
t253 = t427 + t431;
t481 = -m(6) * t253 + mrSges(2,2) - mrSges(3,3);
t341 = t303 * mrSges(3,1) - mrSges(3,2) * t298;
t480 = t522 * t298 + mrSges(2,1) + t341;
t476 = Ifges(6,4) * t472 + Ifges(6,2) * t471 + Ifges(6,6) * t447;
t475 = Ifges(6,1) * t472 + Ifges(6,4) * t471 + Ifges(6,5) * t447;
t473 = m(6) * pkin(4);
t470 = Ifges(5,4) * t464 + Ifges(5,2) * t463 + Ifges(5,6) * t446;
t469 = Ifges(5,1) * t464 + Ifges(5,4) * t463 + Ifges(5,5) * t446;
t467 = t43 / 0.2e1;
t465 = t44 / 0.2e1;
t462 = Ifges(4,1) * t453 + Ifges(4,4) * t452 + Ifges(4,5) * t445;
t415 = Ifges(5,4) * t156;
t84 = Ifges(5,2) * t344 + Ifges(5,6) * t264 + t415;
t461 = -t84 / 0.2e1;
t460 = t84 / 0.2e1;
t150 = Ifges(5,4) * t344;
t85 = Ifges(5,1) * t156 + Ifges(5,5) * t264 + t150;
t459 = -t85 / 0.2e1;
t458 = t85 / 0.2e1;
t456 = t345 / 0.2e1;
t454 = t95 / 0.2e1;
t451 = -t344 / 0.2e1;
t450 = t344 / 0.2e1;
t449 = -t156 / 0.2e1;
t448 = t156 / 0.2e1;
t443 = t239 / 0.2e1;
t441 = t256 / 0.2e1;
t440 = -t264 / 0.2e1;
t439 = t264 / 0.2e1;
t432 = pkin(3) * t239;
t424 = g(3) * t298;
t289 = t298 * pkin(6);
t421 = mrSges(5,3) * t344;
t420 = mrSges(5,3) * t156;
t418 = Ifges(3,4) * t303;
t417 = Ifges(4,4) * t297;
t416 = Ifges(4,4) * t302;
t414 = t161 * mrSges(4,3);
t413 = t162 * mrSges(4,3);
t412 = t239 * Ifges(4,4);
t396 = t297 * t299;
t394 = t297 * t304;
t249 = t342 * qJD(2);
t380 = qJD(2) * t298;
t364 = pkin(6) * t380;
t384 = t302 * t249 + t297 * t364;
t252 = pkin(3) * t397 + t289;
t376 = qJD(3) * t298;
t285 = pkin(6) * t378;
t362 = Ifges(4,5) * t146 + Ifges(4,6) * t147 + Ifges(4,3) * t235;
t181 = pkin(3) * t313 + t285;
t134 = t238 * Ifges(4,2) + t271 * Ifges(4,6) + t412;
t355 = -t297 * t134 / 0.2e1;
t100 = t301 * t160 - t169 * t296;
t206 = -qJDD(2) * pkin(2) + t234;
t338 = mrSges(4,1) * t297 + mrSges(4,2) * t302;
t335 = Ifges(4,1) * t302 - t417;
t334 = Ifges(4,1) * t297 + t416;
t332 = -Ifges(4,2) * t297 + t416;
t331 = Ifges(4,2) * t302 + t417;
t330 = Ifges(3,5) * t303 - Ifges(3,6) * t298;
t329 = Ifges(4,5) * t302 - Ifges(4,6) * t297;
t328 = Ifges(4,5) * t297 + Ifges(4,6) * t302;
t72 = -pkin(4) * t303 + pkin(9) * t204 + t100;
t203 = t241 * t298;
t81 = -pkin(9) * t203 + t101;
t37 = -t295 * t81 + t300 * t72;
t38 = t295 * t72 + t300 * t81;
t127 = -t203 * t300 + t204 * t295;
t128 = -t203 * t295 - t204 * t300;
t157 = -t241 * t295 - t300 * t325;
t158 = t241 * t300 - t295 * t325;
t327 = t245 * t303 - t293 * t298;
t326 = t303 * t279 - t298 * t305;
t216 = -t297 * t389 + t299 * t302;
t214 = t297 * t392 + t302 * t304;
t321 = t260 * t338;
t121 = t297 * t249 + t255 * t375 + (-t298 * t379 - t303 * t377) * pkin(6);
t103 = -pkin(8) * t313 + t121;
t99 = t324 * qJD(2) + (-t273 + (pkin(8) * t298 - t255) * t297) * qJD(3) + t384;
t35 = t301 * t103 + t160 * t373 - t169 * t374 + t296 * t99;
t112 = -pkin(3) * t147 + t206;
t314 = -t297 * t376 + t302 * t378;
t312 = Ifges(4,5) * t298 + t303 * t335;
t311 = Ifges(4,6) * t298 + t303 * t332;
t310 = Ifges(4,3) * t298 + t303 * t329;
t36 = -qJD(4) * t101 - t103 * t296 + t301 * t99;
t258 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t381;
t228 = t338 * t298;
t217 = t302 * t389 + t396;
t215 = -t299 * t390 + t394;
t208 = pkin(3) * t398 + t278 * t295;
t207 = -pkin(3) * t399 + t278 * t300;
t194 = pkin(4) * t325 - t279;
t186 = -pkin(6) * t395 + t237;
t175 = mrSges(4,1) * t271 - mrSges(4,3) * t239;
t174 = -mrSges(4,2) * t271 + mrSges(4,3) * t238;
t173 = -pkin(6) * t358 + t219;
t163 = pkin(4) * t203 + t252;
t126 = mrSges(5,1) * t264 - t420;
t125 = -mrSges(5,2) * t264 + t421;
t122 = -qJD(3) * t187 + t384;
t120 = -t192 * t295 - t193 * t300;
t119 = -t192 * t300 + t193 * t295;
t115 = t432 + t428;
t114 = -mrSges(4,2) * t235 + mrSges(4,3) * t147;
t113 = mrSges(4,1) * t235 - mrSges(4,3) * t146;
t105 = -qJD(2) * t319 + t204 * t487;
t104 = -qJD(2) * t318 - t167 * t298;
t98 = -mrSges(5,1) * t344 + mrSges(5,2) * t156;
t88 = -mrSges(4,1) * t147 + mrSges(4,2) * t146;
t82 = -pkin(4) * t105 + t181;
t78 = mrSges(6,1) * t256 - mrSges(6,3) * t95;
t77 = -mrSges(6,2) * t256 + mrSges(6,3) * t345;
t73 = t146 * Ifges(4,4) + t147 * Ifges(4,2) + t235 * Ifges(4,6);
t66 = -qJD(5) * t158 + t166 * t295 - t167 * t300;
t65 = qJD(5) * t157 - t166 * t300 - t167 * t295;
t55 = -mrSges(5,2) * t226 + mrSges(5,3) * t61;
t54 = mrSges(5,1) * t226 - mrSges(5,3) * t60;
t46 = -mrSges(6,1) * t345 + mrSges(6,2) * t95;
t41 = -pkin(4) * t61 + t112;
t40 = -qJD(5) * t128 - t104 * t295 + t105 * t300;
t39 = qJD(5) * t127 + t104 * t300 + t105 * t295;
t32 = -mrSges(5,1) * t61 + mrSges(5,2) * t60;
t29 = pkin(9) * t105 + t35;
t26 = pkin(4) * t380 - pkin(9) * t104 + t36;
t23 = t300 * t49 - t411;
t22 = -t295 * t49 - t406;
t15 = -mrSges(6,2) * t210 + mrSges(6,3) * t19;
t14 = mrSges(6,1) * t210 - mrSges(6,3) * t18;
t8 = -mrSges(6,1) * t19 + mrSges(6,2) * t18;
t5 = -qJD(5) * t38 + t26 * t300 - t29 * t295;
t4 = qJD(5) * t37 + t26 * t295 + t29 * t300;
t1 = [(-t161 * t314 - t162 * t313 - t393 * t76 - t397 * t75) * mrSges(4,3) + t341 * t404 + t489 * mrSges(3,3) + (-qJDD(2) * mrSges(3,1) + mrSges(3,3) * t251 + t88) * t289 + t418 * t520 + (t419 + t333) * t250 / 0.2e1 + (Ifges(6,4) * t128 + Ifges(6,2) * t127) * t471 + (Ifges(6,4) * t39 + Ifges(6,2) * t40) * t456 - t517 * t370 + t238 * (qJD(2) * t311 - t331 * t376) / 0.2e1 + t271 * (qJD(2) * t310 - t328 * t376) / 0.2e1 + t206 * t228 - t258 * t364 + t128 * t475 + t127 * t476 + (Ifges(5,5) * t104 + Ifges(5,6) * t105) * t439 + (-Ifges(5,5) * t204 - Ifges(5,6) * t203) * t446 + t186 * t113 + t187 * t114 + t121 * t174 + t122 * t175 + t176 * (-mrSges(5,1) * t105 + mrSges(5,2) * t104) + t181 * t98 + t163 * t8 - t73 * t397 / 0.2e1 + (t161 * mrSges(4,1) + t486 / 0.2e1 - t162 * mrSges(4,2) + Ifges(5,5) * t448 + Ifges(5,6) * t450 + Ifges(6,5) * t454 + Ifges(6,6) * t456 + Ifges(5,3) * t439 + Ifges(6,3) * t441 - t510) * t380 + (Ifges(6,5) * t128 + Ifges(6,6) * t127) * t447 + t104 * t458 + t105 * t460 + t393 * t462 + t39 * t465 + t40 * t467 - t204 * t469 - t203 * t470 + m(4) * (t121 * t162 + t122 * t161 + t186 * t76 + t187 * t75) + m(6) * (t109 * t82 + t163 * t41 + t2 * t38 + t20 * t5 + t21 * t4 + t3 * t37) + m(5) * (t100 * t13 + t101 * t12 + t112 * t252 + t176 * t181 + t35 * t63 + t36 * t62) + (t511 + (Ifges(3,2) / 0.2e1 + pkin(6) * mrSges(3,3)) * t250 - Ifges(4,3) * t445 - Ifges(6,3) * t447 - Ifges(4,6) * t452 - Ifges(4,5) * t453 - Ifges(5,6) * t463 - Ifges(5,5) * t464 - Ifges(5,3) * t446 + (-mrSges(3,2) * pkin(6) + Ifges(3,6)) * qJDD(2) - Ifges(6,6) * t471 - Ifges(6,5) * t472 + (-Ifges(3,2) * t298 + t418) * t370 / 0.2e1 - t366 / 0.2e1 - t367 / 0.2e1 - t362 / 0.2e1 + Ifges(3,4) * t520 + t512 + t513) * t303 + (t355 + t260 * t521 + t490 / 0.2e1) * t378 + (Ifges(3,1) * t251 + Ifges(3,5) * qJDD(2) + t206 * t521 + t329 * t445 + t332 * t452 + t335 * t453) * t298 + (-t104 * t62 + t105 * t63 - t12 * t203 + t13 * t204) * mrSges(5,3) - (t302 * t134 + t297 * t135) * t376 / 0.2e1 + (Ifges(6,5) * t39 + Ifges(6,6) * t40) * t441 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t489) - t494 * t285 + (qJD(2) * t312 - t334 * t376) * t443 + (Ifges(5,4) * t104 + Ifges(5,2) * t105) * t450 + (-Ifges(5,4) * t204 - Ifges(5,2) * t203) * t463 + t112 * (mrSges(5,1) * t203 - mrSges(5,2) * t204) + (Ifges(6,1) * t128 + Ifges(6,4) * t127) * t472 + (Ifges(6,1) * t39 + Ifges(6,4) * t40) * t454 + t41 * (-mrSges(6,1) * t127 + mrSges(6,2) * t128) + t35 * t125 + t36 * t126 + t109 * (-mrSges(6,1) * t40 + mrSges(6,2) * t39) + Ifges(2,3) * qJDD(1) + t100 * t54 + t101 * t55 + t5 * t78 + t82 * t46 + t4 * t77 + (Ifges(5,1) * t104 + Ifges(5,4) * t105) * t448 + (-Ifges(5,1) * t204 - Ifges(5,4) * t203) * t464 + qJD(2) ^ 2 * t330 / 0.2e1 + t37 * t14 + t38 * t15 + (t127 * t2 - t128 * t3 - t20 * t39 + t21 * t40) * mrSges(6,3) - pkin(1) * (-mrSges(3,1) * t250 + mrSges(3,2) * t251) + t252 * t32 + (-t396 * t474 - t217 * mrSges(4,1) - t200 * mrSges(5,1) - t185 * mrSges(6,1) - t216 * mrSges(4,2) - t199 * mrSges(5,2) - t184 * mrSges(6,2) + t485 * (t304 * pkin(1) + t299 * pkin(6)) + t481 * t299 + (-m(4) * t343 - m(5) * t326 - m(6) * t327 - t480) * t304) * g(2) + (-t394 * t474 - t215 * mrSges(4,1) - t198 * mrSges(5,1) - t183 * mrSges(6,1) - t214 * mrSges(4,2) - t197 * mrSges(5,2) - t182 * mrSges(6,2) + (-m(5) * (-pkin(1) - t326) - m(4) * t255 - m(6) * (-pkin(1) - t327) + m(3) * pkin(1) + t480) * t299 + (t485 * pkin(6) + t481) * t304) * g(1) + t260 * (mrSges(4,1) * t313 + mrSges(4,2) * t314); -t325 * t470 + (Ifges(5,4) * t241 - Ifges(5,2) * t325) * t463 + (Ifges(5,1) * t241 - Ifges(5,4) * t325) * t464 + (Ifges(5,5) * t241 - Ifges(5,6) * t325) * t446 + t112 * (mrSges(5,1) * t325 + mrSges(5,2) * t241) + (-t12 * t325 - t13 * t241 + t192 * t63 - t193 * t62) * mrSges(5,3) - t176 * (mrSges(5,1) * t192 - mrSges(5,2) * t193) + (-Ifges(5,5) * t193 - Ifges(5,6) * t192) * t440 + (-Ifges(5,4) * t193 - Ifges(5,2) * t192) * t451 + (-Ifges(5,1) * t193 - Ifges(5,4) * t192) * t449 + (-t161 * (mrSges(4,1) * t298 - mrSges(4,3) * t390) - t162 * (-mrSges(4,2) * t298 - mrSges(4,3) * t395)) * qJD(1) - t330 * t370 / 0.2e1 - t321 * t381 + (Ifges(6,1) * t120 + Ifges(6,4) * t119) * t455 + (t304 * g(1) + t299 * g(2)) * (t298 * t515 + t303 * t514 + t340) + (t298 * t514 - t303 * t515 - t341) * g(3) + t517 * qJD(1) ^ 2 + (t355 + t321) * qJD(3) + t495 * t98 + t496 * t46 + (-pkin(2) * t206 - t161 * t172 - t162 * t173 - t260 * t284) * m(4) - t486 * t382 / 0.2e1 + t258 * t283 - t377 * t413 + t158 * t475 + t157 * t476 + t194 * t8 - t173 * t174 - t172 * t175 + t170 * t54 + t171 * t55 + (Ifges(6,5) * t120 + Ifges(6,6) * t119) * t442 + (Ifges(5,5) * t449 + Ifges(6,5) * t455 + Ifges(5,6) * t451 + Ifges(6,6) * t457 + Ifges(5,3) * t440 + Ifges(6,3) * t442 + t510) * t382 - t193 * t459 - t192 * t461 + t297 * t462 + t120 * t466 + t119 * t468 + t241 * t469 + (Ifges(6,4) * t158 + Ifges(6,2) * t157) * t471 + (Ifges(6,1) * t158 + Ifges(6,4) * t157) * t472 + (Ifges(6,1) * t454 + Ifges(6,4) * t456 + Ifges(6,5) * t441 - t435 + t465) * t65 - (Ifges(5,1) * t448 + Ifges(5,4) * t450 + Ifges(5,5) * t439 + t458 + t483) * t166 - (Ifges(5,4) * t448 + Ifges(5,2) * t450 + Ifges(5,6) * t439 + t460 + t484) * t167 + (t238 * t332 + t239 * t335 + t271 * t329) * qJD(3) / 0.2e1 - (t238 * t311 + t239 * t312 + t271 * t310) * qJD(1) / 0.2e1 + (m(4) * ((-t161 * t302 - t162 * t297) * qJD(3) + t488) - t175 * t375 - t174 * t377 - t297 * t113 + t302 * t114) * pkin(7) + t488 * mrSges(4,3) - (-Ifges(3,2) * t382 + t282 + t490) * t381 / 0.2e1 + t494 * t284 + t328 * t445 + (Ifges(6,5) * t158 + Ifges(6,6) * t157) * t447 + t331 * t452 + t334 * t453 + (Ifges(6,4) * t120 + Ifges(6,2) * t119) * t457 + t134 * t359 / 0.2e1 + (t135 / 0.2e1 - t414) * t375 + ((-t120 + t65) * mrSges(6,2) + (t119 - t66) * mrSges(6,1)) * t109 + t41 * (-mrSges(6,1) * t157 + mrSges(6,2) * t158) + (Ifges(6,4) * t454 + Ifges(6,2) * t456 + Ifges(6,6) * t441 + t434 + t467) * t66 + (-t119 * t21 + t120 * t20 + t157 * t2 - t158 * t3) * mrSges(6,3) + t79 * t14 + t80 * t15 - pkin(2) * t88 + t206 * t339 - t233 * mrSges(3,2) - t234 * mrSges(3,1) + t500 * t126 + t501 * t125 + (-t112 * t279 + t12 * t171 + t13 * t170 + t176 * t495 + t500 * t62 + t501 * t63) * m(5) + Ifges(3,6) * t250 + Ifges(3,5) * t251 - t279 * t32 + t503 * t77 + t504 * t78 + (t109 * t496 + t194 * t41 + t2 * t80 + t20 * t504 + t21 * t503 + t3 * t79) * m(6) + t302 * t73 / 0.2e1 + Ifges(3,3) * qJDD(2); -t512 - (Ifges(5,4) * t449 + Ifges(5,2) * t451 + Ifges(5,6) * t440 + t461 - t484) * t156 + t516 + g(3) * t228 - t239 * (Ifges(4,1) * t238 - t412) / 0.2e1 + t208 * t15 + (t12 * t296 + t13 * t301 + (-t296 * t62 + t301 * t63) * qJD(4)) * t474 + t207 * t14 - t161 * t174 + t162 * t175 + t239 * t413 + t238 * t414 + (Ifges(5,1) * t449 + Ifges(5,4) * t451 + Ifges(5,5) * t440 + t459 - t483) * t344 - (-Ifges(4,2) * t239 + t135 + t231) * t238 / 0.2e1 - t98 * t432 - m(5) * (t176 * t432 + t62 * t68 + t63 * t69) + (m(5) * t431 + t493) * t424 + t134 * t443 + t54 * t430 + (t125 * t373 - t126 * t374 + t296 * t55) * pkin(3) + t362 - t69 * t125 - t68 * t126 - t115 * t46 + t498 * t78 + t499 * t77 + (-t109 * t115 + t2 * t208 + t20 * t498 + t207 * t3 + t21 * t499 + t253 * t424) * m(6) + (-mrSges(4,2) * t215 - m(6) * (-t253 * t392 - t254 * t304) + t502 * t214 + t491) * g(2) + (mrSges(4,2) * t217 - m(6) * (-t253 * t389 + t254 * t299) - t502 * t216 + t492) * g(1) - t260 * (mrSges(4,1) * t239 + mrSges(4,2) * t238) - t271 * (Ifges(4,5) * t238 - Ifges(4,6) * t239) / 0.2e1; (t2 * t295 + t3 * t300 + (-t20 * t295 + t21 * t300) * qJD(5)) * t473 - t176 * (mrSges(5,1) * t156 + mrSges(5,2) * t344) - t46 * t428 - m(6) * (t109 * t428 + t20 * t22 + t21 * t23) + t84 * t448 + (Ifges(5,1) * t344 - t415) * t449 + (Ifges(5,5) * t344 - Ifges(5,6) * t156) * t440 - t23 * t77 - t22 * t78 + (t420 + t126) * t63 + (t421 - t125) * t62 + (-Ifges(5,2) * t156 + t150 + t85) * t451 + (m(6) * t427 + t493) * t424 + (t197 * t473 + t491) * g(2) + (-t199 * t473 + t492) * g(1) + (t300 * t14 + t295 * t15 + t371 * t77 - t372 * t78) * pkin(4) + t516; -t109 * (mrSges(6,1) * t95 + mrSges(6,2) * t345) + (Ifges(6,1) * t345 - t433) * t455 + t43 * t454 + (Ifges(6,5) * t345 - Ifges(6,6) * t95) * t442 - t20 * t77 + t21 * t78 - g(1) * t387 - g(2) * t388 - t336 * t424 + (t20 * t345 + t21 * t95) * mrSges(6,3) + t323 + (-Ifges(6,2) * t95 + t44 + t89) * t457;];
tau = t1;
