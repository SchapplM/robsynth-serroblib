% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:51:17
% EndTime: 2018-11-23 17:51:27
% DurationCPUTime: 10.81s
% Computational Cost: add. (22982->685), mult. (58941->933), div. (0->0), fcn. (44710->10), ass. (0->337)
t338 = qJD(2) + qJD(3);
t517 = t338 / 0.2e1;
t342 = sin(qJ(5));
t449 = -t342 / 0.2e1;
t340 = cos(pkin(11));
t347 = cos(qJ(3));
t339 = sin(pkin(11));
t343 = sin(qJ(3));
t402 = t339 * t343;
t435 = pkin(2) * qJD(3);
t287 = (t340 * t347 - t402) * t435;
t346 = cos(qJ(5));
t348 = cos(qJ(2));
t472 = -pkin(8) - pkin(7);
t326 = t472 * t348;
t317 = qJD(1) * t326;
t302 = t347 * t317;
t344 = sin(qJ(2));
t325 = t472 * t344;
t316 = qJD(1) * t325;
t260 = -t316 * t343 + t302;
t395 = t347 * t348;
t310 = -t343 * t344 + t395;
t297 = t310 * qJD(1);
t413 = qJ(4) * t297;
t226 = t260 - t413;
t299 = t343 * t317;
t261 = t347 * t316 + t299;
t312 = t343 * t348 + t347 * t344;
t298 = t312 * qJD(1);
t290 = t298 * qJ(4);
t227 = -t290 + t261;
t158 = t226 * t339 + t227 * t340;
t360 = t297 * t339 + t340 * t298;
t446 = pkin(3) * t298;
t486 = t340 * t297 - t339 * t298;
t163 = pkin(4) * t360 - pkin(9) * t486 + t446;
t393 = qJD(1) * t344;
t335 = pkin(2) * t393;
t159 = t163 + t335;
t87 = t346 * t158 + t342 * t159;
t516 = t287 * t346 - t87;
t86 = -t158 * t342 + t346 * t159;
t515 = -t287 * t342 - t86;
t514 = t486 * Ifges(5,2) / 0.2e1;
t332 = pkin(2) * t347 + pkin(3);
t400 = t340 * t343;
t289 = pkin(2) * t400 + t339 * t332;
t282 = pkin(9) + t289;
t443 = -pkin(10) - t282;
t378 = qJD(5) * t443;
t407 = t486 * t342;
t386 = pkin(10) * t407;
t513 = t342 * t378 + t386 + t516;
t406 = t486 * t346;
t373 = pkin(5) * t360 - pkin(10) * t406;
t512 = t346 * t378 - t373 + t515;
t330 = pkin(3) * t339 + pkin(9);
t442 = -pkin(10) - t330;
t377 = qJD(5) * t442;
t305 = qJD(2) * pkin(2) + t316;
t256 = t305 * t343 - t302;
t219 = t256 + t413;
t210 = t339 * t219;
t255 = t347 * t305 + t299;
t218 = t255 - t290;
t142 = t218 * t340 - t210;
t85 = t346 * t142 + t342 * t163;
t511 = t342 * t377 + t386 - t85;
t84 = -t142 * t342 + t346 * t163;
t510 = t346 * t377 - t373 - t84;
t224 = t338 * t346 - t342 * t360;
t236 = qJD(5) - t486;
t225 = t338 * t342 + t346 * t360;
t428 = t225 * Ifges(6,4);
t118 = t224 * Ifges(6,2) + t236 * Ifges(6,6) + t428;
t209 = pkin(3) * t338 + t218;
t136 = t209 * t340 - t210;
t134 = -pkin(4) * t338 - t136;
t371 = mrSges(6,1) * t342 + mrSges(6,2) * t346;
t509 = t118 * t449 + t134 * t371;
t333 = -pkin(2) * t348 - pkin(1);
t324 = qJD(1) * t333;
t266 = -t297 * pkin(3) + qJD(4) + t324;
t424 = t486 * Ifges(5,4);
t508 = t266 * mrSges(5,2) + t424 / 0.2e1;
t421 = t360 * Ifges(5,4);
t507 = Ifges(5,6) * t517 + t421 / 0.2e1 + t514;
t104 = -pkin(5) * t224 + t134;
t341 = sin(qJ(6));
t345 = cos(qJ(6));
t154 = t224 * t341 + t225 * t345;
t264 = t338 * t310;
t251 = t264 * qJD(1);
t265 = t338 * t312;
t252 = t265 * qJD(1);
t191 = t251 * t339 + t340 * t252;
t185 = Ifges(7,3) * t191;
t401 = t340 * t219;
t137 = t339 * t209 + t401;
t135 = pkin(9) * t338 + t137;
t145 = -pkin(4) * t486 - pkin(9) * t360 + t266;
t76 = t135 * t346 + t145 * t342;
t60 = pkin(10) * t224 + t76;
t417 = t341 * t60;
t75 = -t135 * t342 + t346 * t145;
t59 = -pkin(10) * t225 + t75;
t51 = pkin(5) * t236 + t59;
t20 = t345 * t51 - t417;
t416 = t345 * t60;
t21 = t341 * t51 + t416;
t375 = t345 * t224 - t225 * t341;
t436 = Ifges(7,4) * t154;
t233 = qJD(6) + t236;
t458 = -t233 / 0.2e1;
t464 = -t154 / 0.2e1;
t506 = t185 + (Ifges(7,5) * t375 - Ifges(7,6) * t154) * t458 + (t154 * t21 + t20 * t375) * mrSges(7,3) - t104 * (mrSges(7,1) * t154 + mrSges(7,2) * t375) + (Ifges(7,1) * t375 - t436) * t464;
t359 = t341 * t342 - t345 * t346;
t172 = t359 * t486;
t483 = qJD(5) + qJD(6);
t262 = t483 * t359;
t504 = t262 - t172;
t311 = t341 * t346 + t342 * t345;
t171 = t311 * t486;
t263 = t483 * t311;
t503 = t263 - t171;
t485 = t340 * t226 - t227 * t339 + (t339 * t347 + t400) * t435;
t192 = t251 * t340 - t252 * t339;
t123 = -qJD(5) * t225 - t192 * t342;
t387 = qJD(5) * t346;
t388 = qJD(5) * t342;
t384 = qJD(2) * t472;
t389 = qJD(3) * t347;
t390 = qJD(3) * t343;
t196 = t298 * t384 + t305 * t389 + t317 * t390;
t128 = -qJ(4) * t252 + qJD(4) * t297 + t196;
t313 = t343 * t325;
t351 = (t395 * t472 - t313) * qJD(2) * qJD(1);
t358 = -t251 * qJ(4) - t298 * qJD(4);
t67 = t340 * t128 + (-t305 * t390 + t317 * t389 + t351 + t358) * t339;
t230 = pkin(3) * t252 + qJD(2) * t335;
t90 = pkin(4) * t191 - pkin(9) * t192 + t230;
t16 = -t135 * t388 + t145 * t387 + t342 * t90 + t346 * t67;
t12 = pkin(10) * t123 + t16;
t122 = qJD(5) * t224 + t192 * t346;
t410 = qJD(5) * t76;
t17 = -t342 * t67 + t346 * t90 - t410;
t7 = pkin(5) * t191 - pkin(10) * t122 + t17;
t3 = qJD(6) * t20 + t12 * t345 + t341 * t7;
t4 = -qJD(6) * t21 - t12 * t341 + t345 * t7;
t40 = qJD(6) * t375 + t122 * t345 + t123 * t341;
t41 = -qJD(6) * t154 - t122 * t341 + t123 * t345;
t502 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + Ifges(7,5) * t40 + Ifges(7,6) * t41;
t146 = Ifges(7,4) * t375;
t501 = -Ifges(7,2) * t154 + t146;
t366 = Ifges(6,5) * t346 - Ifges(6,6) * t342;
t437 = Ifges(6,4) * t346;
t368 = -Ifges(6,2) * t342 + t437;
t438 = Ifges(6,4) * t342;
t370 = Ifges(6,1) * t346 - t438;
t220 = Ifges(6,4) * t224;
t119 = t225 * Ifges(6,1) + t236 * Ifges(6,5) + t220;
t397 = t346 * t119;
t459 = t225 / 0.2e1;
t500 = t397 / 0.2e1 + t236 * t366 / 0.2e1 + t370 * t459 + t224 * t368 / 0.2e1 + t509;
t499 = -t266 * mrSges(5,1) - t75 * mrSges(6,1) - t20 * mrSges(7,1) + t76 * mrSges(6,2) + t21 * mrSges(7,2) + t507;
t478 = t40 / 0.2e1;
t477 = t41 / 0.2e1;
t462 = t191 / 0.2e1;
t306 = t442 * t342;
t337 = t346 * pkin(10);
t307 = t330 * t346 + t337;
t249 = t306 * t345 - t307 * t341;
t493 = qJD(6) * t249 + t341 * t510 + t345 * t511;
t250 = t306 * t341 + t307 * t345;
t492 = -qJD(6) * t250 - t341 * t511 + t345 * t510;
t267 = t443 * t342;
t268 = t282 * t346 + t337;
t203 = t267 * t341 + t268 * t345;
t491 = -qJD(6) * t203 - t341 * t513 + t345 * t512;
t202 = t267 * t345 - t268 * t341;
t490 = qJD(6) * t202 + t341 * t512 + t345 * t513;
t259 = t310 * t339 + t312 * t340;
t194 = t359 * t259;
t269 = t347 * t325 + t326 * t343;
t239 = -qJ(4) * t312 + t269;
t270 = -t347 * t326 + t313;
t240 = qJ(4) * t310 + t270;
t175 = t239 * t339 + t240 * t340;
t169 = t346 * t175;
t258 = -t340 * t310 + t312 * t339;
t276 = -t310 * pkin(3) + t333;
t173 = t258 * pkin(4) - t259 * pkin(9) + t276;
t95 = t342 * t173 + t169;
t232 = pkin(5) * t407;
t385 = pkin(5) * t388;
t487 = t385 - t232 + t485;
t484 = t287 - t158;
t481 = t17 * mrSges(6,1) - t16 * mrSges(6,2) + Ifges(6,5) * t122 + Ifges(6,6) * t123 + t502;
t480 = Ifges(7,4) * t478 + Ifges(7,2) * t477 + Ifges(7,6) * t462;
t479 = Ifges(7,1) * t478 + Ifges(7,4) * t477 + Ifges(7,5) * t462;
t69 = Ifges(7,2) * t375 + Ifges(7,6) * t233 + t436;
t476 = -t69 / 0.2e1;
t475 = t69 / 0.2e1;
t70 = Ifges(7,1) * t154 + Ifges(7,5) * t233 + t146;
t474 = -t70 / 0.2e1;
t473 = t70 / 0.2e1;
t471 = pkin(1) * mrSges(3,1);
t470 = pkin(1) * mrSges(3,2);
t468 = t122 / 0.2e1;
t467 = t123 / 0.2e1;
t466 = -t375 / 0.2e1;
t465 = t375 / 0.2e1;
t463 = t154 / 0.2e1;
t461 = -t224 / 0.2e1;
t460 = -t225 / 0.2e1;
t457 = t233 / 0.2e1;
t456 = -t236 / 0.2e1;
t452 = t297 / 0.2e1;
t451 = t298 / 0.2e1;
t450 = -t338 / 0.2e1;
t448 = t346 / 0.2e1;
t447 = m(4) * t324;
t445 = pkin(5) * t346;
t444 = t75 * mrSges(6,3);
t441 = mrSges(4,3) * t297;
t440 = Ifges(3,4) * t344;
t439 = Ifges(4,4) * t298;
t434 = t375 * Ifges(7,6);
t433 = t154 * Ifges(7,5);
t432 = t16 * t346;
t431 = t17 * t342;
t174 = -t340 * t239 + t240 * t339;
t197 = -t256 * qJD(3) + t351;
t66 = t128 * t339 - t340 * (t197 + t358);
t430 = t174 * t66;
t429 = t224 * Ifges(6,6);
t427 = t225 * Ifges(6,5);
t426 = t233 * Ifges(7,3);
t425 = t236 * Ifges(6,3);
t422 = t360 * Ifges(5,1);
t420 = t298 * mrSges(4,3);
t419 = t338 * Ifges(5,5);
t415 = Ifges(3,5) * qJD(2);
t414 = Ifges(3,6) * qJD(2);
t412 = qJD(2) * mrSges(3,1);
t411 = qJD(2) * mrSges(3,2);
t409 = t137 * t360;
t405 = t259 * t342;
t394 = mrSges(5,1) * t338 + mrSges(6,1) * t224 - mrSges(6,2) * t225 - mrSges(5,3) * t360;
t392 = qJD(1) * t348;
t391 = qJD(2) * t344;
t331 = -pkin(3) * t340 - pkin(4);
t383 = t415 / 0.2e1;
t382 = -t414 / 0.2e1;
t246 = pkin(2) * t391 + pkin(3) * t265;
t198 = t264 * t339 + t340 * t265;
t199 = t264 * t340 - t265 * t339;
t101 = pkin(4) * t198 - pkin(9) * t199 + t246;
t318 = t344 * t384;
t319 = t348 * t384;
t205 = t347 * t318 + t343 * t319 + t325 * t389 + t326 * t390;
t155 = -qJ(4) * t265 + qJD(4) * t310 + t205;
t206 = -qJD(3) * t270 - t318 * t343 + t347 * t319;
t156 = -qJ(4) * t264 - qJD(4) * t312 + t206;
t80 = t155 * t340 + t156 * t339;
t379 = t346 * t101 - t342 * t80;
t376 = t191 * mrSges(5,1) + t192 * mrSges(5,2);
t79 = t155 * t339 - t340 * t156;
t94 = t346 * t173 - t175 * t342;
t141 = t218 * t339 + t401;
t288 = -pkin(2) * t402 + t332 * t340;
t281 = -pkin(4) - t288;
t372 = mrSges(6,1) * t346 - mrSges(6,2) * t342;
t369 = Ifges(6,1) * t342 + t437;
t367 = Ifges(6,2) * t346 + t438;
t365 = Ifges(6,5) * t342 + Ifges(6,6) * t346;
t364 = -t16 * t342 - t17 * t346;
t73 = pkin(5) * t258 - t259 * t337 + t94;
t77 = -pkin(10) * t405 + t95;
t30 = -t341 * t77 + t345 * t73;
t31 = t341 * t73 + t345 * t77;
t363 = -t342 * t76 - t346 * t75;
t362 = t342 * t75 - t346 * t76;
t165 = -mrSges(6,2) * t236 + mrSges(6,3) * t224;
t166 = mrSges(6,1) * t236 - mrSges(6,3) * t225;
t361 = t346 * t165 - t342 * t166;
t357 = t199 * t342 + t259 * t387;
t22 = t342 * t101 + t173 * t387 - t175 * t388 + t346 * t80;
t350 = m(6) * (qJD(5) * t363 - t431 + t432);
t117 = t425 + t427 + t429;
t178 = t419 + t422 + t424;
t237 = t297 * Ifges(4,2) + t338 * Ifges(4,6) + t439;
t293 = Ifges(4,4) * t297;
t238 = t298 * Ifges(4,1) + t338 * Ifges(4,5) + t293;
t37 = -pkin(5) * t123 + t66;
t44 = t122 * Ifges(6,4) + t123 * Ifges(6,2) + t191 * Ifges(6,6);
t45 = t122 * Ifges(6,1) + t123 * Ifges(6,4) + t191 * Ifges(6,5);
t68 = t426 + t433 + t434;
t349 = (t406 * t75 + t407 * t76 + t432) * mrSges(6,3) + (t20 * t504 - t21 * t503 - t3 * t359 - t4 * t311) * mrSges(7,3) + (mrSges(7,1) * t503 - mrSges(7,2) * t504) * t104 + (-t372 - mrSges(5,1)) * t66 + (-Ifges(7,5) * t172 - Ifges(7,6) * t171) * t458 + (-Ifges(7,1) * t172 - Ifges(7,4) * t171) * t464 - t324 * (mrSges(4,1) * t298 + mrSges(4,2) * t297) + (-Ifges(7,4) * t172 - Ifges(7,2) * t171) * t466 + (t136 * mrSges(5,3) + Ifges(5,5) * t450 + t366 * t456 + t368 * t461 + t370 * t460 - t508 - t509) * t486 + t500 * qJD(5) + t311 * t479 + t367 * t467 + t369 * t468 - t262 * t473 - t172 * t474 - t263 * t475 - t171 * t476 + t44 * t448 + (Ifges(4,5) * t297 - Ifges(4,6) * t298) * t450 + t237 * t451 + t255 * t441 + (-Ifges(7,5) * t262 - Ifges(7,6) * t263) * t457 + (-Ifges(7,1) * t262 - Ifges(7,4) * t263) * t463 + (-Ifges(7,4) * t262 - Ifges(7,2) * t263) * t465 - (-Ifges(4,2) * t298 + t238 + t293) * t297 / 0.2e1 - (t397 + t178) * t486 / 0.2e1 - (Ifges(5,1) * t486 + t117 - t421 + t68) * t360 / 0.2e1 + (Ifges(7,1) * t311 - Ifges(7,4) * t359) * t478 + (Ifges(7,4) * t311 - Ifges(7,2) * t359) * t477 + t37 * (mrSges(7,1) * t359 + mrSges(7,2) * t311) + (Ifges(7,5) * t311 - Ifges(7,6) * t359 + t365) * t462 - t359 * t480 + t342 * t45 / 0.2e1 - Ifges(4,6) * t252 + Ifges(4,5) * t251 + Ifges(5,5) * t192 - t196 * mrSges(4,2) + t197 * mrSges(4,1) - Ifges(5,6) * t191 - t298 * (Ifges(4,1) * t297 - t439) / 0.2e1 + (Ifges(6,5) * t460 + Ifges(7,5) * t464 - Ifges(5,6) * t450 + Ifges(6,6) * t461 + Ifges(7,6) * t466 + Ifges(6,3) * t456 + Ifges(7,3) * t458 + t499 + t514) * t360 - t67 * mrSges(5,2);
t334 = Ifges(3,4) * t392;
t322 = mrSges(3,3) * t392 - t411;
t321 = -mrSges(3,3) * t393 + t412;
t320 = t331 - t445;
t296 = Ifges(3,1) * t393 + t334 + t415;
t295 = t414 + (t348 * Ifges(3,2) + t440) * qJD(1);
t275 = t281 - t445;
t274 = mrSges(4,1) * t338 - t420;
t273 = -mrSges(4,2) * t338 + t441;
t272 = t335 + t446;
t254 = -mrSges(4,1) * t297 + mrSges(4,2) * t298;
t228 = -mrSges(5,2) * t338 + mrSges(5,3) * t486;
t193 = t311 * t259;
t186 = Ifges(6,3) * t191;
t183 = -mrSges(5,1) * t486 + mrSges(5,2) * t360;
t131 = pkin(5) * t405 + t174;
t111 = mrSges(7,1) * t233 - mrSges(7,3) * t154;
t110 = -mrSges(7,2) * t233 + mrSges(7,3) * t375;
t109 = t141 + t232;
t83 = -mrSges(7,1) * t375 + mrSges(7,2) * t154;
t82 = -mrSges(6,2) * t191 + mrSges(6,3) * t123;
t81 = mrSges(6,1) * t191 - mrSges(6,3) * t122;
t61 = -mrSges(6,1) * t123 + mrSges(6,2) * t122;
t53 = t194 * t483 - t311 * t199;
t52 = -t199 * t359 - t259 * t263;
t50 = pkin(5) * t357 + t79;
t33 = -mrSges(7,2) * t191 + mrSges(7,3) * t41;
t32 = mrSges(7,1) * t191 - mrSges(7,3) * t40;
t25 = t345 * t59 - t417;
t24 = -t341 * t59 - t416;
t23 = -qJD(5) * t95 + t379;
t18 = -pkin(10) * t357 + t22;
t14 = -t199 * t337 + pkin(5) * t198 + (-t169 + (pkin(10) * t259 - t173) * t342) * qJD(5) + t379;
t13 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t6 = -qJD(6) * t31 + t14 * t345 - t18 * t341;
t5 = qJD(6) * t30 + t14 * t341 + t18 * t345;
t1 = [(-pkin(7) * t321 + t296 / 0.2e1 + t383 + (-0.2e1 * t470 + 0.3e1 / 0.2e1 * Ifges(3,4) * t348) * qJD(1)) * t348 * qJD(2) + (t196 * t310 - t197 * t312 - t251 * t269 - t252 * t270 - t255 * t264 - t256 * t265) * mrSges(4,3) + (-t310 * t252 - t265 * t452) * Ifges(4,2) + (t310 * t251 - t252 * t312 + t264 * t452 - t265 * t451) * Ifges(4,4) + t333 * (mrSges(4,1) * t252 + mrSges(4,2) * t251) - t394 * t79 + t276 * t376 + (-t193 * t3 + t194 * t4 - t20 * t52 + t21 * t53) * mrSges(7,3) + (-Ifges(7,4) * t194 - Ifges(7,2) * t193) * t477 + (-Ifges(7,1) * t194 - Ifges(7,4) * t193) * t478 + (-Ifges(7,5) * t194 - Ifges(7,6) * t193) * t462 + t37 * (mrSges(7,1) * t193 - mrSges(7,2) * t194) + m(7) * (t104 * t50 + t131 * t37 + t20 * t6 + t21 * t5 + t3 * t31 + t30 * t4) + (t230 * mrSges(5,1) - Ifges(5,4) * t192 + t186 / 0.2e1 + t185 / 0.2e1 - t67 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) + Ifges(7,3) / 0.2e1) * t191 + t481) * t258 + m(4) * (t196 * t270 + t197 * t269 + t205 * t256 + t206 * t255) + (t419 / 0.2e1 + t422 / 0.2e1 + t178 / 0.2e1 + t363 * mrSges(6,3) + t500 + t508) * t199 + (t426 / 0.2e1 + t68 / 0.2e1 + t117 / 0.2e1 + t434 / 0.2e1 + t433 / 0.2e1 + t425 / 0.2e1 + t429 / 0.2e1 + t427 / 0.2e1 - t499 - t507) * t198 - t194 * t479 - t193 * t480 + (Ifges(7,1) * t52 + Ifges(7,4) * t53) * t463 + (Ifges(7,4) * t52 + Ifges(7,2) * t53) * t465 + t52 * t473 + t53 * t475 + (Ifges(7,5) * t52 + Ifges(7,6) * t53) * t457 + t324 * (mrSges(4,1) * t265 + mrSges(4,2) * t264) + (-t136 * t199 - t137 * t198 + t174 * t192 - t175 * t191) * mrSges(5,3) - t265 * t237 / 0.2e1 + t205 * t273 + t206 * t274 + t264 * t238 / 0.2e1 + t246 * t183 + t80 * t228 + m(6) * (t134 * t79 + t16 * t95 + t17 * t94 + t22 * t76 + t23 * t75 + t430) + m(5) * (-t136 * t79 + t137 * t80 + t175 * t67 + t230 * t276 + t246 * t266 + t430) + (t251 * t312 + t264 * t451) * Ifges(4,1) + (t366 * t462 + t370 * t468 + t368 * t467 + t230 * mrSges(5,2) + Ifges(5,1) * t192 - Ifges(5,4) * t191 + t44 * t449 + t45 * t448 + (mrSges(5,3) + t371) * t66 + t364 * mrSges(6,3) + (-t346 * t118 / 0.2e1 + t119 * t449 + t365 * t456 + t367 * t461 + t369 * t460 + t134 * t372 + t362 * mrSges(6,3)) * qJD(5)) * t259 + (-pkin(7) * t322 - t295 / 0.2e1 + t382 + (-0.2e1 * t471 - 0.3e1 / 0.2e1 * t440 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t348) * qJD(1) + (t254 + qJD(1) * (-mrSges(4,1) * t310 + mrSges(4,2) * t312) + 0.2e1 * t447) * pkin(2)) * t391 + t30 * t32 + t31 * t33 + t50 * t83 + t94 * t81 + t95 * t82 + t104 * (-mrSges(7,1) * t53 + mrSges(7,2) * t52) + t5 * t110 + t6 * t111 + t131 * t13 + t22 * t165 + t23 * t166 + t174 * t61 + (Ifges(4,5) * t264 - Ifges(4,6) * t265) * t517; ((t273 * t347 - t274 * t343) * qJD(3) + (-t251 * t347 - t252 * t343) * mrSges(4,3)) * pkin(2) + ((t196 * t343 + t197 * t347 + (-t255 * t343 + t256 * t347) * qJD(3)) * pkin(2) - t255 * t260 - t256 * t261) * m(4) + (-t289 * t191 - t288 * t192 + t409) * mrSges(5,3) + t349 + (t485 * t134 + t281 * t66 + t515 * t75 + t516 * t76) * m(6) + t490 * t110 + (t104 * t487 + t20 * t491 + t202 * t4 + t203 * t3 + t21 * t490 + t275 * t37) * m(7) + t491 * t111 + t256 * t420 + t282 * t350 + t484 * t228 - t485 * t394 + (-t136 * t485 + t137 * t484 - t266 * t272 - t288 * t66 + t289 * t67) * m(5) + t487 * t83 + t281 * t61 - t272 * t183 - t261 * t273 - t260 * t274 + t275 * t13 + t202 * t32 + t203 * t33 + (-t287 * t166 + (-qJD(5) * t165 - t81) * t282 + (-t17 - t410) * mrSges(6,3)) * t342 + (t287 * t165 + t282 * t82 + (-t166 * t282 - t444) * qJD(5)) * t346 + ((t383 - t296 / 0.2e1 - t334 / 0.2e1 + qJD(1) * t470 + (t321 - t412) * pkin(7)) * t348 + (t382 + t295 / 0.2e1 + (t471 + t440 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t348) * qJD(1) + (t322 + t411) * pkin(7) + (-t254 - t447) * pkin(2)) * t344) * qJD(1) - t87 * t165 - t86 * t166; (t409 + (-t191 * t339 - t192 * t340) * pkin(3)) * mrSges(5,3) + t394 * t141 + t349 + (-t342 * t81 + t346 * t82) * t330 + t330 * t350 + t492 * t111 + t493 * t110 + t331 * t61 + t320 * t13 - t255 * t273 + t249 * t32 + t250 * t33 - t142 * t228 + (t274 + t420) * t256 - mrSges(6,3) * t431 + ((-t166 * t330 - t444) * t346 + (-t76 * mrSges(6,3) + pkin(5) * t83 - t165 * t330) * t342) * qJD(5) - t183 * t446 - t109 * t83 - t85 * t165 - t84 * t166 + (t249 * t4 + t250 * t3 + t320 * t37 + t493 * t21 + t492 * t20 + (-t109 + t385) * t104) * m(7) + (-t134 * t141 + t331 * t66 - t75 * t84 - t76 * t85) * m(6) + ((t339 * t67 - t340 * t66) * pkin(3) + t136 * t141 - t137 * t142 - t266 * t446) * m(5); -t359 * t32 + t311 * t33 + t342 * t82 + t346 * t81 - t503 * t111 - t504 * t110 + t361 * qJD(5) + (-t228 - t361) * t486 + (-t83 + t394) * t360 + t376 + (-t104 * t360 - t20 * t503 - t21 * t504 + t3 * t311 - t359 * t4) * m(7) + (-t134 * t360 - t236 * t362 - t364) * m(6) + (t136 * t360 - t137 * t486 + t230) * m(5); -m(7) * (t20 * t24 + t21 * t25) + t481 + (-Ifges(6,2) * t225 + t119 + t220) * t461 + t375 * t474 + (Ifges(6,1) * t224 - t428) * t460 + t501 * t466 + (Ifges(6,5) * t224 - Ifges(6,6) * t225) * t456 + t118 * t459 - t154 * t476 + (t224 * t75 + t225 * t76) * mrSges(6,3) + t186 - t134 * (mrSges(6,1) * t225 + mrSges(6,2) * t224) - t25 * t110 - t24 * t111 - t75 * t165 + t76 * t166 + (-t225 * t83 + t345 * t32 + t341 * t33 + (t110 * t345 - t111 * t341) * qJD(6) + (-t104 * t225 + t3 * t341 + t345 * t4 + (-t20 * t341 + t21 * t345) * qJD(6)) * m(7)) * pkin(5) + t506; t69 * t463 - t20 * t110 + t21 * t111 + (t501 + t70) * t466 + t502 + t506;];
tauc  = t1(:);
