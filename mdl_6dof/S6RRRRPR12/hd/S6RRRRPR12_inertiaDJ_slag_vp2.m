% Calculate time derivative of joint inertia matrix for
% S6RRRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR12_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR12_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:20
% EndTime: 2019-03-09 23:29:49
% DurationCPUTime: 12.42s
% Computational Cost: add. (28206->952), mult. (84207->1390), div. (0->0), fcn. (86174->14), ass. (0->390)
t340 = sin(pkin(6));
t490 = 0.2e1 * t340;
t350 = cos(qJ(2));
t342 = cos(pkin(6));
t440 = pkin(1) * t342;
t331 = t350 * t440;
t325 = qJD(2) * t331;
t346 = sin(qJ(2));
t341 = cos(pkin(7));
t364 = t340 * (-pkin(10) * t341 - pkin(9));
t357 = t346 * t364;
t243 = pkin(2) * t342 + t331 + t357;
t417 = t243 * t341;
t489 = qJD(2) * t357 + qJD(3) * t417 + t325;
t338 = sin(pkin(13));
t348 = cos(qJ(4));
t344 = sin(qJ(4));
t419 = cos(pkin(13));
t367 = t419 * t344;
t297 = t338 * t348 + t367;
t343 = sin(qJ(6));
t394 = qJD(6) * t343;
t412 = t338 * t344;
t352 = t348 * t419 - t412;
t289 = t352 * qJD(4);
t347 = cos(qJ(6));
t402 = t347 * t289;
t353 = t297 * t394 - t402;
t445 = t343 / 0.2e1;
t444 = t344 / 0.2e1;
t442 = t347 / 0.2e1;
t441 = t348 / 0.2e1;
t370 = -t394 / 0.2e1;
t349 = cos(qJ(3));
t403 = t346 * t349;
t345 = sin(qJ(3));
t404 = t345 * t350;
t355 = t341 * t404 + t403;
t339 = sin(pkin(7));
t398 = qJD(3) * t345;
t380 = t339 * t398;
t182 = t342 * t380 + (t355 * qJD(3) + (t341 * t403 + t404) * qJD(2)) * t340;
t397 = qJD(3) * t349;
t379 = t339 * t397;
t401 = t349 * t350;
t405 = t345 * t346;
t486 = t341 * t401 - t405;
t183 = t342 * t379 + (t486 * qJD(3) + (-t341 * t405 + t401) * qJD(2)) * t340;
t411 = t339 * t345;
t234 = t340 * t355 + t342 * t411;
t409 = t340 * t350;
t287 = -t339 * t409 + t342 * t341;
t189 = t234 * t348 + t287 * t344;
t400 = qJD(2) * t340;
t382 = t346 * t400;
t365 = t339 * t382;
t124 = -qJD(4) * t189 - t183 * t344 + t348 * t365;
t188 = -t234 * t344 + t287 * t348;
t125 = qJD(4) * t188 + t183 * t348 + t344 * t365;
t67 = -t124 * t419 + t125 * t338;
t68 = t338 * t124 + t125 * t419;
t28 = Ifges(6,5) * t68 - Ifges(6,6) * t67 + Ifges(6,3) * t182;
t48 = Ifges(5,5) * t125 + Ifges(5,6) * t124 + Ifges(5,3) * t182;
t488 = t28 + t48;
t290 = t341 * t348 - t344 * t411;
t246 = qJD(4) * t290 + t348 * t379;
t291 = t341 * t344 + t348 * t411;
t247 = -qJD(4) * t291 - t344 * t379;
t173 = t246 * t338 - t247 * t419;
t174 = t246 * t419 + t338 * t247;
t106 = Ifges(6,5) * t174 - Ifges(6,6) * t173 + Ifges(6,3) * t380;
t163 = Ifges(5,5) * t246 + Ifges(5,6) * t247 + Ifges(5,3) * t380;
t487 = -t106 - t163;
t438 = pkin(10) * t339;
t260 = (-pkin(2) * t350 - t346 * t438 - pkin(1)) * t340;
t177 = -t243 * t339 + t341 * t260;
t410 = t339 * t349;
t233 = -t486 * t340 - t342 * t410;
t129 = pkin(3) * t233 - pkin(11) * t234 + t177;
t330 = t346 * t440;
t295 = pkin(9) * t409 + t330;
t229 = (t339 * t342 + t341 * t409) * pkin(10) + t295;
t408 = t341 * t345;
t140 = t349 * t229 + t243 * t408 + t260 * t411;
t134 = pkin(11) * t287 + t140;
t74 = t344 * t129 + t348 * t134;
t294 = pkin(2) * t408 + pkin(10) * t410;
t270 = pkin(11) * t341 + t294;
t271 = (-pkin(3) * t349 - pkin(11) * t345 - pkin(2)) * t339;
t201 = t348 * t270 + t344 * t271;
t326 = pkin(10) * t411;
t407 = t341 * t349;
t292 = pkin(2) * t407 - t326;
t393 = qJD(6) * t347;
t354 = t343 * t289 + t297 * t393;
t335 = -pkin(4) * t348 - pkin(3);
t226 = -pkin(5) * t352 - pkin(12) * t297 + t335;
t435 = -qJ(5) - pkin(11);
t313 = t435 * t348;
t251 = -t313 * t419 + t412 * t435;
t166 = t226 * t347 - t251 * t343;
t368 = qJD(4) * t435;
t286 = qJD(5) * t348 + t344 * t368;
t351 = -qJD(5) * t344 + t348 * t368;
t206 = t286 * t419 + t338 * t351;
t288 = t297 * qJD(4);
t396 = qJD(4) * t344;
t390 = pkin(4) * t396;
t207 = pkin(5) * t288 - pkin(12) * t289 + t390;
t90 = qJD(6) * t166 + t206 * t347 + t207 * t343;
t167 = t226 * t343 + t251 * t347;
t91 = -qJD(6) * t167 - t206 * t343 + t207 * t347;
t485 = -t343 * t91 + t347 * t90;
t399 = qJD(3) * t339;
t277 = (pkin(3) * t345 - pkin(11) * t349) * t399;
t278 = t292 * qJD(3);
t147 = -t201 * qJD(4) + t348 * t277 - t278 * t344;
t102 = pkin(4) * t380 - qJ(5) * t246 - qJD(5) * t291 + t147;
t395 = qJD(4) * t348;
t146 = -t270 * t396 + t271 * t395 + t344 * t277 + t348 * t278;
t112 = qJ(5) * t247 + qJD(5) * t290 + t146;
t57 = t338 * t102 + t419 * t112;
t55 = pkin(12) * t380 + t57;
t200 = -t270 * t344 + t348 * t271;
t162 = -pkin(4) * t410 - qJ(5) * t291 + t200;
t175 = qJ(5) * t290 + t201;
t111 = t338 * t162 + t419 * t175;
t101 = -pkin(12) * t410 + t111;
t210 = -t290 * t419 + t291 * t338;
t211 = t338 * t290 + t291 * t419;
t269 = t326 + (-pkin(2) * t349 - pkin(3)) * t341;
t215 = -pkin(4) * t290 + t269;
t132 = pkin(5) * t210 - pkin(12) * t211 + t215;
t58 = -t101 * t343 + t132 * t347;
t279 = t294 * qJD(3);
t199 = -pkin(4) * t247 + t279;
t85 = pkin(5) * t173 - pkin(12) * t174 + t199;
t18 = qJD(6) * t58 + t343 * t85 + t347 * t55;
t59 = t101 * t347 + t132 * t343;
t19 = -qJD(6) * t59 - t343 * t55 + t347 * t85;
t484 = t18 * t347 - t19 * t343;
t73 = t348 * t129 - t134 * t344;
t45 = pkin(4) * t233 - qJ(5) * t189 + t73;
t53 = qJ(5) * t188 + t74;
t25 = t338 * t45 + t419 * t53;
t23 = pkin(12) * t233 + t25;
t135 = -t188 * t419 + t189 * t338;
t136 = t338 * t188 + t189 * t419;
t139 = -t345 * t229 + t349 * (t260 * t339 + t417);
t133 = -pkin(3) * t287 - t139;
t89 = -pkin(4) * t188 + t133;
t39 = pkin(5) * t135 - pkin(12) * t136 + t89;
t10 = -t23 * t343 + t347 * t39;
t245 = (t350 * t364 - t330) * qJD(2);
t264 = (pkin(2) * t346 - t350 * t438) * t400;
t366 = -t229 * t397 - t260 * t380 - t345 * t489;
t78 = -t245 * t407 + (-pkin(3) * t382 - t264 * t349) * t339 - t366;
t40 = -pkin(4) * t124 + t78;
t15 = pkin(5) * t67 - pkin(12) * t68 + t40;
t81 = -t229 * t398 + t245 * t408 + t260 * t379 + t264 * t411 + t349 * t489;
t77 = pkin(11) * t365 + t81;
t185 = -t245 * t339 + t341 * t264;
t88 = pkin(3) * t182 - pkin(11) * t183 + t185;
t27 = -t74 * qJD(4) - t344 * t77 + t348 * t88;
t14 = pkin(4) * t182 - qJ(5) * t125 - qJD(5) * t189 + t27;
t26 = t129 * t395 - t134 * t396 + t344 * t88 + t348 * t77;
t17 = qJ(5) * t124 + qJD(5) * t188 + t26;
t6 = t338 * t14 + t419 * t17;
t4 = pkin(12) * t182 + t6;
t1 = qJD(6) * t10 + t15 * t343 + t347 * t4;
t11 = t23 * t347 + t343 * t39;
t2 = -qJD(6) * t11 + t15 * t347 - t343 * t4;
t483 = t1 * t347 - t2 * t343;
t482 = -m(5) * pkin(3) - mrSges(5,1) * t348 + mrSges(5,2) * t344;
t481 = 2 * m(4);
t480 = 0.2e1 * m(5);
t479 = 2 * m(6);
t478 = 2 * m(7);
t477 = -2 * mrSges(3,3);
t476 = -2 * mrSges(4,3);
t475 = -2 * mrSges(6,3);
t205 = t286 * t338 - t351 * t419;
t474 = 0.2e1 * t205;
t250 = -t313 * t338 - t367 * t435;
t473 = 0.2e1 * t250;
t471 = m(6) * pkin(4);
t30 = Ifges(6,1) * t68 - Ifges(6,4) * t67 + Ifges(6,5) * t182;
t470 = t30 / 0.2e1;
t71 = Ifges(6,1) * t136 - Ifges(6,4) * t135 + Ifges(6,5) * t233;
t469 = t71 / 0.2e1;
t92 = -t136 * t343 + t233 * t347;
t468 = t92 / 0.2e1;
t93 = t136 * t347 + t233 * t343;
t467 = t93 / 0.2e1;
t108 = Ifges(6,1) * t174 - Ifges(6,4) * t173 + Ifges(6,5) * t380;
t466 = t108 / 0.2e1;
t150 = Ifges(6,1) * t211 - Ifges(6,4) * t210 - Ifges(6,5) * t410;
t465 = t150 / 0.2e1;
t190 = -t211 * t343 - t347 * t410;
t464 = t190 / 0.2e1;
t356 = -t211 * t347 + t343 * t410;
t463 = -t356 / 0.2e1;
t429 = Ifges(7,4) * t347;
t361 = -Ifges(7,2) * t343 + t429;
t195 = -Ifges(7,6) * t352 + t297 * t361;
t462 = t195 / 0.2e1;
t430 = Ifges(7,4) * t343;
t362 = Ifges(7,1) * t347 - t430;
t196 = -Ifges(7,5) * t352 + t297 * t362;
t461 = t196 / 0.2e1;
t219 = Ifges(6,1) * t289 - Ifges(6,4) * t288;
t460 = t219 / 0.2e1;
t238 = Ifges(6,1) * t297 + Ifges(6,4) * t352;
t459 = t238 / 0.2e1;
t336 = Ifges(7,5) * t393;
t458 = Ifges(7,6) * t370 + t336 / 0.2e1;
t305 = t361 * qJD(6);
t457 = t305 / 0.2e1;
t431 = Ifges(5,4) * t348;
t306 = (-Ifges(5,2) * t344 + t431) * qJD(4);
t456 = t306 / 0.2e1;
t307 = t362 * qJD(6);
t455 = t307 / 0.2e1;
t432 = Ifges(5,4) * t344;
t308 = (Ifges(5,1) * t348 - t432) * qJD(4);
t454 = t308 / 0.2e1;
t453 = Ifges(7,5) * t445 + Ifges(7,6) * t442;
t316 = Ifges(7,2) * t347 + t430;
t451 = t316 / 0.2e1;
t317 = Ifges(5,2) * t348 + t432;
t450 = t317 / 0.2e1;
t318 = Ifges(7,1) * t343 + t429;
t449 = t318 / 0.2e1;
t319 = Ifges(5,1) * t344 + t431;
t448 = t319 / 0.2e1;
t447 = t341 / 0.2e1;
t446 = -t343 / 0.2e1;
t443 = -t347 / 0.2e1;
t439 = pkin(4) * t338;
t434 = Ifges(4,4) * t345;
t433 = Ifges(4,4) * t349;
t428 = Ifges(7,6) * t343;
t426 = t185 * mrSges(4,1);
t425 = t185 * mrSges(4,2);
t280 = -pkin(9) * t382 + t325;
t423 = t280 * mrSges(3,2);
t281 = t295 * qJD(2);
t422 = t281 * mrSges(3,1);
t418 = t205 * t250;
t416 = t297 * t343;
t415 = t297 * t347;
t333 = pkin(12) + t439;
t414 = t333 * t343;
t413 = t333 * t347;
t217 = Ifges(6,5) * t289 - Ifges(6,6) * t288;
t34 = qJD(6) * t92 + t182 * t343 + t347 * t68;
t35 = -qJD(6) * t93 + t182 * t347 - t343 * t68;
t7 = Ifges(7,5) * t34 + Ifges(7,6) * t35 + Ifges(7,3) * t67;
t29 = Ifges(6,4) * t68 - Ifges(6,2) * t67 + Ifges(6,6) * t182;
t391 = t7 / 0.2e1 - t29 / 0.2e1;
t36 = Ifges(7,5) * t93 + Ifges(7,6) * t92 + Ifges(7,3) * t135;
t70 = Ifges(6,4) * t136 - Ifges(6,2) * t135 + Ifges(6,6) * t233;
t389 = -t70 / 0.2e1 + t36 / 0.2e1;
t388 = mrSges(7,3) * t394;
t387 = mrSges(7,3) * t393;
t149 = Ifges(6,4) * t211 - Ifges(6,2) * t210 - Ifges(6,6) * t410;
t97 = -Ifges(7,5) * t356 + Ifges(7,6) * t190 + Ifges(7,3) * t210;
t385 = -t149 / 0.2e1 + t97 / 0.2e1;
t107 = Ifges(6,4) * t174 - Ifges(6,2) * t173 + Ifges(6,6) * t380;
t117 = qJD(6) * t190 + t174 * t347 + t343 * t380;
t118 = qJD(6) * t356 - t174 * t343 + t347 * t380;
t42 = Ifges(7,5) * t117 + Ifges(7,6) * t118 + Ifges(7,3) * t173;
t384 = t42 / 0.2e1 - t107 / 0.2e1;
t121 = Ifges(4,5) * t183 - Ifges(4,6) * t182 + Ifges(4,3) * t365;
t383 = t419 * pkin(4);
t378 = t333 * t394;
t377 = t333 * t393;
t141 = -Ifges(7,5) * t353 - Ifges(7,6) * t354 + Ifges(7,3) * t288;
t218 = Ifges(6,4) * t289 - Ifges(6,2) * t288;
t374 = -t218 / 0.2e1 + t141 / 0.2e1;
t194 = -Ifges(7,3) * t352 + (Ifges(7,5) * t347 - t428) * t297;
t237 = Ifges(6,4) * t297 + Ifges(6,2) * t352;
t373 = -t237 / 0.2e1 + t194 / 0.2e1;
t304 = Ifges(5,5) * t395 - Ifges(5,6) * t396;
t372 = t304 / 0.2e1 + t217 / 0.2e1;
t371 = Ifges(5,5) * t444 + Ifges(5,6) * t441 + Ifges(6,5) * t297 / 0.2e1 + Ifges(6,6) * t352 / 0.2e1;
t31 = t67 * mrSges(6,1) + t68 * mrSges(6,2);
t369 = t393 / 0.2e1;
t115 = t173 * mrSges(6,1) + t174 * mrSges(6,2);
t216 = t288 * mrSges(6,1) + t289 * mrSges(6,2);
t311 = -mrSges(7,1) * t347 + mrSges(7,2) * t343;
t363 = mrSges(7,1) * t343 + mrSges(7,2) * t347;
t360 = t26 * t348 - t27 * t344;
t359 = t146 * t348 - t147 * t344;
t5 = t14 * t419 - t338 * t17;
t24 = -t338 * t53 + t419 * t45;
t56 = t102 * t419 - t338 * t112;
t110 = t162 * t419 - t338 * t175;
t334 = -t383 - pkin(5);
t324 = Ifges(3,5) * t350 * t400;
t323 = Ifges(4,5) * t379;
t302 = (mrSges(5,1) * t344 + mrSges(5,2) * t348) * qJD(4);
t301 = t363 * qJD(6);
t300 = -mrSges(4,2) * t341 + mrSges(4,3) * t410;
t299 = mrSges(4,1) * t341 - mrSges(4,3) * t411;
t293 = -pkin(9) * t340 * t346 + t331;
t276 = (Ifges(4,1) * t349 - t434) * t399;
t275 = (-Ifges(4,2) * t345 + t433) * t399;
t274 = -Ifges(4,6) * t380 + t323;
t273 = (mrSges(4,1) * t345 + mrSges(4,2) * t349) * t399;
t268 = Ifges(4,5) * t341 + (Ifges(4,1) * t345 + t433) * t339;
t267 = Ifges(4,6) * t341 + (Ifges(4,2) * t349 + t434) * t339;
t254 = -mrSges(5,1) * t410 - mrSges(5,3) * t291;
t253 = mrSges(5,2) * t410 + mrSges(5,3) * t290;
t235 = -mrSges(6,1) * t352 + mrSges(6,2) * t297;
t228 = -mrSges(7,1) * t352 - mrSges(7,3) * t415;
t227 = mrSges(7,2) * t352 - mrSges(7,3) * t416;
t221 = -mrSges(5,1) * t290 + mrSges(5,2) * t291;
t220 = t363 * t297;
t209 = -mrSges(5,2) * t380 + mrSges(5,3) * t247;
t208 = mrSges(5,1) * t380 - mrSges(5,3) * t246;
t204 = Ifges(5,1) * t291 + Ifges(5,4) * t290 - Ifges(5,5) * t410;
t203 = Ifges(5,4) * t291 + Ifges(5,2) * t290 - Ifges(5,6) * t410;
t202 = Ifges(5,5) * t291 + Ifges(5,6) * t290 - Ifges(5,3) * t410;
t198 = -mrSges(6,1) * t410 - mrSges(6,3) * t211;
t197 = mrSges(6,2) * t410 - mrSges(6,3) * t210;
t193 = mrSges(4,1) * t287 - mrSges(4,3) * t234;
t192 = -mrSges(4,2) * t287 - mrSges(4,3) * t233;
t187 = -mrSges(7,2) * t288 - mrSges(7,3) * t354;
t186 = mrSges(7,1) * t288 + mrSges(7,3) * t353;
t176 = -mrSges(5,1) * t247 + mrSges(5,2) * t246;
t165 = Ifges(5,1) * t246 + Ifges(5,4) * t247 + Ifges(5,5) * t380;
t164 = Ifges(5,4) * t246 + Ifges(5,2) * t247 + Ifges(5,6) * t380;
t161 = mrSges(4,1) * t365 - mrSges(4,3) * t183;
t160 = -mrSges(4,2) * t365 - mrSges(4,3) * t182;
t158 = mrSges(6,1) * t380 - mrSges(6,3) * t174;
t157 = -mrSges(6,2) * t380 - mrSges(6,3) * t173;
t156 = mrSges(7,1) * t354 - mrSges(7,2) * t353;
t155 = mrSges(6,1) * t210 + mrSges(6,2) * t211;
t154 = Ifges(4,1) * t234 - Ifges(4,4) * t233 + Ifges(4,5) * t287;
t153 = Ifges(4,4) * t234 - Ifges(4,2) * t233 + Ifges(4,6) * t287;
t152 = mrSges(5,1) * t233 - mrSges(5,3) * t189;
t151 = -mrSges(5,2) * t233 + mrSges(5,3) * t188;
t148 = Ifges(6,5) * t211 - Ifges(6,6) * t210 - Ifges(6,3) * t410;
t145 = mrSges(7,1) * t210 + mrSges(7,3) * t356;
t144 = -mrSges(7,2) * t210 + mrSges(7,3) * t190;
t143 = -Ifges(7,1) * t353 - Ifges(7,4) * t354 + Ifges(7,5) * t288;
t142 = -Ifges(7,4) * t353 - Ifges(7,2) * t354 + Ifges(7,6) * t288;
t138 = -mrSges(7,1) * t190 - mrSges(7,2) * t356;
t137 = -mrSges(5,1) * t188 + mrSges(5,2) * t189;
t131 = mrSges(4,1) * t182 + mrSges(4,2) * t183;
t123 = Ifges(4,1) * t183 - Ifges(4,4) * t182 + Ifges(4,5) * t365;
t122 = Ifges(4,4) * t183 - Ifges(4,2) * t182 + Ifges(4,6) * t365;
t105 = Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * t233;
t104 = Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * t233;
t103 = Ifges(5,5) * t189 + Ifges(5,6) * t188 + Ifges(5,3) * t233;
t100 = pkin(5) * t410 - t110;
t99 = -Ifges(7,1) * t356 + Ifges(7,4) * t190 + Ifges(7,5) * t210;
t98 = -Ifges(7,4) * t356 + Ifges(7,2) * t190 + Ifges(7,6) * t210;
t96 = mrSges(6,1) * t233 - mrSges(6,3) * t136;
t95 = -mrSges(6,2) * t233 - mrSges(6,3) * t135;
t84 = mrSges(5,1) * t182 - mrSges(5,3) * t125;
t83 = -mrSges(5,2) * t182 + mrSges(5,3) * t124;
t82 = (t245 * t341 + t264 * t339) * t349 + t366;
t80 = -mrSges(7,2) * t173 + mrSges(7,3) * t118;
t79 = mrSges(7,1) * t173 - mrSges(7,3) * t117;
t75 = mrSges(6,1) * t135 + mrSges(6,2) * t136;
t72 = -mrSges(5,1) * t124 + mrSges(5,2) * t125;
t69 = Ifges(6,5) * t136 - Ifges(6,6) * t135 + Ifges(6,3) * t233;
t66 = -mrSges(7,1) * t118 + mrSges(7,2) * t117;
t61 = mrSges(7,1) * t135 - mrSges(7,3) * t93;
t60 = -mrSges(7,2) * t135 + mrSges(7,3) * t92;
t54 = -pkin(5) * t380 - t56;
t51 = -mrSges(7,1) * t92 + mrSges(7,2) * t93;
t50 = Ifges(5,1) * t125 + Ifges(5,4) * t124 + Ifges(5,5) * t182;
t49 = Ifges(5,4) * t125 + Ifges(5,2) * t124 + Ifges(5,6) * t182;
t47 = mrSges(6,1) * t182 - mrSges(6,3) * t68;
t46 = -mrSges(6,2) * t182 - mrSges(6,3) * t67;
t44 = Ifges(7,1) * t117 + Ifges(7,4) * t118 + Ifges(7,5) * t173;
t43 = Ifges(7,4) * t117 + Ifges(7,2) * t118 + Ifges(7,6) * t173;
t38 = Ifges(7,1) * t93 + Ifges(7,4) * t92 + Ifges(7,5) * t135;
t37 = Ifges(7,4) * t93 + Ifges(7,2) * t92 + Ifges(7,6) * t135;
t22 = -t233 * pkin(5) - t24;
t21 = -mrSges(7,2) * t67 + mrSges(7,3) * t35;
t20 = mrSges(7,1) * t67 - mrSges(7,3) * t34;
t12 = -mrSges(7,1) * t35 + mrSges(7,2) * t34;
t9 = Ifges(7,1) * t34 + Ifges(7,4) * t35 + Ifges(7,5) * t67;
t8 = Ifges(7,4) * t34 + Ifges(7,2) * t35 + Ifges(7,6) * t67;
t3 = -t182 * pkin(5) - t5;
t13 = [(0.2e1 * (t280 * t350 + t281 * t346) * mrSges(3,3) + ((t293 * t477 + Ifges(3,5) * t342 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t350) * t490) * t350 + (t339 * (Ifges(4,5) * t234 - Ifges(4,6) * t233 + Ifges(4,3) * t287) + t295 * t477 - 0.2e1 * Ifges(3,6) * t342 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t346 + (Ifges(3,1) - Ifges(3,2)) * t350) * t490) * t346) * qJD(2)) * t340 + (t324 - 0.2e1 * t422 - 0.2e1 * t423) * t342 + 0.2e1 * m(3) * (t280 * t295 - t281 * t293) + (-t122 + 0.2e1 * t426 + t488) * t233 + (t123 + 0.2e1 * t425) * t234 + t287 * t121 + 0.2e1 * t81 * t192 + 0.2e1 * t82 * t193 + t188 * t49 + t189 * t50 + t183 * t154 + 0.2e1 * t177 * t131 + 0.2e1 * t140 * t160 + 0.2e1 * t139 * t161 + (-t70 + t36) * t67 + 0.2e1 * t26 * t151 + 0.2e1 * t27 * t152 + 0.2e1 * t78 * t137 + t136 * t30 + 0.2e1 * t133 * t72 + t124 * t104 + t125 * t105 + 0.2e1 * t6 * t95 + 0.2e1 * t5 * t96 + t92 * t8 + t93 * t9 + 0.2e1 * t89 * t31 + 0.2e1 * t74 * t83 + 0.2e1 * t73 * t84 + 0.2e1 * t40 * t75 + t68 * t71 + 0.2e1 * t1 * t60 + 0.2e1 * t2 * t61 + 0.2e1 * t3 * t51 + 0.2e1 * t25 * t46 + 0.2e1 * t24 * t47 + t35 * t37 + t34 * t38 + 0.2e1 * t22 * t12 + 0.2e1 * t10 * t20 + 0.2e1 * t11 * t21 + (t103 + t69 - t153) * t182 + (t7 - t29) * t135 + (t1 * t11 + t10 * t2 + t22 * t3) * t478 + (t24 * t5 + t25 * t6 + t40 * t89) * t479 + (t133 * t78 + t26 * t74 + t27 * t73) * t480 + (t139 * t82 + t140 * t81 + t177 * t185) * t481; (-t275 / 0.2e1 + t163 / 0.2e1 + t106 / 0.2e1) * t233 + t9 * t463 + t8 * t464 + t68 * t465 + t136 * t466 + t44 * t467 + t43 * t468 + t174 * t469 + t211 * t470 + t121 * t447 + (-t267 / 0.2e1 + t202 / 0.2e1 + t148 / 0.2e1) * t182 + ((t425 + t123 / 0.2e1) * t345 + (-t426 + t122 / 0.2e1 - t28 / 0.2e1 - t48 / 0.2e1) * t349 + (Ifges(4,3) * t447 + (Ifges(4,5) * t345 + Ifges(4,6) * t349) * t339 / 0.2e1) * t382 + (-m(4) * t185 - t131) * pkin(2) + ((-t139 * mrSges(4,3) + t154 / 0.2e1) * t349 + (-t140 * mrSges(4,3) + t69 / 0.2e1 + t103 / 0.2e1 - t153 / 0.2e1) * t345) * qJD(3)) * t339 + t391 * t210 + t389 * t173 + t384 * t135 + t385 * t67 - Ifges(3,6) * t382 + (t137 - t193) * t279 + m(4) * (-t139 * t279 + t140 * t278 + t292 * t82 + t294 * t81) + t82 * t299 + t81 * t300 + t290 * t49 / 0.2e1 + t291 * t50 / 0.2e1 + t292 * t161 + t294 * t160 + t287 * t274 / 0.2e1 + t234 * t276 / 0.2e1 + t278 * t192 + t183 * t268 / 0.2e1 + t269 * t72 + t177 * t273 + t26 * t253 + t27 * t254 + t246 * t105 / 0.2e1 + t247 * t104 / 0.2e1 + t78 * t221 + t215 * t31 + t73 * t208 + t74 * t209 + t124 * t203 / 0.2e1 + t125 * t204 / 0.2e1 + t199 * t75 + t200 * t84 + t201 * t83 + t6 * t197 + t5 * t198 + t188 * t164 / 0.2e1 + t189 * t165 / 0.2e1 + t133 * t176 + t40 * t155 + t25 * t157 + t24 * t158 + t146 * t151 + t147 * t152 + t2 * t145 + t1 * t144 + t3 * t138 + t117 * t38 / 0.2e1 + t118 * t37 / 0.2e1 + t89 * t115 + t110 * t47 + t111 * t46 + t56 * t96 + t35 * t98 / 0.2e1 + t34 * t99 / 0.2e1 + t100 * t12 + t57 * t95 + t10 * t79 + t11 * t80 + t22 * t66 + t18 * t60 + t19 * t61 + t58 * t20 + t59 * t21 + t54 * t51 + t324 + m(5) * (t133 * t279 + t146 * t74 + t147 * t73 + t200 * t27 + t201 * t26 + t269 * t78) + m(6) * (t110 * t5 + t111 * t6 + t199 * t89 + t215 * t40 + t24 * t56 + t25 * t57) + m(7) * (t1 * t59 + t10 * t19 + t100 * t3 + t11 * t18 + t2 * t58 + t22 * t54) - t422 - t423; -t356 * t44 + 0.2e1 * (-t299 + t221) * t279 + (t42 - t107) * t210 + (-0.2e1 * pkin(2) * t273 + t345 * t276 + (t275 + t487) * t349 + ((t292 * t476 + t268) * t349 + (t294 * t476 + t148 + t202 - t267) * t345) * qJD(3)) * t339 + t341 * t274 + 0.2e1 * t278 * t300 + t290 * t164 + t291 * t165 + 0.2e1 * t269 * t176 + 0.2e1 * t146 * t253 + 0.2e1 * t147 * t254 + t246 * t204 + t247 * t203 + t211 * t108 + 0.2e1 * t215 * t115 + 0.2e1 * t200 * t208 + 0.2e1 * t201 * t209 + 0.2e1 * t199 * t155 + 0.2e1 * t57 * t197 + 0.2e1 * t56 * t198 + t190 * t43 + t174 * t150 + 0.2e1 * t111 * t157 + 0.2e1 * t110 * t158 + 0.2e1 * t19 * t145 + (-t149 + t97) * t173 + 0.2e1 * t18 * t144 + 0.2e1 * t54 * t138 + t117 * t99 + t118 * t98 + 0.2e1 * t100 * t66 + 0.2e1 * t58 * t79 + 0.2e1 * t59 * t80 + (t100 * t54 + t18 * t59 + t19 * t58) * t478 + (t110 * t56 + t111 * t57 + t199 * t215) * t479 + (t146 * t201 + t147 * t200 + t269 * t279) * t480 + (t278 * t294 - t279 * t292) * t481; -(-t6 * mrSges(6,3) + t391) * t352 + t189 * t454 + t188 * t456 + t68 * t459 + t136 * t460 + t34 * t461 + t35 * t462 + t143 * t467 + t142 * t468 + t49 * t441 + t50 * t444 + t125 * t448 + t124 * t450 + (-t24 * mrSges(6,3) + t37 * t446 + t38 * t442 + t469) * t289 + (t470 + t9 * t442 + t8 * t446 - t5 * mrSges(6,3) + (t37 * t443 + t38 * t446) * qJD(6)) * t297 + t482 * t78 + (t105 * t441 + (pkin(4) * t75 - t104 / 0.2e1) * t344) * qJD(4) + m(7) * (t1 * t167 + t10 * t91 + t11 * t90 + t166 * t2 + t205 * t22 + t250 * t3) + (-t152 * t395 - t151 * t396 + m(5) * (-t395 * t73 - t396 * t74 + t360) + t348 * t83 - t344 * t84) * pkin(11) + (-t25 * mrSges(6,3) + t389) * t288 + m(6) * (-t205 * t24 + t206 * t25 - t250 * t5 + t251 * t6 + t335 * t40 + t390 * t89) + t373 * t67 + t374 * t135 + t371 * t182 + t372 * t233 + t121 + (t12 - t47) * t250 + (t51 - t96) * t205 + t335 * t31 + t133 * t302 + t251 * t46 + t40 * t235 + t1 * t227 + t2 * t228 + t3 * t220 + t89 * t216 + t206 * t95 + ((-t344 * t74 - t348 * t73) * qJD(4) + t360) * mrSges(5,3) + t10 * t186 + t11 * t187 + t167 * t21 + t166 * t20 + t22 * t156 + t90 * t60 + t91 * t61 - t81 * mrSges(4,2) + t82 * mrSges(4,1) - pkin(3) * t72; -(-t57 * mrSges(6,3) + t384) * t352 + t291 * t454 + t290 * t456 + t174 * t459 + t211 * t460 + t117 * t461 + t118 * t462 + t143 * t463 + t142 * t464 + t164 * t441 + t165 * t444 + t246 * t448 + t247 * t450 + (-mrSges(4,1) + t482) * t279 + (-t110 * mrSges(6,3) + t442 * t99 + t446 * t98 + t465) * t289 + (t466 + t44 * t442 + t43 * t446 - t56 * mrSges(6,3) + (t443 * t98 + t446 * t99) * qJD(6)) * t297 + (t204 * t441 + (pkin(4) * t155 - t203 / 0.2e1) * t344) * qJD(4) + (-t254 * t395 - t253 * t396 + m(5) * (-t200 * t395 - t201 * t396 + t359) + t348 * t209 - t344 * t208) * pkin(11) + (-t372 * t349 + (-Ifges(4,6) + t371) * t398) * t339 + m(6) * (-t110 * t205 + t111 * t206 + t199 * t335 + t215 * t390 - t250 * t56 + t251 * t57) + (-t111 * mrSges(6,3) + t385) * t288 + t373 * t173 + t374 * t210 + t335 * t115 + t269 * t302 - t278 * mrSges(4,2) + t251 * t157 + t199 * t235 + t18 * t227 + t19 * t228 + t54 * t220 + t215 * t216 + t206 * t197 + t58 * t186 + t59 * t187 + ((-t200 * t348 - t201 * t344) * qJD(4) + t359) * mrSges(5,3) - pkin(3) * t176 + t167 * t80 + t166 * t79 + t100 * t156 + t91 * t145 + t90 * t144 + t323 + (-t158 + t66) * t250 + (t138 - t198) * t205 + m(7) * (t100 * t205 + t166 * t19 + t167 * t18 + t250 * t54 + t58 * t91 + t59 * t90); -0.2e1 * pkin(3) * t302 + t156 * t473 + 0.2e1 * t166 * t186 + 0.2e1 * t167 * t187 + t220 * t474 + 0.2e1 * t335 * t216 + 0.2e1 * t90 * t227 + 0.2e1 * t91 * t228 + t348 * t306 + t344 * t308 + (t348 * t319 + (0.2e1 * pkin(4) * t235 - t317) * t344) * qJD(4) + (t206 * t251 + t335 * t390 + t418) * t479 + (t166 * t91 + t167 * t90 + t418) * t478 - (t206 * t475 + t141 - t218) * t352 + (t251 * t475 + t194 - t237) * t288 + (mrSges(6,3) * t473 - t195 * t343 + t196 * t347 + t238) * t289 + (mrSges(6,3) * t474 - t343 * t142 + t347 * t143 + t219 + (-t195 * t347 - t196 * t343) * qJD(6)) * t297; t488 + t67 * t453 + t93 * t455 + t92 * t457 + t135 * t458 + t8 * t442 + t9 * t445 + t34 * t449 + t35 * t451 + t46 * t439 + t21 * t413 + t47 * t383 + m(7) * (t3 * t334 + ((-t10 * t347 - t11 * t343) * qJD(6) + t483) * t333) + t483 * mrSges(7,3) - t20 * t414 - t10 * t387 - t11 * t388 - t61 * t377 - t60 * t378 + t334 * t12 + t3 * t311 + t22 * t301 - t26 * mrSges(5,2) + t27 * mrSges(5,1) - t6 * mrSges(6,2) + t5 * mrSges(6,1) + (t338 * t6 + t419 * t5) * t471 + t38 * t369 + t37 * t370; -t356 * t455 - t487 + t173 * t453 + t190 * t457 + t210 * t458 + t43 * t442 + t44 * t445 + t117 * t449 + t118 * t451 + t157 * t439 + t80 * t413 + t158 * t383 + m(7) * (t334 * t54 + ((-t343 * t59 - t347 * t58) * qJD(6) + t484) * t333) + t484 * mrSges(7,3) - t79 * t414 - t58 * t387 - t59 * t388 - t145 * t377 - t144 * t378 + t334 * t66 + t54 * t311 + t100 * t301 - t146 * mrSges(5,2) + t147 * mrSges(5,1) - t57 * mrSges(6,2) + t56 * mrSges(6,1) + (t338 * t57 + t419 * t56) * t471 + t99 * t369 + t98 * t370; -t352 * t458 + m(7) * ((-t166 * t347 - t167 * t343) * qJD(6) + t485) * t333 + t304 + (t297 * t318 + t195) * t370 + (-t288 * t439 - t289 * t383) * mrSges(6,3) + t288 * t453 + t415 * t455 + t142 * t442 + t143 * t445 + t402 * t449 + t187 * t413 + t485 * mrSges(7,3) - t354 * t316 / 0.2e1 - t186 * t414 - t305 * t416 / 0.2e1 - t166 * t387 - t167 * t388 - t228 * t377 - t227 * t378 + t334 * t156 + t250 * t301 + (-mrSges(5,1) * t395 + mrSges(5,2) * t396) * pkin(11) + (t338 * t471 - mrSges(6,2)) * t206 + (m(7) * t334 - t419 * t471 - mrSges(6,1) + t311) * t205 + t217 + t196 * t369; 0.2e1 * t301 * t334 + t305 * t347 + t307 * t343 + (-t316 * t343 + t318 * t347) * qJD(6); t347 * t20 + t343 * t21 + (-t343 * t61 + t347 * t60) * qJD(6) + m(7) * (t1 * t343 + t2 * t347 + (-t10 * t343 + t11 * t347) * qJD(6)) + m(6) * t40 + t31; t343 * t80 + t347 * t79 + (t144 * t347 - t145 * t343) * qJD(6) + m(7) * (t18 * t343 + t19 * t347 + (-t343 * t58 + t347 * t59) * qJD(6)) + m(6) * t199 + t115; m(7) * (t343 * t90 + t347 * t91 + (-t166 * t343 + t167 * t347) * qJD(6)) + t227 * t393 + t343 * t187 - t228 * t394 + t347 * t186 + m(6) * t390 + t216; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t7; mrSges(7,1) * t19 - mrSges(7,2) * t18 + t42; mrSges(7,1) * t91 - mrSges(7,2) * t90 + t141; t336 + (t311 * t333 - t428) * qJD(6); -t301; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t13(1) t13(2) t13(4) t13(7) t13(11) t13(16); t13(2) t13(3) t13(5) t13(8) t13(12) t13(17); t13(4) t13(5) t13(6) t13(9) t13(13) t13(18); t13(7) t13(8) t13(9) t13(10) t13(14) t13(19); t13(11) t13(12) t13(13) t13(14) t13(15) t13(20); t13(16) t13(17) t13(18) t13(19) t13(20) t13(21);];
Mq  = res;
