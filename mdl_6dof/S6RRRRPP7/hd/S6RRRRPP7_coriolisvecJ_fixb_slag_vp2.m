% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:17:10
% EndTime: 2019-03-09 21:18:08
% DurationCPUTime: 29.45s
% Computational Cost: add. (19167->882), mult. (50568->1192), div. (0->0), fcn. (39181->10), ass. (0->369)
t346 = sin(qJ(2));
t342 = sin(pkin(6));
t404 = qJD(1) * t342;
t390 = t346 * t404;
t343 = cos(pkin(6));
t349 = cos(qJ(2));
t457 = pkin(1) * t349;
t395 = t343 * t457;
t291 = -pkin(8) * t390 + qJD(1) * t395;
t364 = (pkin(2) * t346 - pkin(9) * t349) * t342;
t292 = qJD(1) * t364;
t345 = sin(qJ(3));
t348 = cos(qJ(3));
t211 = t348 * t291 + t345 * t292;
t197 = pkin(10) * t390 + t211;
t331 = t343 * t346 * pkin(1);
t374 = pkin(3) * t345 - pkin(10) * t348;
t418 = t342 * t349;
t219 = (t331 + (pkin(8) + t374) * t418) * qJD(1);
t344 = sin(qJ(4));
t347 = cos(qJ(4));
t130 = -t344 * t197 + t347 * t219;
t414 = t348 * t349;
t259 = (t344 * t346 + t347 * t414) * t404;
t319 = -pkin(3) * t348 - pkin(10) * t345 - pkin(2);
t415 = t347 * t348;
t333 = pkin(9) * t415;
t389 = t349 * t404;
t377 = t345 * t389;
t396 = qJD(5) * t347;
t312 = t374 * qJD(3);
t401 = qJD(3) * t345;
t452 = pkin(9) * t344;
t406 = t347 * t312 + t401 * t452;
t453 = pkin(4) * t345;
t539 = -pkin(4) * t377 + t259 * qJ(5) - t130 - t345 * t396 + (-qJ(5) * t415 + t453) * qJD(3) + (-t333 + (qJ(5) * t345 - t319) * t344) * qJD(4) + t406;
t131 = t347 * t197 + t344 * t219;
t258 = (-t344 * t414 + t346 * t347) * t404;
t397 = qJD(4) * t347;
t407 = t344 * t312 + t319 * t397;
t416 = t345 * t347;
t538 = qJ(5) * t258 + t131 - (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t416 - (-qJD(5) * t345 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t348) * t344 - t407;
t327 = qJD(1) * t343 + qJD(2);
t402 = qJD(2) * t349;
t387 = t345 * t402;
t400 = qJD(3) * t348;
t228 = t327 * t401 + (t346 * t400 + t387) * t404;
t470 = t228 / 0.2e1;
t271 = t327 * t345 + t348 * t390;
t320 = qJD(3) - t389;
t220 = -t271 * t344 + t320 * t347;
t386 = t348 * t402;
t227 = t327 * t400 + (-t346 * t401 + t386) * t404;
t403 = qJD(2) * t342;
t381 = qJD(1) * t403;
t376 = t346 * t381;
t138 = qJD(4) * t220 + t227 * t347 + t344 * t376;
t221 = t271 * t347 + t320 * t344;
t139 = -qJD(4) * t221 - t227 * t344 + t347 * t376;
t341 = sin(pkin(11));
t423 = cos(pkin(11));
t79 = t138 * t423 + t341 * t139;
t487 = t79 / 0.2e1;
t78 = t138 * t341 - t139 * t423;
t488 = t78 / 0.2e1;
t489 = -t78 / 0.2e1;
t508 = Ifges(7,4) + Ifges(6,5);
t509 = Ifges(6,1) + Ifges(7,1);
t529 = Ifges(6,4) * t489 + Ifges(7,5) * t488 + t508 * t470 + t509 * t487;
t507 = Ifges(7,5) - Ifges(6,4);
t189 = -t258 * t423 + t259 * t341;
t379 = t423 * t344;
t308 = t341 * t347 + t379;
t421 = t341 * t344;
t360 = t423 * t347 - t421;
t398 = qJD(4) * t345;
t214 = -t308 * t400 - t360 * t398;
t410 = t189 + t214;
t190 = t341 * t258 + t259 * t423;
t215 = -t308 * t398 + t360 * t400;
t409 = t190 - t215;
t537 = t377 - t401;
t525 = t538 * t341 + t539 * t423;
t524 = t539 * t341 - t538 * t423;
t210 = -t345 * t291 + t292 * t348;
t196 = -pkin(3) * t390 - t210;
t454 = pkin(4) * t344;
t516 = pkin(4) * t258 - t196 + t397 * t453 + (pkin(9) + t454) * t400;
t270 = t327 * t348 - t345 * t390;
t187 = t308 * t270;
t298 = t308 * qJD(4);
t412 = t187 - t298;
t188 = t360 * t270;
t299 = t360 * qJD(4);
t411 = t188 - t299;
t265 = qJD(4) - t270;
t440 = t221 * Ifges(5,4);
t125 = t220 * Ifges(5,2) + t265 * Ifges(5,6) + t440;
t218 = Ifges(5,4) * t220;
t126 = t221 * Ifges(5,1) + t265 * Ifges(5,5) + t218;
t405 = pkin(8) * t418 + t331;
t294 = t405 * qJD(1);
t253 = t327 * pkin(9) + t294;
t286 = (-pkin(2) * t349 - pkin(9) * t346 - pkin(1)) * t342;
t263 = qJD(1) * t286;
t183 = -t345 * t253 + t263 * t348;
t168 = -pkin(3) * t320 - t183;
t252 = -t327 * pkin(2) - t291;
t164 = -t270 * pkin(3) - t271 * pkin(10) + t252;
t184 = t348 * t253 + t345 * t263;
t169 = pkin(10) * t320 + t184;
t98 = t347 * t164 - t169 * t344;
t99 = t164 * t344 + t169 * t347;
t365 = t344 * t99 + t347 * t98;
t367 = Ifges(5,5) * t347 - Ifges(5,6) * t344;
t445 = Ifges(5,4) * t347;
t369 = -Ifges(5,2) * t344 + t445;
t446 = Ifges(5,4) * t344;
t371 = Ifges(5,1) * t347 - t446;
t372 = mrSges(5,1) * t344 + mrSges(5,2) * t347;
t458 = t347 / 0.2e1;
t459 = -t344 / 0.2e1;
t466 = t265 / 0.2e1;
t471 = t221 / 0.2e1;
t473 = t220 / 0.2e1;
t535 = t365 * mrSges(5,3) - t125 * t459 - t126 * t458 - t168 * t372 - t367 * t466 - t369 * t473 - t371 * t471;
t435 = t270 * Ifges(4,2);
t439 = t221 * Ifges(5,5);
t441 = t220 * Ifges(5,6);
t124 = t265 * Ifges(5,3) + t439 + t441;
t145 = -t423 * t220 + t221 * t341;
t82 = qJ(5) * t220 + t99;
t424 = t341 * t82;
t81 = -qJ(5) * t221 + t98;
t60 = pkin(4) * t265 + t81;
t26 = t423 * t60 - t424;
t17 = -t265 * pkin(5) + qJD(6) - t26;
t76 = t423 * t82;
t27 = t341 * t60 + t76;
t18 = qJ(6) * t265 + t27;
t361 = t341 * t220 + t221 * t423;
t375 = -Ifges(5,3) / 0.2e1 - Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t393 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t394 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t533 = -t320 * Ifges(4,6) / 0.2e1 - t271 * Ifges(4,4) / 0.2e1;
t531 = -t435 / 0.2e1 + t533;
t64 = Ifges(6,5) * t361 - t145 * Ifges(6,6) + t265 * Ifges(6,3);
t65 = Ifges(7,4) * t361 + t265 * Ifges(7,2) + t145 * Ifges(7,6);
t495 = t393 * t145 - t375 * t265 + t394 * t361 - t17 * mrSges(7,1) - t184 * mrSges(4,3) - t27 * mrSges(6,2) - t99 * mrSges(5,2) + t124 / 0.2e1 + t531 + t64 / 0.2e1 + t65 / 0.2e1 + t18 * mrSges(7,3) + t441 / 0.2e1 + t439 / 0.2e1 + t252 * mrSges(4,1) + t26 * mrSges(6,1) + t98 * mrSges(5,1) + t533;
t534 = t435 / 0.2e1 - t495;
t293 = qJD(2) * t364;
t280 = qJD(1) * t293;
t419 = t342 * t346;
t328 = pkin(8) * t419;
t304 = -t328 + t395;
t295 = t304 * qJD(2);
t281 = qJD(1) * t295;
t122 = -t253 * t401 + t263 * t400 + t345 * t280 + t348 * t281;
t111 = pkin(10) * t376 + t122;
t296 = t405 * qJD(2);
t282 = qJD(1) * t296;
t135 = t228 * pkin(3) - t227 * pkin(10) + t282;
t34 = -qJD(4) * t99 - t111 * t344 + t347 * t135;
t12 = pkin(4) * t228 - qJ(5) * t138 - qJD(5) * t221 + t34;
t399 = qJD(4) * t344;
t33 = t347 * t111 + t344 * t135 + t164 * t397 - t169 * t399;
t14 = qJ(5) * t139 + qJD(5) * t220 + t33;
t4 = t341 * t12 + t423 * t14;
t1 = qJ(6) * t228 + qJD(6) * t265 + t4;
t511 = -t228 / 0.2e1;
t123 = -t253 * t400 - t263 * t401 + t280 * t348 - t345 * t281;
t112 = -pkin(3) * t376 - t123;
t61 = -pkin(4) * t139 + t112;
t9 = pkin(5) * t78 - qJ(6) * t79 - qJD(6) * t361 + t61;
t493 = mrSges(6,1) * t61 + mrSges(7,1) * t9 - mrSges(7,2) * t1 - mrSges(6,3) * t4 - t79 * Ifges(6,4) / 0.2e1 + Ifges(6,6) * t511 + Ifges(7,6) * t470 + 0.2e1 * Ifges(7,3) * t488 + (-t489 + t488) * Ifges(6,2) + (t507 + Ifges(7,5)) * t487;
t3 = t12 * t423 - t341 * t14;
t2 = -t228 * pkin(5) - t3;
t532 = t61 * mrSges(6,2) + mrSges(7,2) * t2 - mrSges(6,3) * t3 - t9 * mrSges(7,3) + 0.2e1 * t529;
t504 = t507 * t145 + t508 * t265 + t361 * t509;
t530 = t504 / 0.2e1;
t528 = Ifges(6,6) - Ifges(7,6);
t527 = -qJ(6) * t537 - qJD(6) * t348 + t524;
t526 = pkin(5) * t537 - t525;
t63 = Ifges(7,5) * t361 + Ifges(7,6) * t265 + Ifges(7,3) * t145;
t66 = Ifges(6,4) * t361 - Ifges(6,2) * t145 + Ifges(6,6) * t265;
t523 = t66 - t63;
t522 = Ifges(4,5) * t227;
t520 = t281 * mrSges(3,2);
t450 = -qJ(5) - pkin(10);
t380 = qJD(4) * t450;
t297 = t344 * t380 + t396;
t357 = -qJD(5) * t344 + t347 * t380;
t216 = t297 * t341 - t357 * t423;
t205 = pkin(3) * t271 - pkin(10) * t270;
t120 = t347 * t183 + t344 * t205;
t422 = t270 * t344;
t103 = -qJ(5) * t422 + t120;
t119 = -t183 * t344 + t347 * t205;
t92 = -qJ(5) * t270 * t347 + pkin(4) * t271 + t119;
t41 = -t341 * t103 + t423 * t92;
t519 = -t41 - t216;
t217 = t297 * t423 + t341 * t357;
t42 = t423 * t103 + t341 * t92;
t518 = -t42 + t217;
t288 = t360 * t345;
t517 = -pkin(5) * t410 + qJ(6) * t409 - qJD(6) * t288 + t516;
t264 = Ifges(4,4) * t270;
t430 = t320 * Ifges(4,5);
t433 = t271 * Ifges(4,1);
t180 = t264 + t430 + t433;
t356 = t183 * mrSges(4,3) - t180 / 0.2e1 - t252 * mrSges(4,2) - t430 / 0.2e1;
t515 = t356 + t535;
t514 = -t123 * mrSges(4,1) + t122 * mrSges(4,2);
t510 = -Ifges(3,6) * t327 / 0.2e1;
t503 = -qJ(6) * t271 + t518;
t502 = t271 * pkin(5) - t519;
t149 = pkin(4) * t422 + t184;
t501 = pkin(4) * t399 - pkin(5) * t412 + qJ(6) * t411 - qJD(6) * t308 - t149;
t117 = mrSges(6,1) * t265 - mrSges(6,3) * t361;
t118 = -mrSges(7,1) * t265 + mrSges(7,2) * t361;
t500 = t117 - t118;
t284 = t328 + (-pkin(2) - t457) * t343;
t300 = -t343 * t348 + t345 * t419;
t301 = t343 * t345 + t348 * t419;
t192 = t300 * pkin(3) - t301 * pkin(10) + t284;
t285 = pkin(9) * t343 + t405;
t207 = t348 * t285 + t345 * t286;
t194 = -pkin(10) * t418 + t207;
t114 = t344 * t192 + t347 * t194;
t32 = t423 * t81 - t424;
t499 = -t32 + qJD(6);
t279 = t344 * t319 + t333;
t498 = t33 * t347 - t34 * t344;
t497 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t20 = Ifges(6,5) * t79 - Ifges(6,6) * t78 + Ifges(6,3) * t228;
t21 = Ifges(7,4) * t79 + Ifges(7,2) * t228 + Ifges(7,6) * t78;
t136 = Ifges(5,6) * t139;
t137 = Ifges(5,5) * t138;
t54 = Ifges(5,3) * t228 + t136 + t137;
t496 = t54 + t21 + t20;
t121 = -pkin(4) * t220 + qJD(5) + t168;
t45 = pkin(5) * t145 - qJ(6) * t361 + t121;
t494 = t121 * mrSges(6,1) + t45 * mrSges(7,1) + t63 / 0.2e1 - t66 / 0.2e1 - t18 * mrSges(7,2) - t27 * mrSges(6,3);
t56 = t138 * Ifges(5,1) + t139 * Ifges(5,4) + t228 * Ifges(5,5);
t490 = t56 / 0.2e1;
t485 = pkin(1) * mrSges(3,1);
t484 = pkin(1) * mrSges(3,2);
t483 = -t125 / 0.2e1;
t482 = t138 / 0.2e1;
t481 = t139 / 0.2e1;
t480 = -t145 / 0.2e1;
t479 = t145 / 0.2e1;
t477 = -t361 / 0.2e1;
t476 = t361 / 0.2e1;
t474 = -t220 / 0.2e1;
t472 = -t221 / 0.2e1;
t468 = -t264 / 0.2e1;
t467 = -t265 / 0.2e1;
t464 = -t300 / 0.2e1;
t462 = t301 / 0.2e1;
t460 = t343 / 0.2e1;
t456 = pkin(4) * t221;
t455 = pkin(4) * t341;
t451 = pkin(9) * t348;
t240 = -qJD(3) * t300 + t342 * t386;
t241 = -t301 * t344 - t347 * t418;
t388 = t346 * t403;
t163 = qJD(4) * t241 + t240 * t347 + t344 * t388;
t239 = qJD(3) * t301 + t342 * t387;
t363 = -t301 * t347 + t344 * t418;
t140 = -t285 * t401 + t286 * t400 + t345 * t293 + t348 * t295;
t128 = pkin(10) * t388 + t140;
t154 = t239 * pkin(3) - t240 * pkin(10) + t296;
t44 = -qJD(4) * t114 - t128 * t344 + t347 * t154;
t16 = pkin(4) * t239 - qJ(5) * t163 + qJD(5) * t363 + t44;
t162 = qJD(4) * t363 - t240 * t344 + t347 * t388;
t43 = t347 * t128 + t344 * t154 + t192 * t397 - t194 * t399;
t28 = qJ(5) * t162 + qJD(5) * t241 + t43;
t8 = t341 * t16 + t423 * t28;
t100 = qJ(5) * t241 + t114;
t113 = t347 * t192 - t194 * t344;
t90 = pkin(4) * t300 + qJ(5) * t363 + t113;
t40 = t423 * t100 + t341 * t90;
t449 = mrSges(5,3) * t220;
t448 = mrSges(5,3) * t221;
t447 = Ifges(3,4) * t346;
t438 = t227 * Ifges(4,1);
t437 = t227 * Ifges(4,4);
t436 = t228 * Ifges(4,4);
t417 = t344 * t345;
t152 = -mrSges(5,1) * t220 + mrSges(5,2) * t221;
t232 = mrSges(4,1) * t320 - mrSges(4,3) * t271;
t413 = -t152 + t232;
t408 = -mrSges(3,1) * t327 - mrSges(4,1) * t270 + mrSges(4,2) * t271 + mrSges(3,3) * t390;
t310 = t347 * t319;
t234 = -qJ(5) * t416 + t310 + (-pkin(4) - t452) * t348;
t245 = -qJ(5) * t417 + t279;
t171 = t341 * t234 + t423 * t245;
t313 = pkin(4) * t417 + t345 * pkin(9);
t392 = -Ifges(4,6) * t228 + Ifges(4,3) * t376 + t522;
t338 = -pkin(4) * t347 - pkin(3);
t391 = t423 * pkin(4);
t30 = t78 * mrSges(6,1) + t79 * mrSges(6,2);
t29 = t78 * mrSges(7,1) - t79 * mrSges(7,3);
t53 = -t228 * mrSges(7,1) + t79 * mrSges(7,2);
t206 = -t345 * t285 + t286 * t348;
t193 = pkin(3) * t418 - t206;
t373 = mrSges(5,1) * t347 - mrSges(5,2) * t344;
t370 = Ifges(5,1) * t344 + t445;
t368 = Ifges(5,2) * t347 + t446;
t366 = Ifges(5,5) * t344 + Ifges(5,6) * t347;
t141 = -t285 * t400 - t286 * t401 + t293 * t348 - t345 * t295;
t7 = t16 * t423 - t341 * t28;
t39 = -t341 * t100 + t423 * t90;
t170 = t234 * t423 - t341 * t245;
t148 = -pkin(4) * t241 + t193;
t322 = Ifges(3,4) * t389;
t358 = -t291 * mrSges(3,3) + t327 * Ifges(3,5) + Ifges(3,1) * t390 / 0.2e1 + t322 / 0.2e1;
t129 = -pkin(3) * t388 - t141;
t84 = -pkin(4) * t162 + t129;
t355 = -t34 * mrSges(5,1) - t3 * mrSges(6,1) + t2 * mrSges(7,1) + t33 * mrSges(5,2) + t4 * mrSges(6,2) - t1 * mrSges(7,3);
t354 = t183 * mrSges(4,1) + t320 * Ifges(4,3) + t271 * Ifges(4,5) + t270 * Ifges(4,6) + t510 - (Ifges(3,2) * t349 + t447) * t404 / 0.2e1 - t184 * mrSges(4,2) - t294 * mrSges(3,3);
t337 = -t391 - pkin(5);
t334 = qJ(6) + t455;
t321 = t450 * t347;
t318 = Ifges(3,5) * t349 * t381;
t290 = -t327 * mrSges(3,2) + mrSges(3,3) * t389;
t287 = t308 * t345;
t278 = -t344 * t451 + t310;
t251 = -t321 * t423 + t421 * t450;
t250 = -t321 * t341 - t379 * t450;
t231 = -mrSges(4,2) * t320 + mrSges(4,3) * t270;
t229 = -pkin(5) * t360 - qJ(6) * t308 + t338;
t204 = -qJD(4) * t279 + t406;
t203 = (-t347 * t401 - t348 * t399) * pkin(9) + t407;
t201 = -mrSges(4,2) * t376 - mrSges(4,3) * t228;
t200 = mrSges(4,1) * t376 - mrSges(4,3) * t227;
t195 = pkin(5) * t287 - qJ(6) * t288 + t313;
t175 = mrSges(5,1) * t265 - t448;
t174 = -mrSges(5,2) * t265 + t449;
t173 = t341 * t241 - t363 * t423;
t172 = -t241 * t423 - t341 * t363;
t166 = t348 * pkin(5) - t170;
t165 = -qJ(6) * t348 + t171;
t159 = mrSges(4,1) * t228 + mrSges(4,2) * t227;
t143 = Ifges(4,5) * t376 - t436 + t438;
t142 = -t228 * Ifges(4,2) + Ifges(4,6) * t376 + t437;
t116 = -mrSges(6,2) * t265 - mrSges(6,3) * t145;
t115 = -mrSges(7,2) * t145 + mrSges(7,3) * t265;
t106 = -mrSges(5,2) * t228 + mrSges(5,3) * t139;
t105 = mrSges(5,1) * t228 - mrSges(5,3) * t138;
t96 = t341 * t162 + t163 * t423;
t94 = -t162 * t423 + t163 * t341;
t86 = mrSges(6,1) * t145 + mrSges(6,2) * t361;
t85 = mrSges(7,1) * t145 - mrSges(7,3) * t361;
t83 = -mrSges(5,1) * t139 + mrSges(5,2) * t138;
t59 = pkin(5) * t361 + qJ(6) * t145 + t456;
t57 = pkin(5) * t172 - qJ(6) * t173 + t148;
t55 = t138 * Ifges(5,4) + t139 * Ifges(5,2) + t228 * Ifges(5,6);
t52 = mrSges(6,1) * t228 - mrSges(6,3) * t79;
t51 = -mrSges(6,2) * t228 - mrSges(6,3) * t78;
t50 = -mrSges(7,2) * t78 + mrSges(7,3) * t228;
t36 = -t300 * pkin(5) - t39;
t35 = qJ(6) * t300 + t40;
t31 = t341 * t81 + t76;
t10 = pkin(5) * t94 - qJ(6) * t96 - qJD(6) * t173 + t84;
t6 = -t239 * pkin(5) - t7;
t5 = qJ(6) * t239 + qJD(6) * t300 + t8;
t11 = [t227 * (Ifges(4,1) * t301 - Ifges(4,4) * t300) / 0.2e1 + (Ifges(4,4) * t301 - Ifges(4,2) * t300) * t511 + (-Ifges(5,4) * t363 + Ifges(5,2) * t241 + Ifges(5,6) * t300) * t481 + (-Ifges(5,1) * t363 + Ifges(5,4) * t241 + Ifges(5,5) * t300) * t482 + t112 * (-mrSges(5,1) * t241 - mrSges(5,2) * t363) + t34 * (mrSges(5,1) * t300 + mrSges(5,3) * t363) - t363 * t490 + (t358 * t349 + (t510 + t354) * t346) * t403 + (t1 * t300 - t173 * t9 + t18 * t239 - t45 * t96) * mrSges(7,3) + (Ifges(7,5) * t96 + Ifges(7,6) * t239) * t479 + (Ifges(7,5) * t173 + Ifges(7,6) * t300) * t488 + (-t122 * t300 - t123 * t301 - t183 * t240 - t184 * t239) * mrSges(4,3) - t343 * t520 + (-Ifges(5,5) * t363 + Ifges(5,6) * t241 - t172 * t528 + t173 * t508 + t300 * t497) * t470 + (Ifges(5,5) * t163 + Ifges(5,6) * t162 + t239 * t497 + t508 * t96 - t528 * t94) * t466 + (Ifges(6,4) * t96 + Ifges(6,6) * t239) * t480 + (Ifges(6,4) * t173 + Ifges(6,6) * t300) * t489 + t271 * (Ifges(4,1) * t240 - Ifges(4,4) * t239) / 0.2e1 + m(3) * (t281 * t405 - t282 * t304 - t291 * t296 + t294 * t295) + ((Ifges(3,5) * t460 - t304 * mrSges(3,3) + (-0.2e1 * t484 + 0.3e1 / 0.2e1 * Ifges(3,4) * t349) * t342) * t349 + (Ifges(4,5) * t462 + Ifges(4,6) * t464 - Ifges(3,6) * t343 - t405 * mrSges(3,3) + (-0.2e1 * t485 - 0.3e1 / 0.2e1 * t447 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - Ifges(4,3) / 0.2e1) * t349) * t342) * t346) * t381 + (-t392 / 0.2e1 - t522 / 0.2e1 + t281 * mrSges(3,3) - Ifges(4,6) * t511 + t514) * t418 + t270 * (Ifges(4,4) * t240 - Ifges(4,2) * t239) / 0.2e1 + m(4) * (t122 * t207 + t123 * t206 + t140 * t184 + t141 * t183 + t252 * t296 + t282 * t284) + m(5) * (t112 * t193 + t113 * t34 + t114 * t33 + t129 * t168 + t43 * t99 + t44 * t98) + m(6) * (t121 * t84 + t148 * t61 + t26 * t7 + t27 * t8 + t3 * t39 + t4 * t40) + m(7) * (t1 * t35 + t10 * t45 + t17 * t6 + t18 * t5 + t2 * t36 + t57 * t9) + (t124 + t65 + t64) * t239 / 0.2e1 + t320 * (Ifges(4,5) * t240 - Ifges(4,6) * t239) / 0.2e1 + (t173 * t509 + t300 * t508) * t487 + (t239 * t508 + t509 * t96) * t476 + (-Ifges(6,2) * t480 + Ifges(7,3) * t479 + t476 * t507 + t494) * t94 + (t121 * t96 + t173 * t61 - t239 * t27 - t300 * t4) * mrSges(6,2) + t496 * t300 / 0.2e1 + t493 * t172 + (-mrSges(3,1) * t343 + mrSges(4,1) * t300 + mrSges(4,2) * t301 + mrSges(3,3) * t419) * t282 + t408 * t296 + t35 * t50 + t40 * t51 + t39 * t52 + t36 * t53 + t57 * t29 + t10 * t85 + t84 * t86 + t173 * t529 + t96 * t530 + t239 * t531 + t318 * t460 + t143 * t462 + t142 * t464 + (Ifges(5,1) * t163 + Ifges(5,4) * t162 + Ifges(5,5) * t239) * t471 + (Ifges(5,4) * t163 + Ifges(5,2) * t162 + Ifges(5,6) * t239) * t473 + t113 * t105 + t114 * t106 + t5 * t115 + t8 * t116 + t7 * t117 + t6 * t118 + t148 * t30 + t129 * t152 + t162 * t125 / 0.2e1 + t163 * t126 / 0.2e1 + t168 * (-mrSges(5,1) * t162 + mrSges(5,2) * t163) + t43 * t174 + t44 * t175 + t193 * t83 + t206 * t200 + t207 * t201 + t140 * t231 + t141 * t232 + t26 * (mrSges(6,1) * t239 - mrSges(6,3) * t96) + t17 * (-mrSges(7,1) * t239 + mrSges(7,2) * t96) + t99 * (-mrSges(5,2) * t239 + mrSges(5,3) * t162) + t98 * (mrSges(5,1) * t239 - mrSges(5,3) * t163) + t240 * t180 / 0.2e1 + t241 * t55 / 0.2e1 + t252 * (mrSges(4,1) * t239 + mrSges(4,2) * t240) + t284 * t159 + t295 * t290 + t3 * (mrSges(6,1) * t300 - mrSges(6,3) * t173) + t2 * (-mrSges(7,1) * t300 + mrSges(7,2) * t173) + t33 * (-mrSges(5,2) * t300 + mrSges(5,3) * t241); ((t433 / 0.2e1 + (-m(4) * t183 + m(5) * t168 - t413) * pkin(9) + t264 / 0.2e1 - t515) * t348 + ((-m(4) * t184 - t231) * pkin(9) - t534) * t345) * qJD(3) + ((qJD(2) * (Ifges(4,5) * t345 + Ifges(4,6) * t348) / 0.2e1 + (t485 + t447 / 0.2e1) * t404 + (t327 / 0.2e1 - qJD(2)) * Ifges(3,6) - t354) * t346 + (-t322 / 0.2e1 + (t484 + (Ifges(3,2) / 0.2e1 - Ifges(3,1) / 0.2e1) * t346) * t404 + (t468 - t433 / 0.2e1 + t356) * t348 + t534 * t345 - t358) * t349) * t404 + t532 * t288 + (-t131 + t203) * t174 + (-t470 * t528 + t493) * t287 + t318 + (t112 * t372 + t371 * t482 + t369 * t481 + t367 * t470 + t438 / 0.2e1 - t436 / 0.2e1 + t282 * mrSges(4,2) + t55 * t459 + t56 * t458 - t123 * mrSges(4,3) + t143 / 0.2e1 + (-t33 * t344 - t34 * t347) * mrSges(5,3) + (-m(4) * t123 + m(5) * t112 - t200 + t83) * pkin(9) + (t126 * t459 + t347 * t483 + t168 * t373 + t368 * t474 + t370 * t472 + (t344 * t98 - t347 * t99) * mrSges(5,3)) * qJD(4)) * t345 + t516 * t86 + (-t17 * t409 + t18 * t410) * mrSges(7,2) + m(5) * (t203 * t99 + t204 * t98 + t278 * t34 + t279 * t33) + (Ifges(6,4) * t190 + Ifges(7,5) * t215 - Ifges(6,2) * t189 - Ifges(7,3) * t214) * t479 + t504 * (-t190 / 0.2e1 + t215 / 0.2e1) + t523 * (t189 / 0.2e1 + t214 / 0.2e1) + (t26 * t409 + t27 * t410) * mrSges(6,3) + (-t130 + t204) * t175 - m(4) * (t183 * t210 + t184 * t211 + t252 * t294) - m(5) * (t130 * t98 + t131 * t99 + t168 * t196) + t524 * t116 + t525 * t117 + (t121 * t516 + t170 * t3 + t171 * t4 + t26 * t525 + t27 * t524 + t313 * t61) * m(6) + t526 * t118 + (t1 * t165 + t166 * t2 + t17 * t526 + t18 * t527 + t195 * t9 + t45 * t517) * m(7) + t527 * t115 + (-t258 * t99 + t259 * t98) * mrSges(5,3) + (Ifges(6,4) * t215 + Ifges(7,5) * t190 + Ifges(6,2) * t214 + Ifges(7,3) * t189) * t480 + m(4) * (-pkin(2) * t282 + t122 * t451) + (-t54 / 0.2e1 + t355 - t394 * t79 + (-Ifges(4,2) / 0.2e1 + t375) * t228 - t20 / 0.2e1 - t393 * t78 - t21 / 0.2e1 - t137 / 0.2e1 - t136 / 0.2e1 + pkin(9) * t201 + t122 * mrSges(4,3) + t437 / 0.2e1 - t282 * mrSges(4,1) + t142 / 0.2e1) * t348 + (-mrSges(7,1) * t410 + mrSges(7,3) * t409) * t45 + (-mrSges(6,1) * t410 - mrSges(6,2) * t409) * t121 - t408 * t294 - t520 + (t214 * t528 + t215 * t508) * t466 + (Ifges(5,5) * t259 + Ifges(5,6) * t258 - t189 * t528 + t190 * t508 + t366 * t398) * t467 + (t189 * t507 + t190 * t509) * t477 + (-t214 * t507 + t215 * t509) * t476 + (Ifges(5,1) * t259 + Ifges(5,4) * t258) * t472 + (Ifges(5,4) * t259 + Ifges(5,2) * t258) * t474 + t258 * t483 + t517 * t85 - pkin(2) * t159 + t165 * t50 + t166 * t53 + t170 * t52 + t171 * t51 + t195 * t29 - t196 * t152 - t211 * t231 - t210 * t232 - t259 * t126 / 0.2e1 - t168 * (-mrSges(5,1) * t258 + mrSges(5,2) * t259) + t278 * t105 + t279 * t106 - t282 * mrSges(3,1) - t291 * t290 + t313 * t30; ((m(6) * t121 + t86) * t454 - t535) * qJD(4) + t532 * t308 - t493 * t360 - t514 + t392 + (t53 - t52) * t250 + (-pkin(3) * t112 - t119 * t98 - t120 * t99 - t168 * t184) * m(5) + (t26 * t411 + t27 * t412) * mrSges(6,3) + (t468 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t271 + t515) * t270 - t495 * t271 - m(6) * (t121 * t149 + t26 * t41 + t27 * t42) + (Ifges(6,4) * t299 + Ifges(7,5) * t188 - Ifges(6,2) * t298 + Ifges(7,3) * t187) * t480 + (Ifges(6,4) * t188 + Ifges(7,5) * t299 - Ifges(6,2) * t187 + Ifges(7,3) * t298) * t479 + (t360 * t528 + t366) * t470 + (-t187 * t528 + t188 * t508) * t467 + (-t298 * t528 + t299 * t508) * t466 + t523 * (t187 / 0.2e1 - t298 / 0.2e1) + (-t17 * t411 + t18 * t412) * mrSges(7,2) + (t50 + t51) * t251 + (t298 * t507 + t299 * t509) * t476 + (t187 * t507 + t188 * t509) * t477 + t504 * (-t188 / 0.2e1 + t299 / 0.2e1) + t503 * t115 + (t1 * t251 + t17 * t502 + t18 * t503 + t2 * t250 + t229 * t9 + t45 * t501) * m(7) + t501 * t85 + t502 * t118 + (-t105 * t344 + t106 * t347 + (-m(5) * t365 - t344 * t174 - t347 * t175) * qJD(4) + m(5) * t498) * pkin(10) + t498 * mrSges(5,3) + m(6) * (-t216 * t26 + t217 * t27 - t250 * t3 + t251 * t4 + t338 * t61) + (-mrSges(7,1) * t412 + mrSges(7,3) * t411) * t45 + (-mrSges(6,1) * t412 - mrSges(6,2) * t411) * t121 + t413 * t184 - pkin(3) * t83 + t55 * t458 + t368 * t481 + t370 * t482 + t344 * t490 + t518 * t116 + t519 * t117 - t112 * t373 - t149 * t86 - t120 * t174 - t119 * t175 + t229 * t29 - t183 * t231 + t338 * t30; t499 * t115 + (t449 - t174) * t98 + (-Ifges(5,2) * t221 + t126 + t218) * t474 - t355 + (t448 + t175) * t99 + (-t121 * t456 + t26 * t31 - t27 * t32 + (t3 * t423 + t341 * t4) * pkin(4)) * m(6) + t496 + t500 * t31 + (-Ifges(6,2) * t479 + Ifges(7,3) * t480 - t528 * t467 + t477 * t507 - t494) * t361 - t86 * t456 - t59 * t85 + t51 * t455 + (Ifges(5,5) * t220 - Ifges(5,6) * t221) * t467 + t125 * t471 + (Ifges(5,1) * t220 - t440) * t472 + (t1 * t334 - t17 * t31 + t18 * t499 + t2 * t337 - t45 * t59) * m(7) - t32 * t116 - t168 * (mrSges(5,1) * t221 + mrSges(5,2) * t220) + t52 * t391 + t334 * t50 + t337 * t53 + (t121 * mrSges(6,2) + t17 * mrSges(7,2) - t26 * mrSges(6,3) - t45 * mrSges(7,3) - Ifges(6,4) * t479 - Ifges(7,5) * t480 - t508 * t467 - t509 * t477 + t530) * t145; -(-t115 - t116) * t145 + t500 * t361 + t29 + t30 + (t145 * t18 - t17 * t361 + t9) * m(7) + (t145 * t27 + t26 * t361 + t61) * m(6); -t265 * t115 + t361 * t85 + 0.2e1 * (t2 / 0.2e1 + t45 * t476 + t18 * t467) * m(7) + t53;];
tauc  = t11(:);
