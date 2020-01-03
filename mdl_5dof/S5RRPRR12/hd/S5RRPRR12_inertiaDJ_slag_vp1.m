% Calculate time derivative of joint inertia matrix for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR12_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR12_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR12_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:04
% EndTime: 2019-12-31 20:29:31
% DurationCPUTime: 16.73s
% Computational Cost: add. (19037->785), mult. (53634->1108), div. (0->0), fcn. (55978->8), ass. (0->373)
t287 = sin(qJ(2));
t290 = cos(qJ(2));
t439 = Icges(4,5) * t290;
t443 = Icges(3,4) * t290;
t502 = -t439 + t443 + (Icges(3,1) + Icges(4,1)) * t287;
t501 = t502 * qJD(2);
t460 = sin(qJ(4));
t461 = cos(qJ(4));
t238 = t287 * t461 - t290 * t460;
t288 = sin(qJ(1));
t291 = cos(qJ(1));
t314 = -t287 * t460 - t290 * t461;
t221 = t314 * t288;
t286 = sin(qJ(5));
t289 = cos(qJ(5));
t187 = t221 * t286 + t289 * t291;
t188 = -t221 * t289 + t286 * t291;
t220 = t238 * t288;
t110 = Icges(6,5) * t188 + Icges(6,6) * t187 - Icges(6,3) * t220;
t112 = Icges(6,4) * t188 + Icges(6,2) * t187 - Icges(6,6) * t220;
t114 = Icges(6,1) * t188 + Icges(6,4) * t187 - Icges(6,5) * t220;
t223 = t314 * t291;
t189 = t223 * t286 - t288 * t289;
t190 = -t223 * t289 - t286 * t288;
t222 = t238 * t291;
t46 = -t110 * t222 + t112 * t189 + t114 * t190;
t111 = Icges(6,5) * t190 + Icges(6,6) * t189 - Icges(6,3) * t222;
t113 = Icges(6,4) * t190 + Icges(6,2) * t189 - Icges(6,6) * t222;
t115 = Icges(6,1) * t190 + Icges(6,4) * t189 - Icges(6,5) * t222;
t47 = -t111 * t222 + t113 * t189 + t115 * t190;
t27 = t288 * t47 - t46 * t291;
t153 = -Icges(5,5) * t221 + Icges(5,6) * t220 + Icges(5,3) * t291;
t155 = -Icges(5,4) * t221 + Icges(5,2) * t220 + Icges(5,6) * t291;
t157 = -Icges(5,1) * t221 + Icges(5,4) * t220 + Icges(5,5) * t291;
t308 = -t153 * t288 + t155 * t222 - t157 * t223;
t156 = -Icges(5,4) * t223 + Icges(5,2) * t222 - Icges(5,6) * t288;
t158 = -Icges(5,1) * t223 + Icges(5,4) * t222 - Icges(5,5) * t288;
t154 = -Icges(5,5) * t223 + Icges(5,6) * t222 - Icges(5,3) * t288;
t433 = t154 * t288;
t69 = t156 * t222 - t158 * t223 - t433;
t481 = t288 * t69 - t291 * t308 + t27;
t44 = -t110 * t220 + t112 * t187 + t114 * t188;
t45 = -t111 * t220 + t113 * t187 + t115 * t188;
t26 = t288 * t45 - t44 * t291;
t309 = t153 * t291 + t220 * t155 - t221 * t157;
t432 = t154 * t291;
t68 = t220 * t156 - t221 * t158 + t432;
t482 = t288 * t68 - t291 * t309 + t26;
t500 = (t288 * t482 + t291 * t481) * qJD(1);
t180 = Icges(5,4) * t238 + Icges(5,2) * t314;
t491 = -t180 / 0.2e1;
t499 = -t288 / 0.2e1;
t484 = t288 / 0.2e1;
t498 = -t291 / 0.2e1;
t465 = t291 / 0.2e1;
t497 = -qJD(1) / 0.2e1;
t496 = qJD(1) / 0.2e1;
t184 = rSges(5,1) * t238 + rSges(5,2) * t314;
t284 = t288 ^ 2;
t285 = t291 ^ 2;
t410 = t284 + t285;
t495 = t184 * t410;
t377 = qJD(4) * t460;
t378 = qJD(4) * t461;
t379 = qJD(2) * t460;
t380 = qJD(2) * t461;
t177 = (-t377 + t379) * t290 + (t378 - t380) * t287;
t493 = -t177 / 0.2e1;
t179 = Icges(5,5) * t238 + Icges(5,6) * t314;
t492 = -t179 / 0.2e1;
t181 = Icges(5,1) * t238 + Icges(5,4) * t314;
t490 = t181 / 0.2e1;
t489 = t220 / 0.2e1;
t488 = -t222 / 0.2e1;
t487 = t223 / 0.2e1;
t486 = -t314 / 0.2e1;
t485 = -t238 / 0.2e1;
t355 = -Icges(3,2) * t287 + t443;
t204 = Icges(3,6) * t288 + t291 * t355;
t444 = Icges(3,4) * t287;
t359 = Icges(3,1) * t290 - t444;
t208 = Icges(3,5) * t288 + t291 * t359;
t340 = t204 * t287 - t208 * t290;
t323 = t340 * t288;
t203 = -Icges(3,6) * t291 + t288 * t355;
t207 = -Icges(3,5) * t291 + t288 * t359;
t341 = t203 * t287 - t207 * t290;
t324 = t341 * t291;
t351 = Icges(4,3) * t287 + t439;
t198 = Icges(4,6) * t288 + t291 * t351;
t440 = Icges(4,5) * t287;
t357 = Icges(4,1) * t290 + t440;
t206 = Icges(4,4) * t288 + t291 * t357;
t342 = t198 * t287 + t206 * t290;
t325 = t342 * t288;
t197 = -Icges(4,6) * t291 + t288 * t351;
t205 = -Icges(4,4) * t291 + t288 * t357;
t343 = t197 * t287 + t205 * t290;
t326 = t343 * t291;
t363 = -t188 * rSges(6,1) - t187 * rSges(6,2);
t116 = -t220 * rSges(6,3) - t363;
t454 = t221 * pkin(4);
t425 = -t220 * pkin(8) + t116 - t454;
t480 = t425 * t291;
t479 = t238 * qJD(1);
t428 = t287 * t291;
t478 = -rSges(3,2) * t428 + t288 * rSges(3,3);
t352 = Icges(3,5) * t290 - Icges(3,6) * t287;
t199 = -Icges(3,3) * t291 + t288 * t352;
t353 = Icges(4,4) * t290 + Icges(4,6) * t287;
t201 = -Icges(4,2) * t291 + t288 * t353;
t300 = t314 * qJD(2);
t426 = t290 * t291;
t139 = -t288 * t479 - t291 * t300 - t377 * t428 - t378 * t426;
t295 = (-qJD(2) + qJD(4)) * t238;
t406 = qJD(1) * t288;
t140 = t291 * t295 + t314 * t406;
t405 = qJD(1) * t291;
t98 = -qJD(5) * t190 - t140 * t286 - t289 * t405;
t99 = qJD(5) * t189 + t140 * t289 - t286 * t405;
t52 = Icges(6,5) * t99 + Icges(6,6) * t98 - Icges(6,3) * t139;
t53 = Icges(6,4) * t99 + Icges(6,2) * t98 - Icges(6,6) * t139;
t54 = Icges(6,1) * t99 + Icges(6,4) * t98 - Icges(6,5) * t139;
t10 = -t111 * t139 + t113 * t98 + t115 * t99 + t189 * t53 + t190 * t54 - t222 * t52;
t427 = t288 * t290;
t141 = -t287 * t288 * t379 - qJD(4) * t221 - t291 * t479 - t380 * t427;
t147 = -Icges(6,3) * t314 + (Icges(6,5) * t289 - Icges(6,6) * t286) * t238;
t148 = -Icges(6,6) * t314 + (Icges(6,4) * t289 - Icges(6,2) * t286) * t238;
t149 = -Icges(6,5) * t314 + (Icges(6,1) * t289 - Icges(6,4) * t286) * t238;
t398 = qJD(5) * t238;
t383 = t286 * t398;
t178 = qJD(4) * t314 - t300;
t429 = t178 * t289;
t329 = -t383 + t429;
t382 = t289 * t398;
t430 = t178 * t286;
t330 = -t382 - t430;
t79 = Icges(6,5) * t329 + Icges(6,6) * t330 + Icges(6,3) * t177;
t80 = Icges(6,4) * t329 + Icges(6,2) * t330 + Icges(6,6) * t177;
t81 = Icges(6,1) * t329 + Icges(6,4) * t330 + Icges(6,5) * t177;
t16 = -t139 * t147 + t148 * t98 + t149 * t99 + t189 * t80 + t190 * t81 - t222 * t79;
t142 = -qJD(1) * t223 + t288 * t295;
t100 = -qJD(5) * t188 - t142 * t286 - t289 * t406;
t101 = qJD(5) * t187 + t142 * t289 - t286 * t406;
t310 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t141;
t311 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t141;
t312 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t141;
t293 = -t139 * t110 + t98 * t112 + t99 * t114 + t189 * t311 + t190 * t312 - t222 * t310;
t61 = -t147 * t222 + t148 * t189 + t149 * t190;
t1 = -t10 * t222 - t47 * t139 + t46 * t141 - t16 * t314 + t61 * t177 - t220 * t293;
t349 = -t112 * t286 + t114 * t289;
t12 = t177 * t110 - t314 * t310 + t349 * t178 + (t289 * t312 - t286 * t311 + (-t112 * t289 - t114 * t286) * qJD(5)) * t238;
t348 = -t113 * t286 + t115 * t289;
t13 = t111 * t177 - t314 * t52 + t348 * t178 + (-t286 * t53 + t289 * t54 + (-t113 * t289 - t115 * t286) * qJD(5)) * t238;
t60 = -t147 * t220 + t148 * t187 + t149 * t188;
t18 = -t220 * t44 - t222 * t45 - t314 * t60;
t19 = -t220 * t46 - t222 * t47 - t314 * t61;
t11 = t100 * t113 + t101 * t115 + t111 * t141 + t187 * t53 + t188 * t54 - t220 * t52;
t17 = t100 * t148 + t101 * t149 + t141 * t147 + t187 * t80 + t188 * t81 - t220 * t79;
t292 = t100 * t112 + t101 * t114 + t141 * t110 + t187 * t311 + t188 * t312 - t220 * t310;
t2 = -t11 * t222 - t45 * t139 + t44 * t141 - t17 * t314 + t60 * t177 - t220 * t292;
t3 = t10 * t288 - t293 * t291 + (t288 * t46 + t291 * t47) * qJD(1);
t4 = t11 * t288 - t292 * t291 + (t288 * t44 + t291 * t45) * qJD(1);
t50 = -t110 * t314 + t238 * t349;
t51 = -t111 * t314 + t238 * t348;
t477 = -(t18 * t484 + t19 * t465) * qJD(1) + t139 * t27 / 0.2e1 - t141 * t26 / 0.2e1 + (t51 * t288 - t291 * t50) * t493 + t4 * t489 + t222 * t3 / 0.2e1 + t314 * (-t12 * t291 + t13 * t288 + (t50 * t288 + t51 * t291) * qJD(1)) / 0.2e1 + t1 * t499 + t2 * t465;
t117 = t190 * rSges(6,1) + t189 * rSges(6,2) - t222 * rSges(6,3);
t424 = -t223 * pkin(4) - pkin(8) * t222 + t117;
t58 = -t288 * t425 - t291 * t424;
t55 = t99 * rSges(6,1) + t98 * rSges(6,2) - t139 * rSges(6,3);
t452 = t140 * pkin(4) - t139 * pkin(8) + t55;
t455 = t142 * pkin(4);
t364 = -t101 * rSges(6,1) - t100 * rSges(6,2);
t56 = t141 * rSges(6,3) - t364;
t476 = -(t141 * pkin(8) + t455 + t56) * t288 - t452 * t291;
t475 = 2 * m(3);
t474 = 2 * m(4);
t473 = 2 * m(5);
t472 = 2 * m(6);
t471 = m(4) / 0.2e1;
t470 = m(5) / 0.2e1;
t469 = m(6) / 0.2e1;
t468 = -pkin(2) - pkin(3);
t301 = Icges(5,5) * t142 - Icges(5,6) * t141 - Icges(5,3) * t406;
t302 = -Icges(5,4) * t142 + Icges(5,2) * t141 + Icges(5,6) * t406;
t303 = Icges(5,1) * t142 - Icges(5,4) * t141 - Icges(5,5) * t406;
t83 = Icges(5,5) * t140 + Icges(5,6) * t139 - Icges(5,3) * t405;
t84 = Icges(5,4) * t140 + Icges(5,2) * t139 - Icges(5,6) * t405;
t85 = Icges(5,1) * t140 + Icges(5,4) * t139 - Icges(5,5) * t405;
t6 = (t139 * t156 + t140 * t158 + t222 * t84 - t223 * t85 - t288 * t83) * t288 - (t139 * t155 + t140 * t157 - t153 * t405 - t222 * t302 - t223 * t303 - t288 * t301) * t291 + (t69 * t291 + (t308 - t432) * t288) * qJD(1);
t467 = t6 + t3;
t7 = (-t141 * t156 + t142 * t158 + t220 * t84 - t221 * t85 + t291 * t83) * t288 - (-t141 * t155 + t142 * t157 - t153 * t406 - t220 * t302 - t221 * t303 + t291 * t301) * t291 + (t68 * t291 + (t309 - t433) * t288) * qJD(1);
t466 = -t7 - t4;
t464 = -rSges(4,1) - pkin(2);
t463 = -rSges(5,3) - pkin(7);
t462 = -rSges(6,3) - pkin(8);
t257 = rSges(3,1) * t287 + rSges(3,2) * t290;
t459 = m(3) * t257;
t458 = m(5) * t184;
t457 = pkin(7) * t288;
t456 = pkin(7) * t291;
t281 = t288 * pkin(6);
t315 = t238 * t289 * t81 + t177 * t147 - t148 * t430 + t149 * t429 - t314 * t79;
t447 = t286 * t80;
t65 = -t147 * t314 + (-t148 * t286 + t149 * t289) * t238;
t453 = -((-t447 + (-t148 * t289 - t149 * t286) * qJD(5)) * t238 + t315) * t314 + t65 * t177;
t450 = rSges(4,1) * t287;
t449 = rSges(4,2) * t291;
t448 = rSges(3,3) * t291;
t280 = t288 * rSges(4,2);
t446 = -rSges(4,3) - qJ(3);
t82 = rSges(6,1) * t329 + rSges(6,2) * t330 + rSges(6,3) * t177;
t445 = pkin(4) * t178 + pkin(8) * t177 + t82;
t434 = qJ(3) * t287;
t365 = t221 * rSges(5,1) - t220 * rSges(5,2);
t159 = rSges(5,3) * t291 - t365;
t431 = t159 * t291;
t423 = t140 * rSges(5,1) + t139 * rSges(5,2);
t150 = -rSges(6,3) * t314 + (rSges(6,1) * t289 - rSges(6,2) * t286) * t238;
t422 = pkin(4) * t238 - pkin(8) * t314 + t150;
t421 = -t223 * rSges(5,1) + t222 * rSges(5,2);
t362 = pkin(2) * t290 + t434;
t228 = t362 * t288;
t229 = pkin(2) * t426 + qJ(3) * t428;
t420 = t288 * t228 + t291 * t229;
t224 = qJD(2) * t362 - qJD(3) * t290;
t367 = rSges(4,1) * t290 + rSges(4,3) * t287;
t419 = -t367 * qJD(2) - t224;
t273 = pkin(3) * t426;
t242 = t273 - t457;
t418 = -t229 - t242;
t255 = pkin(2) * t287 - qJ(3) * t290;
t230 = t255 * t406;
t388 = t287 * t406;
t417 = pkin(3) * t388 + t230;
t256 = -rSges(4,3) * t290 + t450;
t416 = -t255 - t256;
t400 = qJD(2) * t291;
t384 = t290 * t400;
t399 = qJD(3) * t287;
t415 = qJ(3) * t384 + t291 * t399;
t414 = rSges(4,2) * t405 + rSges(4,3) * t384;
t413 = rSges(3,2) * t388 + rSges(3,3) * t405;
t402 = qJD(2) * t288;
t386 = t287 * t402;
t412 = pkin(3) * t386 + pkin(7) * t406;
t411 = t291 * pkin(1) + t281;
t409 = qJD(1) * t184;
t200 = Icges(3,3) * t288 + t291 * t352;
t408 = qJD(1) * t200;
t202 = Icges(4,2) * t288 + t291 * t353;
t407 = qJD(1) * t202;
t403 = qJD(2) * t287;
t401 = qJD(2) * t290;
t396 = t13 / 0.2e1 + t16 / 0.2e1;
t395 = t17 / 0.2e1 + t12 / 0.2e1;
t394 = -t51 / 0.2e1 - t61 / 0.2e1;
t393 = t60 / 0.2e1 + t50 / 0.2e1;
t264 = pkin(2) * t386;
t385 = t287 * t400;
t313 = -t290 * t406 - t385;
t387 = t290 * t405;
t392 = t288 * (pkin(2) * t387 + t288 * t399 - t264 + (t287 * t405 + t401 * t288) * qJ(3)) + t291 * (pkin(2) * t313 - qJ(3) * t388 + t415) + t228 * t405;
t278 = pkin(6) * t405;
t391 = t278 + t415;
t216 = rSges(4,1) * t426 + rSges(4,3) * t428 + t280;
t381 = (t352 + t353) * qJD(2) / 0.2e1;
t376 = -pkin(3) * t287 - t255;
t194 = t416 * t291;
t375 = qJD(1) * t422;
t241 = pkin(3) * t427 + t456;
t374 = t288 * t241 + t291 * t242 + t420;
t373 = t411 + t229;
t370 = -t184 + t376;
t369 = -pkin(3) * t401 - t224;
t368 = rSges(3,1) * t290 - rSges(3,2) * t287;
t366 = -t142 * rSges(5,1) + t141 * rSges(5,2);
t361 = -t288 * (-rSges(5,3) * t406 - t366) - t291 * (-rSges(5,3) * t405 + t423);
t360 = t273 + t373;
t354 = Icges(3,2) * t290 + t444;
t350 = -Icges(4,3) * t290 + t440;
t160 = -rSges(5,3) * t288 + t421;
t102 = -t288 * t159 - t160 * t291;
t339 = t376 - t422;
t217 = rSges(3,1) * t426 + t478;
t133 = rSges(5,1) * t178 - rSges(5,2) * t177;
t338 = -t133 + t369;
t337 = -pkin(1) - t368;
t130 = Icges(5,5) * t178 - Icges(5,6) * t177;
t131 = Icges(5,4) * t178 - Icges(5,2) * t177;
t132 = Icges(5,1) * t178 - Icges(5,4) * t177;
t336 = t130 * t465 + t131 * t489 - t221 * t132 / 0.2e1 + t141 * t491 + t142 * t490 + t406 * t492 + t155 * t493 + t178 * t157 / 0.2e1 + t302 * t486 + t238 * t303 / 0.2e1 + t395;
t335 = t130 * t484 + t131 * t488 + t132 * t487 + t139 * t491 - t140 * t181 / 0.2e1 + t179 * t405 / 0.2e1 + t156 * t177 / 0.2e1 - t158 * t178 / 0.2e1 + t84 * t486 + t85 * t485 - t396;
t334 = t155 * t486 + t157 * t485 + t220 * t491 + t221 * t490 + t291 * t492 - t393;
t333 = t156 * t486 + t158 * t485 + t179 * t484 + t180 * t488 + t181 * t487 + t394;
t144 = t370 * t291;
t331 = t288 * (pkin(3) * t387 - t412) + t291 * (pkin(3) * t313 - pkin(7) * t405) + t241 * t405 + t392;
t328 = t369 - t445;
t327 = qJD(2) * t257;
t320 = qJD(2) * t354;
t319 = qJD(2) * (-Icges(4,4) * t287 + Icges(4,6) * t290);
t318 = qJD(2) * (-Icges(3,5) * t287 - Icges(3,6) * t290);
t317 = qJD(2) * t350;
t95 = t339 * t291;
t316 = t290 * t468 - pkin(1) - t434;
t307 = t316 * t288;
t306 = t316 * t291;
t305 = t287 * t446 + t290 * t464 - pkin(1);
t304 = t385 * t468 + t391;
t299 = t305 * t288;
t298 = t307 - t456;
t297 = t264 + (-qJ(3) * t401 - t399) * t288 + t412;
t296 = t291 * t463 + t307;
t282 = t291 * pkin(6);
t240 = t368 * qJD(2);
t215 = t288 * t368 - t448;
t214 = t288 * t367 - t449;
t193 = t416 * t288;
t192 = t217 + t411;
t191 = t288 * t337 + t282 + t448;
t168 = t288 * t319 + t407;
t167 = -qJD(1) * t201 + t291 * t319;
t166 = t288 * t318 + t408;
t165 = -qJD(1) * t199 + t291 * t318;
t162 = t373 + t216;
t161 = t282 + t299 + t449;
t146 = t257 * t402 + ((-rSges(3,3) - pkin(6)) * t288 + t337 * t291) * qJD(1);
t145 = rSges(3,1) * t313 - rSges(3,2) * t384 - pkin(1) * t406 + t278 + t413;
t143 = t370 * t288;
t129 = qJD(1) * t194 + t288 * t419;
t128 = t256 * t406 + t291 * t419 + t230;
t127 = t288 * t200 - t291 * t340;
t126 = t288 * t199 - t324;
t125 = t288 * t202 + t291 * t342;
t124 = t288 * t201 + t326;
t123 = -t200 * t291 - t323;
t122 = -t199 * t291 - t288 * t341;
t121 = -t202 * t291 + t325;
t120 = -t201 * t291 + t288 * t343;
t119 = t288 * t463 + t360 + t421;
t118 = t282 + t296 + t365;
t108 = t288 * t214 + t216 * t291 + t420;
t106 = t422 * t291;
t105 = t422 * t288;
t104 = t264 + (-t399 + (t290 * t446 + t450) * qJD(2)) * t288 + ((-rSges(4,2) - pkin(6)) * t288 + t305 * t291) * qJD(1);
t103 = qJD(1) * t299 + t385 * t464 + t391 + t414;
t94 = t339 * t288;
t78 = t360 + t424 - t457;
t77 = -t220 * t462 + t282 + t298 + t363 + t454;
t75 = qJD(1) * t144 + t288 * t338;
t74 = t184 * t406 + t291 * t338 + t417;
t72 = -t117 * t314 + t150 * t222;
t71 = t116 * t314 - t150 * t220;
t70 = -t102 + t374;
t64 = -t116 * t222 + t117 * t220;
t63 = ((rSges(5,3) - pkin(6)) * t288 + t306) * qJD(1) + t297 + t366;
t62 = qJD(1) * t296 + t304 + t423;
t57 = t291 * t414 + (-t256 * t284 - t285 * t450) * qJD(2) + (t291 * t214 + (-t216 - t229 + t280) * t288) * qJD(1) + t392;
t49 = t288 * t445 + t291 * t375;
t48 = t291 * t445 - t406 * t422;
t43 = t374 - t58;
t40 = qJD(1) * t95 + t288 * t328;
t39 = t288 * t375 + t291 * t328 + t417;
t36 = (t160 * t288 - t431) * qJD(1) + t361;
t31 = -t455 + t462 * t141 + (t306 - t281) * qJD(1) + t297 + t364;
t30 = qJD(1) * t298 + t304 + t452;
t28 = (t431 + (-t160 + t418) * t288) * qJD(1) + t331 - t361;
t25 = -t116 * t177 + t141 * t150 - t220 * t82 + t314 * t56;
t24 = t117 * t177 + t139 * t150 + t222 * t82 - t314 * t55;
t20 = -t116 * t139 - t117 * t141 + t220 * t55 - t222 * t56;
t15 = (t288 * t424 - t480) * qJD(1) + t476;
t14 = (t480 + (t418 - t424) * t288) * qJD(1) + t331 - t476;
t5 = [t315 + t314 * t131 - t177 * t180 + t178 * t181 + (t145 * t192 + t146 * t191) * t475 + (t103 * t162 + t104 * t161) * t474 + (t118 * t63 + t119 * t62) * t473 + (t30 * t78 + t31 * t77) * t472 - t148 * t382 - t149 * t383 + (-t447 + t132) * t238 + (t350 - t354 + t357 + t359) * t403 + (t355 - t351 + t502) * t401; m(4) * (t103 * t193 + t104 * t194 + t128 * t161 + t129 * t162) + m(5) * (t118 * t74 + t119 * t75 + t143 * t62 + t144 * t63) + m(6) * (t30 * t94 + t31 * t95 + t39 * t77 + t40 * t78) + (m(3) * (-t146 * t257 - t191 * t240) + t381 * t291 + (t198 * t496 + t204 * t497 + t317 * t499 + t320 * t484) * t290 + ((t206 + t208) * t497 + t501 * t484) * t287 - t336) * t291 + (m(3) * (-t145 * t257 - t192 * t240) + t381 * t288 + (t197 * t496 + t203 * t497 + t317 * t465 + t320 * t498) * t290 + (t501 * t498 + (t205 + t207) * t497) * t287 - t335) * t288 + (-t326 / 0.2e1 + t324 / 0.2e1 + t325 / 0.2e1 - t323 / 0.2e1) * qJD(2) + ((-t192 * t459 + (-t198 / 0.2e1 + t204 / 0.2e1) * t290 + (t206 / 0.2e1 + t208 / 0.2e1) * t287 - t333) * t291 + (t191 * t459 + (-t197 / 0.2e1 + t203 / 0.2e1) * t290 + (t205 / 0.2e1 + t207 / 0.2e1) * t287 - t334) * t288) * qJD(1); -t291 * t4 + t288 * t3 + (t43 * t14 + t39 * t95 + t40 * t94) * t472 - t291 * t7 + (t143 * t75 + t144 * t74 + t28 * t70) * t473 + t288 * t6 + (t108 * t57 + t128 * t194 + t129 * t193) * t474 - t291 * ((t291 * t168 + (t121 - t326) * qJD(1)) * t291 + (t120 * qJD(1) + (t198 * t401 - t206 * t403 + t407) * t288 + (-t167 + (-t197 * t290 + t205 * t287) * qJD(2) + t342 * qJD(1)) * t291) * t288) - t291 * ((t291 * t166 + (t123 + t324) * qJD(1)) * t291 + (t122 * qJD(1) + (-t204 * t401 - t208 * t403 + t408) * t288 + (-t165 + (t203 * t290 + t207 * t287) * qJD(2) - t340 * qJD(1)) * t291) * t288) + ((t288 * t215 + t217 * t291) * ((qJD(1) * t215 - t291 * t327 + t413) * t291 + (-t288 * t327 + (-t217 + t478) * qJD(1)) * t288) + t410 * t257 * t240) * t475 + t288 * ((t288 * t167 + (t124 - t325) * qJD(1)) * t288 + (t125 * qJD(1) + (-t197 * t401 + t205 * t403) * t291 + (-t168 + (t198 * t290 - t206 * t287) * qJD(2) + (t202 + t343) * qJD(1)) * t288) * t291) + t288 * ((t288 * t165 + (t126 + t323) * qJD(1)) * t288 + (t127 * qJD(1) + (t203 * t401 + t207 * t403) * t291 + (-t166 + (-t204 * t290 - t208 * t287) * qJD(2) + (t200 - t341) * qJD(1)) * t288) * t291) + ((-t120 - t122) * t291 + (t121 + t123) * t288 + t482) * t406 + ((-t124 - t126) * t291 + (t125 + t127) * t288 + t481) * t405; 0.2e1 * ((t288 * t78 + t291 * t77) * t469 + (t118 * t291 + t119 * t288) * t470 + (t161 * t291 + t162 * t288) * t471) * t401 + 0.2e1 * ((t288 * t30 + t291 * t31 + t405 * t78 - t406 * t77) * t469 + (-t118 * t406 + t119 * t405 + t288 * t62 + t291 * t63) * t470 + (t103 * t288 + t104 * t291 - t161 * t406 + t162 * t405) * t471) * t287; 0.2e1 * ((t400 * t95 + t402 * t94 - t14) * t469 + (t143 * t402 + t144 * t400 - t28) * t470 + (t193 * t402 + t194 * t400 - t57) * t471) * t290 + 0.2e1 * ((qJD(2) * t43 + t288 * t40 + t291 * t39 + t405 * t94 - t406 * t95) * t469 + (qJD(2) * t70 + t143 * t405 - t144 * t406 + t288 * t75 + t291 * t74) * t470 + (qJD(2) * t108 + t128 * t291 + t129 * t288 + t193 * t405 - t194 * t406) * t471) * t287; 0.4e1 * (t471 + t470 + t469) * (-0.1e1 + t410) * t287 * t401; m(6) * (t105 * t30 + t106 * t31 + t48 * t77 + t49 * t78) + (m(5) * (t118 * t133 + t184 * t63) + (t119 * t458 + t333) * qJD(1) + t336) * t291 + (m(5) * (t119 * t133 + t184 * t62) + (-t118 * t458 + t334) * qJD(1) + t335) * t288; m(6) * (t105 * t40 + t106 * t39 + t58 * t14 + t15 * t43 + t48 * t95 + t49 * t94) + m(5) * (t102 * t28 + t36 * t70) + (m(5) * (t133 * t144 + t143 * t409 + t184 * t74) - t466) * t291 + (m(5) * (t133 * t143 - t144 * t409 + t184 * t75) - t467) * t288 - t500; 0.2e1 * ((qJD(2) * t495 - t36) * t470 + (t105 * t402 + t106 * t400 - t15) * t469) * t290 + 0.2e1 * ((qJD(2) * t102 + t133 * t410) * t470 + (qJD(2) * t58 + t105 * t405 - t106 * t406 + t288 * t49 + t291 * t48) * t469) * t287; t466 * t291 + t467 * t288 + (t102 * t36 + t133 * t495) * t473 + (t105 * t49 + t106 * t48 + t58 * t15) * t472 + t500; m(6) * (t24 * t78 + t25 * t77 + t30 * t72 + t31 * t71) - t396 * t222 - t395 * t220 + t393 * t141 + t394 * t139 + t453; m(6) * (t64 * t14 + t20 * t43 + t24 * t94 + t25 * t95 + t39 * t71 + t40 * t72) - t477; m(6) * ((-t20 + (t288 * t72 + t291 * t71) * qJD(2)) * t290 + (qJD(2) * t64 + t24 * t288 + t25 * t291 + (-t288 * t71 + t291 * t72) * qJD(1)) * t287); m(6) * (t105 * t24 + t106 * t25 + t64 * t15 + t20 * t58 + t48 * t71 + t49 * t72) + t477; (t64 * t20 + t24 * t72 + t25 * t71) * t472 - t139 * t19 - t222 * t1 + t141 * t18 - t220 * t2 + t177 * (-t220 * t50 - t222 * t51 - t314 * t65) - t314 * (-t12 * t220 - t13 * t222 - t51 * t139 + t50 * t141 + t453);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t5(1), t5(2), t5(4), t5(7), t5(11); t5(2), t5(3), t5(5), t5(8), t5(12); t5(4), t5(5), t5(6), t5(9), t5(13); t5(7), t5(8), t5(9), t5(10), t5(14); t5(11), t5(12), t5(13), t5(14), t5(15);];
Mq = res;
