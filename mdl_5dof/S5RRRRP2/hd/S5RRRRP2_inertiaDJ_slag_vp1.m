% Calculate time derivative of joint inertia matrix for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRP2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP2_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:47:34
% EndTime: 2019-12-05 18:47:53
% DurationCPUTime: 11.26s
% Computational Cost: add. (12695->465), mult. (10040->622), div. (0->0), fcn. (7628->8), ass. (0->282)
t511 = Icges(5,1) + Icges(6,1);
t258 = qJ(3) + qJ(4);
t250 = sin(t258);
t510 = (Icges(5,4) + Icges(6,4)) * t250;
t508 = Icges(5,5) + Icges(6,5);
t507 = Icges(5,6) + Icges(6,6);
t252 = cos(t258);
t406 = Icges(6,4) * t252;
t307 = -Icges(6,2) * t250 + t406;
t408 = Icges(5,4) * t252;
t308 = -Icges(5,2) * t250 + t408;
t509 = t307 + t308;
t505 = t511 * t252 - t510;
t259 = qJ(1) + qJ(2);
t251 = sin(t259);
t253 = cos(t259);
t473 = t507 * t251 + t509 * t253;
t471 = t508 * t251 + t505 * t253;
t504 = (Icges(5,2) + Icges(6,2)) * t252 + t510;
t503 = t511 * t250 + t406 + t408;
t506 = Icges(5,3) + Icges(6,3);
t283 = t307 * t251;
t474 = -t251 * t308 + t507 * t253 - t283;
t472 = -t505 * t251 + t508 * t253;
t501 = -t507 * t250 + t508 * t252;
t470 = t508 * t250 + t507 * t252;
t502 = t473 * t250 - t471 * t252;
t256 = qJD(3) + qJD(4);
t500 = t503 * t256;
t499 = t504 * t256;
t498 = -t501 * t251 + t506 * t253;
t497 = t506 * t251 + t501 * t253;
t257 = qJD(1) + qJD(2);
t496 = t505 * t257;
t495 = t474 * t250 - t472 * t252;
t494 = t507 * t257 - t499;
t493 = -t508 * t257 + t500;
t492 = t509 * t256;
t491 = t505 * t256;
t490 = t502 * t251;
t489 = t501 * t257;
t468 = t504 * t250 - t503 * t252;
t488 = t470 * t256 - t506 * t257;
t487 = -t495 * t251 - t498 * t253;
t486 = -t497 * t253 - t490;
t284 = t308 * t257;
t485 = -t251 * t284 + t494 * t253 - t257 * t283;
t376 = t253 * t257;
t484 = t494 * t251 + t253 * t284 + t307 * t376;
t483 = t496 * t251 + t493 * t253;
t482 = t493 * t251 - t496 * t253;
t264 = -pkin(8) - pkin(7);
t255 = -qJ(5) + t264;
t384 = t250 * t253;
t481 = rSges(6,2) * t384 + t251 * t255;
t480 = t495 * t253;
t479 = t498 * t251 - t480;
t478 = t497 * t251 - t502 * t253;
t477 = -t489 * t251 - t488 * t253;
t476 = t488 * t251 - t489 * t253;
t239 = t251 * t264;
t379 = t252 * t253;
t414 = rSges(6,3) * t251;
t290 = -rSges(6,1) * t379 - t414;
t262 = cos(qJ(3));
t246 = t262 * pkin(3) + pkin(2);
t212 = pkin(4) * t252 + t246;
t364 = t212 - t246;
t475 = t253 * t364 + t239 - t290 - t481;
t469 = t501 * t256 + t468 * t257;
t467 = t471 * t256 + t485;
t466 = t472 * t256 - t484;
t465 = -t473 * t256 - t483;
t464 = -t474 * t256 + t482;
t357 = qJD(3) * t262;
t260 = sin(qJ(3));
t375 = t253 * t260;
t463 = t251 * t357 + t257 * t375;
t377 = t253 * t256;
t381 = t251 * t257;
t462 = -t250 * t377 - t252 * t381;
t378 = t252 * t256;
t461 = t250 * t376 + t251 * t378;
t348 = t252 * t377;
t373 = t255 * t257;
t460 = -t462 * rSges(6,1) + rSges(6,2) * t348 + t212 * t381 + t253 * t373;
t459 = -t470 * t257 + (-t491 + t499) * t252 + (t492 + t500) * t250;
t458 = t486 * t251 + t487 * t253;
t410 = Icges(4,4) * t262;
t309 = -Icges(4,2) * t260 + t410;
t148 = Icges(4,6) * t253 - t251 * t309;
t411 = Icges(4,4) * t260;
t312 = Icges(4,1) * t262 - t411;
t150 = Icges(4,5) * t253 - t251 * t312;
t298 = t148 * t260 - t150 * t262;
t457 = t253 * t298;
t385 = t250 * t251;
t362 = rSges(6,2) * t385 + t253 * rSges(6,3);
t454 = rSges(5,2) * t385 + t253 * rSges(5,3);
t358 = qJD(3) * t260;
t383 = t250 * t256;
t172 = -pkin(3) * t358 - pkin(4) * t383;
t353 = t251 * t383;
t453 = qJD(5) * t253 + (-t172 + t373) * t251 + rSges(6,1) * t353 + t461 * rSges(6,2);
t423 = pkin(2) - t246;
t446 = pkin(7) * t251 + t253 * t423;
t445 = t462 * rSges(5,1) - rSges(5,2) * t348 + t257 * t454;
t227 = Icges(4,5) * t260 + Icges(4,6) * t262;
t444 = -Icges(4,3) * t257 + qJD(3) * t227;
t228 = Icges(4,2) * t262 + t411;
t443 = -Icges(4,6) * t257 + qJD(3) * t228;
t229 = Icges(4,1) * t260 + t410;
t442 = -Icges(4,5) * t257 + qJD(3) * t229;
t203 = t309 * qJD(3);
t204 = t312 * qJD(3);
t441 = t203 * t260 - t204 * t262 - t227 * t257 + (t228 * t262 + t229 * t260) * qJD(3);
t266 = t262 * t203 + t260 * t204 - t228 * t358 + t229 * t357 + t491 * t250 + t492 * t252 + t503 * t378 - t504 * t383;
t438 = 2 * m(3);
t437 = 2 * m(4);
t436 = 2 * m(5);
t435 = 2 * m(6);
t434 = t251 / 0.2e1;
t433 = t253 / 0.2e1;
t432 = -rSges(4,3) - pkin(7);
t417 = rSges(4,2) * t260;
t420 = rSges(4,1) * t262;
t209 = (-t417 + t420) * qJD(3);
t431 = m(4) * t209;
t234 = rSges(4,1) * t260 + rSges(4,2) * t262;
t430 = m(4) * t234;
t419 = rSges(5,1) * t252;
t166 = (-rSges(5,2) * t250 + t419) * t256;
t429 = m(5) * t166;
t199 = rSges(5,1) * t250 + rSges(5,2) * t252;
t428 = m(5) * t199;
t263 = cos(qJ(1));
t427 = pkin(1) * t263;
t426 = pkin(3) * t260;
t261 = sin(qJ(1));
t424 = t261 * pkin(1);
t336 = t253 * t358;
t372 = t257 * t264;
t345 = pkin(3) * t336 + t246 * t381 + t253 * t372;
t422 = (qJD(5) * t251 + t172 * t253 + t362 * t257 + t345 - t460) * t253;
t338 = t251 * t358;
t363 = pkin(3) * t338 + t251 * t372;
t421 = -t257 * t290 + t364 * t376 + t363 - t453;
t418 = rSges(6,1) * t252;
t416 = rSges(6,2) * t250;
t415 = rSges(5,3) * t251;
t413 = pkin(1) * qJD(1);
t412 = t251 * rSges(4,3);
t244 = t253 * rSges(4,3);
t392 = t166 * t251;
t382 = t251 * t252;
t380 = t251 * t260;
t374 = t253 * t264;
t371 = t475 * t253;
t370 = -(-t255 + t264) * t253 + t364 * t251 + rSges(6,1) * t382 - t362;
t245 = t253 * pkin(7);
t117 = t251 * t423 - t245 - t374;
t134 = -rSges(5,1) * t382 + t454;
t368 = -t117 - t134;
t198 = rSges(6,1) * t250 + rSges(6,2) * t252;
t219 = pkin(4) * t385;
t367 = t198 * t381 + t257 * t219;
t365 = t463 * pkin(3);
t119 = t251 * t198 + t219;
t167 = rSges(3,1) * t381 + rSges(3,2) * t376;
t235 = rSges(4,2) * t380;
t360 = t235 + t244;
t359 = t251 ^ 2 + t253 ^ 2;
t238 = pkin(3) * t380;
t356 = t251 * t420;
t355 = t263 * t413;
t354 = pkin(3) * t357;
t346 = -t117 + t370;
t343 = rSges(5,1) * t353 + t461 * rSges(5,2);
t340 = -t253 * rSges(4,2) * t357 - rSges(4,1) * t336 - t257 * t356;
t339 = rSges(4,1) * t338 + t463 * rSges(4,2);
t335 = -t381 / 0.2e1;
t334 = t376 / 0.2e1;
t333 = -pkin(2) - t420;
t332 = -pkin(4) * t250 - t198;
t201 = -rSges(3,1) * t253 + t251 * rSges(3,2);
t331 = -t246 - t419;
t330 = -t212 - t418;
t165 = (-t416 + t418) * t256;
t59 = t461 * pkin(4) + t251 * t165 + t198 * t376;
t317 = ((t476 * t253 + (t480 - t486) * t257) * t253 + t479 * t376) * t253 + ((t477 * t251 + (-t479 + t490) * t257) * t251 + t478 * t376 + (t476 * t251 + t477 * t253 + (t471 * t251 - t472 * t253) * t383 + (t473 * t251 - t474 * t253) * t378 + (t478 + t487) * t257 + ((t465 + t483) * t251 + (-t464 + t482) * t253 + (-t472 * t251 - t471 * t253) * t257) * t252 + ((-t467 + t485) * t251 + (t466 + t484) * t253 + (t474 * t251 + t473 * t253) * t257) * t250) * t253) * t251;
t316 = -pkin(4) * t378 - t165;
t168 = -rSges(3,1) * t376 + rSges(3,2) * t381;
t200 = -rSges(3,1) * t251 - rSges(3,2) * t253;
t306 = Icges(4,5) * t262 - Icges(4,6) * t260;
t299 = t148 * t262 + t150 * t260;
t149 = Icges(4,6) * t251 + t253 * t309;
t151 = Icges(4,5) * t251 + t253 * t312;
t297 = t149 * t262 + t151 * t260;
t296 = t149 * t260 - t151 * t262;
t292 = t228 * t260 - t229 * t262;
t291 = -rSges(5,1) * t379 - t415;
t288 = t312 * t257;
t285 = t309 * t257;
t282 = t306 * t257;
t277 = t296 * t251;
t275 = t253 * t331 - t415;
t274 = t253 * t330 - t414;
t271 = t306 * qJD(3) + t257 * t292;
t112 = t251 * t333 + t245 + t360;
t270 = t251 * t432 + t253 * t333;
t269 = t458 * t381 + t317;
t236 = rSges(4,2) * t375;
t113 = t236 + t270;
t218 = rSges(5,2) * t384;
t105 = t218 + t239 + t275;
t100 = t274 + t481;
t99 = t251 * t330 - t253 * t255 + t362;
t104 = t251 * t331 - t374 + t454;
t233 = pkin(2) * t381;
t64 = t233 + (t253 * t432 - t235) * t257 - t340;
t267 = (t465 * t250 + t469 * t251 + t467 * t252 - t459 * t253) * t434 + (t464 * t250 + t459 * t251 + t466 * t252 + t469 * t253) * t433 + (t472 * t250 + t468 * t251 + t474 * t252 + t470 * t253) * t335 + (t471 * t250 + t470 * t251 + t473 * t252 - t468 * t253) * t334;
t42 = t345 - t445;
t65 = t257 * t270 + t339;
t43 = t257 * t275 + t343 + t363;
t30 = (-rSges(6,3) * t257 - t172) * t253 + (-t257 * t416 - qJD(5)) * t251 + t460;
t31 = t257 * t274 + t453;
t265 = t267 + (-qJD(3) * t296 + t271 * t251 - t253 * t441 + t260 * (-t251 * t288 - t253 * t442) + t262 * (-t251 * t285 - t253 * t443)) * t434 + (-qJD(3) * t298 + t251 * t441 + t271 * t253 + t260 * (t251 * t442 - t253 * t288) + t262 * (t251 * t443 - t253 * t285)) * t433 + (t227 * t253 + t251 * t292 + t299) * t335 + (t227 * t251 - t253 * t292 + t297) * t334;
t247 = t261 * t413;
t210 = t257 * t238;
t171 = t201 - t427;
t170 = t200 - t424;
t153 = t253 * t420 - t236 + t412;
t152 = -t356 + t360;
t147 = Icges(4,3) * t251 + t253 * t306;
t146 = Icges(4,3) * t253 - t251 * t306;
t140 = t168 - t355;
t139 = t247 + t167;
t138 = (-t199 - t426) * t253;
t137 = t199 * t251 + t238;
t136 = -t218 - t291;
t120 = t332 * t253;
t118 = -t239 - t446;
t116 = t253 * t136;
t114 = t253 * t118;
t111 = (t332 - t426) * t253;
t110 = t238 + t119;
t109 = t113 - t427;
t108 = t112 - t424;
t102 = t105 - t427;
t101 = t104 - t424;
t98 = t257 * t446 + t363;
t93 = t251 * t444 - t253 * t282;
t92 = -t251 * t282 - t253 * t444;
t91 = t100 - t427;
t90 = t99 - t424;
t89 = t253 * (-pkin(7) * t376 + t233 - t345);
t86 = t257 * t291 + t343;
t82 = -t134 * t251 + t116;
t69 = t199 * t376 + t365 + t392;
t68 = t199 * t381 + t210 + (-t166 - t354) * t253;
t67 = t253 * t445;
t58 = t253 * t316 + t367;
t53 = t65 - t355;
t52 = t247 + t64;
t50 = t147 * t251 - t253 * t296;
t49 = t146 * t251 - t457;
t48 = t147 * t253 + t277;
t47 = t146 * t253 + t251 * t298;
t46 = t59 + t365;
t45 = t210 + (t316 - t354) * t253 + t367;
t41 = t43 - t355;
t40 = t247 + t42;
t29 = t31 - t355;
t28 = t247 + t30;
t27 = t251 * t368 + t114 + t116;
t24 = t251 * t370 + t371;
t21 = -t134 * t376 + t67 + (-t136 * t257 - t86) * t251;
t10 = t251 * t346 + t114 + t371;
t7 = t67 + t89 + t368 * t376 + (-t86 - t98 + (-t118 - t136) * t257) * t251;
t6 = t370 * t376 + (-t257 * t475 + t421) * t251 + t422;
t5 = t89 + t346 * t376 + (-t98 + (-t118 - t475) * t257 + t421) * t251 + t422;
t1 = [(t139 * t171 + t140 * t170) * t438 + (t108 * t53 + t109 * t52) * t437 + (t101 * t41 + t102 * t40) * t436 + (t28 * t91 + t29 * t90) * t435 + t266; t266 + m(3) * (t139 * t201 + t140 * t200 + t167 * t171 + t168 * t170) + m(4) * (t108 * t65 + t109 * t64 + t112 * t53 + t113 * t52) + m(5) * (t101 * t43 + t102 * t42 + t104 * t41 + t105 * t40) + m(6) * (t100 * t28 + t29 * t99 + t30 * t91 + t31 * t90); t266 + (t100 * t30 + t31 * t99) * t435 + (t104 * t43 + t105 * t42) * t436 + (t112 * t65 + t113 * t64) * t437 + (t167 * t201 + t168 * t200) * t438; ((t109 * t257 - t53) * t253 + (t108 * t257 + t52) * t251) * t430 + m(6) * (t110 * t28 + t111 * t29 + t45 * t90 + t46 * t91) + m(5) * (t101 * t68 + t102 * t69 + t137 * t40 + t138 * t41) + (-t108 * t253 + t109 * t251) * t431 + t265; t265 + (-t112 * t253 + t113 * t251) * t431 + ((t113 * t257 - t65) * t253 + (t112 * t257 + t64) * t251) * t430 + m(6) * (t100 * t46 + t110 * t30 + t111 * t31 + t45 * t99) + m(5) * (t104 * t68 + t105 * t69 + t137 * t42 + t138 * t43); (t10 * t5 + t110 * t46 + t111 * t45) * t435 + (t137 * t69 + t138 * t68 + t27 * t7) * t436 + ((-t152 * t251 + t153 * t253) * (t253 * t340 - t251 * t339 + ((-t152 + t244) * t253 + (t412 - t153 + (t417 + t420) * t253) * t251) * t257) + t359 * t234 * t209) * t437 + (t251 * t50 + t253 * t49) * t376 + t251 * ((t251 * t92 + (-t49 + t277) * t257) * t251 + (t50 * t257 + (-t148 * t357 - t150 * t358) * t253 + (-t297 * qJD(3) + t257 * t298 + t93) * t251) * t253) + t253 * ((t253 * t93 + (t48 + t457) * t257) * t253 + (-t47 * t257 + (t149 * t357 + t151 * t358) * t251 + (t299 * qJD(3) + t257 * t296 + t92) * t253) * t251) + t317 + (-t251 * t48 - t253 * t47 + t458) * t381; ((t102 * t257 - t41) * t253 + (t101 * t257 + t40) * t251) * t428 + (-t101 * t253 + t102 * t251) * t429 + m(6) * (t119 * t28 + t120 * t29 + t58 * t90 + t59 * t91) + t267; m(6) * (t100 * t59 + t119 * t30 + t120 * t31 + t58 * t99) + ((t105 * t257 - t43) * t253 + (t104 * t257 + t42) * t251) * t428 + (-t104 * t253 + t105 * t251) * t429 + t267; m(6) * (t10 * t6 + t110 * t59 + t111 * t58 + t119 * t46 + t120 * t45 + t24 * t5) + m(5) * (-t138 * t166 * t253 + t137 * t392 + t21 * t27 + t7 * t82) + ((t137 * t257 - t68) * t253 + (t138 * t257 + t69) * t251) * t428 + t269; (t166 * t199 * t359 + t21 * t82) * t436 + (t119 * t59 + t120 * t58 + t24 * t6) * t435 + t269; m(6) * ((t257 * t90 + t28) * t253 + (-t257 * t91 + t29) * t251); m(6) * ((t257 * t99 + t30) * t253 + (-t100 * t257 + t31) * t251); m(6) * ((t111 * t257 + t46) * t253 + (-t110 * t257 + t45) * t251); m(6) * ((t120 * t257 + t59) * t253 + (-t119 * t257 + t58) * t251); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
