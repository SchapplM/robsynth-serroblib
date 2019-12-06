% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PPRRP2_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRRP2_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:34
% EndTime: 2019-12-05 15:08:49
% DurationCPUTime: 9.36s
% Computational Cost: add. (20265->500), mult. (27804->740), div. (0->0), fcn. (30336->6), ass. (0->304)
t458 = rSges(6,1) + pkin(4);
t453 = rSges(6,3) + qJ(5);
t261 = pkin(8) + qJ(3);
t260 = cos(t261);
t263 = cos(pkin(7));
t264 = sin(qJ(4));
t366 = t263 * t264;
t262 = sin(pkin(7));
t265 = cos(qJ(4));
t367 = t262 * t265;
t246 = t260 * t367 - t366;
t230 = Icges(6,5) * t246;
t365 = t263 * t265;
t368 = t262 * t264;
t245 = t260 * t368 + t365;
t259 = sin(t261);
t381 = t259 * t262;
t134 = Icges(6,6) * t381 + Icges(6,3) * t245 + t230;
t395 = Icges(5,4) * t246;
t140 = -Icges(5,2) * t245 + Icges(5,6) * t381 + t395;
t490 = t134 - t140;
t248 = t260 * t365 + t368;
t231 = Icges(6,5) * t248;
t247 = t260 * t366 - t367;
t380 = t259 * t263;
t135 = Icges(6,6) * t380 + Icges(6,3) * t247 + t231;
t394 = Icges(5,4) * t248;
t141 = -Icges(5,2) * t247 + Icges(5,6) * t380 + t394;
t489 = t135 - t141;
t136 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t381;
t138 = Icges(6,4) * t246 + Icges(6,2) * t381 + Icges(6,6) * t245;
t476 = t136 + t138;
t137 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t380;
t139 = Icges(6,4) * t248 + Icges(6,2) * t380 + Icges(6,6) * t247;
t475 = t137 + t139;
t389 = Icges(6,5) * t245;
t142 = Icges(6,1) * t246 + Icges(6,4) * t381 + t389;
t232 = Icges(5,4) * t245;
t144 = Icges(5,1) * t246 + Icges(5,5) * t381 - t232;
t488 = t142 + t144;
t388 = Icges(6,5) * t247;
t143 = Icges(6,1) * t248 + Icges(6,4) * t380 + t388;
t233 = Icges(5,4) * t247;
t145 = Icges(5,1) * t248 + Icges(5,5) * t380 - t233;
t487 = t143 + t145;
t392 = Icges(5,4) * t265;
t306 = -Icges(5,2) * t264 + t392;
t204 = -Icges(5,6) * t260 + t259 * t306;
t378 = t259 * t265;
t255 = Icges(6,5) * t378;
t379 = t259 * t264;
t384 = Icges(6,6) * t260;
t474 = Icges(6,3) * t379 - t204 + t255 - t384;
t387 = Icges(6,5) * t264;
t309 = Icges(6,1) * t265 + t387;
t206 = -Icges(6,4) * t260 + t259 * t309;
t393 = Icges(5,4) * t264;
t310 = Icges(5,1) * t265 - t393;
t208 = -Icges(5,5) * t260 + t259 * t310;
t472 = t206 + t208;
t162 = -Icges(5,5) * t245 - Icges(5,6) * t246;
t164 = -Icges(6,4) * t245 + Icges(6,6) * t246;
t486 = t162 + t164;
t163 = -Icges(5,5) * t247 - Icges(5,6) * t248;
t165 = -Icges(6,4) * t247 + Icges(6,6) * t248;
t485 = t163 + t165;
t303 = Icges(5,5) * t265 - Icges(5,6) * t264;
t200 = -Icges(5,3) * t260 + t259 * t303;
t305 = Icges(6,4) * t265 + Icges(6,6) * t264;
t202 = -Icges(6,2) * t260 + t259 * t305;
t473 = t200 + t202;
t357 = Icges(5,2) * t248 - t145 + t233;
t359 = Icges(6,3) * t248 - t143 - t388;
t484 = t357 + t359;
t358 = Icges(5,2) * t246 - t144 + t232;
t360 = Icges(6,3) * t246 - t142 - t389;
t483 = t358 + t360;
t361 = -Icges(5,1) * t247 - t141 - t394;
t363 = -Icges(6,1) * t247 + t135 + t231;
t482 = t361 + t363;
t362 = -Icges(5,1) * t245 - t140 - t395;
t364 = -Icges(6,1) * t245 + t134 + t230;
t481 = t362 + t364;
t480 = t489 * t245 + t487 * t246 + t475 * t381;
t479 = t490 * t247 + t488 * t248 + t476 * t380;
t470 = t453 * t264 + t458 * t265;
t447 = t260 * rSges(6,2) - t470 * t259;
t478 = t447 * t262;
t477 = t447 * t263;
t471 = t474 * t264 + t472 * t265;
t424 = t263 ^ 2;
t425 = t262 ^ 2;
t430 = t424 + t425;
t469 = Icges(4,5) * t259 + Icges(4,6) * t260;
t468 = t483 * t245 + t481 * t246 + t486 * t381;
t467 = t484 * t245 + t482 * t246 + t485 * t381;
t466 = t483 * t247 + t481 * t248 + t486 * t380;
t465 = t484 * t247 + t482 * t248 + t485 * t380;
t318 = rSges(5,1) * t265 - rSges(5,2) * t264;
t211 = -t260 * rSges(5,3) + t259 * t318;
t195 = t211 * t262;
t197 = t211 * t263;
t462 = t472 + (t387 - t393 + (-Icges(5,2) - Icges(6,3)) * t265) * t259;
t461 = -(-Icges(5,1) * t264 - t392) * t259 + Icges(6,1) * t379 - t255 - t474;
t235 = (-Icges(5,5) * t264 - Icges(5,6) * t265) * t259;
t236 = (-Icges(6,4) * t264 + Icges(6,6) * t265) * t259;
t460 = (-t235 - t236) * t260;
t442 = t473 * t260;
t302 = Icges(6,5) * t265 + Icges(6,3) * t264;
t272 = -t259 * t302 + t384;
t182 = t272 * t262;
t188 = t204 * t262;
t190 = t206 * t262;
t192 = t208 * t262;
t297 = -t140 * t264 + t144 * t265;
t287 = -t200 * t262 - t297;
t299 = t134 * t264 + t142 * t265;
t289 = t202 * t262 + t299;
t457 = (t289 - t287) * t260 + ((-t190 - t192) * t265 + (t182 + t188) * t264 + t476) * t259;
t183 = t272 * t263;
t189 = t204 * t263;
t191 = t206 * t263;
t193 = t208 * t263;
t296 = -t141 * t264 + t145 * t265;
t286 = -t200 * t263 - t296;
t298 = t135 * t264 + t143 * t265;
t288 = t202 * t263 + t298;
t456 = (-t286 + t288) * t260 + ((-t191 - t193) * t265 + (t183 + t189) * t264 + t475) * t259;
t455 = t490 * t245 + t488 * t246 + t476 * t381;
t454 = t489 * t247 + t487 * t248 + t475 * t380;
t452 = t474 * t245 + t472 * t246 + t473 * t381;
t451 = t474 * t247 + t472 * t248 + t473 * t380;
t450 = t471 * t259 - t442;
t449 = (-t302 + t306) * t260 + (Icges(5,6) - Icges(6,6)) * t259;
t448 = (-t309 - t310) * t260 + (-Icges(6,4) - Icges(5,5)) * t259;
t446 = ((t303 + t305) * t260 + (Icges(6,2) + Icges(5,3)) * t259 - t471) * t260;
t375 = t260 * t138;
t92 = t259 * t299 - t375;
t374 = t260 * t139;
t93 = t259 * t298 - t374;
t377 = t260 * t136;
t94 = t259 * t297 - t377;
t376 = t260 * t137;
t95 = t259 * t296 - t376;
t443 = (t93 + t95) * t263 + (t92 + t94) * t262;
t441 = t479 * t262;
t440 = t480 * t263;
t147 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t381;
t149 = rSges(5,1) * t248 - rSges(5,2) * t247 + rSges(5,3) * t380;
t252 = pkin(3) * t260 + pkin(6) * t259;
t339 = t430 * t252;
t100 = t147 * t262 + t149 * t263 + t339;
t174 = -rSges(5,1) * t245 - rSges(5,2) * t246;
t178 = -rSges(5,1) * t247 - rSges(5,2) * t248;
t118 = t174 * t262 + t178 * t263;
t251 = pkin(3) * t259 - pkin(6) * t260;
t326 = -t251 + t447;
t124 = t326 * t262;
t126 = t326 * t263;
t345 = -t211 - t251;
t153 = t345 * t262;
t155 = t345 * t263;
t338 = (-t458 * t264 + t453 * t265) * t259;
t157 = t338 * t262;
t158 = t338 * t263;
t242 = (-rSges(5,1) * t264 - rSges(5,2) * t265) * t259;
t355 = rSges(6,2) * t380 + t453 * t247 + t458 * t248;
t356 = rSges(6,2) * t381 + t453 * t245 + t458 * t246;
t76 = t262 * t356 + t263 * t355 + t339;
t353 = -t458 * t247 + t453 * t248;
t354 = -t458 * t245 + t453 * t246;
t91 = t262 * t354 + t263 * t353;
t439 = -m(6) * (-t124 * t157 - t126 * t158 + t76 * t91) - m(5) * (t100 * t118 + (-t153 * t262 - t155 * t263) * t242);
t438 = -t259 / 0.2e1;
t437 = t259 / 0.2e1;
t436 = -t260 / 0.2e1;
t435 = t260 / 0.2e1;
t434 = -t262 / 0.2e1;
t409 = t262 / 0.2e1;
t408 = -t263 / 0.2e1;
t407 = t263 / 0.2e1;
t433 = (t462 * t245 + t461 * t246) * t260 + (t467 * t263 + (t460 + t468) * t262) * t259;
t432 = (t462 * t247 + t461 * t248) * t260 + ((t460 + t465) * t263 + t466 * t262) * t259;
t159 = t245 * t262 + t247 * t263;
t369 = t260 * t264;
t123 = (-t159 + t369) * t379;
t323 = t369 / 0.2e1;
t130 = (t323 - t159 / 0.2e1) * m(6);
t404 = m(6) * qJD(5);
t431 = t130 * qJD(1) + t123 * t404;
t428 = t469 * t430;
t426 = t260 ^ 2;
t423 = 2 * qJD(3);
t422 = 2 * qJD(4);
t421 = 4 * qJD(4);
t420 = m(5) / 0.2e1;
t419 = m(6) / 0.2e1;
t98 = -t259 * t478 + t260 * t356;
t99 = t259 * t477 - t260 * t355;
t313 = -t262 * t99 - t263 * t98;
t279 = -t262 * t355 + t263 * t356;
t48 = t279 * t260 + (-t262 * t477 + t263 * t478) * t259;
t344 = t259 * rSges(6,2) + t470 * t260;
t70 = (t262 * t344 - t356) * t259;
t71 = (-t263 * t344 + t355) * t259;
t77 = t279 * t259;
t418 = m(6) * (t245 * t71 + t247 * t70 + (t260 * t77 + (t313 + t48) * t259) * t264);
t417 = m(6) * (t48 * t77 + t70 * t98 + t71 * t99);
t295 = t147 * t263 - t149 * t262;
t110 = t295 * t259;
t116 = t147 * t260 + t211 * t381;
t117 = -t149 * t260 - t211 * t380;
t79 = t295 * t260 + (-t195 * t263 + t197 * t262) * t259;
t213 = t259 * rSges(5,3) + t260 * t318;
t96 = (t213 * t262 - t147) * t259;
t97 = (-t213 * t263 + t149) * t259;
t416 = m(5) * (t110 * t79 + t116 * t96 + t117 * t97);
t133 = (t245 * t263 - t247 * t262) * t259;
t258 = t259 ^ 2;
t180 = t245 * t260 + t258 * t368;
t181 = -t247 * t260 - t258 * t366;
t414 = m(6) * (t124 * t181 + t126 * t180 + t133 * t76 + t159 * t77 + t313 * t379);
t413 = m(6) * (t246 * t124 + t248 * t126 - t245 * t157 - t247 * t158 + (t264 * t91 + t265 * t76) * t259);
t406 = m(6) * qJD(3);
t405 = m(6) * qJD(4);
t343 = -t213 - t252;
t340 = t430 * t251;
t337 = qJD(4) * t259;
t268 = -t259 * t289 + t375;
t50 = t245 * t182 - t246 * t190 + t262 * t268;
t267 = -t259 * t288 + t374;
t51 = t245 * t183 - t246 * t191 + t262 * t267;
t270 = t259 * t287 + t377;
t52 = t245 * t188 - t246 * t192 + t262 * t270;
t269 = t259 * t286 + t376;
t53 = t245 * t189 - t246 * t193 + t262 * t269;
t335 = ((t53 + t51) * t263 + (t52 + t50 - t446) * t262 + t452) * t437 + ((-t442 + t455) * t262 + t448 * t246 + t449 * t245 + t440) * t435;
t54 = t247 * t182 - t248 * t190 + t263 * t268;
t55 = t247 * t183 - t248 * t191 + t263 * t267;
t56 = t247 * t188 - t248 * t192 + t263 * t270;
t57 = t247 * t189 - t248 * t193 + t263 * t269;
t333 = ((t55 + t57 - t446) * t263 + (t54 + t56) * t262 + t451) * t437 + ((-t442 + t454) * t263 + t448 * t248 + t449 * t247 + t441) * t435;
t332 = (t456 * t263 + t457 * t262 + (t264 * t449 + t265 * t448 - t473) * t260 + t450) * t438 + (t443 + t446) * t436;
t331 = t468 * t407 + t467 * t434;
t330 = t466 * t408 + t465 * t409;
t329 = (t455 * t262 + t440) * t437 + t452 * t436;
t328 = (t454 * t263 + t441) * t437 + t451 * t436;
t327 = t435 * t450 + t438 * t443;
t325 = -t252 - t344;
t324 = t378 / 0.2e1;
t249 = rSges(4,1) * t259 + rSges(4,2) * t260;
t266 = -m(5) * (t262 * t96 - t263 * t97) / 0.2e1 - m(6) * (t262 * t70 - t263 * t71) / 0.2e1;
t283 = (t157 * t263 - t158 * t262) * t419;
t18 = t283 + t266;
t19 = 0.2e1 * (t48 / 0.4e1 - t91 / 0.4e1) * m(6) + 0.2e1 * (t79 / 0.4e1 - t118 / 0.4e1) * m(5);
t301 = -t19 * qJD(1) + t18 * qJD(2);
t300 = -t124 * t262 - t126 * t263;
t282 = (t180 * t262 - t181 * t263) * t419;
t285 = m(6) * (-t246 * t263 + t248 * t262);
t104 = t282 - t285 / 0.2e1;
t128 = (t324 - t133 / 0.2e1) * m(6);
t292 = -t128 * qJD(1) + t104 * qJD(2);
t291 = -t331 - t335;
t290 = -t330 + t333;
t225 = t469 * t263;
t224 = t469 * t262;
t156 = t343 * t263;
t154 = t343 * t262;
t152 = t430 * t249;
t131 = m(6) * t323 + t159 * t419;
t129 = m(6) * t324 + t133 * t419;
t127 = t325 * t263;
t125 = t325 * t262;
t122 = -t178 * t260 - t242 * t380;
t121 = t174 * t260 + t242 * t381;
t119 = t258 * t264 * t265 + t245 * t246 + t247 * t248;
t114 = (t174 * t263 - t178 * t262) * t259;
t109 = -t195 * t262 - t197 * t263 - t340;
t103 = t282 + t285 / 0.2e1;
t102 = -t158 * t259 - t260 * t353;
t101 = t260 * t354 + t338 * t381;
t90 = t262 * t478 + t263 * t477 - t340;
t88 = (-t262 * t353 + t263 * t354) * t259;
t75 = -t260 * t163 + (t264 * t357 + t265 * t361) * t259;
t74 = -t260 * t162 + (t264 * t358 + t265 * t362) * t259;
t73 = -t260 * t165 + (t264 * t359 + t265 * t363) * t259;
t72 = -t260 * t164 + (t264 * t360 + t265 * t364) * t259;
t45 = t159 * t76 + t300 * t379;
t32 = t133 * t77 + t180 * t98 + t181 * t99;
t31 = t262 * t57 - t263 * t56;
t30 = t262 * t55 - t263 * t54;
t29 = t262 * t53 - t263 * t52;
t28 = t262 * t51 - t263 * t50;
t26 = t413 / 0.2e1;
t20 = (t118 + t79) * t420 + (t48 + t91) * t419;
t17 = t283 - t266;
t15 = t414 / 0.2e1;
t6 = t418 / 0.2e1;
t5 = t26 + t6 - t414 / 0.2e1;
t4 = t15 + t26 - t418 / 0.2e1;
t3 = t15 + t6 - t413 / 0.2e1;
t2 = t262 * t330 + t263 * t331 - t439;
t1 = t416 + t417 + (t262 * t329 + t263 * t328 + t332) * t260 + (t262 * t335 + t263 * t333 - t327) * t259;
t7 = [0, 0, t20 * qJD(4) + t131 * qJD(5) + (-m(4) * t152 / 0.2e1 + t109 * t420 + t90 * t419) * t423, t20 * qJD(3) + (t114 * t420 + t419 * t88) * t422 + t129 * qJD(5), qJD(3) * t131 + qJD(4) * t129; 0, 0, t17 * qJD(4) + ((-t154 * t263 + t156 * t262) * t420 + (-t125 * t263 + t127 * t262) * t419) * t423, t17 * qJD(3) + ((t121 * t262 - t122 * t263) * t420 + (t101 * t262 - t102 * t263) * t419) * t422 + t103 * qJD(5), t103 * qJD(4); -qJD(4) * t19 - qJD(5) * t130, t18 * qJD(4), t2 * qJD(4) + t45 * t404 + (m(6) * (t124 * t125 + t126 * t127 + t76 * t90) + m(5) * (t100 * t109 + t153 * t154 + t155 * t156) + m(4) * (-t152 + t249) * t430 * (rSges(4,1) * t260 - rSges(4,2) * t259) + (t30 + t31 - t425 * t225 + (t262 * t224 - t428) * t263) * t409 + (t28 + t29 - t424 * t224 + (t263 * t225 - t428) * t262) * t408) * qJD(3), t2 * qJD(3) + t4 * qJD(5) + (-t416 / 0.4e1 - t417 / 0.4e1) * t421 + ((t100 * t114 + t110 * t118 + t121 * t155 + t122 * t153 + (-t116 * t263 - t117 * t262) * t242) * t420 + (t101 * t126 + t102 * t124 - t157 * t99 - t158 * t98 + t76 * t88 + t77 * t91) * t419) * t422 + (t262 * t291 - t263 * t290 + t327) * t337 + t301 + (((t72 / 0.2e1 + t74 / 0.2e1 - t328) * t263 + (-t73 / 0.2e1 - t75 / 0.2e1 - t329) * t262 - t332) * t260 + t432 * t409 + t433 * t408) * qJD(4), t4 * qJD(4) + t406 * t45 - t431; qJD(3) * t19 - qJD(5) * t128, -qJD(3) * t18 + qJD(5) * t104, t1 * qJD(4) + t3 * qJD(5) + ((t124 * t71 + t125 * t99 + t126 * t70 + t127 * t98 + t48 * t76 + t77 * t90) * t419 + (t100 * t79 + t109 * t110 + t116 * t156 + t117 * t154 + t153 * t97 + t155 * t96) * t420) * t423 - t301 + (t291 * t263 + t290 * t262 + (t456 * t434 + (t480 * t262 - t455 * t263) * t409 + (t454 * t262 - t479 * t263 + t457) * t407) * t260 + ((-t94 / 0.2e1 - t92 / 0.2e1 + t31 / 0.2e1 + t30 / 0.2e1) * t263 + (t95 / 0.2e1 + t93 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1) * t262) * t259 + t439) * qJD(3), t1 * qJD(3) + (m(6) * (t98 * t101 + t99 * t102 + t77 * t88) / 0.4e1 + m(5) * (t110 * t114 + t116 * t121 + t117 * t122) / 0.4e1) * t421 + t32 * t404 + (-t235 / 0.2e1 - t236 / 0.2e1) * qJD(4) * t260 * t426 + (((t75 + t73) * t263 + (t74 + t72) * t262 + (t462 * t264 + t461 * t265) * t260) * t436 + t433 * t409 + t432 * t407) * t337, t3 * qJD(3) + t32 * t405 + (t133 * t379 + t180 * t247 + t181 * t245 - t119) * t404 + t292; qJD(3) * t130 + qJD(4) * t128, -t104 * qJD(4), (t125 * t245 + t127 * t247 - t45 + (t260 * t76 + (t300 + t90) * t259) * t264) * t406 + t5 * qJD(4) + t431, t5 * qJD(3) + (t247 * t101 + t245 * t102 + t246 * t99 + t248 * t98 + (t264 * t88 + t265 * t77) * t259 - t32) * t405 + t119 * t404 - t292, 0.4e1 * (t123 * qJD(3) / 0.4e1 + t119 * qJD(4) / 0.4e1) * m(6);];
Cq = t7;
