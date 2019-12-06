% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:41
% EndTime: 2019-12-05 15:33:03
% DurationCPUTime: 10.92s
% Computational Cost: add. (20759->497), mult. (26591->723), div. (0->0), fcn. (28410->8), ass. (0->314)
t272 = qJ(2) + pkin(8);
t269 = cos(t272);
t274 = cos(pkin(7));
t279 = cos(qJ(4));
t389 = t274 * t279;
t273 = sin(pkin(7));
t277 = sin(qJ(4));
t392 = t273 * t277;
t244 = -t269 * t392 - t389;
t390 = t274 * t277;
t391 = t273 * t279;
t245 = t269 * t391 - t390;
t268 = sin(t272);
t396 = t268 * t273;
t133 = Icges(6,5) * t245 + Icges(6,6) * t244 + Icges(6,3) * t396;
t135 = Icges(5,5) * t245 + Icges(5,6) * t244 + Icges(5,3) * t396;
t511 = t133 + t135;
t246 = -t269 * t390 + t391;
t247 = t269 * t389 + t392;
t395 = t268 * t274;
t134 = Icges(6,5) * t247 + Icges(6,6) * t246 + Icges(6,3) * t395;
t136 = Icges(5,5) * t247 + Icges(5,6) * t246 + Icges(5,3) * t395;
t510 = t134 + t136;
t416 = Icges(6,4) * t245;
t137 = Icges(6,2) * t244 + Icges(6,6) * t396 + t416;
t420 = Icges(5,4) * t245;
t139 = Icges(5,2) * t244 + Icges(5,6) * t396 + t420;
t526 = t137 + t139;
t415 = Icges(6,4) * t247;
t138 = Icges(6,2) * t246 + Icges(6,6) * t395 + t415;
t419 = Icges(5,4) * t247;
t140 = Icges(5,2) * t246 + Icges(5,6) * t395 + t419;
t525 = t138 + t140;
t230 = Icges(6,4) * t244;
t141 = Icges(6,1) * t245 + Icges(6,5) * t396 + t230;
t232 = Icges(5,4) * t244;
t143 = Icges(5,1) * t245 + Icges(5,5) * t396 + t232;
t524 = t141 + t143;
t231 = Icges(6,4) * t246;
t142 = Icges(6,1) * t247 + Icges(6,5) * t395 + t231;
t233 = Icges(5,4) * t246;
t144 = Icges(5,1) * t247 + Icges(5,5) * t395 + t233;
t523 = t142 + t144;
t413 = Icges(6,4) * t279;
t324 = -Icges(6,2) * t277 + t413;
t197 = -Icges(6,6) * t269 + t268 * t324;
t417 = Icges(5,4) * t279;
t325 = -Icges(5,2) * t277 + t417;
t199 = -Icges(5,6) * t269 + t268 * t325;
t508 = t197 + t199;
t414 = Icges(6,4) * t277;
t330 = Icges(6,1) * t279 - t414;
t201 = -Icges(6,5) * t269 + t268 * t330;
t418 = Icges(5,4) * t277;
t331 = Icges(5,1) * t279 - t418;
t203 = -Icges(5,5) * t269 + t268 * t331;
t507 = t201 + t203;
t151 = Icges(6,5) * t244 - Icges(6,6) * t245;
t153 = Icges(5,5) * t244 - Icges(5,6) * t245;
t522 = t151 + t153;
t152 = Icges(6,5) * t246 - Icges(6,6) * t247;
t154 = Icges(5,5) * t246 - Icges(5,6) * t247;
t521 = t152 + t154;
t320 = Icges(6,5) * t279 - Icges(6,6) * t277;
t193 = -Icges(6,3) * t269 + t268 * t320;
t321 = Icges(5,5) * t279 - Icges(5,6) * t277;
t195 = -Icges(5,3) * t269 + t268 * t321;
t509 = t193 + t195;
t378 = -Icges(5,2) * t247 + t144 + t233;
t380 = -Icges(6,2) * t247 + t142 + t231;
t520 = t378 + t380;
t379 = -Icges(5,2) * t245 + t143 + t232;
t381 = -Icges(6,2) * t245 + t141 + t230;
t519 = t379 + t381;
t382 = Icges(5,1) * t246 - t140 - t419;
t384 = Icges(6,1) * t246 - t138 - t415;
t518 = t382 + t384;
t383 = Icges(5,1) * t244 - t139 - t420;
t385 = Icges(6,1) * t244 - t137 - t416;
t517 = t383 + t385;
t493 = rSges(6,1) + pkin(4);
t516 = t525 * t244 + t523 * t245 + t510 * t396;
t515 = t526 * t246 + t524 * t247 + t511 * t395;
t338 = rSges(6,1) * t279 - rSges(6,2) * t277;
t436 = pkin(4) * t279;
t481 = (rSges(6,3) + qJ(5)) * t269 + (-t338 - t436) * t268;
t514 = t481 * t273;
t513 = t481 * t274;
t270 = t273 ^ 2;
t271 = t274 ^ 2;
t363 = t270 + t271;
t506 = -t508 * t277 + t507 * t279;
t505 = Icges(4,5) * t268 + Icges(4,6) * t269;
t504 = t519 * t244 + t517 * t245 + t522 * t396;
t503 = t520 * t244 + t518 * t245 + t521 * t396;
t502 = t519 * t246 + t517 * t247 + t522 * t395;
t501 = t520 * t246 + t518 * t247 + t521 * t395;
t339 = rSges(5,1) * t279 - rSges(5,2) * t277;
t206 = -rSges(5,3) * t269 + t268 * t339;
t184 = t206 * t273;
t186 = t206 * t274;
t498 = t507 + (-t414 - t418 + (-Icges(5,2) - Icges(6,2)) * t279) * t268;
t497 = t508 + (t413 + t417 + (Icges(5,1) + Icges(6,1)) * t277) * t268;
t236 = (-Icges(6,5) * t277 - Icges(6,6) * t279) * t268;
t237 = (-Icges(5,5) * t277 - Icges(5,6) * t279) * t268;
t496 = (-t236 - t237) * t269;
t474 = t509 * t269;
t278 = sin(qJ(2));
t280 = cos(qJ(2));
t494 = 0.2e1 * t278 * (Icges(3,1) - Icges(3,2)) * t280 + (-0.2e1 * t278 ^ 2 + 0.2e1 * t280 ^ 2) * Icges(3,4);
t175 = t197 * t273;
t177 = t199 * t273;
t179 = t201 * t273;
t181 = t203 * t273;
t316 = -t139 * t277 + t143 * t279;
t304 = -t195 * t273 - t316;
t318 = -t137 * t277 + t141 * t279;
t306 = -t193 * t273 - t318;
t492 = (-t306 - t304) * t269 + ((-t179 - t181) * t279 + (t175 + t177) * t277 + t511) * t268;
t176 = t197 * t274;
t178 = t199 * t274;
t180 = t201 * t274;
t182 = t203 * t274;
t315 = -t140 * t277 + t144 * t279;
t303 = -t195 * t274 - t315;
t317 = -t138 * t277 + t142 * t279;
t305 = -t193 * t274 - t317;
t491 = (-t305 - t303) * t269 + ((-t180 - t182) * t279 + (t176 + t178) * t277 + t510) * t268;
t490 = t526 * t244 + t524 * t245 + t511 * t396;
t489 = t525 * t246 + t523 * t247 + t510 * t395;
t323 = -Icges(3,5) * t278 - Icges(3,6) * t280;
t249 = t323 * t273;
t250 = t323 * t274;
t488 = t508 * t244 + t507 * t245 + t509 * t396;
t487 = t508 * t246 + t507 * t247 + t509 * t395;
t486 = t506 * t268 - t474;
t485 = (t324 + t325) * t269 + (Icges(5,6) + Icges(6,6)) * t268;
t484 = (-t330 - t331) * t269 + (-Icges(5,5) - Icges(6,5)) * t268;
t483 = -t273 * t505 + t249;
t482 = -t274 * t505 + t250;
t480 = ((t321 + t320) * t269 + (Icges(5,3) + Icges(6,3)) * t268 - t506) * t269;
t406 = t133 * t269;
t93 = t268 * t318 - t406;
t405 = t134 * t269;
t94 = t268 * t317 - t405;
t404 = t135 * t269;
t95 = t268 * t316 - t404;
t403 = t136 * t269;
t96 = t268 * t315 - t403;
t475 = (t94 + t96) * t274 + (t93 + t95) * t273;
t473 = t515 * t273;
t472 = t516 * t274;
t374 = -rSges(6,2) * t247 + t246 * t493;
t375 = -rSges(6,2) * t245 + t244 * t493;
t102 = t273 * t375 + t274 * t374;
t257 = pkin(3) * t268 - pkin(6) * t269;
t440 = pkin(2) * t278;
t350 = -t257 - t440;
t311 = t350 + t481;
t116 = t311 * t273;
t118 = t311 * t274;
t164 = rSges(5,1) * t244 - rSges(5,2) * t245;
t166 = rSges(5,1) * t246 - rSges(5,2) * t247;
t122 = t164 * t273 + t166 * t274;
t341 = -t206 + t350;
t127 = t341 * t273;
t129 = t341 * t274;
t342 = (rSges(6,2) * t279 + t493 * t277) * t268;
t169 = t342 * t273;
t170 = t342 * t274;
t243 = (-rSges(5,1) * t277 - rSges(5,2) * t279) * t268;
t258 = pkin(3) * t269 + pkin(6) * t268;
t439 = pkin(2) * t280;
t371 = t363 * t439;
t344 = t258 * t363 + t371;
t188 = qJ(5) * t268 + t269 * t436;
t386 = rSges(6,1) * t247 + rSges(6,2) * t246 + rSges(6,3) * t395 + pkin(4) * t392 + t188 * t274;
t387 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t396 - pkin(4) * t390 + t188 * t273;
t57 = t273 * t387 + t274 * t386 + t344;
t146 = rSges(5,1) * t245 + rSges(5,2) * t244 + rSges(5,3) * t396;
t148 = rSges(5,1) * t247 + rSges(5,2) * t246 + rSges(5,3) * t395;
t78 = t146 * t273 + t148 * t274 + t344;
t471 = -m(6) * (t102 * t57 + t116 * t169 + t118 * t170) - m(5) * (t122 * t78 + (-t127 * t273 - t129 * t274) * t243);
t470 = -t268 / 0.2e1;
t469 = t268 / 0.2e1;
t468 = -t269 / 0.2e1;
t467 = t269 / 0.2e1;
t466 = -t273 / 0.2e1;
t443 = t273 / 0.2e1;
t442 = -t274 / 0.2e1;
t441 = t274 / 0.2e1;
t465 = (-t244 * t498 + t245 * t497) * t269 + (t503 * t274 + (t496 + t504) * t273) * t268;
t464 = (-t246 * t498 + t247 * t497) * t269 + ((t496 + t501) * t274 + t502 * t273) * t268;
t462 = t505 * t363;
t460 = t269 ^ 2;
t459 = 2 * qJD(2);
t458 = 2 * qJD(4);
t457 = 4 * qJD(4);
t456 = m(4) / 0.2e1;
t455 = m(5) / 0.2e1;
t454 = m(6) / 0.2e1;
t393 = t269 * t274;
t394 = t269 * t273;
t91 = -t268 * t514 + t269 * t387;
t92 = t268 * t513 - t269 * t386;
t433 = t91 * t393 + t92 * t394;
t296 = -t273 * t386 + t274 * t387;
t45 = t296 * t269 + (-t273 * t513 + t274 * t514) * t268;
t372 = rSges(6,3) * t268 + t269 * t338 + t188;
t46 = (t273 * t372 - t387) * t268;
t47 = (-t274 * t372 + t386) * t268;
t74 = t296 * t268;
t453 = m(6) * (-t269 * t45 + (t273 * t47 + t274 * t46 + t74) * t268 + t433);
t452 = m(6) * (t45 * t74 + t46 * t91 + t47 * t92);
t314 = t146 * t274 - t148 * t273;
t107 = t314 * t268;
t120 = t146 * t269 + t206 * t396;
t121 = -t148 * t269 - t206 * t395;
t80 = t314 * t269 + (-t184 * t274 + t186 * t273) * t268;
t208 = rSges(5,3) * t268 + t269 * t339;
t97 = (t273 * t208 - t146) * t268;
t98 = (-t274 * t208 + t148) * t268;
t451 = m(5) * (t107 * t80 + t120 * t97 + t121 * t98);
t248 = t363 * t268;
t448 = m(6) * (t248 * t74 + t433);
t447 = m(6) * (-t102 * t269 + (t169 * t273 + t170 * t274) * t268);
t446 = t248 / 0.2e1;
t432 = m(6) * qJD(2);
t431 = m(6) * qJD(5);
t397 = t268 * t269;
t388 = t116 * t394 + t118 * t393;
t364 = t363 * t397;
t362 = qJD(4) * t268;
t190 = (t446 + t470) * m(6);
t361 = t190 * qJD(1);
t285 = t268 * t306 + t406;
t49 = -t175 * t244 - t179 * t245 + t273 * t285;
t284 = t268 * t305 + t405;
t50 = -t176 * t244 - t180 * t245 + t273 * t284;
t283 = t268 * t304 + t404;
t51 = -t177 * t244 - t181 * t245 + t273 * t283;
t282 = t268 * t303 + t403;
t52 = -t178 * t244 - t182 * t245 + t273 * t282;
t360 = ((t50 + t52) * t274 + (t49 + t51 - t480) * t273 + t488) * t469 + ((-t474 + t490) * t273 + t484 * t245 - t485 * t244 + t472) * t467;
t53 = -t175 * t246 - t179 * t247 + t274 * t285;
t54 = -t176 * t246 - t180 * t247 + t274 * t284;
t55 = -t177 * t246 - t181 * t247 + t274 * t283;
t56 = -t178 * t246 - t182 * t247 + t274 * t282;
t359 = ((t56 + t54 - t480) * t274 + (t55 + t53) * t273 + t487) * t469 + ((-t474 + t489) * t274 + t484 * t247 - t485 * t246 + t473) * t467;
t358 = (t491 * t274 + t492 * t273 + (t277 * t485 + t279 * t484 - t509) * t269 + t486) * t470 + (t475 + t480) * t468;
t357 = t441 * t504 + t466 * t503;
t356 = t442 * t502 + t443 * t501;
t355 = (t273 * t490 + t472) * t469 + t488 * t468;
t354 = (t274 * t489 + t473) * t469 + t487 * t468;
t353 = t467 * t486 + t470 * t475;
t255 = rSges(4,1) * t268 + rSges(4,2) * t269;
t352 = -t255 - t440;
t256 = rSges(4,1) * t269 - rSges(4,2) * t268;
t351 = -t256 - t439;
t349 = -t258 - t439;
t343 = t363 * t440;
t340 = -t208 + t349;
t261 = t278 * rSges(3,1) + rSges(3,2) * t280;
t281 = -m(5) * (t273 * t97 - t274 * t98) / 0.2e1 - m(6) * (t273 * t46 - t274 * t47) / 0.2e1;
t300 = (-t169 * t274 + t170 * t273) * t454;
t16 = t300 + t281;
t17 = 0.2e1 * (t45 / 0.4e1 - t102 / 0.4e1) * m(6) + 0.2e1 * (t80 / 0.4e1 - t122 / 0.4e1) * m(5);
t319 = -t17 * qJD(1) + t16 * qJD(3);
t310 = t349 - t372;
t309 = -t357 - t360;
t308 = -t356 + t359;
t307 = t342 * t268;
t126 = -t363 * t255 - t343;
t302 = t126 * t256;
t301 = -t257 * t363 - t343;
t295 = t273 * t494 + t250;
t294 = -t274 * t494 + t249;
t210 = t351 * t274;
t209 = t351 * t273;
t189 = m(6) * t446 + t268 * t454;
t168 = t364 - t397;
t167 = t363 * t261;
t130 = t340 * t274;
t128 = t340 * t273;
t125 = -t166 * t269 - t243 * t395;
t124 = t164 * t269 + t243 * t396;
t119 = t310 * t274;
t117 = t310 * t273;
t114 = (t164 * t274 - t166 * t273) * t268;
t109 = -t269 * t374 + t274 * t307;
t108 = t269 * t375 - t273 * t307;
t101 = -t184 * t273 - t186 * t274 + t301;
t99 = (-t273 * t374 + t274 * t375) * t268;
t76 = t447 / 0.2e1;
t75 = t273 * t514 + t274 * t513 + t301;
t73 = -t154 * t269 + (-t277 * t378 + t279 * t382) * t268;
t72 = -t153 * t269 + (-t277 * t379 + t279 * t383) * t268;
t71 = -t152 * t269 + (-t277 * t380 + t279 * t384) * t268;
t70 = -t151 * t269 + (-t277 * t381 + t279 * t385) * t268;
t40 = t248 * t57 + t388;
t30 = t448 / 0.2e1;
t29 = t273 * t56 - t274 * t55;
t28 = t273 * t54 - t274 * t53;
t27 = t273 * t52 - t274 * t51;
t26 = t273 * t50 - t274 * t49;
t18 = (t122 + t80) * t455 + (t102 + t45) * t454;
t15 = t300 - t281;
t6 = t453 / 0.2e1;
t5 = t30 + t76 - t453 / 0.2e1;
t4 = t30 + t6 - t447 / 0.2e1;
t3 = t76 + t6 - t448 / 0.2e1;
t2 = t273 * t356 + t274 * t357 - t471;
t1 = t451 + t452 + (t273 * t355 + t274 * t354 + t358) * t269 + (t273 * t360 + t274 * t359 - t353) * t268;
t7 = [0, t18 * qJD(4) + t189 * qJD(5) + (-m(3) * t167 / 0.2e1 + t126 * t456 + t101 * t455 + t75 * t454) * t459, 0, t18 * qJD(2) + (t114 * t455 + t454 * t99) * t458, t189 * qJD(2); -qJD(4) * t17 + qJD(5) * t190, t2 * qJD(4) + t40 * t431 + (m(6) * (t116 * t117 + t118 * t119 + t57 * t75) + m(5) * (t101 * t78 + t127 * t128 + t129 * t130) + m(4) * (t371 * t126 + (t210 * t352 + t274 * t302) * t274 + (t209 * t352 + t273 * t302) * t273) + m(3) * (-t167 + t261) * t363 * (rSges(3,1) * t280 - t278 * rSges(3,2)) + (t29 + t28 + (t295 * t274 - t462 + (t294 - t483) * t273) * t274 + t482 * t270) * t443 + (t27 + t26 + (t294 * t273 - t462 + (t295 - t482) * t274) * t273 + t483 * t271) * t442) * qJD(2), t16 * qJD(4), t2 * qJD(2) + t5 * qJD(5) + (-t451 / 0.4e1 - t452 / 0.4e1) * t457 + ((t107 * t122 + t114 * t78 + t124 * t129 + t125 * t127 + (-t120 * t274 - t121 * t273) * t243) * t455 + (t102 * t74 + t108 * t118 + t109 * t116 + t169 * t92 + t170 * t91 + t57 * t99) * t454) * t458 + (t273 * t309 - t274 * t308 + t353) * t362 + t319 + (((t70 / 0.2e1 + t72 / 0.2e1 - t354) * t274 + (-t71 / 0.2e1 - t73 / 0.2e1 - t355) * t273 - t358) * t269 + t464 * t443 + t465 * t442) * qJD(4), t361 + t40 * t432 + t5 * qJD(4) + (-t248 * t269 - t168 + t364) * t431; 0, t15 * qJD(4) + ((-t117 * t274 + t119 * t273) * t454 + (-t128 * t274 + t130 * t273) * t455 + (-t209 * t274 + t210 * t273) * t456) * t459, 0, t15 * qJD(2) + ((t124 * t273 - t125 * t274) * t455 + (t108 * t273 - t109 * t274) * t454) * t458, 0; t17 * qJD(2), t1 * qJD(4) + t4 * qJD(5) + ((t116 * t47 + t117 * t92 + t118 * t46 + t119 * t91 + t45 * t57 + t74 * t75) * t454 + (t101 * t107 + t120 * t130 + t121 * t128 + t127 * t98 + t129 * t97 + t78 * t80) * t455) * t459 - t319 + (t309 * t274 + t308 * t273 + (t491 * t466 + (t516 * t273 - t490 * t274) * t443 + (t489 * t273 - t515 * t274 + t492) * t441) * t269 + ((-t93 / 0.2e1 - t95 / 0.2e1 + t29 / 0.2e1 + t28 / 0.2e1) * t274 + (t94 / 0.2e1 + t96 / 0.2e1 + t26 / 0.2e1 + t27 / 0.2e1) * t273) * t268 + t471) * qJD(2), -t16 * qJD(2), t1 * qJD(2) + (m(6) * (t108 * t91 + t109 * t92 + t74 * t99) / 0.4e1 + m(5) * (t107 * t114 + t120 * t124 + t121 * t125) / 0.4e1) * t457 + (-t237 / 0.2e1 - t236 / 0.2e1) * qJD(4) * t269 * t460 + (((t73 + t71) * t274 + (t72 + t70) * t273 + (t277 * t498 + t279 * t497) * t269) * t468 + t465 * t443 + t464 * t441) * t362, t4 * qJD(2); -t190 * qJD(2), -t361 + (-t269 * t75 + (t117 * t273 + t119 * t274 + t57) * t268 - t40 + t388) * t432 + t3 * qJD(4) + t168 * t431, 0, t3 * qJD(2) + m(6) * (-t269 * t99 + (t108 * t274 + t109 * t273) * t268) * qJD(4), t168 * t432;];
Cq = t7;
