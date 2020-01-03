% Calculate time derivative of joint inertia matrix for
% S5RRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP6_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP6_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP6_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:56:47
% EndTime: 2019-12-31 19:57:12
% DurationCPUTime: 15.48s
% Computational Cost: add. (19212->830), mult. (28161->1152), div. (0->0), fcn. (26514->8), ass. (0->406)
t277 = qJ(2) + pkin(8);
t270 = sin(t277);
t271 = cos(t277);
t282 = sin(qJ(4));
t285 = cos(qJ(4));
t331 = Icges(6,5) * t285 - Icges(6,6) * t282;
t184 = -Icges(6,3) * t271 + t270 * t331;
t332 = Icges(5,5) * t285 - Icges(5,6) * t282;
t185 = -Icges(5,3) * t271 + t270 * t332;
t524 = t184 + t185;
t457 = Icges(6,4) * t285;
t335 = -Icges(6,2) * t282 + t457;
t186 = -Icges(6,6) * t271 + t270 * t335;
t459 = Icges(5,4) * t285;
t336 = -Icges(5,2) * t282 + t459;
t187 = -Icges(5,6) * t271 + t270 * t336;
t520 = t187 + t186;
t458 = Icges(6,4) * t282;
t341 = Icges(6,1) * t285 - t458;
t188 = -Icges(6,5) * t271 + t270 * t341;
t460 = Icges(5,4) * t282;
t342 = Icges(5,1) * t285 - t460;
t189 = -Icges(5,5) * t271 + t270 * t342;
t519 = -t189 - t188;
t280 = -qJ(5) - pkin(7);
t523 = rSges(6,3) - t280;
t522 = t524 * t271 + (t520 * t282 + t519 * t285) * t270;
t406 = qJD(4) * t270;
t124 = (-Icges(6,5) * t282 - Icges(6,6) * t285) * t406 + (Icges(6,3) * t270 + t271 * t331) * qJD(2);
t125 = (-Icges(5,5) * t282 - Icges(5,6) * t285) * t406 + (Icges(5,3) * t270 + t271 * t332) * qJD(2);
t521 = -t125 - t124;
t126 = (-Icges(6,2) * t285 - t458) * t406 + (Icges(6,6) * t270 + t271 * t335) * qJD(2);
t127 = (-Icges(5,2) * t285 - t460) * t406 + (Icges(5,6) * t270 + t271 * t336) * qJD(2);
t518 = (-t126 - t127) * t282;
t284 = sin(qJ(1));
t287 = cos(qJ(1));
t431 = t284 * t285;
t433 = t282 * t287;
t224 = -t271 * t433 + t431;
t430 = t285 * t287;
t432 = t284 * t282;
t225 = t271 * t430 + t432;
t438 = t270 * t287;
t147 = Icges(6,5) * t225 + Icges(6,6) * t224 + Icges(6,3) * t438;
t151 = Icges(6,4) * t225 + Icges(6,2) * t224 + Icges(6,6) * t438;
t155 = Icges(6,1) * t225 + Icges(6,4) * t224 + Icges(6,5) * t438;
t326 = -t151 * t282 + t155 * t285;
t63 = -t147 * t271 + t270 * t326;
t149 = Icges(5,5) * t225 + Icges(5,6) * t224 + Icges(5,3) * t438;
t153 = Icges(5,4) * t225 + Icges(5,2) * t224 + Icges(5,6) * t438;
t157 = Icges(5,1) * t225 + Icges(5,4) * t224 + Icges(5,5) * t438;
t324 = -t153 * t282 + t157 * t285;
t65 = -t149 * t271 + t270 * t324;
t474 = t63 + t65;
t222 = -t271 * t432 - t430;
t223 = t271 * t431 - t433;
t440 = t270 * t284;
t146 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t440;
t150 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t440;
t154 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t440;
t327 = -t150 * t282 + t154 * t285;
t62 = -t146 * t271 + t270 * t327;
t148 = Icges(5,5) * t223 + Icges(5,6) * t222 + Icges(5,3) * t440;
t152 = Icges(5,4) * t223 + Icges(5,2) * t222 + Icges(5,6) * t440;
t156 = Icges(5,1) * t223 + Icges(5,4) * t222 + Icges(5,5) * t440;
t325 = -t152 * t282 + t156 * t285;
t64 = -t148 * t271 + t270 * t325;
t475 = t62 + t64;
t517 = t284 * t475 + t287 * t474;
t516 = t284 * t474 - t287 * t475;
t368 = qJD(4) * t271 - qJD(1);
t314 = t285 * t368;
t367 = qJD(1) * t271 - qJD(4);
t410 = qJD(2) * t284;
t384 = t270 * t410;
t132 = -t284 * t314 + (-t287 * t367 + t384) * t282;
t313 = t368 * t282;
t409 = qJD(2) * t285;
t133 = t367 * t430 + (-t270 * t409 - t313) * t284;
t381 = t271 * t410;
t415 = qJD(1) * t287;
t294 = t270 * t415 + t381;
t71 = Icges(6,5) * t133 + Icges(6,6) * t132 + Icges(6,3) * t294;
t75 = Icges(6,4) * t133 + Icges(6,2) * t132 + Icges(6,6) * t294;
t79 = Icges(6,1) * t133 + Icges(6,4) * t132 + Icges(6,5) * t294;
t19 = (qJD(2) * t327 - t71) * t271 + (qJD(2) * t146 - t282 * t75 + t285 * t79 + (-t150 * t285 - t154 * t282) * qJD(4)) * t270;
t73 = Icges(5,5) * t133 + Icges(5,6) * t132 + Icges(5,3) * t294;
t77 = Icges(5,4) * t133 + Icges(5,2) * t132 + Icges(5,6) * t294;
t81 = Icges(5,1) * t133 + Icges(5,4) * t132 + Icges(5,5) * t294;
t21 = (qJD(2) * t325 - t73) * t271 + (qJD(2) * t148 - t282 * t77 + t285 * t81 + (-t152 * t285 - t156 * t282) * qJD(4)) * t270;
t515 = -t19 - t21;
t407 = qJD(2) * t287;
t383 = t270 * t407;
t501 = t284 * t367 + t383;
t130 = t282 * t501 - t287 * t314;
t131 = -t285 * t501 - t287 * t313;
t379 = t271 * t407;
t416 = qJD(1) * t284;
t387 = t270 * t416;
t293 = t379 - t387;
t70 = Icges(6,5) * t131 + Icges(6,6) * t130 + Icges(6,3) * t293;
t74 = Icges(6,4) * t131 + Icges(6,2) * t130 + Icges(6,6) * t293;
t78 = Icges(6,1) * t131 + Icges(6,4) * t130 + Icges(6,5) * t293;
t20 = (qJD(2) * t326 - t70) * t271 + (qJD(2) * t147 - t282 * t74 + t285 * t78 + (-t151 * t285 - t155 * t282) * qJD(4)) * t270;
t72 = Icges(5,5) * t131 + Icges(5,6) * t130 + Icges(5,3) * t293;
t76 = Icges(5,4) * t131 + Icges(5,2) * t130 + Icges(5,6) * t293;
t80 = Icges(5,1) * t131 + Icges(5,4) * t130 + Icges(5,5) * t293;
t22 = (qJD(2) * t324 - t72) * t271 + (qJD(2) * t149 - t282 * t76 + t285 * t80 + (-t153 * t285 - t157 * t282) * qJD(4)) * t270;
t514 = t20 + t22;
t54 = t148 * t440 + t152 * t222 + t156 * t223;
t55 = t149 * t440 + t153 * t222 + t157 * t223;
t351 = t284 * t54 + t287 * t55;
t52 = t146 * t440 + t150 * t222 + t154 * t223;
t53 = t147 * t440 + t151 * t222 + t155 * t223;
t352 = t284 * t52 + t287 * t53;
t86 = t184 * t440 + t186 * t222 + t188 * t223;
t87 = t185 * t440 + t187 * t222 + t189 * t223;
t478 = (-t86 - t87) * t271 + (t351 + t352) * t270;
t58 = t148 * t438 + t224 * t152 + t225 * t156;
t59 = t149 * t438 + t224 * t153 + t225 * t157;
t349 = t284 * t58 + t287 * t59;
t56 = t146 * t438 + t224 * t150 + t225 * t154;
t57 = t147 * t438 + t224 * t151 + t225 * t155;
t350 = t284 * t56 + t287 * t57;
t88 = t184 * t438 + t224 * t186 + t225 * t188;
t89 = t185 * t438 + t224 * t187 + t225 * t189;
t513 = (-t88 - t89) * t271 + (t349 + t350) * t270;
t479 = pkin(7) + t280;
t512 = t271 * t479;
t283 = sin(qJ(2));
t286 = cos(qJ(2));
t463 = Icges(3,4) * t286;
t340 = -Icges(3,2) * t283 + t463;
t210 = Icges(3,6) * t284 + t287 * t340;
t464 = Icges(3,4) * t283;
t346 = Icges(3,1) * t286 - t464;
t212 = Icges(3,5) * t284 + t287 * t346;
t316 = t210 * t283 - t212 * t286;
t511 = t284 * t316;
t461 = Icges(4,4) * t271;
t338 = -Icges(4,2) * t270 + t461;
t199 = Icges(4,6) * t284 + t287 * t338;
t462 = Icges(4,4) * t270;
t344 = Icges(4,1) * t271 - t462;
t201 = Icges(4,5) * t284 + t287 * t344;
t318 = t199 * t270 - t201 * t271;
t510 = t284 * t318;
t267 = pkin(2) * t286 + pkin(1);
t482 = pkin(1) - t267;
t509 = t284 * t482;
t209 = -Icges(3,6) * t287 + t284 * t340;
t211 = -Icges(3,5) * t287 + t284 * t346;
t317 = t209 * t283 - t211 * t286;
t508 = t287 * t317;
t198 = -Icges(4,6) * t287 + t284 * t338;
t200 = -Icges(4,5) * t287 + t284 * t344;
t319 = t198 * t270 - t200 * t271;
t507 = t287 * t319;
t266 = pkin(4) * t285 + pkin(3);
t437 = t271 * t287;
t506 = t225 * rSges(6,1) + t224 * rSges(6,2) + pkin(4) * t432 + t266 * t437 + t523 * t438;
t274 = t284 * rSges(4,3);
t505 = -rSges(4,2) * t438 + t274;
t404 = qJD(4) * t285;
t396 = pkin(4) * t404;
t402 = pkin(4) * t433;
t403 = qJD(5) * t270;
t504 = t131 * rSges(6,1) + t130 * rSges(6,2) + rSges(6,3) * t379 + qJD(1) * t402 + t280 * t387 + t284 * t396 + t287 * t403;
t128 = (-Icges(6,1) * t282 - t457) * t406 + (Icges(6,5) * t270 + t271 * t341) * qJD(2);
t129 = (-Icges(5,1) * t282 - t459) * t406 + (Icges(5,5) * t270 + t271 * t342) * qJD(2);
t414 = qJD(2) * t270;
t502 = (t128 + t129) * t270 * t285 + t524 * t414 - t519 * t271 * t409;
t333 = Icges(4,5) * t271 - Icges(4,6) * t270;
t196 = -Icges(4,3) * t287 + t284 * t333;
t334 = Icges(3,5) * t286 - Icges(3,6) * t283;
t207 = -Icges(3,3) * t287 + t284 * t334;
t256 = pkin(3) * t437;
t220 = pkin(7) * t438 + t256;
t427 = -t220 + t506;
t481 = pkin(3) - t266;
t291 = -t270 * t479 - t271 * t481;
t354 = -t223 * rSges(6,1) - t222 * rSges(6,2);
t428 = rSges(6,3) * t440 + t284 * t291 - t354 - t402;
t500 = -t284 * t428 - t287 * t427;
t499 = 2 * m(3);
t498 = 2 * m(4);
t497 = 2 * m(5);
t496 = 2 * m(6);
t278 = t284 ^ 2;
t279 = t287 ^ 2;
t494 = t284 / 0.2e1;
t491 = -rSges(5,3) - pkin(7);
t251 = rSges(3,1) * t283 + rSges(3,2) * t286;
t490 = m(3) * t251;
t489 = pkin(2) * t283;
t488 = pkin(3) * t270;
t487 = pkin(3) * t271;
t486 = pkin(4) * t282;
t485 = t284 * pkin(6);
t276 = t287 * pkin(6);
t281 = -qJ(3) - pkin(6);
t480 = -pkin(6) - t281;
t412 = qJD(2) * t282;
t476 = t502 + (-t412 * t520 + t521) * t271 + ((t282 * t519 - t285 * t520) * qJD(4) + t518) * t270;
t246 = pkin(7) * t379;
t405 = qJD(4) * t282;
t397 = pkin(4) * t405;
t473 = -t246 + (pkin(7) * t416 + t407 * t481) * t270 + ((-qJD(2) * t280 - t397) * t287 + t481 * t416) * t271 - rSges(6,3) * t387 + t504;
t245 = pkin(3) * t384;
t355 = t133 * rSges(6,1) + t132 * rSges(6,2);
t441 = t266 * t270;
t472 = t245 + (qJD(1) * t291 - t396) * t287 + (t403 - pkin(4) * t313 + (-t441 - t512) * qJD(2)) * t284 + rSges(6,3) * t294 + t355;
t471 = t522 * t414;
t470 = rSges(3,1) * t286;
t469 = rSges(4,1) * t271;
t468 = rSges(3,2) * t283;
t467 = rSges(3,3) * t287;
t466 = rSges(6,3) * t270;
t275 = t284 * rSges(3,3);
t357 = -t223 * rSges(5,1) - t222 * rSges(5,2);
t159 = rSges(5,3) * t440 - t357;
t450 = t159 * t287;
t449 = t198 * t271;
t448 = t199 * t271;
t447 = t200 * t270;
t446 = t201 * t270;
t445 = t209 * t286;
t444 = t210 * t286;
t443 = t211 * t283;
t442 = t212 * t283;
t436 = t281 * t287;
t353 = rSges(6,1) * t285 - rSges(6,2) * t282;
t378 = t270 * t405;
t429 = -pkin(4) * t378 - qJD(5) * t271 + (-rSges(6,1) * t282 - rSges(6,2) * t285) * t406 + (t271 * t353 + t291 + t466) * qJD(2);
t426 = -rSges(6,3) * t271 + t512 + (t353 - t481) * t270;
t194 = t276 + t436 - t509;
t257 = t287 * t267;
t195 = -pkin(1) * t287 + t284 * t480 + t257;
t425 = t284 * t194 + t287 * t195;
t424 = -t195 - t220;
t236 = -pkin(7) * t271 + t488;
t386 = t283 * t416;
t261 = pkin(2) * t386;
t423 = t236 * t416 + t261;
t422 = rSges(4,2) * t387 + rSges(4,3) * t415;
t421 = t287 * t470 + t275;
t420 = t278 + t279;
t197 = Icges(4,3) * t284 + t287 * t333;
t419 = qJD(1) * t197;
t208 = Icges(3,3) * t284 + t287 * t334;
t418 = qJD(1) * t208;
t417 = qJD(1) * t281;
t413 = qJD(2) * t271;
t411 = qJD(2) * t283;
t408 = qJD(2) * t286;
t401 = t287 * t468;
t399 = pkin(2) * t411;
t398 = pkin(2) * t408;
t35 = t53 * t284 - t287 * t52;
t36 = t55 * t284 - t287 * t54;
t395 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t57 * t284 - t287 * t56;
t38 = t59 * t284 - t287 * t58;
t394 = t37 / 0.2e1 + t38 / 0.2e1;
t390 = t131 * rSges(5,1) + t130 * rSges(5,2) + rSges(5,3) * t379;
t272 = qJD(3) * t284;
t388 = qJD(3) * t287 + t281 * t416 + t284 * t399;
t389 = t284 * ((-t287 * t482 - t485) * qJD(1) - t388) + t287 * (-t287 * t399 + t272 + (t287 * t480 + t509) * qJD(1)) + t194 * t415;
t161 = t225 * rSges(5,1) + t224 * rSges(5,2) + rSges(5,3) * t438;
t356 = rSges(5,1) * t285 - rSges(5,2) * t282;
t191 = -rSges(5,3) * t271 + t270 * t356;
t385 = t191 * t416;
t234 = rSges(4,1) * t270 + rSges(4,2) * t271;
t376 = -t234 - t489;
t375 = -t236 - t489;
t374 = t428 * t287;
t373 = t426 * t287;
t372 = t426 * t284;
t371 = -t284 * t281 + t257;
t370 = -t266 * t271 - t267;
t369 = qJD(1) * t426;
t362 = pkin(7) * t270 + t487;
t219 = t362 * t284;
t366 = t284 * t219 + t287 * t220 + t425;
t364 = -t191 + t375;
t363 = -t362 * qJD(2) - t398;
t361 = t284 * t369;
t360 = -t468 + t470;
t359 = -rSges(4,2) * t270 + t469;
t358 = t133 * rSges(5,1) + t132 * rSges(5,2);
t60 = t270 * t372 + t271 * t428;
t61 = -t270 * t373 - t271 * t427;
t348 = t284 * t61 + t287 * t60;
t292 = -t270 * t523 + t370;
t95 = (-t281 + t486) * t287 + t292 * t284 + t354;
t96 = t371 + t506;
t347 = t284 * t96 + t287 * t95;
t345 = Icges(3,1) * t283 + t463;
t343 = Icges(4,1) * t270 + t461;
t339 = Icges(3,2) * t286 + t464;
t337 = Icges(4,2) * t271 + t462;
t315 = t375 - t426;
t101 = t315 * t284;
t102 = t315 * t287;
t330 = t101 * t284 + t102 * t287;
t323 = -t161 * t284 + t450;
t322 = -t284 * t159 - t161 * t287;
t203 = rSges(4,1) * t437 + t505;
t135 = (-rSges(5,1) * t282 - rSges(5,2) * t285) * t406 + (rSges(5,3) * t270 + t271 * t356) * qJD(2);
t312 = -t135 + t363;
t311 = -pkin(1) - t360;
t309 = t62 / 0.2e1 + t64 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1;
t308 = t88 / 0.2e1 + t89 / 0.2e1 + t63 / 0.2e1 + t65 / 0.2e1;
t123 = t364 * t287;
t307 = -t267 - t359;
t305 = t284 * (pkin(7) * t294 + qJD(1) * t256 - t245) + t287 * (-pkin(7) * t387 + t246 + (-t271 * t416 - t383) * pkin(3)) + t219 * t415 + t389;
t304 = t363 - t429;
t303 = qJD(2) * t251;
t302 = qJD(2) * t234;
t301 = qJD(2) * t345;
t300 = qJD(2) * t343;
t299 = qJD(2) * t339;
t298 = qJD(2) * t337;
t297 = qJD(2) * (-Icges(3,5) * t283 - Icges(3,6) * t286);
t296 = qJD(2) * (-Icges(4,5) * t270 - Icges(4,6) * t271);
t295 = t270 * t491 - t267 - t487;
t290 = -t284 * t427 + t374;
t289 = rSges(3,2) * t386 + rSges(3,3) * t415 - t287 * t303;
t288 = t284 * t295 - t436;
t241 = t360 * qJD(2);
t229 = t359 * qJD(2);
t218 = -t401 + t421;
t217 = t284 * t360 - t467;
t202 = -rSges(4,3) * t287 + t284 * t359;
t193 = t376 * t287;
t192 = t376 * t284;
t183 = t485 + (pkin(1) - t468) * t287 + t421;
t182 = t284 * t311 + t276 + t467;
t173 = t203 + t371;
t172 = (rSges(4,3) - t281) * t287 + t307 * t284;
t167 = t284 * t297 + t418;
t166 = -qJD(1) * t207 + t287 * t297;
t141 = t284 * t296 + t419;
t140 = -qJD(1) * t196 + t287 * t296;
t137 = t251 * t410 + ((-rSges(3,3) - pkin(6)) * t284 + t311 * t287) * qJD(1);
t136 = (t276 + (-pkin(1) - t470) * t284) * qJD(1) + t289;
t122 = t364 * t284;
t117 = -t234 * t415 - t284 * t229 + (-t283 * t415 - t284 * t408) * pkin(2);
t116 = t234 * t416 + t261 + (-t229 - t398) * t287;
t112 = t284 * t208 - t316 * t287;
t111 = t284 * t207 - t508;
t110 = -t208 * t287 - t511;
t109 = -t207 * t287 - t284 * t317;
t108 = t234 * t410 + (t287 * t307 - t274) * qJD(1) + t388;
t107 = t272 + (-t267 - t469) * t416 + (qJD(2) * t376 - t417) * t287 + t422;
t106 = t371 + t161 + t220;
t105 = t288 + t357;
t104 = -t271 * t161 - t191 * t438;
t103 = t159 * t271 + t191 * t440;
t100 = t284 * t197 - t318 * t287;
t99 = t284 * t196 - t507;
t98 = -t197 * t287 - t510;
t97 = -t196 * t287 - t284 * t319;
t90 = t323 * t270;
t85 = rSges(5,3) * t294 + t358;
t83 = -rSges(5,3) * t387 + t390;
t67 = qJD(1) * t123 + t284 * t312;
t66 = t287 * t312 + t385 + t423;
t51 = t295 * t415 + t381 * t491 + t245 - t358 + t388;
t50 = t246 + t272 + (-t488 - t489) * t407 + t288 * qJD(1) + t390;
t49 = -t322 + t366;
t48 = t290 * t270;
t47 = qJD(1) * t102 + t284 * t304;
t46 = t287 * t304 + t361 + t423;
t45 = (qJD(1) * t292 + t396) * t287 + (-t403 + t368 * t486 + (-t271 * t523 + t441) * qJD(2)) * t284 - t355 + t388;
t44 = t272 + (t370 - t466) * t416 + (-t271 * t397 - t417 + (-t271 * t280 - t441 - t489) * qJD(2)) * t287 + t504;
t43 = t366 - t500;
t42 = (t191 * t410 + t85) * t271 + (-qJD(2) * t159 + t284 * t135 + t191 * t415) * t270;
t41 = (-t191 * t407 - t83) * t271 + (qJD(2) * t161 - t135 * t287 + t385) * t270;
t34 = t125 * t440 + t222 * t127 + t223 * t129 + t132 * t187 + t133 * t189 + t185 * t294;
t33 = t124 * t440 + t222 * t126 + t223 * t128 + t132 * t186 + t133 * t188 + t184 * t294;
t32 = t125 * t438 + t224 * t127 + t225 * t129 + t130 * t187 + t131 * t189 + t185 * t293;
t31 = t124 * t438 + t224 * t126 + t225 * t128 + t130 * t186 + t131 * t188 + t184 * t293;
t30 = t323 * t413 + (qJD(1) * t322 - t284 * t83 + t287 * t85) * t270;
t25 = t284 * t85 + t287 * t83 + (t450 + (-t161 + t424) * t284) * qJD(1) + t305;
t24 = (qJD(2) * t372 + t472) * t271 + (-qJD(2) * t428 + t284 * t429 + t287 * t369) * t270;
t23 = (-qJD(2) * t373 - t473) * t271 + (qJD(2) * t427 - t287 * t429 + t361) * t270;
t18 = t132 * t153 + t133 * t157 + t149 * t294 + t222 * t76 + t223 * t80 + t440 * t72;
t17 = t132 * t152 + t133 * t156 + t148 * t294 + t222 * t77 + t223 * t81 + t440 * t73;
t16 = t132 * t151 + t133 * t155 + t147 * t294 + t222 * t74 + t223 * t78 + t440 * t70;
t15 = t132 * t150 + t133 * t154 + t146 * t294 + t222 * t75 + t223 * t79 + t440 * t71;
t14 = t130 * t153 + t131 * t157 + t149 * t293 + t224 * t76 + t225 * t80 + t438 * t72;
t13 = t130 * t152 + t131 * t156 + t148 * t293 + t224 * t77 + t225 * t81 + t438 * t73;
t12 = t130 * t151 + t131 * t155 + t147 * t293 + t224 * t74 + t225 * t78 + t438 * t70;
t11 = t130 * t150 + t131 * t154 + t146 * t293 + t224 * t75 + t225 * t79 + t438 * t71;
t10 = t290 * t413 + (qJD(1) * t500 - t473 * t284 + t472 * t287) * t270;
t9 = t473 * t287 + t472 * t284 + (t374 + (t424 - t427) * t284) * qJD(1) + t305;
t8 = qJD(1) * t351 - t17 * t287 + t18 * t284;
t7 = qJD(1) * t352 - t15 * t287 + t16 * t284;
t6 = qJD(1) * t349 - t13 * t287 + t14 * t284;
t5 = qJD(1) * t350 - t11 * t287 + t12 * t284;
t4 = (qJD(2) * t351 - t34) * t271 + (-qJD(1) * t36 + qJD(2) * t87 + t17 * t284 + t18 * t287) * t270;
t3 = (qJD(2) * t352 - t33) * t271 + (-qJD(1) * t35 + qJD(2) * t86 + t15 * t284 + t16 * t287) * t270;
t2 = (qJD(2) * t349 - t32) * t271 + (-qJD(1) * t38 + qJD(2) * t89 + t13 * t284 + t14 * t287) * t270;
t1 = (qJD(2) * t350 - t31) * t271 + (-qJD(1) * t37 + qJD(2) * t88 + t11 * t284 + t12 * t287) * t270;
t26 = [(t136 * t183 + t137 * t182) * t499 + (t107 * t173 + t108 * t172) * t498 + (t105 * t51 + t106 * t50) * t497 + (t44 * t96 + t45 * t95) * t496 + (-t337 + t344) * t414 + (t338 + t343) * t413 + (t346 - t339) * t411 + (t345 + t340) * t408 + t519 * t378 + t521 * t271 + t518 * t270 + t502 + t520 * (-t270 * t404 - t271 * t412); m(3) * ((-t136 * t284 - t137 * t287) * t251 + (-t182 * t287 - t183 * t284) * t241) + m(4) * (t107 * t192 + t108 * t193 + t116 * t172 + t117 * t173) + m(5) * (t105 * t66 + t106 * t67 + t122 * t50 + t123 * t51) + m(6) * (t101 * t44 + t102 * t45 + t46 * t95 + t47 * t96) + ((t444 / 0.2e1 + t442 / 0.2e1 - t183 * t490 + t448 / 0.2e1 + t446 / 0.2e1 + t308) * t287 + (t182 * t490 + t445 / 0.2e1 + t443 / 0.2e1 + t449 / 0.2e1 + t447 / 0.2e1 + t309) * t284) * qJD(1) + (t334 + t333) * qJD(2) * (t279 / 0.2e1 + t278 / 0.2e1) + ((-qJD(1) * t198 - t287 * t298) * t271 + (-qJD(1) * t200 - t287 * t300) * t270 + (-qJD(1) * t209 - t287 * t299) * t286 + (-qJD(1) * t211 - t287 * t301) * t283 + t31 + t32 + (-t316 - t318) * qJD(2) + t514) * t494 - ((qJD(1) * t199 - t284 * t298) * t271 + (qJD(1) * t201 - t284 * t300) * t270 + (qJD(1) * t210 - t284 * t299) * t286 + (qJD(1) * t212 - t284 * t301) * t283 + t33 + t34 + (-t317 - t319) * qJD(2) - t515) * t287 / 0.2e1; (t122 * t67 + t123 * t66 + t25 * t49) * t497 + t284 * t5 + (t101 * t47 + t102 * t46 + t43 * t9) * t496 + t284 * t6 - t287 * t8 - t287 * t7 + (t193 * t116 + t192 * t117 + (t284 * t202 + t203 * t287 + t425) * ((qJD(1) * t202 - t287 * t302 + t422) * t287 + (-t284 * t302 + (-t195 - t203 + t505) * qJD(1)) * t284 + t389)) * t498 - t287 * ((t287 * t141 + (t98 + t507) * qJD(1)) * t287 + (t97 * qJD(1) + (-t199 * t413 - t201 * t414 + t419) * t284 + (-t140 + (t447 + t449) * qJD(2) - t318 * qJD(1)) * t287) * t284) + t284 * ((t284 * t166 + (t111 + t511) * qJD(1)) * t284 + (t112 * qJD(1) + (t209 * t408 + t211 * t411) * t287 + (-t167 + (-t442 - t444) * qJD(2) + (t208 - t317) * qJD(1)) * t284) * t287) - t287 * ((t287 * t167 + (t110 + t508) * qJD(1)) * t287 + (t109 * qJD(1) + (-t210 * t408 - t212 * t411 + t418) * t284 + (-t166 + (t443 + t445) * qJD(2) - t316 * qJD(1)) * t287) * t284) + t284 * ((t284 * t140 + (t99 + t510) * qJD(1)) * t284 + (t100 * qJD(1) + (t198 * t413 + t200 * t414) * t287 + (-t141 + (-t446 - t448) * qJD(2) + (t197 - t319) * qJD(1)) * t284) * t287) + ((t284 * t217 + t218 * t287) * ((qJD(1) * t217 + t289) * t287 + (-t284 * t303 + (-t218 - t401 + t275) * qJD(1)) * t284) + t420 * t251 * t241) * t499 + (t35 + t36 + (-t109 - t97) * t287 + (t110 + t98) * t284) * t416 + (t37 + t38 + (-t111 - t99) * t287 + (t100 + t112) * t284) * t415; m(6) * (qJD(1) * t347 + t284 * t45 - t287 * t44) + m(5) * (t284 * t51 - t287 * t50 + (t105 * t287 + t106 * t284) * qJD(1)) + m(4) * (-t107 * t287 + t284 * t108 + (t172 * t287 + t173 * t284) * qJD(1)); m(5) * (t284 * t66 - t287 * t67 + (t122 * t284 + t123 * t287) * qJD(1)) + m(6) * (qJD(1) * t330 + t284 * t46 - t287 * t47) + m(4) * (t284 * t116 - t117 * t287 + (t192 * t284 + t193 * t287) * qJD(1)); 0; m(5) * (t103 * t51 + t104 * t50 + t105 * t42 + t106 * t41) + m(6) * (t23 * t96 + t24 * t95 + t44 * t61 + t45 * t60) + ((t284 * t309 + t287 * t308) * qJD(2) - t476) * t271 + ((t20 / 0.2e1 + t22 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1) * t287 + (t33 / 0.2e1 + t34 / 0.2e1 + t19 / 0.2e1 + t21 / 0.2e1) * t284 + (-t284 * t308 + t287 * t309) * qJD(1)) * t270 - t471; m(5) * (t103 * t66 + t104 * t67 + t122 * t41 + t123 * t42 + t25 * t90 + t30 * t49) + m(6) * (t10 * t43 + t101 * t23 + t102 * t24 + t46 * t60 + t47 * t61 + t48 * t9) + (-t3 / 0.2e1 - t4 / 0.2e1 + t394 * t413) * t287 + (t1 / 0.2e1 + t2 / 0.2e1 + t395 * t413) * t284 + ((-t284 * t394 + t287 * t395) * qJD(1) + (t7 + t8) * t494 + (t5 + t6) * t287 / 0.2e1 + t516 * qJD(2) / 0.2e1) * t270 - (qJD(1) * t517 + t514 * t284 + t515 * t287) * t271 / 0.2e1 + (t478 * t284 + t287 * t513) * qJD(1) / 0.2e1; m(5) * (t42 * t284 - t287 * t41 + (t103 * t287 + t104 * t284) * qJD(1)) + m(6) * (qJD(1) * t348 - t23 * t287 + t24 * t284); (t10 * t48 + t23 * t61 + t24 * t60) * t496 + (t103 * t42 + t104 * t41 + t30 * t90) * t497 + (t476 * t271 + ((-t271 * t474 + t513) * t287 + (-t271 * t475 + t478) * t284) * qJD(2) + t471) * t271 + ((t1 + t2) * t287 + (t3 + t4) * t284 + (t284 * t515 - t514 * t287) * t271 + (t517 * t270 + t522 * t271) * qJD(2) + (t271 * t516 - t284 * t513 + t478 * t287) * qJD(1)) * t270; m(6) * (t347 * t413 + (t284 * t44 + t287 * t45 + (-t284 * t95 + t287 * t96) * qJD(1)) * t270); m(6) * ((qJD(2) * t330 - t9) * t271 + (qJD(2) * t43 + t284 * t47 + t287 * t46 + (t101 * t287 - t102 * t284) * qJD(1)) * t270); 0; m(6) * ((qJD(2) * t348 - t10) * t271 + (qJD(2) * t48 + t23 * t284 + t24 * t287 + (-t284 * t60 + t287 * t61) * qJD(1)) * t270); (-0.1e1 + t420) * t270 * t413 * t496;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t26(1), t26(2), t26(4), t26(7), t26(11); t26(2), t26(3), t26(5), t26(8), t26(12); t26(4), t26(5), t26(6), t26(9), t26(13); t26(7), t26(8), t26(9), t26(10), t26(14); t26(11), t26(12), t26(13), t26(14), t26(15);];
Mq = res;
