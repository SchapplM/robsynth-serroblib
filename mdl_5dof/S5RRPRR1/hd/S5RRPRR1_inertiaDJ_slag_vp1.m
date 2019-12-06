% Calculate time derivative of joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-12-05 18:26
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:23:54
% EndTime: 2019-12-05 18:24:29
% DurationCPUTime: 16.35s
% Computational Cost: add. (16870->704), mult. (24544->1000), div. (0->0), fcn. (22905->8), ass. (0->374)
t296 = cos(qJ(2));
t293 = sin(qJ(2));
t477 = Icges(4,4) * t293;
t479 = Icges(3,4) * t293;
t519 = -t477 - t479 + (-Icges(3,2) - Icges(4,2)) * t296;
t476 = Icges(4,4) * t296;
t478 = Icges(3,4) * t296;
t518 = -t476 - t478 + (-Icges(3,1) - Icges(4,1)) * t293;
t517 = (t518 * t293 + t519 * t296) * qJD(2);
t294 = sin(qJ(1));
t288 = t294 ^ 2;
t297 = cos(qJ(1));
t289 = t297 ^ 2;
t427 = t288 + t289;
t290 = qJ(2) + qJ(4);
t280 = sin(t290);
t287 = qJD(2) + qJD(4);
t295 = cos(qJ(5));
t292 = sin(qJ(5));
t472 = Icges(6,4) * t295;
t350 = -Icges(6,2) * t292 + t472;
t281 = cos(t290);
t448 = t281 * t287;
t473 = Icges(6,4) * t292;
t116 = t350 * t448 + (Icges(6,6) * t287 + (-Icges(6,2) * t295 - t473) * qJD(5)) * t280;
t185 = -Icges(6,6) * t281 + t280 * t350;
t356 = Icges(6,1) * t295 - t473;
t186 = -Icges(6,5) * t281 + t280 * t356;
t516 = -t116 * t292 + (-t185 * t295 - t186 * t292) * qJD(5);
t474 = Icges(5,4) * t281;
t351 = -Icges(5,2) * t280 + t474;
t191 = Icges(5,6) * t294 + t297 * t351;
t475 = Icges(5,4) * t280;
t357 = Icges(5,1) * t281 - t475;
t193 = Icges(5,5) * t294 + t297 * t357;
t336 = t191 * t280 - t193 * t281;
t515 = t294 * t336;
t190 = -Icges(5,6) * t297 + t294 * t351;
t192 = -Icges(5,5) * t297 + t294 * t357;
t337 = t190 * t280 - t192 * t281;
t514 = t297 * t337;
t355 = -Icges(3,2) * t293 + t478;
t209 = Icges(3,6) * t294 + t297 * t355;
t361 = Icges(3,1) * t296 - t479;
t213 = Icges(3,5) * t294 + t297 * t361;
t332 = t209 * t293 - t213 * t296;
t315 = t332 * t294;
t208 = -Icges(3,6) * t297 + t294 * t355;
t212 = -Icges(3,5) * t297 + t294 * t361;
t333 = t208 * t293 - t212 * t296;
t316 = t333 * t297;
t353 = -Icges(4,2) * t293 + t476;
t207 = Icges(4,6) * t294 + t297 * t353;
t359 = Icges(4,1) * t296 - t477;
t211 = Icges(4,5) * t294 + t297 * t359;
t334 = t207 * t293 - t211 * t296;
t317 = t334 * t294;
t206 = -Icges(4,6) * t297 + t294 * t353;
t210 = -Icges(4,5) * t297 + t294 * t359;
t335 = t206 * t293 - t210 * t296;
t318 = t335 * t297;
t346 = Icges(6,5) * t295 - Icges(6,6) * t292;
t115 = t346 * t448 + (Icges(6,3) * t287 + (-Icges(6,5) * t292 - Icges(6,6) * t295) * qJD(5)) * t280;
t458 = t185 * t292;
t513 = -t287 * t458 - t115;
t439 = t295 * t297;
t443 = t292 * t294;
t220 = -t281 * t443 - t439;
t441 = t294 * t295;
t442 = t292 * t297;
t221 = t281 * t441 - t442;
t369 = -t221 * rSges(6,1) - t220 * rSges(6,2);
t450 = t280 * t294;
t145 = rSges(6,3) * t450 - t369;
t222 = -t281 * t442 + t441;
t223 = t281 * t439 + t443;
t449 = t280 * t297;
t146 = t223 * rSges(6,1) + t222 * rSges(6,2) + rSges(6,3) * t449;
t512 = -t294 * t145 - t297 * t146;
t368 = rSges(6,1) * t295 - rSges(6,2) * t292;
t118 = t368 * t448 + (rSges(6,3) * t287 + (-rSges(6,1) * t292 - rSges(6,2) * t295) * qJD(5)) * t280;
t187 = -rSges(6,3) * t281 + t280 * t368;
t423 = qJD(1) * t294;
t178 = t187 * t423;
t511 = -t118 * t297 + t178;
t438 = t296 * t297;
t487 = rSges(4,2) * t293;
t510 = t294 * rSges(4,3) - t297 * t487;
t216 = rSges(4,1) * t438 + t510;
t282 = t294 * qJ(3);
t245 = pkin(1) * t438 + t282;
t175 = t216 + t245;
t419 = qJD(2) * t296;
t422 = qJD(1) * t297;
t509 = -t293 * t422 - t294 * t419;
t440 = t294 * t296;
t480 = t297 * rSges(3,3);
t488 = rSges(3,2) * t293;
t215 = rSges(3,1) * t440 - t294 * t488 - t480;
t373 = rSges(3,1) * t296 - t488;
t484 = rSges(3,3) * t294;
t217 = t297 * t373 + t484;
t508 = t215 * t297 - t217 * t294;
t232 = rSges(5,1) * t280 + rSges(5,2) * t281;
t298 = pkin(2) + pkin(1);
t490 = pkin(1) - t298;
t246 = t490 * t293;
t491 = t293 * pkin(1);
t388 = t246 - t491;
t376 = -t232 + t388;
t163 = t376 * t294;
t164 = t376 * t297;
t507 = t163 * t294 + t164 * t297;
t347 = Icges(5,5) * t281 - Icges(5,6) * t280;
t188 = -Icges(5,3) * t297 + t294 * t347;
t506 = qJD(1) * t188;
t348 = Icges(4,5) * t296 - Icges(4,6) * t293;
t202 = -Icges(4,3) * t297 + t294 * t348;
t349 = Icges(3,5) * t296 - Icges(3,6) * t293;
t204 = -Icges(3,3) * t297 + t294 * t349;
t380 = qJD(1) * t281 - qJD(5);
t444 = t287 * t297;
t410 = t280 * t444;
t505 = t294 * t380 + t410;
t229 = Icges(5,2) * t281 + t475;
t230 = Icges(5,1) * t280 + t474;
t331 = t229 * t280 - t230 * t281;
t503 = qJD(1) * t331 + t347 * t287;
t502 = 2 * m(3);
t501 = 2 * m(4);
t500 = 2 * m(5);
t499 = 2 * m(6);
t498 = t294 / 0.2e1;
t497 = -t297 / 0.2e1;
t496 = -rSges(4,1) - pkin(1);
t495 = -rSges(6,3) - pkin(4);
t257 = rSges(3,1) * t293 + rSges(3,2) * t296;
t494 = m(3) * t257;
t493 = m(5) * t232;
t489 = rSges(5,1) * t281;
t486 = rSges(4,2) * t296;
t485 = rSges(5,2) * t280;
t483 = rSges(4,3) * t297;
t139 = Icges(6,5) * t221 + Icges(6,6) * t220 + Icges(6,3) * t450;
t141 = Icges(6,4) * t221 + Icges(6,2) * t220 + Icges(6,6) * t450;
t143 = Icges(6,1) * t221 + Icges(6,4) * t220 + Icges(6,5) * t450;
t345 = -t141 * t292 + t143 * t295;
t381 = -qJD(5) * t281 + qJD(1);
t446 = t287 * t294;
t127 = t381 * t441 + (t280 * t446 - t297 * t380) * t292;
t445 = t287 * t295;
t128 = t380 * t439 + (-t280 * t445 + t292 * t381) * t294;
t409 = t281 * t446;
t307 = t280 * t422 + t409;
t69 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t307;
t71 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t307;
t73 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t307;
t19 = (t287 * t345 - t69) * t281 + (t139 * t287 - t292 * t71 + t295 * t73 + (-t141 * t295 - t143 * t292) * qJD(5)) * t280;
t482 = t19 * t297;
t140 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t449;
t142 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t449;
t144 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t449;
t344 = -t142 * t292 + t144 * t295;
t330 = t381 * t297;
t125 = t292 * t505 + t295 * t330;
t126 = t292 * t330 - t295 * t505;
t401 = t280 * t423;
t408 = t281 * t444;
t306 = -t401 + t408;
t68 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t306;
t70 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t306;
t72 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t306;
t20 = (t287 * t344 - t68) * t281 + (t140 * t287 - t292 * t70 + t295 * t72 + (-t142 * t295 - t144 * t292) * qJD(5)) * t280;
t481 = t20 * t294;
t284 = t294 * rSges(5,3);
t453 = t229 * t287;
t452 = t230 * t287;
t451 = t280 * t287;
t447 = t281 * t297;
t437 = t296 * t298;
t291 = -pkin(3) - qJ(3);
t436 = t297 * t291;
t382 = -t291 * t294 + t297 * t437;
t181 = t382 - t245;
t435 = -t181 - t245;
t371 = -t485 + t489;
t194 = -t297 * rSges(5,3) + t294 * t371;
t195 = rSges(5,1) * t447 - rSges(5,2) * t449 + t284;
t129 = t294 * t194 + t297 * t195;
t399 = t293 * t423;
t270 = pkin(1) * t399;
t434 = -t246 * t423 + t270;
t283 = t297 * qJ(3);
t244 = pkin(1) * t440 - t283;
t433 = t294 * t244 + t297 * t245;
t432 = rSges(5,2) * t401 + rSges(5,3) * t422;
t421 = qJD(2) * t293;
t394 = t298 * t421;
t431 = -t291 * t423 - t294 * t394;
t430 = rSges(4,2) * t399 + rSges(4,3) * t422;
t420 = qJD(2) * t294;
t396 = t293 * t420;
t269 = pkin(1) * t396;
t279 = qJD(3) * t297;
t429 = -t269 - t279;
t275 = qJ(3) * t422;
t278 = qJD(3) * t294;
t428 = t275 + t278;
t189 = Icges(5,3) * t294 + t297 * t347;
t426 = qJD(1) * t189;
t203 = Icges(4,3) * t294 + t297 * t348;
t425 = qJD(1) * t203;
t205 = Icges(3,3) * t294 + t297 * t349;
t424 = qJD(1) * t205;
t418 = qJD(2) * t297;
t416 = pkin(4) * t448;
t263 = pkin(4) * t447;
t414 = pkin(1) * t419;
t61 = -t139 * t281 + t280 * t345;
t184 = -Icges(6,3) * t281 + t280 * t346;
t78 = t184 * t450 + t185 * t220 + t186 * t221;
t413 = t61 / 0.2e1 + t78 / 0.2e1;
t62 = -t140 * t281 + t280 * t344;
t79 = t184 * t449 + t185 * t222 + t186 * t223;
t412 = t79 / 0.2e1 + t62 / 0.2e1;
t117 = t356 * t448 + (Icges(6,5) * t287 + (-Icges(6,1) * t292 - t472) * qJD(5)) * t280;
t407 = t280 * t295 * t117 + t281 * t186 * t445 + t184 * t451;
t400 = t281 * t423;
t308 = -t400 - t410;
t323 = t232 * t287;
t406 = t294 * (-t294 * t323 + (t297 * t371 + t284) * qJD(1)) + t297 * (rSges(5,1) * t308 - rSges(5,2) * t408 + t432) + t194 * t422;
t405 = t126 * rSges(6,1) + t125 * rSges(6,2) + rSges(6,3) * t408;
t395 = t293 * t418;
t397 = t296 * t423;
t404 = t294 * (qJD(1) * t245 + t429) + t297 * ((-t395 - t397) * pkin(1) + t428) + t244 * t422;
t403 = t279 - t431;
t402 = t187 * t422;
t391 = t296 * t490;
t390 = t423 / 0.2e1;
t389 = t422 / 0.2e1;
t256 = rSges(4,1) * t293 + t486;
t387 = -t256 - t491;
t132 = -qJD(1) * t190 - t297 * t453;
t386 = t193 * t287 + t132;
t133 = qJD(1) * t191 - t294 * t453;
t385 = t192 * t287 + t133;
t134 = -qJD(1) * t192 - t297 * t452;
t384 = -t191 * t287 + t134;
t135 = qJD(1) * t193 - t294 * t452;
t383 = t190 * t287 - t135;
t180 = -t294 * t391 + t283 + t436;
t379 = t294 * t180 + t297 * t181 + t433;
t80 = t427 * pkin(4) * t280 - t512;
t378 = -pkin(4) * t451 - t118;
t377 = -t187 + t388;
t375 = qJD(2) * t391 - t414;
t74 = -rSges(6,3) * t401 + t405;
t370 = t128 * rSges(6,1) + t127 * rSges(6,2);
t75 = rSges(6,3) * t307 + t370;
t374 = t145 * t422 + t294 * t75 + t297 * t74 + t427 * t416;
t372 = rSges(4,1) * t296 - t487;
t50 = t139 * t450 + t141 * t220 + t143 * t221;
t51 = t140 * t450 + t142 * t220 + t144 * t221;
t39 = t294 * t51 - t297 * t50;
t367 = t294 * t50 + t297 * t51;
t52 = t139 * t449 + t141 * t222 + t143 * t223;
t53 = t140 * t449 + t142 * t222 + t144 * t223;
t40 = t294 * t53 - t297 * t52;
t366 = t294 * t52 + t297 * t53;
t365 = t294 * t62 - t297 * t61;
t364 = t294 * t61 + t297 * t62;
t228 = Icges(5,5) * t280 + Icges(5,6) * t281;
t319 = t287 * t228;
t130 = -t297 * t319 - t506;
t131 = -t294 * t319 + t426;
t13 = t125 * t141 + t126 * t143 + t139 * t306 + t222 * t71 + t223 * t73 + t449 * t69;
t14 = t125 * t142 + t126 * t144 + t140 * t306 + t222 * t70 + t223 * t72 + t449 * t68;
t8 = qJD(1) * t366 - t13 * t297 + t14 * t294;
t88 = -t188 * t297 - t294 * t337;
t89 = -t189 * t297 - t515;
t90 = t188 * t294 - t514;
t91 = t189 * t294 - t297 * t336;
t363 = t39 * t423 + t40 * t422 + (-t90 * t422 - t88 * t423) * t297 + (t8 + (t91 * qJD(1) + (t133 * t280 - t135 * t281 + t190 * t448 + t192 * t451 - t506) * t297) * t297 + t89 * t423 + t91 * t422 + ((t90 + t515) * qJD(1) + (-t131 + t384 * t281 - t386 * t280 + (t189 - t337) * qJD(1)) * t297 + t130 * t294) * t294) * t294;
t362 = -t437 - t489;
t343 = t145 * t297 - t146 * t294;
t305 = t362 + t485;
t167 = (rSges(5,3) - t291) * t297 + t305 * t294;
t168 = t195 + t382;
t338 = t167 * t297 + t168 * t294;
t201 = t371 * t287;
t329 = -t201 + t375;
t328 = t377 * t297;
t326 = t294 * (t269 + (-t297 * t391 - t282) * qJD(1) + t431) + t297 * (-t275 + t490 * t395 + (t440 * t490 - t436) * qJD(1)) + t180 * t422 + t404;
t325 = t296 * t496 + t487;
t324 = t280 * t495 - t437;
t322 = qJD(2) * t256;
t310 = qJD(2) * (-Icges(3,5) * t293 - Icges(3,6) * t296);
t309 = qJD(2) * (-Icges(4,5) * t293 - Icges(4,6) * t296);
t12 = (t131 * t297 + (t89 + t514) * qJD(1)) * t297 + (t88 * qJD(1) + (-t132 * t280 + t134 * t281 - t191 * t448 - t193 * t451 + t426) * t294 + (-t130 + t383 * t281 + t385 * t280 + (-t188 - t336) * qJD(1)) * t297) * t294;
t15 = t127 * t141 + t128 * t143 + t139 * t307 + t220 * t71 + t221 * t73 + t450 * t69;
t16 = t127 * t142 + t128 * t144 + t140 * t307 + t220 * t70 + t221 * t72 + t450 * t68;
t9 = qJD(1) * t367 - t15 * t297 + t16 * t294;
t304 = (-t12 - t9) * t297 + t363;
t303 = t375 + t378;
t24 = t280 * t367 - t281 * t78;
t25 = t280 * t366 - t281 * t79;
t29 = t115 * t449 + t116 * t222 + t117 * t223 + t125 * t185 + t126 * t186 + t184 * t306;
t3 = (t287 * t366 - t29) * t281 + (-qJD(1) * t40 + t13 * t294 + t14 * t297 + t287 * t79) * t280;
t30 = t115 * t450 + t116 * t220 + t117 * t221 + t127 * t185 + t128 * t186 + t184 * t307;
t4 = (t287 * t367 - t30) * t281 + (-qJD(1) * t39 + t15 * t294 + t16 * t297 + t287 * t78) * t280;
t302 = t3 * t498 + t4 * t497 + t9 * t450 / 0.2e1 - t281 * (qJD(1) * t364 + t481 - t482) / 0.2e1 + t24 * t390 - t40 * t401 / 0.2e1 + t365 * t451 / 0.2e1 + t8 * t449 / 0.2e1 + (t294 * t39 + t297 * t40) * t448 / 0.2e1 + (t280 * t39 + t25) * t389;
t301 = t294 * t324 - t436;
t197 = t351 * t287;
t198 = t357 * t287;
t300 = qJD(1) * t228 + (t198 - t453) * t281 + (-t197 - t452) * t280;
t299 = -t482 / 0.2e1 + t481 / 0.2e1 + (t280 * t384 + t281 * t386 + t294 * t503 + t300 * t297 + t29) * t498 + (-t280 * t383 + t281 * t385 + t300 * t294 - t297 * t503 + t30) * t497 + (t190 * t281 + t192 * t280 - t228 * t297 - t294 * t331 + t61 + t78) * t390 + (t191 * t281 + t193 * t280 + t228 * t294 - t297 * t331 + t62 + t79) * t389;
t261 = t294 * t281 * pkin(4);
t248 = qJD(1) * t263;
t241 = t373 * qJD(2);
t240 = t372 * qJD(2);
t214 = t294 * t372 - t483;
t200 = t387 * t297;
t199 = t387 * t294;
t174 = t294 * t325 + t283 + t483;
t166 = -t187 * t297 + t263;
t165 = -t187 * t294 + t261;
t162 = -rSges(3,1) * t396 + (rSges(3,1) * t438 + t484) * qJD(1) + t509 * rSges(3,2);
t161 = -t257 * t418 + (-t294 * t373 + t480) * qJD(1);
t152 = t294 * t310 + t424;
t151 = -qJD(1) * t204 + t297 * t310;
t150 = t294 * t309 + t425;
t149 = -qJD(1) * t202 + t297 * t309;
t148 = pkin(1) * t509 - t240 * t294 - t256 * t422;
t147 = t256 * t423 + t270 + (-t240 - t414) * t297;
t114 = t263 + t328;
t113 = t294 * t377 + t261;
t107 = t256 * t420 + ((-rSges(4,3) - qJ(3)) * t294 + t325 * t297) * qJD(1) - t429;
t106 = t496 * t397 + (t293 * t496 - t486) * t418 + t428 + t430;
t105 = pkin(4) * t449 + t146 + t382;
t104 = t301 + t369;
t103 = t205 * t294 - t332 * t297;
t102 = t204 * t294 - t316;
t101 = t203 * t294 - t334 * t297;
t100 = t202 * t294 - t318;
t99 = -t205 * t297 - t315;
t98 = -t204 * t297 - t294 * t333;
t97 = -t203 * t297 - t317;
t96 = -t202 * t297 - t294 * t335;
t95 = t232 * t446 + (t297 * t305 - t284) * qJD(1) + t403;
t94 = t278 + t362 * t423 + (-qJD(1) * t291 - t323 - t394) * t297 + t432;
t93 = -t146 * t281 - t187 * t449;
t92 = t145 * t281 + t187 * t450;
t87 = qJD(1) * t164 + t294 * t329;
t86 = t232 * t423 + t297 * t329 + t434;
t85 = -t184 * t281 + (t186 * t295 - t458) * t280;
t84 = t85 * t451;
t83 = t343 * t280;
t82 = t294 * t378 + t248 - t402;
t81 = pkin(4) * t308 + t511;
t65 = t379 + t129;
t64 = qJD(1) * t328 + t294 * t303 + t248;
t63 = -pkin(4) * t400 + t297 * t303 + t178 + t434;
t60 = -t195 * t423 + t406;
t55 = t324 * t422 + t409 * t495 - t370 + t403;
t54 = t278 + (-t394 + t416) * t297 + t301 * qJD(1) + t405;
t45 = t80 + t379;
t44 = (t187 * t446 + t75) * t281 + (t118 * t294 - t145 * t287 + t402) * t280;
t43 = (-t187 * t444 - t74) * t281 + (t146 * t287 + t511) * t280;
t41 = t516 * t280 + t513 * t281 + t407;
t32 = -t146 * t423 + t374;
t31 = (-t195 + t435) * t423 + t326 + t406;
t26 = t343 * t448 + (qJD(1) * t512 - t294 * t74 + t297 * t75) * t280;
t21 = (-t146 + t435) * t423 + t326 + t374;
t1 = [t407 + (t104 * t55 + t105 * t54) * t499 + (t167 * t95 + t168 * t94) * t500 + (t106 * t175 + t107 * t174) * t501 + (t161 * t217 + t162 * t215) * t502 - t229 * t451 + t230 * t448 + (t197 + t513) * t281 + (t359 + t361 + t519) * t421 + (t353 + t355 - t518) * t419 + (t198 + t516) * t280; m(6) * (t104 * t63 + t105 * t64 + t113 * t54 + t114 * t55) + m(5) * (t163 * t94 + t164 * t95 + t167 * t86 + t168 * t87) + m(4) * (t106 * t199 + t107 * t200 + t147 * t174 + t148 * t175) + t299 + m(3) * ((-t161 * t294 + t162 * t297) * t257 + t508 * t241) + ((-t217 * t494 + (t207 / 0.2e1 + t209 / 0.2e1) * t296 + (t211 / 0.2e1 + t213 / 0.2e1) * t293) * t297 + (-t215 * t494 + (t206 / 0.2e1 + t208 / 0.2e1) * t296 + (t210 / 0.2e1 + t212 / 0.2e1) * t293) * t294) * qJD(1) + (t517 * t297 + ((-t206 - t208) * t296 + (-t210 - t212) * t293) * qJD(1)) * t498 + (t517 * t294 + ((t207 + t209) * t296 + (t211 + t213) * t293) * qJD(1)) * t497 + (t316 / 0.2e1 + t318 / 0.2e1 - t317 / 0.2e1 - t315 / 0.2e1 + (t348 + t349) * (t289 / 0.2e1 + t288 / 0.2e1)) * qJD(2); (t200 * t147 + t199 * t148 + (t214 * t294 + t216 * t297 + t433) * ((qJD(1) * t214 - t297 * t322 + t430) * t297 + (-t294 * t322 + (t510 - t175) * qJD(1)) * t294 + t404)) * t501 + t363 + t294 * ((t294 * t151 + (t102 + t315) * qJD(1)) * t294 + (t103 * qJD(1) + (t208 * t419 + t212 * t421) * t297 + (-t152 + (-t209 * t296 - t213 * t293) * qJD(2) + (t205 - t333) * qJD(1)) * t294) * t297) - t297 * ((t297 * t152 + (t99 + t316) * qJD(1)) * t297 + (t98 * qJD(1) + (-t209 * t419 - t213 * t421 + t424) * t294 + (-t151 + (t208 * t296 + t212 * t293) * qJD(2) - t332 * qJD(1)) * t297) * t294) + t294 * ((t294 * t149 + (t100 + t317) * qJD(1)) * t294 + (t101 * qJD(1) + (t206 * t419 + t210 * t421) * t297 + (-t150 + (-t207 * t296 - t211 * t293) * qJD(2) + (t203 - t335) * qJD(1)) * t294) * t297) - t297 * ((t297 * t150 + (t97 + t318) * qJD(1)) * t297 + (t96 * qJD(1) + (-t207 * t419 - t211 * t421 + t425) * t294 + (-t149 + (t206 * t296 + t210 * t293) * qJD(2) - t334 * qJD(1)) * t297) * t294) + (t113 * t64 + t114 * t63 + t21 * t45) * t499 + (t163 * t87 + t164 * t86 + t31 * t65) * t500 + ((t215 * t294 + t217 * t297) * (qJD(1) * t508 + t297 * t161 + t294 * t162) + t427 * t257 * t241) * t502 - t297 * t12 - t297 * t9 + ((-t96 - t98) * t297 + (t97 + t99) * t294) * t423 + ((-t100 - t102) * t297 + (t101 + t103) * t294) * t422; m(6) * (t294 * t55 - t297 * t54 + (t104 * t297 + t105 * t294) * qJD(1)) + m(5) * (qJD(1) * t338 + t294 * t95 - t297 * t94) + m(4) * (-t106 * t297 + t107 * t294 + (t174 * t297 + t175 * t294) * qJD(1)); m(6) * (t294 * t63 - t297 * t64 + (t113 * t294 + t114 * t297) * qJD(1)) + m(5) * (qJD(1) * t507 + t294 * t86 - t297 * t87) + m(4) * (t147 * t294 - t148 * t297 + (t199 * t294 + t200 * t297) * qJD(1)); 0; t299 + m(6) * (t104 * t81 + t105 * t82 + t165 * t54 + t166 * t55) - m(5) * t338 * t201 + (-t294 * t94 - t297 * t95 + (t167 * t294 - t168 * t297) * qJD(1)) * t493; m(6) * (t113 * t82 + t114 * t81 + t165 * t64 + t166 * t63 + t21 * t80 + t32 * t45) + m(5) * (t129 * t31 - t201 * t507 + t60 * t65) + (-t294 * t87 - t297 * t86 + (-t163 * t297 + t164 * t294) * qJD(1)) * t493 + t304; m(6) * (t294 * t81 - t297 * t82 + (t165 * t294 + t166 * t297) * qJD(1)); (t201 * t232 * t427 + t129 * t60) * t500 + (t165 * t82 + t166 * t81 + t32 * t80) * t499 + t304; t84 + m(6) * (t104 * t44 + t105 * t43 + t54 * t93 + t55 * t92) + (-t41 + (t294 * t413 + t297 * t412) * t287) * t281 + ((t29 / 0.2e1 + t20 / 0.2e1) * t297 + (t19 / 0.2e1 + t30 / 0.2e1) * t294 + (-t294 * t412 + t297 * t413) * qJD(1)) * t280; t302 + m(6) * (t113 * t43 + t114 * t44 + t21 * t83 + t26 * t45 + t63 * t92 + t64 * t93); m(6) * (t294 * t44 - t297 * t43 + (t294 * t93 + t297 * t92) * qJD(1)); t302 + m(6) * (t165 * t43 + t166 * t44 + t26 * t80 + t32 * t83 + t81 * t92 + t82 * t93); (t26 * t83 + t43 * t93 + t44 * t92) * t499 + (t41 * t281 - t84 + (t294 * t24 + t297 * t25 - t281 * t364) * t287) * t281 + (t297 * t3 + t294 * t4 + t364 * t451 + (-t19 * t294 - t20 * t297 - t287 * t85) * t281 + (t297 * t24 - t294 * t25 + t281 * t365) * qJD(1)) * t280;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
