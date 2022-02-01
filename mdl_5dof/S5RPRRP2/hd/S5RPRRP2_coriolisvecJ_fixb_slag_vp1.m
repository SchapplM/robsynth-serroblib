% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m [6x1]
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:27:52
% EndTime: 2022-01-23 09:28:05
% DurationCPUTime: 10.48s
% Computational Cost: add. (12066->401), mult. (8983->501), div. (0->0), fcn. (6890->8), ass. (0->248)
t251 = sin(qJ(4));
t253 = cos(qJ(4));
t203 = Icges(6,5) * t253 - Icges(6,6) * t251;
t205 = Icges(5,5) * t253 - Icges(5,6) * t251;
t499 = t203 + t205;
t516 = Icges(5,5) + Icges(6,5);
t515 = Icges(5,6) + Icges(6,6);
t514 = Icges(5,3) + Icges(6,3);
t249 = qJ(1) + pkin(8);
t244 = qJ(3) + t249;
t240 = cos(t244);
t239 = sin(t244);
t392 = t239 * t253;
t393 = t239 * t251;
t124 = Icges(6,4) * t392 - Icges(6,2) * t393 - Icges(6,6) * t240;
t126 = Icges(5,4) * t392 - Icges(5,2) * t393 - Icges(5,6) * t240;
t504 = t124 + t126;
t196 = Icges(6,4) * t393;
t128 = Icges(6,1) * t392 - Icges(6,5) * t240 - t196;
t197 = Icges(5,4) * t393;
t130 = Icges(5,1) * t392 - Icges(5,5) * t240 - t197;
t511 = t128 + t130;
t411 = Icges(6,4) * t251;
t211 = Icges(6,1) * t253 - t411;
t294 = t211 * t240;
t129 = Icges(6,5) * t239 + t294;
t412 = Icges(5,4) * t251;
t213 = Icges(5,1) * t253 - t412;
t295 = t213 * t240;
t131 = Icges(5,5) * t239 + t295;
t502 = t129 + t131;
t206 = Icges(6,2) * t253 + t411;
t208 = Icges(5,2) * t253 + t412;
t513 = t206 + t208;
t245 = Icges(6,4) * t253;
t210 = Icges(6,1) * t251 + t245;
t246 = Icges(5,4) * t253;
t212 = Icges(5,1) * t251 + t246;
t509 = -t210 - t212;
t512 = t499 * t240;
t483 = t240 * t514 - t392 * t516 + t393 * t515;
t482 = t239 * t514 + t512;
t508 = t211 + t213;
t316 = -Icges(6,2) * t251 + t245;
t317 = -Icges(5,2) * t251 + t246;
t507 = t316 + t317;
t506 = t504 * t251;
t505 = t502 * t392;
t292 = t316 * t240;
t125 = Icges(6,6) * t239 + t292;
t293 = t317 * t240;
t127 = Icges(5,6) * t239 + t293;
t503 = t125 + t127;
t495 = t251 * t513 + t253 * t509;
t478 = -t253 * t511 + t506;
t501 = t507 * qJD(4);
t500 = t508 * qJD(4);
t202 = Icges(6,5) * t251 + Icges(6,6) * t253;
t204 = Icges(5,5) * t251 + Icges(5,6) * t253;
t498 = t204 + t202;
t248 = qJD(1) + qJD(3);
t497 = -qJD(4) * t513 + t248 * t515;
t496 = t509 * qJD(4) + t248 * t516;
t494 = -t240 * t482 + t505;
t389 = t240 * t253;
t493 = t239 * t483 - t389 * t511;
t455 = t239 * t482 + t389 * t502;
t492 = -t239 * t478 + t240 * t483;
t465 = -t393 * t503 + t494;
t390 = t240 * t251;
t464 = -t390 * t504 - t493;
t463 = -t390 * t503 + t455;
t397 = t204 * t240;
t400 = t202 * t240;
t491 = -t239 * t495 - t397 - t400;
t398 = t204 * t239;
t401 = t202 * t239;
t490 = -t240 * t495 + t398 + t401;
t235 = t240 * pkin(7);
t163 = pkin(3) * t239 - t235;
t250 = -qJ(5) - pkin(7);
t214 = t240 * t250;
t427 = pkin(4) * t253;
t241 = pkin(3) + t427;
t287 = -rSges(6,1) * t392 + rSges(6,2) * t393 + rSges(6,3) * t240 - t239 * t241 - t214;
t475 = t163 + t287;
t489 = t475 * t248;
t488 = t503 * t251;
t462 = t251 * t511 + t504 * t253;
t461 = t251 * t502 + t253 * t503;
t395 = t239 * t248;
t487 = t497 * t240 - t395 * t507;
t486 = (t292 + t293) * t248 + t497 * t239;
t485 = t496 * t240 - t395 * t508;
t484 = (-t294 - t295) * t248 - t496 * t239;
t481 = t500 * t253 - t501 * t251 + t498 * t248 + (t251 * t509 - t253 * t513) * qJD(4);
t480 = -t498 * qJD(4) + t248 * t514;
t479 = -t253 * t502 + t488;
t477 = qJD(4) * t499 + t495 * t248;
t476 = t490 * t248;
t474 = (t239 * t463 - t240 * t464) * qJD(4);
t473 = (t239 * t465 - t240 * t492) * qJD(4);
t472 = t491 * t248;
t256 = qJD(1) ^ 2;
t471 = t472 + t473;
t470 = t474 + t476;
t469 = qJD(4) * t478 + t251 * t484 - t253 * t486;
t468 = -qJD(4) * t479 + t251 * t485 + t253 * t487;
t467 = t239 * t477 + t240 * t481;
t466 = t239 * t481 - t240 * t477;
t460 = -qJD(4) * t461 + t248 * t482 - t251 * t487 + t253 * t485;
t459 = qJD(4) * t462 + t248 * t483 + t251 * t486 + t253 * t484;
t458 = t483 + t488;
t359 = qJD(4) * t253;
t387 = t248 * t251;
t456 = t239 * t359 + t240 * t387;
t454 = (t478 + t512) * t248 + t480 * t239;
t453 = t240 * t480 + t248 * t479 - t395 * t499;
t451 = 0.2e1 * qJD(4);
t369 = rSges(5,2) * t393 + rSges(5,3) * t240;
t133 = rSges(5,1) * t392 - t369;
t119 = t248 * t133;
t158 = t248 * t163;
t450 = -t119 - t158;
t368 = rSges(6,1) * t389 + rSges(6,3) * t239;
t449 = rSges(5,1) * t389 + rSges(5,3) * t239;
t221 = qJD(5) * t239;
t419 = t253 * rSges(6,2);
t215 = rSges(6,1) * t251 + t419;
t337 = pkin(4) * t251 + t215;
t361 = qJD(4) * t240;
t285 = -t337 * t361 + t221;
t358 = t239 * t387;
t391 = t240 * t248;
t448 = rSges(6,2) * t358 + rSges(6,3) * t391 + t221;
t234 = t239 * pkin(7);
t164 = pkin(3) * t240 + t234;
t447 = pkin(2) * cos(t249) + cos(qJ(1)) * pkin(1);
t394 = t239 * t250;
t308 = t240 * t241 + t368 - t394;
t360 = qJD(4) * t251;
t350 = t239 * t360;
t222 = qJD(5) * t240;
t370 = pkin(4) * t350 + t222;
t445 = -rSges(6,1) * t350 - rSges(6,2) * t456 - t248 * t394 - t370;
t444 = -t158 + t489;
t216 = rSges(5,1) * t251 + rSges(5,2) * t253;
t362 = qJD(4) * t239;
t160 = t216 * t362;
t135 = -rSges(5,2) * t390 + t449;
t286 = t135 + t164;
t438 = -t248 * t286 + t160;
t381 = -rSges(6,2) * t390 - t164 + t308;
t435 = t215 * t362 - t248 * (t164 + t381) + t370;
t374 = -Icges(5,2) * t392 + t130 - t197;
t378 = t212 * t239 + t126;
t434 = -t251 * t374 - t253 * t378;
t376 = -Icges(6,2) * t392 + t128 - t196;
t380 = t210 * t239 + t124;
t433 = -t251 * t376 - t253 * t380;
t432 = t239 / 0.2e1;
t431 = -t240 / 0.2e1;
t429 = t248 / 0.2e1;
t426 = pkin(3) - t241;
t191 = pkin(7) * t391;
t347 = t240 * t360;
t357 = t248 * t392;
t288 = -t347 - t357;
t346 = t240 * t359;
t425 = -pkin(4) * t347 - t191 + (t239 * t426 - t214) * t248 + rSges(6,1) * t288 - rSges(6,2) * t346 + t448;
t424 = t445 + (-t240 * t426 - t234 + t368) * t248;
t423 = rSges(5,1) * t253;
t422 = rSges(6,1) * t253;
t420 = rSges(6,2) * t251;
t233 = t240 * rSges(4,1);
t325 = -pkin(2) * sin(t249) - sin(qJ(1)) * pkin(1);
t303 = t325 * qJD(1);
t348 = t216 * t361;
t272 = t303 - t348;
t60 = (-t133 - t163) * t248 + t272;
t417 = t60 * t240;
t162 = -rSges(4,2) * t239 + t233;
t302 = t447 * qJD(1);
t117 = t162 * t248 + t302;
t161 = rSges(4,1) * t239 + rSges(4,2) * t240;
t408 = t117 * t161;
t399 = t203 * t248;
t396 = t205 * t248;
t388 = t248 * t161;
t379 = -t210 * t240 - t125;
t377 = -t212 * t240 - t127;
t375 = -t206 * t240 + t129;
t373 = -t208 * t240 + t131;
t367 = -t206 + t211;
t366 = t210 + t316;
t365 = -t208 + t213;
t364 = t212 + t317;
t352 = rSges(5,1) * t350 + rSges(5,2) * t456;
t343 = -pkin(3) - t423;
t342 = -t362 / 0.2e1;
t339 = t361 / 0.2e1;
t218 = -t420 + t422;
t336 = -t218 - t427;
t327 = -pkin(4) * t390 - t215 * t240;
t180 = t218 * qJD(4);
t326 = -pkin(4) * t359 - t180;
t321 = -rSges(5,2) * t251 + t423;
t61 = t302 - t438;
t320 = -t239 * t61 - t417;
t311 = t133 * t239 + t135 * t240;
t306 = (-qJD(4) * t427 - t180) * qJD(4);
t305 = t325 * t256;
t304 = t447 * t256;
t296 = rSges(5,3) * t391 + (-t346 + t358) * rSges(5,2);
t284 = t248 * (-pkin(3) * t395 + t191) + t305;
t283 = -t251 * t375 + t253 * t379;
t282 = -t251 * t373 + t253 * t377;
t281 = t239 * t343 + t235 + t369;
t116 = t303 - t388;
t280 = (-t251 * t366 + t253 * t367) * t248;
t279 = (-t251 * t364 + t253 * t365) * t248;
t92 = rSges(5,1) * t288 + t296;
t94 = t248 * t449 - t352;
t271 = (t92 + t119) * t240 + (-t135 * t248 + t94) * t239;
t266 = t303 + t285;
t259 = ((t455 * t239 + ((t482 + t506) * t240 + t465 + t493 - t505) * t240) * qJD(4) + t476) * t339 + (-qJD(4) * t495 + t251 * t500 + t253 * t501) * t248 + (((t240 * t458 - t455 + t463) * t240 + (t239 * t458 + t464 - t494) * t239) * qJD(4) + t471 - t472) * t342 + (t467 + t468) * t362 / 0.2e1 - (t466 - t469 + t470) * t361 / 0.2e1 + ((t462 + t491) * t239 + (t461 + t490) * t240) * qJD(4) * t429;
t258 = t60 * t352 + t61 * (-rSges(5,1) * t347 + t191 + t296) + (t343 * t417 + (t60 * (-rSges(5,3) - pkin(7)) + t61 * t343) * t239) * t248;
t23 = t306 * t239 + (t425 + t285) * t248 + t284;
t50 = (-t163 + t475) * t248 + t266;
t51 = t302 - t435;
t257 = t50 * (-rSges(6,3) * t395 - t445) + t51 * (-rSges(6,1) * t357 - t241 * t395 + t448) + (-t23 * t420 + (t50 * (-t241 - t422) - t51 * t250) * t248 + t51 * (-t419 + (-rSges(6,1) - pkin(4)) * t251) * qJD(4)) * t240;
t188 = rSges(4,2) * t395;
t181 = t321 * qJD(4);
t157 = t216 * t240;
t155 = t216 * t239;
t154 = t215 * t239;
t139 = t164 * t248;
t138 = rSges(4,1) * t391 - t188;
t103 = -t138 * t248 - t304;
t102 = -t248 * t388 + t305;
t66 = qJD(4) * t311 + qJD(2);
t49 = -t181 * t361 - t304 + (-t139 - t94 + t160) * t248;
t48 = t248 * t92 + (-t181 * t239 - t216 * t391) * qJD(4) + t284;
t42 = qJD(2) + (-t239 * t475 + t240 * t381) * qJD(4);
t25 = t271 * qJD(4);
t24 = -t304 + t306 * t240 + (t337 * t362 - t139 + t222 - t424) * t248;
t5 = ((t425 - t489) * t240 + (-t248 * t381 + t424) * t239) * qJD(4);
t1 = [m(4) * (t103 * (-t161 + t325) + t116 * t188 + t102 * (t162 + t447) + (-t116 * t233 - t408) * t248 + (-t116 * t447 + t117 * t325) * qJD(1)) + t259 + (t24 * (t287 + t325) + t23 * (t308 + t447) + (t325 * t51 - t447 * t50) * qJD(1) + t257 - (-t50 + t266 + t444) * t51) * m(6) + (t49 * (t281 + t325) + t48 * (t286 + t447) + (t325 * t61 - t447 * t60) * qJD(1) + t258 - (-t60 + t272 + t450) * t61) * m(5); m(5) * t25 + m(6) * t5; t259 + (t23 * t308 + t24 * t287 + t257 - t51 * (t285 + t444) - t435 * t50) * m(6) + (t49 * t281 + t48 * t286 + t258 - t60 * t438 - t61 * (-t348 + t450)) * m(5) + (-(-t116 * t162 - t408) * t248 + t102 * t162 - t103 * t161 - t116 * t138 - t117 * t388) * m(4); -(((t364 + t366) * t253 + (t365 + t367) * t251) * t248 + (((-t374 - t376) * t240 + (t373 + t375) * t239) * t253 + ((t378 + t380) * t240 + (t377 + t379) * t239) * t251) * qJD(4)) * t248 / 0.2e1 + ((t248 * t461 + t469) * t240 + (t248 * t462 + t468) * t239) * t429 + ((-t362 * t397 + t396) * t239 + (t279 + (-t434 * t240 + (t398 + t282) * t239) * qJD(4)) * t240 + (-t362 * t400 + t399) * t239 + (t280 + (-t433 * t240 + (t401 + t283) * t239) * qJD(4)) * t240) * t342 + ((-t361 * t398 - t396) * t240 + (t279 + (t282 * t239 + (t397 - t434) * t240) * qJD(4)) * t239 + (-t361 * t401 - t399) * t240 + (t280 + (t283 * t239 + (t400 - t433) * t240) * qJD(4)) * t239) * t339 + (t25 * t311 + t66 * t271 + t320 * t181 + ((-t248 * t61 - t49) * t240 + (t248 * t60 - t48) * t239) * t216 - (t155 * t60 - t157 * t61) * t248 - (t66 * (-t155 * t239 - t157 * t240) + t320 * t321) * qJD(4)) * m(5) + (t467 * t248 + ((t459 * t240 + t463 * t248) * t240 + (t453 * t239 + t464 * t248 + (-t454 + t460) * t240) * t239) * t451) * t432 + (t466 * t248 + ((t454 * t240 + t465 * t248) * t240 + (t460 * t239 + t492 * t248 + (-t453 + t459) * t240) * t239) * t451) * t431 + ((-t23 * t337 + t51 * t326 - t5 * t475 + t42 * t424 + (t215 * t50 - t381 * t42) * t248) * t239 + (-t24 * t337 + t50 * t326 + t5 * t381 + t42 * t425 + (-t337 * t51 - t42 * t475) * t248) * t240 - (t50 * t154 + t327 * t51) * t248 - ((t327 * t42 + t336 * t50) * t240 + (t51 * t336 + (-pkin(4) * t393 - t154) * t42) * t239) * qJD(4)) * m(6) + (t471 + t473) * t395 / 0.2e1 + (t470 + t474) * t391 / 0.2e1; 0.2e1 * (t23 * t431 + t24 * t432) * m(6);];
tauc = t1(:);
