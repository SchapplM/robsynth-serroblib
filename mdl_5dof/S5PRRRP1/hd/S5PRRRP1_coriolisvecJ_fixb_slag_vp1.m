% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:47
% EndTime: 2019-12-05 16:40:07
% DurationCPUTime: 10.68s
% Computational Cost: add. (11958->399), mult. (8743->487), div. (0->0), fcn. (6780->6), ass. (0->249)
t246 = sin(qJ(4));
t247 = cos(qJ(4));
t202 = Icges(6,5) * t247 - Icges(6,6) * t246;
t204 = Icges(5,5) * t247 - Icges(5,6) * t246;
t486 = t202 + t204;
t503 = Icges(5,5) + Icges(6,5);
t502 = Icges(5,6) + Icges(6,6);
t501 = Icges(5,3) + Icges(6,3);
t243 = pkin(8) + qJ(2);
t240 = qJ(3) + t243;
t236 = cos(t240);
t235 = sin(t240);
t383 = t235 * t247;
t384 = t235 * t246;
t122 = Icges(6,4) * t383 - Icges(6,2) * t384 - Icges(6,6) * t236;
t124 = Icges(5,4) * t383 - Icges(5,2) * t384 - Icges(5,6) * t236;
t491 = t122 + t124;
t195 = Icges(6,4) * t384;
t126 = Icges(6,1) * t383 - Icges(6,5) * t236 - t195;
t196 = Icges(5,4) * t384;
t128 = Icges(5,1) * t383 - Icges(5,5) * t236 - t196;
t498 = t126 + t128;
t400 = Icges(6,4) * t246;
t210 = Icges(6,1) * t247 - t400;
t288 = t210 * t236;
t127 = Icges(6,5) * t235 + t288;
t401 = Icges(5,4) * t246;
t212 = Icges(5,1) * t247 - t401;
t289 = t212 * t236;
t129 = Icges(5,5) * t235 + t289;
t489 = t127 + t129;
t205 = Icges(6,2) * t247 + t400;
t207 = Icges(5,2) * t247 + t401;
t500 = t205 + t207;
t241 = Icges(6,4) * t247;
t209 = Icges(6,1) * t246 + t241;
t242 = Icges(5,4) * t247;
t211 = Icges(5,1) * t246 + t242;
t496 = -t209 - t211;
t499 = t486 * t236;
t470 = t236 * t501 - t383 * t503 + t384 * t502;
t469 = t235 * t501 + t499;
t495 = t210 + t212;
t306 = -Icges(6,2) * t246 + t241;
t307 = -Icges(5,2) * t246 + t242;
t494 = t306 + t307;
t493 = t491 * t246;
t492 = t489 * t383;
t286 = t306 * t236;
t123 = Icges(6,6) * t235 + t286;
t287 = t307 * t236;
t125 = Icges(5,6) * t235 + t287;
t490 = t123 + t125;
t482 = t246 * t500 + t247 * t496;
t465 = -t247 * t498 + t493;
t488 = t494 * qJD(4);
t487 = t495 * qJD(4);
t201 = Icges(6,5) * t246 + Icges(6,6) * t247;
t203 = Icges(5,5) * t246 + Icges(5,6) * t247;
t485 = t203 + t201;
t244 = qJD(2) + qJD(3);
t484 = -qJD(4) * t500 + t244 * t502;
t483 = t496 * qJD(4) + t244 * t503;
t481 = -t236 * t469 + t492;
t378 = t236 * t247;
t480 = t235 * t470 - t378 * t498;
t441 = t235 * t469 + t378 * t489;
t479 = -t235 * t465 + t236 * t470;
t452 = -t384 * t490 + t481;
t379 = t236 * t246;
t451 = -t379 * t491 - t480;
t450 = -t379 * t490 + t441;
t388 = t203 * t236;
t391 = t201 * t236;
t478 = -t235 * t482 - t388 - t391;
t389 = t203 * t235;
t392 = t201 * t235;
t477 = -t236 * t482 + t389 + t392;
t232 = t236 * pkin(7);
t162 = pkin(3) * t235 - t232;
t245 = -qJ(5) - pkin(7);
t213 = t236 * t245;
t416 = pkin(4) * t247;
t237 = pkin(3) + t416;
t281 = -rSges(6,1) * t383 + rSges(6,2) * t384 + rSges(6,3) * t236 - t235 * t237 - t213;
t462 = t162 + t281;
t476 = t462 * t244;
t475 = t490 * t246;
t449 = t246 * t498 + t247 * t491;
t448 = t246 * t489 + t247 * t490;
t386 = t235 * t244;
t474 = t236 * t484 - t386 * t494;
t473 = (t286 + t287) * t244 + t484 * t235;
t472 = t236 * t483 - t386 * t495;
t471 = (-t288 - t289) * t244 - t483 * t235;
t468 = t487 * t247 - t488 * t246 + t485 * t244 + (t246 * t496 - t247 * t500) * qJD(4);
t467 = -t485 * qJD(4) + t244 * t501;
t466 = -t247 * t489 + t475;
t464 = qJD(4) * t486 + t482 * t244;
t463 = t477 * t244;
t461 = (t235 * t450 - t236 * t451) * qJD(4);
t460 = (t235 * t452 - t236 * t479) * qJD(4);
t459 = t478 * t244;
t458 = t459 + t460;
t457 = t461 + t463;
t456 = qJD(4) * t465 + t246 * t471 - t247 * t473;
t455 = -qJD(4) * t466 + t246 * t472 + t247 * t474;
t454 = t235 * t464 + t236 * t468;
t453 = t235 * t468 - t236 * t464;
t447 = -qJD(4) * t448 + t244 * t469 - t246 * t474 + t247 * t472;
t446 = qJD(4) * t449 + t244 * t470 + t246 * t473 + t247 * t471;
t359 = rSges(5,2) * t384 + rSges(5,3) * t236;
t131 = rSges(5,1) * t383 - t359;
t117 = t244 * t131;
t156 = t244 * t162;
t380 = t236 * t244;
t190 = pkin(7) * t380;
t350 = qJD(4) * t247;
t333 = t236 * t350;
t377 = t244 * t246;
t345 = t235 * t377;
t290 = rSges(5,3) * t380 + (-t333 + t345) * rSges(5,2);
t351 = qJD(4) * t246;
t334 = t236 * t351;
t445 = -rSges(5,1) * t334 + t117 + t156 + t190 + t290;
t344 = t244 * t383;
t218 = qJD(5) * t235;
t435 = rSges(6,2) * t345 + rSges(6,3) * t380 + t218;
t444 = -rSges(6,1) * t344 - t237 * t386 + t156 + t435 - t476;
t443 = t470 + t475;
t442 = t235 * t350 + t236 * t377;
t238 = sin(t243);
t409 = pkin(2) * qJD(2);
t347 = t238 * t409;
t160 = rSges(4,1) * t235 + rSges(4,2) * t236;
t393 = t160 * t244;
t134 = -t347 - t393;
t440 = (t465 + t499) * t244 + t467 * t235;
t439 = t236 * t467 + t244 * t466 - t386 * t486;
t438 = 0.2e1 * qJD(4);
t358 = rSges(6,1) * t378 + rSges(6,3) * t235;
t436 = rSges(5,1) * t378 + rSges(5,3) * t235;
t407 = t247 * rSges(6,2);
t214 = rSges(6,1) * t246 + t407;
t324 = pkin(4) * t246 + t214;
t352 = qJD(4) * t236;
t279 = -t324 * t352 + t218;
t231 = t235 * pkin(7);
t163 = pkin(3) * t236 + t231;
t385 = t235 * t245;
t298 = t236 * t237 + t358 - t385;
t336 = t235 * t351;
t219 = qJD(5) * t236;
t360 = pkin(4) * t336 + t219;
t434 = -rSges(6,1) * t336 - rSges(6,2) * t442 - t244 * t385 - t360;
t215 = rSges(5,1) * t246 + rSges(5,2) * t247;
t353 = qJD(4) * t235;
t158 = t215 * t353;
t133 = -rSges(5,2) * t379 + t436;
t280 = t133 + t163;
t428 = -t244 * t280 + t158;
t371 = -rSges(6,2) * t379 - t163 + t298;
t425 = t214 * t353 - t244 * (t163 + t371) + t360;
t364 = -Icges(5,2) * t383 + t128 - t196;
t368 = t211 * t235 + t124;
t424 = -t246 * t364 - t247 * t368;
t366 = -Icges(6,2) * t383 + t126 - t195;
t370 = t209 * t235 + t122;
t423 = -t246 * t366 - t247 * t370;
t422 = t235 / 0.2e1;
t421 = -t236 / 0.2e1;
t419 = t244 / 0.2e1;
t418 = pkin(2) * t238;
t417 = pkin(2) * qJD(2) ^ 2;
t415 = pkin(3) - t237;
t282 = -t334 - t344;
t414 = -pkin(4) * t334 - t190 + (t235 * t415 - t213) * t244 + rSges(6,1) * t282 - rSges(6,2) * t333 + t435;
t413 = t434 + (-t236 * t415 - t231 + t358) * t244;
t412 = rSges(5,1) * t247;
t411 = rSges(6,1) * t247;
t410 = rSges(6,2) * t246;
t337 = t215 * t352;
t283 = -t337 - t347;
t64 = (-t131 - t162) * t244 + t283;
t408 = t236 * t64;
t390 = t202 * t244;
t387 = t204 * t244;
t369 = -t209 * t236 - t123;
t367 = -t211 * t236 - t125;
t365 = -t205 * t236 + t127;
t363 = -t207 * t236 + t129;
t357 = -t205 + t210;
t356 = t209 + t306;
t355 = -t207 + t212;
t354 = t211 + t307;
t349 = t238 * t417;
t239 = cos(t243);
t348 = t239 * t417;
t346 = t239 * t409;
t339 = rSges(5,1) * t336 + rSges(5,2) * t442;
t330 = -pkin(3) - t412;
t329 = -t353 / 0.2e1;
t326 = t352 / 0.2e1;
t216 = -t410 + t411;
t323 = -t216 - t416;
t161 = rSges(4,1) * t236 - rSges(4,2) * t235;
t315 = t244 * (-pkin(3) * t386 + t190) - t349;
t314 = -pkin(4) * t379 - t214 * t236;
t138 = rSges(4,1) * t380 - rSges(4,2) * t386;
t180 = t216 * qJD(4);
t313 = -pkin(4) * t350 - t180;
t311 = -rSges(5,2) * t246 + t412;
t65 = t346 - t428;
t310 = -t235 * t65 - t408;
t301 = t131 * t235 + t133 * t236;
t296 = (-qJD(4) * t416 - t180) * qJD(4);
t278 = -t246 * t365 + t247 * t369;
t277 = -t246 * t363 + t247 * t367;
t276 = t235 * t330 + t232 + t359;
t275 = -rSges(6,3) * t386 - t434;
t274 = (-t246 * t356 + t247 * t357) * t244;
t273 = (-t246 * t354 + t247 * t355) * t244;
t266 = t279 - t347;
t92 = rSges(5,1) * t282 + t290;
t94 = t244 * t436 - t339;
t264 = (t92 + t117) * t236 + (-t133 * t244 + t94) * t235;
t252 = (t330 * t408 + (t64 * (-rSges(5,3) - pkin(7)) + t65 * t330) * t235) * t244;
t251 = ((t441 * t235 + ((t469 + t493) * t236 + t452 + t480 - t492) * t236) * qJD(4) + t463) * t326 + (-qJD(4) * t482 + t246 * t487 + t247 * t488) * t244 + (((t236 * t443 - t441 + t450) * t236 + (t235 * t443 + t451 - t481) * t235) * qJD(4) + t458 - t459) * t329 + (t454 + t455) * t353 / 0.2e1 - (t453 - t456 + t457) * t352 / 0.2e1 + ((t449 + t478) * t235 + (t448 + t477) * t236) * qJD(4) * t419;
t23 = t296 * t235 + (t414 + t279) * t244 + t315;
t50 = (-t162 + t462) * t244 + t266;
t51 = t346 - t425;
t250 = (-t23 * t410 + (t50 * (-t237 - t411) - t51 * t245) * t244 + t51 * (-t407 + (-rSges(6,1) - pkin(4)) * t246) * qJD(4)) * t236;
t234 = pkin(2) * t239;
t181 = t311 * qJD(4);
t155 = t215 * t236;
t153 = t215 * t235;
t152 = t214 * t235;
t139 = t163 * t244;
t135 = t161 * t244 + t346;
t115 = -t138 * t244 - t348;
t114 = -t244 * t393 - t349;
t66 = qJD(4) * t301 + qJD(1);
t49 = -t348 - t181 * t352 + (-t139 - t94 + t158) * t244;
t48 = t244 * t92 + (-t181 * t235 - t215 * t380) * qJD(4) + t315;
t42 = qJD(1) + (-t235 * t462 + t236 * t371) * qJD(4);
t25 = t264 * qJD(4);
t24 = -t348 + t296 * t236 + (t324 * t353 - t139 + t219 - t413) * t244;
t5 = ((t414 - t476) * t236 + (-t244 * t371 + t413) * t235) * qJD(4);
t1 = [m(5) * t25 + m(6) * t5; t251 + m(4) * (t115 * (-t160 - t418) + t114 * (t161 + t234) + (-t138 - t346 + t135) * t134) + (t24 * (t281 - t418) + t50 * (t275 - t346) + t23 * (t234 + t298) + t250 + (t50 - t266 - t347 + t444) * t51) * m(6) + (t49 * (t276 - t418) + t64 * (t339 - t346) + t48 * (t234 + t280) + t252 + (-t283 + t64 - t347 + t445) * t65) * m(5); t251 + (t23 * t298 + t24 * t281 + t250 + (-t279 + t444) * t51 + (-t425 + t275) * t50) * m(6) + (t276 * t49 + t280 * t48 + t252 + (t337 + t445) * t65 + (t339 - t428) * t64) * m(5) + (t114 * t161 - t115 * t160 - t134 * t138 - t135 * t393 - (-t134 * t161 - t135 * t160) * t244) * m(4); -(((t354 + t356) * t247 + (t355 + t357) * t246) * t244 + (((-t364 - t366) * t236 + (t363 + t365) * t235) * t247 + ((t368 + t370) * t236 + (t367 + t369) * t235) * t246) * qJD(4)) * t244 / 0.2e1 + ((t244 * t448 + t456) * t236 + (t244 * t449 + t455) * t235) * t419 + ((-t353 * t388 + t387) * t235 + (t273 + (-t424 * t236 + (t389 + t277) * t235) * qJD(4)) * t236 + (-t353 * t391 + t390) * t235 + (t274 + (-t423 * t236 + (t392 + t278) * t235) * qJD(4)) * t236) * t329 + ((-t352 * t389 - t387) * t236 + (t273 + (t277 * t235 + (t388 - t424) * t236) * qJD(4)) * t235 + (-t352 * t392 - t390) * t236 + (t274 + (t278 * t235 + (t391 - t423) * t236) * qJD(4)) * t235) * t326 + (t25 * t301 + t66 * t264 + t310 * t181 + ((-t244 * t65 - t49) * t236 + (t244 * t64 - t48) * t235) * t215 - (t153 * t64 - t155 * t65) * t244 - (t66 * (-t153 * t235 - t155 * t236) + t310 * t311) * qJD(4)) * m(5) + (t454 * t244 + ((t446 * t236 + t450 * t244) * t236 + (t439 * t235 + t451 * t244 + (-t440 + t447) * t236) * t235) * t438) * t422 + (t453 * t244 + ((t440 * t236 + t452 * t244) * t236 + (t447 * t235 + t479 * t244 + (-t439 + t446) * t236) * t235) * t438) * t421 + ((-t23 * t324 + t51 * t313 - t5 * t462 + t42 * t413 + (t214 * t50 - t371 * t42) * t244) * t235 + (-t24 * t324 + t50 * t313 + t5 * t371 + t42 * t414 + (-t324 * t51 - t42 * t462) * t244) * t236 - (t50 * t152 + t314 * t51) * t244 - ((t314 * t42 + t323 * t50) * t236 + (t51 * t323 + (-pkin(4) * t384 - t152) * t42) * t235) * qJD(4)) * m(6) + (t458 + t460) * t386 / 0.2e1 + (t457 + t461) * t380 / 0.2e1; 0.2e1 * (t23 * t421 + t24 * t422) * m(6);];
tauc = t1(:);
