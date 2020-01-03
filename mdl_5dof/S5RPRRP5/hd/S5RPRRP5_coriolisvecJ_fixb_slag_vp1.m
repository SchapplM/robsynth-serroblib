% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP5
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
% Datum: 2019-12-31 18:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:40:42
% EndTime: 2019-12-31 18:40:55
% DurationCPUTime: 10.50s
% Computational Cost: add. (12099->410), mult. (9380->517), div. (0->0), fcn. (7206->8), ass. (0->240)
t251 = cos(qJ(4));
t249 = sin(qJ(4));
t403 = Icges(5,4) * t249;
t209 = Icges(5,2) * t251 + t403;
t244 = Icges(6,5) * t249;
t313 = Icges(6,3) * t251 - t244;
t492 = t209 + t313;
t402 = Icges(6,5) * t251;
t211 = Icges(6,1) * t249 - t402;
t245 = Icges(5,4) * t251;
t213 = Icges(5,1) * t249 + t245;
t496 = t211 + t213;
t248 = qJ(1) + pkin(8);
t243 = qJ(3) + t248;
t239 = sin(t243);
t240 = cos(t243);
t315 = Icges(6,1) * t251 + t244;
t126 = -Icges(6,4) * t240 + t239 * t315;
t387 = t239 * t249;
t197 = Icges(5,4) * t387;
t386 = t239 * t251;
t128 = Icges(5,1) * t386 - Icges(5,5) * t240 - t197;
t495 = t126 + t128;
t289 = t315 * t240;
t127 = Icges(6,4) * t239 + t289;
t214 = Icges(5,1) * t251 - t403;
t290 = t214 * t240;
t129 = Icges(5,5) * t239 + t290;
t488 = t127 + t129;
t204 = Icges(6,3) * t249 + t402;
t314 = -Icges(5,2) * t249 + t245;
t494 = t204 - t314;
t493 = t214 + t315;
t483 = -t249 * t492 + t251 * t496;
t383 = t240 * t251;
t196 = Icges(6,5) * t383;
t384 = t240 * t249;
t119 = Icges(6,6) * t239 + Icges(6,3) * t384 + t196;
t206 = Icges(5,5) * t251 - Icges(5,6) * t249;
t286 = t206 * t240;
t121 = Icges(5,3) * t239 + t286;
t208 = Icges(6,4) * t251 + Icges(6,6) * t249;
t287 = t208 * t240;
t123 = Icges(6,2) * t239 + t287;
t490 = t119 * t384 + t488 * t383 + (t121 + t123) * t239;
t122 = -Icges(6,2) * t240 + t208 * t239;
t112 = t239 * t122;
t118 = -Icges(6,6) * t240 + t204 * t239;
t120 = Icges(5,5) * t386 - Icges(5,6) * t387 - Icges(5,3) * t240;
t489 = -t118 * t384 - t239 * t120 - t383 * t495 - t112;
t487 = t494 * qJD(4);
t486 = t493 * qJD(4);
t205 = Icges(5,5) * t249 + Icges(5,6) * t251;
t207 = Icges(6,4) * t249 - Icges(6,6) * t251;
t485 = t205 + t207;
t484 = -t206 - t208;
t247 = qJD(1) + qJD(3);
t482 = ((Icges(6,4) + Icges(5,5)) * t247 - t496 * qJD(4)) * t249 - ((-Icges(5,6) + Icges(6,6)) * t247 + t492 * qJD(4)) * t251;
t124 = Icges(5,4) * t386 - Icges(5,2) * t387 - Icges(5,6) * t240;
t396 = t124 * t249;
t308 = -t128 * t251 + t396;
t397 = t122 * t240;
t311 = t118 * t249 + t126 * t251;
t440 = t239 * t311;
t50 = -t397 + t440;
t481 = -t120 * t240 - t239 * t308 + t50;
t456 = -t124 * t384 - t489;
t288 = t314 * t240;
t125 = Icges(5,6) * t239 + t288;
t455 = -t125 * t384 + t490;
t390 = t207 * t240;
t393 = t205 * t240;
t480 = t239 * t483 - t390 - t393;
t391 = t207 * t239;
t394 = t205 * t239;
t479 = t240 * t483 + t391 + t394;
t322 = -t119 * t387 + t123 * t240 - t127 * t386;
t103 = t129 * t386;
t328 = t121 * t240 - t103;
t53 = -t125 * t387 - t328;
t478 = -t322 + t53;
t220 = rSges(6,1) * t251 + rSges(6,3) * t249;
t477 = pkin(4) * t251 + qJ(5) * t249 + t220;
t463 = rSges(6,1) + pkin(4);
t454 = rSges(6,3) + qJ(5);
t474 = t486 * t251 + t487 * t249 + t485 * t247 + (-t249 * t496 - t251 * t492) * qJD(4);
t473 = (Icges(6,2) + Icges(5,3)) * t247 - t485 * qJD(4);
t395 = t125 * t249;
t472 = -t119 * t249 - t251 * t488 + t395;
t471 = t308 - t311;
t470 = qJD(4) * t484 + t483 * t247;
t469 = t479 * t247;
t385 = t240 * t247;
t232 = t240 * rSges(6,2);
t370 = t239 * t477 - t232;
t468 = t247 * t370;
t467 = ((t286 + t287 + t471) * t247 + t473 * t239) * t240;
t466 = (t239 * t455 - t240 * t456) * qJD(4);
t465 = (t239 * t478 - t240 * t481) * qJD(4);
t464 = t480 * t247;
t253 = qJD(1) ^ 2;
t462 = t464 + t465;
t461 = t466 + t469;
t460 = t204 * t385 * t251 + t471 * qJD(4) + (-t288 * t251 + (-t289 - t290) * t249) * t247 - t482 * t239;
t388 = t239 * t247;
t459 = -t472 * qJD(4) + (-t249 * t493 + t251 * t494) * t388 + t482 * t240;
t458 = -t239 * t470 + t240 * t474;
t457 = t239 * t474 + t240 * t470;
t453 = (-t119 + t125) * t251 + t488 * t249;
t452 = (-t118 + t124) * t251 + t495 * t249;
t358 = t249 * t463 - t251 * t454;
t329 = t240 * t358;
t451 = rSges(6,2) * t239 + pkin(4) * t383;
t448 = t397 + t490;
t351 = qJD(4) * t251;
t343 = t239 * t351;
t381 = t247 * t249;
t446 = t240 * t381 + t343;
t444 = t473 * t240 + t472 * t247 + t388 * t484;
t442 = 0.2e1 * qJD(4);
t349 = qJD(5) * t251;
t369 = rSges(6,1) * t383 + t384 * t454 + t451;
t47 = -t349 + qJD(2) + (t239 * t370 + t240 * t369) * qJD(4);
t441 = qJD(4) * t47;
t363 = rSges(5,2) * t387 + rSges(5,3) * t240;
t131 = rSges(5,1) * t386 - t363;
t117 = t247 * t131;
t235 = t240 * pkin(7);
t168 = pkin(3) * t239 - t235;
t161 = t247 * t168;
t439 = -t117 - t161;
t169 = pkin(3) * t240 + pkin(7) * t239;
t438 = t169 + t369;
t352 = qJD(4) * t249;
t344 = t239 * t352;
t437 = t463 * t344;
t436 = rSges(5,1) * t383 + rSges(5,3) * t239;
t435 = pkin(2) * cos(t248) + cos(qJ(1)) * pkin(1);
t350 = qJD(5) * t249;
t192 = t240 * t350;
t340 = t240 * t351;
t433 = rSges(6,2) * t385 + t340 * t454 + t192;
t432 = -t161 - t468;
t217 = rSges(5,1) * t249 + rSges(5,2) * t251;
t354 = qJD(4) * t239;
t164 = t217 * t354;
t133 = -rSges(5,2) * t384 + t436;
t282 = t133 + t169;
t426 = -t247 * t282 + t164;
t372 = -Icges(5,2) * t386 + t128 - t197;
t376 = t213 * t239 + t124;
t423 = -t249 * t372 - t251 * t376;
t374 = -t239 * t313 + t126;
t378 = -t211 * t239 + t118;
t422 = -t249 * t374 + t251 * t378;
t339 = t239 * t350;
t413 = rSges(6,3) * t343 + t339 + t446 * qJ(5) - t437 + (t220 * t240 + t451) * t247;
t341 = t240 * t352;
t283 = -t247 * t386 - t341;
t348 = t239 * t381;
t414 = t283 * t463 - t348 * t454 + t433;
t5 = (t350 + (t414 + t468) * t240 + (-t247 * t369 + t413) * t239) * qJD(4);
t421 = m(6) * t5;
t417 = t247 / 0.2e1;
t412 = rSges(5,1) * t251;
t233 = t240 * rSges(4,1);
t324 = -pkin(2) * sin(t248) - sin(qJ(1)) * pkin(1);
t299 = t324 * qJD(1);
t353 = qJD(4) * t240;
t342 = t217 * t353;
t269 = t299 - t342;
t58 = (-t131 - t168) * t247 + t269;
t410 = t240 * t58;
t297 = -qJD(4) * t329 + t192;
t263 = t299 + t297;
t48 = (-t168 - t370) * t247 + t263;
t409 = t247 * t48;
t167 = -rSges(4,2) * t239 + t233;
t298 = t435 * qJD(1);
t115 = t167 * t247 + t298;
t166 = rSges(4,1) * t239 + rSges(4,2) * t240;
t399 = t115 * t166;
t392 = t206 * t247;
t389 = t208 * t247;
t382 = t247 * t166;
t377 = -Icges(6,1) * t384 + t119 + t196;
t375 = -t213 * t240 - t125;
t373 = -t240 * t313 + t127;
t371 = -t209 * t240 + t129;
t368 = t358 * t239;
t366 = -qJD(4) * t477 + t349;
t362 = -t313 + t315;
t361 = t204 - t211;
t360 = -t209 + t214;
t359 = t213 + t314;
t356 = t232 + t235;
t345 = rSges(5,1) * t344 + rSges(5,2) * t446;
t336 = -pkin(3) - t412;
t335 = -t354 / 0.2e1;
t332 = t353 / 0.2e1;
t327 = -t120 + t395;
t321 = t349 + t366;
t318 = -rSges(5,2) * t249 + t412;
t59 = t298 - t426;
t317 = -t239 * t59 - t410;
t306 = t131 * t239 + t133 * t240;
t302 = -t354 * t358 + t339;
t301 = t324 * t253;
t300 = t435 * t253;
t291 = rSges(5,3) * t385 + (-t340 + t348) * rSges(5,2);
t191 = pkin(7) * t385;
t281 = t247 * (-pkin(3) * t388 + t191) + t301;
t280 = -t249 * t373 + t251 * t377;
t279 = -t249 * t371 + t251 * t375;
t278 = t239 * t336 + t235 + t363;
t114 = t299 - t382;
t277 = (t249 * t361 + t251 * t362) * t247;
t276 = (-t249 * t359 + t251 * t360) * t247;
t88 = rSges(5,1) * t283 + t291;
t90 = t247 * t436 - t345;
t268 = (t88 + t117) * t240 + (-t133 * t247 + t90) * t239;
t256 = (((t53 - t103 + (t121 + t396) * t240 + t489) * t240 + (t50 - t440 + t448) * t239) * qJD(4) + t469) * t332 + (qJD(4) * t483 + t249 * t486 - t251 * t487) * t247 + (t458 + t459) * t354 / 0.2e1 + (((t240 * t327 - t448 + t455) * t240 + (t239 * t327 - t112 + t322 + t328 + t456) * t239) * qJD(4) + t462 - t464) * t335 - (t457 - t460 + t461) * t353 / 0.2e1 + ((t452 + t480) * t239 + (t453 + t479) * t240) * qJD(4) * t417;
t255 = t58 * t345 + t59 * (-rSges(5,1) * t341 + t191 + t291) + (t336 * t410 + (t58 * (-rSges(5,3) - pkin(7)) + t59 * t336) * t239) * t247;
t138 = t169 * t247;
t23 = -t300 + t321 * t353 + (-t138 + (qJD(4) * t358 - t350) * t239 - t413) * t247;
t49 = t247 * t438 + t298 + t302;
t254 = t48 * t437 + t49 * (t191 + t433) + (-t49 * t463 * t352 + (-t249 * t454 - t251 * t463 - pkin(3)) * t409) * t240 + (-t23 * pkin(3) + (-qJD(5) * t48 - t23 * t454) * t249 + (-qJD(4) * t454 * t48 - t23 * t463) * t251 + (t48 * (-rSges(6,2) - pkin(7)) + t49 * (-pkin(3) - t477)) * t247) * t239;
t188 = rSges(4,2) * t388;
t180 = t318 * qJD(4);
t159 = t217 * t240;
t155 = t217 * t239;
t137 = rSges(4,1) * t385 - t188;
t100 = -t137 * t247 - t300;
t99 = -t247 * t382 + t301;
t64 = qJD(4) * t306 + qJD(2);
t46 = -t180 * t353 - t300 + (-t138 - t90 + t164) * t247;
t45 = t247 * t88 + (-t180 * t239 - t217 * t385) * qJD(4) + t281;
t25 = t268 * qJD(4);
t22 = (t192 + t414) * t247 + (t239 * t321 - t247 * t329) * qJD(4) + t281;
t1 = [m(4) * (t100 * (-t166 + t324) + t114 * t188 + t99 * (t167 + t435) + (-t114 * t233 - t399) * t247 + (-t114 * t435 + t115 * t324) * qJD(1)) + t256 + (-(-t48 + t263 + t432) * t49 + t23 * (t324 + t356) + t22 * (t438 + t435) + (t324 * t49 - t435 * t48) * qJD(1) + t254) * m(6) + (t46 * (t278 + t324) + t45 * (t282 + t435) + (t324 * t59 - t435 * t58) * qJD(1) + t255 - (-t58 + t269 + t439) * t59) * m(5); m(5) * t25 + t421; t256 + (t48 * t302 - t49 * (t297 + t432) + t23 * t356 + t254 + (t409 + t22) * t438) * m(6) + (t46 * t278 + t45 * t282 + t255 - t58 * t426 - t59 * (-t342 + t439)) * m(5) + (-t100 * t166 - t114 * t137 - t115 * t382 + t167 * t99 - (-t114 * t167 - t399) * t247) * m(4); -(((t359 - t361) * t251 + (t360 + t362) * t249) * t247 + (((-t372 - t374) * t240 + (t371 + t373) * t239) * t251 + ((t376 - t378) * t240 + (t375 + t377) * t239) * t249) * qJD(4)) * t247 / 0.2e1 + ((t247 * t453 + t460) * t240 + (t247 * t452 + t459) * t239) * t417 + ((-t354 * t390 + t389) * t239 + (t277 + (-t422 * t240 + (t391 + t280) * t239) * qJD(4)) * t240 + (-t354 * t393 + t392) * t239 + (t276 + (-t423 * t240 + (t394 + t279) * t239) * qJD(4)) * t240) * t335 + ((-t353 * t391 - t389) * t240 + (t277 + (t280 * t239 + (t390 - t422) * t240) * qJD(4)) * t239 + (-t353 * t394 - t392) * t240 + (t276 + (t279 * t239 + (t393 - t423) * t240) * qJD(4)) * t239) * t332 + (-(t249 * t47 + (t239 * t49 + t240 * t48) * t251) * qJD(5) - (-t329 * t49 + t368 * t48) * t247 - ((-t329 * t47 - t477 * t48) * t240 + (-t368 * t47 - t477 * t49) * t239) * qJD(4) + (-t23 * t358 + t48 * t366 + t5 * t369 + t47 * t414 + (-t358 * t49 + t370 * t47) * t247) * t240 + (-t22 * t358 + t49 * t366 + t5 * t370 + t47 * t413 + (t358 * t48 - t369 * t47) * t247) * t239) * m(6) + (-(t155 * t58 - t159 * t59) * t247 - (t64 * (-t155 * t239 - t159 * t240) + t317 * t318) * qJD(4) + t25 * t306 + t64 * t268 + t317 * t180 + ((-t247 * t59 - t46) * t240 + (t247 * t58 - t45) * t239) * t217) * m(5) + (t458 * t247 + (t455 * t385 + (t444 * t239 + t456 * t247 - t467) * t239) * t442) * t239 / 0.2e1 - (t457 * t247 + ((t247 * t478 + t467) * t240 + (-t444 * t240 + t247 * t481) * t239) * t442) * t240 / 0.2e1 + (t462 + t465) * t388 / 0.2e1 + (t461 + t466) * t385 / 0.2e1; -t251 * t421 + 0.2e1 * (m(6) * (t22 * t239 + t23 * t240 + t441) / 0.2e1 - m(6) * (t239 ^ 2 + t240 ^ 2) * t441 / 0.2e1) * t249;];
tauc = t1(:);
