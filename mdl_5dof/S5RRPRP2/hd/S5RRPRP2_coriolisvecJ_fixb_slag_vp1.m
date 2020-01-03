% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP2
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:40
% DurationCPUTime: 11.15s
% Computational Cost: add. (12297->445), mult. (9480->544), div. (0->0), fcn. (7254->8), ass. (0->267)
t257 = cos(qJ(4));
t255 = sin(qJ(4));
t419 = Icges(5,4) * t255;
t213 = Icges(5,2) * t257 + t419;
t249 = Icges(6,5) * t255;
t321 = Icges(6,3) * t257 - t249;
t508 = t213 + t321;
t418 = Icges(6,5) * t257;
t215 = Icges(6,1) * t255 - t418;
t250 = Icges(5,4) * t257;
t217 = Icges(5,1) * t255 + t250;
t512 = t215 + t217;
t254 = qJ(1) + qJ(2);
t246 = pkin(8) + t254;
t243 = sin(t246);
t244 = cos(t246);
t323 = Icges(6,1) * t257 + t249;
t128 = -Icges(6,4) * t244 + t243 * t323;
t403 = t243 * t255;
t201 = Icges(5,4) * t403;
t402 = t243 * t257;
t130 = Icges(5,1) * t402 - Icges(5,5) * t244 - t201;
t511 = t128 + t130;
t298 = t323 * t244;
t129 = Icges(6,4) * t243 + t298;
t218 = Icges(5,1) * t257 - t419;
t299 = t218 * t244;
t131 = Icges(5,5) * t243 + t299;
t504 = t129 + t131;
t208 = Icges(6,3) * t255 + t418;
t322 = -Icges(5,2) * t255 + t250;
t510 = t208 - t322;
t509 = t218 + t323;
t499 = -t508 * t255 + t512 * t257;
t399 = t244 * t257;
t200 = Icges(6,5) * t399;
t400 = t244 * t255;
t121 = Icges(6,6) * t243 + Icges(6,3) * t400 + t200;
t210 = Icges(5,5) * t257 - Icges(5,6) * t255;
t295 = t210 * t244;
t123 = Icges(5,3) * t243 + t295;
t212 = Icges(6,4) * t257 + Icges(6,6) * t255;
t296 = t212 * t244;
t125 = Icges(6,2) * t243 + t296;
t506 = t121 * t400 + t504 * t399 + (t123 + t125) * t243;
t124 = -Icges(6,2) * t244 + t212 * t243;
t112 = t243 * t124;
t120 = -Icges(6,6) * t244 + t208 * t243;
t122 = Icges(5,5) * t402 - Icges(5,6) * t403 - Icges(5,3) * t244;
t505 = -t120 * t400 - t243 * t122 - t511 * t399 - t112;
t503 = t510 * qJD(4);
t502 = t509 * qJD(4);
t209 = Icges(5,5) * t255 + Icges(5,6) * t257;
t211 = Icges(6,4) * t255 - Icges(6,6) * t257;
t501 = t209 + t211;
t500 = -t210 - t212;
t253 = qJD(1) + qJD(2);
t498 = ((Icges(6,4) + Icges(5,5)) * t253 - t512 * qJD(4)) * t255 - ((-Icges(5,6) + Icges(6,6)) * t253 + t508 * qJD(4)) * t257;
t126 = Icges(5,4) * t402 - Icges(5,2) * t403 - Icges(5,6) * t244;
t413 = t126 * t255;
t316 = -t130 * t257 + t413;
t414 = t124 * t244;
t319 = t120 * t255 + t128 * t257;
t455 = t243 * t319;
t50 = -t414 + t455;
t497 = -t122 * t244 - t243 * t316 + t50;
t472 = -t126 * t400 - t505;
t297 = t322 * t244;
t127 = Icges(5,6) * t243 + t297;
t471 = -t127 * t400 + t506;
t406 = t211 * t244;
t409 = t209 * t244;
t496 = t499 * t243 - t406 - t409;
t407 = t211 * t243;
t410 = t209 * t243;
t495 = t499 * t244 + t407 + t410;
t329 = -t121 * t403 + t125 * t244 - t129 * t402;
t103 = t131 * t402;
t336 = t123 * t244 - t103;
t53 = -t127 * t403 - t336;
t494 = -t329 + t53;
t225 = rSges(6,1) * t257 + rSges(6,3) * t255;
t493 = pkin(4) * t257 + qJ(5) * t255 + t225;
t479 = rSges(6,1) + pkin(4);
t470 = rSges(6,3) + qJ(5);
t490 = t502 * t257 + t503 * t255 + t501 * t253 + (-t255 * t512 - t508 * t257) * qJD(4);
t489 = (Icges(6,2) + Icges(5,3)) * t253 - t501 * qJD(4);
t412 = t127 * t255;
t488 = -t121 * t255 - t257 * t504 + t412;
t487 = t316 - t319;
t486 = qJD(4) * t500 + t253 * t499;
t485 = t495 * t253;
t401 = t244 * t253;
t237 = t244 * rSges(6,2);
t385 = t493 * t243 - t237;
t484 = t385 * t253;
t483 = ((t295 + t296 + t487) * t253 + t489 * t243) * t244;
t482 = (t471 * t243 - t472 * t244) * qJD(4);
t481 = (t494 * t243 - t497 * t244) * qJD(4);
t480 = t496 * t253;
t478 = t480 + t481;
t477 = t482 + t485;
t476 = t208 * t401 * t257 + t487 * qJD(4) + (-t297 * t257 + (-t298 - t299) * t255) * t253 - t498 * t243;
t404 = t243 * t253;
t475 = -t488 * qJD(4) + (-t509 * t255 + t510 * t257) * t404 + t498 * t244;
t474 = -t486 * t243 + t490 * t244;
t473 = t490 * t243 + t486 * t244;
t469 = (-t121 + t127) * t257 + t504 * t255;
t468 = (-t120 + t126) * t257 + t511 * t255;
t373 = t479 * t255 - t470 * t257;
t337 = t244 * t373;
t205 = rSges(5,1) * t399;
t452 = t243 * rSges(5,3) + t205;
t135 = -rSges(5,2) * t400 + t452;
t241 = t244 * pkin(3);
t172 = t243 * pkin(7) + t241;
t248 = cos(t254);
t245 = pkin(2) * t248;
t451 = t245 + t172;
t286 = t135 + t451;
t240 = t244 * pkin(7);
t171 = pkin(3) * t243 - t240;
t195 = pkin(7) * t401;
t467 = t253 * t171 + t195;
t466 = t243 * rSges(6,2) + pkin(4) * t399;
t465 = -t244 * rSges(4,1) - t245;
t462 = t414 + t506;
t368 = qJD(4) * t257;
t353 = t243 * t368;
t396 = t253 * t255;
t461 = t244 * t396 + t353;
t256 = sin(qJ(1));
t425 = pkin(1) * qJD(1);
t360 = t256 * t425;
t247 = sin(t254);
t175 = rSges(3,1) * t247 + rSges(3,2) * t248;
t411 = t175 * t253;
t138 = -t360 - t411;
t398 = t247 * t253;
t364 = pkin(2) * t398;
t460 = t364 - t360;
t458 = t489 * t244 + t488 * t253 + t404 * t500;
t384 = rSges(6,1) * t399 + t470 * t400 + t466;
t290 = t451 + t384;
t457 = 0.2e1 * qJD(4);
t366 = qJD(5) * t257;
t47 = -t366 + qJD(3) + (t385 * t243 + t384 * t244) * qJD(4);
t456 = qJD(4) * t47;
t454 = -rSges(4,2) * t243 - t465;
t369 = qJD(4) * t255;
t354 = t243 * t369;
t453 = t479 * t354;
t367 = qJD(5) * t255;
t196 = t244 * t367;
t351 = t244 * t368;
t450 = rSges(6,2) * t401 + t470 * t351 + t196;
t378 = rSges(5,2) * t403 + t244 * rSges(5,3);
t133 = rSges(5,1) * t402 - t378;
t117 = t253 * t133;
t359 = t243 * t396;
t300 = rSges(5,3) * t401 + (-t351 + t359) * rSges(5,2);
t352 = t244 * t369;
t449 = -rSges(5,1) * t352 + t117 + t300 + t467;
t222 = rSges(5,1) * t255 + rSges(5,2) * t257;
t371 = qJD(4) * t243;
t167 = t222 * t371;
t444 = t253 * t286 - t167;
t350 = t243 * t367;
t441 = t253 * t290 - t373 * t371 + t350;
t440 = -t364 + t450 + t467 + t484;
t387 = -Icges(5,2) * t402 + t130 - t201;
t391 = t217 * t243 + t126;
t439 = -t255 * t387 - t257 * t391;
t389 = -t321 * t243 + t128;
t393 = -t215 * t243 + t120;
t438 = -t255 * t389 + t257 * t393;
t252 = t253 ^ 2;
t427 = rSges(6,3) * t353 + t350 + t461 * qJ(5) - t453 + (t225 * t244 + t466) * t253;
t292 = -t253 * t402 - t352;
t428 = t479 * t292 - t470 * t359 + t450;
t5 = (t367 + (t428 + t484) * t244 + (-t384 * t253 + t427) * t243) * qJD(4);
t437 = m(6) * t5;
t433 = t253 / 0.2e1;
t431 = pkin(1) * t256;
t430 = pkin(2) * t247;
t429 = pkin(2) * t252;
t258 = cos(qJ(1));
t251 = t258 * pkin(1);
t426 = rSges(5,1) * t257;
t408 = t210 * t253;
t405 = t212 * t253;
t397 = t248 * t253;
t392 = -Icges(6,1) * t400 + t121 + t200;
t390 = -t217 * t244 - t127;
t388 = -t321 * t244 + t129;
t386 = -t213 * t244 + t131;
t383 = t373 * t243;
t381 = -qJD(4) * t493 + t366;
t377 = -t321 + t323;
t376 = t208 - t215;
t375 = -t213 + t218;
t374 = t217 + t322;
t370 = qJD(4) * t244;
t361 = t258 * t425;
t59 = t361 + t444;
t365 = t59 * t430;
t259 = qJD(1) ^ 2;
t363 = t259 * t431;
t362 = t259 * t251;
t357 = rSges(5,1) * t354 + t461 * rSges(5,2);
t355 = t222 * t370;
t347 = -pkin(3) - t426;
t346 = -t371 / 0.2e1;
t343 = t370 / 0.2e1;
t169 = rSges(4,1) * t243 + rSges(4,2) * t244;
t293 = -t169 - t430;
t340 = -t171 - t430;
t338 = t240 - t430;
t176 = t248 * rSges(3,1) - rSges(3,2) * t247;
t335 = -t122 + t412;
t331 = t237 + t338;
t163 = rSges(3,1) * t397 - rSges(3,2) * t398;
t328 = t366 + t381;
t326 = -rSges(5,2) * t255 + t426;
t294 = -t355 - t360;
t58 = (-t133 + t340) * t253 + t294;
t325 = -t243 * t59 - t244 * t58;
t314 = t133 * t243 + t135 * t244;
t311 = -pkin(2) * t397 + t453;
t302 = -t247 * t429 - t363;
t301 = -t248 * t429 - t362;
t291 = t253 * (-pkin(3) * t404 + t195) + t302;
t289 = -t255 * t388 + t257 * t392;
t288 = -t255 * t386 + t257 * t390;
t285 = -qJD(4) * t337 + t196 - t360;
t284 = (t255 * t376 + t257 * t377) * t253;
t283 = (-t255 * t374 + t257 * t375) * t253;
t276 = t243 * t347 + t338 + t378;
t88 = rSges(5,1) * t292 + t300;
t90 = t253 * t452 - t357;
t275 = (t88 + t117) * t244 + (-t135 * t253 + t90) * t243;
t114 = t253 * t293 - t360;
t115 = t253 * t454 + t361;
t263 = (t114 * t465 + t115 * t293) * t253;
t262 = (((t53 - t103 + (t123 + t413) * t244 + t505) * t244 + (t50 - t455 + t462) * t243) * qJD(4) + t485) * t343 + (qJD(4) * t499 + t255 * t502 - t257 * t503) * t253 + (t474 + t475) * t371 / 0.2e1 + (((t244 * t335 - t462 + t471) * t244 + (t243 * t335 - t112 + t329 + t336 + t472) * t243) * qJD(4) + t478 - t480) * t346 - (t473 - t476 + t477) * t370 / 0.2e1 + ((t468 + t496) * t243 + (t469 + t495) * t244) * qJD(4) * t433;
t261 = (t58 * (-t205 - t241 - t245) - t365 + (t58 * (-rSges(5,3) - pkin(7)) + t59 * t347) * t243) * t253;
t140 = t172 * t253;
t23 = t328 * t370 + (-t140 + (qJD(4) * t373 - t367) * t243 - t427) * t253 + t301;
t48 = (t340 - t385) * t253 + t285;
t49 = t361 + t441;
t260 = (-t49 * t479 * t369 + t48 * (-t255 * t470 - t257 * t479 - pkin(3)) * t253) * t244 + (-t23 * pkin(3) + (-t48 * qJD(5) - t23 * t470) * t255 + (-qJD(4) * t470 * t48 - t23 * t479) * t257 + (t48 * (-rSges(6,2) - pkin(7)) + t49 * (-pkin(3) - t493)) * t253) * t243;
t192 = rSges(4,2) * t404;
t184 = t326 * qJD(4);
t161 = t253 * t169;
t159 = t222 * t244;
t155 = t222 * t243;
t139 = t176 * t253 + t361;
t119 = -t163 * t253 - t362;
t118 = -t253 * t411 - t363;
t100 = -t253 * (rSges(4,1) * t401 - t192) + t301;
t99 = -t169 * t252 + t302;
t64 = qJD(4) * t314 + qJD(3);
t46 = -t184 * t370 + (-t140 - t90 + t167) * t253 + t301;
t45 = t253 * t88 + (-t184 * t243 - t222 * t401) * qJD(4) + t291;
t25 = t275 * qJD(4);
t22 = (t196 + t428) * t253 + (t243 * t328 - t253 * t337) * qJD(4) + t291;
t1 = [m(3) * (t119 * (-t175 - t431) + t118 * (t176 + t251) + (-t163 - t361 + t139) * t138) + t262 + (t23 * (t331 - t431) + t48 * (t311 - t361) + t22 * (t251 + t290) + t260 + (-t285 + t48 + t440 + t460) * t49) * m(6) + (t46 * (t276 - t431) + t58 * (t357 - t361) + t45 * (t251 + t286) + t261 + (-t294 + t58 + t449 + t460) * t59) * m(5) + (t100 * (t293 - t431) + t114 * (t192 - t361) + t99 * (t251 + t454) + t263 + (t114 + t161 + t364) * t115) * m(4); t262 + (t22 * t290 + t23 * t331 + t260 + (t430 * t253 + t370 * t373 - t196 + t440) * t49 + (t311 + t441) * t48) * m(6) + (t365 * t253 + t46 * t276 + t45 * t286 + t261 + (t355 + t449) * t59 + (t357 + t444) * t58) * m(5) + (t100 * t293 + t114 * t192 + t99 * t454 + t263 + t115 * t161 - (-t114 * t454 - t115 * t430) * t253) * m(4) + (t118 * t176 - t119 * t175 - t138 * t163 - t139 * t411 - (-t138 * t176 - t139 * t175) * t253) * m(3); m(5) * t25 + t437; -(((t374 - t376) * t257 + (t375 + t377) * t255) * t253 + (((-t387 - t389) * t244 + (t386 + t388) * t243) * t257 + ((t391 - t393) * t244 + (t390 + t392) * t243) * t255) * qJD(4)) * t253 / 0.2e1 + ((t469 * t253 + t476) * t244 + (t468 * t253 + t475) * t243) * t433 + ((-t371 * t406 + t405) * t243 + (t284 + (-t438 * t244 + (t407 + t289) * t243) * qJD(4)) * t244 + (-t371 * t409 + t408) * t243 + (t283 + (-t439 * t244 + (t410 + t288) * t243) * qJD(4)) * t244) * t346 + ((-t370 * t410 - t408) * t244 + (t283 + (t288 * t243 + (t409 - t439) * t244) * qJD(4)) * t243 + (-t370 * t407 - t405) * t244 + (t284 + (t289 * t243 + (t406 - t438) * t244) * qJD(4)) * t243) * t343 + ((-t23 * t373 + t48 * t381 + t5 * t384 + t47 * t428 + (-t373 * t49 + t385 * t47) * t253) * t244 + (-t22 * t373 + t49 * t381 + t5 * t385 + t47 * t427 + (t373 * t48 - t384 * t47) * t253) * t243 - (t255 * t47 + (t243 * t49 + t244 * t48) * t257) * qJD(5) - (-t337 * t49 + t383 * t48) * t253 - ((-t337 * t47 - t48 * t493) * t244 + (-t383 * t47 - t49 * t493) * t243) * qJD(4)) * m(6) + (-(t155 * t58 - t159 * t59) * t253 - (t64 * (-t155 * t243 - t159 * t244) + t325 * t326) * qJD(4) + t25 * t314 + t64 * t275 + t325 * t184 + ((-t253 * t59 - t46) * t244 + (t253 * t58 - t45) * t243) * t222) * m(5) + (t474 * t253 + (t471 * t401 + (t458 * t243 + t472 * t253 - t483) * t243) * t457) * t243 / 0.2e1 - (t473 * t253 + ((t494 * t253 + t483) * t244 + (-t458 * t244 + t497 * t253) * t243) * t457) * t244 / 0.2e1 + (t478 + t481) * t404 / 0.2e1 + (t477 + t482) * t401 / 0.2e1; -t257 * t437 + 0.2e1 * (m(6) * (t22 * t243 + t23 * t244 + t456) / 0.2e1 - m(6) * (t243 ^ 2 + t244 ^ 2) * t456 / 0.2e1) * t255;];
tauc = t1(:);
