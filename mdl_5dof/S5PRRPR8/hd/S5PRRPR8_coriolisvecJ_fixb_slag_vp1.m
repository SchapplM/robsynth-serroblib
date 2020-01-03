% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPR8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR8_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:19
% EndTime: 2019-12-31 17:42:36
% DurationCPUTime: 14.35s
% Computational Cost: add. (17569->557), mult. (17847->845), div. (0->0), fcn. (16613->10), ass. (0->311)
t291 = qJ(2) + qJ(3);
t283 = pkin(9) + t291;
t277 = sin(t283);
t278 = cos(t283);
t284 = sin(t291);
t285 = cos(t291);
t551 = -Icges(4,5) * t284 - Icges(5,5) * t277 - Icges(4,6) * t285 - Icges(5,6) * t278;
t295 = sin(qJ(2));
t297 = cos(qJ(2));
t299 = qJD(2) ^ 2;
t292 = sin(pkin(8));
t288 = t292 ^ 2;
t293 = cos(pkin(8));
t289 = t293 ^ 2;
t440 = t288 + t289;
t532 = t299 * t440;
t550 = m(3) * (rSges(3,1) * t295 + rSges(3,2) * t297) * t532;
t549 = -Icges(4,1) + Icges(4,2);
t548 = -Icges(5,1) + Icges(5,2);
t547 = t551 * t292;
t546 = t551 * t293;
t294 = sin(qJ(5));
t454 = t293 * t294;
t296 = cos(qJ(5));
t455 = t292 * t296;
t243 = -t278 * t454 + t455;
t225 = Icges(6,4) * t243;
t453 = t293 * t296;
t456 = t292 * t294;
t244 = t278 * t453 + t456;
t460 = t277 * t293;
t100 = Icges(6,1) * t244 + Icges(6,5) * t460 + t225;
t290 = qJD(2) + qJD(3);
t462 = t277 * t290;
t428 = t294 * t462;
t134 = -qJD(5) * t244 + t293 * t428;
t427 = t296 * t462;
t135 = qJD(5) * t243 - t293 * t427;
t459 = t278 * t290;
t457 = t290 * t293;
t425 = t278 * t457;
t74 = Icges(6,5) * t135 + Icges(6,6) * t134 + Icges(6,3) * t425;
t96 = Icges(6,5) * t244 + Icges(6,6) * t243 + Icges(6,3) * t460;
t348 = t277 * t74 + t459 * t96;
t76 = Icges(6,4) * t135 + Icges(6,2) * t134 + Icges(6,6) * t425;
t78 = Icges(6,1) * t135 + Icges(6,4) * t134 + Icges(6,5) * t425;
t469 = Icges(6,4) * t244;
t98 = Icges(6,2) * t243 + Icges(6,6) * t460 + t469;
t20 = t100 * t135 + t134 * t98 + t243 * t76 + t244 * t78 + t293 * t348;
t282 = qJD(2) * t292;
t272 = qJD(3) * t292 + t282;
t437 = qJD(5) * t277;
t215 = t293 * t437 + t272;
t496 = t215 / 0.2e1;
t545 = t20 * t496;
t241 = -t278 * t456 - t453;
t242 = t278 * t455 - t454;
t461 = t277 * t292;
t101 = rSges(6,1) * t242 + rSges(6,2) * t241 + rSges(6,3) * t461;
t394 = rSges(6,1) * t296 - rSges(6,2) * t294;
t169 = -rSges(6,3) * t278 + t277 * t394;
t434 = qJD(5) * t292;
t216 = t277 * t434 - t457;
t436 = qJD(5) * t278;
t544 = -t101 * t436 - t169 * t216;
t494 = t216 / 0.2e1;
t542 = t547 * t290;
t541 = t546 * t290;
t536 = 0.2e1 * Icges(4,4) * t284 + t285 * t549;
t517 = -Icges(4,5) * t292 + t536 * t293;
t535 = -0.2e1 * Icges(4,4) * t285 + t284 * t549;
t519 = -Icges(4,6) * t292 + t535 * t293;
t538 = 0.2e1 * Icges(5,4) * t277 + t278 * t548;
t521 = -Icges(5,5) * t292 + t538 * t293;
t537 = -0.2e1 * Icges(5,4) * t278 + t277 * t548;
t523 = -Icges(5,6) * t292 + t537 * t293;
t540 = (t277 * t521 + t278 * t523 + t284 * t517 + t285 * t519) * t290;
t518 = Icges(4,5) * t293 + t292 * t536;
t520 = Icges(4,6) * t293 + t292 * t535;
t522 = Icges(5,5) * t293 + t292 * t538;
t524 = Icges(5,6) * t293 + t292 * t537;
t539 = (t277 * t522 + t278 * t524 + t284 * t518 + t285 * t520) * t290;
t132 = -qJD(5) * t242 + t292 * t428;
t133 = qJD(5) * t241 - t292 * t427;
t458 = t290 * t292;
t426 = t278 * t458;
t73 = Icges(6,5) * t133 + Icges(6,6) * t132 + Icges(6,3) * t426;
t95 = Icges(6,5) * t242 + Icges(6,6) * t241 + Icges(6,3) * t461;
t349 = t277 * t73 + t459 * t95;
t75 = Icges(6,4) * t133 + Icges(6,2) * t132 + Icges(6,6) * t426;
t77 = Icges(6,1) * t133 + Icges(6,4) * t132 + Icges(6,5) * t426;
t470 = Icges(6,4) * t242;
t97 = Icges(6,2) * t241 + Icges(6,6) * t461 + t470;
t224 = Icges(6,4) * t241;
t99 = Icges(6,1) * t242 + Icges(6,5) * t461 + t224;
t19 = t134 * t97 + t135 * t99 + t243 * t75 + t244 * t77 + t293 * t349;
t467 = Icges(6,4) * t296;
t371 = -Icges(6,2) * t294 + t467;
t156 = -Icges(6,6) * t278 + t277 * t371;
t468 = Icges(6,4) * t294;
t379 = Icges(6,1) * t296 - t468;
t158 = -Icges(6,5) * t278 + t277 * t379;
t366 = Icges(6,5) * t296 - Icges(6,6) * t294;
t154 = -Icges(6,3) * t278 + t277 * t366;
t334 = t366 * t278;
t365 = -Icges(6,5) * t294 - Icges(6,6) * t296;
t83 = t290 * t334 + (Icges(6,3) * t290 + qJD(5) * t365) * t277;
t346 = t154 * t459 + t277 * t83;
t48 = t243 * t97 + t244 * t99 + t460 * t95;
t479 = t292 * t48;
t49 = t100 * t244 + t243 * t98 + t460 * t96;
t389 = t293 * t49 + t479;
t335 = t371 * t278;
t370 = -Icges(6,2) * t296 - t468;
t84 = t290 * t335 + (Icges(6,6) * t290 + qJD(5) * t370) * t277;
t336 = t379 * t278;
t378 = -Icges(6,1) * t294 - t467;
t85 = t290 * t336 + (Icges(6,5) * t290 + qJD(5) * t378) * t277;
t318 = (-t134 * t156 - t135 * t158 - t243 * t84 - t244 * t85 + t290 * t389 - t293 * t346) * t278;
t64 = t154 * t460 + t156 * t243 + t158 * t244;
t534 = t19 * t494 + t545 + (t462 * t64 + t318) * qJD(5) / 0.2e1;
t533 = t547 * t457 / 0.2e1 - t546 * t272 / 0.2e1;
t531 = (-t539 / 0.2e1 + t541 / 0.2e1) * t293 + (-t542 / 0.2e1 - t540 / 0.2e1) * t292;
t256 = rSges(5,1) * t277 + rSges(5,2) * t278;
t341 = t256 * t292;
t171 = t290 * t341;
t214 = t256 * t293;
t172 = t290 * t214;
t525 = -t171 * t292 - t293 * t172 + t272 * t341;
t516 = (t272 * t519 - t457 * t520) * t285 + (t272 * t517 - t457 * t518) * t284 + (t272 * t523 - t457 * t524) * t278 + (t272 * t521 - t457 * t522) * t277;
t486 = pkin(3) * t284;
t487 = pkin(2) * t295;
t271 = -t486 - t487;
t415 = t271 + t487;
t484 = pkin(3) * t290;
t395 = t284 * t484;
t438 = qJD(4) * t293;
t141 = -t292 * t395 - t438;
t280 = qJD(4) * t292;
t142 = -t293 * t395 + t280;
t450 = t141 * t292 + t293 * t142;
t515 = -t272 * t292 * t415 + t450;
t130 = t169 * t292;
t258 = pkin(4) * t277 - pkin(7) * t278;
t344 = t258 * t292;
t173 = t290 * t344;
t218 = t258 * t293;
t174 = t290 * t218;
t79 = rSges(6,1) * t133 + rSges(6,2) * t132 + rSges(6,3) * t426;
t80 = rSges(6,1) * t135 + rSges(6,2) * t134 + rSges(6,3) * t425;
t513 = t215 * t130 + t272 * t344 + (-t173 + t79) * t292 + (-t174 + t80 + t544) * t293;
t475 = Icges(3,4) * t297;
t476 = Icges(3,4) * t295;
t512 = (-t295 * (-Icges(3,2) * t297 - t476) + t297 * (-Icges(3,1) * t295 - t475)) * qJD(2);
t46 = t241 * t97 + t242 * t99 + t461 * t95;
t47 = t100 * t242 + t241 * t98 + t461 * t96;
t63 = t154 * t461 + t156 * t241 + t158 * t242;
t511 = t293 * (t215 * t49 + t216 * t48 - t436 * t64) + t292 * (t215 * t47 + t216 * t46 - t436 * t63);
t267 = rSges(4,1) * t284 + rSges(4,2) * t285;
t342 = t292 * t267;
t193 = t290 * t342;
t245 = t267 * t293;
t194 = t290 * t245;
t268 = rSges(4,1) * t285 - rSges(4,2) * t284;
t203 = -rSges(4,3) * t293 + t268 * t292;
t204 = rSges(4,3) * t292 + t268 * t293;
t286 = t297 * pkin(2);
t197 = -pkin(6) * t293 + t286 * t292;
t198 = pkin(6) * t292 + t286 * t293;
t439 = qJD(2) * t293;
t420 = t197 * t282 + t198 * t439 + qJD(1);
t510 = (-t193 * t292 - t293 * t194 + t245 * t457 + t272 * t342) * (t203 * t272 + t204 * t457 + t420);
t257 = rSges(5,1) * t278 - rSges(5,2) * t277;
t485 = pkin(3) * t285;
t507 = t272 * (-t257 - t485);
t506 = t272 * t486;
t254 = t457 * t486;
t502 = t293 * t254 + t450;
t386 = t100 * t296 - t294 * t98;
t387 = -t294 * t97 + t296 * t99;
t499 = (t154 * t293 + t386) * t215 - (-t154 * t292 - t387) * t216;
t305 = t215 * (-Icges(6,2) * t244 + t100 + t225) + t216 * (-Icges(6,2) * t242 + t224 + t99) - t436 * (t277 * t370 + t158);
t497 = -t215 / 0.2e1;
t495 = -t216 / 0.2e1;
t480 = pkin(2) * qJD(2);
t478 = t293 * t47;
t463 = t154 * t278;
t122 = -qJ(4) * t293 + t292 * t485;
t123 = qJ(4) * t292 + t293 * t485;
t451 = t122 * t292 + t293 * t123;
t447 = t197 * t292 + t293 * t198;
t446 = t203 * t292 + t293 * t204;
t255 = t457 * t485;
t443 = -t257 * t457 - t255;
t435 = qJD(5) * t290;
t433 = t299 * t286;
t432 = t285 * t484;
t431 = t295 * t480;
t430 = t297 * t480;
t102 = rSges(6,1) * t244 + rSges(6,2) * t243 + rSges(6,3) * t460;
t423 = t102 * t436;
t422 = t278 * t434;
t419 = -t436 / 0.2e1;
t418 = t436 / 0.2e1;
t417 = t435 / 0.2e1;
t416 = -t267 - t487;
t414 = -t256 - t486;
t412 = -t258 - t486;
t410 = t292 * t433;
t409 = t293 * t433;
t408 = t292 * t431;
t407 = t292 * t430;
t406 = t293 * t431;
t405 = t293 * t430;
t185 = -rSges(5,3) * t293 + t257 * t292;
t186 = rSges(5,3) * t292 + t257 * t293;
t404 = t185 * t292 + t293 * t186 + t451;
t403 = t122 * t272 + t420;
t401 = t277 * t417;
t400 = t278 * t417;
t205 = t257 * t290;
t399 = -t205 - t432;
t259 = pkin(4) * t278 + pkin(7) * t277;
t206 = t259 * t290;
t398 = -t206 - t432;
t397 = -t169 + t412;
t231 = t268 * t290;
t396 = -t231 - t430;
t393 = -rSges(6,1) * t294 - rSges(6,2) * t296;
t392 = t215 * t96 + t216 * t95;
t340 = t394 * t278;
t88 = t290 * t340 + (rSges(6,3) * t290 + qJD(5) * t393) * t277;
t41 = -t410 - t215 * t88 + t398 * t272 + (t102 * t462 + (-t169 * t457 - t80) * t278) * qJD(5);
t354 = -t255 * t290 - t409;
t42 = -t206 * t457 + t216 * t88 + (-t101 * t462 + (t169 * t458 + t79) * t278) * qJD(5) + t354;
t391 = t292 * t42 - t293 * t41;
t390 = t292 * t46 + t478;
t55 = t277 * t387 - t278 * t95;
t56 = t277 * t386 - t278 * t96;
t388 = t292 * t55 + t56 * t293;
t385 = Icges(3,1) * t297 - t476;
t377 = -Icges(3,2) * t295 + t475;
t369 = -Icges(3,5) * t295 - Icges(3,6) * t297;
t364 = t101 * t293 - t102 * t292;
t363 = -t156 * t294 + t158 * t296;
t358 = t398 - t88;
t357 = t487 * t532;
t356 = t292 * t400;
t355 = t293 * t400;
t351 = -t256 + t271;
t217 = t259 * t292;
t219 = t259 * t293;
t350 = t451 + (t102 + t219) * t293 + (t101 + t217) * t292;
t343 = Icges(6,3) * t277 + t334 - t363;
t339 = -t254 + t280 - t406;
t337 = -t169 - t258 + t271;
t333 = -t430 - t432;
t330 = qJD(2) * t369;
t329 = t141 * t272 - t357;
t327 = -t408 - t438;
t326 = -t205 + t333;
t323 = -t206 + t333 - t88;
t170 = rSges(6,3) * t277 + t340;
t322 = -t101 * t437 - t130 * t436 + t169 * t422 + t170 * t216 - t259 * t457 - t255;
t321 = -(Icges(6,5) * t241 - Icges(6,6) * t242) * t216 - (Icges(6,5) * t243 - Icges(6,6) * t244) * t215 + t365 * t277 * t436;
t320 = (t290 * t388 - (t290 * t363 - t83) * t278 - (t154 * t290 - t294 * t84 + t296 * t85 + (-t156 * t296 - t158 * t294) * qJD(5)) * t277) * t278;
t319 = (-t132 * t156 - t133 * t158 - t241 * t84 - t242 * t85 + t290 * t390 - t292 * t346) * t278;
t313 = t277 * t321;
t308 = (-(-Icges(3,6) * t293 + t292 * t377) * t297 - (-Icges(3,5) * t293 + t292 * t385) * t295) * qJD(2) + t512 * t292;
t307 = (-(Icges(3,6) * t292 + t293 * t377) * t297 - (Icges(3,5) * t292 + t293 * t385) * t295) * qJD(2) + t512 * t293;
t306 = (Icges(6,1) * t243 - t469 - t98) * t215 + (Icges(6,1) * t241 - t470 - t97) * t216 - (t277 * t378 - t156) * t436;
t304 = -t170 * t215 + t102 * t437 + (-t259 - t485) * t272;
t11 = (t290 * t387 - t73) * t278 + (t290 * t95 - t294 * t75 + t296 * t77 + (-t294 * t99 - t296 * t97) * qJD(5)) * t277;
t12 = (t290 * t386 - t74) * t278 + (t290 * t96 - t294 * t76 + t296 * t78 + (-t100 * t294 - t296 * t98) * qJD(5)) * t277;
t126 = t156 * t292;
t127 = t156 * t293;
t128 = t158 * t292;
t129 = t158 * t293;
t157 = Icges(6,6) * t277 + t335;
t159 = Icges(6,5) * t277 + t336;
t17 = t132 * t97 + t133 * t99 + t241 * t75 + t242 * t77 + t292 * t349;
t18 = t100 * t133 + t132 * t98 + t241 * t76 + t242 * t78 + t292 * t348;
t70 = t277 * t363 - t463;
t23 = t215 * t56 + t216 * t55 - t436 * t70;
t3 = t17 * t216 + t18 * t215 + (t462 * t63 + t319) * qJD(5);
t300 = (-t343 * t436 - t499) * t277;
t301 = ((-t127 * t241 - t129 * t242) * t215 + (-t126 * t241 - t128 * t242) * t216 + (t63 * t277 + (-t157 * t241 - t159 * t242 + t478) * t278) * qJD(5)) * t495 + ((-t127 * t243 - t129 * t244) * t215 + (-t126 * t243 - t128 * t244) * t216 + (t64 * t277 + (-t157 * t243 - t159 * t244 + t479) * t278) * qJD(5)) * t497 - t23 * t437 / 0.2e1 + (((t127 * t294 - t129 * t296 + t96) * t215 + (t126 * t294 - t128 * t296 + t95) * t216 + t70 * qJD(5)) * t277 + ((t343 * t278 + (t157 * t294 - t159 * t296 - t154) * t277 + t388) * qJD(5) + t499) * t278) * t418 + t511 * t419 + (t18 * t494 + (((t46 - t463) * qJD(5) + t392) * t278 + t300) * t495 + t545 + t56 * t401 + t49 * t355 + t47 * t356 + t534 + t12 * t419 + (t516 / 0.2e1 + t531) * t457) * t292 + (-t17 * t494 - t19 * t496 + (((t49 - t463) * qJD(5) + t392) * t278 + t300) * t497 - t55 * t401 - t48 * t355 - t46 * t356 - t3 / 0.2e1 - t11 * t419 + (t292 * t539 - t542 * t293 + t533) * t457) * t293 + ((-t516 / 0.2e1 + t531) * t293 + (t292 * t541 + t540 * t293 + t533) * t292) * t272;
t262 = t369 * t293;
t261 = t369 * t292;
t247 = t293 * t330;
t246 = t292 * t330;
t223 = t393 * t277;
t213 = t415 * t293;
t139 = -t267 * t457 - t406;
t138 = -t267 * t272 - t408;
t120 = -t231 * t457 - t409;
t119 = -t231 * t272 - t410;
t116 = rSges(6,1) * t243 - rSges(6,2) * t244;
t115 = rSges(6,1) * t241 - rSges(6,2) * t242;
t92 = -t256 * t457 + t339;
t91 = t272 * t414 + t327;
t87 = -t205 * t457 + t354;
t86 = t272 * t399 - t410;
t81 = -t193 * t272 - t194 * t457 - t357;
t58 = -t258 * t457 + t339 - t544;
t57 = -t169 * t215 + t272 * t412 + t327 - t423;
t54 = -t171 * t272 - (-t142 + t172) * t457 + t329;
t45 = t185 * t272 - (-t123 - t186) * t457 + t403;
t34 = t101 * t215 - t102 * t216 + t217 * t272 - (-t123 - t219) * t457 + t403;
t21 = -t173 * t272 + t215 * t79 - t216 * t80 - (-t142 + t174) * t457 + t364 * t278 * t435 + t329;
t1 = [m(4) * t81 + m(5) * t54 + m(6) * t21 - t550; -(t262 * qJD(2) * t288 - t292 * t261 * t439) * t282 / 0.2e1 - (t292 * (-t247 * t293 + t292 * t307) - t293 * (-t246 * t293 + t292 * t308)) * t439 + (t292 * (t247 * t292 + t293 * t307) - t293 * (t246 * t292 + t293 * t308)) * t282 + t301 + (t261 * qJD(2) * t289 - t293 * t262 * t282) * t439 / 0.2e1 + (-t58 * (t322 - t405) - t57 * (t304 - t407) + t21 * (t350 + t447) + (t323 * t58 + t337 * t42) * t293 + (t323 * t57 + t337 * t41) * t292 + (t102 * t422 + (-t213 + t218) * t457 + t513 + t515) * t34) * m(6) + (-t92 * (-t405 + t443) - t91 * (-t407 + t507) + t54 * (t404 + t447) + (t326 * t92 + t351 * t87) * t293 + (t326 * t91 + t351 * t86) * t292 + ((-t213 + t214) * t457 + t515 + t525) * t45) * m(5) + (t81 * (t446 + t447) + (t120 * t416 + t139 * t396) * t293 + (t119 * t416 + t138 * t396) * t292 - t139 * (-t268 * t457 - t405) - t138 * (-t268 * t272 - t407) + t510) * m(4) + (-t440 + 0.1e1) * (rSges(3,1) * t297 - rSges(3,2) * t295) * t550; t301 + (-t58 * t322 - t57 * t304 + t21 * t350 + (t358 * t58 + t397 * t42) * t293 + (t57 * t358 + t41 * t397) * t292 + (t218 * t457 - (-t423 - t506) * t292 + t502 + t513) * t34) * m(6) + (-t92 * t443 - t91 * t507 + t54 * t404 + (t399 * t92 + t414 * t87) * t293 + (t399 * t91 + t414 * t86) * t292 + (t214 * t457 + t292 * t506 + t502 + t525) * t45) * m(5) + (-(-t138 * t272 - t139 * t457) * t268 + t81 * t446 + (-t119 * t292 - t120 * t293) * t267 + (-t138 * t292 - t139 * t293) * t231 + t510) * m(4); (-m(5) * t45 - m(6) * t34) * (-t272 * t293 + t292 * t457) + m(5) * (t292 * t87 - t293 * t86) + m(6) * t391; t460 * t534 + (t277 * t389 - t278 * t64) * t355 + (t318 + (t19 * t292 + t20 * t293 + t290 * t64) * t277) * t496 + t3 * t461 / 0.2e1 + (t277 * t390 - t278 * t63) * t356 + (t319 + (t17 * t292 + t18 * t293 + t290 * t63) * t277) * t494 + t23 * t462 / 0.2e1 - t278 * (t11 * t216 + t12 * t215 + (t462 * t70 + t320) * qJD(5)) / 0.2e1 + (t277 * t388 - t278 * t70) * t401 + (t320 + (t11 * t292 + t12 * t293 + t290 * t70) * t277) * t419 + (t243 * t305 + t244 * t306 - t293 * t313) * t497 + (t241 * t305 + t242 * t306 - t292 * t313) * t495 + (t321 * t278 + (-t294 * t305 + t306 * t296) * t277) * t418 + t511 * t459 / 0.2e1 + ((t42 * t101 - t41 * t102 - t57 * t80 + t58 * t79 + (t34 * t364 + (t292 * t58 - t293 * t57) * t169) * t290) * t278 + (t58 * (-t101 * t290 + t292 * t88) + t57 * (t102 * t290 - t293 * t88) + t21 * t364 + t34 * (-t292 * t80 + t293 * t79) + t391 * t169) * t277 - t58 * (t115 * t436 + t216 * t223) - t57 * (-t116 * t436 - t215 * t223) - t34 * (t115 * t215 - t116 * t216)) * m(6);];
tauc = t1(:);
