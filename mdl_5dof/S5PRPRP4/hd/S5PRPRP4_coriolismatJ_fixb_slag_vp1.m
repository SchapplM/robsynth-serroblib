% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRPRP4
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRPRP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:06
% EndTime: 2019-12-05 15:35:32
% DurationCPUTime: 11.86s
% Computational Cost: add. (20782->529), mult. (28564->781), div. (0->0), fcn. (31054->8), ass. (0->327)
t504 = rSges(6,1) + pkin(4);
t499 = rSges(6,3) + qJ(5);
t281 = qJ(2) + pkin(8);
t278 = cos(t281);
t283 = cos(pkin(7));
t285 = sin(qJ(4));
t403 = t283 * t285;
t282 = sin(pkin(7));
t287 = cos(qJ(4));
t404 = t282 * t287;
t255 = t278 * t404 - t403;
t239 = Icges(6,5) * t255;
t402 = t283 * t287;
t405 = t282 * t285;
t254 = t278 * t405 + t402;
t277 = sin(t281);
t410 = t277 * t282;
t139 = Icges(6,6) * t410 + Icges(6,3) * t254 + t239;
t432 = Icges(5,4) * t255;
t145 = -Icges(5,2) * t254 + Icges(5,6) * t410 + t432;
t538 = t139 - t145;
t257 = t278 * t402 + t405;
t240 = Icges(6,5) * t257;
t256 = t278 * t403 - t404;
t409 = t277 * t283;
t140 = Icges(6,6) * t409 + Icges(6,3) * t256 + t240;
t431 = Icges(5,4) * t257;
t146 = -Icges(5,2) * t256 + Icges(5,6) * t409 + t431;
t537 = t140 - t146;
t141 = Icges(5,5) * t255 - Icges(5,6) * t254 + Icges(5,3) * t410;
t143 = Icges(6,4) * t255 + Icges(6,2) * t410 + Icges(6,6) * t254;
t523 = t141 + t143;
t142 = Icges(5,5) * t257 - Icges(5,6) * t256 + Icges(5,3) * t409;
t144 = Icges(6,4) * t257 + Icges(6,2) * t409 + Icges(6,6) * t256;
t522 = t142 + t144;
t426 = Icges(6,5) * t254;
t147 = Icges(6,1) * t255 + Icges(6,4) * t410 + t426;
t241 = Icges(5,4) * t254;
t149 = Icges(5,1) * t255 + Icges(5,5) * t410 - t241;
t536 = t147 + t149;
t425 = Icges(6,5) * t256;
t148 = Icges(6,1) * t257 + Icges(6,4) * t409 + t425;
t242 = Icges(5,4) * t256;
t150 = Icges(5,1) * t257 + Icges(5,5) * t409 - t242;
t535 = t148 + t150;
t429 = Icges(5,4) * t287;
t336 = -Icges(5,2) * t285 + t429;
t207 = -Icges(5,6) * t278 + t277 * t336;
t407 = t277 * t287;
t272 = Icges(6,5) * t407;
t408 = t277 * t285;
t421 = Icges(6,6) * t278;
t521 = Icges(6,3) * t408 - t207 + t272 - t421;
t424 = Icges(6,5) * t285;
t341 = Icges(6,1) * t287 + t424;
t209 = -Icges(6,4) * t278 + t277 * t341;
t430 = Icges(5,4) * t285;
t342 = Icges(5,1) * t287 - t430;
t211 = -Icges(5,5) * t278 + t277 * t342;
t519 = t209 + t211;
t162 = -Icges(5,5) * t254 - Icges(5,6) * t255;
t164 = -Icges(6,4) * t254 + Icges(6,6) * t255;
t534 = t162 + t164;
t163 = -Icges(5,5) * t256 - Icges(5,6) * t257;
t165 = -Icges(6,4) * t256 + Icges(6,6) * t257;
t533 = t163 + t165;
t332 = Icges(5,5) * t287 - Icges(5,6) * t285;
t203 = -Icges(5,3) * t278 + t277 * t332;
t335 = Icges(6,4) * t287 + Icges(6,6) * t285;
t205 = -Icges(6,2) * t278 + t277 * t335;
t520 = t203 + t205;
t394 = Icges(5,2) * t257 - t150 + t242;
t396 = Icges(6,3) * t257 - t148 - t425;
t532 = t394 + t396;
t395 = Icges(5,2) * t255 - t149 + t241;
t397 = Icges(6,3) * t255 - t147 - t426;
t531 = t395 + t397;
t398 = -Icges(5,1) * t256 - t146 - t431;
t400 = -Icges(6,1) * t256 + t140 + t240;
t530 = t398 + t400;
t399 = -Icges(5,1) * t254 - t145 - t432;
t401 = -Icges(6,1) * t254 + t139 + t239;
t529 = t399 + t401;
t528 = t537 * t254 + t535 * t255 + t522 * t410;
t527 = t538 * t256 + t536 * t257 + t523 * t409;
t517 = t499 * t285 + t504 * t287;
t491 = rSges(6,2) * t278 - t277 * t517;
t526 = t491 * t282;
t525 = t491 * t283;
t279 = t282 ^ 2;
t280 = t283 ^ 2;
t473 = t279 + t280;
t518 = t521 * t285 + t519 * t287;
t516 = Icges(4,5) * t277 + Icges(4,6) * t278;
t515 = t531 * t254 + t529 * t255 + t534 * t410;
t514 = t532 * t254 + t530 * t255 + t533 * t410;
t513 = t531 * t256 + t529 * t257 + t534 * t409;
t512 = t532 * t256 + t530 * t257 + t533 * t409;
t352 = rSges(5,1) * t287 - rSges(5,2) * t285;
t214 = -rSges(5,3) * t278 + t277 * t352;
t196 = t214 * t282;
t198 = t214 * t283;
t509 = t519 + (t424 - t430 + (-Icges(5,2) - Icges(6,3)) * t287) * t277;
t508 = -(-Icges(5,1) * t285 - t429) * t277 + Icges(6,1) * t408 - t272 - t521;
t244 = (-Icges(5,5) * t285 - Icges(5,6) * t287) * t277;
t245 = (-Icges(6,4) * t285 + Icges(6,6) * t287) * t277;
t507 = (-t244 - t245) * t278;
t484 = t520 * t278;
t286 = sin(qJ(2));
t288 = cos(qJ(2));
t505 = 0.2e1 * t288 * (Icges(3,1) - Icges(3,2)) * t286 + (-0.2e1 * t286 ^ 2 + 0.2e1 * t288 ^ 2) * Icges(3,4);
t331 = Icges(6,5) * t287 + Icges(6,3) * t285;
t295 = -t277 * t331 + t421;
t183 = t295 * t282;
t189 = t207 * t282;
t191 = t209 * t282;
t193 = t211 * t282;
t326 = -t145 * t285 + t149 * t287;
t314 = -t203 * t282 - t326;
t328 = t139 * t285 + t147 * t287;
t316 = t205 * t282 + t328;
t503 = (t316 - t314) * t278 + ((-t191 - t193) * t287 + (t183 + t189) * t285 + t523) * t277;
t184 = t295 * t283;
t190 = t207 * t283;
t192 = t209 * t283;
t194 = t211 * t283;
t325 = -t146 * t285 + t150 * t287;
t313 = -t203 * t283 - t325;
t327 = t140 * t285 + t148 * t287;
t315 = t205 * t283 + t327;
t502 = (-t313 + t315) * t278 + ((-t192 - t194) * t287 + (t184 + t190) * t285 + t522) * t277;
t501 = t538 * t254 + t536 * t255 + t523 * t410;
t500 = t537 * t256 + t535 * t257 + t522 * t409;
t334 = -Icges(3,5) * t286 - Icges(3,6) * t288;
t258 = t334 * t282;
t259 = t334 * t283;
t498 = t254 * t521 + t519 * t255 + t410 * t520;
t497 = t256 * t521 + t257 * t519 + t409 * t520;
t496 = t277 * t518 - t484;
t495 = (-t331 + t336) * t278 + (Icges(5,6) - Icges(6,6)) * t277;
t494 = (-t341 - t342) * t278 + (-Icges(6,4) - Icges(5,5)) * t277;
t493 = -t282 * t516 + t258;
t492 = -t283 * t516 + t259;
t490 = ((t332 + t335) * t278 + (Icges(6,2) + Icges(5,3)) * t277 - t518) * t278;
t417 = t143 * t278;
t93 = t277 * t328 - t417;
t416 = t144 * t278;
t94 = t277 * t327 - t416;
t419 = t141 * t278;
t95 = t277 * t326 - t419;
t418 = t142 * t278;
t96 = t277 * t325 - t418;
t485 = (t94 + t96) * t283 + (t93 + t95) * t282;
t483 = t527 * t282;
t482 = t528 * t283;
t174 = -rSges(5,1) * t254 - rSges(5,2) * t255;
t178 = -rSges(5,1) * t256 - rSges(5,2) * t257;
t118 = t174 * t282 + t178 * t283;
t266 = pkin(3) * t277 - pkin(6) * t278;
t449 = pkin(2) * t286;
t362 = -t266 - t449;
t321 = t362 + t491;
t121 = t321 * t282;
t123 = t321 * t283;
t354 = -t214 + t362;
t134 = t354 * t282;
t136 = t354 * t283;
t378 = (-t504 * t285 + t499 * t287) * t277;
t157 = t378 * t282;
t158 = t378 * t283;
t251 = (-rSges(5,1) * t285 - rSges(5,2) * t287) * t277;
t267 = pkin(3) * t278 + pkin(6) * t277;
t448 = pkin(2) * t288;
t387 = t473 * t448;
t356 = t473 * t267 + t387;
t392 = rSges(6,2) * t409 + t499 * t256 + t504 * t257;
t393 = rSges(6,2) * t410 + t499 * t254 + t504 * t255;
t70 = t282 * t393 + t283 * t392 + t356;
t152 = rSges(5,1) * t255 - rSges(5,2) * t254 + rSges(5,3) * t410;
t154 = rSges(5,1) * t257 - rSges(5,2) * t256 + rSges(5,3) * t409;
t78 = t152 * t282 + t154 * t283 + t356;
t390 = -t504 * t256 + t499 * t257;
t391 = -t504 * t254 + t499 * t255;
t92 = t282 * t391 + t283 * t390;
t481 = -m(6) * (-t121 * t157 - t123 * t158 + t70 * t92) - m(5) * (t118 * t78 + (-t134 * t282 - t136 * t283) * t251);
t480 = -t277 / 0.2e1;
t479 = t277 / 0.2e1;
t478 = -t278 / 0.2e1;
t477 = t278 / 0.2e1;
t476 = -t282 / 0.2e1;
t452 = t282 / 0.2e1;
t451 = -t283 / 0.2e1;
t450 = t283 / 0.2e1;
t475 = (t254 * t509 + t255 * t508) * t278 + (t514 * t283 + (t507 + t515) * t282) * t277;
t474 = (t256 * t509 + t257 * t508) * t278 + ((t507 + t512) * t283 + t513 * t282) * t277;
t159 = t254 * t282 + t256 * t283;
t406 = t278 * t285;
t127 = (-t159 + t406) * t408;
t365 = t406 / 0.2e1;
t131 = (t365 - t159 / 0.2e1) * m(6);
t443 = m(6) * qJD(5);
t472 = t131 * qJD(1) + t127 * t443;
t470 = t516 * t473;
t468 = t278 ^ 2;
t467 = 2 * qJD(2);
t466 = 2 * qJD(4);
t465 = 4 * qJD(4);
t464 = m(4) / 0.2e1;
t463 = m(5) / 0.2e1;
t462 = m(6) / 0.2e1;
t100 = t277 * t525 - t278 * t392;
t99 = -t277 * t526 + t278 * t393;
t347 = -t100 * t282 - t283 * t99;
t304 = -t282 * t392 + t283 * t393;
t48 = t304 * t278 + (-t282 * t525 + t283 * t526) * t277;
t381 = rSges(6,2) * t277 + t278 * t517;
t71 = (t282 * t381 - t393) * t277;
t72 = (-t283 * t381 + t392) * t277;
t77 = t304 * t277;
t461 = m(6) * (t254 * t72 + t256 * t71 + (t278 * t77 + (t347 + t48) * t277) * t285);
t460 = m(6) * (t100 * t72 + t48 * t77 + t71 * t99);
t324 = t152 * t283 - t154 * t282;
t110 = t324 * t277;
t116 = t152 * t278 + t214 * t410;
t117 = -t154 * t278 - t214 * t409;
t80 = t324 * t278 + (-t196 * t283 + t198 * t282) * t277;
t216 = rSges(5,3) * t277 + t278 * t352;
t97 = (t282 * t216 - t152) * t277;
t98 = (-t283 * t216 + t154) * t277;
t459 = m(5) * (t110 * t80 + t116 * t97 + t117 * t98);
t138 = (t254 * t283 - t256 * t282) * t277;
t276 = t277 ^ 2;
t181 = t254 * t278 + t276 * t405;
t182 = -t256 * t278 - t276 * t403;
t457 = m(6) * (t121 * t182 + t123 * t181 + t138 * t70 + t159 * t77 + t347 * t408);
t456 = m(6) * (t121 * t255 + t123 * t257 - t157 * t254 - t158 * t256 + (t285 * t92 + t287 * t70) * t277);
t445 = m(6) * qJD(2);
t444 = m(6) * qJD(4);
t377 = qJD(4) * t277;
t291 = -t277 * t316 + t417;
t50 = t183 * t254 - t191 * t255 + t282 * t291;
t290 = -t277 * t315 + t416;
t51 = t184 * t254 - t192 * t255 + t282 * t290;
t293 = t277 * t314 + t419;
t52 = t189 * t254 - t193 * t255 + t282 * t293;
t292 = t277 * t313 + t418;
t53 = t190 * t254 - t194 * t255 + t282 * t292;
t375 = ((t51 + t53) * t283 + (t50 + t52 - t490) * t282 + t498) * t479 + ((-t484 + t501) * t282 + t494 * t255 + t495 * t254 + t482) * t477;
t54 = t183 * t256 - t191 * t257 + t283 * t291;
t55 = t184 * t256 - t192 * t257 + t283 * t290;
t56 = t189 * t256 - t193 * t257 + t283 * t293;
t57 = t190 * t256 - t194 * t257 + t283 * t292;
t373 = ((t55 + t57 - t490) * t283 + (t54 + t56) * t282 + t497) * t479 + ((-t484 + t500) * t283 + t494 * t257 + t495 * t256 + t483) * t477;
t372 = (t502 * t283 + t503 * t282 + (t495 * t285 + t494 * t287 - t520) * t278 + t496) * t480 + (t485 + t490) * t478;
t371 = t515 * t450 + t514 * t476;
t370 = t513 * t451 + t512 * t452;
t369 = (t501 * t282 + t482) * t479 + t498 * t478;
t368 = (t500 * t283 + t483) * t479 + t497 * t478;
t367 = t496 * t477 + t485 * t480;
t366 = t407 / 0.2e1;
t264 = rSges(4,1) * t277 + rSges(4,2) * t278;
t364 = -t264 - t449;
t265 = rSges(4,1) * t278 - rSges(4,2) * t277;
t363 = -t265 - t448;
t361 = -t267 - t448;
t355 = t473 * t449;
t353 = -t216 + t361;
t270 = t286 * rSges(3,1) + rSges(3,2) * t288;
t289 = -m(5) * (t282 * t97 - t283 * t98) / 0.2e1 - m(6) * (t282 * t71 - t283 * t72) / 0.2e1;
t308 = (t157 * t283 - t158 * t282) * t462;
t18 = t308 + t289;
t19 = 0.2e1 * (t48 / 0.4e1 - t92 / 0.4e1) * m(6) + 0.2e1 * (t80 / 0.4e1 - t118 / 0.4e1) * m(5);
t330 = -t19 * qJD(1) + t18 * qJD(3);
t329 = -t121 * t282 - t123 * t283;
t320 = t361 - t381;
t307 = (t181 * t282 - t182 * t283) * t462;
t312 = m(6) * (-t255 * t283 + t257 * t282);
t105 = t307 - t312 / 0.2e1;
t128 = (t366 - t138 / 0.2e1) * m(6);
t319 = -t128 * qJD(1) + t105 * qJD(3);
t318 = -t371 - t375;
t317 = -t370 + t373;
t130 = -t473 * t264 - t355;
t311 = t130 * t265;
t309 = -t473 * t266 - t355;
t303 = t505 * t282 + t259;
t302 = -t283 * t505 + t258;
t218 = t363 * t283;
t217 = t363 * t282;
t180 = t473 * t270;
t137 = t353 * t283;
t135 = t353 * t282;
t132 = m(6) * t365 + t159 * t462;
t129 = m(6) * t366 + t138 * t462;
t126 = -t178 * t278 - t251 * t409;
t125 = t174 * t278 + t251 * t410;
t124 = t320 * t283;
t122 = t320 * t282;
t119 = t276 * t285 * t287 + t254 * t255 + t256 * t257;
t114 = (t174 * t283 - t178 * t282) * t277;
t104 = t307 + t312 / 0.2e1;
t103 = -t196 * t282 - t198 * t283 + t309;
t102 = -t158 * t277 - t278 * t390;
t101 = t278 * t391 + t378 * t410;
t90 = t282 * t526 + t283 * t525 + t309;
t89 = (-t282 * t390 + t283 * t391) * t277;
t76 = -t163 * t278 + (t285 * t394 + t287 * t398) * t277;
t75 = -t162 * t278 + (t285 * t395 + t287 * t399) * t277;
t74 = -t165 * t278 + (t285 * t396 + t287 * t400) * t277;
t73 = -t164 * t278 + (t285 * t397 + t287 * t401) * t277;
t43 = t159 * t70 + t329 * t408;
t33 = t100 * t182 + t138 * t77 + t181 * t99;
t32 = t282 * t57 - t283 * t56;
t31 = t282 * t55 - t283 * t54;
t30 = t282 * t53 - t283 * t52;
t29 = t282 * t51 - t283 * t50;
t26 = t456 / 0.2e1;
t20 = (t118 + t80) * t463 + (t48 + t92) * t462;
t17 = t308 - t289;
t15 = t457 / 0.2e1;
t6 = t461 / 0.2e1;
t5 = t26 + t6 - t457 / 0.2e1;
t4 = t15 + t26 - t461 / 0.2e1;
t3 = t15 + t6 - t456 / 0.2e1;
t2 = t282 * t370 + t283 * t371 - t481;
t1 = t459 + t460 + (t282 * t369 + t283 * t368 + t372) * t278 + (t282 * t375 + t283 * t373 - t367) * t277;
t7 = [0, t20 * qJD(4) + t132 * qJD(5) + (-m(3) * t180 / 0.2e1 + t130 * t464 + t103 * t463 + t90 * t462) * t467, 0, t20 * qJD(2) + (t114 * t463 + t462 * t89) * t466 + t129 * qJD(5), qJD(2) * t132 + qJD(4) * t129; -qJD(4) * t19 - qJD(5) * t131, t2 * qJD(4) + t43 * t443 + (m(4) * (t387 * t130 + (t218 * t364 + t283 * t311) * t283 + (t217 * t364 + t282 * t311) * t282) + m(6) * (t121 * t122 + t123 * t124 + t70 * t90) + m(5) * (t103 * t78 + t134 * t135 + t136 * t137) + m(3) * (-t180 + t270) * t473 * (rSges(3,1) * t288 - t286 * rSges(3,2)) + (t32 + t31 + (t303 * t283 - t470 + (t302 - t493) * t282) * t283 + t492 * t279) * t452 + (t30 + t29 + (t302 * t282 - t470 + (t303 - t492) * t283) * t282 + t493 * t280) * t451) * qJD(2), t18 * qJD(4), t2 * qJD(2) + t4 * qJD(5) + (-t459 / 0.4e1 - t460 / 0.4e1) * t465 + ((t110 * t118 + t114 * t78 + t125 * t136 + t126 * t134 + (-t116 * t283 - t117 * t282) * t251) * t463 + (-t100 * t157 + t101 * t123 + t102 * t121 - t158 * t99 + t70 * t89 + t77 * t92) * t462) * t466 + (t282 * t318 - t283 * t317 + t367) * t377 + t330 + (((t73 / 0.2e1 + t75 / 0.2e1 - t368) * t283 + (-t74 / 0.2e1 - t76 / 0.2e1 - t369) * t282 - t372) * t278 + t474 * t452 + t475 * t451) * qJD(4), t4 * qJD(4) + t43 * t445 - t472; 0, t17 * qJD(4) + ((-t122 * t283 + t124 * t282) * t462 + (-t135 * t283 + t137 * t282) * t463 + (-t217 * t283 + t218 * t282) * t464) * t467, 0, t17 * qJD(2) + ((t125 * t282 - t126 * t283) * t463 + (t101 * t282 - t102 * t283) * t462) * t466 + t104 * qJD(5), t104 * qJD(4); qJD(2) * t19 - qJD(5) * t128, t1 * qJD(4) + t3 * qJD(5) + ((t100 * t122 + t121 * t72 + t123 * t71 + t124 * t99 + t48 * t70 + t77 * t90) * t462 + (t103 * t110 + t116 * t137 + t117 * t135 + t134 * t98 + t136 * t97 + t78 * t80) * t463) * t467 - t330 + (t318 * t283 + t317 * t282 + (t502 * t476 + (t528 * t282 - t501 * t283) * t452 + (t500 * t282 - t527 * t283 + t503) * t450) * t278 + ((t32 / 0.2e1 + t31 / 0.2e1 - t95 / 0.2e1 - t93 / 0.2e1) * t283 + (t30 / 0.2e1 + t29 / 0.2e1 + t96 / 0.2e1 + t94 / 0.2e1) * t282) * t277 + t481) * qJD(2), -qJD(2) * t18 + qJD(5) * t105, t1 * qJD(2) + (m(6) * (t100 * t102 + t101 * t99 + t77 * t89) / 0.4e1 + m(5) * (t110 * t114 + t116 * t125 + t117 * t126) / 0.4e1) * t465 + t33 * t443 + (-t245 / 0.2e1 - t244 / 0.2e1) * qJD(4) * t278 * t468 + (((t74 + t76) * t283 + (t73 + t75) * t282 + (t285 * t509 + t287 * t508) * t278) * t478 + t475 * t452 + t474 * t450) * t377, t3 * qJD(2) + t33 * t444 + (t138 * t408 + t181 * t256 + t182 * t254 - t119) * t443 + t319; qJD(2) * t131 + qJD(4) * t128, (t122 * t254 + t124 * t256 - t43 + (t278 * t70 + (t329 + t90) * t277) * t285) * t445 + t5 * qJD(4) + t472, -t105 * qJD(4), t5 * qJD(2) + (t100 * t255 + t101 * t256 + t102 * t254 + t257 * t99 + (t285 * t89 + t287 * t77) * t277 - t33) * t444 + t119 * t443 - t319, 0.4e1 * (t127 * qJD(2) / 0.4e1 + t119 * qJD(4) / 0.4e1) * m(6);];
Cq = t7;
