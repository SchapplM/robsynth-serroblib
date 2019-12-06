% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:06
% EndTime: 2019-12-05 18:16:23
% DurationCPUTime: 8.84s
% Computational Cost: add. (17538->530), mult. (11058->686), div. (0->0), fcn. (8580->10), ass. (0->335)
t276 = qJD(1) ^ 2;
t269 = qJ(4) + qJ(5);
t262 = sin(t269);
t263 = cos(t269);
t208 = rSges(6,1) * t262 + rSges(6,2) * t263;
t268 = qJ(1) + pkin(9);
t261 = qJ(3) + t268;
t256 = cos(t261);
t272 = cos(qJ(4));
t476 = pkin(4) * t272;
t257 = pkin(3) + t476;
t219 = t256 * t257;
t255 = sin(t261);
t274 = -pkin(8) - pkin(7);
t472 = pkin(7) + t274;
t477 = pkin(3) * t256;
t127 = t255 * t472 - t219 + t477;
t426 = t256 * t262;
t221 = rSges(6,2) * t426;
t425 = t256 * t263;
t384 = rSges(6,1) * t425;
t316 = rSges(6,3) * t255 + t384;
t143 = -t221 + t316;
t414 = t127 - t143;
t432 = t255 * t262;
t217 = Icges(6,4) * t432;
t431 = t255 * t263;
t139 = -Icges(6,1) * t431 + Icges(6,5) * t256 + t217;
t218 = Icges(6,4) * t426;
t140 = Icges(6,1) * t425 + Icges(6,5) * t255 - t218;
t266 = qJD(4) + qJD(5);
t185 = t255 * t266;
t186 = t256 * t266;
t267 = qJD(1) + qJD(3);
t254 = Icges(6,4) * t263;
t329 = -Icges(6,2) * t262 + t254;
t515 = Icges(6,1) * t262 + t254;
t528 = t515 + t329;
t284 = t185 * (-Icges(6,2) * t425 + t140 - t218) + t186 * (Icges(6,2) * t431 + t139 + t217) + t267 * t528;
t305 = t329 * t255;
t137 = Icges(6,6) * t256 - t305;
t138 = Icges(6,4) * t425 - Icges(6,2) * t426 + Icges(6,6) * t255;
t457 = Icges(6,4) * t262;
t204 = Icges(6,2) * t263 + t457;
t207 = Icges(6,1) * t263 - t457;
t493 = t185 * (t256 * t515 + t138) + t186 * (-t255 * t515 + t137) + t267 * (t204 - t207);
t534 = t284 * t262 + t493 * t263;
t533 = 2 * qJD(4);
t160 = t208 * t255;
t437 = t208 * t256;
t388 = rSges(6,1) * t431;
t529 = t208 * t266;
t378 = -t529 * t256 - t267 * t388;
t517 = rSges(6,2) * t432 + t256 * rSges(6,3);
t89 = t267 * t517 + t378;
t532 = t160 * t185 + t186 * t437 + t256 * t89;
t270 = sin(qJ(4));
t423 = t256 * t270;
t264 = Icges(5,4) * t272;
t330 = -Icges(5,2) * t270 + t264;
t514 = Icges(5,1) * t270 + t264;
t397 = t514 + t330;
t458 = Icges(5,4) * t270;
t240 = Icges(5,2) * t272 + t458;
t243 = Icges(5,1) * t272 - t458;
t398 = t240 - t243;
t531 = (t270 * t397 + t272 * t398) * t267;
t467 = rSges(6,1) * t263;
t209 = -rSges(6,2) * t262 + t467;
t430 = t255 * t267;
t530 = t267 * t160 - t186 * t209 - t208 * t430;
t473 = pkin(3) - t257;
t519 = t472 * t256;
t527 = -t255 * t473 + t519;
t260 = cos(t268);
t273 = cos(qJ(1));
t480 = pkin(1) * t273;
t340 = pkin(2) * t260 + t480;
t509 = t340 * qJD(1);
t259 = sin(t268);
t478 = pkin(2) * t259;
t271 = sin(qJ(1));
t481 = pkin(1) * t271;
t341 = t478 + t481;
t313 = t341 * qJD(1);
t181 = t209 * t266;
t424 = t256 * t267;
t475 = pkin(7) * t255;
t188 = t475 + t477;
t393 = qJD(4) * t270;
t374 = t255 * t393;
t228 = pkin(4) * t374;
t347 = t185 * t208 + t228;
t58 = -t509 + (-t188 + t414) * t267 + t347;
t526 = (t181 * t255 - t185 * t209 + t208 * t424 - t267 * t437) * t58;
t136 = Icges(6,5) * t425 - Icges(6,6) * t426 + Icges(6,3) * t255;
t63 = t256 * t136 + t138 * t432 - t140 * t431;
t325 = t204 * t262 - t263 * t515;
t202 = Icges(6,5) * t262 + Icges(6,6) * t263;
t440 = t202 * t256;
t92 = t255 * t325 + t440;
t525 = t185 * t63 + t92 * t267;
t470 = rSges(4,1) * t256;
t187 = -t255 * rSges(4,2) + t470;
t421 = t267 * t187;
t134 = -t509 - t421;
t334 = rSges(4,1) * t255 + rSges(4,2) * t256;
t524 = t134 * t334;
t422 = t256 * t272;
t148 = Icges(5,4) * t422 - Icges(5,2) * t423 + Icges(5,6) * t255;
t235 = Icges(5,4) * t423;
t150 = Icges(5,1) * t422 + Icges(5,5) * t255 - t235;
t326 = t148 * t270 - t150 * t272;
t521 = t326 * t256;
t328 = t138 * t262 - t140 * t263;
t520 = t328 * t256;
t164 = t334 * t267;
t237 = rSges(5,2) * t423;
t389 = rSges(5,1) * t422;
t317 = -rSges(5,3) * t255 - t389;
t152 = -t237 - t317;
t142 = t267 * t152;
t244 = rSges(5,1) * t270 + rSges(5,2) * t272;
t395 = qJD(4) * t255;
t184 = t244 * t395;
t518 = -t142 + t184;
t429 = t255 * t270;
t399 = rSges(5,2) * t429 + t256 * rSges(5,3);
t253 = t259 * rSges(3,2);
t471 = rSges(3,1) * t260;
t338 = -t471 - t480;
t516 = t253 + t338;
t466 = rSges(5,2) * t270;
t469 = rSges(5,1) * t272;
t333 = -t466 + t469;
t222 = t333 * qJD(4);
t512 = -t222 * t256 + t244 * t430;
t438 = t204 * t266;
t510 = -Icges(6,6) * t267 + t438;
t508 = t414 * t267 + t347;
t306 = t330 * t267;
t502 = -Icges(5,6) * t267 + qJD(4) * t240;
t102 = t502 * t255 - t256 * t306;
t308 = t267 * t243;
t500 = -Icges(5,5) * t267 + qJD(4) * t514;
t104 = t500 * t255 - t256 * t308;
t239 = Icges(5,5) * t272 - Icges(5,6) * t270;
t304 = t239 * t255;
t145 = Icges(5,3) * t256 - t304;
t147 = Icges(5,6) * t256 - t255 * t330;
t234 = Icges(5,4) * t429;
t428 = t255 * t272;
t149 = -Icges(5,1) * t428 + Icges(5,5) * t256 + t234;
t96 = t147 * t272 + t149 * t270;
t507 = qJD(4) * t96 + t102 * t270 - t104 * t272 - t145 * t267;
t101 = -t255 * t306 - t502 * t256;
t103 = -t255 * t308 - t500 * t256;
t146 = Icges(5,5) * t422 - Icges(5,6) * t423 + Icges(5,3) * t255;
t97 = t148 * t272 + t150 * t270;
t506 = qJD(4) * t97 + t101 * t270 - t103 * t272 - t146 * t267;
t505 = -Icges(6,3) * t267 + t202 * t266;
t504 = -Icges(6,5) * t267 + t266 * t515;
t238 = Icges(5,5) * t270 + Icges(5,6) * t272;
t503 = -Icges(5,3) * t267 + qJD(4) * t238;
t213 = t330 * qJD(4);
t214 = t243 * qJD(4);
t501 = qJD(4) * (t240 * t272 + t270 * t514) + t213 * t270 - t214 * t272 - t238 * t267;
t406 = -Icges(5,2) * t422 + t150 - t235;
t408 = t256 * t514 + t148;
t498 = t270 * t406 + t272 * t408;
t407 = Icges(5,2) * t428 + t149 + t234;
t409 = -t255 * t514 + t147;
t497 = -t270 * t407 - t272 * t409;
t348 = -t207 * t266 + t438;
t349 = t528 * t266;
t496 = -t202 * t267 + t262 * t349 + t263 * t348;
t354 = t140 * t266 - t510 * t256 - t267 * t305;
t307 = t267 * t207;
t356 = t138 * t266 + t255 * t307 + t504 * t256;
t495 = -t136 * t267 + t262 * t354 + t263 * t356;
t203 = Icges(6,5) * t263 - Icges(6,6) * t262;
t303 = t203 * t255;
t135 = Icges(6,3) * t256 - t303;
t355 = t139 * t266 + t510 * t255 - t329 * t424;
t357 = t137 * t266 - t504 * t255 + t256 * t307;
t494 = -t135 * t267 + t262 * t355 + t263 * t357;
t162 = t266 * t430;
t492 = -t162 / 0.2e1;
t163 = t267 * t186;
t491 = t163 / 0.2e1;
t490 = -t185 / 0.2e1;
t489 = t185 / 0.2e1;
t488 = -t186 / 0.2e1;
t487 = t186 / 0.2e1;
t486 = t255 / 0.2e1;
t485 = t256 / 0.2e1;
t484 = -t267 / 0.2e1;
t483 = t267 / 0.2e1;
t482 = -rSges(5,3) - pkin(7);
t479 = pkin(1) * t276;
t474 = t256 * pkin(7);
t76 = t184 + (-t152 - t188) * t267 - t509;
t463 = t255 * t76;
t154 = t202 * t255;
t93 = -t256 * t325 + t154;
t462 = t93 * t267;
t367 = t256 * t473;
t419 = t267 * t274;
t401 = t255 * t419 + t228;
t108 = (t367 + t475) * t267 + t401;
t379 = t267 * t221 + t529 * t255;
t90 = -t267 * t316 + t379;
t461 = -t108 - t90;
t169 = t238 * t255;
t323 = t270 * t240 - t272 * t514;
t111 = -t256 * t323 + t169;
t448 = t111 * t267;
t445 = t137 * t262;
t444 = t139 * t263;
t443 = t147 * t270;
t442 = t149 * t272;
t435 = t238 * t256;
t434 = t239 * t267;
t175 = t244 * t255;
t433 = t244 * t256;
t427 = t256 * t127;
t418 = t256 * t135 + t137 * t432;
t417 = -t255 * t135 - t139 * t425;
t416 = t256 * t145 + t147 * t429;
t415 = t255 * t145 + t149 * t422;
t405 = -t256 * t419 - t257 * t430;
t258 = t271 * t479;
t396 = t276 * t478 + t258;
t394 = qJD(4) * t256;
t392 = qJD(4) * t272;
t391 = (qJD(4) ^ 2) * t476;
t390 = rSges(5,1) * t428;
t385 = pkin(4) * t393;
t371 = t256 * t392;
t372 = t256 * t393;
t377 = -rSges(5,1) * t372 - rSges(5,2) * t371 - t267 * t390;
t376 = rSges(5,1) * t374 + (t255 * t392 + t267 * t423) * rSges(5,2);
t375 = t244 * t394;
t370 = -t430 / 0.2e1;
t369 = t424 / 0.2e1;
t368 = -pkin(3) - t469;
t364 = -t395 / 0.2e1;
t362 = -t394 / 0.2e1;
t361 = t394 / 0.2e1;
t358 = -t257 - t467;
t351 = -t146 - t442;
t346 = pkin(4) * t372;
t345 = -pkin(3) - t358;
t342 = t399 + t474;
t339 = -t256 * t274 + t517;
t335 = rSges(3,1) * t259 + rSges(3,2) * t260;
t81 = t138 * t263 + t140 * t262;
t327 = -t442 + t443;
t322 = t390 - t399;
t321 = -t219 + t221 - t384;
t320 = t388 - t517;
t64 = -t137 * t426 - t417;
t67 = t256 * t146 + t148 * t429 - t150 * t428;
t141 = t267 * t322;
t182 = t267 * (-pkin(3) * t255 + t474);
t319 = t141 - t182 + t375;
t318 = t519 - t517;
t314 = t340 * t276;
t66 = -t149 * t428 + t416;
t310 = (t255 * t67 + t256 * t66) * qJD(4);
t68 = -t147 * t423 + t415;
t69 = t146 * t255 - t521;
t309 = (t255 * t69 + t256 * t68) * qJD(4);
t302 = t237 - t389 - t477;
t300 = -t267 ^ 2 * t188 - t314;
t297 = t154 * t186 - t185 * t440 + t203 * t267;
t183 = t267 * t188;
t295 = t183 + t509;
t294 = -t505 * t256 + (-t303 + t328) * t267;
t293 = -t203 * t424 + t505 * t255 + (-t444 + t445) * t267;
t292 = -t503 * t256 + (-t304 + t326) * t267;
t291 = -t239 * t424 + t503 * t255 + t327 * t267;
t290 = t203 * t266 + t267 * t325;
t289 = t239 * qJD(4) + t267 * t323;
t13 = t293 * t255 - t494 * t256;
t14 = t294 * t255 - t495 * t256;
t15 = t494 * t255 + t293 * t256;
t16 = t495 * t255 + t294 * t256;
t62 = -t139 * t431 + t418;
t28 = t186 * t62 + t525;
t65 = t136 * t255 - t520;
t29 = t185 * t65 + t186 * t64 + t462;
t40 = -t262 * t357 + t263 * t355;
t41 = -t262 * t356 + t263 * t354;
t44 = t290 * t255 - t496 * t256;
t45 = t496 * t255 + t290 * t256;
t80 = t137 * t263 + t139 * t262;
t288 = (t13 * t186 + t14 * t185 - t162 * t64 + t163 * t65 + t267 * t44) * t486 + (t297 * t255 - t534 * t256) * t490 + (t534 * t255 + t297 * t256) * t488 + (t15 * t186 + t16 * t185 - t162 * t62 + t163 * t63 + t267 * t45) * t485 + (-t262 * t493 + t263 * t284) * t484 + t28 * t370 + t29 * t369 + ((t267 * t65 + t13) * t256 + (-t267 * t64 + t14) * t255) * t489 + (t255 * t63 + t256 * t62) * t492 + (t255 * t65 + t256 * t64) * t491 + ((t267 * t63 + t15) * t256 + (-t267 * t62 + t16) * t255) * t487 + ((t267 * t81 + t40) * t256 + (-t267 * t80 + t41) * t255) * t483;
t287 = t186 * t208 - t182 + t346 + (t320 + t527) * t267;
t286 = t256 * t152 + t255 * t322;
t281 = (-t136 - t444) * t255 + t520 + t418;
t105 = t399 * t267 + t377;
t106 = t267 * t317 + t376;
t280 = (-t106 - t142) * t255 + (t105 + t141) * t256;
t110 = t255 * t323 + t435;
t109 = t110 * t267;
t32 = t109 + t310;
t33 = t309 + t448;
t49 = -qJD(4) * t327 + t102 * t272 + t104 * t270;
t50 = -qJD(4) * t326 + t101 * t272 + t103 * t270;
t54 = t289 * t255 - t501 * t256;
t55 = t501 * t255 + t289 * t256;
t279 = ((t65 + t281) * t186 + t525) * t490 + (t109 + ((t416 + t69 + t521) * t256 + (-t68 + (t351 - t443) * t256 + t67 + t415) * t255) * qJD(4)) * t364 + (t80 + t92) * t492 + (t81 + t93) * t491 + (-t462 + (t63 + (-t136 + t445) * t256 - t328 * t255 + t417) * t186 + (-t62 + t281) * t185 + t29) * t488 + (t40 + t45) * t487 + (t33 - t448 + ((t67 + (-t146 + t443) * t256 - t415) * t256 + (t255 * t351 + t416 - t66) * t255) * qJD(4)) * t362 + (t49 + t55) * t361 + (t28 + t41 + t44) * t489 + (t50 + t54 + t32) * t395 / 0.2e1 + ((t110 + t96) * t364 + (t111 + t97) * t361 - qJD(4) * t323 + t213 * t272 + t214 * t270 - t262 * t348 + t263 * t349) * t267;
t231 = pkin(3) * t430;
t166 = pkin(7) * t424 - t231;
t60 = t222 * t395 + (-t105 - t166 + t375) * t267 + t396;
t61 = t512 * qJD(4) + t106 * t267 + t300;
t75 = t313 + t319;
t278 = t76 * (t231 - t377) + (t61 * t368 + t60 * t482) * t255 + ((-t466 * t76 - t482 * t75) * t255 + (-t368 * t75 + t482 * t76) * t256) * t267 - t75 * t376;
t107 = t231 + (-pkin(7) * t267 - t385) * t256 + t405;
t36 = t255 * t391 + t163 * t208 + t181 * t185 + (-t107 - t166 - t89 + t346) * t267 + t396;
t37 = -t256 * t391 + t162 * t208 - t181 * t186 + (t228 - t461) * t267 + t300;
t57 = t313 + t287;
t277 = t58 * (t346 - t378 - t405) + (t36 * (-rSges(6,3) + t274) + t37 * t358) * t255 + (-t58 * t517 - t57 * (-t316 - t219)) * t267 - t57 * (t379 + t401);
t230 = rSges(4,2) * t430;
t165 = -rSges(4,1) * t424 + t230;
t133 = t313 + t164;
t123 = t256 * t143;
t119 = t165 * t267 - t314;
t118 = t164 * t267 + t396;
t82 = qJD(4) * t286 + qJD(2);
t51 = qJD(2) + t186 * t143 + t185 * t320 + (t527 * t255 - t427) * qJD(4);
t46 = t280 * qJD(4);
t12 = -t162 * t143 + t186 * t89 + t163 * t320 - t185 * t90 + ((t267 * t519 + t107) * t256 + (-t108 + (t127 - t367) * t267) * t255) * qJD(4);
t1 = [m(3) * ((t276 * t335 + t258) * t516 + (t273 * t479 + (-0.2e1 * t253 - t338 + t471 + t516) * t276) * (t335 + t481)) + t279 + (-(t295 + t58 - t508) * t57 + t36 * (t321 - t340) + t37 * (t339 - t341) + (t340 * t57 + t341 * t58) * qJD(1) + t277) * m(6) + (-(t76 + t295 - t518) * t75 + t60 * (t302 - t340) + t61 * (-t341 + t342) + (t340 * t75 + t341 * t76) * qJD(1) + t278) * m(5) + (t118 * (-t187 - t340) + t119 * (-t334 - t341) + t524 * t267 + t134 * t313 + (t470 * t267 - t134 - t230 - t421) * t133) * m(4); m(5) * t46 + m(6) * t12; t279 + (t36 * t321 + t37 * t339 + t277 - t58 * t287 + t57 * (-t183 + t508)) * m(6) + (t60 * t302 + t61 * t342 + t278 - t76 * t319 + t75 * (-t183 + t518)) * m(5) + (-(t133 * t187 + t524) * t267 - t118 * t187 - t119 * t334 - t133 * t165 + t134 * t164) * m(4); ((-t270 * t398 + t272 * t397) * t267 + ((t255 * t406 + t256 * t407) * t272 + (-t255 * t408 - t256 * t409) * t270) * qJD(4)) * t484 + ((t267 * t97 + t49) * t256 + (-t267 * t96 + t50) * t255) * t483 + t288 + ((t169 * t394 + t434) * t256 + (t531 + (t498 * t255 + (-t435 - t497) * t256) * qJD(4)) * t255) * t362 + ((-t395 * t435 + t434) * t255 + (-t531 + (t497 * t256 + (t169 - t498) * t255) * qJD(4)) * t256) * t364 + (t267 * t54 + ((t291 * t255 - t507 * t256 + t267 * t69) * t256 + (t292 * t255 - t506 * t256 - t267 * t68) * t255) * t533) * t486 + (t267 * t55 + ((t507 * t255 + t291 * t256 + t267 * t67) * t256 + (t506 * t255 + t292 * t256 - t267 * t66) * t255) * t533) * t485 + (t310 + t32) * t370 + (t309 + t33) * t369 + (t12 * (-t427 + t123 + (t255 * t345 + t318) * t255) + (t255 * t36 - t256 * t37) * (pkin(4) * t270 + t208) + (-(-t255 ^ 2 - t256 ^ 2) * t385 + t256 * t107 + t461 * t255 + (t318 * t256 + (t256 * t345 + t414) * t255) * t267 + t532) * t51 + t526 + (-pkin(4) * t371 - (-pkin(4) * t392 - t181) * t256 + t530) * t57) * m(6) + (t76 * t244 * t424 + t60 * t175 + t222 * t463 + t82 * t280 + t46 * t286 - t61 * t433 - t512 * t75 - (-t175 * t75 + t433 * t76) * t267 - (t82 * (-t175 * t255 - t256 * t433) + (t256 * t75 + t463) * t333) * qJD(4)) * m(5); t288 + (t12 * (t255 * t320 + t123) + t36 * t160 - t37 * t437 + (-t255 * t90 + (-t255 * t143 + t256 * t320) * t267 + t532) * t51 + t526 + (t256 * t181 + t530) * t57) * m(6);];
tauc = t1(:);
