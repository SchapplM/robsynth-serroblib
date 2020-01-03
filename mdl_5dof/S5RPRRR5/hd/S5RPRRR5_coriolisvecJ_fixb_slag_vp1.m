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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:54:03
% DurationCPUTime: 9.80s
% Computational Cost: add. (17538->543), mult. (11058->726), div. (0->0), fcn. (8580->10), ass. (0->333)
t265 = qJ(1) + pkin(9);
t256 = qJ(3) + t265;
t250 = cos(t256);
t269 = cos(qJ(4));
t458 = pkin(4) * t269;
t251 = pkin(3) + t458;
t208 = t250 * t251;
t249 = sin(t256);
t271 = -pkin(8) - pkin(7);
t459 = pkin(3) * t250;
t126 = t459 - t208 + (pkin(7) + t271) * t249;
t264 = qJD(1) + qJD(3);
t118 = t264 * t126;
t266 = qJ(4) + qJ(5);
t257 = sin(t266);
t418 = t250 * t257;
t210 = rSges(6,2) * t418;
t258 = cos(t266);
t417 = t250 * t258;
t372 = rSges(6,1) * t417;
t138 = rSges(6,3) * t249 - t210 + t372;
t127 = t264 * t138;
t186 = pkin(7) * t249 + t459;
t179 = t264 * t186;
t416 = t250 * t264;
t187 = t251 * t416;
t263 = qJD(4) + qJD(5);
t181 = t249 * t263;
t199 = rSges(6,1) * t257 + rSges(6,2) * t258;
t267 = sin(qJ(4));
t376 = qJD(4) * t267;
t358 = t249 * t376;
t336 = pkin(4) * t358;
t298 = -t181 * t199 - t336;
t421 = t249 * t264;
t396 = rSges(6,3) * t421 + t264 * t372;
t504 = t118 - t127 - t179 - t298 + t187 + t396;
t415 = t250 * t267;
t226 = rSges(5,2) * t415;
t414 = t250 * t269;
t373 = rSges(5,1) * t414;
t149 = rSges(5,3) * t249 - t226 + t373;
t136 = t264 * t149;
t389 = pkin(3) * t416 + pkin(7) * t421;
t394 = rSges(5,3) * t421 + t264 * t373;
t503 = t389 + t394 - t136 - t179;
t259 = Icges(5,4) * t269;
t323 = -Icges(5,2) * t267 + t259;
t303 = t323 * t264;
t446 = Icges(5,4) * t267;
t229 = Icges(5,2) * t269 + t446;
t479 = -Icges(5,6) * t264 + qJD(4) * t229;
t102 = -t249 * t479 + t250 * t303;
t232 = Icges(5,1) * t269 - t446;
t489 = Icges(5,1) * t267 + t259;
t476 = -Icges(5,5) * t264 + qJD(4) * t489;
t104 = t232 * t416 - t249 * t476;
t144 = Icges(5,4) * t414 - Icges(5,2) * t415 + Icges(5,6) * t249;
t224 = Icges(5,4) * t415;
t146 = Icges(5,1) * t414 + Icges(5,5) * t249 - t224;
t319 = t144 * t269 + t146 * t267;
t305 = t232 * t249;
t145 = -Icges(5,5) * t250 + t305;
t428 = t145 * t269;
t143 = -Icges(5,6) * t250 + t249 * t323;
t430 = t143 * t267;
t320 = t428 - t430;
t502 = -qJD(4) * t320 - t102 * t269 - t104 * t267 + t264 * t319;
t492 = t249 * rSges(4,1) + t250 * rSges(4,2);
t162 = t492 * t264;
t254 = sin(t265);
t268 = sin(qJ(1));
t261 = t268 * pkin(1);
t491 = pkin(2) * t254 + t261;
t309 = t491 * qJD(1);
t128 = t309 + t162;
t445 = Icges(6,4) * t257;
t198 = Icges(6,1) * t258 - t445;
t304 = t198 * t249;
t134 = -Icges(6,5) * t250 + t304;
t206 = Icges(6,4) * t418;
t135 = Icges(6,1) * t417 + Icges(6,5) * t249 - t206;
t182 = t250 * t263;
t195 = Icges(6,2) * t258 + t445;
t248 = Icges(6,4) * t258;
t322 = -Icges(6,2) * t257 + t248;
t490 = Icges(6,1) * t257 + t248;
t497 = t490 + t322;
t279 = t181 * (-Icges(6,2) * t417 + t135 - t206) - t182 * (-t195 * t249 + t134) + t264 * t497;
t302 = t322 * t249;
t132 = -Icges(6,6) * t250 + t302;
t133 = Icges(6,4) * t417 - Icges(6,2) * t418 + Icges(6,6) * t249;
t471 = t181 * (t250 * t490 + t133) - t182 * (t249 * t490 + t132) + t264 * (t195 - t198);
t501 = t279 * t257 + t258 * t471;
t386 = t489 + t323;
t387 = t229 - t232;
t499 = (t267 * t386 + t269 * t387) * t264;
t496 = 0.2e1 * qJD(4);
t194 = Icges(6,5) * t258 - Icges(6,6) * t257;
t130 = -Icges(6,3) * t250 + t194 * t249;
t432 = t133 * t257;
t321 = -t135 * t258 + t432;
t310 = -t130 + t321;
t495 = t182 * t310;
t425 = t195 * t263;
t486 = -Icges(6,6) * t264 + t425;
t157 = t199 * t249;
t158 = t199 * t250;
t454 = rSges(6,2) * t257;
t456 = rSges(6,1) * t258;
t200 = -t454 + t456;
t243 = t249 * pkin(3);
t185 = -pkin(7) * t250 + t243;
t237 = t250 * t271;
t391 = t249 * t251 + t237;
t125 = -t185 + t391;
t422 = t249 * t258;
t209 = rSges(6,1) * t422;
t423 = t249 * t257;
t137 = -rSges(6,2) * t423 - rSges(6,3) * t250 + t209;
t51 = t137 * t181 + t138 * t182 + qJD(2) + (t125 * t249 - t126 * t250) * qJD(4);
t356 = t250 * t376;
t335 = pkin(4) * t356;
t313 = -t182 * t199 - t335;
t365 = t125 + t137 + t185;
t57 = t264 * t365 + t309 - t313;
t270 = cos(qJ(1));
t255 = cos(t265);
t379 = qJD(1) * t255;
t453 = pkin(1) * qJD(1);
t384 = pkin(2) * t379 + t270 * t453;
t407 = -t126 + t138;
t58 = (t186 + t407) * t264 + t298 + t384;
t485 = -t57 * (-t264 * t157 + t182 * t200) - t51 * (-t157 * t181 - t182 * t158) - t58 * (-t158 * t264 - t181 * t200);
t228 = Icges(5,5) * t269 - Icges(5,6) * t267;
t141 = -Icges(5,3) * t250 + t228 * t249;
t96 = t143 * t269 + t145 * t267;
t484 = qJD(4) * t96 + t102 * t267 - t104 * t269 - t141 * t264;
t193 = Icges(6,5) * t257 + Icges(6,6) * t258;
t483 = -Icges(6,3) * t264 + t193 * t263;
t482 = -Icges(6,5) * t264 + t263 * t490;
t419 = t249 * t269;
t225 = rSges(5,1) * t419;
t420 = t249 * t267;
t148 = -rSges(5,2) * t420 - rSges(5,3) * t250 + t225;
t233 = rSges(5,1) * t267 + rSges(5,2) * t269;
t377 = qJD(4) * t250;
t359 = t233 * t377;
t481 = t264 * (t148 + t185) + t359;
t227 = Icges(5,5) * t267 + Icges(5,6) * t269;
t480 = -Icges(5,3) * t264 + qJD(4) * t227;
t202 = t323 * qJD(4);
t203 = t232 * qJD(4);
t478 = qJD(4) * (t229 * t269 + t267 * t489) + t202 * t267 - t203 * t269 - t227 * t264;
t101 = t249 * t303 + t250 * t479;
t103 = t250 * t476 + t264 * t305;
t142 = Icges(5,5) * t414 - Icges(5,6) * t415 + Icges(5,3) * t249;
t477 = qJD(4) * t319 - t101 * t267 + t103 * t269 - t142 * t264;
t337 = -t198 * t263 + t425;
t338 = t497 * t263;
t474 = -t193 * t264 + t257 * t338 + t258 * t337;
t131 = Icges(6,5) * t417 - Icges(6,6) * t418 + Icges(6,3) * t249;
t342 = t135 * t263 - t250 * t486 - t264 * t302;
t344 = t133 * t263 + t250 * t482 + t264 * t304;
t473 = -t131 * t264 + t257 * t342 + t258 * t344;
t343 = t134 * t263 - t249 * t486 + t322 * t416;
t345 = t132 * t263 - t198 * t416 + t249 * t482;
t472 = -t130 * t264 + t257 * t343 + t258 * t345;
t273 = qJD(1) ^ 2;
t160 = t264 * t181;
t470 = t160 / 0.2e1;
t161 = t264 * t182;
t469 = -t161 / 0.2e1;
t468 = -t181 / 0.2e1;
t467 = t181 / 0.2e1;
t466 = t182 / 0.2e1;
t465 = -t182 / 0.2e1;
t464 = -t249 / 0.2e1;
t463 = -t250 / 0.2e1;
t462 = -t264 / 0.2e1;
t461 = t264 / 0.2e1;
t460 = rSges(5,3) + pkin(7);
t247 = pkin(2) * t255;
t262 = t270 * pkin(1);
t457 = rSges(5,1) * t269;
t455 = rSges(5,2) * t267;
t452 = t264 * t58;
t378 = qJD(4) * t249;
t360 = t233 * t378;
t312 = -t360 + t384;
t76 = (t149 + t186) * t264 + t312;
t451 = t264 * t76;
t427 = t193 * t249;
t93 = -t195 * t418 + t417 * t490 + t427;
t450 = t93 * t264;
t220 = pkin(7) * t416;
t107 = t335 + t220 + (t237 + (-pkin(3) + t251) * t249) * t264;
t411 = t258 * t263;
t368 = rSges(6,2) * t411;
t412 = t257 * t264;
t369 = rSges(6,2) * t412;
t397 = rSges(6,3) * t416 + t249 * t369;
t413 = t257 * t263;
t89 = t250 * t368 + (t250 * t413 + t258 * t421) * rSges(6,1) - t397;
t449 = -t107 - t89;
t408 = t267 * t229;
t424 = t227 * t249;
t111 = -t250 * t408 + t414 * t489 + t424;
t436 = t111 * t264;
t433 = t132 * t257;
t431 = t134 * t258;
t429 = t144 * t267;
t152 = t193 * t250;
t300 = t194 * t264;
t168 = t227 * t250;
t301 = t228 * t264;
t410 = t264 * t267;
t409 = t264 * t271;
t402 = t249 * t489 + t143;
t401 = t250 * t489 + t144;
t400 = -t229 * t249 + t145;
t399 = -Icges(5,2) * t414 + t146 - t224;
t370 = rSges(5,2) * t410;
t395 = rSges(5,3) * t416 + t249 * t370;
t390 = t208 - t210;
t388 = t225 + t243;
t252 = t273 * t262;
t385 = t273 * t247 + t252;
t381 = t247 + t262;
t380 = qJD(1) * t254;
t375 = qJD(4) * t269;
t374 = qJD(4) ^ 2 * t458;
t371 = rSges(6,1) * t413;
t367 = pkin(4) * t376;
t90 = -t249 * t371 + (-t249 * t411 - t250 * t412) * rSges(6,2) + t396;
t366 = t137 * t416 + (-t127 + t90) * t249;
t364 = t264 * t389 + t385;
t362 = t220 + t395;
t361 = t209 + t391;
t357 = t249 * t375;
t354 = t421 / 0.2e1;
t353 = -t416 / 0.2e1;
t352 = pkin(3) + t457;
t350 = t378 / 0.2e1;
t348 = t377 / 0.2e1;
t346 = pkin(4) * t267 + t199;
t341 = -t131 - t433;
t340 = -t131 + t431;
t163 = rSges(4,1) * t416 - rSges(4,2) * t421;
t184 = rSges(4,1) * t250 - rSges(4,2) * t249;
t129 = t184 * t264 + t384;
t190 = rSges(3,1) * t255 - rSges(3,2) * t254;
t327 = -t455 + t457;
t75 = t309 + t481;
t326 = -t249 * t76 + t250 * t75;
t81 = -t133 * t258 - t135 * t257;
t318 = -t146 * t269 + t429;
t317 = t148 * t249 + t149 * t250;
t316 = -t195 * t257 + t258 * t490;
t314 = t269 * t489 - t408;
t311 = t491 * t273;
t308 = qJD(4) * t233;
t122 = t145 * t419;
t66 = -t141 * t250 - t143 * t420 + t122;
t123 = t146 * t419;
t67 = t142 * t250 + t144 * t420 - t123;
t307 = (-t249 * t67 - t250 * t66) * qJD(4);
t124 = t143 * t415;
t68 = -t141 * t249 - t145 * t414 + t124;
t69 = t249 * t142 - t318 * t250;
t306 = (-t249 * t69 - t250 * t68) * qJD(4);
t299 = -pkin(2) * t380 - t268 * t453;
t295 = t152 * t181 - t182 * t427 - t300;
t293 = -t249 * t300 - t250 * t483 + t264 * t321;
t292 = -t250 * t300 + t483 * t249 + (t431 - t433) * t264;
t291 = -t249 * t301 - t250 * t480 + t264 * t318;
t290 = t249 * t480 - t250 * t301 + t264 * t320;
t289 = -t194 * t263 + t264 * t316;
t288 = -t228 * qJD(4) + t264 * t314;
t287 = t267 * t400 + t269 * t402;
t286 = t267 * t399 + t269 * t401;
t13 = t292 * t249 + t250 * t472;
t14 = t293 * t249 - t250 * t473;
t15 = -t249 * t472 + t292 * t250;
t16 = t249 * t473 + t293 * t250;
t112 = t134 * t422;
t62 = -t130 * t250 - t132 * t423 + t112;
t113 = t135 * t422;
t63 = t131 * t250 + t133 * t423 - t113;
t92 = t249 * t316 - t152;
t91 = t92 * t264;
t28 = -t181 * t63 - t182 * t62 + t91;
t114 = t132 * t418;
t64 = -t130 * t249 - t134 * t417 + t114;
t65 = t131 * t249 - t250 * t321;
t29 = -t181 * t65 - t182 * t64 - t450;
t40 = -t257 * t345 + t258 * t343;
t41 = t257 * t344 - t258 * t342;
t44 = t289 * t249 + t250 * t474;
t45 = -t249 * t474 + t289 * t250;
t80 = t132 * t258 + t134 * t257;
t285 = (-t13 * t182 - t14 * t181 + t160 * t64 - t161 * t65 + t264 * t44) * t464 + (t295 * t249 + t250 * t501) * t467 + (-t249 * t501 + t295 * t250) * t466 + (-t15 * t182 - t16 * t181 + t160 * t62 - t161 * t63 + t264 * t45) * t463 + (-t257 * t471 + t258 * t279) * t462 + t28 * t354 + t29 * t353 + ((-t264 * t65 - t13) * t250 + (t264 * t64 - t14) * t249) * t468 + (-t249 * t63 - t250 * t62) * t470 + (-t249 * t65 - t250 * t64) * t469 + ((-t264 * t63 - t15) * t250 + (t264 * t62 - t16) * t249) * t465 + ((-t264 * t81 - t40) * t250 + (t264 * t80 - t41) * t249) * t461;
t283 = -t367 - t368 - t371;
t105 = rSges(5,2) * t250 * t375 + (t264 * t419 + t356) * rSges(5,1) - t395;
t106 = -rSges(5,1) * t358 + (-t250 * t410 - t357) * rSges(5,2) + t394;
t281 = (t148 * t264 - t105) * t250 + (t106 - t136) * t249;
t110 = t249 * t314 - t168;
t109 = t110 * t264;
t32 = t109 + t307;
t33 = t306 - t436;
t50 = qJD(4) * t318 + t101 * t269 + t103 * t267;
t54 = t288 * t249 + t250 * t478;
t55 = -t249 * t478 + t288 * t250;
t276 = (t91 - (t249 * t341 + t112 + t65) * t182 + (t113 - t114 + t64 + (t130 - t432) * t249) * t181 + (t181 * t340 - t495) * t250) * t467 + t81 * t469 + t161 * t93 / 0.2e1 + (t109 + ((t68 + t123 - t124 + (t141 - t429) * t249) * t249 + (-t122 - t69 + (t141 - t318) * t250 + (t428 + t430) * t249) * t250) * qJD(4)) * t350 + (t80 + t92) * t470 + (t450 - (-t114 + t63) * t182 + (-t112 + t62) * t181 + (-t181 * t310 - t182 * t340) * t250 + (-t181 * t341 + t495) * t249 + t29) * t466 + (t40 + t45) * t465 + (t33 + t436 + ((t124 - t67 + (t142 - t428) * t250) * t250 + (-t122 + t66 + (t142 + t430) * t249) * t249) * qJD(4)) * t348 + (qJD(4) * t314 + t202 * t269 + t203 * t267 - t257 * t337 + t258 * t338) * t264 + (t41 + t44 + t28) * t468 - (t50 + t54 + t32) * t378 / 0.2e1 - (t55 - t502) * t377 / 0.2e1 + (t250 * t111 + (t110 + t96) * t249) * qJD(4) * t461;
t164 = pkin(3) * t421 - t220;
t211 = t327 * qJD(4);
t60 = -t211 * t378 - t311 + (-t105 - t164 - t359) * t264;
t61 = t106 * t264 + (t211 * t250 - t233 * t421) * qJD(4) + t364;
t275 = (-t308 * t75 - t352 * t451 - t455 * t61 + t460 * t60) * t249 + (-t308 * t76 + t352 * t60 - t370 * t75 - t460 * t61) * t250;
t178 = t200 * t263;
t36 = -t249 * t374 - t161 * t199 - t178 * t181 - t311 + (-t164 - t335 + t449) * t264;
t108 = t187 + (-t367 - t409) * t249 - t389;
t37 = t250 * t374 - t160 * t199 + t178 * t182 + (t108 + t90 - t336) * t264 + t364;
t274 = (t36 * (rSges(6,3) - t271) - t37 * t454 + t57 * t283 + (t58 * (-t251 - t456) - t57 * t271) * t264) * t249 + (t36 * t456 + t58 * (t283 - t409) - t37 * rSges(6,3) - t57 * t369) * t250;
t174 = t233 * t250;
t173 = t233 * t249;
t121 = t249 * t137;
t120 = t163 * t264 + t385;
t119 = -t162 * t264 - t311;
t82 = qJD(4) * t317 + qJD(2);
t46 = t281 * qJD(4);
t12 = t137 * t161 - t138 * t160 + t181 * t90 - t182 * t89 + ((t125 * t264 - t107) * t250 + (t108 + t118) * t249) * qJD(4);
t1 = [m(4) * (t119 * (t184 + t381) + t120 * (t491 + t492) + (-t129 + t163 + t384) * t128) + m(3) * (-t252 + t273 * (t190 + t262) + (-0.2e1 * rSges(3,1) * t379 + 0.2e1 * rSges(3,2) * t380 + qJD(1) * t190) * qJD(1)) * (-t254 * rSges(3,1) - t255 * rSges(3,2) - t261) + t276 + (t36 * (t381 + t390) + t58 * (t299 + t397) + t37 * (t361 + t491) + t274 + (t58 + t504) * t57) * m(6) + (t60 * (-t226 + t381) + t76 * (t299 + t362) + t61 * (t491 + t388) + t275 + (t384 - t312 + t76 + t503) * t75) * m(5); m(5) * t46 + m(6) * t12; t276 + (t36 * t390 + t37 * t361 + t365 * t452 + t274 + (-t313 + t397) * t58 + t504 * t57) * m(6) + (-t60 * t226 + t61 * t388 + t275 + (t362 + t481) * t76 + (t360 + t503) * t75) * m(5) + (t119 * t184 + t120 * t492 + t128 * t163 - t129 * t162 - (t128 * t184 - t129 * t492) * t264) * m(4); (t502 * t250 + (t264 * t96 - t50) * t249) * t461 + t285 + ((t168 * t378 - t301) * t249 + (t499 + (-t287 * t250 + (-t424 + t286) * t249) * qJD(4)) * t250) * t350 + ((-t377 * t424 - t301) * t250 + (-t499 + (-t286 * t249 + (t168 + t287) * t250) * qJD(4)) * t249) * t348 + ((-t267 * t387 + t269 * t386) * t264 + ((t249 * t399 - t250 * t400) * t269 + (-t249 * t401 + t250 * t402) * t267) * qJD(4)) * t462 + (t264 * t54 + ((-t290 * t249 - t250 * t484 - t264 * t69) * t250 + (-t291 * t249 + t250 * t477 + t264 * t68) * t249) * t496) * t464 + (t264 * t55 + ((t249 * t484 - t290 * t250 - t264 * t67) * t250 + (-t249 * t477 - t291 * t250 + t264 * t66) * t249) * t496) * t463 + (t307 + t32) * t354 + (t306 + t33) * t353 + (t12 * t121 + t51 * t366 + (t12 * t125 + t51 * t108 - t36 * t346 + t58 * (-pkin(4) * t375 - t178) + (t51 * t126 - t346 * t57) * t264) * t249 + (t12 * t407 + t51 * t449 + t37 * t346 + t57 * t178 + (t51 * t125 - t346 * t58) * t264) * t250 - (-t58 * t357 + ((-t249 * t57 - t250 * t58) * t264 + t51 * (-t249 ^ 2 - t250 ^ 2) * qJD(4)) * t267) * pkin(4) + t485) * m(6) + (-(-t173 * t75 - t174 * t76) * t264 - (t82 * (-t173 * t249 - t174 * t250) + t326 * t327) * qJD(4) + t46 * t317 + t82 * t281 + t326 * t211 + ((t61 - t451) * t250 + (-t264 * t75 - t60) * t249) * t233) * m(5); t285 + (t12 * (t138 * t250 + t121) + t51 * (-t250 * t89 + t366) + (-t249 * t58 + t250 * t57) * t178 + ((t37 - t452) * t250 + (-t264 * t57 - t36) * t249) * t199 + t485) * m(6);];
tauc = t1(:);
