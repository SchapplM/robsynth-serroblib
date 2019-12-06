% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:33
% EndTime: 2019-12-05 17:39:50
% DurationCPUTime: 8.51s
% Computational Cost: add. (9549->547), mult. (10482->714), div. (0->0), fcn. (7976->8), ass. (0->315)
t260 = pkin(8) + qJ(4);
t240 = qJ(5) + t260;
t231 = sin(t240);
t232 = cos(t240);
t265 = sin(qJ(1));
t402 = t232 * t265;
t201 = Icges(6,4) * t402;
t403 = t231 * t265;
t266 = cos(qJ(1));
t409 = Icges(6,5) * t266;
t123 = Icges(6,1) * t403 + t201 + t409;
t261 = qJD(4) + qJD(5);
t209 = t261 * t265;
t210 = t261 * t266;
t412 = Icges(6,4) * t231;
t176 = Icges(6,1) * t232 - t412;
t308 = Icges(6,2) * t232 + t412;
t475 = -t308 + t176;
t411 = Icges(6,4) * t232;
t310 = Icges(6,1) * t231 + t411;
t124 = -Icges(6,5) * t265 + t266 * t310;
t174 = -Icges(6,2) * t231 + t411;
t477 = t174 * t266 + t124;
t450 = qJD(1) * t475 - t209 * t477 + t210 * (-Icges(6,2) * t403 + t123 + t201);
t476 = t174 + t310;
t122 = -Icges(6,6) * t265 + t266 * t308;
t478 = -t176 * t266 + t122;
t121 = Icges(6,6) * t266 + t265 * t308;
t479 = -t176 * t265 + t121;
t451 = qJD(1) * t476 - t209 * t478 + t210 * t479;
t481 = t231 * t451 - t232 * t450;
t239 = cos(t260);
t480 = pkin(4) * t239;
t313 = rSges(6,1) * t231 + rSges(6,2) * t232;
t238 = sin(t260);
t314 = rSges(5,1) * t238 + rSges(5,2) * t239;
t400 = t239 * t265;
t216 = Icges(5,4) * t400;
t401 = t238 * t265;
t410 = Icges(5,5) * t266;
t132 = Icges(5,1) * t401 + t216 + t410;
t413 = Icges(5,4) * t239;
t311 = Icges(5,1) * t238 + t413;
t133 = -Icges(5,5) * t265 + t266 * t311;
t185 = -Icges(5,2) * t238 + t413;
t156 = t185 * t266;
t277 = t265 * (t133 + t156) - t266 * (-Icges(5,2) * t401 + t132 + t216);
t414 = Icges(5,4) * t238;
t309 = Icges(5,2) * t239 + t414;
t130 = Icges(5,6) * t266 + t265 * t309;
t131 = -Icges(5,6) * t265 + t266 * t309;
t187 = Icges(5,1) * t239 - t414;
t157 = t187 * t265;
t158 = t187 * t266;
t278 = t265 * (t131 - t158) - t266 * (t130 - t157);
t474 = -t278 * t238 + t277 * t239;
t379 = t185 + t311;
t380 = -t309 + t187;
t473 = (t238 * t379 - t239 * t380) * qJD(1);
t472 = 2 * qJD(4);
t471 = rSges(5,2) * t238;
t306 = Icges(6,5) * t231 + Icges(6,6) * t232;
t120 = -Icges(6,3) * t265 + t266 * t306;
t45 = -t266 * t120 - t122 * t402 - t124 * t403;
t300 = t174 * t232 + t176 * t231;
t172 = Icges(6,5) * t232 - Icges(6,6) * t231;
t405 = t172 * t266;
t58 = t265 * t300 + t405;
t469 = t58 * qJD(1) + t209 * t45;
t426 = rSges(6,2) * t231;
t429 = rSges(6,1) * t232;
t178 = -t426 + t429;
t148 = t178 * t265;
t149 = t178 * t266;
t255 = t266 * rSges(6,3);
t76 = -t261 * t149 + (t265 * t313 + t255) * qJD(1);
t468 = t148 * t209 + t210 * t149 + t266 * t76;
t304 = t122 * t232 + t124 * t231;
t467 = t266 * t304;
t301 = t131 * t239 + t133 * t238;
t464 = t301 * t266;
t214 = t266 * pkin(1) + t265 * qJ(2);
t262 = sin(pkin(8));
t399 = t262 * t265;
t227 = pkin(3) * t399;
t463 = t214 + t227;
t407 = qJ(3) * t266;
t427 = rSges(4,2) * cos(pkin(8));
t460 = -rSges(4,1) * t399 - t266 * rSges(4,3) - t265 * t427;
t462 = t214 + t407 - t460;
t332 = -rSges(3,2) * t266 + t265 * rSges(3,3);
t461 = t214 + t332;
t307 = Icges(5,5) * t238 + Icges(5,6) * t239;
t129 = -Icges(5,3) * t265 + t266 * t307;
t368 = qJD(1) * t129;
t69 = t131 * t238 - t133 * t239;
t86 = qJD(1) * t130 - qJD(4) * t156;
t88 = -qJD(4) * t158 + (t265 * t311 + t410) * qJD(1);
t459 = qJD(4) * t69 + t238 * t88 + t239 * t86 + t368;
t167 = t309 * qJD(4);
t168 = t311 * qJD(4);
t183 = Icges(5,5) * t239 - Icges(5,6) * t238;
t458 = qJD(1) * t183 + qJD(4) * (t185 * t238 - t187 * t239) + t167 * t239 + t168 * t238;
t302 = t130 * t238 - t132 * t239;
t128 = Icges(5,3) * t266 + t265 * t307;
t369 = qJD(1) * t128;
t362 = qJD(4) * t265;
t87 = qJD(1) * t131 + t185 * t362;
t89 = qJD(1) * t133 + qJD(4) * t157;
t457 = qJD(4) * t302 - t238 * t89 - t239 * t87 + t369;
t341 = t239 * t362;
t202 = pkin(4) * t341;
t241 = qJD(3) * t266;
t242 = qJD(2) * t265;
t373 = t241 + t242;
t294 = t209 * t178 + t202 + t373;
t398 = t262 * t266;
t356 = pkin(3) * t398;
t264 = -pkin(6) - qJ(3);
t394 = qJ(3) + t264;
t163 = t265 * t394 + t356;
t245 = t266 * qJ(2);
t211 = pkin(1) * t265 - t245;
t408 = qJ(3) * t265;
t327 = -t211 - t408;
t319 = t163 + t327;
t252 = t265 * rSges(6,3);
t126 = t313 * t266 - t252;
t436 = pkin(4) * t238;
t437 = pkin(3) * t262;
t196 = t436 + t437;
t259 = -pkin(7) + t264;
t381 = -t266 * t196 - t265 * t259;
t96 = t264 * t265 + t356 + t381;
t415 = t126 - t96;
t42 = (t319 + t415) * qJD(1) + t294;
t361 = qJD(4) * t266;
t204 = t361 * t480;
t243 = qJD(2) * t266;
t321 = qJD(3) * t265 - t243;
t317 = -t204 + t321;
t318 = -t266 * t394 + t407 + t463;
t125 = rSges(6,1) * t403 + rSges(6,2) * t402 + t255;
t180 = t265 * t196;
t97 = t180 - t227 + (-t259 + t264) * t266;
t416 = -t125 - t97;
t43 = -t178 * t210 + (t318 - t416) * qJD(1) + t317;
t456 = -t42 * (qJD(1) * t149 - t209 * t313) - t43 * (qJD(1) * t148 + t210 * t313);
t324 = t476 * t261;
t325 = t475 * t261;
t454 = qJD(1) * t172 + t231 * t324 - t232 * t325;
t328 = -qJD(1) * t121 + t477 * t261;
t330 = (t265 * t310 + t409) * qJD(1) + t478 * t261;
t370 = qJD(1) * t120;
t453 = t231 * t330 - t232 * t328 + t370;
t329 = qJD(1) * t122 + t123 * t261 + t174 * t209;
t331 = -qJD(1) * t124 + t479 * t261;
t119 = Icges(6,3) * t266 + t265 * t306;
t371 = qJD(1) * t119;
t452 = t231 * t331 - t232 * t329 + t371;
t449 = t265 ^ 2;
t366 = qJD(1) * t265;
t191 = t261 * t366;
t448 = -t191 / 0.2e1;
t192 = qJD(1) * t210;
t447 = t192 / 0.2e1;
t446 = -t209 / 0.2e1;
t445 = t209 / 0.2e1;
t444 = -t210 / 0.2e1;
t443 = t210 / 0.2e1;
t442 = t265 / 0.2e1;
t441 = -t266 / 0.2e1;
t440 = t266 / 0.2e1;
t439 = rSges(3,2) - pkin(1);
t438 = -rSges(6,3) - pkin(1);
t434 = -qJD(1) / 0.2e1;
t433 = qJD(1) / 0.2e1;
t365 = qJD(1) * t266;
t348 = t209 * t429 + t313 * t365;
t353 = t261 * t426;
t77 = (-rSges(6,3) * qJD(1) - t353) * t265 + t348;
t349 = t196 * t365 + t259 * t366 + t202;
t343 = t262 * t365;
t376 = pkin(3) * t343 + t264 * t366;
t82 = t349 - t376;
t432 = -t77 - t82;
t424 = rSges(3,3) * t266;
t140 = t313 * t261;
t358 = qJD(1) * qJD(3);
t359 = qJD(1) * qJD(2);
t235 = qJ(2) * t365;
t375 = t235 + t242;
t384 = qJD(1) * (-pkin(1) * t366 + t375) + t265 * t359;
t406 = qJ(3) * qJD(1) ^ 2;
t289 = -t265 * t406 + 0.2e1 * t266 * t358 + t384;
t285 = qJD(1) * (qJ(3) * t366 + t376) + t289;
t355 = (qJD(4) ^ 2) * t436;
t22 = t266 * t355 + t140 * t210 + t178 * t191 + (t202 - t432) * qJD(1) + t285;
t422 = t22 * t266;
t234 = t266 * t359;
t296 = -0.2e1 * t265 * t358 - t266 * t406 + t234;
t169 = qJD(1) * t214 - t243;
t224 = t264 * t365;
t385 = t224 - (t227 - t407) * qJD(1) - t169;
t83 = -t204 + t224 + (-t259 * t266 + (t196 - t437) * t265) * qJD(1);
t23 = -t265 * t355 - t140 * t209 + t178 * t192 + (-t76 - t83 + t204 + t385) * qJD(1) + t296;
t421 = t23 * t265;
t170 = t314 * qJD(4);
t189 = rSges(5,1) * t239 - t471;
t342 = t189 * t361;
t160 = t189 * t266;
t256 = t266 * rSges(5,3);
t92 = -qJD(4) * t160 + (t265 * t314 + t256) * qJD(1);
t41 = -t170 * t362 + (-t92 + t342 + t385) * qJD(1) + t296;
t420 = t265 * t41;
t162 = t189 * t362;
t350 = qJD(4) * t471;
t293 = -rSges(5,3) * qJD(1) - t350;
t347 = rSges(5,1) * t341 + t314 * t365;
t93 = t265 * t293 + t347;
t40 = t170 * t361 + (t93 + t162) * qJD(1) + t285;
t419 = t266 * t40;
t418 = t266 * t42;
t417 = t42 * t178;
t404 = t183 * t266;
t397 = t265 * t120;
t142 = t265 * t172;
t153 = t265 * t183;
t59 = t266 * t300 - t142;
t396 = t59 * qJD(1);
t299 = t239 * t185 + t238 * t187;
t66 = t266 * t299 - t153;
t395 = t66 * qJD(1);
t354 = t266 * t427;
t377 = rSges(4,1) * t343 + qJD(1) * t354;
t374 = rSges(3,2) * t366 + rSges(3,3) * t365;
t197 = qJD(1) * t211;
t372 = t242 - t197;
t363 = qJD(4) * t239;
t360 = t307 * qJD(1);
t357 = -rSges(4,3) - pkin(1) - qJ(3);
t44 = t266 * t119 + t121 * t402 + t123 * t403;
t51 = t266 * t128 + t130 * t400 + t132 * t401;
t52 = -t266 * t129 - t131 * t400 - t133 * t401;
t138 = rSges(5,1) * t401 + rSges(5,2) * t400 + t256;
t346 = t235 + t373;
t345 = t241 + t372;
t340 = -t366 / 0.2e1;
t339 = t365 / 0.2e1;
t338 = -t362 / 0.2e1;
t336 = -t361 / 0.2e1;
t334 = t178 + t480;
t315 = rSges(4,1) * t262 + t427;
t253 = t265 * rSges(5,3);
t139 = t266 * t314 - t253;
t55 = t162 + (t139 + t319) * qJD(1) + t373;
t56 = -t342 + (t138 + t318) * qJD(1) + t321;
t312 = t265 * t55 - t266 * t56;
t305 = t121 * t232 + t123 * t231;
t61 = t122 * t231 - t124 * t232;
t303 = t130 * t239 + t132 * t238;
t297 = t44 + t397;
t291 = (t265 * t52 + t266 * t51) * qJD(4);
t115 = t265 * t128;
t53 = -t303 * t266 + t115;
t54 = -t265 * t129 + t464;
t290 = (t265 * t54 + t266 * t53) * qJD(4);
t78 = (-t138 * t265 - t139 * t266) * qJD(4);
t288 = t314 + t437;
t284 = -qJD(1) * t306 + t142 * t210 - t209 * t405;
t109 = t265 * t119;
t46 = -t266 * t305 + t109;
t282 = -qJD(1) * t304 - t261 * t405 + t371;
t281 = qJD(1) * t305 + t142 * t261 + t370;
t280 = -qJD(1) * t301 - qJD(4) * t404 + t369;
t279 = qJD(1) * t303 + qJD(4) * t153 + t368;
t276 = qJD(1) * t300 - t306 * t261;
t275 = t299 * qJD(1) - t307 * qJD(4);
t10 = t282 * t265 - t266 * t453;
t11 = -t265 * t452 + t281 * t266;
t12 = t265 * t453 + t282 * t266;
t20 = t210 * t44 + t469;
t47 = -t397 + t467;
t21 = t209 * t47 + t210 * t46 - t396;
t26 = t276 * t265 + t266 * t454;
t27 = -t265 * t454 + t276 * t266;
t30 = -t231 * t329 - t232 * t331;
t31 = t231 * t328 + t232 * t330;
t60 = -t121 * t231 + t123 * t232;
t9 = t281 * t265 + t266 * t452;
t274 = (qJD(1) * t26 + t10 * t209 - t191 * t46 + t192 * t47 + t210 * t9) * t442 + (-t231 * t450 - t232 * t451) * t434 + t20 * t340 + t21 * t339 + (qJD(1) * t27 + t11 * t210 + t12 * t209 - t191 * t44 + t192 * t45) * t440 + (t265 * t45 + t266 * t44) * t448 + (t265 * t47 + t266 * t46) * t447 + (t10 * t265 + t266 * t9 + (-t265 * t46 + t266 * t47) * qJD(1)) * t445 + (t11 * t266 + t12 * t265 + (-t265 * t44 + t266 * t45) * qJD(1)) * t443 + (t265 * t31 + t266 * t30 + (-t265 * t60 + t266 * t61) * qJD(1)) * t433 + (t284 * t265 + t481 * t266) * t446 + (-t481 * t265 + t284 * t266) * t444;
t212 = rSges(3,2) * t265 + t424;
t159 = t189 * t265;
t152 = -rSges(4,3) * t265 + t266 * t315;
t150 = qJD(1) * t163;
t114 = qJD(1) * t461 - t243;
t113 = t242 + (-t211 + t212) * qJD(1);
t112 = t266 * t126;
t95 = t234 + (-qJD(1) * t332 - t169) * qJD(1);
t94 = qJD(1) * t374 + t384;
t91 = qJD(1) * t462 + t321;
t90 = (t152 + t327) * qJD(1) + t373;
t65 = t265 * t299 + t404;
t64 = (qJD(1) * t460 - t169) * qJD(1) + t296;
t63 = qJD(1) * (-rSges(4,3) * t366 + t377) + t289;
t62 = t65 * qJD(1);
t39 = -t125 * t209 - t126 * t210 + (-t265 * t97 + t266 * t96) * qJD(4);
t35 = t301 * qJD(4) - t238 * t86 + t239 * t88;
t34 = -qJD(4) * t303 - t238 * t87 + t239 * t89;
t33 = -t265 * t458 + t275 * t266;
t32 = t275 * t265 + t266 * t458;
t25 = t290 - t395;
t24 = t62 + t291;
t8 = -t125 * t192 + t126 * t191 - t209 * t77 + t210 * t76 + (-t265 * t82 + t266 * t83 + (-t265 * t96 - t266 * t97) * qJD(1)) * qJD(4);
t1 = [(t62 + ((-t53 + t115 + t52) * t265 + (t54 - t464 + (t129 - t303) * t265 + t51) * t266) * qJD(4)) * t338 + ((t297 + t47 - t467) * t210 + t469) * t446 + t61 * t447 - t192 * t59 / 0.2e1 + (t60 + t58) * t448 + (t396 + (t304 * t265 - t109 + t45) * t210 + (t297 - t44) * t209 + ((t120 + t305) * t210 - t304 * t209) * t266 + t21) * t444 + (t30 + t27) * t443 + (t395 + (t449 * t129 + (-t115 + t52 + (t129 + t303) * t266) * t266) * qJD(4) + t25) * t336 + (-qJD(4) * t299 + t167 * t238 - t168 * t239 - t231 * t325 - t232 * t324) * qJD(1) + (t23 * (-t211 - t252 - t381) - t42 * t317 + t22 * (t180 + t125 + t214) + t43 * (-t265 * t353 + t346 + t348 + t349) + (-t22 * t259 + t23 * t313 + t261 * t417) * t266 + ((t259 + t438) * t418 + (t42 * (-qJ(2) - t196 - t313) + t43 * t438) * t265) * qJD(1) - (t150 - t197 - t42 + (-t408 + t415) * qJD(1) + t294) * t43) * m(6) + ((-t264 * t266 + t138 + t463) * t40 + (t245 - t253 + t288 * t266 + (-pkin(1) + t264) * t265) * t41 + (t224 + t243 + (rSges(5,1) * t363 - pkin(1) * qJD(1) + t293) * t266 + (-qJD(3) + (-qJ(2) - t288) * qJD(1)) * t265) * t55 + (t346 + t347 + t376 + (-t350 + (-rSges(5,3) - pkin(1)) * qJD(1)) * t265 - t150 - t162 + t55 - (t139 - t408) * qJD(1) - t345) * t56) * m(5) + (t64 * (rSges(4,1) * t398 + t245 + t354) + t90 * t243 + t63 * t462 + t91 * (t346 + t377) + (-t90 * qJD(3) + t64 * t357) * t265 + (t90 * t357 * t266 + (t90 * (-qJ(2) - t315) + t91 * t357) * t265) * qJD(1) - (-t90 + (t152 - t408) * qJD(1) + t345) * t91) * m(4) + (t95 * (t265 * t439 + t245 + t424) + t113 * t243 + t94 * t461 + t114 * (t374 + t375) + (t113 * t439 * t266 + (t113 * (-rSges(3,3) - qJ(2)) - t114 * pkin(1)) * t265) * qJD(1) - (qJD(1) * t212 - t113 + t372) * t114) * m(3) + (t31 + t26 + t20) * t445 + (t35 + t32 + t24) * t362 / 0.2e1 + (qJD(1) * t69 + t33 + t34) * t361 / 0.2e1 + (t266 * t66 + (-t302 + t65) * t265) * qJD(4) * t434; 0.2e1 * (-t422 / 0.2e1 + t421 / 0.2e1) * m(6) + 0.2e1 * (t420 / 0.2e1 - t419 / 0.2e1) * m(5) + 0.2e1 * (t441 * t63 + t442 * t64) * m(4) + 0.2e1 * (t441 * t94 + t442 * t95) * m(3); m(4) * (t265 * t63 + t266 * t64) + m(5) * (t265 * t40 + t266 * t41) + m(6) * (t22 * t265 + t23 * t266); ((t153 * t361 - t360) * t266 + (-t473 + (-t266 * t404 - t474) * qJD(4)) * t265) * t336 + ((-t362 * t404 - t360) * t265 + (t473 + (t265 * t153 + t474) * qJD(4)) * t266) * t338 + t274 + (t265 * t35 + t266 * t34 + (t265 * t302 + t69 * t266) * qJD(1)) * t433 + ((-t238 * t380 - t239 * t379) * qJD(1) + (t238 * t277 + t239 * t278) * qJD(4)) * t434 + (qJD(1) * t32 + ((t279 * t265 + t266 * t457) * t266 + (t280 * t265 - t266 * t459) * t265 + (-t53 * t265 + t54 * t266) * qJD(1)) * t472) * t442 + (qJD(1) * t33 + ((-t265 * t457 + t279 * t266) * t266 + (t265 * t459 + t280 * t266) * t265 + (-t51 * t265 + t52 * t266) * qJD(1)) * t472) * t440 + (t291 + t24) * t340 + (t25 + t290) * t339 + (-t8 * t112 + (t417 * qJD(1) + t43 * t140 - t22 * t334 + t8 * t96) * t266 + (t43 * qJD(1) * t178 - t42 * t140 + t23 * t334 + t8 * t416) * t265 + ((qJD(1) * t416 + t83) * t266 + (t415 * qJD(1) + t432) * t265 - (-t266 ^ 2 - t449) * pkin(4) * t363 + t468) * t39 + t456) * m(6) + (-(t159 * t56 + t160 * t55) * qJD(1) - (t78 * (-t159 * t265 - t160 * t266) - t312 * t314) * qJD(4) + 0.2e1 * t78 * (-t265 * t93 + t266 * t92 + (-t138 * t266 + t139 * t265) * qJD(1)) - t312 * t170 + (t420 - t419 + (t265 * t56 + t266 * t55) * qJD(1)) * t189) * m(5); t274 + (t8 * (-t125 * t265 - t112) - (t265 * t42 - t266 * t43) * t140 + (-t422 + t421 + (t265 * t43 + t418) * qJD(1)) * t178 + (-t265 * t77 + (-t125 * t266 + t126 * t265) * qJD(1) + t468) * t39 + t456) * m(6);];
tauc = t1(:);
