% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPPR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPPR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:58
% EndTime: 2019-12-05 15:24:27
% DurationCPUTime: 18.11s
% Computational Cost: add. (11900->571), mult. (13646->890), div. (0->0), fcn. (12665->10), ass. (0->323)
t251 = qJ(2) + pkin(8);
t245 = sin(t251);
t247 = cos(t251);
t258 = sin(qJ(2));
t259 = cos(qJ(2));
t474 = -Icges(3,5) * t258 - Icges(4,5) * t245 - Icges(3,6) * t259 - Icges(4,6) * t247;
t250 = pkin(9) + qJ(5);
t244 = sin(t250);
t246 = cos(t250);
t347 = rSges(6,1) * t246 - rSges(6,2) * t244;
t144 = -rSges(6,3) * t247 + t245 * t347;
t255 = cos(pkin(7));
t384 = qJD(5) * t245;
t253 = sin(pkin(7));
t391 = qJD(2) * t253;
t221 = t255 * t384 + t391;
t477 = t144 * t221;
t390 = qJD(2) * t255;
t222 = t253 * t384 - t390;
t476 = t144 * t222;
t248 = t253 ^ 2;
t249 = t255 ^ 2;
t394 = t248 + t249;
t475 = t474 * qJD(2);
t252 = sin(pkin(9));
t254 = cos(pkin(9));
t424 = Icges(3,4) * t258;
t328 = -Icges(3,2) * t259 - t424;
t423 = Icges(3,4) * t259;
t333 = -Icges(3,1) * t258 - t423;
t392 = qJD(2) * t247;
t393 = qJD(2) * t245;
t421 = Icges(4,4) * t247;
t422 = Icges(4,4) * t245;
t473 = -(-t258 * t328 + t259 * t333) * qJD(2) - (-Icges(4,1) * t245 - t421) * t392 + (-(-Icges(5,5) * t254 + Icges(5,6) * t252) * t245 - t422 + (-Icges(5,3) - Icges(4,2)) * t247) * t393;
t453 = -rSges(5,1) * t254 + rSges(5,2) * t252;
t451 = rSges(5,3) * t247 + t245 * t453;
t472 = t475 * t253;
t471 = t475 * t255;
t470 = t474 * t253;
t469 = t474 * t255;
t406 = t253 * t254;
t407 = t252 * t255;
t202 = -t247 * t407 + t406;
t405 = t254 * t255;
t408 = t252 * t253;
t203 = t247 * t405 + t408;
t411 = t245 * t255;
t339 = t252 * (Icges(5,4) * t203 + Icges(5,2) * t202 + Icges(5,6) * t411) - t254 * (Icges(5,1) * t203 + Icges(5,4) * t202 + Icges(5,5) * t411);
t200 = -t247 * t408 - t405;
t201 = t247 * t406 - t407;
t412 = t245 * t253;
t340 = t252 * (Icges(5,4) * t201 + Icges(5,2) * t200 + Icges(5,6) * t412) - t254 * (Icges(5,1) * t201 + Icges(5,4) * t200 + Icges(5,5) * t412);
t327 = -Icges(4,2) * t245 + t421;
t332 = Icges(4,1) * t247 - t422;
t460 = -(Icges(4,6) * t253 + t255 * t327) * t247 - (Icges(4,5) * t253 + t255 * t332) * t245;
t461 = (-Icges(4,6) * t255 + t253 * t327) * t247 + (-Icges(4,5) * t255 + t253 * t332) * t245;
t77 = Icges(5,5) * t201 + Icges(5,6) * t200 + Icges(5,3) * t412;
t78 = Icges(5,5) * t203 + Icges(5,6) * t202 + Icges(5,3) * t411;
t468 = (t253 * t339 - t255 * t340) * t245 + (t253 * t78 - t255 * t77) * t247 + t253 * t460 + t255 * t461;
t329 = -Icges(3,2) * t258 + t423;
t184 = Icges(3,6) * t253 + t255 * t329;
t334 = Icges(3,1) * t259 - t424;
t186 = Icges(3,5) * t253 + t255 * t334;
t465 = t473 * t255 + (t184 * t259 + t186 * t258 - t245 * t339 - t247 * t78 - t460) * qJD(2);
t183 = -Icges(3,6) * t255 + t253 * t329;
t185 = -Icges(3,5) * t255 + t253 * t334;
t464 = t473 * t253 + (t183 * t259 + t185 * t258 - t245 * t340 - t247 * t77 + t461) * qJD(2);
t386 = qJD(4) * t253;
t232 = t245 * t386;
t225 = pkin(3) * t245 - qJ(4) * t247;
t300 = qJD(2) * t225;
t131 = -t253 * t300 + t232;
t385 = qJD(4) * t255;
t234 = t245 * t385;
t132 = -t255 * t300 + t234;
t440 = pkin(2) * t258;
t357 = t394 * t440;
t429 = pkin(2) * qJD(2);
t379 = t258 * t429;
t388 = qJD(3) * t255;
t223 = -t253 * t379 - t388;
t243 = qJD(3) * t253;
t224 = -t255 * t379 + t243;
t396 = t253 * t223 + t255 * t224;
t463 = t357 * qJD(2) - qJD(4) * t245 + t253 * t131 + t255 * t132 + t300 * t394 + t396;
t462 = t394 * qJD(2) * t451;
t444 = t221 / 0.2e1;
t442 = t222 / 0.2e1;
t459 = qJD(5) / 0.2e1;
t236 = rSges(3,1) * t258 + rSges(3,2) * t259;
t302 = qJD(2) * t236;
t322 = Icges(6,5) * t246 - Icges(6,6) * t244;
t135 = -Icges(6,3) * t247 + t245 * t322;
t417 = Icges(6,4) * t246;
t325 = -Icges(6,2) * t244 + t417;
t137 = -Icges(6,6) * t247 + t245 * t325;
t418 = Icges(6,4) * t244;
t330 = Icges(6,1) * t246 - t418;
t139 = -Icges(6,5) * t247 + t245 * t330;
t404 = t255 * t244;
t181 = t246 * t253 - t247 * t404;
t409 = t247 * t255;
t182 = t244 * t253 + t246 * t409;
t419 = Icges(6,4) * t182;
t69 = Icges(6,2) * t181 + Icges(6,6) * t411 + t419;
t163 = Icges(6,4) * t181;
t71 = Icges(6,1) * t182 + Icges(6,5) * t411 + t163;
t342 = -t244 * t69 + t246 * t71;
t410 = t247 * t253;
t179 = -t244 * t410 - t246 * t255;
t180 = t246 * t410 - t404;
t420 = Icges(6,4) * t180;
t68 = Icges(6,2) * t179 + Icges(6,6) * t412 + t420;
t162 = Icges(6,4) * t179;
t70 = Icges(6,1) * t180 + Icges(6,5) * t412 + t162;
t343 = -t244 * t68 + t246 * t70;
t450 = -(-t135 * t255 - t342) * t221 - (-t135 * t253 - t343) * t222;
t257 = -pkin(6) - qJ(4);
t403 = qJ(4) + t257;
t241 = pkin(4) * t254 + pkin(3);
t437 = pkin(3) - t241;
t449 = t245 * t437 - t247 * t403;
t448 = -t258 * (t328 * t253 + t185) - t259 * (-t333 * t253 + t183);
t165 = (-Icges(6,2) * t246 - t418) * t245;
t383 = qJD(5) * t247;
t269 = t221 * (-Icges(6,2) * t182 + t163 + t71) + t222 * (-Icges(6,2) * t180 + t162 + t70) - t383 * (t139 + t165);
t260 = qJD(2) ^ 2;
t372 = t247 * t390;
t374 = t245 * t390;
t98 = -qJD(5) * t182 + t244 * t374;
t99 = qJD(5) * t181 - t246 * t374;
t48 = Icges(6,5) * t99 + Icges(6,6) * t98 + Icges(6,3) * t372;
t67 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t411;
t308 = t245 * t48 + t392 * t67;
t50 = Icges(6,4) * t99 + Icges(6,2) * t98 + Icges(6,6) * t372;
t52 = Icges(6,1) * t99 + Icges(6,4) * t98 + Icges(6,5) * t372;
t10 = t181 * t50 + t182 * t52 + t255 * t308 + t69 * t98 + t71 * t99;
t136 = Icges(6,3) * t245 + t247 * t322;
t164 = (-Icges(6,5) * t244 - Icges(6,6) * t246) * t245;
t62 = qJD(2) * t136 + qJD(5) * t164;
t307 = t135 * t392 + t245 * t62;
t26 = t181 * t69 + t182 * t71 + t411 * t67;
t66 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t412;
t25 = t181 * t68 + t182 * t70 + t411 * t66;
t427 = t25 * t253;
t341 = t255 * t26 + t427;
t40 = t135 * t411 + t137 * t181 + t139 * t182;
t425 = t40 * t245;
t138 = Icges(6,6) * t245 + t247 * t325;
t63 = qJD(2) * t138 + qJD(5) * t165;
t140 = Icges(6,5) * t245 + t247 * t330;
t166 = (-Icges(6,1) * t244 - t417) * t245;
t64 = qJD(2) * t140 + qJD(5) * t166;
t267 = -(t137 * t98 + t139 * t99 + t181 * t63 + t182 * t64 + t255 * t307) * t247 + (t247 * t341 + t425) * qJD(2);
t373 = t247 * t391;
t375 = t245 * t391;
t96 = -qJD(5) * t180 + t244 * t375;
t97 = qJD(5) * t179 - t246 * t375;
t47 = Icges(6,5) * t97 + Icges(6,6) * t96 + Icges(6,3) * t373;
t309 = t245 * t47 + t392 * t66;
t49 = Icges(6,4) * t97 + Icges(6,2) * t96 + Icges(6,6) * t373;
t51 = Icges(6,1) * t97 + Icges(6,4) * t96 + Icges(6,5) * t373;
t9 = t181 * t49 + t182 * t51 + t255 * t309 + t68 * t98 + t70 * t99;
t447 = t10 * t444 + t267 * t459 + t442 * t9;
t446 = -t394 * t393 / 0.2e1;
t445 = -t221 / 0.2e1;
t443 = -t222 / 0.2e1;
t441 = -t247 / 0.2e1;
t439 = pkin(2) * t259;
t24 = t179 * t69 + t180 * t71 + t412 * t67;
t428 = t24 * t255;
t39 = t135 * t412 + t137 * t179 + t139 * t180;
t426 = t39 * t245;
t413 = t135 * t247;
t128 = -t245 * t403 - t247 * t437;
t112 = t128 * qJD(2);
t227 = pkin(3) * t247 + qJ(4) * t245;
t387 = qJD(4) * t247;
t167 = qJD(2) * t227 - t387;
t402 = -t112 - t167;
t150 = -qJ(3) * t255 + t253 * t439;
t151 = qJ(3) * t253 + t255 * t439;
t400 = t253 * t150 + t255 * t151;
t397 = t223 * t391 + t224 * t390;
t395 = t234 + t243;
t381 = qJD(2) * qJD(4);
t380 = t260 * t439;
t378 = t259 * t429;
t371 = t150 * t391 + t151 * t390 + qJD(1);
t370 = t247 * t381;
t368 = t390 / 0.2e1;
t367 = -t383 / 0.2e1;
t366 = t383 / 0.2e1;
t365 = -t225 - t440;
t226 = rSges(4,1) * t245 + rSges(4,2) * t247;
t364 = -t226 - t440;
t363 = -t227 - t439;
t228 = rSges(4,1) * t247 - rSges(4,2) * t245;
t362 = -t228 - t439;
t361 = qJD(2) * t459;
t360 = t232 - t388;
t196 = t227 * t253;
t197 = t227 * t255;
t358 = t253 * t196 + t255 * t197 + t400;
t356 = t365 + t449;
t355 = -t128 + t363;
t354 = t365 + t451;
t149 = rSges(5,3) * t245 - t247 * t453;
t353 = -t149 + t363;
t352 = t245 * t361;
t351 = t247 * t361;
t350 = -t167 - t378;
t214 = t228 * qJD(2);
t349 = -t214 - t378;
t237 = rSges(3,1) * t259 - rSges(3,2) * t258;
t230 = t253 * t370;
t54 = rSges(6,1) * t99 + rSges(6,2) * t98 + rSges(6,3) * t372;
t145 = rSges(6,3) * t245 + t247 * t347;
t168 = (-rSges(6,1) * t244 - rSges(6,2) * t246) * t245;
t65 = qJD(2) * t145 + qJD(5) * t168;
t73 = rSges(6,1) * t182 + rSges(6,2) * t181 + rSges(6,3) * t411;
t17 = -t253 * t380 - t54 * t383 - t221 * t65 + t230 + (t402 * t253 + (-t144 * t409 + t245 * t73) * qJD(5)) * qJD(2);
t231 = t255 * t370;
t53 = rSges(6,1) * t97 + rSges(6,2) * t96 + rSges(6,3) * t373;
t72 = rSges(6,1) * t180 + rSges(6,2) * t179 + rSges(6,3) * t412;
t18 = -t255 * t380 + t53 * t383 + t222 * t65 + t231 + (t402 * t255 + (t144 * t410 - t245 * t72) * qJD(5)) * qJD(2);
t346 = -t17 * t255 + t18 * t253;
t345 = t67 * t221 + t66 * t222;
t23 = t179 * t68 + t180 * t70 + t412 * t66;
t344 = t23 * t253 + t428;
t29 = t245 * t343 - t247 * t66;
t30 = t245 * t342 - t247 * t67;
t338 = t29 * t253 + t30 * t255;
t337 = -t253 * t73 + t255 * t72;
t336 = qJD(2) * t364;
t335 = t131 * t391 + t132 * t390 + t245 * t381 + t397;
t321 = -t137 * t244 + t139 * t246;
t320 = t394 * t237;
t319 = t394 * t302;
t318 = -t144 + t356;
t317 = t253 * t351;
t316 = t255 * t351;
t141 = t149 * qJD(2);
t315 = -t141 + t350;
t312 = qJD(2) * t356;
t311 = qJD(2) * t354;
t310 = -qJD(2) * t214 - t380;
t306 = t136 - t321;
t292 = t196 * t391 + t197 * t390 + t371 - t387;
t74 = -pkin(4) * t407 + t128 * t253;
t75 = pkin(4) * t408 + t128 * t255;
t16 = t221 * t72 - t222 * t73 + (t253 * t74 + t255 * t75) * qJD(2) + t292;
t305 = t16 * t337;
t160 = -rSges(4,3) * t255 + t228 * t253;
t161 = rSges(4,3) * t253 + t228 * t255;
t45 = (t160 * t253 + t161 * t255) * qJD(2) + t371;
t304 = t45 * t226;
t303 = -t112 + t350 - t65;
t301 = qJD(2) * t226;
t293 = -t380 + (-t141 - t167) * qJD(2);
t291 = (t333 * t255 - t184) * t259 + (-t328 * t255 - t186) * t258;
t290 = t164 * t383 - t221 * (Icges(6,5) * t181 - Icges(6,6) * t182) - t222 * (Icges(6,5) * t179 - Icges(6,6) * t180);
t288 = qJD(2) * t449;
t287 = Icges(5,5) * t247 + (-Icges(5,1) * t254 + Icges(5,4) * t252) * t245;
t285 = Icges(5,6) * t247 + (-Icges(5,4) * t254 + Icges(5,2) * t252) * t245;
t281 = t245 * t290;
t85 = rSges(5,1) * t201 + rSges(5,2) * t200 + rSges(5,3) * t412;
t86 = rSges(5,1) * t203 + rSges(5,2) * t202 + rSges(5,3) * t411;
t27 = (t253 * t85 + t255 * t86) * qJD(2) + t292;
t280 = t27 * t451;
t278 = t16 * (-t241 * t245 - t247 * t257 + t225);
t277 = qJD(2) * t287;
t276 = qJD(2) * t285;
t270 = (Icges(6,1) * t181 - t419 - t69) * t221 + (Icges(6,1) * t179 - t420 - t68) * t222 - (-t137 + t166) * t383;
t268 = -(t137 * t96 + t139 * t97 + t179 * t63 + t180 * t64 + t253 * t307) * t247 + (t247 * t344 + t426) * qJD(2);
t46 = t245 * t321 - t413;
t266 = -((qJD(2) * t321 - t62) * t247 + (qJD(2) * t135 - t244 * t63 + t246 * t64 + (-t137 * t246 - t139 * t244) * qJD(5)) * t245) * t247 + (t46 * t245 + t247 * t338) * qJD(2);
t261 = (-t306 * t383 - t450) * t245;
t235 = t247 * t385;
t233 = t247 * t386;
t178 = t255 * t301;
t177 = t253 * t301;
t134 = t310 * t255;
t133 = t310 * t253;
t130 = t255 * t336 + t243;
t129 = t253 * t336 - t388;
t126 = t287 * t255;
t125 = t287 * t253;
t124 = t285 * t255;
t123 = t285 * t253;
t111 = t139 * t255;
t110 = t139 * t253;
t109 = t137 * t255;
t108 = t137 * t253;
t105 = t255 * t277;
t104 = t253 * t277;
t103 = t255 * t276;
t102 = t253 * t276;
t95 = t319 * qJD(2);
t94 = rSges(6,1) * t181 - rSges(6,2) * t182;
t93 = rSges(6,1) * t179 - rSges(6,2) * t180;
t84 = t255 * t288;
t83 = t253 * t288;
t76 = qJD(2) * t320 + qJD(1);
t61 = t255 * t311 + t395;
t60 = t253 * t311 + t360;
t57 = t255 * t293 + t231;
t56 = t253 * t293 + t230;
t55 = (-t177 * t253 - t178 * t255) * qJD(2) + t397;
t34 = t255 * t312 + t383 * t72 + t395 + t476;
t33 = t253 * t312 - t383 * t73 + t360 - t477;
t28 = qJD(2) * t462 + t335;
t12 = t221 * t30 + t222 * t29 - t383 * t46;
t11 = t221 * t53 - t222 * t54 + (t253 * t83 + t255 * t84 + t337 * t383) * qJD(2) + t335;
t8 = t179 * t50 + t180 * t52 + t253 * t308 + t69 * t96 + t71 * t97;
t7 = t179 * t49 + t180 * t51 + t253 * t309 + t68 * t96 + t70 * t97;
t6 = t221 * t26 + t222 * t25 - t383 * t40;
t5 = t221 * t24 + t222 * t23 - t383 * t39;
t4 = (qJD(2) * t342 - t48) * t247 + (qJD(2) * t67 - t244 * t50 + t246 * t52 + (-t244 * t71 - t246 * t69) * qJD(5)) * t245;
t3 = (qJD(2) * t343 - t47) * t247 + (qJD(2) * t66 - t244 * t49 + t246 * t51 + (-t244 * t70 - t246 * t68) * qJD(5)) * t245;
t1 = qJD(5) * t268 + t221 * t8 + t222 * t7;
t2 = [-m(3) * t95 + m(4) * t55 + m(5) * t28 + m(6) * t11; -t12 * t384 / 0.2e1 + t253 * t447 + (t253 * t8 - t255 * t7) * t442 + ((-t109 * t179 - t111 * t180) * t221 + (-t108 * t179 - t110 * t180) * t222 + (t426 + (-t138 * t179 - t140 * t180 + t428) * t247) * qJD(5) + (((t23 - t413) * qJD(5) + t345) * t247 + t261) * t253) * t443 + (t10 * t253 - t255 * t9) * t444 + ((-t109 * t181 - t111 * t182) * t221 + (-t108 * t181 - t110 * t182) * t222 + (t425 + (-t138 * t181 - t140 * t182 + t427) * t247) * qJD(5) + (((t26 - t413) * qJD(5) + t345) * t247 + t261) * t255) * t445 + (-t25 * t255 + t253 * t26) * t316 + (-t23 * t255 + t24 * t253) * t317 + (t253 * t30 - t255 * t29) * t352 - t255 * t1 / 0.2e1 + (((t109 * t244 - t111 * t246 + t67) * t221 + (t108 * t244 - t110 * t246 + t66) * t222 + t46 * qJD(5)) * t245 + ((t306 * t247 + (t138 * t244 - t140 * t246 - t135) * t245 + t338) * qJD(5) + t450) * t247) * t366 + (t11 * t358 + (t18 * t318 + t34 * t303 + t11 * (t73 + t75)) * t255 + (t17 * t318 + t33 * t303 + t11 * (t72 + t74)) * t253 - t34 * (t145 * t222 + t235) - t33 * (-t145 * t221 + t233) - ((t33 * t73 - t34 * t72) * t245 + t305 * t247) * qJD(5) - ((t255 * t278 + t34 * t355) * t255 + (t253 * t278 + t33 * t355) * t253) * qJD(2) + (t463 + (t54 + t84 - t476) * t255 + (t53 + t83 + t477) * t253) * t16) * m(6) + (t28 * t358 + (t28 * t86 + t315 * t61 + t354 * t57) * t255 + (t28 * t85 + t315 * t60 + t354 * t56) * t253 - t61 * t235 - t60 * t233 - ((t255 * t280 + t353 * t61) * t255 + (t253 * t280 + t353 * t60) * t253) * qJD(2) + (t462 + t463) * t27) * m(5) + (-(-t45 * t357 + (t130 * t362 - t255 * t304) * t255 + (t129 * t362 - t253 * t304) * t253) * qJD(2) + t55 * t400 + t45 * t396 + (t130 * t349 + t134 * t364 + t55 * t161 - t45 * t178) * t255 + (t129 * t349 + t133 * t364 + t55 * t160 - t45 * t177) * t253) * m(4) + ((t103 * t202 + t105 * t203 + t471 * t253) * t253 + (-t102 * t202 - t104 * t203 + t464 * t255 + (-t465 - t472) * t253) * t255) * t391 + ((t102 * t200 + t104 * t201 - t472 * t255) * t255 + (-t103 * t200 - t105 * t201 + t465 * t253 + (-t464 + t471) * t255) * t253) * t390 - ((t124 * t202 + t126 * t203) * t391 + t469 * t248 * qJD(2) + (-t202 * t123 - t203 * t125 - t448 * t255 + (t291 - t470) * t253 + t468) * t390) * t391 / 0.2e1 + (-(t123 * t200 + t125 * t201) * t390 + t470 * qJD(2) * t249 + (t200 * t124 + t201 * t126 + t291 * t253 + (-t448 - t469) * t255 + t468) * t391) * t368 + ((-t3 + t6) * t255 + (t4 + t5) * t253) * t367 + (-t76 * t319 - t95 * t320 + (t236 * t237 * t260 + t302 * t76) * t394) * m(3); m(4) * (-t133 * t255 + t134 * t253) + m(5) * (t253 * t57 - t255 * t56) + m(6) * t346; 0.2e1 * (t11 * t441 + t16 * t446) * m(6) + 0.2e1 * (t27 * t446 + t28 * t441) * m(5) + 0.2e1 * (m(5) * (qJD(2) * t27 + t253 * t56 + t255 * t57) / 0.2e1 + m(6) * (qJD(2) * t16 + t17 * t253 + t18 * t255) / 0.2e1) * t245; t247 * t6 * t368 + t411 * t447 + (t245 * t341 - t247 * t40) * t316 + ((t10 * t255 + t253 * t9) * t245 + t267) * t444 + t5 * t373 / 0.2e1 + t1 * t412 / 0.2e1 + (t245 * t344 - t247 * t39) * t317 + ((t253 * t7 + t255 * t8) * t245 + t268) * t442 + t12 * t393 / 0.2e1 + (qJD(5) * t266 + t221 * t4 + t222 * t3) * t441 + (t245 * t338 - t247 * t46) * t352 + ((t253 * t3 + t255 * t4) * t245 + t266) * t367 + (t181 * t269 + t182 * t270 - t255 * t281) * t445 + (t179 * t269 + t180 * t270 - t253 * t281) * t443 + (t290 * t247 + (-t244 * t269 + t270 * t246) * t245) * t366 + ((-t17 * t73 + t18 * t72 - t33 * t54 + t34 * t53 + (t305 + (t253 * t34 - t255 * t33) * t144) * qJD(2)) * t247 + (t34 * (-qJD(2) * t72 + t253 * t65) + t33 * (qJD(2) * t73 - t255 * t65) + t11 * t337 + t16 * (-t253 * t54 + t255 * t53) + t346 * t144) * t245 - t34 * (t168 * t222 + t383 * t93) - t33 * (-t168 * t221 - t383 * t94) - t16 * (t221 * t93 - t222 * t94)) * m(6);];
tauc = t2(:);
