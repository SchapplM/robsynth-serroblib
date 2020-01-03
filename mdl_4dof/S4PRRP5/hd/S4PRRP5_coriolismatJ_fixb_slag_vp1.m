% Calculate matrix of centrifugal and coriolis load on the joints for
% S4PRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4PRRP5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP5_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:28:52
% EndTime: 2019-12-31 16:29:02
% DurationCPUTime: 8.24s
% Computational Cost: add. (9994->466), mult. (24411->680), div. (0->0), fcn. (26314->6), ass. (0->292)
t247 = sin(pkin(6));
t248 = cos(pkin(6));
t252 = cos(qJ(3));
t250 = sin(qJ(3));
t253 = cos(qJ(2));
t349 = t250 * t253;
t216 = -t247 * t349 - t248 * t252;
t341 = t252 * t253;
t352 = t248 * t250;
t217 = t247 * t341 - t352;
t251 = sin(qJ(2));
t354 = t247 * t251;
t121 = Icges(5,5) * t217 + Icges(5,6) * t216 + Icges(5,3) * t354;
t123 = Icges(4,5) * t217 + Icges(4,6) * t216 + Icges(4,3) * t354;
t469 = t121 + t123;
t218 = t247 * t252 - t248 * t349;
t355 = t247 * t250;
t219 = t248 * t341 + t355;
t351 = t248 * t251;
t122 = Icges(5,5) * t219 + Icges(5,6) * t218 + Icges(5,3) * t351;
t124 = Icges(4,5) * t219 + Icges(4,6) * t218 + Icges(4,3) * t351;
t468 = t122 + t124;
t285 = Icges(5,5) * t252 - Icges(5,6) * t250;
t189 = -Icges(5,3) * t253 + t251 * t285;
t286 = Icges(4,5) * t252 - Icges(4,6) * t250;
t191 = -Icges(4,3) * t253 + t251 * t286;
t467 = t189 + t191;
t374 = Icges(5,4) * t252;
t288 = -Icges(5,2) * t250 + t374;
t193 = -Icges(5,6) * t253 + t251 * t288;
t378 = Icges(4,4) * t252;
t289 = -Icges(4,2) * t250 + t378;
t195 = -Icges(4,6) * t253 + t251 * t289;
t454 = t193 + t195;
t375 = Icges(5,4) * t250;
t292 = Icges(5,1) * t252 - t375;
t197 = -Icges(5,5) * t253 + t251 * t292;
t379 = Icges(4,4) * t250;
t293 = Icges(4,1) * t252 - t379;
t199 = -Icges(4,5) * t253 + t251 * t293;
t453 = t197 + t199;
t410 = t253 ^ 2;
t342 = t251 * t253;
t140 = Icges(5,5) * t216 - Icges(5,6) * t217;
t142 = Icges(4,5) * t216 - Icges(4,6) * t217;
t466 = t140 + t142;
t141 = Icges(5,5) * t218 - Icges(5,6) * t219;
t143 = Icges(4,5) * t218 - Icges(4,6) * t219;
t465 = t141 + t143;
t213 = Icges(4,4) * t218;
t132 = Icges(4,1) * t219 + Icges(4,5) * t351 + t213;
t332 = -Icges(4,2) * t219 + t132 + t213;
t211 = Icges(5,4) * t218;
t130 = Icges(5,1) * t219 + Icges(5,5) * t351 + t211;
t334 = -Icges(5,2) * t219 + t130 + t211;
t463 = t332 + t334;
t212 = Icges(4,4) * t216;
t131 = Icges(4,1) * t217 + Icges(4,5) * t354 + t212;
t333 = -Icges(4,2) * t217 + t131 + t212;
t210 = Icges(5,4) * t216;
t129 = Icges(5,1) * t217 + Icges(5,5) * t354 + t210;
t335 = -Icges(5,2) * t217 + t129 + t210;
t462 = t333 + t335;
t380 = Icges(4,4) * t219;
t128 = Icges(4,2) * t218 + Icges(4,6) * t351 + t380;
t336 = Icges(4,1) * t218 - t128 - t380;
t376 = Icges(5,4) * t219;
t126 = Icges(5,2) * t218 + Icges(5,6) * t351 + t376;
t338 = Icges(5,1) * t218 - t126 - t376;
t461 = t336 + t338;
t381 = Icges(4,4) * t217;
t127 = Icges(4,2) * t216 + Icges(4,6) * t354 + t381;
t337 = Icges(4,1) * t216 - t127 - t381;
t377 = Icges(5,4) * t217;
t125 = Icges(5,2) * t216 + Icges(5,6) * t354 + t377;
t339 = Icges(5,1) * t216 - t125 - t377;
t460 = t337 + t339;
t440 = rSges(5,1) + pkin(3);
t458 = t125 + t127;
t457 = t126 + t128;
t456 = t129 + t131;
t455 = t130 + t132;
t452 = -t454 * t250 + t453 * t252;
t451 = t467 * t251;
t450 = t468 * t251;
t449 = t469 * t251;
t448 = t462 * t216 + t460 * t217 + t354 * t466;
t447 = t216 * t463 + t217 * t461 + t354 * t465;
t446 = t462 * t218 + t460 * t219 + t351 * t466;
t445 = t218 * t463 + t219 * t461 + t351 * t465;
t444 = t453 + (-t375 - t379 + (-Icges(4,2) - Icges(5,2)) * t252) * t251;
t443 = t454 + (t374 + t378 + (Icges(4,1) + Icges(5,1)) * t250) * t251;
t226 = (-Icges(5,5) * t250 - Icges(5,6) * t252) * t251;
t227 = (-Icges(4,5) * t250 - Icges(4,6) * t252) * t251;
t442 = (-t226 - t227) * t253;
t423 = t467 * t253;
t441 = 0.2e1 * (Icges(3,1) - Icges(3,2)) * t342 + (-0.2e1 * t251 ^ 2 + 0.2e1 * t410) * Icges(3,4);
t168 = t193 * t247;
t170 = t195 * t247;
t172 = t197 * t247;
t174 = t199 * t247;
t282 = -t127 * t250 + t131 * t252;
t272 = -t191 * t247 - t282;
t284 = -t125 * t250 + t129 * t252;
t274 = -t189 * t247 - t284;
t439 = (-t274 - t272) * t253 + ((-t172 - t174) * t252 + (t168 + t170) * t250 + t469) * t251;
t169 = t193 * t248;
t171 = t195 * t248;
t173 = t197 * t248;
t175 = t199 * t248;
t281 = -t128 * t250 + t132 * t252;
t271 = -t191 * t248 - t281;
t283 = -t126 * t250 + t130 * t252;
t273 = -t189 * t248 - t283;
t438 = (-t271 - t273) * t253 + ((-t173 - t175) * t252 + (t169 + t171) * t250 + t468) * t251;
t437 = t458 * t216 + t456 * t217 + t449 * t247;
t436 = t457 * t216 + t455 * t217 + t450 * t247;
t435 = t458 * t218 + t456 * t219 + t449 * t248;
t434 = t457 * t218 + t455 * t219 + t450 * t248;
t287 = -Icges(3,5) * t251 - Icges(3,6) * t253;
t220 = t287 * t247;
t221 = t287 * t248;
t433 = t454 * t216 + t453 * t217 + t451 * t247;
t432 = t454 * t218 + t453 * t219 + t451 * t248;
t431 = t452 * t251 - t423;
t430 = (t288 + t289) * t253 + (Icges(4,6) + Icges(5,6)) * t251;
t429 = (t292 + t293) * t253 + (Icges(4,5) + Icges(5,5)) * t251;
t244 = t247 ^ 2;
t245 = t248 ^ 2;
t414 = t244 + t245;
t298 = rSges(5,1) * t252 - rSges(5,2) * t250;
t390 = pkin(3) * t252;
t428 = (rSges(5,3) + qJ(4)) * t253 + (-t298 - t390) * t251;
t427 = ((t286 + t285) * t253 + (Icges(4,3) + Icges(5,3)) * t251 - t452) * t253;
t367 = t121 * t253;
t88 = t251 * t284 - t367;
t366 = t122 * t253;
t89 = t251 * t283 - t366;
t365 = t123 * t253;
t90 = t251 * t282 - t365;
t364 = t124 * t253;
t91 = t251 * t281 - t364;
t424 = (t89 + t91) * t248 + (t88 + t90) * t247;
t159 = rSges(4,1) * t216 - rSges(4,2) * t217;
t161 = rSges(4,1) * t218 - rSges(4,2) * t219;
t108 = t159 * t247 + t161 * t248;
t237 = t251 * pkin(2) - pkin(5) * t253;
t304 = -t237 + t428;
t115 = t304 * t247;
t117 = t304 * t248;
t299 = rSges(4,1) * t252 - rSges(4,2) * t250;
t202 = -rSges(4,3) * t253 + t251 * t299;
t319 = -t202 - t237;
t154 = t319 * t247;
t156 = t319 * t248;
t300 = (rSges(5,2) * t252 + t440 * t250) * t251;
t162 = t300 * t247;
t163 = t300 * t248;
t233 = (-rSges(4,1) * t250 - rSges(4,2) * t252) * t251;
t238 = pkin(2) * t253 + t251 * pkin(5);
t316 = t414 * t238;
t182 = qJ(4) * t251 + t390 * t253;
t330 = rSges(5,1) * t219 + rSges(5,2) * t218 + rSges(5,3) * t351 + pkin(3) * t355 + t182 * t248;
t331 = rSges(5,1) * t217 + rSges(5,2) * t216 + rSges(5,3) * t354 - pkin(3) * t352 + t182 * t247;
t65 = t331 * t247 + t330 * t248 + t316;
t136 = rSges(4,1) * t217 + rSges(4,2) * t216 + rSges(4,3) * t354;
t138 = rSges(4,1) * t219 + rSges(4,2) * t218 + rSges(4,3) * t351;
t93 = t136 * t247 + t138 * t248 + t316;
t326 = -rSges(5,2) * t219 + t440 * t218;
t327 = -rSges(5,2) * t217 + t440 * t216;
t97 = t327 * t247 + t326 * t248;
t422 = -m(5) * (t115 * t162 + t117 * t163 + t65 * t97) - m(4) * (t108 * t93 + (-t154 * t247 - t156 * t248) * t233);
t421 = -t247 / 0.2e1;
t395 = t247 / 0.2e1;
t394 = -t248 / 0.2e1;
t393 = t248 / 0.2e1;
t420 = -t251 / 0.2e1;
t419 = t251 / 0.2e1;
t418 = -t253 / 0.2e1;
t416 = (-t444 * t216 + t443 * t217) * t253 + (t447 * t248 + (t442 + t448) * t247) * t251;
t415 = (-t444 * t218 + t443 * t219) * t253 + ((t442 + t445) * t248 + t446 * t247) * t251;
t324 = t251 * rSges(5,3) + t253 * t298 + t182;
t412 = t324 * t251 - t253 * t428;
t409 = 2 * qJD(2);
t408 = 2 * qJD(3);
t407 = 4 * qJD(3);
t406 = m(4) / 0.2e1;
t405 = m(5) / 0.2e1;
t350 = t248 * t253;
t353 = t247 * t253;
t302 = t331 * t253;
t86 = -t354 * t428 + t302;
t301 = t330 * t253;
t87 = t351 * t428 - t301;
t388 = t86 * t350 + t87 * t353;
t328 = t428 * t248;
t329 = t428 * t247;
t41 = (t329 * t251 + t302) * t248 + (-t328 * t251 - t301) * t247;
t51 = t247 * t412 - t331 * t251 + t329 * t253;
t52 = -t248 * t412 + t330 * t251 - t328 * t253;
t71 = (-t330 * t247 + t331 * t248) * t251;
t404 = m(5) * (-t253 * t41 + (t247 * t52 + t248 * t51 + t71) * t251 + t388);
t403 = m(5) * (t41 * t71 + t51 * t86 + t52 * t87);
t363 = t136 * t253;
t111 = t202 * t354 + t363;
t362 = t138 * t253;
t112 = -t202 * t351 - t362;
t177 = t202 * t247;
t179 = t202 * t248;
t83 = (-t177 * t251 + t363) * t248 + (t179 * t251 - t362) * t247;
t204 = t251 * rSges(4,3) + t253 * t299;
t278 = t202 * t253 + t204 * t251;
t94 = -t251 * t136 - t177 * t253 + t247 * t278;
t95 = t251 * t138 + t179 * t253 - t248 * t278;
t99 = (t136 * t248 - t138 * t247) * t251;
t402 = m(4) * (t111 * t94 + t112 * t95 + t83 * t99);
t234 = t414 * t251;
t399 = m(5) * (t234 * t71 + t388);
t398 = m(5) * (-t253 * t97 + (t162 * t247 + t163 * t248) * t251);
t397 = t234 / 0.2e1;
t387 = m(5) * qJD(2);
t386 = m(5) * qJD(4);
t15 = 0.2e1 * (t41 / 0.4e1 - t97 / 0.4e1) * m(5) + 0.2e1 * (t83 / 0.4e1 - t108 / 0.4e1) * m(4);
t361 = t15 * qJD(1);
t340 = t115 * t353 + t117 * t350;
t318 = -t204 - t238;
t317 = t414 * t237;
t315 = t414 * t342;
t314 = qJD(3) * t251;
t184 = (t397 + t420) * m(5);
t313 = t184 * qJD(1);
t257 = t251 * t274 + t367;
t43 = -t216 * t168 - t217 * t172 + t247 * t257;
t256 = t251 * t273 + t366;
t44 = -t216 * t169 - t217 * t173 + t247 * t256;
t255 = t251 * t272 + t365;
t45 = -t216 * t170 - t217 * t174 + t247 * t255;
t254 = t251 * t271 + t364;
t46 = -t216 * t171 - t217 * t175 + t247 * t254;
t312 = t433 * t419 + (t216 * t430 + t217 * t429) * t418 + ((-t423 + t437) * t253 + (t43 + t45 - t427) * t251) * t395 + (t436 * t253 + (t44 + t46) * t251) * t393;
t47 = -t218 * t168 - t219 * t172 + t248 * t257;
t48 = -t218 * t169 - t219 * t173 + t248 * t256;
t49 = -t218 * t170 - t219 * t174 + t248 * t255;
t50 = -t218 * t171 - t219 * t175 + t248 * t254;
t311 = t432 * t419 + (t218 * t430 + t219 * t429) * t418 + (t435 * t253 + (t47 + t49) * t251) * t395 + ((-t423 + t434) * t253 + (t50 + t48 - t427) * t251) * t393;
t310 = ((t250 * t430 - t252 * t429 - t467) * t253 + t438 * t248 + t439 * t247 + t431) * t420 + (t424 + t427) * t418;
t309 = t448 * t393 + t447 * t421;
t308 = t446 * t394 + t445 * t395;
t307 = (t437 * t247 + t436 * t248) * t419 + t433 * t418;
t306 = (t435 * t247 + t434 * t248) * t419 + t432 * t418;
t305 = t424 * t420 + t431 * t253 / 0.2e1;
t303 = -t238 - t324;
t235 = t251 * rSges(3,1) + rSges(3,2) * t253;
t277 = -t309 - t312;
t276 = -t308 + t311;
t275 = t300 * t251;
t267 = t441 * t247 + t221;
t266 = -t441 * t248 + t220;
t183 = m(5) * t397 + t251 * t405;
t180 = t315 - t342;
t157 = t318 * t248;
t155 = t318 * t247;
t139 = t414 * t235;
t120 = -t161 * t253 - t233 * t351;
t119 = t159 * t253 + t233 * t354;
t118 = t303 * t248;
t116 = t303 * t247;
t106 = (t159 * t248 - t161 * t247) * t251;
t105 = t275 * t248 - t326 * t253;
t104 = -t275 * t247 + t327 * t253;
t98 = -t177 * t247 - t179 * t248 - t317;
t92 = (-t326 * t247 + t327 * t248) * t251;
t73 = t329 * t247 + t328 * t248 - t317;
t70 = t398 / 0.2e1;
t69 = -t143 * t253 + (-t332 * t250 + t336 * t252) * t251;
t68 = -t142 * t253 + (-t333 * t250 + t337 * t252) * t251;
t67 = -t141 * t253 + (-t334 * t250 + t338 * t252) * t251;
t66 = -t140 * t253 + (-t335 * t250 + t339 * t252) * t251;
t39 = t234 * t65 + t340;
t31 = t399 / 0.2e1;
t25 = t247 * t50 - t248 * t49;
t24 = t247 * t48 - t248 * t47;
t23 = t247 * t46 - t248 * t45;
t22 = t247 * t44 - t248 * t43;
t16 = (t108 + t83) * t406 + (t41 + t97) * t405;
t6 = t404 / 0.2e1;
t5 = t31 + t70 - t404 / 0.2e1;
t4 = t31 + t6 - t398 / 0.2e1;
t3 = t70 + t6 - t399 / 0.2e1;
t2 = t308 * t247 + t309 * t248 - t422;
t1 = t402 + t403 + (t307 * t247 + t306 * t248 + t310) * t253 + (t312 * t247 + t311 * t248 - t305) * t251;
t7 = [0, t16 * qJD(3) + t183 * qJD(4) + (-m(3) * t139 / 0.2e1 + t98 * t406 + t73 * t405) * t409, t16 * qJD(2) + (t106 * t406 + t92 * t405) * t408, t183 * qJD(2); -qJD(3) * t15 + qJD(4) * t184, t2 * qJD(3) + t39 * t386 + (m(5) * (t115 * t116 + t117 * t118 + t65 * t73) + m(4) * (t154 * t155 + t156 * t157 + t93 * t98) + m(3) * (-t139 + t235) * t414 * (rSges(3,1) * t253 - t251 * rSges(3,2)) + (t244 * t221 + (t267 * t248 + (-t220 + t266) * t247) * t248 + t25 + t24) * t395 + (t245 * t220 + (t266 * t247 + (-t221 + t267) * t248) * t247 + t23 + t22) * t394) * qJD(2), -t361 + t2 * qJD(2) + t5 * qJD(4) + (-t403 / 0.4e1 - t402 / 0.4e1) * t407 + ((t104 * t117 + t105 * t115 + t162 * t87 + t163 * t86 + t65 * t92 + t71 * t97) * t405 + (t106 * t93 + t108 * t99 + t119 * t156 + t120 * t154 + (-t111 * t248 - t112 * t247) * t233) * t406) * t408 + (t277 * t247 - t276 * t248 + t305) * t314 + (((t66 / 0.2e1 + t68 / 0.2e1 - t306) * t248 + (-t67 / 0.2e1 - t69 / 0.2e1 - t307) * t247 - t310) * t253 + t415 * t395 + t416 * t394) * qJD(3), t313 + t39 * t387 + t5 * qJD(3) + (-t234 * t253 - t180 + t315) * t386; t15 * qJD(2), t361 + t1 * qJD(3) + t4 * qJD(4) + ((t115 * t52 + t116 * t87 + t117 * t51 + t118 * t86 + t41 * t65 + t71 * t73) * t405 + (t111 * t157 + t112 * t155 + t154 * t95 + t156 * t94 + t83 * t93 + t98 * t99) * t406) * t409 + (t277 * t248 + t276 * t247 + (t438 * t421 + (t436 * t247 - t437 * t248) * t395 + (t434 * t247 - t435 * t248 + t439) * t393) * t253 + ((-t90 / 0.2e1 - t88 / 0.2e1 + t25 / 0.2e1 + t24 / 0.2e1) * t248 + (t91 / 0.2e1 + t89 / 0.2e1 + t22 / 0.2e1 + t23 / 0.2e1) * t247) * t251 + t422) * qJD(2), t1 * qJD(2) + (m(5) * (t104 * t86 + t105 * t87 + t71 * t92) / 0.4e1 + m(4) * (t106 * t99 + t111 * t119 + t112 * t120) / 0.4e1) * t407 + (-t227 / 0.2e1 - t226 / 0.2e1) * qJD(3) * t253 * t410 + (t416 * t395 + t415 * t393 + ((t444 * t250 + t443 * t252) * t253 + (t69 + t67) * t248 + (t68 + t66) * t247) * t418) * t314, t4 * qJD(2); -t184 * qJD(2), -t313 + (-t253 * t73 + (t116 * t247 + t118 * t248 + t65) * t251 - t39 + t340) * t387 + t3 * qJD(3) + t180 * t386, t3 * qJD(2) + m(5) * (-t253 * t92 + (t104 * t248 + t105 * t247) * t251) * qJD(3), t180 * t387;];
Cq = t7;
