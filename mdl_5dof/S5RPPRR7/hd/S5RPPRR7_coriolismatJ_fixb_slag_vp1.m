% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR7_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:32
% EndTime: 2019-12-31 17:59:45
% DurationCPUTime: 9.74s
% Computational Cost: add. (30139->502), mult. (36917->751), div. (0->0), fcn. (40335->8), ass. (0->294)
t296 = qJ(1) + pkin(8);
t294 = sin(t296);
t295 = cos(t296);
t298 = sin(qJ(4));
t427 = pkin(4) * t298;
t428 = cos(qJ(1)) * pkin(1);
t453 = pkin(2) + pkin(6);
t300 = cos(qJ(5));
t297 = sin(qJ(5));
t392 = t297 * t298;
t239 = -t294 * t392 + t295 * t300;
t388 = t298 * t300;
t240 = t294 * t388 + t295 * t297;
t301 = cos(qJ(4));
t397 = t294 * t301;
t184 = t240 * rSges(6,1) + t239 * rSges(6,2) - rSges(6,3) * t397;
t470 = -pkin(7) * t397 + t184;
t150 = t428 + t453 * t295 + (qJ(3) + t427) * t294 + t470;
t202 = rSges(6,1) * t239 - rSges(6,2) * t240;
t241 = t294 * t300 + t295 * t392;
t242 = -t294 * t297 + t295 * t388;
t203 = rSges(6,1) * t241 + rSges(6,2) * t242;
t395 = t295 * t301;
t284 = pkin(7) * t395;
t344 = -sin(qJ(1)) * pkin(1) + t295 * qJ(3);
t309 = -t453 * t294 + t344;
t417 = rSges(6,3) * t301;
t469 = t242 * rSges(6,1) - t241 * rSges(6,2);
t480 = (-t417 + t427) * t295 - t284 + t309 + t469;
t489 = m(6) * (t150 * t202 - t203 * t480);
t458 = m(6) / 0.2e1;
t421 = m(6) * qJD(5);
t175 = Icges(6,5) * t240 + Icges(6,6) * t239 - Icges(6,3) * t397;
t410 = Icges(6,4) * t240;
t178 = Icges(6,2) * t239 - Icges(6,6) * t397 + t410;
t234 = Icges(6,4) * t239;
t181 = Icges(6,1) * t240 - Icges(6,5) * t397 + t234;
t94 = t175 * t395 + t241 * t178 - t242 * t181;
t488 = t294 * t94;
t487 = t295 * t94;
t339 = rSges(6,1) * t300 - rSges(6,2) * t297;
t258 = rSges(6,3) * t298 + t339 * t301;
t186 = rSges(6,3) * t395 - t469;
t400 = t186 * t298;
t486 = t258 * t395 - t400;
t177 = -Icges(6,5) * t242 + Icges(6,6) * t241 + Icges(6,3) * t395;
t402 = t177 * t298;
t236 = Icges(6,4) * t242;
t180 = Icges(6,2) * t241 + Icges(6,6) * t395 - t236;
t235 = Icges(6,4) * t241;
t182 = Icges(6,1) * t242 - Icges(6,5) * t395 - t235;
t484 = t180 * t297 + t182 * t300;
t109 = t484 * t301 - t402;
t326 = -t180 * t241 - t182 * t242;
t377 = t239 * t178 + t240 * t181;
t485 = t326 + t377 + (-t175 * t294 - t177 * t295) * t301;
t93 = -t177 * t397 + t239 * t180 - t182 * t240;
t459 = m(5) / 0.2e1;
t408 = Icges(6,4) * t300;
t333 = -Icges(6,2) * t297 + t408;
t246 = Icges(6,6) * t298 + t333 * t301;
t409 = Icges(6,4) * t297;
t335 = Icges(6,1) * t300 - t409;
t248 = Icges(6,5) * t298 + t335 * t301;
t330 = Icges(6,5) * t300 - Icges(6,6) * t297;
t244 = Icges(6,3) * t298 + t330 * t301;
t383 = t301 * t244;
t134 = t241 * t246 - t242 * t248 + t295 * t383;
t482 = t298 * t134;
t412 = Icges(5,4) * t298;
t334 = Icges(5,2) * t301 + t412;
t229 = Icges(5,6) * t295 + t334 * t294;
t276 = Icges(5,4) * t397;
t398 = t294 * t298;
t231 = Icges(5,1) * t398 + Icges(5,5) * t295 + t276;
t478 = t301 * t229 + t298 * t231;
t481 = t478 * t295;
t230 = -Icges(5,6) * t294 + t334 * t295;
t411 = Icges(5,4) * t301;
t336 = Icges(5,1) * t298 + t411;
t232 = -Icges(5,5) * t294 + t336 * t295;
t479 = t301 * t230 + t298 * t232;
t443 = -t294 / 0.2e1;
t442 = t294 / 0.2e1;
t440 = t295 / 0.2e1;
t475 = t294 * t93;
t474 = t295 * t93;
t472 = t479 * t295;
t189 = ((rSges(6,3) + pkin(7)) * t298 + (pkin(4) + t339) * t301) * t295;
t274 = rSges(5,1) * t301 - rSges(5,2) * t298;
t255 = t274 * t294;
t256 = t274 * t295;
t219 = rSges(6,3) * t398 + t339 * t397;
t351 = -pkin(4) * t397 - pkin(7) * t398 - t219;
t380 = (t189 * t294 + t295 * t351) * t458 + (-t255 * t295 + t256 * t294) * t459;
t243 = Icges(6,3) * t301 - t330 * t298;
t386 = t300 * t248;
t393 = t297 * t246;
t322 = t386 - t393;
t313 = t243 - t322;
t399 = t244 * t298;
t467 = t313 * t301 - t399;
t314 = -t244 * t295 + t484;
t466 = t314 * t301 - t402;
t327 = -t178 * t297 + t181 * t300;
t315 = t244 * t294 - t327;
t403 = t175 * t298;
t465 = t315 * t301 - t403;
t263 = (-Icges(6,2) * t300 - t409) * t301;
t264 = (-Icges(6,1) * t297 - t408) * t301;
t464 = -(t248 / 0.2e1 + t263 / 0.2e1) * t297 + (t264 / 0.2e1 - t246 / 0.2e1) * t300;
t269 = -Icges(5,2) * t298 + t411;
t252 = t269 * t295;
t271 = Icges(5,1) * t301 - t412;
t254 = t271 * t295;
t463 = t298 * (t230 - t254) - t301 * (t232 + t252);
t251 = -Icges(5,2) * t398 + t276;
t253 = t271 * t294;
t462 = t298 * (t229 - t253) - (t231 + t251) * t301;
t292 = t294 ^ 2;
t293 = t295 ^ 2;
t461 = 4 * qJD(1);
t460 = 2 * qJD(4);
t132 = t239 * t246 + t240 * t248 - t294 * t383;
t130 = t132 * t298;
t92 = -t175 * t397 + t377;
t338 = -t294 * t92 + t474;
t39 = t338 * t301 + t130;
t457 = -t39 / 0.2e1;
t215 = t246 * t294;
t217 = t248 * t294;
t80 = (-t215 * t297 + t217 * t300 + t175) * t301 + t315 * t298;
t456 = t80 / 0.2e1;
t216 = t246 * t295;
t218 = t248 * t295;
t81 = (t216 * t297 - t218 * t300 + t177) * t301 + t314 * t298;
t455 = t81 / 0.2e1;
t197 = Icges(6,5) * t241 + Icges(6,6) * t242;
t373 = Icges(6,2) * t242 - t182 + t235;
t375 = -Icges(6,1) * t241 + t180 - t236;
t86 = t197 * t298 + (-t373 * t297 - t375 * t300) * t301;
t454 = t86 / 0.2e1;
t257 = -t339 * t298 + t417;
t321 = t257 * t301 - t258 * t298;
t118 = t184 * t301 + t219 * t298 + t321 * t294;
t220 = t258 * t295;
t119 = -t186 * t301 + t220 * t298 + t321 * t295;
t123 = (-t184 * t295 - t186 * t294) * t301;
t401 = t184 * t298;
t155 = t258 * t397 + t401;
t96 = (-t219 * t301 + t401) * t295 + (t220 * t301 + t400) * t294;
t451 = m(6) * (t118 * t155 + t119 * t486 + t123 * t96);
t450 = m(6) * (t118 * t150 + t119 * t480 - t155 * t351 + t189 * t486);
t280 = pkin(4) * t301 + pkin(7) * t298;
t363 = t258 + t280;
t205 = t363 * t294;
t207 = t363 * t295;
t265 = (-rSges(6,1) * t297 - rSges(6,2) * t300) * t301;
t448 = m(6) * (-t202 * t207 - t203 * t205 + (-t150 * t295 + t294 * t480) * t265);
t446 = m(6) * (-t150 * t351 + t189 * t480);
t445 = m(6) * (t150 * t294 + t295 * t480);
t444 = -t109 / 0.2e1;
t441 = t294 / 0.4e1;
t439 = t295 / 0.4e1;
t438 = t298 / 0.2e1;
t437 = m(4) * ((rSges(4,3) * t295 + t344) * t295 + (t428 + (rSges(4,3) + qJ(3)) * t294) * t294);
t341 = rSges(5,1) * t298 + rSges(5,2) * t301;
t308 = -t294 * rSges(5,3) + t341 * t295;
t192 = t308 + t309;
t193 = t428 + (rSges(5,3) + t453) * t295 + (qJ(3) + t341) * t294;
t436 = m(5) * (t192 * t256 + t193 * t255);
t435 = m(5) * (t192 * t295 + t193 * t294);
t433 = m(6) * (t155 * t294 + t295 * t486);
t144 = -t202 * t295 - t203 * t294;
t431 = m(6) * t144;
t430 = m(6) * (t205 * t295 - t207 * t294);
t95 = t177 * t395 - t326;
t426 = t95 + t485;
t423 = -t92 + t485;
t416 = t294 * t39;
t337 = t295 * t95 - t488;
t40 = t337 * t301 + t482;
t415 = t295 * t40;
t396 = t295 * t298;
t245 = Icges(6,6) * t301 - t333 * t298;
t394 = t297 * t245;
t262 = (-Icges(6,5) * t297 - Icges(6,6) * t300) * t301;
t389 = t298 * t262;
t247 = Icges(6,5) * t301 - t335 * t298;
t387 = t300 * t247;
t376 = -Icges(6,1) * t239 + t178 + t410;
t374 = -Icges(6,2) * t240 + t181 + t234;
t366 = t246 - t264;
t365 = t248 + t263;
t364 = pkin(7) * t301 + t257 - t427;
t362 = t292 + t293;
t361 = qJD(1) * t298;
t360 = qJD(1) * t301;
t359 = qJD(5) * t301;
t196 = Icges(6,5) * t239 - Icges(6,6) * t240;
t70 = -t196 * t397 + t374 * t239 - t376 * t240;
t71 = -t197 * t397 + t373 * t239 - t375 * t240;
t28 = t294 * t71 + t295 * t70;
t100 = t239 * t245 + t240 * t247 - t294 * t467;
t66 = t215 * t239 + t217 * t240 - t294 * t465;
t67 = -t216 * t239 - t218 * t240 - t294 * t466;
t8 = (-t294 * t66 + t295 * t67 + t132) * t301 + (t100 - t338) * t298;
t357 = -t28 / 0.2e1 + t8 / 0.2e1;
t72 = t196 * t395 + t374 * t241 + t376 * t242;
t73 = t197 * t395 + t373 * t241 + t375 * t242;
t29 = t294 * t73 + t295 * t72;
t101 = t241 * t245 - t242 * t247 + t295 * t467;
t68 = t215 * t241 - t217 * t242 + t295 * t465;
t69 = -t216 * t241 + t218 * t242 + t295 * t466;
t9 = (-t294 * t68 + t295 * t69 + t134) * t301 + (t101 - t337) * t298;
t356 = t29 / 0.2e1 - t9 / 0.2e1;
t13 = -t482 + (t423 * t295 + t488) * t301;
t353 = -t13 / 0.2e1 - t40 / 0.2e1;
t12 = t130 + (-t426 * t294 + t474) * t301;
t352 = t457 + t12 / 0.2e1;
t331 = Icges(5,5) * t298 + Icges(5,6) * t301;
t138 = t295 * (Icges(5,3) * t295 + t331 * t294) + t478 * t294;
t228 = -Icges(5,3) * t294 + t331 * t295;
t139 = -t295 * t228 - t479 * t294;
t349 = -t397 / 0.4e1;
t347 = t395 / 0.4e1;
t343 = t362 * t341;
t111 = t365 * t239 - t366 * t240 - t262 * t397;
t112 = t365 * t241 + t366 * t242 + t262 * t395;
t85 = t196 * t298 + (-t374 * t297 - t376 * t300) * t301;
t342 = t448 / 0.2e1 + (t112 + t86) * t441 + (t111 + t85) * t439;
t332 = Icges(5,5) * t301 - Icges(5,6) * t298;
t311 = m(6) * (-t118 * t295 + t119 * t294);
t320 = t362 * t265 * t458;
t65 = -t311 / 0.2e1 + t320;
t145 = -t202 * t294 + t203 * t295;
t74 = 0.2e1 * (t96 / 0.4e1 - t145 / 0.4e1) * m(6);
t329 = t74 * qJD(2) - t65 * qJD(3);
t108 = t327 * t301 + t403;
t328 = -t108 * t294 - t109 * t295;
t141 = -t294 * t228 + t472;
t17 = t423 * t294 - t487;
t46 = t294 * t95 + t487;
t319 = t292 * t228 / 0.2e1 + t141 * t442 + t17 / 0.2e1 + t46 / 0.2e1 + (t139 + (t228 + t478) * t295 - t481) * t440;
t16 = t426 * t295 + t475;
t45 = t295 * t92 + t475;
t318 = -t138 * t295 / 0.2e1 + t139 * t443 + (t139 + t481) * t442 + (t141 - t472 + (t228 - t478) * t294 + t138) * t440 - t45 / 0.2e1 + t16 / 0.2e1;
t312 = t109 / 0.2e1 + t444;
t307 = t12 * t441 + t13 * t439 + t16 * t347 - t416 / 0.4e1 + t415 / 0.4e1 - t45 * t395 / 0.4e1 + (t17 + t46) * t349;
t117 = (t244 + t387 - t394) * t301 + t313 * t298;
t160 = t322 * t301 + t399;
t306 = t117 * t438 + t160 * t301 / 0.2e1 + t450 / 0.2e1 + (t108 + t132) * t398 / 0.4e1 + (t80 + t100) * t349 - (-t109 + t134) * t396 / 0.4e1 + (t81 + t101) * t347;
t305 = t387 / 0.2e1 - t394 / 0.2e1 + t244 / 0.2e1 - t336 / 0.2e1 - t269 / 0.2e1;
t304 = t386 / 0.2e1 - t393 / 0.2e1 + t271 / 0.2e1 - t334 / 0.2e1 - t243 / 0.2e1;
t250 = t332 * t295;
t249 = t294 * t332;
t206 = t364 * t295;
t204 = t364 * t294;
t191 = -t255 * t294 - t256 * t295;
t168 = -t203 * t298 + t265 * t395;
t167 = t202 * t298 + t265 * t397;
t146 = t430 / 0.2e1;
t136 = t431 / 0.2e1;
t135 = t144 * t301;
t129 = (t389 + (-t365 * t297 - t366 * t300) * t301) * t298;
t122 = t351 * t294 + (-t280 * t295 - t220) * t295;
t115 = (-pkin(4) * t396 + t186 + t284) * t295 + (-pkin(4) * t398 - t470) * t294;
t104 = t433 / 0.2e1;
t75 = (t145 + t96) * t458;
t64 = t311 / 0.2e1 + t320;
t60 = t146 - t380;
t59 = -t430 / 0.2e1 + t380;
t58 = t146 + t380;
t55 = t115 * t145 + (t205 * t294 + t207 * t295) * t265;
t54 = t389 / 0.2e1 + t489 + t464 * t301;
t52 = t435 + t437 + t445;
t42 = t160 * t298 + t328 * t301;
t33 = -t304 * t298 + t305 * t301 + t436 + t446;
t32 = t104 - t431 / 0.2e1;
t31 = t136 + t104;
t30 = t136 - t433 / 0.2e1;
t27 = t294 * t69 + t295 * t68;
t26 = t294 * t67 + t295 * t66;
t22 = t112 * t298 + (-t294 * t72 + t295 * t73) * t301;
t21 = t111 * t298 + (-t294 * t70 + t295 * t71) * t301;
t18 = (-t80 * t294 + t81 * t295 + t160) * t301 + (t117 - t328) * t298;
t7 = m(6) * t55 + t28 * t440 + t29 * t442;
t6 = (t353 * t294 + t352 * t295) * t301;
t5 = t318 * t294 + t319 * t295;
t4 = t451 + (t8 * t443 + t9 * t440 + t42 / 0.2e1) * t301 + (t416 / 0.2e1 - t415 / 0.2e1 + t18 / 0.2e1) * t298;
t3 = t307 - t450 / 0.2e1 + (-t160 / 0.2e1 + (-t101 / 0.4e1 - t81 / 0.4e1) * t295 + (t100 / 0.4e1 + t80 / 0.4e1) * t294) * t301 + (-t117 / 0.2e1 + (t134 / 0.4e1 - t109 / 0.4e1) * t295 + (-t132 / 0.4e1 - t108 / 0.4e1) * t294) * t298 + t342;
t2 = t306 + (-t112 / 0.4e1 - t86 / 0.4e1) * t294 + (-t111 / 0.4e1 - t85 / 0.4e1) * t295 + t307 - t448 / 0.2e1;
t1 = (t39 / 0.4e1 - t12 / 0.4e1 + (t46 / 0.4e1 + t17 / 0.4e1) * t301) * t294 + t306 + (-t40 / 0.4e1 - t13 / 0.4e1 + (t45 / 0.4e1 - t16 / 0.4e1) * t301) * t295 + t342;
t10 = [t52 * qJD(3) + t33 * qJD(4) + t54 * qJD(5), 0, qJD(1) * t52 + qJD(4) * t58 + qJD(5) * t31, t33 * qJD(1) + t58 * qJD(3) + t1 * qJD(5) + (m(6) * (-t150 * t206 + t189 * t205 + t204 * t480 + t207 * t351) + (t456 + t100 / 0.2e1 + m(5) * (t193 * t341 - t255 * t274) - t331 * t440 + (-t229 / 0.2e1 + t253 / 0.2e1) * t301 + (-t231 / 0.2e1 - t251 / 0.2e1) * t298 - t319) * t295 + (t101 / 0.2e1 + t455 + m(5) * (-t192 * t341 + t256 * t274) - t331 * t442 + (t230 / 0.2e1 - t254 / 0.2e1) * t301 + (t232 / 0.2e1 + t252 / 0.2e1) * t298 - t318) * t294) * qJD(4), t54 * qJD(1) + t31 * qJD(3) + t1 * qJD(4) + t129 * qJD(5) + (t150 * t167 + t155 * t202 + t168 * t480 - t203 * t486) * t421 + ((t112 / 0.2e1 + t454 - t352) * t295 + (-t85 / 0.2e1 - t111 / 0.2e1 - t353) * t294) * t359; 0, 0, 0, (t122 * t458 + t191 * t459) * t460 + t75 * qJD(5), t75 * qJD(4) + t135 * t421; t59 * qJD(4) + t30 * qJD(5) + (-t445 / 0.4e1 - t435 / 0.4e1 - t437 / 0.4e1) * t461, 0, 0, t59 * qJD(1) + ((t204 * t294 + t206 * t295) * t458 - t343 * t459) * t460 + t64 * qJD(5), t30 * qJD(1) + t64 * qJD(4) + (-t167 * t295 + t168 * t294) * t421; t60 * qJD(3) + t5 * qJD(4) + t3 * qJD(5) + (-t446 / 0.4e1 - t436 / 0.4e1) * t461 + t312 * qJD(1) * t295 - t305 * t360 + t304 * t361, -qJD(5) * t74, qJD(1) * t60 + qJD(5) * t65, t5 * qJD(1) + (m(5) * ((-t295 * t308 + (-t295 * rSges(5,3) - t341 * t294) * t294) * t191 - t274 * t343) + m(6) * (t115 * t122 + t204 * t205 + t206 * t207) + (-t292 * t250 + (t462 * t295 + (t249 - t463) * t294) * t295 + t27) * t442 + (t293 * t249 + (t463 * t294 + (-t250 - t462) * t295) * t294 + t26) * t440) * qJD(4) + t7 * qJD(5), t3 * qJD(1) + t7 * qJD(4) + (-t42 / 0.2e1 + t356 * t295 + t357 * t294) * t359 - t329 + (m(6) * (t115 * t135 + t123 * t145 - t167 * t207 + t168 * t205 + (-t155 * t295 + t294 * t486) * t265) + t21 * t440 + t22 * t442 - t451 + (-t18 / 0.2e1 + (t85 / 0.2e1 + t40 / 0.2e1) * t295 + (t454 + t457) * t294) * t298) * qJD(5); -t262 * t361 / 0.2e1 + t32 * qJD(3) + t2 * qJD(4) + t6 * qJD(5) - qJD(1) * t489 + (-t294 * t312 - t464) * t360, qJD(4) * t74, qJD(1) * t32 - qJD(4) * t65, t2 * qJD(1) + (((t27 / 0.2e1 + t108 / 0.2e1) * t301 + (-t46 / 0.2e1 + t456) * t298 + t357) * t295 + ((-t26 / 0.2e1 + t444) * t301 + (t45 / 0.2e1 + t455) * t298 - t356) * t294 + (t115 * t96 - t118 * t207 + t119 * t205 + t122 * t123 - t155 * t206 + t204 * t486 - t55) * m(6)) * qJD(4) + t4 * qJD(5) + t329, t6 * qJD(1) + t4 * qJD(4) + (m(6) * (t123 * t135 + t155 * t167 + t168 * t486) + t129 * t438 + (t21 * t443 + t22 * t440 + (-t294 * t85 + t295 * t86) * t438) * t301) * qJD(5);];
Cq = t10;
