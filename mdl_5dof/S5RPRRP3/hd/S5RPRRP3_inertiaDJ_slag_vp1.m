% Calculate time derivative of joint inertia matrix for
% S5RPRRP3
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:10
% EndTime: 2019-12-05 18:03:32
% DurationCPUTime: 8.70s
% Computational Cost: add. (9208->383), mult. (7890->513), div. (0->0), fcn. (6082->8), ass. (0->232)
t419 = Icges(5,4) + Icges(6,4);
t418 = Icges(5,1) + Icges(6,1);
t417 = Icges(5,2) + Icges(6,2);
t202 = qJ(3) + qJ(4);
t197 = cos(t202);
t416 = t419 * t197;
t196 = sin(t202);
t415 = t419 * t196;
t414 = Icges(5,5) + Icges(6,5);
t413 = Icges(5,6) + Icges(6,6);
t412 = -t417 * t196 + t416;
t411 = t418 * t197 - t415;
t201 = qJ(1) + pkin(8);
t194 = sin(t201);
t195 = cos(t201);
t382 = t413 * t194 + t412 * t195;
t380 = t414 * t194 + t411 * t195;
t410 = Icges(5,3) + Icges(6,3);
t383 = -t412 * t194 + t413 * t195;
t381 = -t411 * t194 + t414 * t195;
t409 = -t413 * t196 + t414 * t197;
t402 = t417 * t197 + t415;
t401 = t418 * t196 + t416;
t408 = t382 * t196 - t380 * t197;
t406 = -t409 * t194 + t410 * t195;
t405 = t410 * t194 + t409 * t195;
t407 = t383 * t196 - t381 * t197;
t200 = qJD(3) + qJD(4);
t404 = t412 * t200;
t403 = t411 * t200;
t379 = t414 * t196 + t413 * t197;
t400 = t408 * t194;
t377 = t402 * t196 - t401 * t197;
t398 = -t407 * t194 - t406 * t195;
t397 = -t405 * t195 - t400;
t303 = t195 * t200;
t396 = t383 * qJD(1) - t402 * t303;
t307 = t194 * t200;
t395 = t382 * qJD(1) - t402 * t307;
t394 = -t381 * qJD(1) + t401 * t303;
t393 = -t380 * qJD(1) + t401 * t307;
t392 = t406 * qJD(1);
t391 = t405 * qJD(1);
t390 = t407 * t195;
t389 = t406 * t194 - t390;
t388 = t405 * t194 - t195 * t408;
t387 = -t379 * t303 + t392;
t386 = t379 * t307 - t391;
t207 = -pkin(7) - pkin(6);
t184 = t194 * t207;
t199 = -qJ(5) + t207;
t337 = rSges(6,1) * t197;
t242 = -rSges(6,2) * t196 + t337;
t205 = cos(qJ(3));
t190 = t205 * pkin(3) + pkin(2);
t163 = pkin(4) * t197 + t190;
t299 = t163 - t190;
t334 = rSges(6,3) * t194;
t385 = -t194 * t199 + t184 + t334 + (t242 + t299) * t195;
t219 = t163 + t242;
t384 = t219 * t194;
t378 = t377 * qJD(1) + t409 * t200;
t376 = t380 * t200 + t396;
t375 = t381 * t200 - t395;
t374 = -t382 * t200 - t394;
t373 = -t383 * t200 + t393;
t203 = sin(qJ(3));
t282 = qJD(3) * t205;
t286 = qJD(1) * t195;
t372 = t194 * t282 + t203 * t286;
t301 = t197 * t200;
t371 = t194 * t301 + t196 * t286;
t370 = (t402 * t200 - t403) * t197 + (t401 * t200 + t404) * t196 - t379 * qJD(1);
t153 = rSges(5,1) * t196 + rSges(5,2) * t197;
t188 = t195 * rSges(5,3);
t338 = rSges(5,1) * t197;
t243 = -rSges(5,2) * t196 + t338;
t369 = -(t243 * t194 - t188) * qJD(1) - t153 * t303;
t368 = t397 * t194 + t398 * t195;
t179 = rSges(4,1) * t203 + rSges(4,2) * t205;
t284 = qJD(3) * t195;
t366 = t179 * t284;
t328 = Icges(4,4) * t205;
t235 = -Icges(4,2) * t203 + t328;
t122 = Icges(4,6) * t195 - t235 * t194;
t329 = Icges(4,4) * t203;
t239 = Icges(4,1) * t205 - t329;
t124 = Icges(4,5) * t195 - t239 * t194;
t223 = t122 * t203 - t124 * t205;
t365 = t223 * t195;
t306 = t194 * t203;
t294 = rSges(4,2) * t306 + t195 * rSges(4,3);
t283 = qJD(3) * t203;
t279 = pkin(3) * t283;
t302 = t196 * t200;
t140 = -pkin(4) * t302 - t279;
t277 = t194 * t302;
t287 = qJD(1) * t194;
t362 = rSges(6,1) * t277 + t371 * rSges(6,2) + qJD(5) * t195 - t140 * t194 + t199 * t287;
t344 = sin(qJ(1)) * pkin(1);
t191 = qJD(1) * t344;
t231 = Icges(4,5) * t205 - Icges(4,6) * t203;
t121 = Icges(4,3) * t194 + t231 * t195;
t123 = Icges(4,6) * t194 + t235 * t195;
t125 = Icges(4,5) * t194 + t239 * t195;
t343 = pkin(2) - t190;
t358 = pkin(6) * t194 + t343 * t195;
t355 = 2 * m(4);
t354 = 2 * m(5);
t353 = 2 * m(6);
t192 = t194 ^ 2;
t193 = t195 ^ 2;
t352 = t194 / 0.2e1;
t351 = t195 / 0.2e1;
t350 = -rSges(4,3) - pkin(6);
t349 = m(4) * t179;
t348 = m(5) * t153;
t347 = cos(qJ(1)) * pkin(1);
t346 = pkin(3) * t203;
t152 = rSges(6,1) * t196 + rSges(6,2) * t197;
t187 = t195 * rSges(6,3);
t285 = qJD(1) * t207;
t273 = t190 * t287 + (t279 + t285) * t195;
t281 = qJD(5) * t194;
t304 = t195 * t199;
t342 = (t140 * t195 - t152 * t303 + t273 + t281 + (t187 - t304 - t384) * qJD(1)) * t195;
t270 = t194 * t283;
t296 = pkin(3) * t270 + t194 * t285;
t305 = t195 * t197;
t341 = t299 * t286 + t296 - (-rSges(6,1) * t305 - t334) * qJD(1) - t362;
t340 = t385 * t195;
t339 = rSges(4,1) * t205;
t336 = rSges(4,3) * t194;
t335 = rSges(5,3) * t194;
t333 = -rSges(6,3) + t199;
t309 = t194 * t196;
t298 = rSges(6,2) * t309 + t187;
t308 = t194 * t197;
t332 = rSges(6,1) * t308 - t298 - (-t199 + t207) * t195 + t299 * t194;
t297 = rSges(5,2) * t309 + t188;
t115 = -rSges(5,1) * t308 + t297;
t98 = (-pkin(6) - t207) * t195 + t343 * t194;
t331 = -t115 - t98;
t314 = t122 * t205;
t313 = t123 * t205;
t312 = t124 * t203;
t311 = t125 * t203;
t136 = t243 * t200;
t310 = t136 * t194;
t166 = pkin(4) * t309;
t300 = qJD(1) * t166 + t152 * t287;
t100 = t194 * t152 + t166;
t295 = t372 * pkin(3);
t293 = t192 + t193;
t120 = Icges(4,3) * t195 - t231 * t194;
t290 = qJD(1) * t120;
t183 = pkin(3) * t306;
t280 = -t98 + t332;
t278 = pkin(3) * t282;
t274 = rSges(5,1) * t277 + t371 * rSges(5,2);
t272 = -rSges(4,1) * t270 - t372 * rSges(4,2);
t267 = -pkin(2) - t339;
t264 = -pkin(4) * t196 - t152;
t263 = -t190 - t338;
t262 = -t163 - t337;
t135 = t242 * t200;
t54 = t371 * pkin(4) + t194 * t135 + t152 * t286;
t249 = -pkin(4) * t301 - t135;
t248 = ((t386 * t195 + (t390 - t397) * qJD(1)) * t195 + t389 * t286) * t195 + ((t387 * t194 + (-t389 + t400) * qJD(1)) * t194 + t388 * t286 + ((t386 - t391) * t194 + (t387 + t392) * t195 + (t380 * t194 - t381 * t195) * t302 + (t382 * t194 - t383 * t195) * t301 + ((t374 + t394) * t194 + (-t373 + t393) * t195) * t197 + ((-t376 + t396) * t194 + (t375 + t395) * t195) * t196 + ((t405 + t407) * t194 + (t408 - t406) * t195 + t388 + t398) * qJD(1)) * t195) * t194;
t246 = -t335 - t347;
t244 = -rSges(4,2) * t203 + t339;
t238 = Icges(4,1) * t203 + t328;
t234 = Icges(4,2) * t205 + t329;
t230 = Icges(4,5) * t203 + Icges(4,6) * t205;
t222 = t123 * t203 - t125 * t205;
t218 = t350 * t194 - t347;
t215 = t222 * t194;
t214 = qJD(3) * t238;
t213 = qJD(3) * t234;
t209 = t368 * t287 + t248;
t208 = (t378 * t194 - t370 * t195 + t374 * t196 + t376 * t197) * t352 + (t370 * t194 + t378 * t195 + t373 * t196 + t375 * t197) * t351 - (t377 * t194 + t379 * t195 + t381 * t196 + t383 * t197) * t287 / 0.2e1 + (t379 * t194 - t377 * t195 + t380 * t196 + t382 * t197) * t286 / 0.2e1;
t185 = pkin(2) * t287;
t174 = qJD(1) * t183;
t162 = t244 * qJD(3);
t127 = t244 * t195 + t336;
t126 = -t194 * t339 + t294;
t119 = (-t153 - t346) * t195;
t118 = t153 * t194 + t183;
t117 = t243 * t195 + t335;
t101 = t264 * t195;
t99 = -t184 - t358;
t97 = t195 * t117;
t95 = t195 * t99;
t94 = (t264 - t346) * t195;
t93 = t183 + t100;
t92 = (-pkin(2) - t244) * t195 + t218;
t91 = t195 * pkin(6) + t267 * t194 + t294 - t344;
t88 = qJD(1) * t358 + t296;
t83 = t230 * t194 * qJD(3) - qJD(1) * t121;
t82 = -t230 * t284 + t290;
t80 = t184 + (-t190 - t243) * t195 + t246;
t79 = t263 * t194 - t195 * t207 + t297 - t344;
t78 = t195 * (-pkin(6) * t286 + t185 - t273);
t77 = t333 * t194 - t219 * t195 - t347;
t76 = t262 * t194 + t298 - t304 - t344;
t75 = (-rSges(5,1) * t305 - t335) * qJD(1) + t274;
t61 = t153 * t286 + t295 + t310;
t60 = t153 * t287 + t174 + (-t136 - t278) * t195;
t59 = t195 * t369;
t57 = -t115 * t194 + t97;
t56 = (t267 * t195 + t218) * qJD(1) - t272;
t55 = t185 + t191 + t366 + (t244 * t194 + t350 * t195) * qJD(1);
t53 = t249 * t195 + t300;
t43 = t54 + t295;
t42 = t174 + (t249 - t278) * t195 + t300;
t40 = (t263 * t195 + t246) * qJD(1) + t274 + t296;
t39 = t191 + t273 - t369;
t38 = t121 * t194 - t222 * t195;
t37 = t120 * t194 - t365;
t36 = t121 * t195 + t215;
t35 = t120 * t195 + t223 * t194;
t26 = (t262 * t195 - t334 - t347) * qJD(1) + t362;
t25 = -t281 + t191 + (t152 * t200 - t140) * t195 + (t333 * t195 + t384) * qJD(1);
t24 = t331 * t194 + t95 + t97;
t23 = t332 * t194 + t340;
t22 = ((-t127 + t336) * qJD(1) + t272) * t194 + (-t366 + (-t126 + t294) * qJD(1)) * t195;
t21 = -t194 * t75 + t59 + (-t115 * t195 - t117 * t194) * qJD(1);
t8 = t280 * t194 + t340 + t95;
t7 = t59 + t78 + (-t75 - t88) * t194 + (t331 * t195 + (-t117 - t99) * t194) * qJD(1);
t6 = t341 * t194 + (-t194 * t385 + t332 * t195) * qJD(1) + t342;
t5 = t78 + (-t88 + t341) * t194 + (t280 * t195 + (-t99 - t385) * t194) * qJD(1) + t342;
t1 = [t205 * t214 + t239 * t283 - t203 * t213 + t235 * t282 + (t55 * t92 + t56 * t91) * t355 + (t39 * t80 + t40 * t79) * t354 + (t25 * t77 + t26 * t76) * t353 - t402 * t302 + t401 * t301 + t404 * t197 + t403 * t196; 0; 0; m(4) * ((t194 * t55 - t195 * t56) * t179 + (t194 * t92 - t195 * t91) * t162) + (-t222 * qJD(3) + t203 * (t124 * qJD(1) - t238 * t284) + t205 * (t122 * qJD(1) - t234 * t284)) * t352 + m(6) * (t25 * t93 + t26 * t94 + t42 * t76 + t43 * t77) + m(5) * (t118 * t39 + t119 * t40 + t60 * t79 + t61 * t80) + (-t223 * qJD(3) + t203 * (-qJD(1) * t125 + t194 * t214) + t205 * (-qJD(1) * t123 + t194 * t213)) * t351 + t208 + ((t91 * t349 - t314 / 0.2e1 - t312 / 0.2e1) * t194 + (t92 * t349 + t313 / 0.2e1 + t311 / 0.2e1) * t195) * qJD(1) + (t193 / 0.2e1 + t192 / 0.2e1) * t231 * qJD(3); m(4) * t22 + m(5) * t7 + m(6) * t5; (t42 * t94 + t43 * t93 + t5 * t8) * t353 + (t118 * t61 + t119 * t60 + t24 * t7) * t354 + ((-t126 * t194 + t127 * t195) * t22 + t293 * t179 * t162) * t355 + (t194 * t38 + t195 * t37) * t286 + t194 * ((t194 * t82 + (-t37 + t215) * qJD(1)) * t194 + (t38 * qJD(1) + (-t122 * t282 - t124 * t283 + t290) * t195 + (t83 + (-t311 - t313) * qJD(3) + t223 * qJD(1)) * t194) * t195) + t195 * ((t195 * t83 + (t36 + t365) * qJD(1)) * t195 + (-t35 * qJD(1) + (t123 * t282 + t125 * t283) * t194 + (t82 + (t312 + t314) * qJD(3) + (-t120 + t222) * qJD(1)) * t195) * t194) + t248 + (-t194 * t36 - t195 * t35 + t368) * t287; m(6) * (t100 * t25 + t101 * t26 + t53 * t76 + t54 * t77) + t208 + (t194 * t39 - t195 * t40 + (t194 * t79 + t195 * t80) * qJD(1)) * t348 + m(5) * (t194 * t80 - t195 * t79) * t136; m(5) * t21 + m(6) * t6; m(6) * (t100 * t43 + t101 * t42 + t23 * t5 + t53 * t94 + t54 * t93 + t6 * t8) + m(5) * (-t119 * t136 * t195 + t118 * t310 + t21 * t24 + t57 * t7) + (t194 * t61 - t195 * t60 + (t118 * t195 + t119 * t194) * qJD(1)) * t348 + t209; (t293 * t153 * t136 + t21 * t57) * t354 + (t100 * t54 + t101 * t53 + t23 * t6) * t353 + t209; m(6) * (t194 * t26 + t195 * t25 + (-t194 * t77 + t195 * t76) * qJD(1)); 0; m(6) * (t194 * t42 + t195 * t43 + (-t194 * t93 + t195 * t94) * qJD(1)); m(6) * (t194 * t53 + t195 * t54 + (-t100 * t194 + t101 * t195) * qJD(1)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
