% Calculate time derivative of joint inertia matrix for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRRP3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:32
% EndTime: 2019-12-05 16:43:53
% DurationCPUTime: 9.48s
% Computational Cost: add. (9158->367), mult. (7808->508), div. (0->0), fcn. (6030->6), ass. (0->224)
t409 = Icges(5,4) + Icges(6,4);
t408 = Icges(5,1) + Icges(6,1);
t407 = -Icges(5,2) - Icges(6,2);
t195 = qJ(3) + qJ(4);
t190 = cos(t195);
t406 = t409 * t190;
t189 = sin(t195);
t405 = t409 * t189;
t404 = Icges(5,5) + Icges(6,5);
t403 = Icges(5,6) + Icges(6,6);
t402 = t407 * t189 + t406;
t401 = t408 * t190 - t405;
t400 = Icges(5,3) + Icges(6,3);
t193 = pkin(8) + qJ(2);
t187 = sin(t193);
t188 = cos(t193);
t372 = t403 * t187 + t402 * t188;
t370 = t404 * t187 + t401 * t188;
t398 = t407 * t190 - t405;
t397 = t408 * t189 + t406;
t399 = -t403 * t189 + t404 * t190;
t373 = t402 * t187 - t403 * t188;
t371 = t401 * t187 - t404 * t188;
t358 = t404 * t189 + t403 * t190;
t395 = t399 * t187 - t400 * t188;
t387 = t400 * t187 + t399 * t188;
t396 = t372 * t189 - t370 * t190;
t194 = qJD(3) + qJD(4);
t392 = t397 * t194;
t391 = t398 * t194;
t394 = t358 * t194;
t393 = t373 * t189 - t371 * t190;
t390 = t372 * qJD(2) + t391 * t187;
t389 = -t370 * qJD(2) + t392 * t187;
t386 = t402 * t194;
t385 = t401 * t194;
t384 = -t398 * t189 - t397 * t190;
t382 = t387 * qJD(2);
t381 = t395 * qJD(2);
t380 = t396 * t187;
t379 = t393 * t187 + t188 * t395;
t378 = -t188 * t387 - t380;
t377 = -t188 * t394 - t381;
t376 = t187 * t394 - t382;
t375 = -qJD(2) * t373 + t188 * t391;
t374 = qJD(2) * t371 + t188 * t392;
t197 = cos(qJ(3));
t184 = t197 * pkin(3) + pkin(2);
t163 = pkin(4) * t190 + t184;
t180 = t187 * rSges(6,3);
t292 = t188 * t190;
t293 = t188 * t189;
t369 = rSges(6,1) * t292 - rSges(6,2) * t293 + t188 * t163 + t180;
t368 = t371 * t194 + t390;
t367 = -t373 * t194 - t389;
t366 = t393 * t188;
t196 = sin(qJ(3));
t272 = qJD(3) * t196;
t266 = pkin(3) * t272;
t291 = t189 * t194;
t147 = -pkin(4) * t291 - t266;
t154 = rSges(6,1) * t189 + rSges(6,2) * t190;
t216 = t154 * t194;
t365 = qJD(5) * t188 - (t147 - t216) * t187;
t364 = -t187 * t395 + t366;
t363 = t387 * t187 - t188 * t396;
t323 = rSges(6,1) * t190;
t244 = -rSges(6,2) * t189 + t323;
t198 = -pkin(7) - pkin(6);
t192 = -qJ(5) + t198;
t281 = t192 - t198;
t285 = t163 - t184;
t91 = t285 * t187 + t281 * t188;
t362 = -rSges(6,3) * t188 + t244 * t187 + t91;
t164 = t188 * t184;
t361 = -t281 * t187 - t164 + t369;
t155 = rSges(5,1) * t189 + rSges(5,2) * t190;
t217 = t155 * t194;
t360 = t187 * t217;
t357 = (t385 + t391) * t190 + (-t386 - t392) * t189 + t358 * qJD(2);
t355 = -t370 * t194 - t375;
t354 = -t372 * t194 - t374;
t353 = qJD(2) * t198 + t266;
t275 = qJD(2) * t189;
t263 = t187 * t275;
t276 = qJD(2) * t188;
t352 = rSges(6,2) * t263 + rSges(6,3) * t276 + qJD(5) * t187 + t188 * t147;
t351 = t384 * qJD(2) + t399 * t194;
t290 = t190 * t194;
t350 = (t376 * t188 + (-t366 - t378) * qJD(2)) * t188 + ((t375 * t189 + t374 * t190 + t372 * t290 + t370 * t291 - t382) * t187 + (-t368 * t189 + t367 * t190 + t377) * t188 + ((t396 + t395) * t188 + t379) * qJD(2)) * t187;
t322 = rSges(4,2) * t196;
t325 = rSges(4,1) * t197;
t246 = -t322 + t325;
t321 = rSges(4,3) * t188;
t134 = t246 * t187 - t321;
t267 = t188 * t322;
t182 = t187 * rSges(4,3);
t283 = t188 * t325 + t182;
t135 = -t267 + t283;
t173 = rSges(4,1) * t196 + rSges(4,2) * t197;
t215 = qJD(3) * t173;
t274 = qJD(2) * t196;
t262 = t187 * t274;
t202 = rSges(4,2) * t262 + rSges(4,3) * t276 - t188 * t215;
t348 = t187 * t215;
t22 = (qJD(2) * t134 + t202) * t188 + (-t348 + (-t135 - t267 + t182) * qJD(2)) * t187;
t337 = 2 * m(4);
t349 = t22 * t337;
t317 = Icges(4,4) * t197;
t237 = -Icges(4,2) * t196 + t317;
t131 = Icges(4,6) * t187 + t237 * t188;
t318 = Icges(4,4) * t196;
t241 = Icges(4,1) * t197 - t318;
t133 = Icges(4,5) * t187 + t241 * t188;
t225 = t131 * t196 - t133 * t197;
t347 = t225 * t188;
t328 = pkin(2) - t184;
t342 = t328 * t187;
t233 = Icges(4,5) * t197 - Icges(4,6) * t196;
t128 = -Icges(4,3) * t188 + t233 * t187;
t130 = -Icges(4,6) * t188 + t237 * t187;
t132 = -Icges(4,5) * t188 + t241 * t187;
t336 = 2 * m(5);
t335 = 2 * m(6);
t185 = t187 ^ 2;
t186 = t188 ^ 2;
t334 = t187 / 0.2e1;
t333 = -t188 / 0.2e1;
t332 = m(4) * t173;
t331 = m(5) * t155;
t330 = pkin(3) * t196;
t329 = t187 * pkin(6);
t183 = t188 * pkin(6);
t327 = -pkin(6) - t198;
t106 = t188 * t198 + t183 - t342;
t107 = -t188 * pkin(2) + t327 * t187 + t164;
t326 = t187 * t106 + t188 * t107;
t324 = rSges(5,1) * t190;
t181 = t187 * rSges(5,3);
t320 = rSges(6,3) - t192;
t303 = t130 * t197;
t302 = t131 * t197;
t301 = t132 * t196;
t300 = t133 * t196;
t245 = -rSges(5,2) * t189 + t324;
t144 = t245 * t194;
t299 = t144 * t187;
t294 = t187 * t192;
t123 = -t188 * rSges(5,3) + t245 * t187;
t125 = rSges(5,1) * t292 - rSges(5,2) * t293 + t181;
t55 = t187 * t123 + t188 * t125;
t277 = qJD(2) * t187;
t288 = pkin(4) * t263 + t154 * t277;
t286 = rSges(5,2) * t263 + rSges(5,3) * t276;
t284 = t353 * t187;
t282 = t185 + t186;
t129 = Icges(4,3) * t187 + t233 * t188;
t278 = qJD(2) * t129;
t271 = qJD(3) * t197;
t249 = t188 * t266;
t269 = t187 * ((-t328 * t188 - t329) * qJD(2) - t284) + t188 * (-t249 + (t327 * t188 + t342) * qJD(2)) + t106 * t276;
t205 = -t188 * t291 - t190 * t277;
t264 = t188 * t290;
t268 = t123 * t276 + t187 * (-t360 + (t245 * t188 + t181) * qJD(2)) + t188 * (rSges(5,1) * t205 - rSges(5,2) * t264 + t286);
t265 = pkin(3) * t271;
t259 = -t155 - t330;
t258 = -pkin(4) * t189 - t154;
t23 = t362 * t187 + t361 * t188;
t143 = t244 * t194;
t248 = -pkin(4) * t290 - t143;
t247 = (t364 * t276 + t379 * t277) * t188 + ((t377 * t187 + (-t364 + t380) * qJD(2)) * t187 + t378 * t277 + t363 * t276 + (t376 * t187 + (t373 * t290 + t371 * t291 - t381) * t188 + (t354 * t187 + t389 * t188) * t190 + (t355 * t187 + t390 * t188) * t189 + ((t387 - t393) * t187 + t363) * qJD(2)) * t188) * t187;
t240 = Icges(4,1) * t196 + t317;
t236 = Icges(4,2) * t197 + t318;
t226 = t130 * t196 - t132 * t197;
t222 = t362 * t276 + (rSges(6,1) * t205 - rSges(6,2) * t264 - t91 * qJD(2) + t249 + t352) * t188 + (t284 + (t180 - t294 + (t244 + t285) * t188) * qJD(2) - t365) * t187;
t221 = -pkin(2) - t246;
t220 = t258 - t330;
t219 = -t184 - t245;
t218 = -t163 - t244;
t208 = qJD(3) * t240;
t207 = qJD(3) * t236;
t206 = qJD(3) * (-Icges(4,5) * t196 - Icges(4,6) * t197);
t95 = t220 * t188;
t204 = t248 - t265;
t203 = t350 * t188 + t247;
t201 = (t351 * t187 + t357 * t188 + t354 * t189 - t355 * t190) * t334 + (t357 * t187 - t351 * t188 + t367 * t189 + t368 * t190) * t333 + (-t187 * t384 - t358 * t188 + t371 * t189 + t373 * t190) * t277 / 0.2e1 + (t358 * t187 - t188 * t384 + t370 * t189 + t372 * t190) * t276 / 0.2e1;
t169 = pkin(3) * t262;
t162 = t246 * qJD(3);
t127 = t259 * t188;
t126 = t259 * t187;
t109 = t258 * t188;
t108 = t258 * t187;
t97 = t329 + (pkin(2) - t322) * t188 + t283;
t96 = t221 * t187 + t183 + t321;
t94 = t220 * t187;
t90 = -t187 * t198 + t125 + t164;
t89 = (rSges(5,3) - t198) * t188 + t219 * t187;
t84 = t187 * t206 + t278;
t83 = -qJD(2) * t128 + t188 * t206;
t77 = -t294 + t369;
t76 = t218 * t187 + t320 * t188;
t63 = t348 + ((-rSges(4,3) - pkin(6)) * t187 + t221 * t188) * qJD(2);
t62 = (t183 + (-pkin(2) - t325) * t187) * qJD(2) + t202;
t61 = -t155 * t276 - t299 + (-t187 * t271 - t188 * t274) * pkin(3);
t60 = t155 * t277 + t169 + (-t144 - t265) * t188;
t54 = -t154 * t276 - t143 * t187 + (-t187 * t290 - t188 * t275) * pkin(4);
t53 = t248 * t188 + t288;
t44 = qJD(2) * t95 + t187 * t204;
t43 = t188 * t204 + t169 + t288;
t42 = t360 + (t188 * t219 - t181) * qJD(2) + t284;
t41 = (-t184 - t324) * t277 + (-t217 - t353) * t188 + t286;
t38 = t129 * t187 - t347;
t37 = t128 * t187 - t226 * t188;
t36 = -t129 * t188 - t225 * t187;
t35 = -t128 * t188 - t226 * t187;
t26 = (-t320 * t187 + t218 * t188) * qJD(2) + t365;
t25 = -t188 * t216 + (-t188 * t192 + (-t163 - t323) * t187) * qJD(2) + t352;
t24 = t55 + t326;
t21 = -t125 * t277 + t268;
t8 = t23 + t326;
t7 = (-t107 - t125) * t277 + t268 + t269;
t6 = -t277 * t361 + t222;
t5 = (-t107 - t361) * t277 + t222 + t269;
t1 = [0; 0; (t25 * t77 + t26 * t76) * t335 + (t41 * t90 + t42 * t89) * t336 + (t62 * t97 + t63 * t96) * t337 + t398 * t291 + t397 * t290 + (t241 - t236) * t272 + (t240 + t237) * t271 + t386 * t190 + t385 * t189; m(4) * t22 + m(5) * t7 + m(6) * t5; (-t225 * qJD(3) + t196 * (-t132 * qJD(2) - t188 * t208) + t197 * (-t130 * qJD(2) - t188 * t207)) * t334 + (-t226 * qJD(3) + t196 * (qJD(2) * t133 - t187 * t208) + t197 * (qJD(2) * t131 - t187 * t207)) * t333 + m(5) * (t126 * t41 + t127 * t42 + t60 * t89 + t61 * t90) + m(6) * (t25 * t94 + t26 * t95 + t43 * t76 + t44 * t77) + m(4) * ((-t187 * t62 - t188 * t63) * t173 + (-t187 * t97 - t188 * t96) * t162) + (t186 / 0.2e1 + t185 / 0.2e1) * t233 * qJD(3) + t201 + ((-t97 * t332 + t302 / 0.2e1 + t300 / 0.2e1) * t188 + (t96 * t332 + t303 / 0.2e1 + t301 / 0.2e1) * t187) * qJD(2); (t43 * t95 + t44 * t94 + t5 * t8) * t335 + (t126 * t61 + t127 * t60 + t24 * t7) * t336 + t282 * t173 * t162 * t337 + t247 + (t134 * t349 + t38 * t276 + t36 * t277 + (t37 * qJD(2) + (t225 * qJD(2) + t83) * t187) * t187) * t187 + (-t35 * t277 + t135 * t349 - t37 * t276 + (-t36 * qJD(2) + (-t226 * qJD(2) - t84) * t188) * t188 + ((t83 - (t301 + t303) * qJD(3) + t130 * t271 + t132 * t272) * t188 + (t131 * t271 + t133 * t272 - t278 - t84 + (-t300 - t302) * qJD(3)) * t187 + (-t35 + t38 + t347 + (t129 - t226) * t187) * qJD(2)) * t187 + t350) * t188; m(5) * t21 + m(6) * t6; m(5) * (-t187 * t90 - t188 * t89) * t144 + t201 + (-t187 * t41 - t188 * t42 + (t187 * t89 - t188 * t90) * qJD(2)) * t331 + m(6) * (t108 * t25 + t109 * t26 + t53 * t76 + t54 * t77); m(6) * (t108 * t44 + t109 * t43 + t23 * t5 + t53 * t95 + t54 * t94 + t6 * t8) + m(5) * (-t127 * t144 * t188 - t126 * t299 + t21 * t24 + t55 * t7) + (-t187 * t61 - t188 * t60 + (-t126 * t188 + t127 * t187) * qJD(2)) * t331 + t203; (t144 * t155 * t282 + t21 * t55) * t336 + (t108 * t54 + t109 * t53 + t23 * t6) * t335 + t203; 0; m(6) * (t187 * t26 - t188 * t25 + (t187 * t77 + t188 * t76) * qJD(2)); m(6) * (t187 * t43 - t188 * t44 + (t187 * t94 + t188 * t95) * qJD(2)); m(6) * (t187 * t53 - t188 * t54 + (t108 * t187 + t109 * t188) * qJD(2)); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
