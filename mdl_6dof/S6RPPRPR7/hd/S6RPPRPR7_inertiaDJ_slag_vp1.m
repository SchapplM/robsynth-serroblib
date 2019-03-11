% Calculate time derivative of joint inertia matrix for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR7_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:54
% EndTime: 2019-03-09 01:53:11
% DurationCPUTime: 10.16s
% Computational Cost: add. (16049->697), mult. (17822->1001), div. (0->0), fcn. (16551->10), ass. (0->332)
t425 = Icges(6,5) / 0.2e1;
t424 = Icges(6,6) / 0.2e1;
t423 = Icges(6,3) / 0.2e1;
t237 = sin(qJ(1));
t397 = -t237 / 0.2e1;
t422 = -qJD(1) / 0.2e1;
t421 = rSges(6,3) + qJ(5);
t228 = pkin(9) + qJ(4);
t217 = sin(t228);
t219 = cos(t228);
t336 = qJD(4) * t237;
t308 = t219 * t336;
t238 = cos(qJ(1));
t339 = qJD(1) * t238;
t420 = t217 * t339 + t308;
t335 = qJD(4) * t238;
t309 = t217 * t335;
t340 = qJD(1) * t237;
t419 = t219 * t340 + t309;
t310 = t217 * t336;
t311 = t219 * t339;
t242 = t310 - t311;
t227 = pkin(10) + qJ(6);
t216 = sin(t227);
t218 = cos(t227);
t374 = Icges(7,4) * t218;
t265 = -Icges(7,2) * t216 + t374;
t137 = Icges(7,6) * t217 + t219 * t265;
t375 = Icges(7,4) * t216;
t269 = Icges(7,1) * t218 - t375;
t138 = Icges(7,5) * t217 + t219 * t269;
t418 = -t137 * t216 + t138 * t218;
t231 = sin(pkin(10));
t233 = cos(pkin(10));
t270 = Icges(6,1) * t233 - Icges(6,4) * t231;
t145 = Icges(6,5) * t217 + t219 * t270;
t417 = -t145 / 0.2e1;
t213 = pkin(5) * t233 + pkin(4);
t235 = -pkin(8) - qJ(5);
t300 = qJD(6) * t217 + qJD(1);
t255 = t237 * t300;
t299 = qJD(1) * t217 + qJD(6);
t405 = t238 * t299 + t308;
t82 = -t216 * t405 - t218 * t255;
t83 = -t216 * t255 + t218 * t405;
t47 = t83 * rSges(7,1) + t82 * rSges(7,2) + rSges(7,3) * t242;
t416 = t213 * t420 + t235 * t311 + t47;
t355 = t231 * t238;
t211 = pkin(5) * t355;
t357 = t219 * t237;
t360 = t217 * t237;
t353 = t237 * t216;
t358 = t218 * t238;
t160 = -t217 * t353 + t358;
t352 = t237 * t218;
t161 = t216 * t238 + t217 * t352;
t94 = t161 * rSges(7,1) + t160 * rSges(7,2) - rSges(7,3) * t357;
t329 = -t213 * t360 - t235 * t357 - t211 - t94;
t295 = rSges(5,1) * t217 + rSges(5,2) * t219;
t248 = t238 * t295;
t377 = Icges(5,4) * t217;
t267 = Icges(5,2) * t219 + t377;
t150 = Icges(5,6) * t238 + t237 * t267;
t376 = Icges(5,4) * t219;
t271 = Icges(5,1) * t217 + t376;
t152 = Icges(5,5) * t238 + t237 * t271;
t258 = t150 * t219 + t152 * t217;
t247 = t258 * t238;
t389 = pkin(4) - t213;
t415 = t389 * t217;
t351 = t237 * t231;
t354 = t233 * t238;
t173 = -t217 * t351 + t354;
t350 = t237 * t233;
t174 = t217 * t350 + t355;
t205 = pkin(4) * t360;
t414 = -t174 * rSges(6,1) - t173 * rSges(6,2) - t205;
t413 = -t137 * t218 - t138 * t216;
t232 = sin(pkin(9));
t391 = pkin(3) * t232;
t212 = t238 * t391;
t224 = t238 * qJ(2);
t236 = -pkin(7) - qJ(3);
t316 = t237 * t236 + t212 + t224;
t393 = -rSges(5,3) - pkin(1);
t118 = t237 * t393 + t248 + t316;
t225 = t238 * rSges(5,3);
t157 = rSges(5,1) * t360 + rSges(5,2) * t357 + t225;
t226 = t238 * pkin(1);
t343 = t237 * qJ(2) + t226;
t253 = -t236 * t238 + t237 * t391 + t343;
t119 = t253 + t157;
t412 = -t118 * t237 + t119 * t238;
t344 = qJ(2) * t339 + qJD(2) * t237;
t315 = qJD(3) * t238 + t344;
t273 = qJD(1) * t212 + t236 * t340 + t315;
t317 = rSges(5,1) * t420 + rSges(5,2) * t311;
t338 = qJD(4) * t217;
t74 = (-rSges(5,2) * t338 + qJD(1) * t393) * t237 + t273 + t317;
t386 = rSges(5,2) * t217;
t189 = rSges(5,1) * t219 - t386;
t222 = qJD(2) * t238;
t301 = -qJD(3) * t237 + t222;
t297 = t236 * t339 + t301;
t306 = -qJ(2) - t391;
t75 = t189 * t335 + (t393 * t238 + (-t295 + t306) * t237) * qJD(1) + t297;
t411 = t237 * t75 - t238 * t74;
t364 = t213 * t217;
t379 = rSges(7,3) - t235;
t410 = -t219 * t379 + t364;
t307 = t219 * t335;
t409 = t237 * t299 - t307;
t263 = Icges(5,5) * t217 + Icges(5,6) * t219;
t408 = -Icges(5,3) * t237 + t238 * t263;
t407 = -Icges(5,6) * t237 + t238 * t267;
t406 = -Icges(5,5) * t237 + t238 * t271;
t290 = rSges(7,1) * t218 - rSges(7,2) * t216;
t140 = rSges(7,3) * t217 + t219 * t290;
t254 = t300 * t238;
t80 = -t216 * t409 + t218 * t254;
t81 = t216 * t254 + t218 * t409;
t298 = t81 * rSges(7,1) + t80 * rSges(7,2);
t46 = -rSges(7,3) * t419 + t298;
t332 = qJD(6) * t219;
t87 = (-rSges(7,1) * t216 - rSges(7,2) * t218) * t332 + (rSges(7,3) * t219 - t217 * t290) * qJD(4);
t359 = t217 * t238;
t162 = t216 * t359 + t352;
t163 = -t217 * t358 + t353;
t291 = -t163 * rSges(7,1) - t162 * rSges(7,2);
t356 = t219 * t238;
t95 = rSges(7,3) * t356 - t291;
t21 = (-t140 * t335 - t46) * t217 + (-qJD(4) * t95 - t140 * t340 + t238 * t87) * t219;
t22 = (-t140 * t336 + t47) * t217 + (qJD(4) * t94 + t140 * t339 + t237 * t87) * t219;
t58 = t140 * t357 + t217 * t94;
t59 = t140 * t356 - t217 * t95;
t404 = qJD(1) * (t237 * t58 + t238 * t59) + t21 * t237 - t22 * t238;
t403 = 2 * m(5);
t402 = 2 * m(6);
t401 = 2 * m(7);
t229 = t237 ^ 2;
t230 = t238 ^ 2;
t400 = m(6) / 0.2e1;
t399 = m(7) / 0.2e1;
t398 = t217 / 0.2e1;
t395 = t238 / 0.2e1;
t394 = rSges(3,2) - pkin(1);
t392 = m(5) * t189;
t390 = pkin(4) * t217;
t261 = Icges(7,5) * t218 - Icges(7,6) * t216;
t136 = Icges(7,3) * t217 + t219 * t261;
t337 = qJD(4) * t219;
t84 = (-Icges(7,5) * t216 - Icges(7,6) * t218) * t332 + (Icges(7,3) * t219 - t217 * t261) * qJD(4);
t86 = (-Icges(7,1) * t216 - t374) * t332 + (Icges(7,5) * t219 - t217 * t269) * qJD(4);
t240 = t219 * t218 * t86 + t136 * t337 + t217 * t84 - t338 * t418;
t85 = (-Icges(7,2) * t218 - t375) * t332 + (Icges(7,6) * t219 - t217 * t265) * qJD(4);
t384 = t216 * t85;
t53 = t136 * t217 + t219 * t418;
t388 = ((qJD(6) * t413 - t384) * t219 + t240) * t217 + t53 * t337;
t206 = pkin(4) * t359;
t349 = -qJ(5) - t235;
t303 = t349 * t219;
t330 = pkin(5) * t351;
t387 = t95 + t330 + t206 + (t303 - t364) * t238;
t385 = rSges(6,3) * t219;
t383 = t237 * rSges(5,3);
t380 = rSges(4,3) + qJ(3);
t378 = (t303 + t415) * qJD(4) + t87;
t363 = t213 * t219;
t362 = t217 * t231;
t361 = t217 * t233;
t348 = t357 * t421 + t414;
t347 = t217 * t349 - t219 * t389 + t140;
t154 = qJD(5) * t217 + (qJ(5) * t219 - t390) * qJD(4);
t188 = pkin(4) * t219 + qJ(5) * t217;
t346 = t237 * t154 + t188 * t339;
t342 = t229 + t230;
t148 = Icges(5,3) * t238 + t237 * t263;
t341 = qJD(1) * t148;
t334 = qJD(5) * t219;
t333 = qJD(5) * t237;
t331 = -pkin(1) - t380;
t90 = Icges(7,4) * t161 + Icges(7,2) * t160 - Icges(7,6) * t357;
t92 = Icges(7,1) * t161 + Icges(7,4) * t160 - Icges(7,5) * t357;
t288 = t216 * t90 - t218 * t92;
t88 = Icges(7,5) * t161 + Icges(7,6) * t160 - Icges(7,3) * t357;
t32 = t217 * t88 - t219 * t288;
t49 = -t136 * t357 + t137 * t160 + t138 * t161;
t328 = t32 / 0.2e1 + t49 / 0.2e1;
t91 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t356;
t93 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t356;
t287 = t216 * t91 - t218 * t93;
t89 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t356;
t33 = t217 * t89 - t219 * t287;
t50 = t136 * t356 + t162 * t137 + t163 * t138;
t327 = -t50 / 0.2e1 - t33 / 0.2e1;
t102 = Icges(6,5) * t174 + Icges(6,6) * t173 - Icges(6,3) * t357;
t326 = t102 * t357;
t325 = t102 * t356;
t175 = t217 * t355 + t350;
t249 = t217 * t354 - t351;
t103 = -Icges(6,5) * t249 + Icges(6,6) * t175 + Icges(6,3) * t356;
t324 = t103 * t357;
t323 = t103 * t356;
t129 = -qJD(1) * t175 - t231 * t308;
t130 = qJD(1) * t249 + t233 * t308;
t322 = -t130 * rSges(6,1) - t129 * rSges(6,2) - rSges(6,3) * t310;
t319 = pkin(4) * t420 + qJ(5) * t310;
t318 = pkin(4) * t307 + qJ(5) * t419;
t314 = -pkin(5) * t231 - pkin(1);
t305 = t219 * t421;
t181 = t295 * qJD(4);
t304 = t181 * t342;
t302 = qJD(1) * t347;
t296 = rSges(4,1) * t232 + rSges(4,2) * cos(pkin(9));
t127 = qJD(1) * t173 + t231 * t307;
t128 = qJD(1) * t174 - t233 * t307;
t294 = -t128 * rSges(6,1) - t127 * rSges(6,2);
t293 = rSges(6,1) * t249 - t175 * rSges(6,2);
t292 = rSges(6,1) * t233 - rSges(6,2) * t231;
t24 = (t314 * qJD(1) - t235 * t338 - t334) * t237 + t273 + t416;
t25 = (-t334 + (t217 * t379 + t363) * qJD(4)) * t238 + (t314 * t238 + (t306 - t410) * t237) * qJD(1) + t297 - t298;
t286 = t237 * t25 - t238 * t24;
t26 = t160 * t90 + t161 * t92 - t357 * t88;
t27 = t160 * t91 + t161 * t93 - t357 * t89;
t17 = t27 * t237 + t238 * t26;
t285 = t237 * t26 - t238 * t27;
t28 = t162 * t90 + t163 * t92 + t356 * t88;
t29 = t162 * t91 + t163 * t93 + t356 * t89;
t18 = t29 * t237 + t238 * t28;
t284 = t237 * t28 - t238 * t29;
t169 = t188 * t340;
t30 = t169 + t237 * t302 + (-t154 - t378) * t238;
t31 = t237 * t378 + t238 * t302 + t346;
t283 = t31 * t237 - t238 * t30;
t282 = t33 * t237 + t32 * t238;
t281 = t237 * t32 - t238 * t33;
t243 = -t237 * pkin(1) - t238 * t305;
t38 = qJD(1) * t243 - t219 * t333 + t273 + t319 - t322;
t39 = (rSges(6,3) * t338 - t334) * t238 + (-t226 + (t306 + t385 - t390) * t237) * qJD(1) + t294 + t297 + t318;
t280 = t237 * t39 - t238 * t38;
t54 = t314 * t237 + t238 * t410 + t291 + t316;
t55 = t253 - t329;
t279 = t237 * t54 - t238 * t55;
t139 = (-t217 * t292 + t385) * qJD(4);
t147 = rSges(6,3) * t217 + t219 * t292;
t56 = t147 * t340 + t169 + (-t139 - t154) * t238;
t57 = t237 * t139 + t147 * t339 + t346;
t278 = t57 * t237 - t238 * t56;
t277 = t237 * t59 - t238 * t58;
t70 = t206 + t243 + t293 + t316;
t71 = -t237 * t305 + t253 - t414;
t275 = t237 * t70 - t238 * t71;
t274 = t237 * t95 + t238 * t94;
t272 = Icges(5,1) * t219 - t377;
t268 = -Icges(5,2) * t217 + t376;
t266 = Icges(6,4) * t233 - Icges(6,2) * t231;
t264 = Icges(5,5) * t219 - Icges(5,6) * t217;
t262 = Icges(6,5) * t233 - Icges(6,6) * t231;
t257 = -t217 * t406 - t219 * t407;
t256 = (t400 + t399) * t338;
t252 = rSges(3,3) * t238 + t237 * t394;
t246 = t257 * t237;
t245 = qJD(4) * t272;
t244 = qJD(4) * t268;
t239 = t237 * t331 + t238 * t296;
t177 = t237 * t188;
t172 = -rSges(3,2) * t238 + t237 * rSges(3,3) + t343;
t171 = t224 + t252;
t168 = qJ(5) * t356 - t206;
t159 = t238 * t168;
t158 = t383 - t248;
t144 = Icges(6,6) * t217 + t219 * t266;
t142 = t222 + (t394 * t238 + (-rSges(3,3) - qJ(2)) * t237) * qJD(1);
t141 = qJD(1) * t252 + t344;
t135 = (Icges(6,5) * t219 - t217 * t270) * qJD(4);
t134 = (Icges(6,6) * t219 - t217 * t266) * qJD(4);
t132 = t237 * t296 + t238 * t380 + t343;
t131 = t224 + t239;
t117 = (-t147 - t188) * t238;
t116 = t147 * t237 + t177;
t115 = rSges(6,3) * t356 - t293;
t109 = qJD(1) * t408 + t264 * t336;
t108 = -t264 * t335 + t341;
t107 = -Icges(6,1) * t249 + Icges(6,4) * t175 + Icges(6,5) * t356;
t106 = Icges(6,1) * t174 + Icges(6,4) * t173 - Icges(6,5) * t357;
t105 = -Icges(6,4) * t249 + Icges(6,2) * t175 + Icges(6,6) * t356;
t104 = Icges(6,4) * t174 + Icges(6,2) * t173 - Icges(6,6) * t357;
t101 = (t331 * t238 + (-qJ(2) - t296) * t237) * qJD(1) + t301;
t100 = qJD(1) * t239 + t315;
t97 = (-qJ(5) * t339 - t333) * t219 + t319;
t96 = t238 * (qJD(1) * t205 + t238 * t334 - t318);
t73 = (-t188 - t347) * t238;
t72 = t237 * t347 + t177;
t69 = -t237 * t408 - t238 * t257;
t68 = t237 * t148 - t247;
t67 = -t238 * t408 + t246;
t66 = t148 * t238 + t237 * t258;
t65 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t242;
t64 = Icges(6,1) * t128 + Icges(6,4) * t127 - Icges(6,5) * t419;
t63 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t242;
t62 = Icges(6,4) * t128 + Icges(6,2) * t127 - Icges(6,6) * t419;
t51 = t274 * t219;
t48 = t115 * t238 + t237 * t348 + t159;
t45 = Icges(7,1) * t83 + Icges(7,4) * t82 + Icges(7,5) * t242;
t44 = Icges(7,1) * t81 + Icges(7,4) * t80 - Icges(7,5) * t419;
t43 = Icges(7,4) * t83 + Icges(7,2) * t82 + Icges(7,6) * t242;
t42 = Icges(7,4) * t81 + Icges(7,2) * t80 - Icges(7,6) * t419;
t41 = Icges(7,5) * t83 + Icges(7,6) * t82 + Icges(7,3) * t242;
t40 = Icges(7,5) * t81 + Icges(7,6) * t80 - Icges(7,3) * t419;
t37 = t175 * t105 - t107 * t249 + t323;
t36 = t175 * t104 - t106 * t249 + t325;
t35 = t105 * t173 + t107 * t174 - t324;
t34 = t104 * t173 + t106 * t174 - t326;
t23 = t237 * t329 + t238 * t387 + t159;
t19 = t96 + t238 * (-rSges(6,3) * t309 - t294) + (-t97 + t322) * t237 + (t348 * t238 + (-t115 - t168) * t237) * qJD(1);
t16 = t136 * t242 + t82 * t137 + t83 * t138 + t160 * t85 + t161 * t86 - t357 * t84;
t15 = -t136 * t419 + t80 * t137 + t81 * t138 + t162 * t85 + t163 * t86 + t356 * t84;
t14 = t274 * t338 + (-t237 * t46 - t238 * t47 + (t237 * t94 - t238 * t95) * qJD(1)) * t219;
t13 = t50 * t217 - t219 * t284;
t12 = t49 * t217 - t219 * t285;
t11 = (qJD(4) * t287 + t40) * t217 + (qJD(4) * t89 - t216 * t42 + t218 * t44 + (-t216 * t93 - t218 * t91) * qJD(6)) * t219;
t10 = (qJD(4) * t288 + t41) * t217 + (qJD(4) * t88 - t216 * t43 + t218 * t45 + (-t216 * t92 - t218 * t90) * qJD(6)) * t219;
t9 = t89 * t310 + t160 * t42 + t161 * t44 + t82 * t91 + t83 * t93 + (-t237 * t40 - t339 * t89) * t219;
t8 = t88 * t310 + t160 * t43 + t161 * t45 + t82 * t90 + t83 * t92 + (-t237 * t41 - t339 * t88) * t219;
t7 = -t89 * t309 + t162 * t42 + t163 * t44 + t80 * t91 + t81 * t93 + (t238 * t40 - t340 * t89) * t219;
t6 = -t88 * t309 + t162 * t43 + t163 * t45 + t80 * t90 + t81 * t92 + (t238 * t41 - t340 * t88) * t219;
t5 = t96 + (t46 + (t217 * t235 - t363) * t335 + t318) * t238 + (t235 * t310 + t319 - t416 - t97) * t237 + ((t329 + t211) * t238 + (t330 - t168 + ((-qJ(5) + t235) * t219 - t415) * t238 - t387) * t237) * qJD(1);
t4 = -qJD(1) * t285 + t9 * t237 + t238 * t8;
t3 = -qJD(1) * t284 + t7 * t237 + t238 * t6;
t2 = (qJD(4) * t285 + t16) * t217 + (-qJD(1) * t17 + qJD(4) * t49 - t237 * t8 + t238 * t9) * t219;
t1 = (qJD(4) * t284 + t15) * t217 + (-qJD(1) * t18 + qJD(4) * t50 - t237 * t6 + t238 * t7) * t219;
t20 = [t240 - t271 * t337 + (t24 * t55 + t25 * t54) * t401 + (t38 * t71 + t39 * t70) * t402 + (t118 * t75 + t119 * t74) * t403 + 0.2e1 * m(4) * (t100 * t132 + t101 * t131) + 0.2e1 * m(3) * (t141 * t172 + t142 * t171) + t413 * t332 + (Icges(6,3) * t337 - t245) * t217 + (t144 * t231 - t145 * t233 - t217 * t262 + t267) * t338 + (Icges(6,3) * t338 - t134 * t231 + t135 * t233 + t262 * t337 - t244 - t384) * t219; m(7) * ((t237 * t55 + t238 * t54) * qJD(1) + t286) + m(6) * ((t237 * t71 + t238 * t70) * qJD(1) + t280) + m(5) * ((t118 * t238 + t119 * t237) * qJD(1) + t411) + m(4) * (-t100 * t238 + t237 * t101 + (t131 * t238 + t132 * t237) * qJD(1)) + m(3) * (-t141 * t238 + t237 * t142 + (t171 * t238 + t172 * t237) * qJD(1)); 0; m(7) * (-qJD(1) * t279 + t237 * t24 + t238 * t25) + m(6) * (-qJD(1) * t275 + t237 * t38 + t238 * t39) + m(5) * (qJD(1) * t412 + t237 * t74 + t238 * t75) + m(4) * (t237 * t100 + t101 * t238 + (-t131 * t237 + t132 * t238) * qJD(1)); 0; 0; m(6) * (t116 * t39 + t117 * t38 + t56 * t71 + t57 * t70) + m(7) * (t24 * t73 + t25 * t72 + t30 * t55 + t31 * t54) + m(5) * (t181 * t412 + t189 * t411) + ((t129 * t424 + t130 * t425 + t242 * t423 + t244 * t397 + t407 * t422) * t238 + (t150 * t422 + t268 * t335 / 0.2e1 + t128 * t425 + t127 * t424 - t419 * t423) * t237) * t217 + ((-t144 * t173 / 0.2e1 + t174 * t417 + t119 * t392 + (t150 / 0.2e1 - t102 / 0.2e1) * t217 + (t104 * t231 / 0.2e1 - t106 * t233 / 0.2e1) * t219 - t328) * t237 + (t175 * t144 / 0.2e1 + t249 * t417 + t118 * t392 + (t407 / 0.2e1 + t103 / 0.2e1) * t217 + (-t406 / 0.2e1 - t105 * t231 / 0.2e1 + t107 * t233 / 0.2e1) * t219 - t327) * t238) * qJD(1) + (-(t229 / 0.2e1 + t230 / 0.2e1) * t263 - t247 / 0.2e1 + t257 * t397) * qJD(4) + (t127 * t144 + t128 * t145 + t175 * t134 - t249 * t135 + t11 + t15 + (-t231 * t62 + t233 * t64 - t272 * t335) * t219 + (t103 * t219 + t105 * t362 - t107 * t361) * qJD(4)) * t237 / 0.2e1 + (t129 * t144 + t130 * t145 + t173 * t134 + t174 * t135 + t16 + t10 + (qJD(1) * t406 - t231 * t63 + t233 * t65 + t237 * t245) * t219 + (t102 * t219 + t104 * t362 - t106 * t361) * qJD(4)) * t395; m(6) * ((t116 * t238 + t117 * t237) * qJD(1) + t278) + m(7) * ((t237 * t73 + t238 * t72) * qJD(1) + t283) - m(5) * t304; m(6) * (t56 * t237 + t238 * t57 + (-t116 * t237 + t117 * t238) * qJD(1)) + m(7) * (t30 * t237 + t238 * t31 + (-t237 * t72 + t238 * t73) * qJD(1)); (t23 * t5 + t30 * t73 + t31 * t72) * t401 + t238 * t4 + t237 * t3 + (t116 * t57 + t117 * t56 + t48 * t19) * t402 + t237 * ((t237 * t108 + (-t68 + t246) * qJD(1)) * t237 + (t69 * qJD(1) + (t150 * t338 - t152 * t337 + t341) * t238 + (t109 + (-t217 * t407 + t219 * t406) * qJD(4) + t258 * qJD(1)) * t237) * t238) + ((-t237 * t157 + t158 * t238) * (-t237 * t317 + (-t189 * t230 + t229 * t386) * qJD(4) + ((-t157 + t225) * t238 + (-t158 + t248 + t383) * t237) * qJD(1)) - t189 * t304) * t403 + t238 * ((t109 * t238 + (t67 + t247) * qJD(1)) * t238 + (-t66 * qJD(1) + (-t337 * t406 + t338 * t407) * t237 + (t108 + (-t150 * t217 + t152 * t219) * qJD(4) + (-t148 + t257) * qJD(1)) * t238) * t237) + t238 * ((t129 * t104 + t130 * t106 + t173 * t63 + t174 * t65 + (t35 - t325) * qJD(1)) * t238 + (t129 * t105 + t130 * t107 + t173 * t62 + t174 * t64 + (-t34 - t323) * qJD(1)) * t237) + t237 * ((t127 * t105 + t128 * t107 + t175 * t62 - t249 * t64 + (-t36 - t324) * qJD(1)) * t237 + (t127 * t104 + t128 * t106 + t175 * t63 - t249 * t65 + (t37 - t326) * qJD(1)) * t238) + (-t17 + (-t34 - t66) * t238 + (-t35 - t67) * t237) * t340 + (t18 + (t36 + t68) * t238 + (t37 + t69) * t237) * t339; 0.2e1 * (t275 * t400 + t279 * t399) * t338 + 0.2e1 * ((-t339 * t54 - t340 * t55 - t286) * t399 + (-t339 * t70 - t340 * t71 - t280) * t400) * t219; 0.2e1 * t342 * t256; 0; 0.2e1 * ((-t335 * t73 + t336 * t72 + t5) * t399 + (t116 * t336 - t117 * t335 + t19) * t400) * t217 + 0.2e1 * ((qJD(4) * t23 - t339 * t72 - t340 * t73 - t283) * t399 + (qJD(4) * t48 - t116 * t339 - t117 * t340 - t278) * t400) * t219; 0.4e1 * (0.1e1 - t342) * t219 * t256; m(7) * (t21 * t54 + t22 * t55 + t24 * t58 + t25 * t59) + (t237 * t328 + t238 * t327) * t338 + ((t15 / 0.2e1 + t11 / 0.2e1) * t238 + (-t10 / 0.2e1 - t16 / 0.2e1) * t237 + (t237 * t327 - t238 * t328) * qJD(1)) * t219 + t388; m(7) * t404; m(7) * (-qJD(1) * t277 + t21 * t238 + t22 * t237); m(7) * (t14 * t23 + t21 * t72 + t22 * t73 + t30 * t58 + t31 * t59 - t5 * t51) + ((qJD(1) * t33 + t10) * t398 + t2 / 0.2e1 + qJD(1) * t13 / 0.2e1 - t18 * t338 / 0.2e1) * t238 + ((-qJD(1) * t32 + t11) * t398 + t12 * t422 + t1 / 0.2e1 + t17 * t338 / 0.2e1) * t237 + (qJD(4) * t282 / 0.2e1 + t3 * t395 + t4 * t397 + (t18 * t397 - t238 * t17 / 0.2e1) * qJD(1)) * t219; m(7) * ((qJD(4) * t277 + t14) * t217 + (-qJD(4) * t51 - t404) * t219); (-t14 * t51 + t21 * t59 + t22 * t58) * t401 + ((t237 * t12 - t238 * t13 + t217 * t281) * qJD(4) + t388) * t217 + (-t237 * t2 + t238 * t1 + t217 * (-t10 * t237 + t11 * t238) + (t53 * t217 - t219 * t281) * qJD(4) + (-t238 * t12 - t237 * t13 - t217 * t282) * qJD(1)) * t219;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
