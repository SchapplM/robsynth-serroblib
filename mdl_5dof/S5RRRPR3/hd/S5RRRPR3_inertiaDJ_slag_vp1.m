% Calculate time derivative of joint inertia matrix for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:50
% EndTime: 2020-01-03 12:09:06
% DurationCPUTime: 9.81s
% Computational Cost: add. (10947->521), mult. (7952->704), div. (0->0), fcn. (5984->10), ass. (0->304)
t237 = qJ(3) + pkin(9);
t227 = cos(t237);
t226 = sin(t237);
t371 = Icges(5,4) * t226;
t171 = Icges(5,2) * t227 + t371;
t242 = cos(qJ(3));
t240 = sin(qJ(3));
t373 = Icges(4,4) * t240;
t192 = Icges(4,2) * t242 + t373;
t418 = -t171 * t226 - t192 * t240;
t370 = Icges(5,4) * t227;
t172 = Icges(5,1) * t226 + t370;
t372 = Icges(4,4) * t242;
t193 = Icges(4,1) * t240 + t372;
t415 = -t172 * t227 - t193 * t242 - t418;
t170 = Icges(5,5) * t226 + Icges(5,6) * t227;
t191 = Icges(4,5) * t240 + Icges(4,6) * t242;
t417 = -t170 - t191;
t236 = qJD(1) + qJD(2);
t285 = Icges(5,5) * t227 - Icges(5,6) * t226;
t286 = Icges(4,5) * t242 - Icges(4,6) * t240;
t416 = -t415 * t236 + (-t285 - t286) * qJD(3);
t228 = qJ(5) + t237;
t217 = sin(t228);
t377 = rSges(6,2) * t217;
t218 = cos(t228);
t380 = rSges(6,1) * t218;
t414 = -t377 + t380;
t378 = rSges(5,2) * t226;
t382 = rSges(5,1) * t227;
t413 = -t378 + t382;
t288 = -Icges(5,2) * t226 + t370;
t155 = t288 * qJD(3);
t291 = Icges(5,1) * t227 - t371;
t156 = t291 * qJD(3);
t289 = -Icges(4,2) * t240 + t372;
t178 = t289 * qJD(3);
t292 = Icges(4,1) * t242 - t373;
t179 = t292 * qJD(3);
t412 = t155 * t226 - t156 * t227 + t178 * t240 - t179 * t242 + (-0.2e1 * t170 - t191) * t236 + (0.2e1 * t171 * t227 + 0.2e1 * t172 * t226 + t192 * t242 + t193 * t240) * qJD(3);
t235 = qJD(3) + qJD(5);
t141 = t414 * t235;
t238 = qJ(1) + qJ(2);
t230 = cos(t238);
t323 = qJD(3) * t242;
t312 = t230 * t323;
t190 = pkin(3) * t312;
t325 = qJD(3) * t227;
t229 = sin(t238);
t376 = rSges(6,2) * t218;
t381 = rSges(6,1) * t217;
t161 = t376 + t381;
t385 = pkin(4) * t226;
t386 = pkin(3) * t240;
t295 = -t385 - t386;
t270 = -t161 + t295;
t93 = t270 * t229;
t374 = t236 * t93;
t36 = t190 + (pkin(4) * t325 + t141) * t230 + t374;
t411 = t36 - t374;
t133 = -Icges(4,6) * t229 - t230 * t289;
t135 = -Icges(4,5) * t229 - t230 * t292;
t276 = t133 * t240 - t135 * t242;
t410 = t229 * t276;
t120 = -Icges(5,6) * t229 - t230 * t288;
t122 = -Icges(5,5) * t229 - t230 * t291;
t280 = t120 * t226 - t122 * t227;
t409 = t229 * t280;
t368 = Icges(6,4) * t218;
t287 = -Icges(6,2) * t217 + t368;
t108 = -Icges(6,6) * t229 - t230 * t287;
t369 = Icges(6,4) * t217;
t290 = Icges(6,1) * t218 - t369;
t110 = -Icges(6,5) * t229 - t230 * t290;
t282 = t108 * t217 - t110 * t218;
t408 = t229 * t282;
t239 = -qJ(4) - pkin(7);
t234 = -pkin(8) + t239;
t327 = t234 - t239;
t407 = t229 * t327;
t175 = rSges(3,1) * t229 + rSges(3,2) * t230;
t160 = Icges(6,1) * t217 + t368;
t351 = t160 * t235;
t406 = -Icges(6,5) * t236 + t351;
t159 = Icges(6,2) * t218 + t369;
t352 = t159 * t235;
t405 = -Icges(6,6) * t236 + t352;
t158 = Icges(6,5) * t217 + Icges(6,6) * t218;
t404 = -Icges(6,3) * t236 + t158 * t235;
t403 = -Icges(5,3) * t236 + qJD(3) * t170;
t400 = -Icges(4,3) * t236 + qJD(3) * t191;
t399 = -Icges(4,6) * t236 + qJD(3) * t192;
t398 = -Icges(4,5) * t236 + qJD(3) * t193;
t139 = t287 * t235;
t140 = t290 * t235;
t395 = t217 * (t139 + t351) + t218 * (-t140 + t352) - t158 * t236;
t394 = 2 * m(3);
t393 = 2 * m(4);
t392 = 2 * m(5);
t391 = 2 * m(6);
t390 = -t229 / 0.2e1;
t389 = -t230 / 0.2e1;
t379 = rSges(4,2) * t240;
t383 = rSges(4,1) * t242;
t185 = (-t379 + t383) * qJD(3);
t388 = m(4) * t185;
t204 = rSges(4,1) * t240 + rSges(4,2) * t242;
t387 = m(4) * t204;
t384 = pkin(4) * t227;
t232 = t242 * pkin(3);
t222 = t232 + pkin(2);
t375 = pkin(1) * qJD(1);
t284 = Icges(6,5) * t218 - Icges(6,6) * t217;
t105 = -Icges(6,3) * t230 + t229 * t284;
t358 = t105 * t236;
t106 = -Icges(6,3) * t229 - t230 * t284;
t357 = t106 * t236;
t174 = rSges(5,1) * t226 + rSges(5,2) * t227;
t309 = -t174 - t386;
t125 = t309 * t229;
t354 = t125 * t236;
t348 = t217 * t235;
t347 = t218 * t235;
t346 = t229 * t235;
t345 = t229 * t236;
t344 = t229 * t242;
t343 = t230 * t235;
t342 = t230 * t236;
t209 = t230 * t239;
t341 = t230 * t240;
t340 = t236 * t239;
t173 = t295 * qJD(3);
t184 = t222 + t384;
t339 = t173 * t229 + t184 * t342;
t338 = t229 * t184 + t230 * t234;
t318 = t236 * t377;
t337 = -rSges(6,3) * t342 - t229 * t318;
t183 = t230 * t380;
t336 = rSges(6,3) * t345 + t183 * t236;
t319 = t236 * t378;
t335 = rSges(5,3) * t342 + t229 * t319;
t187 = t230 * t382;
t334 = rSges(5,3) * t345 + t187 * t236;
t320 = t229 * t379;
t333 = rSges(4,3) * t342 + t236 * t320;
t206 = t230 * t383;
t332 = rSges(4,3) * t345 + t206 * t236;
t331 = t229 * t222 + t209;
t330 = -pkin(2) * t342 - pkin(7) * t345;
t329 = pkin(2) * t230 + pkin(7) * t229;
t328 = t229 ^ 2 + t230 ^ 2;
t326 = qJD(3) * t226;
t324 = qJD(3) * t240;
t220 = t229 * pkin(2);
t111 = pkin(7) * t230 - t220 + t331;
t189 = t230 * t222;
t298 = t229 * t239 - t189;
t112 = t298 + t329;
t168 = t222 * t342;
t322 = t229 * (-qJD(4) * t230 + t168 + (-pkin(3) * t324 - t340) * t229 + t330) + t112 * t345 + t111 * t342;
t113 = -rSges(6,3) * t230 + t229 * t414;
t114 = -rSges(6,3) * t229 + t230 * t377 - t183;
t321 = t113 * t342 + t229 * (-t346 * t381 + (-t217 * t342 - t218 * t346) * rSges(6,2) + t336) + t114 * t345;
t241 = sin(qJ(1));
t317 = t241 * t375;
t107 = -Icges(6,6) * t230 + t229 * t287;
t109 = -Icges(6,5) * t230 + t229 * t290;
t283 = t107 * t217 - t109 * t218;
t24 = -t105 * t230 - t229 * t283;
t25 = -t106 * t230 - t408;
t257 = t283 * t230;
t26 = -t105 * t229 + t257;
t27 = -t106 * t229 + t230 * t282;
t261 = t287 * t236;
t57 = t229 * t261 + t230 * t405;
t301 = t110 * t235 + t57;
t264 = t290 * t236;
t59 = t229 * t264 + t230 * t406;
t303 = -t108 * t235 + t59;
t258 = t284 * t236;
t55 = t229 * t258 + t230 * t404;
t56 = -t229 * t404 + t230 * t258;
t58 = -t229 * t405 + t230 * t261;
t60 = -t229 * t406 + t230 * t264;
t316 = -t229 * ((t229 * t55 + (t26 + t408) * t236) * t229 + (-t27 * t236 + (-t107 * t347 - t109 * t348 - t217 * t58 + t218 * t60 + t358) * t230 + (t357 + t56 + (-t109 * t236 + t303) * t218 + (t107 * t236 - t301) * t217) * t229) * t230) + (-t229 * t25 - t230 * t24) * t345;
t302 = t109 * t235 + t58;
t304 = -t107 * t235 + t60;
t2 = (t230 * t56 + (-t25 + t257) * t236) * t230 + (t24 * t236 + (t108 * t347 + t110 * t348 + t217 * t57 - t218 * t59 - t357) * t229 + (-t358 + t55 + (-t110 * t236 - t304) * t218 + (t108 * t236 + t302) * t217) * t230) * t229;
t5 = -t229 * t27 - t230 * t26;
t315 = -t236 * t5 - t2;
t314 = t229 * t324;
t313 = t230 * t324;
t311 = t345 / 0.2e1;
t310 = -t342 / 0.2e1;
t149 = t230 * t173;
t210 = qJD(4) * t229;
t248 = -t161 * t235 - t234 * t236;
t22 = t149 + t210 + (-t184 - t380) * t345 + t248 * t230 - t337;
t20 = t22 - t317;
t231 = t241 * pkin(1);
t79 = t113 + t338;
t75 = t231 + t79;
t308 = -t236 * t75 - t20;
t243 = cos(qJ(1));
t223 = t243 * t375;
t23 = (-qJD(4) - t318) * t230 + t248 * t229 + t336 + t339;
t21 = t223 + t23;
t233 = t243 * pkin(1);
t153 = t230 * t184;
t80 = -t229 * t234 - t114 + t153;
t76 = t233 + t80;
t307 = t236 * t76 - t21;
t306 = -t236 * t79 - t22;
t305 = t236 * t80 - t23;
t176 = rSges(3,1) * t230 - rSges(3,2) * t229;
t297 = pkin(3) * t313;
t147 = rSges(3,1) * t342 - rSges(3,2) * t345;
t296 = rSges(4,1) * t344 - t320;
t262 = t288 * t229;
t119 = -Icges(5,6) * t230 + t262;
t265 = t291 * t229;
t121 = -Icges(5,5) * t230 + t265;
t281 = t119 * t226 - t121 * t227;
t132 = -Icges(4,6) * t230 + t229 * t289;
t134 = -Icges(4,5) * t230 + t229 * t292;
t279 = t132 * t242 + t134 * t240;
t278 = t132 * t240 - t134 * t242;
t277 = t133 * t242 + t135 * t240;
t275 = t159 * t217 - t160 * t218;
t137 = rSges(4,2) * t341 - rSges(4,3) * t229 - t206;
t124 = -rSges(5,3) * t229 + t230 * t378 - t187;
t251 = -t235 * t284 - t236 * t275;
t269 = (t217 * t303 + t218 * t301 + t229 * t251 + t230 * t395) * t390 + (t217 * t304 + t218 * t302 - t229 * t395 + t230 * t251) * t389 + (t107 * t218 + t109 * t217 - t158 * t230 - t229 * t275) * t311 + (t108 * t218 + t110 * t217 - t158 * t229 + t230 * t275) * t310;
t146 = t175 * t236;
t268 = qJD(3) * t204;
t267 = qJD(3) * t174;
t266 = t292 * t236;
t263 = t289 * t236;
t260 = t286 * t229;
t259 = t285 * t229;
t256 = t281 * t230;
t255 = t278 * t230;
t123 = -rSges(5,3) * t230 + t229 * t413;
t254 = -t229 * t323 - t236 * t341;
t103 = -t137 + t329;
t89 = t123 + t331;
t90 = -t124 - t298;
t102 = t220 + (-rSges(4,3) - pkin(7)) * t230 + t296;
t247 = qJD(3) * t309 - t340;
t246 = t139 * t218 + t140 * t217 + t155 * t227 + t156 * t226 - t159 * t348 + t160 * t347 + t172 * t325 + t178 * t242 + t179 * t240 + t193 * t323;
t54 = -rSges(4,1) * t314 + rSges(4,2) * t254 - t330 + t332;
t245 = qJD(3) * t418 + t246;
t244 = t269 + (t240 * (t229 * t266 + t230 * t398) + t242 * (t229 * t263 + t230 * t399) + (t226 * t265 + t227 * t262) * t236 + t412 * t230 + t416 * t229 + (-t276 - t280) * qJD(3)) * t390 + (t240 * (-t229 * t398 + t230 * t266) + t242 * (-t229 * t399 + t230 * t263) + (t226 * t291 + t227 * t288) * t342 + t416 * t230 - t412 * t229 + (-t278 - t281) * qJD(3)) * t389 + (t119 * t227 + t121 * t226 - t229 * t415 + t230 * t417 + t279) * t311 + (t120 * t227 + t122 * t226 + t229 * t417 + t230 * t415 + t277) * t310;
t202 = pkin(7) * t342;
t53 = -rSges(4,2) * t312 - pkin(2) * t345 + t202 + (-t236 * t344 - t313) * rSges(4,1) + t333;
t34 = t210 + (-t222 - t382) * t345 + t247 * t230 + t335;
t35 = t168 + (-qJD(4) - t319) * t230 + t247 * t229 + t334;
t208 = pkin(3) * t341;
t157 = t413 * qJD(3);
t151 = t176 + t233;
t150 = t231 + t175;
t136 = -rSges(4,3) * t230 + t296;
t131 = -Icges(4,3) * t229 - t230 * t286;
t130 = -Icges(4,3) * t230 + t260;
t129 = t147 + t223;
t128 = -t146 - t317;
t126 = t174 * t230 + t208;
t118 = -Icges(5,3) * t229 - t230 * t285;
t117 = -Icges(5,3) * t230 + t259;
t104 = t229 * t113;
t101 = t229 * t111;
t96 = t103 + t233;
t95 = t231 + t102;
t94 = t208 + (t161 + t385) * t230;
t92 = -t153 + t189 + t407;
t91 = -t331 + t338;
t88 = t233 + t90;
t87 = t231 + t89;
t82 = -t229 * t400 + t286 * t342;
t81 = t230 * t400 + t236 * t260;
t74 = t297 + t202 - t210 + (t209 + (-pkin(2) + t222) * t229) * t236;
t69 = -t229 * t403 + t285 * t342;
t68 = t230 * t403 + t236 * t259;
t67 = pkin(3) * t254 - t157 * t229 - t174 * t342;
t66 = t157 * t230 + t190 + t354;
t62 = t343 * t376 + (t217 * t343 + t218 * t345) * rSges(6,1) + t337;
t61 = -t114 * t230 + t104;
t49 = t223 + t54;
t48 = t53 - t317;
t41 = -t131 * t229 + t230 * t276;
t40 = -t130 * t229 + t255;
t39 = -t131 * t230 - t410;
t38 = -t130 * t230 - t229 * t278;
t37 = t270 * t342 + (-t141 + (-t232 - t384) * qJD(3)) * t229;
t33 = t223 + t35;
t32 = t34 - t317;
t31 = -t118 * t229 + t230 * t280;
t30 = -t117 * t229 + t256;
t29 = -t118 * t230 - t409;
t28 = -t117 * t230 - t229 * t281;
t11 = -t230 * t62 + t321;
t10 = t229 * t91 + t101 + t104 + (-t112 - t114 - t92) * t230;
t3 = (pkin(3) * t314 - t168 + (t92 - t407) * t236 + t339) * t229 + (t297 + t149 - t62 - t74 + (t91 - t327 * t230 + (-t184 + t222) * t229) * t236) * t230 + t321 + t322;
t1 = [(t128 * t151 + t129 * t150) * t394 + (t48 * t96 + t49 * t95) * t393 + (t32 * t88 + t33 * t87) * t392 + (t20 * t76 + t21 * t75) * t391 + t246 - t171 * t326 - t192 * t324; m(3) * (t128 * t176 + t129 * t175 - t146 * t151 + t147 * t150) + m(4) * (t102 * t49 + t103 * t48 + t53 * t96 + t54 * t95) + m(5) * (t32 * t90 + t33 * t89 + t34 * t88 + t35 * t87) + m(6) * (t20 * t80 + t21 * t79 + t22 * t76 + t23 * t75) + t245; (t22 * t80 + t23 * t79) * t391 + (t34 * t90 + t35 * t89) * t392 + (t102 * t54 + t103 * t53) * t393 + (-t146 * t176 + t147 * t175) * t394 + t245; m(6) * (t20 * t93 + t21 * t94 + t36 * t75 + t37 * t76) + m(5) * (t125 * t32 + t126 * t33 + t66 * t87 + t67 * t88) + t244 + ((-t236 * t96 + t49) * t230 + (-t236 * t95 - t48) * t229) * t387 + (-t229 * t96 + t230 * t95) * t388; ((-t103 * t236 + t54) * t230 + (-t102 * t236 - t53) * t229) * t387 + (t102 * t230 - t103 * t229) * t388 + t244 + m(6) * (t22 * t93 + t23 * t94 + t36 * t79 + t37 * t80) + m(5) * (t125 * t34 + t126 * t35 + t66 * t89 + t67 * t90); (t10 * t3 + t36 * t94 + t37 * t93) * t391 - t230 * t2 + ((t123 * t229 + t101 + (-t112 - t124) * t230) * ((t236 * t124 - t229 * t267 + t334) * t229 + (-t74 - t230 * t267 + (t123 + (-t378 - t382) * t229) * t236 + t335) * t230 + t322) + t125 * t67 + t126 * t66) * t392 + ((t136 * t229 - t137 * t230) * ((t236 * t136 - t230 * t268 + t333) * t230 + (-t229 * t268 + (t137 + (-t379 - t383) * t230) * t236 + t332) * t229) + t328 * t204 * t185) * t393 - t230 * ((t230 * t69 + (-t29 + t256) * t236) * t230 + (t28 * t236 + (t120 * t325 + t122 * t326) * t229 + (t68 + (qJD(3) * t119 - t122 * t236) * t227 + (qJD(3) * t121 + t120 * t236) * t226) * t230) * t229) - t230 * ((t230 * t82 + (-t39 + t255) * t236) * t230 + (t38 * t236 + (t133 * t323 + t135 * t324) * t229 + (t279 * qJD(3) + t236 * t276 + t81) * t230) * t229) - t229 * ((t229 * t81 + (t40 + t410) * t236) * t229 + (-t41 * t236 + (-t132 * t323 - t134 * t324) * t230 + (-t277 * qJD(3) + t236 * t278 + t82) * t229) * t230) - t229 * ((t229 * t68 + (t30 + t409) * t236) * t229 + (-t31 * t236 + (-t119 * t325 - t121 * t326) * t230 + (t69 + (-qJD(3) * t120 - t121 * t236) * t227 + (-qJD(3) * t122 + t119 * t236) * t226) * t229) * t230) + t316 + ((-t28 - t38) * t230 + (-t29 - t39) * t229) * t345 + (-t5 + (t30 + t40) * t230 + (t31 + t41) * t229) * t342; m(5) * ((-t236 * t87 - t32) * t230 + (t236 * t88 - t33) * t229) + m(6) * (t229 * t307 + t230 * t308); m(6) * (t229 * t305 + t230 * t306) + m(5) * ((-t236 * t89 - t34) * t230 + (t236 * t90 - t35) * t229); m(6) * ((-t236 * t94 - t37) * t230 - t411 * t229) + m(5) * ((-t126 * t236 - t67) * t230 + (-t66 + t354) * t229); 0; m(6) * ((-t229 * t76 + t230 * t75) * t141 + (t229 * t308 - t230 * t307) * t161) + t269; m(6) * ((-t229 * t80 + t230 * t79) * t141 + (t229 * t306 - t230 * t305) * t161) + t269; m(6) * (-t141 * t229 * t93 + t11 * t10 + t61 * t3 + (-t229 * t37 - t345 * t94) * t161) + (m(6) * (t141 * t94 + t161 * t411) + t315) * t230 + t316; 0; (t141 * t161 * t328 + t11 * t61) * t391 + t315 * t230 + t316;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
