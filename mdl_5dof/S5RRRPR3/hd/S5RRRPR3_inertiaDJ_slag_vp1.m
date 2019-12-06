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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:42:19
% EndTime: 2019-12-05 18:42:34
% DurationCPUTime: 8.35s
% Computational Cost: add. (10947->510), mult. (7952->688), div. (0->0), fcn. (5984->10), ass. (0->304)
t235 = qJ(3) + pkin(9);
t227 = cos(t235);
t226 = sin(t235);
t378 = Icges(5,4) * t226;
t167 = Icges(5,2) * t227 + t378;
t240 = cos(qJ(3));
t238 = sin(qJ(3));
t380 = Icges(4,4) * t238;
t202 = Icges(4,2) * t240 + t380;
t431 = -t167 * t226 - t202 * t238;
t377 = Icges(5,4) * t227;
t168 = Icges(5,1) * t226 + t377;
t379 = Icges(4,4) * t240;
t203 = Icges(4,1) * t238 + t379;
t429 = -t168 * t227 - t203 * t240 - t431;
t166 = Icges(5,5) * t226 + Icges(5,6) * t227;
t201 = Icges(4,5) * t238 + Icges(4,6) * t240;
t430 = t166 + t201;
t236 = qJ(1) + qJ(2);
t229 = sin(t236);
t234 = qJD(1) + qJD(2);
t335 = qJD(3) * t240;
t230 = cos(t236);
t350 = t230 * t238;
t428 = t229 * t335 + t234 * t350;
t337 = qJD(3) * t227;
t352 = t230 * t234;
t427 = t226 * t352 + t229 * t337;
t228 = qJ(5) + t235;
t220 = cos(t228);
t233 = qJD(3) + qJD(5);
t355 = t220 * t233;
t219 = sin(t228);
t356 = t219 * t233;
t426 = rSges(6,1) * t356 + rSges(6,2) * t355;
t282 = Icges(5,5) * t227 - Icges(5,6) * t226;
t283 = Icges(4,5) * t240 - Icges(4,6) * t238;
t285 = -Icges(5,2) * t226 + t377;
t288 = Icges(5,1) * t227 - t378;
t425 = (t282 + t283) * qJD(3) + (-t226 * t288 - t227 * t285 + t429) * t234;
t147 = t285 * qJD(3);
t148 = t288 * qJD(3);
t286 = -Icges(4,2) * t238 + t379;
t174 = t286 * qJD(3);
t289 = Icges(4,1) * t240 - t380;
t175 = t289 * qJD(3);
t424 = t147 * t226 - t148 * t227 + t174 * t238 - t175 * t240 + (-0.2e1 * t166 - t201) * t234 + (0.2e1 * t167 * t227 + 0.2e1 * t168 * t226 + t202 * t240 + t203 * t238) * qJD(3);
t261 = t286 * t229;
t128 = Icges(4,6) * t230 - t261;
t130 = Icges(4,5) * t230 - t229 * t289;
t275 = t128 * t238 - t130 * t240;
t423 = t230 * t275;
t115 = Icges(5,6) * t230 - t229 * t285;
t117 = Icges(5,5) * t230 - t229 * t288;
t278 = t115 * t226 - t117 * t227;
t422 = t230 * t278;
t375 = Icges(6,4) * t220;
t284 = -Icges(6,2) * t219 + t375;
t259 = t284 * t229;
t103 = Icges(6,6) * t230 - t259;
t376 = Icges(6,4) * t219;
t287 = Icges(6,1) * t220 - t376;
t105 = Icges(6,5) * t230 - t229 * t287;
t280 = t103 * t219 - t105 * t220;
t421 = t230 * t280;
t385 = rSges(6,2) * t219;
t345 = t230 * rSges(6,3) + t229 * t385;
t152 = Icges(6,1) * t219 + t375;
t359 = t152 * t233;
t420 = -Icges(6,5) * t234 + t359;
t151 = Icges(6,2) * t220 + t376;
t360 = t151 * t233;
t419 = -Icges(6,6) * t234 + t360;
t150 = Icges(6,5) * t219 + Icges(6,6) * t220;
t418 = -Icges(6,3) * t234 + t150 * t233;
t231 = t240 * pkin(3);
t222 = t231 + pkin(2);
t392 = pkin(2) - t222;
t417 = t229 * pkin(7) + t230 * t392;
t416 = -Icges(5,3) * t234 + qJD(3) * t166;
t413 = -Icges(4,3) * t234 + qJD(3) * t201;
t412 = -Icges(4,6) * t234 + qJD(3) * t202;
t411 = -Icges(4,5) * t234 + qJD(3) * t203;
t135 = t284 * t233;
t136 = t287 * t233;
t408 = t219 * (t135 + t359) + t220 * (-t136 + t360) - t150 * t234;
t407 = 2 * m(3);
t406 = 2 * m(4);
t405 = 2 * m(5);
t404 = 2 * m(6);
t403 = t229 / 0.2e1;
t402 = t230 / 0.2e1;
t401 = -rSges(4,3) - pkin(7);
t388 = rSges(4,2) * t238;
t391 = rSges(4,1) * t240;
t187 = (-t388 + t391) * qJD(3);
t400 = m(4) * t187;
t208 = rSges(4,1) * t238 + rSges(4,2) * t240;
t399 = m(4) * t208;
t239 = sin(qJ(1));
t398 = pkin(1) * t239;
t241 = cos(qJ(1));
t397 = pkin(1) * t241;
t396 = pkin(3) * t238;
t395 = pkin(4) * t226;
t394 = pkin(4) * t227;
t237 = -qJ(4) - pkin(7);
t390 = rSges(5,1) * t227;
t389 = rSges(6,1) * t220;
t387 = rSges(5,2) * t226;
t386 = rSges(5,2) * t230;
t384 = pkin(1) * qJD(1);
t383 = t229 * rSges(4,3);
t382 = t229 * rSges(5,3);
t381 = t229 * rSges(6,3);
t218 = t230 * rSges(4,3);
t217 = t230 * rSges(5,3);
t354 = t229 * t234;
t353 = t229 * t238;
t351 = t230 * t237;
t232 = -pkin(8) + t237;
t349 = t232 * t234;
t348 = t234 * t237;
t221 = t230 * pkin(7);
t107 = t229 * t392 - t221 - t351;
t333 = t229 * t390;
t341 = t229 * t387 + t217;
t347 = -t107 + t333 - t341;
t186 = t222 + t394;
t346 = t186 * t354 + t230 * t349;
t344 = t186 - t222;
t343 = t428 * pkin(3);
t336 = qJD(3) * t238;
t315 = t229 * t336;
t342 = pkin(3) * t315 + t229 * t348;
t142 = rSges(3,1) * t354 + rSges(3,2) * t352;
t209 = rSges(4,2) * t353;
t340 = t209 + t218;
t339 = t229 ^ 2 + t230 ^ 2;
t338 = qJD(3) * t226;
t212 = pkin(3) * t353;
t334 = t229 * t391;
t331 = t229 * t389;
t329 = t241 * t384;
t183 = t230 * t385;
t109 = -t331 + t345;
t297 = t229 * t344;
t328 = -t107 - t109 - (-t232 + t237) * t230 + t297;
t327 = rSges(5,1) * t338;
t324 = t234 * t183 + t426 * t229;
t323 = -t426 * t230 - t234 * t331;
t322 = -t230 * t327 - t234 * t333 - t337 * t386;
t321 = -t427 * rSges(5,2) - t229 * t327;
t313 = t230 * t336;
t320 = pkin(3) * t313 + t222 * t354 + t230 * t348;
t319 = -t230 * rSges(4,2) * t335 - rSges(4,1) * t313 - t234 * t334;
t318 = rSges(4,1) * t315 + t428 * rSges(4,2);
t214 = qJD(4) * t230;
t317 = t214 + t342;
t312 = -t354 / 0.2e1;
t311 = t352 / 0.2e1;
t310 = -pkin(2) - t391;
t293 = -t395 - t396;
t169 = t293 * qJD(3);
t22 = (-rSges(6,3) * t234 - t169) * t230 + (-t234 * t385 - qJD(4)) * t229 - t323 + t346;
t223 = t239 * t384;
t20 = t223 + t22;
t302 = -t186 - t389;
t79 = t229 * t302 - t230 * t232 + t345;
t75 = t79 - t398;
t309 = t234 * t75 + t20;
t249 = t230 * t302 - t381;
t294 = (-t169 + t349) * t229;
t23 = t234 * t249 + t214 + t294 + t324;
t21 = t23 - t329;
t211 = t229 * t232;
t80 = t183 + t211 + t249;
t76 = t80 - t397;
t308 = t234 * t76 - t21;
t307 = t234 * t79 + t22;
t306 = t234 * t80 - t23;
t137 = (-t385 + t389) * t233;
t188 = t234 * t212;
t153 = rSges(6,1) * t219 + rSges(6,2) * t220;
t292 = (t153 + t395) * t229;
t36 = t188 + t234 * t292 + (-t137 + (-t231 - t394) * qJD(3)) * t230;
t93 = t212 + t292;
t305 = t234 * t93 - t36;
t37 = t427 * pkin(4) + t137 * t229 + t153 * t352 + t343;
t94 = (-t153 + t293) * t230;
t304 = t234 * t94 + t37;
t172 = -t230 * rSges(3,1) + t229 * rSges(3,2);
t303 = -t222 - t390;
t262 = t287 * t234;
t60 = t229 * t420 - t230 * t262;
t301 = -t103 * t233 + t60;
t104 = Icges(6,6) * t229 + t230 * t284;
t59 = -t229 * t262 - t230 * t420;
t300 = -t104 * t233 + t59;
t58 = t229 * t419 - t284 * t352;
t299 = t105 * t233 + t58;
t106 = Icges(6,5) * t229 + t230 * t287;
t57 = -t230 * t419 - t234 * t259;
t298 = t106 * t233 + t57;
t143 = -rSges(3,1) * t352 + rSges(3,2) * t354;
t171 = -rSges(3,1) * t229 - rSges(3,2) * t230;
t281 = Icges(6,5) * t220 - Icges(6,6) * t219;
t279 = t104 * t219 - t106 * t220;
t116 = Icges(5,6) * t229 + t230 * t285;
t118 = Icges(5,5) * t229 + t230 * t288;
t277 = t116 * t226 - t118 * t227;
t276 = t128 * t240 + t130 * t238;
t129 = Icges(4,6) * t229 + t230 * t286;
t131 = Icges(4,5) * t229 + t230 * t289;
t274 = t129 * t240 + t131 * t238;
t273 = t129 * t238 - t131 * t240;
t272 = t151 * t219 - t152 * t220;
t101 = Icges(6,3) * t230 - t229 * t281;
t24 = t101 * t230 + t229 * t280;
t102 = Icges(6,3) * t229 + t230 * t281;
t255 = t279 * t229;
t25 = t102 * t230 + t255;
t26 = t101 * t229 - t421;
t27 = t102 * t229 - t230 * t279;
t256 = t281 * t234;
t55 = -t229 * t256 - t230 * t418;
t56 = t229 * t418 - t230 * t256;
t267 = -t354 * (t229 * t25 + t230 * t24) + t229 * ((t229 * t55 + (-t26 + t255) * t234) * t229 + (t27 * t234 + (-t103 * t355 - t105 * t356 - t219 * t58 + t220 * t60) * t230 + (t56 + (-t105 * t234 + t300) * t220 + (t103 * t234 - t298) * t219) * t229) * t230) + t230 * ((t230 * t56 + (t25 + t421) * t234) * t230 + (-t24 * t234 + (t104 * t355 + t106 * t356 + t219 * t57 - t220 * t59) * t229 + (t55 + (-t106 * t234 - t301) * t220 + (t104 * t234 + t299) * t219) * t230) * t229) + (t229 * t27 + t230 * t26) * t352;
t266 = -t230 * t389 - t381;
t248 = t281 * t233 + t234 * t272;
t265 = (t219 * t300 + t220 * t298 + t248 * t229 - t230 * t408) * t403 + (t219 * t301 + t220 * t299 + t229 * t408 + t248 * t230) * t402 + (t103 * t220 + t105 * t219 + t150 * t230 + t229 * t272) * t312 + (t104 * t220 + t106 * t219 + t150 * t229 - t230 * t272) * t311;
t264 = t289 * t234;
t258 = t283 * t229;
t257 = t282 * t229;
t254 = t277 * t229;
t253 = t273 * t229;
t250 = t230 * t303 - t382;
t98 = t229 * t310 + t221 + t340;
t245 = t229 * t401 + t230 * t310;
t210 = rSges(4,2) * t350;
t99 = t210 + t245;
t193 = t226 * t386;
t213 = t229 * t237;
t90 = t193 + t213 + t250;
t89 = t229 * t303 + t341 - t351;
t244 = t220 * t135 + t219 * t136 + t227 * t147 + t226 * t148 - t151 * t356 + t152 * t355 + t168 * t337 + t240 * t174 + t238 * t175 + t203 * t335;
t207 = pkin(2) * t354;
t53 = t207 + (t230 * t401 - t209) * t234 - t319;
t54 = t234 * t245 + t318;
t35 = t234 * t250 + t317 - t321;
t34 = -rSges(5,3) * t352 + (-t234 * t387 - qJD(4)) * t229 + t320 - t322;
t243 = qJD(3) * t431 + t244;
t242 = t265 + (t238 * (-t229 * t264 - t230 * t411) + t240 * (-t230 * t412 - t234 * t261) - t424 * t230 + t425 * t229 + (-t273 - t277) * qJD(3)) * t403 + (t238 * (t229 * t411 - t230 * t264) + t240 * (t229 * t412 - t286 * t352) + t425 * t230 + t424 * t229 + (-t275 - t278) * qJD(3)) * t402 + (t115 * t227 + t117 * t226 + t429 * t229 + t430 * t230 + t276) * t312 + (t116 * t227 + t118 * t226 + t430 * t229 - t429 * t230 + t274) * t311;
t170 = rSges(5,1) * t226 + rSges(5,2) * t227;
t149 = (-t387 + t390) * qJD(3);
t145 = t172 - t397;
t144 = t171 - t398;
t133 = t230 * t391 - t210 + t383;
t132 = -t334 + t340;
t127 = Icges(4,3) * t229 + t230 * t283;
t126 = Icges(4,3) * t230 - t258;
t125 = t143 - t329;
t124 = t223 + t142;
t122 = (-t170 - t396) * t230;
t121 = t170 * t229 + t212;
t120 = t230 * t390 - t193 + t382;
t114 = Icges(5,3) * t229 + t230 * t282;
t113 = Icges(5,3) * t230 - t257;
t110 = -t183 - t266;
t108 = -t213 - t417;
t100 = t230 * t110;
t97 = t230 * t108;
t96 = t99 - t397;
t95 = t98 - t398;
t92 = t230 * t344 - t211 + t213;
t88 = t90 - t397;
t87 = t89 - t398;
t82 = t229 * t413 - t283 * t352;
t81 = -t230 * t413 - t234 * t258;
t74 = t234 * t417 + t317;
t69 = t229 * t416 - t282 * t352;
t68 = -t230 * t416 - t234 * t257;
t67 = t149 * t229 + t170 * t352 + t343;
t66 = t170 * t354 + t188 + (-pkin(3) * t335 - t149) * t230;
t65 = t230 * (-pkin(7) * t352 + qJD(4) * t229 + t207 - t320);
t62 = t234 * t266 + t324;
t61 = -t109 * t229 + t100;
t52 = t230 * (t234 * t345 + t323);
t49 = t54 - t329;
t48 = t223 + t53;
t41 = t127 * t229 - t273 * t230;
t40 = t126 * t229 - t423;
t39 = t127 * t230 + t253;
t38 = t126 * t230 + t275 * t229;
t33 = t35 - t329;
t32 = t223 + t34;
t31 = t114 * t229 - t230 * t277;
t30 = t113 * t229 - t422;
t29 = t114 * t230 + t254;
t28 = t113 * t230 + t229 * t278;
t11 = -t109 * t352 + t52 + (-t110 * t234 - t62) * t229;
t10 = t229 * t328 + t230 * t92 + t100 + t97;
t3 = t65 + t230 * (t169 * t230 + t320 - t346) + t52 + (-t294 - t62 - t74 + t342) * t229 + ((-t108 - t110 - t92) * t229 + (t297 + t328) * t230) * t234;
t1 = [(t124 * t145 + t125 * t144) * t407 + (t48 * t96 + t49 * t95) * t406 + (t32 * t88 + t33 * t87) * t405 + (t20 * t76 + t21 * t75) * t404 + t244 - t167 * t338 - t202 * t336; m(3) * (t124 * t172 + t125 * t171 + t142 * t145 + t143 * t144) + m(4) * (t48 * t99 + t49 * t98 + t53 * t96 + t54 * t95) + m(5) * (t32 * t90 + t33 * t89 + t34 * t88 + t35 * t87) + m(6) * (t20 * t80 + t21 * t79 + t22 * t76 + t23 * t75) + t243; t243 + (t22 * t80 + t23 * t79) * t404 + (t34 * t90 + t35 * t89) * t405 + (t53 * t99 + t54 * t98) * t406 + (t142 * t172 + t143 * t171) * t407; ((t234 * t96 - t49) * t230 + (t234 * t95 + t48) * t229) * t399 + m(6) * (t20 * t93 + t21 * t94 + t36 * t75 + t37 * t76) + m(5) * (t121 * t32 + t122 * t33 + t66 * t87 + t67 * t88) + t242 + (t229 * t96 - t230 * t95) * t400; ((t234 * t99 - t54) * t230 + (t234 * t98 + t53) * t229) * t399 + (t229 * t99 - t230 * t98) * t400 + m(6) * (t22 * t93 + t23 * t94 + t36 * t79 + t37 * t80) + m(5) * (t121 * t34 + t122 * t35 + t66 * t89 + t67 * t90) + t242; (t10 * t3 + t36 * t94 + t37 * t93) * t404 + ((t120 * t230 + t229 * t347 + t97) * (t65 + t230 * t322 + (-t74 + t321) * t229 + ((t347 + t217) * t230 + (t382 - t108 - t120 + (t387 + t390) * t230) * t229) * t234) + t121 * t67 + t122 * t66) * t405 + t229 * ((t229 * t68 + (-t30 + t254) * t234) * t229 + (t31 * t234 + (-t115 * t337 - t117 * t338) * t230 + (t69 + (-qJD(3) * t116 - t117 * t234) * t227 + (-qJD(3) * t118 + t115 * t234) * t226) * t229) * t230) + t229 * ((t229 * t81 + (-t40 + t253) * t234) * t229 + (t41 * t234 + (-t128 * t335 - t130 * t336) * t230 + (-t274 * qJD(3) + t234 * t275 + t82) * t229) * t230) + t230 * ((t230 * t69 + (t29 + t422) * t234) * t230 + (-t28 * t234 + (t116 * t337 + t118 * t338) * t229 + (t68 + (qJD(3) * t115 - t118 * t234) * t227 + (qJD(3) * t117 + t116 * t234) * t226) * t230) * t229) + t230 * ((t230 * t82 + (t39 + t423) * t234) * t230 + (-t38 * t234 + (t129 * t335 + t131 * t336) * t229 + (t276 * qJD(3) + t234 * t273 + t81) * t230) * t229) + ((-t132 * t229 + t133 * t230) * (t230 * t319 - t229 * t318 + ((-t132 + t218) * t230 + (t383 - t133 + (t388 + t391) * t230) * t229) * t234) + t339 * t208 * t187) * t406 + t267 + ((-t28 - t38) * t230 + (-t29 - t39) * t229) * t354 + ((t30 + t40) * t230 + (t31 + t41) * t229) * t352; m(5) * ((t234 * t87 + t32) * t230 + (-t234 * t88 + t33) * t229) + m(6) * (-t229 * t308 + t230 * t309); m(6) * (-t229 * t306 + t230 * t307) + m(5) * ((t234 * t89 + t34) * t230 + (-t234 * t90 + t35) * t229); m(6) * (-t229 * t305 + t230 * t304) + m(5) * ((t122 * t234 + t67) * t230 + (-t121 * t234 + t66) * t229); 0; m(6) * ((t229 * t76 - t230 * t75) * t137 + (t229 * t309 + t230 * t308) * t153) + t265; m(6) * ((t229 * t80 - t230 * t79) * t137 + (t229 * t307 + t230 * t306) * t153) + t265; m(6) * (t11 * t10 + t61 * t3 + (t229 * t93 - t230 * t94) * t137 + (t229 * t304 + t230 * t305) * t153) + t267; 0; (t137 * t153 * t339 + t61 * t11) * t404 + t267;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
