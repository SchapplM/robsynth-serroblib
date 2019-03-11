% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:31:35
% EndTime: 2019-03-09 05:32:06
% DurationCPUTime: 14.39s
% Computational Cost: add. (37250->646), mult. (122888->903), div. (0->0), fcn. (105921->14), ass. (0->276)
t264 = sin(pkin(6));
t262 = sin(pkin(12));
t265 = cos(pkin(12));
t268 = sin(qJ(3));
t377 = cos(pkin(7));
t326 = t268 * t377;
t400 = cos(qJ(3));
t286 = t262 * t326 - t400 * t265;
t283 = t264 * t286;
t221 = qJD(1) * t283;
t263 = sin(pkin(7));
t334 = t263 * t400;
t420 = qJD(3) * t334 + t221;
t308 = t377 * t400;
t378 = cos(pkin(6));
t325 = t378 * t263;
t358 = t264 * t265;
t275 = t308 * t358 + t400 * t325;
t347 = qJD(1) * t264;
t331 = t262 * t347;
t199 = -t275 * qJD(1) + t268 * t331;
t196 = qJD(4) + t199;
t267 = sin(qJ(4));
t270 = cos(qJ(4));
t398 = -qJ(5) - pkin(10);
t327 = qJD(4) * t398;
t366 = t199 * t267;
t210 = t264 * (t400 * t262 + t265 * t326) + t268 * t325;
t203 = t210 * qJD(1);
t156 = pkin(3) * t203 + pkin(10) * t199;
t321 = qJD(1) * t378;
t313 = pkin(1) * t321;
t330 = t265 * t347;
t228 = qJ(2) * t330 + t262 * t313;
t281 = (t377 * t358 + t325) * pkin(9);
t189 = qJD(1) * t281 + t228;
t248 = t265 * t313;
t360 = t262 * t264;
t276 = t378 * pkin(2) + (-t377 * pkin(9) - qJ(2)) * t360;
t197 = qJD(1) * t276 + t248;
t361 = t262 * t263;
t222 = (-pkin(2) * t265 - pkin(9) * t361 - pkin(1)) * t264;
t215 = qJD(1) * t222 + qJD(2);
t405 = -t268 * t189 + t197 * t308 + t215 * t334;
t88 = t267 * t156 + t270 * t405;
t419 = -qJ(5) * t366 + t270 * qJD(5) + t267 * t327 - t88;
t87 = t270 * t156 - t267 * t405;
t418 = -pkin(4) * t203 - t267 * qJD(5) - t87 + (-qJ(5) * t199 + t327) * t270;
t261 = sin(pkin(13));
t376 = cos(pkin(13));
t323 = t376 * t267;
t238 = t261 * t270 + t323;
t353 = t196 * t238;
t362 = t261 * t267;
t291 = t376 * t270 - t362;
t153 = t291 * t199;
t232 = t291 * qJD(4);
t417 = t153 + t232;
t359 = t263 * t268;
t234 = t267 * t377 + t270 * t359;
t316 = t263 * t331;
t351 = qJD(4) * t234 + t420 * t267 + t270 * t316;
t233 = -t267 * t359 + t270 * t377;
t350 = -qJD(4) * t233 + t267 * t316 - t420 * t270;
t379 = t418 * t261 + t419 * t376;
t133 = t268 * (t377 * t197 + t215 * t263) + t400 * t189;
t344 = qJD(4) * t267;
t311 = -t133 + (t344 + t366) * pkin(4);
t305 = t378 * t377;
t223 = -qJD(1) * t305 + t263 * t330 - qJD(3);
t166 = t203 * t270 - t223 * t267;
t157 = -t197 * t263 + t377 * t215;
t106 = pkin(3) * t199 - pkin(10) * t203 + t157;
t109 = -t223 * pkin(10) + t133;
t70 = t270 * t106 - t109 * t267;
t62 = -qJ(5) * t166 + t70;
t59 = pkin(4) * t196 + t62;
t164 = t203 * t267 + t223 * t270;
t71 = t106 * t267 + t109 * t270;
t63 = -qJ(5) * t164 + t71;
t60 = t376 * t63;
t26 = t261 * t59 + t60;
t22 = pkin(11) * t196 + t26;
t266 = sin(qJ(6));
t269 = cos(qJ(6));
t113 = t376 * t164 + t166 * t261;
t292 = -t261 * t164 + t376 * t166;
t108 = t223 * pkin(3) - t405;
t85 = t164 * pkin(4) + qJD(5) + t108;
t51 = t113 * pkin(5) - pkin(11) * t292 + t85;
t14 = t22 * t269 + t266 * t51;
t284 = t262 * t308 + t265 * t268;
t282 = t264 * t284;
t277 = qJD(2) * t282;
t403 = qJD(3) * t133;
t103 = qJD(1) * t277 + t403;
t343 = qJD(4) * t270;
t336 = t268 * t360;
t315 = qJD(3) * t336;
t338 = qJD(1) * qJD(3);
t401 = qJD(1) * t315 - t275 * t338;
t138 = t203 * t343 - t223 * t344 - t267 * t401;
t75 = t138 * pkin(4) + t103;
t137 = t203 * t344 + t223 * t343 + t270 * t401;
t83 = -t137 * t261 + t376 * t138;
t84 = -t376 * t137 - t261 * t138;
t33 = t83 * pkin(5) - t84 * pkin(11) + t75;
t202 = t210 * qJD(3);
t185 = qJD(1) * t202;
t278 = qJD(2) * t283;
t402 = qJD(3) * t405;
t102 = -qJD(1) * t278 + t402;
t346 = qJD(2) * t264;
t314 = t346 * t361;
t301 = qJD(1) * t314;
t146 = t185 * pkin(3) + pkin(10) * t401 + t301;
t46 = -qJD(4) * t71 - t267 * t102 + t270 * t146;
t27 = t185 * pkin(4) + t137 * qJ(5) - t166 * qJD(5) + t46;
t45 = t270 * t102 + t106 * t343 - t109 * t344 + t267 * t146;
t31 = -qJ(5) * t138 - qJD(5) * t164 + t45;
t8 = t261 * t27 + t376 * t31;
t6 = pkin(11) * t185 + t8;
t2 = -qJD(6) * t14 - t266 * t6 + t269 * t33;
t410 = qJD(6) + t113;
t416 = t14 * t410 + t2;
t304 = t22 * t266 - t269 * t51;
t1 = -qJD(6) * t304 + t266 * t33 + t269 * t6;
t415 = t304 * t410 + t1;
t414 = -pkin(11) * t203 + t379;
t413 = -t353 * pkin(5) + t417 * pkin(11) - t311;
t412 = t113 * t292;
t354 = t351 * t261 + t350 * t376;
t220 = qJD(1) * t282;
t345 = qJD(3) * t268;
t307 = t263 * t345 - t220;
t318 = t269 * t410;
t411 = -t266 * t83 - t318 * t410;
t380 = t419 * t261 - t418 * t376;
t409 = -t196 * t70 + t45;
t408 = -t71 * t196 - t46;
t406 = t264 ^ 2 * (t262 ^ 2 + t265 ^ 2);
t335 = pkin(1) * t378;
t252 = t265 * t335;
t211 = t252 + t276;
t167 = -t211 * t263 + t377 * t222;
t209 = -t275 + t336;
t123 = pkin(3) * t209 - pkin(10) * t210 + t167;
t349 = qJ(2) * t358 + t262 * t335;
t206 = t281 + t349;
t195 = t400 * t206;
t324 = t377 * t211;
t140 = t222 * t359 + t268 * t324 + t195;
t230 = t263 * t358 - t305;
t131 = -pkin(10) * t230 + t140;
t77 = t267 * t123 + t270 * t131;
t139 = -t268 * t206 + t211 * t308 + t222 * t334;
t93 = -t269 * t196 + t266 * t292;
t95 = t196 * t266 + t269 * t292;
t399 = t95 * t93;
t257 = -pkin(4) * t270 - pkin(3);
t204 = -pkin(5) * t291 - pkin(11) * t238 + t257;
t244 = t398 * t270;
t217 = -t376 * t244 + t398 * t362;
t159 = t204 * t266 + t217 * t269;
t397 = qJD(6) * t159 + t414 * t266 + t413 * t269;
t158 = t204 * t269 - t217 * t266;
t396 = -qJD(6) * t158 + t413 * t266 - t414 * t269;
t169 = t210 * t267 + t230 * t270;
t201 = -qJD(3) * t275 + t315;
t145 = -qJD(4) * t169 - t201 * t270;
t170 = t210 * t270 - t230 * t267;
t118 = qJD(3) * t139 - t278;
t151 = pkin(3) * t202 + pkin(10) * t201 + t314;
t50 = -qJD(4) * t77 - t267 * t118 + t270 * t151;
t37 = t202 * pkin(4) - t145 * qJ(5) - t170 * qJD(5) + t50;
t144 = qJD(4) * t170 - t201 * t267;
t49 = t270 * t118 + t123 * t343 - t131 * t344 + t267 * t151;
t41 = -qJ(5) * t144 - qJD(5) * t169 + t49;
t12 = t261 * t37 + t376 * t41;
t341 = qJD(6) * t269;
t322 = -t269 * t185 + t266 * t84;
t54 = qJD(6) * t95 + t322;
t395 = -t266 * t54 - t93 * t341;
t76 = t270 * t123 - t131 * t267;
t65 = pkin(4) * t209 - qJ(5) * t170 + t76;
t69 = -qJ(5) * t169 + t77;
t39 = t261 * t65 + t376 * t69;
t394 = t292 * t93;
t134 = t376 * t169 + t170 * t261;
t393 = t134 * t83;
t391 = t261 * t63;
t389 = t266 * t93;
t388 = t266 * t95;
t387 = t269 * t93;
t342 = qJD(6) * t266;
t53 = -t266 * t185 - t196 * t341 - t269 * t84 + t292 * t342;
t386 = t53 * t266;
t385 = t54 * t269;
t383 = t83 * t291;
t382 = t95 * t292;
t381 = t203 * pkin(5) + t380;
t375 = t292 ^ 2;
t374 = t292 * t196;
t373 = t113 ^ 2;
t372 = t113 * t196;
t371 = t164 * t196;
t370 = t166 * t164;
t369 = t166 * t196;
t368 = t185 * t209;
t367 = t196 * t203;
t365 = t203 * t199;
t364 = t238 * t266;
t363 = t238 * t269;
t187 = t261 * t233 + t376 * t234;
t171 = -t266 * t187 - t269 * t334;
t357 = qJD(6) * t171 + t307 * t266 - t354 * t269;
t290 = -t269 * t187 + t266 * t334;
t356 = qJD(6) * t290 + t354 * t266 + t307 * t269;
t355 = t350 * t261 - t351 * t376;
t339 = qJD(1) * qJD(2);
t333 = t400 * t103;
t328 = t264 * t339;
t7 = -t261 * t31 + t376 * t27;
t124 = -t153 * t266 - t269 * t203;
t320 = -t232 * t266 + t124;
t125 = -t153 * t269 + t203 * t266;
t319 = -t232 * t269 + t125;
t317 = t270 * t196;
t271 = qJD(1) ^ 2;
t310 = t264 * t271 * t378;
t35 = pkin(11) * t209 + t39;
t135 = -t261 * t169 + t376 * t170;
t130 = t230 * pkin(3) - t139;
t92 = t169 * pkin(4) + t130;
t55 = t134 * pkin(5) - t135 * pkin(11) + t92;
t16 = t266 * t55 + t269 * t35;
t15 = -t266 * t35 + t269 * t55;
t303 = t387 + t388;
t100 = t135 * t269 + t209 * t266;
t302 = (-qJ(2) * t331 + t248) * t262 - t228 * t265;
t300 = t269 * t83 + (-t113 * t266 - t342) * t410;
t298 = -t269 * t53 - t342 * t95;
t297 = -0.2e1 * t321 * t346;
t11 = -t261 * t41 + t376 * t37;
t25 = t376 * t59 - t391;
t38 = -t261 * t69 + t376 * t65;
t295 = -pkin(10) * t185 + t108 * t196;
t21 = -t196 * pkin(5) - t25;
t255 = pkin(4) * t261 + pkin(11);
t294 = t21 * t410 - t255 * t83;
t289 = t238 * t341 - t320;
t288 = -t238 * t342 - t319;
t119 = t277 + (t195 + (t222 * t263 + t324) * t268) * qJD(3);
t80 = t144 * pkin(4) + t119;
t256 = -t376 * pkin(4) - pkin(5);
t216 = -t244 * t261 - t398 * t323;
t186 = -t376 * t233 + t234 * t261;
t107 = t196 * t202 + t368;
t99 = t135 * t266 - t209 * t269;
t90 = -t261 * t144 + t376 * t145;
t89 = t376 * t144 + t145 * t261;
t66 = pkin(4) * t166 + pkin(5) * t292 + pkin(11) * t113;
t58 = qJD(6) * t100 - t202 * t269 + t266 * t90;
t57 = t135 * t342 - t202 * t266 - t209 * t341 - t269 * t90;
t42 = t89 * pkin(5) - t90 * pkin(11) + t80;
t34 = -t209 * pkin(5) - t38;
t30 = t376 * t62 - t391;
t29 = t261 * t62 + t60;
t18 = t266 * t66 + t269 * t30;
t17 = -t266 * t30 + t269 * t66;
t10 = pkin(11) * t202 + t12;
t9 = -t202 * pkin(5) - t11;
t5 = -pkin(5) * t185 - t7;
t4 = -qJD(6) * t16 - t266 * t10 + t269 * t42;
t3 = qJD(6) * t15 + t269 * t10 + t266 * t42;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262 * t297, t265 * t297, 0.2e1 * t339 * t406 ((t265 * t349 + (qJ(2) * t360 - t252) * t262) * qJD(1) - t302) * t346, -t203 * t201 - t210 * t401, -t210 * t185 + t201 * t199 - t203 * t202 + t209 * t401, t201 * t223 + t230 * t401, t199 * t202 + t368, t185 * t230 + t202 * t223, 0, t103 * t230 + t119 * t223 + t157 * t202 + t167 * t185 + (qJD(1) * t209 + t199) * t314, t102 * t230 + t118 * t223 - t157 * t201 - t167 * t401 + t203 * t314 + t210 * t301, -t102 * t209 + t103 * t210 - t118 * t199 + t119 * t203 - t133 * t202 + t139 * t401 - t140 * t185 + t201 * t405, t102 * t140 - t103 * t139 + t133 * t118 - t405 * t119 + (qJD(1) * t167 + t157) * t314, -t137 * t170 + t145 * t166, t137 * t169 - t138 * t170 - t144 * t166 - t145 * t164, -t137 * t209 + t145 * t196 + t166 * t202 + t170 * t185, t138 * t169 + t144 * t164, -t138 * t209 - t144 * t196 - t164 * t202 - t169 * t185, t107, t103 * t169 + t108 * t144 + t119 * t164 + t130 * t138 + t185 * t76 + t196 * t50 + t202 * t70 + t209 * t46, t103 * t170 + t108 * t145 + t119 * t166 - t130 * t137 - t185 * t77 - t196 * t49 - t202 * t71 - t209 * t45, t137 * t76 - t138 * t77 - t144 * t71 - t145 * t70 - t164 * t49 - t166 * t50 - t169 * t45 - t170 * t46, t103 * t130 + t108 * t119 + t45 * t77 + t46 * t76 + t49 * t71 + t50 * t70, t135 * t84 + t292 * t90, -t113 * t90 - t134 * t84 - t135 * t83 - t292 * t89, t135 * t185 + t196 * t90 + t202 * t292 + t209 * t84, t113 * t89 + t393, -t113 * t202 - t134 * t185 - t196 * t89 - t209 * t83, t107, t11 * t196 + t113 * t80 + t134 * t75 + t185 * t38 + t202 * t25 + t209 * t7 + t83 * t92 + t85 * t89, -t12 * t196 + t135 * t75 - t185 * t39 - t202 * t26 - t209 * t8 + t292 * t80 + t84 * t92 + t85 * t90, -t11 * t292 - t113 * t12 - t134 * t8 - t135 * t7 - t25 * t90 - t26 * t89 - t38 * t84 - t39 * t83, t11 * t25 + t12 * t26 + t38 * t7 + t39 * t8 + t75 * t92 + t80 * t85, -t100 * t53 - t57 * t95, -t100 * t54 + t53 * t99 + t57 * t93 - t58 * t95, t100 * t83 - t134 * t53 - t410 * t57 + t89 * t95, t54 * t99 + t58 * t93, -t134 * t54 - t410 * t58 - t83 * t99 - t89 * t93, t410 * t89 + t393, t134 * t2 + t15 * t83 + t21 * t58 - t304 * t89 + t34 * t54 + t4 * t410 + t5 * t99 + t9 * t93, -t1 * t134 + t100 * t5 - t14 * t89 - t16 * t83 - t21 * t57 - t3 * t410 - t34 * t53 + t9 * t95, -t1 * t99 - t100 * t2 - t14 * t58 + t15 * t53 - t16 * t54 - t3 * t93 - t304 * t57 - t4 * t95, t1 * t16 + t14 * t3 + t15 * t2 + t21 * t9 - t304 * t4 + t34 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262 * t310, t265 * t310, -t271 * t406, t302 * t347, 0, 0, 0, 0, 0, 0, t377 * t185 - t199 * t316 + t307 * t223, -t203 * t316 + t223 * t420 - t377 * t401, -t185 * t359 - t199 * t420 + t203 * t307 + t334 * t401, t405 * t220 + t133 * t221 + (-t333 + t102 * t268 + (t400 * t133 - t268 * t405) * qJD(3) + (t377 * qJD(2) - t157) * t331) * t263, 0, 0, 0, 0, 0, 0, -t138 * t334 + t164 * t307 + t233 * t185 - t196 * t351, t137 * t334 + t166 * t307 - t234 * t185 + t196 * t350, t233 * t137 - t234 * t138 + t164 * t350 + t166 * t351, t108 * t307 + t46 * t233 + t45 * t234 - t263 * t333 - t350 * t71 - t351 * t70, 0, 0, 0, 0, 0, 0, t113 * t307 - t186 * t185 + t196 * t355 - t334 * t83, -t187 * t185 + t196 * t354 + t292 * t307 - t334 * t84, t113 * t354 + t186 * t84 - t187 * t83 - t292 * t355, -t7 * t186 + t8 * t187 - t85 * t220 + (t85 * t345 - t400 * t75) * t263 - t354 * t26 + t355 * t25, 0, 0, 0, 0, 0, 0, t171 * t83 + t186 * t54 - t355 * t93 + t356 * t410, -t186 * t53 + t290 * t83 - t355 * t95 - t357 * t410, t171 * t53 + t290 * t54 - t356 * t95 - t357 * t93, -t1 * t290 + t14 * t357 + t2 * t171 + t5 * t186 - t21 * t355 - t304 * t356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t365, -t199 ^ 2 + t203 ^ 2, -t199 * t223 - t401, -t365, -t203 * t223 - t210 * t338, 0, -t133 * t223 - t157 * t203 - t284 * t328 - t403, t157 * t199 - t223 * t405 + t286 * t328 - t402, 0, 0, -t137 * t267 + t166 * t317 (-t137 - t371) * t270 + (-t138 - t369) * t267, -t166 * t203 + t267 * t185 + t196 * t317, -t138 * t270 + t267 * t371, -t196 ^ 2 * t267 + t164 * t203 + t270 * t185, -t367, -pkin(3) * t138 - t103 * t270 - t133 * t164 - t70 * t203 + (-pkin(10) * t343 - t87) * t196 + t295 * t267, pkin(3) * t137 + t103 * t267 - t133 * t166 + t71 * t203 + (pkin(10) * t344 + t88) * t196 + t295 * t270, t164 * t88 + t166 * t87 + ((qJD(4) * t166 - t138) * pkin(10) + t409) * t270 + ((qJD(4) * t164 - t137) * pkin(10) + t408) * t267, -t103 * pkin(3) - t108 * t133 - t70 * t87 - t71 * t88 + (-t46 * t267 + t45 * t270 + (-t267 * t71 - t270 * t70) * qJD(4)) * pkin(10), t84 * t238 + t292 * t417, -t113 * t417 - t238 * t83 + t291 * t84 - t292 * t353, t238 * t185 + t196 * t417 - t203 * t292, t113 * t353 - t383, t113 * t203 + t185 * t291 - t196 * t353, -t367, t113 * t311 - t185 * t216 - t196 * t380 - t203 * t25 + t257 * t83 - t291 * t75 + t353 * t85, -t185 * t217 - t196 * t379 + t203 * t26 + t238 * t75 + t257 * t84 + t292 * t311 + t417 * t85, -t113 * t379 + t216 * t84 - t217 * t83 - t238 * t7 - t25 * t417 - t26 * t353 + t291 * t8 + t292 * t380, -t216 * t7 + t217 * t8 - t25 * t380 + t257 * t75 + t26 * t379 + t311 * t85, t288 * t95 - t363 * t53, t95 * t124 + t125 * t93 - t303 * t232 + (t386 - t385 + (-t269 * t95 + t389) * qJD(6)) * t238, t288 * t410 + t291 * t53 + t353 * t95 + t363 * t83, t289 * t93 + t364 * t54, -t289 * t410 + t291 * t54 - t353 * t93 - t364 * t83, t353 * t410 - t383, t158 * t83 - t2 * t291 + t21 * t289 + t216 * t54 - t304 * t353 + t364 * t5 + t381 * t93 - t397 * t410, t1 * t291 - t14 * t353 - t159 * t83 + t21 * t288 - t216 * t53 + t363 * t5 + t381 * t95 + t396 * t410, t158 * t53 - t159 * t54 + t397 * t95 + t396 * t93 + t320 * t14 - t319 * t304 + (-t1 * t266 - t2 * t269 + (-t14 * t269 - t266 * t304) * qJD(6)) * t238, t1 * t159 - t14 * t396 + t2 * t158 + t21 * t381 + t5 * t216 + t304 * t397; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t370, -t164 ^ 2 + t166 ^ 2, -t137 + t371, -t370, -t138 + t369, t185, -t108 * t166 - t408, t108 * t164 - t409, 0, 0, t412, -t373 + t375, t84 + t372, -t412, -t83 + t374, t185, -t85 * t292 + t29 * t196 + (-t113 * t166 + t185 * t376) * pkin(4) + t7, t85 * t113 + t30 * t196 + (-t166 * t292 - t185 * t261) * pkin(4) - t8 (-t261 * t83 - t376 * t84) * pkin(4) + (t26 - t29) * t292 + (t30 - t25) * t113, t25 * t29 - t26 * t30 + (-t166 * t85 + t261 * t8 + t376 * t7) * pkin(4), t318 * t95 - t386, -t113 * t303 + t298 + t395, -t382 - t411, t389 * t410 - t385, t300 + t394, -t410 * t292, t304 * t292 + t256 * t54 - t5 * t269 - t29 * t93 + (-t255 * t341 - t17) * t410 + t294 * t266, t292 * t14 - t256 * t53 + t266 * t5 - t29 * t95 + (t255 * t342 + t18) * t410 + t294 * t269, t17 * t95 + t18 * t93 + (t113 * t304 - t255 * t54 + t1 + (t255 * t95 + t304) * qJD(6)) * t269 + (-t113 * t14 - t255 * t53 - t2 + (t255 * t93 - t14) * qJD(6)) * t266, t304 * t17 - t14 * t18 - t21 * t29 + t5 * t256 + (t1 * t269 - t2 * t266 + (-t14 * t266 + t269 * t304) * qJD(6)) * t255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 + t374, t84 - t372, -t373 - t375, t113 * t26 + t25 * t292 + t75, 0, 0, 0, 0, 0, 0, t300 - t394, -t382 + t411 -(t387 - t388) * t113 - t298 + t395, -t21 * t292 + t415 * t266 + t416 * t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t399, -t93 ^ 2 + t95 ^ 2, t410 * t93 - t53, -t399, -t322 + (-qJD(6) + t410) * t95, t83, -t21 * t95 + t416, t21 * t93 - t415, 0, 0;];
tauc_reg  = t13;
