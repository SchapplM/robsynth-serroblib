% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRRR7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:59:46
% EndTime: 2019-03-09 14:00:03
% DurationCPUTime: 7.13s
% Computational Cost: add. (14025->572), mult. (30851->747), div. (0->0), fcn. (20801->8), ass. (0->280)
t235 = sin(qJ(4));
t239 = cos(qJ(4));
t240 = cos(qJ(2));
t318 = qJD(1) * t240;
t236 = sin(qJ(2));
t319 = qJD(1) * t236;
t162 = -t235 * t318 + t239 * t319;
t227 = qJD(2) - qJD(4);
t234 = sin(qJ(5));
t238 = cos(qJ(5));
t126 = t162 * t234 + t227 * t238;
t330 = t239 * t240;
t332 = t235 * t236;
t171 = t330 + t332;
t159 = t171 * qJD(1);
t409 = qJD(5) + t159;
t388 = t409 * t126;
t128 = t162 * t238 - t227 * t234;
t387 = t409 * t128;
t216 = pkin(7) * t319;
t178 = pkin(8) * t319 - t216;
t419 = qJD(3) - t178;
t164 = -qJD(1) * pkin(1) - pkin(2) * t318 - qJ(3) * t319;
t138 = pkin(3) * t318 - t164;
t81 = pkin(4) * t159 - pkin(9) * t162 + t138;
t241 = -pkin(2) - pkin(3);
t299 = t241 * qJD(2);
t144 = t299 + t419;
t229 = qJD(2) * qJ(3);
t217 = pkin(7) * t318;
t273 = -pkin(8) * t318 + t217;
t163 = t229 + t273;
t102 = t144 * t235 + t163 * t239;
t89 = -pkin(9) * t227 + t102;
t39 = -t234 * t89 + t238 * t81;
t311 = qJD(5) * t238;
t312 = qJD(5) * t234;
t313 = qJD(4) * t239;
t315 = qJD(2) * t240;
t252 = t235 * t315 + t236 * t313;
t306 = qJD(1) * qJD(4);
t293 = t240 * t306;
t307 = qJD(1) * qJD(2);
t295 = t236 * t307;
t112 = qJD(1) * t252 - t235 * t293 - t239 * t295;
t220 = t236 * qJD(3);
t294 = t240 * t307;
t322 = qJ(3) * t294 + qJD(1) * t220;
t323 = t239 * t293 + t306 * t332;
t45 = t323 * pkin(9) + t112 * pkin(4) + (-pkin(9) * t330 + (-pkin(9) * t235 + t241) * t236) * t307 + t322;
t317 = qJD(2) * t236;
t381 = pkin(7) - pkin(8);
t179 = t381 * t317;
t228 = qJD(2) * qJD(3);
t149 = -qJD(1) * t179 + t228;
t208 = pkin(7) * t294;
t167 = -pkin(8) * t294 + t208;
t314 = qJD(4) * t235;
t53 = t144 * t313 + t239 * t149 - t163 * t314 + t235 * t167;
t9 = t234 * t45 + t238 * t53 + t81 * t311 - t312 * t89;
t418 = -t409 * t39 + t9;
t40 = t234 * t81 + t238 * t89;
t10 = -qJD(5) * t40 - t234 * t53 + t238 * t45;
t417 = -t409 * t40 - t10;
t253 = t171 * qJD(2);
t248 = qJD(1) * t253 - t323;
t301 = t162 * t311 - t227 * t312 + t234 * t248;
t63 = t162 * t312 + t227 * t311 - t238 * t248;
t416 = (t387 + t301) * t234 + (t63 + t388) * t238;
t183 = -t235 * qJ(3) + t239 * t241;
t145 = t239 * qJD(3) + qJD(4) * t183;
t258 = pkin(10) * t159 * t238 + pkin(5) * t162;
t184 = t239 * qJ(3) + t235 * t241;
t177 = -pkin(9) + t184;
t377 = pkin(10) - t177;
t289 = qJD(5) * t377;
t116 = t239 * t178 + t235 * t273;
t113 = pkin(4) * t162 + pkin(9) * t159;
t211 = qJ(3) * t318;
t152 = t241 * t319 + t211;
t82 = -t113 + t152;
t47 = -t116 * t234 + t238 * t82;
t415 = t234 * t145 - t238 * t289 - t258 + t47;
t342 = t159 * t234;
t303 = pkin(10) * t342;
t48 = t238 * t116 + t234 * t82;
t414 = -t238 * t145 - t234 * t289 - t303 + t48;
t233 = sin(qJ(6));
t237 = cos(qJ(6));
t259 = t126 * t233 - t237 * t128;
t69 = t237 * t126 + t128 * t233;
t378 = t69 * t259;
t380 = -pkin(10) - pkin(9);
t300 = qJD(5) * t380;
t101 = t144 * t239 - t235 * t163;
t56 = t238 * t101 + t234 * t113;
t413 = t234 * t300 - t303 - t56;
t55 = -t101 * t234 + t238 * t113;
t412 = -t238 * t300 + t258 + t55;
t408 = t259 ^ 2 - t69 ^ 2;
t33 = -pkin(10) * t128 + t39;
t26 = pkin(5) * t409 + t33;
t34 = -pkin(10) * t126 + t40;
t368 = t237 * t34;
t12 = t233 * t26 + t368;
t6 = pkin(5) * t112 + pkin(10) * t63 + t10;
t7 = -pkin(10) * t301 + t9;
t2 = -qJD(6) * t12 - t233 * t7 + t237 * t6;
t88 = pkin(4) * t227 - t101;
t59 = pkin(5) * t126 + t88;
t407 = t59 * t259 + t2;
t148 = qJD(6) + t409;
t309 = qJD(6) * t237;
t310 = qJD(6) * t233;
t24 = t126 * t309 + t128 * t310 + t233 * t301 + t237 * t63;
t406 = t148 * t69 - t24;
t1 = (qJD(6) * t26 + t7) * t237 + t233 * t6 - t34 * t310;
t405 = t59 * t69 - t1;
t335 = t233 * t238;
t173 = t234 * t237 + t335;
t54 = qJD(4) * t102 + t235 * t149 - t239 * t167;
t27 = pkin(5) * t301 + t54;
t336 = t233 * t234;
t170 = -t237 * t238 + t336;
t305 = qJD(5) + qJD(6);
t383 = t170 * t159 - t237 * t311 - t238 * t309 + t305 * t336;
t404 = t12 * t162 + t173 * t27 - t383 * t59;
t370 = t233 * t34;
t11 = t237 * t26 - t370;
t119 = t305 * t173;
t382 = t173 * t159 + t119;
t403 = -t11 * t162 + t170 * t27 + t382 * t59;
t246 = qJD(6) * t259 + t233 * t63 - t237 * t301;
t402 = -t170 * t246 + t382 * t69;
t401 = -t148 * t259 + t246;
t400 = -t112 * t173 + t148 * t383 - t162 * t259;
t399 = t112 * t170 + t148 * t382 - t162 * t69;
t354 = t112 * t238;
t389 = t409 ^ 2;
t398 = -t126 * t162 + t234 * t389 - t354;
t272 = t301 * t238;
t397 = t234 * t388 - t272;
t334 = t234 * t112;
t396 = -t128 * t162 + t238 * t389 + t334;
t364 = t63 * t234;
t395 = t238 * t387 - t364;
t394 = -t24 * t173 + t259 * t383;
t393 = -t1 * t170 + t11 * t383 - t12 * t382 - t173 * t2;
t392 = t170 * t24 + t173 * t246 + t259 * t382 + t383 * t69;
t391 = t409 * t88;
t192 = t381 * t236;
t194 = t381 * t240;
t132 = t192 * t235 + t194 * t239;
t124 = t238 * t132;
t186 = -t240 * pkin(2) - t236 * qJ(3) - pkin(1);
t168 = t240 * pkin(3) - t186;
t172 = -t235 * t240 + t236 * t239;
t99 = pkin(4) * t171 - pkin(9) * t172 + t168;
t58 = t234 * t99 + t124;
t390 = -t102 * t227 - t54;
t154 = t170 * t235;
t153 = t173 * t235;
t386 = t162 * t227 + t112;
t325 = qJD(4) * t184 + t419 * t235 + t239 * t273;
t385 = t239 * t192 - t194 * t235;
t384 = (t312 + t342) * pkin(5);
t379 = pkin(5) * t238;
t139 = t377 * t234;
t140 = t377 * t238;
t86 = t139 * t233 - t140 * t237;
t376 = qJD(6) * t86 - t414 * t233 + t415 * t237;
t85 = t139 * t237 + t140 * t233;
t375 = -qJD(6) * t85 + t415 * t233 + t414 * t237;
t191 = t380 * t234;
t193 = t380 * t238;
t131 = t191 * t233 - t193 * t237;
t374 = qJD(6) * t131 + t413 * t233 + t412 * t237;
t129 = t191 * t237 + t193 * t233;
t373 = -qJD(6) * t129 + t412 * t233 - t413 * t237;
t372 = qJD(2) * pkin(2);
t371 = t385 * t54;
t369 = t234 * t54;
t365 = t54 * t238;
t316 = qJD(2) * t239;
t158 = -t234 * t316 + t238 * t319;
t161 = t234 * t319 + t238 * t316;
t363 = -t154 * t305 + t158 * t237 - t161 * t233 + t173 * t313;
t362 = t153 * t305 + t158 * t233 + t161 * t237 + t170 * t313;
t357 = -t384 + t325;
t356 = t101 * t227;
t84 = t112 * t171;
t121 = -qJD(4) * t171 + t253;
t353 = t121 * t234;
t352 = t121 * t238;
t351 = t126 * t234;
t350 = t126 * t238;
t349 = t128 * t126;
t348 = t128 * t234;
t347 = t128 * t238;
t346 = t138 * t162;
t345 = t148 * t162;
t344 = t409 * t162;
t343 = t159 * t227;
t341 = t162 * t159;
t339 = t172 * t234;
t338 = t172 * t238;
t333 = t235 * t162;
t243 = qJD(1) ^ 2;
t328 = t240 * t243;
t242 = qJD(2) ^ 2;
t327 = t242 * t236;
t221 = t242 * t240;
t326 = -t116 + t145;
t321 = qJ(3) * t315 + t220;
t230 = t236 ^ 2;
t320 = t240 ^ 2 - t230;
t298 = t172 * t312;
t297 = t159 ^ 2 - t162 ^ 2;
t120 = -t236 * t316 - t240 * t314 + t252;
t276 = t236 * t299;
t135 = t276 + t321;
t52 = pkin(4) * t120 - pkin(9) * t121 + t135;
t277 = qJD(2) * t194;
t76 = qJD(4) * t385 - t239 * t179 + t235 * t277;
t290 = -t234 * t76 + t238 * t52;
t57 = -t132 * t234 + t238 * t99;
t288 = -0.2e1 * pkin(1) * t307;
t287 = qJD(3) - t372;
t279 = qJD(1) * t186 + t164;
t278 = t227 * t235;
t275 = t236 * t294;
t274 = -t102 + t384;
t176 = pkin(4) - t183;
t268 = t162 * t40 + t369;
t267 = -t39 * t162 - t365;
t41 = pkin(5) * t171 - pkin(10) * t338 + t57;
t46 = -pkin(10) * t339 + t58;
t22 = -t233 * t46 + t237 * t41;
t23 = t233 * t41 + t237 * t46;
t265 = -t234 * t40 - t238 * t39;
t264 = t234 * t39 - t238 * t40;
t137 = pkin(2) * t295 - t322;
t155 = pkin(2) * t317 - t321;
t257 = -pkin(7) * t242 - qJD(1) * t155 - t137;
t256 = t172 * t311 + t353;
t255 = -t298 + t352;
t20 = -t132 * t312 + t234 * t52 + t238 * t76 + t99 * t311;
t254 = -pkin(9) * t112 + t391;
t250 = -t177 * t112 - t145 * t409 - t391;
t249 = -t138 * t159 + t53;
t247 = -qJD(2) * t159 + t323;
t77 = qJD(4) * t132 - t179 * t235 - t239 * t277;
t245 = qJD(5) * t265 - t10 * t234 + t9 * t238;
t182 = -pkin(7) * t295 + t228;
t185 = t216 + t287;
t188 = t217 + t229;
t244 = t182 * t240 + (t185 * t240 + (-t188 + t217) * t236) * qJD(2);
t214 = -pkin(4) - t379;
t203 = t236 * t328;
t190 = -0.2e1 * t275;
t189 = 0.2e1 * t275;
t187 = t320 * t243;
t175 = pkin(2) * t319 - t211;
t169 = t320 * t307;
t157 = t176 + t379;
t125 = qJD(1) * t276 + t322;
t110 = t170 * t172;
t109 = t173 * t172;
t100 = pkin(5) * t339 - t385;
t37 = pkin(5) * t256 + t77;
t29 = t121 * t335 - t233 * t298 - t310 * t339 + (t305 * t338 + t353) * t237;
t28 = t119 * t172 + t121 * t170;
t21 = -qJD(5) * t58 + t290;
t15 = t237 * t33 - t370;
t14 = -t233 * t33 - t368;
t13 = -pkin(10) * t256 + t20;
t8 = -pkin(10) * t352 + pkin(5) * t120 + (-t124 + (pkin(10) * t172 - t99) * t234) * qJD(5) + t290;
t4 = -qJD(6) * t23 - t13 * t233 + t237 * t8;
t3 = qJD(6) * t22 + t13 * t237 + t233 * t8;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t189, 0.2e1 * t169, t221, t190, -t327, 0, -pkin(7) * t221 + t236 * t288, pkin(7) * t327 + t240 * t288, 0, 0, t189, t221, -0.2e1 * t169, 0, t327, t190, t240 * t257 + t279 * t317, t244, t236 * t257 - t279 * t315, pkin(7) * t244 + t137 * t186 + t155 * t164, t162 * t121 + t172 * t248, -t172 * t112 - t162 * t120 - t121 * t159 - t171 * t248, -t121 * t227, t120 * t159 + t84, t120 * t227, 0, t112 * t168 + t120 * t138 + t125 * t171 + t135 * t159 + t227 * t77, t138 * t121 + t125 * t172 + t135 * t162 + t168 * t248 + t76 * t227, -t101 * t121 - t102 * t120 - t132 * t112 - t76 * t159 + t77 * t162 - t53 * t171 + t54 * t172 + t247 * t385, -t101 * t77 + t102 * t76 + t125 * t168 + t132 * t53 + t135 * t138 - t371, t128 * t255 - t338 * t63 (-t348 - t350) * t121 + (-t272 + t364 + (-t347 + t351) * qJD(5)) * t172, t112 * t338 + t120 * t128 - t171 * t63 + t255 * t409, t126 * t256 + t301 * t339, -t126 * t120 - t171 * t301 - t172 * t334 - t256 * t409, t120 * t409 + t84, t21 * t409 + t57 * t112 + t10 * t171 + t39 * t120 + t77 * t126 - t385 * t301 + t88 * t353 + (t311 * t88 + t369) * t172, t88 * t352 - t112 * t58 - t120 * t40 + t128 * t77 + t385 * t63 - t409 * t20 - t171 * t9 + (-t312 * t88 + t365) * t172, -t20 * t126 - t58 * t301 - t21 * t128 + t57 * t63 + t265 * t121 + (qJD(5) * t264 - t10 * t238 - t9 * t234) * t172, t10 * t57 + t20 * t40 + t21 * t39 + t58 * t9 + t77 * t88 - t371, t110 * t24 + t259 * t28, t109 * t24 - t110 * t246 + t259 * t29 + t28 * t69, -t110 * t112 - t120 * t259 - t148 * t28 - t171 * t24, -t109 * t246 + t29 * t69, -t109 * t112 - t120 * t69 - t148 * t29 + t171 * t246, t120 * t148 + t84, -t100 * t246 + t109 * t27 + t11 * t120 + t112 * t22 + t148 * t4 + t171 * t2 + t29 * t59 + t37 * t69, -t1 * t171 - t100 * t24 - t110 * t27 - t112 * t23 - t12 * t120 - t148 * t3 - t259 * t37 - t28 * t59, -t1 * t109 + t11 * t28 + t110 * t2 - t12 * t29 + t22 * t24 + t23 * t246 + t259 * t4 - t3 * t69, t1 * t23 + t100 * t27 + t11 * t4 + t12 * t3 + t2 * t22 + t37 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203, -t187, 0, t203, 0, 0, t243 * pkin(1) * t236, pkin(1) * t328, 0, 0, -t203, 0, t187, 0, 0, t203 (-t164 * t236 + t175 * t240) * qJD(1) ((t188 - t229) * t236 + (-t185 + t287) * t240) * qJD(1), 0.2e1 * t228 + (t164 * t240 + t175 * t236) * qJD(1), qJ(3) * t182 + qJD(3) * t188 - t164 * t175 + (t188 * t236 + (-t185 - t372) * t240) * qJD(1) * pkin(7), -t341, t297, t247 + t343, t341, t386, 0, -t152 * t159 + t227 * t325 + t346 + t54, -t152 * t162 + t227 * t326 + t249, -t184 * t112 + t183 * t247 + (-t102 + t325) * t162 + (t101 - t326) * t159, -t101 * t325 + t102 * t326 - t138 * t152 - t183 * t54 + t184 * t53, -t395, t416, -t396, -t397, t398, t344, t176 * t301 + (-t177 * t311 - t47) * t409 + t325 * t126 + t250 * t234 - t267, -t176 * t63 + (t177 * t312 + t48) * t409 + t325 * t128 + t250 * t238 - t268, t48 * t126 + t47 * t128 + (-t145 * t126 - t177 * t301 - t9 + t39 * t159 + (t128 * t177 + t39) * qJD(5)) * t238 + (t145 * t128 + t40 * t159 - t177 * t63 + t10 + (t126 * t177 + t40) * qJD(5)) * t234, -t145 * t264 + t176 * t54 + t177 * t245 + t325 * t88 - t39 * t47 - t40 * t48, -t394, -t392, t400, -t402, t399, t345, t112 * t85 - t148 * t376 - t157 * t246 + t357 * t69 - t403, -t112 * t86 + t148 * t375 - t157 * t24 - t259 * t357 - t404, t24 * t85 + t246 * t86 - t259 * t376 + t375 * t69 - t393, t1 * t86 - t11 * t376 - t12 * t375 + t157 * t27 + t2 * t85 + t357 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203, 0, -t230 * t243 - t242, -qJD(2) * t188 + t164 * t319 + t208, 0, 0, 0, 0, 0, 0, -t159 * t319 - t227 * t278, -t227 ^ 2 * t239 - t162 * t319, -t235 * t112 + t239 * t323 + (-t239 * t159 + t333) * qJD(4) - qJD(2) * t333, -t138 * t319 + t390 * t239 + (t53 + t356) * t235, 0, 0, 0, 0, 0, 0, -t239 * t301 + (-t234 * t313 - t158) * t409 + (-t126 * t227 - t311 * t409 - t334) * t235, t239 * t63 + (-t238 * t313 + t161) * t409 + (-t128 * t227 + t312 * t409 - t354) * t235, t161 * t126 + t158 * t128 + (t348 - t350) * t313 + (-t272 - t364 + (t347 + t351) * qJD(5)) * t235, -t158 * t39 - t161 * t40 + (-qJD(4) * t264 - t54) * t239 + (-t227 * t88 + t245) * t235, 0, 0, 0, 0, 0, 0, -t112 * t153 - t148 * t363 + t239 * t246 - t278 * t69, t112 * t154 + t148 * t362 + t239 * t24 + t259 * t278, -t153 * t24 - t154 * t246 - t259 * t363 + t362 * t69, -t1 * t154 - t11 * t363 - t12 * t362 - t153 * t2 - t239 * t27 - t278 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341, -t297, t248 - t343, -t341, -t386, 0, -t346 + t390, -t249 - t356, 0, 0, t395, -t416, t396, t397, -t398, -t344, -pkin(4) * t301 - t102 * t126 + (-pkin(9) * t311 - t55) * t409 + t254 * t234 + t267, pkin(4) * t63 - t102 * t128 + (pkin(9) * t312 + t56) * t409 + t254 * t238 + t268, t56 * t126 + t55 * t128 + ((qJD(5) * t128 - t301) * pkin(9) + t418) * t238 + ((qJD(5) * t126 - t63) * pkin(9) + t417) * t234, -pkin(4) * t54 + pkin(9) * t245 - t102 * t88 - t39 * t55 - t40 * t56, t394, t392, -t400, t402, -t399, -t345, t112 * t129 - t148 * t374 - t214 * t246 + t274 * t69 + t403, -t112 * t131 + t148 * t373 - t214 * t24 - t259 * t274 + t404, t129 * t24 + t131 * t246 - t259 * t374 + t373 * t69 + t393, t1 * t131 - t11 * t374 - t12 * t373 + t129 * t2 + t214 * t27 + t274 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t349, -t126 ^ 2 + t128 ^ 2, -t63 + t388, -t349, -t301 + t387, t112, -t128 * t88 - t417, t126 * t88 - t418, 0, 0, -t378, t408, t406, t378, t401, t112, -t14 * t148 + (t112 * t237 - t128 * t69 - t148 * t310) * pkin(5) + t407, t148 * t15 + (-t112 * t233 + t128 * t259 - t148 * t309) * pkin(5) + t405, -t11 * t69 - t12 * t259 - t14 * t259 + t15 * t69 + (t233 * t246 + t237 * t24 + (-t233 * t259 - t237 * t69) * qJD(6)) * pkin(5), -t11 * t14 - t12 * t15 + (t1 * t233 - t128 * t59 + t2 * t237 + (-t11 * t233 + t12 * t237) * qJD(6)) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t378, t408, t406, t378, t401, t112, t12 * t148 + t407, t11 * t148 + t405, 0, 0;];
tauc_reg  = t5;
