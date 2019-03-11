% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:27:38
% EndTime: 2019-03-08 22:27:57
% DurationCPUTime: 7.00s
% Computational Cost: add. (6609->444), mult. (18354->669), div. (0->0), fcn. (15653->14), ass. (0->229)
t200 = cos(pkin(7));
t208 = cos(qJ(3));
t209 = cos(qJ(2));
t298 = t208 * t209;
t204 = sin(qJ(3));
t205 = sin(qJ(2));
t302 = t204 * t205;
t226 = -t200 * t302 + t298;
t284 = qJD(3) * t208;
t266 = t200 * t284;
t198 = sin(pkin(6));
t292 = qJD(1) * t198;
t342 = -pkin(2) * t266 + t226 * t292;
t197 = sin(pkin(7));
t241 = pkin(3) * t204 - qJ(4) * t208;
t215 = qJD(3) * t241 - qJD(4) * t204;
t265 = t205 * t292;
t341 = (t215 - t265) * t197;
t285 = qJD(3) * t204;
t269 = t197 * t285;
t340 = -pkin(9) * t269 + qJD(4) * t200 - t342;
t196 = sin(pkin(13));
t199 = cos(pkin(13));
t311 = t340 * t196 - t341 * t199;
t310 = t341 * t196 + t340 * t199;
t202 = sin(qJ(6));
t206 = cos(qJ(6));
t289 = qJD(2) * t200;
t188 = qJD(3) + t289;
t290 = qJD(2) * t197;
t271 = t204 * t290;
t137 = t188 * t199 - t196 * t271;
t138 = t188 * t196 + t199 * t271;
t203 = sin(qJ(5));
t207 = cos(qJ(5));
t91 = -t207 * t137 + t138 * t203;
t89 = qJD(6) + t91;
t336 = t89 ^ 2;
t164 = t196 * t207 + t199 * t203;
t305 = t197 * t208;
t220 = t164 * t305;
t216 = qJD(3) * t220;
t234 = t137 * t203 + t138 * t207;
t324 = qJD(5) * t234;
t56 = qJD(2) * t216 + t324;
t339 = -t202 * t336 + t206 * t56;
t303 = t199 * t208;
t222 = (pkin(4) * t204 - pkin(10) * t303) * t197;
t218 = qJD(3) * t222;
t338 = -t218 + t311;
t268 = t197 * t284;
t251 = t196 * t268;
t337 = -pkin(10) * t251 + t310;
t288 = qJD(2) * t208;
t270 = t197 * t288;
t183 = -qJD(5) + t270;
t335 = t91 * t183;
t299 = t207 * t199;
t163 = t196 * t203 - t299;
t219 = t163 * t305;
t125 = qJD(2) * t219;
t157 = t163 * qJD(5);
t295 = t125 - t157;
t300 = t205 * t208;
t301 = t204 * t209;
t228 = t200 * t300 + t301;
t267 = t200 * t285;
t296 = -pkin(2) * t267 - pkin(9) * t268 + t228 * t292;
t294 = -qJD(2) * t220 + t164 * qJD(5);
t273 = t209 * t292;
t201 = cos(pkin(6));
t291 = qJD(1) * t201;
t274 = t197 * t291;
t334 = qJD(2) * t273 + qJD(3) * t274;
t315 = qJD(2) * pkin(2);
t175 = t273 + t315;
t333 = t175 * t200 + t274;
t72 = t206 * t183 + t202 * t234;
t332 = t234 * t72;
t235 = t183 * t202 - t206 * t234;
t331 = t234 * t235;
t282 = qJD(5) * t207;
t283 = qJD(5) * t203;
t306 = t197 * t204;
t156 = t196 * t200 + t199 * t306;
t142 = pkin(9) * t305 + (pkin(2) * t204 + qJ(4)) * t200;
t242 = -pkin(3) * t208 - qJ(4) * t204;
t143 = (-pkin(2) + t242) * t197;
t96 = -t142 * t196 + t199 * t143;
t65 = -pkin(4) * t305 - pkin(10) * t156 + t96;
t154 = t196 * t306 - t199 * t200;
t97 = t199 * t142 + t196 * t143;
t76 = -pkin(10) * t154 + t97;
t330 = -t338 * t203 + t337 * t207 + t65 * t282 - t76 * t283;
t316 = t203 * t65 + t207 * t76;
t329 = t316 * qJD(5) + t337 * t203 + t338 * t207;
t328 = t183 * t234;
t327 = t204 * t208;
t318 = pkin(10) + qJ(4);
t181 = t318 * t196;
t182 = t318 * t199;
t233 = -t181 * t207 - t182 * t203;
t149 = t241 * t290;
t162 = pkin(9) * t290 + t265;
t152 = t204 * t162;
t98 = t208 * t333 - t152;
t66 = t199 * t149 - t196 * t98;
t45 = qJD(2) * t222 + t66;
t253 = t196 * t270;
t67 = t196 * t149 + t199 * t98;
t54 = -pkin(10) * t253 + t67;
t326 = -qJD(4) * t163 + qJD(5) * t233 - t203 * t45 - t207 * t54;
t123 = -t181 * t203 + t182 * t207;
t325 = qJD(4) * t164 + qJD(5) * t123 - t203 * t54 + t207 * t45;
t297 = pkin(4) * t251 - t296;
t187 = t200 * t291;
t111 = t187 + (qJD(2) * t242 - t175) * t197;
t99 = t208 * t162 + t333 * t204;
t84 = qJ(4) * t188 + t99;
t42 = t199 * t111 - t196 * t84;
t33 = -pkin(4) * t270 - pkin(10) * t138 + t42;
t43 = t196 * t111 + t199 * t84;
t36 = pkin(10) * t137 + t43;
t14 = t203 * t33 + t207 * t36;
t108 = (t215 + t265) * t290;
t254 = t200 * t265;
t232 = qJD(2) * t254;
t275 = t175 * t266 + t334 * t208;
t53 = t188 * qJD(4) + (-qJD(3) * t162 - t232) * t204 + t275;
t34 = t199 * t108 - t196 * t53;
t25 = qJD(2) * t218 + t34;
t278 = qJD(2) * qJD(3);
t264 = t197 * t278;
t247 = t208 * t264;
t231 = t196 * t247;
t35 = t196 * t108 + t199 * t53;
t29 = -pkin(10) * t231 + t35;
t213 = -qJD(5) * t14 - t203 * t29 + t207 * t25;
t246 = t204 * t264;
t4 = -pkin(5) * t246 - t213;
t323 = (pkin(5) * t234 + t89 * pkin(11)) * t89 + t4;
t192 = -pkin(4) * t199 - pkin(3);
t116 = pkin(5) * t163 - pkin(11) * t164 + t192;
t80 = pkin(4) * t253 + t99;
t13 = -t203 * t36 + t207 * t33;
t9 = pkin(5) * t183 - t13;
t322 = -t9 * qJD(6) * t164 - t116 * t56 + (-t294 * pkin(5) + t295 * pkin(11) + qJD(6) * t123 + t80) * t89;
t55 = t137 * t282 + t247 * t299 + (-qJD(5) * t138 - t231) * t203;
t22 = -qJD(6) * t235 + t202 * t55 - t206 * t246;
t60 = t162 * t284 + t175 * t267 + t334 * t204 + t208 * t232;
t48 = pkin(4) * t231 + t60;
t16 = pkin(5) * t56 - pkin(11) * t55 + t48;
t10 = -pkin(11) * t183 + t14;
t81 = -t188 * pkin(3) + qJD(4) - t98;
t57 = -t137 * pkin(4) + t81;
t19 = t91 * pkin(5) - pkin(11) * t234 + t57;
t240 = t10 * t202 - t19 * t206;
t225 = t203 * t25 + t207 * t29 + t33 * t282 - t36 * t283;
t3 = pkin(11) * t246 + t225;
t1 = -qJD(6) * t240 + t202 * t16 + t206 * t3;
t321 = -pkin(5) * t269 + t329;
t320 = t72 * t89;
t319 = t235 * t89;
t314 = t202 * t56;
t280 = qJD(6) * t206;
t281 = qJD(6) * t202;
t21 = -t183 * t280 + t202 * t246 + t206 * t55 - t234 * t281;
t313 = t21 * t202;
t312 = pkin(5) * t271 + t325;
t309 = t164 * t206;
t193 = t197 ^ 2;
t210 = qJD(2) ^ 2;
t307 = t193 * t210;
t304 = t198 * t210;
t293 = t204 ^ 2 - t208 ^ 2;
t287 = qJD(3) * t196;
t286 = qJD(3) * t199;
t279 = qJD(3) - t188;
t277 = t202 * t305;
t276 = t205 * t304;
t272 = t193 * t288;
t259 = t206 * t89;
t258 = t188 + t289;
t257 = 0.2e1 * t193 * t278;
t256 = t193 * t276;
t252 = t198 * t205 * t290;
t112 = t207 * t154 + t156 * t203;
t113 = -t154 * t203 + t156 * t207;
t145 = pkin(9) * t306 + (-pkin(2) * t208 - pkin(3)) * t200;
t114 = t154 * pkin(4) + t145;
t37 = t112 * pkin(5) - t113 * pkin(11) + t114;
t250 = -pkin(11) * t269 - qJD(6) * t37 - t330;
t27 = -pkin(11) * t305 + t316;
t69 = -qJD(3) * t219 - qJD(5) * t112;
t70 = qJD(5) * t113 + t216;
t245 = -pkin(5) * t70 + pkin(11) * t69 + qJD(6) * t27 - t297;
t6 = t10 * t206 + t19 * t202;
t238 = -t203 * t76 + t207 * t65;
t227 = t200 * t301 + t300;
t118 = t198 * t227 + t201 * t306;
t155 = -t197 * t198 * t209 + t200 * t201;
t85 = -t118 * t196 + t155 * t199;
t86 = t118 * t199 + t155 * t196;
t38 = t203 * t86 - t207 * t85;
t39 = t203 * t85 + t207 * t86;
t229 = t200 * t298 - t302;
t117 = -t198 * t229 - t201 * t305;
t237 = t117 * t206 - t202 * t39;
t236 = t117 * t202 + t206 * t39;
t87 = t113 * t202 + t206 * t305;
t110 = -t125 * t206 + t202 * t271;
t221 = -t157 * t206 - t164 * t281 - t110;
t217 = -pkin(11) * t56 + (t13 + t9) * t89;
t214 = -qJ(4) * t285 + (-pkin(3) * qJD(3) + qJD(4) - t81) * t208;
t2 = -qJD(6) * t6 + t206 * t16 - t202 * t3;
t211 = -t123 * t56 - t9 * t157 + t4 * t164 + (pkin(11) * t271 - qJD(6) * t116 - t326) * t89;
t132 = -t197 * t175 + t187;
t109 = -t125 * t202 - t206 * t271;
t88 = t113 * t206 - t277;
t79 = t201 * t268 + (qJD(2) * t226 + qJD(3) * t229) * t198;
t78 = t201 * t269 + (qJD(2) * t228 + qJD(3) * t227) * t198;
t63 = t196 * t252 + t199 * t79;
t62 = -t196 * t79 + t199 * t252;
t32 = -qJD(6) * t277 + t113 * t280 + t202 * t69 - t206 * t269;
t31 = -qJD(6) * t87 + t202 * t269 + t206 * t69;
t26 = pkin(5) * t305 - t238;
t12 = qJD(5) * t39 + t203 * t63 - t207 * t62;
t11 = -qJD(5) * t38 + t203 * t62 + t207 * t63;
t5 = [0, 0, -t276, -t209 * t304, 0, 0, 0, 0, 0, t155 * t246 - t188 * t78 - t208 * t256, t155 * t247 - t188 * t79 + t204 * t256, -t78 * t137 + (-t208 * t62 + (t117 * t196 * t208 + t204 * t85) * qJD(3)) * t290, t78 * t138 + (t208 * t63 + (t117 * t303 - t204 * t86) * qJD(3)) * t290, t63 * t137 - t62 * t138 + (-t196 * t86 - t199 * t85) * t247, t117 * t60 + t34 * t85 + t35 * t86 + t42 * t62 + t43 * t63 + t78 * t81, 0, 0, 0, 0, 0, t117 * t56 + t12 * t183 - t246 * t38 + t78 * t91, t11 * t183 + t117 * t55 + t234 * t78 - t246 * t39, 0, 0, 0, 0, 0 (-qJD(6) * t236 - t202 * t11 + t78 * t206) * t89 + t237 * t56 + t12 * t72 + t38 * t22 -(qJD(6) * t237 + t206 * t11 + t78 * t202) * t89 - t236 * t56 - t12 * t235 + t38 * t21; 0, 0, 0, 0, t257 * t327, -t293 * t257, t258 * t268, -t258 * t269, 0, -t60 * t200 + t296 * t188 + (t132 * t197 - t193 * t315) * t285 -(-t204 * t232 + t275) * t200 + t342 * t188 + (-pkin(2) * t272 + t200 * t152 + (pkin(9) * t188 * t204 + t132 * t208) * t197) * qJD(3), t60 * t154 + t296 * t137 + ((qJD(2) * t96 + t42) * t285 + (t81 * t287 - t34 + (t145 * t287 + t311) * qJD(2)) * t208) * t197, t60 * t156 - t296 * t138 + ((-qJD(2) * t97 - t43) * t285 + (t81 * t286 + t35 + (t145 * t286 + t310) * qJD(2)) * t208) * t197, -t35 * t154 - t34 * t156 + t311 * t138 + t310 * t137 + (-t196 * t43 - t199 * t42 + (-t196 * t97 - t199 * t96) * qJD(2)) * t268, t60 * t145 - t296 * t81 + t310 * t43 - t311 * t42 + t34 * t96 + t35 * t97, t113 * t55 + t234 * t69, -t112 * t55 - t113 * t56 - t234 * t70 - t69 * t91, -t69 * t183 + (-t208 * t55 + (qJD(2) * t113 + t234) * t285) * t197, t70 * t183 + (t208 * t56 + (-qJD(2) * t112 - t91) * t285) * t197 (-t183 * t197 - t272) * t285, t48 * t112 + t114 * t56 + t57 * t70 + t297 * t91 + t329 * t183 + (-t213 * t208 + (qJD(2) * t238 + t13) * t285) * t197, t48 * t113 + t114 * t55 + t57 * t69 + t297 * t234 + t330 * t183 + (t225 * t208 + (-qJD(2) * t316 - t14) * t285) * t197, t21 * t88 - t235 * t31, -t21 * t87 - t22 * t88 + t235 * t32 - t31 * t72, t112 * t21 - t235 * t70 + t31 * t89 + t56 * t88, -t112 * t22 - t32 * t89 - t56 * t87 - t70 * t72, t112 * t56 + t70 * t89 (-t202 * t27 + t206 * t37) * t56 + t2 * t112 - t240 * t70 + t26 * t22 + t4 * t87 + t9 * t32 + (t202 * t250 - t206 * t245) * t89 + t321 * t72 -(t202 * t37 + t206 * t27) * t56 - t1 * t112 - t6 * t70 + t26 * t21 + t4 * t88 + t9 * t31 + (t202 * t245 + t206 * t250) * t89 - t321 * t235; 0, 0, 0, 0, -t307 * t327, t293 * t307, t279 * t270, -t279 * t271, 0, -t132 * t271 + t188 * t99 - t60, t162 * t285 + t98 * t188 + (-t132 * t305 + t204 * t254) * qJD(2) - t275, t99 * t137 - t60 * t199 + (t196 * t214 - t204 * t42 + t208 * t66) * t290, -t99 * t138 + t60 * t196 + (t199 * t214 + t204 * t43 - t208 * t67) * t290, -t67 * t137 + t66 * t138 + (qJD(4) * t137 + t270 * t42 + t35) * t199 + (qJD(4) * t138 + t270 * t43 - t34) * t196, -t60 * pkin(3) - t42 * t66 - t43 * t67 - t81 * t99 + (-t196 * t42 + t199 * t43) * qJD(4) + (-t34 * t196 + t35 * t199) * qJ(4), t55 * t164 + t234 * t295, -t55 * t163 - t164 * t56 - t234 * t294 - t295 * t91, -t295 * t183 + (qJD(3) * t164 - t234) * t271, t294 * t183 + (-qJD(3) * t163 + t91) * t271, t183 * t271, t48 * t163 + t192 * t56 - t80 * t91 + t294 * t57 + t325 * t183 + (qJD(3) * t233 - t13) * t271, t48 * t164 + t192 * t55 - t80 * t234 + t295 * t57 + t326 * t183 + (-qJD(3) * t123 + t14) * t271, t21 * t309 - t221 * t235, -t235 * t109 + t110 * t72 - (t202 * t235 - t206 * t72) * t157 + (-t313 - t206 * t22 + (t202 * t72 + t206 * t235) * qJD(6)) * t164, t21 * t163 + t221 * t89 - t235 * t294 + t309 * t56, -t164 * t314 - t22 * t163 - t294 * t72 + (t157 * t202 - t164 * t280 + t109) * t89, t56 * t163 + t294 * t89, -t9 * t109 + t2 * t163 + t211 * t202 - t322 * t206 - t22 * t233 - t240 * t294 + t312 * t72, -t1 * t163 - t9 * t110 + t322 * t202 + t211 * t206 - t21 * t233 - t235 * t312 - t294 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t138 + t287) * t270 (-t137 + t286) * t270, -t137 ^ 2 - t138 ^ 2, -t137 * t43 + t138 * t42 + t60, 0, 0, 0, 0, 0, t56 - t328, t55 + t335, 0, 0, 0, 0, 0, -t332 + t339, -t206 * t336 - t314 + t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234 * t91, t234 ^ 2 - t91 ^ 2, t55 - t335, -t164 * t247 - t324 - t328, t246, -t14 * t183 - t234 * t57 + t213, -t13 * t183 + t57 * t91 - t225, -t235 * t259 + t313 (t21 - t320) * t206 + (-t22 + t319) * t202, t259 * t89 + t314 + t331, t332 + t339, -t89 * t234, -pkin(5) * t22 - t14 * t72 + t217 * t202 - t323 * t206 + t234 * t240, -pkin(5) * t21 + t14 * t235 + t323 * t202 + t217 * t206 + t6 * t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235 * t72, t235 ^ 2 - t72 ^ 2, t21 + t320, -t22 - t319, t56, t235 * t9 + t6 * t89 + t2, -t240 * t89 + t9 * t72 - t1;];
tauc_reg  = t5;
