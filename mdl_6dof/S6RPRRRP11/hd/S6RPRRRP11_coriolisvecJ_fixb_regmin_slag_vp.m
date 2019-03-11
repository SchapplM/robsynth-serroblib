% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP11_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:13
% EndTime: 2019-03-09 06:41:31
% DurationCPUTime: 7.17s
% Computational Cost: add. (14481->480), mult. (47814->691), div. (0->0), fcn. (40859->12), ass. (0->231)
t199 = sin(pkin(12));
t205 = sin(qJ(3));
t201 = sin(pkin(6));
t289 = qJD(1) * t201;
t273 = t205 * t289;
t202 = cos(pkin(12));
t312 = cos(pkin(7));
t332 = cos(qJ(3));
t252 = t312 * t332;
t236 = t252 * t289;
t200 = sin(pkin(7));
t313 = cos(pkin(6));
t265 = t313 * t200;
t245 = t332 * t265;
t293 = -qJD(1) * t245 - t202 * t236;
t139 = t199 * t273 + t293;
t234 = qJD(4) + t139;
t207 = cos(qJ(4));
t266 = t205 * t312;
t150 = (t332 * t199 + t202 * t266) * t201 + t205 * t265;
t142 = t150 * qJD(3);
t349 = qJD(1) * t142;
t352 = t207 * t349;
t233 = t199 * t266 - t332 * t202;
t226 = t201 * t233;
t160 = qJD(1) * t226;
t270 = qJD(3) * t332;
t351 = t200 * t270 + t160;
t204 = sin(qJ(4));
t261 = qJD(1) * t313;
t256 = pkin(1) * t261;
t189 = t202 * t256;
t307 = t199 * t201;
t217 = t313 * pkin(2) + (-t312 * pkin(9) - qJ(2)) * t307;
t138 = qJD(1) * t217 + t189;
t264 = t312 * t138;
t126 = t205 * t264;
t274 = t202 * t289;
t165 = qJ(2) * t274 + t199 * t256;
t304 = t201 * t202;
t224 = (t312 * t304 + t265) * pkin(9);
t134 = qJD(1) * t224 + t165;
t128 = t332 * t134;
t308 = t199 * t200;
t162 = (-pkin(2) * t202 - pkin(9) * t308 - pkin(1)) * t201;
t154 = qJD(1) * t162 + qJD(2);
t305 = t200 * t205;
t89 = t154 * t305 + t126 + t128;
t350 = -t89 + t234 * (pkin(4) * t204 - pkin(11) * t207);
t169 = t204 * t305 - t207 * t312;
t275 = t199 * t289;
t259 = t200 * t275;
t348 = -t169 * qJD(4) - t204 * t259 + t351 * t207;
t231 = t199 * t252 + t202 * t205;
t225 = t201 * t231;
t272 = qJD(3) * t305;
t347 = qJD(1) * t225 - t272;
t143 = qJD(1) * t150;
t249 = t313 * t312;
t181 = t200 * t274;
t281 = qJD(3) - t181;
t230 = -qJD(1) * t249 - t281;
t118 = t204 * t143 + t207 * t230;
t133 = t139 * qJD(3);
t299 = t207 * t133;
t212 = -qJD(4) * t118 - t299;
t346 = qJD(5) * t234 + t212;
t330 = -qJ(6) - pkin(11);
t345 = -qJ(6) * t118 + qJD(5) * t330;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t287 = qJD(4) * t204;
t280 = pkin(10) * t287;
t344 = t203 * t280 + t206 * t350;
t343 = t201 ^ 2 * (t199 ^ 2 + t202 ^ 2);
t276 = t200 * t332;
t342 = -t205 * t134 + t138 * t252 + t154 * t276;
t170 = t204 * t312 + t207 * t305;
t296 = t170 * qJD(4) + t351 * t204 + t207 * t259;
t228 = qJD(4) * t234;
t340 = -t204 * t349 - t207 * t228;
t277 = pkin(1) * t313;
t292 = qJ(2) * t304 + t199 * t277;
t146 = t224 + t292;
t193 = t202 * t277;
t151 = t193 + t217;
t339 = -t205 * t146 + t151 * t252 + t162 * t276;
t285 = qJD(4) * t207;
t309 = t139 * t207;
t338 = t285 + t309;
t302 = t204 * t139;
t337 = t287 + t302;
t183 = -pkin(4) * t207 - pkin(11) * t204 - pkin(3);
t283 = qJD(5) * t206;
t110 = pkin(3) * t143 + pkin(10) * t139;
t314 = t204 * t110 + t207 * t342;
t48 = pkin(11) * t143 + t314;
t336 = -t183 * t283 - t203 * t350 + t206 * t48;
t335 = (t154 * t200 + t264) * t205 + t128;
t120 = t207 * t143 - t204 * t230;
t95 = t206 * t120 + t203 * t234;
t334 = t95 ^ 2;
t117 = qJD(5) + t118;
t115 = -t138 * t200 + t312 * t154;
t65 = pkin(3) * t139 - pkin(10) * t143 + t115;
t69 = -pkin(10) * t230 + t89;
t36 = t204 * t65 + t207 * t69;
t30 = pkin(11) * t234 + t36;
t68 = pkin(3) * t230 - t342;
t42 = t118 * pkin(4) - t120 * pkin(11) + t68;
t13 = -t203 * t30 + t206 * t42;
t9 = -qJ(6) * t95 + t13;
t8 = pkin(5) * t117 + t9;
t333 = t8 - t9;
t331 = pkin(5) * t203;
t35 = -t204 * t69 + t207 * t65;
t78 = pkin(4) * t120 + pkin(11) * t118;
t329 = t203 * t78 + t206 * t35;
t244 = t202 * t252;
t306 = t199 * t205;
t279 = t201 * t306;
t149 = -t201 * t244 - t245 + t279;
t121 = -t151 * t200 + t312 * t162;
t80 = pkin(3) * t149 - pkin(10) * t150 + t121;
t168 = t200 * t304 - t249;
t216 = t332 * t146 + (t312 * t151 + t162 * t200) * t205;
t87 = -t168 * pkin(10) + t216;
t326 = t204 * t80 + t207 * t87;
t39 = pkin(11) * t149 + t326;
t122 = t150 * t204 + t168 * t207;
t123 = t150 * t207 - t168 * t204;
t86 = t168 * pkin(3) - t339;
t51 = t122 * pkin(4) - t123 * pkin(11) + t86;
t328 = t203 * t51 + t206 * t39;
t93 = t120 * t203 - t206 * t234;
t325 = t117 * t93;
t288 = qJD(2) * t201;
t257 = t288 * t308;
t246 = qJD(1) * t257;
t101 = pkin(3) * t349 + t133 * pkin(10) + t246;
t219 = qJD(2) * t226;
t59 = -qJD(1) * t219 + qJD(3) * t342;
t263 = -t207 * t101 + t204 * t59 + t69 * t285 + t65 * t287;
t17 = -pkin(4) * t349 + t263;
t324 = t17 * t203;
t323 = t17 * t206;
t284 = qJD(5) * t203;
t43 = t120 * t284 - t203 * t349 - t346 * t206;
t322 = t203 * t43;
t303 = t204 * t133;
t91 = t120 * qJD(4) - t303;
t321 = t203 * t91;
t320 = t206 * t91;
t319 = t95 * t117;
t300 = t206 * t207;
t107 = -t139 * t300 + t143 * t203;
t194 = pkin(10) * t300;
t282 = qJD(6) * t206;
t310 = qJ(6) * t204;
t318 = pkin(5) * t302 + qJ(6) * t107 + t203 * t48 - t204 * t282 + (pkin(5) * t204 - qJ(6) * t300) * qJD(4) + (-t194 + (-t183 + t310) * t203) * qJD(5) + t344;
t106 = -t206 * t143 - t203 * t309;
t301 = t204 * t206;
t317 = qJ(6) * t106 + (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t301 + (-qJD(6) * t204 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t207) * t203 - t336;
t316 = t345 * t203 + t282 - t329;
t73 = t206 * t78;
t315 = -pkin(5) * t120 - t73 + t345 * t206 + (-qJD(6) + t35) * t203;
t235 = -t206 * t170 + t203 * t276;
t298 = -qJD(5) * t235 + t348 * t203 + t347 * t206;
t155 = -t203 * t170 - t206 * t276;
t297 = -qJD(5) * t155 + t347 * t203 - t348 * t206;
t291 = t203 * t183 + t194;
t286 = qJD(4) * t206;
t278 = pkin(10) + t331;
t271 = t117 * t284;
t269 = -t203 * t39 + t206 * t51;
t268 = -t204 * t87 + t207 * t80;
t262 = t110 * t207 - t204 * t342;
t260 = t117 * t206;
t208 = qJD(1) ^ 2;
t254 = t201 * t208 * t313;
t250 = t206 * t285 - t107;
t14 = t203 * t42 + t206 * t30;
t103 = t123 * t206 + t149 * t203;
t102 = t123 * t203 - t149 * t206;
t248 = (-qJ(2) * t275 + t189) * t199 - t165 * t202;
t141 = (t245 + (t244 - t306) * t201) * qJD(3);
t105 = pkin(3) * t142 - pkin(10) * t141 + t257;
t74 = t339 * qJD(3) - t219;
t247 = t105 * t207 - t204 * t74 - t87 * t285 - t80 * t287;
t38 = -pkin(4) * t149 - t268;
t47 = -pkin(4) * t143 - t262;
t243 = -t117 * t283 - t321;
t242 = -t204 * t101 - t207 * t59 - t65 * t285 + t69 * t287;
t16 = pkin(11) * t349 - t242;
t28 = t91 * pkin(4) - t212 * pkin(11) + qJD(3) * t126 + t134 * t270 + t154 * t272 + (t199 * t236 + t202 * t273) * qJD(2);
t5 = t206 * t16 + t203 * t28 + t42 * t283 - t30 * t284;
t238 = t204 * t105 + t207 * t74 + t80 * t285 - t87 * t287;
t22 = pkin(11) * t142 + t238;
t100 = -t122 * qJD(4) + t141 * t207;
t218 = qJD(2) * t225;
t75 = qJD(3) * t216 + t218;
t99 = t123 * qJD(4) + t141 * t204;
t33 = t99 * pkin(4) - t100 * pkin(11) + t75;
t241 = t203 * t33 + t206 * t22 + t51 * t283 - t39 * t284;
t240 = -0.2e1 * t261 * t288;
t29 = -pkin(4) * t234 - t35;
t239 = -pkin(11) * t91 + t117 * t29;
t23 = -pkin(4) * t142 - t247;
t229 = t139 * t234;
t6 = -t14 * qJD(5) - t16 * t203 + t206 * t28;
t223 = -t328 * qJD(5) - t203 * t22 + t206 * t33;
t44 = t120 * t283 + t346 * t203 - t206 * t349;
t7 = t44 * pkin(5) + t17;
t185 = t330 * t206;
t184 = t330 * t203;
t177 = t206 * t183;
t157 = -t203 * t310 + t291;
t148 = -qJ(6) * t301 + t177 + (-pkin(10) * t203 - pkin(5)) * t207;
t92 = t93 ^ 2;
t60 = qJD(1) * t218 + t335 * qJD(3);
t53 = -t102 * qJD(5) + t100 * t206 + t142 * t203;
t52 = t103 * qJD(5) + t100 * t203 - t142 * t206;
t25 = t93 * pkin(5) + qJD(6) + t29;
t12 = -qJ(6) * t102 + t328;
t11 = pkin(5) * t122 - qJ(6) * t103 + t269;
t10 = -qJ(6) * t93 + t14;
t4 = -qJ(6) * t52 - qJD(6) * t102 + t241;
t3 = pkin(5) * t99 - qJ(6) * t53 - qJD(6) * t103 + t223;
t2 = -qJ(6) * t44 - qJD(6) * t93 + t5;
t1 = pkin(5) * t91 + qJ(6) * t43 - qJD(6) * t95 + t6;
t15 = [0, 0, 0, t199 * t240, t202 * t240, 0.2e1 * qJD(2) * qJD(1) * t343 ((t202 * t292 + (qJ(2) * t307 - t193) * t199) * qJD(1) - t248) * t288, -t133 * t150 + t141 * t143, t133 * t149 - t141 * t139 - t143 * t142 - t150 * t349, t133 * t168 - t141 * t230 (-t281 + (t168 - t249) * qJD(1)) * t142, 0, t115 * t142 + t121 * t349 + t139 * t257 + t149 * t246 + t60 * t168 + t230 * t75, t115 * t141 - t121 * t133 + 0.2e1 * t143 * t257 + t59 * t168 + t74 * t230, t120 * t100 + t123 * t212, -t100 * t118 - t120 * t99 - t122 * t212 - t123 * t91, t100 * t234 + t120 * t142 + t123 * t349 + t149 * t212, -t118 * t142 - t122 * t349 - t91 * t149 - t234 * t99 (qJD(4) + t293 + (t149 + t279) * qJD(1)) * t142, t75 * t118 + t60 * t122 + t35 * t142 - t263 * t149 + t247 * t234 + t268 * t349 + t68 * t99 + t86 * t91, t68 * t100 + t75 * t120 + t60 * t123 - t36 * t142 + t242 * t149 + t86 * t212 - t238 * t234 - t326 * t349, -t103 * t43 + t53 * t95, t102 * t43 - t103 * t44 - t52 * t95 - t53 * t93, t103 * t91 + t117 * t53 - t122 * t43 + t95 * t99, -t102 * t91 - t117 * t52 - t122 * t44 - t93 * t99, t117 * t99 + t122 * t91, t17 * t102 + t117 * t223 + t6 * t122 + t13 * t99 + t23 * t93 + t269 * t91 + t29 * t52 + t38 * t44, t17 * t103 - t241 * t117 - t5 * t122 - t14 * t99 + t23 * t95 + t29 * t53 - t328 * t91 - t38 * t43, -t1 * t103 - t10 * t52 - t102 * t2 + t11 * t43 - t12 * t44 - t3 * t95 - t4 * t93 - t53 * t8, t2 * t12 + t10 * t4 + t1 * t11 + t8 * t3 + t7 * (pkin(5) * t102 + t38) + t25 * (pkin(5) * t52 + t23); 0, 0, 0, t199 * t254, t202 * t254, -t208 * t343, t248 * t289, 0, 0, 0, 0, 0, -t139 * t259 - t230 * t347 + t312 * t349, -t312 * t133 + t160 * t230 + (-t143 * t275 + t230 * t270) * t200, 0, 0, 0, 0, 0, -t118 * t347 - t169 * t349 - t296 * t234 - t91 * t276, -t120 * t347 - t170 * t349 - t212 * t276 - t348 * t234, 0, 0, 0, 0, 0, -t117 * t298 + t155 * t91 + t169 * t44 + t296 * t93, t117 * t297 - t169 * t43 + t235 * t91 + t296 * t95, t155 * t43 + t235 * t44 + t297 * t93 + t298 * t95, t1 * t155 - t10 * t297 + t169 * t7 - t2 * t235 + t25 * t296 - t298 * t8; 0, 0, 0, 0, 0, 0, 0, t143 * t139, -t139 ^ 2 + t143 ^ 2, -t139 * t230 - t133, t143 * t281 + (t143 * t249 - t142) * qJD(1), 0, -t115 * t143 - t89 * t181 + (-t231 * t288 + t89 * t249) * qJD(1) + (-t335 + t89) * qJD(3), t115 * t139 - t342 * t181 + (t233 * t288 + t249 * t342) * qJD(1), -qJD(4) * t204 ^ 2 * t143 + ((-qJD(4) * t230 - t133) * t204 + t234 * t120) * t207, -t338 * t118 - t337 * t120 - t204 * t91 + t207 * t212, -t120 * t143 + t207 * t229 - t340, t118 * t143 + t352 + (-t228 - t229) * t204, -t234 * t143, -pkin(3) * t91 + t340 * pkin(10) - t89 * t118 - t35 * t143 - t60 * t207 - t234 * t262 + t337 * t68, -pkin(10) * t352 - pkin(3) * t212 - t89 * t120 + t36 * t143 + t60 * t204 + t338 * t68 + (t280 + t314) * t234, -t43 * t301 + (-t204 * t284 + t250) * t95, t106 * t95 + t107 * t93 + (-t203 * t95 - t206 * t93) * t285 + (t322 - t206 * t44 + (t203 * t93 - t206 * t95) * qJD(5)) * t204, t207 * t43 + t250 * t117 + (t234 * t95 - t271 + t320) * t204, t207 * t44 + (-t203 * t285 + t106) * t117 + (-t234 * t93 + t243) * t204, t117 * t204 * t234 - t207 * t91, -t29 * t106 + t177 * t91 - t47 * t93 + ((-qJD(5) * t183 + t48) * t203 + t344) * t117 + (t29 * t203 * qJD(4) - t6 + (qJD(4) * t93 + t243) * pkin(10)) * t207 + (pkin(10) * t44 + t13 * t234 + t283 * t29 + t324) * t204, -t291 * t91 - t47 * t95 - t29 * t107 + t336 * t117 + (t29 * t286 + t5 + (qJD(4) * t95 + t271) * pkin(10)) * t207 + (-t29 * t284 + t323 - t234 * t14 + (t117 * t286 - t43) * pkin(10)) * t204, t10 * t106 + t107 * t8 + t148 * t43 - t157 * t44 - t318 * t95 - t317 * t93 + (-t10 * t203 - t206 * t8) * t285 + (-t1 * t206 - t2 * t203 + (-t10 * t206 + t203 * t8) * qJD(5)) * t204, t278 * t7 * t204 + t1 * t148 + t317 * t10 + t2 * t157 + t318 * t8 + (t278 * t285 - t47 + (t283 * t204 - t106) * pkin(5)) * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120 * t118, -t118 ^ 2 + t120 ^ 2, t139 * t118 - t299, t120 * t139 + t303, t349, -t68 * t120 + t234 * t36 - t263, t68 * t118 + t234 * t35 + t242, t260 * t95 - t322 (-t43 - t325) * t206 + (-t44 - t319) * t203, t117 * t260 - t120 * t95 + t321, -t117 ^ 2 * t203 + t120 * t93 + t320, -t117 * t120, -pkin(4) * t44 - t13 * t120 - t323 - t36 * t93 + (-pkin(11) * t283 - t73) * t117 + (t35 * t117 + t239) * t203, pkin(4) * t43 + t14 * t120 + t324 - t36 * t95 + (pkin(11) * t284 + t329) * t117 + t239 * t206, t184 * t43 + t185 * t44 - t315 * t95 - t316 * t93 + (-t117 * t8 + t2) * t206 + (-t10 * t117 - t1) * t203, -t2 * t185 + t1 * t184 + t7 * (-pkin(5) * t206 - pkin(4)) + t315 * t8 + (t117 * t331 - t36) * t25 + t316 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t93, -t92 + t334, -t43 + t325, -t44 + t319, t91, t117 * t14 - t29 * t95 + t6, t117 * t13 + t29 * t93 - t5, pkin(5) * t43 - t333 * t93, t333 * t10 + (-t25 * t95 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92 - t334, t10 * t93 + t8 * t95 + t7;];
tauc_reg  = t15;
