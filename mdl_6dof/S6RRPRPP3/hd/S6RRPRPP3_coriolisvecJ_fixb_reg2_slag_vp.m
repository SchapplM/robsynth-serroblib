% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:09
% EndTime: 2019-03-09 09:57:21
% DurationCPUTime: 4.75s
% Computational Cost: add. (6702->470), mult. (16916->588), div. (0->0), fcn. (11940->6), ass. (0->237)
t201 = sin(pkin(9));
t205 = sin(qJ(2));
t291 = qJD(1) * t205;
t274 = t201 * t291;
t202 = cos(pkin(9));
t284 = t202 * qJD(2);
t150 = -t274 + t284;
t273 = t202 * t291;
t285 = t201 * qJD(2);
t151 = t273 + t285;
t204 = sin(qJ(4));
t332 = cos(qJ(4));
t100 = -t332 * t150 + t204 * t151;
t206 = cos(qJ(2));
t283 = t206 * qJD(1);
t187 = -qJD(4) + t283;
t343 = t187 * t100;
t279 = t332 * t202;
t256 = t206 * t279;
t240 = qJD(1) * t256;
t282 = qJD(1) * qJD(2);
t269 = t206 * t282;
t254 = t201 * t269;
t271 = qJD(4) * t332;
t60 = t204 * (qJD(4) * t151 + t254) - qJD(2) * t240 - t150 * t271;
t212 = t60 - t343;
t344 = t60 + t343;
t278 = t201 * t283;
t286 = qJD(4) * t204;
t336 = -t201 * t286 + t202 * t271;
t296 = -t204 * t278 + t240 - t336;
t97 = t100 ^ 2;
t246 = pkin(2) * t205 - qJ(3) * t206;
t157 = t246 * qJD(1);
t141 = t201 * t157;
t306 = t202 * t205;
t307 = t201 * t206;
t226 = -pkin(7) * t306 - pkin(8) * t307;
t105 = qJD(1) * t226 + t141;
t329 = pkin(8) + qJ(3);
t166 = t329 * t201;
t167 = t329 * t202;
t118 = pkin(7) * t274 + t202 * t157;
t305 = t202 * t206;
t235 = pkin(3) * t205 - pkin(8) * t305;
t88 = qJD(1) * t235 + t118;
t322 = -qJD(3) * t279 + t332 * t105 + t166 * t271 + (qJD(3) * t201 + qJD(4) * t167 + t88) * t204;
t155 = t332 * t201 + t204 * t202;
t140 = t155 * qJD(4);
t222 = t206 * t155;
t297 = qJD(1) * t222 - t140;
t229 = t204 * t150 + t332 * t151;
t162 = -t206 * pkin(2) - t205 * qJ(3) - pkin(1);
t143 = t162 * qJD(1);
t193 = pkin(7) * t283;
t169 = qJD(2) * qJ(3) + t193;
t108 = t202 * t143 - t201 * t169;
t219 = -pkin(3) * t283 - t151 * pkin(8) + t108;
t65 = t332 * t219;
t109 = t201 * t143 + t202 * t169;
t74 = t150 * pkin(8) + t109;
t35 = t204 * t74 - t65;
t233 = pkin(5) * t229 + t35;
t299 = qJD(5) + t233;
t98 = t229 ^ 2;
t342 = -0.2e1 * t282;
t326 = qJ(5) * t291 + t322;
t272 = -t98 + t97;
t116 = -t204 * t166 + t332 * t167;
t341 = -qJD(3) * t155 - qJD(4) * t116 + t204 * t105 - t332 * t88;
t181 = t187 ^ 2;
t340 = -t98 - t181;
t339 = t229 * t187;
t172 = qJD(5) * t187;
t191 = t205 * t282;
t182 = qJ(5) * t191;
t338 = t182 - t172;
t144 = pkin(3) * t278 + t193;
t337 = -t296 * qJ(5) + t155 * qJD(5) + t144;
t335 = -t204 * t201 + t279;
t199 = t205 ^ 2;
t200 = t206 ^ 2;
t334 = qJD(1) * (t199 - 0.2e1 * t200);
t333 = qJD(4) * t229;
t132 = t155 * t205;
t288 = qJD(2) * t205;
t220 = qJD(2) * t222;
t61 = qJD(1) * t220 + t333;
t83 = t336 * t205 + t220;
t14 = -t187 * t83 - t206 * t61 + (qJD(1) * t132 + t100) * t288;
t133 = t335 * t205;
t277 = t206 * t285;
t82 = -qJD(2) * t256 + t140 * t205 + t204 * t277;
t15 = t187 * t82 + t206 * t60 + (qJD(1) * t133 + t229) * t288;
t32 = t296 * t187 + (qJD(2) * t155 - t229) * t291;
t31 = t297 * t187 + (-qJD(2) * t335 - t100) * t291;
t213 = t296 * t100 - t61 * t155 + t229 * t297 - t335 * t60;
t331 = t100 * pkin(5);
t330 = pkin(4) + qJ(6);
t270 = t205 * t330;
t328 = -t296 * pkin(5) + qJD(1) * t270 - t341;
t327 = -t297 * pkin(5) + t326;
t325 = -pkin(4) * t291 + t341;
t324 = -qJD(6) * t335 - t297 * t330 - t337;
t321 = t297 * pkin(4) + t337;
t320 = qJ(5) * t61;
t319 = qJD(2) * pkin(2);
t192 = pkin(7) * t291;
t265 = qJD(3) - t319;
t161 = t192 + t265;
t117 = -t150 * pkin(3) + t161;
t218 = -qJ(5) * t229 + t117;
t40 = t100 * pkin(4) + t218;
t318 = t229 * t40;
t26 = t330 * t100 + t218;
t317 = t26 * t229;
t149 = t202 * t162;
t107 = -pkin(8) * t306 + t149 + (-pkin(7) * t201 - pkin(3)) * t206;
t185 = pkin(7) * t305;
t123 = t201 * t162 + t185;
t308 = t201 * t205;
t114 = -pkin(8) * t308 + t123;
t53 = t204 * t107 + t332 * t114;
t315 = t100 * qJ(5);
t314 = t100 * t229;
t312 = t150 * t202;
t208 = qJD(1) ^ 2;
t309 = t200 * t208;
t303 = t206 * t208;
t207 = qJD(2) ^ 2;
t302 = t207 * t205;
t301 = t207 * t206;
t300 = -qJD(5) - t35;
t216 = t204 * t219;
t36 = t332 * t74 + t216;
t25 = t36 - t331;
t298 = -qJD(6) - t25;
t136 = qJD(2) * t246 - t205 * qJD(3);
t126 = t136 * qJD(1);
t160 = (qJD(3) - t192) * qJD(2);
t87 = t201 * t126 + t202 * t160;
t280 = pkin(7) * t288;
t112 = t202 * t136 + t201 * t280;
t186 = pkin(7) * t269;
t135 = pkin(3) * t254 + t186;
t287 = qJD(2) * t206;
t194 = pkin(7) * t287;
t145 = pkin(3) * t277 + t194;
t158 = pkin(3) * t308 + t205 * pkin(7);
t292 = t199 - t200;
t115 = t332 * t166 + t204 * t167;
t290 = qJD(2) * t115;
t289 = qJD(2) * t116;
t281 = pkin(7) * t307;
t190 = -t202 * pkin(3) - pkin(2);
t275 = t187 * t291;
t224 = t235 * qJD(2);
t86 = t202 * t126 - t201 * t160;
t64 = qJD(1) * t224 + t86;
t73 = -pkin(8) * t254 + t87;
t268 = -qJD(4) * t65 - t204 * t64 + t74 * t286 - t332 * t73;
t267 = qJD(4) * t216 + t204 * t73 + t74 * t271 - t332 * t64;
t266 = pkin(1) * t342;
t264 = -t191 + t314;
t259 = -t150 + t284;
t258 = -t151 + t285;
t257 = pkin(4) * t191;
t253 = t206 * t191;
t252 = qJD(2) * t270;
t250 = -t172 - t268;
t48 = t206 * qJ(5) - t53;
t249 = -t161 + t265;
t247 = -t133 * qJ(5) + t158;
t52 = t332 * t107 - t204 * t114;
t22 = t100 * t83 + t61 * t132;
t21 = -t60 * t133 - t229 * t82;
t245 = -t115 * t60 - t116 * t61;
t241 = -t98 - t97;
t238 = qJD(1) * t259;
t237 = qJD(1) * t258;
t234 = -t155 * qJ(5) + t190;
t232 = -t61 * pkin(5) - t268;
t231 = -t60 * pkin(5) + t267;
t49 = t206 * pkin(4) - t52;
t230 = -t187 * t36 - t267;
t78 = t224 + t112;
t128 = t201 * t136;
t89 = qJD(2) * t226 + t128;
t228 = t107 * t286 + t114 * t271 + t204 * t89 - t332 * t78;
t19 = t107 * t271 - t114 * t286 + t204 * t78 + t332 * t89;
t227 = t206 * t237;
t11 = -t60 * t155 - t229 * t296;
t12 = -t297 * t100 - t335 * t61;
t225 = t82 * qJ(5) - t133 * qJD(5) + t145;
t5 = -t257 + t267;
t221 = -t100 * t26 + t232;
t13 = t61 * pkin(4) + t60 * qJ(5) - qJD(5) * t229 + t135;
t217 = t82 * t100 + t60 * t132 - t61 * t133 - t229 * t83;
t3 = t61 * qJ(6) + t100 * qJD(6) + t13;
t214 = -qJD(1) * t252 + t231;
t16 = -qJ(5) * t288 + t206 * qJD(5) - t19;
t29 = t187 * qJ(5) - t36;
t211 = -t155 * t269 - t333;
t210 = t61 - t339;
t209 = t211 - t339;
t198 = t202 ^ 2;
t197 = t201 ^ 2;
t183 = t205 * t303;
t177 = 0.2e1 * t182;
t170 = -0.2e1 * t253;
t127 = (-t187 - t283) * t288;
t122 = t149 - t281;
t119 = -pkin(7) * t273 + t141;
t113 = -t202 * t280 + t128;
t96 = -pkin(4) * t335 + t234;
t81 = pkin(5) * t335 + t116;
t80 = t155 * pkin(5) + t115;
t72 = -t330 * t335 + t234;
t68 = t132 * pkin(4) + t247;
t51 = pkin(4) * t229 + t315;
t50 = t330 * t132 + t247;
t41 = -t132 * pkin(5) - t48;
t39 = t133 * pkin(5) + t206 * qJ(6) + t49;
t38 = t61 + t339;
t34 = t229 * t330 + t315;
t28 = t187 * pkin(4) - t300;
t27 = t83 * pkin(4) + t225;
t23 = qJD(6) - t29 - t331;
t18 = t330 * t187 + t299;
t17 = -pkin(4) * t288 + t228;
t10 = t132 * qJD(6) + t330 * t83 + t225;
t9 = -t83 * pkin(5) - t16;
t8 = -t82 * pkin(5) + t206 * qJD(6) + t228 - t252;
t4 = -t182 - t250;
t2 = t232 + t338;
t1 = qJD(6) * t187 + t214;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t253, t292 * t342, t301, t170, -t302, 0, -pkin(7) * t301 + t205 * t266, pkin(7) * t302 + t206 * t266, 0, 0 (t151 * t202 + t198 * t291) * t287 (t312 + (-t151 - 0.2e1 * t273) * t201) * t287 (t151 * t205 + t202 * t334) * qJD(2) (-t150 * t201 + t197 * t291) * t287 (t150 * t205 - t201 * t334) * qJD(2), t170 (-qJD(1) * t112 - t86) * t206 + ((-pkin(7) * t150 + t161 * t201) * t206 + (t108 + (t122 + 0.2e1 * t281) * qJD(1)) * t205) * qJD(2) (qJD(1) * t113 + t87) * t206 + ((pkin(7) * t151 + t161 * t202) * t206 + (-t109 + (-t123 + 0.2e1 * t185) * qJD(1)) * t205) * qJD(2), -t112 * t151 + t113 * t150 + (-t201 * t87 - t202 * t86) * t205 + (-t108 * t202 - t109 * t201 + (-t122 * t202 - t123 * t201) * qJD(1)) * t287, t108 * t112 + t109 * t113 + t86 * t122 + t87 * t123 + (t161 + t192) * t194, t21, t217, t15, t22, -t14, t127, t145 * t100 + t117 * t83 + t135 * t132 + t158 * t61 + t228 * t187 + t267 * t206 + (qJD(1) * t52 - t35) * t288, t145 * t229 - t117 * t82 + t135 * t133 - t158 * t60 + t19 * t187 - t268 * t206 + (-qJD(1) * t53 - t36) * t288, -t19 * t100 + t132 * t268 + t133 * t267 + t228 * t229 - t35 * t82 - t36 * t83 + t52 * t60 - t53 * t61, t117 * t145 + t135 * t158 + t36 * t19 + t228 * t35 - t267 * t52 - t268 * t53, t127, -t15, t14, t21, t217, t22, t16 * t100 + t4 * t132 + t5 * t133 + t17 * t229 - t28 * t82 + t29 * t83 + t48 * t61 - t49 * t60, -t27 * t100 - t13 * t132 - t17 * t187 - t5 * t206 - t40 * t83 - t68 * t61 + (qJD(1) * t49 + t28) * t288, -t27 * t229 - t13 * t133 + t16 * t187 + t4 * t206 + t40 * t82 + t68 * t60 + (-qJD(1) * t48 - t29) * t288, t13 * t68 + t29 * t16 + t28 * t17 + t40 * t27 + t4 * t48 + t5 * t49, t127, t14, t15, t22, -t217, t21, t1 * t133 - t100 * t9 - t132 * t2 - t18 * t82 + t229 * t8 - t23 * t83 - t39 * t60 - t41 * t61, -t10 * t229 - t3 * t133 - t9 * t187 - t2 * t206 + t26 * t82 + t50 * t60 + (qJD(1) * t41 + t23) * t288, t1 * t206 + t10 * t100 + t3 * t132 + t8 * t187 + t26 * t83 + t50 * t61 + (-qJD(1) * t39 - t18) * t288, t1 * t39 + t10 * t26 + t18 * t8 + t2 * t41 + t23 * t9 + t3 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183, t292 * t208, 0, t183, 0, 0, t208 * pkin(1) * t205, pkin(1) * t303, 0, 0, t202 * t227 (-t312 + t151 * t201 + (-t197 + t198) * qJD(2)) * t283, t202 * t309 + t205 * t237, -t259 * t278, -t201 * t309 + t205 * t238, t183 ((-qJ(3) * t285 - t108) * t205 + (-pkin(7) * t259 + t201 * t249 + t118) * t206) * qJD(1) ((-qJ(3) * t284 + t109) * t205 + (pkin(7) * t258 + t202 * t249 - t119) * t206) * qJD(1), t118 * t151 - t119 * t150 + (qJD(3) * t150 + t108 * t283 + t87) * t202 + (qJD(3) * t151 + t109 * t283 - t86) * t201, -t108 * t118 - t109 * t119 + (-t108 * t201 + t109 * t202) * qJD(3) + (-t86 * t201 + t87 * t202) * qJ(3) + (-t161 - t319) * t193, t11, t213, t32, t12, -t31, t275, -t144 * t100 - t135 * t335 + t190 * t61 - t341 * t187 - t297 * t117 + (t35 - t290) * t291, -t144 * t229 + t135 * t155 - t190 * t60 - t322 * t187 - t296 * t117 + (t36 - t289) * t291, t100 * t322 + t155 * t267 - t229 * t341 - t268 * t335 - t296 * t35 + t297 * t36 + t245, t115 * t267 - t116 * t268 - t117 * t144 + t135 * t190 - t322 * t36 - t341 * t35, t275, -t32, t31, t11, t213, t12, t100 * t326 + t5 * t155 - t229 * t325 - t28 * t296 - t29 * t297 - t335 * t4 + t245, t13 * t335 - t96 * t61 + t297 * t40 + t325 * t187 + t321 * t100 + (-t28 + t290) * t291, -t13 * t155 + t96 * t60 + t296 * t40 + t326 * t187 + t321 * t229 + (t29 + t289) * t291, t5 * t115 - t4 * t116 + t13 * t96 - t28 * t325 + t29 * t326 - t321 * t40, t275, t31, t32, t12, -t213, t11, t1 * t155 + t100 * t327 - t18 * t296 + t2 * t335 + t229 * t328 + t23 * t297 - t80 * t60 - t81 * t61, -t3 * t155 + t72 * t60 + t296 * t26 + t327 * t187 - t324 * t229 + (qJD(2) * t81 - t23) * t291, -t3 * t335 + t72 * t61 - t297 * t26 + t328 * t187 + t324 * t100 + (-qJD(2) * t80 + t18) * t291, t1 * t80 + t18 * t328 + t2 * t81 - t23 * t327 + t26 * t324 + t3 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t227, t206 * t238, -t150 ^ 2 - t151 ^ 2, t108 * t151 - t109 * t150 + t186, 0, 0, 0, 0, 0, 0, t210, -t212, t241, t100 * t36 - t229 * t35 + t135, 0, 0, 0, 0, 0, 0, t241, t211 + t339, t212, -t100 * t29 - t229 * t28 + t13, 0, 0, 0, 0, 0, 0, t241, t212, t210, t100 * t23 - t18 * t229 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, -t272, -t344, -t314, t209, t191, -t117 * t229 + t230, t100 * t117 + t187 * t35 + t268, 0, 0, t191, t344, t38, t314, -t272, -t314, pkin(4) * t60 - t320 + (-t29 - t36) * t229 + (t28 + t300) * t100, t100 * t51 - t230 - 0.2e1 * t257 + t318, -t40 * t100 + t187 * t300 + t229 * t51 + t177 + t250, -t5 * pkin(4) - t4 * qJ(5) - t28 * t36 + t29 * t300 - t40 * t51, t191, t38, -t344, -t314, t272, t314, -t320 + t330 * t60 + (t23 + t298) * t229 + (t18 - t299) * t100, -t187 * t233 + t229 * t34 - 0.2e1 * t172 + t177 + t221, -t34 * t100 - t317 + (-0.2e1 * qJD(6) - t25) * t187 + 0.2e1 * t330 * t191 - t231, t2 * qJ(5) - t1 * t330 + t18 * t298 + t23 * t299 - t26 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t344, -t264, t340, -t187 * t29 + t318 + t5, 0, 0, 0, 0, 0, 0, -t344, t340, t264, t317 + (qJD(6) + t23) * t187 + t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t209, t191 + t314, -t97 - t181, -t18 * t187 + t221 + t338;];
tauc_reg  = t6;
