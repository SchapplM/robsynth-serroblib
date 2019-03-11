% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:21:00
% EndTime: 2019-03-09 06:21:10
% DurationCPUTime: 4.03s
% Computational Cost: add. (8683->401), mult. (22067->515), div. (0->0), fcn. (16944->8), ass. (0->193)
t183 = sin(pkin(10));
t184 = cos(pkin(10));
t187 = sin(qJ(3));
t189 = cos(qJ(3));
t154 = t183 * t189 + t184 * t187;
t301 = -t183 * t187 + t189 * t184;
t310 = t301 * qJD(1);
t186 = sin(qJ(4));
t282 = -pkin(9) - pkin(8);
t232 = qJD(4) * t282;
t258 = t310 * t186;
t290 = t154 * qJD(1);
t103 = pkin(3) * t290 - pkin(8) * t310;
t188 = cos(qJ(4));
t279 = pkin(7) + qJ(2);
t162 = t279 * t183;
t155 = qJD(1) * t162;
t163 = t279 * t184;
t156 = qJD(1) * t163;
t294 = -t155 * t189 - t187 * t156;
t268 = t186 * t103 + t188 * t294;
t309 = -pkin(9) * t258 - t186 * t232 + t268;
t257 = t310 * t188;
t93 = t188 * t103;
t308 = -pkin(4) * t290 + pkin(9) * t257 + t186 * t294 + t188 * t232 - t93;
t141 = qJD(4) - t310;
t136 = qJD(5) + t141;
t149 = t154 * qJD(3);
t137 = qJD(1) * t149;
t185 = sin(qJ(5));
t281 = cos(qJ(5));
t231 = t281 * t186;
t158 = t185 * t188 + t231;
t247 = t185 * t186;
t202 = t281 * t188 - t247;
t289 = qJD(4) + qJD(5);
t223 = t281 * qJD(5);
t292 = t281 * qJD(4) + t223;
t266 = -t292 * t188 + t202 * t310 + t289 * t247;
t215 = -t266 * t136 + t158 * t137;
t122 = qJD(3) * t188 - t186 * t290;
t123 = qJD(3) * t186 + t188 * t290;
t204 = t185 * t122 + t281 * t123;
t269 = t204 * t290;
t307 = t215 + t269;
t114 = t289 * t158;
t265 = -t158 * t310 + t114;
t216 = -t265 * t136 + t202 * t137;
t65 = -t281 * t122 + t123 * t185;
t273 = t290 * t65;
t306 = t216 - t273;
t101 = -qJD(3) * pkin(3) - t294;
t60 = -pkin(4) * t122 + t101;
t28 = pkin(5) * t65 - qJ(6) * t204 + t60;
t305 = t28 * t65;
t304 = t60 * t65;
t280 = t204 * t65;
t303 = t154 * qJD(2);
t148 = t301 * qJD(3);
t238 = qJD(4) * t188;
t302 = t148 * t186 + t154 * t238;
t300 = t290 * qJD(3);
t283 = t204 ^ 2;
t299 = -t65 ^ 2 + t283;
t193 = qJD(3) * (qJD(4) + t310);
t239 = qJD(4) * t186;
t192 = -t188 * t193 + t239 * t290;
t194 = qJD(1) * t148;
t233 = qJD(3) * t239 + t186 * t194 + t238 * t290;
t237 = qJD(5) * t185;
t30 = -t122 * t223 + t123 * t237 + t185 * t233 + t281 * t192;
t17 = t136 * t65 - t30;
t42 = pkin(5) * t204 + qJ(6) * t65;
t164 = t282 * t186;
t165 = t282 * t188;
t203 = t281 * t164 + t185 * t165;
t297 = -t203 * qJD(5) - t308 * t185 + t309 * t281;
t121 = t185 * t164 - t281 * t165;
t296 = -t121 * qJD(5) + t309 * t185 + t308 * t281;
t108 = -t187 * t155 + t189 * t156;
t218 = -t108 + (t239 - t258) * pkin(4);
t246 = t186 * t137;
t125 = t188 * t137;
t295 = t141 * t239 - t125;
t118 = t162 * t189 + t187 * t163;
t293 = -t123 * t239 - t192 * t188;
t228 = t154 * t239;
t255 = t148 * t188;
t291 = -t228 + t255;
t132 = t137 * pkin(5);
t102 = qJD(3) * pkin(8) + t108;
t89 = t137 * pkin(3) - pkin(8) * t194;
t86 = t188 * t89;
t214 = -t102 * t238 + t86;
t195 = t301 * qJD(2);
t70 = qJD(1) * t195 + qJD(3) * t294;
t270 = t186 * t70;
t176 = -pkin(2) * t184 - pkin(1);
t161 = qJD(1) * t176 + qJD(2);
t81 = -pkin(3) * t310 - pkin(8) * t290 + t161;
t11 = t137 * pkin(4) + pkin(9) * t192 - t81 * t239 + t214 - t270;
t200 = -t102 * t239 + t186 * t89 + t188 * t70 + t81 * t238;
t16 = -t233 * pkin(9) + t200;
t54 = -t102 * t186 + t188 * t81;
t44 = -pkin(9) * t123 + t54;
t36 = pkin(4) * t141 + t44;
t55 = t188 * t102 + t186 * t81;
t45 = pkin(9) * t122 + t55;
t221 = -t281 * t11 + t185 * t16 + t45 * t223 + t36 * t237;
t2 = -t132 + t221;
t198 = t204 * t28 + t2;
t31 = t122 * t237 + t123 * t223 - t185 * t192 + t281 * t233;
t288 = t136 * t204 - t31;
t287 = -t60 * t204 - t221;
t286 = -t141 ^ 2 * t188 - t246;
t285 = -t202 * t30 - t204 * t265;
t119 = -t162 * t187 + t163 * t189;
t253 = t154 * t188;
t105 = -pkin(3) * t301 - pkin(8) * t154 + t176;
t98 = t188 * t105;
t50 = -pkin(4) * t301 - pkin(9) * t253 - t119 * t186 + t98;
t254 = t154 * t186;
t112 = t188 * t119;
t267 = t186 * t105 + t112;
t56 = -pkin(9) * t254 + t267;
t206 = t185 * t50 + t281 * t56;
t82 = -t118 * qJD(3) + t195;
t104 = pkin(3) * t149 - pkin(8) * t148;
t94 = t188 * t104;
t23 = -pkin(9) * t255 + t149 * pkin(4) - t186 * t82 + t94 + (-t112 + (pkin(9) * t154 - t105) * t186) * qJD(4);
t199 = t186 * t104 + t105 * t238 - t119 * t239 + t188 * t82;
t27 = -pkin(9) * t302 + t199;
t284 = -qJD(5) * t206 - t185 * t27 + t281 * t23;
t278 = qJ(6) * t290 + t297;
t277 = -pkin(5) * t290 + t296;
t276 = -t265 * pkin(5) - t266 * qJ(6) + qJD(6) * t158 - t218;
t234 = t281 * t45;
t13 = t185 * t36 + t234;
t274 = t13 * t136;
t271 = t185 * t45;
t19 = t281 * t44 - t271;
t264 = pkin(4) * t223 + qJD(6) - t19;
t263 = t203 * t137;
t262 = t121 * t137;
t261 = t122 * t290;
t260 = t123 * t141;
t259 = t123 * t290;
t245 = t188 * t122;
t12 = t281 * t36 - t271;
t243 = qJD(6) - t12;
t242 = t183 ^ 2 + t184 ^ 2;
t241 = qJD(3) * t187;
t240 = qJD(3) * t189;
t236 = qJD(1) * qJD(2);
t180 = -pkin(4) * t188 - pkin(3);
t226 = t101 * t238;
t222 = -t185 * t11 - t281 * t16 - t36 * t223 + t45 * t237;
t220 = t242 * qJD(1) ^ 2;
t71 = t303 * qJD(1) - t155 * t241 + t156 * t240;
t83 = -t162 * t241 + t163 * t240 + t303;
t219 = -t158 * t31 + t266 * t65;
t18 = t185 * t44 + t234;
t217 = pkin(4) * t237 - t18;
t84 = pkin(4) * t254 + t118;
t127 = t136 * qJD(6);
t129 = t137 * qJ(6);
t1 = t129 + t127 - t222;
t212 = 0.2e1 * t242 * t236;
t211 = t141 * t258 - t295;
t210 = t12 * t136 + t222;
t208 = -t185 * t56 + t281 * t50;
t58 = t302 * pkin(4) + t83;
t205 = t185 * t23 + t50 * t223 - t56 * t237 + t281 * t27;
t49 = t233 * pkin(4) + t71;
t179 = -t281 * pkin(4) - pkin(5);
t175 = pkin(4) * t185 + qJ(6);
t109 = t137 * t301;
t106 = -pkin(5) * t202 - qJ(6) * t158 + t180;
t96 = t202 * t154;
t95 = t158 * t154;
t39 = pkin(5) * t95 - qJ(6) * t96 + t84;
t38 = t148 * t231 - t185 * t228 - t237 * t254 + (t148 * t185 + t292 * t154) * t188;
t37 = t114 * t154 - t202 * t148;
t33 = pkin(4) * t123 + t42;
t25 = pkin(5) * t301 - t208;
t24 = -qJ(6) * t301 + t206;
t8 = t136 * qJ(6) + t13;
t7 = -t136 * pkin(5) + t243;
t6 = pkin(5) * t38 + qJ(6) * t37 - qJD(6) * t96 + t58;
t5 = t31 * pkin(5) + t30 * qJ(6) - qJD(6) * t204 + t49;
t4 = -t149 * pkin(5) - t284;
t3 = qJ(6) * t149 - qJD(6) * t301 + t205;
t9 = [0, 0, 0, 0, 0, t212, qJ(2) * t212, t148 * t290 + t154 * t194, -t154 * t137 + t148 * t310 - t149 * t290 + t194 * t301, t148 * qJD(3), -t149 * qJD(3), 0, -qJD(3) * t83 + t137 * t176 + t149 * t161, t161 * t148 + (t176 * t310 - t82) * qJD(3), t123 * t255 + t293 * t154, -t122 * t228 - t123 * t302 + t148 * t245 + t192 * t254 - t233 * t253, t123 * t149 + t154 * t125 + t291 * t141 + t192 * t301, t122 * t149 - t141 * t302 - t154 * t246 + t233 * t301, t141 * t149 - t109 (-t119 * t238 + t94) * t141 + t98 * t137 - t214 * t301 + t54 * t149 - t83 * t122 + t118 * t233 + t154 * t226 + ((-qJD(4) * t105 - t82) * t141 - t119 * t137 - (-qJD(4) * t81 - t70) * t301 + t71 * t154 + t101 * t148) * t186, t291 * t101 - t118 * t192 + t83 * t123 - t267 * t137 - t199 * t141 - t55 * t149 + t200 * t301 + t71 * t253, -t204 * t37 - t30 * t96, -t204 * t38 + t30 * t95 - t31 * t96 + t37 * t65, -t136 * t37 + t137 * t96 + t149 * t204 + t30 * t301, -t136 * t38 - t137 * t95 - t149 * t65 + t301 * t31, t136 * t149 - t109, t12 * t149 + t136 * t284 + t208 * t137 + t221 * t301 + t84 * t31 + t60 * t38 + t49 * t95 + t58 * t65, -t13 * t149 - t205 * t136 - t206 * t137 + t204 * t58 - t222 * t301 - t84 * t30 - t60 * t37 + t49 * t96, -t136 * t4 - t137 * t25 - t149 * t7 + t2 * t301 + t28 * t38 + t31 * t39 + t5 * t95 + t6 * t65, -t1 * t95 + t2 * t96 + t204 * t4 - t24 * t31 - t25 * t30 - t3 * t65 - t37 * t7 - t38 * t8, -t1 * t301 + t136 * t3 + t137 * t24 + t149 * t8 - t204 * t6 + t28 * t37 + t30 * t39 - t5 * t96, t1 * t24 + t2 * t25 + t28 * t6 + t3 * t8 + t39 * t5 + t4 * t7; 0, 0, 0, 0, 0, -t220, -qJ(2) * t220, 0, 0, 0, 0, 0, 0.2e1 * t300, 0.2e1 * t310 * qJD(3), 0, 0, 0, 0, 0, t211 + t261, -t259 + t286, 0, 0, 0, 0, 0, t306, -t307, t306, t219 - t285, t307, t1 * t158 - t2 * t202 + t265 * t7 - t266 * t8 - t28 * t290; 0, 0, 0, 0, 0, 0, 0, -t290 * t310, t290 ^ 2 - t310 ^ 2, 0, 0, 0, qJD(3) * t108 - t161 * t290 - t71, -t161 * t310 - t236 * t301, -qJD(4) * t186 ^ 2 * t290 + (t186 * t193 + t260) * t188, t122 * t238 + t123 * t258 - t186 * t233 - t245 * t310 + t293, -t259 - t286, t211 - t261, -t141 * t290, -pkin(8) * t246 - pkin(3) * t233 + t108 * t122 - t54 * t290 - t71 * t188 + (-pkin(8) * t238 - t93 + (t101 + t294) * t186) * t141, pkin(3) * t192 + t295 * pkin(8) - t101 * t257 - t108 * t123 + t268 * t141 + t71 * t186 + t290 * t55 + t226, -t30 * t158 - t204 * t266, t219 + t285, t215 - t269, t216 + t273, -t136 * t290, -t12 * t290 + t296 * t136 + t180 * t31 - t202 * t49 + t218 * t65 + t265 * t60 + t263, t13 * t290 + t297 * t136 + t49 * t158 - t180 * t30 + t218 * t204 - t266 * t60 - t262, t106 * t31 + t277 * t136 - t202 * t5 + t265 * t28 - t276 * t65 + t290 * t7 + t263, t1 * t202 - t121 * t31 + t2 * t158 + t203 * t30 - t204 * t277 - t265 * t8 - t266 * t7 + t278 * t65, t106 * t30 - t278 * t136 - t5 * t158 + t204 * t276 + t266 * t28 - t290 * t8 + t262, t1 * t121 + t5 * t106 - t2 * t203 - t276 * t28 - t277 * t7 - t278 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123 * t122, -t122 ^ 2 + t123 ^ 2, -t122 * t141 - t192, -t233 + t260, t137, -t101 * t123 - t270 + t86 + (-qJD(4) + t141) * t55, -t101 * t122 + t141 * t54 - t200, t280, t299, t17, t288, t137, t18 * t136 + (-t123 * t65 - t136 * t237 + t281 * t137) * pkin(4) + t287, t19 * t136 + t304 + (-t123 * t204 - t136 * t223 - t185 * t137) * pkin(4) + t222, -t136 * t217 - t179 * t137 - t33 * t65 - t198, -t175 * t31 - t179 * t30 + (t217 + t8) * t204 + (-t264 + t7) * t65, t264 * t136 + t175 * t137 + t204 * t33 + t1 - t305, t1 * t175 + t2 * t179 + t217 * t7 + t264 * t8 - t28 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, t299, t17, t288, t137, t274 + t287, t210 + t304, -t42 * t65 + t132 - t198 + t274, pkin(5) * t30 - t31 * qJ(6) + (-t13 + t8) * t204 + (t7 - t243) * t65, t204 * t42 + 0.2e1 * t127 + 0.2e1 * t129 - t210 - t305, -t2 * pkin(5) + t1 * qJ(6) - t7 * t13 + t243 * t8 - t28 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280 - t300, t17, -t136 ^ 2 - t283, -t136 * t8 + t198;];
tauc_reg  = t9;
