% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:43:05
% EndTime: 2019-03-09 10:43:14
% DurationCPUTime: 4.32s
% Computational Cost: add. (7065->455), mult. (21463->603), div. (0->0), fcn. (17167->10), ass. (0->233)
t181 = sin(pkin(11));
t185 = sin(qJ(2));
t188 = cos(qJ(2));
t288 = cos(pkin(11));
t211 = t181 * t188 + t185 * t288;
t182 = sin(pkin(6));
t271 = qJD(1) * t182;
t153 = t211 * t271;
t184 = sin(qJ(4));
t187 = cos(qJ(4));
t289 = cos(pkin(6));
t247 = t289 * qJD(1);
t221 = t247 + qJD(2);
t194 = -t187 * t153 - t184 * t221;
t316 = qJD(6) - t194;
t208 = t187 * t221;
t118 = t184 * t153 - t208;
t252 = t288 * t188;
t237 = t182 * t252;
t167 = qJD(1) * t237;
t259 = t185 * t271;
t150 = -t181 * t259 + t167;
t146 = qJD(4) - t150;
t183 = sin(qJ(6));
t186 = cos(qJ(6));
t72 = -t186 * t118 + t146 * t183;
t319 = t316 * t72;
t269 = qJD(4) * t184;
t280 = t150 * t184;
t318 = t269 - t280;
t285 = t118 * t146;
t262 = qJD(1) * qJD(2);
t254 = t182 * t262;
t238 = t185 * t254;
t203 = qJD(2) * t167 - t181 * t238;
t70 = -qJD(4) * t208 + t153 * t269 - t187 * t203;
t317 = t70 - t285;
t260 = pkin(1) * t289;
t240 = t188 * t260;
t172 = qJD(1) * t240;
t302 = pkin(8) + qJ(3);
t256 = t302 * t185;
t239 = t182 * t256;
t199 = pkin(2) * t289 - t239;
t122 = qJD(2) * pkin(2) + qJD(1) * t199 + t172;
t241 = t185 * t260;
t276 = t182 * t188;
t149 = t302 * t276 + t241;
t138 = t149 * qJD(1);
t253 = t288 * t138;
t67 = t181 * t122 + t253;
t57 = pkin(9) * t221 + t67;
t236 = (-pkin(2) * t188 - pkin(1)) * t182;
t217 = qJD(1) * t236;
t161 = qJD(3) + t217;
t84 = -t150 * pkin(3) - t153 * pkin(9) + t161;
t34 = t184 * t57 - t187 * t84;
t273 = qJD(5) + t34;
t270 = qJD(2) * t182;
t258 = t185 * t270;
t315 = pkin(2) * t258;
t267 = qJD(4) * t187;
t229 = qJD(2) * t247;
t219 = pkin(1) * t229;
t170 = t188 * t219;
t196 = (-qJD(2) * t256 + qJD(3) * t188) * t182;
t112 = qJD(1) * t196 + t170;
t277 = t182 * t185;
t124 = -t149 * qJD(2) - qJD(3) * t277;
t113 = t124 * qJD(1);
t50 = t112 * t288 + t181 * t113;
t156 = t211 * t182;
t151 = qJD(2) * t156;
t144 = qJD(1) * t151;
t169 = pkin(2) * t238;
t85 = t144 * pkin(3) - pkin(9) * t203 + t169;
t250 = t184 * t50 - t187 * t85 + t57 * t267 + t84 * t269;
t305 = pkin(4) * t144;
t11 = t250 - t305;
t35 = t184 * t84 + t187 * t57;
t28 = -qJ(5) * t146 - t35;
t314 = -t146 * t28 - t11;
t282 = t194 * t146;
t71 = -qJD(4) * t194 + t184 * t203;
t313 = -t71 - t282;
t106 = t153 * t186 + t183 * t280;
t264 = qJD(6) * t187;
t200 = t183 * t269 - t186 * t264 - t106;
t312 = t200 * t316;
t137 = -qJD(1) * t239 + t172;
t88 = t181 * t137 + t253;
t311 = pkin(4) * t318 - t184 * qJD(5) - t88;
t274 = -pkin(5) * t194 + t273;
t310 = -qJD(6) + t316;
t304 = pkin(5) * t118;
t19 = -t28 - t304;
t306 = pkin(4) + pkin(10);
t309 = t306 * t70 + (t19 - t35 + t304) * t316;
t308 = t194 ^ 2;
t307 = t146 ^ 2;
t190 = qJD(1) ^ 2;
t175 = pkin(2) * t181 + pkin(9);
t303 = pkin(5) + t175;
t265 = qJD(6) * t186;
t266 = qJD(6) * t183;
t30 = t118 * t265 + t186 * t144 - t146 * t266 + t183 * t71;
t74 = t118 * t183 + t146 * t186;
t301 = t30 * t184 + t74 * t267;
t155 = t181 * t277 - t237;
t104 = t155 * pkin(3) - t156 * pkin(9) + t236;
t136 = t240 + t199;
t98 = t181 * t136 + t288 * t149;
t87 = pkin(9) * t289 + t98;
t300 = t184 * t104 + t187 * t87;
t100 = pkin(2) * t259 + pkin(3) * t153 - pkin(9) * t150;
t127 = t181 * t138;
t89 = t137 * t288 - t127;
t299 = t184 * t100 + t187 * t89;
t66 = t122 * t288 - t127;
t56 = -pkin(3) * t221 - t66;
t191 = qJ(5) * t194 + t56;
t33 = t118 * pkin(4) + t191;
t298 = t194 * t33;
t296 = t146 * t72;
t295 = t150 * t74;
t294 = t183 * t70;
t248 = t183 * t144 - t186 * t71;
t31 = qJD(6) * t74 + t248;
t293 = t184 * t31;
t63 = t186 * t70;
t292 = t30 * t186;
t36 = -qJ(5) * t153 - t299;
t291 = pkin(5) * t280 - t303 * t269 + t36;
t287 = qJ(5) * t187;
t290 = -qJ(5) * t267 + t150 * t287 + t311;
t286 = t118 * qJ(5);
t284 = t118 * t153;
t283 = t194 * t118;
t281 = t194 * t153;
t139 = t144 * qJ(5);
t279 = t175 * t144;
t178 = t182 ^ 2;
t278 = t178 * t190;
t275 = t184 * t144;
t272 = t185 ^ 2 - t188 ^ 2;
t268 = qJD(4) * t186;
t261 = t188 * t278;
t257 = t175 * t269;
t164 = t303 * t187;
t255 = t178 * t262;
t251 = -t184 * t85 - t187 * t50 - t84 * t267 + t57 * t269;
t249 = t187 * t104 - t184 * t87;
t49 = t181 * t112 - t288 * t113;
t173 = qJD(2) * t240;
t123 = t173 + t196;
t59 = t181 * t123 - t288 * t124;
t245 = t146 * t187;
t244 = t316 * t183;
t242 = 0.2e1 * t255;
t40 = -qJ(5) * t155 - t300;
t176 = -pkin(2) * t288 - pkin(3);
t235 = t182 * t190 * t289;
t234 = -0.2e1 * pkin(1) * t255;
t163 = t303 * t184;
t233 = -qJD(6) * t163 - t311 - t146 * (pkin(10) * t184 - t287);
t207 = -t184 * qJ(5) + t176;
t157 = -t306 * t187 + t207;
t82 = t184 * t89;
t232 = qJD(6) * t157 - qJD(4) * t164 + t82 + (pkin(5) * t150 - t100) * t187 - t306 * t153;
t105 = t153 * t183 - t186 * t280;
t231 = (t183 * t264 + t184 * t268 + t105) * t316;
t198 = t70 * qJ(5) + qJD(5) * t194 + t49;
t12 = t306 * t71 + t198;
t6 = -t70 * pkin(5) - t306 * t144 + t250;
t228 = t186 * t12 + t183 * t6;
t17 = -t306 * t146 + t274;
t22 = t306 * t118 + t191;
t3 = t17 * t186 - t183 * t22;
t4 = t17 * t183 + t186 * t22;
t131 = t156 * t187 + t184 * t289;
t21 = t131 * pkin(5) - t306 * t155 - t249;
t130 = t156 * t184 - t187 * t289;
t97 = t136 * t288 - t181 * t149;
t86 = -pkin(3) * t289 - t97;
t193 = -t131 * qJ(5) + t86;
t32 = t306 * t130 + t193;
t225 = -t183 * t32 + t186 * t21;
t224 = t183 * t21 + t186 * t32;
t223 = t130 * t186 - t155 * t183;
t103 = t130 * t183 + t155 * t186;
t152 = (-t181 * t185 + t252) * t270;
t101 = pkin(3) * t151 - pkin(9) * t152 + t315;
t60 = t123 * t288 + t181 * t124;
t222 = t187 * t101 - t104 * t269 - t184 * t60 - t87 * t267;
t220 = 0.2e1 * t247 + qJD(2);
t10 = -t146 * qJD(5) - t139 + t251;
t216 = t187 * t144 - t146 * t318;
t215 = t146 * t35 - t250;
t214 = -t244 * t316 - t63;
t213 = t184 * t101 + t104 * t267 + t187 * t60 - t269 * t87;
t7 = -pkin(5) * t71 - t10;
t212 = t7 + (t306 * t316 + t286) * t316;
t210 = t146 * t56 - t279;
t209 = -t146 * t33 + t279;
t206 = -t186 * t316 ^ 2 + t294;
t205 = t146 * t245 + t275;
t204 = -pkin(8) * t276 - t241;
t202 = -t216 + t284;
t96 = -t156 * t269 + (qJD(4) * t289 + t152) * t187;
t197 = -t96 * qJ(5) - t131 * qJD(5) + t59;
t2 = -qJD(6) * t4 - t183 * t12 + t186 * t6;
t13 = -qJ(5) * t151 - qJD(5) * t155 - t213;
t192 = t221 * t204;
t162 = -t187 * pkin(4) + t207;
t95 = qJD(4) * t131 + t152 * t184;
t65 = t70 * t184;
t58 = -pkin(4) * t194 + t286;
t46 = t70 * t131;
t43 = t130 * pkin(4) + t193;
t41 = -pkin(4) * t155 - t249;
t39 = qJD(6) * t223 + t151 * t186 + t95 * t183;
t38 = qJD(6) * t103 + t151 * t183 - t95 * t186;
t37 = -t153 * pkin(4) - t187 * t100 + t82;
t26 = -pkin(5) * t130 - t40;
t25 = -pkin(4) * t146 + t273;
t18 = pkin(4) * t95 + t197;
t16 = pkin(4) * t71 + t198;
t15 = t306 * t95 + t197;
t14 = -pkin(4) * t151 - t222;
t9 = -pkin(5) * t95 - t13;
t8 = t96 * pkin(5) - t306 * t151 - t222;
t1 = qJD(6) * t3 + t228;
t5 = [0, 0, 0, t185 * t188 * t242, -t272 * t242, t220 * t188 * t270, -t220 * t258, 0, qJD(2) * t192 + t185 * t234 + t204 * t229, t188 * t234 - (-pkin(8) * t258 + t173) * t221 - (-pkin(8) * t238 + t170) * t289, -t98 * t144 + t60 * t150 - t67 * t151 - t66 * t152 + t59 * t153 - t50 * t155 + t49 * t156 - t203 * t97, -t49 * t97 + t50 * t98 - t66 * t59 + t67 * t60 + (t161 + t217) * t315, -t194 * t96 - t46, -t118 * t96 + t130 * t70 - t131 * t71 + t194 * t95, t131 * t144 + t146 * t96 - t151 * t194 - t155 * t70, -t118 * t151 - t130 * t144 - t146 * t95 - t155 * t71, t144 * t155 + t146 * t151, t59 * t118 + t49 * t130 + t144 * t249 + t146 * t222 - t34 * t151 - t155 * t250 + t56 * t95 + t86 * t71, t49 * t131 - t144 * t300 - t146 * t213 - t35 * t151 + t155 * t251 - t194 * t59 + t56 * t96 - t86 * t70, t10 * t130 + t11 * t131 + t118 * t13 - t14 * t194 + t25 * t96 + t28 * t95 + t40 * t71 - t41 * t70, t11 * t155 - t118 * t18 - t130 * t16 + t14 * t146 + t144 * t41 + t151 * t25 - t33 * t95 - t43 * t71, -t10 * t155 - t13 * t146 - t131 * t16 - t144 * t40 - t151 * t28 + t18 * t194 - t33 * t96 + t43 * t70, t10 * t40 + t11 * t41 + t13 * t28 + t14 * t25 + t16 * t43 + t18 * t33, t103 * t30 + t39 * t74, -t103 * t31 + t223 * t30 - t38 * t74 - t39 * t72, -t103 * t70 + t131 * t30 + t316 * t39 + t74 * t96, -t131 * t31 - t223 * t70 - t316 * t38 - t72 * t96, t316 * t96 - t46 (-qJD(6) * t224 - t183 * t15 + t186 * t8) * t316 - t225 * t70 + t2 * t131 + t3 * t96 + t9 * t72 + t26 * t31 - t7 * t223 + t19 * t38 -(qJD(6) * t225 + t186 * t15 + t183 * t8) * t316 + t224 * t70 - t1 * t131 - t4 * t96 + t9 * t74 + t26 * t30 + t7 * t103 + t19 * t39; 0, 0, 0, -t185 * t261, t272 * t278, -t188 * t235, t185 * t235, 0, -pkin(8) * t188 * t254 - qJD(1) * t192 + (pkin(1) * t278 - t219) * t185, -t170 + pkin(1) * t261 + (-pkin(8) * t259 + t172) * t247 + t172 * qJD(2) (t67 - t88) * t153 + (-t89 + t66) * t150 + (-t181 * t144 - t203 * t288) * pkin(2), t66 * t88 - t67 * t89 + (-t161 * t259 + t181 * t50 - t288 * t49) * pkin(2), -t194 * t245 - t65 (-t70 - t285) * t187 + (-t71 + t282) * t184, t205 + t281, t216 + t284, -t146 * t153, -t88 * t118 + t34 * t153 + t176 * t71 - t49 * t187 + (t82 + (-qJD(4) * t175 - t100) * t187) * t146 + t210 * t184, t88 * t194 + t35 * t153 - t176 * t70 + t49 * t184 + (t257 + t299) * t146 + t210 * t187, -t36 * t118 + t37 * t194 + (-t150 * t25 - t175 * t71 - t10 + (-t175 * t194 + t25) * qJD(4)) * t187 + (-t150 * t28 - t175 * t70 + t11 + (t118 * t175 + t28) * qJD(4)) * t184, -t25 * t153 + t16 * t187 - t162 * t71 + (t175 * t267 - t37) * t146 - t290 * t118 + t209 * t184, t28 * t153 - t16 * t184 + t162 * t70 + (t36 - t257) * t146 + t290 * t194 + t209 * t187, t16 * t162 - t25 * t37 - t28 * t36 + t290 * t33 + (-t10 * t187 + t11 * t184 + (t184 * t28 + t187 * t25) * qJD(4)) * t175, -t30 * t183 * t187 + t200 * t74, t74 * t105 + t106 * t72 + (-t183 * t72 + t186 * t74) * t269 + (t183 * t31 - t292 + (t183 * t74 + t186 * t72) * qJD(6)) * t187 (t294 - t295) * t187 + t312 + t301, -t293 + (t63 - t296) * t187 + t231, t245 * t316 - t65 -(-t157 * t183 + t163 * t186) * t70 + t164 * t31 - t19 * t105 + t291 * t72 + (-t19 * t268 + t2) * t184 + (t183 * t233 - t186 * t232) * t316 + (t146 * t3 + t7 * t186 - t19 * t266) * t187 (t157 * t186 + t163 * t183) * t70 + t164 * t30 - t19 * t106 + t291 * t74 + (qJD(4) * t183 * t19 - t1) * t184 + (t183 * t232 + t186 * t233) * t316 + (-t146 * t4 - t7 * t183 - t19 * t265) * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150 ^ 2 - t153 ^ 2, -t150 * t67 + t153 * t66 + t169, 0, 0, 0, 0, 0, -t202, -t307 * t187 - t275 + t281, t313 * t184 + t187 * t317, t202, t205 - t281, -t33 * t153 + t314 * t187 + (t146 * t25 - t10) * t184, 0, 0, 0, 0, 0, t293 + (t63 + t296) * t187 + t231 (-t294 - t295) * t187 - t312 + t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t283, -t118 ^ 2 + t308, -t317, t313, t144, t194 * t56 + t215, t118 * t56 - t146 * t34 + t251, pkin(4) * t70 - qJ(5) * t71 - (-t28 - t35) * t194 + (t25 - t273) * t118, t118 * t58 - t215 - t298 - 0.2e1 * t305, -t33 * t118 + t146 * t273 - t194 * t58 - t10 + t139, -t11 * pkin(4) - t10 * qJ(5) - t25 * t35 - t273 * t28 - t33 * t58, -t244 * t74 + t292 (-t316 * t74 - t31) * t186 + (-t30 + t319) * t183, t74 * t118 + t214, -t72 * t118 + t206, t316 * t118, qJ(5) * t31 + t3 * t118 + t212 * t183 + t309 * t186 + t274 * t72, qJ(5) * t30 - t4 * t118 - t309 * t183 + t212 * t186 + t274 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t317, t144 + t283, -t307 - t308, -t298 - t314, 0, 0, 0, 0, 0, t214 - t296, -t146 * t74 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t72, -t72 ^ 2 + t74 ^ 2, t30 + t319, t310 * t74 - t248, -t70, -t19 * t74 + t316 * t4 + t2, t19 * t72 + t310 * t3 - t228;];
tauc_reg  = t5;
