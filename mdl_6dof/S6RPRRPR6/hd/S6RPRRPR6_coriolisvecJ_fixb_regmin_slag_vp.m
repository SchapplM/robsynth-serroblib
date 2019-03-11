% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:17:42
% EndTime: 2019-03-09 05:17:55
% DurationCPUTime: 4.61s
% Computational Cost: add. (7637->385), mult. (19979->527), div. (0->0), fcn. (15875->10), ass. (0->203)
t217 = cos(pkin(10));
t223 = cos(qJ(3));
t272 = t223 * t217;
t206 = qJD(1) * t272;
t215 = sin(pkin(10));
t220 = sin(qJ(3));
t278 = t215 * t220;
t253 = qJD(1) * t278;
t179 = t206 - t253;
t170 = qJD(4) - t179;
t218 = sin(qJ(6));
t221 = cos(qJ(6));
t261 = qJD(6) * t218;
t190 = t215 * t223 + t217 * t220;
t181 = t190 * qJD(1);
t219 = sin(qJ(4));
t222 = cos(qJ(4));
t260 = t222 * qJD(3);
t150 = t181 * t219 - t260;
t152 = qJD(3) * t219 + t181 * t222;
t214 = sin(pkin(11));
t216 = cos(pkin(11));
t240 = -t150 * t216 - t152 * t214;
t297 = t221 * t240;
t203 = qJD(3) * t206;
t165 = -qJD(3) * t253 + t203;
t263 = qJD(4) * t219;
t106 = qJD(4) * t260 + t222 * t165 - t181 * t263;
t107 = t152 * qJD(4) + t219 * t165;
t66 = -t106 * t214 - t107 * t216;
t67 = t106 * t216 - t107 * t214;
t97 = t150 * t214 - t152 * t216;
t12 = qJD(6) * t297 + t218 * t66 + t221 * t67 + t261 * t97;
t164 = qJD(6) + t170;
t52 = t218 * t97 + t297;
t301 = t164 * t52;
t327 = t12 - t301;
t312 = t218 * t240 - t221 * t97;
t326 = t312 * t52;
t189 = t214 * t222 + t216 * t219;
t316 = t170 * t189;
t237 = t214 * t219 - t216 * t222;
t315 = t170 * t237;
t325 = t312 ^ 2 - t52 ^ 2;
t209 = -pkin(2) * t217 - pkin(1);
t196 = t209 * qJD(1) + qJD(2);
t111 = -t179 * pkin(3) - t181 * pkin(8) + t196;
t305 = pkin(7) + qJ(2);
t197 = t305 * t215;
t194 = qJD(1) * t197;
t198 = t305 * t217;
t195 = qJD(1) * t198;
t143 = -t220 * t194 + t223 * t195;
t136 = qJD(3) * pkin(8) + t143;
t81 = t111 * t219 + t136 * t222;
t61 = -qJ(5) * t150 + t81;
t299 = t216 * t61;
t80 = t222 * t111 - t136 * t219;
t60 = -qJ(5) * t152 + t80;
t47 = pkin(4) * t170 + t60;
t29 = t214 * t47 + t299;
t311 = pkin(9) * t240;
t16 = t29 + t311;
t15 = t16 * t261;
t183 = t190 * qJD(3);
t166 = qJD(1) * t183;
t188 = -t272 + t278;
t227 = t188 * qJD(2);
t309 = -t194 * t223 - t220 * t195;
t102 = -qJD(1) * t227 + qJD(3) * t309;
t123 = pkin(3) * t166 - pkin(8) * t165;
t117 = t222 * t123;
t225 = -t81 * qJD(4) - t219 * t102 + t117;
t21 = t166 * pkin(4) - t106 * qJ(5) - t152 * qJD(5) + t225;
t262 = qJD(4) * t222;
t228 = t222 * t102 + t111 * t262 + t219 * t123 - t136 * t263;
t24 = -qJ(5) * t107 - qJD(5) * t150 + t228;
t6 = t216 * t21 - t214 * t24;
t2 = pkin(5) * t166 - pkin(9) * t67 + t6;
t135 = -qJD(3) * pkin(3) - t309;
t91 = pkin(4) * t150 + qJD(5) + t135;
t48 = -pkin(5) * t240 + t91;
t324 = -t218 * t2 - t48 * t52 + t15;
t13 = qJD(6) * t312 + t218 * t67 - t221 * t66;
t296 = t312 * t164;
t322 = -t13 + t296;
t304 = -qJ(5) - pkin(8);
t248 = qJD(4) * t304;
t137 = pkin(3) * t181 - pkin(8) * t179;
t269 = t219 * t137 + t222 * t309;
t285 = t179 * t219;
t321 = qJ(5) * t285 + t222 * qJD(5) + t219 * t248 - t269;
t129 = t222 * t137;
t320 = -pkin(4) * t181 - t129 + (qJ(5) * t179 + t248) * t222 + (-qJD(5) + t309) * t219;
t319 = t263 - t285;
t7 = t214 * t21 + t216 * t24;
t3 = pkin(9) * t66 + t7;
t254 = t221 * t2 - t218 * t3;
t318 = -t48 * t312 + t254;
t317 = pkin(9) * t97;
t141 = t189 * t221 - t218 * t237;
t302 = t141 * qJD(6) - t315 * t218 + t316 * t221;
t314 = t190 * qJD(2);
t182 = t188 * qJD(3);
t252 = t190 * t262;
t313 = -t182 * t219 + t252;
t294 = -t321 * t214 + t320 * t216;
t293 = t320 * t214 + t321 * t216;
t146 = t197 * t223 + t220 * t198;
t308 = t319 * pkin(4) - t143;
t239 = -t189 * t218 - t221 * t237;
t303 = t239 * qJD(6) - t218 * t316 - t221 * t315;
t307 = -t141 * t166 - t303 * t164;
t306 = pkin(4) * t214;
t112 = -t146 * qJD(3) - t227;
t138 = pkin(3) * t183 + pkin(8) * t182;
t130 = t222 * t138;
t139 = pkin(3) * t188 - pkin(8) * t190 + t209;
t147 = -t197 * t220 + t198 * t223;
t145 = t222 * t147;
t236 = qJ(5) * t182 - qJD(5) * t190;
t36 = t183 * pkin(4) - t219 * t112 + t130 + t236 * t222 + (-t145 + (qJ(5) * t190 - t139) * t219) * qJD(4);
t257 = t222 * t112 + t219 * t138 + t139 * t262;
t40 = -qJ(5) * t252 + (-qJD(4) * t147 + t236) * t219 + t257;
t11 = t214 * t36 + t216 * t40;
t54 = t214 * t61;
t33 = t216 * t60 - t54;
t132 = t222 * t139;
t282 = t190 * t222;
t71 = pkin(4) * t188 - qJ(5) * t282 - t147 * t219 + t132;
t268 = t219 * t139 + t145;
t283 = t190 * t219;
t82 = -qJ(5) * t283 + t268;
t43 = t214 * t71 + t216 * t82;
t300 = t181 * t52;
t28 = t216 * t47 - t54;
t14 = pkin(5) * t170 + t28 + t317;
t298 = t221 * t14;
t295 = t312 * t181;
t292 = pkin(5) * t316 + t308;
t291 = t106 * t219;
t289 = t150 * t170;
t288 = t150 * t181;
t287 = t152 * t170;
t286 = t152 * t181;
t273 = t219 * t166;
t156 = t222 * t166;
t199 = t304 * t219;
t200 = t304 * t222;
t149 = t214 * t199 - t216 * t200;
t267 = t215 ^ 2 + t217 ^ 2;
t266 = qJD(3) * t220;
t265 = qJD(3) * t223;
t264 = qJD(4) * t190;
t259 = qJD(1) * qJD(2);
t255 = -pkin(4) * t222 - pkin(3);
t250 = qJD(6) * t14 + t3;
t10 = -t214 * t40 + t216 * t36;
t32 = -t214 * t60 - t299;
t42 = -t214 * t82 + t216 * t71;
t247 = t267 * qJD(1) ^ 2;
t148 = t216 * t199 + t200 * t214;
t246 = t170 * t222;
t103 = qJD(1) * t314 - t194 * t266 + t195 * t265;
t113 = -t197 * t266 + t198 * t265 + t314;
t245 = -t164 * t302 + t239 * t166;
t121 = -pkin(9) * t189 + t148;
t244 = pkin(9) * t316 - qJD(6) * t121 - t293;
t122 = -pkin(9) * t237 + t149;
t243 = pkin(5) * t181 - pkin(9) * t315 + qJD(6) * t122 - t294;
t242 = pkin(4) * t283 + t146;
t5 = t218 * t14 + t221 * t16;
t125 = t189 * t190;
t126 = t237 * t190;
t241 = -t221 * t125 + t126 * t218;
t86 = -t125 * t218 - t126 * t221;
t235 = 0.2e1 * t267 * t259;
t234 = -t319 * t170 + t156;
t208 = pkin(4) * t216 + pkin(5);
t233 = t208 * t218 + t221 * t306;
t232 = t208 * t221 - t218 * t306;
t231 = pkin(4) * t313 + t113;
t230 = -t182 * t222 - t190 * t263;
t70 = pkin(4) * t107 + t103;
t229 = -pkin(8) * t166 + t170 * t135;
t157 = pkin(5) * t237 + t255;
t144 = t166 * t188;
t88 = pkin(4) * t152 - pkin(5) * t97;
t87 = pkin(5) * t125 + t242;
t84 = -t237 * t182 + t189 * t264;
t83 = t189 * t182 + t237 * t264;
t44 = -pkin(5) * t83 + t231;
t41 = -pkin(5) * t66 + t70;
t31 = -pkin(9) * t125 + t43;
t30 = pkin(5) * t188 + pkin(9) * t126 + t42;
t26 = t86 * qJD(6) - t218 * t84 - t221 * t83;
t25 = t241 * qJD(6) + t218 * t83 - t221 * t84;
t18 = t33 + t317;
t17 = t32 - t311;
t9 = pkin(9) * t83 + t11;
t8 = pkin(5) * t183 + pkin(9) * t84 + t10;
t4 = -t16 * t218 + t298;
t1 = [0, 0, 0, 0, 0, t235, qJ(2) * t235, t165 * t190 - t181 * t182, -t165 * t188 - t166 * t190 - t179 * t182 - t181 * t183, -t182 * qJD(3), -t183 * qJD(3), 0, -qJD(3) * t113 + t166 * t209 + t183 * t196, -qJD(3) * t112 + t165 * t209 - t182 * t196, t106 * t282 + t230 * t152 -(-t150 * t222 - t152 * t219) * t182 + (-t291 - t107 * t222 + (t150 * t219 - t152 * t222) * qJD(4)) * t190, t106 * t188 + t152 * t183 + t190 * t156 + t230 * t170, -t107 * t188 - t150 * t183 - t170 * t313 - t190 * t273, t170 * t183 + t144 (-t147 * t262 + t130) * t170 + t132 * t166 + (-t136 * t262 + t117) * t188 + t80 * t183 + t113 * t150 + t146 * t107 + t135 * t252 + ((-qJD(4) * t139 - t112) * t170 - t147 * t166 + (-qJD(4) * t111 - t102) * t188 + t103 * t190 - t135 * t182) * t219 -(-t147 * t263 + t257) * t170 - t268 * t166 - t228 * t188 - t81 * t183 + t113 * t152 + t146 * t106 + t103 * t282 + t230 * t135, t10 * t97 + t11 * t240 - t125 * t7 + t126 * t6 + t28 * t84 + t29 * t83 - t42 * t67 + t43 * t66, t28 * t10 + t29 * t11 + t231 * t91 + t242 * t70 + t6 * t42 + t7 * t43, t12 * t86 + t25 * t312, t12 * t241 - t13 * t86 + t25 * t52 - t26 * t312, t12 * t188 + t164 * t25 + t166 * t86 + t183 * t312, -t13 * t188 - t164 * t26 + t166 * t241 + t183 * t52, t164 * t183 + t144 (-t218 * t9 + t221 * t8) * t164 + (-t218 * t31 + t221 * t30) * t166 + t254 * t188 + t4 * t183 - t44 * t52 + t87 * t13 - t41 * t241 + t48 * t26 + ((-t218 * t30 - t221 * t31) * t164 - t5 * t188) * qJD(6), t87 * t12 + t15 * t188 - t5 * t183 + t48 * t25 + t41 * t86 + t44 * t312 + (-(-qJD(6) * t31 + t8) * t164 - t30 * t166 - t2 * t188) * t218 + (-(qJD(6) * t30 + t9) * t164 - t31 * t166 - t250 * t188) * t221; 0, 0, 0, 0, 0, -t247, -qJ(2) * t247, 0, 0, 0, 0, 0, 0.2e1 * t181 * qJD(3), t203 + (t179 - t253) * qJD(3), 0, 0, 0, 0, 0, t234 - t288, -t170 ^ 2 * t222 - t273 - t286, t189 * t66 + t237 * t67 - t240 * t315 - t316 * t97, -t91 * t181 + t7 * t189 - t237 * t6 - t28 * t316 - t29 * t315, 0, 0, 0, 0, 0, t245 + t300, -t295 + t307; 0, 0, 0, 0, 0, 0, 0, -t181 * t179, -t179 ^ 2 + t181 ^ 2, t203 + (-t179 - t253) * qJD(3), 0, 0, qJD(3) * t143 - t181 * t196 - t103, -t196 * t179 + t188 * t259, t152 * t246 + t291 (t106 - t289) * t222 + (-t107 - t287) * t219, t170 * t246 + t273 - t286, t234 + t288, -t170 * t181, -pkin(3) * t107 - t103 * t222 - t143 * t150 - t80 * t181 + (-pkin(8) * t262 - t129) * t170 + (t170 * t309 + t229) * t219, -pkin(3) * t106 + t103 * t219 - t143 * t152 + t81 * t181 + (pkin(8) * t263 + t269) * t170 + t229 * t222, -t148 * t67 + t149 * t66 - t6 * t189 - t237 * t7 + t240 * t293 + t28 * t315 - t29 * t316 + t294 * t97, t6 * t148 + t7 * t149 + t70 * t255 + t294 * t28 + t293 * t29 + t308 * t91, t12 * t141 + t303 * t312, t12 * t239 - t141 * t13 - t302 * t312 + t303 * t52, -t295 - t307, t245 - t300, -t164 * t181 (t121 * t221 - t122 * t218) * t166 + t157 * t13 - t41 * t239 - t4 * t181 - t292 * t52 + t302 * t48 + (t218 * t244 - t221 * t243) * t164 -(t121 * t218 + t122 * t221) * t166 + t157 * t12 + t41 * t141 + t5 * t181 + t292 * t312 + t303 * t48 + (t218 * t243 + t221 * t244) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152 * t150, -t150 ^ 2 + t152 ^ 2, t106 + t289, -t107 + t287, t166, -t135 * t152 + t81 * t170 + t225, t135 * t150 + t170 * t80 - t228 (t214 * t66 - t216 * t67) * pkin(4) + (t28 - t33) * t240 + (-t32 - t29) * t97, -t28 * t32 - t29 * t33 + (-t152 * t91 + t214 * t7 + t216 * t6) * pkin(4), -t326, t325, t327, t322, t166, t232 * t166 - (t17 * t221 - t18 * t218) * t164 + t88 * t52 + (-t164 * t233 - t5) * qJD(6) + t318, -t233 * t166 - t221 * t3 + (t17 * t218 + t18 * t221) * t164 - t88 * t312 + (-t164 * t232 - t298) * qJD(6) + t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240 ^ 2 - t97 ^ 2, -t240 * t29 - t28 * t97 + t70, 0, 0, 0, 0, 0, t13 + t296, t12 + t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t326, t325, t327, t322, t166 (-qJD(6) + t164) * t5 + t318, t4 * t164 - t221 * t250 + t324;];
tauc_reg  = t1;
