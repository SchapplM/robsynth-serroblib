% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:20:59
% EndTime: 2019-03-09 09:21:11
% DurationCPUTime: 4.23s
% Computational Cost: add. (3480->438), mult. (8996->597), div. (0->0), fcn. (6481->8), ass. (0->214)
t156 = sin(pkin(6));
t161 = sin(qJ(2));
t252 = qJD(1) * t161;
t221 = t156 * t252;
t164 = cos(qJ(2));
t157 = cos(pkin(6));
t240 = t157 * qJD(1);
t238 = pkin(1) * t240;
t256 = pkin(8) * t221 - t164 * t238;
t290 = qJD(3) + t256;
t292 = -qJ(4) * t221 + t290;
t138 = qJD(2) + t240;
t160 = sin(qJ(5));
t163 = cos(qJ(5));
t253 = qJD(1) * t156;
t227 = t164 * t253;
t84 = -t163 * t138 + t160 * t227;
t259 = qJD(6) - t84;
t117 = qJD(5) + t221;
t159 = sin(qJ(6));
t162 = cos(qJ(6));
t85 = t160 * t138 + t163 * t227;
t50 = -t162 * t117 - t159 * t85;
t291 = t259 * t50;
t237 = pkin(1) * qJD(2) * t157;
t134 = t164 * t237;
t147 = t157 * qJD(3);
t250 = qJD(2) * t161;
t226 = t156 * t250;
t288 = (pkin(8) * t250 + qJD(4) * t164) * t156;
t47 = -qJ(4) * t226 - t134 - t147 + t288;
t208 = qJD(1) * t237;
t115 = t164 * t208;
t118 = t138 * qJD(3);
t239 = qJD(1) * qJD(2);
t218 = t156 * t239;
t203 = t161 * t218;
t231 = qJ(4) * t203 + t115 + t118;
t34 = qJD(1) * t288 - t231;
t52 = t159 * t117 - t162 * t85;
t289 = t117 * t52;
t283 = pkin(1) * t161;
t144 = t157 * t283;
t264 = t156 * t164;
t287 = pkin(8) * t264 + t144;
t165 = -pkin(2) - pkin(3);
t286 = qJD(2) * t165;
t120 = t164 * t218;
t153 = -pkin(9) + t165;
t31 = t153 * t138 + t292;
t182 = (pkin(4) * t161 + pkin(9) * t164) * t156;
t79 = -pkin(1) * t253 - pkin(2) * t227 - qJ(3) * t221;
t60 = pkin(3) * t227 + qJD(4) - t79;
t33 = qJD(1) * t182 + t60;
t14 = t160 * t33 + t163 * t31;
t172 = t156 * (pkin(4) * t164 + t153 * t161);
t170 = qJD(2) * t172;
t265 = t156 * t161;
t135 = qJD(3) * t265;
t258 = qJ(3) * t120 + qJD(1) * t135;
t28 = qJD(1) * t170 + t258;
t248 = qJD(2) * t164;
t228 = qJ(4) * t248;
t247 = qJD(4) * t161;
t86 = pkin(8) * t120 + t161 * t208;
t46 = (-t228 - t247) * t253 + t86;
t169 = -t14 * qJD(5) - t160 * t46 + t163 * t28;
t4 = -pkin(5) * t120 - t169;
t285 = (-t85 * pkin(5) + t259 * pkin(10)) * t259 + t4;
t88 = -t156 * pkin(1) - pkin(2) * t264 - qJ(3) * t265;
t72 = pkin(3) * t264 - t88;
t44 = t182 + t72;
t140 = pkin(8) * t265;
t199 = -qJ(4) * t265 + t140;
t229 = -pkin(1) * t164 - pkin(2);
t211 = -pkin(3) + t229;
t54 = (-pkin(9) + t211) * t157 + t199;
t190 = t160 * t44 + t163 * t54;
t225 = t156 * t248;
t257 = qJ(3) * t225 + t135;
t32 = t170 + t257;
t59 = -t156 * t247 + (t144 + (pkin(8) - qJ(4)) * t264) * qJD(2);
t284 = -qJD(5) * t190 - t160 * t59 + t163 * t32;
t224 = qJD(5) * t264;
t204 = qJD(1) * t224;
t243 = qJD(5) * t163;
t55 = -t138 * t243 + t160 * t204 + t163 * t203;
t245 = qJD(5) * t160;
t56 = t138 * t245 - t160 * t203 + t163 * t204;
t10 = -t56 * pkin(5) - t55 * pkin(10) - t34;
t12 = t117 * pkin(10) + t14;
t121 = t138 * qJ(3);
t92 = pkin(8) * t227 + t161 * t238;
t76 = -qJ(4) * t227 + t92;
t57 = -t121 - t76;
t41 = t138 * pkin(4) - t57;
t20 = -t84 * pkin(5) + t85 * pkin(10) + t41;
t197 = t159 * t12 - t162 * t20;
t179 = t160 * t28 + t163 * t46 + t33 * t243 - t31 * t245;
t3 = pkin(10) * t120 + t179;
t1 = -t197 * qJD(6) + t159 * t10 + t162 * t3;
t80 = t159 * t163 * t221 - t162 * t227;
t282 = t80 * t259;
t263 = t161 * t163;
t81 = (t159 * t164 + t162 * t263) * t253;
t281 = t81 * t259;
t158 = qJ(3) + pkin(4);
t280 = -t292 + t117 * (pkin(5) * t160 - pkin(10) * t163);
t124 = qJ(3) * t227;
t42 = qJD(1) * t172 + t124;
t279 = t160 * t42 + t163 * t76;
t278 = pkin(8) * qJD(2);
t277 = t117 * t85;
t276 = t159 * t56;
t275 = t159 * t259;
t241 = qJD(6) * t162;
t242 = qJD(6) * t159;
t16 = t117 * t241 + t159 * t120 + t162 * t55 + t85 * t242;
t274 = t16 * t159;
t273 = t162 * t52;
t272 = t162 * t56;
t271 = t162 * t259;
t270 = t84 * t117;
t269 = t256 * t138;
t268 = t117 * t160;
t267 = t117 * t163;
t152 = t156 ^ 2;
t266 = t152 * qJD(1) ^ 2;
t260 = -qJD(4) - t60;
t154 = t161 ^ 2;
t155 = t164 ^ 2;
t255 = t154 - t155;
t254 = qJ(3) * qJD(2);
t251 = qJD(2) * t160;
t249 = qJD(2) * t163;
t246 = qJD(5) * t159;
t244 = qJD(5) * t162;
t236 = t161 * t268;
t235 = t117 * t263;
t234 = t164 * t266;
t233 = t259 * t246;
t232 = t259 * t244;
t87 = t157 * qJ(3) + t287;
t223 = t117 * t245;
t222 = t117 * t243;
t219 = t152 * t239;
t215 = -t162 * t120 + t159 * t55;
t214 = t260 * t161;
t107 = t163 * pkin(5) + t160 * pkin(10) + t158;
t213 = -pkin(10) * t227 + qJD(6) * t107 - t279;
t212 = t138 + t240;
t210 = 0.2e1 * t219;
t209 = t165 * t265;
t207 = t161 * t234;
t206 = t138 * t227;
t202 = t92 * t138 - t86;
t201 = -0.2e1 * pkin(1) * t219;
t6 = t162 * t12 + t159 * t20;
t93 = t287 * qJD(2);
t196 = -t93 * t138 - t86 * t157;
t19 = pkin(10) * t265 + t190;
t71 = qJ(4) * t264 - t87;
t61 = t157 * pkin(4) - t71;
t94 = -t157 * t163 + t160 * t264;
t95 = -t157 * t160 - t163 * t264;
t25 = -t94 * pkin(5) - t95 * pkin(10) + t61;
t195 = t159 * t25 + t162 * t19;
t194 = -t159 * t19 + t162 * t25;
t13 = -t160 * t31 + t163 * t33;
t192 = -t160 * t76 + t163 * t42;
t191 = -t160 * t54 + t163 * t44;
t189 = qJD(2) * t209;
t188 = -t138 ^ 2 - t154 * t266;
t187 = -pkin(8) * t226 + t134;
t186 = -t120 + t206;
t185 = -t241 * t259 + t276;
t184 = t242 * t259 + t272;
t183 = -t95 * t159 + t162 * t265;
t69 = t159 * t265 + t95 * t162;
t180 = -t153 * t248 - t161 * t41;
t178 = t160 * t32 + t163 * t59 + t44 * t243 - t54 * t245;
t177 = t259 * t117;
t176 = t117 * t50;
t175 = -pkin(8) * t203 + t115;
t11 = -t117 * pkin(5) - t13;
t171 = pkin(10) * t56 + (t11 + t13) * t259;
t2 = -t6 * qJD(6) + t162 * t10 - t159 * t3;
t168 = t184 + t289;
t167 = t176 + t185;
t96 = t138 * t221;
t90 = pkin(2) * t221 - t124;
t89 = t229 * t157 + t140;
t83 = t147 + t187;
t77 = pkin(2) * t226 - t257;
t75 = qJD(1) * t209 + t124;
t73 = t121 + t92;
t70 = -t138 * pkin(2) + t290;
t67 = qJD(5) * t94 + t163 * t226;
t66 = -t157 * t245 + t160 * t226 - t163 * t224;
t65 = t118 + t175;
t63 = pkin(2) * t203 - t258;
t62 = t211 * t157 + t199;
t58 = t189 + t257;
t45 = qJD(1) * t189 + t258;
t39 = t165 * t138 + t292;
t24 = qJD(6) * t183 + t159 * t225 + t67 * t162;
t23 = qJD(6) * t69 + t67 * t159 - t162 * t225;
t21 = -pkin(5) * t227 - t192;
t18 = -pkin(5) * t265 - t191;
t17 = qJD(6) * t52 + t215;
t15 = t66 * pkin(5) - t67 * pkin(10) - t47;
t8 = -pkin(5) * t225 - t284;
t7 = pkin(10) * t225 + t178;
t5 = [0, 0, 0, t161 * t164 * t210, -t255 * t210, t212 * t225, -t212 * t226, 0, t161 * t201 + t196, -t138 * t187 - t157 * t175 + t164 * t201 (t79 * t250 - t164 * t63 + (-t164 * t77 + t88 * t250) * qJD(1)) * t156 + t196 (t161 * t86 + t164 * t65 + (-t161 * t73 + t164 * t70) * qJD(2) + (t161 * t93 + t164 * t83 + (-t161 * t87 + t164 * t89) * qJD(2)) * qJD(1)) * t156, t83 * t138 + t65 * t157 + (-t79 * t248 - t161 * t63 + (-t161 * t77 - t88 * t248) * qJD(1)) * t156, t63 * t88 + t65 * t87 + t70 * t93 + t73 * t83 + t79 * t77 + t86 * t89, -t47 * t138 - t34 * t157 + (t60 * t248 + t161 * t45 + (t161 * t58 + t72 * t248) * qJD(1)) * t156, t59 * t138 + t46 * t157 + (t60 * t250 - t164 * t45 + (-t164 * t58 + t72 * t250) * qJD(1)) * t156 (-t161 * t46 + t164 * t34 + (-t161 * t57 - t164 * t39) * qJD(2) + (-t161 * t59 + t164 * t47 + (-t161 * t71 - t164 * t62) * qJD(2)) * qJD(1)) * t156, t34 * t71 + t39 * t59 + t45 * t72 + t46 * t62 + t57 * t47 + t60 * t58, t55 * t95 - t85 * t67, t55 * t94 + t95 * t56 + t85 * t66 + t67 * t84, t67 * t117 + (t161 * t55 + (qJD(1) * t95 - t85) * t248) * t156, -t66 * t117 + (t161 * t56 + (qJD(1) * t94 + t84) * t248) * t156 (t117 * t156 + t152 * t252) * t248, t284 * t117 + t47 * t84 - t61 * t56 + t34 * t94 + t41 * t66 + (t169 * t161 + (qJD(1) * t191 + t13) * t248) * t156, -t178 * t117 + t47 * t85 + t61 * t55 - t34 * t95 + t41 * t67 + (-t179 * t161 + (-t190 * qJD(1) - t14) * t248) * t156, t16 * t69 + t52 * t24, t16 * t183 - t69 * t17 - t52 * t23 - t24 * t50, -t16 * t94 + t24 * t259 + t52 * t66 - t69 * t56, t17 * t94 - t183 * t56 - t23 * t259 - t50 * t66, t259 * t66 + t56 * t94 (-qJD(6) * t195 + t162 * t15 - t159 * t7) * t259 - t194 * t56 - t2 * t94 - t197 * t66 + t8 * t50 + t18 * t17 - t4 * t183 + t11 * t23 -(qJD(6) * t194 + t159 * t15 + t162 * t7) * t259 + t195 * t56 + t1 * t94 - t6 * t66 + t8 * t52 + t18 * t16 + t4 * t69 + t11 * t24; 0, 0, 0, -t207, t255 * t266, -t186, t96 - t203, 0, t266 * t283 + t202, pkin(1) * t234 - t175 - t269 (-t161 * t79 + t164 * t90) * t253 + t202 ((t73 - t92 - t254) * t161 + (-pkin(2) * qJD(2) + t290 - t70) * t164) * t253, t269 + t115 + 0.2e1 * t118 + (t164 * t79 + (t90 - t278) * t161) * t253, -t86 * pkin(2) + t65 * qJ(3) + t290 * t73 - t70 * t92 - t79 * t90, t292 * t138 + (t260 * t164 + (-t75 - t278) * t161) * t253 + t231, -t76 * t138 + ((-qJ(4) * qJD(2) + t75) * t164 + t214) * t253 + t86 ((t57 + t76 + t254) * t161 + (-t292 + t39 - t286) * t164) * t253, -t34 * qJ(3) + t46 * t165 - t292 * t57 - t39 * t76 - t60 * t75, -t55 * t160 + t85 * t267 (-t55 - t270) * t163 + (-t56 - t277) * t160, -t222 + (-t235 + (t85 - t251) * t164) * t253, t223 + (t236 + (-t84 - t249) * t164) * t253, -t117 * t227, -t158 * t56 - t34 * t163 - t192 * t117 - t292 * t84 + (-t153 * t267 - t41 * t160) * qJD(5) + (-t13 * t164 + t160 * t180) * t253, t158 * t55 + t34 * t160 + t279 * t117 - t292 * t85 + (t153 * t268 - t41 * t163) * qJD(5) + (t14 * t164 + t163 * t180) * t253, -t16 * t162 * t160 + (t160 * t242 - t162 * t243 - t81) * t52, t81 * t50 + t52 * t80 + (t159 * t52 + t162 * t50) * t243 + (t274 + t162 * t17 + (-t159 * t50 + t273) * qJD(6)) * t160, -t281 + (t16 - t232) * t163 + (t184 - t289) * t160, t282 + (-t17 + t233) * t163 + (t176 - t185) * t160, -t56 * t163 - t259 * t268, -t107 * t272 - t11 * t80 - t21 * t50 - (t159 * t213 + t162 * t280) * t259 + (-t11 * t246 + t2 + (qJD(5) * t50 + t185) * t153) * t163 + (t197 * t221 - t11 * t241 + t153 * t17 - t4 * t159 + (t153 * t275 + t197) * qJD(5)) * t160, t107 * t276 - t11 * t81 - t21 * t52 - (-t159 * t280 + t162 * t213) * t259 + (-t11 * t244 - t1 + (qJD(5) * t52 + t184) * t153) * t163 + (t6 * t221 + t11 * t242 + t153 * t16 - t4 * t162 + (t153 * t271 + t6) * qJD(5)) * t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t207, -t186, t188, -t73 * t138 + t79 * t221 + t86, t188, t207, t186, t57 * t138 + (t214 - t228) * t253 + t86, 0, 0, 0, 0, 0, -t222 + t138 * t84 + (-t160 * t248 - t235) * t253, t223 + t138 * t85 + (-t163 * t248 + t236) * t253, 0, 0, 0, 0, 0, -t138 * t271 + (t159 * t177 + t17) * t160 + t167 * t163, t138 * t275 + (t162 * t177 + t16) * t160 + t168 * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120 + t206, t96 + t203 (-t154 - t155) * t266 (-t164 * t57 + (t39 + t286) * t161) * t253 + t258, 0, 0, 0, 0, 0, -t223 + (-t236 + (-t84 + t249) * t164) * t253, -t222 + (-t235 + (-t85 - t251) * t164) * t253, 0, 0, 0, 0, 0, -t282 + (-t17 - t233) * t163 + t167 * t160, -t281 + (-t16 - t232) * t163 + t168 * t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t84, -t84 ^ 2 + t85 ^ 2, t55 - t270, t56 - t277, t120, t14 * t117 + t41 * t85 + t169, t117 * t13 - t41 * t84 - t179, t259 * t273 + t274 (t16 - t291) * t162 + (-t259 * t52 - t17) * t159, t259 * t271 + t52 * t85 - t276, -t259 * t275 - t50 * t85 - t272, t259 * t85, -pkin(5) * t17 - t14 * t50 + t171 * t159 - t162 * t285 - t197 * t85, -pkin(5) * t16 - t14 * t52 + t159 * t285 + t171 * t162 - t6 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t50, -t50 ^ 2 + t52 ^ 2, t16 + t291, -t215 + (-qJD(6) + t259) * t52, -t56, -t11 * t52 + t259 * t6 + t2, t11 * t50 - t197 * t259 - t1;];
tauc_reg  = t5;
