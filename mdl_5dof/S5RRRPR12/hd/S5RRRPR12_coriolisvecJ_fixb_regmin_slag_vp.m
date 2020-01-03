% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:40:35
% EndTime: 2019-12-31 21:40:52
% DurationCPUTime: 5.81s
% Computational Cost: add. (5594->430), mult. (15308->633), div. (0->0), fcn. (11868->10), ass. (0->195)
t201 = sin(pkin(5));
t209 = cos(qJ(2));
t256 = qJD(1) * t209;
t240 = t201 * t256;
t297 = qJD(3) - t240;
t205 = sin(qJ(3));
t208 = cos(qJ(3));
t203 = cos(pkin(5));
t257 = qJD(1) * t203;
t233 = qJD(2) + t257;
t206 = sin(qJ(2));
t258 = qJD(1) * t201;
t244 = t206 * t258;
t293 = -t205 * t244 + t208 * t233;
t133 = qJD(5) - t293;
t207 = cos(qJ(5));
t141 = t205 * t233 + t208 * t244;
t200 = sin(pkin(10));
t202 = cos(pkin(10));
t100 = t141 * t202 + t200 * t297;
t204 = sin(qJ(5));
t275 = t100 * t204;
t98 = t141 * t200 - t202 * t297;
t49 = t207 * t98 + t275;
t296 = t133 * t49;
t247 = pkin(1) * t257;
t161 = pkin(7) * t240 + t206 * t247;
t295 = -qJD(4) * t205 - t161 + t297 * (pkin(3) * t205 - qJ(4) * t208);
t51 = t100 * t207 - t204 * t98;
t253 = qJD(3) * t205;
t246 = pkin(8) * t253;
t158 = -pkin(7) * t244 + t209 * t247;
t218 = (pkin(2) * t206 - pkin(8) * t209) * t201;
t159 = qJD(1) * t218;
t260 = t158 * t208 + t159 * t205;
t82 = qJ(4) * t244 + t260;
t281 = t295 * t202 + (t246 + t82) * t200;
t292 = -t200 * t295 + t202 * t82;
t291 = t206 * t209;
t290 = -qJD(5) + t133;
t249 = qJD(1) * qJD(2);
t238 = t201 * t249;
t228 = t209 * t238;
t104 = qJD(3) * t141 + t205 * t228;
t210 = qJD(1) ^ 2;
t289 = pkin(1) * t206;
t288 = pkin(1) * t209;
t287 = pkin(4) * t205;
t286 = pkin(9) + qJ(4);
t122 = pkin(8) * t233 + t161;
t155 = (-pkin(2) * t209 - pkin(8) * t206 - pkin(1)) * t201;
t132 = qJD(1) * t155;
t160 = qJD(2) * t218;
t148 = qJD(1) * t160;
t267 = t201 * t206;
t191 = pkin(7) * t267;
t162 = (t203 * t288 - t191) * qJD(2);
t149 = qJD(1) * t162;
t252 = qJD(3) * t208;
t213 = t122 * t253 - t132 * t252 - t148 * t205 - t149 * t208;
t229 = t206 * t238;
t30 = qJ(4) * t229 + qJD(4) * t297 - t213;
t103 = qJD(3) * t293 + t208 * t228;
t266 = t201 * t209;
t248 = pkin(7) * t266;
t163 = (t203 * t289 + t248) * qJD(2);
t150 = qJD(1) * t163;
t34 = pkin(3) * t104 - qJ(4) * t103 - qJD(4) * t141 + t150;
t11 = t200 * t34 + t202 * t30;
t154 = t248 + (pkin(8) + t289) * t203;
t212 = -t154 * t253 + t155 * t252 + t160 * t205 + t162 * t208;
t255 = qJD(2) * t206;
t42 = (qJ(4) * t255 - qJD(4) * t209) * t201 + t212;
t168 = t203 * t205 + t208 * t267;
t242 = qJD(2) * t266;
t114 = qJD(3) * t168 + t205 * t242;
t167 = -t203 * t208 + t205 * t267;
t115 = -qJD(3) * t167 + t208 * t242;
t48 = pkin(3) * t114 - qJ(4) * t115 - qJD(4) * t168 + t163;
t14 = t200 * t48 + t202 * t42;
t121 = -pkin(2) * t233 - t158;
t57 = -pkin(3) * t293 - t141 * qJ(4) + t121;
t67 = t122 * t208 + t132 * t205;
t60 = qJ(4) * t297 + t67;
t24 = t200 * t57 + t202 * t60;
t66 = -t122 * t205 + t132 * t208;
t87 = pkin(3) * t141 - qJ(4) * t293;
t39 = t200 * t87 + t202 * t66;
t153 = t191 + (-pkin(2) - t288) * t203;
t77 = pkin(3) * t167 - qJ(4) * t168 + t153;
t261 = t154 * t208 + t155 * t205;
t78 = -qJ(4) * t266 + t261;
t37 = t200 * t77 + t202 * t78;
t268 = t200 * t208;
t123 = -t202 * t244 + t240 * t268;
t263 = t208 * t209;
t124 = (t200 * t206 + t202 * t263) * t258;
t269 = t200 * t204;
t173 = -t202 * t207 + t269;
t174 = t200 * t207 + t202 * t204;
t251 = qJD(5) * t205;
t285 = t123 * t204 - t124 * t207 - t173 * t252 - t174 * t251;
t250 = qJD(5) * t207;
t265 = t202 * t205;
t284 = -t123 * t207 - t124 * t204 + t174 * t252 + t250 * t265 - t251 * t269;
t231 = t122 * t252 + t132 * t253 - t148 * t208 + t149 * t205;
t35 = -pkin(3) * t229 + t231;
t283 = t200 * t35;
t282 = t202 * t35;
t280 = -t202 * t246 - t292;
t279 = t133 * t173;
t278 = t133 * t174;
t245 = pkin(4) * t200 + pkin(8);
t235 = -t158 * t205 + t159 * t208;
t83 = -pkin(3) * t244 - t235;
t277 = -pkin(4) * t123 + t245 * t252 - t83;
t276 = qJ(4) * t104;
t274 = t293 * t297;
t273 = t293 * t200;
t272 = t141 * t297;
t215 = t297 * t205;
t271 = t297 * t208;
t197 = t201 ^ 2;
t270 = t197 * t210;
t264 = t202 * t208;
t59 = -pkin(3) * t297 + qJD(4) - t66;
t262 = -qJD(4) + t59;
t182 = -pkin(3) * t208 - qJ(4) * t205 - pkin(2);
t138 = pkin(8) * t264 + t182 * t200;
t259 = t206 ^ 2 - t209 ^ 2;
t254 = qJD(2) * t208;
t243 = t201 * t255;
t241 = t201 * t203 * t210;
t239 = t197 * t249;
t10 = -t200 * t30 + t202 * t34;
t13 = -t200 * t42 + t202 * t48;
t23 = -t200 * t60 + t202 * t57;
t36 = -t200 * t78 + t202 * t77;
t38 = -t200 * t66 + t202 * t87;
t80 = t103 * t200 - t202 * t229;
t81 = t103 * t202 + t200 * t229;
t237 = t204 * t81 + t207 * t80;
t236 = -t154 * t205 + t155 * t208;
t234 = 0.2e1 * t239;
t232 = qJD(2) + 0.2e1 * t257;
t116 = -pkin(9) * t200 * t205 + t138;
t227 = -pkin(9) * t124 + qJD(5) * t116 + t240 * t287 - (-pkin(9) * t264 + t287) * qJD(3) - t281;
t226 = -0.2e1 * pkin(1) * t239;
t6 = pkin(4) * t104 - pkin(9) * t81 + t10;
t7 = -pkin(9) * t80 + t11;
t225 = t204 * t6 + t207 * t7;
t171 = t202 * t182;
t105 = -pkin(9) * t265 + t171 + (-pkin(8) * t200 - pkin(4)) * t208;
t224 = -pkin(9) * t123 - qJD(5) * t105 - (-pkin(8) * t265 - pkin(9) * t268) * qJD(3) + t292;
t79 = pkin(3) * t266 - t236;
t12 = -pkin(4) * t293 - pkin(9) * t100 + t23;
t17 = -pkin(9) * t98 + t24;
t3 = t12 * t207 - t17 * t204;
t4 = t12 * t204 + t17 * t207;
t113 = t168 * t202 - t200 * t266;
t19 = pkin(4) * t167 - pkin(9) * t113 + t36;
t112 = t168 * t200 + t202 * t266;
t25 = -pkin(9) * t112 + t37;
t222 = t19 * t207 - t204 * t25;
t221 = t19 * t204 + t207 * t25;
t61 = t112 * t207 + t113 * t204;
t62 = -t112 * t204 + t113 * t207;
t219 = -t154 * t252 - t155 * t253 + t160 * t208 - t162 * t205;
t184 = t286 * t200;
t217 = -pkin(9) * t273 - qJD(4) * t202 + qJD(5) * t184 + t39;
t185 = t286 * t202;
t216 = -pkin(9) * t202 * t293 + pkin(4) * t141 + qJD(4) * t200 + qJD(5) * t185 + t38;
t15 = -qJD(5) * t275 - t204 * t80 + t207 * t81 - t250 * t98;
t211 = pkin(1) * (-t203 * t249 + t270);
t47 = -pkin(3) * t243 - t219;
t2 = -qJD(5) * t4 - t204 * t7 + t207 * t6;
t16 = qJD(5) * t51 + t237;
t196 = -pkin(4) * t202 - pkin(3);
t175 = t245 * t205;
t157 = t173 * t205;
t156 = t174 * t205;
t137 = -pkin(8) * t268 + t171;
t92 = t115 * t202 + t200 * t243;
t91 = t115 * t200 - t202 * t243;
t53 = pkin(4) * t273 + t67;
t52 = pkin(4) * t112 + t79;
t41 = pkin(4) * t98 + t59;
t26 = pkin(4) * t91 + t47;
t21 = qJD(5) * t62 + t204 * t92 + t207 * t91;
t20 = -qJD(5) * t61 - t204 * t91 + t207 * t92;
t18 = pkin(4) * t80 + t35;
t9 = -pkin(9) * t91 + t14;
t8 = pkin(4) * t114 - pkin(9) * t92 + t13;
t1 = qJD(5) * t3 + t225;
t5 = [0, 0, 0, t234 * t291, -t259 * t234, t232 * t242, -t232 * t243, 0, -t150 * t203 - t163 * t233 + t206 * t226, -t149 * t203 - t162 * t233 + t209 * t226, t103 * t168 + t115 * t141, -t103 * t167 - t104 * t168 - t114 * t141 + t115 * t293, t115 * t297 + (-t103 * t209 + (qJD(1) * t168 + t141) * t255) * t201, -t114 * t297 + (t104 * t209 + (-qJD(1) * t167 + t293) * t255) * t201, (-t197 * t256 + t201 * t297) * t255, t219 * t297 - t163 * t293 + t153 * t104 + t150 * t167 + t121 * t114 + (t231 * t209 + (qJD(1) * t236 + t66) * t255) * t201, -t212 * t297 + t163 * t141 + t153 * t103 + t150 * t168 + t121 * t115 + (-t213 * t209 + (-qJD(1) * t261 - t67) * t255) * t201, t10 * t167 + t104 * t36 + t112 * t35 + t114 * t23 - t13 * t293 + t47 * t98 + t59 * t91 + t79 * t80, t100 * t47 - t104 * t37 - t11 * t167 + t113 * t35 - t114 * t24 + t14 * t293 + t59 * t92 + t79 * t81, -t10 * t113 - t100 * t13 - t11 * t112 - t14 * t98 - t23 * t92 - t24 * t91 - t36 * t81 - t37 * t80, t10 * t36 + t11 * t37 + t13 * t23 + t14 * t24 + t35 * t79 + t47 * t59, t15 * t62 + t20 * t51, -t15 * t61 - t16 * t62 - t20 * t49 - t21 * t51, t104 * t62 + t114 * t51 + t133 * t20 + t15 * t167, -t104 * t61 - t114 * t49 - t133 * t21 - t16 * t167, t104 * t167 + t114 * t133, (-qJD(5) * t221 - t204 * t9 + t207 * t8) * t133 + t222 * t104 + t2 * t167 + t3 * t114 + t26 * t49 + t52 * t16 + t18 * t61 + t41 * t21, t26 * t51 + t52 * t15 + t18 * t62 + t41 * t20 - (qJD(5) * t222 + t204 * t8 + t207 * t9) * t133 - t221 * t104 - t1 * t167 - t4 * t114; 0, 0, 0, -t270 * t291, t259 * t270, -t209 * t241, t206 * t241, 0, -pkin(7) * t228 + t161 * t233 + t206 * t211, pkin(7) * t229 + t158 * t233 + t209 * t211, t103 * t205 + t141 * t271, (t103 + t274) * t208 + (-t104 - t272) * t205, t297 * t252 + (-t297 * t263 + (t205 * qJD(2) - t141) * t206) * t258, -t297 * t253 + (t209 * t215 + (-t293 + t254) * t206) * t258, -t297 * t244, -pkin(2) * t104 - t150 * t208 - t235 * t297 + t161 * t293 + (-pkin(8) * t271 + t121 * t205) * qJD(3) + (-t66 * t206 + (-pkin(8) * t255 - t121 * t209) * t205) * t258, -pkin(2) * t103 + t150 * t205 + t260 * t297 - t161 * t141 + (pkin(8) * t215 + t121 * t208) * qJD(3) + (-t121 * t263 + (-pkin(8) * t254 + t67) * t206) * t258, t104 * t137 - t123 * t59 - t83 * t98 - t281 * t293 + (-t10 + (pkin(8) * t98 + t200 * t59) * qJD(3)) * t208 + (pkin(8) * t80 + t23 * t297 + t283) * t205, -t100 * t83 - t104 * t138 - t124 * t59 + t280 * t293 + (t11 + (pkin(8) * t100 + t202 * t59) * qJD(3)) * t208 + (pkin(8) * t81 - t24 * t297 + t282) * t205, t123 * t24 + t124 * t23 - t137 * t81 - t138 * t80 - t280 * t98 + (-t10 * t202 - t11 * t200) * t205 - t281 * t100 + (-t200 * t24 - t202 * t23) * t252, t10 * t137 + t11 * t138 - t59 * t83 + t280 * t24 + t281 * t23 + (t205 * t35 + t252 * t59) * pkin(8), -t15 * t157 + t285 * t51, -t15 * t156 + t157 * t16 - t284 * t51 - t285 * t49, -t104 * t157 + t133 * t285 - t15 * t208 + t215 * t51, -t104 * t156 - t133 * t284 + t16 * t208 - t215 * t49, -t104 * t208 + t133 * t215, (t105 * t207 - t116 * t204) * t104 - t2 * t208 + t175 * t16 + t18 * t156 + t277 * t49 + t284 * t41 + t3 * t215 + (t204 * t224 - t207 * t227) * t133, t175 * t15 - t18 * t157 - (t105 * t204 + t116 * t207) * t104 + t1 * t208 + t277 * t51 + t285 * t41 - t4 * t215 + (t204 * t227 + t207 * t224) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t141 * t293, t141 ^ 2 - t293 ^ 2, t103 - t274, -t104 + t272, t229, -t121 * t141 + t297 * t67 - t231, -t121 * t293 + t297 * t66 + t213, -t200 * t276 - pkin(3) * t80 - t141 * t23 - t282 - t67 * t98 - (t200 * t262 - t38) * t293, -t202 * t276 - pkin(3) * t81 - t100 * t67 + t141 * t24 + t283 - (t202 * t262 + t39) * t293, t100 * t38 + t39 * t98 + (-qJ(4) * t80 - qJD(4) * t98 + t23 * t293 + t11) * t202 + (qJ(4) * t81 + qJD(4) * t100 + t24 * t293 - t10) * t200, -pkin(3) * t35 - t23 * t38 - t24 * t39 - t59 * t67 + (-t200 * t23 + t202 * t24) * qJD(4) + (-t10 * t200 + t11 * t202) * qJ(4), t15 * t174 - t279 * t51, -t15 * t173 - t16 * t174 - t278 * t51 + t279 * t49, t104 * t174 - t133 * t279 - t141 * t51, -t104 * t173 - t133 * t278 + t141 * t49, -t133 * t141, (-t184 * t207 - t185 * t204) * t104 + t196 * t16 + t18 * t173 - t3 * t141 - t53 * t49 + t278 * t41 + (t204 * t217 - t207 * t216) * t133, t196 * t15 + t18 * t174 - (-t184 * t204 + t185 * t207) * t104 - t53 * t51 + t4 * t141 - t279 * t41 + (t204 * t216 + t207 * t217) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 * t293 + t80, t293 * t98 + t81, -t100 ^ 2 - t98 ^ 2, t100 * t23 + t24 * t98 + t35, 0, 0, 0, 0, 0, t133 * t51 + t16, t15 - t296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51 * t49, -t49 ^ 2 + t51 ^ 2, t15 + t296, t290 * t51 - t237, t104, t133 * t4 - t41 * t51 + t2, t290 * t3 + t41 * t49 - t225;];
tauc_reg = t5;
