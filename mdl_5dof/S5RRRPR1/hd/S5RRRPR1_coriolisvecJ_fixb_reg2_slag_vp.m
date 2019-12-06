% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:38:51
% EndTime: 2019-12-05 18:39:01
% DurationCPUTime: 3.53s
% Computational Cost: add. (7593->326), mult. (20178->443), div. (0->0), fcn. (14944->8), ass. (0->194)
t207 = cos(qJ(2));
t278 = cos(qJ(3));
t234 = t278 * t207;
t223 = qJD(1) * t234;
t204 = sin(qJ(3));
t205 = sin(qJ(2));
t241 = qJD(1) * t205;
t233 = t204 * t241;
t159 = -t223 + t233;
t173 = t204 * t207 + t205 * t278;
t161 = qJD(1) * t173;
t201 = sin(pkin(9));
t202 = cos(pkin(9));
t119 = t159 * t202 + t161 * t201;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t216 = -t159 * t201 + t161 * t202;
t238 = qJD(5) * t206;
t239 = qJD(5) * t203;
t198 = qJD(2) + qJD(3);
t247 = t204 * t205;
t219 = t198 * t247;
t243 = t198 * t223;
t125 = qJD(1) * t219 - t243;
t135 = t198 * t173;
t126 = t135 * qJD(1);
t71 = -t125 * t201 + t126 * t202;
t72 = -t125 * t202 - t126 * t201;
t21 = t119 * t238 + t203 * t71 - t206 * t72 + t216 * t239;
t197 = qJD(5) + t198;
t64 = t119 * t206 + t203 * t216;
t264 = t64 * t197;
t8 = -t21 + t264;
t217 = t119 * t203 - t206 * t216;
t212 = qJD(5) * t217 - t203 * t72 - t206 * t71;
t265 = t197 * t217;
t9 = t212 - t265;
t272 = t217 ^ 2;
t273 = t64 ^ 2;
t10 = t272 - t273;
t287 = t217 * t64;
t114 = t216 * pkin(8);
t279 = -pkin(7) - pkin(6);
t184 = t279 * t207;
t179 = qJD(1) * t184;
t162 = t204 * t179;
t183 = t279 * t205;
t177 = qJD(1) * t183;
t266 = qJD(2) * pkin(2);
t168 = t177 + t266;
t127 = t168 * t278 + t162;
t156 = t161 * qJ(4);
t100 = t127 - t156;
t89 = pkin(3) * t198 + t100;
t166 = t278 * t179;
t128 = t168 * t204 - t166;
t253 = t159 * qJ(4);
t101 = t128 - t253;
t92 = t201 * t101;
t47 = t202 * t89 - t92;
t29 = pkin(4) * t198 - t114 + t47;
t276 = t119 * pkin(8);
t249 = t202 * t101;
t48 = t201 * t89 + t249;
t32 = t48 - t276;
t235 = qJD(2) * t279;
t224 = qJD(1) * t235;
t169 = t205 * t224;
t170 = t207 * t224;
t231 = t278 * qJD(3);
t240 = qJD(3) * t204;
t225 = -t168 * t231 - t169 * t278 - t170 * t204 - t179 * t240;
t40 = -qJ(4) * t126 - qJD(4) * t159 - t225;
t226 = -t204 * t169 + t170 * t278;
t76 = -qJD(3) * t128 + t226;
t41 = t125 * qJ(4) - t161 * qJD(4) + t76;
t19 = -t201 * t40 + t202 * t41;
t6 = -pkin(8) * t72 + t19;
t20 = t201 * t41 + t202 * t40;
t7 = -pkin(8) * t71 + t20;
t1 = (qJD(5) * t29 + t7) * t206 + t203 * t6 - t32 * t239;
t194 = -pkin(2) * t207 - pkin(1);
t182 = qJD(1) * t194;
t136 = t159 * pkin(3) + qJD(4) + t182;
t82 = t119 * pkin(4) + t136;
t211 = t82 * t64 - t1;
t12 = t203 * t29 + t206 * t32;
t2 = -qJD(5) * t12 - t203 * t7 + t206 * t6;
t210 = t217 * t82 + t2;
t132 = -t177 * t204 + t166;
t102 = t132 + t253;
t133 = t177 * t278 + t162;
t103 = -t156 + t133;
t248 = t202 * t204;
t267 = pkin(2) * qJD(3);
t261 = -t202 * t102 + t201 * t103 + (-t201 * t278 - t248) * t267;
t250 = t201 * t204;
t260 = -t201 * t102 - t202 * t103 + (t202 * t278 - t250) * t267;
t286 = t276 - t261;
t285 = t114 + t260;
t284 = t119 * t216;
t196 = t205 * t266;
t283 = 0.2e1 * t196;
t237 = qJD(1) * qJD(2);
t282 = -0.2e1 * t237;
t138 = t183 * t204 - t184 * t278;
t277 = pkin(3) * t201;
t275 = t161 * pkin(3);
t193 = pkin(2) * t278 + pkin(3);
t152 = -pkin(2) * t250 + t193 * t202;
t145 = pkin(4) + t152;
t154 = pkin(2) * t248 + t193 * t201;
t111 = t145 * t203 + t154 * t206;
t269 = qJD(5) * t111 + t203 * t285 + t206 * t286;
t110 = t145 * t206 - t154 * t203;
t268 = -qJD(5) * t110 + t203 * t286 - t206 * t285;
t172 = -t234 + t247;
t178 = t205 * t235;
t180 = t207 * t235;
t86 = t178 * t278 + t180 * t204 + t183 * t231 + t184 * t240;
t55 = -qJ(4) * t135 - qJD(4) * t172 + t86;
t134 = -qJD(2) * t234 - t207 * t231 + t219;
t87 = -qJD(3) * t138 - t204 * t178 + t180 * t278;
t56 = t134 * qJ(4) - t173 * qJD(4) + t87;
t28 = t201 * t56 + t202 * t55;
t50 = t100 * t202 - t92;
t191 = pkin(3) * t202 + pkin(4);
t155 = t191 * t203 + t206 * t277;
t49 = -t100 * t201 - t249;
t33 = t49 + t276;
t34 = -t114 + t50;
t263 = qJD(5) * t155 - t203 * t34 + t206 * t33;
t153 = t191 * t206 - t203 * t277;
t262 = -qJD(5) * t153 + t203 * t33 + t206 * t34;
t259 = t216 * t198;
t257 = t119 ^ 2;
t256 = t119 * t198;
t255 = t216 ^ 2;
t252 = t161 * t159;
t251 = t182 * t161;
t209 = qJD(1) ^ 2;
t246 = t207 * t209;
t208 = qJD(2) ^ 2;
t245 = t208 * t205;
t244 = t208 * t207;
t137 = t183 * t278 + t184 * t204;
t115 = -qJ(4) * t173 + t137;
t116 = -qJ(4) * t172 + t138;
t60 = t115 * t201 + t116 * t202;
t242 = t205 ^ 2 - t207 ^ 2;
t195 = pkin(2) * t241;
t236 = t205 * t246;
t230 = t205 * t237;
t105 = pkin(2) * t230 + pkin(3) * t126;
t122 = pkin(3) * t135 + t196;
t27 = -t201 * t55 + t202 * t56;
t227 = pkin(1) * t282;
t59 = t115 * t202 - t116 * t201;
t222 = t207 * t230;
t221 = t119 * t136 - t20;
t88 = pkin(4) * t216 + t275;
t11 = -t203 * t32 + t206 * t29;
t220 = -t11 * t64 - t12 * t217;
t218 = -t119 * t47 + t216 * t48;
t131 = -t172 * t201 + t173 * t202;
t43 = -pkin(8) * t131 + t59;
t130 = t172 * t202 + t173 * t201;
t44 = -pkin(8) * t130 + t60;
t23 = -t203 * t44 + t206 * t43;
t24 = t203 * t43 + t206 * t44;
t74 = -t130 * t203 + t131 * t206;
t142 = pkin(3) * t172 + t194;
t42 = pkin(4) * t71 + t105;
t215 = -t136 * t216 + t19;
t214 = t182 * t159 + t225;
t139 = t195 + t275;
t104 = -t159 ^ 2 + t161 ^ 2;
t96 = pkin(4) * t130 + t142;
t90 = t243 + (t159 - t233) * t198;
t83 = t195 + t88;
t81 = -t134 * t202 - t135 * t201;
t80 = -t134 * t201 + t135 * t202;
t73 = t130 * t206 + t131 * t203;
t54 = pkin(4) * t80 + t122;
t46 = t72 + t256;
t45 = -t71 + t259;
t31 = t255 - t257;
t26 = qJD(5) * t74 + t203 * t81 + t206 * t80;
t25 = t130 * t238 + t131 * t239 + t203 * t80 - t206 * t81;
t18 = -pkin(8) * t80 + t28;
t17 = -pkin(8) * t81 + t27;
t4 = -qJD(5) * t24 + t206 * t17 - t203 * t18;
t3 = qJD(5) * t23 + t203 * t17 + t206 * t18;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t222, t242 * t282, t244, -0.2e1 * t222, -t245, 0, -pkin(6) * t244 + t205 * t227, pkin(6) * t245 + t207 * t227, 0, 0, -t125 * t173 - t134 * t161, t125 * t172 - t126 * t173 + t134 * t159 - t135 * t161, -t134 * t198, t126 * t172 + t135 * t159, -t135 * t198, 0, t194 * t126 + t182 * t135 + t87 * t198 + (qJD(1) * t172 + t159) * t196, -t194 * t125 - t182 * t134 + t161 * t283 - t86 * t198, t125 * t137 - t126 * t138 + t127 * t134 - t128 * t135 - t159 * t86 - t161 * t87 + t172 * t225 - t173 * t76, t127 * t87 + t128 * t86 + t76 * t137 - t138 * t225 + t182 * t283, t131 * t72 + t216 * t81, -t119 * t81 - t130 * t72 - t131 * t71 - t216 * t80, t81 * t198, t119 * t80 + t130 * t71, -t80 * t198, 0, t105 * t130 + t119 * t122 + t136 * t80 + t142 * t71 + t198 * t27, t105 * t131 + t122 * t216 + t136 * t81 + t142 * t72 - t198 * t28, -t119 * t28 - t130 * t20 - t131 * t19 - t216 * t27 - t47 * t81 - t48 * t80 - t59 * t72 - t60 * t71, t105 * t142 + t122 * t136 + t19 * t59 + t20 * t60 + t27 * t47 + t28 * t48, -t21 * t74 + t217 * t25, t21 * t73 + t212 * t74 + t217 * t26 + t25 * t64, -t25 * t197, -t212 * t73 + t26 * t64, -t26 * t197, 0, t197 * t4 - t212 * t96 + t26 * t82 + t42 * t73 + t54 * t64, -t197 * t3 - t21 * t96 - t217 * t54 - t25 * t82 + t42 * t74, -t1 * t73 + t11 * t25 - t12 * t26 - t2 * t74 + t21 * t23 + t212 * t24 + t217 * t4 - t3 * t64, t1 * t24 + t11 * t4 + t12 * t3 + t2 * t23 + t42 * t96 + t54 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236, t242 * t209, 0, t236, 0, 0, t209 * pkin(1) * t205, pkin(1) * t246, 0, 0, t252, t104, t90, -t252, 0, 0, -t159 * t195 - t132 * t198 - t251 + (t166 + (-pkin(2) * t198 - t168) * t204) * qJD(3) + t226, t133 * t198 + (-t161 * t241 - t198 * t231) * pkin(2) + t214, (t128 + t132) * t161 + (-t127 + t133) * t159 + (t278 * t125 - t126 * t204 + (-t159 * t278 + t161 * t204) * qJD(3)) * pkin(2), -t127 * t132 - t128 * t133 + (-t182 * t241 + t278 * t76 - t204 * t225 + (-t127 * t204 + t128 * t278) * qJD(3)) * pkin(2), t284, t31, t46, -t284, t45, 0, -t139 * t119 + t198 * t261 + t215, -t139 * t216 - t198 * t260 + t221, -t119 * t260 - t152 * t72 - t154 * t71 - t216 * t261 + t218, -t136 * t139 + t19 * t152 + t20 * t154 + t260 * t48 + t261 * t47, -t287, t10, t8, t287, t9, 0, -t197 * t269 - t83 * t64 + t210, t197 * t268 + t217 * t83 + t211, t110 * t21 + t111 * t212 - t217 * t269 + t268 * t64 + t220, t1 * t111 - t11 * t269 + t2 * t110 - t12 * t268 - t82 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t104, t90, -t252, 0, 0, t128 * t198 - t251 + t76, t127 * t198 + t214, 0, 0, t284, t31, t46, -t284, t45, 0, -t119 * t275 - t198 * t49 + t215, t198 * t50 - t216 * t275 + t221, t50 * t119 + t49 * t216 + (-t201 * t71 - t202 * t72) * pkin(3) + t218, -t47 * t49 - t48 * t50 + (-t136 * t161 + t19 * t202 + t20 * t201) * pkin(3), -t287, t10, t8, t287, t9, 0, -t197 * t263 - t88 * t64 + t210, t197 * t262 + t217 * t88 + t211, t153 * t21 + t155 * t212 - t217 * t263 + t262 * t64 + t220, t1 * t155 - t11 * t263 - t12 * t262 + t2 * t153 - t82 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 + t259, t72 - t256, -t255 - t257, t119 * t48 + t216 * t47 + t105, 0, 0, 0, 0, 0, 0, -t212 - t265, -t21 - t264, -t272 - t273, -t11 * t217 + t12 * t64 + t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, t10, t8, t287, t9, 0, t12 * t197 + t210, t11 * t197 + t211, 0, 0;];
tauc_reg = t5;
