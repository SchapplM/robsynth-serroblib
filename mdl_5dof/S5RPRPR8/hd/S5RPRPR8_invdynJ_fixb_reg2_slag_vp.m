% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:22:16
% EndTime: 2019-12-31 18:22:23
% DurationCPUTime: 3.88s
% Computational Cost: add. (3634->384), mult. (7820->520), div. (0->0), fcn. (5235->14), ass. (0->200)
t163 = cos(qJ(3));
t150 = g(3) * t163;
t161 = sin(qJ(3));
t152 = qJ(1) + pkin(8);
t143 = sin(t152);
t145 = cos(t152);
t199 = g(1) * t145 + g(2) * t143;
t174 = t199 * t161 - t150;
t156 = sin(pkin(8));
t137 = t156 * pkin(1) + pkin(6);
t123 = t137 * qJD(1);
t231 = qJD(3) * t163;
t121 = t137 * qJDD(1);
t271 = qJD(2) * qJD(3) + t121;
t42 = t163 * qJDD(2) - t123 * t231 - t271 * t161;
t35 = -qJDD(3) * pkin(3) + qJDD(4) - t42;
t168 = t35 - t174;
t157 = cos(pkin(9));
t148 = t163 * qJDD(1);
t207 = t161 * t148;
t225 = qJD(1) * qJD(3);
t153 = t161 ^ 2;
t154 = t163 ^ 2;
t235 = t153 - t154;
t266 = t225 * t235 - t207;
t275 = t266 * t157;
t155 = sin(pkin(9));
t228 = t157 * qJD(3);
t233 = qJD(1) * t161;
t103 = -t155 * t233 + t228;
t229 = t155 * qJD(3);
t104 = t157 * t233 + t229;
t160 = sin(qJ(5));
t264 = cos(qJ(5));
t51 = -t264 * t103 + t160 * t104;
t274 = t51 ^ 2;
t54 = t160 * t103 + t264 * t104;
t273 = t54 ^ 2;
t226 = t163 * qJD(1);
t132 = -qJD(5) + t226;
t272 = t51 * t132;
t222 = t161 * qJDD(1);
t270 = t163 * t225 + t222;
t246 = pkin(1) * qJDD(1);
t223 = t155 * qJDD(3);
t73 = t270 * t157 + t223;
t181 = t104 * t231 + t161 * t73;
t269 = t155 * t181;
t212 = qJD(5) * t264;
t230 = qJD(5) * t160;
t268 = -t155 * t230 + t157 * t212;
t208 = t161 * t225;
t267 = t208 - t148;
t138 = t157 * pkin(4) + pkin(3);
t255 = pkin(7) + qJ(4);
t190 = t163 * t138 + t161 * t255;
t227 = t161 * qJD(2);
t89 = t163 * t123 + t227;
t76 = qJD(3) * qJ(4) + t89;
t195 = t163 * pkin(3) + t161 * qJ(4);
t186 = -pkin(2) - t195;
t158 = cos(pkin(8));
t257 = t158 * pkin(1);
t99 = t186 - t257;
t79 = t99 * qJD(1);
t24 = -t155 * t76 + t157 * t79;
t20 = -pkin(4) * t226 - t104 * pkin(7) + t24;
t25 = t155 * t79 + t157 * t76;
t22 = t103 * pkin(7) + t25;
t184 = t160 * t22 - t264 * t20;
t111 = t161 * t123;
t221 = t161 * qJDD(2) + t271 * t163;
t30 = qJDD(3) * qJ(4) + (qJD(4) - t111) * qJD(3) + t221;
t194 = pkin(3) * t161 - qJ(4) * t163;
t93 = qJD(3) * t194 - t161 * qJD(4);
t38 = qJD(1) * t93 + qJDD(1) * t99;
t13 = -t155 * t30 + t157 * t38;
t8 = t267 * pkin(4) - t73 * pkin(7) + t13;
t14 = t155 * t38 + t157 * t30;
t72 = -t157 * qJDD(3) + t270 * t155;
t9 = -t72 * pkin(7) + t14;
t1 = -qJD(5) * t184 + t160 * t8 + t264 * t9;
t265 = t72 * pkin(4);
t263 = pkin(4) * t155;
t262 = g(1) * t143;
t259 = g(2) * t145;
t258 = g(3) * t161;
t256 = t54 * t51;
t203 = t160 * t73 + t264 * t72;
t16 = qJD(5) * t54 + t203;
t219 = t264 * t157;
t108 = t264 * t155 + t160 * t157;
t96 = t108 * qJD(5);
t39 = t160 * t163 * t229 + t161 * t96 - t219 * t231;
t180 = -t160 * t155 + t219;
t85 = t180 * t161;
t254 = -t85 * t16 + t39 * t51;
t106 = qJDD(5) + t267;
t178 = t163 * t108;
t40 = qJD(3) * t178 + t268 * t161;
t84 = t108 * t161;
t253 = -t84 * t106 + t40 * t132;
t239 = t157 * t163;
t187 = pkin(4) * t161 - pkin(7) * t239;
t113 = t194 * qJD(1);
t88 = t163 * qJD(2) - t111;
t45 = t157 * t113 - t155 * t88;
t23 = qJD(1) * t187 + t45;
t217 = t155 * t226;
t46 = t155 * t113 + t157 * t88;
t29 = -pkin(7) * t217 + t46;
t119 = t255 * t155;
t120 = t255 * t157;
t59 = -t264 * t119 - t160 * t120;
t252 = qJD(4) * t180 + qJD(5) * t59 - t160 * t23 - t264 * t29;
t60 = -t160 * t119 + t264 * t120;
t251 = -qJD(4) * t108 - qJD(5) * t60 + t160 * t29 - t264 * t23;
t250 = qJD(1) * t178 - t96;
t249 = t180 * t226 - t268;
t215 = t103 * t231;
t240 = t157 * t161;
t248 = t157 * t215 - t72 * t240;
t247 = t163 * t72;
t232 = qJD(3) * t161;
t216 = t137 * t232;
t48 = t155 * t216 + t157 * t93;
t57 = t137 * t239 + t155 * t99;
t245 = t103 * t161;
t244 = t143 * t163;
t243 = t145 * t163;
t166 = qJD(1) ^ 2;
t242 = t154 * t166;
t241 = t155 * t163;
t236 = (t154 * t225 + t207) * t155;
t139 = -pkin(2) - t257;
t124 = qJD(1) * t139;
t234 = qJD(1) * t155;
t122 = qJDD(1) * t139;
t164 = cos(qJ(1));
t220 = t164 * pkin(1) + t145 * pkin(2) + t143 * pkin(6);
t218 = t153 * t234;
t213 = t104 * t226;
t210 = t155 * t148;
t209 = t157 * t148;
t162 = sin(qJ(1));
t206 = -t162 * pkin(1) + t145 * pkin(6);
t205 = t137 + t263;
t15 = -t103 * t212 + t104 * t230 + t160 * t72 - t264 * t73;
t204 = t15 * t163 + t54 * t232;
t202 = t104 * t232 - t163 * t73;
t200 = t163 * t208;
t198 = -t259 + t262;
t197 = g(1) * t162 - g(2) * t164;
t196 = -t84 * t15 + t54 * t40;
t193 = -t85 * t106 - t39 * t132;
t192 = -t13 * t155 + t14 * t157;
t191 = -t24 * t155 + t25 * t157;
t188 = qJD(1) * (-t103 + t228);
t87 = t157 * t99;
t37 = -pkin(7) * t240 + t87 + (-t137 * t155 - pkin(4)) * t163;
t44 = -t155 * t161 * pkin(7) + t57;
t17 = -t160 * t44 + t264 * t37;
t6 = t160 * t20 + t264 * t22;
t18 = t160 * t37 + t264 * t44;
t183 = t163 * t16 - t232 * t51;
t179 = -qJD(1) * t124 + t199;
t74 = -qJD(3) * pkin(3) + qJD(4) - t88;
t177 = -qJ(4) * t232 + (qJD(4) - t74) * t163;
t165 = qJD(3) ^ 2;
t176 = t137 * t165 + 0.2e1 * t122 + t259;
t175 = 0.2e1 * qJD(3) * t124 - qJDD(3) * t137;
t172 = -t163 * t199 - t258;
t2 = -qJD(5) * t6 - t160 * t9 + t264 * t8;
t41 = -t123 * t232 + t221;
t169 = -t42 * t161 + t41 * t163 + (-t161 * t89 - t163 * t88) * qJD(3);
t151 = pkin(9) + qJ(5);
t144 = cos(t151);
t142 = sin(t151);
t131 = t161 * t166 * t163;
t127 = t161 * t262;
t116 = qJDD(3) * t163 - t165 * t161;
t115 = qJDD(3) * t161 + t165 * t163;
t100 = t154 * qJDD(1) - 0.2e1 * t200;
t92 = t205 * t161;
t83 = t205 * t231;
t81 = t155 * t93;
t71 = t143 * t142 + t144 * t243;
t70 = -t142 * t243 + t143 * t144;
t69 = t145 * t142 - t144 * t244;
t68 = t142 * t244 + t145 * t144;
t61 = t227 + (pkin(4) * t234 + t123) * t163;
t56 = -t137 * t241 + t87;
t49 = -t157 * t216 + t81;
t47 = -t103 * pkin(4) + t74;
t36 = t81 + (-pkin(7) * t241 - t137 * t240) * qJD(3);
t28 = qJD(3) * t187 + t48;
t21 = t35 + t265;
t4 = -t18 * qJD(5) - t160 * t36 + t264 * t28;
t3 = t17 * qJD(5) + t160 * t28 + t264 * t36;
t5 = [0, 0, 0, 0, 0, qJDD(1), t197, g(1) * t164 + g(2) * t162, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t158 * t246 + t198, -0.2e1 * t156 * t246 + t199, 0, (t197 + (t156 ^ 2 + t158 ^ 2) * t246) * pkin(1), t153 * qJDD(1) + 0.2e1 * t200, -0.2e1 * t266, t115, t100, t116, 0, t175 * t161 + (-t176 + t262) * t163, t161 * t176 + t163 * t175 - t127, (t153 + t154) * t121 + t169 - t199, t122 * t139 - g(1) * (-t143 * pkin(2) + t206) - g(2) * t220 + t169 * t137, t181 * t157, t248 - t269, t202 + t275, (t161 * t72 - t215) * t155, t247 + (-t218 + t245) * qJD(3) + t236, t100, -t199 * t155 + (t137 * t72 + t35 * t155 + (qJD(1) * t56 + t24) * qJD(3)) * t161 + (-t48 * qJD(1) - t56 * qJDD(1) - t13 + t198 * t157 + (-t103 * t137 + t155 * t74) * qJD(3)) * t163, -t199 * t157 + (t137 * t73 + t35 * t157 + (-qJD(1) * t57 - t25) * qJD(3)) * t161 + (t49 * qJD(1) + t57 * qJDD(1) + t14 - t198 * t155 + (t104 * t137 + t157 * t74) * qJD(3)) * t163, t49 * t103 - t48 * t104 - t56 * t73 - t57 * t72 + t127 + (-t155 * t25 - t157 * t24) * t231 + (-t13 * t157 - t14 * t155 - t259) * t161, t14 * t57 + t25 * t49 + t13 * t56 + t24 * t48 - g(1) * t206 - g(2) * (t145 * t195 + t220) - t186 * t262 + (t35 * t161 + t231 * t74) * t137, -t15 * t85 - t54 * t39, -t196 + t254, -t193 + t204, t16 * t84 + t51 * t40, t183 + t253, -t106 * t163 - t132 * t232, -g(1) * t69 - g(2) * t71 + t17 * t106 - t4 * t132 + t92 * t16 - t2 * t163 - t184 * t232 + t21 * t84 + t47 * t40 + t83 * t51, -g(1) * t68 - g(2) * t70 + t1 * t163 - t18 * t106 + t3 * t132 - t92 * t15 + t21 * t85 - t232 * t6 - t47 * t39 + t83 * t54, -t1 * t84 + t17 * t15 - t18 * t16 - t161 * t259 - t184 * t39 - t2 * t85 - t3 * t51 - t4 * t54 - t6 * t40 + t127, t1 * t18 + t6 * t3 + t2 * t17 - t184 * t4 + t21 * t92 + t47 * t83 - g(1) * (t145 * t263 + t206) - g(2) * (t190 * t145 + t220) + (-g(1) * (-pkin(2) - t190) - g(2) * t263) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t116, -t115, 0, t41 * t161 + t42 * t163 - g(3) + (-t161 * t88 + t163 * t89) * qJD(3), 0, 0, 0, 0, 0, 0, -t247 + (-t218 - t245) * qJD(3) + t236, t202 - t275, t248 + t269, -t35 * t163 - g(3) + t192 * t161 + (t161 * t74 + t163 * t191) * qJD(3), 0, 0, 0, 0, 0, 0, -t183 + t253, t193 + t204, t196 + t254, t1 * t85 - t21 * t163 + t184 * t40 - t2 * t84 + t232 * t47 - t6 * t39 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, t235 * t166, t222, t131, t148, qJDD(3), t89 * qJD(3) + t161 * t179 - t150 + t42, t258 + (t88 + t111) * qJD(3) + t179 * t163 - t221, 0, 0, t73 * t155 - t157 * t213, -t155 * t72 + t73 * t157 + (-t103 * t157 + t104 * t155) * t226, -t210 + t157 * t242 + (-t104 + t229) * t233, t103 * t217 - t72 * t157, -t155 * t242 + t161 * t188 - t209, t131, qJ(4) * t210 - pkin(3) * t72 + t89 * t103 - t168 * t157 + (t155 * t177 - t24 * t161 + t163 * t45) * qJD(1), qJ(4) * t209 - pkin(3) * t73 - t89 * t104 + t168 * t155 + (t157 * t177 + t25 * t161 - t163 * t46) * qJD(1), -t46 * t103 + t45 * t104 + (-qJ(4) * t72 + qJD(4) * t103 + t226 * t24 + t14) * t157 + (qJ(4) * t73 + qJD(4) * t104 + t226 * t25 - t13) * t155 + t172, -t24 * t45 - t25 * t46 - t74 * t89 + t191 * qJD(4) - t168 * pkin(3) + (t172 + t192) * qJ(4), -t15 * t108 - t249 * t54, -t108 * t16 - t15 * t180 + t249 * t51 + t250 * t54, t108 * t106 + t132 * t249 - t233 * t54, -t16 * t180 - t250 * t51, t106 * t180 - t132 * t250 + t233 * t51, t132 * t233, t59 * t106 - t132 * t251 - t138 * t16 + t144 * t174 - t180 * t21 + t184 * t233 - t250 * t47 - t61 * t51, -t60 * t106 + t21 * t108 + t132 * t252 + t138 * t15 - t142 * t174 + t233 * t6 - t249 * t47 - t61 * t54, t1 * t180 - t2 * t108 + t59 * t15 - t60 * t16 - t184 * t249 + t250 * t6 - t251 * t54 - t252 * t51 + t172, -g(3) * t190 + t1 * t60 - t21 * t138 + t2 * t59 - t251 * t184 + t252 * t6 - t47 * t61 + t199 * (t138 * t161 - t163 * t255); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 - t213, t157 * t222 + t163 * t188 + t223, -t103 ^ 2 - t104 ^ 2, -t25 * t103 + t24 * t104 + t168, 0, 0, 0, 0, 0, 0, -t54 * t132 + t16, -t15 + t272, -t273 - t274, -t184 * t54 + t51 * t6 + t168 + t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t273 - t274, -t15 - t272, -t256, -t203 + (-qJD(5) - t132) * t54, t106, -g(1) * t70 + g(2) * t68 - t6 * t132 + t142 * t258 - t47 * t54 + t2, g(1) * t71 - g(2) * t69 + t132 * t184 + t144 * t258 + t47 * t51 - t1, 0, 0;];
tau_reg = t5;
