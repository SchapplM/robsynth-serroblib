% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR12
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:13:13
% EndTime: 2019-12-31 19:13:20
% DurationCPUTime: 3.14s
% Computational Cost: add. (5066->399), mult. (9678->505), div. (0->0), fcn. (6238->10), ass. (0->217)
t149 = sin(qJ(3));
t152 = cos(qJ(4));
t148 = sin(qJ(4));
t153 = cos(qJ(3));
t246 = t148 * t153;
t84 = t149 * t152 + t246;
t79 = t84 * qJD(1);
t280 = qJD(5) + t79;
t147 = sin(qJ(5));
t151 = cos(qJ(5));
t222 = qJD(3) + qJD(4);
t201 = t151 * t222;
t234 = qJD(1) * t153;
t235 = qJD(1) * t149;
t81 = -t148 * t235 + t152 * t234;
t62 = t147 * t81 - t201;
t283 = t280 * t62;
t158 = qJD(1) ^ 2;
t154 = cos(qJ(1));
t138 = g(2) * t154;
t150 = sin(qJ(1));
t239 = g(1) * t150 - t138;
t282 = -t158 * qJ(2) - t239;
t156 = -pkin(1) - pkin(6);
t107 = t156 * qJD(1) + qJD(2);
t233 = qJD(3) * t149;
t223 = t153 * qJDD(1);
t227 = qJD(1) * qJD(3);
t278 = t149 * t227 - t223;
t105 = t156 * qJDD(1) + qJDD(2);
t86 = t153 * t105;
t52 = qJDD(3) * pkin(3) + t278 * pkin(7) - t107 * t233 + t86;
t232 = qJD(3) * t153;
t209 = t153 * t227;
t224 = t149 * qJDD(1);
t279 = t209 + t224;
t54 = -t279 * pkin(7) + t149 * t105 + t107 * t232;
t206 = t148 * t54 - t152 * t52;
t70 = (-pkin(7) * qJD(1) + t107) * t149;
t256 = t152 * t70;
t71 = -pkin(7) * t234 + t153 * t107;
t68 = qJD(3) * pkin(3) + t71;
t45 = t148 * t68 + t256;
t15 = -t45 * qJD(4) - t206;
t221 = qJDD(3) + qJDD(4);
t13 = -t221 * pkin(4) - t15;
t146 = qJ(3) + qJ(4);
t133 = cos(t146);
t249 = t133 * t150;
t219 = g(1) * t249;
t281 = t13 + t219;
t132 = sin(t146);
t172 = g(3) * t132 + t133 * t138;
t42 = t222 * pkin(8) + t45;
t92 = pkin(3) * t235 + qJD(1) * qJ(2);
t48 = pkin(4) * t79 - pkin(8) * t81 + t92;
t20 = -t147 * t42 + t151 * t48;
t21 = t147 * t48 + t151 * t42;
t178 = t147 * t21 + t151 * t20;
t141 = qJDD(1) * qJ(2);
t226 = qJD(1) * qJD(4);
t277 = -t149 * t226 - t278;
t189 = g(1) * t154 + g(2) * t150;
t276 = t189 * t133;
t144 = t149 ^ 2;
t145 = t153 ^ 2;
t237 = t144 + t145;
t203 = t237 * t105;
t120 = t133 * pkin(8);
t190 = pkin(4) * t132 - t120;
t142 = qJD(1) * qJD(2);
t217 = 0.2e1 * t142;
t275 = 0.2e1 * t141 + t217 - t189;
t200 = t222 * t153;
t231 = qJD(4) * t148;
t59 = -t148 * t233 - t149 * t231 + t152 * t200;
t274 = -t221 * t84 - t222 * t59;
t122 = g(3) * t133;
t271 = g(3) * t149;
t136 = t149 * pkin(3);
t230 = qJD(4) * t152;
t14 = t148 * t52 + t152 * t54 + t68 * t230 - t70 * t231;
t12 = t221 * pkin(8) + t14;
t166 = -t279 * t148 + t277 * t152 - t226 * t246;
t171 = t277 * t148;
t40 = (qJD(1) * t200 + t224) * t152 + t171;
t69 = t279 * pkin(3) + t141 + t142;
t17 = pkin(4) * t40 - pkin(8) * t166 + t69;
t2 = qJD(5) * t20 + t151 * t12 + t147 * t17;
t1 = t2 * t151;
t260 = t148 * t70;
t44 = t152 * t68 - t260;
t41 = -t222 * pkin(4) - t44;
t270 = t41 * t79;
t64 = t147 * t222 + t151 * t81;
t269 = t64 * t62;
t268 = t280 * t81;
t267 = t81 * t79;
t266 = pkin(7) - t156;
t49 = t148 * t71 + t256;
t265 = t45 - t49;
t58 = t222 * t84;
t85 = -t148 * t149 + t152 * t153;
t264 = t85 * t221 - t58 * t222;
t38 = qJDD(5) + t40;
t262 = t147 * t38;
t261 = t147 * t62;
t258 = t151 * t38;
t257 = t151 * t64;
t229 = qJD(5) * t147;
t22 = -qJD(5) * t201 - t147 * t221 - t151 * t166 + t81 * t229;
t255 = t22 * t147;
t205 = t147 * t166 - t151 * t221;
t23 = t64 * qJD(5) + t205;
t254 = t23 * t151;
t253 = pkin(1) * qJDD(1);
t252 = qJD(1) * t92;
t251 = qJD(5) * t280;
t250 = t132 * t150;
t248 = t147 * t150;
t247 = t147 * t154;
t245 = t150 * t151;
t244 = t151 * t154;
t119 = qJ(2) + t136;
t241 = (t217 + t141) * qJ(2);
t240 = t154 * pkin(1) + t150 * qJ(2);
t238 = t144 - t145;
t157 = qJD(3) ^ 2;
t236 = -t157 - t158;
t228 = qJD(5) * t151;
t108 = pkin(3) * t232 + qJD(2);
t225 = qJDD(3) * t149;
t220 = pkin(8) * t251;
t102 = t132 * t138;
t218 = pkin(3) * t230;
t124 = pkin(3) * t148 + pkin(8);
t216 = t124 * t251;
t215 = t85 * t229;
t214 = t85 * t228;
t213 = t153 * t158 * t149;
t35 = t41 * t229;
t36 = t41 * t228;
t212 = g(1) * t250 - t102 + t122;
t211 = g(1) * (pkin(4) * t249 + pkin(8) * t250);
t91 = t266 * t153;
t135 = t154 * qJ(2);
t207 = -t150 * pkin(1) + t135;
t204 = t151 * t280;
t202 = qJD(5) * t84 + qJD(1);
t199 = t237 * qJDD(1);
t198 = qJDD(2) - t253;
t196 = t149 * t209;
t195 = pkin(3) * t231 - t49;
t55 = pkin(4) * t81 + pkin(8) * t79;
t194 = t266 * t233;
t193 = -pkin(8) * t38 + t270;
t192 = t2 * t84 + t21 * t59;
t16 = t151 * t17;
t3 = -qJD(5) * t21 - t147 * t12 + t16;
t191 = -t20 * t59 - t3 * t84;
t187 = t13 * t85 - t41 * t58;
t186 = -t85 * t22 - t58 * t64;
t185 = -t22 * t84 + t64 * t59;
t184 = t85 * t23 - t58 * t62;
t183 = -t23 * t84 - t62 * t59;
t182 = -t280 * t58 + t38 * t85;
t181 = t166 * t85 - t58 * t81;
t180 = t40 * t84 + t59 * t79;
t177 = t147 * t20 - t151 * t21;
t56 = pkin(4) * t84 - pkin(8) * t85 + t119;
t90 = t266 * t149;
t61 = -t148 * t91 - t152 * t90;
t29 = -t147 * t61 + t151 * t56;
t30 = t147 * t56 + t151 * t61;
t60 = -t148 * t90 + t152 * t91;
t176 = t281 * t147 + t21 * t81 + t36;
t175 = t172 * t151 - t20 * t81 + t35;
t155 = -pkin(7) - pkin(6);
t174 = t154 * t136 + t150 * t155 + t207;
t173 = -t178 * t79 + t1 - t212;
t170 = t150 * t136 - t154 * t155 + t240;
t169 = -qJD(5) * t48 - t12 + t122;
t168 = 0.2e1 * qJ(2) * t227 + qJDD(3) * t156;
t165 = t92 * t79 - t14 + t212;
t164 = -t124 * t38 - t218 * t280 + t270;
t163 = -t178 * qJD(5) - t3 * t147;
t162 = -t92 * t81 + t172 - t206 - t219;
t161 = t14 * t84 + t15 * t85 - t44 * t58 + t45 * t59 - t239;
t160 = t163 + t1;
t159 = -t156 * t157 + t275;
t131 = qJDD(3) * t153;
t125 = -pkin(3) * t152 - pkin(4);
t83 = qJD(3) * t91;
t75 = t132 * t244 - t248;
t74 = t132 * t247 + t245;
t73 = t132 * t245 + t247;
t72 = -t132 * t248 + t244;
t53 = pkin(3) * t234 + t55;
t50 = t152 * t71 - t260;
t43 = -t79 ^ 2 + t81 ^ 2;
t34 = t61 * qJD(4) - t148 * t83 - t152 * t194;
t33 = -t60 * qJD(4) + t148 * t194 - t152 * t83;
t32 = t81 * t222 + (-t222 * t234 - t224) * t152 - t171;
t31 = t79 * t222 + t166;
t28 = pkin(4) * t59 + pkin(8) * t58 + t108;
t27 = t147 * t55 + t151 * t44;
t26 = -t147 * t44 + t151 * t55;
t25 = t147 * t53 + t151 * t50;
t24 = -t147 * t50 + t151 * t53;
t10 = t204 * t280 - t64 * t81 + t262;
t9 = -t147 * t280 ^ 2 + t62 * t81 + t258;
t8 = t261 * t280 - t254;
t7 = t64 * t204 - t255;
t6 = -t30 * qJD(5) - t147 * t33 + t151 * t28;
t5 = t29 * qJD(5) + t147 * t28 + t151 * t33;
t4 = (-t22 - t283) * t151 + (-t280 * t64 - t23) * t147;
t11 = [0, 0, 0, 0, 0, qJDD(1), t239, t189, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t239 - 0.2e1 * t253, t275, -t198 * pkin(1) - g(1) * t207 - g(2) * t240 + t241, qJDD(1) * t145 - 0.2e1 * t196, -0.2e1 * t149 * t223 + 0.2e1 * t238 * t227, -t149 * t157 + t131, qJDD(1) * t144 + 0.2e1 * t196, -t153 * t157 - t225, 0, t149 * t159 + t153 * t168, -t149 * t168 + t153 * t159, -t156 * t199 - t203 + t239, -g(1) * (t156 * t150 + t135) - g(2) * (pkin(6) * t154 + t240) + t156 * t203 + t241, t181, -t166 * t84 - t40 * t85 + t58 * t79 - t59 * t81, t264, t180, t274, 0, t108 * t79 + t119 * t40 - t132 * t189 - t221 * t60 - t222 * t34 + t92 * t59 + t69 * t84, t108 * t81 + t119 * t166 - t221 * t61 - t222 * t33 - t92 * t58 + t69 * t85 - t276, t166 * t60 - t33 * t79 + t34 * t81 - t40 * t61 - t161, -g(1) * t174 - g(2) * t170 + t92 * t108 + t69 * t119 + t14 * t61 - t15 * t60 + t45 * t33 - t44 * t34, t151 * t186 - t215 * t64, (t147 * t64 + t151 * t62) * t58 + (t255 - t254 + (-t257 + t261) * qJD(5)) * t85, t151 * t182 - t215 * t280 + t185, t147 * t184 + t214 * t62, -t147 * t182 - t214 * t280 + t183, t280 * t59 + t38 * t84, -g(1) * t75 - g(2) * t73 + t147 * t187 + t60 * t23 + t280 * t6 + t29 * t38 + t34 * t62 + t36 * t85 - t191, g(1) * t74 - g(2) * t72 + t151 * t187 - t60 * t22 - t280 * t5 - t30 * t38 + t34 * t64 - t35 * t85 - t192, t29 * t22 - t30 * t23 - t5 * t62 - t6 * t64 + t178 * t58 + t276 + (qJD(5) * t177 - t147 * t2 - t151 * t3) * t85, t2 * t30 + t21 * t5 + t3 * t29 + t20 * t6 + t13 * t60 + t41 * t34 - g(1) * (t154 * t190 + t174) - g(2) * (t150 * t190 + t170); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t158, t282 + t198, 0, 0, 0, 0, 0, 0, t236 * t149 + t131, t236 * t153 - t225, -t199, t203 + t282, 0, 0, 0, 0, 0, 0, -qJD(1) * t79 + t264, -qJD(1) * t81 + t274, -t180 - t181, t161 - t252, 0, 0, 0, 0, 0, 0, -t84 * t262 + (-t147 * t59 - t151 * t202) * t280 - t184, -t84 * t258 + (t147 * t202 - t151 * t59) * t280 - t186, (t202 * t64 + t183) * t151 + (t202 * t62 + t185) * t147, (-t20 * t202 + t192) * t151 + (-t202 * t21 + t191) * t147 - t187 - t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t238 * t158, t223, -t213, -t224, qJDD(3), t153 * t282 + t271 + t86, g(3) * t153 + (-t105 - t282) * t149, 0, 0, t267, t43, t31, -t267, t32, t221, t49 * qJD(3) - t265 * qJD(4) + (t152 * t221 - t222 * t231 - t234 * t79) * pkin(3) + t162, t50 * t222 + (-t148 * t221 - t222 * t230 - t234 * t81) * pkin(3) + t165, t265 * t81 + (-t44 + t50) * t79 + (-t148 * t40 - t152 * t166 + (t148 * t81 - t152 * t79) * qJD(4)) * pkin(3), t44 * t49 - t45 * t50 + (t271 + t14 * t148 + t15 * t152 + (-t148 * t44 + t152 * t45) * qJD(4) + (-t239 - t252) * t153) * pkin(3), t7, t4, t10, t8, t9, -t268, t125 * t23 - t24 * t280 + t195 * t62 + (-t281 - t216) * t151 + t164 * t147 + t175, -t125 * t22 + t25 * t280 + t195 * t64 + t164 * t151 + (-t172 + t216) * t147 + t176, t24 * t64 + t25 * t62 + (-t62 * t218 - t124 * t23 + (t124 * t64 - t20) * qJD(5)) * t151 + (t64 * t218 - t124 * t22 - t3 + (t124 * t62 - t21) * qJD(5)) * t147 + t173, t13 * t125 - t21 * t25 - t20 * t24 - t41 * t49 - t211 - g(3) * (-t136 - t190) - (-pkin(4) * t133 - pkin(8) * t132) * t138 + (-t239 * t153 + (t148 * t41 - t152 * t177) * qJD(4)) * pkin(3) + t160 * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, t43, t31, -t267, t32, t221, t45 * qJD(3) + t162, t222 * t44 + t165, 0, 0, t7, t4, t10, t8, t9, -t268, -pkin(4) * t23 - t26 * t280 - t45 * t62 + t193 * t147 + (-t281 - t220) * t151 + t175, pkin(4) * t22 + t27 * t280 - t45 * t64 + t193 * t151 + (-t172 + t220) * t147 + t176, t26 * t64 + t27 * t62 + (-t255 - t254 + (t257 + t261) * qJD(5)) * pkin(8) + t163 + t173, -t21 * t27 - t20 * t26 - t41 * t45 - t211 - g(3) * t120 + (-t13 + t172) * pkin(4) + (t160 + t102) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t269, -t62 ^ 2 + t64 ^ 2, -t22 + t283, -t269, -t205 + (-qJD(5) + t280) * t64, t38, -g(1) * t72 - g(2) * t74 + t147 * t169 + t21 * t280 - t228 * t42 - t41 * t64 + t16, g(1) * t73 - g(2) * t75 + t20 * t280 + t41 * t62 + (qJD(5) * t42 - t17) * t147 + t169 * t151, 0, 0;];
tau_reg = t11;
