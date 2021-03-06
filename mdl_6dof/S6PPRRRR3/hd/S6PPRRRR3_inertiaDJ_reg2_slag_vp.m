% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PPRRRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_inertiaDJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:11:41
% EndTime: 2019-03-08 19:11:55
% DurationCPUTime: 5.68s
% Computational Cost: add. (6943->458), mult. (21914->829), div. (0->0), fcn. (24067->16), ass. (0->219)
t103 = sin(qJ(5));
t104 = sin(qJ(4));
t108 = cos(qJ(4));
t225 = cos(pkin(8));
t180 = t108 * t225;
t97 = sin(pkin(8));
t234 = t104 * t97;
t130 = -pkin(3) * t180 + pkin(10) * t234;
t64 = -t225 * pkin(4) + t130;
t107 = cos(qJ(5));
t201 = t103 * t234;
t69 = -t107 * t225 + t201;
t229 = t107 * t97;
t70 = t103 * t225 + t104 * t229;
t113 = t69 * pkin(5) - t70 * pkin(12) + t64;
t181 = t104 * t225;
t169 = pkin(3) * t181;
t128 = t225 * pkin(11) + t169;
t126 = qJD(5) * t128;
t65 = t130 * qJD(4);
t115 = -t103 * t126 - t107 * t65;
t163 = pkin(4) * t104 - pkin(11) * t108;
t193 = -pkin(11) * t104 - pkin(3);
t153 = pkin(4) * t108 - t193;
t251 = pkin(10) * t108;
t123 = t103 * t251 + t107 * t153;
t258 = t123 * qJD(5);
t261 = -(-t258 + (t104 * pkin(12) + t103 * t163) * qJD(4)) * t97 - t115 - qJD(6) * t113;
t105 = sin(qJ(3));
t220 = qJD(3) * t105;
t98 = sin(pkin(7));
t194 = t98 * t220;
t173 = t97 * t194;
t252 = cos(qJ(3));
t165 = t225 * t252;
t151 = t98 * t165;
t137 = t108 * t151;
t226 = cos(pkin(7));
t184 = t226 * t97;
t167 = qJD(4) * t184;
t199 = t98 * t252;
t170 = qJD(3) * t199;
t200 = t104 * t105 * t98;
t33 = -qJD(4) * t137 + (t225 * qJD(3) + qJD(4)) * t200 + (-t167 - t170) * t108;
t129 = t103 * t173 - t107 * t33;
t224 = t105 * t108;
t53 = t98 * t224 + (t151 + t184) * t104;
t68 = -t97 * t199 + t225 * t226;
t36 = t103 * t53 - t107 * t68;
t114 = -t36 * qJD(5) + t129;
t213 = qJD(5) * t108;
t219 = qJD(4) * t104;
t260 = t103 * t213 + t107 * t219;
t216 = qJD(5) * t103;
t218 = qJD(4) * t108;
t195 = t97 * t218;
t221 = t70 * qJD(5);
t54 = t103 * t195 + t221;
t147 = -t107 * t54 + t69 * t216;
t102 = sin(qJ(6));
t92 = t102 ^ 2;
t106 = cos(qJ(6));
t94 = t106 ^ 2;
t239 = t92 - t94;
t183 = qJD(6) * t239;
t136 = t103 * t153;
t259 = (t107 * t251 - t136) * qJD(5);
t177 = t225 * qJD(5);
t131 = t177 + t195;
t86 = qJD(5) * t201;
t257 = -t107 * t131 + t86;
t101 = cos(pkin(6));
t100 = cos(pkin(14));
t182 = t100 * t226;
t96 = sin(pkin(14));
t99 = sin(pkin(6));
t240 = t96 * t99;
t51 = (t101 * t98 + t99 * t182) * t105 + t252 * t240;
t256 = 0.2e1 * t97;
t255 = t108 ^ 2;
t254 = pkin(11) * t97;
t150 = t252 * t182;
t50 = t101 * t199 + (-t105 * t96 + t150) * t99;
t67 = -t99 * t100 * t98 + t101 * t226;
t139 = t225 * t50 + t67 * t97;
t22 = t104 * t139 + t51 * t108;
t35 = t67 * t225 - t50 * t97;
t13 = t103 * t22 - t35 * t107;
t235 = t104 * t51;
t46 = -t99 * qJD(3) * t150 - t101 * t170 + t220 * t240;
t47 = t51 * qJD(3);
t12 = -t47 * t181 - t46 * t108 + (t108 * t139 - t235) * qJD(4);
t14 = t103 * t35 + t107 * t22;
t3 = qJD(5) * t14 + t103 * t12 - t47 * t229;
t253 = t13 * t3;
t250 = pkin(11) * t102;
t249 = t103 * pkin(11);
t248 = t103 * t3;
t11 = t22 * qJD(4) - t46 * t104 + t47 * t180;
t228 = t108 * t97;
t21 = -t50 * t180 - t67 * t228 + t235;
t247 = t11 * t21;
t206 = pkin(10) * t228;
t73 = t169 + t206;
t66 = t73 * qJD(4);
t246 = t21 * t66;
t37 = t103 * t68 + t107 * t53;
t16 = qJD(5) * t37 - t103 * t33 - t107 * t173;
t245 = t36 * t16;
t91 = t97 ^ 2;
t244 = t47 * t91;
t243 = t47 * t97;
t34 = t104 * t167 + ((t104 * t165 + t224) * qJD(4) + (t252 * t104 + t105 * t180) * qJD(3)) * t98;
t52 = -t108 * t184 - t137 + t200;
t242 = t52 * t34;
t241 = t52 * t66;
t43 = -t97 * t136 + t107 * (t128 + t206);
t93 = t103 ^ 2;
t238 = -t107 ^ 2 + t93;
t196 = t97 * t219;
t55 = t102 * t70 + t106 * t228;
t28 = qJD(6) * t55 - t102 * t196 + t106 * t257;
t237 = t102 * t28;
t236 = t102 * t55;
t233 = t106 * t28;
t202 = t102 * t228;
t211 = qJD(6) * t106;
t29 = -qJD(6) * t202 - t102 * t257 - t106 * t196 + t70 * t211;
t232 = t106 * t29;
t56 = t106 * t70 - t202;
t231 = t106 * t56;
t227 = t16 * t103;
t223 = t106 * t107;
t217 = qJD(5) * t102;
t215 = qJD(5) * t106;
t214 = qJD(5) * t107;
t212 = qJD(6) * t102;
t210 = qJD(6) * t107;
t209 = 0.2e1 * t69 * t54;
t208 = -0.2e1 * pkin(4) * qJD(5);
t207 = -0.2e1 * pkin(5) * qJD(6);
t205 = t107 * t250;
t204 = pkin(11) * t223;
t203 = pkin(11) * t214;
t198 = t91 * t218;
t191 = t106 * t214;
t190 = t102 * t210;
t189 = t103 * t211;
t188 = t106 * t210;
t187 = t102 * t211;
t186 = t103 * t214;
t179 = t238 * qJD(5);
t178 = qJD(4) * t225;
t176 = 0.2e1 * t186;
t174 = t91 * t194;
t172 = t93 * t187;
t171 = t104 * t198;
t168 = t106 * t186;
t166 = t13 * t16 + t3 * t36;
t164 = t104 * t178;
t162 = -pkin(5) * t107 - pkin(12) * t103;
t161 = pkin(5) * t103 - pkin(12) * t107;
t160 = t11 * t52 + t21 * t34;
t10 = t102 * t21 + t106 * t14;
t9 = -t102 * t14 + t106 * t21;
t159 = -t10 * t102 - t106 * t9;
t39 = -pkin(12) * t228 + t43;
t17 = -t102 * t39 + t106 * t113;
t18 = t102 * t113 + t106 * t39;
t158 = -t102 * t18 - t106 * t17;
t23 = -t102 * t37 + t106 * t52;
t24 = t102 * t52 + t106 * t37;
t157 = -t102 * t24 - t106 * t23;
t156 = -t102 * t56 - t106 * t55;
t152 = pkin(4) - t162;
t135 = t106 * t152;
t61 = -t135 - t205;
t62 = -t102 * t152 + t204;
t155 = -t102 * t62 - t106 * t61;
t138 = qJD(4) * t163;
t25 = (-t103 * t138 + t258) * t97 - t115;
t116 = t103 * t65 - t107 * t126;
t26 = (t107 * t138 - t259) * t97 + t116;
t154 = -t26 * t103 - t25 * t107;
t149 = t102 * t3 + t13 * t211;
t148 = -t106 * t3 + t13 * t212;
t146 = t102 * t16 + t36 * t211;
t145 = -t106 * t16 + t36 * t212;
t20 = (t259 + (-t104 * pkin(5) - t107 * t163) * qJD(4)) * t97 - t116;
t127 = t103 * t128;
t38 = t127 + (-t107 * t193 + (t107 * pkin(4) + pkin(10) * t103 + pkin(5)) * t108) * t97;
t144 = t102 * t20 + t38 * t211;
t143 = -t106 * t20 + t38 * t212;
t142 = t102 * t54 + t69 * t211;
t141 = -t106 * t54 + t69 * t212;
t140 = t161 * t102;
t133 = t103 * t219 - t107 * t213;
t132 = t103 * t215 + t190;
t4 = -qJD(5) * t13 + t103 * t243 + t107 * t12;
t1 = -qJD(6) * t10 - t102 * t4 + t106 * t11;
t2 = qJD(6) * t9 + t102 * t11 + t106 * t4;
t121 = qJD(6) * t159 - t1 * t102 + t106 * t2;
t111 = pkin(3) * t164 + t54 * pkin(5) + pkin(10) * t195 + pkin(12) * t257;
t5 = -t102 * t111 + t261 * t106 + t39 * t212;
t6 = t261 * t102 + t106 * t111 - t39 * t211;
t120 = qJD(6) * t158 - t102 * t6 - t106 * t5;
t7 = -t102 * t34 - t106 * t114 - t52 * t211 + t37 * t212;
t8 = -qJD(6) * t24 - t102 * t114 + t106 * t34;
t119 = qJD(6) * t157 - t102 * t8 - t106 * t7;
t118 = t248 + t107 * t4 + (-t103 * t14 + t107 * t13) * qJD(5);
t40 = pkin(11) * t132 - qJD(5) * t140 + qJD(6) * t135;
t41 = -t62 * qJD(6) + (t102 * t249 + t106 * t161) * qJD(5);
t117 = qJD(6) * t155 - t102 * t41 - t106 * t40;
t110 = t107 * t129 - t37 * t216 + t227;
t89 = -0.2e1 * t186;
t79 = -0.2e1 * t171;
t57 = -t102 * t191 + t103 * t183;
t42 = -t123 * t97 - t127;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t46 * t51 - 0.2e1 * t47 * t50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t12 * t22 + 0.2e1 * t243 * t35 + 0.2e1 * t247, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t14 * t4 + 0.2e1 * t247 + 0.2e1 * t253, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t1 * t9 + 0.2e1 * t10 * t2 + 0.2e1 * t253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-t252 * t47 - t105 * t46 + (-t105 * t50 + t252 * t51) * qJD(3)) * t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t53 - t22 * t33 + (t194 * t35 + t47 * t68) * t97 + t160, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t14 + t4 * t37 + t160 + t166, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t23 - t10 * t7 + t2 * t24 + t8 * t9 + t166; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t173 * t68 - 0.2e1 * t33 * t53 + 0.2e1 * t242, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t114 * t37 + 0.2e1 * t242 + 0.2e1 * t245, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t23 * t8 - 0.2e1 * t24 * t7 + 0.2e1 * t245; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t46, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t244 - t11 * t225 + t35 * t196, t104 * t244 - t12 * t225 + t195 * t35 (t104 * t11 + t108 * t12 + (-t104 * t22 + t108 * t21) * qJD(4)) * t97, -pkin(3) * t244 + t11 * t130 + t12 * t73 - t22 * t65 + t246, 0, 0, 0, 0, 0, 0, t11 * t69 + t21 * t54 + (t108 * t3 - t13 * t219) * t97, t11 * t70 - t14 * t196 - t21 * t257 + t228 * t4, -t13 * t257 - t14 * t54 + t3 * t70 - t4 * t69, t11 * t64 - t13 * t26 - t14 * t25 - t3 * t42 + t4 * t43 + t246, 0, 0, 0, 0, 0, 0, t1 * t69 + t13 * t29 + t3 * t55 + t54 * t9, -t10 * t54 - t13 * t28 - t2 * t69 + t3 * t56, -t1 * t56 - t10 * t29 - t2 * t55 + t28 * t9, t1 * t17 - t10 * t5 + t13 * t20 + t18 * t2 + t3 * t38 + t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t194, -t170, 0, 0, 0, 0, 0, 0, 0, 0, -t108 * t174 + t68 * t196 - t34 * t225, t104 * t174 + t195 * t68 + t225 * t33 (t104 * t34 - t108 * t33 + (-t104 * t53 + t108 * t52) * qJD(4)) * t97, -pkin(3) * t174 + t130 * t34 - t33 * t73 - t53 * t65 + t241, 0, 0, 0, 0, 0, 0, t34 * t69 + t52 * t54 + (t108 * t16 - t219 * t36) * t97, t114 * t228 - t196 * t37 - t257 * t52 + t34 * t70, -t114 * t69 + t16 * t70 - t257 * t36 - t37 * t54, t114 * t43 - t16 * t42 - t37 * t25 - t36 * t26 + t34 * t64 + t241, 0, 0, 0, 0, 0, 0, t16 * t55 + t23 * t54 + t29 * t36 + t69 * t8, t16 * t56 - t24 * t54 - t28 * t36 + t69 * t7, t23 * t28 - t24 * t29 + t55 * t7 - t56 * t8, t16 * t38 + t17 * t8 - t18 * t7 + t20 * t36 + t23 * t6 - t24 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t171, 0.2e1 * (-t104 ^ 2 + t255) * t91 * qJD(4), 0.2e1 * t178 * t228, t79, -0.2e1 * t97 * t164, 0, -0.2e1 * t91 * pkin(3) * t219 - 0.2e1 * t66 * t225, -0.2e1 * pkin(3) * t198 + 0.2e1 * t225 * t65 (t104 * t66 - t108 * t65 + (-t104 * t73 + t108 * t130) * qJD(4)) * t256, 0.2e1 * t130 * t66 - 0.2e1 * t65 * t73, -0.2e1 * t70 * t257, 0.2e1 * t257 * t69 - 0.2e1 * t70 * t54 (-(t107 * t177 - t86) * t108 + (t104 * t70 - t255 * t229) * qJD(4)) * t256, t209 (t108 * t54 - t219 * t69) * t256, t79, 0.2e1 * t64 * t54 + 0.2e1 * t66 * t69 + 0.2e1 * (-t108 * t26 + t219 * t42) * t97, -0.2e1 * t196 * t43 - 0.2e1 * t228 * t25 - 0.2e1 * t257 * t64 + 0.2e1 * t66 * t70, 0.2e1 * t25 * t69 + 0.2e1 * t257 * t42 - 0.2e1 * t26 * t70 - 0.2e1 * t43 * t54, -0.2e1 * t25 * t43 + 0.2e1 * t26 * t42 + 0.2e1 * t64 * t66, -0.2e1 * t56 * t28, 0.2e1 * t28 * t55 - 0.2e1 * t29 * t56, -0.2e1 * t28 * t69 + 0.2e1 * t54 * t56, 0.2e1 * t55 * t29, -0.2e1 * t29 * t69 - 0.2e1 * t54 * t55, t209, 0.2e1 * t17 * t54 + 0.2e1 * t20 * t55 + 0.2e1 * t29 * t38 + 0.2e1 * t6 * t69, -0.2e1 * t18 * t54 + 0.2e1 * t20 * t56 - 0.2e1 * t28 * t38 + 0.2e1 * t5 * t69, 0.2e1 * t17 * t28 - 0.2e1 * t18 * t29 + 0.2e1 * t5 * t55 - 0.2e1 * t56 * t6, 0.2e1 * t17 * t6 - 0.2e1 * t18 * t5 + 0.2e1 * t20 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t11 + t21 * t216, t103 * t11 + t21 * t214, t118, -pkin(4) * t11 + pkin(11) * t118, 0, 0, 0, 0, 0, 0 (t13 * t217 - t1) * t107 + (qJD(5) * t9 + t149) * t103 (t13 * t215 + t2) * t107 + (-qJD(5) * t10 - t148) * t103, t159 * t214 + (-t1 * t106 - t102 * t2 + (-t10 * t106 + t102 * t9) * qJD(6)) * t103, t1 * t61 - t10 * t40 + t2 * t62 + t41 * t9 + (t13 * t214 + t248) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t33, 0, 0, 0, 0, 0, 0, 0, 0, -t107 * t34 + t216 * t52, t103 * t34 + t214 * t52, t110, -t34 * pkin(4) + pkin(11) * t110, 0, 0, 0, 0, 0, 0 (t217 * t36 - t8) * t107 + (qJD(5) * t23 + t146) * t103 (t215 * t36 - t7) * t107 + (-qJD(5) * t24 - t145) * t103, t157 * t214 + (t102 * t7 - t106 * t8 + (t102 * t23 - t106 * t24) * qJD(6)) * t103, t23 * t41 - t24 * t40 + t61 * t8 - t62 * t7 + (t214 * t36 + t227) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195, 0, -t196, 0, -t66, t65, 0, 0, -t86 * t103 + (t103 * t131 + t221) * t107, -t103 * t54 - t107 * t257 - t214 * t69 - t216 * t70, t133 * t97, t147, t260 * t97, 0, -pkin(4) * t54 - t66 * t107 - t133 * t254 + t216 * t64, pkin(4) * t257 + t66 * t103 + t214 * t64 - t260 * t254, t147 * pkin(11) + t203 * t70 - t214 * t42 - t216 * t43 - t249 * t257 + t154, -pkin(4) * t66 + ((-t43 * t103 - t42 * t107) * qJD(5) + t154) * pkin(11), t56 * t191 + (-t212 * t56 - t233) * t103, t156 * t214 + (t237 - t232 + (-t231 + t236) * qJD(6)) * t103 (t215 * t69 + t28) * t107 + (qJD(5) * t56 - t141) * t103, t55 * t189 + (t103 * t29 + t214 * t55) * t102 (-t217 * t69 + t29) * t107 + (-qJD(5) * t55 - t142) * t103, t147, t41 * t69 + t54 * t61 + (-t6 + (pkin(11) * t55 + t102 * t38) * qJD(5)) * t107 + (pkin(11) * t29 + qJD(5) * t17 + t144) * t103, t40 * t69 - t54 * t62 + (-t5 + (pkin(11) * t56 + t106 * t38) * qJD(5)) * t107 + (-pkin(11) * t28 - qJD(5) * t18 - t143) * t103, t28 * t61 - t29 * t62 + t40 * t55 - t41 * t56 + t158 * t214 + (t102 * t5 - t106 * t6 + (t102 * t17 - t106 * t18) * qJD(6)) * t103, t17 * t41 - t18 * t40 - t5 * t62 + t6 * t61 + (t103 * t20 + t214 * t38) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, -0.2e1 * t179, 0, t89, 0, 0, t103 * t208, t107 * t208, 0, 0, 0.2e1 * t186 * t94 - 0.2e1 * t172, -0.4e1 * t102 * t168 + 0.2e1 * t183 * t93, 0.2e1 * t103 * t190 + 0.2e1 * t215 * t238, 0.2e1 * t186 * t92 + 0.2e1 * t172, -0.2e1 * t102 * t179 + 0.2e1 * t103 * t188, t89, 0.2e1 * t61 * t216 - 0.2e1 * t107 * t41 + 0.2e1 * (t102 * t176 + t211 * t93) * pkin(11), -0.2e1 * t62 * t216 - 0.2e1 * t107 * t40 + 0.2e1 * (-t212 * t93 + 0.2e1 * t168) * pkin(11), 0.2e1 * t155 * t214 + 0.2e1 * (t102 * t40 - t106 * t41 + (t102 * t61 - t106 * t62) * qJD(6)) * t103, 0.2e1 * pkin(11) ^ 2 * t186 - 0.2e1 * t40 * t62 + 0.2e1 * t41 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, 0, 0, 0, t148, t149, t121, -pkin(5) * t3 + pkin(12) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t114, 0, 0, 0, 0, 0, 0, 0, 0, t145, t146, t119, -pkin(5) * t16 + pkin(12) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257, 0, -t54, t196, t26, t25, 0, 0, t211 * t56 - t237, qJD(6) * t156 - t102 * t29 - t233, t142, t212 * t55 - t232, -t141, 0, -pkin(5) * t29 - pkin(12) * t142 + t143, pkin(5) * t28 + pkin(12) * t141 + t144 (-t237 - t232 + (t231 + t236) * qJD(6)) * pkin(12) + t120, -pkin(5) * t20 + pkin(12) * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, 0, -t216, 0, -t203, pkin(11) * t216, 0, 0, -t57, -0.4e1 * t103 * t187 - t214 * t239, t102 * t216 - t188, t57, t132, 0 (pkin(12) * t223 + (-pkin(5) * t106 + t250) * t103) * qJD(6) + (t102 * t162 - t204) * qJD(5) (t106 * t249 + t140) * qJD(6) + (t106 * t162 + t205) * qJD(5), t117, -pkin(5) * t203 + pkin(12) * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t187, -0.2e1 * t183, 0, -0.2e1 * t187, 0, 0, t102 * t207, t106 * t207, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, 0, -t29, t54, t6, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103 * t212 + t191, 0, -t102 * t214 - t189, t216, t41, t40, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, -t212, 0, -pkin(12) * t211, pkin(12) * t212, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t15;
