% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR13_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:58
% EndTime: 2019-12-31 18:33:04
% DurationCPUTime: 2.94s
% Computational Cost: add. (3630->361), mult. (8503->439), div. (0->0), fcn. (6233->10), ass. (0->188)
t123 = sin(pkin(8));
t124 = cos(pkin(8));
t127 = sin(qJ(3));
t230 = cos(qJ(3));
t88 = t123 * t230 + t124 * t127;
t241 = t88 * qJD(1);
t247 = qJD(5) + t241;
t129 = cos(qJ(5));
t126 = sin(qJ(5));
t192 = t126 * qJD(3);
t183 = t230 * t124;
t168 = qJD(1) * t183;
t206 = t127 * t123;
t182 = qJD(1) * t206;
t79 = -t168 + t182;
t54 = -t129 * t79 + t192;
t178 = t247 * t54;
t193 = qJD(5) * t129;
t179 = qJDD(1) * t230;
t189 = t123 * qJDD(1);
t157 = -t124 * t179 + t127 * t189;
t84 = t88 * qJD(3);
t49 = qJD(1) * t84 + t157;
t21 = qJD(5) * t192 - qJDD(3) * t129 - t126 * t49 - t193 * t79;
t249 = t21 - t178;
t114 = t124 * pkin(2) + pkin(1);
t92 = -qJD(1) * t114 + qJD(2);
t143 = -qJ(4) * t241 + t92;
t233 = pkin(3) + pkin(7);
t24 = t233 * t79 + t143;
t219 = pkin(6) + qJ(2);
t93 = t219 * t123;
t89 = qJD(1) * t93;
t94 = t219 * t124;
t90 = qJD(1) * t94;
t50 = t127 * t90 + t230 * t89;
t199 = -qJD(4) - t50;
t198 = pkin(4) * t241 - t199;
t27 = -qJD(3) * t233 + t198;
t6 = -t126 * t24 + t129 * t27;
t251 = t6 * t247;
t188 = t124 * qJDD(1);
t187 = qJD(3) * t168 + t123 * t179 + t127 * t188;
t48 = qJD(3) * t182 - t187;
t215 = t48 * qJ(4);
t91 = -qJDD(1) * t114 + qJDD(2);
t136 = -qJD(4) * t241 + t215 + t91;
t5 = t233 * t49 + t136;
t7 = t126 * t27 + t129 * t24;
t181 = qJD(3) * t230;
t194 = qJD(3) * t127;
t191 = qJD(1) * qJD(2);
t235 = qJDD(1) * t219 + t191;
t62 = t235 * t123;
t63 = t235 * t124;
t176 = t127 * t63 + t181 * t90 - t194 * t89 + t230 * t62;
t158 = qJDD(4) + t176;
t9 = -t48 * pkin(4) - qJDD(3) * t233 + t158;
t2 = -qJD(5) * t7 - t126 * t5 + t129 * t9;
t238 = t247 * t7 + t2;
t242 = t129 * t247;
t56 = qJD(3) * t129 + t126 * t79;
t250 = t56 * t242;
t175 = t126 * t247;
t46 = -qJDD(5) + t48;
t149 = -t129 * t46 - t175 * t247;
t212 = qJDD(1) * pkin(1);
t128 = sin(qJ(1));
t130 = cos(qJ(1));
t240 = g(1) * t128 - g(2) * t130;
t150 = -qJDD(2) + t212 + t240;
t122 = pkin(8) + qJ(3);
t116 = sin(t122);
t117 = cos(t122);
t165 = g(1) * t130 + g(2) * t128;
t142 = -g(3) * t116 - t117 * t165;
t177 = t127 * t62 + t181 * t89 + t194 * t90 - t230 * t63;
t248 = t142 - t177;
t246 = 0.2e1 * qJD(3) * t241 + t157;
t234 = t79 ^ 2;
t77 = t241 ^ 2;
t244 = -t234 - t77;
t243 = -t234 + t77;
t196 = pkin(3) * t117 + t116 * qJ(4);
t239 = qJ(2) * qJDD(1);
t231 = t79 * pkin(4);
t51 = -t127 * t89 + t230 * t90;
t44 = -qJD(3) * qJ(4) - t51;
t28 = -t44 - t231;
t237 = t233 * t46 + t247 * t28;
t36 = (qJD(2) * t123 + qJD(3) * t94) * t127 - qJD(2) * t183 + t93 * t181;
t53 = -t127 * t93 + t230 * t94;
t236 = t36 * qJD(3) - t53 * qJDD(3) - t116 * t240;
t232 = t49 * pkin(3);
t112 = g(3) * t117;
t225 = t117 * pkin(7);
t224 = t56 * t54;
t223 = t56 * t79;
t222 = t79 * t54;
t221 = t79 * t241;
t220 = pkin(4) + t219;
t218 = t129 * t21;
t171 = t126 * qJDD(3) - t129 * t49;
t22 = qJD(5) * t56 + t171;
t216 = t22 * t126;
t214 = t79 * qJ(4);
t211 = qJDD(3) * pkin(3);
t210 = t116 * t128;
t209 = t116 * t130;
t208 = t117 * t128;
t207 = t117 * t130;
t205 = t128 * t126;
t204 = t128 * t129;
t203 = t130 * t219;
t202 = t130 * t126;
t201 = t130 * t129;
t200 = t51 * qJD(3);
t120 = t123 ^ 2;
t121 = t124 ^ 2;
t195 = t120 + t121;
t190 = qJDD(3) * qJ(4);
t87 = -t183 + t206;
t186 = qJD(5) * t126 * t87;
t185 = t87 * t193;
t184 = -g(1) * t209 - g(2) * t210 + t112;
t52 = t127 * t94 + t230 * t93;
t172 = t195 * qJD(1) ^ 2;
t98 = t130 * t114;
t170 = g(2) * (pkin(3) * t207 + qJ(4) * t209 + t98);
t169 = 0.2e1 * t195;
t167 = g(1) * t208 - g(2) * t207;
t1 = qJD(5) * t6 + t126 * t9 + t129 * t5;
t166 = t1 - t251;
t16 = -qJD(3) * qJD(4) + t177 - t190;
t10 = -pkin(4) * t49 - t16;
t163 = t10 * t87 + t28 * t84;
t162 = -t126 * t6 + t129 * t7;
t161 = t247 * t84 - t46 * t87;
t83 = t123 * t194 - t124 * t181;
t160 = -t241 * t83 - t48 * t88;
t159 = t49 * t87 + t79 * t84;
t151 = -t88 * qJ(4) - t114;
t32 = t233 * t87 + t151;
t38 = pkin(4) * t88 + t52;
t15 = t126 * t38 + t129 * t32;
t14 = -t126 * t32 + t129 * t38;
t156 = t83 * qJ(4) - t88 * qJD(4);
t153 = qJD(3) * t83 - qJDD(3) * t88;
t152 = qJD(3) * t84 + qJDD(3) * t87;
t148 = -t114 - t196;
t147 = -t176 - t184;
t146 = t126 * t46 - t242 * t247;
t144 = t150 + t212;
t141 = -t241 * t84 + t48 * t87 - t49 * t88 + t79 * t83;
t37 = qJD(2) * t88 + qJD(3) * t53;
t140 = -qJD(3) * t37 - qJDD(3) * t52 + t167;
t139 = -t240 + t91;
t35 = t79 * pkin(3) + t143;
t138 = t241 * t35 + qJDD(4) - t147;
t137 = t169 * t191 - t165;
t135 = qJD(5) * t233 * t247 + t10 + t142;
t134 = t241 * t37 + t36 * t79 - t48 * t52 - t49 * t53 - t165;
t97 = qJ(4) * t207;
t95 = qJ(4) * t208;
t74 = -t116 * t205 + t201;
t73 = t116 * t204 + t202;
t72 = t116 * t202 + t204;
t71 = t116 * t201 - t205;
t65 = qJD(3) * t79;
t47 = pkin(3) * t87 + t151;
t45 = pkin(3) * t241 + t214;
t43 = -qJD(3) * pkin(3) - t199;
t39 = -pkin(4) * t87 + t53;
t34 = t51 - t231;
t31 = pkin(3) * t84 + t156;
t30 = t48 - t65;
t29 = t233 * t241 + t214;
t26 = -t83 * pkin(4) + t37;
t25 = -t84 * pkin(4) - t36;
t23 = t233 * t84 + t156;
t20 = t129 * t22;
t17 = t158 - t211;
t13 = t136 + t232;
t12 = t126 * t34 + t129 * t29;
t11 = -t126 * t29 + t129 * t34;
t4 = -qJD(5) * t15 - t126 * t23 + t129 * t26;
t3 = qJD(5) * t14 + t126 * t26 + t129 * t23;
t8 = [0, 0, 0, 0, 0, qJDD(1), t240, t165, 0, 0, t120 * qJDD(1), 0.2e1 * t123 * t188, 0, t121 * qJDD(1), 0, 0, t144 * t124, -t144 * t123, t169 * t239 + t137, t150 * pkin(1) + (t195 * t239 + t137) * qJ(2), t160, t141, -t153, t159, -t152, 0, -t114 * t49 + t84 * t92 + t87 * t91 + t140, t114 * t48 - t92 * t83 + t91 * t88 + t236, t176 * t88 + t177 * t87 - t50 * t83 - t51 * t84 + t134, -t177 * t53 - t51 * t36 + t176 * t52 + t50 * t37 - t91 * t114 - g(1) * (-t128 * t114 + t203) - g(2) * (t128 * t219 + t98), 0, t153, t152, t160, t141, t159, t16 * t87 + t17 * t88 - t43 * t83 + t44 * t84 + t134, -t13 * t87 - t31 * t79 - t35 * t84 - t47 * t49 - t140, -t13 * t88 - t241 * t31 + t35 * t83 + t47 * t48 - t236, t13 * t47 + t35 * t31 - t16 * t53 + t44 * t36 + t17 * t52 + t43 * t37 - g(1) * t203 - t170 + (-g(1) * t148 - g(2) * t219) * t128, t56 * t185 + (-t21 * t87 + t56 * t84) * t126, (-t126 * t54 + t129 * t56) * t84 + (-t216 - t218 + (-t126 * t56 - t129 * t54) * qJD(5)) * t87, t126 * t161 + t185 * t247 - t21 * t88 - t56 * t83, t54 * t186 + (-t22 * t87 - t54 * t84) * t129, t129 * t161 - t186 * t247 - t22 * t88 + t54 * t83, -t247 * t83 - t46 * t88, -g(1) * t74 - g(2) * t72 - t129 * t163 - t14 * t46 + t186 * t28 + t2 * t88 + t39 * t22 + t247 * t4 + t25 * t54 - t6 * t83, g(1) * t73 - g(2) * t71 - t1 * t88 + t126 * t163 + t15 * t46 + t185 * t28 - t39 * t21 - t247 * t3 + t25 * t56 + t7 * t83, t14 * t21 - t15 * t22 - t3 * t54 - t4 * t56 + t162 * t84 + (t1 * t129 - t2 * t126 + (-t126 * t7 - t129 * t6) * qJD(5)) * t87 + t167, t1 * t15 + t7 * t3 + t2 * t14 + t6 * t4 + t10 * t39 + t28 * t25 - t170 + (-g(1) * t220 - g(2) * t225) * t130 + (-g(1) * (t148 - t225) - g(2) * t220) * t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, t189, -t172, -qJ(2) * t172 - t150, 0, 0, 0, 0, 0, 0, t246, (-t79 - t182) * qJD(3) + t187, t244, -t241 * t50 + t51 * t79 + t139, 0, 0, 0, 0, 0, 0, t244, -t246, t48 + t65, t232 + t215 - t44 * t79 + (-qJD(4) - t43) * t241 + t139, 0, 0, 0, 0, 0, 0, t146 + t222, t223 - t149, -t126 * t249 - t20 + t250, -t126 * t238 + t129 * t166 + t28 * t79 - t240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, t243, (t79 - t182) * qJD(3) + t187, -t221, -t157, qJDD(3), -t241 * t92 + t147 + t200, -t50 * qJD(3) + t92 * t79 - t248, 0, 0, qJDD(3), t30, t157, t221, t243, -t221, pkin(3) * t48 - qJ(4) * t49 + (-t44 - t51) * t241 + (t43 + t199) * t79, t45 * t79 + t138 - t200 - 0.2e1 * t211, 0.2e1 * t190 - t35 * t79 + t45 * t241 + (0.2e1 * qJD(4) + t50) * qJD(3) + t248, -t16 * qJ(4) - t17 * pkin(3) - t35 * t45 - t43 * t51 - g(1) * (-pkin(3) * t209 + t97) - g(2) * (-pkin(3) * t210 + t95) - g(3) * t196 + t199 * t44, -t175 * t56 - t218, -t20 - t250 + (t21 + t178) * t126, t149 + t223, t242 * t54 + t216, t146 - t222, t247 * t79, qJ(4) * t22 - t11 * t247 + t126 * t135 + t129 * t237 + t198 * t54 + t6 * t79, -qJ(4) * t21 + t12 * t247 - t126 * t237 + t129 * t135 + t198 * t56 - t7 * t79, t11 * t56 + t12 * t54 + (-t233 * t21 - t7 * t241 - t2 + (t233 * t54 - t7) * qJD(5)) * t129 + (t233 * t22 + t6 * t241 - t1 + (-t233 * t56 + t6) * qJD(5)) * t126 - t184, t10 * qJ(4) - t7 * t12 - t6 * t11 - g(1) * t97 - g(2) * t95 - g(3) * (t196 + t225) + t198 * t28 + (-qJD(5) * t162 - t1 * t126 + t116 * t165 - t2 * t129) * t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, qJDD(3) - t221, -qJD(3) ^ 2 - t77, qJD(3) * t44 + t138 - t211, 0, 0, 0, 0, 0, 0, -qJD(3) * t54 + t149, -qJD(3) * t56 + t146, t249 * t129 + (t247 * t56 - t22) * t126, -t28 * qJD(3) + t126 * t166 + t129 * t238 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, -t54 ^ 2 + t56 ^ 2, -t249, -t224, -t171 + (-qJD(5) + t247) * t56, -t46, -g(1) * t71 - g(2) * t73 + t112 * t129 - t28 * t56 + t238, g(1) * t72 - g(2) * t74 + t28 * t54 + t251 + (-qJD(5) * t27 - t5) * t129 + (qJD(5) * t24 - t112 - t9) * t126, 0, 0;];
tau_reg = t8;
