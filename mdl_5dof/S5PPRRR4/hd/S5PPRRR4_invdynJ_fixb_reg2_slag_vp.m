% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:50
% EndTime: 2019-12-05 15:19:57
% DurationCPUTime: 3.50s
% Computational Cost: add. (3516->361), mult. (9311->537), div. (0->0), fcn. (8672->14), ass. (0->174)
t122 = sin(qJ(4));
t240 = pkin(8) * t122;
t213 = cos(pkin(5));
t102 = t213 * qJDD(1) + qJDD(2);
t105 = t213 * qJD(1) + qJD(2);
t123 = sin(qJ(3));
t126 = cos(qJ(3));
t119 = sin(pkin(5));
t120 = cos(pkin(11));
t212 = cos(pkin(6));
t173 = t120 * t212;
t167 = t119 * t173;
t154 = qJDD(1) * t167;
t117 = sin(pkin(11));
t207 = t117 * t119;
t184 = qJD(1) * t207;
t171 = qJD(3) * t184;
t118 = sin(pkin(6));
t198 = qJD(3) * t123;
t183 = t118 * t198;
t204 = t118 * t126;
t155 = qJD(1) * t167;
t238 = qJD(3) * t155 + qJDD(1) * t207;
t157 = -t102 * t204 + t105 * t183 + (-t154 + t171) * t126 + t238 * t123;
t125 = cos(qJ(4));
t163 = pkin(4) * t122 - pkin(9) * t125;
t95 = t163 * qJD(4);
t230 = -t125 * pkin(4) - t122 * pkin(9);
t99 = -pkin(3) + t230;
t12 = qJD(3) * t95 + t99 * qJDD(3) + t157;
t121 = sin(qJ(5));
t124 = cos(qJ(5));
t138 = (t117 * t126 + t123 * t173) * t119;
t205 = t118 * t123;
t50 = qJD(1) * t138 + t105 * t205;
t48 = qJD(3) * pkin(8) + t50;
t186 = t119 * t120 * t118;
t70 = -qJD(1) * t186 + t212 * t105;
t35 = t122 * t70 + t125 * t48;
t33 = qJD(4) * pkin(9) + t35;
t98 = t123 * t184;
t49 = (t105 * t118 + t155) * t126 - t98;
t44 = t99 * qJD(3) - t49;
t158 = t121 * t33 - t124 * t44;
t196 = qJD(4) * t125;
t208 = qJDD(3) * pkin(8);
t182 = qJD(3) * t204;
t169 = -t102 * t205 - t105 * t182 - t123 * t154 - t126 * t238;
t28 = -t123 * t171 - t169;
t26 = t28 + t208;
t69 = -qJDD(1) * t186 + t212 * t102;
t187 = -t122 * t69 - t125 * t26 - t70 * t196;
t197 = qJD(4) * t122;
t7 = -t48 * t197 - t187;
t5 = qJDD(4) * pkin(9) + t7;
t1 = -t158 * qJD(5) + t121 * t12 + t124 * t5;
t190 = t125 * qJD(3);
t106 = -qJD(5) + t190;
t239 = -t158 * t106 + t1;
t199 = qJD(3) * t122;
t237 = qJD(5) * t199 - qJDD(4);
t176 = t118 * t213;
t236 = t119 * (-t117 * t123 + t126 * t173) + t126 * t176;
t188 = t122 * qJDD(3);
t52 = t121 * ((qJD(5) + t190) * qJD(4) + t188) + t237 * t124;
t235 = -pkin(8) * qJD(5) * t125 - t50 + t95;
t10 = t121 * t44 + t124 * t33;
t2 = -qJD(5) * t10 + t124 * t12 - t121 * t5;
t234 = t10 * t106 - t2;
t211 = cos(pkin(10));
t160 = t213 * t211;
t210 = sin(pkin(10));
t135 = t210 * t117 - t120 * t160;
t175 = t119 * t211;
t233 = t118 * t175 + t135 * t212;
t159 = t213 * t210;
t136 = t211 * t117 + t120 * t159;
t174 = t119 * t210;
t232 = -t118 * t174 + t136 * t212;
t66 = t125 * t69;
t8 = -t35 * qJD(4) - t122 * t26 + t66;
t191 = t124 * qJD(4);
t87 = t121 * t199 - t191;
t192 = t121 * qJD(4);
t89 = t124 * t199 + t192;
t226 = t89 * t87;
t194 = qJD(5) * t124;
t202 = t124 * t125;
t224 = t121 * t235 - t191 * t240 + t99 * t194 - t49 * t202;
t195 = qJD(5) * t121;
t203 = t121 * t125;
t223 = t124 * t235 + t192 * t240 - t99 * t195 + t49 * t203;
t222 = qJD(3) * pkin(3);
t220 = t106 * t87;
t219 = t121 * t87;
t218 = t122 * t48;
t217 = t122 * t87;
t189 = qJD(3) * qJD(4);
t181 = t125 * t189;
t51 = -qJD(5) * t191 + (-t181 - t188) * t124 + t237 * t121;
t216 = t51 * t121;
t215 = t52 * t124;
t214 = t89 * t106;
t209 = qJDD(3) * pkin(3);
t115 = t122 ^ 2;
t116 = t125 ^ 2;
t201 = t115 - t116;
t112 = t125 * qJDD(3);
t128 = qJD(3) ^ 2;
t185 = t122 * t128 * t125;
t178 = t122 * t189;
t170 = t125 * t178;
t162 = -t10 * t121 + t124 * t158;
t140 = t213 * t212 - t186;
t57 = t123 * t176 + t138;
t43 = t140 * t122 + t57 * t125;
t20 = -t121 * t236 + t43 * t124;
t19 = -t43 * t121 - t124 * t236;
t34 = t125 * t70 - t218;
t156 = qJDD(3) * t126 - t123 * t128;
t77 = t122 * t212 + t125 * t205;
t62 = -t121 * t77 - t124 * t204;
t152 = t121 * t204 - t124 * t77;
t86 = qJDD(5) - t112 + t178;
t151 = -t106 * t194 + t121 * t86;
t150 = t106 * t195 + t124 * t86;
t74 = t117 * t160 + t210 * t120;
t39 = -t123 * t233 + t74 * t126;
t75 = -t117 * t159 + t211 * t120;
t41 = -t123 * t232 + t75 * t126;
t42 = t57 * t122 - t140 * t125;
t58 = t135 * t118 - t212 * t175;
t59 = t136 * t118 + t212 * t174;
t149 = -g(1) * (-t41 * t122 + t59 * t125) - g(2) * (-t39 * t122 + t58 * t125) + g(3) * t42;
t16 = t58 * t122 + t39 * t125;
t18 = t59 * t122 + t41 * t125;
t148 = -g(1) * t18 - g(2) * t16 - g(3) * t43;
t38 = t74 * t123 + t126 * t233;
t40 = t75 * t123 + t126 * t232;
t147 = g(1) * t40 + g(2) * t38 - g(3) * t236;
t146 = g(1) * t41 + g(2) * t39 + g(3) * t57;
t6 = -qJDD(4) * pkin(4) - t8;
t144 = t149 - t6;
t76 = t122 * t205 - t125 * t212;
t32 = -qJD(4) * pkin(4) - t34;
t143 = -pkin(9) * t86 - t106 * t32;
t47 = -t49 - t222;
t141 = -pkin(8) * qJDD(4) + (t47 + t49 - t222) * qJD(4);
t139 = t50 * qJD(3) + t147;
t137 = -g(1) * t174 + g(2) * t175 - g(3) * t213;
t134 = pkin(9) * qJD(5) * t106 + t144;
t127 = qJD(4) ^ 2;
t27 = t157 - t209;
t130 = -pkin(8) * t127 + t139 + t209 - t27;
t129 = -t8 * t122 + t7 * t125 + (-t122 * t35 - t125 * t34) * qJD(4) - t146;
t94 = t163 * qJD(3);
t72 = pkin(8) * t202 + t121 * t99;
t71 = -pkin(8) * t203 + t124 * t99;
t61 = t77 * qJD(4) + t122 * t182;
t60 = -t76 * qJD(4) + t125 * t182;
t55 = t236 * pkin(3);
t54 = t57 * qJD(3);
t53 = t236 * qJD(3);
t37 = t40 * pkin(3);
t36 = t38 * pkin(3);
t31 = t152 * qJD(5) - t121 * t60 + t124 * t183;
t30 = t62 * qJD(5) + t121 * t183 + t124 * t60;
t22 = t121 * t94 + t124 * t34;
t21 = -t121 * t34 + t124 * t94;
t14 = -t42 * qJD(4) + t53 * t125;
t13 = t43 * qJD(4) + t53 * t122;
t4 = t19 * qJD(5) + t54 * t121 + t14 * t124;
t3 = -t20 * qJD(5) - t14 * t121 + t54 * t124;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t213 - g(3) + (t117 ^ 2 + t120 ^ 2) * t119 ^ 2 * qJDD(1), 0, 0, 0, 0, 0, 0, -t54 * qJD(3) + qJDD(3) * t236, -t53 * qJD(3) - t57 * qJDD(3), 0, t69 * t140 - t157 * t236 + t28 * t57 - t49 * t54 + t50 * t53 - g(3), 0, 0, 0, 0, 0, 0, t236 * t112 - t13 * qJD(4) - t42 * qJDD(4) + (-t125 * t54 - t197 * t236) * qJD(3), -t236 * t188 - t14 * qJD(4) - t43 * qJDD(4) + (t122 * t54 - t196 * t236) * qJD(3), (t122 * t42 + t125 * t43) * qJDD(3) + (t122 * t13 + t125 * t14 + (-t122 * t43 + t125 * t42) * qJD(4)) * qJD(3), -t34 * t13 + t35 * t14 - t236 * t27 - t8 * t42 + t7 * t43 + t47 * t54 - g(3), 0, 0, 0, 0, 0, 0, -t3 * t106 + t13 * t87 + t19 * t86 + t42 * t52, t4 * t106 + t13 * t89 - t20 * t86 - t42 * t51, t19 * t51 - t20 * t52 - t3 * t89 - t4 * t87, t1 * t20 + t10 * t4 + t32 * t13 - t158 * t3 + t2 * t19 + t6 * t42 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137 + t102, 0, 0, 0, 0, 0, 0, t156 * t118, (-qJDD(3) * t123 - t126 * t128) * t118, 0, t69 * t212 + (t123 * t28 - t126 * t157 + (-t123 * t49 + t126 * t50) * qJD(3)) * t118 + t137, 0, 0, 0, 0, 0, 0, -t61 * qJD(4) - t76 * qJDD(4) + (t125 * t156 - t126 * t178) * t118, -t60 * qJD(4) - t77 * qJDD(4) + (-t122 * t156 - t126 * t181) * t118, (t122 * t76 + t125 * t77) * qJDD(3) + (t122 * t61 + t125 * t60 + (-t122 * t77 + t125 * t76) * qJD(4)) * qJD(3), -t34 * t61 + t35 * t60 + t7 * t77 - t8 * t76 + (-t126 * t27 + t198 * t47) * t118 + t137, 0, 0, 0, 0, 0, 0, -t31 * t106 + t76 * t52 + t61 * t87 + t62 * t86, t30 * t106 + t152 * t86 - t76 * t51 + t61 * t89, t152 * t52 - t30 * t87 - t31 * t89 + t62 * t51, -t1 * t152 + t10 * t30 - t158 * t31 + t2 * t62 + t32 * t61 + t6 * t76 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t139 - t157, (t49 + t98) * qJD(3) + t146 + t169, 0, 0, t115 * qJDD(3) + 0.2e1 * t170, 0.2e1 * t122 * t112 - 0.2e1 * t201 * t189, qJDD(4) * t122 + t127 * t125, t116 * qJDD(3) - 0.2e1 * t170, qJDD(4) * t125 - t127 * t122, 0, t122 * t141 + t125 * t130, -t122 * t130 + t125 * t141, t129 + (-t49 * qJD(3) + t208) * (t115 + t116), -t27 * pkin(3) + g(1) * t37 + g(2) * t36 - g(3) * t55 - t47 * t50 + (t34 * t122 - t35 * t125) * t49 + t129 * pkin(8), t89 * t125 * t191 + (-t51 * t124 - t195 * t89) * t122, (-t121 * t89 - t124 * t87) * t196 + (t216 - t215 + (-t124 * t89 + t219) * qJD(5)) * t122, (-t106 * t191 + t51) * t125 + (qJD(4) * t89 + t150) * t122, t194 * t217 + (t122 * t52 + t196 * t87) * t121, (t106 * t192 + t52) * t125 + (-qJD(4) * t87 - t151) * t122, -t106 * t197 - t86 * t125, t71 * t86 - t223 * t106 - t146 * t121 + (-t2 + (pkin(8) * t87 + t121 * t32) * qJD(4) + t147 * t124) * t125 + (pkin(8) * t52 - qJD(4) * t158 + t6 * t121 + t194 * t32 - t49 * t87) * t122, -t72 * t86 + t224 * t106 - t146 * t124 + (t1 + (pkin(8) * t89 + t124 * t32) * qJD(4) - t147 * t121) * t125 + (-pkin(8) * t51 - t10 * qJD(4) + t6 * t124 - t195 * t32 - t49 * t89) * t122, t71 * t51 - t72 * t52 - t223 * t89 - t224 * t87 + t162 * t196 + (-t1 * t121 - t124 * t2 + (-t10 * t124 - t121 * t158) * qJD(5) + t147) * t122, t1 * t72 + t2 * t71 - t32 * t122 * t49 - g(1) * (t230 * t40 - t37) - g(2) * (t230 * t38 - t36) - g(3) * (-t230 * t236 + t55) - t223 * t158 + t224 * t10 + (t122 * t6 + t196 * t32 - t146) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, t201 * t128, t188, t185, t112, qJDD(4), t66 + (-qJD(3) * t47 - t26) * t122 + t149, -t47 * t190 + (t34 + t218) * qJD(4) - t148 + t187, 0, 0, -t124 * t214 - t216, (-t51 + t220) * t124 + (-t52 + t214) * t121, (t106 * t202 - t122 * t89) * qJD(3) + t151, -t106 * t219 - t215, (-t106 * t203 + t217) * qJD(3) + t150, t106 * t199, -pkin(4) * t52 + t21 * t106 + t121 * t143 + t124 * t134 + t158 * t199 - t35 * t87, pkin(4) * t51 + t10 * t199 - t22 * t106 - t121 * t134 + t124 * t143 - t35 * t89, t21 * t89 + t22 * t87 + ((qJD(5) * t89 - t52) * pkin(9) + t239) * t124 + ((qJD(5) * t87 - t51) * pkin(9) + t234) * t121 + t148, -t10 * t22 + t158 * t21 - t32 * t35 + t144 * pkin(4) + (qJD(5) * t162 + t1 * t124 - t2 * t121 + t148) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t226, -t87 ^ 2 + t89 ^ 2, -t51 - t220, -t226, -t214 - t52, t86, -t32 * t89 - g(1) * (-t121 * t18 + t124 * t40) - g(2) * (-t121 * t16 + t124 * t38) - g(3) * t19 - t234, t32 * t87 - g(1) * (-t121 * t40 - t124 * t18) - g(2) * (-t121 * t38 - t124 * t16) + g(3) * t20 - t239, 0, 0;];
tau_reg = t9;
