% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRPP2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:52:04
% EndTime: 2019-12-31 20:52:08
% DurationCPUTime: 2.07s
% Computational Cost: add. (1848->349), mult. (2689->364), div. (0->0), fcn. (1332->8), ass. (0->196)
t127 = sin(qJ(3));
t124 = t127 ^ 2;
t130 = cos(qJ(3));
t125 = t130 ^ 2;
t201 = t124 + t125;
t119 = qJDD(1) + qJDD(2);
t128 = sin(qJ(2));
t193 = qJDD(1) * t128;
t131 = cos(qJ(2));
t198 = qJD(2) * t131;
t37 = t119 * pkin(7) + (qJD(1) * t198 + t193) * pkin(1);
t250 = t201 * t37;
t219 = qJ(4) * t130;
t236 = pkin(3) + pkin(4);
t239 = t236 * t127;
t249 = t219 - t239;
t226 = pkin(1) * qJD(1);
t241 = t201 * t131;
t248 = t241 * t226;
t195 = qJD(3) * t130;
t247 = qJ(4) * t195 + t127 * qJD(4);
t181 = t236 * qJD(3);
t120 = qJD(1) + qJD(2);
t208 = t120 * t127;
t189 = t128 * t226;
t63 = pkin(7) * t120 + t189;
t49 = t127 * t63;
t26 = qJ(5) * t208 - t49;
t205 = qJD(4) - t26;
t15 = -t181 + t205;
t209 = t119 * t125;
t210 = t119 * t124;
t246 = t209 + t210;
t176 = t236 * qJDD(3);
t200 = qJD(1) * t131;
t186 = pkin(1) * t200;
t231 = t120 * pkin(2);
t64 = -t186 - t231;
t245 = t128 * t64 + t241 * t63;
t31 = t127 * t37;
t46 = t63 * t195;
t183 = qJDD(4) + t31 + t46;
t216 = qJDD(3) * pkin(3);
t10 = t183 - t216;
t121 = qJDD(3) * qJ(4);
t122 = qJD(3) * qJD(4);
t196 = qJD(3) * t127;
t33 = t130 * t37;
t9 = -t63 * t196 + t121 + t122 + t33;
t244 = t10 * t127 + t9 * t130;
t126 = qJ(1) + qJ(2);
t109 = sin(t126);
t110 = cos(t126);
t227 = g(1) * t110 + g(2) * t109;
t123 = qJD(3) * qJ(4);
t207 = t120 * t130;
t50 = t130 * t63;
t27 = -qJ(5) * t207 + t50;
t23 = t123 + t27;
t232 = g(2) * t110;
t98 = g(1) * t109;
t242 = t98 - t232;
t111 = t127 * qJ(4);
t174 = pkin(2) + t111;
t12 = t186 + qJD(5) + (t236 * t130 + t174) * t120;
t204 = qJD(5) + t12;
t240 = t204 * t127;
t166 = -qJD(3) * pkin(3) + qJD(4);
t35 = t166 + t49;
t41 = t50 + t123;
t238 = t35 * t195 - t41 * t196 + t244;
t197 = qJD(3) * t120;
t202 = -t124 + t125;
t92 = t130 * t119;
t237 = 0.2e1 * t127 * t92 + 0.2e1 * t202 * t197;
t235 = pkin(1) * t131;
t134 = qJD(3) ^ 2;
t234 = pkin(7) * t134;
t129 = sin(qJ(1));
t233 = g(1) * t129;
t106 = t119 * pkin(2);
t115 = t130 * pkin(3);
t114 = t130 * pkin(4);
t230 = pkin(7) - qJ(5);
t212 = t110 * t127;
t214 = t109 * t127;
t229 = -g(1) * t214 + g(2) * t212;
t228 = t110 * pkin(2) + t109 * pkin(7);
t203 = t115 + t111;
t182 = t114 + t203;
t52 = pkin(2) + t182;
t225 = t120 * t52;
t223 = t130 * t23;
t222 = t35 * t130;
t199 = qJD(2) * t128;
t188 = pkin(1) * t199;
t221 = qJD(1) * t188 - qJDD(1) * t235;
t220 = pkin(7) * qJDD(3);
t218 = qJ(5) * t119;
t102 = pkin(1) * t128 + pkin(7);
t206 = -qJ(5) + t102;
t55 = t206 * t130;
t217 = qJD(3) * t55;
t215 = t102 * t134;
t213 = t109 * t130;
t211 = t110 * t130;
t91 = t127 * t119;
t194 = t127 * qJD(5);
t192 = qJDD(3) * t102;
t191 = t23 * t196 + t227;
t165 = t120 * t189;
t83 = g(1) * t213;
t190 = t130 * t165 + t186 * t196 + t83;
t187 = pkin(1) * t198;
t185 = 0.2e1 * t121 + 0.2e1 * t122 + t33;
t36 = -t106 + t221;
t44 = pkin(3) * t196 - t247;
t103 = -pkin(2) - t235;
t180 = t120 * t199;
t179 = t120 * t196;
t178 = t120 * t195;
t100 = qJ(5) * t196;
t177 = -t36 - t232;
t71 = t230 * t130;
t146 = pkin(3) * t92 + qJ(4) * t91 + t247 * t120 - t36;
t143 = pkin(4) * t92 + qJDD(5) + t146;
t3 = -t236 * t179 + t143;
t175 = t12 * t195 + t3 * t127 - t229;
t172 = t36 * t127 + t64 * t195 + t229;
t171 = pkin(3) * t211 + t110 * t111 + t228;
t95 = t110 * pkin(7);
t170 = -t110 * qJ(5) + t95;
t167 = g(1) * t212 + g(2) * t214 - g(3) * t130 - t31;
t163 = t127 * t178;
t162 = g(1) * (-t109 * pkin(2) + t95);
t132 = cos(qJ(1));
t116 = t132 * pkin(1);
t161 = t116 + t171;
t160 = -qJD(5) + t187;
t159 = -t106 + t234;
t157 = -g(2) * t132 + t233;
t156 = -qJDD(4) + t167;
t155 = pkin(3) * t127 - t219;
t28 = -pkin(4) * t196 - t44;
t20 = t28 - t188;
t51 = t103 - t203;
t38 = t114 - t51;
t154 = t119 * t38 + t120 * t20;
t153 = t119 * t52 + t120 * t28;
t152 = t127 * t35 + t130 * t41;
t151 = -t221 + t242;
t150 = t201 * t120 * t187 + t246 * t102 - t227;
t149 = -t174 - t115;
t148 = -t46 + t156;
t147 = t120 * t51 - t187;
t6 = pkin(3) * t179 - t146;
t65 = -pkin(2) - t203;
t145 = -t119 * t65 - t120 * t44 - t234 - t6;
t29 = t44 + t188;
t144 = -t119 * t51 - t120 * t29 - t215 - t6;
t142 = -qJ(5) * t91 - t148;
t141 = pkin(1) * t180 + t103 * t119 + t215;
t140 = t246 * pkin(7) - t120 * t248 - t227;
t139 = -g(1) * t95 - t149 * t98;
t138 = -t192 + (t103 * t120 - t187) * qJD(3);
t137 = (-t41 * t127 + t222) * qJD(3) + t244;
t136 = (-g(1) * (t149 - t114) + g(2) * qJ(5)) * t109;
t118 = t120 ^ 2;
t78 = pkin(4) * t211;
t77 = qJ(4) * t211;
t75 = qJ(4) * t213;
t74 = t127 * t118 * t130;
t72 = t120 * t100;
t70 = t230 * t127;
t68 = -t118 * t124 - t134;
t67 = qJDD(3) * t130 - t134 * t127;
t66 = qJDD(3) * t127 + t130 * t134;
t60 = t127 * t165;
t58 = -qJDD(3) - t74;
t54 = t202 * t118;
t53 = t206 * t127;
t47 = t64 * t196;
t45 = t155 * t120;
t43 = qJD(3) * t71 - t194;
t42 = -pkin(7) * t196 - qJD(5) * t130 + t100;
t40 = -0.2e1 * t163 + t209;
t39 = 0.2e1 * t163 + t210;
t25 = t249 * t120;
t21 = t149 * t120 - t186;
t17 = t160 * t127 + t217;
t16 = -t102 * t196 + t160 * t130 + t100;
t13 = t21 * t196;
t5 = t72 + (-qJD(5) * t120 - t218) * t130 + t9;
t4 = -t120 * t194 - t176 + (-t178 - t91) * qJ(5) + t183;
t2 = t3 * t130;
t1 = [0, 0, 0, 0, 0, qJDD(1), t157, g(1) * t132 + g(2) * t129, 0, 0, 0, 0, 0, 0, 0, t119, (t119 * t131 - t180) * pkin(1) + t151, ((-qJDD(1) - t119) * t128 + (-qJD(1) - t120) * t198) * pkin(1) + t227, 0, (t157 + (t128 ^ 2 + t131 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t39, t237, t66, t40, t67, 0, t47 + t83 + t138 * t127 + (-t141 + t177) * t130, t141 * t127 + t138 * t130 + t172, t150 + t250, t36 * t103 - t162 - g(2) * (t116 + t228) + t102 * t250 + (t245 * qJD(2) + t233) * pkin(1), t39, t66, -t237, 0, -t67, t40, t13 + t83 + (t147 * qJD(3) - t192) * t127 + (t144 - t232) * t130, t150 + t238, (t192 + (-t147 - t21) * qJD(3)) * t130 + t144 * t127 - t229, t6 * t51 + t21 * t29 - g(2) * t161 + (t152 * t198 + t233) * pkin(1) + t137 * t102 + t139, t39, -t237, -t66, t40, t67, 0, -qJDD(3) * t53 + t2 + t83 + (t154 - t232) * t130 + (-t17 + (-t120 * t38 - t12) * t127) * qJD(3), qJDD(3) * t55 + t154 * t127 + (t38 * t207 + t16) * qJD(3) + t175, (-t119 * t53 - t4 + (-t17 + t217) * t120) * t127 + (-t119 * t55 - t120 * t16 - t5 + (-t120 * t53 - t15) * qJD(3)) * t130 + t191, t5 * t55 + t23 * t16 + t4 * t53 + t15 * t17 + t3 * t38 + t12 * t20 - g(1) * (-t129 * pkin(1) + t170) - g(2) * (t78 + t161) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t151 + t165, (-t193 + (-qJD(2) + t120) * t200) * pkin(1) + t227, 0, 0, t39, t237, t66, t40, t67, 0, t47 + (-pkin(2) * t197 - t220) * t127 + (-t159 + t177) * t130 + t190, -t60 + t159 * t127 + (-t220 + (t186 - t231) * qJD(3)) * t130 + t172, t140 + t250, -t36 * pkin(2) + pkin(7) * t250 - g(2) * t228 - t245 * t226 - t162, t39, t66, -t237, 0, -t67, t40, t13 + (t65 * t197 - t220) * t127 + (t145 - t232) * t130 + t190, t140 + t238, t60 + (t220 + (-t120 * t65 - t186 - t21) * qJD(3)) * t130 + t145 * t127 - t229, t6 * t65 + t21 * t44 - g(2) * t171 + (-t128 * t21 - t152 * t131) * t226 + t137 * pkin(7) + t139, t39, -t237, -t66, t40, t67, 0, -qJDD(3) * t70 + t2 + (t153 - t232) * t130 + (-t43 + (-t12 - t225) * t127) * qJD(3) + t190, qJDD(3) * t71 + t60 + t153 * t127 + (t42 + (-t186 + t225) * t130) * qJD(3) + t175, (-t119 * t70 - t4) * t127 + (-qJD(3) * t15 - t119 * t71 - t5) * t130 + (-t127 * t43 - t130 * t42 + (t127 * t71 - t130 * t70) * qJD(3) + t248) * t120 + t191, t5 * t71 + t23 * t42 + t4 * t70 + t15 * t43 + t3 * t52 + t12 * t28 - g(1) * t170 - g(2) * (t78 + t171) + t136 + (t12 * t128 + (-t127 * t15 - t223) * t131) * t226; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, -t54, t91, t74, t92, qJDD(3), -t64 * t208 + t167, g(3) * t127 - t33 + (-t120 * t64 + t227) * t130, 0, 0, -t74, t91, t54, qJDD(3), -t92, t74, 0.2e1 * t216 + (-t127 * t21 + t130 * t45) * t120 + t156, -t155 * t119 + ((t41 - t123) * t127 + (t166 - t35) * t130) * t120, (t120 * t45 - g(3)) * t127 + (t120 * t21 - t227) * t130 + t185, t9 * qJ(4) - t10 * pkin(3) - t21 * t45 - t63 * t222 - g(1) * (-pkin(3) * t212 + t77) - g(2) * (-pkin(3) * t214 + t75) - g(3) * t203 + (qJD(4) + t49) * t41, -t74, t54, -t91, t74, t92, qJDD(3), qJD(3) * t27 + 0.2e1 * t176 + ((qJ(5) * qJD(3) - t25) * t130 + t240) * t120 - t142, -qJD(3) * t26 + t72 + (-qJD(3) * t63 - t120 * t25 - g(3)) * t127 + (-t204 * t120 - t218 - t227) * t130 + t185, -t249 * t119, -g(1) * t77 - g(2) * t75 - g(3) * t182 + t5 * qJ(4) - t12 * t25 - t15 * t27 + t205 * t23 + t227 * t239 - t236 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t91, t68, -qJD(3) * t41 + t21 * t208 - t148 - t216, 0, 0, 0, 0, 0, 0, t58, t68, -t91, -qJD(3) * t23 - t176 + (-qJ(5) * t195 - t240) * t120 + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92 - 0.2e1 * t179, t91 + 0.2e1 * t178, -t201 * t118, (t223 + (t15 - t181) * t127) * t120 + t143 + t242;];
tau_reg = t1;
