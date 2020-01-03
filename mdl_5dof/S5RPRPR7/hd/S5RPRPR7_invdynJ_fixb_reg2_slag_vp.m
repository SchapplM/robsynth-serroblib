% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR7
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:19:56
% EndTime: 2019-12-31 18:20:00
% DurationCPUTime: 2.38s
% Computational Cost: add. (3645->366), mult. (7873->466), div. (0->0), fcn. (5376->14), ass. (0->193)
t242 = 2 * qJD(3);
t126 = sin(pkin(8));
t103 = t126 * pkin(1) + pkin(6);
t92 = t103 * qJDD(1);
t168 = -(qJD(2) * qJD(3)) - t92;
t129 = sin(qJ(5));
t125 = sin(pkin(9));
t130 = sin(qJ(3));
t194 = qJD(1) * t130;
t133 = cos(qJ(3));
t203 = cos(pkin(9));
t175 = t203 * t133;
t98 = qJD(1) * t175;
t78 = t125 * t194 - t98;
t75 = qJD(5) + t78;
t173 = t129 * t75;
t132 = cos(qJ(5));
t86 = t125 * t133 + t203 * t130;
t81 = t86 * qJD(1);
t63 = t129 * qJD(3) + t132 * t81;
t241 = t63 * t173;
t202 = pkin(1) * qJDD(1);
t121 = qJ(3) + pkin(9);
t109 = sin(t121);
t122 = qJ(1) + pkin(8);
t110 = sin(t122);
t112 = cos(t122);
t176 = -g(1) * t110 + g(2) * t112;
t240 = t176 * t109;
t94 = t103 * qJD(1);
t170 = qJ(4) * qJD(1) + t94;
t239 = t170 * t133;
t127 = cos(pkin(8));
t105 = -t127 * pkin(1) - pkin(2);
t119 = t133 * pkin(3);
t238 = t105 - t119;
t111 = cos(t121);
t237 = t111 * pkin(4) + t109 * pkin(7);
t162 = g(1) * t112 + g(2) * t110;
t102 = t125 * pkin(3) + pkin(7);
t143 = -g(3) * t111 + t162 * t109;
t116 = t133 * qJDD(2);
t187 = qJD(1) * qJD(4);
t27 = qJDD(3) * pkin(3) + t116 - qJD(3) * t239 + (-qJ(4) * qJDD(1) + t168 - t187) * t130;
t188 = qJD(1) * qJD(3);
t177 = t130 * t188;
t184 = t133 * qJDD(1);
t182 = -t130 * qJDD(2) + t168 * t133;
t193 = qJD(3) * t130;
t43 = -t94 * t193 - t182;
t30 = t133 * t187 + (-t177 + t184) * qJ(4) + t43;
t7 = -t125 * t30 + t203 * t27;
t5 = -qJDD(3) * pkin(4) - t7;
t236 = -qJD(5) * t102 * t75 + t143 - t5;
t227 = g(3) * t109;
t235 = -t162 * t111 - t227;
t185 = t130 * qJDD(1);
t156 = -qJDD(1) * t175 + t125 * t185;
t80 = t86 * qJD(3);
t50 = qJD(1) * t80 + t156;
t49 = qJDD(5) + t50;
t148 = -t125 * t130 + t175;
t83 = t148 * qJD(3);
t158 = t49 * t86 + t75 * t83;
t192 = qJD(5) * t129;
t180 = t86 * t192;
t234 = -t158 * t132 + t75 * t180;
t233 = t81 ^ 2;
t190 = t130 * qJD(2);
t66 = t190 + t239;
t57 = t203 * t66;
t118 = t133 * qJD(2);
t65 = -t170 * t130 + t118;
t59 = qJD(3) * pkin(3) + t65;
t29 = t125 * t59 + t57;
t25 = qJD(3) * pkin(7) + t29;
t77 = qJD(1) * t238 + qJD(4);
t35 = t78 * pkin(4) - t81 * pkin(7) + t77;
t9 = -t129 * t25 + t132 * t35;
t232 = t9 * t75;
t231 = pkin(3) * t130;
t225 = g(3) * t133;
t10 = t129 * t35 + t132 * t25;
t224 = t10 * t75;
t131 = sin(qJ(1));
t221 = t131 * pkin(1);
t189 = t132 * qJD(3);
t61 = t129 * t81 - t189;
t220 = t61 * t78;
t219 = t61 * t83;
t218 = t63 * t61;
t217 = t63 * t81;
t216 = t63 * t83;
t215 = t81 * t61;
t214 = t81 * t78;
t142 = t86 * qJDD(1) - t125 * t177;
t51 = qJD(3) * t98 + t142;
t174 = -t132 * qJDD(3) + t129 * t51;
t21 = t63 * qJD(5) + t174;
t205 = t21 * t132;
t213 = -t132 * t219 - t86 * t205;
t191 = qJD(5) * t132;
t212 = -t129 * t21 - t61 * t191;
t8 = t125 * t27 + t203 * t30;
t211 = -t86 * t50 - t83 * t78;
t20 = -qJD(5) * t189 - t129 * qJDD(3) - t132 * t51 + t81 * t192;
t210 = t148 * t20 + t63 * t80;
t209 = t125 * t66;
t208 = t129 * t49;
t207 = t130 * t94;
t206 = t133 * t94;
t108 = t119 + pkin(2);
t134 = cos(qJ(1));
t120 = t134 * pkin(1);
t204 = t112 * t108 + t120;
t201 = qJD(5) * t61;
t200 = t110 * t129;
t199 = t110 * t132;
t198 = t112 * t129;
t197 = t112 * t132;
t196 = qJ(4) + t103;
t123 = t130 ^ 2;
t124 = t133 ^ 2;
t195 = t123 - t124;
t95 = qJD(1) * t105;
t93 = qJDD(1) * t105;
t183 = pkin(3) * t193;
t179 = t86 * t191;
t136 = qJD(1) ^ 2;
t178 = t130 * t136 * t133;
t172 = t132 * t75;
t171 = t196 * t130;
t169 = qJD(3) * t196;
t166 = t63 * t179;
t165 = t133 * t177;
t28 = t203 * t59 - t209;
t24 = -qJD(3) * pkin(4) - t28;
t164 = t24 * t83 + t5 * t86;
t160 = g(1) * t131 - g(2) * t134;
t159 = t148 * t21 - t80 * t61;
t157 = t148 * t51 - t81 * t80;
t155 = t10 * t132 - t129 * t9;
t154 = -t10 * t129 - t132 * t9;
t128 = -qJ(4) - pkin(6);
t153 = -t112 * t128 - t221;
t41 = -pkin(4) * t148 - t86 * pkin(7) + t238;
t84 = t196 * t133;
t47 = -t125 * t171 + t203 * t84;
t16 = -t129 * t47 + t132 * t41;
t17 = t129 * t41 + t132 * t47;
t73 = t190 + t206;
t151 = t132 * t49 - t78 * t173 - t75 * t192;
t6 = qJDD(3) * pkin(7) + t8;
t150 = -qJD(5) * t35 + t227 - t6;
t149 = -t102 * t49 + t75 * t24;
t146 = -qJD(1) * t95 + t162;
t145 = -t130 * qJD(4) - t133 * t169;
t144 = -qJDD(3) * t103 + t95 * t242;
t64 = pkin(3) * t177 + qJDD(1) * t238 + qJDD(4);
t135 = qJD(3) ^ 2;
t141 = -t103 * t135 - t176 - 0.2e1 * t93;
t140 = -t158 * t129 - t75 * t179;
t14 = t50 * pkin(4) - t51 * pkin(7) + t64;
t1 = qJD(5) * t9 + t129 * t14 + t132 * t6;
t13 = t132 * t14;
t2 = -qJD(5) * t10 - t129 * t6 + t13;
t139 = t154 * qJD(5) + t1 * t132 - t2 * t129;
t44 = -t73 * qJD(3) - t130 * t92 + t116;
t72 = t118 - t207;
t138 = -t44 * t130 + t43 * t133 + (-t130 * t73 - t133 * t72) * qJD(3);
t104 = -t203 * pkin(3) - pkin(4);
t91 = qJDD(3) * t133 - t135 * t130;
t90 = qJDD(3) * t130 + t135 * t133;
t76 = t78 ^ 2;
t71 = t111 * t197 + t200;
t70 = -t111 * t198 + t199;
t69 = -t111 * t199 + t198;
t68 = t111 * t200 + t197;
t67 = t133 * qJD(4) - t130 * t169;
t53 = t83 * qJD(3) + t86 * qJDD(3);
t52 = -t80 * qJD(3) + qJDD(3) * t148;
t46 = t125 * t84 + t203 * t171;
t39 = t80 * pkin(4) - t83 * pkin(7) + t183;
t38 = pkin(3) * t194 + t81 * pkin(4) + t78 * pkin(7);
t34 = t125 * t145 + t203 * t67;
t33 = t125 * t67 - t203 * t145;
t32 = t203 * t65 - t209;
t31 = t125 * t65 + t57;
t12 = t129 * t38 + t132 * t32;
t11 = -t129 * t32 + t132 * t38;
t4 = -t17 * qJD(5) - t129 * t34 + t132 * t39;
t3 = t16 * qJD(5) + t129 * t39 + t132 * t34;
t15 = [0, 0, 0, 0, 0, qJDD(1), t160, g(1) * t134 + g(2) * t131, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t127 * t202 - t176, -0.2e1 * t126 * t202 + t162, 0, (t160 + (t126 ^ 2 + t127 ^ 2) * t202) * pkin(1), t123 * qJDD(1) + 0.2e1 * t165, 0.2e1 * t130 * t184 - 0.2e1 * t195 * t188, t90, t124 * qJDD(1) - 0.2e1 * t165, t91, 0, t144 * t130 + t141 * t133, -t141 * t130 + t144 * t133, (t123 + t124) * t92 + t138 - t162, t93 * t105 - g(1) * (-t110 * pkin(2) + t112 * pkin(6) - t221) - g(2) * (t112 * pkin(2) + t110 * pkin(6) + t120) + t138 * t103, t51 * t86 + t81 * t83, t157 + t211, t53, -t148 * t50 + t78 * t80, t52, 0, -t46 * qJDD(3) + t238 * t50 - t64 * t148 + t77 * t80 - t176 * t111 + (t78 * t231 - t33) * qJD(3), -t47 * qJDD(3) + t238 * t51 + t64 * t86 + t77 * t83 + t240 + (t81 * t231 - t34) * qJD(3), t148 * t8 - t28 * t83 - t29 * t80 + t33 * t81 - t34 * t78 + t46 * t51 - t47 * t50 - t7 * t86 - t162, t8 * t47 + t29 * t34 - t7 * t46 - t28 * t33 + t64 * t238 + t77 * t183 - g(1) * (-t110 * t108 + t153) - g(2) * (-t110 * t128 + t204), -t63 * t180 + (-t20 * t86 + t216) * t132, -t166 + (-t216 + (t20 + t201) * t86) * t129 + t213, t210 - t234, t61 * t179 + (t21 * t86 + t219) * t129, t140 + t159, -t148 * t49 + t75 * t80, -g(1) * t69 - g(2) * t71 + t129 * t164 - t148 * t2 + t16 * t49 + t24 * t179 + t46 * t21 + t33 * t61 + t4 * t75 + t9 * t80, -g(1) * t68 - g(2) * t70 + t1 * t148 - t10 * t80 + t164 * t132 - t17 * t49 - t24 * t180 - t46 * t20 - t3 * t75 + t33 * t63, t16 * t20 - t17 * t21 - t3 * t61 - t4 * t63 + t154 * t83 - t240 + (-qJD(5) * t155 - t1 * t129 - t2 * t132) * t86, t1 * t17 + t10 * t3 + t2 * t16 + t9 * t4 + t5 * t46 + t24 * t33 - g(1) * t153 - g(2) * (t237 * t112 + t204) + (-g(1) * (-t108 - t237) + g(2) * t128) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t91, -t90, 0, t43 * t130 + t44 * t133 - g(3) + (-t130 * t72 + t133 * t73) * qJD(3), 0, 0, 0, 0, 0, 0, t52, -t53, -t157 + t211, t148 * t7 - t28 * t80 + t29 * t83 + t8 * t86 - g(3), 0, 0, 0, 0, 0, 0, t140 - t159, t210 + t234, t166 + (t216 + (-t20 + t201) * t86) * t129 + t213, t139 * t86 - t148 * t5 + t155 * t83 + t24 * t80 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t195 * t136, t185, t178, t184, qJDD(3), -t225 + t116 + (t73 - t206) * qJD(3) + (t146 + t168) * t130, g(3) * t130 + (t72 + t207) * qJD(3) + t146 * t133 + t182, 0, 0, t214, -t76 + t233, (t98 + t78) * qJD(3) + t142, -t214, -t156, qJDD(3), t31 * qJD(3) - t77 * t81 + (t203 * qJDD(3) - t78 * t194) * pkin(3) + t143 + t7, t32 * qJD(3) + t77 * t78 + (-qJDD(3) * t125 - t81 * t194) * pkin(3) - t8 - t235, (t29 - t31) * t81 + (-t28 + t32) * t78 + (-t125 * t50 - t203 * t51) * pkin(3), t28 * t31 - t29 * t32 + (t203 * t7 - t225 + t125 * t8 + (-qJD(1) * t77 + t162) * t130) * pkin(3), -t20 * t129 + t172 * t63, (-t20 - t220) * t132 - t241 + t212, t172 * t75 + t208 - t217, t61 * t173 - t205, t151 + t215, -t75 * t81, t104 * t21 - t11 * t75 + t149 * t129 + t236 * t132 - t31 * t61 - t9 * t81, t10 * t81 - t104 * t20 + t12 * t75 - t236 * t129 + t149 * t132 - t31 * t63, t11 * t63 + t12 * t61 + (-t102 * t21 - t78 * t9 + t1 + (t102 * t63 - t9) * qJD(5)) * t132 + (-t10 * t78 - t102 * t20 - t2 + (t102 * t61 - t10) * qJD(5)) * t129 + t235, t5 * t104 - t10 * t12 - t9 * t11 - t24 * t31 - g(3) * (t119 + t237) + t139 * t102 + t162 * (pkin(4) * t109 - pkin(7) * t111 + t231); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81 * t242 + t156, (t98 - t78) * qJD(3) + t142, -t76 - t233, t28 * t81 + t29 * t78 + t176 + t64, 0, 0, 0, 0, 0, 0, t151 - t215, -t75 ^ 2 * t132 - t208 - t217, (t20 - t220) * t132 + t241 + t212, -t24 * t81 + (t2 + t224) * t132 + (t1 - t232) * t129 + t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t61 ^ 2 + t63 ^ 2, t61 * t75 - t20, -t218, -t174 + (-qJD(5) + t75) * t63, t49, -g(1) * t70 + g(2) * t68 + t129 * t150 - t25 * t191 - t24 * t63 + t13 + t224, g(1) * t71 - g(2) * t69 + t24 * t61 + t232 + (qJD(5) * t25 - t14) * t129 + t150 * t132, 0, 0;];
tau_reg = t15;
