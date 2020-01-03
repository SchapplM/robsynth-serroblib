% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:46
% EndTime: 2019-12-31 20:15:48
% DurationCPUTime: 1.60s
% Computational Cost: add. (3063->285), mult. (3934->336), div. (0->0), fcn. (2272->12), ass. (0->173)
t128 = qJ(1) + qJ(2);
t115 = sin(t128);
t117 = cos(t128);
t227 = -g(1) * t115 + g(2) * t117;
t129 = sin(qJ(5));
t130 = sin(qJ(4));
t192 = qJD(5) * t129;
t194 = qJD(4) * t130;
t232 = -t129 * t194 - t130 * t192;
t124 = qJD(1) + qJD(2);
t131 = sin(qJ(2));
t217 = pkin(1) * qJD(1);
t186 = t131 * t217;
t206 = qJ(3) * t124;
t73 = t186 + t206;
t231 = -t124 * t73 + t227;
t118 = t130 * pkin(4);
t137 = -pkin(8) - pkin(7);
t230 = t115 * t137 + t117 * t118;
t125 = t130 ^ 2;
t134 = cos(qJ(4));
t126 = t134 ^ 2;
t197 = t125 + t126;
t138 = -pkin(2) - pkin(7);
t135 = cos(qJ(2));
t188 = t135 * t217;
t164 = qJD(3) - t188;
t62 = t138 * t124 + t164;
t229 = t197 * t62;
t133 = cos(qJ(5));
t71 = t129 * t134 + t130 * t133;
t56 = t71 * t124;
t208 = pkin(1) * qJDD(1);
t122 = qJDD(1) + qJDD(2);
t228 = t122 * t138;
t99 = qJ(3) + t118;
t163 = -g(1) * t117 - g(2) * t115;
t72 = -t129 * t130 + t133 * t134;
t123 = qJD(4) + qJD(5);
t201 = t134 * t122;
t184 = qJD(2) * t217;
t209 = t131 * t184 - t135 * t208;
t178 = qJDD(3) + t209;
t48 = t178 + t228;
t43 = t134 * t48;
t15 = -t62 * t194 + qJDD(4) * pkin(4) + t43 + (t124 * t194 - t201) * pkin(8);
t193 = qJD(4) * t134;
t181 = t124 * t193;
t203 = t130 * t122;
t151 = t181 + t203;
t16 = -pkin(8) * t151 + t130 * t48 + t62 * t193;
t205 = t124 * t134;
t41 = -pkin(8) * t205 + t134 * t62;
t37 = qJD(4) * pkin(4) + t41;
t40 = (-pkin(8) * t124 + t62) * t130;
t4 = (qJD(5) * t37 + t16) * t133 + t129 * t15 - t40 * t192;
t214 = t129 * t40;
t17 = t133 * t37 - t214;
t213 = t133 * t40;
t18 = t129 * t37 + t213;
t38 = t123 * t71;
t173 = t123 * t134;
t39 = t133 * t173 + t232;
t5 = -qJD(5) * t18 - t129 * t16 + t133 * t15;
t159 = t17 * t38 - t18 * t39 - t4 * t71 - t5 * t72 - t227;
t120 = t124 ^ 2;
t108 = -pkin(1) * t135 - pkin(2);
t95 = -pkin(7) + t108;
t226 = -pkin(8) + t95;
t225 = pkin(2) * t122;
t132 = sin(qJ(1));
t224 = g(1) * t132;
t223 = g(3) * t130;
t58 = t72 * t124;
t222 = t58 * t56;
t221 = -pkin(8) + t138;
t76 = t221 * t130;
t77 = t221 * t134;
t44 = -t129 * t76 + t133 * t77;
t111 = pkin(8) * t194;
t67 = -t138 * t194 + t111;
t68 = qJD(4) * t77;
t220 = qJD(5) * t44 + t129 * t67 + t133 * t68 - t71 * t186;
t45 = t129 * t77 + t133 * t76;
t219 = -qJD(5) * t45 - t129 * t68 + t133 * t67 - t72 * t186;
t121 = qJDD(4) + qJDD(5);
t28 = t72 * t121 - t38 * t123;
t207 = qJ(3) * t122;
t210 = t131 * t208 + t135 * t184;
t51 = qJD(3) * t124 + t207 + t210;
t218 = t51 * t130 + t73 * t193;
t59 = t99 * t124 + t186;
t216 = t124 * t59;
t212 = t135 * t73;
t211 = t117 * pkin(2) + t115 * qJ(3);
t139 = qJD(4) ^ 2;
t199 = -t120 - t139;
t198 = t125 - t126;
t196 = qJD(2) * t131;
t195 = qJD(2) * t135;
t191 = qJDD(4) * t95;
t190 = qJDD(4) * t130;
t189 = qJDD(4) * t138;
t187 = pkin(1) * t196;
t136 = cos(qJ(1));
t119 = t136 * pkin(1);
t185 = t119 + t211;
t183 = t134 * t120 * t130;
t66 = t226 * t134;
t182 = t124 * t196;
t97 = t117 * qJ(3);
t177 = -pkin(2) * t115 + t97;
t98 = pkin(1) * t131 + qJ(3);
t175 = t197 * t48;
t172 = t124 * t186;
t171 = pkin(1) * t182;
t170 = -t210 - t163;
t169 = t227 + t209;
t168 = t130 * t181;
t112 = pkin(4) * t193;
t167 = t112 + t164;
t166 = t129 * t203 - t133 * t201;
t87 = pkin(1) * t195 + qJD(3);
t161 = -g(2) * t136 + t224;
t23 = t124 * t38 + t166;
t6 = -t23 * t72 - t38 * t58;
t157 = t232 * t124 + t129 * t201;
t24 = (t124 * t173 + t203) * t133 + t157;
t7 = t24 * t71 + t39 * t56;
t160 = t51 * t98 + t73 * t87;
t29 = -t121 * t71 - t123 * t39;
t65 = t226 * t130;
t34 = t129 * t66 + t133 * t65;
t33 = -t129 * t65 + t133 * t66;
t158 = t51 * qJ(3) + t73 * qJD(3);
t156 = t115 * t118 - t117 * t137 + t211;
t155 = -t132 * pkin(1) + t177;
t55 = t178 - t225;
t154 = g(1) * (t138 * t115 + t97);
t153 = t124 * t98 + t187;
t152 = -t186 + t206;
t147 = -t169 + t172;
t127 = qJ(4) + qJ(5);
t114 = sin(t127);
t32 = pkin(4) * t151 + t51;
t146 = t163 * t114 + t32 * t71 + t59 * t39;
t116 = cos(t127);
t145 = t163 * t116 + t32 * t72 - t59 * t38;
t144 = t122 * t98 + t124 * t87 - t139 * t95 + t163;
t143 = g(3) * t116 - t227 * t114 + t56 * t59 - t4;
t142 = t164 * t124 - t138 * t139 + t163 + t207;
t141 = g(3) * t114 + t227 * t116 - t58 * t59 + t5;
t113 = qJDD(4) * t134;
t100 = t117 * pkin(7);
t82 = -t130 * t139 + t113;
t81 = -t134 * t139 - t190;
t80 = t118 + t98;
t70 = -pkin(2) * t124 + t164;
t69 = t112 + t87;
t61 = t122 * t126 - 0.2e1 * t168;
t60 = t122 * t125 + 0.2e1 * t168;
t50 = qJD(4) * t66 + t130 * t187;
t49 = t134 * t187 - t95 * t194 + t111;
t47 = t51 * t134;
t42 = 0.2e1 * t198 * t124 * qJD(4) - 0.2e1 * t130 * t201;
t25 = -t56 ^ 2 + t58 ^ 2;
t20 = t133 * t41 - t214;
t19 = -t129 * t41 - t213;
t13 = t123 * t58 + (-t123 * t205 - t203) * t133 - t157;
t11 = -qJD(5) * t34 - t129 * t50 + t133 * t49;
t10 = qJD(5) * t33 + t129 * t49 + t133 * t50;
t1 = t23 * t71 - t24 * t72 + t38 * t56 - t39 * t58;
t2 = [0, 0, 0, 0, 0, qJDD(1), t161, g(1) * t136 + g(2) * t132, 0, 0, 0, 0, 0, 0, 0, t122, (t122 * t135 - t182) * pkin(1) - t169, (-t122 * t131 - t124 * t195) * pkin(1) + t170, 0, (t161 + (t131 ^ 2 + t135 ^ 2) * t208) * pkin(1), t122, 0, 0, 0, 0, 0, 0, t171 + qJDD(3) + (-pkin(2) + t108) * t122 + t169, (qJD(3) + t87) * t124 + (qJ(3) + t98) * t122 - t170, -g(1) * t155 - g(2) * t185 + t55 * t108 + t70 * t187 + t160, t61, t42, t82, t60, t81, 0, (t153 * qJD(4) + t191) * t134 + t144 * t130 + t218, t47 + (-t191 + (-t153 - t73) * qJD(4)) * t130 + t144 * t134, -t227 + t197 * (-t122 * t95 - t171 - t48), -t154 - g(2) * (t100 + t185) + t95 * t175 + (t196 * t229 + t224) * pkin(1) + t160, t6, t1, t28, t7, t29, 0, t11 * t123 + t121 * t33 + t24 * t80 + t56 * t69 + t146, -t10 * t123 - t121 * t34 - t23 * t80 + t58 * t69 + t145, -t10 * t56 - t11 * t58 + t23 * t33 - t24 * t34 + t159, t4 * t34 + t18 * t10 + t5 * t33 + t17 * t11 + t32 * t80 + t59 * t69 - g(1) * (t155 + t230) - g(2) * (t119 + t156); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t147, t124 * t188 + t170, 0, 0, t122, 0, 0, 0, 0, 0, 0, qJDD(3) - t147 - 0.2e1 * t225, 0.2e1 * t207 + (0.2e1 * qJD(3) - t188) * t124 - t170, -t55 * pkin(2) - g(1) * t177 - g(2) * t211 + (-t131 * t70 - t212) * t217 + t158, t61, t42, t82, t60, t81, 0, (t152 * qJD(4) + t189) * t134 + t142 * t130 + t218, t47 + (-t189 + (-t152 - t73) * qJD(4)) * t130 + t142 * t134, -t227 + t197 * (t172 - t48 - t228), -t154 - g(2) * (t100 + t211) + t138 * t175 + (-t131 * t229 - t212) * t217 + t158, t6, t1, t28, t7, t29, 0, t121 * t44 + t219 * t123 + t167 * t56 + t24 * t99 + t146, -t121 * t45 - t220 * t123 + t167 * t58 - t23 * t99 + t145, -t219 * t58 - t220 * t56 + t23 * t44 - t24 * t45 + t159, t4 * t45 + t5 * t44 + t32 * t99 - g(1) * (t177 + t230) - g(2) * t156 + t167 * t59 + t220 * t18 + t219 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, -t120, t231 + t55, 0, 0, 0, 0, 0, 0, t199 * t130 + t113, t199 * t134 - t190, -t197 * t122, t175 + t231, 0, 0, 0, 0, 0, 0, -t124 * t56 + t28, -t124 * t58 + t29, -t6 - t7, -t216 - t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t183, -t198 * t120, t201, -t183, -t203, qJDD(4), t134 * t231 + t223 + t43, g(3) * t134 + (-t231 - t48) * t130, 0, 0, t222, t25, -t166, -t222, t13, t121, -t123 * t19 + (t121 * t133 - t123 * t192 - t56 * t205) * pkin(4) + t141, t123 * t20 + (-qJD(5) * t123 * t133 - t121 * t129 - t58 * t205) * pkin(4) + t143, (t18 + t19) * t58 + (-t17 + t20) * t56 + (-t129 * t24 + t133 * t23 + (t129 * t58 - t133 * t56) * qJD(5)) * pkin(4), -t17 * t19 - t18 * t20 + (t223 + t129 * t4 + t133 * t5 + (-t129 * t17 + t133 * t18) * qJD(5) + (t227 - t216) * t134) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t25, -t166, -t222, t13, t121, t123 * t18 + t141, t123 * t17 + t143, 0, 0;];
tau_reg = t2;
