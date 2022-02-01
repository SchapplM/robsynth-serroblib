% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR2
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
% Datum: 2022-01-23 09:19
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:19:07
% EndTime: 2022-01-23 09:19:10
% DurationCPUTime: 1.24s
% Computational Cost: add. (2291->215), mult. (3904->268), div. (0->0), fcn. (2533->16), ass. (0->144)
t138 = cos(pkin(8));
t115 = t138 * pkin(1) + pkin(2);
t136 = sin(pkin(8));
t200 = pkin(1) * t136;
t176 = qJD(3) * t200;
t206 = qJD(1) * t176 - t115 * qJDD(1);
t133 = qJ(1) + pkin(8);
t124 = qJ(3) + t133;
t113 = cos(t124);
t103 = g(2) * t113;
t141 = sin(qJ(3));
t144 = cos(qJ(3));
t188 = pkin(1) * qJDD(1);
t175 = t136 * t188;
t93 = t115 * qJD(1);
t192 = t141 * t93;
t167 = qJD(3) * t192 + t141 * t175 + t206 * t144;
t157 = qJDD(4) + t167;
t128 = qJDD(1) + qJDD(3);
t199 = t128 * pkin(3);
t31 = t157 - t199;
t205 = t31 + t103;
t177 = qJD(1) * t200;
t54 = -t141 * t177 + t144 * t93;
t162 = qJD(4) - t54;
t112 = sin(t124);
t104 = g(1) * t112;
t202 = t104 - t103;
t137 = cos(pkin(9));
t143 = cos(qJ(5));
t183 = t143 * t137;
t135 = sin(pkin(9));
t140 = sin(qJ(5));
t186 = t140 * t135;
t72 = -t183 + t186;
t73 = t143 * t135 + t140 * t137;
t132 = qJD(1) + qJD(3);
t60 = t73 * t132;
t129 = t135 ^ 2;
t130 = t137 ^ 2;
t180 = t129 + t130;
t204 = t132 * t180;
t203 = g(1) * t113 + g(2) * t112;
t201 = t60 ^ 2;
t198 = t137 * pkin(4);
t125 = t137 * pkin(7);
t174 = t132 * t186;
t58 = -t132 * t183 + t174;
t197 = t60 * t58;
t139 = -pkin(7) - qJ(4);
t196 = t205 * t135;
t84 = t139 * t135;
t85 = t137 * qJ(4) + t125;
t49 = -t140 * t85 + t143 * t84;
t195 = t49 * qJD(5) - t162 * t72;
t50 = t140 * t84 + t143 * t85;
t194 = -t50 * qJD(5) - t162 * t73;
t163 = t72 * t128;
t69 = t73 * qJD(5);
t39 = t132 * t69 + t163;
t172 = qJD(5) * t183;
t173 = qJD(5) * t186;
t68 = -t172 + t173;
t193 = -t73 * t39 + t68 * t58;
t55 = t144 * t177 + t192;
t191 = t55 * t132;
t67 = t141 * t115 + t144 * t200;
t63 = t67 * qJD(3);
t190 = t63 * t132;
t189 = t113 * pkin(3) + t112 * qJ(4);
t179 = qJD(3) * t144;
t168 = -t206 * t141 + t144 * t175 + t93 * t179;
t26 = t128 * qJ(4) + t132 * qJD(4) + t168;
t21 = t135 * qJDD(2) + t137 * t26;
t48 = t132 * qJ(4) + t55;
t41 = t135 * qJD(2) + t137 * t48;
t187 = t137 * t128;
t123 = cos(t133);
t145 = cos(qJ(1));
t181 = t145 * pkin(1) + pkin(2) * t123;
t178 = t73 * t128 + t132 * t172;
t114 = pkin(3) + t198;
t171 = -t112 * pkin(3) + t113 * qJ(4);
t117 = t137 * qJDD(2);
t15 = t117 + (-pkin(7) * t128 - t26) * t135;
t16 = pkin(7) * t187 + t21;
t169 = -t140 * t16 + t143 * t15;
t166 = -t112 * t139 + t113 * t114;
t165 = t180 * t128;
t66 = t144 * t115 - t141 * t200;
t65 = -pkin(3) - t66;
t121 = sin(t133);
t142 = sin(qJ(1));
t161 = -t142 * pkin(1) - pkin(2) * t121;
t159 = g(1) * t142 - g(2) * t145;
t38 = t132 * t173 - t178;
t158 = -t72 * t38 + t60 * t69;
t156 = -t191 - t199;
t119 = t137 * qJD(2);
t36 = t119 + (-pkin(7) * t132 - t48) * t135;
t37 = t132 * t125 + t41;
t10 = -t140 * t37 + t143 * t36;
t11 = t140 * t36 + t143 * t37;
t153 = t140 * t15 + t143 * t16;
t3 = t10 * qJD(5) + t153;
t4 = -t11 * qJD(5) + t169;
t155 = t10 * t68 - t11 * t69 - t3 * t72 - t4 * t73 - t203;
t154 = t128 * t65 + t190;
t64 = qJ(4) + t67;
t51 = (-pkin(7) - t64) * t135;
t52 = t137 * t64 + t125;
t23 = -t140 * t52 + t143 * t51;
t24 = t140 * t51 + t143 * t52;
t152 = -t112 * t114 - t113 * t139;
t62 = t115 * t179 - t141 * t176;
t20 = -t135 * t26 + t117;
t151 = -t20 * t135 + t21 * t137 - t203;
t150 = -t168 + t203;
t149 = -t167 + t202;
t131 = pkin(9) + qJ(5);
t120 = sin(t131);
t22 = -t114 * t128 + t157;
t43 = -t114 * t132 + t162;
t148 = -t202 * t120 + t22 * t73 - t43 * t68;
t122 = cos(t131);
t147 = t202 * t122 + t22 * t72 + t43 * t69;
t134 = qJDD(2) - g(3);
t107 = t130 * t128;
t106 = t129 * t128;
t83 = t137 * t104;
t79 = 0.2e1 * t135 * t187;
t57 = t58 ^ 2;
t56 = qJD(4) + t62;
t53 = t65 - t198;
t47 = -t132 * pkin(3) + t162;
t45 = -t69 * qJD(5) - t72 * qJDD(5);
t44 = -t68 * qJD(5) + t73 * qJDD(5);
t40 = -t135 * t48 + t119;
t13 = t39 * t72 + t58 * t69;
t12 = -t38 * t73 - t60 * t68;
t7 = -t24 * qJD(5) - t73 * t56;
t6 = t23 * qJD(5) - t72 * t56;
t5 = -t158 + t193;
t1 = [0, 0, 0, 0, 0, qJDD(1), t159, g(1) * t145 + g(2) * t142, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t121 - g(2) * t123 + 0.2e1 * t138 * t188, g(1) * t123 + g(2) * t121 - 0.2e1 * t175, 0, (t159 + (t136 ^ 2 + t138 ^ 2) * t188) * pkin(1), 0, 0, 0, 0, 0, t128, t66 * t128 + t149 - t190, -t67 * t128 - t62 * t132 + t150, 0, -g(1) * t161 - g(2) * t181 - t167 * t66 + t168 * t67 - t54 * t63 + t55 * t62, t106, t79, 0, t107, 0, 0, t83 + (-t154 - t205) * t137, (t154 - t104) * t135 + t196, t64 * t165 + t56 * t204 + t151, t31 * t65 + t47 * t63 - g(1) * (t161 + t171) - g(2) * (t181 + t189) + (t21 * t64 + t41 * t56) * t137 + (-t20 * t64 - t40 * t56) * t135, t12, t5, t44, t13, t45, 0, t7 * qJD(5) + t23 * qJDD(5) + t53 * t39 + t63 * t58 + t147, -t6 * qJD(5) - t24 * qJDD(5) - t53 * t38 + t63 * t60 + t148, t23 * t38 - t24 * t39 - t6 * t58 - t7 * t60 + t155, t3 * t24 + t11 * t6 + t4 * t23 + t10 * t7 + t22 * t53 + t43 * t63 - g(1) * (t152 + t161) - g(2) * (t166 + t181); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t135 + t20 * t137 - g(3), 0, 0, 0, 0, 0, 0, t45, -t44, t158 + t193, -t10 * t69 - t11 * t68 + t3 * t73 - t4 * t72 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t149 + t191, t54 * t132 + t150, 0, 0, t106, t79, 0, t107, 0, 0, t83 + (-t156 - t205) * t137, (t156 - t104) * t135 + t196, qJ(4) * t165 + t162 * t204 + t151, -t31 * pkin(3) - t47 * t55 - g(1) * t171 - g(2) * t189 + (t21 * qJ(4) + t162 * t41) * t137 + (-t20 * qJ(4) - t162 * t40) * t135, t12, t5, t44, t13, t45, 0, t194 * qJD(5) + t49 * qJDD(5) - t114 * t39 - t55 * t58 + t147, -qJD(5) * t195 - t50 * qJDD(5) + t114 * t38 - t55 * t60 + t148, -t194 * t60 - t195 * t58 + t49 * t38 - t50 * t39 + t155, -g(1) * t152 - g(2) * t166 + t10 * t194 + t11 * t195 - t22 * t114 + t3 * t50 + t4 * t49 - t43 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, t135 * t128, -t180 * t132 ^ 2, (t40 * t135 - t41 * t137) * t132 + t31 - t202, 0, 0, 0, 0, 0, 0, 0.2e1 * t60 * qJD(5) + t163, (-t58 - t174) * qJD(5) + t178, -t57 - t201, t10 * t60 + t11 * t58 - t202 + t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t197, -t57 + t201, (t58 - t174) * qJD(5) + t178, -t197, -t163, qJDD(5), -g(3) * t122 + t120 * t203 - t43 * t60 + t169, g(3) * t120 + t122 * t203 + t43 * t58 - t153, 0, 0;];
tau_reg = t1;
