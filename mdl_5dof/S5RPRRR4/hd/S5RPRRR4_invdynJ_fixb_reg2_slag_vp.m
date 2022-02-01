% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:46
% EndTime: 2022-01-23 09:34:48
% DurationCPUTime: 1.08s
% Computational Cost: add. (2634->199), mult. (4959->247), div. (0->0), fcn. (2884->16), ass. (0->131)
t86 = qJ(1) + pkin(9);
t82 = qJ(3) + t86;
t73 = qJ(4) + t82;
t68 = cos(t73);
t163 = g(2) * t68;
t84 = qJDD(1) + qJDD(3);
t77 = qJDD(4) + t84;
t161 = t77 * pkin(4);
t91 = cos(pkin(9));
t72 = t91 * pkin(1) + pkin(2);
t55 = t72 * qJD(1);
t143 = qJD(3) * t55;
t90 = sin(pkin(9));
t164 = pkin(1) * t90;
t53 = t72 * qJDD(1);
t94 = sin(qJ(3));
t98 = cos(qJ(3));
t103 = -t94 * t143 + t98 * t53 + (-qJD(1) * qJD(3) * t98 - qJDD(1) * t94) * t164;
t21 = t84 * pkin(3) + t103;
t134 = qJD(1) * t164;
t125 = t94 * t134;
t144 = pkin(1) * qJDD(1);
t132 = t90 * t144;
t25 = (t132 + t143) * t98 - qJD(3) * t125 + t94 * t53;
t93 = sin(qJ(4));
t97 = cos(qJ(4));
t128 = -t97 * t21 + t93 * t25;
t37 = t98 * t134 + t94 * t55;
t151 = t97 * t37;
t36 = t98 * t55 - t125;
t85 = qJD(1) + qJD(3);
t32 = t85 * pkin(3) + t36;
t18 = t93 * t32 + t151;
t8 = -t18 * qJD(4) - t128;
t6 = -t161 - t8;
t169 = t6 + t163;
t80 = qJD(4) + t85;
t16 = t80 * pkin(8) + t18;
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t11 = t96 * qJD(2) - t92 * t16;
t153 = t96 * t16;
t12 = t92 * qJD(2) + t153;
t139 = t11 * qJD(5);
t141 = qJD(4) * t97;
t142 = qJD(4) * t93;
t126 = -t32 * t141 + t37 * t142 - t93 * t21 - t97 * t25;
t5 = t77 * pkin(8) - t126;
t2 = t92 * qJDD(2) + t96 * t5 + t139;
t81 = t96 * qJDD(2);
t3 = -t12 * qJD(5) - t92 * t5 + t81;
t108 = t2 * t96 - t3 * t92 + (-t11 * t96 - t12 * t92) * qJD(5);
t45 = -t94 * t164 + t98 * t72;
t44 = pkin(3) + t45;
t46 = t98 * t164 + t94 * t72;
t29 = t93 * t44 + t97 * t46;
t67 = sin(t73);
t148 = -g(1) * t68 - g(2) * t67;
t154 = t93 * t37;
t23 = t97 * t36 - t154;
t168 = pkin(3) * t141 - t23;
t61 = g(1) * t67;
t167 = t61 - t163;
t70 = sin(t82);
t64 = g(1) * t70;
t71 = cos(t82);
t166 = -g(2) * t71 + t64;
t100 = qJD(5) ^ 2;
t22 = t93 * t36 + t151;
t74 = t93 * pkin(3) + pkin(8);
t158 = t97 * pkin(3);
t75 = -pkin(4) - t158;
t165 = -t100 * t74 - (pkin(3) * t142 - t22) * t80 - t75 * t77;
t160 = t80 * pkin(4);
t28 = t97 * t44 - t93 * t46;
t42 = t45 * qJD(3);
t43 = t46 * qJD(3);
t9 = t28 * qJD(4) + t97 * t42 - t93 * t43;
t159 = t9 * t80;
t10 = t29 * qJD(4) + t93 * t42 + t97 * t43;
t157 = t10 * t80;
t17 = t97 * t32 - t154;
t156 = t17 * t80;
t155 = t18 * t80;
t152 = t96 * t77;
t15 = -t17 - t160;
t150 = t15 * qJD(5) * t92 + t96 * t61;
t149 = t68 * pkin(4) + t67 * pkin(8);
t79 = cos(t86);
t99 = cos(qJ(1));
t147 = t99 * pkin(1) + pkin(2) * t79;
t87 = t92 ^ 2;
t88 = t96 ^ 2;
t146 = t87 - t88;
t145 = t87 + t88;
t140 = qJD(5) * t96;
t89 = qJDD(2) - g(3);
t137 = t15 * t140 + t169 * t92;
t76 = t80 ^ 2;
t136 = t92 * t76 * t96;
t66 = pkin(3) * t71;
t135 = t66 + t147;
t130 = -t67 * pkin(4) + t68 * pkin(8);
t129 = t145 * t77;
t124 = t92 * t80 * t140;
t122 = g(1) * t130;
t78 = sin(t86);
t95 = sin(qJ(1));
t121 = -t95 * pkin(1) - pkin(2) * t78;
t120 = g(1) * t95 - g(2) * t99;
t119 = t11 * t92 - t12 * t96;
t117 = t126 - t148;
t115 = -pkin(3) * t70 + t121;
t114 = pkin(8) * t100 - t155 - t161;
t26 = -pkin(4) - t28;
t27 = pkin(8) + t29;
t113 = t100 * t27 + t26 * t77 + t157;
t112 = -pkin(8) * qJDD(5) + (t17 - t160) * qJD(5);
t111 = -qJDD(5) * t27 + (t26 * t80 - t9) * qJD(5);
t109 = -qJD(2) * qJD(5) - t15 * t80 - t148 - t5;
t107 = -qJDD(5) * t74 + (t75 * t80 - t168) * qJD(5);
t106 = t148 + t108;
t105 = g(1) * t71 + g(2) * t70 - t25;
t104 = t8 + t167;
t102 = t103 + t166;
t52 = qJDD(5) * t96 - t100 * t92;
t51 = qJDD(5) * t92 + t100 * t96;
t39 = t88 * t77 - 0.2e1 * t124;
t38 = t87 * t77 + 0.2e1 * t124;
t30 = -0.2e1 * t146 * t80 * qJD(5) + 0.2e1 * t92 * t152;
t1 = [0, 0, 0, 0, 0, qJDD(1), t120, g(1) * t99 + g(2) * t95, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(1) * t78 - g(2) * t79 + 0.2e1 * t91 * t144, g(1) * t79 + g(2) * t78 - 0.2e1 * t132, 0, (t120 + (t90 ^ 2 + t91 ^ 2) * t144) * pkin(1), 0, 0, 0, 0, 0, t84, -t43 * t85 + t45 * t84 + t102, -t42 * t85 - t46 * t84 + t105, 0, -g(1) * t121 - g(2) * t147 + t103 * t45 + t25 * t46 - t36 * t43 + t37 * t42, 0, 0, 0, 0, 0, t77, t28 * t77 + t104 - t157, -t29 * t77 + t117 - t159, 0, -g(1) * t115 - g(2) * t135 - t17 * t10 - t126 * t29 + t18 * t9 + t8 * t28, t38, t30, t51, t39, t52, 0, t111 * t92 + (-t113 - t169) * t96 + t150, t111 * t96 + (t113 - t61) * t92 + t137, t27 * t129 + t145 * t159 + t106, t6 * t26 + t15 * t10 - g(1) * (t115 + t130) - g(2) * (t135 + t149) - t119 * t9 + t108 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, t52, -t51, 0, -t119 * qJD(5) + t2 * t92 + t3 * t96 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t37 * t85 + t102, t36 * t85 + t105, 0, 0, 0, 0, 0, 0, 0, t77, t77 * t158 + t22 * t80 + (-t151 + (-pkin(3) * t80 - t32) * t93) * qJD(4) - t128 + t167, t23 * t80 + (-t80 * t141 - t77 * t93) * pkin(3) + t117, 0, t17 * t22 - t18 * t23 + (-t126 * t93 + t8 * t97 + (-t17 * t93 + t18 * t97) * qJD(4) + t166) * pkin(3), t38, t30, t51, t39, t52, 0, t107 * t92 + (-t169 + t165) * t96 + t150, t107 * t96 + (-t165 - t61) * t92 + t137, t168 * t80 * t145 + t74 * t129 + t106, t6 * t75 - t15 * t22 - t122 - g(2) * (t66 + t149) + t119 * t23 + (t64 + (-t119 * t97 + t15 * t93) * qJD(4)) * pkin(3) + t108 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t104 + t155, t117 + t156, 0, 0, t38, t30, t51, t39, t52, 0, t112 * t92 + (-t114 - t169) * t96 + t150, t112 * t96 + (t114 - t61) * t92 + t137, pkin(8) * t129 - t145 * t156 + t106, -t6 * pkin(4) + t108 * pkin(8) - g(2) * t149 + t119 * t17 - t15 * t18 - t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t146 * t76, t92 * t77, t136, t152, qJDD(5), -g(3) * t96 + t81 + (t12 - t153) * qJD(5) + t109 * t92, t139 + (qJD(5) * t16 - t89) * t92 + t109 * t96, 0, 0;];
tau_reg = t1;
