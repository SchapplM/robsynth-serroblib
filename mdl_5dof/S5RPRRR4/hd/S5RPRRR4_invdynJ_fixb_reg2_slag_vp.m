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
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:52:27
% EndTime: 2020-01-03 11:52:29
% DurationCPUTime: 1.00s
% Computational Cost: add. (2634->196), mult. (4959->246), div. (0->0), fcn. (2884->16), ass. (0->129)
t86 = qJ(1) + pkin(9);
t81 = qJ(3) + t86;
t72 = qJ(4) + t81;
t65 = sin(t72);
t66 = cos(t72);
t122 = -g(2) * t66 - g(3) * t65;
t84 = qJDD(1) + qJDD(3);
t76 = qJDD(4) + t84;
t158 = t76 * pkin(4);
t91 = cos(pkin(9));
t71 = pkin(1) * t91 + pkin(2);
t55 = t71 * qJD(1);
t142 = qJD(3) * t55;
t90 = sin(pkin(9));
t163 = pkin(1) * t90;
t53 = t71 * qJDD(1);
t94 = sin(qJ(3));
t98 = cos(qJ(3));
t105 = -t94 * t142 + t98 * t53 + (-qJD(1) * qJD(3) * t98 - qJDD(1) * t94) * t163;
t21 = t84 * pkin(3) + t105;
t134 = qJD(1) * t163;
t126 = t94 * t134;
t143 = pkin(1) * qJDD(1);
t132 = t90 * t143;
t25 = (t132 + t142) * t98 - qJD(3) * t126 + t94 * t53;
t93 = sin(qJ(4));
t97 = cos(qJ(4));
t130 = -t97 * t21 + t93 * t25;
t37 = t134 * t98 + t55 * t94;
t150 = t97 * t37;
t36 = t98 * t55 - t126;
t85 = qJD(1) + qJD(3);
t32 = pkin(3) * t85 + t36;
t18 = t32 * t93 + t150;
t8 = -qJD(4) * t18 - t130;
t6 = -t158 - t8;
t117 = t122 - t6;
t79 = qJD(4) + t85;
t16 = pkin(8) * t79 + t18;
t92 = sin(qJ(5));
t96 = cos(qJ(5));
t11 = qJD(2) * t96 - t16 * t92;
t155 = t16 * t96;
t12 = qJD(2) * t92 + t155;
t138 = t11 * qJD(5);
t140 = qJD(4) * t97;
t141 = qJD(4) * t93;
t127 = -t32 * t140 + t37 * t141 - t93 * t21 - t97 * t25;
t5 = pkin(8) * t76 - t127;
t2 = t92 * qJDD(2) + t96 * t5 + t138;
t80 = t96 * qJDD(2);
t3 = -t12 * qJD(5) - t92 * t5 + t80;
t108 = t2 * t96 - t3 * t92 + (-t11 * t96 - t12 * t92) * qJD(5);
t45 = -t94 * t163 + t98 * t71;
t44 = pkin(3) + t45;
t46 = t163 * t98 + t71 * t94;
t29 = t93 * t44 + t97 * t46;
t60 = g(3) * t66;
t148 = -g(2) * t65 + t60;
t152 = t37 * t93;
t23 = t36 * t97 - t152;
t165 = pkin(3) * t140 - t23;
t100 = qJD(5) ^ 2;
t22 = t36 * t93 + t150;
t73 = pkin(3) * t93 + pkin(8);
t162 = pkin(3) * t97;
t74 = -pkin(4) - t162;
t164 = t100 * t73 + (pkin(3) * t141 - t22) * t79 + t74 * t76;
t161 = pkin(4) * t79;
t28 = t44 * t97 - t46 * t93;
t42 = t45 * qJD(3);
t43 = t46 * qJD(3);
t9 = qJD(4) * t28 + t97 * t42 - t93 * t43;
t157 = t79 * t9;
t10 = t29 * qJD(4) + t93 * t42 + t97 * t43;
t156 = t10 * t79;
t17 = t32 * t97 - t152;
t154 = t17 * t79;
t153 = t18 * t79;
t151 = t96 * t76;
t149 = t66 * pkin(4) + t65 * pkin(8);
t77 = sin(t86);
t95 = sin(qJ(1));
t147 = t95 * pkin(1) + pkin(2) * t77;
t78 = cos(t86);
t99 = cos(qJ(1));
t146 = t99 * pkin(1) + pkin(2) * t78;
t87 = t92 ^ 2;
t88 = t96 ^ 2;
t145 = t87 - t88;
t144 = t87 + t88;
t139 = qJD(5) * t96;
t89 = qJDD(2) - g(1);
t75 = t79 ^ 2;
t136 = t92 * t75 * t96;
t70 = cos(t81);
t64 = pkin(3) * t70;
t135 = t64 + t146;
t131 = t144 * t76;
t15 = -t17 - t161;
t128 = -t117 * t92 + t15 * t139;
t125 = t79 * t92 * t139;
t58 = t65 * pkin(4);
t69 = sin(t81);
t63 = pkin(3) * t69;
t124 = -pkin(8) * t66 + t58 + t63;
t121 = -g(2) * t70 - g(3) * t69;
t120 = -g(2) * t99 - g(3) * t95;
t119 = t11 * t92 - t12 * t96;
t116 = t127 - t148;
t114 = pkin(8) * t100 - t153 - t158;
t26 = -pkin(4) - t28;
t27 = pkin(8) + t29;
t113 = t100 * t27 + t26 * t76 + t156;
t112 = -pkin(8) * qJDD(5) + (t17 - t161) * qJD(5);
t111 = -qJDD(5) * t27 + (t26 * t79 - t9) * qJD(5);
t109 = -qJD(2) * qJD(5) - t15 * t79 - t148 - t5;
t107 = -qJDD(5) * t73 + (t74 * t79 - t165) * qJD(5);
t106 = t148 + t108;
t104 = g(2) * t69 - g(3) * t70 - t25;
t103 = t122 + t8;
t102 = t105 + t121;
t52 = qJDD(5) * t96 - t100 * t92;
t51 = qJDD(5) * t92 + t100 * t96;
t39 = t76 * t88 - 0.2e1 * t125;
t38 = t76 * t87 + 0.2e1 * t125;
t30 = -0.2e1 * qJD(5) * t145 * t79 + 0.2e1 * t151 * t92;
t13 = t15 * qJD(5) * t92;
t1 = [0, 0, 0, 0, 0, qJDD(1), t120, g(2) * t95 - g(3) * t99, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -g(2) * t78 - g(3) * t77 + 0.2e1 * t143 * t91, g(2) * t77 - g(3) * t78 - 0.2e1 * t132, 0, (t120 + (t90 ^ 2 + t91 ^ 2) * t143) * pkin(1), 0, 0, 0, 0, 0, t84, -t43 * t85 + t45 * t84 + t102, -t42 * t85 - t46 * t84 + t104, 0, -g(2) * t146 - g(3) * t147 + t105 * t45 + t25 * t46 - t36 * t43 + t37 * t42, 0, 0, 0, 0, 0, t76, t28 * t76 + t103 - t156, -t29 * t76 + t116 - t157, 0, -t127 * t29 + t18 * t9 + t8 * t28 - t17 * t10 - g(2) * t135 - g(3) * (t63 + t147), t38, t30, t51, t39, t52, 0, t13 + t111 * t92 + (-t113 + t117) * t96, t111 * t96 + t113 * t92 + t128, t131 * t27 + t144 * t157 + t106, t6 * t26 + t15 * t10 - g(2) * (t135 + t149) - g(3) * (t124 + t147) - t119 * t9 + t108 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, 0, 0, 0, 0, 0, 0, t52, -t51, 0, -qJD(5) * t119 + t2 * t92 + t3 * t96 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t37 * t85 + t102, t36 * t85 + t104, 0, 0, 0, 0, 0, 0, 0, t76, t76 * t162 + t22 * t79 + (-t150 + (-pkin(3) * t79 - t32) * t93) * qJD(4) + t122 - t130, t23 * t79 + (-t140 * t79 - t76 * t93) * pkin(3) + t116, 0, t17 * t22 - t18 * t23 + (-t127 * t93 + t8 * t97 + (-t17 * t93 + t18 * t97) * qJD(4) + t121) * pkin(3), t38, t30, t51, t39, t52, 0, t13 + t107 * t92 + (t117 - t164) * t96, t107 * t96 + t164 * t92 + t128, t165 * t79 * t144 + t73 * t131 + t106, t6 * t74 - t15 * t22 - g(2) * (t64 + t149) - g(3) * t124 + t119 * t23 + (-t119 * t97 + t15 * t93) * qJD(4) * pkin(3) + t108 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t103 + t153, t116 + t154, 0, 0, t38, t30, t51, t39, t52, 0, t13 + t112 * t92 + (-t114 + t117) * t96, t112 * t96 + t114 * t92 + t128, pkin(8) * t131 - t144 * t154 + t106, -t6 * pkin(4) - t15 * t18 - g(2) * t149 - g(3) * t58 + t119 * t17 + (t108 + t60) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, t145 * t75, t92 * t76, t136, t151, qJDD(5), -g(1) * t96 + t80 + (t12 - t155) * qJD(5) + t109 * t92, t138 + (qJD(5) * t16 - t89) * t92 + t109 * t96, 0, 0;];
tau_reg = t1;
