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
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:14:42
% EndTime: 2019-12-05 18:14:44
% DurationCPUTime: 0.99s
% Computational Cost: add. (2634->198), mult. (4959->247), div. (0->0), fcn. (2884->16), ass. (0->127)
t82 = qJ(1) + pkin(9);
t79 = qJ(3) + t82;
t70 = qJ(4) + t79;
t65 = sin(t70);
t66 = cos(t70);
t119 = g(2) * t66 + g(3) * t65;
t86 = sin(pkin(9));
t158 = pkin(1) * t86;
t133 = qJD(1) * t158;
t87 = cos(pkin(9));
t69 = pkin(1) * t87 + pkin(2);
t55 = t69 * qJD(1);
t90 = sin(qJ(3));
t94 = cos(qJ(3));
t37 = t94 * t133 + t55 * t90;
t93 = cos(qJ(4));
t146 = t93 * t37;
t125 = t90 * t133;
t36 = t94 * t55 - t125;
t81 = qJD(1) + qJD(3);
t32 = pkin(3) * t81 + t36;
t89 = sin(qJ(4));
t18 = t32 * t89 + t146;
t77 = qJD(4) + t81;
t16 = pkin(8) * t77 + t18;
t88 = sin(qJ(5));
t92 = cos(qJ(5));
t11 = qJD(2) * t92 - t16 * t88;
t151 = t16 * t92;
t12 = qJD(2) * t88 + t151;
t138 = qJD(5) * t11;
t139 = qJD(4) * t93;
t140 = qJD(4) * t89;
t141 = qJD(3) * t55;
t53 = t69 * qJDD(1);
t100 = -t90 * t141 + t94 * t53 + (-qJD(1) * qJD(3) * t94 - qJDD(1) * t90) * t158;
t80 = qJDD(1) + qJDD(3);
t21 = t80 * pkin(3) + t100;
t142 = pkin(1) * qJDD(1);
t131 = t86 * t142;
t25 = (t131 + t141) * t94 - qJD(3) * t125 + t90 * t53;
t126 = -t32 * t139 + t37 * t140 - t89 * t21 - t93 * t25;
t74 = qJDD(4) + t80;
t5 = pkin(8) * t74 - t126;
t2 = t88 * qJDD(2) + t92 * t5 + t138;
t78 = t92 * qJDD(2);
t3 = -qJD(5) * t12 - t88 * t5 + t78;
t104 = t2 * t92 - t3 * t88 + (-t11 * t92 - t12 * t88) * qJD(5);
t45 = -t90 * t158 + t94 * t69;
t44 = pkin(3) + t45;
t46 = t94 * t158 + t69 * t90;
t29 = t89 * t44 + t93 * t46;
t60 = g(2) * t65;
t145 = -g(3) * t66 + t60;
t67 = sin(t79);
t68 = cos(t79);
t118 = g(2) * t68 + g(3) * t67;
t148 = t37 * t89;
t23 = t36 * t93 - t148;
t160 = pkin(3) * t139 - t23;
t22 = t36 * t89 + t146;
t71 = pkin(3) * t89 + pkin(8);
t157 = pkin(3) * t93;
t72 = -pkin(4) - t157;
t96 = qJD(5) ^ 2;
t159 = (pkin(3) * t140 - t22) * t77 + t71 * t96 + t72 * t74;
t128 = -t93 * t21 + t89 * t25;
t8 = -qJD(4) * t18 - t128;
t156 = pkin(4) * t77;
t155 = t74 * pkin(4);
t28 = t44 * t93 - t46 * t89;
t42 = t45 * qJD(3);
t43 = t46 * qJD(3);
t9 = t28 * qJD(4) + t93 * t42 - t89 * t43;
t154 = t77 * t9;
t137 = qJD(5) * t92;
t17 = t32 * t93 - t148;
t15 = -t17 - t156;
t6 = -t155 - t8;
t153 = t15 * t137 + t6 * t88;
t10 = t29 * qJD(4) + t89 * t42 + t93 * t43;
t152 = t10 * t77;
t150 = t17 * t77;
t149 = t18 * t77;
t147 = t92 * t74;
t83 = t88 ^ 2;
t84 = t92 ^ 2;
t144 = t83 - t84;
t143 = t83 + t84;
t85 = qJDD(2) - g(1);
t73 = t77 ^ 2;
t135 = t88 * t73 * t92;
t134 = t15 * qJD(5) * t88 + t119 * t92;
t57 = t66 * pkin(8);
t130 = -pkin(4) * t65 + t57;
t129 = t143 * t74;
t124 = t77 * t88 * t137;
t75 = sin(t82);
t91 = sin(qJ(1));
t122 = -pkin(1) * t91 - pkin(2) * t75;
t76 = cos(t82);
t95 = cos(qJ(1));
t121 = -pkin(1) * t95 - pkin(2) * t76;
t120 = -pkin(4) * t66 - pkin(8) * t65;
t117 = g(2) * t95 + g(3) * t91;
t116 = t11 * t88 - t12 * t92;
t114 = t126 - t145;
t112 = -pkin(3) * t67 + t122;
t111 = -pkin(3) * t68 + t121;
t110 = -pkin(8) * t96 + t149 + t155;
t26 = -pkin(4) - t28;
t27 = pkin(8) + t29;
t109 = -t26 * t74 - t27 * t96 - t152;
t108 = -pkin(8) * qJDD(5) + (t17 - t156) * qJD(5);
t107 = -qJDD(5) * t27 + (t26 * t77 - t9) * qJD(5);
t105 = -qJD(2) * qJD(5) - t15 * t77 - t145 - t5;
t103 = -qJDD(5) * t71 + (t72 * t77 - t160) * qJD(5);
t102 = t145 + t104;
t101 = t8 + t119;
t99 = -g(2) * t67 + g(3) * t68 - t25;
t98 = t100 + t118;
t52 = qJDD(5) * t92 - t88 * t96;
t51 = qJDD(5) * t88 + t92 * t96;
t39 = t74 * t84 - 0.2e1 * t124;
t38 = t74 * t83 + 0.2e1 * t124;
t30 = -0.2e1 * t144 * t77 * qJD(5) + 0.2e1 * t88 * t147;
t1 = [0, 0, 0, 0, 0, qJDD(1), t117, -g(2) * t91 + g(3) * t95, 0, 0, 0, 0, 0, 0, 0, qJDD(1), g(2) * t76 + g(3) * t75 + 0.2e1 * t87 * t142, -g(2) * t75 + g(3) * t76 - 0.2e1 * t131, 0, (t117 + (t86 ^ 2 + t87 ^ 2) * t142) * pkin(1), 0, 0, 0, 0, 0, t80, -t43 * t81 + t45 * t80 + t98, -t42 * t81 - t46 * t80 + t99, 0, -g(2) * t121 - g(3) * t122 + t100 * t45 + t25 * t46 - t36 * t43 + t37 * t42, 0, 0, 0, 0, 0, t74, t28 * t74 + t101 - t152, -t29 * t74 + t114 - t154, 0, -g(2) * t111 - g(3) * t112 - t17 * t10 - t126 * t29 + t18 * t9 + t8 * t28, t38, t30, t51, t39, t52, 0, t107 * t88 + (t109 - t6) * t92 + t134, t107 * t92 + (-t109 - t119) * t88 + t153, t27 * t129 + t143 * t154 + t102, t6 * t26 + t15 * t10 - g(2) * (t111 + t120) - g(3) * (t112 + t130) - t116 * t9 + t104 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, 0, 0, 0, 0, 0, t52, -t51, 0, -t116 * qJD(5) + t2 * t88 + t3 * t92 - g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t37 * t81 + t98, t36 * t81 + t99, 0, 0, 0, 0, 0, 0, 0, t74, t74 * t157 + t22 * t77 + (-t146 + (-pkin(3) * t77 - t32) * t89) * qJD(4) - t128 + t119, t23 * t77 + (-t77 * t139 - t74 * t89) * pkin(3) + t114, 0, t17 * t22 - t18 * t23 + (-t126 * t89 + t8 * t93 + (-t17 * t89 + t18 * t93) * qJD(4) + t118) * pkin(3), t38, t30, t51, t39, t52, 0, t103 * t88 + (-t159 - t6) * t92 + t134, t103 * t92 + (-t119 + t159) * t88 + t153, t160 * t77 * t143 + t71 * t129 + t102, t6 * t72 - t15 * t22 - g(2) * t120 - g(3) * t130 + t116 * t23 + ((-t116 * t93 + t15 * t89) * qJD(4) + t118) * pkin(3) + t104 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t101 + t149, t114 + t150, 0, 0, t38, t30, t51, t39, t52, 0, t108 * t88 + (t110 - t6) * t92 + t134, t108 * t92 + (-t110 - t119) * t88 + t153, pkin(8) * t129 - t143 * t150 + t102, -g(3) * t57 - t15 * t18 + t116 * t17 + (t119 - t6) * pkin(4) + (t104 + t60) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t135, t144 * t73, t88 * t74, t135, t147, qJDD(5), -g(1) * t92 + t78 + (t12 - t151) * qJD(5) + t105 * t88, t138 + (qJD(5) * t16 - t85) * t88 + t105 * t92, 0, 0;];
tau_reg = t1;
