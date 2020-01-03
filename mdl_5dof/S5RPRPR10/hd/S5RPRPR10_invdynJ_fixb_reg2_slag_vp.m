% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR10
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR10_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:09
% EndTime: 2019-12-31 18:26:11
% DurationCPUTime: 1.10s
% Computational Cost: add. (2285->209), mult. (3337->248), div. (0->0), fcn. (1813->10), ass. (0->129)
t134 = qJD(1) - qJD(3);
t92 = sin(pkin(8));
t93 = cos(pkin(8));
t95 = sin(qJ(3));
t98 = cos(qJ(3));
t48 = t92 * t95 - t93 * t98;
t174 = t134 * t48;
t181 = t134 * t174;
t180 = t134 ^ 2;
t138 = qJ(2) * qJD(1);
t100 = -pkin(1) - pkin(2);
t67 = t100 * qJD(1) + qJD(2);
t38 = -t95 * t138 + t98 * t67;
t179 = t134 * t38;
t178 = -qJD(3) * t138 + t100 * qJDD(1) + qJDD(2);
t30 = -pkin(3) * t134 + t38;
t39 = t98 * t138 + t67 * t95;
t32 = t93 * t39;
t15 = t92 * t30 + t32;
t10 = -pkin(7) * t134 + t15;
t94 = sin(qJ(5));
t97 = cos(qJ(5));
t7 = qJD(4) * t97 - t10 * t94;
t143 = t7 * qJD(5);
t88 = qJDD(1) - qJDD(3);
t165 = t88 * pkin(3);
t136 = qJD(1) * qJD(2);
t137 = qJ(2) * qJDD(1);
t176 = qJD(3) * t67 + t136 + t137;
t21 = -t176 * t95 + t178 * t98;
t13 = t21 - t165;
t20 = t176 * t98 + t178 * t95;
t6 = t92 * t13 + t93 * t20;
t4 = -pkin(7) * t88 + t6;
t1 = t94 * qJDD(4) + t97 * t4 + t143;
t164 = t10 * t97;
t8 = qJD(4) * t94 + t164;
t81 = t97 * qJDD(4);
t2 = -t8 * qJD(5) - t94 * t4 + t81;
t106 = -(t7 * t97 + t8 * t94) * qJD(5) + t1 * t97 - t2 * t94;
t135 = qJ(3) + pkin(8);
t128 = sin(t135);
t129 = cos(t135);
t96 = sin(qJ(1));
t99 = cos(qJ(1));
t40 = -t96 * t128 - t99 * t129;
t41 = t99 * t128 - t96 * t129;
t123 = g(1) * t40 + g(2) * t41;
t177 = t106 + t123;
t175 = -t39 * t134 + t21;
t149 = t99 * pkin(1) + t96 * qJ(2);
t79 = pkin(3) * t98 + pkin(2);
t173 = t99 * t79 + t149;
t58 = -qJ(2) * t95 + t98 * t100;
t101 = qJD(5) ^ 2;
t49 = t92 * t98 + t93 * t95;
t152 = t134 * t49;
t115 = -t134 * t152 + t48 * t88;
t172 = -t101 * t49 + t115;
t171 = t41 * pkin(4) + t40 * pkin(7);
t154 = t99 * t98;
t157 = t96 * t95;
t51 = -t154 - t157;
t156 = t96 * t98;
t158 = t95 * t99;
t52 = -t156 + t158;
t170 = g(1) * t52 - g(2) * t51;
t36 = t98 * qJD(2) + t58 * qJD(3);
t59 = qJ(2) * t98 + t100 * t95;
t37 = -t95 * qJD(2) - t59 * qJD(3);
t19 = t36 * t93 + t37 * t92;
t163 = t19 * t134;
t161 = t39 * t92;
t23 = t38 * t93 - t161;
t162 = t23 * t134;
t159 = t94 * t97;
t155 = t97 * t88;
t153 = -t93 * t13 + t92 * t20;
t53 = -pkin(3) + t58;
t27 = t92 * t53 + t93 * t59;
t70 = pkin(3) * t158;
t150 = pkin(3) * t156 - t70;
t148 = g(1) * t96 - g(2) * t99;
t90 = t94 ^ 2;
t91 = t97 ^ 2;
t147 = t90 - t91;
t146 = t90 + t91;
t142 = pkin(1) * qJDD(1);
t140 = qJD(5) * t134;
t139 = qJDD(4) + g(3);
t133 = t180 * t159;
t132 = 0.2e1 * t136;
t131 = t146 * t88;
t126 = qJDD(2) - t142;
t125 = -0.2e1 * t140 * t159;
t124 = g(1) * t41 - g(2) * t40;
t122 = g(1) * t99 + g(2) * t96;
t119 = t7 * t94 - t8 * t97;
t14 = t30 * t93 - t161;
t26 = t53 * t93 - t59 * t92;
t68 = pkin(3) * t157;
t118 = pkin(4) * t40 - pkin(7) * t41 - t68;
t83 = t99 * qJ(2);
t117 = t70 + t83 + (-pkin(1) - t79) * t96;
t114 = t123 + t6;
t18 = t36 * t92 - t93 * t37;
t113 = t134 * t18 - t124;
t22 = t38 * t92 + t32;
t112 = -t134 * t22 + t124;
t24 = pkin(4) - t26;
t25 = -pkin(7) + t27;
t9 = pkin(4) * t134 - t14;
t110 = -qJDD(5) * t25 + (-t134 * t24 - t19 - t9) * qJD(5);
t72 = pkin(3) * t92 + pkin(7);
t73 = -pkin(3) * t93 - pkin(4);
t109 = -qJDD(5) * t72 + (-t134 * t73 + t23 + t9) * qJD(5);
t108 = -0.2e1 * t174 * qJD(5) - qJDD(5) * t49;
t107 = -qJD(4) * qJD(5) + t134 * t9 - t123 - t4;
t3 = pkin(4) * t88 + t153;
t105 = t101 * t25 - t24 * t88 - t113 - t3;
t104 = -t101 * t72 + t73 * t88 + t112 - t3;
t103 = -g(1) * t51 - g(2) * t52 - t20;
t102 = qJD(1) ^ 2;
t61 = qJDD(5) * t97 - t101 * t94;
t60 = qJDD(5) * t94 + t101 * t97;
t35 = t88 * t91 + t125;
t34 = -t88 * t90 + t125;
t28 = t147 * t140 - t94 * t155;
t5 = [0, 0, 0, 0, 0, qJDD(1), t148, t122, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + 0.2e1 * t142 + t148, 0, -t122 + t132 + 0.2e1 * t137, -t126 * pkin(1) - g(1) * (-pkin(1) * t96 + t83) - g(2) * t149 + (t132 + t137) * qJ(2), 0, 0, 0, 0, 0, t88, -t134 * t37 - t58 * t88 - t170 - t21, t134 * t36 + t59 * t88 - t103, 0, t20 * t59 + t39 * t36 + t21 * t58 + t38 * t37 - g(1) * (t100 * t96 + t83) - g(2) * (pkin(2) * t99 + t149), 0, 0, 0, 0, 0, t88, -t26 * t88 + t113 + t153, t27 * t88 + t114 + t163, 0, t6 * t27 + t15 * t19 - t153 * t26 - t14 * t18 - g(1) * t117 - g(2) * (t68 + t173), -t34, -0.2e1 * t28, -t60, t35, -t61, 0, -t105 * t97 + t110 * t94, t105 * t94 + t110 * t97, -t131 * t25 - t146 * t163 - t177, t3 * t24 + t9 * t18 - g(1) * (t117 + t171) - g(2) * (-t118 + t173) - t119 * t19 + t106 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t102, -qJ(2) * t102 + t126 - t148, 0, 0, 0, 0, 0, 0, -t180 * t95 - t98 * t88, -t180 * t98 + t95 * t88, 0, t175 * t98 + (t20 + t179) * t95 - t148, 0, 0, 0, 0, 0, 0, t115, t49 * t88 + t181, 0, t152 * t14 + t15 * t174 + t153 * t48 + t6 * t49 - t148, 0, 0, 0, 0, 0, 0, t108 * t94 + t172 * t97, t108 * t97 - t172 * t94, -t49 * t131 - t146 * t181, t106 * t49 - t174 * t119 - t152 * t9 + t3 * t48 - t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, t170 + t175, t103 - t179, 0, 0, 0, 0, 0, 0, 0, -t88, -t93 * t165 + t112 - t153, t92 * t165 - t114 - t162, 0, -t15 * t23 + t14 * t22 - g(1) * t150 + g(2) * t68 + (g(2) * t154 - t153 * t93 + t6 * t92) * pkin(3), t34, 0.2e1 * t28, t60, -t35, t61, 0, t104 * t97 + t109 * t94, -t104 * t94 + t109 * t97, -t131 * t72 + t146 * t162 + t177, t3 * t73 - t9 * t22 - g(1) * (t150 - t171) - g(2) * (-pkin(3) * t154 + t118) + t119 * t23 + t106 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, 0, 0, 0, 0, 0, t61, -t60, 0, -qJD(5) * t119 + t1 * t94 + t2 * t97 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, t147 * t180, -t94 * t88, t133, -t155, qJDD(5), g(3) * t97 + t81 + (t8 - t164) * qJD(5) + t107 * t94, t143 + (qJD(5) * t10 - t139) * t94 + t107 * t97, 0, 0;];
tau_reg = t5;
