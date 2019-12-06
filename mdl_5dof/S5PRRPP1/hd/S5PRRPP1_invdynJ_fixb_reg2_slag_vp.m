% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:07:01
% EndTime: 2019-12-05 16:07:04
% DurationCPUTime: 1.31s
% Computational Cost: add. (1445->253), mult. (3148->293), div. (0->0), fcn. (2124->8), ass. (0->131)
t101 = sin(qJ(3));
t102 = cos(qJ(3));
t98 = sin(pkin(8));
t99 = cos(pkin(8));
t61 = t99 * t101 + t98 * t102;
t170 = t61 * qJD(2);
t141 = qJD(2) * qJD(3);
t130 = t101 * t141;
t111 = t61 * qJDD(2) - t98 * t130;
t129 = t102 * t141;
t32 = t99 * t129 + t111;
t56 = t61 * qJD(3);
t153 = t99 * t102;
t60 = t98 * t101 - t153;
t119 = -t170 * t56 - t32 * t60;
t139 = t102 * qJDD(2);
t140 = t101 * qJDD(2);
t121 = -t99 * t139 + t98 * t140;
t31 = qJD(2) * t56 + t121;
t146 = qJD(2) * t102;
t133 = t99 * t146;
t147 = qJD(2) * t101;
t54 = t98 * t147 - t133;
t145 = qJD(3) * t101;
t59 = qJD(3) * t153 - t98 * t145;
t155 = -t61 * t31 - t59 * t54;
t175 = -t119 + t155;
t174 = t119 + t155;
t92 = pkin(7) + qJ(2);
t85 = sin(t92);
t87 = cos(t92);
t171 = g(1) * t85 - g(2) * t87;
t125 = g(1) * t87 + g(2) * t85;
t167 = t170 ^ 2;
t51 = t54 ^ 2;
t173 = -t51 - t167;
t172 = -t51 + t167;
t157 = qJ(4) + pkin(6);
t131 = t157 * t101;
t72 = t157 * t102;
t39 = -t98 * t131 + t99 * t72;
t93 = qJ(3) + pkin(8);
t86 = sin(t93);
t169 = -t39 * qJDD(3) - t171 * t86;
t128 = qJD(3) * t157;
t123 = t102 * t128;
t142 = qJD(1) * qJD(3);
t90 = t102 * qJDD(1);
t20 = qJDD(3) * pkin(3) + t90 - qJD(2) * t123 + (-qJD(2) * qJD(4) - t157 * qJDD(2) - t142) * t101;
t138 = pkin(6) * t139 + t101 * qJDD(1) + t102 * t142;
t49 = t102 * qJD(4) - t101 * t128;
t24 = qJ(4) * t139 + t49 * qJD(2) + t138;
t5 = t98 * t20 + t99 * t24;
t88 = cos(t93);
t168 = g(3) * t86 + t125 * t88 - t5;
t163 = pkin(3) * t101;
t162 = g(3) * t102;
t161 = t102 * pkin(3);
t160 = t170 * t54;
t144 = t101 * qJD(1);
t50 = qJD(2) * t72 + t144;
t158 = t98 * t50;
t43 = t99 * t50;
t4 = t99 * t20 - t98 * t24;
t91 = t102 * qJD(1);
t48 = -qJD(2) * t131 + t91;
t46 = qJD(3) * pkin(3) + t48;
t23 = t98 * t46 + t43;
t95 = t101 ^ 2;
t96 = t102 ^ 2;
t154 = t95 - t96;
t152 = pkin(6) * qJDD(2);
t151 = qJDD(2) * pkin(2);
t150 = qJDD(3) * pkin(4);
t26 = t98 * t48 + t43;
t149 = t26 * qJD(3);
t28 = t99 * t48 - t158;
t148 = qJD(5) - t28;
t137 = pkin(3) * t145;
t136 = pkin(6) * t147;
t135 = pkin(6) * t146;
t104 = qJD(2) ^ 2;
t134 = t101 * t104 * t102;
t83 = pkin(2) + t161;
t126 = t101 * t129;
t122 = t88 * pkin(4) + t86 * qJ(5);
t120 = t31 * t60 + t54 * t56;
t22 = t99 * t46 - t158;
t35 = qJD(3) * t56 + qJDD(3) * t60;
t118 = -g(3) * t88 + t125 * t86 + t4;
t117 = -0.2e1 * pkin(2) * t141 - pkin(6) * qJDD(3);
t69 = -t83 * qJD(2) + qJD(4);
t115 = pkin(2) * t104 + t125;
t38 = t99 * t131 + t98 * t72;
t114 = -t38 * qJDD(3) + t171 * t88;
t47 = pkin(3) * t130 - t83 * qJDD(2) + qJDD(4);
t112 = -t101 * qJD(4) - t123;
t103 = qJD(3) ^ 2;
t110 = -pkin(6) * t103 + 0.2e1 * t151 + t171;
t21 = t54 * pkin(4) - qJ(5) * t170 + t69;
t109 = -t170 * t21 - qJDD(5) + t118;
t27 = -t99 * t112 + t98 * t49;
t29 = t98 * t112 + t99 * t49;
t108 = t170 * t27 - t29 * t54 - t39 * t31 + t38 * t32 - t125;
t107 = t31 * pkin(4) - t32 * qJ(5) + t47;
t106 = 0.2e1 * t170 * qJD(3) + t121;
t40 = -pkin(6) * t130 + t138;
t41 = -t101 * t142 + t90 + (-t129 - t140) * pkin(6);
t67 = t91 - t136;
t68 = t135 + t144;
t105 = -t41 * t101 + t40 * t102 + (-t101 * t68 - t102 * t67) * qJD(3) - t125;
t97 = qJDD(1) - g(3);
t94 = qJDD(3) * qJ(5);
t81 = -t99 * pkin(3) - pkin(4);
t78 = t98 * pkin(3) + qJ(5);
t71 = qJDD(3) * t102 - t103 * t101;
t70 = qJDD(3) * t101 + t103 * t102;
t66 = t87 * t83;
t34 = t59 * qJD(3) + t61 * qJDD(3);
t30 = t60 * pkin(4) - t61 * qJ(5) - t83;
t25 = pkin(3) * t147 + pkin(4) * t170 + t54 * qJ(5);
t17 = qJD(3) * qJ(5) + t23;
t14 = -qJD(3) * pkin(4) + qJD(5) - t22;
t11 = (t54 + t133) * qJD(3) + t111;
t10 = (-t54 + t133) * qJD(3) + t111;
t9 = t56 * pkin(4) - t59 * qJ(5) - t61 * qJD(5) + t137;
t6 = t170 * t59 + t32 * t61;
t3 = qJDD(5) - t150 - t4;
t2 = qJD(3) * qJD(5) + t5 + t94;
t1 = -qJD(5) * t170 + t107;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, t71, -t70, 0, t40 * t101 + t41 * t102 - g(3) + (-t101 * t67 + t102 * t68) * qJD(3), 0, 0, 0, 0, 0, 0, -t35, -t34, t175, -t22 * t56 + t23 * t59 - t4 * t60 + t5 * t61 - g(3), 0, 0, 0, 0, 0, 0, -t35, t175, t34, t14 * t56 + t17 * t59 + t2 * t61 + t3 * t60 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t171, t125, 0, 0, t95 * qJDD(2) + 0.2e1 * t126, 0.2e1 * t101 * t139 - 0.2e1 * t154 * t141, t70, t96 * qJDD(2) - 0.2e1 * t126, t71, 0, t117 * t101 + t110 * t102, -t110 * t101 + t117 * t102, (t95 + t96) * t152 + t105, (t171 + t151) * pkin(2) + t105 * pkin(6), t6, t174, t34, t120, -t35, 0, -t83 * t31 + t47 * t60 + t69 * t56 + (t54 * t163 - t27) * qJD(3) + t114, -t83 * t32 + t47 * t61 + t69 * t59 + (t163 * t170 - t29) * qJD(3) + t169, -t22 * t59 - t23 * t56 - t4 * t61 - t5 * t60 + t108, t5 * t39 + t23 * t29 - t4 * t38 - t22 * t27 - t47 * t83 + t69 * t137 - g(1) * (t157 * t87 - t85 * t83) - g(2) * (t157 * t85 + t66), t6, t34, -t174, 0, t35, t120, -t27 * qJD(3) + t1 * t60 + t21 * t56 + t30 * t31 + t9 * t54 + t114, t14 * t59 - t17 * t56 - t2 * t60 + t3 * t61 + t108, t29 * qJD(3) - t1 * t61 - t170 * t9 - t21 * t59 - t30 * t32 - t169, -g(2) * t66 + t1 * t30 + t14 * t27 + t17 * t29 + t2 * t39 + t21 * t9 + t3 * t38 + (-g(1) * t157 - g(2) * t122) * t87 + (-g(1) * (-t122 - t83) - g(2) * t157) * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t154 * t104, t140, t134, t139, qJDD(3), -t162 + t90 + (t68 - t135) * qJD(3) + (t115 - t142 - t152) * t101, g(3) * t101 + (t67 + t136) * qJD(3) + t115 * t102 - t138, 0, 0, t160, t172, t11, -t160, -t121, qJDD(3), t149 - t69 * t170 + (qJDD(3) * t99 - t54 * t147) * pkin(3) + t118, t28 * qJD(3) + t69 * t54 + (-qJDD(3) * t98 - t147 * t170) * pkin(3) + t168, (t23 - t26) * t170 + (-t22 + t28) * t54 + (-t31 * t98 - t32 * t99) * pkin(3), t22 * t26 - t23 * t28 + (-t162 + t4 * t99 + t5 * t98 + (-qJD(2) * t69 + t125) * t101) * pkin(3), t160, t11, -t172, qJDD(3), t121, -t160, t149 - t25 * t54 + (pkin(4) - t81) * qJDD(3) + t109, -t78 * t31 + t81 * t32 + (t17 - t26) * t170 + (t14 - t148) * t54, t78 * qJDD(3) - t21 * t54 + t25 * t170 + t94 + (0.2e1 * qJD(5) - t28) * qJD(3) - t168, t2 * t78 + t3 * t81 - t21 * t25 - t14 * t26 - g(3) * (t122 + t161) + t148 * t17 + t125 * (pkin(4) * t86 - qJ(5) * t88 + t163); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t10, t173, t170 * t22 + t23 * t54 - t171 + t47, 0, 0, 0, 0, 0, 0, t106, t173, -t10, t17 * t54 + (-qJD(5) - t14) * t170 + t107 - t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) + t160, t11, -t167 - t103, -t17 * qJD(3) - t109 - t150;];
tau_reg = t7;
