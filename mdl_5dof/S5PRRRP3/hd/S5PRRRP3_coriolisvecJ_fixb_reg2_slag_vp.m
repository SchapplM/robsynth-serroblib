% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:44:13
% EndTime: 2019-12-05 16:44:17
% DurationCPUTime: 1.04s
% Computational Cost: add. (1589->195), mult. (4070->250), div. (0->0), fcn. (2706->4), ass. (0->124)
t100 = cos(qJ(3));
t155 = cos(qJ(4));
t98 = sin(qJ(4));
t99 = sin(qJ(3));
t72 = t100 * t98 + t155 * t99;
t138 = qJD(2) * t72;
t139 = t138 * qJ(5);
t133 = t99 * qJD(1);
t156 = -pkin(7) - pkin(6);
t82 = t156 * t100;
t59 = -qJD(2) * t82 + t133;
t53 = t98 * t59;
t143 = qJD(3) * pkin(3);
t123 = qJD(2) * t156;
t94 = t100 * qJD(1);
t58 = t123 * t99 + t94;
t56 = t58 + t143;
t24 = t155 * t56 - t53;
t12 = t24 - t139;
t129 = t99 * t143;
t159 = 0.2e1 * t129;
t95 = qJD(3) + qJD(4);
t131 = qJD(2) * qJD(3);
t158 = -0.2e1 * t131;
t81 = t156 * t99;
t41 = -t155 * t82 + t98 * t81;
t157 = t138 ^ 2;
t149 = t98 * t99;
t111 = t95 * t149;
t118 = t155 * qJD(4);
t121 = t155 * t100;
t37 = -qJD(3) * t121 - t100 * t118 + t111;
t154 = t37 * t95;
t112 = qJD(2) * t121;
t137 = qJD(2) * t99;
t125 = t98 * t137;
t61 = -t112 + t125;
t153 = t61 * t95;
t152 = t138 * t61;
t93 = -pkin(3) * t100 - pkin(2);
t80 = qJD(2) * t93;
t150 = t80 * t138;
t11 = pkin(4) * t95 + t12;
t148 = t11 - t12;
t38 = t95 * t72;
t33 = t38 * qJD(2);
t147 = -t33 * t72 + t37 * t61;
t114 = pkin(3) * t118;
t146 = -pkin(3) * t33 * t98 - t114 * t61;
t27 = t155 * t58 - t53;
t145 = t95 * t112;
t144 = -t100 ^ 2 + t99 ^ 2;
t101 = qJD(3) ^ 2;
t142 = t101 * t99;
t32 = qJD(2) * t111 - t145;
t141 = t32 * qJ(5);
t140 = t61 * qJ(5);
t136 = qJD(4) * t98;
t102 = qJD(2) ^ 2;
t135 = t100 * t102;
t134 = t101 * t100;
t117 = t61 * pkin(4) + qJD(5);
t39 = t117 + t80;
t132 = qJD(5) + t39;
t130 = pkin(3) * t136;
t128 = pkin(3) * t137;
t127 = pkin(6) * t137;
t126 = t155 * pkin(3);
t124 = t99 * t135;
t55 = t155 * t59;
t120 = t99 * t131;
t23 = pkin(3) * t120 + pkin(4) * t33;
t122 = qJD(3) * t156;
t76 = t99 * t122;
t90 = qJD(3) * t94;
t50 = qJD(2) * t76 + t90;
t51 = (t100 * t123 - t133) * qJD(3);
t119 = t155 * t51 - t98 * t50;
t26 = -t58 * t98 - t55;
t40 = t155 * t81 + t82 * t98;
t116 = -t118 * t56 + t136 * t59 - t155 * t50 - t51 * t98;
t115 = pkin(2) * t158;
t113 = t100 * t120;
t71 = -t121 + t149;
t110 = t138 * t38 - t32 * t71;
t109 = t61 * t80 + t116;
t25 = t56 * t98 + t55;
t108 = t33 * qJ(5) + t116;
t79 = pkin(6) * qJD(2) * t100 + t133;
t77 = t100 * t122;
t14 = t118 * t81 + t136 * t82 + t155 * t76 + t77 * t98;
t107 = -t125 * t95 + t145;
t106 = t132 * t61 + t108;
t10 = -qJD(4) * t25 + t119;
t15 = -qJD(4) * t41 + t155 * t77 - t98 * t76;
t68 = -pkin(6) * t120 + t90;
t69 = t79 * qJD(3);
t78 = t94 - t127;
t105 = t68 * t100 + t69 * t99 + (-t100 * t78 - t79 * t99) * qJD(3);
t104 = t10 + t141;
t103 = (-t55 + (-pkin(3) * t95 - t56) * t98) * qJD(4) + t119;
t92 = t126 + pkin(4);
t83 = t95 * t114;
t60 = t61 ^ 2;
t52 = pkin(4) * t71 + t93;
t46 = pkin(4) * t138 + t128;
t36 = t38 * t95;
t31 = pkin(4) * t38 + t129;
t29 = -qJ(5) * t71 + t41;
t28 = -qJ(5) * t72 + t40;
t21 = -t60 + t157;
t18 = t107 + t153;
t17 = -t139 + t27;
t16 = t26 + t140;
t13 = t25 - t140;
t8 = t33 * t71 + t38 * t61;
t7 = -t138 * t37 - t32 * t72;
t6 = t37 * qJ(5) - t72 * qJD(5) + t15;
t5 = -qJ(5) * t38 - qJD(5) * t71 + t14;
t4 = -qJD(5) * t138 + t104;
t3 = -qJD(5) * t61 - t108;
t2 = t110 + t147;
t1 = -t110 + t147;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, -t134, 0, -t69 * t100 + t68 * t99 + (t100 * t79 - t78 * t99) * qJD(3), 0, 0, 0, 0, 0, 0, -t36, t154, t2, -t10 * t71 - t116 * t72 - t24 * t38 - t25 * t37, 0, 0, 0, 0, 0, 0, -t36, t154, t2, -t11 * t38 - t13 * t37 + t3 * t72 - t4 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t113, t144 * t158, t134, -0.2e1 * t113, -t142, 0, -pkin(6) * t134 + t115 * t99, pkin(6) * t142 + t100 * t115, t105, t105 * pkin(6), t7, t1, -t154, t8, -t36, 0, t15 * t95 + t93 * t33 + t80 * t38 + (qJD(2) * t71 + t61) * t129, t138 * t159 - t14 * t95 - t93 * t32 - t80 * t37, -t10 * t72 + t116 * t71 - t138 * t15 - t14 * t61 + t24 * t37 - t25 * t38 + t32 * t40 - t33 * t41, t10 * t40 - t116 * t41 + t25 * t14 + t24 * t15 + t159 * t80, t7, t1, -t154, t8, -t36, 0, t23 * t71 + t31 * t61 + t33 * t52 + t38 * t39 + t6 * t95, t138 * t31 + t23 * t72 - t32 * t52 - t37 * t39 - t5 * t95, t11 * t37 - t13 * t38 - t138 * t6 + t28 * t32 - t29 * t33 - t3 * t71 - t4 * t72 - t5 * t61, t11 * t6 + t13 * t5 + t23 * t52 + t28 * t4 + t29 * t3 + t31 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t144 * t102, 0, t124, 0, 0, t102 * pkin(2) * t99, pkin(2) * t135 - t90 + (t78 + t127) * qJD(3), 0, 0, t152, t21, t18, -t152, 0, 0, -t128 * t61 - t26 * t95 + t103 - t150, -t128 * t138 + t27 * t95 + t109 - t83, t32 * t126 + (-t24 + t27) * t61 + (t25 + t26 + t130) * t138 + t146, -t24 * t26 - t25 * t27 + (-t80 * t137 + t155 * t10 - t116 * t98 + (t155 * t25 - t24 * t98) * qJD(4)) * pkin(3), t152, t21, t18, -t152, 0, 0, -t132 * t138 - t16 * t95 - t46 * t61 + t103 + t141, -t138 * t46 + t17 * t95 + t106 - t83, t92 * t32 + (-t11 + t17) * t61 + (t13 + t16 + t130) * t138 + t146, -t11 * t16 - t13 * t17 - t39 * t46 + t4 * t92 + (t3 * t98 + (-t11 * t98 + t13 * t155) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, t21, t18, -t152, 0, 0, t25 * t95 + t10 - t150, t24 * t95 + t109, 0, 0, t152, t21, t18, -t152, 0, 0, t13 * t95 + (-t117 - t39) * t138 + t104, -pkin(4) * t157 + t12 * t95 + t106, t32 * pkin(4) - t148 * t61, t148 * t13 + (-t138 * t39 + t4) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95 * t138 + t33, t107 - t153, -t60 - t157, t11 * t138 + t13 * t61 + t23;];
tauc_reg = t9;
