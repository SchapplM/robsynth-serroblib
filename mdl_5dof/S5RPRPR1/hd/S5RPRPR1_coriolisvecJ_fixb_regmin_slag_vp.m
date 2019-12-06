% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:43
% EndTime: 2019-12-05 17:47:48
% DurationCPUTime: 0.87s
% Computational Cost: add. (1015->152), mult. (2283->228), div. (0->0), fcn. (1574->6), ass. (0->105)
t100 = sin(qJ(5));
t102 = cos(qJ(5));
t120 = qJD(5) * t100;
t101 = sin(qJ(3));
t103 = cos(qJ(3));
t98 = sin(pkin(8));
t99 = cos(pkin(8));
t110 = t99 * t101 + t98 * t103;
t68 = t110 * qJD(1);
t58 = t102 * t68;
t109 = t98 * t101 - t99 * t103;
t117 = qJD(1) * qJD(3);
t62 = t109 * t117;
t63 = t110 * t117;
t123 = qJD(1) * t103;
t124 = qJD(1) * t101;
t71 = t123 * t99 - t124 * t98;
t1 = -qJD(5) * t58 + t100 * t62 - t102 * t63 - t120 * t71;
t26 = t100 * t71 + t58;
t93 = qJD(3) + qJD(5);
t136 = t26 * t93;
t151 = t1 + t136;
t111 = -t100 * t68 + t102 * t71;
t150 = t111 * t26;
t137 = t111 * t93;
t2 = qJD(5) * t111 - t100 * t63 - t102 * t62;
t149 = -t2 + t137;
t148 = t111 ^ 2 - t26 ^ 2;
t140 = t68 * pkin(7);
t104 = -pkin(1) - pkin(6);
t82 = t104 * qJD(1) + qJD(2);
t66 = -qJ(4) * t124 + t101 * t82;
t135 = t99 * t66;
t67 = -qJ(4) * t123 + t103 * t82;
t55 = qJD(3) * pkin(3) + t67;
t18 = t98 * t55 + t135;
t10 = t18 - t140;
t79 = pkin(3) * t124 + qJD(1) * qJ(2) + qJD(4);
t38 = t68 * pkin(4) + t79;
t147 = t10 * t120 + t38 * t26;
t145 = qJD(5) - t93;
t118 = t103 * qJD(4);
t122 = qJD(3) * t101;
t39 = -t82 * t122 + (qJ(4) * t122 - t118) * qJD(1);
t119 = t101 * qJD(4);
t121 = qJD(3) * t103;
t40 = t82 * t121 + (-qJ(4) * t121 - t119) * qJD(1);
t13 = t99 * t39 - t98 * t40;
t4 = t63 * pkin(7) + t13;
t14 = t98 * t39 + t99 * t40;
t5 = t62 * pkin(7) + t14;
t144 = -t100 * t5 + t102 * t4 - t38 * t111;
t94 = qJD(1) * qJD(2);
t143 = 0.2e1 * t94;
t142 = pkin(3) * t98;
t30 = -t100 * t109 + t102 * t110;
t69 = -t121 * t99 + t122 * t98;
t70 = t110 * qJD(3);
t6 = -qJD(5) * t30 + t100 * t69 - t102 * t70;
t141 = t6 * t93;
t31 = -t100 * t110 - t102 * t109;
t7 = qJD(5) * t31 - t100 * t70 - t102 * t69;
t139 = t7 * t93;
t138 = t71 * pkin(7);
t51 = t98 * t66;
t127 = qJ(4) - t104;
t64 = t122 * t127 - t118;
t81 = t127 * t103;
t65 = -qJD(3) * t81 - t119;
t20 = t98 * t64 + t99 * t65;
t24 = t99 * t67 - t51;
t80 = t127 * t101;
t37 = -t99 * t80 - t98 * t81;
t114 = t103 * t117;
t134 = pkin(3) * t114 + t94;
t133 = t101 ^ 2 - t103 ^ 2;
t132 = t101 * pkin(3) + qJ(2);
t105 = qJD(3) ^ 2;
t131 = t105 * t101;
t130 = t105 * t103;
t106 = qJD(1) ^ 2;
t129 = t106 * qJ(2);
t128 = t106 * t103;
t126 = pkin(3) * t121 + qJD(2);
t125 = -t105 - t106;
t116 = 0.2e1 * qJD(1);
t17 = t99 * t55 - t51;
t19 = t99 * t64 - t98 * t65;
t23 = -t98 * t67 - t135;
t36 = t98 * t80 - t99 * t81;
t9 = qJD(3) * pkin(4) - t138 + t17;
t112 = -t102 * t10 - t100 * t9;
t108 = -t109 * t13 + t110 * t14 - t17 * t70 - t18 * t69;
t88 = t99 * pkin(3) + pkin(4);
t56 = pkin(4) * t110 + t132;
t44 = pkin(3) * t123 + t71 * pkin(4);
t41 = -t69 * pkin(4) + t126;
t32 = -t62 * pkin(4) + t134;
t22 = -pkin(7) * t110 + t37;
t21 = pkin(7) * t109 + t36;
t16 = t24 - t138;
t15 = t23 + t140;
t12 = t69 * pkin(7) + t20;
t11 = t70 * pkin(7) + t19;
t3 = [0, 0, 0, 0, t143, qJ(2) * t143, -0.2e1 * t101 * t114, 0.2e1 * t133 * t117, -t131, -t130, 0, -t104 * t131 + (qJ(2) * t121 + qJD(2) * t101) * t116, -t104 * t130 + (-qJ(2) * t122 + qJD(2) * t103) * t116, -t19 * t71 - t20 * t68 + t36 * t63 + t37 * t62 - t108, t79 * t126 + t13 * t36 + t134 * t132 + t14 * t37 + t17 * t19 + t18 * t20, t1 * t31 + t111 * t6, -t1 * t30 - t111 * t7 - t31 * t2 - t6 * t26, t141, -t139, 0, t41 * t26 + t56 * t2 + t32 * t30 + t38 * t7 + (-t100 * t12 + t102 * t11 + (-t100 * t21 - t102 * t22) * qJD(5)) * t93, t41 * t111 + t56 * t1 + t32 * t31 + t38 * t6 - (t100 * t11 + t102 * t12 + (-t100 * t22 + t102 * t21) * qJD(5)) * t93; 0, 0, 0, 0, -t106, -t129, 0, 0, 0, 0, 0, t125 * t101, t125 * t103, -t109 * t63 + t110 * t62 + t69 * t68 + t70 * t71, -t79 * qJD(1) + t108, 0, 0, 0, 0, 0, -qJD(1) * t26 + t141, -qJD(1) * t111 - t139; 0, 0, 0, 0, 0, 0, t101 * t128, -t133 * t106, 0, 0, 0, -qJ(2) * t128, t101 * t129, (t18 + t23) * t71 - (t17 - t24) * t68 + (t62 * t98 + t63 * t99) * pkin(3), -t17 * t23 - t18 * t24 + (-t123 * t79 + t13 * t99 + t14 * t98) * pkin(3), t150, t148, t151, t149, 0, -t44 * t26 - (-t100 * t16 + t102 * t15) * t93 + ((-t100 * t88 - t102 * t142) * t93 + t112) * qJD(5) + t144, -t102 * t5 - t100 * t4 - t44 * t111 + (t100 * t15 + t102 * t16) * t93 + (-(-t100 * t142 + t102 * t88) * t93 - t102 * t9) * qJD(5) + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68 ^ 2 - t71 ^ 2, t17 * t71 + t18 * t68 + t134, 0, 0, 0, 0, 0, t2 + t137, t1 - t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t148, t151, t149, 0, t145 * t112 + t144, (-t10 * t93 - t4) * t100 + (-t145 * t9 - t5) * t102 + t147;];
tauc_reg = t3;
