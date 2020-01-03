% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:29
% DurationCPUTime: 1.07s
% Computational Cost: add. (1718->182), mult. (3271->214), div. (0->0), fcn. (1683->6), ass. (0->126)
t75 = cos(qJ(3));
t121 = qJD(3) * t75;
t114 = pkin(2) * t121;
t123 = pkin(1) * qJD(1);
t76 = cos(qJ(2));
t112 = t76 * t123;
t73 = sin(qJ(2));
t115 = t73 * t123;
t72 = sin(qJ(3));
t57 = t72 * t115;
t43 = t112 * t75 - t57;
t152 = t114 - t43;
t71 = sin(qJ(4));
t69 = t71 ^ 2;
t74 = cos(qJ(4));
t70 = t74 ^ 2;
t124 = t69 + t70;
t127 = (qJD(2) + qJD(3)) * t57;
t68 = qJD(1) + qJD(2);
t51 = pkin(2) * t68 + t112;
t99 = qJD(2) * t112;
t14 = (qJD(3) * t51 + t99) * t75 - t127;
t156 = t124 * t14;
t65 = qJD(3) + t68;
t106 = t124 * t65;
t155 = 2 * qJD(4);
t11 = t74 * t14;
t33 = t115 * t75 + t51 * t72;
t28 = pkin(8) * t65 + t33;
t137 = t71 * t28;
t2 = t11 + (qJD(5) - t137) * qJD(4);
t119 = qJD(4) * t74;
t9 = t71 * t14;
t4 = t119 * t28 + t9;
t154 = t2 * t74 + t4 * t71;
t153 = t124 * t28;
t135 = t73 * t75;
t91 = t72 * t76 + t135;
t151 = qJD(2) * t91 + t73 * t121;
t125 = -t69 + t70;
t150 = t125 * t65 * t155;
t77 = qJD(4) ^ 2;
t149 = pkin(8) * t77;
t148 = t65 * pkin(3);
t147 = t75 * pkin(2);
t122 = qJD(3) * t72;
t136 = t72 * t73;
t62 = pkin(1) * t76 + pkin(2);
t21 = t62 * t121 + (-t73 * t122 + (t75 * t76 - t136) * qJD(2)) * pkin(1);
t145 = t21 * t65;
t80 = t151 * pkin(1);
t22 = t122 * t62 + t80;
t144 = t22 * t65;
t143 = t33 * t65;
t126 = pkin(1) * t135 + t62 * t72;
t41 = pkin(8) + t126;
t142 = t41 * t77;
t53 = -pkin(4) * t74 - qJ(5) * t71 - pkin(3);
t141 = t53 * t65;
t60 = pkin(2) * t72 + pkin(8);
t140 = t60 * t77;
t139 = t65 * t71;
t138 = t65 * t74;
t134 = t74 * t28;
t111 = t51 * t122;
t15 = qJD(1) * t80 + t111;
t32 = t51 * t75 - t57;
t27 = -t32 - t148;
t133 = t119 * t27 + t15 * t71;
t132 = t124 * t145;
t120 = qJD(4) * t71;
t131 = t120 * t32 + t138 * t33;
t42 = t91 * t123;
t130 = t120 * t43 + t138 * t42;
t116 = qJD(4) * qJ(5);
t117 = t71 * qJD(5);
t45 = pkin(4) * t120 - t116 * t74 - t117;
t129 = t33 - t45;
t113 = pkin(2) * t122;
t34 = t45 + t113;
t128 = t34 - t42;
t20 = t116 + t134;
t118 = t20 * qJD(4);
t109 = t65 * t120;
t94 = pkin(4) * t71 - qJ(5) * t74;
t5 = (qJD(4) * t94 - t117) * t65 + t15;
t108 = -t5 - t149;
t107 = -t5 - t140;
t103 = -pkin(1) * t136 + t62 * t75;
t30 = -t103 + t53;
t104 = t30 * t65 - t21;
t101 = -qJD(4) * pkin(4) + qJD(5);
t100 = t74 * t109;
t98 = t32 * t106;
t97 = -t42 + t113;
t96 = (-qJD(2) + t68) * t123;
t95 = pkin(1) * qJD(2) * (-qJD(1) - t68);
t19 = t101 + t137;
t93 = t19 * t71 + t20 * t74;
t92 = t142 + t144;
t90 = (-pkin(2) * t65 - t51) * qJD(3);
t40 = -pkin(3) - t103;
t89 = qJD(4) * (t40 * t65 - t21);
t7 = t22 + t45;
t88 = -t65 * t7 - t142 - t5;
t46 = t53 - t147;
t87 = t46 * t65 - t114;
t61 = -pkin(3) - t147;
t86 = t61 * t65 - t114;
t85 = -t118 * t71 + t119 * t19 + t154;
t84 = t152 * t106;
t81 = (t19 * t74 - t20 * t71) * qJD(4) + t154;
t79 = t151 * t123;
t78 = -t79 - t111;
t67 = t77 * t74;
t66 = t77 * t71;
t64 = t65 ^ 2;
t52 = t71 * t64 * t74;
t50 = -0.2e1 * t100;
t49 = 0.2e1 * t100;
t44 = t125 * t64;
t39 = t94 * t65;
t24 = t27 * t120;
t13 = -t32 + t141;
t6 = t13 * t120;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t95, t76 * t95, 0, 0, 0, 0, 0, 0, 0, 0, t78 - t144, -t14 - t145, 0, -t103 * t15 + t126 * t14 + t21 * t33 - t22 * t32, t49, t150, t67, t50, -t66, 0, t24 + t71 * t89 + (-t15 - t92) * t74, t71 * t92 + t74 * t89 + t133, t132 + t156, t15 * t40 + t153 * t21 + t156 * t41 + t27 * t22, t49, t67, -t150, 0, t66, t50, t104 * t120 + t74 * t88 + t6, t85 + t132, t88 * t71 + (-t104 - t13) * t119, t13 * t7 + t21 * t93 + t5 * t30 + t41 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t96, t76 * t96, 0, 0, 0, 0, 0, 0, 0, 0, t42 * t65 + t72 * t90 - t79, t43 * t65 + (t90 - t99) * t75 + t127, 0, t32 * t42 - t33 * t43 + (t14 * t72 - t15 * t75 + (-t32 * t72 + t33 * t75) * qJD(3)) * pkin(2), t49, t150, t67, t50, -t66, 0, t24 + t86 * t120 + (-t113 * t65 - t140 - t15) * t74 + t130, (t65 * t97 + t140) * t71 + (t43 + t86) * t119 + t133, t84 + t156, t15 * t61 + t152 * t153 + t156 * t60 + t97 * t27, t49, t67, -t150, 0, t66, t50, t6 + t87 * t120 + (-t34 * t65 + t107) * t74 + t130, t84 + t85, (-t128 * t65 + t107) * t71 + (-t13 - t43 - t87) * t119, t128 * t13 + t152 * t93 + t5 * t46 + t81 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 + t143, t32 * t65 - t14, 0, 0, t49, t150, t67, t50, -t66, 0, -pkin(3) * t109 + t24 + (-t15 - t149) * t74 + t131, (-t143 + t149) * t71 + (t32 - t148) * t119 + t133, -t98 + t156, -t15 * pkin(3) + pkin(8) * t156 - t153 * t32 - t27 * t33, t49, t67, -t150, 0, t66, t50, t53 * t109 + t6 + (-t45 * t65 + t108) * t74 + t131, -t98 + t85, (-t13 - t32 - t141) * t119 + (t129 * t65 + t108) * t71, pkin(8) * t81 - t129 * t13 - t32 * t93 + t5 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t44, 0, t52, 0, 0, -t139 * t27 - t9, -t138 * t27 - t11, 0, 0, -t52, 0, t44, 0, 0, t52, -t9 + (-t13 * t71 + t39 * t74) * t65, ((t20 - t116) * t71 + (t101 - t19) * t74) * t65, qJD(5) * t155 + t11 + (t13 * t74 + t39 * t71) * t65, -t19 * t134 - t4 * pkin(4) + t2 * qJ(5) - t13 * t39 + (qJD(5) + t137) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, 0, -t64 * t69 - t77, t13 * t139 - t118 + t4;];
tauc_reg = t1;
