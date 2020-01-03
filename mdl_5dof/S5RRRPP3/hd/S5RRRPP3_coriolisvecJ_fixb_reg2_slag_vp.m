% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPP3_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:41
% EndTime: 2019-12-31 20:53:44
% DurationCPUTime: 0.94s
% Computational Cost: add. (1080->213), mult. (1932->242), div. (0->0), fcn. (862->4), ass. (0->136)
t92 = sin(qJ(3));
t89 = t92 ^ 2;
t94 = cos(qJ(3));
t90 = t94 ^ 2;
t139 = t89 + t90;
t87 = qJD(1) + qJD(2);
t148 = t87 * t94;
t95 = cos(qJ(2));
t160 = t139 * t95;
t127 = qJD(3) * qJ(4);
t111 = -qJD(5) - t127;
t138 = pkin(1) * qJD(1);
t93 = sin(qJ(2));
t124 = t93 * t138;
t55 = t87 * pkin(7) + t124;
t42 = t94 * t55;
t27 = pkin(4) * t148 + t42;
t17 = -t111 + t27;
t134 = qJD(3) * t92;
t137 = pkin(1) * qJD(2);
t119 = qJD(1) * t137;
t112 = t95 * t119;
t63 = t94 * t112;
t88 = qJD(3) * qJD(4);
t141 = t63 + t88;
t15 = t55 * t134 - t141;
t133 = qJD(3) * t94;
t61 = t92 * t112;
t20 = t55 * t133 + t61;
t159 = -t15 * t94 + t20 * t92;
t91 = -pkin(3) - qJ(5);
t158 = qJD(3) * t91;
t41 = t92 * t55;
t157 = -qJD(4) - t41;
t156 = 0.2e1 * (-t89 + t90) * t87 * qJD(3);
t96 = qJD(3) ^ 2;
t155 = pkin(7) * t96;
t154 = t87 * pkin(2);
t153 = t95 * pkin(1);
t135 = t92 * qJ(4);
t40 = t91 * t94 - pkin(2) - t135;
t151 = t40 * t87;
t108 = -t94 * pkin(3) - t135;
t57 = -pkin(2) + t108;
t150 = t57 * t87;
t149 = t87 * t92;
t86 = t87 ^ 2;
t147 = t90 * t86;
t81 = t96 * t92;
t82 = t96 * t94;
t123 = t95 * t138;
t56 = -t123 - t154;
t73 = t93 * t119;
t146 = t56 * t133 + t92 * t73;
t125 = t95 * t137;
t145 = t139 * t87 * t125;
t114 = t87 * t124;
t132 = qJD(3) * t95;
t120 = t92 * t132;
t144 = t94 * t114 + t120 * t138;
t143 = t139 * t112;
t142 = t63 + 0.2e1 * t88;
t122 = t87 * t134;
t140 = pkin(3) * t122 + t73;
t136 = qJ(4) * t94;
t33 = -t42 - t127;
t131 = t33 * qJD(3);
t130 = t92 * qJD(4);
t129 = -qJD(1) - t87;
t26 = -pkin(4) * t149 - t41;
t128 = qJD(4) - t26;
t121 = t87 * t133;
t126 = pkin(4) * t121 + t20;
t80 = t93 * t137;
t104 = qJ(5) * t92 - t136;
t1 = (t104 * qJD(3) - t94 * qJD(5) - t130) * t87 + t140;
t78 = pkin(3) * t134;
t18 = qJ(5) * t134 + t111 * t94 - t130 + t78;
t13 = t80 + t18;
t118 = -t13 * t87 - t1;
t117 = -t18 * t87 - t1;
t19 = -t123 + t150;
t116 = -t19 - t150;
t113 = t92 * t121;
t110 = (-pkin(4) * t87 - t55) * t92;
t99 = -t94 * t127 - t130;
t35 = t78 + t99;
t109 = t35 * t87 + t155;
t28 = t35 + t80;
t76 = t93 * pkin(1) + pkin(7);
t107 = t28 * t87 + t76 * t96;
t29 = -qJD(3) * pkin(3) - t157;
t106 = t29 * t94 + t33 * t92;
t105 = t29 * t92 - t33 * t94;
t103 = t92 * t131 + t29 * t133 + t159;
t102 = t139 * t123;
t14 = t128 + t158;
t6 = qJD(3) * t110 + t141;
t7 = -qJD(3) * qJD(5) + t126;
t101 = t14 * t133 - t17 * t134 + t6 * t94 + t7 * t92;
t100 = t87 * t102;
t43 = t57 - t153;
t98 = qJD(3) * (-t43 * t87 + t125 - t19);
t97 = t106 * qJD(3) + t159;
t85 = t94 * pkin(4);
t84 = t92 * pkin(4);
t79 = pkin(4) * t133;
t77 = -pkin(2) - t153;
t74 = t89 * t86;
t72 = pkin(3) * t149;
t70 = qJD(4) * t148;
t67 = t94 * pkin(7) + t85;
t66 = t92 * pkin(7) + t84;
t65 = t92 * t86 * t94;
t58 = -t74 - t96;
t54 = -0.2e1 * t113;
t53 = 0.2e1 * t113;
t50 = t92 * t114;
t48 = pkin(7) * t133 + t79;
t47 = (-pkin(4) - pkin(7)) * t134;
t46 = t94 * t76 + t85;
t45 = t92 * t76 + t84;
t44 = -t74 + t147;
t38 = t56 * t134;
t36 = -t87 * t136 + t72;
t32 = t40 - t153;
t23 = t104 * t87 + t72;
t22 = t92 * t125 + t76 * t133 + t79;
t21 = t94 * t125 + (-pkin(4) - t76) * t134;
t12 = t19 * t149;
t11 = t99 * t87 + t140;
t10 = -t123 + t151;
t8 = t11 * t94;
t5 = t10 * t134;
t4 = t10 * t148;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87 * t80 - t73, t129 * t125, 0, 0, t53, t156, t82, t54, -t81, 0, t77 * t122 - t76 * t82 + t38 + (t129 * t94 * t93 - t120) * t137, t77 * t121 + t76 * t81 + (-t94 * t132 + t93 * t149) * t137 + t146, t143 + t145, ((qJD(1) * t77 + t56) * t93 + (qJD(1) * t76 + t55) * t160) * t137, 0, -t82, t81, t53, t156, t54, t103 + t145, t107 * t94 + t92 * t98 + t8, (-t107 - t11) * t92 + t94 * t98, t105 * t125 + t11 * t43 + t19 * t28 + t97 * t76, 0, t81, t82, t54, -t156, t53, (t21 * t94 + t22 * t92 + (t45 * t94 - t46 * t92) * qJD(3)) * t87 + t101, t118 * t92 + (t21 + (-t32 * t87 - t10) * t94) * qJD(3), t5 + t118 * t94 + (t32 * t149 - t22) * qJD(3), t1 * t32 + t10 * t13 + t14 * t22 + t17 * t21 + t7 * t45 + t6 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73 + t114, (-qJD(2) + t87) * t123, 0, 0, t53, t156, t82, t54, -t81, 0, -pkin(2) * t122 + t38 + (-t73 - t155) * t94 + t144, pkin(7) * t81 - t50 + (t123 - t154) * t133 + t146, -t100 + t143, ((-pkin(2) * qJD(2) - t56) * t93 + (pkin(7) * qJD(2) - t55) * t160) * t138, 0, -t82, t81, t53, t156, t54, -t100 + t103, t109 * t94 + t116 * t134 - t144 + t8, t50 + (-t109 - t11) * t92 + (t116 - t123) * t133, t11 * t57 + t19 * t35 + (-t105 * t95 - t19 * t93) * t138 + t97 * pkin(7), 0, t81, t82, t54, -t156, t53, (t47 * t94 + t48 * t92 + (t66 * t94 - t67 * t92) * qJD(3) - t102) * t87 + t101, t50 + t117 * t92 + (t47 + (-t10 - t123 - t151) * t94) * qJD(3), t5 + t117 * t94 + (t40 * t149 - t48) * qJD(3) + t144, t1 * t40 + t10 * t18 + t14 * t48 + t17 * t47 + t6 * t67 + t7 * t66 + (-t10 * t93 + (-t14 * t92 - t17 * t94) * t95) * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t44, 0, t65, 0, 0, -t56 * t149 - t61, -t56 * t148 - t63, 0, 0, 0, 0, 0, -t65, -t44, t65, t70 + (t108 * qJD(3) - t106) * t87, -t36 * t148 + t12 + t61, (t19 * t94 + t36 * t92) * t87 + t142, -t20 * pkin(3) - t15 * qJ(4) + t157 * t33 - t19 * t36 - t29 * t42, 0, 0, 0, t65, t44, -t65, t70 + (-t14 - t26 + t158) * t148, t23 * t149 + t4 + (-t26 + t110) * qJD(3) + t142, (-t10 * t92 + t23 * t94) * t87 + (0.2e1 * qJD(5) + t27) * qJD(3) - t126, t6 * qJ(4) - t10 * t23 + t7 * t91 + t128 * t17 + (-qJD(5) - t27) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t58, t12 + t131 + t20, 0, 0, 0, 0, 0, 0, 0, t58, -t65, t10 * t149 + (-qJD(5) - t17) * qJD(3) + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t96 - t147, t4 + (t14 + t110) * qJD(3) + t141;];
tauc_reg = t2;
