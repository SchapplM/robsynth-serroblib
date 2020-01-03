% Calculate inertial parameters regressor of coriolis joint torque vector for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x(4*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:21
% DurationCPUTime: 0.78s
% Computational Cost: add. (1474->160), mult. (2710->239), div. (0->0), fcn. (1633->6), ass. (0->128)
t88 = sin(qJ(3));
t85 = t88 ^ 2;
t91 = cos(qJ(3));
t86 = t91 ^ 2;
t134 = t85 + t86;
t130 = qJD(3) * t91;
t90 = cos(qJ(4));
t144 = t90 * t91;
t156 = -qJD(4) * t144 - t90 * t130;
t87 = sin(qJ(4));
t57 = t87 * t91 + t90 * t88;
t84 = qJD(1) + qJD(2);
t49 = t57 * t84;
t132 = pkin(1) * qJD(2);
t114 = qJD(1) * t132;
t92 = cos(qJ(2));
t102 = t92 * t114;
t133 = pkin(1) * qJD(1);
t89 = sin(qJ(2));
t121 = t89 * t133;
t63 = t84 * pkin(6) + t121;
t112 = pkin(7) * t84 + t63;
t98 = qJD(3) * t112;
t26 = t91 * t102 - t88 * t98;
t42 = t112 * t88;
t37 = qJD(3) * pkin(3) - t42;
t155 = (qJD(4) * t37 + t26) * t90;
t83 = qJD(3) + qJD(4);
t154 = -pkin(7) - pkin(6);
t153 = t92 * pkin(1);
t77 = t89 * pkin(1) + pkin(6);
t152 = -pkin(7) - t77;
t147 = t87 * t88;
t125 = t84 * t147;
t47 = -t84 * t144 + t125;
t151 = t49 * t47;
t120 = t92 * t133;
t79 = -t91 * pkin(3) - pkin(2);
t50 = t79 * t84 - t120;
t150 = t50 * t49;
t149 = t84 * t88;
t43 = t112 * t91;
t148 = t87 * t43;
t146 = t89 * t91;
t145 = t90 * t43;
t93 = qJD(3) ^ 2;
t143 = t93 * t88;
t80 = t93 * t91;
t72 = t154 * t88;
t81 = t91 * pkin(7);
t73 = t91 * pkin(6) + t81;
t38 = t90 * t72 - t87 * t73;
t56 = -t144 + t147;
t111 = qJD(3) * t154;
t58 = t88 * t111;
t59 = t91 * t111;
t142 = t38 * qJD(4) + t56 * t120 + t90 * t58 + t87 * t59;
t39 = t87 * t72 + t90 * t73;
t141 = -t39 * qJD(4) + t57 * t120 - t87 * t58 + t90 * t59;
t33 = t83 * t57;
t103 = t89 * t114;
t131 = qJD(3) * t88;
t119 = t84 * t131;
t51 = pkin(3) * t119 + t103;
t140 = t50 * t33 + t51 * t56;
t99 = t83 * t147;
t32 = t99 + t156;
t139 = -t50 * t32 + t51 * t57;
t64 = -t84 * pkin(2) - t120;
t138 = t88 * t103 + t64 * t130;
t137 = t156 * t84;
t136 = t134 * t102;
t135 = t85 - t86;
t129 = qJD(3) * t92;
t128 = -qJD(1) - t84;
t127 = -qJD(2) + t84;
t126 = pkin(3) * t149;
t82 = t84 ^ 2;
t124 = t88 * t82 * t91;
t123 = t92 * t132;
t122 = pkin(3) * t131;
t118 = t84 * t130;
t117 = t88 * t129;
t113 = -pkin(3) * t83 - t37;
t17 = t90 * t37 - t148;
t18 = t87 * t37 + t145;
t27 = -t88 * t102 - t91 * t98;
t108 = -qJD(4) * t148 + t87 * t27;
t4 = t108 + t155;
t109 = -t87 * t26 + t90 * t27;
t5 = -t18 * qJD(4) + t109;
t110 = t17 * t32 - t18 * t33 - t4 * t56 - t5 * t57;
t107 = qJD(3) * t152;
t106 = t134 * qJD(2);
t104 = t88 * t118;
t101 = t127 * t133;
t100 = t128 * t132;
t54 = t152 * t88;
t55 = t91 * t77 + t81;
t28 = t90 * t54 - t87 * t55;
t29 = t87 * t54 + t90 * t55;
t97 = t50 * t47 - t108;
t96 = -t121 + t122;
t95 = t91 * t129 - t89 * t149;
t94 = -t64 * t84 - t102;
t78 = -pkin(2) - t153;
t67 = t79 - t153;
t62 = -0.2e1 * t104;
t61 = 0.2e1 * t104;
t60 = t89 * t132 + t122;
t52 = t64 * t131;
t46 = -0.2e1 * t135 * t84 * qJD(3);
t41 = t91 * t107 - t88 * t123;
t40 = t88 * t107 + t91 * t123;
t31 = t33 * t83;
t30 = t32 * t83;
t24 = t33 * t84;
t23 = t84 * t99 + t137;
t20 = -t90 * t42 - t148;
t19 = t87 * t42 - t145;
t14 = -t47 ^ 2 + t49 ^ 2;
t12 = -t137 + (-t125 + t47) * t83;
t9 = -t29 * qJD(4) - t87 * t40 + t90 * t41;
t8 = t28 * qJD(4) + t90 * t40 + t87 * t41;
t7 = t24 * t56 + t47 * t33;
t6 = -t23 * t57 - t49 * t32;
t1 = t23 * t56 - t57 * t24 + t32 * t47 - t49 * t33;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t100, t92 * t100, 0, 0, t61, t46, t80, t62, -t143, 0, t78 * t119 - t77 * t80 + t52 + (t128 * t146 - t117) * t132, t78 * t118 - t95 * t132 + t77 * t143 + t138, t84 * t106 * t153 + t136, ((qJD(1) * t78 + t64) * t89 + (qJD(1) * t77 + t63) * t92 * t134) * t132, t6, t1, -t30, t7, -t31, 0, t67 * t24 + t60 * t47 + t9 * t83 + t140, -t67 * t23 + t60 * t49 - t8 * t83 + t139, t28 * t23 - t29 * t24 - t8 * t47 - t9 * t49 + t110, t17 * t9 + t18 * t8 + t5 * t28 + t4 * t29 + t50 * t60 + t51 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t101, t92 * t101, 0, 0, t61, t46, t80, t62, -t143, 0, -pkin(2) * t119 - pkin(6) * t80 + t52 + (t127 * t146 + t117) * t133, -pkin(2) * t118 + pkin(6) * t143 + t95 * t133 + t138, -t134 * t84 * t120 + t136, ((-pkin(2) * qJD(2) - t64) * t89 + (pkin(6) * t106 - t134 * t63) * t92) * t133, t6, t1, -t30, t7, -t31, 0, t141 * t83 + t79 * t24 + t96 * t47 + t140, -t142 * t83 - t79 * t23 + t96 * t49 + t139, -t141 * t49 - t142 * t47 + t38 * t23 - t39 * t24 + t110, t141 * t17 + t142 * t18 + t5 * t38 + t4 * t39 + t96 * t50 + t51 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t124, t135 * t82, 0, t124, 0, 0, t94 * t88, t94 * t91, 0, 0, t151, t14, t12, -t151, 0, 0, -t47 * t126 - t19 * t83 - t150 + (t113 * t87 - t145) * qJD(4) + t109, -t49 * t126 + t20 * t83 + (t113 * qJD(4) - t26) * t90 + t97, (t18 + t19) * t49 + (-t17 + t20) * t47 + (t23 * t90 - t24 * t87 + (-t47 * t90 + t49 * t87) * qJD(4)) * pkin(3), -t17 * t19 - t18 * t20 + (-t50 * t149 + t4 * t87 + t5 * t90 + (-t17 * t87 + t18 * t90) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t14, t12, -t151, 0, 0, t18 * t83 - t150 + t5, t17 * t83 - t155 + t97, 0, 0;];
tauc_reg = t2;
