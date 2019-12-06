% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:32
% DurationCPUTime: 1.19s
% Computational Cost: add. (1177->221), mult. (3068->334), div. (0->0), fcn. (2136->8), ass. (0->128)
t79 = sin(qJ(4));
t120 = t79 * qJD(3);
t77 = sin(pkin(5));
t129 = qJD(1) * t77;
t83 = cos(qJ(3));
t84 = cos(qJ(2));
t143 = t83 * t84;
t80 = sin(qJ(3));
t95 = pkin(3) * t80 - pkin(8) * t83;
t60 = t95 * qJD(3);
t81 = sin(qJ(2));
t82 = cos(qJ(4));
t165 = t80 * pkin(7) * t120 + t82 * t60 - (-t79 * t143 + t81 * t82) * t129;
t122 = qJD(4) * t82;
t64 = -t83 * pkin(3) - t80 * pkin(8) - pkin(2);
t164 = -(t82 * t143 + t79 * t81) * t129 + t64 * t122 + t79 * t60;
t128 = qJD(1) * t83;
t108 = t81 * t129;
t62 = qJD(2) * pkin(7) + t108;
t78 = cos(pkin(5));
t37 = t78 * t128 - t80 * t62;
t115 = qJD(3) * qJD(4);
t33 = (t83 * t120 + t80 * t122) * qJD(2) + t79 * t115;
t126 = qJD(2) * t80;
t57 = t82 * t126 + t120;
t163 = t57 ^ 2;
t117 = t83 * qJD(2);
t71 = -qJD(4) + t117;
t148 = t78 * t80;
t70 = qJD(1) * t148;
t38 = t83 * t62 + t70;
t30 = qJD(3) * pkin(8) + t38;
t107 = t84 * t129;
t39 = t64 * qJD(2) - t107;
t12 = -t79 * t30 + t82 * t39;
t8 = -t57 * qJ(5) + t12;
t3 = -t71 * pkin(4) + t8;
t162 = t3 - t8;
t118 = t82 * qJD(5);
t132 = qJ(5) * t80;
t144 = t82 * t83;
t72 = pkin(7) * t144;
t131 = qJ(5) * t82;
t92 = pkin(4) * t80 - t83 * t131;
t161 = -t80 * t118 + t92 * qJD(3) + (-t72 + (-t64 + t132) * t79) * qJD(4) + t165;
t145 = t80 * t82;
t160 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t145 + (-qJD(5) * t80 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t83) * t79 + t164;
t124 = qJD(3) * t83;
t127 = qJD(2) * t77;
t112 = t84 * t127;
t97 = t80 * t112;
t18 = qJD(1) * t97 + qJD(3) * t70 + t62 * t124;
t159 = t18 * t79;
t158 = t18 * t82;
t29 = -qJD(3) * pkin(3) - t37;
t157 = t29 * t79;
t119 = t82 * qJD(3);
t105 = t83 * t119;
t106 = t79 * t126;
t32 = -qJD(2) * t105 + qJD(4) * t106 - t82 * t115;
t156 = t32 * t79;
t55 = t106 - t119;
t155 = t55 * t71;
t154 = t57 * t71;
t153 = t71 * t79;
t152 = t71 * t82;
t151 = t77 * t81;
t150 = t77 * t84;
t86 = qJD(2) ^ 2;
t149 = t77 * t86;
t147 = t79 * t39;
t146 = t79 * t83;
t85 = qJD(3) ^ 2;
t142 = t85 * t80;
t141 = t85 * t83;
t140 = -qJ(5) - pkin(8);
t59 = t95 * qJD(2);
t139 = t82 * t37 + t79 * t59;
t98 = qJD(4) * t140;
t138 = t118 - t139 + (qJ(5) * t117 + t98) * t79;
t99 = -t79 * t37 + t82 * t59;
t137 = -t92 * qJD(2) - t79 * qJD(5) + t82 * t98 - t99;
t134 = t79 * t64 + t72;
t75 = t80 ^ 2;
t133 = -t83 ^ 2 + t75;
t130 = qJD(2) * pkin(2);
t125 = qJD(3) * t80;
t123 = qJD(4) * t79;
t121 = t29 * qJD(4);
t116 = qJD(2) * qJD(3);
t114 = t81 * t149;
t113 = t81 * t127;
t111 = t71 * t123;
t110 = t80 * t123;
t109 = t71 * t122;
t103 = pkin(4) * t79 + pkin(7);
t101 = t80 * t116;
t17 = -t62 * t125 + (qJD(3) * t78 + t112) * t128;
t36 = (t60 + t108) * qJD(2);
t100 = t79 * t17 - t82 * t36;
t96 = t83 * t112;
t63 = -t107 - t130;
t94 = -t63 - t107;
t13 = t82 * t30 + t147;
t10 = t33 * pkin(4) + t18;
t93 = qJD(2) * t75 - t71 * t83;
t44 = t83 * t151 + t148;
t22 = -t82 * t150 - t44 * t79;
t91 = t79 * t150 - t44 * t82;
t43 = t80 * t151 - t78 * t83;
t90 = t39 * t122 - t30 * t123 + t82 * t17 + t79 * t36;
t88 = qJD(3) * (-t94 - t130);
t87 = -t13 * qJD(4) - t100;
t66 = t140 * t82;
t65 = t140 * t79;
t54 = t82 * t64;
t52 = t55 ^ 2;
t24 = -t79 * t132 + t134;
t21 = t44 * qJD(3) + t97;
t20 = -t43 * qJD(3) + t96;
t19 = -t80 * t131 + t54 + (-pkin(7) * t79 - pkin(4)) * t83;
t15 = t55 * pkin(4) + qJD(5) + t29;
t9 = -t55 * qJ(5) + t13;
t6 = t22 * qJD(4) + t79 * t113 + t20 * t82;
t5 = t91 * qJD(4) + t82 * t113 - t20 * t79;
t2 = -t33 * qJ(5) - t55 * qJD(5) + t90;
t1 = pkin(4) * t101 + t32 * qJ(5) - t57 * qJD(5) + t87;
t4 = [0, 0, -t114, -t84 * t149, 0, 0, 0, 0, 0, -t83 * t114 + (-t21 - t97) * qJD(3), t80 * t114 + (-t20 - t96) * qJD(3), 0, 0, 0, 0, 0, t101 * t22 + t21 * t55 + t43 * t33 - t5 * t71, t101 * t91 + t21 * t57 - t43 * t32 + t6 * t71, t22 * t32 + t33 * t91 - t5 * t57 - t6 * t55, t1 * t22 + t10 * t43 + t15 * t21 - t2 * t91 + t3 * t5 + t9 * t6; 0, 0, 0, 0, 0.2e1 * t83 * t101, -0.2e1 * t133 * t116, t141, -t142, 0, -pkin(7) * t141 + t80 * t88, pkin(7) * t142 + t83 * t88, -t32 * t145 + (t105 - t110) * t57, (-t55 * t82 - t57 * t79) * t124 + (t156 - t33 * t82 + (t55 * t79 - t57 * t82) * qJD(4)) * t80, t71 * t110 + t32 * t83 + (t57 * t80 + t93 * t82) * qJD(3), t80 * t109 + t33 * t83 + (-t55 * t80 - t93 * t79) * qJD(3), (-t71 - t117) * t125, (t64 * t123 - t165) * t71 + ((pkin(7) * t55 + t157) * qJD(3) + (t147 + (pkin(7) * t71 + t30) * t82) * qJD(4) + t100) * t83 + (-t55 * t107 + t82 * t121 + pkin(7) * t33 + t159 + ((-pkin(7) * t146 + t54) * qJD(2) + t12) * qJD(3)) * t80, t164 * t71 + (t29 * t119 + (qJD(3) * t57 - t111) * pkin(7) + t90) * t83 + (-t57 * t107 - t79 * t121 - pkin(7) * t32 + t158 + (-pkin(7) * t152 - t134 * qJD(2) - t13) * qJD(3)) * t80, t19 * t32 - t24 * t33 - t161 * t57 - t160 * t55 + (-t3 * t82 - t79 * t9) * t124 + (-t1 * t82 - t2 * t79 + (t3 * t79 - t82 * t9) * qJD(4)) * t80, t1 * t19 + t2 * t24 + t160 * t9 + t161 * t3 + t15 * t103 * t124 + (t10 * t103 + (pkin(4) * t122 - t107) * t15) * t80; 0, 0, 0, 0, -t80 * t86 * t83, t133 * t86, 0, 0, 0, t38 * qJD(3) - t63 * t126 - t18, t94 * t117, -t57 * t152 - t156, (-t32 + t155) * t82 + (-t33 + t154) * t79, -t109 + (t71 * t144 + (-t57 + t120) * t80) * qJD(2), t111 + (-t71 * t146 + (t55 + t119) * t80) * qJD(2), t71 * t126, -pkin(3) * t33 - t158 + t99 * t71 - t38 * t55 + (pkin(8) * t152 + t157) * qJD(4) + (-t12 * t80 + (-pkin(8) * t125 - t29 * t83) * t79) * qJD(2), pkin(3) * t32 + t159 - t139 * t71 - t38 * t57 + (-pkin(8) * t153 + t29 * t82) * qJD(4) + (-t29 * t144 + (-pkin(8) * t119 + t13) * t80) * qJD(2), t65 * t32 + t66 * t33 - t137 * t57 - t138 * t55 + (t3 * t71 + t2) * t82 + (t71 * t9 - t1) * t79, -t2 * t66 + t1 * t65 + t10 * (-t82 * pkin(4) - pkin(3)) + t138 * t9 + t137 * t3 + (-pkin(4) * t153 - t38) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t55, -t52 + t163, -t32 - t155, -t154 - t33, t101, -t13 * t71 - t29 * t57 + t87, -t12 * t71 + t29 * t55 - t90, pkin(4) * t32 - t162 * t55, t162 * t9 + (-t15 * t57 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 - t163, t3 * t57 + t9 * t55 + t10;];
tauc_reg = t4;
