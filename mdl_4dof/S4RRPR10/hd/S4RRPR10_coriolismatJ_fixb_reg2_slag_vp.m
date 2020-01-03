% Calculate inertial parameters regressor of coriolis matrix for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRPR10_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR10_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:58
% EndTime: 2019-12-31 17:12:00
% DurationCPUTime: 1.08s
% Computational Cost: add. (970->149), mult. (1946->218), div. (0->0), fcn. (1606->4), ass. (0->130)
t165 = pkin(3) + pkin(5);
t90 = sin(qJ(2));
t64 = t165 * t90;
t91 = cos(qJ(4));
t157 = t91 * t64;
t151 = t90 * qJ(3);
t92 = cos(qJ(2));
t93 = -pkin(2) - pkin(6);
t104 = -t92 * t93 + t151;
t46 = -pkin(1) - t104;
t89 = sin(qJ(4));
t30 = t89 * t46 - t157;
t160 = t89 * t64;
t31 = t91 * t46 + t160;
t11 = (t30 * t89 + t31 * t91) * t90;
t127 = t92 * qJD(1);
t121 = t91 * t127;
t85 = t89 ^ 2;
t87 = t91 ^ 2;
t68 = t85 - t87;
t96 = t68 * qJD(2) + 0.2e1 * t89 * t121;
t126 = t92 * qJD(3);
t166 = qJD(2) * t104 - t126;
t150 = t92 * qJ(3);
t60 = t90 * pkin(2) - t150;
t48 = t90 * pkin(6) + t60;
t161 = t89 * t48;
t65 = t165 * t92;
t56 = t65 * t91;
t32 = t56 - t161;
t164 = t32 * t91;
t158 = t91 * t48;
t162 = t65 * t89;
t33 = t158 + t162;
t163 = t33 * t89;
t159 = t89 * t92;
t156 = t91 * t92;
t86 = t90 ^ 2;
t88 = t92 ^ 2;
t69 = t88 - t86;
t1 = t33 * t156 - t32 * t159 - t11;
t155 = t1 * qJD(1);
t2 = -t30 * t32 + t31 * t33 - t65 * t64;
t154 = t2 * qJD(1);
t4 = (-t30 - t157) * t92 + (t32 - t56) * t90;
t153 = t4 * qJD(1);
t5 = (t31 - t160) * t92 + (t33 - t162) * t90;
t152 = t5 * qJD(1);
t149 = qJD(1) * t90;
t148 = qJD(1) * t91;
t147 = qJD(3) * t90;
t146 = qJD(3) * t91;
t145 = qJD(4) * t90;
t144 = qJD(4) * t91;
t143 = qJD(4) * t93;
t142 = t11 * qJD(1);
t17 = t65 * t156 - t30 * t90;
t141 = t17 * qJD(1);
t18 = -t65 * t159 - t31 * t90;
t140 = t18 * qJD(1);
t107 = -t92 * pkin(2) - t151;
t57 = -pkin(1) + t107;
t34 = t57 * t92 + t60 * t90;
t139 = t34 * qJD(1);
t35 = -t57 * t90 + t60 * t92;
t138 = t35 * qJD(1);
t42 = (t85 + t87) * t92 * t90;
t137 = t42 * qJD(1);
t47 = (t87 / 0.2e1 - t85 / 0.2e1) * t92;
t136 = t47 * qJD(4);
t52 = t69 * t89;
t135 = t52 * qJD(1);
t54 = t69 * t91;
t134 = t54 * qJD(1);
t133 = t69 * qJD(1);
t132 = t86 * qJD(1);
t131 = t86 * qJD(3);
t130 = t89 * qJD(2);
t129 = t90 * qJD(2);
t128 = t91 * qJD(2);
t81 = t92 * qJD(2);
t125 = qJ(3) * qJD(4);
t83 = qJD(2) * qJ(3);
t124 = pkin(1) * t149;
t123 = pkin(1) * t127;
t122 = pkin(5) * t129;
t120 = t89 * t145;
t119 = t90 * t144;
t118 = t57 * t60 * qJD(1);
t117 = t57 * t149;
t116 = t89 * t81;
t74 = t90 * t81;
t73 = t90 * t127;
t115 = t89 * t144;
t114 = t89 * t128;
t102 = -t163 / 0.2e1 - t164 / 0.2e1;
t108 = (pkin(5) / 0.2e1 + pkin(3) / 0.2e1) * t92;
t7 = t108 + t102;
t113 = t7 * qJD(1) + t83;
t112 = t88 * t115;
t110 = t92 * t114;
t109 = qJD(4) + t149;
t105 = t163 + t164;
t103 = -t132 - t145;
t101 = -t93 * t90 / 0.2e1 - t150 / 0.2e1;
t100 = t109 * t159;
t38 = t47 * qJD(1) + t114;
t95 = t48 / 0.2e1 + t101;
t26 = t95 * t89;
t99 = -t26 * qJD(1) - t91 * t83;
t27 = t95 * t91;
t98 = -t27 * qJD(1) + t89 * t83;
t36 = t89 * t88 * t148 - t47 * qJD(2);
t53 = t68 * t88;
t97 = -t53 * qJD(1) + 0.2e1 * t110;
t94 = t107 * qJD(2) + t126;
t84 = qJ(3) * qJD(3);
t80 = pkin(5) * t81;
t77 = t81 / 0.2e1;
t76 = t91 * t81;
t75 = t90 * t148;
t72 = t89 * t149;
t61 = t69 * qJD(2);
t51 = -t75 - t144;
t50 = -qJD(4) * t89 - t72;
t49 = t73 + t92 * qJD(4) / 0.2e1;
t20 = -t162 - t158 / 0.2e1 + t101 * t91;
t19 = t56 - t161 / 0.2e1 + t101 * t89;
t6 = t108 - t102;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t61, 0, -t74, 0, 0, -pkin(1) * t129, -pkin(1) * t81, 0, 0, 0, 0, 0, t74, t61, -t74, 0, t35 * qJD(2) - t90 * t126, -t34 * qJD(2) + t131, (qJD(2) * t60 - t147) * t57, -t85 * t74 + t112, -t53 * qJD(4) - 0.2e1 * t110 * t90, -t52 * qJD(2) - t119 * t92, -t74 * t87 - t112, -t54 * qJD(2) + t120 * t92, t74, t4 * qJD(2) + t18 * qJD(4) + t131 * t89, -t5 * qJD(2) - t17 * qJD(4) + t131 * t91, -t1 * qJD(2) + t42 * qJD(3), t2 * qJD(2) - t11 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t133, t81, -t73, -t129, 0, -t80 - t124, t122 - t123, 0, 0, 0, -t81, t129, t73, t133, -t73, t94, t80 + t138, -t122 - t139, t94 * pkin(5) + t118, -t136 + (-t85 * t127 + t114) * t90, 0.2e1 * t92 * t115 - t96 * t90, t76 - t135, t136 + (-t127 * t87 - t114) * t90, -t116 - t134, t49, t19 * qJD(4) - t64 * t130 - t166 * t91 + t153, t20 * qJD(4) - t64 * t128 + t166 * t89 - t152, -qJD(2) * t105 - t155, t154 + (-t64 * qJ(3) + t105 * t93) * qJD(2) + t6 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t73, t132, t80 - t117, 0, 0, 0, 0, 0, 0, t132 * t89 + t76, t132 * t91 - t116, t137, t6 * qJD(2) - t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t97, -t109 * t156, -t36, t100, t77, t19 * qJD(2) - t31 * qJD(4) + t140, t20 * qJD(2) + t30 * qJD(4) - t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, -t133, 0, t73, 0, 0, t124, t123, 0, 0, 0, 0, 0, -t73, -t133, t73, 0, -t138, t139, -t118, t85 * t73 - t136, 0.2e1 * t91 * t100, -t120 + t135, t73 * t87 + t136, -t119 + t134, -t49, t26 * qJD(4) - t153, t27 * qJD(4) + t152, t155, t7 * qJD(3) - t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t84, -t115, t68 * qJD(4), 0, t115, 0, 0, qJD(3) * t89 + t125 * t91, -t125 * t89 + t146, 0, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t83, 0, 0, 0, 0, 0, 0, t130, t128, 0, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t96, t50, t38, t51, -t127 / 0.2e1, -t143 * t89 - t99, -t143 * t91 - t98, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t132, t117, 0, 0, 0, 0, 0, 0, t103 * t89, t103 * t91, -t137, -t7 * qJD(2) + t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t83, 0, 0, 0, 0, 0, 0, -t130, -t128, 0, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t97, (t121 + t130) * t90, t36, (-t127 * t89 + t128) * t90, t77, -t26 * qJD(2) + t89 * t147 - t140, -t27 * qJD(2) + t90 * t146 + t141, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t96, t72, -t38, t75, t127 / 0.2e1, t99, t98, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t3;
