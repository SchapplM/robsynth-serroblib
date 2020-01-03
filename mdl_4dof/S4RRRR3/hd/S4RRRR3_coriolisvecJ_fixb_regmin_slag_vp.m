% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRR3
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
% tauc_reg [4x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:36
% EndTime: 2019-12-31 17:24:40
% DurationCPUTime: 0.96s
% Computational Cost: add. (1063->156), mult. (2896->232), div. (0->0), fcn. (2076->6), ass. (0->115)
t85 = sin(qJ(4));
t123 = qJD(4) * t85;
t143 = pkin(5) + pkin(6);
t90 = cos(qJ(2));
t72 = t143 * t90;
t68 = qJD(1) * t72;
t89 = cos(qJ(3));
t58 = t89 * t68;
t128 = qJD(2) * pkin(2);
t87 = sin(qJ(2));
t71 = t143 * t87;
t66 = qJD(1) * t71;
t60 = -t66 + t128;
t86 = sin(qJ(3));
t101 = -t86 * t60 - t58;
t126 = qJD(1) * t90;
t117 = t89 * t126;
t127 = qJD(1) * t87;
t118 = t86 * t127;
t51 = -t117 + t118;
t141 = t51 * pkin(7);
t18 = -t101 - t141;
t88 = cos(qJ(4));
t45 = t88 * t51;
t53 = -t86 * t126 - t89 * t127;
t31 = t85 * t53 - t45;
t80 = -t90 * pkin(2) - pkin(1);
t70 = t80 * qJD(1);
t41 = t51 * pkin(3) + t70;
t149 = t18 * t123 - t41 * t31;
t122 = qJD(1) * qJD(2);
t112 = t90 * t122;
t82 = qJD(2) + qJD(3);
t34 = qJD(3) * t117 + t89 * t112 - t82 * t118;
t113 = qJD(2) * t143;
t104 = qJD(1) * t113;
t61 = t87 * t104;
t62 = t90 * t104;
t109 = t86 * t61 - t89 * t62;
t95 = t101 * qJD(3) + t109;
t8 = -t34 * pkin(7) + t95;
t121 = -qJD(3) - qJD(4);
t81 = qJD(2) - t121;
t148 = (-t18 * t81 - t8) * t85 + t149;
t102 = t85 * t51 + t88 * t53;
t139 = t102 * t31;
t3 = t102 ^ 2 - t31 ^ 2;
t65 = t86 * t90 + t89 * t87;
t40 = t82 * t65;
t35 = t40 * qJD(1);
t4 = -qJD(4) * t45 + t53 * t123 + t88 * t34 - t85 * t35;
t1 = -t31 * t81 + t4;
t93 = t102 * qJD(4) - t85 * t34 - t88 * t35;
t2 = -t102 * t81 + t93;
t125 = qJD(3) * t86;
t108 = -t68 * t125 - t86 * t62;
t146 = (qJD(3) * t60 - t61) * t89;
t7 = -t35 * pkin(7) + t108 + t146;
t99 = t102 * t41 - t85 * t7 + t88 * t8;
t147 = -0.2e1 * t122;
t54 = t86 * t68;
t110 = t89 * t60 - t54;
t49 = t53 * pkin(7);
t17 = t110 + t49;
t145 = qJD(1) * t65;
t144 = qJD(4) - t81;
t142 = pkin(2) * t81;
t140 = t53 * pkin(3);
t137 = t53 * t51;
t136 = t70 * t53;
t135 = t86 * t88;
t134 = t88 * t18;
t92 = qJD(1) ^ 2;
t133 = t90 * t92;
t91 = qJD(2) ^ 2;
t132 = t91 * t87;
t131 = t91 * t90;
t130 = -t89 * t66 - t54;
t129 = t87 ^ 2 - t90 ^ 2;
t124 = qJD(3) * t89;
t120 = t87 * t128;
t119 = pkin(2) * t127;
t115 = -pkin(2) * t82 - t60;
t13 = t82 * pkin(3) + t17;
t114 = -pkin(3) * t81 - t13;
t107 = t86 * t66 - t58;
t105 = pkin(1) * t147;
t103 = -t85 * t13 - t134;
t64 = t86 * t87 - t89 * t90;
t37 = t88 * t64 + t85 * t65;
t38 = -t85 * t64 + t88 * t65;
t100 = t86 * t71 - t89 * t72;
t98 = t70 * t51 - t108;
t67 = t87 * t113;
t69 = t90 * t113;
t96 = -t71 * t124 - t72 * t125 - t89 * t67 - t86 * t69;
t94 = t100 * qJD(3) + t86 * t67 - t89 * t69;
t79 = t89 * pkin(2) + pkin(3);
t44 = t64 * pkin(3) + t80;
t42 = t119 - t140;
t39 = t82 * t64;
t33 = t40 * pkin(3) + t120;
t24 = -t64 * pkin(7) - t100;
t23 = -t65 * pkin(7) - t89 * t71 - t86 * t72;
t22 = t35 * pkin(3) + qJD(2) * t119;
t21 = t49 + t130;
t20 = t107 + t141;
t19 = -t51 ^ 2 + t53 ^ 2;
t16 = (-t145 - t53) * t82;
t15 = t51 * t82 + t34;
t12 = t39 * pkin(7) + t94;
t11 = -t40 * pkin(7) + t96;
t10 = t38 * qJD(4) - t85 * t39 + t88 * t40;
t9 = -t37 * qJD(4) - t88 * t39 - t85 * t40;
t5 = [0, 0, 0, 0.2e1 * t87 * t112, t129 * t147, t131, -t132, 0, -pkin(5) * t131 + t87 * t105, pkin(5) * t132 + t90 * t105, t34 * t65 + t53 * t39, -t34 * t64 - t65 * t35 + t39 * t51 + t53 * t40, -t39 * t82, -t40 * t82, 0, t80 * t35 + t70 * t40 + t94 * t82 + (qJD(1) * t64 + t51) * t120, t80 * t34 - t70 * t39 - t96 * t82 + (-t53 + t145) * t120, -t102 * t9 + t4 * t38, t10 * t102 + t31 * t9 - t4 * t37 + t38 * t93, t9 * t81, -t10 * t81, 0, -t33 * t31 - t44 * t93 + t22 * t37 + t41 * t10 + (-t85 * t11 + t88 * t12 + (-t23 * t85 - t24 * t88) * qJD(4)) * t81, -t33 * t102 + t44 * t4 + t22 * t38 + t41 * t9 - (t88 * t11 + t85 * t12 + (t23 * t88 - t24 * t85) * qJD(4)) * t81; 0, 0, 0, -t87 * t133, t129 * t92, 0, 0, 0, t92 * pkin(1) * t87, pkin(1) * t133, -t137, t19, t15, t16, 0, -t51 * t119 + t136 - t107 * t82 + (t115 * t86 - t58) * qJD(3) + t109, t53 * t119 + t130 * t82 + (t115 * qJD(3) + t61) * t89 + t98, t139, t3, t1, t2, 0, t42 * t31 - (t88 * t20 - t85 * t21) * t81 + (-t85 * t89 - t135) * qJD(3) * t142 + ((-pkin(2) * t135 - t79 * t85) * t81 + t103) * qJD(4) + t99, t42 * t102 + (-t121 * t86 * t142 + t20 * t81 - t8) * t85 + (-qJD(4) * t13 - t7 + (-pkin(2) * t124 - qJD(4) * t79 + t21) * t81) * t88 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, t19, t15, t16, 0, -t101 * t82 + t136 + t95, t110 * t82 - t146 + t98, t139, t3, t1, t2, 0, -t31 * t140 - (-t85 * t17 - t134) * t81 + (t114 * t85 - t134) * qJD(4) + t99, -t102 * t140 + (t114 * qJD(4) + t17 * t81 - t7) * t88 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t3, t1, t2, 0, t144 * t103 + t99, (-t144 * t13 - t7) * t88 + t148;];
tauc_reg = t5;
