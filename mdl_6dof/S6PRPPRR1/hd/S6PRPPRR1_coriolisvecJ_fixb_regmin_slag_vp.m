% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:14
% EndTime: 2019-03-08 19:16:17
% DurationCPUTime: 0.98s
% Computational Cost: add. (1532->200), mult. (4108->291), div. (0->0), fcn. (3419->12), ass. (0->126)
t86 = sin(pkin(12));
t128 = qJD(2) * t86;
t93 = sin(qJ(5));
t119 = t93 * t128;
t89 = cos(pkin(12));
t127 = qJD(2) * t89;
t96 = cos(qJ(5));
t78 = t96 * t127;
t60 = -t78 + t119;
t59 = qJD(6) + t60;
t88 = sin(pkin(6));
t129 = qJD(1) * t88;
t97 = cos(qJ(2));
t120 = t97 * t129;
t90 = cos(pkin(11));
t112 = t90 * t120;
t94 = sin(qJ(2));
t121 = t94 * t129;
t87 = sin(pkin(11));
t75 = t87 * t121;
t53 = -t75 + t112;
t151 = t53 - qJD(4);
t70 = qJD(2) * t112;
t40 = t70 + (qJD(4) - t75) * qJD(2);
t68 = t93 * t86 - t96 * t89;
t154 = t68 * t40;
t69 = t96 * t86 + t93 * t89;
t153 = t69 * t40;
t54 = (t87 * t94 - t90 * t97) * t88;
t62 = t69 * qJD(2);
t152 = -qJD(6) + t59;
t130 = t86 ^ 2 + t89 ^ 2;
t150 = qJD(2) * t130;
t71 = qJD(2) * pkin(2) + t120;
t76 = t90 * t121;
t42 = t87 * t71 + t76;
t39 = qJD(2) * qJ(4) + t42;
t91 = cos(pkin(6));
t79 = t91 * qJD(1) + qJD(3);
t73 = t89 * t79;
t21 = t73 + (-pkin(8) * qJD(2) - t39) * t86;
t26 = t89 * t39 + t86 * t79;
t22 = pkin(8) * t127 + t26;
t10 = t93 * t21 + t96 * t22;
t4 = t10 * qJD(5) + t153;
t149 = (t62 * pkin(5) + t59 * pkin(9)) * t59 + t4;
t41 = t90 * t71 - t75;
t111 = qJD(4) - t41;
t118 = -t89 * pkin(4) - pkin(3);
t32 = t118 * qJD(2) + t111;
t11 = t60 * pkin(5) - t62 * pkin(9) + t32;
t80 = t87 * pkin(2) + qJ(4);
t145 = pkin(8) + t80;
t65 = t145 * t86;
t66 = t145 * t89;
t105 = -t96 * t65 - t93 * t66;
t133 = -t105 * qJD(5) - t151 * t68;
t146 = t90 * pkin(2);
t74 = t118 - t146;
t28 = t68 * pkin(5) - t69 * pkin(9) + t74;
t9 = t96 * t21 - t93 * t22;
t3 = t9 * qJD(5) - t154;
t30 = -t93 * t65 + t96 * t66;
t64 = t69 * qJD(5);
t57 = qJD(2) * t64;
t63 = t68 * qJD(5);
t7 = -qJD(5) * pkin(5) - t9;
t148 = -(qJD(6) * t11 + t3) * t68 + t4 * t69 - t7 * t63 + (-qJD(6) * t28 + t133) * t59 - t30 * t57;
t92 = sin(qJ(6));
t125 = qJD(6) * t92;
t95 = cos(qJ(6));
t102 = t69 * t125 + t95 * t63;
t49 = t95 * t57;
t147 = -t102 * t59 + t69 * t49;
t47 = t92 * qJD(5) + t95 * t62;
t77 = qJD(5) * t78;
t56 = -qJD(5) * t119 + t77;
t24 = t47 * qJD(6) + t92 * t56;
t123 = t95 * qJD(5);
t23 = qJD(6) * t123 - t62 * t125 + t95 * t56;
t144 = t23 * t92;
t143 = t28 * t57;
t55 = (t87 * t97 + t90 * t94) * t88;
t51 = qJD(2) * t55;
t43 = qJD(1) * t51;
t142 = t43 * t54;
t45 = t92 * t62 - t123;
t141 = t45 * t59;
t140 = t47 * t59;
t139 = t47 * t62;
t138 = t62 * t45;
t98 = qJD(2) ^ 2;
t137 = t88 * t98;
t135 = t92 * t57;
t134 = t23 * t68 + t47 * t64;
t132 = t30 * qJD(5) - t151 * t69;
t50 = t87 * t120 + t76;
t131 = t64 * pkin(5) + t63 * pkin(9) - t50;
t126 = qJD(6) * t69;
t124 = t63 * qJD(5);
t117 = t130 * t40;
t114 = t59 * t95;
t113 = qJD(2) * t50 - t43;
t8 = qJD(5) * pkin(9) + t10;
t1 = t95 * t11 - t92 * t8;
t2 = t92 * t11 + t95 * t8;
t37 = -t55 * t86 + t91 * t89;
t38 = t55 * t89 + t91 * t86;
t14 = t93 * t37 + t96 * t38;
t110 = t95 * t14 + t54 * t92;
t109 = -t92 * t14 + t54 * t95;
t108 = -t68 * t24 - t64 * t45;
t107 = (-t86 * t39 + t73) * t86 - t26 * t89;
t106 = t96 * t37 - t93 * t38;
t103 = t49 + (-t60 * t92 - t125) * t59;
t101 = -pkin(9) * t57 + (t7 + t9) * t59;
t100 = (-t95 * t126 + t92 * t63) * t59 - t69 * t135;
t58 = t64 * qJD(5);
t52 = qJD(2) * t54;
t44 = -qJD(2) * t75 + t70;
t36 = -qJD(2) * pkin(3) + t111;
t15 = t57 * pkin(5) - t56 * pkin(9) + t43;
t12 = t95 * t15;
t6 = t14 * qJD(5) - t69 * t52;
t5 = t106 * qJD(5) + t68 * t52;
t13 = [0, 0, -t94 * t137, -t97 * t137, -t41 * t51 - t42 * t52 + t44 * t55 + t142, -t51 * t127, t51 * t128, -t52 * t150, t36 * t51 + t142 + t107 * t52 + (-t37 * t86 + t38 * t89) * t40, 0, 0, 0, 0, 0, -t6 * qJD(5) + t51 * t60 + t54 * t57, -t5 * qJD(5) + t51 * t62 + t54 * t56, 0, 0, 0, 0, 0 (-qJD(6) * t110 - t92 * t5 + t51 * t95) * t59 + t109 * t57 + t6 * t45 - t106 * t24 -(qJD(6) * t109 + t95 * t5 + t51 * t92) * t59 - t110 * t57 + t6 * t47 - t106 * t23; 0, 0, 0, 0, t41 * t50 - t42 * t53 + (-t43 * t90 + t44 * t87) * pkin(2), t113 * t89, -t113 * t86, -t151 * t150 + t117, t43 * (-pkin(3) - t146) - t36 * t50 + t80 * t117 + t151 * t107, t56 * t69 - t62 * t63, -t56 * t68 - t69 * t57 + t63 * t60 - t62 * t64, -t124, -t58, 0, -t132 * qJD(5) + t32 * t64 + t43 * t68 - t50 * t60 + t74 * t57, qJD(5) * t133 - t32 * t63 + t43 * t69 - t50 * t62 + t74 * t56, t23 * t95 * t69 - t102 * t47 -(-t45 * t95 - t47 * t92) * t63 + (-t144 - t24 * t95 + (t45 * t92 - t47 * t95) * qJD(6)) * t69, t134 + t147, t100 + t108, t57 * t68 + t59 * t64, t1 * t64 + t12 * t68 - t105 * t24 + t132 * t45 + (t143 + t131 * t59 + (-t30 * t59 - t8 * t68 + t7 * t69) * qJD(6)) * t95 + t148 * t92, -t2 * t64 - t105 * t23 + t132 * t47 + (-t143 - (-qJD(6) * t8 + t15) * t68 - t7 * t126 + (qJD(6) * t30 - t131) * t59) * t92 + t148 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, t124, 0, 0, 0, 0, 0, t100 - t108, t134 - t147; 0, 0, 0, 0, 0, 0, 0, -t130 * t98 (qJD(1) * t55 + t107) * qJD(2), 0, 0, 0, 0, 0, 0.2e1 * t62 * qJD(5), t77 + (-t60 - t119) * qJD(5), 0, 0, 0, 0, 0, t103 - t138, -t59 ^ 2 * t95 - t135 - t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t60, -t60 ^ 2 + t62 ^ 2, t77 + (t60 - t119) * qJD(5), 0, 0, -t32 * t62 - t153, t32 * t60 + t154, t114 * t47 + t144 (t23 - t141) * t95 + (-t24 - t140) * t92, t114 * t59 + t135 - t139, t103 + t138, -t59 * t62, -pkin(5) * t24 - t1 * t62 - t10 * t45 + t101 * t92 - t149 * t95, -pkin(5) * t23 - t10 * t47 + t101 * t95 + t149 * t92 + t2 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t45, -t45 ^ 2 + t47 ^ 2, t23 + t141, t140 - t24, t57, t152 * t2 - t92 * t3 - t7 * t47 + t12, t152 * t1 - t92 * t15 - t95 * t3 + t7 * t45;];
tauc_reg  = t13;
