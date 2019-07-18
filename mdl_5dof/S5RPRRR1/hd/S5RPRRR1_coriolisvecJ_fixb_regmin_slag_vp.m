% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% tauc_reg [5x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:02
% EndTime: 2019-07-18 13:26:07
% DurationCPUTime: 1.32s
% Computational Cost: add. (860->214), mult. (2456->334), div. (0->0), fcn. (1816->6), ass. (0->127)
t48 = sin(qJ(5));
t114 = qJD(5) * t48;
t50 = sin(qJ(3));
t133 = t48 * t50;
t49 = sin(qJ(4));
t112 = t49 * qJD(2);
t52 = cos(qJ(4));
t108 = qJ(2) * qJD(1);
t53 = cos(qJ(3));
t86 = t53 * t108;
t25 = t52 * t86 + t112;
t51 = cos(qJ(5));
t113 = qJD(5) * t51;
t118 = qJD(3) * t53;
t57 = t50 * t113 + t48 * t118;
t44 = t52 * qJD(3);
t105 = qJD(2) * qJD(4);
t83 = t52 * t105;
t115 = qJD(4) * t53;
t96 = t49 * t115;
t110 = t52 * qJD(2);
t99 = t53 * t110;
t7 = t83 + (t99 + (-t50 * t44 - t96) * qJ(2)) * qJD(1);
t1 = (t57 * qJ(2) + qJD(2) * t133) * qJD(1) - t25 * t114 + t51 * t7;
t111 = t49 * qJD(3);
t121 = qJD(1) * t50;
t28 = t52 * t121 + t111;
t109 = t53 * qJD(1);
t42 = -qJD(4) + t109;
t11 = t48 * t28 + t51 * t42;
t153 = t11 * t42;
t92 = t49 * t121;
t26 = -t44 + t92;
t22 = qJD(5) + t26;
t116 = qJD(4) * t52;
t126 = t53 * t51;
t63 = t52 * t126 + t133;
t19 = t63 * qJD(1);
t69 = t51 * t116 - t19;
t56 = -t49 * t114 + t69;
t152 = t56 * t22;
t104 = qJD(3) * qJD(4);
t94 = t50 * t116;
t58 = t53 * t111 + t94;
t17 = qJD(1) * t58 + t49 * t104;
t13 = t51 * t28 - t48 * t42;
t107 = qJD(1) * qJD(3);
t81 = t53 * t107;
t16 = -qJD(4) * t92 + (t104 + t81) * t52;
t82 = t50 * t107;
t4 = qJD(5) * t13 + t48 * t16 - t51 * t82;
t3 = -t42 * t113 - t28 * t114 + t51 * t16 + t48 * t82;
t150 = t3 * t48;
t149 = t49 * t3;
t148 = t52 * t4;
t132 = t49 * t50;
t79 = qJ(2) * t107;
t68 = t49 * t105 - t79 * t132;
t8 = (qJ(2) * t116 + t112) * t109 + t68;
t147 = t8 * t48;
t146 = t8 * t51;
t145 = t11 * t22;
t144 = t13 * t22;
t143 = t16 * t49;
t142 = t17 * t52;
t77 = t22 * t51;
t141 = t26 * t42;
t140 = t26 * t50;
t139 = t28 * t42;
t138 = t28 * t50;
t137 = t28 * t52;
t136 = t42 * t49;
t135 = t42 * t53;
t134 = t48 * t17;
t131 = t49 * t53;
t130 = t50 * t17;
t129 = t50 * t52;
t128 = t51 * t17;
t127 = t53 * t48;
t46 = t50 ^ 2;
t125 = -t53 ^ 2 + t46;
t124 = qJ(2) * t50;
t54 = qJD(3) ^ 2;
t123 = qJ(2) * t54;
t122 = qJD(1) * t46;
t120 = qJD(3) * t50;
t119 = qJD(3) * t51;
t117 = qJD(4) * t49;
t106 = qJD(2) * qJD(1);
t103 = t52 * t135;
t102 = t52 * t127;
t55 = qJD(1) ^ 2;
t101 = t50 * t55 * t53;
t100 = t13 * t109;
t98 = t50 * t117;
t97 = t51 * t117;
t95 = t42 * t116;
t91 = t26 * t121;
t90 = 0.2e1 * t106;
t89 = t13 * t117 - t3 * t52;
t88 = t22 * t108;
t87 = t50 * t108;
t85 = t50 * t106;
t78 = (qJD(4) + t42) * t50;
t76 = t26 + t44;
t75 = -t28 + t111;
t74 = qJ(2) * t101;
t73 = 0.2e1 * t81;
t72 = -qJD(5) + t44;
t18 = qJD(1) * t102 - t51 * t121;
t70 = t48 * t116 - t18;
t67 = (-t42 + t109) * t50;
t66 = t122 - t135;
t65 = t122 + t135;
t64 = -t109 * t136 + t42 * t117 + t52 * t82;
t21 = t51 * t129 - t127;
t20 = t48 * t129 + t126;
t61 = t70 * t22;
t60 = -t22 * t113 - t134;
t59 = t65 * t52;
t15 = t51 * t25 + t48 * t87;
t2 = -qJD(5) * t15 + t79 * t126 - t48 * t7 + t51 * t85;
t24 = t49 * t86 - t110;
t14 = -t48 * t25 + t51 * t87;
t6 = t72 * t126 + (-t97 + (-qJD(5) * t52 + qJD(3)) * t48) * t50;
t5 = -t53 * t114 - t50 * t119 - t48 * t98 + t52 * t57;
t9 = [0, 0, 0, 0, t90, qJ(2) * t90, t50 * t73, -0.2e1 * t125 * t107, t54 * t53, -t54 * t50, 0, -t53 * t123, t50 * t123, t16 * t129 + (t44 * t53 - t98) * t28, (-t26 * t52 - t28 * t49) * t118 + (-t143 - t142 + (t26 * t49 - t137) * qJD(4)) * t50, t42 * t98 - t16 * t53 + (t52 * t66 + t138) * qJD(3), t42 * t94 + t17 * t53 + (-t49 * t66 - t140) * qJD(3), (-t42 - t109) * t120, -t24 * t120 + t8 * t53 + (t49 * t65 + t140) * qJD(2) + (t130 + qJD(4) * t59 + (t26 * t53 + t49 * t67) * qJD(3)) * qJ(2), -t25 * t120 + t7 * t53 + (t59 + t138) * qJD(2) + (t16 * t50 - t65 * t117 + (t28 * t53 + t52 * t67) * qJD(3)) * qJ(2), t13 * t6 + t3 * t21, -t6 * t11 - t13 * t5 - t3 * t20 - t21 * t4, t13 * t58 + t132 * t3 + t21 * t17 + t6 * t22, -t11 * t58 - t132 * t4 - t20 * t17 - t5 * t22, t130 * t49 + t22 * t58, t2 * t132 + t8 * t20 + t24 * t5 + t58 * t14 + ((t50 * t51 - t102) * t22 + t11 * t131) * qJD(2) + ((t22 * t48 * t72 - t11 * t111 + t128) * t50 + ((t117 * t48 + t119) * t22 + t49 * t4 + (qJD(4) * t11 + t60) * t52) * t53) * qJ(2), -t1 * t132 + t8 * t21 + t24 * t6 - t58 * t15 + (t13 * t131 - t22 * t63) * qJD(2) + ((-t111 * t13 + t72 * t77 - t134) * t50 + (-(qJD(3) * t48 - t97) * t22 + t149 + (qJD(4) * t13 + t22 * t114 - t128) * t52) * t53) * qJ(2); 0, 0, 0, 0, -t55, -t55 * qJ(2), 0, 0, 0, 0, 0, 0.2e1 * t82, t73, 0, 0, 0, 0, 0, t64 - t91, t95 + (-t103 + (-t28 - t111) * t50) * qJD(1), 0, 0, 0, 0, 0, -t148 - t61 + (t60 - t153) * t49, (-t100 - t128) * t49 - t152 + t89; 0, 0, 0, 0, 0, 0, -t101, t125 * t55, 0, 0, 0, -0.2e1 * t85, -0.2e1 * t53 * t106, -t137 * t42 + t143, (t16 + t141) * t52 + (-t17 + t139) * t49, -t95 + (t50 * t75 + t103) * qJD(1), t64 + t91, t42 * t121, -t49 * t74 + ((t24 - t110) * t50 + (t49 * t78 - t53 * t76) * qJ(2)) * qJD(1), -t52 * t74 + ((t25 + t112) * t50 + (t52 * t78 + t53 * t75) * qJ(2)) * qJD(1), t13 * t56 + t51 * t149, t19 * t11 + t13 * t18 + (-t11 * t51 - t13 * t48) * t116 + (-t150 - t4 * t51 + (t11 * t48 - t13 * t51) * qJD(5)) * t49, (-t100 + t128) * t49 + t152 + t89, t148 - t61 + (t60 + t153) * t49, -t22 * t136 - t142, -t2 * t52 + t70 * t24 - t20 * t88 + (t24 * t113 + t14 * qJD(4) + t147 + (t11 * t124 - t14 * t53) * qJD(1)) * t49, t1 * t52 + t69 * t24 - t21 * t88 + (-t24 * t114 - t15 * qJD(4) + t146 + (t124 * t13 + t15 * t53) * qJD(1)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t26, -t26 ^ 2 + t28 ^ 2, t16 - t141, -t139 - t17, t82, -t25 * t42 + (-t53 * t112 + (-t115 * t52 - t138) * qJ(2)) * qJD(1) - t68, -t83 + t24 * t42 + (-t99 + (t50 * t76 + t96) * qJ(2)) * qJD(1), t13 * t77 + t150, (t3 - t145) * t51 + (-t4 - t144) * t48, -t13 * t28 + t22 * t77 + t134, -t22 ^ 2 * t48 + t11 * t28 + t128, -t22 * t28, -t25 * t11 - t14 * t28 - t146, -t25 * t13 + t15 * t28 + t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t11, -t11 ^ 2 + t13 ^ 2, t3 + t145, t144 - t4, t17, -t24 * t13 + t15 * t22 + t2, t24 * t11 + t14 * t22 - t1;];
tauc_reg  = t9;
