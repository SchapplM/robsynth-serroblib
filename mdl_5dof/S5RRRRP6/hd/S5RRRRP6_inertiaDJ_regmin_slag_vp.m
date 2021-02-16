% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP6
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
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:11
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:10:24
% EndTime: 2021-01-16 00:10:31
% DurationCPUTime: 1.54s
% Computational Cost: add. (1863->209), mult. (4292->346), div. (0->0), fcn. (3734->6), ass. (0->118)
t82 = sin(qJ(3));
t120 = qJD(3) * t82;
t143 = pkin(7) + pkin(6);
t104 = qJD(2) * t143;
t83 = sin(qJ(2));
t53 = t83 * t104;
t86 = cos(qJ(2));
t54 = t86 * t104;
t63 = t143 * t83;
t64 = t143 * t86;
t85 = cos(qJ(3));
t90 = (-qJD(3) * t63 - t53) * t85 - t64 * t120 - t82 * t54;
t50 = t82 * t83 - t85 * t86;
t51 = t82 * t86 + t85 * t83;
t72 = -t86 * pkin(2) - pkin(1);
t94 = -t50 * pkin(3) + t51 * pkin(8) - t72;
t146 = qJD(4) * t94 - t90;
t36 = -t82 * t63 + t85 * t64;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t127 = t84 * t36 - t81 * t94;
t111 = pkin(2) * t120;
t76 = qJD(4) * t81;
t74 = pkin(4) * t76;
t55 = t74 + t111;
t139 = t85 * pkin(2);
t140 = t84 * pkin(4);
t71 = -pkin(3) - t140;
t60 = t71 - t139;
t77 = qJD(4) * t84;
t126 = t55 * t81 + t60 * t77;
t43 = t60 * t76;
t145 = -t55 * t84 + t43;
t79 = t81 ^ 2;
t80 = t84 ^ 2;
t123 = t79 - t80;
t100 = t123 * qJD(4);
t144 = qJD(2) + qJD(3);
t142 = pkin(4) * t79;
t33 = t144 * t51;
t141 = t33 * pkin(4);
t136 = t51 * t81;
t35 = t85 * t63 + t82 * t64;
t23 = pkin(4) * t136 + t35;
t16 = t36 * qJD(3) - t82 * t53 + t85 * t54;
t32 = t144 * t50;
t133 = t81 * t32;
t20 = t51 * t77 - t133;
t7 = t20 * pkin(4) + t16;
t138 = t23 * t77 + t7 * t81;
t137 = t51 * t32;
t135 = t51 * t84;
t132 = t81 * t33;
t131 = t84 * t32;
t130 = t84 * t33;
t129 = -qJ(5) - pkin(8);
t128 = t16 * t81 + t35 * t77;
t70 = -pkin(3) - t139;
t125 = t81 * t111 + t70 * t77;
t124 = qJD(4) * t142 + t71 * t77;
t122 = pkin(2) * qJD(3);
t78 = t84 * qJ(5);
t69 = t82 * pkin(2) + pkin(8);
t121 = qJ(5) + t69;
t119 = qJD(4) * t51;
t118 = t83 * qJD(2);
t117 = t86 * qJD(2);
t116 = -0.2e1 * pkin(1) * qJD(2);
t112 = pkin(2) * t118;
t92 = t33 * pkin(3) + t32 * pkin(8) + t112;
t115 = t146 * t84 - t81 * t92;
t114 = pkin(3) * t76;
t113 = pkin(3) * t77;
t110 = t85 * t122;
t109 = pkin(4) * t77;
t108 = t51 * t76;
t21 = t23 * t76;
t58 = t71 * t76;
t107 = t81 * t77;
t106 = qJ(5) * t119;
t105 = -t7 * t84 + t21;
t103 = -0.4e1 * t81 * t135;
t102 = -t81 * t36 - t84 * t94;
t98 = t50 * t69 - t51 * t70;
t97 = -t84 * t111 + t70 * t76;
t96 = t108 + t131;
t95 = t50 * t76 - t130;
t93 = -t32 * qJ(5) + qJD(4) * t36 + t51 * qJD(5);
t88 = -t32 * t70 - t33 * t69 + (-t50 * t85 + t51 * t82) * t122;
t12 = t84 * t92;
t87 = -t93 * t84 + t12 + (t106 + t146) * t81;
t75 = t84 * qJD(5);
t67 = 0.2e1 * t107;
t65 = t84 * t110;
t62 = t84 * pkin(8) + t78;
t61 = t129 * t81;
t49 = -0.2e1 * t100;
t48 = t51 ^ 2;
t46 = t84 * t69 + t78;
t45 = t121 * t81;
t40 = -t81 * qJD(5) + t129 * t77;
t39 = -t129 * t76 - t75;
t37 = t39 * t84;
t31 = (-qJD(5) - t110) * t81 - t121 * t77;
t30 = t121 * t76 - t65 - t75;
t28 = t35 * t76;
t27 = t30 * t84;
t19 = t50 * t77 + t132;
t13 = -t51 * t100 - t81 * t131;
t10 = -qJ(5) * t136 + t127;
t9 = qJD(4) * t103 + t123 * t32;
t8 = t50 * pkin(4) - t51 * t78 + t102;
t5 = -t127 * qJD(4) - t81 * t90 + t12;
t4 = t36 * t76 + t115;
t3 = t84 * t106 + t93 * t81 + t115;
t2 = t3 * t84;
t1 = t87 + t141;
t6 = [0, 0, 0, 0.2e1 * t83 * t117, 0.2e1 * (-t83 ^ 2 + t86 ^ 2) * qJD(2), 0, 0, 0, t83 * t116, t86 * t116, -0.2e1 * t137, 0.2e1 * t32 * t50 - 0.2e1 * t51 * t33, 0, 0, 0, 0.2e1 * t50 * t112 + 0.2e1 * t72 * t33, 0.2e1 * t51 * t112 - 0.2e1 * t72 * t32, -0.2e1 * t107 * t48 - 0.2e1 * t80 * t137, 0.2e1 * t48 * t100 - t32 * t103, 0.2e1 * t51 * t130 - 0.2e1 * t50 * t96, -0.2e1 * t51 * t132 - 0.2e1 * t20 * t50, 0.2e1 * t50 * t33, 0.2e1 * t102 * t33 + 0.2e1 * t16 * t136 + 0.2e1 * t20 * t35 + 0.2e1 * t5 * t50, -0.2e1 * t127 * t33 + 0.2e1 * t16 * t135 - 0.2e1 * t96 * t35 + 0.2e1 * t4 * t50, 0.2e1 * t1 * t50 + 0.2e1 * t7 * t136 + 0.2e1 * t20 * t23 + 0.2e1 * t8 * t33, -0.2e1 * t10 * t33 + 0.2e1 * t7 * t135 - 0.2e1 * t23 * t96 + 0.2e1 * t3 * t50, -0.2e1 * (-t10 * t81 - t8 * t84) * t32 + 0.2e1 * (-t1 * t84 + t3 * t81 + (-t10 * t84 + t8 * t81) * qJD(4)) * t51, 0.2e1 * t8 * t1 - 0.2e1 * t10 * t3 + 0.2e1 * t23 * t7; 0, 0, 0, 0, 0, t117, -t118, 0, -pkin(6) * t117, pkin(6) * t118, 0, 0, -t32, -t33, 0, -t16, -t90, t13, t9, t19, -t95, 0, t28 + (-qJD(4) * t98 - t16) * t84 + t88 * t81, t98 * t76 + t84 * t88 + t128, t126 * t51 - t60 * t133 + t31 * t50 - t45 * t33 + t105, -t60 * t131 - t145 * t51 + t30 * t50 - t46 * t33 + t138, -t2 + (-t31 * t51 - t32 * t45 + (-t46 * t51 - t8) * qJD(4)) * t84 + (t30 * t51 + t32 * t46 - t1 + (-t45 * t51 - t10) * qJD(4)) * t81, -t1 * t45 - t10 * t30 + t23 * t55 - t3 * t46 + t8 * t31 + t7 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t111, -0.2e1 * t110, t67, t49, 0, 0, 0, 0.2e1 * t97, 0.2e1 * t125, 0.2e1 * t145, 0.2e1 * t126, -0.2e1 * t31 * t81 - 0.2e1 * t27 + 0.2e1 * (t45 * t84 - t46 * t81) * qJD(4), -0.2e1 * t46 * t30 - 0.2e1 * t45 * t31 + 0.2e1 * t60 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t33, 0, -t16, -t90, t13, t9, t19, -t95, 0, t28 + (pkin(3) * t32 - pkin(8) * t33) * t81 + (-t16 + (-pkin(3) * t51 - pkin(8) * t50) * qJD(4)) * t84, pkin(3) * t96 + pkin(8) * t95 + t128, -t71 * t133 + t61 * t33 + t40 * t50 + (t71 * t84 + t142) * t119 + t105, -t71 * t131 - t62 * t33 + t39 * t50 + (-t71 + t140) * t108 + t138, -t2 + (t32 * t61 - t40 * t51 + (-t51 * t62 - t8) * qJD(4)) * t84 + (t32 * t62 + t39 * t51 - t1 + (t51 * t61 - t10) * qJD(4)) * t81, pkin(4) * t21 + t1 * t61 - t10 * t39 - t3 * t62 + t8 * t40 + t7 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t110, t67, t49, 0, 0, 0, t97 - t114, -t113 + t125, t43 + t58 + (-t55 - t74) * t84, t124 + t126, -t27 - t37 + (-t31 - t40) * t81 + ((t45 - t61) * t84 + (-t46 - t62) * t81) * qJD(4), pkin(4) * t43 - t30 * t62 + t31 * t61 - t46 * t39 - t45 * t40 + t55 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t49, 0, 0, 0, -0.2e1 * t114, -0.2e1 * t113, -0.2e1 * pkin(4) * t107 + 0.2e1 * t58, 0.2e1 * t124, -0.2e1 * t40 * t81 - 0.2e1 * t37 + 0.2e1 * (-t61 * t84 - t62 * t81) * qJD(4), 0.2e1 * pkin(4) * t58 - 0.2e1 * t62 * t39 + 0.2e1 * t61 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t20, t33, t5, t4, t87 + 0.2e1 * t141, t3, t96 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t76, 0, -t81 * t110 - t69 * t77, t69 * t76 - t65, t31, t30, -t109, t31 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t76, 0, -pkin(8) * t77, pkin(8) * t76, t40, t39, -t109, t40 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t96, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t77, 0, t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t77, 0, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
