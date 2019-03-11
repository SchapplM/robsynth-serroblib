% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:15:16
% EndTime: 2019-03-09 09:15:19
% DurationCPUTime: 1.17s
% Computational Cost: add. (1732->162), mult. (3756->303), div. (0->0), fcn. (3697->8), ass. (0->113)
t77 = cos(qJ(6));
t74 = t77 ^ 2;
t75 = sin(qJ(6));
t124 = t75 ^ 2 - t74;
t139 = qJD(6) * t124;
t142 = 0.2e1 * t139;
t76 = sin(qJ(2));
t78 = cos(qJ(2));
t141 = -t78 * pkin(2) - t76 * qJ(3);
t122 = sin(pkin(10));
t123 = cos(pkin(10));
t92 = t76 * t122 + t78 * t123;
t140 = t92 * qJD(2);
t134 = sin(qJ(5));
t135 = cos(qJ(5));
t49 = t134 * t122 - t135 * t123;
t108 = t78 * t122;
t109 = t76 * t123;
t50 = t109 - t108;
t86 = t135 * t92;
t29 = t134 * t50 + t86;
t85 = t134 * t92;
t30 = t135 * t50 - t85;
t79 = -pkin(2) - pkin(3);
t57 = t123 * qJ(3) + t122 * t79;
t56 = -t122 * qJ(3) + t123 * t79;
t91 = -pkin(4) + t56;
t84 = t135 * t91;
t31 = t134 * t57 + pkin(5) - t84;
t82 = t134 * t91 + t135 * t57;
t32 = -pkin(9) + t82;
t111 = qJD(5) * t134;
t112 = qJD(5) * t135;
t126 = pkin(7) - qJ(4);
t60 = t126 * t76;
t61 = t126 * t78;
t34 = -t122 * t61 + t123 * t60;
t23 = -t50 * pkin(8) + t34;
t35 = t122 * t60 + t123 * t61;
t24 = -t92 * pkin(8) + t35;
t119 = t76 * qJD(2);
t43 = -t78 * qJD(4) - t126 * t119;
t93 = qJD(2) * t61 - t76 * qJD(4);
t22 = t122 * t93 + t123 * t43;
t65 = qJD(2) * t108;
t94 = -qJD(2) * t109 + t65;
t80 = -t94 * pkin(8) + t22;
t21 = t122 * t43 - t123 * t93;
t81 = pkin(8) * t140 + t21;
t4 = t23 * t111 + t24 * t112 + t134 * t80 + t135 * t81;
t138 = -t4 + (t29 * t32 - t30 * t31) * qJD(6);
t137 = 2 * qJD(3);
t136 = t4 * t75;
t51 = t135 * t122 + t134 * t123;
t19 = t51 * qJD(3) + t82 * qJD(5);
t133 = t19 * t75;
t132 = t19 * t77;
t14 = qJD(5) * t86 + t50 * t111 + t134 * t94 - t135 * t140;
t131 = t30 * t14;
t130 = t30 * t77;
t15 = -qJD(5) * t85 + t50 * t112 + t134 * t140 + t135 * t94;
t129 = t75 * t15;
t128 = t77 * t14;
t127 = t77 * t15;
t69 = t78 * qJD(2);
t125 = qJ(3) * t69 + t76 * qJD(3);
t121 = qJD(6) * t75;
t120 = qJD(6) * t77;
t118 = -0.2e1 * pkin(1) * qJD(2);
t117 = -0.2e1 * pkin(5) * qJD(6);
t58 = -pkin(1) + t141;
t116 = pkin(7) * t119;
t115 = pkin(7) * t69;
t114 = t75 * t120;
t113 = 0.4e1 * t75 * t130;
t110 = qJD(6) * (pkin(5) + t31);
t48 = t78 * pkin(3) - t58;
t106 = t123 * qJD(3);
t105 = t122 * qJD(3);
t104 = pkin(5) * t14 - pkin(9) * t15;
t103 = pkin(5) * t30 + pkin(9) * t29;
t33 = t92 * pkin(4) + t48;
t10 = t29 * pkin(5) - t30 * pkin(9) + t33;
t12 = t134 * t23 + t135 * t24;
t100 = t77 * t10 - t75 * t12;
t99 = t75 * t10 + t77 * t12;
t97 = t29 * t51 - t30 * t49;
t96 = t30 * t120 - t75 * t14;
t95 = -t30 * t121 - t128;
t9 = t29 * t120 + t129;
t46 = t49 * qJD(5);
t47 = t51 * qJD(5);
t90 = -t14 * t49 - t15 * t51 + t29 * t46 + t30 * t47;
t87 = t141 * qJD(2) + t78 * qJD(3);
t11 = t134 * t24 - t135 * t23;
t18 = t49 * qJD(3) - qJD(5) * t84 + t57 * t111;
t83 = -qJD(6) * t11 - t14 * t31 - t15 * t32 + t18 * t29 + t19 * t30;
t25 = t65 * pkin(4) + (-t123 * pkin(4) + t79) * t119 + t125;
t62 = 0.2e1 * t114;
t55 = -0.2e1 * t139;
t44 = pkin(2) * t119 - t125;
t36 = t79 * t119 + t125;
t28 = t30 ^ 2;
t27 = -t49 * t121 + t47 * t77;
t26 = t49 * t120 + t47 * t75;
t8 = t29 * t121 - t127;
t7 = t75 * t128 + t139 * t30;
t6 = t15 * pkin(5) + t14 * pkin(9) + t25;
t5 = qJD(6) * t113 - t124 * t14;
t3 = t24 * t111 - t23 * t112 + t134 * t81 - t135 * t80;
t2 = -t99 * qJD(6) + t75 * t3 + t77 * t6;
t1 = -t100 * qJD(6) + t77 * t3 - t75 * t6;
t13 = [0, 0, 0, 0.2e1 * t76 * t69, 0.2e1 * (-t76 ^ 2 + t78 ^ 2) * qJD(2), 0, 0, 0, t76 * t118, t78 * t118, 0.2e1 * t58 * t119 - 0.2e1 * t44 * t78, 0, -0.2e1 * t44 * t76 - 0.2e1 * t58 * t69, 0.2e1 * t58 * t44, 0.2e1 * t36 * t92 + 0.2e1 * t48 * t94, 0.2e1 * t140 * t48 + 0.2e1 * t36 * t50, -0.2e1 * t140 * t34 + 0.2e1 * t21 * t50 - 0.2e1 * t22 * t92 - 0.2e1 * t35 * t94, -0.2e1 * t34 * t21 + 0.2e1 * t35 * t22 + 0.2e1 * t48 * t36, -0.2e1 * t131, 0.2e1 * t14 * t29 - 0.2e1 * t30 * t15, 0, 0, 0, 0.2e1 * t33 * t15 + 0.2e1 * t25 * t29, -0.2e1 * t33 * t14 + 0.2e1 * t25 * t30, -0.2e1 * t28 * t114 - 0.2e1 * t74 * t131, t14 * t113 + t142 * t28, 0.2e1 * t30 * t127 + 0.2e1 * t95 * t29, -0.2e1 * t30 * t129 - 0.2e1 * t96 * t29, 0.2e1 * t29 * t15, 0.2e1 * t100 * t15 + 0.2e1 * t96 * t11 + 0.2e1 * t30 * t136 + 0.2e1 * t2 * t29, 0.2e1 * t1 * t29 + 0.2e1 * t95 * t11 + 0.2e1 * t4 * t130 - 0.2e1 * t99 * t15; 0, 0, 0, 0, 0, t69, -t119, 0, -t115, t116, -t115, t87, -t116, t87 * pkin(7), t21, t22, t50 * t105 - t92 * t106 - t140 * t56 - t57 * t94, -t21 * t56 + t22 * t57 + (-t122 * t34 + t123 * t35) * qJD(3), 0, 0, t14, t15, 0, t4, -t3, t7, t5, -t9, t8, 0, -t138 * t77 + t83 * t75, t138 * t75 + t83 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, qJ(3) * t137, 0.2e1 * t105, 0.2e1 * t106, 0 (-t122 * t56 + t123 * t57) * t137, 0, 0, 0, 0, 0, 0.2e1 * t19, -0.2e1 * t18, t62, t55, 0, 0, 0, -0.2e1 * t31 * t121 + 0.2e1 * t132, -0.2e1 * t31 * t120 - 0.2e1 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t115, 0, 0, -t123 ^ 2 * t69 - t122 * t65, t22 * t122 - t21 * t123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t120 + t90 * t75, t97 * t121 + t90 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t46, 0, 0, 0, 0, 0, t27, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t140, 0, t36, 0, 0, 0, 0, 0, t15, -t14, 0, 0, 0, 0, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, -t4, t3, -t7, -t5, t9, -t8, 0, -t4 * t77 + t104 * t75 + (-t103 * t77 + t11 * t75) * qJD(6), t136 + t104 * t77 + (t103 * t75 + t11 * t77) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, t18, -0.2e1 * t114, t142, 0, 0, 0, t75 * t110 - t132, t77 * t110 + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t46, 0, 0, 0, 0, 0, -t27, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t55, 0, 0, 0, t75 * t117, t77 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t96, t15, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, t121, 0, -t32 * t120 + t75 * t18, t121 * t32 + t77 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t120 + t75 * t46, t121 * t51 + t77 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, -t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t121, 0, -pkin(9) * t120, pkin(9) * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t13;
