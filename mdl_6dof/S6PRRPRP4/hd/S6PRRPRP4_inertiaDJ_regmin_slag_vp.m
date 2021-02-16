% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRP4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:47
% EndTime: 2021-01-16 03:12:54
% DurationCPUTime: 1.39s
% Computational Cost: add. (993->192), mult. (2486->340), div. (0->0), fcn. (2085->8), ass. (0->116)
t130 = pkin(4) + pkin(8);
t131 = pkin(3) + pkin(9);
t117 = qJ(6) + t131;
t66 = cos(qJ(3));
t119 = qJ(4) * t66;
t63 = sin(qJ(3));
t135 = t131 * t63 - t119;
t120 = qJ(4) * t63;
t134 = t131 * t66 + t120;
t59 = t66 ^ 2;
t86 = qJD(3) * (t63 ^ 2 - t59);
t62 = sin(qJ(5));
t56 = t62 ^ 2;
t65 = cos(qJ(5));
t123 = -t65 ^ 2 + t56;
t85 = t123 * qJD(5);
t45 = t130 * t66;
t106 = t45 * qJD(5);
t111 = qJD(4) * t66;
t133 = t134 * qJD(3) - t106 - t111;
t132 = 0.2e1 * qJD(4);
t129 = pkin(5) * t65;
t60 = sin(pkin(6));
t64 = sin(qJ(2));
t127 = t60 * t64;
t100 = t63 * t127;
t61 = cos(pkin(6));
t67 = cos(qJ(2));
t115 = qJD(2) * t67;
t92 = t60 * t115;
t15 = -qJD(3) * t100 + (qJD(3) * t61 + t92) * t66;
t29 = t66 * t127 + t61 * t63;
t128 = t29 * t15;
t126 = t60 * t67;
t125 = t62 * t63;
t36 = -pkin(2) - t134;
t96 = t130 * t63;
t124 = t65 * t36 + t62 * t96;
t121 = qJ(4) * t62;
t118 = qJ(6) * t66;
t116 = qJD(2) * t64;
t114 = qJD(3) * t62;
t113 = qJD(3) * t65;
t112 = qJD(4) * t63;
t110 = qJD(5) * t62;
t109 = qJD(5) * t65;
t108 = qJD(5) * t66;
t107 = qJD(5) * t131;
t105 = t62 * qJD(6);
t104 = t63 * qJD(3);
t103 = t65 * qJD(6);
t53 = t66 * qJD(3);
t102 = qJ(4) * qJD(5);
t101 = -0.2e1 * pkin(2) * qJD(3);
t99 = pkin(5) * t110;
t98 = pkin(8) * t104;
t97 = pkin(8) * t53;
t95 = t62 * t108;
t94 = t65 * t108;
t93 = t60 * t116;
t91 = t62 * t53;
t90 = t63 * t53;
t89 = t65 * t104;
t88 = t62 * t109;
t87 = qJ(6) * t108;
t42 = t117 * t65;
t39 = t65 * t96;
t84 = t62 * t89;
t83 = t130 * t53;
t19 = -pkin(5) * t95 + (-t129 - t130) * t104;
t52 = t62 * pkin(5) + qJ(4);
t82 = -t52 * t108 + t19;
t81 = -t36 * t109 + t65 * t83;
t80 = -t66 * pkin(3) - t120;
t13 = t63 * pkin(5) + t39 + (-t36 + t118) * t62;
t14 = -t65 * t118 + t124;
t79 = -t13 * t62 + t14 * t65;
t28 = -t61 * t66 + t100;
t17 = t62 * t126 + t28 * t65;
t74 = t65 * t126 - t28 * t62;
t78 = -t17 * t62 - t65 * t74;
t76 = -qJD(5) * t130 + qJD(4);
t11 = t29 * t109 + t15 * t62;
t12 = -t29 * t110 + t15 * t65;
t7 = t36 * t110 - t65 * (t135 * qJD(3) - t112) - t62 * t83 - qJD(5) * t39;
t33 = t62 * t104 - t94;
t40 = t130 * t104;
t73 = t135 * qJD(5) - t40;
t27 = t66 * t129 + t45;
t50 = pkin(5) * t109 + qJD(4);
t72 = -qJD(5) * t27 + t52 * t104 - t50 * t66;
t71 = t80 * qJD(3) + t111;
t16 = qJD(3) * t29 + t63 * t92;
t70 = t15 * t66 + t16 * t63 + (t28 * t66 - t29 * t63) * qJD(3);
t69 = t65 * t87 + t66 * t105 + (-t117 * qJD(3) + t76) * t125 + t81;
t49 = 0.2e1 * t90;
t43 = -pkin(2) + t80;
t41 = t117 * t62;
t32 = -t63 * t109 - t91;
t31 = t89 + t95;
t30 = -t63 * t110 + t65 * t53;
t26 = -t112 + (pkin(3) * t63 - t119) * qJD(3);
t24 = -qJD(5) * t42 - t105;
t23 = t117 * t110 - t103;
t22 = (t67 * t104 + t66 * t116) * t60;
t21 = (t63 * t116 - t67 * t53) * t60;
t10 = t17 * qJD(5) + t16 * t62 + t65 * t93;
t9 = t74 * qJD(5) + t16 * t65 - t62 * t93;
t8 = qJ(4) * t91 + (-t131 * qJD(3) + t76) * t125 + t81;
t6 = t23 * t65 + t24 * t62 + (-t41 * t65 + t42 * t62) * qJD(5);
t5 = -qJ(6) * t89 + t66 * t103 - t62 * t87 + t7;
t4 = (pkin(5) + t121) * t53 + t69;
t3 = (t29 * t114 - t10) * t63 + (qJD(3) * t74 - t11) * t66;
t2 = (-t29 * t113 + t9) * t63 + (qJD(3) * t17 + t12) * t66;
t1 = t78 * qJD(5) + t10 * t62 + t9 * t65;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t60 ^ 2 * t64 * t115 + 0.2e1 * t28 * t16 + 0.2e1 * t128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t10 * t74 + 0.2e1 * t17 * t9 + 0.2e1 * t128; 0, 0, -t93, -t92, 0, 0, 0, 0, 0, -t22, t21, t70, t22, -t21, (t43 * t116 - t26 * t67) * t60 + t70 * pkin(8), 0, 0, 0, 0, 0, t2, t3, t2, t3, t78 * t104 + (-t10 * t65 + t62 * t9 + (t17 * t65 - t62 * t74) * qJD(5)) * t66, t10 * t14 + t9 * t13 + t15 * t27 + t17 * t4 + t29 * t19 + t5 * t74; 0, 0, 0, 0, t49, -0.2e1 * t86, 0, 0, 0, t63 * t101, t66 * t101, 0, -0.2e1 * t43 * t104 + 0.2e1 * t26 * t66, -0.2e1 * t26 * t63 - 0.2e1 * t43 * t53, 0.2e1 * t43 * t26, -0.2e1 * t56 * t90 + 0.2e1 * t59 * t88, -0.2e1 * t59 * t85 - 0.4e1 * t66 * t84, 0.2e1 * t62 * t86 - 0.2e1 * t63 * t94, 0.2e1 * t63 * t95 + 0.2e1 * t65 * t86, t49, 0.2e1 * (-t45 * t113 + t8) * t63 + 0.2e1 * ((-t62 * t36 + t39) * qJD(3) - t40 * t65 - t62 * t106) * t66, 0.2e1 * (t45 * t114 + t7) * t63 + 0.2e1 * (-t124 * qJD(3) - t65 * t106 + t40 * t62) * t66, 0.2e1 * (-t27 * t113 + t4) * t63 + 0.2e1 * (qJD(3) * t13 - t27 * t110 + t19 * t65) * t66, 0.2e1 * (t27 * t114 + t5) * t63 + 0.2e1 * (-qJD(3) * t14 - t27 * t109 - t19 * t62) * t66, 0.2e1 * t79 * t104 + 0.2e1 * (t4 * t62 + t5 * t65 + (t13 * t65 + t14 * t62) * qJD(5)) * t66, 0.2e1 * t13 * t4 - 0.2e1 * t14 * t5 + 0.2e1 * t27 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t15, 0, t16, t15, -t16 * pkin(3) + t15 * qJ(4) + t29 * qJD(4), 0, 0, 0, 0, 0, t11, t12, t11, t12, -t1, -t10 * t41 + t15 * t52 + t17 * t23 - t24 * t74 + t29 * t50 - t9 * t42; 0, 0, 0, 0, 0, 0, t53, -t104, 0, -t97, t98, t71, t97, -t98, t71 * pkin(8), t66 * t85 + t84, -t123 * t104 + 0.4e1 * t66 * t88, t30, t32, 0, -t133 * t65 + t73 * t62, t133 * t62 + t73 * t65, t23 * t63 - t42 * t53 + t82 * t62 - t72 * t65, -t24 * t63 + t41 * t53 + t72 * t62 + t82 * t65, (-t41 * t104 - t24 * t66 - t4 + (-t42 * t66 - t14) * qJD(5)) * t65 + (t42 * t104 + t23 * t66 + t5 + (-t41 * t66 + t13) * qJD(5)) * t62, t13 * t23 + t14 * t24 + t19 * t52 + t27 * t50 - t4 * t42 + t5 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, qJ(4) * t132, -0.2e1 * t88, 0.2e1 * t85, 0, 0, 0, 0.2e1 * qJD(4) * t62 + 0.2e1 * t65 * t102, 0.2e1 * qJD(4) * t65 - 0.2e1 * t62 * t102, 0.2e1 * t52 * t109 + 0.2e1 * t50 * t62, -0.2e1 * t52 * t110 + 0.2e1 * t50 * t65, -0.2e1 * t6, -0.2e1 * t42 * t23 - 0.2e1 * t41 * t24 + 0.2e1 * t52 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, 0, t97, 0, 0, 0, 0, 0, t30, t32, t30, t32, 0, qJD(5) * t79 + t4 * t65 - t5 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t10, t9, -t10, 0, t9 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t31, t53, t8, t7, (0.2e1 * pkin(5) + t121) * t53 + t69, t5, -t33 * pkin(5), t4 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109, 0, t62 * t107, t65 * t107, t23, -t24, t99, t23 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t109, -t110, -t109, 0, -t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t33, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t110, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t18;
