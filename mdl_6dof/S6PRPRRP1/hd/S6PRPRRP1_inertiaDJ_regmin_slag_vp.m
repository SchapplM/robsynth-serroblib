% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x23]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRPRRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:18
% EndTime: 2021-01-16 01:24:24
% DurationCPUTime: 1.29s
% Computational Cost: add. (1080->180), mult. (2980->330), div. (0->0), fcn. (2789->10), ass. (0->122)
t56 = sin(qJ(4));
t138 = -0.4e1 * t56;
t59 = cos(qJ(4));
t109 = t59 * qJD(4);
t58 = cos(qJ(5));
t47 = qJD(5) * t58;
t55 = sin(qJ(5));
t30 = t55 * t109 + t56 * t47;
t137 = t30 * pkin(5);
t53 = sin(pkin(6));
t52 = sin(pkin(11));
t54 = cos(pkin(11));
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t79 = t52 * t60 + t54 * t57;
t71 = qJD(4) * t79;
t78 = t52 * t57 - t54 * t60;
t72 = qJD(2) * t78;
t116 = cos(pkin(6));
t89 = qJD(4) * t116;
t135 = (t56 * t71 + t59 * t72) * t53 - t59 * t89;
t74 = t78 * t53;
t136 = -qJD(5) * t74 + t135;
t42 = t52 * pkin(2) + pkin(8);
t134 = t56 * qJD(6) + (qJ(6) * qJD(4) + qJD(5) * t42) * t59;
t133 = 0.2e1 * qJD(5);
t132 = t58 * pkin(5);
t93 = t42 * t109;
t19 = t93 + t137;
t131 = t19 * t55;
t130 = t19 * t58;
t123 = -qJ(6) - pkin(9);
t37 = t123 * t55;
t129 = t37 * t56;
t38 = t123 * t58;
t128 = t38 * t56;
t127 = t42 * t55;
t126 = t55 * t59;
t125 = t56 * t58;
t124 = t58 * t59;
t43 = -t54 * pkin(2) - pkin(3);
t87 = -t59 * pkin(4) - t56 * pkin(9);
t67 = t43 + t87;
t66 = qJD(5) * t67;
t86 = pkin(4) * t56 - pkin(9) * t59;
t76 = t86 * qJD(4);
t122 = -t55 * t76 - t58 * t66;
t45 = t56 * qJD(4);
t95 = t42 * t45;
t121 = t55 * t95 + t58 * t76;
t23 = t55 * t67;
t35 = t42 * t124;
t120 = t35 + t23;
t48 = t55 ^ 2;
t50 = t58 ^ 2;
t119 = t48 - t50;
t49 = t56 ^ 2;
t118 = -t59 ^ 2 + t49;
t117 = qJ(6) * t56;
t115 = qJD(2) * t53;
t114 = qJD(4) * t55;
t113 = qJD(4) * t58;
t112 = qJD(5) * t55;
t111 = qJD(5) * t59;
t108 = -0.2e1 * pkin(4) * qJD(5);
t107 = 0.2e1 * qJD(4) * t43;
t106 = pkin(5) * t45;
t105 = pkin(5) * t112;
t104 = t56 * t112;
t103 = t55 * t111;
t102 = t58 * t111;
t75 = t79 * t53;
t17 = -t116 * t59 + t56 * t75;
t101 = t17 * t112;
t27 = (pkin(5) * t55 + t42) * t56;
t100 = t27 * t112;
t99 = t48 * t109;
t98 = t57 * t115;
t97 = t55 * t47;
t96 = t56 * t109;
t94 = t58 * t109;
t44 = -pkin(4) - t132;
t92 = -t44 + t132;
t91 = t119 * qJD(5);
t90 = t118 * qJD(4);
t88 = t55 * t94;
t84 = pkin(5) * t48 + t44 * t58;
t18 = t116 * t56 + t59 * t75;
t13 = -t55 * t18 + t58 * t74;
t14 = t58 * t18 + t55 * t74;
t83 = -t13 * t58 - t14 * t55;
t82 = t13 * t55 - t14 * t58;
t24 = t58 * t67;
t15 = -t58 * t117 + t24 + (-pkin(5) - t127) * t59;
t16 = -t55 * t117 + t120;
t81 = -t15 * t58 - t16 * t55;
t80 = t15 * t55 - t16 * t58;
t12 = t56 * t89 + (-t56 * t72 + t59 * t71) * t53;
t5 = t12 * t55 + t17 * t47;
t6 = -t12 * t58 + t101;
t73 = qJD(2) * t79;
t70 = qJD(4) * t78;
t28 = -t94 + t104;
t29 = t58 * t45 + t103;
t69 = t53 * t73;
t3 = t18 * t112 + t136 * t58 - t55 * t69;
t4 = t136 * t55 - t18 * t47 + t58 * t69;
t65 = t83 * qJD(5) - t3 * t58 - t4 * t55;
t25 = -t58 * qJD(6) - t123 * t112;
t26 = -t55 * qJD(6) + t123 * t47;
t64 = -t25 * t58 - t26 * t55 + (-t37 * t58 + t38 * t55) * qJD(5);
t61 = qJ(6) * t104 - t134 * t58 - t55 * t66 + t121;
t40 = t50 * t109;
t36 = t50 * t96;
t31 = t55 * t45 - t102;
t10 = -t120 * qJD(5) + t121;
t9 = t29 * t42 + t122;
t8 = (qJ(6) * qJD(5) + qJD(4) * t42) * t125 + t134 * t55 + t122;
t7 = t61 + t106;
t2 = (t17 * t114 - t4) * t59 + (qJD(4) * t13 + t5) * t56;
t1 = (t17 * t113 - t3) * t59 + (-qJD(4) * t14 - t6) * t56;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t17 * t12 + 0.2e1 * t13 * t4 - 0.2e1 * t14 * t3; 0, 0, -t98, -t60 * t115, (-t52 ^ 2 - t54 ^ 2) * pkin(2) * t98, 0, 0, 0, 0, 0, (t56 * t70 - t59 * t73) * t53, (t56 * t73 + t59 * t70) * t53, 0, 0, 0, 0, 0, t2, t1, t2, t1, t83 * t109 + (t82 * qJD(5) + t3 * t55 - t4 * t58) * t56, t12 * t27 + t13 * t7 - t14 * t8 + t4 * t15 - t3 * t16 + t17 * t19; 0, 0, 0, 0, 0, 0.2e1 * t96, -0.2e1 * t90, 0, 0, 0, t56 * t107, t59 * t107, -0.2e1 * t49 * t97 + 0.2e1 * t36, t119 * t49 * t133 + t88 * t138, 0.2e1 * t56 * t103 + 0.2e1 * t118 * t113, 0.2e1 * t56 * t102 - 0.2e1 * t55 * t90, -0.2e1 * t96, 0.2e1 * t24 * t45 - 0.2e1 * t10 * t59 + 0.2e1 * (t49 * t47 + t55 * t96) * t42, -0.2e1 * t49 * t42 * t112 - 0.2e1 * t9 * t59 + 0.2e1 * (-t23 + t35) * t45, 0.2e1 * (t27 * t114 - t7) * t59 + 0.2e1 * (qJD(4) * t15 + t27 * t47 + t131) * t56, 0.2e1 * (t27 * t113 - t8) * t59 + 0.2e1 * (-qJD(4) * t16 - t100 + t130) * t56, 0.2e1 * t81 * t109 + 0.2e1 * (t80 * qJD(5) + t55 * t8 - t58 * t7) * t56, 0.2e1 * t15 * t7 - 0.2e1 * t16 * t8 + 0.2e1 * t27 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t82 * qJD(4) - t12) * t59 + (qJD(4) * t17 + t65) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t80 * qJD(4) - t19) * t59 + (qJD(4) * t27 + qJD(5) * t81 - t7 * t55 - t8 * t58) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t36 + 0.2e1 * (t48 - 0.1e1) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t135, 0, 0, 0, 0, 0, t6, t5, t6, t5, t65, pkin(5) * t101 + t12 * t44 + t13 * t26 - t14 * t25 + t3 * t38 + t4 * t37; 0, 0, 0, 0, 0, 0, 0, t109, -t45, 0, -t93, t95, -t56 * t91 + t88, t97 * t138 + t40 - t99, t31, t29, 0, (pkin(9) * t124 + (-pkin(4) * t58 + t127) * t56) * qJD(5) + (t87 * t55 - t35) * qJD(4), (t42 * t125 + t86 * t55) * qJD(5) + (t42 * t126 + t87 * t58) * qJD(4), -t130 - t26 * t59 + (t44 * t126 + t129) * qJD(4) + (t27 * t55 + t84 * t56) * qJD(5), t131 - t25 * t59 + (t44 * t124 + t128) * qJD(4) + (t92 * t56 * t55 + t27 * t58) * qJD(5), (-t37 * t109 - t26 * t56 - t8 + (-t15 + t128) * qJD(5)) * t58 + (t38 * t109 + t25 * t56 - t7 + (-t16 + t129) * qJD(5)) * t55, pkin(5) * t100 + t15 * t26 - t16 * t25 + t19 * t44 + t7 * t37 + t8 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, -t109, 0, 0, 0, 0, 0, -t29, t31, -t29, t31, t40 + t99, (-t105 + (-t37 * t55 - t38 * t58) * qJD(4)) * t59 + (qJD(4) * t44 + t64) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t97, -0.2e1 * t91, 0, 0, 0, t55 * t108, t58 * t108, -0.2e1 * t92 * t112, t84 * t133, 0.2e1 * t64, 0.2e1 * t44 * t105 + 0.2e1 * t38 * t25 + 0.2e1 * t37 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, t4, t3, 0, t4 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t30, t45, t10, t9, t61 + 0.2e1 * t106, t8, t28 * pkin(5), t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t28, -t30, t28, 0, -t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t112, 0, -pkin(9) * t47, pkin(9) * t112, t26, t25, -pkin(5) * t47, t26 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t47, 0, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
