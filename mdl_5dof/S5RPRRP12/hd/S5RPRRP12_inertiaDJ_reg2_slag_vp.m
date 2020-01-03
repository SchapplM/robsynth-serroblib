% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaDJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:25
% DurationCPUTime: 1.26s
% Computational Cost: add. (749->158), mult. (1729->273), div. (0->0), fcn. (1203->4), ass. (0->109)
t59 = cos(qJ(3));
t57 = sin(qJ(3));
t98 = t57 * qJD(3);
t123 = -qJ(5) * t98 + t59 * qJD(5);
t56 = sin(qJ(4));
t52 = t56 ^ 2;
t58 = cos(qJ(4));
t54 = t58 ^ 2;
t108 = t52 - t54;
t122 = 2 * qJD(2);
t121 = 2 * qJD(4);
t120 = t57 * pkin(3);
t119 = t58 * pkin(4);
t118 = t59 * pkin(7);
t100 = qJD(4) * t59;
t88 = t58 * t100;
t28 = -t56 * t98 + t88;
t60 = -pkin(1) - pkin(6);
t85 = t60 * t98;
t13 = t28 * pkin(4) + t85;
t117 = t13 * t56;
t116 = t13 * t58;
t109 = -qJ(5) - pkin(7);
t35 = t109 * t56;
t115 = t35 * t59;
t36 = t109 * t58;
t114 = t36 * t59;
t47 = -pkin(3) - t119;
t113 = t47 * t58;
t112 = t56 * t60;
t111 = t57 * t60;
t110 = t59 * t60;
t40 = t58 * t111;
t73 = -t118 + t120;
t67 = qJ(2) + t73;
t17 = t56 * t67 + t40;
t107 = t52 + t54;
t53 = t57 ^ 2;
t55 = t59 ^ 2;
t106 = t53 - t55;
t105 = t53 + t55;
t104 = qJ(5) * t59;
t32 = (pkin(4) * t56 - t60) * t59;
t103 = qJD(3) * t32;
t102 = qJD(3) * t58;
t101 = qJD(4) * t56;
t50 = qJD(4) * t58;
t99 = qJD(4) * t60;
t49 = t59 * qJD(3);
t96 = qJ(2) * qJD(3);
t95 = t56 * t111;
t94 = -0.2e1 * t101;
t31 = t58 * t67;
t74 = pkin(3) * t59 + pkin(7) * t57;
t65 = t74 * qJD(3) + qJD(2);
t82 = t60 * t49;
t93 = -qJD(4) * t31 - t56 * t65 - t58 * t82;
t92 = pkin(4) * t101;
t91 = t58 * t104;
t90 = t56 * t100;
t89 = t56 * t99;
t87 = t32 * t101;
t86 = t56 * t50;
t84 = t58 * t98;
t83 = t57 * t49;
t80 = -t47 + t119;
t79 = pkin(4) - t112;
t78 = t107 * t59;
t77 = t106 * qJD(3);
t43 = 0.2e1 * t83;
t76 = t56 * t84;
t75 = t55 * t86;
t5 = t79 * t57 + t31 - t91;
t8 = -t56 * t104 + t17;
t72 = t5 * t58 + t56 * t8;
t71 = t5 * t56 - t58 * t8;
t70 = pkin(4) * t52 + t113;
t16 = t31 - t95;
t69 = t16 * t58 + t17 * t56;
t68 = t16 * t56 - t17 * t58;
t66 = t84 + t90;
t27 = t56 * t49 + t57 * t50;
t3 = t57 * t89 + t93;
t62 = -t17 * qJD(4) + t58 * t65;
t4 = -t56 * t82 + t62;
t64 = -t69 * qJD(4) - t3 * t58 - t4 * t56;
t18 = -t58 * qJD(5) - t109 * t101;
t19 = -t56 * qJD(5) + t109 * t50;
t63 = -t18 * t58 - t19 * t56 + (-t35 * t58 + t36 * t56) * qJD(4);
t61 = qJ(5) * t90 - t123 * t58 + t62;
t51 = qJ(2) * t122;
t42 = -0.2e1 * t86;
t41 = 0.2e1 * t86;
t33 = -0.2e1 * t108 * qJD(4);
t26 = t105 * t50;
t24 = t57 * t101 - t58 * t49;
t23 = t105 * t101;
t22 = qJD(3) * t78;
t15 = -0.2e1 * t54 * t83 - 0.2e1 * t75;
t14 = -0.2e1 * t52 * t83 + 0.2e1 * t75;
t12 = t108 * t100 + t76;
t11 = t108 * t98 - 0.4e1 * t59 * t86;
t10 = 0.2e1 * t56 * t77 - 0.2e1 * t57 * t88;
t9 = -0.2e1 * t106 * t102 - 0.2e1 * t57 * t90;
t7 = t108 * t55 * t121 + 0.4e1 * t59 * t76;
t6 = (-0.1e1 + t107) * t43;
t2 = (t91 + t95) * qJD(4) + t93 + t123 * t56;
t1 = t79 * t49 + t61;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, t51, -0.2e1 * t83, 0.2e1 * t77, 0, t43, 0, 0, 0.2e1 * qJD(2) * t57 + 0.2e1 * t59 * t96, 0.2e1 * qJD(2) * t59 - 0.2e1 * t57 * t96, 0, t51, t15, t7, t9, t14, t10, t43, -0.2e1 * t55 * t58 * t99 + 0.2e1 * t4 * t57 + 0.2e1 * (t16 + 0.2e1 * t95) * t49, 0.2e1 * t55 * t89 + 0.2e1 * t3 * t57 + 0.2e1 * (-t17 + 0.2e1 * t40) * t49, 0.2e1 * t69 * t98 + 0.2e1 * (t68 * qJD(4) + t3 * t56 - t4 * t58) * t59, -0.2e1 * t60 ^ 2 * t83 + 0.2e1 * t16 * t4 - 0.2e1 * t17 * t3, t15, t7, t9, t14, t10, t43, 0.2e1 * (-t56 * t103 + t1) * t57 + 0.2e1 * (qJD(3) * t5 + t32 * t50 + t117) * t59, 0.2e1 * (-t32 * t102 + t2) * t57 + 0.2e1 * (-qJD(3) * t8 + t116 - t87) * t59, 0.2e1 * t72 * t98 + 0.2e1 * (t71 * qJD(4) - t1 * t58 + t2 * t56) * t59, 0.2e1 * t5 * t1 + 0.2e1 * t32 * t13 - 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t23, 0, -t68 * t49 + (t64 - 0.2e1 * t82) * t57, 0, 0, 0, 0, 0, 0, -t26, t23, 0, (-t71 * qJD(3) - t13) * t59 + (-t72 * qJD(4) - t1 * t56 - t2 * t58 + t103) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, 0, -t49, 0, -t85, -t82, 0, 0, -t12, t11, t27, t12, -t24, 0, (-t56 * t110 - t74 * t58) * qJD(4) + (t73 * t56 - t40) * qJD(3), (-t58 * t110 + t74 * t56) * qJD(4) + (-t58 * t118 + (pkin(3) * t58 + t112) * t57) * qJD(3), t64, -pkin(3) * t85 + t64 * pkin(7), -t12, t11, t27, t12, -t24, 0, -t116 + t19 * t57 + (-t47 * t56 * t57 + t115) * qJD(3) + (t32 * t56 + t70 * t59) * qJD(4), t117 + t18 * t57 + (-t57 * t113 + t114) * qJD(3) + (t80 * t59 * t56 + t32 * t58) * qJD(4), (t35 * t98 - t19 * t59 - t2 + (-t5 + t114) * qJD(4)) * t58 + (-t36 * t98 + t18 * t59 - t1 + (-t8 + t115) * qJD(4)) * t56, pkin(4) * t87 + t1 * t35 + t13 * t47 - t8 * t18 + t5 * t19 + t2 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t49, 0, 0, 0, 0, 0, 0, 0, 0, -t66, -t28, t22, (pkin(7) * t78 - t120) * qJD(3), 0, 0, 0, 0, 0, 0, -t66, -t28, t22, (-t92 + (-t35 * t56 - t36 * t58) * qJD(3)) * t59 + (qJD(3) * t47 + t63) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t33, 0, t42, 0, 0, pkin(3) * t94, -0.2e1 * pkin(3) * t50, 0, 0, t41, t33, 0, t42, 0, 0, t80 * t94, t70 * t121, 0.2e1 * t63, 0.2e1 * t36 * t18 + 0.2e1 * t35 * t19 + 0.2e1 * t47 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, 0, -t28, t49, t4, t3, 0, 0, 0, 0, -t66, 0, -t28, t49, (0.2e1 * pkin(4) - t112) * t49 + t61, t2, t66 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t24, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t24, 0, -t27 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t101, 0, -pkin(7) * t50, pkin(7) * t101, 0, 0, 0, 0, t50, 0, -t101, 0, t19, t18, -pkin(4) * t50, t19 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t66, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t50, 0, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t20;
