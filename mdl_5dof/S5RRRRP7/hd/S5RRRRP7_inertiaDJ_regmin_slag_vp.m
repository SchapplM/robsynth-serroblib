% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRRP7
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:44
% EndTime: 2019-12-31 21:57:48
% DurationCPUTime: 1.18s
% Computational Cost: add. (1842->169), mult. (4204->295), div. (0->0), fcn. (3680->6), ass. (0->110)
t67 = sin(qJ(3));
t68 = sin(qJ(2));
t70 = cos(qJ(3));
t71 = cos(qJ(2));
t48 = t67 * t68 - t70 * t71;
t49 = t67 * t71 + t70 * t68;
t61 = -t71 * pkin(2) - pkin(1);
t25 = t48 * pkin(3) - t49 * pkin(8) + t61;
t130 = -pkin(7) - pkin(6);
t54 = t130 * t68;
t55 = t130 * t71;
t32 = t67 * t54 - t70 * t55;
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t134 = t66 * t25 + t69 * t32;
t64 = t66 ^ 2;
t65 = t69 ^ 2;
t113 = t64 - t65;
t90 = t113 * qJD(4);
t133 = qJD(2) + qJD(3);
t28 = t133 * t48;
t29 = t133 * t49;
t105 = t68 * qJD(2);
t99 = pkin(2) * t105;
t13 = t29 * pkin(3) + t28 * pkin(8) + t99;
t110 = qJD(3) * t70;
t111 = qJD(3) * t67;
t93 = qJD(2) * t130;
t88 = t71 * t93;
t89 = t68 * t93;
t18 = -t54 * t110 - t55 * t111 - t67 * t88 - t70 * t89;
t5 = -qJD(4) * t134 + t69 * t13 + t66 * t18;
t85 = t69 * pkin(4) + t66 * qJ(5);
t132 = t85 * qJD(4) - t69 * qJD(5);
t131 = 0.2e1 * qJD(5);
t129 = pkin(8) * t29;
t128 = pkin(8) * t48;
t127 = t29 * pkin(4);
t126 = t70 * pkin(2);
t59 = t67 * pkin(2) + pkin(8);
t125 = t29 * t59;
t124 = t48 * t59;
t123 = t49 * t28;
t122 = t49 * t66;
t121 = t49 * t69;
t120 = t66 * t29;
t119 = t69 * t28;
t118 = t69 * t29;
t19 = t32 * qJD(3) + t67 * t89 - t70 * t88;
t31 = -t70 * t54 - t67 * t55;
t63 = qJD(4) * t69;
t117 = t19 * t66 + t31 * t63;
t106 = t66 * qJD(5);
t108 = qJD(4) * t66;
t39 = pkin(4) * t108 - qJ(5) * t63 - t106;
t98 = pkin(2) * t111;
t33 = t39 + t98;
t115 = -t33 - t39;
t60 = -pkin(3) - t126;
t114 = t60 * t63 + t66 * t98;
t112 = t29 * qJ(5);
t84 = t66 * pkin(4) - t69 * qJ(5);
t17 = t84 * t49 + t31;
t109 = qJD(4) * t17;
t107 = t48 * qJD(5);
t103 = t71 * qJD(2);
t102 = -0.2e1 * pkin(1) * qJD(2);
t101 = pkin(3) * t108;
t100 = pkin(3) * t63;
t97 = pkin(2) * t110;
t96 = pkin(8) * t108;
t95 = pkin(8) * t63;
t94 = t66 * t63;
t92 = -0.4e1 * t66 * t121;
t8 = t48 * qJ(5) + t134;
t83 = t69 * t25 - t66 * t32;
t9 = -t48 * pkin(4) - t83;
t87 = t66 * t9 + t69 * t8;
t86 = -t66 * t8 + t69 * t9;
t81 = -t49 * t60 + t124;
t37 = (t64 + t65) * t97;
t80 = t60 * t108 - t69 * t98;
t53 = -pkin(3) - t85;
t79 = -t66 * t28 + t49 * t63;
t78 = t49 * t108 + t119;
t77 = t48 * t108 - t118;
t4 = t32 * t108 - t66 * t13 + t69 * t18 - t25 * t63;
t76 = -t28 * t53 + t39 * t49 - t129;
t6 = t132 * t49 - t84 * t28 + t19;
t75 = -t6 + (t49 * t53 - t128) * qJD(4);
t44 = t53 - t126;
t74 = -t6 + (t44 * t49 - t124) * qJD(4);
t2 = t107 - t4 + t112;
t3 = -t127 - t5;
t1 = t86 * qJD(4) + t2 * t69 + t3 * t66;
t73 = -t28 * t44 + t33 * t49 - t48 * t97 - t125;
t72 = -t28 * t60 - t125 + (-t48 * t70 + t49 * t67) * qJD(3) * pkin(2);
t57 = 0.2e1 * t94;
t47 = -0.2e1 * t90;
t46 = t49 ^ 2;
t43 = t53 * t108;
t36 = t44 * t108;
t35 = t59 * t63 + t66 * t97;
t34 = t59 * t108 - t69 * t97;
t26 = t31 * t108;
t22 = t48 * t63 + t120;
t14 = t17 * t108;
t12 = -t66 * t119 - t49 * t90;
t7 = qJD(4) * t92 + t113 * t28;
t10 = [0, 0, 0, 0.2e1 * t68 * t103, 0.2e1 * (-t68 ^ 2 + t71 ^ 2) * qJD(2), 0, 0, 0, t68 * t102, t71 * t102, -0.2e1 * t123, 0.2e1 * t28 * t48 - 0.2e1 * t49 * t29, 0, 0, 0, 0.2e1 * t61 * t29 + 0.2e1 * t48 * t99, -0.2e1 * t61 * t28 + 0.2e1 * t49 * t99, -0.2e1 * t65 * t123 - 0.2e1 * t46 * t94, -t28 * t92 + 0.2e1 * t46 * t90, 0.2e1 * t49 * t118 - 0.2e1 * t78 * t48, -0.2e1 * t49 * t120 - 0.2e1 * t79 * t48, 0.2e1 * t48 * t29, 0.2e1 * t19 * t122 + 0.2e1 * t83 * t29 + 0.2e1 * t79 * t31 + 0.2e1 * t5 * t48, 0.2e1 * t19 * t121 - 0.2e1 * t134 * t29 - 0.2e1 * t78 * t31 + 0.2e1 * t4 * t48, 0.2e1 * t6 * t122 + 0.2e1 * t79 * t17 - 0.2e1 * t9 * t29 - 0.2e1 * t3 * t48, -0.2e1 * t86 * t28 + 0.2e1 * (-t87 * qJD(4) - t2 * t66 + t3 * t69) * t49, -0.2e1 * t6 * t121 + 0.2e1 * t78 * t17 + 0.2e1 * t2 * t48 + 0.2e1 * t8 * t29, 0.2e1 * t17 * t6 + 0.2e1 * t8 * t2 + 0.2e1 * t9 * t3; 0, 0, 0, 0, 0, t103, -t105, 0, -pkin(6) * t103, pkin(6) * t105, 0, 0, -t28, -t29, 0, -t19, t18, t12, t7, t22, -t77, 0, t26 + (-t81 * qJD(4) - t19) * t69 + t72 * t66, t81 * t108 + t72 * t69 + t117, t73 * t66 + t74 * t69 + t14, t1, t74 * t66 + (-t73 - t109) * t69, t1 * t59 + t17 * t33 + t6 * t44 + t87 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t98, -0.2e1 * t97, t57, t47, 0, 0, 0, 0.2e1 * t80, 0.2e1 * t114, -0.2e1 * t33 * t69 + 0.2e1 * t36, 0.2e1 * t37, -0.2e1 * t33 * t66 - 0.2e1 * t44 * t63, 0.2e1 * t44 * t33 + 0.2e1 * t59 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, -t19, t18, t12, t7, t22, -t77, 0, t26 + (pkin(3) * t28 - t129) * t66 + (-t19 + (-pkin(3) * t49 - t128) * qJD(4)) * t69, t78 * pkin(3) + t77 * pkin(8) + t117, t76 * t66 + t75 * t69 + t14, t1, t75 * t66 + (-t76 - t109) * t69, t1 * pkin(8) + t17 * t39 + t6 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t98, -t97, t57, t47, 0, 0, 0, t80 - t101, -t100 + t114, t115 * t69 + t36 + t43, t37, t115 * t66 + (-t44 - t53) * t63, pkin(8) * t37 + t33 * t53 + t44 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t47, 0, 0, 0, -0.2e1 * t101, -0.2e1 * t100, -0.2e1 * t39 * t69 + 0.2e1 * t43, 0, -0.2e1 * t39 * t66 - 0.2e1 * t53 * t63, 0.2e1 * t53 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, -t79, t29, t5, t4, t5 + 0.2e1 * t127, t85 * t28 + (t84 * qJD(4) - t106) * t49, 0.2e1 * t107 - t4 + 0.2e1 * t112, -t3 * pkin(4) + t2 * qJ(5) + t8 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t108, 0, -t35, t34, -t35, -t132, -t34, -t132 * t59 - t84 * t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t108, 0, -t95, t96, -t95, -t132, -t96, -t132 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, qJ(5) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29, -t78, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t10;
