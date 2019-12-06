% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRPR6_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaDJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:47
% EndTime: 2019-12-05 16:32:55
% DurationCPUTime: 1.44s
% Computational Cost: add. (1246->192), mult. (3486->380), div. (0->0), fcn. (3276->10), ass. (0->107)
t65 = sin(qJ(5));
t108 = qJD(5) * t65;
t61 = sin(pkin(10));
t63 = cos(pkin(10));
t122 = cos(qJ(5));
t89 = qJD(5) * t122;
t35 = t61 * t108 - t63 * t89;
t91 = t122 * t63;
t42 = t65 * t61 - t91;
t66 = sin(qJ(3));
t68 = cos(qJ(3));
t127 = (t66 ^ 2 - t68 ^ 2) * qJD(3);
t126 = -0.2e1 * t35;
t43 = t122 * t61 + t65 * t63;
t36 = t43 * qJD(5);
t125 = 0.2e1 * t36;
t124 = pkin(3) * t66;
t123 = t68 * pkin(3);
t62 = sin(pkin(5));
t67 = sin(qJ(2));
t118 = t62 * t67;
t64 = cos(pkin(5));
t38 = t68 * t118 + t64 * t66;
t69 = cos(qJ(2));
t110 = qJD(2) * t69;
t96 = t62 * t110;
t27 = t38 * qJD(3) + t66 * t96;
t121 = t27 * t66;
t37 = t66 * t118 - t64 * t68;
t14 = t37 * t27;
t120 = t61 * t66;
t119 = t61 * t68;
t117 = t62 * t69;
t116 = t63 * t66;
t115 = t63 * t68;
t113 = pkin(8) + qJ(4);
t52 = pkin(7) * t115;
t82 = -t66 * qJ(4) - t123;
t77 = -pkin(2) + t82;
t32 = t61 * t77 + t52;
t111 = qJD(2) * t67;
t109 = qJD(4) * t61;
t107 = t63 * qJD(4);
t106 = t66 * qJD(3);
t105 = t66 * qJD(4);
t104 = t68 * qJD(3);
t103 = pkin(7) * t119;
t102 = pkin(7) * t116;
t101 = -0.2e1 * pkin(2) * qJD(3);
t56 = pkin(7) * t104;
t99 = t61 * t104;
t98 = t61 * t105;
t97 = t62 * t111;
t95 = t63 * t104;
t94 = t63 * t105;
t93 = t66 * t104;
t92 = t61 * pkin(7) + pkin(4);
t90 = t113 * t61;
t88 = t122 * qJD(4);
t87 = 0.2e1 * t93;
t86 = t61 * t95;
t57 = t61 ^ 2;
t58 = t63 ^ 2;
t85 = 0.2e1 * (t57 + t58) * qJD(4);
t84 = -0.2e1 * t127;
t81 = -qJ(4) * t68 + t124;
t26 = -t37 * qJD(3) + t68 * t96;
t15 = -t26 * t61 + t63 * t97;
t16 = t26 * t63 + t61 * t97;
t80 = -t15 * t61 + t16 * t63;
t22 = -t94 + (pkin(7) * t120 + t63 * t81) * qJD(3);
t23 = -t98 + (t61 * t81 - t102) * qJD(3);
t79 = -t22 * t61 + t23 * t63;
t78 = t122 * t90;
t76 = -t113 * t68 + t124;
t24 = -t63 * t117 - t38 * t61;
t25 = -t61 * t117 + t38 * t63;
t8 = t122 * t25 + t65 * t24;
t75 = t26 * t68 + t121 + (t37 * t68 - t38 * t66) * qJD(3);
t74 = -t92 * t68 + (-t113 * t66 - pkin(2) - t123) * t63;
t73 = t65 * t74;
t72 = t122 * t74;
t71 = t98 - (t76 * t61 - t102) * qJD(3);
t70 = -t94 + (t76 * t63 + t92 * t66) * qJD(3);
t55 = -t63 * pkin(4) - pkin(3);
t49 = -0.2e1 * t93;
t47 = t113 * t63;
t44 = (pkin(4) * t61 + pkin(7)) * t66;
t39 = pkin(4) * t99 + t56;
t34 = t42 * t66;
t33 = t43 * t66;
t31 = t63 * t77 - t103;
t30 = t122 * t47 - t65 * t90;
t29 = -t65 * t47 - t78;
t28 = -pkin(8) * t120 + t32;
t18 = t43 * t104 - t35 * t66;
t17 = -t91 * t104 + t66 * t36 + t65 * t99;
t13 = -t47 * t89 - t65 * t107 + (t113 * t108 - t88) * t61;
t12 = qJD(5) * t78 - t63 * t88 + (qJD(5) * t47 + t109) * t65;
t7 = t122 * t24 - t65 * t25;
t6 = t122 * t28 + t73;
t5 = -t65 * t28 + t72;
t4 = -qJD(5) * t73 + t122 * t70 - t28 * t89 + t65 * t71;
t3 = -qJD(5) * t72 + t28 * t108 + t122 * t71 - t65 * t70;
t2 = -t8 * qJD(5) + t122 * t15 - t65 * t16;
t1 = t25 * t108 - t122 * t16 - t65 * t15 - t24 * t89;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t62 ^ 2 * t67 * t110 + 0.2e1 * t38 * t26 + 0.2e1 * t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t24 * t15 + 0.2e1 * t25 * t16 + 0.2e1 * t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t8 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, -t96, 0, 0, 0, 0, 0, 0, 0, 0, (-t69 * t106 - t68 * t111) * t62, (-t69 * t104 + t66 * t111) * t62, t75, -pkin(2) * t97 + t75 * pkin(7), 0, 0, 0, 0, 0, 0, t27 * t120 - t15 * t68 + (t37 * t119 + t24 * t66) * qJD(3), t27 * t116 + t16 * t68 + (t37 * t115 - t25 * t66) * qJD(3), (-t15 * t63 - t16 * t61) * t66 + (-t24 * t63 - t25 * t61) * t104, t15 * t31 + t16 * t32 + t24 * t22 + t25 * t23 + (t37 * t104 + t121) * pkin(7), 0, 0, 0, 0, 0, 0, t7 * t106 + t37 * t18 - t2 * t68 + t27 * t33, -t1 * t68 - t8 * t106 - t37 * t17 - t27 * t34, t1 * t33 + t7 * t17 - t8 * t18 + t2 * t34, -t1 * t6 + t2 * t5 + t27 * t44 - t8 * t3 + t37 * t39 + t7 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t84, 0, t49, 0, 0, t66 * t101, t68 * t101, 0, 0, t58 * t87, -0.4e1 * t66 * t86, 0.2e1 * t63 * t127, t57 * t87, t61 * t84, t49, -0.2e1 * t22 * t68 + 0.2e1 * (t31 + 0.2e1 * t103) * t106, 0.2e1 * t23 * t68 + 0.2e1 * (-t32 + 0.2e1 * t52) * t106, 0.2e1 * (-t22 * t63 - t23 * t61) * t66 + 0.2e1 * (-t31 * t63 - t32 * t61) * t104, 0.2e1 * pkin(7) ^ 2 * t93 + 0.2e1 * t31 * t22 + 0.2e1 * t32 * t23, 0.2e1 * t34 * t17, 0.2e1 * t17 * t33 + 0.2e1 * t34 * t18, -0.2e1 * t34 * t106 + 0.2e1 * t17 * t68, 0.2e1 * t33 * t18, -0.2e1 * t33 * t106 + 0.2e1 * t18 * t68, t49, 0.2e1 * t5 * t106 + 0.2e1 * t44 * t18 + 0.2e1 * t39 * t33 - 0.2e1 * t4 * t68, -0.2e1 * t6 * t106 - 0.2e1 * t44 * t17 - 0.2e1 * t3 * t68 - 0.2e1 * t39 * t34, 0.2e1 * t5 * t17 - 0.2e1 * t6 * t18 + 0.2e1 * t3 * t33 + 0.2e1 * t4 * t34, -0.2e1 * t6 * t3 + 0.2e1 * t44 * t39 + 0.2e1 * t5 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, 0, 0, 0, 0, 0, 0, 0, 0, -t27 * t63, t27 * t61, t80, -t27 * pkin(3) + (-t24 * t61 + t25 * t63) * qJD(4) + t80 * qJ(4), 0, 0, 0, 0, 0, 0, t27 * t42 + t37 * t36, t27 * t43 - t37 * t35, t1 * t42 - t2 * t43 + t7 * t35 - t8 * t36, -t1 * t30 - t8 * t12 + t7 * t13 + t2 * t29 + t27 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, -t106, 0, -t56, pkin(7) * t106, 0, 0, t86, (-t57 + t58) * t104, t61 * t106, -t86, t63 * t106, 0, t68 * t109 + (t82 * t61 - t52) * qJD(3), t68 * t107 + (t82 * t63 + t103) * qJD(3), t79, -pkin(3) * t56 + (-t31 * t61 + t32 * t63) * qJD(4) + t79 * qJ(4), -t17 * t43 + t34 * t35, t17 * t42 - t43 * t18 + t35 * t33 + t34 * t36, t43 * t106 + t35 * t68, t18 * t42 + t33 * t36, -t42 * t106 + t36 * t68, 0, t29 * t106 - t13 * t68 + t55 * t18 + t44 * t36 + t39 * t42, -t30 * t106 - t12 * t68 - t55 * t17 - t44 * t35 + t39 * t43, t12 * t33 + t13 * t34 + t29 * t17 - t30 * t18 + t3 * t42 + t5 * t35 - t6 * t36 - t4 * t43, -t6 * t12 + t5 * t13 + t4 * t29 - t3 * t30 + t39 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, qJ(4) * t85, t43 * t126, 0.2e1 * t35 * t42 - 0.2e1 * t43 * t36, 0, t42 * t125, 0, 0, t55 * t125, t55 * t126, 0.2e1 * t12 * t42 - 0.2e1 * t13 * t43 + 0.2e1 * t29 * t35 - 0.2e1 * t30 * t36, -0.2e1 * t30 * t12 + 0.2e1 * t29 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t95, 0, t56, 0, 0, 0, 0, 0, 0, t18, -t17, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, -t18, t106, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, -t36, 0, t13, t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
