% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP10_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:52:04
% EndTime: 2019-12-31 18:52:08
% DurationCPUTime: 1.21s
% Computational Cost: add. (1529->156), mult. (3584->288), div. (0->0), fcn. (3316->6), ass. (0->98)
t65 = sin(qJ(3));
t102 = qJD(3) * t65;
t107 = pkin(6) + qJ(2);
t63 = cos(pkin(8));
t47 = t107 * t63;
t62 = sin(pkin(8));
t116 = cos(qJ(3));
t79 = t107 * t116;
t88 = t116 * qJD(2);
t99 = t65 * qJD(2);
t68 = (-qJD(3) * t79 - t99) * t62 - t47 * t102 + t63 * t88;
t44 = -t116 * t63 + t65 * t62;
t45 = t116 * t62 + t65 * t63;
t54 = -t63 * pkin(2) - pkin(1);
t74 = t44 * pkin(3) - t45 * pkin(7) + t54;
t123 = -qJD(4) * t74 - t68;
t64 = sin(qJ(4));
t60 = t64 ^ 2;
t66 = cos(qJ(4));
t61 = t66 ^ 2;
t105 = t60 - t61;
t90 = t65 * t107;
t34 = t116 * t47 - t62 * t90;
t11 = -t64 * t34 + t66 * t74;
t12 = t66 * t34 + t64 * t74;
t100 = qJD(4) * t64;
t89 = qJD(3) * t116;
t40 = t62 * t102 - t63 * t89;
t41 = t45 * qJD(3);
t85 = t41 * pkin(3) + t40 * pkin(7);
t95 = t123 * t66 - t64 * t85;
t3 = t34 * t100 + t95;
t57 = qJD(4) * t66;
t4 = t123 * t64 - t34 * t57 + t66 * t85;
t122 = t3 * t64 - t4 * t66 + (t11 * t64 - t12 * t66) * qJD(4);
t118 = t41 * pkin(4);
t103 = t40 * qJ(5);
t101 = qJD(4) * t45;
t92 = qJ(5) * t101;
t98 = t66 * qJD(5);
t67 = t66 * t103 - t45 * t98 + t64 * t92 + t4;
t1 = t67 + t118;
t2 = t66 * t92 + (qJD(4) * t34 + t45 * qJD(5) - t103) * t64 + t95;
t104 = qJ(5) * t45;
t5 = t44 * pkin(4) - t66 * t104 + t11;
t8 = -t64 * t104 + t12;
t121 = -t1 * t66 + t2 * t64 + (t5 * t64 - t66 * t8) * qJD(4);
t120 = 0.2e1 * t41;
t119 = 0.2e1 * qJD(4);
t117 = t66 * pkin(4);
t18 = t47 * t89 + t63 * t99 + (-qJD(3) * t90 + t88) * t62;
t33 = t65 * t47 + t62 * t79;
t115 = t33 * t18;
t114 = t45 * t40;
t113 = t45 * t64;
t112 = t45 * t66;
t111 = t64 * t40;
t110 = t64 * t41;
t109 = t66 * t40;
t108 = t66 * t41;
t106 = -qJ(5) - pkin(7);
t32 = t44 * t120;
t97 = t64 * t109;
t96 = -0.2e1 * t100;
t94 = pkin(4) * t100;
t93 = t64 * t57;
t55 = -pkin(3) - t117;
t91 = -t55 + t117;
t42 = t45 ^ 2;
t87 = t42 * t93;
t86 = 0.2e1 * (t62 ^ 2 + t63 ^ 2) * qJD(2);
t84 = pkin(3) * t40 - pkin(7) * t41;
t83 = pkin(3) * t45 + pkin(7) * t44;
t78 = pkin(4) * t60 + t55 * t66;
t77 = t11 * t66 + t12 * t64;
t30 = t45 * t57 - t111;
t75 = t45 * t100 + t109;
t28 = t44 * t57 + t110;
t69 = -t77 * qJD(4) - t3 * t66 - t4 * t64;
t51 = -0.2e1 * t93;
t50 = 0.2e1 * t93;
t49 = t106 * t66;
t48 = t106 * t64;
t46 = -0.2e1 * t105 * qJD(4);
t39 = -t64 * qJD(5) + t106 * t57;
t38 = -t106 * t100 - t98;
t26 = -t44 * t100 + t108;
t21 = (t60 + t61) * t40;
t19 = pkin(4) * t113 + t33;
t16 = -0.2e1 * t61 * t114 - 0.2e1 * t87;
t15 = -0.2e1 * t60 * t114 + 0.2e1 * t87;
t14 = t105 * t101 + t97;
t13 = t105 * t40 - 0.4e1 * t45 * t93;
t10 = t105 * t42 * t119 + 0.4e1 * t45 * t97;
t9 = t30 * pkin(4) + t18;
t7 = -0.2e1 * t45 * t110 - 0.2e1 * t30 * t44;
t6 = 0.2e1 * t45 * t108 - 0.2e1 * t75 * t44;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, qJ(2) * t86, -0.2e1 * t114, 0.2e1 * t40 * t44 - 0.2e1 * t45 * t41, 0, t32, 0, 0, t54 * t120, -0.2e1 * t54 * t40, 0.2e1 * t18 * t45 - 0.2e1 * t33 * t40 - 0.2e1 * t34 * t41 - 0.2e1 * t68 * t44, 0.2e1 * t68 * t34 + 0.2e1 * t115, t16, t10, t6, t15, t7, t32, 0.2e1 * t11 * t41 + 0.2e1 * t18 * t113 + 0.2e1 * t30 * t33 + 0.2e1 * t4 * t44, 0.2e1 * t18 * t112 - 0.2e1 * t12 * t41 + 0.2e1 * t3 * t44 - 0.2e1 * t75 * t33, 0.2e1 * t122 * t45 + 0.2e1 * t77 * t40, 0.2e1 * t11 * t4 - 0.2e1 * t12 * t3 + 0.2e1 * t115, t16, t10, t6, t15, t7, t32, 0.2e1 * t1 * t44 + 0.2e1 * t9 * t113 + 0.2e1 * t30 * t19 + 0.2e1 * t5 * t41, 0.2e1 * t9 * t112 - 0.2e1 * t75 * t19 + 0.2e1 * t2 * t44 - 0.2e1 * t8 * t41, 0.2e1 * (t5 * t66 + t64 * t8) * t40 + 0.2e1 * t121 * t45, 0.2e1 * t5 * t1 + 0.2e1 * t19 * t9 - 0.2e1 * t8 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t28, t21, -t122, 0, 0, 0, 0, 0, 0, t26, -t28, t21, -t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, -t41, 0, -t18, -t68, 0, 0, -t14, t13, t28, t14, t26, 0, -t18 * t66 + t84 * t64 + (t33 * t64 - t83 * t66) * qJD(4), t18 * t64 + t84 * t66 + (t33 * t66 + t83 * t64) * qJD(4), t69, -t18 * pkin(3) + t69 * pkin(7), -t14, t13, t28, t14, t26, 0, -t55 * t111 + t39 * t44 + t48 * t41 - t9 * t66 + (t19 * t64 + t78 * t45) * qJD(4), -t55 * t109 + t38 * t44 + t49 * t41 + t9 * t64 + (t91 * t113 + t19 * t66) * qJD(4), (-t39 * t45 + t40 * t48 - t2 + (t45 * t49 - t5) * qJD(4)) * t66 + (t38 * t45 - t40 * t49 - t1 + (t45 * t48 - t8) * qJD(4)) * t64, t1 * t48 + t19 * t94 + t2 * t49 - t8 * t38 + t5 * t39 + t9 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 * t38 + t66 * t39 + (-t48 * t64 - t49 * t66) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t46, 0, t51, 0, 0, pkin(3) * t96, -0.2e1 * pkin(3) * t57, 0, 0, t50, t46, 0, t51, 0, 0, t91 * t96, t78 * t119, -0.2e1 * t38 * t66 - 0.2e1 * t39 * t64 + 0.2e1 * (-t48 * t66 + t49 * t64) * qJD(4), 0.2e1 * t49 * t38 + 0.2e1 * t48 * t39 + 0.2e1 * t55 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75, 0, -t30, t41, t4, t3, 0, 0, 0, 0, -t75, 0, -t30, t41, t67 + 0.2e1 * t118, t2, t75 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t57, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t57, 0, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, -t100, 0, -pkin(7) * t57, pkin(7) * t100, 0, 0, 0, 0, t57, 0, -t100, 0, t39, t38, -pkin(4) * t57, t39 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t75, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t57, 0, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t17;
