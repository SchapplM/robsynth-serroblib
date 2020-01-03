% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:28
% DurationCPUTime: 0.67s
% Computational Cost: add. (848->127), mult. (2022->225), div. (0->0), fcn. (1729->6), ass. (0->89)
t56 = sin(qJ(5));
t54 = t56 ^ 2;
t59 = cos(qJ(5));
t94 = -t59 ^ 2 + t54;
t79 = t94 * qJD(5);
t105 = qJD(2) + qJD(3);
t62 = 2 * qJD(4);
t104 = pkin(3) + pkin(8);
t103 = -pkin(7) - pkin(6);
t101 = cos(qJ(3));
t58 = sin(qJ(2));
t81 = qJD(2) * t103;
t34 = t58 * t81;
t38 = t103 * t58;
t60 = cos(qJ(2));
t39 = t103 * t60;
t80 = t101 * qJD(3);
t57 = sin(qJ(3));
t93 = qJD(3) * t57;
t98 = t57 * t60;
t11 = -t101 * t34 - t38 * t80 - t39 * t93 - t81 * t98;
t32 = t101 * t58 + t98;
t22 = t105 * t32;
t6 = -t22 * pkin(4) - t11;
t4 = t6 * t56;
t5 = t6 * t59;
t82 = t101 * t60;
t99 = t57 * t58;
t31 = -t82 + t99;
t69 = -t101 * t39 + t57 * t38;
t19 = -t31 * pkin(4) + t69;
t90 = qJD(5) * t59;
t102 = t19 * t90 + t4;
t100 = t56 * t22;
t97 = t59 * t22;
t78 = pkin(2) * t80;
t41 = t78 + qJD(4);
t45 = t57 * pkin(2) + qJ(4);
t96 = t41 * t56 + t45 * t90;
t86 = qJ(4) * qJD(5);
t95 = qJD(4) * t56 + t59 * t86;
t92 = qJD(5) * t19;
t91 = qJD(5) * t56;
t89 = qJD(5) * t104;
t88 = t58 * qJD(2);
t87 = t60 * qJD(2);
t85 = -0.2e1 * pkin(1) * qJD(2);
t84 = t56 * t97;
t51 = pkin(2) * t88;
t50 = pkin(2) * t93;
t83 = t56 * t90;
t49 = -t60 * pkin(2) - pkin(1);
t77 = qJD(2) * t82;
t48 = -t101 * pkin(2) - pkin(3);
t23 = -t101 * t38 - t57 * t39;
t73 = -t32 * qJ(4) + t49;
t15 = t104 * t31 + t73;
t18 = t32 * pkin(4) + t23;
t76 = t59 * t15 + t56 * t18;
t75 = t56 * t15 - t59 * t18;
t74 = -qJ(4) * t22 - qJD(4) * t31;
t21 = t105 * t99 - t60 * t80 - t77;
t72 = -t56 * t21 + t32 * t90;
t13 = -t59 * t21 - t32 * t91;
t71 = t31 * t90 + t100;
t70 = t31 * t91 - t97;
t68 = t21 * qJ(4) - t32 * qJD(4) + t51;
t44 = -pkin(8) + t48;
t67 = qJD(5) * (t31 * t45 - t32 * t44);
t66 = qJD(5) * (qJ(4) * t31 + t104 * t32);
t65 = t104 * t21 + t74;
t64 = -t45 * t22 - t41 * t31 + t32 * t50;
t63 = -t21 * t44 + t64;
t12 = t69 * qJD(3) - t103 * t77 + t57 * t34;
t53 = qJD(4) * t59;
t40 = -0.2e1 * t83;
t36 = t41 * t59;
t30 = 0.2e1 * t79;
t29 = t31 ^ 2;
t20 = t31 * pkin(3) + t73;
t17 = -0.2e1 * t32 * t21;
t10 = -t31 * t79 + t84;
t9 = t22 * pkin(3) + t68;
t8 = -t94 * t22 - 0.4e1 * t31 * t83;
t7 = -t21 * pkin(4) + t12;
t3 = t104 * t22 + t68;
t2 = -t76 * qJD(5) - t56 * t3 + t59 * t7;
t1 = t75 * qJD(5) - t59 * t3 - t56 * t7;
t14 = [0, 0, 0, 0.2e1 * t58 * t87, 0.2e1 * (-t58 ^ 2 + t60 ^ 2) * qJD(2), 0, 0, 0, t58 * t85, t60 * t85, t17, 0.2e1 * t21 * t31 - 0.2e1 * t32 * t22, 0, 0, 0, 0.2e1 * t49 * t22 + 0.2e1 * t31 * t51, -0.2e1 * t49 * t21 + 0.2e1 * t32 * t51, 0.2e1 * t11 * t31 + 0.2e1 * t12 * t32 - 0.2e1 * t23 * t21 - 0.2e1 * t22 * t69, -0.2e1 * t20 * t22 - 0.2e1 * t9 * t31, 0.2e1 * t20 * t21 - 0.2e1 * t9 * t32, -0.2e1 * t11 * t69 + 0.2e1 * t23 * t12 + 0.2e1 * t20 * t9, 0.2e1 * t54 * t31 * t22 + 0.2e1 * t29 * t83, -0.2e1 * t29 * t79 + 0.4e1 * t31 * t84, 0.2e1 * t32 * t100 + 0.2e1 * t72 * t31, 0.2e1 * t13 * t31 + 0.2e1 * t32 * t97, t17, 0.2e1 * t70 * t19 + 0.2e1 * t2 * t32 + 0.2e1 * t75 * t21 - 0.2e1 * t31 * t5, 0.2e1 * t1 * t32 + 0.2e1 * t71 * t19 + 0.2e1 * t76 * t21 + 0.2e1 * t31 * t4; 0, 0, 0, 0, 0, t87, -t88, 0, -pkin(6) * t87, pkin(6) * t88, 0, 0, -t21, -t22, 0, -t12, t11, -t48 * t21 + t64, t12, -t11, -t11 * t45 + t12 * t48 + t23 * t50 + t41 * t69, t10, t8, t13, -t72, 0, t56 * t67 + t63 * t59 + t102, t5 + t59 * t67 + (-t63 - t92) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t50, -0.2e1 * t78, 0, 0.2e1 * t50, 0.2e1 * t41, 0.2e1 * t45 * t41 + 0.2e1 * t48 * t50, t40, t30, 0, 0, 0, 0.2e1 * t96, -0.2e1 * t45 * t91 + 0.2e1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t22, 0, -t12, t11, pkin(3) * t21 + t74, t12, -t11, -t12 * pkin(3) - t11 * qJ(4) + qJD(4) * t69, t10, t8, t13, -t72, 0, t56 * t66 + t65 * t59 + t102, t5 + t59 * t66 + (-t65 - t92) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t78, 0, t50, t62 + t78, -pkin(3) * t50 + t41 * qJ(4) + t45 * qJD(4), t40, t30, 0, 0, 0, t95 + t96, t36 + t53 + (-qJ(4) - t45) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, qJ(4) * t62, t40, t30, 0, 0, 0, 0.2e1 * t95, -0.2e1 * t56 * t86 + 0.2e1 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, t12, 0, 0, 0, 0, 0, t13, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t70, -t21, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t90, 0, -t44 * t91 + t59 * t50, -t44 * t90 - t56 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t90, 0, t56 * t89, t59 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, -t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t14;
