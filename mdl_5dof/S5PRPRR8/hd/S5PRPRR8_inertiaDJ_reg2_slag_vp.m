% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR8_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:04:27
% EndTime: 2019-12-05 16:04:31
% DurationCPUTime: 0.90s
% Computational Cost: add. (552->121), mult. (1473->238), div. (0->0), fcn. (1298->8), ass. (0->91)
t39 = cos(qJ(5));
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t54 = pkin(4) * t40 + pkin(8) * t37;
t98 = t39 * t54;
t36 = sin(qJ(5));
t30 = t36 ^ 2;
t32 = t39 ^ 2;
t87 = t30 - t32;
t58 = qJD(5) * t87;
t97 = 2 * qJD(3);
t35 = cos(pkin(5));
t34 = sin(pkin(5));
t41 = cos(qJ(2));
t91 = t34 * t41;
t14 = t35 * t37 + t40 * t91;
t38 = sin(qJ(2));
t92 = t34 * t38;
t28 = qJD(2) * t92;
t72 = t37 * t91;
t75 = t40 * qJD(4);
t7 = -qJD(4) * t72 - t28 * t40 + t35 * t75;
t96 = t14 * t7;
t95 = t37 * pkin(4);
t94 = t40 * pkin(8);
t93 = t7 * t40;
t42 = -pkin(2) - pkin(7);
t90 = t37 * t42;
t89 = t40 * t42;
t83 = qJD(2) * t41;
t65 = t34 * t83;
t88 = qJ(3) * t65 + qJD(3) * t92;
t86 = t30 + t32;
t31 = t37 ^ 2;
t33 = t40 ^ 2;
t85 = t31 - t33;
t84 = t31 + t33;
t82 = qJD(4) * t14;
t81 = qJD(4) * t39;
t80 = qJD(5) * t36;
t79 = qJD(5) * t39;
t78 = qJD(5) * t40;
t77 = qJD(5) * t42;
t76 = t37 * qJD(4);
t74 = qJ(3) * qJD(4);
t73 = -0.2e1 * pkin(4) * qJD(5);
t71 = t36 * t90;
t70 = t36 * t89;
t69 = t39 * t90;
t68 = t36 * t78;
t67 = t36 * t77;
t66 = t39 * t78;
t64 = t36 * t79;
t63 = t42 * t76;
t62 = t39 * t76;
t61 = t37 * t75;
t60 = t42 * t75;
t59 = t86 * t40;
t57 = t85 * qJD(4);
t26 = 0.2e1 * t61;
t56 = t36 * t62;
t55 = t33 * t64;
t53 = -t94 + t95;
t15 = t35 * t40 - t72;
t8 = -t15 * t36 + t39 * t92;
t9 = t15 * t39 + t36 * t92;
t52 = t36 * t9 + t39 * t8;
t51 = t36 * t8 - t39 * t9;
t48 = qJ(3) + t53;
t45 = t39 * t48;
t11 = t45 - t71;
t12 = t36 * t48 + t69;
t50 = t11 * t39 + t12 * t36;
t49 = t11 * t36 - t12 * t39;
t47 = t14 * t79 + t36 * t7;
t46 = t14 * t80 - t39 * t7;
t6 = -t28 * t37 + t82;
t2 = -qJD(5) * t9 + t6 * t36 + t39 * t65;
t3 = qJD(5) * t8 + t36 * t65 - t6 * t39;
t44 = -qJD(5) * t52 - t2 * t36 + t3 * t39;
t4 = t37 * t67 - t36 * (qJD(4) * t54 + qJD(3)) - qJD(5) * t45 - t39 * t60;
t5 = t39 * qJD(3) - t12 * qJD(5) + (-t70 + t98) * qJD(4);
t43 = -qJD(5) * t50 - t5 * t36 - t4 * t39;
t1 = -t6 * t37 - t93 + (t14 * t37 + t15 * t40) * qJD(4);
t29 = qJ(3) * t97;
t19 = t36 * t76 - t66;
t18 = t36 * t75 + t37 * t79;
t17 = -t62 - t68;
t16 = t37 * t80 - t39 * t75;
t10 = t40 * t58 + t56;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t34 ^ 2 * t38 * t83 - 0.2e1 * t15 * t6 + 0.2e1 * t96, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t2 * t8 + 0.2e1 * t3 * t9 + 0.2e1 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t65, -pkin(2) * t28 + t88, 0, 0, 0, 0, 0, 0, (t37 * t83 + t38 * t75) * t34, (-t38 * t76 + t40 * t83) * t34, -t1, t1 * t42 + t88, 0, 0, 0, 0, 0, 0, (-t36 * t82 + t2) * t37 + (qJD(4) * t8 + t47) * t40, (-t14 * t81 - t3) * t37 + (-qJD(4) * t9 - t46) * t40, t52 * t76 + (qJD(5) * t51 - t2 * t39 - t3 * t36) * t40, t2 * t11 + t3 * t12 - t9 * t4 + t8 * t5 + (t14 * t76 - t93) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t29, -0.2e1 * t61, 0.2e1 * t57, 0, t26, 0, 0, 0.2e1 * qJD(3) * t37 + 0.2e1 * t40 * t74, 0.2e1 * qJD(3) * t40 - 0.2e1 * t37 * t74, 0, t29, -0.2e1 * t32 * t61 - 0.2e1 * t55, 0.2e1 * t33 * t58 + 0.4e1 * t40 * t56, -0.2e1 * t37 * t68 - 0.2e1 * t81 * t85, -0.2e1 * t30 * t61 + 0.2e1 * t55, 0.2e1 * t36 * t57 - 0.2e1 * t37 * t66, t26, -0.2e1 * t33 * t39 * t77 + 0.2e1 * t5 * t37 + 0.2e1 * (t11 + 0.2e1 * t71) * t75, 0.2e1 * t33 * t67 + 0.2e1 * t4 * t37 + 0.2e1 * (-t12 + 0.2e1 * t69) * t75, 0.2e1 * t50 * t76 + 0.2e1 * (qJD(5) * t49 + t36 * t4 - t39 * t5) * t40, -0.2e1 * t42 ^ 2 * t61 + 0.2e1 * t11 * t5 - 0.2e1 * t12 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-qJD(4) * t51 - t7) * t40 + (t44 + t82) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84 * t79, t84 * t80, 0, -t49 * t75 + (t43 - 0.2e1 * t60) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-0.1e1 + t86) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, 0, 0, 0, 0, 0, 0, 0, 0, t46, t47, t44, -t7 * pkin(4) + pkin(8) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, 0, -t75, 0, -t63, -t60, 0, 0, -t10, -0.4e1 * t40 * t64 + t76 * t87, t18, t10, -t16, 0, (-t70 - t98) * qJD(5) + (t36 * t53 - t69) * qJD(4), (t36 * t54 - t39 * t89) * qJD(5) + (-t39 * t94 + (pkin(4) * t39 + t36 * t42) * t37) * qJD(4), t43, -pkin(4) * t63 + pkin(8) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19, qJD(4) * t59, (pkin(8) * t59 - t95) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t64, -0.2e1 * t58, 0, -0.2e1 * t64, 0, 0, t36 * t73, t39 * t73, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, t19, t75, t5, t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, t16, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, -t80, 0, -pkin(8) * t79, pkin(8) * t80, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t13;
