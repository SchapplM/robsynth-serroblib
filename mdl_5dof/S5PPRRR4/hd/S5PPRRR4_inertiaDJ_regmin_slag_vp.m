% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PPRRR4_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_inertiaDJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:43
% EndTime: 2019-12-05 15:19:46
% DurationCPUTime: 0.54s
% Computational Cost: add. (303->91), mult. (1060->196), div. (0->0), fcn. (1073->12), ass. (0->80)
t39 = cos(qJ(5));
t28 = t39 ^ 2;
t36 = sin(qJ(5));
t80 = t36 ^ 2 - t28;
t56 = t80 * qJD(5);
t85 = pkin(8) * t36;
t31 = sin(pkin(6));
t38 = sin(qJ(3));
t84 = t31 * t38;
t41 = cos(qJ(3));
t83 = t31 * t41;
t33 = cos(pkin(11));
t34 = cos(pkin(6));
t82 = t33 * t34;
t40 = cos(qJ(4));
t81 = t39 * t40;
t37 = sin(qJ(4));
t27 = t37 ^ 2;
t79 = -t40 ^ 2 + t27;
t78 = qJD(3) * t38;
t77 = qJD(4) * t36;
t76 = qJD(4) * t39;
t75 = qJD(5) * t36;
t74 = qJD(5) * t39;
t73 = qJD(5) * t40;
t66 = t37 * t84;
t18 = -t40 * t34 + t66;
t72 = t18 * qJD(5);
t71 = t37 * qJD(4);
t70 = t40 * qJD(4);
t69 = -0.2e1 * pkin(3) * qJD(4);
t68 = -0.2e1 * pkin(4) * qJD(5);
t67 = t41 * t82;
t65 = t36 * t73;
t64 = t39 * t73;
t63 = t31 * t78;
t62 = qJD(3) * t83;
t61 = t36 * t71;
t60 = t36 * t74;
t59 = t37 * t70;
t58 = t39 * t71;
t57 = t39 * t70;
t55 = t79 * qJD(4);
t54 = t37 * t57;
t53 = -t40 * pkin(4) - t37 * pkin(9);
t52 = pkin(4) * t37 - pkin(9) * t40;
t30 = sin(pkin(11));
t32 = sin(pkin(5));
t35 = cos(pkin(5));
t13 = -t35 * t83 + (t30 * t38 - t67) * t32;
t14 = t35 * t84 + (t30 * t41 + t38 * t82) * t32;
t17 = -t32 * t33 * t31 + t35 * t34;
t8 = t14 * t40 + t17 * t37;
t51 = t13 * t39 - t8 * t36;
t50 = t13 * t36 + t8 * t39;
t7 = t14 * t37 - t17 * t40;
t11 = -t35 * t62 + (-qJD(3) * t67 + t30 * t78) * t32;
t3 = t8 * qJD(4) - t11 * t37;
t49 = t3 * t36 + t7 * t74;
t48 = -t3 * t39 + t7 * t75;
t19 = t37 * t34 + t40 * t84;
t47 = t36 * t19 + t39 * t83;
t46 = -t39 * t19 + t36 * t83;
t16 = t19 * qJD(4) + t37 * t62;
t45 = t16 * t36 + t39 * t72;
t44 = -t16 * t39 + t36 * t72;
t43 = t58 + t65;
t42 = t61 - t64;
t24 = -pkin(3) + t53;
t21 = t52 * qJD(4);
t15 = qJD(4) * t66 - t34 * t70 - t40 * t62;
t12 = t14 * qJD(3);
t10 = t42 * pkin(8) + t39 * t21 - t24 * t75;
t9 = t43 * pkin(8) - t36 * t21 - t24 * t74;
t6 = t46 * qJD(5) + t36 * t15 + t39 * t63;
t5 = t47 * qJD(5) + t39 * t15 - t36 * t63;
t4 = -qJD(4) * t7 - t11 * t40;
t2 = t51 * qJD(5) + t12 * t36 + t4 * t39;
t1 = -t50 * qJD(5) + t12 * t39 - t4 * t36;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, -t12 * t40 + t13 * t71, t12 * t37 + t13 * t70, 0, 0, 0, 0, 0, (t7 * t77 - t1) * t40 + (t51 * qJD(4) + t49) * t37, (t7 * t76 + t2) * t40 + (-t50 * qJD(4) - t48) * t37; 0, 0, 0, -t63, -t62, 0, 0, 0, 0, 0, (-t40 * t78 - t41 * t71) * t31, (t37 * t78 - t41 * t70) * t31, 0, 0, 0, 0, 0, (t18 * t77 - t6) * t40 + (-t47 * qJD(4) + t45) * t37, (t18 * t76 - t5) * t40 + (t46 * qJD(4) - t44) * t37; 0, 0, 0, 0, 0, 0.2e1 * t59, -0.2e1 * t55, 0, 0, 0, t37 * t69, t40 * t69, -0.2e1 * t27 * t60 + 0.2e1 * t28 * t59, 0.2e1 * t27 * t56 - 0.4e1 * t36 * t54, 0.2e1 * t37 * t65 + 0.2e1 * t79 * t76, -0.2e1 * t36 * t55 + 0.2e1 * t37 * t64, -0.2e1 * t59, 0.2e1 * t24 * t58 - 0.2e1 * t10 * t40 + 0.2e1 * (t27 * t74 + t36 * t59) * pkin(8), -0.2e1 * t24 * t61 - 0.2e1 * t9 * t40 + 0.2e1 * (-t27 * t75 + t54) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0, 0, 0, 0, t48, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, 0, 0, 0, 0, 0, t44, t45; 0, 0, 0, 0, 0, 0, 0, t70, -t71, 0, -pkin(8) * t70, pkin(8) * t71, t36 * t57 - t37 * t56, -0.4e1 * t37 * t60 - t80 * t70, t42, t43, 0, (pkin(9) * t81 + (-pkin(4) * t39 + t85) * t37) * qJD(5) + (-pkin(8) * t81 + t53 * t36) * qJD(4), (pkin(8) * t37 * t39 + t52 * t36) * qJD(5) + (t53 * t39 + t40 * t85) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t60, -0.2e1 * t56, 0, 0, 0, t36 * t68, t39 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 * t75 + t57, -t36 * t70 - t37 * t74, t71, t10, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t75, 0, -pkin(9) * t74, pkin(9) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t20;
