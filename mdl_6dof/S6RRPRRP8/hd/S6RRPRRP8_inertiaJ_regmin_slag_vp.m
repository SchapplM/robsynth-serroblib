% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRP8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:22:47
% EndTime: 2019-05-06 18:22:50
% DurationCPUTime: 0.87s
% Computational Cost: add. (1535->137), mult. (3102->239), div. (0->0), fcn. (3499->8), ass. (0->72)
t53 = sin(pkin(10));
t54 = cos(pkin(10));
t56 = sin(qJ(4));
t58 = cos(qJ(4));
t37 = t56 * t53 - t58 * t54;
t38 = t58 * t53 + t56 * t54;
t55 = sin(qJ(5));
t76 = cos(qJ(5));
t22 = -t55 * t37 + t76 * t38;
t85 = -0.2e1 * t22;
t44 = -t54 * pkin(3) - pkin(2);
t31 = t37 * pkin(4) + t44;
t84 = 0.2e1 * t31;
t83 = 0.2e1 * t44;
t59 = cos(qJ(2));
t82 = -0.2e1 * t59;
t81 = 0.2e1 * t59;
t57 = sin(qJ(2));
t40 = -t59 * pkin(2) - t57 * qJ(3) - pkin(1);
t36 = t54 * t40;
t72 = t54 * t57;
t80 = pkin(7) * t53;
t23 = -pkin(8) * t72 + t36 + (-pkin(3) - t80) * t59;
t77 = t59 * pkin(7);
t29 = t53 * t40 + t54 * t77;
t73 = t53 * t57;
t25 = -pkin(8) * t73 + t29;
t15 = t56 * t23 + t58 * t25;
t32 = t38 * t57;
t13 = -t32 * pkin(9) + t15;
t14 = t58 * t23 - t56 * t25;
t33 = t37 * t57;
t79 = t59 * pkin(4);
t8 = t33 * pkin(9) + t14 - t79;
t4 = t76 * t13 + t55 * t8;
t48 = t57 * pkin(7);
t78 = t59 * pkin(5);
t71 = pkin(8) + qJ(3);
t66 = t71 * t54;
t67 = t71 * t53;
t27 = -t56 * t67 + t58 * t66;
t19 = -t37 * pkin(9) + t27;
t26 = -t56 * t66 - t58 * t67;
t63 = -t38 * pkin(9) + t26;
t11 = t55 * t19 - t76 * t63;
t75 = t11 * t59;
t12 = t76 * t19 + t55 * t63;
t74 = t12 * t59;
t39 = pkin(3) * t73 + t48;
t70 = t53 ^ 2 + t54 ^ 2;
t69 = t59 * qJ(6);
t68 = t76 * pkin(4);
t3 = -t55 * t13 + t76 * t8;
t24 = t32 * pkin(4) + t39;
t65 = -pkin(2) * t57 + qJ(3) * t59;
t28 = -t53 * t77 + t36;
t64 = -t28 * t53 + t29 * t54;
t62 = 0.2e1 * pkin(5);
t60 = 0.2e1 * qJ(6);
t52 = t59 ^ 2;
t51 = t57 ^ 2;
t47 = t55 * pkin(4);
t45 = t68 + pkin(5);
t43 = t47 + qJ(6);
t21 = t76 * t37 + t55 * t38;
t18 = -t55 * t32 - t76 * t33;
t17 = t76 * t32 - t55 * t33;
t10 = t21 * pkin(5) - t22 * qJ(6) + t31;
t5 = t17 * pkin(5) - t18 * qJ(6) + t24;
t2 = -t3 + t78;
t1 = -t69 + t4;
t6 = [1, 0, 0, t51, t57 * t81, 0, 0, 0, pkin(1) * t81, -0.2e1 * pkin(1) * t57, -0.2e1 * t28 * t59 + 0.2e1 * t51 * t80, 0.2e1 * t51 * pkin(7) * t54 + 0.2e1 * t29 * t59, 0.2e1 * (-t28 * t54 - t29 * t53) * t57, t51 * pkin(7) ^ 2 + t28 ^ 2 + t29 ^ 2, t33 ^ 2, 0.2e1 * t33 * t32, -t33 * t82, t32 * t81, t52, -0.2e1 * t14 * t59 + 0.2e1 * t39 * t32, 0.2e1 * t15 * t59 - 0.2e1 * t39 * t33, t18 ^ 2, -0.2e1 * t18 * t17, t18 * t82, t17 * t81, t52, 0.2e1 * t24 * t17 - 0.2e1 * t3 * t59, 0.2e1 * t24 * t18 + 0.2e1 * t4 * t59, 0.2e1 * t5 * t17 + 0.2e1 * t2 * t59, -0.2e1 * t1 * t17 + 0.2e1 * t2 * t18, -0.2e1 * t1 * t59 - 0.2e1 * t5 * t18, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t57, t59, 0, -t48, -t77, -pkin(7) * t72 + t65 * t53, pkin(7) * t73 + t65 * t54, t64, -pkin(2) * t48 + t64 * qJ(3), -t33 * t38, -t38 * t32 + t33 * t37, -t38 * t59, t37 * t59, 0, -t26 * t59 + t44 * t32 + t39 * t37, t27 * t59 - t44 * t33 + t39 * t38, t18 * t22, -t22 * t17 - t18 * t21, -t22 * t59, t21 * t59, 0, t31 * t17 + t24 * t21 + t75, t31 * t18 + t24 * t22 + t74, t10 * t17 + t5 * t21 + t75, -t1 * t21 + t11 * t18 - t12 * t17 + t2 * t22, -t10 * t18 - t5 * t22 - t74, t1 * t12 + t5 * t10 + t2 * t11; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t54, -0.2e1 * pkin(2) * t53, 0.2e1 * t70 * qJ(3), t70 * qJ(3) ^ 2 + pkin(2) ^ 2, t38 ^ 2, -0.2e1 * t38 * t37, 0, 0, 0, t37 * t83, t38 * t83, t22 ^ 2, t21 * t85, 0, 0, 0, t21 * t84, t22 * t84, 0.2e1 * t10 * t21, 0.2e1 * t11 * t22 - 0.2e1 * t12 * t21, t10 * t85, t10 ^ 2 + t11 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, t48, 0, 0, 0, 0, 0, t32, -t33, 0, 0, 0, 0, 0, t17, t18, t17, 0, -t18, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t53, 0, -pkin(2), 0, 0, 0, 0, 0, t37, t38, 0, 0, 0, 0, 0, t21, t22, t21, 0, -t22, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t32, -t59, t14, -t15, 0, 0, t18, -t17, -t59, -t59 * t68 + t3, t55 * t79 - t4 (-pkin(5) - t45) * t59 + t3, -t43 * t17 - t45 * t18 (-qJ(6) - t43) * t59 + t4, t1 * t43 - t2 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t37, 0, t26, -t27, 0, 0, t22, -t21, 0, -t11, -t12, -t11, -t43 * t21 - t45 * t22, t12, -t11 * t45 + t12 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, -0.2e1 * t47, 0.2e1 * t45, 0, 0.2e1 * t43, t43 ^ 2 + t45 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, -t59, t3, -t4, t3 - 0.2e1 * t78, -pkin(5) * t18 - t17 * qJ(6), -0.2e1 * t69 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, 0, -t11, -t12, -t11, -pkin(5) * t22 - t21 * qJ(6), t12, -t11 * pkin(5) + t12 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t68, -t47, t62 + t68, 0, t60 + t47, t45 * pkin(5) + t43 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t62, 0, t60, pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, t18, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
