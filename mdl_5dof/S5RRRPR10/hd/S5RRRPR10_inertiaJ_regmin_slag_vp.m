% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:17
% EndTime: 2021-01-15 23:41:20
% DurationCPUTime: 0.70s
% Computational Cost: add. (998->114), mult. (2374->257), div. (0->0), fcn. (2712->10), ass. (0->80)
t54 = cos(pkin(5));
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t52 = sin(pkin(5));
t57 = sin(qJ(2));
t80 = t52 * t57;
t33 = -t54 * t59 + t56 * t80;
t34 = t54 * t56 + t59 * t80;
t51 = sin(pkin(10));
t53 = cos(pkin(10));
t21 = -t51 * t33 + t53 * t34;
t55 = sin(qJ(5));
t58 = cos(qJ(5));
t60 = cos(qJ(2));
t79 = t52 * t60;
t15 = t55 * t21 + t58 * t79;
t90 = -0.2e1 * t15;
t89 = -0.2e1 * t34;
t47 = -t59 * pkin(3) - pkin(2);
t88 = 0.2e1 * t47;
t87 = 0.2e1 * t59;
t86 = pkin(1) * t57;
t85 = pkin(1) * t60;
t84 = t51 * pkin(3);
t83 = t53 * pkin(3);
t67 = pkin(7) * t79;
t29 = t67 + (pkin(8) + t86) * t54;
t30 = (-pkin(2) * t60 - pkin(8) * t57 - pkin(1)) * t52;
t17 = -t56 * t29 + t59 * t30;
t68 = pkin(3) * t79;
t10 = -t34 * qJ(4) + t17 - t68;
t18 = t59 * t29 + t56 * t30;
t14 = -t33 * qJ(4) + t18;
t6 = t51 * t10 + t53 * t14;
t16 = t58 * t21 - t55 * t79;
t82 = t16 * t55;
t48 = t52 ^ 2;
t81 = t48 * t60;
t78 = t54 * t57;
t20 = t53 * t33 + t51 * t34;
t77 = t55 * t20;
t38 = t51 * t56 - t53 * t59;
t76 = t55 * t38;
t39 = t51 * t59 + t53 * t56;
t75 = t55 * t39;
t45 = pkin(9) + t84;
t74 = t55 * t45;
t73 = t55 * t58;
t72 = t58 * t39;
t71 = t58 * t45;
t70 = qJ(4) + pkin(8);
t69 = 0.2e1 * t79;
t66 = t56 * t79;
t65 = t59 * t79;
t64 = -t53 * t10 + t51 * t14;
t63 = t70 * t56;
t46 = -pkin(4) - t83;
t62 = -t38 * t45 + t39 * t46;
t42 = pkin(7) * t80;
t28 = t42 + (-pkin(2) - t85) * t54;
t22 = t33 * pkin(3) + t28;
t50 = t58 ^ 2;
t49 = t55 ^ 2;
t41 = t70 * t59;
t37 = t39 ^ 2;
t36 = pkin(1) * t78 + t67;
t35 = t54 * t85 - t42;
t32 = t58 * t38;
t26 = t53 * t41 - t51 * t63;
t24 = t51 * t41 + t53 * t63;
t23 = t38 * pkin(4) - t39 * pkin(9) + t47;
t19 = t58 * t20;
t12 = t55 * t23 + t58 * t26;
t11 = t58 * t23 - t55 * t26;
t7 = t20 * pkin(4) - t21 * pkin(9) + t22;
t4 = -pkin(9) * t79 + t6;
t3 = pkin(4) * t79 + t64;
t2 = t58 * t4 + t55 * t7;
t1 = -t55 * t4 + t58 * t7;
t5 = [1, 0, 0, t48 * t57 ^ 2, 0.2e1 * t57 * t81, 0.2e1 * t52 * t78, t54 * t69, t54 ^ 2, 0.2e1 * pkin(1) * t81 + 0.2e1 * t35 * t54, -0.2e1 * t36 * t54 - 0.2e1 * t48 * t86, t34 ^ 2, t33 * t89, t79 * t89, t33 * t69, t48 * t60 ^ 2, -0.2e1 * t17 * t79 + 0.2e1 * t28 * t33, 0.2e1 * t18 * t79 + 0.2e1 * t28 * t34, 0.2e1 * t22 * t20 + 0.2e1 * t64 * t79, 0.2e1 * t22 * t21 + 0.2e1 * t6 * t79, -0.2e1 * t6 * t20 + 0.2e1 * t21 * t64, t22 ^ 2 + t6 ^ 2 + t64 ^ 2, t16 ^ 2, t16 * t90, 0.2e1 * t16 * t20, t20 * t90, t20 ^ 2, 0.2e1 * t1 * t20 + 0.2e1 * t3 * t15, 0.2e1 * t3 * t16 - 0.2e1 * t2 * t20; 0, 0, 0, 0, 0, t80, t79, t54, t35, -t36, t34 * t56, -t56 * t33 + t34 * t59, -t66, -t65, 0, -pkin(2) * t33 + pkin(8) * t66 - t28 * t59, -pkin(2) * t34 + pkin(8) * t65 + t28 * t56, t47 * t20 + t22 * t38 + t24 * t79, t47 * t21 + t22 * t39 + t26 * t79, -t26 * t20 + t24 * t21 - t6 * t38 + t39 * t64, t22 * t47 + t24 * t64 + t6 * t26, t16 * t72, (-t15 * t58 - t82) * t39, t16 * t38 + t20 * t72, -t15 * t38 - t20 * t75, t20 * t38, t1 * t38 + t11 * t20 + t24 * t15 + t3 * t75, -t12 * t20 + t24 * t16 - t2 * t38 + t3 * t72; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56 ^ 2, t56 * t87, 0, 0, 0, pkin(2) * t87, -0.2e1 * pkin(2) * t56, t38 * t88, t39 * t88, 0.2e1 * t24 * t39 - 0.2e1 * t26 * t38, t24 ^ 2 + t26 ^ 2 + t47 ^ 2, t50 * t37, -0.2e1 * t37 * t73, 0.2e1 * t38 * t72, -0.2e1 * t38 * t75, t38 ^ 2, 0.2e1 * t11 * t38 + 0.2e1 * t24 * t75, -0.2e1 * t12 * t38 + 0.2e1 * t24 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t33, -t79, t17, -t18, -t53 * t68 - t64, t51 * t68 - t6, (-t20 * t51 - t21 * t53) * pkin(3), (t51 * t6 - t53 * t64) * pkin(3), t82, -t55 * t15 + t16 * t58, t77, t19, 0, t46 * t15 - t20 * t74 - t3 * t58, t46 * t16 - t20 * t71 + t3 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t59, 0, -t56 * pkin(8), -t59 * pkin(8), -t24, -t26, (-t38 * t51 - t39 * t53) * pkin(3), (-t24 * t53 + t26 * t51) * pkin(3), t55 * t72, (-t49 + t50) * t39, t76, t32, 0, -t24 * t58 + t62 * t55, t24 * t55 + t62 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * t83, -0.2e1 * t84, 0, (t51 ^ 2 + t53 ^ 2) * pkin(3) ^ 2, t49, 0.2e1 * t73, 0, 0, 0, -0.2e1 * t46 * t58, 0.2e1 * t46 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t21, 0, t22, 0, 0, 0, 0, 0, t19, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t39, 0, t47, 0, 0, 0, 0, 0, t32, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t20, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t75, t38, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t58, 0, -t74, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
