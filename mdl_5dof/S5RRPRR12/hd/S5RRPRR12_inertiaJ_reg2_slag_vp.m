% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t43 = sin(qJ(2));
t37 = t43 ^ 2;
t46 = cos(qJ(2));
t40 = t46 ^ 2;
t81 = t37 + t40;
t41 = sin(qJ(5));
t35 = t41 ^ 2;
t44 = cos(qJ(5));
t38 = t44 ^ 2;
t57 = t35 + t38;
t80 = t57 * pkin(8);
t42 = sin(qJ(4));
t45 = cos(qJ(4));
t72 = -pkin(2) - pkin(3);
t23 = t45 * qJ(3) + t42 * t72;
t19 = -pkin(8) + t23;
t53 = t57 * t19;
t32 = t46 * pkin(6);
t25 = -t46 * pkin(7) + t32;
t54 = (pkin(6) - pkin(7)) * t43;
t6 = t42 * t25 - t45 * t54;
t79 = t6 ^ 2;
t12 = t43 * t42 + t46 * t45;
t78 = t12 ^ 2;
t77 = 0.2e1 * t12;
t14 = -t46 * t42 + t43 * t45;
t76 = 0.2e1 * t14;
t75 = -0.2e1 * t41;
t74 = -0.2e1 * t43;
t73 = 0.2e1 * t44;
t71 = t43 * pkin(6);
t70 = t6 * t41;
t69 = t6 * t44;
t68 = t6 * t45;
t21 = t42 * qJ(3) - t45 * t72;
t18 = pkin(4) + t21;
t67 = pkin(4) + t18;
t66 = t41 * t12;
t65 = t41 * t14;
t64 = t41 * t44;
t63 = t43 * t46;
t62 = t44 * t12;
t61 = t44 * t14;
t60 = t45 * t41;
t59 = t45 * t44;
t58 = t81 * pkin(6) ^ 2;
t56 = -0.2e1 * t14 * t12;
t55 = -0.2e1 * t64;
t24 = -t46 * pkin(2) - t43 * qJ(3) - pkin(1);
t15 = t57 * t42;
t10 = t46 * pkin(3) - t24;
t52 = -pkin(4) * t14 - pkin(8) * t12;
t4 = t12 * pkin(4) - t14 * pkin(8) + t10;
t8 = t45 * t25 + t42 * t54;
t2 = t44 * t4 - t41 * t8;
t3 = t41 * t4 + t44 * t8;
t1 = -t2 * t41 + t3 * t44;
t51 = -t43 * pkin(2) + t46 * qJ(3);
t50 = -t12 * t19 + t14 * t18;
t49 = -t42 * t12 - t45 * t14;
t39 = t45 ^ 2;
t36 = t42 ^ 2;
t26 = 0.2e1 * t64;
t20 = 0.2e1 * t81 * pkin(6);
t11 = t14 ^ 2;
t9 = t41 * t61;
t5 = (t35 - t38) * t14;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, 0.2e1 * t63, 0, t40, 0, 0, 0.2e1 * pkin(1) * t46, pkin(1) * t74, t20, pkin(1) ^ 2 + t58, t37, 0, -0.2e1 * t63, 0, 0, t40, -0.2e1 * t24 * t46, t20, t24 * t74, t24 ^ 2 + t58, t11, t56, 0, t78, 0, 0, t10 * t77, t10 * t76, -0.2e1 * t8 * t12 + 0.2e1 * t6 * t14, t10 ^ 2 + t8 ^ 2 + t79, t38 * t11, t11 * t55, t61 * t77, t35 * t11, t41 * t56, t78, 0.2e1 * t2 * t12 + 0.2e1 * t6 * t65, -0.2e1 * t3 * t12 + 0.2e1 * t6 * t61, (-t2 * t44 - t3 * t41) * t76, t2 ^ 2 + t3 ^ 2 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t46, 0, -t71, -t32, 0, 0, 0, t43, 0, 0, -t46, 0, -t71, t51, t32, t51 * pkin(6), 0, 0, -t14, 0, t12, 0, t6, t8, -t23 * t12 + t21 * t14, t6 * t21 + t8 * t23, -t9, t5, -t66, t9, -t62, 0, t50 * t41 + t69, t50 * t44 - t70, -t1, t1 * t19 + t6 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t21, 0.2e1 * t23, 0, t21 ^ 2 + t23 ^ 2, t35, t26, 0, t38, 0, 0, t18 * t73, t18 * t75, -0.2e1 * t53, t57 * t19 ^ 2 + t18 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t71, 0, 0, 0, 0, 0, 0, 0, 0, t49, t8 * t42 - t68, 0, 0, 0, 0, 0, 0, t49 * t41, t49 * t44, 0, t1 * t42 - t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t45, t42, 0, -t21 * t45 + t23 * t42, 0, 0, 0, 0, 0, 0, -t59, t60, -t15, t19 * t15 - t18 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 + t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t36 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, -t12, 0, -t6, -t8, 0, 0, t9, -t5, t66, -t9, t62, 0, t52 * t41 - t69, t52 * t44 + t70, t1, -t6 * pkin(4) + t1 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t21, -t23, 0, 0, -t35, t55, 0, -t38, 0, 0, -t67 * t44, t67 * t41, t53 - t80, -t18 * pkin(4) + pkin(8) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t60, t15, t45 * pkin(4) + pkin(8) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t35, t26, 0, t38, 0, 0, pkin(4) * t73, pkin(4) * t75, 0.2e1 * t80, t57 * pkin(8) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t65, t12, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, -t44, 0, -t41 * t19, -t44 * t19, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 * t42, -t44 * t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t44, 0, -t41 * pkin(8), -t44 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t7;
