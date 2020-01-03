% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t42 = sin(qJ(3));
t43 = sin(qJ(2));
t45 = cos(qJ(3));
t46 = cos(qJ(2));
t20 = t42 * t43 - t45 * t46;
t22 = t42 * t46 + t45 * t43;
t35 = -t46 * pkin(2) - pkin(1);
t54 = -t22 * qJ(4) + t35;
t9 = t20 * pkin(3) + t54;
t79 = -0.2e1 * t9;
t18 = t20 ^ 2;
t19 = t22 ^ 2;
t72 = t42 * pkin(2);
t31 = qJ(4) + t72;
t78 = t31 ^ 2;
t77 = 0.2e1 * t31;
t76 = 0.2e1 * t35;
t75 = 0.2e1 * t46;
t48 = 0.2e1 * qJ(4);
t74 = pkin(3) + pkin(8);
t73 = -pkin(7) - pkin(6);
t71 = t45 * pkin(2);
t70 = t22 * t20;
t69 = t31 * t20;
t41 = sin(qJ(5));
t68 = t41 * t20;
t67 = t41 * t22;
t44 = cos(qJ(5));
t66 = t44 * t20;
t65 = t44 * t41;
t36 = t41 ^ 2;
t38 = t44 ^ 2;
t27 = t36 + t38;
t37 = t43 ^ 2;
t39 = t46 ^ 2;
t64 = t37 + t39;
t63 = qJ(4) * t20;
t62 = t31 * qJ(4);
t61 = qJ(4) + t31;
t60 = -0.2e1 * t70;
t59 = 0.2e1 * t70;
t58 = t73 * t43;
t26 = t73 * t46;
t11 = -t42 * t26 - t45 * t58;
t13 = -t45 * t26 + t42 * t58;
t57 = t11 ^ 2 + t13 ^ 2;
t34 = -pkin(3) - t71;
t25 = t27 * t74;
t4 = t74 * t20 + t54;
t7 = t22 * pkin(4) + t11;
t2 = -t41 * t4 + t44 * t7;
t3 = t44 * t4 + t41 * t7;
t1 = t2 * t44 + t3 * t41;
t29 = -pkin(8) + t34;
t56 = -t22 * t29 + t69;
t55 = t22 * t74 + t63;
t53 = 0.2e1 * t11 * t22 - 0.2e1 * t13 * t20;
t50 = -0.2e1 * pkin(3);
t49 = qJ(4) ^ 2;
t30 = -0.2e1 * t65;
t17 = t44 * t22;
t16 = t20 * t65;
t15 = t27 * t29;
t10 = (-t36 + t38) * t20;
t8 = -t20 * pkin(4) + t13;
t6 = t8 * t44;
t5 = t8 * t41;
t12 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t43 * t75, 0, t39, 0, 0, pkin(1) * t75, -0.2e1 * pkin(1) * t43, 0.2e1 * t64 * pkin(6), t64 * pkin(6) ^ 2 + pkin(1) ^ 2, t19, t60, 0, t18, 0, 0, t20 * t76, t22 * t76, t53, t35 ^ 2 + t57, 0, 0, 0, t19, t60, t18, t53, t20 * t79, t22 * t79, t9 ^ 2 + t57, t36 * t18, 0.2e1 * t18 * t65, t41 * t59, t38 * t18, t44 * t59, t19, 0.2e1 * t2 * t22 - 0.2e1 * t8 * t66, -0.2e1 * t3 * t22 + 0.2e1 * t8 * t68, 0.2e1 * (-t2 * t41 + t3 * t44) * t20, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t46, 0, -t43 * pkin(6), -t46 * pkin(6), 0, 0, 0, 0, t22, 0, -t20, 0, -t11, -t13, (-t20 * t42 - t22 * t45) * pkin(2), (-t11 * t45 + t13 * t42) * pkin(2), 0, -t22, t20, 0, 0, 0, t34 * t22 - t69, t11, t13, t11 * t34 + t13 * t31, t16, t10, t17, -t16, -t67, 0, -t56 * t44 + t5, t56 * t41 + t6, -t1, t1 * t29 + t8 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t71, -0.2e1 * t72, 0, (t42 ^ 2 + t45 ^ 2) * pkin(2) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t34, t77, t34 ^ 2 + t78, t38, t30, 0, t36, 0, 0, t41 * t77, t44 * t77, -0.2e1 * t15, t27 * t29 ^ 2 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, -t20, 0, -t11, -t13, 0, 0, 0, -t22, t20, 0, 0, 0, -pkin(3) * t22 - t63, t11, t13, -t11 * pkin(3) + t13 * qJ(4), t16, t10, t17, -t16, -t67, 0, -t55 * t44 + t5, t55 * t41 + t6, -t1, t8 * qJ(4) - t1 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t71, -t72, 0, 0, 1, 0, 0, 0, 0, 0, 0, t50 - t71, t48 + t72, -t34 * pkin(3) + t62, t38, t30, 0, t36, 0, 0, t61 * t41, t61 * t44, (-t29 + t74) * t27, -t29 * t25 + t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, t50, t48, pkin(3) ^ 2 + t49, t38, t30, 0, t36, 0, 0, t41 * t48, t44 * t48, 0.2e1 * t25, t27 * t74 ^ 2 + t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, 0, t11, 0, 0, 0, 0, 0, 0, t17, -t67, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, t66, t22, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t41, 0, t44 * t29, -t41 * t29, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t41, 0, -t44 * t74, t41 * t74, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t41, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t12;
