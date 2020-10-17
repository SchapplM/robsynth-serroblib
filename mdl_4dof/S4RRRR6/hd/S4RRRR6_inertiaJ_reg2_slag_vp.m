% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:51
% EndTime: 2019-12-31 17:30:53
% DurationCPUTime: 0.62s
% Computational Cost: add. (599->109), mult. (1535->235), div. (0->0), fcn. (1624->8), ass. (0->76)
t33 = cos(pkin(4));
t35 = sin(qJ(3));
t38 = cos(qJ(3));
t32 = sin(pkin(4));
t36 = sin(qJ(2));
t62 = t32 * t36;
t18 = t33 * t35 + t38 * t62;
t34 = sin(qJ(4));
t37 = cos(qJ(4));
t39 = cos(qJ(2));
t61 = t32 * t39;
t8 = t18 * t34 + t37 * t61;
t80 = -0.2e1 * t8;
t16 = -t33 * t38 + t35 * t62;
t79 = t16 ^ 2;
t78 = -0.2e1 * t18;
t77 = 0.2e1 * t32;
t76 = -0.2e1 * t35;
t75 = 0.2e1 * t38;
t74 = pkin(1) * t36;
t73 = pkin(1) * t39;
t72 = pkin(3) * t37;
t71 = pkin(7) * t34;
t70 = t35 * pkin(7);
t48 = pkin(6) * t61;
t14 = t48 + (pkin(7) + t74) * t33;
t15 = (-pkin(2) * t39 - pkin(7) * t36 - pkin(1)) * t32;
t6 = -t35 * t14 + t38 * t15;
t4 = pkin(3) * t61 - t6;
t69 = t4 * t34;
t68 = t4 * t37;
t67 = t8 * t37;
t10 = t18 * t37 - t34 * t61;
t66 = t10 * t34;
t65 = t16 * t38;
t64 = t18 * t35;
t27 = t32 ^ 2;
t63 = t27 * t39;
t60 = t33 * t36;
t59 = t34 * t16;
t58 = t34 * t35;
t57 = t34 * t37;
t56 = t34 * t38;
t55 = t35 * t16;
t54 = t37 * t16;
t53 = t37 * t35;
t52 = t37 * t38;
t28 = t34 ^ 2;
t30 = t37 ^ 2;
t51 = t28 + t30;
t50 = 0.2e1 * t61;
t49 = t35 * t75;
t47 = t35 * t61;
t46 = t38 * t61;
t45 = t34 * t53;
t23 = pkin(6) * t62;
t13 = t23 + (-pkin(2) - t73) * t33;
t3 = t16 * pkin(3) - t18 * pkin(8) + t13;
t7 = t38 * t14 + t35 * t15;
t5 = -pkin(8) * t61 + t7;
t1 = t37 * t3 - t34 * t5;
t2 = t34 * t3 + t37 * t5;
t44 = -t1 * t34 + t2 * t37;
t43 = -t6 * t35 + t7 * t38;
t21 = -t38 * pkin(3) - t35 * pkin(8) - pkin(2);
t11 = -pkin(7) * t56 + t37 * t21;
t12 = pkin(7) * t52 + t34 * t21;
t42 = -t11 * t34 + t12 * t37;
t41 = pkin(7) ^ 2;
t31 = t38 ^ 2;
t29 = t35 ^ 2;
t26 = t29 * t41;
t24 = t27 * t39 ^ 2;
t20 = pkin(1) * t60 + t48;
t19 = t33 * t73 - t23;
t9 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t27 * t36 ^ 2, 0.2e1 * t36 * t63, t60 * t77, t24, t33 * t50, t33 ^ 2, 0.2e1 * pkin(1) * t63 + 0.2e1 * t19 * t33, -0.2e1 * t20 * t33 - 0.2e1 * t27 * t74, (-t19 * t36 + t20 * t39) * t77, t27 * pkin(1) ^ 2 + t19 ^ 2 + t20 ^ 2, t18 ^ 2, t16 * t78, t61 * t78, t79, t16 * t50, t24, 0.2e1 * t13 * t16 - 0.2e1 * t6 * t61, 0.2e1 * t13 * t18 + 0.2e1 * t7 * t61, -0.2e1 * t7 * t16 - 0.2e1 * t6 * t18, t13 ^ 2 + t6 ^ 2 + t7 ^ 2, t10 ^ 2, t10 * t80, 0.2e1 * t10 * t16, t8 ^ 2, t16 * t80, t79, 0.2e1 * t1 * t16 + 0.2e1 * t4 * t8, 0.2e1 * t4 * t10 - 0.2e1 * t2 * t16, -0.2e1 * t1 * t10 - 0.2e1 * t2 * t8, t1 ^ 2 + t2 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t61, t33, t19, -t20, 0, 0, t64, t18 * t38 - t55, -t47, -t65, -t46, 0, -pkin(2) * t16 + pkin(7) * t47 - t13 * t38, -pkin(2) * t18 + pkin(7) * t46 + t13 * t35, (t64 - t65) * pkin(7) + t43, -t13 * pkin(2) + t43 * pkin(7), t10 * t53, (-t66 - t67) * t35, -t10 * t38 + t16 * t53, t8 * t58, -t34 * t55 + t8 * t38, -t65, -t1 * t38 + t11 * t16 + (pkin(7) * t8 + t69) * t35, -t12 * t16 + t2 * t38 + (pkin(7) * t10 + t68) * t35, -t11 * t10 - t12 * t8 + (-t1 * t37 - t2 * t34) * t35, t1 * t11 + t2 * t12 + t4 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t29, t49, 0, t31, 0, 0, pkin(2) * t75, pkin(2) * t76, 0.2e1 * (t29 + t31) * pkin(7), pkin(2) ^ 2 + t31 * t41 + t26, t30 * t29, -0.2e1 * t29 * t57, t52 * t76, t28 * t29, t34 * t49, t31, -0.2e1 * t11 * t38 + 0.2e1 * t29 * t71, 0.2e1 * t29 * pkin(7) * t37 + 0.2e1 * t12 * t38, 0.2e1 * (-t11 * t37 - t12 * t34) * t35, t11 ^ 2 + t12 ^ 2 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, -t61, t6, -t7, 0, 0, t66, t10 * t37 - t34 * t8, t59, -t67, t54, 0, -pkin(3) * t8 - pkin(8) * t59 - t68, -pkin(3) * t10 - pkin(8) * t54 + t69, (t66 - t67) * pkin(8) + t44, -t4 * pkin(3) + t44 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t38, 0, -t70, -t38 * pkin(7), 0, 0, t45, (-t28 + t30) * t35, -t56, -t45, -t52, 0, -pkin(7) * t53 + (-pkin(3) * t35 + pkin(8) * t38) * t34, pkin(8) * t52 + (t71 - t72) * t35, t42, -pkin(3) * t70 + t42 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t28, 0.2e1 * t57, 0, t30, 0, 0, 0.2e1 * t72, -0.2e1 * pkin(3) * t34, 0.2e1 * t51 * pkin(8), t51 * pkin(8) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, -t8, t16, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, 0, -t58, -t38, t11, -t12, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t37, 0, -t34 * pkin(8), -t37 * pkin(8), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t9;
