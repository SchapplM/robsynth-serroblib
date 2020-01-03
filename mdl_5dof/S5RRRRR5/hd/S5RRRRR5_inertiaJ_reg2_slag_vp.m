% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t69 = cos(qJ(2));
t58 = t69 * pkin(1);
t52 = t58 + pkin(2);
t64 = sin(qJ(3));
t68 = cos(qJ(3));
t65 = sin(qJ(2));
t87 = t65 * pkin(1);
t75 = t68 * t87;
t31 = t64 * t52 + t75;
t29 = pkin(8) + t31;
t63 = sin(qJ(4));
t60 = t63 ^ 2;
t67 = cos(qJ(4));
t61 = t67 ^ 2;
t76 = t60 + t61;
t80 = t76 * t29;
t88 = t64 * pkin(2);
t50 = pkin(8) + t88;
t94 = t76 * t50;
t15 = (-pkin(9) - t29) * t63;
t56 = t67 * pkin(9);
t82 = t67 * t29;
t16 = t56 + t82;
t62 = sin(qJ(5));
t66 = cos(qJ(5));
t3 = t66 * t15 - t62 * t16;
t36 = t62 * t63 - t66 * t67;
t38 = t62 * t67 + t66 * t63;
t4 = t62 * t15 + t66 * t16;
t93 = -t3 * t38 - t4 * t36;
t32 = (-pkin(9) - t50) * t63;
t81 = t67 * t50;
t33 = t56 + t81;
t10 = t62 * t32 + t66 * t33;
t9 = t66 * t32 - t62 * t33;
t92 = -t10 * t36 - t9 * t38;
t40 = (-pkin(9) - pkin(8)) * t63;
t85 = t67 * pkin(8);
t41 = t56 + t85;
t17 = t66 * t40 - t62 * t41;
t18 = t62 * t40 + t66 * t41;
t91 = -t17 * t38 - t18 * t36;
t90 = pkin(3) * t63;
t89 = t62 * pkin(4);
t86 = t66 * pkin(4);
t78 = -t68 * t52 + t64 * t87;
t28 = -pkin(3) + t78;
t84 = t28 * t67;
t57 = t68 * pkin(2);
t51 = -t57 - pkin(3);
t83 = t51 * t67;
t77 = pkin(8) * t76;
t53 = -t67 * pkin(4) - pkin(3);
t59 = pkin(3) * t67;
t47 = 0.2e1 * t63 * t67;
t45 = t51 * t63;
t39 = t53 - t57;
t35 = t38 ^ 2;
t34 = t36 ^ 2;
t26 = t53 * t38;
t25 = t53 * t36;
t24 = t28 * t63;
t21 = t53 + t78;
t20 = t39 * t38;
t19 = t39 * t36;
t14 = -0.2e1 * t38 * t36;
t13 = (-t36 * t62 - t38 * t66) * pkin(4);
t12 = t21 * t38;
t11 = t21 * t36;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t58, -0.2e1 * t87, 0, (t65 ^ 2 + t69 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t78, -0.2e1 * t31, 0, t31 ^ 2 + t78 ^ 2, t60, t47, 0, t61, 0, 0, -0.2e1 * t84, 0.2e1 * t24, 0.2e1 * t80, t76 * t29 ^ 2 + t28 ^ 2, t35, t14, 0, t34, 0, 0, 0.2e1 * t11, 0.2e1 * t12, 0.2e1 * t93, t21 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t58, -t87, 0, 0, 0, 0, 0, 0, 0, 1, t57 - t78, -t75 + (-pkin(2) - t52) * t64, 0, (t31 * t64 - t68 * t78) * pkin(2), t60, t47, 0, t61, 0, 0, (-t28 - t51) * t67, t45 + t24, t94 + t80, t28 * t51 + t29 * t94, t35, t14, 0, t34, 0, 0, t19 + t11, t20 + t12, t92 + t93, t4 * t10 + t21 * t39 + t3 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t57, -0.2e1 * t88, 0, (t64 ^ 2 + t68 ^ 2) * pkin(2) ^ 2, t60, t47, 0, t61, 0, 0, -0.2e1 * t83, 0.2e1 * t45, 0.2e1 * t94, t76 * t50 ^ 2 + t51 ^ 2, t35, t14, 0, t34, 0, 0, 0.2e1 * t19, 0.2e1 * t20, 0.2e1 * t92, t10 ^ 2 + t39 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t78, -t31, 0, 0, t60, t47, 0, t61, 0, 0, t59 - t84, t24 - t90, t77 + t80, -t28 * pkin(3) + pkin(8) * t80, t35, t14, 0, t34, 0, 0, t25 + t11, t26 + t12, t91 + t93, t3 * t17 + t4 * t18 + t21 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t57, -t88, 0, 0, t60, t47, 0, t61, 0, 0, t59 - t83, t45 - t90, t77 + t94, -t51 * pkin(3) + pkin(8) * t94, t35, t14, 0, t34, 0, 0, t25 + t19, t26 + t20, t91 + t92, t10 * t18 + t9 * t17 + t39 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t60, t47, 0, t61, 0, 0, 0.2e1 * t59, -0.2e1 * t90, 0.2e1 * t77, t76 * pkin(8) ^ 2 + pkin(3) ^ 2, t35, t14, 0, t34, 0, 0, 0.2e1 * t25, 0.2e1 * t26, 0.2e1 * t91, t17 ^ 2 + t18 ^ 2 + t53 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t67, 0, -t63 * t29, -t82, 0, 0, 0, 0, t38, 0, -t36, 0, t3, -t4, t13, (t3 * t66 + t4 * t62) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t67, 0, -t63 * t50, -t81, 0, 0, 0, 0, t38, 0, -t36, 0, t9, -t10, t13, (t10 * t62 + t66 * t9) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t67, 0, -t63 * pkin(8), -t85, 0, 0, 0, 0, t38, 0, -t36, 0, t17, -t18, t13, (t17 * t66 + t18 * t62) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t86, -0.2e1 * t89, 0, (t62 ^ 2 + t66 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t36, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t36, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t36, 0, t17, -t18, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t86, -t89, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
