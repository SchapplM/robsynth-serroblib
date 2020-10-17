% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR9
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
% Datum: 2019-12-31 22:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR9_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR9_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:29:52
% EndTime: 2019-12-31 22:29:57
% DurationCPUTime: 1.13s
% Computational Cost: add. (1387->155), mult. (2856->296), div. (0->0), fcn. (3102->8), ass. (0->79)
t58 = sin(qJ(4));
t59 = sin(qJ(3));
t62 = cos(qJ(4));
t63 = cos(qJ(3));
t35 = t58 * t59 - t62 * t63;
t48 = -t63 * pkin(3) - pkin(2);
t27 = t35 * pkin(4) + t48;
t94 = 0.2e1 * t27;
t93 = 0.2e1 * t48;
t60 = sin(qJ(2));
t92 = -0.2e1 * t60;
t64 = cos(qJ(2));
t91 = -0.2e1 * t64;
t90 = 0.2e1 * t64;
t89 = -pkin(8) - pkin(7);
t88 = pkin(2) * t63;
t87 = pkin(6) * t59;
t54 = t60 ^ 2;
t86 = t54 * pkin(6);
t57 = sin(qJ(5));
t85 = t57 * pkin(4);
t84 = t58 * pkin(3);
t50 = t60 * pkin(6);
t61 = cos(qJ(5));
t37 = t58 * t63 + t62 * t59;
t28 = t37 * t60;
t40 = -t64 * pkin(2) - t60 * pkin(7) - pkin(1);
t34 = t63 * t40;
t76 = t63 * t60;
t19 = -pkin(8) * t76 + t34 + (-pkin(3) - t87) * t64;
t75 = t63 * t64;
t71 = pkin(6) * t75;
t21 = t71 + (-pkin(8) * t60 + t40) * t59;
t77 = t62 * t21;
t9 = t58 * t19 + t77;
t7 = -t28 * pkin(9) + t9;
t83 = t61 * t7;
t82 = t64 * pkin(3);
t81 = t64 * pkin(4);
t80 = t59 * t60;
t79 = t59 * t63;
t78 = t59 * t64;
t39 = pkin(3) * t80 + t50;
t53 = t59 ^ 2;
t55 = t63 ^ 2;
t74 = t53 + t55;
t73 = t60 * t90;
t72 = t61 * t84;
t70 = t59 * t76;
t30 = -t58 * t80 + t62 * t76;
t8 = t62 * t19 - t58 * t21;
t4 = -t30 * pkin(9) + t8 - t81;
t1 = t61 * t4 - t57 * t7;
t41 = t89 * t59;
t42 = t89 * t63;
t22 = t62 * t41 + t58 * t42;
t52 = t62 * pkin(3);
t47 = t52 + pkin(4);
t31 = t61 * t47 - t57 * t84;
t2 = t57 * t4 + t83;
t25 = -pkin(6) * t78 + t34;
t26 = t59 * t40 + t71;
t69 = -t25 * t59 + t26 * t63;
t23 = t58 * t41 - t62 * t42;
t66 = pkin(6) ^ 2;
t56 = t64 ^ 2;
t51 = t61 * pkin(4);
t49 = t54 * t66;
t32 = t57 * t47 + t72;
t20 = t28 * pkin(4) + t39;
t18 = -t57 * t35 + t61 * t37;
t16 = t61 * t35 + t57 * t37;
t14 = -t35 * pkin(9) + t23;
t13 = -t37 * pkin(9) + t22;
t12 = -t57 * t28 + t61 * t30;
t10 = t61 * t28 + t57 * t30;
t6 = t57 * t13 + t61 * t14;
t5 = t61 * t13 - t57 * t14;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t54, t73, 0, t56, 0, 0, pkin(1) * t90, pkin(1) * t92, 0.2e1 * (t54 + t56) * pkin(6), pkin(1) ^ 2 + t56 * t66 + t49, t55 * t54, -0.2e1 * t54 * t79, t75 * t92, t53 * t54, t59 * t73, t56, -0.2e1 * t25 * t64 + 0.2e1 * t59 * t86, 0.2e1 * t26 * t64 + 0.2e1 * t63 * t86, 0.2e1 * (-t25 * t63 - t26 * t59) * t60, t25 ^ 2 + t26 ^ 2 + t49, t30 ^ 2, -0.2e1 * t30 * t28, t30 * t91, t28 ^ 2, -t28 * t91, t56, 0.2e1 * t39 * t28 - 0.2e1 * t8 * t64, 0.2e1 * t39 * t30 + 0.2e1 * t9 * t64, -0.2e1 * t9 * t28 - 0.2e1 * t8 * t30, t39 ^ 2 + t8 ^ 2 + t9 ^ 2, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t91, t10 ^ 2, t10 * t90, t56, -0.2e1 * t1 * t64 + 0.2e1 * t20 * t10, 0.2e1 * t20 * t12 + 0.2e1 * t2 * t64, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t10, t1 ^ 2 + t2 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, 0, t64, 0, -t50, -t64 * pkin(6), 0, 0, t70, (-t53 + t55) * t60, -t78, -t70, -t75, 0, -pkin(6) * t76 + (-pkin(2) * t60 + pkin(7) * t64) * t59, pkin(7) * t75 + (t87 - t88) * t60, t69, -pkin(2) * t50 + t69 * pkin(7), t30 * t37, -t37 * t28 - t30 * t35, -t64 * t37, t28 * t35, t64 * t35, 0, -t22 * t64 + t48 * t28 + t39 * t35, t23 * t64 + t48 * t30 + t39 * t37, -t22 * t30 - t23 * t28 - t9 * t35 - t8 * t37, t8 * t22 + t9 * t23 + t39 * t48, t12 * t18, -t18 * t10 - t12 * t16, -t18 * t64, t10 * t16, t16 * t64, 0, t27 * t10 + t20 * t16 - t5 * t64, t27 * t12 + t20 * t18 + t6 * t64, -t1 * t18 - t6 * t10 - t5 * t12 - t2 * t16, t1 * t5 + t2 * t6 + t20 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t53, 0.2e1 * t79, 0, t55, 0, 0, 0.2e1 * t88, -0.2e1 * pkin(2) * t59, 0.2e1 * t74 * pkin(7), t74 * pkin(7) ^ 2 + pkin(2) ^ 2, t37 ^ 2, -0.2e1 * t37 * t35, 0, t35 ^ 2, 0, 0, t35 * t93, t37 * t93, -0.2e1 * t22 * t37 - 0.2e1 * t23 * t35, t22 ^ 2 + t23 ^ 2 + t48 ^ 2, t18 ^ 2, -0.2e1 * t18 * t16, 0, t16 ^ 2, 0, 0, t16 * t94, t18 * t94, -0.2e1 * t6 * t16 - 0.2e1 * t5 * t18, t27 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, -t80, -t64, t25, -t26, 0, 0, 0, 0, t30, 0, -t28, -t64, -t62 * t82 + t8, -t77 + (-t19 + t82) * t58, (-t28 * t58 - t30 * t62) * pkin(3), (t58 * t9 + t62 * t8) * pkin(3), 0, 0, t12, 0, -t10, -t64, -t31 * t64 + t1, t32 * t64 - t2, -t32 * t10 - t31 * t12, t1 * t31 + t2 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t63, 0, -t59 * pkin(7), -t63 * pkin(7), 0, 0, 0, 0, t37, 0, -t35, 0, t22, -t23, (-t35 * t58 - t37 * t62) * pkin(3), (t22 * t62 + t23 * t58) * pkin(3), 0, 0, t18, 0, -t16, 0, t5, -t6, -t32 * t16 - t31 * t18, t5 * t31 + t6 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t84, 0, (t58 ^ 2 + t62 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t31, -0.2e1 * t32, 0, t31 ^ 2 + t32 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t28, -t64, t8, -t9, 0, 0, 0, 0, t12, 0, -t10, -t64, -t61 * t81 + t1, -t83 + (-t4 + t81) * t57, (-t10 * t57 - t12 * t61) * pkin(4), (t1 * t61 + t2 * t57) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t35, 0, t22, -t23, 0, 0, 0, 0, t18, 0, -t16, 0, t5, -t6, (-t16 * t57 - t18 * t61) * pkin(4), (t5 * t61 + t57 * t6) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t52, -t84, 0, 0, 0, 0, 0, 0, 0, 1, t31 + t51, -t72 + (-pkin(4) - t47) * t57, 0, (t31 * t61 + t32 * t57) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t51, -0.2e1 * t85, 0, (t57 ^ 2 + t61 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, -t64, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t31, -t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t51, -t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
