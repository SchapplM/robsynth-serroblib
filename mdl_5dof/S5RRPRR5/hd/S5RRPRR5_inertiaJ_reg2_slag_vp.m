% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:02:59
% DurationCPUTime: 0.80s
% Computational Cost: add. (747->75), mult. (1386->135), div. (0->0), fcn. (1567->8), ass. (0->67)
t58 = sin(pkin(9));
t56 = t58 ^ 2;
t59 = cos(pkin(9));
t57 = t59 ^ 2;
t69 = t56 + t57;
t70 = t69 * qJ(3);
t61 = sin(qJ(4));
t64 = cos(qJ(4));
t39 = t58 * t61 - t59 * t64;
t51 = -pkin(3) * t59 - pkin(2);
t30 = pkin(4) * t39 + t51;
t65 = cos(qJ(2));
t77 = t65 * pkin(1);
t29 = t30 - t77;
t88 = 0.2e1 * t29;
t87 = 0.2e1 * t30;
t42 = t51 - t77;
t86 = 0.2e1 * t42;
t85 = 0.2e1 * t51;
t84 = 0.2e1 * t59;
t41 = t58 * t64 + t59 * t61;
t60 = sin(qJ(5));
t63 = cos(qJ(5));
t23 = t39 * t63 + t41 * t60;
t25 = -t39 * t60 + t41 * t63;
t62 = sin(qJ(2));
t79 = t62 * pkin(1);
t50 = qJ(3) + t79;
t33 = (-pkin(7) - t50) * t58;
t55 = t59 * pkin(7);
t34 = t50 * t59 + t55;
t19 = t64 * t33 - t34 * t61;
t81 = t41 * pkin(8);
t11 = t19 - t81;
t20 = t33 * t61 + t34 * t64;
t35 = t39 * pkin(8);
t12 = t20 - t35;
t5 = t11 * t63 - t12 * t60;
t6 = t11 * t60 + t12 * t63;
t83 = -t6 * t23 - t5 * t25;
t43 = (-pkin(7) - qJ(3)) * t58;
t44 = qJ(3) * t59 + t55;
t27 = t64 * t43 - t44 * t61;
t15 = t27 - t81;
t28 = t43 * t61 + t44 * t64;
t16 = t28 - t35;
t7 = t15 * t63 - t16 * t60;
t8 = t15 * t60 + t16 * t63;
t82 = -t8 * t23 - t7 * t25;
t80 = t60 * pkin(4);
t78 = t63 * pkin(4);
t52 = -pkin(2) - t77;
t76 = pkin(2) - t52;
t75 = -t19 * t41 - t20 * t39;
t74 = -t27 * t41 - t28 * t39;
t73 = t29 + t30;
t72 = t42 + t51;
t71 = t69 * t50;
t47 = t58 * t84;
t37 = t41 ^ 2;
t36 = t39 ^ 2;
t26 = -0.2e1 * t41 * t39;
t22 = t25 ^ 2;
t21 = t23 ^ 2;
t10 = -0.2e1 * t25 * t23;
t9 = (-t23 * t60 - t25 * t63) * pkin(4);
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t79, 0, (t62 ^ 2 + t65 ^ 2) * pkin(1) ^ 2, t56, t47, 0, t57, 0, 0, -0.2e1 * t52 * t59, 0.2e1 * t52 * t58, 0.2e1 * t71, t50 ^ 2 * t69 + t52 ^ 2, t37, t26, 0, t36, 0, 0, t39 * t86, t41 * t86, 0.2e1 * t75, t19 ^ 2 + t20 ^ 2 + t42 ^ 2, t22, t10, 0, t21, 0, 0, t23 * t88, t25 * t88, 0.2e1 * t83, t29 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t77, -t79, 0, 0, t56, t47, 0, t57, 0, 0, t76 * t59, -t76 * t58, t70 + t71, -t52 * pkin(2) + t50 * t70, t37, t26, 0, t36, 0, 0, t72 * t39, t72 * t41, t74 + t75, t19 * t27 + t20 * t28 + t42 * t51, t22, t10, 0, t21, 0, 0, t73 * t23, t73 * t25, t82 + t83, t29 * t30 + t5 * t7 + t6 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t56, t47, 0, t57, 0, 0, pkin(2) * t84, -0.2e1 * pkin(2) * t58, 0.2e1 * t70, qJ(3) ^ 2 * t69 + pkin(2) ^ 2, t37, t26, 0, t36, 0, 0, t39 * t85, t41 * t85, 0.2e1 * t74, t27 ^ 2 + t28 ^ 2 + t51 ^ 2, t22, t10, 0, t21, 0, 0, t23 * t87, t25 * t87, 0.2e1 * t82, t30 ^ 2 + t7 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0, t52, 0, 0, 0, 0, 0, 0, t39, t41, 0, t42, 0, 0, 0, 0, 0, 0, t23, t25, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, t58, 0, -pkin(2), 0, 0, 0, 0, 0, 0, t39, t41, 0, t51, 0, 0, 0, 0, 0, 0, t23, t25, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t39, 0, t19, -t20, 0, 0, 0, 0, t25, 0, -t23, 0, t5, -t6, t9, (t5 * t63 + t6 * t60) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t39, 0, t27, -t28, 0, 0, 0, 0, t25, 0, -t23, 0, t7, -t8, t9, (t60 * t8 + t63 * t7) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t78, -0.2e1 * t80, 0, (t60 ^ 2 + t63 ^ 2) * pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, -t23, 0, t7, -t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t78, -t80, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
