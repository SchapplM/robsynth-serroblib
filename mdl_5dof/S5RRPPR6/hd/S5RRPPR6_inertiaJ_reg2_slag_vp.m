% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:05
% EndTime: 2019-12-31 19:33:09
% DurationCPUTime: 0.77s
% Computational Cost: add. (862->85), mult. (1689->188), div. (0->0), fcn. (1952->8), ass. (0->69)
t50 = sin(pkin(8));
t52 = cos(pkin(8));
t54 = sin(qJ(2));
t55 = cos(qJ(2));
t34 = t50 * t55 + t52 * t54;
t85 = -0.2e1 * t34;
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t53 = sin(qJ(5));
t75 = cos(qJ(5));
t84 = -t53 * t49 + t75 * t51;
t67 = -qJ(3) - pkin(6);
t38 = t67 * t55;
t61 = t67 * t54;
t19 = -t50 * t38 - t52 * t61;
t83 = t19 ^ 2;
t30 = t50 * t54 - t52 * t55;
t26 = t30 ^ 2;
t82 = 0.2e1 * t30;
t77 = t52 * pkin(2);
t43 = -pkin(3) - t77;
t37 = -t51 * pkin(4) + t43;
t81 = 0.2e1 * t37;
t44 = -t55 * pkin(2) - pkin(1);
t80 = 0.2e1 * t44;
t79 = 0.2e1 * t55;
t78 = t50 * pkin(2);
t40 = qJ(4) + t78;
t76 = pkin(7) + t40;
t11 = t84 * t34;
t74 = t11 * t84;
t35 = t75 * t49 + t53 * t51;
t73 = t35 * t30;
t72 = t49 * t30;
t71 = t49 * t34;
t70 = t49 * t51;
t69 = t51 * t34;
t16 = t30 * pkin(3) - t34 * qJ(4) + t44;
t21 = -t52 * t38 + t50 * t61;
t6 = t49 * t16 + t51 * t21;
t45 = t49 ^ 2;
t46 = t51 ^ 2;
t66 = t45 + t46;
t47 = t54 ^ 2;
t48 = t55 ^ 2;
t65 = t47 + t48;
t64 = t30 * t85;
t63 = t49 * t69;
t5 = t51 * t16 - t49 * t21;
t60 = t6 * t49 + t5 * t51;
t59 = -t5 * t49 + t6 * t51;
t58 = -t30 * t40 + t34 * t43;
t29 = t35 ^ 2;
t28 = t34 ^ 2;
t27 = t84 ^ 2;
t25 = t76 * t51;
t24 = t76 * t49;
t23 = t51 * t30;
t18 = t84 * t30;
t15 = -t53 * t24 + t75 * t25;
t14 = -t75 * t24 - t53 * t25;
t9 = t35 * t34;
t8 = pkin(4) * t71 + t19;
t7 = t35 * t9;
t4 = -pkin(7) * t71 + t6;
t3 = t30 * pkin(4) - pkin(7) * t69 + t5;
t2 = t53 * t3 + t75 * t4;
t1 = t75 * t3 - t53 * t4;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t47, t54 * t79, 0, t48, 0, 0, pkin(1) * t79, -0.2e1 * pkin(1) * t54, 0.2e1 * t65 * pkin(6), t65 * pkin(6) ^ 2 + pkin(1) ^ 2, t28, t64, 0, t26, 0, 0, t30 * t80, t34 * t80, 0.2e1 * t19 * t34 - 0.2e1 * t21 * t30, t21 ^ 2 + t44 ^ 2 + t83, t46 * t28, -0.2e1 * t28 * t70, t69 * t82, t45 * t28, t49 * t64, t26, 0.2e1 * t19 * t71 + 0.2e1 * t5 * t30, 0.2e1 * t19 * t69 - 0.2e1 * t6 * t30, t60 * t85, t5 ^ 2 + t6 ^ 2 + t83, t11 ^ 2, -0.2e1 * t11 * t9, t11 * t82, t9 ^ 2, -t9 * t82, t26, 0.2e1 * t1 * t30 + 0.2e1 * t8 * t9, 0.2e1 * t8 * t11 - 0.2e1 * t2 * t30, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t55, 0, -t54 * pkin(6), -t55 * pkin(6), 0, 0, 0, 0, t34, 0, -t30, 0, -t19, -t21, (-t30 * t50 - t34 * t52) * pkin(2), (-t19 * t52 + t21 * t50) * pkin(2), t63, (-t45 + t46) * t34, t72, -t63, t23, 0, -t19 * t51 + t49 * t58, t19 * t49 + t51 * t58, t59, t19 * t43 + t40 * t59, t11 * t35, -t7 + t74, t73, -t9 * t84, t18, 0, t14 * t30 + t37 * t9 - t8 * t84, t37 * t11 - t15 * t30 + t8 * t35, -t1 * t35 - t14 * t11 - t15 * t9 + t2 * t84, t1 * t14 + t2 * t15 + t8 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t77, -0.2e1 * t78, 0, (t50 ^ 2 + t52 ^ 2) * pkin(2) ^ 2, t45, 0.2e1 * t70, 0, t46, 0, 0, -0.2e1 * t43 * t51, 0.2e1 * t43 * t49, 0.2e1 * t66 * t40, t40 ^ 2 * t66 + t43 ^ 2, t29, 0.2e1 * t35 * t84, 0, t27, 0, 0, -t84 * t81, t35 * t81, -0.2e1 * t14 * t35 + 0.2e1 * t15 * t84, t14 ^ 2 + t15 ^ 2 + t37 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t34, 0, t44, 0, 0, 0, 0, 0, 0, t23, -t72, -t66 * t34, t60, 0, 0, 0, 0, 0, 0, t18, -t73, -t7 - t74, t1 * t84 + t2 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t84 + t15 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 + t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t69, 0, t19, 0, 0, 0, 0, 0, 0, t9, t11, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t49, 0, t43, 0, 0, 0, 0, 0, 0, -t84, t35, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, t30, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, t84, 0, t14, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
