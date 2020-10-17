% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:34
% EndTime: 2019-12-05 16:37:38
% DurationCPUTime: 0.83s
% Computational Cost: add. (493->117), mult. (1202->244), div. (0->0), fcn. (1351->10), ass. (0->79)
t43 = cos(pkin(5));
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t41 = sin(pkin(5));
t46 = sin(qJ(2));
t76 = t41 * t46;
t17 = t43 * t45 + t48 * t76;
t40 = sin(pkin(10));
t42 = cos(pkin(10));
t49 = cos(qJ(2));
t75 = t41 * t49;
t5 = t17 * t40 + t42 * t75;
t85 = t5 ^ 2;
t15 = -t43 * t48 + t45 * t76;
t84 = t15 ^ 2;
t47 = cos(qJ(5));
t44 = sin(qJ(5));
t70 = t44 * t45;
t18 = t42 * t70 + t47 * t48;
t83 = -0.2e1 * t18;
t82 = -0.2e1 * t40;
t37 = t45 ^ 2;
t81 = t37 * pkin(7);
t32 = t45 * pkin(7);
t80 = t15 * t45;
t33 = t40 ^ 2;
t79 = t33 * t47;
t78 = t40 * t42;
t26 = t40 * t45;
t77 = t40 * t48;
t23 = -t48 * pkin(3) - t45 * qJ(4) - pkin(2);
t74 = t42 * t23;
t73 = t42 * t45;
t72 = t42 * t48;
t71 = t44 * t40;
t69 = t45 * t48;
t68 = t47 * t40;
t67 = t47 * t42;
t13 = pkin(7) * t72 + t40 * t23;
t36 = t44 ^ 2;
t38 = t47 ^ 2;
t66 = t36 + t38;
t65 = qJ(4) * t40;
t64 = qJ(4) * t42;
t63 = t33 * qJ(4);
t62 = -0.2e1 * t78;
t61 = 0.2e1 * t78;
t60 = 0.2e1 * t69;
t59 = t40 * t73;
t7 = t17 * t42 - t40 * t75;
t1 = t15 * t47 - t7 * t44;
t2 = t15 * t44 + t7 * t47;
t58 = t1 * t47 + t2 * t44;
t14 = t32 + (pkin(4) * t40 - pkin(8) * t42) * t45;
t9 = -t48 * pkin(8) + t13;
t3 = t47 * t14 - t44 * t9;
t4 = t44 * t14 + t47 * t9;
t57 = t3 * t47 + t4 * t44;
t56 = t5 * t40 + t7 * t42;
t55 = -pkin(3) * t45 + qJ(4) * t48;
t22 = -t42 * pkin(4) - t40 * pkin(8) - pkin(3);
t10 = t47 * t22 - t44 * t64;
t11 = t44 * t22 + t47 * t64;
t54 = t10 * t47 + t11 * t44;
t12 = -pkin(7) * t77 + t74;
t53 = -t12 * t40 + t13 * t42;
t52 = t17 * t48 + t80;
t51 = pkin(7) ^ 2;
t50 = qJ(4) ^ 2;
t39 = t48 ^ 2;
t35 = t42 ^ 2;
t34 = t41 ^ 2;
t31 = t37 * t51;
t30 = t33 * t50;
t27 = t34 * t49 ^ 2;
t25 = t33 * t37;
t20 = -t44 * t48 + t45 * t67;
t8 = -t74 + (pkin(7) * t40 + pkin(4)) * t48;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t46 ^ 2 + t43 ^ 2 + t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 + t27 + t84, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t84 + t85, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t76, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t75, -t45 * t75, t52, pkin(2) * t75 + pkin(7) * t52, 0, 0, 0, 0, 0, 0, t15 * t26 + t5 * t48, t15 * t73 + t7 * t48, (-t40 * t7 + t42 * t5) * t45, pkin(7) * t80 - t5 * t12 + t7 * t13, 0, 0, 0, 0, 0, 0, t1 * t26 + t5 * t18, -t2 * t26 + t5 * t20, -t1 * t20 - t2 * t18, t1 * t3 + t2 * t4 + t5 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, t60, 0, t39, 0, 0, 0.2e1 * pkin(2) * t48, -0.2e1 * pkin(2) * t45, 0.2e1 * (t37 + t39) * pkin(7), pkin(2) ^ 2 + t39 * t51 + t31, t35 * t37, t37 * t62, -0.2e1 * t42 * t69, t25, t40 * t60, t39, -0.2e1 * t12 * t48 + 0.2e1 * t40 * t81, 0.2e1 * t13 * t48 + 0.2e1 * t42 * t81, 0.2e1 * (-t12 * t42 - t13 * t40) * t45, t12 ^ 2 + t13 ^ 2 + t31, t20 ^ 2, t20 * t83, 0.2e1 * t20 * t26, t18 ^ 2, t26 * t83, t25, 0.2e1 * t8 * t18 + 0.2e1 * t26 * t3, 0.2e1 * t8 * t20 - 0.2e1 * t26 * t4, -0.2e1 * t4 * t18 - 0.2e1 * t3 * t20, t3 ^ 2 + t4 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t15 * t42, t15 * t40, t56, -t15 * pkin(3) + qJ(4) * t56, 0, 0, 0, 0, 0, 0, -t1 * t42 + t5 * t71, t2 * t42 + t5 * t68, -t58 * t40, t1 * t10 + t2 * t11 + t5 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, t48, 0, -t32, -t48 * pkin(7), 0, 0, t59, (-t33 + t35) * t45, -t77, -t59, -t72, 0, -pkin(7) * t73 + t55 * t40, pkin(7) * t26 + t42 * t55, t53, -pkin(3) * t32 + qJ(4) * t53, t20 * t68, (-t18 * t47 - t20 * t44) * t40, -t20 * t42 + t45 * t79, t18 * t71, t18 * t42 - t33 * t70, -t59, -t3 * t42 + (qJ(4) * t18 + t10 * t45 + t44 * t8) * t40, t4 * t42 + (qJ(4) * t20 - t11 * t45 + t47 * t8) * t40, -t10 * t20 - t11 * t18 - t57 * t40, t3 * t10 + t4 * t11 + t65 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t33, t61, 0, t35, 0, 0, 0.2e1 * pkin(3) * t42, pkin(3) * t82, 0.2e1 * (t33 + t35) * qJ(4), pkin(3) ^ 2 + t35 * t50 + t30, t38 * t33, -0.2e1 * t44 * t79, t47 * t62, t36 * t33, t44 * t61, t35, -0.2e1 * t10 * t42 + 0.2e1 * t44 * t63, 0.2e1 * t11 * t42 + 0.2e1 * t47 * t63, t54 * t82, t10 ^ 2 + t11 ^ 2 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t73, 0, t32, 0, 0, 0, 0, 0, 0, t45 * t68, -t40 * t70, -t44 * t18 - t47 * t20, t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, t40, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t67, t44 * t42, -t66 * t40, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, t26, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, -t71, -t42, t10, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
