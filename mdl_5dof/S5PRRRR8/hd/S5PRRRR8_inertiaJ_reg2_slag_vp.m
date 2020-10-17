% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:16:41
% EndTime: 2019-12-05 17:16:45
% DurationCPUTime: 0.86s
% Computational Cost: add. (550->99), mult. (1250->189), div. (0->0), fcn. (1467->10), ass. (0->70)
t48 = sin(qJ(4));
t77 = t48 * pkin(3);
t35 = pkin(9) + t77;
t47 = sin(qJ(5));
t41 = t47 ^ 2;
t51 = cos(qJ(5));
t43 = t51 ^ 2;
t64 = t41 + t43;
t66 = t64 * t35;
t46 = cos(pkin(5));
t49 = sin(qJ(3));
t53 = cos(qJ(3));
t45 = sin(pkin(5));
t50 = sin(qJ(2));
t71 = t45 * t50;
t19 = t46 * t53 - t49 * t71;
t20 = t46 * t49 + t53 * t71;
t52 = cos(qJ(4));
t8 = -t52 * t19 + t48 * t20;
t83 = t8 ^ 2;
t78 = -pkin(8) - pkin(7);
t28 = t78 * t53;
t61 = t78 * t49;
t14 = -t48 * t28 - t52 * t61;
t82 = t14 ^ 2;
t24 = t48 * t49 - t52 * t53;
t81 = t24 ^ 2;
t37 = -t53 * pkin(3) - pkin(2);
t80 = 0.2e1 * t37;
t79 = 0.2e1 * t53;
t76 = t52 * pkin(3);
t75 = t8 * t14;
t74 = t8 * t51;
t36 = -pkin(4) - t76;
t73 = pkin(4) - t36;
t72 = t14 * t51;
t54 = cos(qJ(2));
t70 = t45 * t54;
t26 = t48 * t53 + t52 * t49;
t69 = t47 * t26;
t68 = t47 * t51;
t67 = t51 * t26;
t65 = t64 * pkin(9);
t42 = t49 ^ 2;
t44 = t53 ^ 2;
t63 = t42 + t44;
t62 = -0.2e1 * t26 * t24;
t60 = -pkin(4) * t26 - pkin(9) * t24;
t11 = t24 * pkin(4) - t26 * pkin(9) + t37;
t16 = -t52 * t28 + t48 * t61;
t3 = t51 * t11 - t47 * t16;
t4 = t47 * t11 + t51 * t16;
t1 = -t3 * t47 + t4 * t51;
t10 = t48 * t19 + t52 * t20;
t5 = -t47 * t10 - t51 * t70;
t6 = t51 * t10 - t47 * t70;
t2 = -t5 * t47 + t6 * t51;
t59 = -t19 * t49 + t20 * t53;
t58 = -t24 * t35 + t26 * t36;
t40 = t45 ^ 2;
t33 = t40 * t54 ^ 2;
t31 = 0.2e1 * t68;
t23 = t26 ^ 2;
t22 = t51 * t24;
t21 = t47 * t24;
t18 = t47 * t67;
t13 = t14 * t47;
t12 = (-t41 + t43) * t26;
t7 = t8 * t47;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t50 ^ 2 + t46 ^ 2 + t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 + t20 ^ 2 + t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t33 + t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t6 ^ 2 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t71, 0, 0, 0, 0, 0, 0, 0, 0, t53 * t70, -t49 * t70, t59, pkin(2) * t70 + t59 * pkin(7), 0, 0, 0, 0, 0, 0, -t24 * t70, -t26 * t70, -t10 * t24 + t8 * t26, t10 * t16 - t37 * t70 + t75, 0, 0, 0, 0, 0, 0, t5 * t24 + t8 * t69, -t6 * t24 + t8 * t67, (-t47 * t6 - t5 * t51) * t26, t5 * t3 + t6 * t4 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t42, t49 * t79, 0, t44, 0, 0, pkin(2) * t79, -0.2e1 * pkin(2) * t49, 0.2e1 * t63 * pkin(7), t63 * pkin(7) ^ 2 + pkin(2) ^ 2, t23, t62, 0, t81, 0, 0, t24 * t80, t26 * t80, 0.2e1 * t14 * t26 - 0.2e1 * t16 * t24, t16 ^ 2 + t37 ^ 2 + t82, t43 * t23, -0.2e1 * t23 * t68, 0.2e1 * t24 * t67, t41 * t23, t47 * t62, t81, 0.2e1 * t14 * t69 + 0.2e1 * t3 * t24, 0.2e1 * t14 * t67 - 0.2e1 * t4 * t24, 0.2e1 * (-t3 * t51 - t4 * t47) * t26, t3 ^ 2 + t4 ^ 2 + t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t20, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t10, 0, (t10 * t48 - t52 * t8) * pkin(3), 0, 0, 0, 0, 0, 0, -t74, t7, t2, t2 * t35 + t8 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t53, 0, -t49 * pkin(7), -t53 * pkin(7), 0, 0, 0, 0, t26, 0, -t24, 0, -t14, -t16, (-t24 * t48 - t26 * t52) * pkin(3), (-t14 * t52 + t16 * t48) * pkin(3), t18, t12, t21, -t18, t22, 0, t58 * t47 - t72, t58 * t51 + t13, t1, t1 * t35 + t14 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t76, -0.2e1 * t77, 0, (t48 ^ 2 + t52 ^ 2) * pkin(3) ^ 2, t41, t31, 0, t43, 0, 0, -0.2e1 * t36 * t51, 0.2e1 * t36 * t47, 0.2e1 * t66, t64 * t35 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t10, 0, 0, 0, 0, 0, 0, 0, 0, -t74, t7, t2, -t8 * pkin(4) + t2 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, 0, -t14, -t16, 0, 0, t18, t12, t21, -t18, t22, 0, t60 * t47 - t72, t60 * t51 + t13, t1, -t14 * pkin(4) + t1 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t76, -t77, 0, 0, t41, t31, 0, t43, 0, 0, t73 * t51, -t73 * t47, t65 + t66, -t36 * pkin(4) + pkin(9) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t41, t31, 0, t43, 0, 0, 0.2e1 * pkin(4) * t51, -0.2e1 * pkin(4) * t47, 0.2e1 * t65, t64 * pkin(9) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, -t69, t24, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t51, 0, -t47 * t35, -t51 * t35, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t51, 0, -t47 * pkin(9), -t51 * pkin(9), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t9;
