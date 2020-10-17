% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:01:13
% EndTime: 2019-05-05 07:01:15
% DurationCPUTime: 0.57s
% Computational Cost: add. (707->85), mult. (1496->171), div. (0->0), fcn. (1888->12), ass. (0->74)
t52 = sin(qJ(4));
t53 = sin(qJ(3));
t56 = cos(qJ(4));
t57 = cos(qJ(3));
t35 = t52 * t57 + t56 * t53;
t82 = 0.2e1 * t35;
t51 = sin(qJ(6));
t81 = 0.2e1 * t51;
t55 = cos(qJ(6));
t80 = -0.2e1 * t55;
t79 = 0.2e1 * t57;
t78 = pkin(8) + pkin(9);
t77 = t52 * pkin(3);
t41 = t52 * t53;
t67 = t78 * t53;
t64 = t52 * t67;
t69 = t56 * t57;
t17 = -t41 * qJ(5) + (qJ(5) + t78) * t69 - t64;
t47 = sin(pkin(12));
t49 = cos(pkin(12));
t66 = t78 * t57;
t25 = -t52 * t66 - t56 * t67;
t60 = -t35 * qJ(5) + t25;
t8 = t47 * t17 - t49 * t60;
t76 = t8 * t55;
t50 = cos(pkin(6));
t48 = sin(pkin(6));
t74 = t48 * sin(qJ(2));
t33 = t50 * t57 - t53 * t74;
t34 = t50 * t53 + t57 * t74;
t18 = t52 * t33 + t56 * t34;
t61 = t56 * t33 - t52 * t34;
t11 = t47 * t18 - t49 * t61;
t75 = t11 * t55;
t58 = cos(qJ(2));
t73 = t48 * t58;
t65 = -t41 + t69;
t23 = t47 * t35 - t49 * t65;
t20 = t51 * t23;
t24 = t49 * t35 + t47 * t65;
t72 = t51 * t24;
t71 = t51 * t55;
t70 = t55 * t24;
t44 = t56 * pkin(3);
t42 = t44 + pkin(4);
t30 = t49 * t42 - t47 * t77;
t28 = -pkin(5) - t30;
t40 = -t49 * pkin(4) - pkin(5);
t68 = t28 + t40;
t31 = t47 * t42 + t49 * t77;
t43 = -t57 * pkin(3) - pkin(2);
t29 = pkin(10) + t31;
t63 = -t23 * t29 + t24 * t28;
t39 = t47 * pkin(4) + pkin(10);
t62 = -t23 * t39 + t24 * t40;
t27 = -t65 * pkin(4) + t43;
t46 = t55 ^ 2;
t45 = t51 ^ 2;
t38 = 0.2e1 * t71;
t26 = -t56 * t66 + t64;
t22 = t24 ^ 2;
t21 = t55 * t23;
t19 = t51 * t70;
t14 = (-t45 + t46) * t24;
t13 = t49 * t18 + t47 * t61;
t10 = t49 * t17 + t47 * t60;
t7 = t23 * pkin(5) - t24 * pkin(10) + t27;
t6 = t11 * t51;
t5 = t8 * t51;
t4 = t55 * t13 - t51 * t73;
t3 = -t51 * t13 - t55 * t73;
t2 = t55 * t10 + t51 * t7;
t1 = -t51 * t10 + t55 * t7;
t9 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 ^ 2 * t58 ^ 2 + t11 ^ 2 + t13 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, t73, -t74, 0, 0, 0, 0, 0, t57 * t73, -t53 * t73, 0, 0, 0, 0, 0, t65 * t73, -t35 * t73, t11 * t24 - t13 * t23, t13 * t10 + t11 * t8 - t27 * t73, 0, 0, 0, 0, 0, t11 * t72 + t3 * t23, t11 * t70 - t4 * t23; 0, 1, 0, 0, t53 ^ 2, t53 * t79, 0, 0, 0, pkin(2) * t79, -0.2e1 * pkin(2) * t53, t35 ^ 2, t65 * t82, 0, 0, 0, -0.2e1 * t43 * t65, t43 * t82, -0.2e1 * t10 * t23 + 0.2e1 * t8 * t24, t10 ^ 2 + t27 ^ 2 + t8 ^ 2, t46 * t22, -0.2e1 * t22 * t71, 0.2e1 * t23 * t70, -0.2e1 * t23 * t72, t23 ^ 2, 0.2e1 * t1 * t23 + 0.2e1 * t8 * t72, -0.2e1 * t2 * t23 + 0.2e1 * t8 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t34, 0, 0, 0, 0, 0, t61, -t18, 0, -t11 * t30 + t13 * t31, 0, 0, 0, 0, 0, -t75, t6; 0, 0, 0, 0, 0, 0, t53, t57, 0, -t53 * pkin(8), -t57 * pkin(8), 0, 0, t35, t65, 0, t25, t26, -t31 * t23 - t30 * t24, t10 * t31 - t8 * t30, t19, t14, t20, t21, 0, t63 * t51 - t76, t63 * t55 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t77, 0, t30 ^ 2 + t31 ^ 2, t45, t38, 0, 0, 0, t28 * t80, t28 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t18, 0 (-t11 * t49 + t13 * t47) * pkin(4), 0, 0, 0, 0, 0, -t75, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, t65, 0, t25, t26 (-t23 * t47 - t24 * t49) * pkin(4) (t10 * t47 - t49 * t8) * pkin(4), t19, t14, t20, t21, 0, t62 * t51 - t76, t62 * t55 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t44, -t77, 0 (t30 * t49 + t31 * t47) * pkin(4), t45, t38, 0, 0, 0, -t68 * t55, t68 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t47 ^ 2 + t49 ^ 2) * pkin(4) ^ 2, t45, t38, 0, 0, 0, t40 * t80, t40 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t73, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, t21, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t72, t23, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t55, 0, -t51 * t29, -t55 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t55, 0, -t51 * t39, -t55 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
