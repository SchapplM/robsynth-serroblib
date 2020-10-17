% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:32:41
% EndTime: 2019-05-05 04:32:43
% DurationCPUTime: 0.58s
% Computational Cost: add. (584->101), mult. (1282->187), div. (0->0), fcn. (1618->12), ass. (0->73)
t47 = sin(pkin(12));
t49 = cos(pkin(12));
t52 = sin(qJ(3));
t73 = cos(qJ(3));
t33 = t47 * t52 - t49 * t73;
t81 = -0.2e1 * t33;
t80 = 0.2e1 * t33;
t42 = -t49 * pkin(3) - pkin(4);
t55 = cos(qJ(5));
t38 = -t55 * pkin(5) + t42;
t79 = 0.2e1 * t38;
t78 = t33 * pkin(5);
t50 = sin(qJ(6));
t77 = t50 * pkin(5);
t54 = cos(qJ(6));
t76 = t54 * pkin(5);
t34 = t47 * t73 + t49 * t52;
t44 = -t73 * pkin(3) - pkin(2);
t21 = t33 * pkin(4) - t34 * pkin(9) + t44;
t51 = sin(qJ(5));
t62 = t73 * pkin(8);
t39 = t73 * qJ(4) + t62;
t60 = (-qJ(4) - pkin(8)) * t52;
t24 = t49 * t39 + t47 * t60;
t66 = t55 * t24;
t7 = t66 + (-pkin(10) * t34 + t21) * t51;
t75 = t54 * t7;
t41 = t47 * pkin(3) + pkin(9);
t74 = pkin(10) + t41;
t37 = t50 * t55 + t54 * t51;
t72 = t37 * t33;
t48 = sin(pkin(6));
t53 = sin(qJ(2));
t71 = t48 * t53;
t56 = cos(qJ(2));
t70 = t48 * t56;
t69 = t51 * t33;
t68 = t51 * t34;
t67 = t51 * t55;
t65 = t55 * t34;
t64 = cos(pkin(6));
t63 = 0.2e1 * t73;
t8 = t55 * t21 - t51 * t24;
t6 = -pkin(10) * t65 + t78 + t8;
t1 = -t50 * t7 + t54 * t6;
t61 = t48 * t73;
t22 = t47 * t39 - t49 * t60;
t59 = -t33 * t41 + t34 * t42;
t36 = t50 * t51 - t54 * t55;
t58 = -t52 * t71 + t64 * t73;
t46 = t55 ^ 2;
t45 = t51 ^ 2;
t32 = t34 ^ 2;
t31 = t33 ^ 2;
t30 = t74 * t55;
t29 = t74 * t51;
t28 = t64 * t52 + t53 * t61;
t27 = t55 * t33;
t25 = t36 * t33;
t20 = -t50 * t29 + t54 * t30;
t19 = -t54 * t29 - t50 * t30;
t17 = t36 * t34;
t16 = t37 * t34;
t15 = t49 * t28 + t47 * t58;
t13 = t47 * t28 - t49 * t58;
t12 = pkin(5) * t68 + t22;
t11 = t55 * t15 - t51 * t70;
t10 = -t51 * t15 - t55 * t70;
t9 = t51 * t21 + t66;
t4 = t50 * t10 + t54 * t11;
t3 = t54 * t10 - t50 * t11;
t2 = t50 * t6 + t75;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 ^ 2 * t56 ^ 2 + t13 ^ 2 + t15 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t70, -t71, 0, 0, 0, 0, 0, t56 * t61, -t52 * t70, t13 * t34 - t15 * t33, t13 * t22 + t15 * t24 - t44 * t70, 0, 0, 0, 0, 0, t10 * t33 + t13 * t68, -t11 * t33 + t13 * t65, 0, 0, 0, 0, 0, t13 * t16 + t3 * t33, -t13 * t17 - t4 * t33; 0, 1, 0, 0, t52 ^ 2, t52 * t63, 0, 0, 0, pkin(2) * t63, -0.2e1 * pkin(2) * t52, 0.2e1 * t22 * t34 - 0.2e1 * t24 * t33, t22 ^ 2 + t24 ^ 2 + t44 ^ 2, t46 * t32, -0.2e1 * t32 * t67, t65 * t80, t68 * t81, t31, 0.2e1 * t22 * t68 + 0.2e1 * t8 * t33, 0.2e1 * t22 * t65 - 0.2e1 * t9 * t33, t17 ^ 2, 0.2e1 * t17 * t16, -t17 * t80, t16 * t81, t31, 0.2e1 * t1 * t33 + 0.2e1 * t12 * t16, -0.2e1 * t12 * t17 - 0.2e1 * t2 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, -t28, 0 (-t13 * t49 + t15 * t47) * pkin(3), 0, 0, 0, 0, 0, -t13 * t55, t13 * t51, 0, 0, 0, 0, 0, t13 * t36, t13 * t37; 0, 0, 0, 0, 0, 0, t52, t73, 0, -t52 * pkin(8), -t62 (-t33 * t47 - t34 * t49) * pkin(3) (-t22 * t49 + t24 * t47) * pkin(3), t51 * t65 (-t45 + t46) * t34, t69, t27, 0, -t22 * t55 + t59 * t51, t22 * t51 + t59 * t55, -t17 * t37, -t37 * t16 + t17 * t36, t72, -t25, 0, t12 * t36 + t38 * t16 + t19 * t33, t12 * t37 - t38 * t17 - t20 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t47 ^ 2 + t49 ^ 2) * pkin(3) ^ 2, t45, 0.2e1 * t67, 0, 0, 0, -0.2e1 * t42 * t55, 0.2e1 * t42 * t51, t37 ^ 2, -0.2e1 * t37 * t36, 0, 0, 0, t36 * t79, t37 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, t27, -t69, 0, 0, 0, 0, 0, -t25, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t11, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, -t68, t33, t8, -t9, 0, 0, -t17, -t16, t33, t33 * t76 + t1, -t75 + (-t6 - t78) * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, t55, 0, -t51 * t41, -t55 * t41, 0, 0, t37, -t36, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, -t51, 0, 0, 0, 0, 0, -t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t76, -0.2e1 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t16, t33, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, 0, t19, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t76, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
