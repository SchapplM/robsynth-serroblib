% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRPRR5
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
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRPRR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:26:27
% EndTime: 2019-05-05 05:26:28
% DurationCPUTime: 0.63s
% Computational Cost: add. (686->115), mult. (1602->218), div. (0->0), fcn. (1950->12), ass. (0->73)
t55 = sin(pkin(12));
t57 = cos(pkin(12));
t60 = sin(qJ(5));
t64 = cos(qJ(5));
t40 = t60 * t55 - t64 * t57;
t49 = -t57 * pkin(4) - pkin(3);
t31 = t40 * pkin(5) + t49;
t86 = 0.2e1 * t31;
t85 = 0.2e1 * t49;
t65 = cos(qJ(3));
t84 = -0.2e1 * t65;
t83 = 0.2e1 * t65;
t82 = pkin(8) * t55;
t59 = sin(qJ(6));
t81 = t59 * pkin(5);
t61 = sin(qJ(3));
t50 = t61 * pkin(8);
t63 = cos(qJ(6));
t80 = t63 * pkin(5);
t43 = -t65 * pkin(3) - t61 * qJ(4) - pkin(2);
t38 = t57 * t43;
t73 = t57 * t61;
t21 = -pkin(9) * t73 + t38 + (-pkin(4) - t82) * t65;
t77 = t65 * pkin(8);
t29 = t55 * t43 + t57 * t77;
t76 = t55 * t61;
t25 = -pkin(9) * t76 + t29;
t11 = t60 * t21 + t64 * t25;
t41 = t64 * t55 + t60 * t57;
t32 = t41 * t61;
t9 = -t32 * pkin(10) + t11;
t79 = t63 * t9;
t78 = t65 * pkin(5);
t56 = sin(pkin(6));
t75 = t56 * sin(qJ(2));
t74 = t56 * cos(qJ(2));
t72 = pkin(9) + qJ(4);
t42 = pkin(4) * t76 + t50;
t71 = t55 ^ 2 + t57 ^ 2;
t10 = t64 * t21 - t60 * t25;
t33 = t40 * t61;
t6 = t33 * pkin(10) + t10 - t78;
t1 = -t59 * t9 + t63 * t6;
t44 = t72 * t55;
t45 = t72 * t57;
t26 = -t64 * t44 - t60 * t45;
t70 = -pkin(3) * t61 + qJ(4) * t65;
t58 = cos(pkin(6));
t36 = t58 * t61 + t65 * t75;
t23 = -t36 * t55 - t57 * t74;
t24 = t36 * t57 - t55 * t74;
t69 = -t23 * t55 + t24 * t57;
t28 = -t55 * t77 + t38;
t68 = -t28 * t55 + t29 * t57;
t27 = -t60 * t44 + t64 * t45;
t54 = t65 ^ 2;
t53 = t61 ^ 2;
t35 = -t58 * t65 + t61 * t75;
t22 = t32 * pkin(5) + t42;
t20 = -t59 * t40 + t63 * t41;
t19 = t63 * t40 + t59 * t41;
t17 = -t40 * pkin(10) + t27;
t16 = -t41 * pkin(10) + t26;
t15 = -t59 * t32 - t63 * t33;
t14 = t63 * t32 - t59 * t33;
t13 = t60 * t23 + t64 * t24;
t12 = t64 * t23 - t60 * t24;
t8 = t59 * t16 + t63 * t17;
t7 = t63 * t16 - t59 * t17;
t4 = t59 * t12 + t63 * t13;
t3 = t63 * t12 - t59 * t13;
t2 = t59 * t6 + t79;
t5 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 ^ 2 + t24 ^ 2 + t35 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t74, -t75, 0, 0, 0, 0, 0, t65 * t74, -t61 * t74, -t23 * t65 + t35 * t76, t24 * t65 + t35 * t73 (-t23 * t57 - t24 * t55) * t61, t23 * t28 + t24 * t29 + t35 * t50, 0, 0, 0, 0, 0, -t12 * t65 + t35 * t32, t13 * t65 - t35 * t33, 0, 0, 0, 0, 0, t35 * t14 - t3 * t65, t35 * t15 + t4 * t65; 0, 1, 0, 0, t53, t61 * t83, 0, 0, 0, pkin(2) * t83, -0.2e1 * pkin(2) * t61, -0.2e1 * t28 * t65 + 0.2e1 * t53 * t82, 0.2e1 * t53 * pkin(8) * t57 + 0.2e1 * t29 * t65, 0.2e1 * (-t28 * t57 - t29 * t55) * t61, t53 * pkin(8) ^ 2 + t28 ^ 2 + t29 ^ 2, t33 ^ 2, 0.2e1 * t33 * t32, -t33 * t84, t32 * t83, t54, -0.2e1 * t10 * t65 + 0.2e1 * t42 * t32, 0.2e1 * t11 * t65 - 0.2e1 * t42 * t33, t15 ^ 2, -0.2e1 * t15 * t14, t15 * t84, t14 * t83, t54, -0.2e1 * t1 * t65 + 0.2e1 * t22 * t14, 0.2e1 * t22 * t15 + 0.2e1 * t2 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, -t35 * t57, t35 * t55, t69, -t35 * pkin(3) + t69 * qJ(4), 0, 0, 0, 0, 0, t35 * t40, t35 * t41, 0, 0, 0, 0, 0, t35 * t19, t35 * t20; 0, 0, 0, 0, 0, 0, t61, t65, 0, -t50, -t77, -pkin(8) * t73 + t70 * t55, pkin(8) * t76 + t70 * t57, t68, -pkin(3) * t50 + t68 * qJ(4), -t33 * t41, -t41 * t32 + t33 * t40, -t41 * t65, t40 * t65, 0, -t26 * t65 + t49 * t32 + t42 * t40, t27 * t65 - t49 * t33 + t42 * t41, t15 * t20, -t20 * t14 - t15 * t19, -t20 * t65, t19 * t65, 0, t31 * t14 + t22 * t19 - t7 * t65, t31 * t15 + t22 * t20 + t8 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t57, -0.2e1 * pkin(3) * t55, 0.2e1 * t71 * qJ(4), t71 * qJ(4) ^ 2 + pkin(3) ^ 2, t41 ^ 2, -0.2e1 * t41 * t40, 0, 0, 0, t40 * t85, t41 * t85, t20 ^ 2, -0.2e1 * t20 * t19, 0, 0, 0, t19 * t86, t20 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t73, 0, t50, 0, 0, 0, 0, 0, t32, -t33, 0, 0, 0, 0, 0, t14, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t55, 0, -pkin(3), 0, 0, 0, 0, 0, t40, t41, 0, 0, 0, 0, 0, t19, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t13, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t32, -t65, t10, -t11, 0, 0, t15, -t14, -t65, -t63 * t78 + t1, -t79 + (-t6 + t78) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t40, 0, t26, -t27, 0, 0, t20, -t19, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t80, -0.2e1 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, -t65, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, 0, t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t80, -t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t5;
