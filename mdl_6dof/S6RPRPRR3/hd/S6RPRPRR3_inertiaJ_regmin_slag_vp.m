% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:34:56
% EndTime: 2019-05-05 18:34:57
% DurationCPUTime: 0.59s
% Computational Cost: add. (680->104), mult. (1392->178), div. (0->0), fcn. (1602->10), ass. (0->73)
t49 = sin(pkin(11));
t51 = cos(pkin(11));
t54 = sin(qJ(5));
t56 = cos(qJ(5));
t34 = t54 * t49 - t56 * t51;
t44 = -t51 * pkin(4) - pkin(3);
t25 = t34 * pkin(5) + t44;
t85 = 0.2e1 * t25;
t84 = 0.2e1 * t44;
t55 = sin(qJ(3));
t83 = 0.2e1 * t55;
t57 = cos(qJ(3));
t82 = -0.2e1 * t57;
t81 = 0.2e1 * t57;
t53 = sin(qJ(6));
t80 = t53 * pkin(5);
t79 = t57 * pkin(3);
t78 = t57 * pkin(5);
t77 = cos(qJ(6));
t35 = t56 * t49 + t54 * t51;
t19 = -t53 * t34 + t77 * t35;
t76 = t19 * t57;
t75 = t35 * t57;
t50 = sin(pkin(10));
t42 = t50 * pkin(1) + pkin(7);
t74 = t42 * t49;
t73 = t49 * t55;
t72 = t51 * t55;
t39 = t55 * t42;
t18 = t77 * t34 + t53 * t35;
t71 = t57 * t18;
t70 = t57 * t34;
t69 = t57 * t42;
t68 = t57 * t49;
t67 = t57 * t51;
t66 = pkin(8) + qJ(4);
t52 = cos(pkin(10));
t43 = -t52 * pkin(1) - pkin(2);
t32 = -t55 * qJ(4) + t43 - t79;
t21 = t49 * t32 + t42 * t67;
t30 = pkin(4) * t73 + t39;
t65 = t49 ^ 2 + t51 ^ 2;
t64 = t77 * pkin(5);
t26 = t35 * t55;
t29 = t51 * t32;
t13 = -pkin(8) * t72 + t29 + (-pkin(4) - t74) * t57;
t16 = -pkin(8) * t73 + t21;
t9 = t54 * t13 + t56 * t16;
t5 = -t26 * pkin(9) + t9;
t63 = t77 * t5;
t27 = t34 * t55;
t8 = t56 * t13 - t54 * t16;
t4 = t27 * pkin(9) - t78 + t8;
t1 = t77 * t4 - t53 * t5;
t37 = t66 * t49;
t38 = t66 * t51;
t22 = -t56 * t37 - t54 * t38;
t62 = t65 * qJ(4);
t61 = -pkin(3) * t55 + qJ(4) * t57;
t20 = -t42 * t68 + t29;
t60 = -t20 * t49 + t21 * t51;
t23 = -t54 * t37 + t56 * t38;
t48 = t57 ^ 2;
t47 = t55 ^ 2;
t17 = t26 * pkin(5) + t30;
t15 = -t34 * pkin(9) + t23;
t14 = -t35 * pkin(9) + t22;
t12 = -t53 * t26 - t77 * t27;
t11 = t77 * t26 - t53 * t27;
t7 = t53 * t14 + t77 * t15;
t6 = t77 * t14 - t53 * t15;
t2 = t53 * t4 + t63;
t3 = [1, 0, 0 (t50 ^ 2 + t52 ^ 2) * pkin(1) ^ 2, t47, t55 * t81, 0, 0, 0, t43 * t82, t43 * t83, -0.2e1 * t20 * t57 + 0.2e1 * t47 * t74, 0.2e1 * t47 * t42 * t51 + 0.2e1 * t21 * t57 (-t20 * t51 - t21 * t49) * t83, t47 * t42 ^ 2 + t20 ^ 2 + t21 ^ 2, t27 ^ 2, 0.2e1 * t27 * t26, -t27 * t82, t26 * t81, t48, 0.2e1 * t30 * t26 - 0.2e1 * t8 * t57, -0.2e1 * t30 * t27 + 0.2e1 * t9 * t57, t12 ^ 2, -0.2e1 * t12 * t11, t12 * t82, t11 * t81, t48, -0.2e1 * t1 * t57 + 0.2e1 * t17 * t11, 0.2e1 * t17 * t12 + 0.2e1 * t2 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t60 - t69) * t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t47 + t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t55, t57, 0, -t39, -t69, -t51 * t39 + t61 * t49, t49 * t39 + t61 * t51, t60, -pkin(3) * t39 + t60 * qJ(4), -t27 * t35, -t35 * t26 + t27 * t34, -t75, t70, 0, -t22 * t57 + t44 * t26 + t30 * t34, t23 * t57 - t44 * t27 + t30 * t35, t12 * t19, -t19 * t11 - t12 * t18, -t76, t71, 0, t25 * t11 + t17 * t18 - t6 * t57, t25 * t12 + t17 * t19 + t7 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t55, t67, -t68, t65 * t55, t55 * t62 + t79, 0, 0, 0, 0, 0, -t70, -t75, 0, 0, 0, 0, 0, -t71, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t51, -0.2e1 * pkin(3) * t49, 0.2e1 * t62, t65 * qJ(4) ^ 2 + pkin(3) ^ 2, t35 ^ 2, -0.2e1 * t35 * t34, 0, 0, 0, t34 * t84, t35 * t84, t19 ^ 2, -0.2e1 * t19 * t18, 0, 0, 0, t18 * t85, t19 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t72, 0, t39, 0, 0, 0, 0, 0, t26, -t27, 0, 0, 0, 0, 0, t11, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t49, 0, -pkin(3), 0, 0, 0, 0, 0, t34, t35, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t26, -t57, t8, -t9, 0, 0, t12, -t11, -t57, -t57 * t64 + t1, -t63 + (-t4 + t78) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t27, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, t22, -t23, 0, 0, t19, -t18, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t64, -0.2e1 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, -t57, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t64, -t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
