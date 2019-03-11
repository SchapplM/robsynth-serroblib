% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t58 = cos(qJ(2));
t54 = sin(qJ(3));
t55 = sin(qJ(2));
t74 = t54 * t55;
t80 = cos(qJ(3));
t24 = -t80 * t58 + t74;
t66 = t80 * t55;
t25 = t54 * t58 + t66;
t53 = sin(qJ(5));
t57 = cos(qJ(5));
t13 = -t57 * t24 + t53 * t25;
t87 = 0.2e1 * t13;
t86 = -0.2e1 * t25;
t43 = -t58 * pkin(2) - pkin(1);
t85 = 0.2e1 * t43;
t56 = cos(qJ(6));
t84 = 0.2e1 * t56;
t83 = 0.2e1 * t58;
t82 = -pkin(8) - pkin(7);
t34 = t82 * t58;
t16 = -t80 * t34 + t82 * t74;
t10 = t24 * pkin(9) + t16;
t15 = -t54 * t34 - t82 * t66;
t9 = -t25 * pkin(9) + t15;
t5 = t53 * t10 - t57 * t9;
t52 = sin(qJ(6));
t81 = t5 * t52;
t47 = t54 * pkin(2);
t38 = t47 + qJ(4);
t32 = t53 * t38;
t67 = t80 * pkin(2);
t41 = t67 + pkin(3);
t35 = -pkin(4) - t41;
t19 = -t57 * t35 + t32;
t17 = pkin(5) + t19;
t79 = t17 * t52;
t45 = t53 * qJ(4);
t59 = -pkin(3) - pkin(4);
t29 = -t57 * t59 + t45;
t27 = pkin(5) + t29;
t78 = t27 * t52;
t77 = t52 * t13;
t14 = t53 * t24 + t57 * t25;
t76 = t52 * t14;
t75 = t52 * t56;
t73 = t56 * t13;
t72 = t56 * t14;
t71 = t57 * t56;
t70 = t17 + t27;
t20 = t53 * t35 + t57 * t38;
t30 = t57 * qJ(4) + t53 * t59;
t69 = -0.2e1 * t14 * t13;
t37 = -0.2e1 * t75;
t68 = t52 * t72;
t11 = t24 * pkin(3) - t25 * qJ(4) + t43;
t65 = -pkin(5) * t14 - pkin(10) * t13;
t18 = -pkin(10) + t20;
t64 = -t13 * t18 + t14 * t17;
t28 = -pkin(10) + t30;
t63 = -t13 * t28 + t14 * t27;
t62 = -t13 * t53 - t14 * t57;
t8 = -t24 * pkin(4) - t11;
t61 = 0.2e1 * pkin(3);
t60 = 0.2e1 * qJ(4);
t51 = t56 ^ 2;
t50 = t52 ^ 2;
t49 = pkin(5) * t52;
t39 = t57 * t52;
t36 = 0.2e1 * t75;
t12 = t14 ^ 2;
t7 = (t50 - t51) * t14;
t6 = t57 * t10 + t53 * t9;
t4 = t5 * t56;
t3 = t13 * pkin(5) - t14 * pkin(10) + t8;
t2 = t52 * t3 + t56 * t6;
t1 = t56 * t3 - t52 * t6;
t21 = [1, 0, 0, t55 ^ 2, t55 * t83, 0, 0, 0, pkin(1) * t83, -0.2e1 * pkin(1) * t55, t25 ^ 2, t24 * t86, 0, 0, 0, t24 * t85, t25 * t85, 0.2e1 * t11 * t24, 0.2e1 * t15 * t25 - 0.2e1 * t16 * t24, t11 * t86, t11 ^ 2 + t15 ^ 2 + t16 ^ 2, t12, t69, 0, 0, 0, t8 * t87, 0.2e1 * t8 * t14, t51 * t12, t12 * t37, t72 * t87, t52 * t69, t13 ^ 2, 0.2e1 * t1 * t13 + 0.2e1 * t5 * t76, -0.2e1 * t2 * t13 + 0.2e1 * t5 * t72; 0, 0, 0, 0, 0, t55, t58, 0, -t55 * pkin(7), -t58 * pkin(7), 0, 0, t25, -t24, 0, -t15, -t16, -t15, -t38 * t24 - t41 * t25, t16, -t15 * t41 + t16 * t38, 0, 0, -t14, t13, 0, t5, t6, -t68, t7, -t77, -t73, 0, t64 * t52 + t4, t64 * t56 - t81; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t67, -0.2e1 * t47, 0.2e1 * t41, 0, 0.2e1 * t38, t38 ^ 2 + t41 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t19, 0.2e1 * t20, t50, t36, 0, 0, 0, t17 * t84, -0.2e1 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, -t15, -t16, -t15, -pkin(3) * t25 - t24 * qJ(4), t16, -t15 * pkin(3) + t16 * qJ(4), 0, 0, -t14, t13, 0, t5, t6, -t68, t7, -t77, -t73, 0, t63 * t52 + t4, t63 * t56 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t67, -t47, t61 + t67, 0, t60 + t47, t41 * pkin(3) + t38 * qJ(4), 0, 0, 0, 0, 1, t32 + t45 + (-t35 - t59) * t57, t30 + t20, t50, t36, 0, 0, 0, t70 * t56, -t70 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t61, 0, t60, pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t29, 0.2e1 * t30, t50, t36, 0, 0, 0, t27 * t84, -0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t52, t62 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t41, 0, 0, 0, 0, 0, -t57, t53, 0, 0, 0, 0, 0, -t71, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, -t57, t53, 0, 0, 0, 0, 0, -t71, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, -t5, -t6, t68, -t7, t77, t73, 0, t65 * t52 - t4, t65 * t56 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t19, -t20, -t50, t37, 0, 0, 0 (-pkin(5) - t17) * t56, t49 + t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t29, -t30, -t50, t37, 0, 0, 0 (-pkin(5) - t27) * t56, t49 + t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t53, 0, 0, 0, 0, 0, t71, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t50, t36, 0, 0, 0, pkin(5) * t84, -0.2e1 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, -t76, t13, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t56, 0, -t52 * t18, -t56 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t56, 0, -t52 * t28, -t56 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t53, -t56 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t56, 0, -t52 * pkin(10), -t56 * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t21;
