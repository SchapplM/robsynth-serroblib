% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t33 = g(1) * t57 + g(2) * t54;
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t74 = t57 * t55;
t56 = cos(qJ(2));
t85 = t54 * t56;
t23 = t52 * t85 + t74;
t75 = t57 * t52;
t25 = t54 * t55 - t56 * t75;
t53 = sin(qJ(2));
t90 = g(3) * t53;
t101 = -g(1) * t25 + g(2) * t23 + t52 * t90;
t51 = qJ(3) + qJ(4);
t42 = sin(t51);
t43 = cos(t51);
t76 = t57 * t43;
t16 = t42 * t85 + t76;
t77 = t57 * t42;
t18 = t54 * t43 - t56 * t77;
t5 = -g(1) * t18 + g(2) * t16 + t42 * t90;
t44 = qJ(5) + t51;
t38 = sin(t44);
t39 = cos(t44);
t78 = t57 * t39;
t11 = t38 * t85 + t78;
t79 = t57 * t38;
t13 = t54 * t39 - t56 * t79;
t3 = -g(1) * t13 + g(2) * t11 + t38 * t90;
t60 = -g(3) * t56 + t33 * t53;
t58 = -pkin(9) - pkin(8);
t100 = pkin(4) * t42;
t96 = g(1) * t54;
t88 = t52 * pkin(3);
t87 = t53 * t57;
t86 = t53 * t58;
t84 = t56 * t57;
t28 = -pkin(5) * t38 - t100;
t21 = -t28 + t88;
t83 = t57 * t21;
t31 = t88 + t100;
t82 = t57 * t31;
t41 = qJ(6) + t44;
t35 = sin(t41);
t81 = t57 * t35;
t36 = cos(t41);
t80 = t57 * t36;
t37 = pkin(4) * t43;
t47 = t55 * pkin(3);
t32 = t37 + t47;
t73 = t57 * pkin(1) + t54 * pkin(7);
t50 = -pkin(10) + t58;
t34 = pkin(5) * t39;
t22 = t34 + t32;
t69 = t56 * pkin(2) + t53 * pkin(8);
t67 = -g(2) * t57 + t96;
t15 = pkin(2) + t22;
t45 = -pkin(11) + t50;
t66 = t56 * t15 - t53 * t45;
t30 = pkin(2) + t32;
t64 = t56 * t30 - t53 * t50;
t40 = t47 + pkin(2);
t62 = t56 * t40 - t86;
t48 = t57 * pkin(7);
t29 = t34 + t37;
t27 = t67 * t53;
t26 = t54 * t52 + t56 * t74;
t24 = -t55 * t85 + t75;
t20 = t33 * t56 + t90;
t19 = t54 * t42 + t56 * t76;
t17 = -t43 * t85 + t77;
t14 = t54 * t38 + t56 * t78;
t12 = -t39 * t85 + t79;
t10 = t54 * t35 + t56 * t80;
t9 = t54 * t36 - t56 * t81;
t8 = -t36 * t85 + t81;
t7 = t35 * t85 + t80;
t6 = g(1) * t19 - g(2) * t17 + t43 * t90;
t4 = g(1) * t14 - g(2) * t12 + t39 * t90;
t2 = g(1) * t10 - g(2) * t8 + t36 * t90;
t1 = -g(1) * t9 + g(2) * t7 + t35 * t90;
t46 = [0, 0, 0, 0, 0, 0, t67, t33, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t56, -t27, -t33, -g(1) * (-t54 * pkin(1) + t48) - g(2) * t73, 0, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t26, -g(1) * t23 - g(2) * t25, t27, -g(1) * t48 - g(2) * (t69 * t57 + t73) - (-pkin(1) - t69) * t96, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t19, -g(1) * t16 - g(2) * t18, t27, -g(1) * (pkin(3) * t75 + t48) - g(2) * (t40 * t84 - t57 * t86 + t73) + (-g(1) * (-pkin(1) - t62) - g(2) * t88) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t27, -g(1) * (t48 + t82) - g(2) * (t30 * t84 - t50 * t87 + t73) + (-g(1) * (-pkin(1) - t64) - g(2) * t31) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t27, -g(1) * (t48 + t83) - g(2) * (t15 * t84 - t45 * t87 + t73) + (-g(1) * (-pkin(1) - t66) - g(2) * t21) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t20, 0, 0, 0, 0, 0, 0, 0, 0, t60 * t55, -t60 * t52, -t20, -g(3) * t69 + t33 * (pkin(2) * t53 - pkin(8) * t56) 0, 0, 0, 0, 0, 0, t60 * t43, -t60 * t42, -t20, -g(3) * t62 + t33 * (t40 * t53 + t56 * t58) 0, 0, 0, 0, 0, 0, t60 * t39, -t60 * t38, -t20, -g(3) * t64 + t33 * (t30 * t53 + t50 * t56) 0, 0, 0, 0, 0, 0, t60 * t36, -t60 * t35, -t20, -g(3) * t66 + t33 * (t15 * t53 + t45 * t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, g(1) * t26 - g(2) * t24 + t55 * t90, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t101 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (t54 * t32 - t56 * t82) - g(2) * (-t31 * t85 - t57 * t32) + t31 * t90, 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t54 * t22 - t56 * t83) - g(2) * (-t21 * t85 - t57 * t22) + t21 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t28 * t84 + t54 * t29) - g(2) * (t28 * t85 - t57 * t29) - t28 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t46;
