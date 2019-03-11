% Calculate inertial parameters regressor of gravitation load for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t49 = qJ(4) + pkin(11);
t42 = cos(t49);
t55 = cos(qJ(4));
t46 = t55 * pkin(4);
t27 = pkin(5) * t42 + t46;
t25 = pkin(3) + t27;
t50 = qJ(2) + qJ(3);
t44 = sin(t50);
t45 = cos(t50);
t51 = -qJ(5) - pkin(9);
t48 = -pkin(10) + t51;
t105 = t45 * t25 - t44 * t48;
t39 = t46 + pkin(3);
t104 = t45 * t39 - t44 * t51;
t103 = t45 * pkin(3) + t44 * pkin(9);
t54 = sin(qJ(1));
t57 = cos(qJ(1));
t29 = g(1) * t57 + g(2) * t54;
t76 = t57 * t55;
t52 = sin(qJ(4));
t83 = t54 * t52;
t21 = t45 * t83 + t76;
t77 = t57 * t52;
t82 = t54 * t55;
t23 = -t45 * t77 + t82;
t94 = g(3) * t44;
t102 = -g(1) * t23 + g(2) * t21 + t52 * t94;
t9 = -g(3) * t45 + t29 * t44;
t53 = sin(qJ(2));
t101 = pkin(2) * t53;
t100 = pkin(3) * t44;
t56 = cos(qJ(2));
t47 = t56 * pkin(2);
t40 = t47 + pkin(1);
t30 = t57 * t40;
t96 = g(2) * t30;
t92 = t52 * pkin(4);
t89 = t45 * t54;
t88 = t45 * t57;
t43 = qJ(6) + t49;
t34 = sin(t43);
t87 = t54 * t34;
t35 = cos(t43);
t86 = t54 * t35;
t41 = sin(t49);
t85 = t54 * t41;
t84 = t54 * t42;
t81 = t57 * t34;
t80 = t57 * t35;
t79 = t57 * t41;
t78 = t57 * t42;
t26 = pkin(5) * t41 + t92;
t58 = -pkin(8) - pkin(7);
t75 = t26 - t58;
t72 = -t58 + t92;
t69 = -t100 - t101;
t67 = g(1) * t54 - g(2) * t57;
t65 = t25 * t44 + t45 * t48;
t63 = t39 * t44 + t45 * t51;
t59 = -g(3) * t56 + t29 * t53;
t32 = pkin(9) * t88;
t31 = pkin(9) * t89;
t24 = t45 * t76 + t83;
t22 = -t45 * t82 + t77;
t19 = t67 * t44;
t18 = t45 * t78 + t85;
t17 = -t45 * t79 + t84;
t16 = -t45 * t84 + t79;
t15 = t45 * t85 + t78;
t14 = t45 * t80 + t87;
t13 = -t45 * t81 + t86;
t12 = -t45 * t86 + t81;
t11 = t45 * t87 + t80;
t10 = t29 * t45 + t94;
t8 = t9 * t55;
t7 = t9 * t52;
t6 = t9 * t42;
t5 = t9 * t41;
t4 = t9 * t35;
t3 = t9 * t34;
t2 = g(1) * t14 - g(2) * t12 + t35 * t94;
t1 = -g(1) * t13 + g(2) * t11 + t34 * t94;
t20 = [0, 0, 0, 0, 0, 0, t67, t29, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t56, -t67 * t53, -t29, -g(1) * (-t54 * pkin(1) + t57 * pkin(7)) - g(2) * (t57 * pkin(1) + t54 * pkin(7)) 0, 0, 0, 0, 0, 0, t67 * t45, -t19, -t29, -g(1) * (-t54 * t40 - t57 * t58) - g(2) * (-t54 * t58 + t30) 0, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, -g(1) * t21 - g(2) * t23, t19, -t96 + (g(1) * t58 - g(2) * t103) * t57 + (-g(1) * (-t40 - t103) + g(2) * t58) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t19, -t96 + (-g(1) * t72 - g(2) * t104) * t57 + (-g(1) * (-t40 - t104) - g(2) * t72) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t19, -t96 + (-g(1) * t75 - g(2) * t105) * t57 + (-g(1) * (-t40 - t105) - g(2) * t75) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, g(3) * t53 + t29 * t56, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t59 * pkin(2), 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (t69 * t57 + t32) - g(2) * (t69 * t54 + t31) - g(3) * (t47 + t103) 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * (t47 + t104) + t29 * (t63 + t101) 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * (t47 + t105) + t29 * (t65 + t101); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t10, -g(1) * (-t57 * t100 + t32) - g(2) * (-t54 * t100 + t31) - g(3) * t103, 0, 0, 0, 0, 0, 0, t6, -t5, -t10, -g(3) * t104 + t29 * t63, 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * t105 + t29 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, g(1) * t24 - g(2) * t22 + t55 * t94, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + g(2) * t15 + t41 * t94, g(1) * t18 - g(2) * t16 + t42 * t94, 0, t102 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t26 * t88 + t54 * t27) - g(2) * (-t26 * t89 - t57 * t27) + t26 * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t20;
