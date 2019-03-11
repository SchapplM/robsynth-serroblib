% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t44 = sin(qJ(3));
t43 = sin(qJ(4));
t48 = cos(qJ(1));
t72 = t48 * t43;
t45 = sin(qJ(1));
t46 = cos(qJ(4));
t78 = t45 * t46;
t16 = t44 * t72 + t78;
t71 = t48 * t46;
t79 = t45 * t43;
t14 = -t44 * t79 + t71;
t88 = g(2) * t48;
t89 = g(1) * t45;
t92 = t88 - t89;
t47 = cos(qJ(3));
t51 = -g(3) * t44 - t47 * t92;
t91 = -pkin(1) - pkin(7);
t90 = pkin(4) * t43;
t86 = g(3) * t47;
t42 = qJ(4) + qJ(5);
t35 = sin(t42);
t85 = t35 * t47;
t36 = cos(t42);
t84 = t36 * t47;
t83 = t44 * t45;
t82 = t44 * t48;
t81 = t45 * t35;
t80 = t45 * t36;
t49 = -pkin(9) - pkin(8);
t77 = t45 * t49;
t76 = t47 * t48;
t75 = t47 * t49;
t74 = t48 * t35;
t73 = t48 * t36;
t70 = t48 * t49;
t69 = t16 * pkin(4);
t68 = t48 * pkin(1) + t45 * qJ(2);
t67 = t43 * t86;
t34 = t46 * pkin(4) + pkin(3);
t38 = t48 * qJ(2);
t64 = t34 * t82 + t47 * t70 + t38;
t63 = t48 * pkin(7) + t68;
t10 = t44 * t80 + t74;
t9 = t44 * t81 - t73;
t62 = -t9 * pkin(5) + t10 * qJ(6);
t11 = t44 * t74 + t80;
t12 = t44 * t73 - t81;
t61 = t11 * pkin(5) - t12 * qJ(6);
t60 = g(2) * t63;
t59 = g(1) * t11 + g(2) * t9;
t57 = t44 * pkin(3) - t47 * pkin(8);
t24 = g(1) * t48 + g(2) * t45;
t56 = pkin(5) * t36 + qJ(6) * t35;
t55 = t14 * pkin(4);
t54 = pkin(4) * t72 + t34 * t83 + t45 * t75 + t63;
t53 = t34 + t56;
t52 = (-t90 + t91) * t89;
t1 = g(1) * t9 - g(2) * t11 + g(3) * t85;
t3 = g(1) * t10 - g(2) * t12 + g(3) * t84;
t27 = t44 * t70;
t25 = qJ(6) * t84;
t20 = t45 * t47 * t34;
t18 = t24 * t47;
t17 = t44 * t71 - t79;
t15 = t44 * t78 + t72;
t13 = g(1) * t83 - g(2) * t82 + t86;
t6 = t51 * t36;
t5 = t51 * t35;
t4 = -g(1) * t12 - g(2) * t10;
t2 = [0, 0, 0, 0, 0, 0, -t92, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t24, -g(1) * (-t45 * pkin(1) + t38) - g(2) * t68, 0, 0, 0, 0, 0, 0, -t24 * t44, -t18, -t92, -g(1) * (t45 * t91 + t38) - t60, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15, g(1) * t16 - g(2) * t14, t18, -g(1) * (pkin(3) * t82 - pkin(8) * t76 + t38) - t60 + (-g(1) * t91 - g(2) * t57) * t45, 0, 0, 0, 0, 0, 0, t4, t59, t18, -g(1) * t64 - g(2) * t54 - t52, 0, 0, 0, 0, 0, 0, t4, t18, -t59, -g(1) * (t12 * pkin(5) + t11 * qJ(6) + t64) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t54) - t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t13, 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t46, t51 * t43, -t13, g(3) * t57 + t92 * (pkin(3) * t47 + pkin(8) * t44) 0, 0, 0, 0, 0, 0, -t6, t5, -t13, -g(1) * (-t44 * t77 + t20) - g(2) * (-t34 * t76 + t27) - g(3) * (-t44 * t34 - t75) 0, 0, 0, 0, 0, 0, -t6, -t13, -t5, -g(1) * t20 - g(2) * t27 + (g(1) * t77 + g(3) * t53) * t44 + (g(3) * t49 + t53 * t88 - t56 * t89) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16 + t67, g(1) * t15 - g(2) * t17 + t46 * t86, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, pkin(4) * t67 - g(1) * t55 - g(2) * t69, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (t55 + t62) - g(2) * (t61 + t69) - g(3) * (t25 + (-pkin(5) * t35 - t90) * t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t62 - g(2) * t61 - g(3) * (-pkin(5) * t85 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
