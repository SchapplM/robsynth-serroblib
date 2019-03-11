% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = sin(qJ(5));
t46 = cos(qJ(5));
t80 = pkin(5) * t46 + qJ(6) * t43;
t45 = sin(qJ(1));
t48 = cos(qJ(1));
t27 = g(1) * t48 + g(2) * t45;
t42 = qJ(2) + qJ(3);
t36 = pkin(10) + t42;
t30 = sin(t36);
t54 = t27 * t30;
t31 = cos(t36);
t79 = -t31 * pkin(4) - t30 * pkin(9);
t49 = -pkin(8) - pkin(7);
t37 = sin(t42);
t78 = pkin(3) * t37;
t77 = pkin(4) * t30;
t73 = g(3) * t30;
t72 = g(3) * t43;
t71 = t30 * t48;
t70 = t31 * t48;
t69 = t45 * t43;
t68 = t45 * t46;
t41 = -qJ(4) + t49;
t67 = t48 * t41;
t66 = t48 * t43;
t65 = t48 * t46;
t44 = sin(qJ(2));
t18 = -t44 * pkin(2) - t78;
t20 = t45 * t31 * pkin(9);
t64 = t45 * t18 + t20;
t23 = pkin(9) * t70;
t63 = t48 * t18 + t23;
t38 = cos(t42);
t32 = pkin(3) * t38;
t47 = cos(qJ(2));
t39 = t47 * pkin(2);
t62 = t32 + t39;
t17 = pkin(1) + t62;
t14 = t48 * t17;
t60 = pkin(4) * t70 + pkin(9) * t71 + t14;
t59 = t32 - t79;
t58 = t80 * t31 + t59;
t57 = -t77 - t78;
t10 = t31 * t69 + t65;
t12 = t31 * t66 - t68;
t56 = g(1) * t10 - g(2) * t12;
t26 = g(1) * t45 - g(2) * t48;
t1 = g(1) * t12 + g(2) * t10 + t30 * t72;
t11 = t31 * t68 - t66;
t13 = t31 * t65 + t69;
t53 = g(1) * t13 + g(2) * t11 + t46 * t73;
t5 = -g(3) * t31 + t54;
t7 = -g(3) * t38 + t27 * t37;
t52 = -g(3) * t47 + t27 * t44;
t51 = (-g(1) * (-t17 + t79) + g(2) * t41) * t45;
t50 = (pkin(4) + t80) * t54;
t35 = t39 + pkin(1);
t9 = t26 * t30;
t8 = g(3) * t37 + t27 * t38;
t6 = t27 * t31 + t73;
t4 = t5 * t46;
t3 = -t31 * t72 + t43 * t54;
t2 = g(1) * t11 - g(2) * t13;
t15 = [0, 0, 0, 0, 0, 0, t26, t27, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t47, -t26 * t44, -t27, -g(1) * (-t45 * pkin(1) + t48 * pkin(7)) - g(2) * (t48 * pkin(1) + t45 * pkin(7)) 0, 0, 0, 0, 0, 0, t26 * t38, -t26 * t37, -t27, -g(1) * (-t45 * t35 - t48 * t49) - g(2) * (t48 * t35 - t45 * t49) 0, 0, 0, 0, 0, 0, t26 * t31, -t9, -t27, -g(1) * (-t45 * t17 - t67) - g(2) * (-t45 * t41 + t14) 0, 0, 0, 0, 0, 0, t2, -t56, t9, g(1) * t67 - g(2) * t60 + t51, 0, 0, 0, 0, 0, 0, t2, t9, t56, -g(1) * (-t11 * pkin(5) - t10 * qJ(6) - t67) - g(2) * (t13 * pkin(5) + t12 * qJ(6) + t60) + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, g(3) * t44 + t27 * t47, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t52 * pkin(2), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(3) * t62 - t27 * t18, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-pkin(4) * t71 + t63) - g(2) * (-t45 * t77 + t64) - g(3) * (t39 + t59) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t63 - g(2) * t64 - g(3) * (t39 + t58) + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t7 * pkin(3), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t57 * t48 + t23) - g(2) * (t57 * t45 + t20) - g(3) * t59, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (-t48 * t78 + t23) - g(2) * (-t45 * t78 + t20) - g(3) * t58 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t53, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t53, -g(1) * (-t12 * pkin(5) + t13 * qJ(6)) - g(2) * (-t10 * pkin(5) + t11 * qJ(6)) - (-pkin(5) * t43 + qJ(6) * t46) * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t15;
