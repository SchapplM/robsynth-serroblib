% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t52 = g(1) * t45 + g(2) * t42;
t39 = qJ(3) + qJ(4);
t34 = sin(t39);
t35 = cos(t39);
t59 = t45 * t35;
t44 = cos(qJ(2));
t63 = t42 * t44;
t15 = t34 * t63 + t59;
t60 = t45 * t34;
t17 = t42 * t35 - t44 * t60;
t41 = sin(qJ(2));
t68 = g(3) * t41;
t4 = -g(1) * t17 + g(2) * t15 + t34 * t68;
t48 = -g(3) * t44 + t52 * t41;
t40 = sin(qJ(3));
t26 = t40 * pkin(3) + pkin(4) * t34;
t66 = pkin(7) + t26;
t36 = qJ(5) + t39;
t32 = sin(t36);
t65 = t32 * t41;
t33 = cos(t36);
t64 = t33 * t41;
t62 = t45 * t32;
t61 = t45 * t33;
t58 = t45 * t40;
t43 = cos(qJ(3));
t57 = t45 * t43;
t27 = t43 * pkin(3) + pkin(4) * t35;
t11 = t32 * t63 + t61;
t12 = t33 * t63 - t62;
t55 = -t11 * pkin(5) + t12 * qJ(6);
t13 = -t42 * t33 + t44 * t62;
t14 = t42 * t32 + t44 * t61;
t54 = -t13 * pkin(5) + t14 * qJ(6);
t53 = g(1) * t11 - g(2) * t13;
t51 = g(1) * t42 - g(2) * t45;
t25 = pkin(2) + t27;
t38 = -pkin(10) - pkin(9) - pkin(8);
t50 = t44 * t25 - t41 * t38 + pkin(1);
t49 = pkin(5) * t33 + qJ(6) * t32 + t25;
t1 = g(1) * t13 + g(2) * t11 + g(3) * t65;
t3 = g(1) * t14 + g(2) * t12 + g(3) * t64;
t28 = qJ(6) * t64;
t46 = -g(1) * t54 - g(2) * t55 - g(3) * (-pkin(5) * t65 + t28);
t24 = t51 * t41;
t23 = t42 * t40 + t44 * t57;
t22 = t42 * t43 - t44 * t58;
t21 = -t43 * t63 + t58;
t20 = t40 * t63 + t57;
t19 = t52 * t44 + t68;
t18 = t42 * t34 + t44 * t59;
t16 = -t35 * t63 + t60;
t8 = t48 * t33;
t7 = t48 * t32;
t6 = g(1) * t12 - g(2) * t14;
t5 = g(1) * t18 - g(2) * t16 + t35 * t68;
t2 = [0, t51, t52, 0, 0, 0, 0, 0, t51 * t44, -t24, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t23, -g(1) * t20 - g(2) * t22, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, 0, 0, 0, 0, 0, t6, -t53, t6, t24, t53, -g(1) * (-t12 * pkin(5) - t11 * qJ(6)) - g(2) * (t14 * pkin(5) + t13 * qJ(6)) + (-g(1) * t66 - g(2) * t50) * t45 + (g(1) * t50 - g(2) * t66) * t42; 0, 0, 0, 0, 0, 0, 0, 0, t48, t19, 0, 0, 0, 0, 0, t48 * t43, -t48 * t40, 0, 0, 0, 0, 0, t48 * t35, -t48 * t34, 0, 0, 0, 0, 0, t8, -t7, t8, -t19, t7 (-g(3) * t49 + t52 * t38) * t44 + (g(3) * t38 + t52 * t49) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t22 + g(2) * t20 + t40 * t68, g(1) * t23 - g(2) * t21 + t43 * t68, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, -g(1) * (-t45 * t44 * t26 + t42 * t27 + t54) - g(2) * (-t26 * t63 - t45 * t27 + t55) - g(3) * (t28 + (-pkin(5) * t32 - t26) * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t4 * pkin(4) + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
