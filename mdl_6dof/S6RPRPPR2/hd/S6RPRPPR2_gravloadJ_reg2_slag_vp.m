% Calculate inertial parameters regressor of gravitation load for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = qJ(1) + pkin(9);
t22 = sin(t28);
t24 = cos(t28);
t10 = g(1) * t24 + g(2) * t22;
t27 = qJ(3) + pkin(10);
t21 = sin(t27);
t60 = t10 * t21;
t16 = t21 * qJ(5);
t23 = cos(t27);
t59 = -t23 * pkin(4) - t16;
t2 = g(3) * t21 + t10 * t23;
t31 = sin(qJ(3));
t57 = pkin(3) * t31;
t53 = g(3) * t23;
t52 = t23 * pkin(8);
t32 = sin(qJ(1));
t51 = t32 * pkin(1);
t30 = sin(qJ(6));
t50 = t22 * t30;
t33 = cos(qJ(6));
t49 = t22 * t33;
t48 = t23 * t24;
t47 = t24 * t30;
t46 = t24 * t33;
t34 = cos(qJ(3));
t25 = t34 * pkin(3);
t20 = t25 + pkin(2);
t35 = cos(qJ(1));
t26 = t35 * pkin(1);
t45 = t24 * t20 + t26;
t44 = qJ(5) * t23;
t43 = t25 - t59;
t42 = pkin(4) * t48 + t24 * t16 + t45;
t41 = -pkin(4) * t21 - t57;
t9 = g(1) * t22 - g(2) * t24;
t40 = g(1) * t32 - g(2) * t35;
t29 = -qJ(4) - pkin(7);
t39 = -t24 * t29 - t51;
t38 = -t20 + t59;
t36 = -g(3) * t34 + t10 * t31;
t13 = t24 * t44;
t11 = t22 * t44;
t8 = -t21 * t50 + t46;
t7 = t21 * t49 + t47;
t6 = t21 * t47 + t49;
t5 = t21 * t46 - t50;
t4 = t9 * t23;
t3 = t9 * t21;
t1 = -t53 + t60;
t12 = [0, 0, 0, 0, 0, 0, t40, g(1) * t35 + g(2) * t32, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t40 * pkin(1), 0, 0, 0, 0, 0, 0, t9 * t34, -t9 * t31, -t10, -g(1) * (-t22 * pkin(2) + t24 * pkin(7) - t51) - g(2) * (t24 * pkin(2) + t22 * pkin(7) + t26) 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(1) * (-t22 * t20 + t39) - g(2) * (-t22 * t29 + t45) 0, 0, 0, 0, 0, 0, -t10, -t4, t3, -g(1) * t39 - g(2) * t42 + (-g(1) * t38 + g(2) * t29) * t22, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5, t4, -g(1) * (t24 * pkin(5) + t39) - g(2) * (pkin(8) * t48 + t42) + (-g(1) * (t38 - t52) - g(2) * (pkin(5) - t29)) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t31 + t10 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t36 * pkin(3), 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (t41 * t24 + t13) - g(2) * (t41 * t22 + t11) - g(3) * t43, 0, 0, 0, 0, 0, 0, -t2 * t30, -t2 * t33, t1, -g(1) * (-t24 * t57 + t13) - g(2) * (-t22 * t57 + t11) - g(3) * (t43 + t52) + (pkin(4) + pkin(8)) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t33 * t53, g(1) * t6 - g(2) * t8 - t30 * t53, 0, 0;];
taug_reg  = t12;
