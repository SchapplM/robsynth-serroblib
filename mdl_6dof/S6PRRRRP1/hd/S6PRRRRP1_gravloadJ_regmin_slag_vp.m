% Calculate minimal parameter regressor of gravitation load for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x27]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = sin(pkin(11));
t38 = sin(qJ(2));
t41 = cos(qJ(2));
t52 = cos(pkin(11));
t53 = cos(pkin(6));
t48 = t53 * t52;
t19 = t33 * t38 - t41 * t48;
t63 = g(2) * t19;
t20 = t33 * t41 + t38 * t48;
t62 = g(2) * t20;
t34 = sin(pkin(6));
t61 = g(3) * t34;
t32 = qJ(3) + qJ(4);
t31 = cos(t32);
t36 = sin(qJ(5));
t60 = t31 * t36;
t39 = cos(qJ(5));
t59 = t31 * t39;
t58 = t33 * t34;
t57 = t34 * t38;
t40 = cos(qJ(3));
t56 = t34 * t40;
t55 = t36 * t41;
t54 = t39 * t41;
t51 = pkin(5) * t36 + pkin(8) + pkin(9);
t50 = t33 * t53;
t49 = t34 * t52;
t28 = t39 * pkin(5) + pkin(4);
t30 = sin(t32);
t35 = -qJ(6) - pkin(10);
t47 = t40 * pkin(3) + t28 * t31 - t30 * t35 + pkin(2);
t11 = t20 * t30 + t31 * t49;
t22 = -t38 * t50 + t52 * t41;
t13 = t22 * t30 - t31 * t58;
t17 = t30 * t57 - t53 * t31;
t3 = g(1) * t13 + g(2) * t11 + g(3) * t17;
t12 = t20 * t31 - t30 * t49;
t14 = t22 * t31 + t30 * t58;
t18 = t53 * t30 + t31 * t57;
t5 = g(1) * t14 + g(2) * t12 + g(3) * t18;
t21 = t52 * t38 + t41 * t50;
t46 = -g(1) * t21 + t41 * t61 - t63;
t45 = -g(1) * (-t13 * t28 - t14 * t35) - g(2) * (-t11 * t28 - t12 * t35) - g(3) * (-t17 * t28 - t18 * t35);
t44 = -g(1) * (-t14 * t36 + t21 * t39) - g(2) * (-t12 * t36 + t19 * t39) - g(3) * (-t18 * t36 - t34 * t54);
t37 = sin(qJ(3));
t43 = -g(1) * (-t22 * t37 + t33 * t56) - g(2) * (-t20 * t37 - t40 * t49) - g(3) * (-t37 * t57 + t53 * t40);
t6 = t46 * t30;
t2 = t3 * t39;
t1 = t3 * t36;
t4 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t46, g(1) * t22 + g(3) * t57 + t62, 0, 0, 0, 0, 0, -t46 * t40, t46 * t37, 0, 0, 0, 0, 0, -t46 * t31, t6, 0, 0, 0, 0, 0, -g(1) * (-t21 * t59 + t22 * t36) - g(2) * (-t19 * t59 + t20 * t36) - (t31 * t54 + t36 * t38) * t61, -g(1) * (t21 * t60 + t22 * t39) - g(2) * (t19 * t60 + t20 * t39) - (-t31 * t55 + t38 * t39) * t61, -t6, -g(1) * (-t47 * t21 + t51 * t22) - t51 * t62 + t47 * t63 - (t51 * t38 + t47 * t41) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -g(1) * (-t22 * t40 - t37 * t58) - g(2) * (-t20 * t40 + t37 * t49) - g(3) * (-t53 * t37 - t38 * t56) 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, t2, -t1, -t5, t43 * pkin(3) + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t5, 0, 0, 0, 0, 0, t2, -t1, -t5, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -g(1) * (-t14 * t39 - t21 * t36) - g(2) * (-t12 * t39 - t19 * t36) - g(3) * (-t18 * t39 + t34 * t55) 0, t44 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3;];
taug_reg  = t4;
