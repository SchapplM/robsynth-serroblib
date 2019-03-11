% Calculate inertial parameters regressor of gravitation load for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPPRPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = sin(pkin(9));
t46 = cos(pkin(9));
t50 = sin(qJ(1));
t51 = cos(qJ(1));
t11 = -t50 * t45 - t51 * t46;
t12 = t51 * t45 - t50 * t46;
t39 = g(1) * t11 + g(2) * t12;
t26 = qJ(4) + pkin(10);
t19 = sin(t26);
t20 = cos(t26);
t57 = -g(3) * t20 + t39 * t19;
t54 = g(3) * t19;
t31 = cos(qJ(4));
t52 = t31 * pkin(4);
t28 = sin(qJ(6));
t49 = t20 * t28;
t30 = cos(qJ(6));
t48 = t20 * t30;
t47 = t51 * pkin(1) + t50 * qJ(2);
t44 = t51 * pkin(2) + t47;
t18 = pkin(3) + t52;
t27 = -qJ(5) - pkin(7);
t43 = -t11 * t18 - t12 * t27 + t44;
t42 = -t50 * pkin(1) + t51 * qJ(2);
t41 = t20 * pkin(5) + t19 * pkin(8);
t40 = g(1) * t12 - g(2) * t11;
t38 = t11 * t28 + t12 * t48;
t37 = -t11 * t30 + t12 * t49;
t35 = -t50 * pkin(2) + t42;
t29 = sin(qJ(4));
t33 = g(3) * t31 - t39 * t29;
t32 = -t11 * t27 + t12 * t18 + t35;
t14 = g(1) * t51 + g(2) * t50;
t13 = g(1) * t50 - g(2) * t51;
t4 = -t11 * t48 + t12 * t28;
t3 = t11 * t49 + t12 * t30;
t2 = t40 * t19;
t1 = -t39 * t20 - t54;
t5 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, -t14, -g(1) * t42 - g(2) * t47, 0, 0, 0, 0, 0, 0, -t40, t39, 0, -g(1) * t35 - g(2) * t44, 0, 0, 0, 0, 0, 0, -t40 * t31, t40 * t29, -t39, -g(1) * (t12 * pkin(3) + t11 * pkin(7) + t35) - g(2) * (-t11 * pkin(3) + t12 * pkin(7) + t44) 0, 0, 0, 0, 0, 0, -t40 * t20, t2, -t39, -g(1) * t32 - g(2) * t43, 0, 0, 0, 0, 0, 0, -g(1) * t38 - g(2) * t4, g(1) * t37 - g(2) * t3, -t2, -g(1) * (t41 * t12 + t32) - g(2) * (-t41 * t11 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -g(3) * t29 - t39 * t31, 0, 0, 0, 0, 0, 0, 0, 0, -t57, t1, 0, t33 * pkin(4), 0, 0, 0, 0, 0, 0, -t57 * t30, t57 * t28, -t1, -g(3) * (-t41 - t52) - t39 * (pkin(4) * t29 + pkin(5) * t19 - pkin(8) * t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t37 - t28 * t54, g(1) * t4 - g(2) * t38 - t30 * t54, 0, 0;];
taug_reg  = t5;
