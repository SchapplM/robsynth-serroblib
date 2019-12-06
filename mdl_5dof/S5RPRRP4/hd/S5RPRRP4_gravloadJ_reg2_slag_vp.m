% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = cos(pkin(8));
t29 = sin(qJ(3));
t32 = cos(qJ(1));
t40 = t32 * t29;
t30 = sin(qJ(1));
t31 = cos(qJ(3));
t43 = t30 * t31;
t11 = t28 * t40 - t43;
t27 = sin(pkin(8));
t52 = g(1) * t27;
t39 = t32 * t31;
t44 = t30 * t29;
t9 = t28 * t44 + t39;
t56 = -g(2) * t9 + g(3) * t11 + t29 * t52;
t26 = qJ(3) + qJ(4);
t21 = sin(t26);
t22 = cos(t26);
t41 = t32 * t22;
t46 = t30 * t21;
t5 = t28 * t46 + t41;
t42 = t32 * t21;
t45 = t30 * t22;
t7 = t28 * t42 - t45;
t1 = -g(2) * t5 + g(3) * t7 + t21 * t52;
t33 = -pkin(7) - pkin(6);
t51 = g(2) * t32;
t23 = t32 * qJ(2);
t49 = g(3) * t23;
t48 = t29 * pkin(3);
t15 = pkin(4) * t21 + t48;
t47 = t15 * t28;
t24 = t31 * pkin(3);
t16 = pkin(4) * t22 + t24;
t18 = g(3) * t30 + t51;
t17 = g(2) * t30 - g(3) * t32;
t36 = pkin(2) * t28 + pkin(6) * t27 + pkin(1);
t35 = (pkin(2) + t16) * t28 - (-qJ(5) + t33) * t27 + pkin(1);
t34 = (t24 + pkin(2)) * t28 - t27 * t33 + pkin(1);
t13 = t18 * t27;
t12 = -t28 * t39 - t44;
t10 = t28 * t43 - t40;
t8 = -t28 * t41 - t46;
t6 = t28 * t45 - t42;
t4 = -g(2) * t8 + g(3) * t6;
t3 = -g(2) * t7 - g(3) * t5;
t2 = -g(2) * t6 - g(3) * t8 + t22 * t52;
t14 = [0, 0, 0, 0, 0, 0, t18, -t17, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t28, -t13, t17, -g(2) * (-t32 * pkin(1) - t30 * qJ(2)) - g(3) * (-t30 * pkin(1) + t23), 0, 0, 0, 0, 0, 0, -g(2) * t12 + g(3) * t10, -g(2) * t11 - g(3) * t9, t13, -t49 + t36 * t51 + (g(2) * qJ(2) + g(3) * t36) * t30, 0, 0, 0, 0, 0, 0, t4, t3, t13, -t49 + (g(2) * t34 - g(3) * t48) * t32 + (-g(2) * (-qJ(2) - t48) + g(3) * t34) * t30, 0, 0, 0, 0, 0, 0, t4, t3, t13, -t49 + (g(2) * t35 - g(3) * t15) * t32 + (-g(2) * (-qJ(2) - t15) + g(3) * t35) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -g(2) * t10 - g(3) * t12 + t31 * t52, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t56 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, t15 * t52 - g(2) * (t32 * t16 + t30 * t47) - g(3) * (t30 * t16 - t32 * t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t28 + t17 * t27;];
taug_reg = t14;
