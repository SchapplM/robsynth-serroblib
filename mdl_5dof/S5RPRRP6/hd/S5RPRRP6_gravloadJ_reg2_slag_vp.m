% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = cos(qJ(4));
t18 = t27 * pkin(4) + pkin(3);
t23 = -qJ(5) - pkin(7);
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t31 = t28 * t18 - t25 * t23;
t22 = qJ(1) + pkin(8);
t19 = sin(t22);
t20 = cos(t22);
t14 = g(1) * t20 + g(2) * t19;
t24 = sin(qJ(4));
t42 = t24 * t28;
t11 = t19 * t27 - t20 * t42;
t45 = g(3) * t25;
t9 = t19 * t42 + t20 * t27;
t1 = -g(1) * t11 + g(2) * t9 + t24 * t45;
t7 = -g(3) * t28 + t14 * t25;
t48 = g(1) * t19;
t43 = t20 * t24;
t40 = t27 * t28;
t29 = cos(qJ(1));
t37 = t29 * pkin(1) + t20 * pkin(2) + t19 * pkin(6);
t26 = sin(qJ(1));
t36 = -t26 * pkin(1) + t20 * pkin(6);
t35 = t28 * pkin(3) + t25 * pkin(7);
t33 = -g(2) * t20 + t48;
t32 = g(1) * t26 - g(2) * t29;
t13 = t33 * t25;
t12 = t19 * t24 + t20 * t40;
t10 = -t19 * t40 + t43;
t8 = t14 * t28 + t45;
t6 = t7 * t27;
t5 = t7 * t24;
t4 = -g(1) * t10 - g(2) * t12;
t3 = -g(1) * t9 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t10 + t27 * t45;
t15 = [0, 0, 0, 0, 0, 0, t32, g(1) * t29 + g(2) * t26, 0, 0, 0, 0, 0, 0, 0, 0, t33, t14, 0, t32 * pkin(1), 0, 0, 0, 0, 0, 0, t33 * t28, -t13, -t14, -g(1) * (-t19 * pkin(2) + t36) - g(2) * t37, 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * t36 - g(2) * (t35 * t20 + t37) - (-pkin(2) - t35) * t48, 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * (pkin(4) * t43 + t36) - g(2) * (t31 * t20 + t37) + (-g(1) * (-pkin(2) - t31) - g(2) * pkin(4) * t24) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t35 + t14 * (pkin(3) * t25 - pkin(7) * t28), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t31 + t14 * (t18 * t25 + t23 * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t15;
