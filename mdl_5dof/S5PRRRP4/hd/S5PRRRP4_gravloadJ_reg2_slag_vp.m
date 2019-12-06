% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t25 = sin(qJ(4));
t27 = cos(qJ(4));
t51 = pkin(4) * t27 + qJ(5) * t25;
t22 = qJ(2) + qJ(3);
t19 = sin(t22);
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t34 = g(1) * t24 + g(2) * t23;
t32 = t34 * t19;
t26 = sin(qJ(2));
t50 = pkin(2) * t26;
t49 = pkin(3) * t19;
t20 = cos(t22);
t47 = pkin(7) * t20;
t44 = g(3) * t19;
t43 = g(3) * t25;
t42 = t23 * t25;
t41 = t23 * t27;
t40 = t24 * t25;
t39 = t24 * t27;
t38 = t20 * pkin(3) + t19 * pkin(7);
t36 = t51 * t20 + t38;
t35 = -t49 - t50;
t6 = t20 * t42 + t39;
t8 = t20 * t40 - t41;
t1 = g(1) * t8 + g(2) * t6 + t19 * t43;
t7 = t20 * t41 - t40;
t9 = t20 * t39 + t42;
t31 = g(1) * t9 + g(2) * t7 + t27 * t44;
t4 = -g(3) * t20 + t32;
t28 = cos(qJ(2));
t30 = -g(3) * t28 + t34 * t26;
t29 = (pkin(3) + t51) * t32;
t21 = t28 * pkin(2);
t12 = t24 * t47;
t10 = t23 * t47;
t5 = t34 * t20 + t44;
t3 = t4 * t27;
t2 = -t20 * t43 + t25 * t32;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(3) * t26 + t34 * t28, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, t30 * pkin(2), 0, 0, 0, 0, 0, 0, t3, -t2, -t5, -g(1) * (t35 * t24 + t12) - g(2) * (t35 * t23 + t10) - g(3) * (t21 + t38), 0, 0, 0, 0, 0, 0, t3, -t5, t2, -g(1) * (-t24 * t50 + t12) - g(2) * (-t23 * t50 + t10) - g(3) * (t21 + t36) + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, -t5, -g(1) * (-t24 * t49 + t12) - g(2) * (-t23 * t49 + t10) - g(3) * t38, 0, 0, 0, 0, 0, 0, t3, -t5, t2, -g(1) * t12 - g(2) * t10 - g(3) * t36 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t31, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t31, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (-t6 * pkin(4) + t7 * qJ(5)) - (-pkin(4) * t25 + qJ(5) * t27) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
