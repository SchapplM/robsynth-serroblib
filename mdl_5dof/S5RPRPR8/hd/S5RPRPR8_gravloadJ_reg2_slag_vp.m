% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t22 = cos(pkin(9));
t12 = t22 * pkin(4) + pkin(3);
t23 = -pkin(7) - qJ(4);
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t30 = t26 * t12 - t24 * t23;
t20 = qJ(1) + pkin(8);
t14 = sin(t20);
t16 = cos(t20);
t8 = g(1) * t16 + g(2) * t14;
t5 = -g(3) * t26 + t8 * t24;
t48 = g(1) * t14;
t45 = g(3) * t24;
t43 = t14 * t26;
t21 = sin(pkin(9));
t42 = t16 * t21;
t41 = t16 * t26;
t40 = t21 * t26;
t39 = t22 * t26;
t27 = cos(qJ(1));
t36 = t27 * pkin(1) + t16 * pkin(2) + t14 * pkin(6);
t25 = sin(qJ(1));
t35 = -t25 * pkin(1) + t16 * pkin(6);
t34 = -g(2) * t16 + t48;
t33 = g(1) * t25 - g(2) * t27;
t32 = t26 * pkin(3) + t24 * qJ(4);
t19 = pkin(9) + qJ(5);
t15 = cos(t19);
t13 = sin(t19);
t7 = t34 * t24;
t6 = t8 * t26 + t45;
t4 = t14 * t13 + t15 * t41;
t3 = -t13 * t41 + t14 * t15;
t2 = t16 * t13 - t15 * t43;
t1 = t13 * t43 + t16 * t15;
t9 = [0, 0, 0, 0, 0, 0, t33, g(1) * t27 + g(2) * t25, 0, 0, 0, 0, 0, 0, 0, 0, t34, t8, 0, t33 * pkin(1), 0, 0, 0, 0, 0, 0, t34 * t26, -t7, -t8, -g(1) * (-t14 * pkin(2) + t35) - g(2) * t36, 0, 0, 0, 0, 0, 0, -g(1) * (-t14 * t39 + t42) - g(2) * (t14 * t21 + t16 * t39), -g(1) * (t14 * t40 + t16 * t22) - g(2) * (t14 * t22 - t16 * t40), t7, -g(1) * t35 - g(2) * (t32 * t16 + t36) - (-pkin(2) - t32) * t48, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t4, -g(1) * t1 - g(2) * t3, t7, -g(1) * (pkin(4) * t42 + t35) - g(2) * (t30 * t16 + t36) + (-g(1) * (-pkin(2) - t30) - g(2) * pkin(4) * t21) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t22, -t5 * t21, -t6, -g(3) * t32 + t8 * (pkin(3) * t24 - qJ(4) * t26), 0, 0, 0, 0, 0, 0, t5 * t15, -t5 * t13, -t6, -g(3) * t30 + t8 * (t12 * t24 + t23 * t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t1 + t13 * t45, g(1) * t4 - g(2) * t2 + t15 * t45, 0, 0;];
taug_reg = t9;
