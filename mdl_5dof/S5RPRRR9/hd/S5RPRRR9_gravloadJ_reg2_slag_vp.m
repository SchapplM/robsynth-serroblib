% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = pkin(9) + qJ(3);
t24 = qJ(4) + t28;
t19 = sin(t24);
t20 = cos(t24);
t48 = t20 * pkin(4) + t19 * pkin(8);
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t15 = g(1) * t35 + g(2) * t33;
t3 = -g(3) * t20 + t15 * t19;
t47 = pkin(4) * t19;
t46 = pkin(8) * t20;
t45 = g(3) * t19;
t30 = cos(pkin(9));
t21 = t30 * pkin(2) + pkin(1);
t32 = sin(qJ(5));
t43 = t33 * t32;
t34 = cos(qJ(5));
t42 = t33 * t34;
t41 = t35 * t32;
t40 = t35 * t34;
t31 = -pkin(6) - qJ(2);
t22 = sin(t28);
t38 = -pkin(3) * t22 - t47;
t14 = g(1) * t33 - g(2) * t35;
t23 = cos(t28);
t36 = -g(3) * t23 + t15 * t22;
t27 = -pkin(7) + t31;
t18 = pkin(3) * t23;
t13 = t35 * t46;
t12 = t33 * t46;
t11 = t18 + t21;
t10 = t35 * t11;
t9 = t20 * t40 + t43;
t8 = -t20 * t41 + t42;
t7 = -t20 * t42 + t41;
t6 = t20 * t43 + t40;
t5 = t14 * t19;
t4 = t15 * t20 + t45;
t2 = t3 * t34;
t1 = t3 * t32;
t16 = [0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t30, -t14 * sin(pkin(9)), -t15, -g(1) * (-t33 * pkin(1) + t35 * qJ(2)) - g(2) * (t35 * pkin(1) + t33 * qJ(2)), 0, 0, 0, 0, 0, 0, t14 * t23, -t14 * t22, -t15, -g(1) * (-t33 * t21 - t35 * t31) - g(2) * (t35 * t21 - t33 * t31), 0, 0, 0, 0, 0, 0, t14 * t20, -t5, -t15, -g(1) * (-t33 * t11 - t35 * t27) - g(2) * (-t33 * t27 + t10), 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, -g(1) * t6 - g(2) * t8, t5, -g(2) * t10 + (g(1) * t27 - g(2) * t48) * t35 + (-g(1) * (-t11 - t48) + g(2) * t27) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, g(3) * t22 + t15 * t23, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t36 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t38 * t35 + t13) - g(2) * (t38 * t33 + t12) - g(3) * (t18 + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t35 * t47 + t13) - g(2) * (-t33 * t47 + t12) - g(3) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t6 + t32 * t45, g(1) * t9 - g(2) * t7 + t34 * t45, 0, 0;];
taug_reg = t16;
