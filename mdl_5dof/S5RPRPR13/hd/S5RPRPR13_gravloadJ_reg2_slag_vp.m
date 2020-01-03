% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = pkin(8) + qJ(3);
t20 = sin(t23);
t21 = cos(t23);
t35 = t21 * pkin(3) + t20 * qJ(4);
t28 = sin(qJ(1));
t30 = cos(qJ(1));
t14 = g(1) * t30 + g(2) * t28;
t49 = t14 * t20;
t2 = g(3) * t20 + t14 * t21;
t47 = pkin(3) * t20;
t43 = g(3) * t21;
t42 = t21 * pkin(7);
t26 = -pkin(6) - qJ(2);
t41 = pkin(4) - t26;
t27 = sin(qJ(5));
t40 = t28 * t27;
t29 = cos(qJ(5));
t39 = t28 * t29;
t38 = t30 * t26;
t37 = t30 * t27;
t36 = t30 * t29;
t34 = qJ(4) * t21;
t25 = cos(pkin(8));
t19 = t25 * pkin(2) + pkin(1);
t12 = t30 * t19;
t33 = g(2) * (t35 * t30 + t12);
t13 = g(1) * t28 - g(2) * t30;
t32 = -t19 - t35;
t11 = t30 * t34;
t9 = t28 * t34;
t8 = -t20 * t40 + t36;
t7 = t20 * t39 + t37;
t6 = t20 * t37 + t39;
t5 = t20 * t36 - t40;
t4 = t13 * t21;
t3 = t13 * t20;
t1 = -t43 + t49;
t10 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t25, -t13 * sin(pkin(8)), -t14, -g(1) * (-t28 * pkin(1) + t30 * qJ(2)) - g(2) * (t30 * pkin(1) + t28 * qJ(2)), 0, 0, 0, 0, 0, 0, t4, -t3, -t14, -g(1) * (-t28 * t19 - t38) - g(2) * (-t28 * t26 + t12), 0, 0, 0, 0, 0, 0, -t14, -t4, t3, g(1) * t38 - t33 + (-g(1) * t32 + g(2) * t26) * t28, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, g(1) * t7 - g(2) * t5, t4, -t33 + (-g(1) * t41 - g(2) * t42) * t30 + (-g(1) * (t32 - t42) - g(2) * t41) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -g(1) * (-t30 * t47 + t11) - g(2) * (-t28 * t47 + t9) - g(3) * t35, 0, 0, 0, 0, 0, 0, -t2 * t27, -t2 * t29, t1, -g(1) * t11 - g(2) * t9 - g(3) * (t35 + t42) + (pkin(3) + pkin(7)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7 + t29 * t43, g(1) * t6 - g(2) * t8 - t27 * t43, 0, 0;];
taug_reg = t10;
