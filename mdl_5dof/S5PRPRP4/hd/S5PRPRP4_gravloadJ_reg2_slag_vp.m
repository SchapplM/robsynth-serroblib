% Calculate inertial parameters regressor of gravitation load for
% S5PRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t19 = qJ(2) + pkin(8);
t16 = sin(t19);
t20 = sin(pkin(7));
t21 = cos(pkin(7));
t31 = g(1) * t21 + g(2) * t20;
t43 = t31 * t16;
t23 = sin(qJ(2));
t42 = pkin(2) * t23;
t17 = cos(t19);
t41 = pkin(6) * t17;
t38 = g(3) * t16;
t22 = sin(qJ(4));
t37 = t20 * t22;
t24 = cos(qJ(4));
t36 = t20 * t24;
t35 = t21 * t22;
t34 = t21 * t24;
t25 = cos(qJ(2));
t33 = t25 * pkin(2) + t17 * pkin(3) + t16 * pkin(6);
t32 = -pkin(3) * t16 - t42;
t30 = pkin(4) * t24 + qJ(5) * t22;
t5 = t17 * t37 + t34;
t7 = t17 * t35 - t36;
t1 = g(1) * t7 + g(2) * t5 + t22 * t38;
t6 = t17 * t36 - t35;
t8 = t17 * t34 + t37;
t28 = g(1) * t8 + g(2) * t6 + t24 * t38;
t27 = -g(3) * t17 + t43;
t26 = -g(3) * t25 + t31 * t23;
t11 = t21 * t41;
t10 = t20 * t41;
t9 = -g(1) * t20 + g(2) * t21;
t4 = t31 * t17 + t38;
t3 = t27 * t24;
t2 = t27 * t22;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, g(3) * t23 + t31 * t25, 0, 0, 0, 0, 0, 0, 0, 0, t27, t4, 0, t26 * pkin(2), 0, 0, 0, 0, 0, 0, t3, -t2, -t4, -g(1) * (t32 * t21 + t11) - g(2) * (t32 * t20 + t10) - g(3) * t33, 0, 0, 0, 0, 0, 0, t3, -t4, t2, -g(1) * (-t21 * t42 + t11) - g(2) * (-t20 * t42 + t10) - g(3) * (t30 * t17 + t33) + (pkin(3) + t30) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t28, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t28, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - (-pkin(4) * t22 + qJ(5) * t24) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
