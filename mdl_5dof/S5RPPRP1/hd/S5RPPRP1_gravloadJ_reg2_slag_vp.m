% Calculate inertial parameters regressor of gravitation load for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t21 = sin(pkin(8));
t22 = cos(pkin(8));
t26 = cos(qJ(4));
t42 = -t21 * (-qJ(5) - pkin(6)) + (t26 * pkin(4) + pkin(3)) * t22;
t24 = sin(qJ(4));
t39 = g(1) * t21;
t20 = qJ(1) + pkin(7);
t16 = sin(t20);
t17 = cos(t20);
t35 = t22 * t24;
t5 = -t16 * t35 - t17 * t26;
t7 = -t16 * t26 + t17 * t35;
t1 = -g(2) * t5 - g(3) * t7 + t24 * t39;
t37 = t16 * t24;
t34 = t22 * t26;
t25 = sin(qJ(1));
t33 = t25 * pkin(1) + t16 * pkin(2);
t27 = cos(qJ(1));
t31 = t27 * pkin(1) + t17 * pkin(2) + t16 * qJ(3);
t30 = -t17 * qJ(3) + t33;
t29 = pkin(3) * t22 + pkin(6) * t21;
t11 = g(2) * t17 + g(3) * t16;
t10 = g(2) * t16 - g(3) * t17;
t28 = -g(2) * t27 - g(3) * t25;
t9 = t11 * t21;
t8 = t17 * t34 + t37;
t6 = t16 * t34 - t17 * t24;
t4 = -g(2) * t8 - g(3) * t6;
t3 = g(2) * t7 - g(3) * t5;
t2 = g(2) * t6 - g(3) * t8 + t26 * t39;
t12 = [0, 0, 0, 0, 0, 0, t28, g(2) * t25 - g(3) * t27, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, t28 * pkin(1), 0, 0, 0, 0, 0, 0, -t11 * t22, t9, -t10, -g(2) * t31 - g(3) * t30, 0, 0, 0, 0, 0, 0, t4, t3, -t9, -g(2) * (t29 * t17 + t31) - g(3) * (t29 * t16 + t30), 0, 0, 0, 0, 0, 0, t4, t3, -t9, -g(2) * (pkin(4) * t37 + t31) - g(3) * (t42 * t16 + t33) + (-g(2) * t42 - g(3) * (-pkin(4) * t24 - qJ(3))) * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t22 - t10 * t21;];
taug_reg = t12;
