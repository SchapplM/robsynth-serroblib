% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t29 = qJ(1) + qJ(2);
t23 = pkin(8) + t29;
t16 = sin(t23);
t17 = cos(t23);
t31 = cos(pkin(9));
t18 = t31 * pkin(4) + pkin(3);
t24 = sin(t29);
t19 = pkin(2) * t24;
t32 = -pkin(7) - qJ(4);
t39 = t16 * t18 + t17 * t32 + t19;
t25 = cos(t29);
t20 = pkin(2) * t25;
t38 = t17 * pkin(3) + t16 * qJ(4) + t20;
t37 = -t16 * t32 + t17 * t18 + t20;
t36 = t16 * pkin(3) - t17 * qJ(4) + t19;
t6 = g(2) * t17 + g(3) * t16;
t5 = g(2) * t16 - g(3) * t17;
t8 = -g(2) * t25 - g(3) * t24;
t33 = sin(qJ(1));
t34 = cos(qJ(1));
t35 = -g(2) * t34 - g(3) * t33;
t28 = pkin(9) + qJ(5);
t27 = t34 * pkin(1);
t26 = t33 * pkin(1);
t22 = cos(t28);
t21 = sin(t28);
t7 = g(2) * t24 - g(3) * t25;
t4 = t6 * t31;
t3 = t6 * sin(pkin(9));
t2 = t6 * t22;
t1 = t6 * t21;
t9 = [0, 0, 0, 0, 0, 0, t35, g(2) * t33 - g(3) * t34, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t35 * pkin(1), 0, 0, 0, 0, 0, 0, -t6, t5, 0, -g(2) * (t20 + t27) - g(3) * (t19 + t26), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * (t27 + t38) - g(3) * (t26 + t36), 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(2) * (t27 + t37) - g(3) * (t26 + t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, t8 * pkin(2), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(2) * t38 - g(3) * t36, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(2) * t37 - g(3) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t22 + t5 * t21, g(1) * t21 + t5 * t22, 0, 0;];
taug_reg = t9;
