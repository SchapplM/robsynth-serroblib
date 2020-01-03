% Calculate inertial parameters regressor of gravitation load for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t24 = cos(qJ(3));
t14 = t24 * pkin(3) + pkin(2);
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t27 = -pkin(7) - pkin(6);
t31 = t25 * t14 - t22 * t27;
t23 = sin(qJ(1));
t26 = cos(qJ(1));
t13 = g(1) * t26 + g(2) * t23;
t21 = sin(qJ(3));
t38 = t26 * t21;
t10 = t23 * t24 - t25 * t38;
t45 = g(3) * t22;
t37 = t26 * t24;
t42 = t23 * t25;
t8 = t21 * t42 + t37;
t51 = -g(1) * t10 + g(2) * t8 + t21 * t45;
t29 = -g(3) * t25 + t13 * t22;
t48 = g(1) * t23;
t20 = qJ(3) + qJ(4);
t15 = sin(t20);
t40 = t26 * t15;
t16 = cos(t20);
t39 = t26 * t16;
t36 = t26 * pkin(1) + t23 * pkin(5);
t34 = t25 * pkin(2) + t22 * pkin(6);
t32 = -g(2) * t26 + t48;
t18 = t26 * pkin(5);
t12 = t32 * t22;
t11 = t23 * t21 + t25 * t37;
t9 = -t24 * t42 + t38;
t7 = t13 * t25 + t45;
t6 = t23 * t15 + t25 * t39;
t5 = t23 * t16 - t25 * t40;
t4 = -t16 * t42 + t40;
t3 = t15 * t42 + t39;
t2 = g(1) * t6 - g(2) * t4 + t16 * t45;
t1 = -g(1) * t5 + g(2) * t3 + t15 * t45;
t17 = [0, 0, 0, 0, 0, 0, t32, t13, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t25, -t12, -t13, -g(1) * (-t23 * pkin(1) + t18) - g(2) * t36, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t12, -g(1) * t18 - g(2) * (t34 * t26 + t36) - (-pkin(1) - t34) * t48, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t12, -g(1) * (pkin(3) * t38 + t18) - g(2) * (t31 * t26 + t36) + (-g(1) * (-pkin(1) - t31) - g(2) * pkin(3) * t21) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t7, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t24, -t29 * t21, -t7, -g(3) * t34 + t13 * (pkin(2) * t22 - pkin(6) * t25), 0, 0, 0, 0, 0, 0, t29 * t16, -t29 * t15, -t7, -g(3) * t31 + t13 * (t14 * t22 + t25 * t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, g(1) * t11 - g(2) * t9 + t24 * t45, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t51 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t17;
