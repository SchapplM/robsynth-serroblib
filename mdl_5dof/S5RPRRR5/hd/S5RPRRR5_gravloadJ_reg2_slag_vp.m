% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR5
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = qJ(1) + pkin(9);
t23 = qJ(3) + t28;
t18 = sin(t23);
t19 = cos(t23);
t32 = cos(qJ(4));
t20 = t32 * pkin(4) + pkin(3);
t34 = -pkin(8) - pkin(7);
t43 = t18 * t20 + t19 * t34;
t42 = t19 * pkin(3) + t18 * pkin(7);
t21 = sin(t28);
t31 = sin(qJ(1));
t41 = t31 * pkin(1) + pkin(2) * t21;
t22 = cos(t28);
t33 = cos(qJ(1));
t40 = t33 * pkin(1) + pkin(2) * t22;
t39 = t18 * pkin(3) - t19 * pkin(7);
t38 = -t18 * t34 + t19 * t20;
t37 = g(2) * t19 + g(3) * t18;
t7 = g(2) * t18 - g(3) * t19;
t36 = -g(2) * t33 - g(3) * t31;
t30 = sin(qJ(4));
t35 = -g(1) * t32 + t7 * t30;
t29 = qJ(4) + qJ(5);
t25 = cos(t29);
t24 = sin(t29);
t6 = t37 * t32;
t5 = t37 * t30;
t4 = t37 * t25;
t3 = t37 * t24;
t2 = -g(1) * t25 + t7 * t24;
t1 = g(1) * t24 + t7 * t25;
t8 = [0, 0, 0, 0, 0, 0, t36, g(2) * t31 - g(3) * t33, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t22 - g(3) * t21, g(2) * t21 - g(3) * t22, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, -t37, t7, 0, -g(2) * t40 - g(3) * t41, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * (t40 + t42) - g(3) * (t39 + t41), 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(2) * (t38 + t40) - g(3) * (t41 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * t42 - g(3) * t39, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(2) * t38 - g(3) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, g(1) * t30 + t7 * t32, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t35 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t8;
