% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = qJ(1) + qJ(2);
t23 = pkin(9) + t31;
t18 = sin(t23);
t19 = cos(t23);
t25 = sin(t31);
t20 = pkin(2) * t25;
t34 = cos(qJ(4));
t22 = pkin(4) * t34 + pkin(3);
t36 = -pkin(8) - pkin(7);
t43 = t18 * t22 + t19 * t36 + t20;
t27 = cos(t31);
t21 = pkin(2) * t27;
t42 = t19 * pkin(3) + t18 * pkin(7) + t21;
t41 = t18 * pkin(3) - pkin(7) * t19 + t20;
t40 = -t18 * t36 + t19 * t22 + t21;
t39 = g(2) * t19 + g(3) * t18;
t7 = g(2) * t18 - g(3) * t19;
t10 = -g(2) * t27 - g(3) * t25;
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t38 = -g(2) * t35 - g(3) * t33;
t32 = sin(qJ(4));
t37 = -g(1) * t34 + t7 * t32;
t30 = qJ(4) + qJ(5);
t29 = t35 * pkin(1);
t28 = t33 * pkin(1);
t26 = cos(t30);
t24 = sin(t30);
t9 = g(2) * t25 - g(3) * t27;
t6 = t39 * t34;
t5 = t39 * t32;
t4 = t39 * t26;
t3 = t39 * t24;
t2 = -g(1) * t26 + t7 * t24;
t1 = g(1) * t24 + t7 * t26;
t8 = [0, 0, 0, 0, 0, 0, t38, g(2) * t33 - g(3) * t35, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, t38 * pkin(1), 0, 0, 0, 0, 0, 0, -t39, t7, 0, -g(2) * (t21 + t29) - g(3) * (t20 + t28), 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * (t29 + t42) - g(3) * (t28 + t41), 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(2) * (t29 + t40) - g(3) * (t28 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t7, 0, t10 * pkin(2), 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * t42 - g(3) * t41, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(2) * t40 - g(3) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(1) * t32 + t7 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t37 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t8;
