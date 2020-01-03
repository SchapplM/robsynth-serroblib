% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR3
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t27 = qJ(1) + qJ(2);
t22 = pkin(9) + t27;
t21 = qJ(4) + t22;
t15 = sin(t21);
t16 = cos(t21);
t41 = t16 * pkin(4) + t15 * pkin(8);
t17 = sin(t22);
t13 = pkin(3) * t17;
t23 = sin(t27);
t19 = pkin(2) * t23;
t40 = t13 + t19;
t18 = cos(t22);
t14 = pkin(3) * t18;
t24 = cos(t27);
t20 = pkin(2) * t24;
t39 = t14 + t20;
t29 = sin(qJ(1));
t25 = t29 * pkin(1);
t38 = t19 + t25;
t31 = cos(qJ(1));
t26 = t31 * pkin(1);
t37 = t20 + t26;
t36 = t15 * pkin(4) - t16 * pkin(8);
t35 = t39 + t41;
t34 = g(2) * t16 + g(3) * t15;
t3 = g(2) * t15 - g(3) * t16;
t8 = -g(2) * t24 - g(3) * t23;
t33 = -g(2) * t31 - g(3) * t29;
t32 = t36 + t40;
t30 = cos(qJ(5));
t28 = sin(qJ(5));
t7 = g(2) * t23 - g(3) * t24;
t6 = -g(2) * t18 - g(3) * t17;
t5 = g(2) * t17 - g(3) * t18;
t2 = t34 * t30;
t1 = t34 * t28;
t4 = [0, 0, 0, 0, 0, 0, t33, g(2) * t29 - g(3) * t31, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t33 * pkin(1), 0, 0, 0, 0, 0, 0, t6, t5, 0, -g(2) * t37 - g(3) * t38, 0, 0, 0, 0, 0, 0, -t34, t3, 0, -g(2) * (t14 + t37) - g(3) * (t13 + t38), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * (t26 + t35) - g(3) * (t25 + t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t8 * pkin(2), 0, 0, 0, 0, 0, 0, -t34, t3, 0, -g(2) * t39 - g(3) * t40, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * t35 - g(3) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * t41 - g(3) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t30 + t3 * t28, g(1) * t28 + t3 * t30, 0, 0;];
taug_reg = t4;
