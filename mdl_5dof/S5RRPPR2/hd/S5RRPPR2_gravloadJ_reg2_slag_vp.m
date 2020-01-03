% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR2
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
% Datum: 2020-01-03 11:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = sin(pkin(9));
t32 = cos(pkin(9));
t47 = pkin(4) * t32 + pkin(7) * t31;
t44 = g(1) * t31;
t33 = sin(qJ(5));
t43 = t32 * t33;
t35 = cos(qJ(5));
t42 = t32 * t35;
t30 = qJ(1) + qJ(2);
t25 = pkin(8) + t30;
t21 = sin(t25);
t22 = cos(t25);
t27 = cos(t30);
t24 = pkin(2) * t27;
t41 = t22 * pkin(3) + t21 * qJ(4) + t24;
t26 = sin(t30);
t23 = pkin(2) * t26;
t40 = t21 * pkin(3) - t22 * qJ(4) + t23;
t39 = t47 * t22 + t41;
t10 = g(2) * t22 + g(3) * t21;
t12 = -g(2) * t27 - g(3) * t26;
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t38 = -g(2) * t36 - g(3) * t34;
t37 = t47 * t21 + t40;
t29 = t36 * pkin(1);
t28 = t34 * pkin(1);
t11 = g(2) * t26 - g(3) * t27;
t9 = g(2) * t21 - g(3) * t22;
t8 = t10 * t32;
t7 = t10 * t31;
t6 = t21 * t33 + t22 * t42;
t5 = -t21 * t35 + t22 * t43;
t4 = t21 * t42 - t22 * t33;
t3 = -t21 * t43 - t22 * t35;
t2 = -g(2) * t6 - g(3) * t4;
t1 = g(2) * t5 - g(3) * t3;
t13 = [0, 0, 0, 0, 0, 0, t38, g(2) * t34 - g(3) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, t38 * pkin(1), 0, 0, 0, 0, 0, 0, -t10, t9, 0, -g(2) * (t24 + t29) - g(3) * (t23 + t28), 0, 0, 0, 0, 0, 0, -t8, t7, -t9, -g(2) * (t29 + t41) - g(3) * (t28 + t40), 0, 0, 0, 0, 0, 0, t2, t1, -t7, -g(2) * (t29 + t39) - g(3) * (t28 + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, t12 * pkin(2), 0, 0, 0, 0, 0, 0, -t8, t7, -t9, -g(2) * t41 - g(3) * t40, 0, 0, 0, 0, 0, 0, t2, t1, -t7, -g(2) * t39 - g(3) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t3 - g(3) * t5 + t33 * t44, g(2) * t4 - g(3) * t6 + t35 * t44, 0, 0;];
taug_reg = t13;
