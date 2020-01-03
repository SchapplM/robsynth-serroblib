% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR7
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t28 = cos(qJ(4));
t17 = t28 * pkin(4) + pkin(3);
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t31 = -pkin(8) - pkin(7);
t35 = t29 * t17 - t26 * t31;
t23 = qJ(1) + pkin(9);
t18 = sin(t23);
t19 = cos(t23);
t13 = g(1) * t19 + g(2) * t18;
t25 = sin(qJ(4));
t46 = t25 * t29;
t10 = t18 * t28 - t19 * t46;
t51 = g(3) * t26;
t8 = t18 * t46 + t19 * t28;
t57 = -g(1) * t10 + g(2) * t8 + t25 * t51;
t33 = -g(3) * t29 + t13 * t26;
t54 = g(1) * t18;
t49 = t19 * t25;
t24 = qJ(4) + qJ(5);
t20 = sin(t24);
t48 = t20 * t29;
t21 = cos(t24);
t47 = t21 * t29;
t44 = t28 * t29;
t30 = cos(qJ(1));
t41 = t30 * pkin(1) + t19 * pkin(2) + t18 * pkin(6);
t27 = sin(qJ(1));
t40 = -t27 * pkin(1) + t19 * pkin(6);
t39 = t29 * pkin(3) + t26 * pkin(7);
t37 = -g(2) * t19 + t54;
t36 = g(1) * t27 - g(2) * t30;
t12 = t37 * t26;
t11 = t18 * t25 + t19 * t44;
t9 = -t18 * t44 + t49;
t7 = t13 * t29 + t51;
t6 = t18 * t20 + t19 * t47;
t5 = t18 * t21 - t19 * t48;
t4 = -t18 * t47 + t19 * t20;
t3 = t18 * t48 + t19 * t21;
t2 = g(1) * t6 - g(2) * t4 + t21 * t51;
t1 = -g(1) * t5 + g(2) * t3 + t20 * t51;
t14 = [0, 0, 0, 0, 0, 0, t36, g(1) * t30 + g(2) * t27, 0, 0, 0, 0, 0, 0, 0, 0, t37, t13, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, t37 * t29, -t12, -t13, -g(1) * (-t18 * pkin(2) + t40) - g(2) * t41, 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t12, -g(1) * t40 - g(2) * (t39 * t19 + t41) - (-pkin(2) - t39) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t12, -g(1) * (pkin(4) * t49 + t40) - g(2) * (t35 * t19 + t41) + (-g(1) * (-pkin(2) - t35) - g(2) * pkin(4) * t25) * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t7, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t28, -t33 * t25, -t7, -g(3) * t39 + t13 * (pkin(3) * t26 - pkin(7) * t29), 0, 0, 0, 0, 0, 0, t33 * t21, -t33 * t20, -t7, -g(3) * t35 + t13 * (t17 * t26 + t29 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, g(1) * t11 - g(2) * t9 + t28 * t51, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t57 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t14;
