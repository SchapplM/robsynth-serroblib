% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = cos(qJ(1));
t27 = sin(qJ(1));
t57 = g(1) * t27;
t61 = g(2) * t30 - t57;
t26 = sin(qJ(3));
t25 = sin(qJ(4));
t43 = t30 * t25;
t28 = cos(qJ(4));
t47 = t27 * t28;
t10 = t26 * t43 + t47;
t29 = cos(qJ(3));
t53 = g(3) * t29;
t42 = t30 * t28;
t48 = t27 * t25;
t8 = -t26 * t48 + t42;
t60 = -g(1) * t8 - g(2) * t10 + t25 * t53;
t33 = -g(3) * t26 - t29 * t61;
t58 = -pkin(1) - pkin(6);
t52 = t29 * pkin(7);
t51 = t26 * t30;
t24 = qJ(4) + qJ(5);
t17 = sin(t24);
t50 = t27 * t17;
t18 = cos(t24);
t49 = t27 * t18;
t31 = -pkin(8) - pkin(7);
t46 = t29 * t31;
t45 = t30 * t17;
t44 = t30 * t18;
t41 = t30 * pkin(1) + t27 * qJ(2);
t39 = t30 * pkin(6) + t41;
t38 = g(2) * t39;
t36 = t26 * pkin(3) - t52;
t14 = g(1) * t30 + g(2) * t27;
t16 = t28 * pkin(4) + pkin(3);
t34 = t26 * t16 + t46;
t20 = t30 * qJ(2);
t12 = t14 * t29;
t11 = t26 * t42 - t48;
t9 = t26 * t47 + t43;
t7 = -g(2) * t51 + t26 * t57 + t53;
t6 = t26 * t44 - t50;
t5 = t26 * t45 + t49;
t4 = t26 * t49 + t45;
t3 = -t26 * t50 + t44;
t2 = g(1) * t4 - g(2) * t6 + t18 * t53;
t1 = -g(1) * t3 - g(2) * t5 + t17 * t53;
t13 = [0, 0, 0, 0, 0, 0, -t61, t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t14, -g(1) * (-t27 * pkin(1) + t20) - g(2) * t41, 0, 0, 0, 0, 0, 0, -t14 * t26, -t12, -t61, -g(1) * (t58 * t27 + t20) - t38, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t9, g(1) * t10 - g(2) * t8, t12, -g(1) * (pkin(3) * t51 - t30 * t52 + t20) - t38 + (-g(1) * t58 - g(2) * t36) * t27, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t12, -g(1) * (t16 * t51 + t30 * t46 + t20) - g(2) * (pkin(4) * t43 + t39) + (-g(1) * (-pkin(4) * t25 + t58) - g(2) * t34) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t33 * t28, t33 * t25, -t7, g(3) * t36 + t61 * (pkin(3) * t29 + pkin(7) * t26), 0, 0, 0, 0, 0, 0, -t33 * t18, t33 * t17, -t7, g(3) * t34 + t61 * (t16 * t29 - t26 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, g(1) * t9 - g(2) * t11 + t28 * t53, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t60 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t13;
