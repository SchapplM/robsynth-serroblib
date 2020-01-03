% Calculate inertial parameters regressor of gravitation load for
% S5RPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t31 = cos(qJ(1));
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t36 = t26 * t30 - t27 * t28;
t10 = t36 * t31;
t44 = t26 * t28;
t35 = t27 * t30 + t44;
t48 = -g(1) * t10 + g(3) * t35;
t29 = sin(qJ(1));
t46 = g(1) * t29;
t45 = g(2) * t29;
t43 = t27 * t31;
t42 = t31 * pkin(1) + t29 * qJ(2);
t41 = qJ(3) * t26;
t40 = pkin(4) * t44;
t39 = pkin(2) * t43 + t31 * t41 + t42;
t15 = g(1) * t31 + t45;
t14 = -g(2) * t31 + t46;
t25 = qJ(4) + qJ(5);
t19 = sin(t25);
t20 = cos(t25);
t38 = t27 * t19 - t26 * t20;
t37 = t26 * t19 + t27 * t20;
t34 = -pkin(2) * t27 - pkin(1) - t41;
t32 = -pkin(7) - pkin(6);
t22 = t31 * qJ(2);
t18 = t30 * pkin(4) + pkin(3);
t13 = t14 * t27;
t12 = t14 * t26;
t11 = t35 * t31;
t9 = t35 * t29;
t8 = t36 * t29;
t7 = g(3) * t27 - t15 * t26;
t6 = t37 * t31;
t5 = t38 * t31;
t4 = t37 * t29;
t3 = t38 * t29;
t2 = g(1) * t6 + g(2) * t4 - g(3) * t38;
t1 = g(1) * t5 + g(2) * t3 + g(3) * t37;
t16 = [0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, -t15, -g(1) * (-t29 * pkin(1) + t22) - g(2) * t42, 0, 0, 0, 0, 0, 0, t13, -t15, t12, -g(1) * t22 - g(2) * t39 - t34 * t46, 0, 0, 0, 0, 0, 0, g(1) * t9 - g(2) * t11, g(1) * t8 - g(2) * t10, t15, -g(1) * (-t31 * pkin(6) + t22) - g(2) * (pkin(3) * t43 + t39) + (-g(1) * (-pkin(3) * t27 + t34) + g(2) * pkin(6)) * t29, 0, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t6, -g(1) * t3 + g(2) * t5, t15, -g(1) * (t31 * t32 + t22) - g(2) * (t18 * t43 + t31 * t40 + t39) + (-g(1) * (-t18 * t27 + t34 - t40) - g(2) * t32) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * t8 + t48, g(1) * t11 + g(2) * t9 + g(3) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, (-t36 * t45 + t48) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t16;
