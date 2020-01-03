% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t29 = cos(qJ(1));
t26 = sin(qJ(1));
t48 = g(1) * t26;
t51 = g(2) * t29 - t48;
t25 = sin(qJ(3));
t24 = sin(qJ(4));
t39 = t29 * t24;
t27 = cos(qJ(4));
t41 = t26 * t27;
t11 = t25 * t39 + t41;
t28 = cos(qJ(3));
t44 = g(3) * t28;
t38 = t29 * t27;
t42 = t26 * t24;
t9 = -t25 * t42 + t38;
t1 = -g(1) * t9 - g(2) * t11 + t24 * t44;
t8 = -g(3) * t25 - t28 * t51;
t49 = -pkin(1) - pkin(6);
t43 = t25 * t29;
t40 = t28 * t29;
t37 = t29 * pkin(1) + t26 * qJ(2);
t35 = t29 * pkin(6) + t37;
t34 = g(2) * t35;
t32 = t25 * pkin(3) - t28 * pkin(7);
t15 = g(1) * t29 + g(2) * t26;
t17 = t27 * pkin(4) + pkin(3);
t23 = -qJ(5) - pkin(7);
t30 = t25 * t17 + t28 * t23;
t19 = t29 * qJ(2);
t13 = t15 * t28;
t12 = t25 * t38 - t42;
t10 = t25 * t41 + t39;
t7 = -g(2) * t43 + t25 * t48 + t44;
t6 = t8 * t27;
t5 = t8 * t24;
t4 = -g(1) * t12 - g(2) * t10;
t3 = g(1) * t11 - g(2) * t9;
t2 = g(1) * t10 - g(2) * t12 + t27 * t44;
t14 = [0, 0, 0, 0, 0, 0, -t51, t15, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t15, -g(1) * (-t26 * pkin(1) + t19) - g(2) * t37, 0, 0, 0, 0, 0, 0, -t15 * t25, -t13, -t51, -g(1) * (t26 * t49 + t19) - t34, 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * (pkin(3) * t43 - pkin(7) * t40 + t19) - t34 + (-g(1) * t49 - g(2) * t32) * t26, 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * (t17 * t43 + t23 * t40 + t19) - g(2) * (pkin(4) * t39 + t35) + (-g(1) * (-pkin(4) * t24 + t49) - g(2) * t30) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, g(3) * t32 + t51 * (pkin(3) * t28 + pkin(7) * t25), 0, 0, 0, 0, 0, 0, -t6, t5, -t7, g(3) * t30 + t51 * (t17 * t28 - t23 * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8;];
taug_reg = t14;
