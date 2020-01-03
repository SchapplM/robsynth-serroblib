% Calculate inertial parameters regressor of gravitation load for
% S4RRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRRP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_gravloadJ_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = cos(qJ(3));
t15 = t23 * pkin(3) + pkin(2);
t19 = -qJ(4) - pkin(6);
t21 = sin(qJ(2));
t24 = cos(qJ(2));
t27 = t24 * t15 - t21 * t19;
t22 = sin(qJ(1));
t25 = cos(qJ(1));
t14 = g(1) * t25 + g(2) * t22;
t20 = sin(qJ(3));
t34 = t25 * t20;
t11 = t22 * t23 - t24 * t34;
t39 = g(3) * t21;
t33 = t25 * t23;
t36 = t22 * t24;
t9 = t20 * t36 + t33;
t1 = -g(1) * t11 + g(2) * t9 + t20 * t39;
t7 = -g(3) * t24 + t14 * t21;
t42 = g(1) * t22;
t32 = t25 * pkin(1) + t22 * pkin(5);
t30 = t24 * pkin(2) + t21 * pkin(6);
t28 = -g(2) * t25 + t42;
t17 = t25 * pkin(5);
t13 = t28 * t21;
t12 = t22 * t20 + t24 * t33;
t10 = -t23 * t36 + t34;
t8 = t14 * t24 + t39;
t6 = t7 * t23;
t5 = t7 * t20;
t4 = -g(1) * t10 - g(2) * t12;
t3 = -g(1) * t9 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t10 + t23 * t39;
t16 = [0, 0, 0, 0, 0, 0, t28, t14, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t24, -t13, -t14, -g(1) * (-t22 * pkin(1) + t17) - g(2) * t32, 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * t17 - g(2) * (t30 * t25 + t32) - (-pkin(1) - t30) * t42, 0, 0, 0, 0, 0, 0, t4, t3, t13, -g(1) * (pkin(3) * t34 + t17) - g(2) * (t27 * t25 + t32) + (-g(1) * (-pkin(1) - t27) - g(2) * pkin(3) * t20) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t30 + t14 * (pkin(2) * t21 - pkin(6) * t24), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t27 + t14 * (t15 * t21 + t19 * t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg = t16;
