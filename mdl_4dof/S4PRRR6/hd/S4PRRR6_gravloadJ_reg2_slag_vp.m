% Calculate inertial parameters regressor of gravitation load for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t8 = sin(pkin(7));
t9 = cos(pkin(7));
t20 = g(1) * t9 + g(2) * t8;
t11 = sin(qJ(2));
t13 = cos(qJ(2));
t17 = -g(3) * t13 + t20 * t11;
t26 = g(3) * t11;
t24 = t13 * t8;
t23 = t13 * t9;
t10 = sin(qJ(3));
t22 = t10 * t13;
t12 = cos(qJ(3));
t21 = t12 * t13;
t15 = -g(1) * (t8 * t12 - t9 * t22) - g(2) * (-t9 * t12 - t8 * t22) + t10 * t26;
t14 = -pkin(6) - pkin(5);
t7 = qJ(3) + qJ(4);
t6 = cos(t7);
t5 = sin(t7);
t4 = t12 * pkin(3) + pkin(2);
t3 = t20 * t13 + t26;
t2 = -g(1) * (-t6 * t23 - t8 * t5) - g(2) * (-t6 * t24 + t9 * t5) + t6 * t26;
t1 = -g(1) * (-t5 * t23 + t8 * t6) - g(2) * (-t5 * t24 - t9 * t6) + t5 * t26;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t3, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t12, -t17 * t10, -t3, -g(3) * (t13 * pkin(2) + t11 * pkin(5)) + t20 * (pkin(2) * t11 - pkin(5) * t13), 0, 0, 0, 0, 0, 0, t17 * t6, -t17 * t5, -t3, -g(3) * (-t11 * t14 + t13 * t4) + t20 * (t11 * t4 + t13 * t14); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -g(1) * (-t8 * t10 - t9 * t21) - g(2) * (t9 * t10 - t8 * t21) + t12 * t26, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t15 * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t16;
