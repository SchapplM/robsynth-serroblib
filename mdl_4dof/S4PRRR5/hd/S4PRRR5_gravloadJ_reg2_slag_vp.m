% Calculate inertial parameters regressor of gravitation load for
% S4PRRR5
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t11 = qJ(2) + qJ(3);
t10 = cos(t11);
t12 = sin(pkin(7));
t13 = cos(pkin(7));
t19 = g(1) * t13 + g(2) * t12;
t9 = sin(t11);
t3 = -g(3) * t10 + t19 * t9;
t29 = pkin(3) * t9;
t28 = g(3) * t9;
t27 = t10 * pkin(3) + t9 * pkin(6);
t26 = pkin(6) * t10;
t14 = sin(qJ(4));
t24 = t12 * t14;
t16 = cos(qJ(4));
t23 = t12 * t16;
t22 = t13 * t14;
t21 = t13 * t16;
t15 = sin(qJ(2));
t20 = -pkin(2) * t15 - t29;
t17 = cos(qJ(2));
t18 = -g(3) * t17 + t19 * t15;
t6 = t13 * t26;
t5 = t12 * t26;
t4 = t19 * t10 + t28;
t2 = t3 * t16;
t1 = t3 * t14;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, g(3) * t15 + t19 * t17, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t18 * pkin(2), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t20 * t13 + t6) - g(2) * (t20 * t12 + t5) - g(3) * (t17 * pkin(2) + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t13 * t29 + t6) - g(2) * (-t12 * t29 + t5) - g(3) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t22 + t23) - g(2) * (-t10 * t24 - t21) + t14 * t28, -g(1) * (-t10 * t21 - t24) - g(2) * (-t10 * t23 + t22) + t16 * t28, 0, 0;];
taug_reg = t7;
