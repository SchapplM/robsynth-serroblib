% Calculate inertial parameters regressor of gravitation load for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t17 = qJ(1) + pkin(8);
t15 = qJ(3) + t17;
t11 = sin(t15);
t12 = cos(t15);
t29 = t12 * pkin(3) + t11 * qJ(4);
t28 = t11 * pkin(3);
t14 = cos(t17);
t21 = cos(qJ(1));
t27 = t21 * pkin(1) + pkin(2) * t14;
t26 = (-pkin(3) - pkin(7)) * t11;
t25 = t27 + t29;
t13 = sin(t17);
t19 = sin(qJ(1));
t24 = -t19 * pkin(1) - pkin(2) * t13;
t4 = g(1) * t12 + g(2) * t11;
t3 = g(1) * t11 - g(2) * t12;
t23 = g(1) * t19 - g(2) * t21;
t6 = t12 * qJ(4);
t22 = t24 + t6;
t20 = cos(qJ(5));
t18 = sin(qJ(5));
t7 = t12 * pkin(7);
t2 = t4 * t20;
t1 = t4 * t18;
t5 = [0, 0, 0, 0, 0, 0, t23, g(1) * t21 + g(2) * t19, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t14, g(1) * t14 + g(2) * t13, 0, t23 * pkin(1), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t24 - g(2) * t27, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t22 - t28) - g(2) * t25, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (t26 + t22) - g(2) * (t7 + t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t6 - t28) - g(2) * t29, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (t6 + t26) - g(2) * (t7 + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t18 - t3 * t20, g(3) * t20 + t3 * t18, 0, 0;];
taug_reg = t5;
