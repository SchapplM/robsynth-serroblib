% Calculate inertial parameters regressor of gravitation load for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t14 = qJ(1) + pkin(7);
t11 = sin(t14);
t12 = cos(t14);
t7 = g(1) * t12 + g(2) * t11;
t16 = sin(qJ(3));
t19 = cos(qJ(3));
t21 = -g(3) * t19 + t7 * t16;
t34 = g(1) * t11;
t31 = g(3) * t16;
t15 = sin(qJ(4));
t29 = t15 * t19;
t18 = cos(qJ(4));
t28 = t18 * t19;
t20 = cos(qJ(1));
t27 = t20 * pkin(1) + t12 * pkin(2) + t11 * pkin(5);
t17 = sin(qJ(1));
t26 = -t17 * pkin(1) + t12 * pkin(5);
t25 = t19 * pkin(3) + t16 * pkin(6);
t23 = -g(2) * t12 + t34;
t22 = g(1) * t17 - g(2) * t20;
t6 = t23 * t16;
t5 = t11 * t15 + t12 * t28;
t4 = t11 * t18 - t12 * t29;
t3 = -t11 * t28 + t12 * t15;
t2 = t11 * t29 + t12 * t18;
t1 = t7 * t19 + t31;
t8 = [0, 0, 0, 0, 0, 0, t22, g(1) * t20 + g(2) * t17, 0, 0, 0, 0, 0, 0, 0, 0, t23, t7, 0, t22 * pkin(1), 0, 0, 0, 0, 0, 0, t23 * t19, -t6, -t7, -g(1) * (-t11 * pkin(2) + t26) - g(2) * t27, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5, -g(1) * t2 - g(2) * t4, t6, -g(1) * t26 - g(2) * (t25 * t12 + t27) - (-pkin(2) - t25) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t1, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t18, -t21 * t15, -t1, -g(3) * t25 + t7 * (pkin(3) * t16 - pkin(6) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t2 + t15 * t31, g(1) * t5 - g(2) * t3 + t18 * t31, 0, 0;];
taug_reg = t8;
