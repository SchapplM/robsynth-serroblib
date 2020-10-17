% Calculate inertial parameters regressor of gravitation load for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x(4*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_gravloadJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:06:44
% EndTime: 2019-12-31 17:06:45
% DurationCPUTime: 0.15s
% Computational Cost: add. (111->43), mult. (165->64), div. (0->0), fcn. (161->8), ass. (0->30)
t19 = sin(qJ(1));
t22 = cos(qJ(1));
t8 = g(1) * t22 + g(2) * t19;
t15 = qJ(2) + pkin(7);
t11 = sin(t15);
t12 = cos(t15);
t24 = -g(3) * t12 + t8 * t11;
t32 = g(3) * t11;
t17 = sin(qJ(4));
t30 = t19 * t17;
t20 = cos(qJ(4));
t29 = t19 * t20;
t28 = t22 * t17;
t27 = t22 * t20;
t26 = t12 * pkin(3) + t11 * pkin(6);
t7 = g(1) * t19 - g(2) * t22;
t18 = sin(qJ(2));
t21 = cos(qJ(2));
t23 = -g(3) * t21 + t8 * t18;
t16 = -qJ(3) - pkin(5);
t13 = t21 * pkin(2);
t10 = t13 + pkin(1);
t9 = t22 * t10;
t6 = t12 * t27 + t30;
t5 = -t12 * t28 + t29;
t4 = -t12 * t29 + t28;
t3 = t12 * t30 + t27;
t2 = t7 * t11;
t1 = t8 * t12 + t32;
t14 = [0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t21, -t7 * t18, -t8, -g(1) * (-t19 * pkin(1) + t22 * pkin(5)) - g(2) * (t22 * pkin(1) + t19 * pkin(5)), 0, 0, 0, 0, 0, 0, t7 * t12, -t2, -t8, -g(1) * (-t19 * t10 - t22 * t16) - g(2) * (-t19 * t16 + t9), 0, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5, t2, -g(2) * t9 + (g(1) * t16 - g(2) * t26) * t22 + (-g(1) * (-t10 - t26) + g(2) * t16) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, g(3) * t18 + t8 * t21, 0, 0, 0, 0, 0, 0, 0, 0, t24, t1, 0, t23 * pkin(2), 0, 0, 0, 0, 0, 0, t24 * t20, -t24 * t17, -t1, -g(3) * (t13 + t26) + t8 * (pkin(2) * t18 + pkin(3) * t11 - pkin(6) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 + t17 * t32, g(1) * t6 - g(2) * t4 + t20 * t32, 0, 0;];
taug_reg = t14;
