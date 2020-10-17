% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:36
% EndTime: 2019-12-31 19:26:36
% DurationCPUTime: 0.13s
% Computational Cost: add. (194->41), mult. (116->40), div. (0->0), fcn. (102->8), ass. (0->28)
t19 = qJ(1) + qJ(2);
t16 = sin(t19);
t31 = pkin(2) * t16;
t21 = sin(qJ(1));
t30 = t21 * pkin(1);
t15 = pkin(8) + t19;
t12 = sin(t15);
t13 = cos(t15);
t17 = cos(t19);
t14 = pkin(2) * t17;
t29 = t13 * pkin(3) + t12 * qJ(4) + t14;
t28 = t13 * qJ(4) - t31;
t23 = cos(qJ(1));
t18 = t23 * pkin(1);
t27 = t18 + t29;
t4 = g(1) * t13 + g(2) * t12;
t3 = g(1) * t12 - g(2) * t13;
t5 = g(1) * t16 - g(2) * t17;
t26 = g(1) * t21 - g(2) * t23;
t25 = -t12 * pkin(3) + t28;
t24 = (-pkin(3) - pkin(7)) * t12 + t28;
t22 = cos(qJ(5));
t20 = sin(qJ(5));
t9 = t13 * pkin(7);
t6 = g(1) * t17 + g(2) * t16;
t2 = t4 * t22;
t1 = t4 * t20;
t7 = [0, 0, 0, 0, 0, 0, t26, g(1) * t23 + g(2) * t21, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t26 * pkin(1), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * (-t30 - t31) - g(2) * (t14 + t18), 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * (t25 - t30) - g(2) * t27, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * (t24 - t30) - g(2) * (t9 + t27); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(2), 0, 0, 0, 0, 0, 0, 0, -t3, -t4, -g(1) * t25 - g(2) * t29, 0, 0, 0, 0, 0, 0, -t1, -t2, t3, -g(1) * t24 - g(2) * (t9 + t29); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t20 - t3 * t22, g(3) * t22 + t3 * t20, 0, 0;];
taug_reg = t7;
