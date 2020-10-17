% Calculate inertial parameters regressor of gravitation load for
% S5PPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:09:10
% EndTime: 2019-12-05 15:09:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (133->39), mult. (175->51), div. (0->0), fcn. (180->6), ass. (0->32)
t18 = pkin(8) + qJ(3);
t16 = sin(t18);
t19 = sin(pkin(7));
t20 = cos(pkin(7));
t27 = g(1) * t20 + g(2) * t19;
t38 = t27 * t16;
t37 = pkin(3) * t16;
t17 = cos(t18);
t36 = pkin(6) * t17;
t33 = g(3) * t16;
t21 = sin(qJ(4));
t32 = t19 * t21;
t22 = cos(qJ(4));
t31 = t19 * t22;
t30 = t20 * t21;
t29 = t20 * t22;
t28 = t17 * pkin(3) + t16 * pkin(6);
t26 = pkin(4) * t22 + qJ(5) * t21;
t5 = t17 * t32 + t29;
t7 = t17 * t30 - t31;
t1 = g(1) * t7 + g(2) * t5 + t21 * t33;
t6 = t17 * t31 - t30;
t8 = t17 * t29 + t32;
t24 = g(1) * t8 + g(2) * t6 + t22 * t33;
t23 = -g(3) * t17 + t38;
t11 = t20 * t36;
t10 = t19 * t36;
t9 = -g(1) * t19 + g(2) * t20;
t4 = t27 * t17 + t33;
t3 = t23 * t22;
t2 = t23 * t21;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, -t4, -g(1) * (-t20 * t37 + t11) - g(2) * (-t19 * t37 + t10) - g(3) * t28, 0, 0, 0, 0, 0, 0, t3, -t4, t2, -g(1) * t11 - g(2) * t10 - g(3) * (t26 * t17 + t28) + (pkin(3) + t26) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t24, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t24, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - (-pkin(4) * t21 + qJ(5) * t22) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t12;
