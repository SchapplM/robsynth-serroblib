% Calculate inertial parameters regressor of gravitation load for
% S5PPRRP1
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:07:14
% EndTime: 2019-12-05 15:07:15
% DurationCPUTime: 0.14s
% Computational Cost: add. (117->31), mult. (149->38), div. (0->0), fcn. (148->6), ass. (0->23)
t12 = sin(pkin(7));
t13 = cos(pkin(7));
t18 = g(1) * t13 + g(2) * t12;
t11 = pkin(8) + qJ(3);
t10 = cos(t11);
t9 = sin(t11);
t5 = -g(3) * t10 + t18 * t9;
t27 = g(3) * t9;
t15 = sin(qJ(4));
t23 = t12 * t15;
t16 = cos(qJ(4));
t22 = t12 * t16;
t21 = t13 * t15;
t20 = t13 * t16;
t1 = -g(1) * (-t10 * t21 + t22) - g(2) * (-t10 * t23 - t20) + t15 * t27;
t14 = -qJ(5) - pkin(6);
t8 = t16 * pkin(4) + pkin(3);
t7 = -g(1) * t12 + g(2) * t13;
t6 = t18 * t10 + t27;
t4 = t5 * t16;
t3 = t5 * t15;
t2 = -g(1) * (-t10 * t20 - t23) - g(2) * (-t10 * t22 + t21) + t16 * t27;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (t10 * pkin(3) + t9 * pkin(6)) + t18 * (pkin(3) * t9 - pkin(6) * t10), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * (t10 * t8 - t9 * t14) + t18 * (t10 * t14 + t8 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t17;
