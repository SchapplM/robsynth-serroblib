% Calculate inertial parameters regressor of gravitation load for
% S5PPRRP3
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PPRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:11:16
% EndTime: 2019-12-05 15:11:17
% DurationCPUTime: 0.18s
% Computational Cost: add. (119->46), mult. (301->66), div. (0->0), fcn. (355->8), ass. (0->39)
t23 = sin(pkin(8));
t27 = sin(qJ(4));
t45 = t23 * t27;
t28 = sin(qJ(3));
t44 = t23 * t28;
t29 = cos(qJ(4));
t43 = t23 * t29;
t30 = cos(qJ(3));
t42 = t23 * t30;
t24 = sin(pkin(7));
t41 = t24 * t28;
t40 = t24 * t30;
t26 = cos(pkin(7));
t39 = t26 * t28;
t38 = t26 * t30;
t37 = g(3) * t44;
t25 = cos(pkin(8));
t11 = -t25 * t41 - t38;
t12 = t25 * t40 - t39;
t36 = t11 * pkin(3) + t12 * pkin(6);
t13 = -t25 * t39 + t40;
t14 = t25 * t38 + t41;
t35 = t13 * pkin(3) + t14 * pkin(6);
t34 = pkin(4) * t29 + qJ(5) * t27;
t15 = t25 * t29 + t27 * t42;
t5 = t12 * t27 - t24 * t43;
t7 = t14 * t27 - t26 * t43;
t1 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t16 = -t25 * t27 + t29 * t42;
t6 = t12 * t29 + t24 * t45;
t8 = t14 * t29 + t26 * t45;
t33 = g(1) * t8 + g(2) * t6 + g(3) * t16;
t32 = -g(1) * t13 - g(2) * t11 + t37;
t31 = g(1) * t14 + g(2) * t12 + g(3) * t42;
t20 = pkin(6) * t42;
t17 = -g(1) * t24 + g(2) * t26;
t3 = t32 * t29;
t2 = t32 * t27;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t31, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, -t31, -g(1) * t35 - g(2) * t36 - g(3) * (-pkin(3) * t44 + t20), 0, 0, 0, 0, 0, 0, t3, -t31, t2, -g(1) * (t34 * t13 + t35) - g(2) * (t34 * t11 + t36) - g(3) * t20 - (-pkin(3) - t34) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t33, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t33, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - g(3) * (-t15 * pkin(4) + t16 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t4;
