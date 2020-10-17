% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:37
% EndTime: 2020-01-03 12:07:37
% DurationCPUTime: 0.16s
% Computational Cost: add. (274->44), mult. (138->47), div. (0->0), fcn. (120->10), ass. (0->35)
t27 = qJ(1) + qJ(2);
t22 = sin(t27);
t17 = pkin(2) * t22;
t29 = sin(qJ(1));
t25 = t29 * pkin(1);
t39 = t17 + t25;
t23 = cos(t27);
t18 = pkin(2) * t23;
t31 = cos(qJ(1));
t26 = t31 * pkin(1);
t38 = t18 + t26;
t24 = qJ(3) + t27;
t19 = pkin(9) + t24;
t13 = sin(t19);
t14 = cos(t19);
t21 = cos(t24);
t16 = pkin(3) * t21;
t37 = t14 * pkin(4) + t13 * pkin(8) + t16;
t36 = t18 + t37;
t20 = sin(t24);
t15 = pkin(3) * t20;
t35 = t13 * pkin(4) - t14 * pkin(8) + t15;
t34 = g(2) * t14 + g(3) * t13;
t3 = g(2) * t13 - g(3) * t14;
t6 = -g(2) * t21 - g(3) * t20;
t8 = -g(2) * t23 - g(3) * t22;
t33 = -g(2) * t31 - g(3) * t29;
t32 = t17 + t35;
t30 = cos(qJ(5));
t28 = sin(qJ(5));
t7 = g(2) * t22 - g(3) * t23;
t5 = g(2) * t20 - g(3) * t21;
t2 = t34 * t30;
t1 = t34 * t28;
t4 = [0, 0, 0, 0, 0, 0, t33, g(2) * t29 - g(3) * t31, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, t33 * pkin(1), 0, 0, 0, 0, 0, 0, t6, t5, 0, -g(2) * t38 - g(3) * t39, 0, 0, 0, 0, 0, 0, -t34, t3, 0, -g(2) * (t16 + t38) - g(3) * (t15 + t39), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * (t26 + t36) - g(3) * (t25 + t32); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t7, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t8 * pkin(2), 0, 0, 0, 0, 0, 0, -t34, t3, 0, -g(2) * (t16 + t18) - g(3) * (t15 + t17), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * t36 - g(3) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t34, t3, 0, t6 * pkin(3), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(2) * t37 - g(3) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t30 + t3 * t28, g(1) * t28 + t3 * t30, 0, 0;];
taug_reg = t4;
