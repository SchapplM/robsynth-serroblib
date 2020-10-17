% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:23
% EndTime: 2020-01-03 12:13:23
% DurationCPUTime: 0.20s
% Computational Cost: add. (307->51), mult. (190->54), div. (0->0), fcn. (169->10), ass. (0->40)
t31 = qJ(1) + qJ(2);
t27 = qJ(3) + t31;
t20 = sin(t27);
t21 = cos(t27);
t34 = cos(qJ(4));
t22 = pkin(4) * t34 + pkin(3);
t36 = -pkin(9) - pkin(8);
t47 = t20 * t22 + t21 * t36;
t46 = t21 * pkin(3) + t20 * pkin(8);
t24 = sin(t31);
t18 = pkin(2) * t24;
t45 = t18 + t47;
t26 = cos(t31);
t19 = pkin(2) * t26;
t44 = t19 + t46;
t43 = t20 * pkin(3) - pkin(8) * t21;
t42 = -t20 * t36 + t21 * t22;
t41 = t18 + t43;
t40 = t19 + t42;
t39 = g(2) * t21 + g(3) * t20;
t7 = g(2) * t20 - g(3) * t21;
t10 = -g(2) * t26 - g(3) * t24;
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t38 = -g(2) * t35 - g(3) * t33;
t32 = sin(qJ(4));
t37 = -g(1) * t34 + t7 * t32;
t30 = qJ(4) + qJ(5);
t29 = t35 * pkin(1);
t28 = t33 * pkin(1);
t25 = cos(t30);
t23 = sin(t30);
t9 = g(2) * t24 - g(3) * t26;
t6 = t39 * t34;
t5 = t39 * t32;
t4 = t39 * t25;
t3 = t39 * t23;
t2 = -g(1) * t25 + t7 * t23;
t1 = g(1) * t23 + t7 * t25;
t8 = [0, 0, 0, 0, 0, 0, t38, g(2) * t33 - g(3) * t35, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, t38 * pkin(1), 0, 0, 0, 0, 0, 0, -t39, t7, 0, -g(2) * (t19 + t29) - g(3) * (t18 + t28), 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * (t29 + t44) - g(3) * (t28 + t41), 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(2) * (t29 + t40) - g(3) * (t28 + t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t7, 0, t10 * pkin(2), 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * t44 - g(3) * t41, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(2) * t40 - g(3) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * t46 - g(3) * t43, 0, 0, 0, 0, 0, 0, -t4, t3, -t7, -g(2) * t42 - g(3) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(1) * t32 + t7 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t37 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0;];
taug_reg = t8;
