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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:02:06
% EndTime: 2022-01-20 12:02:07
% DurationCPUTime: 0.20s
% Computational Cost: add. (307->48), mult. (190->54), div. (0->0), fcn. (169->10), ass. (0->40)
t26 = qJ(1) + qJ(2);
t20 = sin(t26);
t43 = pkin(2) * t20;
t28 = sin(qJ(1));
t42 = t28 * pkin(1);
t23 = qJ(3) + t26;
t16 = sin(t23);
t17 = cos(t23);
t41 = t17 * pkin(3) + t16 * pkin(8);
t22 = cos(t26);
t15 = pkin(2) * t22;
t40 = t15 + t41;
t39 = -t16 * pkin(3) + t17 * pkin(8);
t29 = cos(qJ(4));
t18 = t29 * pkin(4) + pkin(3);
t31 = -pkin(9) - pkin(8);
t38 = -t16 * t31 + t17 * t18;
t37 = t15 + t38;
t8 = g(1) * t17 + g(2) * t16;
t7 = g(1) * t16 - g(2) * t17;
t9 = g(1) * t20 - g(2) * t22;
t30 = cos(qJ(1));
t36 = g(1) * t28 - g(2) * t30;
t35 = -t16 * t18 - t17 * t31;
t34 = t39 - t43;
t33 = t35 - t43;
t27 = sin(qJ(4));
t32 = -g(3) * t29 + t8 * t27;
t25 = qJ(4) + qJ(5);
t24 = t30 * pkin(1);
t21 = cos(t25);
t19 = sin(t25);
t10 = g(1) * t22 + g(2) * t20;
t6 = t7 * t29;
t5 = t7 * t27;
t4 = t7 * t21;
t3 = t7 * t19;
t2 = g(3) * t19 + t8 * t21;
t1 = -g(3) * t21 + t8 * t19;
t11 = [0, 0, 0, 0, 0, 0, t36, g(1) * t30 + g(2) * t28, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t36 * pkin(1), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(1) * (-t42 - t43) - g(2) * (t15 + t24), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t34 - t42) - g(2) * (t24 + t40), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t33 - t42) - g(2) * (t24 + t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t9 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t34 - g(2) * t40, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t33 - g(2) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t39 - g(2) * t41, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t35 - g(2) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t27 + t8 * t29, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t32 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t11;
