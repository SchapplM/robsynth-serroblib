% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:55
% EndTime: 2022-01-20 11:02:57
% DurationCPUTime: 0.19s
% Computational Cost: add. (239->48), mult. (174->52), div. (0->0), fcn. (157->10), ass. (0->37)
t35 = sin(qJ(1));
t45 = pkin(1) * t35;
t33 = cos(pkin(9));
t19 = t33 * pkin(3) + pkin(2);
t34 = -pkin(7) - qJ(3);
t31 = qJ(1) + qJ(2);
t25 = sin(t31);
t26 = cos(t31);
t44 = t26 * pkin(2) + t25 * qJ(3);
t30 = pkin(9) + qJ(4);
t43 = -pkin(2) * t25 + t26 * qJ(3);
t23 = cos(t30);
t12 = pkin(4) * t23 + t19;
t29 = pkin(8) - t34;
t42 = t26 * t12 + t25 * t29;
t41 = -t12 * t25 + t29 * t26;
t40 = t26 * t19 - t25 * t34;
t11 = g(1) * t26 + g(2) * t25;
t10 = g(1) * t25 - g(2) * t26;
t36 = cos(qJ(1));
t39 = g(1) * t35 - g(2) * t36;
t38 = -t19 * t25 - t26 * t34;
t22 = sin(t30);
t37 = -g(3) * t23 + t11 * t22;
t28 = t36 * pkin(1);
t24 = qJ(5) + t30;
t18 = cos(t24);
t17 = sin(t24);
t8 = t10 * t33;
t7 = t10 * sin(pkin(9));
t6 = t10 * t23;
t5 = t10 * t22;
t4 = t10 * t18;
t3 = t10 * t17;
t2 = g(3) * t17 + t11 * t18;
t1 = -g(3) * t18 + t11 * t17;
t9 = [0, 0, 0, 0, 0, 0, t39, g(1) * t36 + g(2) * t35, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, t39 * pkin(1), 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * (t43 - t45) - g(2) * (t28 + t44), 0, 0, 0, 0, 0, 0, t6, -t5, -t11, -g(1) * (t38 - t45) - g(2) * (t28 + t40), 0, 0, 0, 0, 0, 0, t4, -t3, -t11, -g(1) * (t41 - t45) - g(2) * (t28 + t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t11, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t11, -g(1) * t43 - g(2) * t44, 0, 0, 0, 0, 0, 0, t6, -t5, -t11, -g(1) * t38 - g(2) * t40, 0, 0, 0, 0, 0, 0, t4, -t3, -t11, -g(1) * t41 - g(2) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(3) * t22 + t11 * t23, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t37 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
