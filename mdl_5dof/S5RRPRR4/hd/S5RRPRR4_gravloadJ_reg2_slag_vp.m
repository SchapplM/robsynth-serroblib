% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR4
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:32
% EndTime: 2022-01-20 10:48:33
% DurationCPUTime: 0.18s
% Computational Cost: add. (240->45), mult. (152->50), div. (0->0), fcn. (135->10), ass. (0->36)
t26 = qJ(1) + qJ(2);
t21 = sin(t26);
t39 = pkin(2) * t21;
t28 = sin(qJ(1));
t38 = t28 * pkin(1);
t19 = pkin(9) + t26;
t15 = sin(t19);
t16 = cos(t19);
t23 = cos(t26);
t17 = pkin(2) * t23;
t37 = t16 * pkin(3) + t15 * pkin(7) + t17;
t29 = cos(qJ(4));
t18 = t29 * pkin(4) + pkin(3);
t31 = -pkin(8) - pkin(7);
t36 = -t15 * t31 + t16 * t18 + t17;
t8 = g(1) * t16 + g(2) * t15;
t7 = g(1) * t15 - g(2) * t16;
t9 = g(1) * t21 - g(2) * t23;
t30 = cos(qJ(1));
t35 = g(1) * t28 - g(2) * t30;
t34 = -t15 * pkin(3) + t16 * pkin(7) - t39;
t33 = -t15 * t18 - t16 * t31 - t39;
t27 = sin(qJ(4));
t32 = -g(3) * t29 + t27 * t8;
t25 = qJ(4) + qJ(5);
t24 = t30 * pkin(1);
t22 = cos(t25);
t20 = sin(t25);
t10 = g(1) * t23 + g(2) * t21;
t6 = t7 * t29;
t5 = t7 * t27;
t4 = t7 * t22;
t3 = t7 * t20;
t2 = g(3) * t20 + t22 * t8;
t1 = -g(3) * t22 + t20 * t8;
t11 = [0, 0, 0, 0, 0, 0, t35, g(1) * t30 + g(2) * t28, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, t35 * pkin(1), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(1) * (-t38 - t39) - g(2) * (t17 + t24), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t34 - t38) - g(2) * (t24 + t37), 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * (t33 - t38) - g(2) * (t24 + t36); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t9 * pkin(2), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t34 - g(2) * t37, 0, 0, 0, 0, 0, 0, t4, -t3, -t8, -g(1) * t33 - g(2) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, g(3) * t27 + t29 * t8, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t32 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t11;
