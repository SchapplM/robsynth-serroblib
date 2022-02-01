% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR6
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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:08:23
% EndTime: 2022-01-20 12:08:24
% DurationCPUTime: 0.20s
% Computational Cost: add. (286->50), mult. (217->61), div. (0->0), fcn. (194->10), ass. (0->42)
t39 = -pkin(8) - pkin(7);
t36 = sin(qJ(1));
t49 = t36 * pkin(1);
t34 = qJ(1) + qJ(2);
t26 = sin(t34);
t28 = cos(t34);
t48 = t28 * pkin(2) + t26 * pkin(7);
t33 = qJ(3) + qJ(4);
t27 = cos(t33);
t37 = cos(qJ(3));
t30 = t37 * pkin(3);
t47 = pkin(4) * t27 + t30;
t46 = -t26 * pkin(2) + t28 * pkin(7);
t14 = pkin(2) + t47;
t32 = pkin(9) - t39;
t45 = t28 * t14 + t26 * t32;
t44 = -t26 * t14 + t32 * t28;
t24 = t30 + pkin(2);
t43 = t28 * t24 - t26 * t39;
t13 = g(1) * t28 + g(2) * t26;
t12 = g(1) * t26 - g(2) * t28;
t38 = cos(qJ(1));
t42 = g(1) * t36 - g(2) * t38;
t41 = -t26 * t24 - t28 * t39;
t25 = sin(t33);
t3 = -g(3) * t27 + t13 * t25;
t35 = sin(qJ(3));
t40 = -g(3) * t37 + t13 * t35;
t31 = t38 * pkin(1);
t29 = qJ(5) + t33;
t23 = cos(t29);
t22 = sin(t29);
t10 = t12 * t37;
t9 = t12 * t35;
t8 = t12 * t27;
t7 = t12 * t25;
t6 = t12 * t23;
t5 = t12 * t22;
t4 = g(3) * t25 + t13 * t27;
t2 = g(3) * t22 + t13 * t23;
t1 = -g(3) * t23 + t13 * t22;
t11 = [0, 0, 0, 0, 0, 0, t42, g(1) * t38 + g(2) * t36, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, 0, t42 * pkin(1), 0, 0, 0, 0, 0, 0, t10, -t9, -t13, -g(1) * (t46 - t49) - g(2) * (t31 + t48), 0, 0, 0, 0, 0, 0, t8, -t7, -t13, -g(1) * (t41 - t49) - g(2) * (t31 + t43), 0, 0, 0, 0, 0, 0, t6, -t5, -t13, -g(1) * (t44 - t49) - g(2) * (t31 + t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t13, -g(1) * t46 - g(2) * t48, 0, 0, 0, 0, 0, 0, t8, -t7, -t13, -g(1) * t41 - g(2) * t43, 0, 0, 0, 0, 0, 0, t6, -t5, -t13, -g(1) * t44 - g(2) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(3) * t35 + t13 * t37, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t40 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t47 - t13 * (-t35 * pkin(3) - pkin(4) * t25); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t11;
