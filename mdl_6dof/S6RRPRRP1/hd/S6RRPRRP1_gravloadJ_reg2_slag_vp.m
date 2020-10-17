% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP1_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:19:17
% EndTime: 2019-05-06 17:19:18
% DurationCPUTime: 0.38s
% Computational Cost: add. (452->86), mult. (425->109), div. (0->0), fcn. (410->10), ass. (0->56)
t37 = qJ(2) + pkin(10);
t33 = qJ(4) + t37;
t27 = sin(t33);
t28 = cos(t33);
t73 = -pkin(4) * t27 + pkin(9) * t28;
t43 = cos(qJ(5));
t29 = t43 * pkin(5) + pkin(4);
t38 = -qJ(6) - pkin(9);
t72 = -t27 * t38 + t28 * t29;
t71 = t28 * pkin(4) + t27 * pkin(9);
t45 = cos(qJ(1));
t59 = t45 * t43;
t40 = sin(qJ(5));
t42 = sin(qJ(1));
t62 = t42 * t40;
t10 = t28 * t62 + t59;
t60 = t45 * t40;
t61 = t42 * t43;
t12 = -t28 * t60 + t61;
t65 = g(3) * t27;
t1 = -g(1) * t12 + g(2) * t10 + t40 * t65;
t23 = g(1) * t45 + g(2) * t42;
t7 = -g(3) * t28 + t23 * t27;
t32 = cos(t37);
t44 = cos(qJ(2));
t34 = t44 * pkin(2);
t57 = pkin(3) * t32 + t34;
t18 = pkin(1) + t57;
t14 = t45 * t18;
t66 = g(2) * t14;
t39 = -qJ(3) - pkin(7);
t36 = -pkin(8) + t39;
t55 = pkin(5) * t40 - t36;
t53 = t73 * t42;
t52 = t73 * t45;
t22 = g(1) * t42 - g(2) * t45;
t49 = t27 * t29 + t28 * t38;
t48 = t49 * t42;
t47 = t49 * t45;
t41 = sin(qJ(2));
t46 = -g(3) * t44 + t23 * t41;
t31 = sin(t37);
t30 = t34 + pkin(1);
t19 = -t41 * pkin(2) - pkin(3) * t31;
t16 = t45 * t19;
t15 = t42 * t19;
t13 = t28 * t59 + t62;
t11 = -t28 * t61 + t60;
t9 = t22 * t27;
t8 = t23 * t28 + t65;
t6 = t7 * t43;
t5 = t7 * t40;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t43 * t65;
t17 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t44, -t22 * t41, -t23, -g(1) * (-t42 * pkin(1) + t45 * pkin(7)) - g(2) * (t45 * pkin(1) + t42 * pkin(7)) 0, 0, 0, 0, 0, 0, t22 * t32, -t22 * t31, -t23, -g(1) * (-t42 * t30 - t45 * t39) - g(2) * (t45 * t30 - t42 * t39) 0, 0, 0, 0, 0, 0, t22 * t28, -t9, -t23, -g(1) * (-t42 * t18 - t45 * t36) - g(2) * (-t42 * t36 + t14) 0, 0, 0, 0, 0, 0, t4, t3, t9, -t66 + (g(1) * t36 - g(2) * t71) * t45 + (-g(1) * (-t18 - t71) + g(2) * t36) * t42, 0, 0, 0, 0, 0, 0, t4, t3, t9, -t66 + (-g(1) * t55 - g(2) * t72) * t45 + (-g(1) * (-t18 - t72) - g(2) * t55) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, g(3) * t41 + t23 * t44, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t32 + t23 * t31, g(3) * t31 + t23 * t32, 0, t46 * pkin(2), 0, 0, 0, 0, 0, 0, t7, t8, 0, -g(3) * t57 - t23 * t19, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t16 + t52) - g(2) * (t15 + t53) - g(3) * (t57 + t71) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t16 - t47) - g(2) * (t15 - t48) - g(3) * (t72 + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * t52 - g(2) * t53 - g(3) * t71, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, g(1) * t47 + g(2) * t48 - g(3) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t17;
