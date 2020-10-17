% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRP3
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:33:03
% EndTime: 2019-05-06 17:33:05
% DurationCPUTime: 0.45s
% Computational Cost: add. (435->93), mult. (479->122), div. (0->0), fcn. (480->10), ass. (0->63)
t40 = sin(qJ(1));
t43 = cos(qJ(1));
t22 = g(1) * t43 + g(2) * t40;
t35 = qJ(2) + pkin(10);
t28 = cos(t35);
t41 = cos(qJ(4));
t57 = t43 * t41;
t38 = sin(qJ(4));
t62 = t40 * t38;
t14 = t28 * t62 + t57;
t58 = t43 * t38;
t61 = t40 * t41;
t16 = -t28 * t58 + t61;
t27 = sin(t35);
t68 = g(3) * t27;
t77 = -g(1) * t16 + g(2) * t14 + t38 * t68;
t36 = qJ(4) + qJ(5);
t29 = sin(t36);
t60 = t43 * t29;
t30 = cos(t36);
t63 = t40 * t30;
t11 = -t28 * t60 + t63;
t59 = t43 * t30;
t64 = t40 * t29;
t9 = t28 * t64 + t59;
t1 = -g(1) * t11 + g(2) * t9 + t29 * t68;
t7 = -g(3) * t28 + t22 * t27;
t44 = -pkin(9) - pkin(8);
t39 = sin(qJ(2));
t75 = pkin(2) * t39;
t42 = cos(qJ(2));
t32 = t42 * pkin(2);
t26 = t32 + pkin(1);
t23 = t43 * t26;
t70 = g(2) * t23;
t66 = t38 * pkin(4);
t19 = pkin(5) * t29 + t66;
t65 = t19 * t28;
t37 = -qJ(3) - pkin(7);
t56 = t19 - t37;
t31 = t41 * pkin(4);
t20 = pkin(5) * t30 + t31;
t53 = -t37 + t66;
t52 = t28 * pkin(3) + t27 * pkin(8);
t21 = g(1) * t40 - g(2) * t43;
t18 = pkin(3) + t20;
t34 = -qJ(6) + t44;
t51 = t28 * t18 - t27 * t34;
t25 = t31 + pkin(3);
t50 = t28 * t25 - t27 * t44;
t45 = -g(3) * t42 + t22 * t39;
t17 = t28 * t57 + t62;
t15 = -t28 * t61 + t58;
t13 = t21 * t27;
t12 = t28 * t59 + t64;
t10 = -t28 * t63 + t60;
t8 = t22 * t28 + t68;
t6 = t7 * t30;
t5 = t7 * t29;
t4 = -g(1) * t10 - g(2) * t12;
t3 = -g(1) * t9 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t10 + t30 * t68;
t24 = [0, 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t42, -t21 * t39, -t22, -g(1) * (-t40 * pkin(1) + t43 * pkin(7)) - g(2) * (t43 * pkin(1) + t40 * pkin(7)) 0, 0, 0, 0, 0, 0, t21 * t28, -t13, -t22, -g(1) * (-t40 * t26 - t43 * t37) - g(2) * (-t40 * t37 + t23) 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t13, -t70 + (g(1) * t37 - g(2) * t52) * t43 + (-g(1) * (-t26 - t52) + g(2) * t37) * t40, 0, 0, 0, 0, 0, 0, t4, t3, t13, -t70 + (-g(1) * t53 - g(2) * t50) * t43 + (-g(1) * (-t26 - t50) - g(2) * t53) * t40, 0, 0, 0, 0, 0, 0, t4, t3, t13, -t70 + (-g(1) * t56 - g(2) * t51) * t43 + (-g(1) * (-t26 - t51) - g(2) * t56) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(3) * t39 + t22 * t42, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t45 * pkin(2), 0, 0, 0, 0, 0, 0, t7 * t41, -t7 * t38, -t8, -g(3) * (t32 + t52) + t22 * (pkin(3) * t27 - pkin(8) * t28 + t75) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t32 + t50) + t22 * (t25 * t27 + t28 * t44 + t75) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t32 + t51) + t22 * (t18 * t27 + t28 * t34 + t75); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, g(1) * t17 - g(2) * t15 + t41 * t68, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t77 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t40 * t20 - t43 * t65) - g(2) * (-t43 * t20 - t40 * t65) + t19 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t24;
