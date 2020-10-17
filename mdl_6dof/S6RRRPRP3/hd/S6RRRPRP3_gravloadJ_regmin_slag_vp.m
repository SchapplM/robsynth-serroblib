% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:39:45
% EndTime: 2019-05-07 07:39:46
% DurationCPUTime: 0.35s
% Computational Cost: add. (439->86), mult. (459->113), div. (0->0), fcn. (470->10), ass. (0->60)
t34 = qJ(2) + qJ(3);
t30 = sin(t34);
t31 = cos(t34);
t39 = sin(qJ(1));
t41 = cos(qJ(1));
t51 = g(1) * t41 + g(2) * t39;
t7 = -g(3) * t31 + t51 * t30;
t33 = pkin(10) + qJ(5);
t28 = sin(t33);
t29 = cos(t33);
t76 = pkin(5) * t29 + qJ(6) * t28;
t75 = t31 * pkin(3) + t30 * qJ(4);
t38 = sin(qJ(2));
t74 = pkin(2) * t38;
t73 = pkin(3) * t30;
t69 = g(3) * t30;
t37 = -pkin(9) - qJ(4);
t67 = t30 * t37;
t36 = cos(pkin(10));
t24 = t36 * pkin(4) + pkin(3);
t15 = t31 * t24;
t66 = t39 * t28;
t65 = t39 * t29;
t35 = sin(pkin(10));
t64 = t39 * t35;
t63 = t39 * t36;
t62 = t41 * t28;
t61 = t41 * t29;
t60 = t41 * t35;
t59 = t41 * t36;
t57 = qJ(4) * t31;
t55 = t76 * t31 + t15;
t42 = -pkin(8) - pkin(7);
t54 = pkin(4) * t35 - t42;
t11 = t31 * t62 - t65;
t9 = t31 * t66 + t61;
t53 = g(1) * t9 - g(2) * t11;
t52 = -t73 - t74;
t50 = g(1) * t39 - g(2) * t41;
t48 = t15 - t67;
t47 = t24 + t76;
t45 = t51 * t31;
t1 = g(1) * t11 + g(2) * t9 + t28 * t69;
t10 = t31 * t65 - t62;
t12 = t31 * t61 + t66;
t44 = g(1) * t12 + g(2) * t10 + t29 * t69;
t40 = cos(qJ(2));
t32 = t40 * pkin(2);
t27 = t32 + pkin(1);
t20 = t41 * t27;
t19 = t41 * t57;
t18 = t39 * t57;
t13 = t50 * t30;
t8 = t45 + t69;
t6 = t7 * t36;
t5 = t7 * t35;
t4 = t7 * t29;
t3 = t7 * t28;
t2 = g(1) * t10 - g(2) * t12;
t14 = [0, t50, t51, 0, 0, 0, 0, 0, t50 * t40, -t50 * t38, 0, 0, 0, 0, 0, t50 * t31, -t13, -g(1) * (-t31 * t63 + t60) - g(2) * (t31 * t59 + t64) -g(1) * (t31 * t64 + t59) - g(2) * (-t31 * t60 + t63) t13, -g(2) * t20 + (g(1) * t42 - g(2) * t75) * t41 + (-g(1) * (-t27 - t75) + g(2) * t42) * t39, 0, 0, 0, 0, 0, t2, -t53, t2, t13, t53, -g(1) * (-t10 * pkin(5) - t9 * qJ(6)) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t20) + (-g(1) * t54 - g(2) * t48) * t41 + (-g(1) * (-t27 - t48) - g(2) * t54) * t39; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t40 + t51 * t38, g(3) * t38 + t51 * t40, 0, 0, 0, 0, 0, t7, t8, t6, -t5, -t8, -g(1) * (t52 * t41 + t19) - g(2) * (t52 * t39 + t18) - g(3) * (t32 + t75) 0, 0, 0, 0, 0, t4, -t3, t4, -t8, t3, -g(3) * (t32 + t55 - t67) + t51 * (t47 * t30 + t31 * t37 + t74); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t6, -t5, -t8, -g(1) * (-t41 * t73 + t19) - g(2) * (-t39 * t73 + t18) - g(3) * t75, 0, 0, 0, 0, 0, t4, -t3, t4, -t8, t3, -g(3) * t55 + t37 * t45 + (g(3) * t37 + t51 * t47) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t44, t1, 0, -t44, -g(1) * (-t11 * pkin(5) + t12 * qJ(6)) - g(2) * (-t9 * pkin(5) + t10 * qJ(6)) - (-pkin(5) * t28 + qJ(6) * t29) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t14;
