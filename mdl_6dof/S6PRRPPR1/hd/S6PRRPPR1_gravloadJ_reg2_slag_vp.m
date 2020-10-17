% Calculate inertial parameters regressor of gravitation load for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:33:59
% EndTime: 2019-05-05 02:34:01
% DurationCPUTime: 0.64s
% Computational Cost: add. (508->128), mult. (858->190), div. (0->0), fcn. (1023->14), ass. (0->63)
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t61 = sin(pkin(10));
t63 = cos(pkin(6));
t50 = t63 * t61;
t62 = cos(pkin(10));
t19 = -t42 * t50 + t62 * t44;
t41 = sin(qJ(3));
t43 = cos(qJ(3));
t37 = sin(pkin(6));
t58 = t37 * t61;
t80 = -t19 * t41 + t43 * t58;
t68 = t37 * t42;
t79 = -t41 * t68 + t63 * t43;
t29 = t43 * pkin(3) + pkin(2);
t67 = t37 * t44;
t20 = t29 * t67;
t78 = g(3) * t20;
t77 = g(3) * t37;
t51 = t63 * t62;
t17 = t42 * t51 + t61 * t44;
t36 = sin(pkin(12));
t76 = t17 * t36;
t75 = t19 * t36;
t34 = pkin(12) + qJ(6);
t30 = sin(t34);
t35 = qJ(3) + pkin(11);
t33 = cos(t35);
t73 = t30 * t33;
t32 = cos(t34);
t72 = t32 * t33;
t71 = t33 * t36;
t38 = cos(pkin(12));
t70 = t33 * t38;
t69 = t33 * t44;
t39 = -qJ(4) - pkin(8);
t66 = t39 * t42;
t16 = t61 * t42 - t44 * t51;
t65 = -t16 * t29 - t17 * t39;
t18 = t62 * t42 + t44 * t50;
t64 = -t18 * t29 - t19 * t39;
t59 = t37 * t62;
t56 = t80 * pkin(3);
t31 = sin(t35);
t54 = pkin(4) * t33 + qJ(5) * t31;
t28 = t38 * pkin(5) + pkin(4);
t40 = -pkin(9) - qJ(5);
t53 = t28 * t33 - t31 * t40;
t52 = t79 * pkin(3);
t12 = t31 * t68 - t63 * t33;
t6 = t17 * t31 + t33 * t59;
t8 = t19 * t31 - t33 * t58;
t49 = g(1) * t8 + g(2) * t6 + g(3) * t12;
t13 = t63 * t31 + t33 * t68;
t7 = t17 * t33 - t31 * t59;
t9 = t19 * t33 + t31 * t58;
t48 = g(1) * t9 + g(2) * t7 + g(3) * t13;
t47 = -t17 * t41 - t43 * t59;
t4 = -g(1) * t18 - g(2) * t16 + g(3) * t67;
t46 = g(1) * t19 + g(2) * t17 + g(3) * t68;
t45 = t47 * pkin(3);
t3 = t4 * t31;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t46, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t43, t4 * t41, -t46, -g(1) * (-t18 * pkin(2) + t19 * pkin(8)) - g(2) * (-t16 * pkin(2) + t17 * pkin(8)) - (pkin(2) * t44 + pkin(8) * t42) * t77, 0, 0, 0, 0, 0, 0, -t4 * t33, t3, -t46, -g(1) * t64 - g(2) * t65 - g(3) * (-t37 * t66 + t20) 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t70 + t75) - g(2) * (-t16 * t70 + t76) - (t36 * t42 + t38 * t69) * t77, -g(1) * (t18 * t71 + t19 * t38) - g(2) * (t16 * t71 + t17 * t38) - (-t36 * t69 + t38 * t42) * t77, -t3, -g(1) * (-t18 * t54 + t64) - g(2) * (-t16 * t54 + t65) - t78 - (t44 * t54 - t66) * t77, 0, 0, 0, 0, 0, 0, -g(1) * (-t18 * t72 + t19 * t30) - g(2) * (-t16 * t72 + t17 * t30) - (t30 * t42 + t32 * t69) * t77, -g(1) * (t18 * t73 + t19 * t32) - g(2) * (t16 * t73 + t17 * t32) - (-t30 * t69 + t32 * t42) * t77, -t3, -g(1) * (pkin(5) * t75 - t18 * t53 + t64) - g(2) * (pkin(5) * t76 - t16 * t53 + t65) - t78 - (t53 * t44 + (pkin(5) * t36 - t39) * t42) * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t47 - g(3) * t79, -g(1) * (-t19 * t43 - t41 * t58) - g(2) * (-t17 * t43 + t41 * t59) - g(3) * (-t63 * t41 - t43 * t68) 0, 0, 0, 0, 0, 0, 0, 0, t49, t48, 0, -g(1) * t56 - g(2) * t45 - g(3) * t52, 0, 0, 0, 0, 0, 0, t49 * t38, -t49 * t36, -t48, -g(1) * (-t8 * pkin(4) + t9 * qJ(5) + t56) - g(2) * (-t6 * pkin(4) + t7 * qJ(5) + t45) - g(3) * (-t12 * pkin(4) + t13 * qJ(5) + t52) 0, 0, 0, 0, 0, 0, t49 * t32, -t49 * t30, -t48, -g(1) * (-t8 * t28 - t9 * t40 + t56) - g(2) * (-t6 * t28 - t7 * t40 + t45) - g(3) * (-t12 * t28 - t13 * t40 + t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t18 * t32 - t9 * t30) - g(2) * (t16 * t32 - t7 * t30) - g(3) * (-t13 * t30 - t32 * t67) -g(1) * (-t18 * t30 - t9 * t32) - g(2) * (-t16 * t30 - t7 * t32) - g(3) * (-t13 * t32 + t30 * t67) 0, 0;];
taug_reg  = t1;
