% Calculate inertial parameters regressor of gravitation load for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:10:55
% EndTime: 2019-05-04 22:10:57
% DurationCPUTime: 0.56s
% Computational Cost: add. (520->115), mult. (1141->181), div. (0->0), fcn. (1438->14), ass. (0->69)
t41 = sin(pkin(11));
t49 = sin(qJ(2));
t52 = cos(qJ(2));
t70 = cos(pkin(11));
t29 = -t52 * t41 - t49 * t70;
t45 = cos(pkin(6));
t27 = t29 * t45;
t42 = sin(pkin(10));
t44 = cos(pkin(10));
t59 = -t49 * t41 + t52 * t70;
t14 = -t44 * t27 + t42 * t59;
t15 = -t42 * t27 - t44 * t59;
t43 = sin(pkin(6));
t26 = t29 * t43;
t2 = g(1) * t15 - g(2) * t14 + g(3) * t26;
t48 = sin(qJ(4));
t84 = t15 * t48;
t83 = t26 * t48;
t40 = qJ(4) + pkin(12);
t39 = cos(t40);
t47 = sin(qJ(6));
t82 = t39 * t47;
t50 = cos(qJ(6));
t81 = t39 * t50;
t80 = t42 * t43;
t79 = t42 * t49;
t78 = t43 * t44;
t77 = t43 * t48;
t51 = cos(qJ(4));
t76 = t43 * t51;
t75 = t43 * t52;
t74 = t45 * t49;
t73 = t45 * t51;
t72 = t45 * t52;
t69 = t42 * t76;
t68 = t44 * t76;
t67 = t44 * t72;
t25 = t59 * t43;
t34 = pkin(2) * t75;
t37 = t51 * pkin(4) + pkin(3);
t46 = -qJ(5) - pkin(8);
t66 = t25 * t37 + t26 * t46 + t34;
t32 = pkin(2) * t67;
t64 = -pkin(2) * t79 + t32;
t38 = sin(t40);
t63 = pkin(5) * t39 + pkin(9) * t38;
t62 = -t14 * t48 - t68;
t61 = -t42 * t72 - t44 * t49;
t55 = t59 * t45;
t13 = t42 * t29 + t44 * t55;
t60 = t13 * t37 - t14 * t46 + t64;
t18 = t26 * t38 + t45 * t39;
t4 = -t14 * t38 - t39 * t78;
t6 = t15 * t38 + t39 * t80;
t58 = g(1) * t6 + g(2) * t4 + g(3) * t18;
t19 = -t26 * t39 + t45 * t38;
t5 = t14 * t39 - t38 * t78;
t7 = -t15 * t39 + t38 * t80;
t57 = g(1) * t7 + g(2) * t5 + g(3) * t19;
t16 = t44 * t29 - t42 * t55;
t3 = g(1) * t16 + g(2) * t13 + g(3) * t25;
t56 = t61 * pkin(2);
t54 = t15 * t46 + t16 * t37 + t56;
t53 = -g(1) * t61 - g(3) * t75;
t35 = pkin(4) * t73;
t31 = pkin(4) * t69;
t24 = -g(3) * t45 + (-g(1) * t42 + g(2) * t44) * t43;
t1 = t3 * t38;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2) * (t67 - t79) + t53, -g(1) * (t42 * t74 - t44 * t52) - g(2) * (-t42 * t52 - t44 * t74) + g(3) * t43 * t49, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t2, 0, -g(2) * t32 + (g(2) * t79 + t53) * pkin(2), 0, 0, 0, 0, 0, 0, -t3 * t51, t3 * t48, t2, -g(1) * (t16 * pkin(3) - t15 * pkin(8) + t56) - g(2) * (t13 * pkin(3) + pkin(8) * t14 + t64) - g(3) * (t25 * pkin(3) - t26 * pkin(8) + t34) 0, 0, 0, 0, 0, 0, -t3 * t39, t1, t2, -g(1) * t54 - g(2) * t60 - g(3) * t66, 0, 0, 0, 0, 0, 0, -g(1) * (-t15 * t47 + t16 * t81) - g(2) * (t13 * t81 + t14 * t47) - g(3) * (t25 * t81 - t26 * t47) -g(1) * (-t15 * t50 - t16 * t82) - g(2) * (-t13 * t82 + t14 * t50) - g(3) * (-t25 * t82 - t26 * t50) -t1, -g(1) * (t63 * t16 + t54) - g(2) * (t63 * t13 + t60) - g(3) * (t63 * t25 + t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t69 + t84) - g(2) * t62 - g(3) * (t73 + t83) -g(1) * (t15 * t51 - t42 * t77) - g(2) * (-t14 * t51 + t44 * t77) - g(3) * (t26 * t51 - t45 * t48) 0, 0, 0, 0, 0, 0, 0, 0, -t58, t57, 0, -g(1) * t31 - g(3) * t35 + (g(2) * t68 - t2 * t48) * pkin(4), 0, 0, 0, 0, 0, 0, -t58 * t50, t58 * t47, -t57, -g(1) * (pkin(4) * t84 + t6 * pkin(5) + t7 * pkin(9) + t31) - g(2) * (t62 * pkin(4) + t4 * pkin(5) + t5 * pkin(9)) - g(3) * (pkin(4) * t83 + t18 * pkin(5) + t19 * pkin(9) + t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t50 - t7 * t47) - g(2) * (-t13 * t50 - t5 * t47) - g(3) * (-t19 * t47 - t25 * t50) -g(1) * (t16 * t47 - t7 * t50) - g(2) * (t13 * t47 - t5 * t50) - g(3) * (-t19 * t50 + t25 * t47) 0, 0;];
taug_reg  = t8;
