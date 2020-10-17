% Calculate inertial parameters regressor of gravitation load for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PPRRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:05:13
% EndTime: 2019-05-04 20:05:15
% DurationCPUTime: 0.50s
% Computational Cost: add. (730->110), mult. (1947->176), div. (0->0), fcn. (2518->16), ass. (0->69)
t70 = sin(pkin(12));
t71 = sin(pkin(11));
t55 = t71 * t70;
t74 = cos(pkin(12));
t75 = cos(pkin(11));
t62 = t75 * t74;
t77 = cos(pkin(6));
t46 = -t62 * t77 + t55;
t72 = sin(pkin(7));
t73 = sin(pkin(6));
t59 = t73 * t72;
t76 = cos(pkin(7));
t85 = t46 * t76 + t75 * t59;
t57 = t71 * t74;
t60 = t75 * t70;
t47 = t57 * t77 + t60;
t56 = t71 * t73;
t84 = t47 * t76 - t72 * t56;
t83 = t74 * t76 * t73 + t77 * t72;
t82 = cos(qJ(3));
t34 = pkin(13) + qJ(6);
t32 = sin(t34);
t40 = cos(qJ(4));
t81 = t32 * t40;
t33 = cos(t34);
t80 = t33 * t40;
t35 = sin(pkin(13));
t79 = t35 * t40;
t36 = cos(pkin(13));
t78 = t36 * t40;
t69 = pkin(5) * t35 + pkin(9);
t25 = t60 * t77 + t57;
t39 = sin(qJ(3));
t11 = t25 * t82 - t85 * t39;
t10 = t25 * t39 + t85 * t82;
t8 = t10 * pkin(3);
t68 = t11 * pkin(9) - t8;
t26 = -t55 * t77 + t62;
t13 = t26 * t82 - t84 * t39;
t12 = t26 * t39 + t84 * t82;
t9 = t12 * pkin(3);
t67 = t13 * pkin(9) - t9;
t58 = t73 * t70;
t19 = t39 * t58 - t83 * t82;
t18 = t19 * pkin(3);
t20 = t83 * t39 + t82 * t58;
t66 = t20 * pkin(9) - t18;
t38 = sin(qJ(4));
t65 = -pkin(4) * t40 - qJ(5) * t38;
t31 = t36 * pkin(5) + pkin(4);
t37 = -pkin(10) - qJ(5);
t64 = -t31 * t40 + t37 * t38;
t61 = t75 * t73;
t45 = -t59 * t74 + t76 * t77;
t14 = t20 * t38 - t40 * t45;
t41 = t46 * t72 - t61 * t76;
t4 = t11 * t38 - t40 * t41;
t42 = t47 * t72 + t56 * t76;
t6 = t13 * t38 - t40 * t42;
t54 = g(1) * t6 + g(2) * t4 + g(3) * t14;
t15 = t20 * t40 + t38 * t45;
t5 = t11 * t40 + t38 * t41;
t7 = t13 * t40 + t38 * t42;
t53 = g(1) * t7 + g(2) * t5 + g(3) * t15;
t52 = g(1) * t12 + g(2) * t10 + g(3) * t19;
t51 = g(1) * t13 + g(2) * t11 + g(3) * t20;
t23 = -g(1) * t56 + g(2) * t61 - g(3) * t77;
t3 = t52 * t38;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, t51, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t40, -t3, -t51, -g(1) * t67 - g(2) * t68 - g(3) * t66, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t78 + t13 * t35) - g(2) * (-t10 * t78 + t11 * t35) - g(3) * (-t19 * t78 + t20 * t35) -g(1) * (t12 * t79 + t13 * t36) - g(2) * (t10 * t79 + t11 * t36) - g(3) * (t19 * t79 + t20 * t36) t3, -g(1) * (t12 * t65 + t67) - g(2) * (t10 * t65 + t68) - g(3) * (t19 * t65 + t66) 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t80 + t13 * t32) - g(2) * (-t10 * t80 + t11 * t32) - g(3) * (-t19 * t80 + t20 * t32) -g(1) * (t12 * t81 + t13 * t33) - g(2) * (t10 * t81 + t11 * t33) - g(3) * (t19 * t81 + t20 * t33) t3, -g(1) * (t12 * t64 + t13 * t69 - t9) - g(2) * (t10 * t64 + t11 * t69 - t8) - g(3) * (t19 * t64 + t20 * t69 - t18); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t53, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t36, -t54 * t35, -t53, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (-t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t14 * pkin(4) + t15 * qJ(5)) 0, 0, 0, 0, 0, 0, t54 * t33, -t54 * t32, -t53, -g(1) * (-t6 * t31 - t7 * t37) - g(2) * (-t4 * t31 - t5 * t37) - g(3) * (-t14 * t31 - t15 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t12 * t33 - t7 * t32) - g(2) * (t10 * t33 - t5 * t32) - g(3) * (-t15 * t32 + t19 * t33) -g(1) * (-t12 * t32 - t7 * t33) - g(2) * (-t10 * t32 - t5 * t33) - g(3) * (-t15 * t33 - t19 * t32) 0, 0;];
taug_reg  = t1;
