% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:18:29
% EndTime: 2019-05-05 01:18:30
% DurationCPUTime: 0.50s
% Computational Cost: add. (431->111), mult. (804->159), div. (0->0), fcn. (961->12), ass. (0->67)
t41 = sin(pkin(11));
t43 = cos(pkin(11));
t47 = sin(qJ(2));
t44 = cos(pkin(6));
t50 = cos(qJ(2));
t71 = t44 * t50;
t27 = t41 * t47 - t43 * t71;
t49 = cos(qJ(4));
t42 = sin(pkin(6));
t46 = sin(qJ(4));
t76 = t42 * t46;
t90 = t27 * t49 + t43 * t76;
t72 = t44 * t47;
t28 = t41 * t50 + t43 * t72;
t51 = -pkin(9) - pkin(8);
t86 = pkin(4) * t46;
t89 = t27 * t51 + t28 * t86;
t29 = t41 * t71 + t43 * t47;
t30 = -t41 * t72 + t43 * t50;
t88 = t29 * t51 + t30 * t86;
t87 = -g(1) * t30 - g(2) * t28;
t83 = g(3) * t42;
t81 = t29 * t49;
t40 = qJ(4) + qJ(5);
t38 = sin(t40);
t45 = sin(qJ(6));
t80 = t38 * t45;
t48 = cos(qJ(6));
t79 = t38 * t48;
t78 = t41 * t42;
t77 = t42 * t43;
t75 = t42 * t47;
t74 = t42 * t49;
t73 = t42 * t50;
t70 = t45 * t47;
t69 = t47 * t48;
t68 = t50 * t51;
t67 = t90 * pkin(4);
t66 = pkin(2) * t73 + qJ(3) * t75;
t65 = t41 * t76;
t63 = t75 * t86 + t66;
t39 = cos(t40);
t11 = t29 * t39 - t38 * t78;
t12 = t29 * t38 + t39 * t78;
t62 = t11 * pkin(5) + t12 * pkin(10);
t13 = t27 * t39 + t38 * t77;
t14 = -t27 * t38 + t39 * t77;
t61 = t13 * pkin(5) - t14 * pkin(10);
t21 = -t44 * t38 - t39 * t73;
t22 = -t38 * t73 + t44 * t39;
t60 = t21 * pkin(5) + t22 * pkin(10);
t25 = t27 * pkin(2);
t59 = t28 * qJ(3) - t25;
t26 = t29 * pkin(2);
t58 = t30 * qJ(3) - t26;
t57 = pkin(5) * t38 - pkin(10) * t39;
t55 = -t44 * t46 - t49 * t73;
t54 = g(1) * t11 + g(2) * t13 + g(3) * t21;
t5 = g(1) * t12 - g(2) * t14 + g(3) * t22;
t53 = g(3) * t55;
t7 = -g(1) * t29 - g(2) * t27 + g(3) * t73;
t52 = g(3) * t75 - t87;
t19 = pkin(4) * t81;
t6 = t52 * t39;
t2 = t54 * t48;
t1 = t54 * t45;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t52, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t52, -g(1) * t58 - g(2) * t59 - g(3) * t66, 0, 0, 0, 0, 0, 0, -t52 * t46, -t52 * t49, -t7, -g(1) * (-t29 * pkin(8) + t58) - g(2) * (-t27 * pkin(8) + t59) - g(3) * (pkin(8) * t73 + t66) 0, 0, 0, 0, 0, 0, -t52 * t38, -t6, -t7, -g(1) * (t58 + t88) - g(2) * (t59 + t89) - g(3) * (-t42 * t68 + t63) 0, 0, 0, 0, 0, 0, -g(1) * (-t29 * t45 + t30 * t79) - g(2) * (-t27 * t45 + t28 * t79) - (t38 * t69 + t45 * t50) * t83, -g(1) * (-t29 * t48 - t30 * t80) - g(2) * (-t27 * t48 - t28 * t80) - (-t38 * t70 + t48 * t50) * t83, t6, -g(1) * (-t26 + t88) - g(2) * (-t25 + t89) - g(3) * t63 - (t47 * t57 - t68) * t83 + t87 * (qJ(3) + t57); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t65 + t81) - g(2) * t90 - t53, -g(1) * (-t29 * t46 - t41 * t74) - g(2) * (-t27 * t46 + t43 * t74) - g(3) * (-t44 * t49 + t46 * t73) 0, 0, 0, 0, 0, 0, 0, 0, -t54, t5, 0, -g(1) * t19 - g(2) * t67 + (g(1) * t65 - t53) * pkin(4), 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * (-pkin(4) * t65 + t19 + t62) - g(2) * (t61 + t67) - g(3) * (pkin(4) * t55 + t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t5, -g(1) * t62 - g(2) * t61 - g(3) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t45 + t30 * t48) - g(2) * (t14 * t45 + t28 * t48) - g(3) * (-t22 * t45 + t42 * t69) -g(1) * (-t12 * t48 - t30 * t45) - g(2) * (t14 * t48 - t28 * t45) - g(3) * (-t22 * t48 - t42 * t70) 0, 0;];
taug_reg  = t3;
