% Calculate minimal parameter regressor of gravitation load for
% S6RRRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRPR13_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR13_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR13_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR13_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 01:35:37
% EndTime: 2019-05-08 01:35:41
% DurationCPUTime: 0.74s
% Computational Cost: add. (568->121), mult. (1521->192), div. (0->0), fcn. (1959->12), ass. (0->65)
t48 = sin(qJ(6));
t52 = cos(qJ(6));
t51 = sin(qJ(2));
t55 = cos(qJ(2));
t77 = cos(pkin(6));
t86 = cos(qJ(1));
t68 = t77 * t86;
t85 = sin(qJ(1));
t37 = t51 * t68 + t85 * t55;
t50 = sin(qJ(3));
t54 = cos(qJ(3));
t47 = sin(pkin(6));
t75 = t47 * t86;
t22 = t37 * t54 - t50 * t75;
t36 = t85 * t51 - t55 * t68;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t8 = t22 * t49 - t36 * t53;
t9 = t22 * t53 + t36 * t49;
t99 = t9 * t48 - t8 * t52;
t98 = t8 * t48 + t9 * t52;
t67 = t77 * t85;
t38 = t86 * t51 + t55 * t67;
t97 = g(1) * t38 + g(2) * t36;
t39 = -t51 * t67 + t86 * t55;
t74 = t47 * t85;
t25 = t39 * t50 - t54 * t74;
t73 = -t37 * t50 - t54 * t75;
t82 = t47 * t51;
t96 = g(3) * (-t50 * t82 + t77 * t54) + g(2) * t73 - g(1) * t25;
t35 = t77 * t50 + t54 * t82;
t78 = t53 * t55;
t19 = t35 * t49 + t47 * t78;
t81 = t47 * t55;
t76 = t49 * t81;
t20 = t35 * t53 - t76;
t26 = t39 * t54 + t50 * t74;
t12 = t26 * t49 - t38 * t53;
t13 = t26 * t53 + t38 * t49;
t3 = t12 * t48 + t13 * t52;
t95 = g(2) * t98 + g(3) * (t19 * t48 + t20 * t52) + g(1) * t3;
t2 = t12 * t52 - t13 * t48;
t94 = g(2) * t99 - g(3) * (t19 * t52 - t20 * t48) - g(1) * t2;
t80 = t49 * t54;
t79 = t53 * t54;
t72 = -g(1) * t8 + g(2) * t12;
t71 = g(1) * t73 + g(2) * t25;
t62 = pkin(3) * t54 + pkin(10) * t50 + pkin(2);
t1 = g(1) * t12 + g(2) * t8 + g(3) * t19;
t60 = g(1) * t13 + g(2) * t9 + g(3) * t20;
t15 = -t36 * t80 - t37 * t53;
t17 = -t38 * t80 - t39 * t53;
t27 = -t53 * t82 + t54 * t76;
t59 = g(1) * t17 + g(2) * t15 + g(3) * t27;
t57 = g(1) * t26 + g(2) * t22 + g(3) * t35;
t56 = g(3) * t81 - t97;
t28 = (t49 * t51 + t54 * t78) * t47;
t18 = -t38 * t79 + t39 * t49;
t16 = -t36 * t79 + t37 * t49;
t14 = t56 * t50;
t7 = t96 * t53;
t6 = t96 * t49;
t5 = g(1) * t9 - g(2) * t13;
t4 = -g(1) * t18 - g(2) * t16 - g(3) * t28;
t10 = [0, g(1) * t85 - g(2) * t86, g(1) * t86 + g(2) * t85, 0, 0, 0, 0, 0, g(1) * t37 - g(2) * t39, -g(1) * t36 + g(2) * t38, 0, 0, 0, 0, 0, g(1) * t22 - g(2) * t26, t71, 0, 0, 0, 0, 0, t5, t72, t5, -t71, -t72, -g(1) * (-t85 * pkin(1) - t37 * pkin(2) - pkin(3) * t22 - pkin(4) * t9 + pkin(8) * t75 - t36 * pkin(9) + pkin(10) * t73 - qJ(5) * t8) - g(2) * (t86 * pkin(1) + t39 * pkin(2) + t26 * pkin(3) + t13 * pkin(4) + pkin(8) * t74 + t38 * pkin(9) + t25 * pkin(10) + t12 * qJ(5)) 0, 0, 0, 0, 0, g(1) * t98 - g(2) * t3, -g(1) * t99 - g(2) * t2; 0, 0, 0, 0, 0, 0, 0, 0, -t56, g(1) * t39 + g(2) * t37 + g(3) * t82, 0, 0, 0, 0, 0, -t56 * t54, t14, 0, 0, 0, 0, 0, t4, t59, t4, -t14, -t59, -g(1) * (t18 * pkin(4) + t39 * pkin(9) + t17 * qJ(5)) - g(2) * (t16 * pkin(4) + t37 * pkin(9) + t15 * qJ(5)) + t97 * t62 + (-t28 * pkin(4) - t27 * qJ(5) - (pkin(9) * t51 + t62 * t55) * t47) * g(3), 0, 0, 0, 0, 0, -g(1) * (t17 * t48 + t18 * t52) - g(2) * (t15 * t48 + t16 * t52) - g(3) * (t27 * t48 + t28 * t52) -g(1) * (t17 * t52 - t18 * t48) - g(2) * (t15 * t52 - t16 * t48) - g(3) * (t27 * t52 - t28 * t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t57, 0, 0, 0, 0, 0, -t7, t6, -t7, -t57, -t6, -t57 * pkin(10) - t96 * (pkin(4) * t53 + qJ(5) * t49 + pkin(3)) 0, 0, 0, 0, 0, -t96 * (t48 * t49 + t52 * t53) t96 * (t48 * t53 - t49 * t52); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t60, t1, 0, -t60, -g(1) * (-t12 * pkin(4) + t13 * qJ(5)) - g(2) * (-t8 * pkin(4) + t9 * qJ(5)) - g(3) * (-t19 * pkin(4) + t20 * qJ(5)) 0, 0, 0, 0, 0, -t94, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t95;];
taug_reg  = t10;
