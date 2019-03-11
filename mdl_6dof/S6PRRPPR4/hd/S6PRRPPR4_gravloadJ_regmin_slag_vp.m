% Calculate minimal parameter regressor of gravitation load for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t54 = sin(qJ(3));
t91 = qJ(4) * t54 + pkin(2);
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t82 = cos(pkin(10));
t83 = cos(pkin(6));
t63 = t83 * t82;
t81 = sin(pkin(10));
t35 = t55 * t63 + t81 * t57;
t51 = sin(pkin(6));
t75 = t51 * t82;
t87 = cos(qJ(3));
t19 = t35 * t54 + t87 * t75;
t62 = t83 * t81;
t37 = -t55 * t62 + t82 * t57;
t74 = t51 * t81;
t21 = t37 * t54 - t87 * t74;
t86 = t51 * t55;
t38 = t54 * t86 - t83 * t87;
t60 = g(1) * t21 + g(2) * t19 + g(3) * t38;
t85 = t51 * t57;
t34 = t81 * t55 - t57 * t63;
t80 = t34 * t87;
t36 = t82 * t55 + t57 * t62;
t79 = t36 * t87;
t50 = sin(pkin(11));
t78 = t50 * t87;
t52 = cos(pkin(11));
t77 = t52 * t87;
t76 = t57 * t87;
t20 = t35 * t87 - t54 * t75;
t73 = -t19 * pkin(3) + t20 * qJ(4);
t22 = t37 * t87 + t54 * t74;
t72 = -t21 * pkin(3) + t22 * qJ(4);
t39 = t83 * t54 + t87 * t86;
t71 = -t38 * pkin(3) + t39 * qJ(4);
t69 = t51 * t76;
t70 = pkin(3) * t69 + pkin(8) * t86 + t91 * t85;
t68 = -pkin(4) * t52 - qJ(5) * t50;
t65 = -pkin(3) * t80 + t35 * pkin(8) - t91 * t34;
t64 = -pkin(3) * t79 + t37 * pkin(8) - t91 * t36;
t11 = -t34 * t78 - t35 * t52;
t13 = -t36 * t78 - t37 * t52;
t23 = t50 * t69 - t52 * t86;
t61 = g(1) * t13 + g(2) * t11 + g(3) * t23;
t59 = g(1) * t22 + g(2) * t20 + g(3) * t39;
t58 = -g(1) * t36 - g(2) * t34 + g(3) * t85;
t56 = cos(qJ(6));
t53 = sin(qJ(6));
t24 = (t50 * t55 + t52 * t76) * t51;
t18 = t39 * t52 - t50 * t85;
t17 = t39 * t50 + t52 * t85;
t14 = -t36 * t77 + t37 * t50;
t12 = -t34 * t77 + t35 * t50;
t10 = t58 * t54;
t9 = t22 * t52 + t36 * t50;
t8 = t22 * t50 - t36 * t52;
t7 = t20 * t52 + t34 * t50;
t6 = t20 * t50 - t34 * t52;
t3 = t60 * t52;
t2 = t60 * t50;
t1 = -g(1) * t14 - g(2) * t12 - g(3) * t24;
t4 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -t58, g(1) * t37 + g(2) * t35 + g(3) * t86, 0, 0, 0, 0, 0, g(1) * t79 + g(2) * t80 - g(3) * t69, t10, t1, t61, -t10, -g(1) * t64 - g(2) * t65 - g(3) * t70, t1, -t10, -t61, -g(1) * (t14 * pkin(4) + t13 * qJ(5) + t64) - g(2) * (t12 * pkin(4) + t11 * qJ(5) + t65) - g(3) * (t24 * pkin(4) + t23 * qJ(5) + t70) 0, 0, 0, 0, 0, -g(1) * (t13 * t53 + t14 * t56) - g(2) * (t11 * t53 + t12 * t56) - g(3) * (t23 * t53 + t24 * t56) -g(1) * (t13 * t56 - t14 * t53) - g(2) * (t11 * t56 - t12 * t53) - g(3) * (t23 * t56 - t24 * t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t59, t3, -t2, -t59, -g(1) * t72 - g(2) * t73 - g(3) * t71, t3, -t59, t2, -g(1) * (t68 * t21 + t72) - g(2) * (t68 * t19 + t73) - g(3) * (t68 * t38 + t71) 0, 0, 0, 0, 0, t60 * (t50 * t53 + t52 * t56) t60 * (t50 * t56 - t52 * t53); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, 0, 0, 0, -t60, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6 - g(3) * t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t9 * t53 + t8 * t56) - g(2) * (-t7 * t53 + t6 * t56) - g(3) * (t17 * t56 - t18 * t53) -g(1) * (-t8 * t53 - t9 * t56) - g(2) * (-t6 * t53 - t7 * t56) - g(3) * (-t17 * t53 - t18 * t56);];
taug_reg  = t4;
