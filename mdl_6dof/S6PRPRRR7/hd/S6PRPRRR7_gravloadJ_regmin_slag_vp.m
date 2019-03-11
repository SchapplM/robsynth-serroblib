% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_gravloadJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t49 = sin(pkin(8));
t54 = cos(pkin(8));
t48 = sin(pkin(14));
t51 = sin(pkin(6));
t52 = cos(pkin(14));
t59 = sin(qJ(2));
t50 = sin(pkin(7));
t88 = cos(pkin(6));
t82 = t88 * t50;
t55 = cos(pkin(7));
t62 = cos(qJ(2));
t92 = t55 * t62;
t65 = t52 * t82 + (-t48 * t59 + t52 * t92) * t51;
t95 = t51 * t62;
t73 = t50 * t95 - t88 * t55;
t104 = -t73 * t49 + t65 * t54;
t53 = cos(pkin(13));
t87 = sin(pkin(13));
t80 = t88 * t87;
t46 = t53 * t62 - t59 * t80;
t45 = -t53 * t59 - t62 * t80;
t84 = t51 * t87;
t72 = t45 * t55 + t50 * t84;
t66 = t46 * t48 - t72 * t52;
t74 = t45 * t50 - t55 * t84;
t103 = t74 * t49 + t66 * t54;
t83 = t53 * t88;
t44 = t59 * t83 + t87 * t62;
t43 = -t87 * t59 + t62 * t83;
t97 = t51 * t53;
t78 = t43 * t55 - t50 * t97;
t68 = t44 * t48 - t78 * t52;
t79 = t43 * t50 + t55 * t97;
t102 = t79 * t49 + t68 * t54;
t101 = cos(qJ(4));
t100 = t48 * t55;
t99 = t50 * t49;
t98 = t50 * t54;
t96 = t51 * t59;
t94 = t52 * t55;
t93 = t55 * t59;
t56 = sin(qJ(6));
t61 = cos(qJ(5));
t91 = t56 * t61;
t60 = cos(qJ(6));
t90 = t60 * t61;
t89 = qJ(3) * t50;
t86 = t50 * t96;
t85 = t54 * t101;
t81 = t101 * t99;
t29 = t44 * t52 + t78 * t48;
t58 = sin(qJ(4));
t10 = t29 * t101 - t102 * t58;
t30 = t46 * t52 + t72 * t48;
t12 = t30 * t101 - t103 * t58;
t39 = t52 * t96 + (t51 * t92 + t82) * t48;
t19 = t39 * t101 + t104 * t58;
t20 = t68 * t49 - t79 * t54;
t21 = t66 * t49 - t74 * t54;
t28 = -t65 * t49 - t73 * t54;
t57 = sin(qJ(5));
t77 = g(1) * (-t12 * t57 + t21 * t61) + g(2) * (-t10 * t57 + t20 * t61) + g(3) * (-t19 * t57 + t28 * t61);
t11 = t103 * t101 + t30 * t58;
t18 = -t104 * t101 + t39 * t58;
t9 = t102 * t101 + t29 * t58;
t76 = g(1) * t11 + g(2) * t9 + g(3) * t18;
t71 = g(1) * t46 + g(2) * t44 + g(3) * t96;
t42 = (-t48 * t93 + t52 * t62) * t51;
t41 = (-t48 * t62 - t52 * t93) * t51;
t35 = -t41 * t49 + t54 * t86;
t34 = -t46 * t100 + t45 * t52;
t33 = -t45 * t48 - t46 * t94;
t32 = -t44 * t100 + t43 * t52;
t31 = -t43 * t48 - t44 * t94;
t25 = t42 * t101 + (t41 * t54 + t49 * t86) * t58;
t24 = -t41 * t85 + t42 * t58 - t81 * t96;
t23 = -t33 * t49 + t46 * t98;
t22 = -t31 * t49 + t44 * t98;
t17 = t25 * t61 + t35 * t57;
t16 = t34 * t101 + (t33 * t54 + t46 * t99) * t58;
t15 = -t33 * t85 + t34 * t58 - t46 * t81;
t14 = t32 * t101 + (t31 * t54 + t44 * t99) * t58;
t13 = -t31 * t85 + t32 * t58 - t44 * t81;
t8 = t19 * t61 + t28 * t57;
t6 = t16 * t61 + t23 * t57;
t5 = t14 * t61 + t22 * t57;
t4 = t12 * t61 + t21 * t57;
t2 = t10 * t61 + t20 * t57;
t1 = [-g(3), 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(1) * t45 - g(2) * t43 - g(3) * t95, t71, -g(1) * t34 - g(2) * t32 - g(3) * t42, -g(1) * t33 - g(2) * t31 - g(3) * t41, -t71 * t50, -g(1) * (t45 * pkin(2) + t46 * t89) - g(2) * (t43 * pkin(2) + t44 * t89) - g(3) * (pkin(2) * t62 + t59 * t89) * t51, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14 - g(3) * t25, g(1) * t15 + g(2) * t13 + g(3) * t24, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t5 - g(3) * t17, -g(1) * (-t16 * t57 + t23 * t61) - g(2) * (-t14 * t57 + t22 * t61) - g(3) * (-t25 * t57 + t35 * t61) 0, 0, 0, 0, 0, -g(1) * (t15 * t56 + t6 * t60) - g(2) * (t13 * t56 + t5 * t60) - g(3) * (t17 * t60 + t24 * t56) -g(1) * (t15 * t60 - t6 * t56) - g(2) * (t13 * t60 - t5 * t56) - g(3) * (-t17 * t56 + t24 * t60); 0, 0, 0, 0, 0, 0, 0, g(1) * t74 + g(2) * t79 + g(3) * t73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, g(1) * t12 + g(2) * t10 + g(3) * t19, 0, 0, 0, 0, 0, t76 * t61, -t76 * t57, 0, 0, 0, 0, 0, -g(1) * (-t11 * t90 + t12 * t56) - g(2) * (t10 * t56 - t9 * t90) - g(3) * (-t18 * t90 + t19 * t56) -g(1) * (t11 * t91 + t12 * t60) - g(2) * (t10 * t60 + t9 * t91) - g(3) * (t18 * t91 + t19 * t60); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, g(1) * t4 + g(2) * t2 + g(3) * t8, 0, 0, 0, 0, 0, -t77 * t60, t77 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t11 * t60 - t4 * t56) - g(2) * (-t2 * t56 + t9 * t60) - g(3) * (t18 * t60 - t8 * t56) -g(1) * (-t11 * t56 - t4 * t60) - g(2) * (-t2 * t60 - t9 * t56) - g(3) * (-t18 * t56 - t8 * t60);];
taug_reg  = t1;
