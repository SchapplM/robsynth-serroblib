% Calculate inertial parameters regressor of gravitation load for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRRRP1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t52 = sin(qJ(3));
t55 = cos(qJ(3));
t76 = cos(pkin(6));
t49 = sin(pkin(6));
t53 = sin(qJ(2));
t87 = t49 * t53;
t97 = -t52 * t87 + t76 * t55;
t56 = cos(qJ(2));
t48 = sin(pkin(11));
t69 = t48 * t76;
t75 = cos(pkin(11));
t34 = -t53 * t69 + t75 * t56;
t86 = t49 * t55;
t96 = -t34 * t52 + t48 * t86;
t44 = t55 * pkin(3) + pkin(2);
t85 = t49 * t56;
t35 = t44 * t85;
t95 = g(3) * t35;
t94 = g(3) * t49;
t62 = t76 * t75;
t32 = t48 * t56 + t53 * t62;
t51 = sin(qJ(5));
t93 = t32 * t51;
t92 = t34 * t51;
t47 = qJ(3) + qJ(4);
t46 = cos(t47);
t90 = t46 * t51;
t54 = cos(qJ(5));
t89 = t46 * t54;
t88 = t48 * t49;
t84 = t51 * t56;
t57 = -pkin(9) - pkin(8);
t83 = t53 * t57;
t82 = t54 * t56;
t45 = sin(t47);
t68 = t49 * t75;
t18 = t32 * t45 + t46 * t68;
t19 = t32 * t46 - t45 * t68;
t43 = t54 * pkin(5) + pkin(4);
t50 = -qJ(6) - pkin(10);
t81 = -t18 * t43 - t19 * t50;
t20 = t34 * t45 - t46 * t88;
t21 = t34 * t46 + t45 * t88;
t80 = -t20 * t43 - t21 * t50;
t27 = t45 * t87 - t76 * t46;
t28 = t76 * t45 + t46 * t87;
t79 = -t27 * t43 - t28 * t50;
t31 = t48 * t53 - t56 * t62;
t78 = -t31 * t44 - t32 * t57;
t33 = t75 * t53 + t56 * t69;
t77 = -t33 * t44 - t34 * t57;
t72 = -t18 * pkin(4) + t19 * pkin(10);
t71 = -t20 * pkin(4) + t21 * pkin(10);
t70 = -t27 * pkin(4) + t28 * pkin(10);
t66 = t96 * pkin(3);
t65 = pkin(4) * t46 + pkin(10) * t45;
t64 = t43 * t46 - t45 * t50;
t63 = t97 * pkin(3);
t7 = g(1) * t20 + g(2) * t18 + g(3) * t27;
t9 = g(1) * t21 + g(2) * t19 + g(3) * t28;
t61 = -t32 * t52 - t55 * t68;
t60 = -g(1) * t33 - g(2) * t31 + g(3) * t85;
t59 = g(1) * t34 + g(2) * t32 + g(3) * t87;
t58 = t61 * pkin(3);
t1 = -g(1) * (-t21 * t51 + t33 * t54) - g(2) * (-t19 * t51 + t31 * t54) - g(3) * (-t28 * t51 - t49 * t82);
t10 = t60 * t45;
t6 = t7 * t54;
t5 = t7 * t51;
t4 = -g(1) * (-t33 * t89 + t92) - g(2) * (-t31 * t89 + t93) - (t46 * t82 + t51 * t53) * t94;
t3 = -g(1) * (t33 * t90 + t34 * t54) - g(2) * (t31 * t90 + t32 * t54) - (-t46 * t84 + t53 * t54) * t94;
t2 = -g(1) * (-t21 * t54 - t33 * t51) - g(2) * (-t19 * t54 - t31 * t51) - g(3) * (-t28 * t54 + t49 * t84);
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t59, 0, 0, 0, 0, 0, 0, 0, 0, -t60 * t55, t60 * t52, -t59, -g(1) * (-t33 * pkin(2) + t34 * pkin(8)) - g(2) * (-t31 * pkin(2) + t32 * pkin(8)) - (pkin(2) * t56 + pkin(8) * t53) * t94, 0, 0, 0, 0, 0, 0, -t60 * t46, t10, -t59, -g(1) * t77 - g(2) * t78 - g(3) * (-t49 * t83 + t35) 0, 0, 0, 0, 0, 0, t4, t3, -t10, -g(1) * (-t65 * t33 + t77) - g(2) * (-t65 * t31 + t78) - t95 - (t65 * t56 - t83) * t94, 0, 0, 0, 0, 0, 0, t4, t3, -t10, -g(1) * (pkin(5) * t92 - t64 * t33 + t77) - g(2) * (pkin(5) * t93 - t64 * t31 + t78) - t95 - (t64 * t56 + (pkin(5) * t51 - t57) * t53) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t61 - g(3) * t97, -g(1) * (-t34 * t55 - t52 * t88) - g(2) * (-t32 * t55 + t52 * t68) - g(3) * (-t76 * t52 - t53 * t86) 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, -g(1) * t66 - g(2) * t58 - g(3) * t63, 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * (t66 + t71) - g(2) * (t58 + t72) - g(3) * (t63 + t70) 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * (t66 + t80) - g(2) * (t58 + t81) - g(3) * (t63 + t79); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t9, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * t71 - g(2) * t72 - g(3) * t70, 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * t80 - g(2) * t81 - g(3) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t8;
