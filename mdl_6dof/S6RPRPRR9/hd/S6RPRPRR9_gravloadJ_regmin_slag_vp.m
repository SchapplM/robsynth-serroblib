% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_gravloadJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t50 = sin(pkin(7));
t48 = sin(pkin(13));
t52 = cos(pkin(13));
t58 = sin(qJ(3));
t62 = cos(qJ(3));
t74 = t62 * t48 + t58 * t52;
t29 = t74 * t50;
t54 = cos(pkin(7));
t31 = t74 * t54;
t55 = cos(pkin(6));
t53 = cos(pkin(12));
t63 = cos(qJ(1));
t82 = t63 * t53;
t49 = sin(pkin(12));
t59 = sin(qJ(1));
t86 = t59 * t49;
t35 = -t55 * t82 + t86;
t83 = t63 * t49;
t85 = t59 * t53;
t36 = t55 * t83 + t85;
t39 = t58 * t48 - t62 * t52;
t51 = sin(pkin(6));
t89 = t51 * t63;
t14 = t29 * t89 + t35 * t31 + t36 * t39;
t23 = -t35 * t50 + t54 * t89;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t5 = t14 * t61 + t23 * t57;
t56 = sin(qJ(6));
t60 = cos(qJ(6));
t28 = t39 * t50;
t30 = t39 * t54;
t70 = t28 * t89 + t35 * t30 - t36 * t74;
t97 = t5 * t56 - t60 * t70;
t96 = t5 * t60 + t56 * t70;
t95 = t14 * t57 - t23 * t61;
t92 = t50 * t58;
t91 = t51 * t53;
t90 = t51 * t59;
t88 = t54 * t58;
t87 = t56 * t61;
t84 = t60 * t61;
t81 = pkin(9) + qJ(4);
t79 = qJ(2) * t51;
t80 = t63 * pkin(1) + t59 * t79;
t78 = t50 * t89;
t77 = -t59 * pkin(1) + t63 * t79;
t76 = g(1) * t63 + g(2) * t59;
t75 = g(1) * t59 - g(2) * t63;
t73 = g(3) * t49 * t51 + g(2) * t36;
t72 = t35 * t54 + t78;
t37 = -t55 * t85 - t83;
t71 = t37 * t54 + t50 * t90;
t34 = -t50 * t91 + t55 * t54;
t25 = -t37 * t50 + t54 * t90;
t38 = -t55 * t86 + t82;
t66 = t29 * t90 + t37 * t31 - t38 * t39;
t6 = t25 * t61 - t57 * t66;
t65 = t55 * t29 + (t31 * t53 - t39 * t49) * t51;
t69 = g(1) * t6 + g(2) * t95 + g(3) * (t34 * t61 - t57 * t65);
t68 = t35 * t88 - t36 * t62 + t58 * t78;
t16 = -t28 * t90 - t37 * t30 - t38 * t74;
t19 = -t55 * t28 + (-t30 * t53 - t49 * t74) * t51;
t67 = g(1) * t16 + g(2) * t70 + g(3) * t19;
t64 = g(2) * t72 - g(3) * (t50 * t55 + t54 * t91);
t46 = t62 * pkin(3) + pkin(2);
t33 = pkin(3) * t88 - t50 * t81;
t32 = pkin(3) * t92 + t54 * t81;
t27 = -g(3) * t55 - t75 * t51;
t22 = t38 * t62 + t58 * t71;
t21 = -t38 * t58 + t62 * t71;
t9 = t34 * t57 + t61 * t65;
t7 = t25 * t57 + t61 * t66;
t2 = -t16 * t56 + t7 * t60;
t1 = -t16 * t60 - t7 * t56;
t3 = [0, t75, t76, g(1) * t36 - g(2) * t38, -g(1) * t35 - g(2) * t37, -t76 * t51, -g(1) * t77 - g(2) * t80, 0, 0, 0, 0, 0, -g(1) * t68 - g(2) * t22, -g(1) * (t36 * t58 + t62 * t72) - g(2) * t21, -g(1) * t23 - g(2) * t25, -g(1) * (t32 * t89 + t35 * t33 - t36 * t46 + t77) - g(2) * (t32 * t90 + t37 * t33 + t38 * t46 + t80) 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t7, g(1) * t95 - g(2) * t6, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t2, g(1) * t97 - g(2) * t1; 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t21 + t73 * t58 + t64 * t62, g(1) * t22 - g(2) * t68 - g(3) * (-t55 * t92 + (-t49 * t62 - t53 * t88) * t51) 0 ((g(1) * t38 + t73) * t58 + (-g(1) * t71 + t64) * t62) * pkin(3), 0, 0, 0, 0, 0, -t67 * t61, t67 * t57, 0, 0, 0, 0, 0, -g(1) * (t16 * t84 + t56 * t66) - g(2) * (-t14 * t56 + t70 * t84) - g(3) * (t19 * t84 + t56 * t65) -g(1) * (-t16 * t87 + t60 * t66) - g(2) * (-t14 * t60 - t70 * t87) - g(3) * (-t19 * t87 + t60 * t65); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t25 + g(2) * t23 - g(3) * t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, g(1) * t7 - g(2) * t5 + g(3) * t9, 0, 0, 0, 0, 0, -t69 * t60, t69 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t97 - g(3) * (-t19 * t60 - t9 * t56) g(1) * t2 - g(2) * t96 - g(3) * (t19 * t56 - t9 * t60);];
taug_reg  = t3;
