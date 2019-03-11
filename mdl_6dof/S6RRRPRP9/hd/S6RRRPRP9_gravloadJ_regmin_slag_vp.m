% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t59 = sin(qJ(3));
t63 = cos(qJ(3));
t65 = cos(qJ(1));
t90 = t65 * t63;
t61 = sin(qJ(1));
t64 = cos(qJ(2));
t94 = t61 * t64;
t36 = t59 * t94 + t90;
t91 = t65 * t59;
t93 = t63 * t64;
t37 = t61 * t93 - t91;
t58 = sin(qJ(5));
t62 = cos(qJ(5));
t107 = -t36 * t62 + t37 * t58;
t60 = sin(qJ(2));
t96 = t60 * t63;
t97 = t59 * t62;
t28 = t58 * t96 - t60 * t97;
t38 = -t61 * t63 + t64 * t91;
t39 = t59 * t61 + t64 * t90;
t77 = -t38 * t62 + t39 * t58;
t2 = g(1) * t77 + g(2) * t107 + g(3) * t28;
t16 = t38 * t58 + t39 * t62;
t112 = -pkin(5) * t77 + qJ(6) * t16;
t10 = t36 * t58 + t37 * t62;
t111 = -pkin(5) * t107 + qJ(6) * t10;
t76 = t58 * t59 + t62 * t63;
t108 = t76 * t60;
t71 = g(1) * t16 + g(2) * t10 + g(3) * t108;
t80 = g(1) * t65 + g(2) * t61;
t109 = t60 * t80;
t105 = -pkin(3) - pkin(4);
t104 = g(1) * t61;
t53 = t60 * pkin(8);
t55 = t64 * pkin(2);
t98 = t59 * t60;
t95 = t60 * t65;
t92 = t64 * t65;
t89 = qJ(4) * t59;
t88 = -pkin(1) - t55;
t87 = -pkin(2) - t89;
t86 = -pkin(3) * t36 + qJ(4) * t37;
t85 = -pkin(3) * t38 + qJ(4) * t39;
t84 = pkin(3) * t93 + t64 * t89 + t53 + t55;
t83 = -pkin(3) * t37 + pkin(7) * t65 - t36 * qJ(4);
t82 = -g(1) * t107 + g(2) * t77;
t81 = g(1) * t36 - g(2) * t38;
t79 = -g(2) * t65 + t104;
t78 = -t28 * pkin(5) + qJ(6) * t108;
t74 = t60 * (t58 * t63 - t97);
t20 = t61 * t74;
t22 = t65 * t74;
t30 = t58 * t93 - t64 * t97;
t70 = g(1) * t22 + g(2) * t20 - g(3) * t30;
t68 = pkin(1) * t65 + pkin(2) * t92 + pkin(3) * t39 + pkin(7) * t61 + pkin(8) * t95 + qJ(4) * t38;
t6 = g(1) * t38 + g(2) * t36 + g(3) * t98;
t67 = g(1) * t39 + g(2) * t37 + g(3) * t96;
t66 = -g(3) * t64 + t109;
t48 = pkin(8) * t92;
t45 = pkin(8) * t94;
t43 = qJ(4) * t96;
t40 = -g(2) * t95 + t104 * t60;
t31 = t76 * t64;
t24 = g(3) * t60 + t64 * t80;
t23 = t65 * t108;
t21 = t61 * t108;
t19 = t66 * t63;
t18 = t66 * t59;
t17 = g(1) * t37 - g(2) * t39;
t4 = g(1) * t23 + g(2) * t21 - g(3) * t31;
t3 = g(1) * t10 - g(2) * t16;
t1 = [0, t79, t80, 0, 0, 0, 0, 0, t79 * t64, -t40, 0, 0, 0, 0, 0, t17, -t81, t17, t40, t81, -g(1) * t83 - g(2) * t68 - (t88 - t53) * t104, 0, 0, 0, 0, 0, t3, t82, t3, -t40, -t82, -g(1) * (-t37 * pkin(4) - pkin(5) * t10 - qJ(6) * t107 + t83) - g(2) * (pkin(4) * t39 + pkin(5) * t16 - pkin(9) * t95 + qJ(6) * t77 + t68) - ((-pkin(8) + pkin(9)) * t60 + t88) * t104; 0, 0, 0, 0, 0, 0, 0, 0, t66, t24, 0, 0, 0, 0, 0, t19, -t18, t19, -t24, t18, -g(1) * t48 - g(2) * t45 - g(3) * t84 + (pkin(3) * t63 - t87) * t109, 0, 0, 0, 0, 0, t4, -t70, t4, t24, t70, -g(1) * (-t23 * pkin(5) - pkin(9) * t92 - t22 * qJ(6) + t48) - g(2) * (-t21 * pkin(5) - pkin(9) * t94 - t20 * qJ(6) + t45) - g(3) * (pkin(4) * t93 + pkin(5) * t31 + qJ(6) * t30 + t84) + (g(3) * pkin(9) + t80 * (-t105 * t63 - t87)) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t67, t6, 0, -t67, -g(1) * t85 - g(2) * t86 - g(3) * (-pkin(3) * t98 + t43) 0, 0, 0, 0, 0, -t2, -t71, -t2, 0, t71, -g(1) * (-pkin(4) * t38 - t112 + t85) - g(2) * (-pkin(4) * t36 - t111 + t86) - g(3) * (t105 * t98 + t43 - t78); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t71, t2, 0, -t71, -g(1) * t112 - g(2) * t111 - g(3) * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg  = t1;
