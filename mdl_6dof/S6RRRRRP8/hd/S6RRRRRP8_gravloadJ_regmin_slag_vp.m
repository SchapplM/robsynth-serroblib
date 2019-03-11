% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t67 = sin(qJ(2));
t71 = cos(qJ(2));
t72 = cos(qJ(1));
t68 = sin(qJ(1));
t89 = cos(pkin(6));
t84 = t68 * t89;
t51 = -t67 * t84 + t72 * t71;
t63 = qJ(3) + qJ(4);
t61 = sin(t63);
t62 = cos(t63);
t64 = sin(pkin(6));
t95 = t64 * t68;
t32 = t51 * t61 - t62 * t95;
t83 = t72 * t89;
t49 = t67 * t83 + t68 * t71;
t92 = t64 * t72;
t86 = -t49 * t61 - t62 * t92;
t96 = t64 * t67;
t6 = -g(3) * (-t61 * t96 + t89 * t62) - g(2) * t86 + g(1) * t32;
t29 = t49 * t62 - t61 * t92;
t48 = t68 * t67 - t71 * t83;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t10 = t29 * t65 - t48 * t69;
t11 = t29 * t69 + t48 * t65;
t33 = t51 * t62 + t61 * t95;
t42 = t89 * t61 + t62 * t96;
t8 = g(1) * t33 + g(2) * t29 + g(3) * t42;
t66 = sin(qJ(3));
t70 = cos(qJ(3));
t94 = t64 * t70;
t34 = -t51 * t66 + t68 * t94;
t79 = t49 * t66 + t70 * t92;
t109 = g(2) * t79 - g(3) * (-t66 * t96 + t89 * t70) - g(1) * t34;
t108 = g(1) * t72 + g(2) * t68;
t50 = t72 * t67 + t71 * t84;
t107 = g(1) * t50 + g(2) * t48;
t98 = t62 * t65;
t97 = t62 * t69;
t93 = t64 * t71;
t91 = t69 * t71;
t88 = t65 * t93;
t85 = t49 * t70 - t66 * t92;
t14 = t33 * t65 - t50 * t69;
t82 = -g(1) * t10 + g(2) * t14;
t81 = g(1) * t86 + g(2) * t32;
t60 = t70 * pkin(3) + pkin(2);
t80 = pkin(4) * t62 + pkin(11) * t61 + t60;
t26 = t42 * t65 + t64 * t91;
t1 = g(1) * t14 + g(2) * t10 + g(3) * t26;
t15 = t33 * t69 + t50 * t65;
t27 = t42 * t69 - t88;
t78 = g(1) * t15 + g(2) * t11 + g(3) * t27;
t16 = -t48 * t98 - t49 * t69;
t18 = -t50 * t98 - t51 * t69;
t38 = t62 * t88 - t69 * t96;
t77 = g(1) * t18 + g(2) * t16 + g(3) * t38;
t76 = g(3) * t93 - t107;
t74 = -t8 * pkin(11) + t6 * (pkin(5) * t69 + qJ(6) * t65 + pkin(4));
t73 = -pkin(10) - pkin(9);
t39 = (t62 * t91 + t65 * t67) * t64;
t35 = t51 * t70 + t66 * t95;
t19 = -t50 * t97 + t51 * t65;
t17 = -t48 * t97 + t49 * t65;
t9 = t76 * t61;
t5 = t6 * t69;
t4 = t6 * t65;
t3 = g(1) * t11 - g(2) * t15;
t2 = -g(1) * t19 - g(2) * t17 - g(3) * t39;
t7 = [0, g(1) * t68 - g(2) * t72, t108, 0, 0, 0, 0, 0, g(1) * t49 - g(2) * t51, -g(1) * t48 + g(2) * t50, 0, 0, 0, 0, 0, g(1) * t85 - g(2) * t35, -g(1) * t79 - g(2) * t34, 0, 0, 0, 0, 0, g(1) * t29 - g(2) * t33, t81, 0, 0, 0, 0, 0, t3, t82, t3, -t81, -t82, -g(1) * (-t68 * pkin(1) - pkin(4) * t29 - pkin(5) * t11 + pkin(11) * t86 - qJ(6) * t10 + t48 * t73 - t49 * t60) - g(2) * (t72 * pkin(1) + t33 * pkin(4) + t15 * pkin(5) + t32 * pkin(11) + t14 * qJ(6) - t50 * t73 + t51 * t60) - t108 * t64 * (pkin(3) * t66 + pkin(8)); 0, 0, 0, 0, 0, 0, 0, 0, -t76, g(1) * t51 + g(2) * t49 + g(3) * t96, 0, 0, 0, 0, 0, -t76 * t70, t76 * t66, 0, 0, 0, 0, 0, -t76 * t62, t9, 0, 0, 0, 0, 0, t2, t77, t2, -t9, -t77, -g(1) * (t19 * pkin(5) + t18 * qJ(6) - t51 * t73) - g(2) * (t17 * pkin(5) + t16 * qJ(6) - t49 * t73) + t107 * t80 + (-t39 * pkin(5) - t38 * qJ(6) - (-t67 * t73 + t71 * t80) * t64) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, g(1) * t35 + g(2) * t85 - g(3) * (-t89 * t66 - t67 * t94) 0, 0, 0, 0, 0, t6, t8, 0, 0, 0, 0, 0, t5, -t4, t5, -t8, t4, t109 * pkin(3) + t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t8, 0, 0, 0, 0, 0, t5, -t4, t5, -t8, t4, t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t78, t1, 0, -t78, -g(1) * (-t14 * pkin(5) + t15 * qJ(6)) - g(2) * (-t10 * pkin(5) + t11 * qJ(6)) - g(3) * (-t26 * pkin(5) + t27 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t7;
