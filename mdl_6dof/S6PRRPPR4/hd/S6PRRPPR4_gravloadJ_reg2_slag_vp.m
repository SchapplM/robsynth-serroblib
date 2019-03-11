% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPPR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t61 = sin(pkin(11));
t100 = qJ(5) * t61;
t109 = cos(qJ(3));
t66 = sin(qJ(2));
t68 = cos(qJ(2));
t98 = cos(pkin(10));
t99 = cos(pkin(6));
t77 = t99 * t98;
t97 = sin(pkin(10));
t46 = t66 * t77 + t97 * t68;
t65 = sin(qJ(3));
t62 = sin(pkin(6));
t88 = t62 * t98;
t27 = t109 * t88 + t46 * t65;
t63 = cos(pkin(11));
t108 = t27 * t63;
t116 = -pkin(4) * t108 - t27 * t100;
t76 = t99 * t97;
t48 = -t66 * t76 + t98 * t68;
t87 = t62 * t97;
t29 = -t109 * t87 + t48 * t65;
t107 = t29 * t63;
t115 = -pkin(4) * t107 - t29 * t100;
t105 = t62 * t66;
t49 = t65 * t105 - t99 * t109;
t106 = t49 * t63;
t114 = -pkin(4) * t106 - t49 * t100;
t74 = g(1) * t29 + g(2) * t27 + g(3) * t49;
t113 = pkin(9) * t65;
t104 = t62 * t68;
t103 = -pkin(9) + qJ(4);
t102 = pkin(2) * t104 + pkin(8) * t105;
t101 = qJ(4) * t65;
t96 = t65 * t104;
t45 = t97 * t66 - t68 * t77;
t95 = t45 * t109;
t47 = t98 * t66 + t68 * t76;
t94 = t47 * t109;
t93 = t61 * t109;
t92 = t63 * t109;
t91 = t68 * t109;
t90 = -t45 * pkin(2) + t46 * pkin(8);
t89 = -t47 * pkin(2) + t48 * pkin(8);
t23 = t27 * pkin(3);
t28 = t46 * t109 - t65 * t88;
t86 = t28 * qJ(4) - t23;
t24 = t29 * pkin(3);
t30 = t48 * t109 + t65 * t87;
t85 = t30 * qJ(4) - t24;
t44 = t49 * pkin(3);
t50 = t109 * t105 + t99 * t65;
t84 = t50 * qJ(4) - t44;
t82 = t62 * t91;
t83 = pkin(3) * t82 + qJ(4) * t96 + t102;
t79 = -pkin(3) * t95 - t45 * t101 + t90;
t78 = -pkin(3) * t94 - t47 * t101 + t89;
t15 = -t45 * t93 - t46 * t63;
t17 = -t47 * t93 - t48 * t63;
t32 = -t63 * t105 + t61 * t82;
t75 = g(1) * t17 + g(2) * t15 + g(3) * t32;
t7 = g(1) * t30 + g(2) * t28 + g(3) * t50;
t33 = (t61 * t66 + t63 * t91) * t62;
t73 = t33 * pkin(4) + t32 * qJ(5) + t83;
t72 = -g(1) * t47 - g(2) * t45 + g(3) * t104;
t71 = g(1) * t48 + g(2) * t46 + g(3) * t105;
t16 = -t45 * t92 + t46 * t61;
t70 = t16 * pkin(4) + t15 * qJ(5) + t79;
t18 = -t47 * t92 + t48 * t61;
t69 = t18 * pkin(4) + t17 * qJ(5) + t78;
t67 = cos(qJ(6));
t64 = sin(qJ(6));
t26 = -t61 * t104 + t50 * t63;
t25 = t63 * t104 + t50 * t61;
t12 = t72 * t65;
t11 = t30 * t63 + t47 * t61;
t10 = t30 * t61 - t47 * t63;
t9 = t28 * t63 + t45 * t61;
t8 = t28 * t61 - t45 * t63;
t4 = t74 * t63;
t3 = t74 * t61;
t2 = -g(1) * t18 - g(2) * t16 - g(3) * t33;
t1 = -g(1) * t10 - g(2) * t8 - g(3) * t25;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t71, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t94 + g(2) * t95 - g(3) * t82, t12, -t71, -g(1) * t89 - g(2) * t90 - g(3) * t102, 0, 0, 0, 0, 0, 0, t2, t75, -t12, -g(1) * t78 - g(2) * t79 - g(3) * t83, 0, 0, 0, 0, 0, 0, t2, -t12, -t75, -g(1) * t69 - g(2) * t70 - g(3) * t73, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t64 + t18 * t67) - g(2) * (t15 * t64 + t16 * t67) - g(3) * (t32 * t64 + t33 * t67) -g(1) * (t17 * t67 - t18 * t64) - g(2) * (t15 * t67 - t16 * t64) - g(3) * (t32 * t67 - t33 * t64) t12, -g(1) * (t18 * pkin(5) + t47 * t113 + t69) - g(2) * (t16 * pkin(5) + t45 * t113 + t70) - g(3) * (t33 * pkin(5) - pkin(9) * t96 + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t7, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t7, -g(1) * t85 - g(2) * t86 - g(3) * t84, 0, 0, 0, 0, 0, 0, t4, -t7, t3, -g(1) * (t85 + t115) - g(2) * (t86 + t116) - g(3) * (t84 + t114) 0, 0, 0, 0, 0, 0, t74 * (t61 * t64 + t63 * t67) t74 * (t61 * t67 - t63 * t64) t7, -g(1) * (-pkin(5) * t107 + t103 * t30 + t115 - t24) - g(2) * (-pkin(5) * t108 + t103 * t28 + t116 - t23) - g(3) * (-pkin(5) * t106 + t103 * t50 + t114 - t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t10 * t67 - t11 * t64) - g(2) * (-t9 * t64 + t8 * t67) - g(3) * (t25 * t67 - t26 * t64) -g(1) * (-t10 * t64 - t11 * t67) - g(2) * (-t8 * t64 - t9 * t67) - g(3) * (-t25 * t64 - t26 * t67) 0, 0;];
taug_reg  = t5;
