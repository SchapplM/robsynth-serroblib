% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRR11
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_gravloadJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t119 = cos(qJ(3));
t109 = cos(pkin(7));
t68 = cos(qJ(1));
t105 = sin(pkin(12));
t118 = sin(qJ(1));
t108 = cos(pkin(12));
t110 = cos(pkin(6));
t97 = t110 * t108;
t83 = t118 * t105 - t68 * t97;
t106 = sin(pkin(7));
t107 = sin(pkin(6));
t94 = t107 * t106;
t127 = t83 * t109 + t68 * t94;
t96 = t110 * t105;
t46 = t118 * t108 + t68 * t96;
t66 = sin(qJ(3));
t26 = -t46 * t119 + t127 * t66;
t95 = t107 * t109;
t37 = -t83 * t106 + t68 * t95;
t61 = pkin(13) + qJ(5);
t58 = sin(t61);
t59 = cos(t61);
t10 = t26 * t59 + t37 * t58;
t23 = t127 * t119 + t46 * t66;
t65 = sin(qJ(6));
t67 = cos(qJ(6));
t131 = t10 * t65 + t23 * t67;
t130 = t10 * t67 - t23 * t65;
t126 = t26 * t58 - t37 * t59;
t77 = t68 * t105 + t118 * t97;
t123 = t77 * t106 + t118 * t95;
t122 = t77 * t109 - t118 * t94;
t121 = t106 * t110 + t108 * t95;
t62 = sin(pkin(13));
t117 = t37 * t62;
t116 = t59 * t65;
t115 = t59 * t67;
t63 = cos(pkin(13));
t57 = t63 * pkin(4) + pkin(3);
t64 = -pkin(10) - qJ(4);
t114 = -t23 * t57 + t26 * t64;
t47 = t68 * t108 - t118 * t96;
t27 = t119 * t122 + t47 * t66;
t28 = t47 * t119 - t122 * t66;
t113 = -t27 * t57 - t28 * t64;
t93 = t107 * t105;
t35 = -t119 * t121 + t66 * t93;
t36 = t119 * t93 + t121 * t66;
t112 = -t35 * t57 - t36 * t64;
t98 = t118 * t107;
t111 = t68 * pkin(1) + qJ(2) * t98;
t104 = t68 * t107;
t11 = -t123 * t59 + t28 * t58;
t103 = g(1) * t126 + g(2) * t11;
t102 = -t118 * pkin(1) + qJ(2) * t104;
t101 = -pkin(5) * t59 - pkin(11) * t58;
t100 = -g(1) * t23 + g(2) * t27;
t45 = -t108 * t94 + t110 * t109;
t17 = -t36 * t58 + t45 * t59;
t89 = g(1) * t11 - g(2) * t126 - g(3) * t17;
t12 = t123 * t58 + t28 * t59;
t18 = t36 * t59 + t45 * t58;
t88 = g(1) * t12 - g(2) * t10 + g(3) * t18;
t87 = g(1) * t27 + g(2) * t23 + g(3) * t35;
t86 = g(1) * t28 - g(2) * t26 + g(3) * t36;
t74 = -t46 * pkin(2) + t37 * pkin(9) + t102;
t73 = t47 * pkin(2) + t123 * pkin(9) + t111;
t72 = pkin(4) * t117 + t23 * t64 + t26 * t57 + t74;
t69 = t123 * t62;
t71 = pkin(4) * t69 - t27 * t64 + t28 * t57 + t73;
t42 = -g(1) * t98 + g(2) * t104 - g(3) * t110;
t3 = t12 * t67 + t27 * t65;
t2 = -t12 * t65 + t27 * t67;
t1 = t87 * t58;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t118 - g(2) * t68, g(1) * t68 + g(2) * t118, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t47, -g(1) * t83 + g(2) * t77, -g(1) * t104 - g(2) * t98, -g(1) * t102 - g(2) * t111, 0, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t28, t100, -g(1) * t37 - g(2) * t123, -g(1) * t74 - g(2) * t73, 0, 0, 0, 0, 0, 0, -g(1) * (t26 * t63 + t117) - g(2) * (t28 * t63 + t69) -g(1) * (-t26 * t62 + t37 * t63) - g(2) * (t123 * t63 - t28 * t62) -t100, -g(1) * (t26 * pkin(3) - qJ(4) * t23 + t74) - g(2) * (t28 * pkin(3) + t27 * qJ(4) + t73) 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, t103, -t100, -g(1) * t72 - g(2) * t71, 0, 0, 0, 0, 0, 0, -g(1) * t130 - g(2) * t3, g(1) * t131 - g(2) * t2, -t103, -g(1) * (t10 * pkin(5) + pkin(11) * t126 + t72) - g(2) * (t12 * pkin(5) + t11 * pkin(11) + t71); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, t86, 0, 0, 0, 0, 0, 0, 0, 0, t87 * t63, -t87 * t62, -t86, -g(1) * (-t27 * pkin(3) + t28 * qJ(4)) - g(2) * (-t23 * pkin(3) - qJ(4) * t26) - g(3) * (-t35 * pkin(3) + t36 * qJ(4)) 0, 0, 0, 0, 0, 0, t87 * t59, -t1, -t86, -g(1) * t113 - g(2) * t114 - g(3) * t112, 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t115 + t28 * t65) - g(2) * (-t23 * t115 - t26 * t65) - g(3) * (-t35 * t115 + t36 * t65) -g(1) * (t27 * t116 + t28 * t67) - g(2) * (t23 * t116 - t26 * t67) - g(3) * (t35 * t116 + t36 * t67) t1, -g(1) * (t101 * t27 + t113) - g(2) * (t101 * t23 + t114) - g(3) * (t101 * t35 + t112); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t88, 0, 0, 0, 0, 0, 0, 0, 0, t89 * t67, -t89 * t65, -t88, -g(1) * (-t11 * pkin(5) + t12 * pkin(11)) - g(2) * (pkin(5) * t126 - pkin(11) * t10) - g(3) * (t17 * pkin(5) + t18 * pkin(11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t131 - g(3) * (-t18 * t65 + t35 * t67) g(1) * t3 - g(2) * t130 - g(3) * (-t18 * t67 - t35 * t65) 0, 0;];
taug_reg  = t4;
