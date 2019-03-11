% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t122 = sin(pkin(7));
t84 = sin(pkin(6));
t117 = t84 * t122;
t90 = cos(qJ(1));
t113 = t90 * t117;
t123 = cos(pkin(7));
t87 = sin(qJ(3));
t115 = t87 * t123;
t131 = cos(qJ(3));
t124 = cos(pkin(6));
t132 = cos(qJ(2));
t110 = t124 * t132;
t129 = sin(qJ(2));
t130 = sin(qJ(1));
t75 = -t110 * t90 + t129 * t130;
t109 = t124 * t129;
t76 = t109 * t90 + t130 * t132;
t46 = -t87 * t113 - t115 * t75 + t131 * t76;
t118 = t84 * t123;
t65 = -t118 * t90 + t122 * t75;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t25 = t46 * t89 + t65 * t86;
t108 = t123 * t131;
t45 = t108 * t75 + t113 * t131 + t76 * t87;
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t8 = t25 * t85 - t45 * t88;
t9 = t25 * t88 + t45 * t85;
t26 = -t46 * t86 + t65 * t89;
t93 = t110 * t130 + t129 * t90;
t140 = t130 * t117 - t93 * t123;
t77 = -t109 * t130 + t132 * t90;
t49 = -t140 * t131 + t77 * t87;
t105 = t124 * t122;
t120 = t84 * t129;
t121 = t84 * t132;
t63 = -t105 * t131 - t108 * t121 + t120 * t87;
t94 = g(1) * t49 + g(2) * t45 + g(3) * t63;
t50 = t77 * t131 + t140 * t87;
t91 = -t118 * t130 - t122 * t93;
t28 = t50 * t86 + t89 * t91;
t64 = t87 * t105 + (t115 * t132 + t129 * t131) * t84;
t74 = -t117 * t132 + t123 * t124;
t97 = -g(3) * (-t64 * t86 + t74 * t89) - g(2) * t26 + g(1) * t28;
t139 = pkin(9) * t84;
t126 = t85 * t89;
t125 = t88 * t89;
t119 = pkin(10) * t122;
t116 = t86 * t122;
t114 = t89 * t122;
t29 = t50 * t89 - t86 * t91;
t12 = t29 * t85 - t49 * t88;
t112 = -g(1) * t8 + g(2) * t12;
t111 = g(1) * t26 + g(2) * t28;
t107 = t123 * t129;
t106 = t122 * t129;
t102 = t84 * t106;
t44 = t64 * t89 + t74 * t86;
t22 = t44 * t85 - t63 * t88;
t1 = g(1) * t12 + g(2) * t8 + g(3) * t22;
t13 = t29 * t88 + t49 * t85;
t23 = t44 * t88 + t63 * t85;
t100 = g(1) * t13 + g(2) * t9 + g(3) * t23;
t14 = -t126 * t45 - t46 * t88;
t16 = -t126 * t49 - t50 * t88;
t30 = -t126 * t63 - t64 * t88;
t99 = g(1) * t16 + g(2) * t14 + g(3) * t30;
t54 = -t115 * t76 - t131 * t75;
t33 = t116 * t76 + t54 * t89;
t53 = t108 * t76 - t75 * t87;
t18 = t33 * t85 - t53 * t88;
t56 = -t115 * t77 - t131 * t93;
t35 = t116 * t77 + t56 * t89;
t55 = t108 * t77 - t87 * t93;
t20 = t35 * t85 - t55 * t88;
t73 = (-t107 * t87 + t131 * t132) * t84;
t58 = t102 * t86 + t73 * t89;
t72 = (t107 * t131 + t132 * t87) * t84;
t36 = t58 * t85 - t72 * t88;
t98 = g(1) * t20 + g(2) * t18 + g(3) * t36;
t96 = g(1) * t29 + g(2) * t25 + g(3) * t44;
t32 = -t114 * t76 + t54 * t86;
t34 = -t114 * t77 + t56 * t86;
t57 = -t102 * t89 + t73 * t86;
t95 = g(1) * t34 + g(2) * t32 + g(3) * t57;
t37 = t58 * t88 + t72 * t85;
t31 = -t125 * t63 + t64 * t85;
t21 = t35 * t88 + t55 * t85;
t19 = t33 * t88 + t53 * t85;
t17 = -t125 * t49 + t50 * t85;
t15 = -t125 * t45 + t46 * t85;
t7 = t94 * t86;
t6 = t97 * t88;
t5 = t97 * t85;
t4 = -g(1) * t21 - g(2) * t19 - g(3) * t37;
t3 = g(1) * t9 - g(2) * t13;
t2 = -g(1) * t17 - g(2) * t15 - g(3) * t31;
t10 = [0, g(1) * t130 - g(2) * t90, g(1) * t90 + g(2) * t130, 0, 0, 0, 0, 0, g(1) * t76 - g(2) * t77, -g(1) * t75 + g(2) * t93, 0, 0, 0, 0, 0, g(1) * t46 - g(2) * t50, -g(1) * t45 + g(2) * t49, 0, 0, 0, 0, 0, g(1) * t25 - g(2) * t29, t111, 0, 0, 0, 0, 0, t3, t112, t3, -t111, -t112, -g(1) * (-pkin(1) * t130 - t76 * pkin(2) - pkin(3) * t46 - pkin(4) * t25 - pkin(5) * t9 - pkin(11) * t45 + t26 * pkin(12) - qJ(6) * t8 + t139 * t90) - g(2) * (t90 * pkin(1) + t77 * pkin(2) + t50 * pkin(3) + t29 * pkin(4) + t13 * pkin(5) + t49 * pkin(11) + t28 * pkin(12) + t12 * qJ(6) + t130 * t139) + (g(1) * t65 + g(2) * t91) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t93 + g(2) * t75 - g(3) * t121, g(1) * t77 + g(2) * t76 + g(3) * t120, 0, 0, 0, 0, 0, -g(1) * t56 - g(2) * t54 - g(3) * t73, g(1) * t55 + g(2) * t53 + g(3) * t72, 0, 0, 0, 0, 0, -g(1) * t35 - g(2) * t33 - g(3) * t58, t95, 0, 0, 0, 0, 0, t4, t98, t4, -t95, -t98, -g(1) * (-pkin(2) * t93 + t56 * pkin(3) + t35 * pkin(4) + t21 * pkin(5) + t55 * pkin(11) + t34 * pkin(12) + t20 * qJ(6) + t119 * t77) - g(2) * (-t75 * pkin(2) + t54 * pkin(3) + t33 * pkin(4) + t19 * pkin(5) + t53 * pkin(11) + t32 * pkin(12) + t18 * qJ(6) + t119 * t76) - g(3) * (t73 * pkin(3) + t58 * pkin(4) + t37 * pkin(5) + t72 * pkin(11) + t57 * pkin(12) + t36 * qJ(6) + (pkin(2) * t132 + pkin(10) * t106) * t84); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, g(1) * t50 + g(2) * t46 + g(3) * t64, 0, 0, 0, 0, 0, t94 * t89, -t7, 0, 0, 0, 0, 0, t2, t99, t2, t7, -t99, -g(1) * (t17 * pkin(5) + t50 * pkin(11) + t16 * qJ(6)) - g(2) * (t15 * pkin(5) + t46 * pkin(11) + t14 * qJ(6)) - g(3) * (t31 * pkin(5) + t64 * pkin(11) + t30 * qJ(6)) + t94 * (pkin(4) * t89 + pkin(12) * t86 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t96, 0, 0, 0, 0, 0, t6, -t5, t6, -t96, t5, -t96 * pkin(12) + t97 * (pkin(5) * t88 + qJ(6) * t85 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t100, t1, 0, -t100, -g(1) * (-pkin(5) * t12 + qJ(6) * t13) - g(2) * (-pkin(5) * t8 + qJ(6) * t9) - g(3) * (-pkin(5) * t22 + qJ(6) * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t10;
