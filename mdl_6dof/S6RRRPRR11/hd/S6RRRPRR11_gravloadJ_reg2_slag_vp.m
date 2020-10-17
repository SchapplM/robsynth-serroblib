% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 14:32:31
% EndTime: 2019-05-07 14:32:35
% DurationCPUTime: 1.20s
% Computational Cost: add. (788->163), mult. (2057->232), div. (0->0), fcn. (2585->12), ass. (0->99)
t128 = cos(pkin(6));
t141 = cos(qJ(1));
t110 = t128 * t141;
t140 = sin(qJ(1));
t86 = sin(qJ(2));
t90 = cos(qJ(2));
t64 = -t90 * t110 + t140 * t86;
t82 = sin(pkin(6));
t124 = t82 * t141;
t65 = t86 * t110 + t140 * t90;
t85 = sin(qJ(3));
t89 = cos(qJ(3));
t40 = t89 * t124 + t65 * t85;
t41 = -t85 * t124 + t65 * t89;
t84 = sin(qJ(5));
t88 = cos(qJ(5));
t8 = t40 * t84 + t41 * t88;
t83 = sin(qJ(6));
t87 = cos(qJ(6));
t158 = t64 * t87 + t8 * t83;
t157 = -t64 * t83 + t8 * t87;
t123 = t82 * t140;
t109 = t128 * t140;
t67 = -t86 * t109 + t141 * t90;
t44 = -t89 * t123 + t67 * t85;
t45 = t85 * t123 + t67 * t89;
t108 = -t44 * t88 + t45 * t84;
t151 = t40 * t88 - t41 * t84;
t135 = t82 * t86;
t62 = -t128 * t89 + t85 * t135;
t63 = t128 * t85 + t89 * t135;
t26 = t62 * t88 - t63 * t84;
t150 = g(1) * t108 - g(2) * t151 - g(3) * t26;
t156 = t150 * t83;
t155 = t150 * t87;
t14 = t44 * t84 + t45 * t88;
t154 = -pkin(5) * t108 + t14 * pkin(11);
t153 = pkin(5) * t151 + t8 * pkin(11);
t27 = t62 * t84 + t63 * t88;
t152 = t26 * pkin(5) + t27 * pkin(11);
t101 = g(1) * t14 + g(2) * t8 + g(3) * t27;
t129 = qJ(4) * t85;
t137 = t64 * t89;
t147 = -pkin(3) * t137 - t64 * t129;
t66 = t90 * t109 + t141 * t86;
t136 = t66 * t89;
t146 = -pkin(3) * t136 - t66 * t129;
t145 = -t84 * t89 + t85 * t88;
t144 = pkin(9) - pkin(10);
t143 = t64 * pkin(9);
t142 = t66 * pkin(9);
t134 = t82 * t90;
t131 = pkin(2) * t134 + pkin(9) * t135;
t130 = t141 * pkin(1) + pkin(8) * t123;
t127 = t85 * t134;
t126 = t89 * t134;
t125 = t67 * pkin(2) + t130;
t58 = t64 * pkin(2);
t122 = t65 * pkin(9) - t58;
t60 = t66 * pkin(2);
t121 = t67 * pkin(9) - t60;
t120 = -t40 * pkin(3) + t41 * qJ(4);
t119 = -t44 * pkin(3) + t45 * qJ(4);
t118 = -t62 * pkin(3) + t63 * qJ(4);
t117 = pkin(3) * t126 + qJ(4) * t127 + t131;
t116 = g(1) * t151 + g(2) * t108;
t115 = -t40 * pkin(4) + t120;
t114 = -t44 * pkin(4) + t119;
t113 = -t62 * pkin(4) + t118;
t112 = -t140 * pkin(1) + pkin(8) * t124;
t111 = -g(1) * t40 + g(2) * t44;
t29 = g(1) * t64 - g(2) * t66;
t107 = t84 * t85 + t88 * t89;
t106 = -t65 * pkin(2) + t112;
t105 = t45 * pkin(3) + t44 * qJ(4) + t125;
t18 = t145 * t64;
t20 = t145 * t66;
t46 = t84 * t126 - t88 * t127;
t100 = g(1) * t20 + g(2) * t18 + g(3) * t46;
t4 = g(1) * t44 + g(2) * t40 + g(3) * t62;
t99 = g(1) * t45 + g(2) * t41 + g(3) * t63;
t98 = g(1) * t141 + g(2) * t140;
t97 = -pkin(4) * t137 + t144 * t65 + t147 - t58;
t96 = -pkin(4) * t136 + t144 * t67 + t146 - t60;
t95 = -g(1) * t66 - g(2) * t64 + g(3) * t134;
t23 = g(1) * t67 + g(2) * t65 + g(3) * t135;
t94 = pkin(4) * t126 - pkin(10) * t135 + t117;
t93 = -pkin(3) * t41 - qJ(4) * t40 + t106;
t92 = t45 * pkin(4) + t144 * t66 + t105;
t91 = -pkin(4) * t41 - t144 * t64 + t93;
t47 = t107 * t134;
t21 = t107 * t66;
t19 = t107 * t64;
t17 = t95 * t89;
t16 = t95 * t85;
t15 = g(1) * t41 - g(2) * t45;
t2 = t14 * t87 - t66 * t83;
t1 = -t14 * t83 - t66 * t87;
t3 = [0, 0, 0, 0, 0, 0, g(1) * t140 - g(2) * t141, t98, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t65 - g(2) * t67, -t29, -t98 * t82, -g(1) * t112 - g(2) * t130, 0, 0, 0, 0, 0, 0, t15, t111, t29, -g(1) * (t106 - t143) - g(2) * (t125 + t142) 0, 0, 0, 0, 0, 0, t15, t29, -t111, -g(1) * (t93 - t143) - g(2) * (t105 + t142) 0, 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t14, t116, -t29, -g(1) * t91 - g(2) * t92, 0, 0, 0, 0, 0, 0, g(1) * t157 - g(2) * t2, -g(1) * t158 - g(2) * t1, -t116, -g(1) * (-pkin(5) * t8 + pkin(11) * t151 + t91) - g(2) * (t14 * pkin(5) + pkin(11) * t108 + t92); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t95, t23, 0, 0, 0, 0, 0, 0, 0, 0, -t17, t16, -t23, -g(1) * t121 - g(2) * t122 - g(3) * t131, 0, 0, 0, 0, 0, 0, -t17, -t23, -t16, -g(1) * (t121 + t146) - g(2) * (t122 + t147) - g(3) * t117, 0, 0, 0, 0, 0, 0, g(1) * t21 + g(2) * t19 - g(3) * t47, t100, t23, -g(1) * t96 - g(2) * t97 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (-t21 * t87 - t67 * t83) - g(2) * (-t19 * t87 - t65 * t83) - g(3) * (-t83 * t135 + t47 * t87) -g(1) * (t21 * t83 - t67 * t87) - g(2) * (t19 * t83 - t65 * t87) - g(3) * (-t87 * t135 - t47 * t83) -t100, -g(1) * (-t21 * pkin(5) + t20 * pkin(11) + t96) - g(2) * (-t19 * pkin(5) + t18 * pkin(11) + t97) - g(3) * (t47 * pkin(5) + t46 * pkin(11) + t94); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t99, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, -t99, -g(1) * t119 - g(2) * t120 - g(3) * t118, 0, 0, 0, 0, 0, 0, -t150, -t101, 0, -g(1) * t114 - g(2) * t115 - g(3) * t113, 0, 0, 0, 0, 0, 0, -t155, t156, t101, -g(1) * (t114 - t154) - g(2) * (t115 - t153) - g(3) * (t113 - t152); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t150, t101, 0, 0, 0, 0, 0, 0, 0, 0, t155, -t156, -t101, -g(1) * t154 - g(2) * t153 - g(3) * t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t158 - g(3) * (t87 * t134 - t27 * t83) g(1) * t2 + g(2) * t157 - g(3) * (-t83 * t134 - t27 * t87) 0, 0;];
taug_reg  = t3;
