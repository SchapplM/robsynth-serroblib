% Calculate minimal parameter regressor of potential energy for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:51
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:50:10
% EndTime: 2021-01-16 00:50:10
% DurationCPUTime: 0.34s
% Computational Cost: add. (327->88), mult. (814->148), div. (0->0), fcn. (1051->14), ass. (0->63)
t116 = sin(qJ(3));
t119 = cos(qJ(3));
t106 = sin(pkin(11));
t110 = cos(pkin(11));
t105 = sin(pkin(12));
t111 = cos(pkin(7));
t131 = t105 * t111;
t109 = cos(pkin(12));
t112 = cos(pkin(6));
t125 = t112 * t111;
t107 = sin(pkin(7));
t108 = sin(pkin(6));
t129 = t108 * t107;
t89 = t109 * t125 - t129;
t83 = t89 * t106 + t110 * t131;
t130 = t105 * t112;
t92 = -t106 * t130 + t110 * t109;
t142 = t83 * t116 - t92 * t119;
t86 = t106 * t131 - t89 * t110;
t91 = t106 * t109 + t110 * t130;
t141 = t86 * t116 - t91 * t119;
t74 = t92 * t116 + t83 * t119;
t77 = t91 * t116 + t86 * t119;
t114 = sin(qJ(5));
t140 = t74 * t114;
t139 = t77 * t114;
t126 = t112 * t107;
t128 = t109 * t111;
t80 = -t119 * t126 + (t105 * t116 - t119 * t128) * t108;
t138 = t80 * t114;
t133 = qJ(2) * t108;
t132 = t105 * t107;
t127 = t111 * t108;
t124 = (t109 * t127 + t126) * t116 + t108 * t105 * t119;
t101 = t111 * pkin(8) + qJ(2);
t94 = t109 * t107 * pkin(8) - t105 * pkin(2);
t123 = t108 * t101 + t94 * t112;
t98 = pkin(3) * t128 + t105 * pkin(9);
t122 = pkin(3) * t129 - t98 * t112;
t97 = -t105 * pkin(3) + pkin(9) * t128;
t121 = pkin(9) * t129 - t97 * t112;
t115 = sin(qJ(4));
t118 = cos(qJ(4));
t87 = t109 * t126 + t127;
t81 = t87 * t106 + t110 * t132;
t70 = t142 * t115 + t81 * t118;
t82 = -t106 * t132 + t87 * t110;
t73 = t141 * t115 - t82 * t118;
t88 = t109 * t129 - t125;
t79 = -t124 * t115 - t118 * t88;
t120 = g(1) * t70 + g(2) * t73 + g(3) * t79;
t117 = cos(qJ(5));
t113 = -qJ(6) - pkin(10);
t102 = t117 * pkin(5) + pkin(4);
t96 = pkin(3) * t131 - t109 * pkin(9);
t95 = t109 * pkin(3) + pkin(9) * t131;
t93 = t109 * pkin(2) + pkin(8) * t132 + pkin(1);
t78 = -t115 * t88 + t124 * t118;
t72 = t81 * t115 - t118 * t142;
t71 = -t82 * t115 - t118 * t141;
t69 = -g(1) * (t72 * t117 + t140) - g(2) * (t71 * t117 + t139) - g(3) * (t78 * t117 + t138);
t68 = -g(1) * (-t72 * t114 + t74 * t117) - g(2) * (-t71 * t114 + t77 * t117) - g(3) * (-t78 * t114 + t80 * t117);
t1 = [-g(3) * qJ(1), -g(1) * (t110 * pkin(1) + t106 * t133) - g(2) * (t106 * pkin(1) - t110 * t133) - g(3) * (t112 * qJ(2) + qJ(1)), 0, g(1) * t142 + g(2) * t141 - g(3) * t124, g(1) * t74 + g(2) * t77 + g(3) * t80, 0, 0, 0, 0, 0, -g(1) * t72 - g(2) * t71 - g(3) * t78, -t120, 0, 0, 0, 0, 0, t69, t68, t69, t68, t120, -g(1) * (t72 * t102 + t70 * t113 + pkin(5) * t140 + (-t121 * t106 + t110 * t95) * t119 + (t122 * t106 - t110 * t96) * t116 + t123 * t106 + t93 * t110) - g(2) * (t71 * t102 + t73 * t113 + pkin(5) * t139 + (t106 * t95 + t121 * t110) * t119 + (-t106 * t96 - t122 * t110) * t116 - t123 * t110 + t93 * t106) - g(3) * (t78 * t102 + t79 * t113 + pkin(5) * t138 + (-pkin(9) * t126 - t97 * t108) * t119 + (pkin(3) * t126 + t98 * t108) * t116 - t94 * t108 + t101 * t112 + qJ(1));];
U_reg = t1;
