% Calculate inertial parameters regressor of potential energy for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6PPRRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:01:29
% EndTime: 2019-03-08 19:01:29
% DurationCPUTime: 0.37s
% Computational Cost: add. (502->110), mult. (1179->171), div. (0->0), fcn. (1486->16), ass. (0->64)
t118 = sin(pkin(7));
t122 = cos(pkin(7));
t123 = cos(pkin(6));
t119 = sin(pkin(6));
t120 = cos(pkin(13));
t151 = t119 * t120;
t97 = -t118 * t151 + t123 * t122;
t116 = sin(pkin(13));
t121 = cos(pkin(12));
t117 = sin(pkin(12));
t152 = t117 * t123;
t100 = -t121 * t116 - t120 * t152;
t149 = t119 * t122;
t90 = -t100 * t118 + t117 * t149;
t159 = cos(qJ(3));
t125 = sin(qJ(4));
t148 = t121 * t123;
t98 = -t117 * t116 + t120 * t148;
t89 = -t98 * t118 - t121 * t149;
t158 = t89 * t125;
t157 = t90 * t125;
t156 = t97 * t125;
t154 = t116 * t119;
t153 = t117 * t119;
t150 = t119 * t121;
t146 = t123 * qJ(2) + qJ(1);
t145 = t121 * pkin(1) + qJ(2) * t153;
t142 = t118 * t159;
t141 = t122 * t159;
t140 = t119 * t142;
t139 = g(1) * t117 - g(2) * t121;
t138 = t117 * pkin(1) - qJ(2) * t150;
t115 = qJ(4) + qJ(5);
t111 = sin(t115);
t112 = cos(t115);
t126 = sin(qJ(3));
t99 = t116 * t148 + t117 * t120;
t78 = t99 * t159 + (-t118 * t150 + t122 * t98) * t126;
t67 = t78 * t111 - t89 * t112;
t101 = -t116 * t152 + t121 * t120;
t80 = t101 * t159 + (t100 * t122 + t118 * t153) * t126;
t69 = t80 * t111 - t90 * t112;
t88 = t123 * t118 * t126 + (t120 * t122 * t126 + t159 * t116) * t119;
t73 = t88 * t111 - t97 * t112;
t137 = g(1) * t69 + g(2) * t67 + g(3) * t73;
t77 = t121 * t140 + t99 * t126 - t98 * t141;
t79 = -t100 * t141 + t101 * t126 - t117 * t140;
t87 = -t123 * t142 + t126 * t154 - t141 * t151;
t136 = g(1) * t79 + g(2) * t77 + g(3) * t87;
t135 = t101 * pkin(2) + t90 * pkin(8) + t145;
t134 = pkin(2) * t154 + t97 * pkin(8) + t146;
t128 = cos(qJ(4));
t109 = t128 * pkin(4) + pkin(3);
t129 = -pkin(10) - pkin(9);
t133 = pkin(4) * t157 + t80 * t109 - t79 * t129 + t135;
t132 = pkin(4) * t156 + t88 * t109 - t87 * t129 + t134;
t131 = t99 * pkin(2) + t89 * pkin(8) + t138;
t130 = pkin(4) * t158 + t78 * t109 - t77 * t129 + t131;
t127 = cos(qJ(6));
t124 = sin(qJ(6));
t74 = t97 * t111 + t88 * t112;
t70 = t90 * t111 + t80 * t112;
t68 = t89 * t111 + t78 * t112;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t121 - g(2) * t117, t139, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t101 - g(2) * t99 - g(3) * t154, -g(1) * t100 - g(2) * t98 - g(3) * t151, -g(3) * t123 - t139 * t119, -g(1) * t145 - g(2) * t138 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * t80 - g(2) * t78 - g(3) * t88, t136, -g(1) * t90 - g(2) * t89 - g(3) * t97, -g(1) * t135 - g(2) * t131 - g(3) * t134, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t128 + t157) - g(2) * (t78 * t128 + t158) - g(3) * (t88 * t128 + t156) -g(1) * (-t80 * t125 + t90 * t128) - g(2) * (-t78 * t125 + t89 * t128) - g(3) * (-t88 * t125 + t97 * t128) -t136, -g(1) * (t80 * pkin(3) + t79 * pkin(9) + t135) - g(2) * (t78 * pkin(3) + t77 * pkin(9) + t131) - g(3) * (t88 * pkin(3) + t87 * pkin(9) + t134) 0, 0, 0, 0, 0, 0, -g(1) * t70 - g(2) * t68 - g(3) * t74, t137, -t136, -g(1) * t133 - g(2) * t130 - g(3) * t132, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t124 + t70 * t127) - g(2) * (t77 * t124 + t68 * t127) - g(3) * (t87 * t124 + t74 * t127) -g(1) * (-t70 * t124 + t79 * t127) - g(2) * (-t68 * t124 + t77 * t127) - g(3) * (-t74 * t124 + t87 * t127) -t137, -g(1) * (t70 * pkin(5) + t69 * pkin(11) + t133) - g(2) * (t68 * pkin(5) + t67 * pkin(11) + t130) - g(3) * (t74 * pkin(5) + t73 * pkin(11) + t132);];
U_reg  = t1;
