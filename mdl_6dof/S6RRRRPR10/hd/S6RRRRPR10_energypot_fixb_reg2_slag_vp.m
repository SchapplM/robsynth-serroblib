% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:08:30
% EndTime: 2019-03-09 23:08:30
% DurationCPUTime: 0.21s
% Computational Cost: add. (308->91), mult. (511->130), div. (0->0), fcn. (603->12), ass. (0->51)
t114 = cos(pkin(6));
t146 = t114 * pkin(8) + pkin(7);
t113 = sin(pkin(6));
t117 = sin(qJ(2));
t145 = t113 * t117;
t118 = sin(qJ(1));
t144 = t113 * t118;
t120 = cos(qJ(3));
t143 = t113 * t120;
t121 = cos(qJ(2));
t142 = t113 * t121;
t122 = cos(qJ(1));
t141 = t113 * t122;
t116 = sin(qJ(3));
t140 = t114 * t116;
t139 = t118 * t117;
t138 = t118 * t121;
t137 = t122 * t117;
t136 = t122 * t121;
t135 = t122 * pkin(1) + pkin(8) * t144;
t134 = t116 * t144;
t110 = t118 * pkin(1);
t133 = -pkin(8) * t141 + t110;
t106 = t120 * pkin(3) + pkin(2);
t123 = -pkin(10) - pkin(9);
t132 = pkin(3) * t140 + t106 * t145 + t123 * t142 + t146;
t131 = g(1) * t118 - g(2) * t122;
t95 = t114 * t138 + t137;
t96 = -t114 * t139 + t136;
t130 = pkin(3) * t134 + t96 * t106 - t95 * t123 + t135;
t112 = qJ(3) + qJ(4);
t107 = sin(t112);
t108 = cos(t112);
t94 = t114 * t137 + t138;
t82 = t94 * t107 + t108 * t141;
t84 = t96 * t107 - t108 * t144;
t89 = t107 * t145 - t114 * t108;
t129 = g(1) * t84 + g(2) * t82 + g(3) * t89;
t83 = -t107 * t141 + t94 * t108;
t85 = t107 * t144 + t96 * t108;
t90 = t114 * t107 + t108 * t145;
t128 = g(1) * t85 + g(2) * t83 + g(3) * t90;
t93 = -t114 * t136 + t139;
t79 = -g(1) * t95 - g(2) * t93 + g(3) * t142;
t127 = t90 * pkin(4) + t89 * qJ(5) + t132;
t126 = t85 * pkin(4) + t84 * qJ(5) + t130;
t125 = t110 + t94 * t106 - t93 * t123 + (-pkin(3) * t116 - pkin(8)) * t141;
t124 = t83 * pkin(4) + t82 * qJ(5) + t125;
t119 = cos(qJ(6));
t115 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t122 - g(2) * t118, t131, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t145, -t79, -g(3) * t114 - t131 * t113, -g(1) * t135 - g(2) * t133 - g(3) * t146, 0, 0, 0, 0, 0, 0, -g(1) * (t96 * t120 + t134) - g(2) * (-t116 * t141 + t94 * t120) - g(3) * (t117 * t143 + t140) -g(1) * (-t96 * t116 + t118 * t143) - g(2) * (-t94 * t116 - t120 * t141) - g(3) * (t114 * t120 - t116 * t145) t79, -g(1) * (t96 * pkin(2) + t95 * pkin(9) + t135) - g(2) * (t94 * pkin(2) + t93 * pkin(9) + t133) - g(3) * ((pkin(2) * t117 - pkin(9) * t121) * t113 + t146) 0, 0, 0, 0, 0, 0, -t128, t129, t79, -g(1) * t130 - g(2) * t125 - g(3) * t132, 0, 0, 0, 0, 0, 0, t79, t128, -t129, -g(1) * t126 - g(2) * t124 - g(3) * t127, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t115 + t95 * t119) - g(2) * (t82 * t115 + t93 * t119) - g(3) * (t89 * t115 - t119 * t142) -g(1) * (-t95 * t115 + t84 * t119) - g(2) * (-t93 * t115 + t82 * t119) - g(3) * (t115 * t142 + t89 * t119) -t128, -g(1) * (t95 * pkin(5) + t85 * pkin(11) + t126) - g(2) * (t93 * pkin(5) + t83 * pkin(11) + t124) - g(3) * (-pkin(5) * t142 + t90 * pkin(11) + t127);];
U_reg  = t1;
