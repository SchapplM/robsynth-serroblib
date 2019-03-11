% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:09:23
% EndTime: 2019-03-09 11:09:24
% DurationCPUTime: 0.20s
% Computational Cost: add. (308->91), mult. (511->129), div. (0->0), fcn. (603->12), ass. (0->50)
t115 = cos(pkin(6));
t144 = t115 * pkin(8) + pkin(7);
t113 = sin(pkin(6));
t118 = sin(qJ(2));
t143 = t113 * t118;
t119 = sin(qJ(1));
t142 = t113 * t119;
t121 = cos(qJ(2));
t141 = t113 * t121;
t122 = cos(qJ(1));
t140 = t113 * t122;
t112 = sin(pkin(11));
t139 = t115 * t112;
t138 = t119 * t118;
t137 = t119 * t121;
t136 = t122 * t118;
t135 = t122 * t121;
t134 = t122 * pkin(1) + pkin(8) * t142;
t133 = t112 * t142;
t109 = t119 * pkin(1);
t132 = -pkin(8) * t140 + t109;
t114 = cos(pkin(11));
t105 = t114 * pkin(3) + pkin(2);
t116 = -pkin(9) - qJ(3);
t94 = t115 * t137 + t136;
t95 = -t115 * t138 + t135;
t131 = pkin(3) * t133 + t95 * t105 - t94 * t116 + t134;
t130 = pkin(3) * t139 + t105 * t143 + t116 * t141 + t144;
t129 = g(1) * t119 - g(2) * t122;
t111 = pkin(11) + qJ(4);
t106 = sin(t111);
t107 = cos(t111);
t93 = t115 * t136 + t137;
t81 = t93 * t106 + t107 * t140;
t83 = t95 * t106 - t107 * t142;
t88 = t106 * t143 - t115 * t107;
t128 = g(1) * t83 + g(2) * t81 + g(3) * t88;
t82 = -t106 * t140 + t93 * t107;
t84 = t106 * t142 + t95 * t107;
t89 = t115 * t106 + t107 * t143;
t127 = g(1) * t84 + g(2) * t82 + g(3) * t89;
t92 = -t115 * t135 + t138;
t78 = -g(1) * t94 - g(2) * t92 + g(3) * t141;
t126 = t89 * pkin(4) + t88 * qJ(5) + t130;
t125 = t84 * pkin(4) + t83 * qJ(5) + t131;
t124 = t109 + t93 * t105 - t92 * t116 + (-pkin(3) * t112 - pkin(8)) * t140;
t123 = t82 * pkin(4) + t81 * qJ(5) + t124;
t120 = cos(qJ(6));
t117 = sin(qJ(6));
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t122 - g(2) * t119, t129, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t95 - g(2) * t93 - g(3) * t143, -t78, -g(3) * t115 - t129 * t113, -g(1) * t134 - g(2) * t132 - g(3) * t144, 0, 0, 0, 0, 0, 0, -g(1) * (t95 * t114 + t133) - g(2) * (-t112 * t140 + t93 * t114) - g(3) * (t114 * t143 + t139) -g(1) * (-t95 * t112 + t114 * t142) - g(2) * (-t93 * t112 - t114 * t140) - g(3) * (-t112 * t143 + t115 * t114) t78, -g(1) * (t95 * pkin(2) + t94 * qJ(3) + t134) - g(2) * (t93 * pkin(2) + t92 * qJ(3) + t132) - g(3) * ((pkin(2) * t118 - qJ(3) * t121) * t113 + t144) 0, 0, 0, 0, 0, 0, -t127, t128, t78, -g(1) * t131 - g(2) * t124 - g(3) * t130, 0, 0, 0, 0, 0, 0, t78, t127, -t128, -g(1) * t125 - g(2) * t123 - g(3) * t126, 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t117 + t94 * t120) - g(2) * (t81 * t117 + t92 * t120) - g(3) * (t88 * t117 - t120 * t141) -g(1) * (-t94 * t117 + t83 * t120) - g(2) * (-t92 * t117 + t81 * t120) - g(3) * (t117 * t141 + t88 * t120) -t127, -g(1) * (t94 * pkin(5) + t84 * pkin(10) + t125) - g(2) * (t92 * pkin(5) + t82 * pkin(10) + t123) - g(3) * (-pkin(5) * t141 + t89 * pkin(10) + t126);];
U_reg  = t1;
