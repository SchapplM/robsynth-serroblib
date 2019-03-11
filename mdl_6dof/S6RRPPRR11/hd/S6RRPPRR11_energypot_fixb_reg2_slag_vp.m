% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:42:43
% EndTime: 2019-03-09 09:42:43
% DurationCPUTime: 0.22s
% Computational Cost: add. (260->97), mult. (470->133), div. (0->0), fcn. (547->12), ass. (0->54)
t112 = sin(qJ(1));
t141 = g(1) * t112;
t115 = cos(qJ(1));
t140 = g(2) * t115;
t108 = cos(pkin(6));
t139 = t108 * pkin(8) + pkin(7);
t105 = sin(pkin(11));
t114 = cos(qJ(2));
t128 = t115 * t114;
t111 = sin(qJ(2));
t131 = t112 * t111;
t88 = -t108 * t128 + t131;
t138 = t88 * t105;
t129 = t115 * t111;
t130 = t112 * t114;
t90 = t108 * t130 + t129;
t137 = t90 * t105;
t106 = sin(pkin(6));
t134 = t106 * t112;
t136 = t115 * pkin(1) + pkin(8) * t134;
t135 = t106 * t111;
t133 = t106 * t114;
t132 = t106 * t115;
t127 = pkin(2) * t135 + t139;
t126 = pkin(8) * t132;
t107 = cos(pkin(11));
t98 = t107 * pkin(4) + pkin(3);
t125 = t108 * t98 + t127;
t102 = t112 * pkin(1);
t89 = t108 * t129 + t130;
t124 = t89 * pkin(2) + t88 * qJ(3) + t102;
t123 = -t140 + t141;
t91 = -t108 * t131 + t128;
t122 = t91 * pkin(2) + t90 * qJ(3) + t136;
t109 = -pkin(9) - qJ(4);
t121 = pkin(4) * t138 - t89 * t109 + t124;
t104 = pkin(11) + qJ(5);
t100 = cos(t104);
t99 = sin(t104);
t74 = -t90 * t100 + t99 * t134;
t76 = t88 * t100 + t99 * t132;
t79 = t100 * t133 + t108 * t99;
t120 = g(1) * t74 - g(2) * t76 + g(3) * t79;
t119 = g(1) * t91 + g(2) * t89 + g(3) * t135;
t118 = -g(1) * t90 - g(2) * t88 + g(3) * t133;
t117 = pkin(4) * t137 - t91 * t109 + t98 * t134 + t122;
t116 = (-g(3) * (-t109 * t111 + (-pkin(4) * t105 - qJ(3)) * t114) - (-pkin(8) - t98) * t140) * t106;
t113 = cos(qJ(6));
t110 = sin(qJ(6));
t83 = -g(3) * t108 - t123 * t106;
t80 = t108 * t100 - t99 * t133;
t77 = -t100 * t132 + t88 * t99;
t75 = t100 * t134 + t90 * t99;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t115 - g(2) * t112, t123, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t119, -t118, t83, -g(1) * t136 - g(2) * (t102 - t126) - g(3) * t139, 0, 0, 0, 0, 0, 0, t83, t119, t118, -g(1) * t122 - g(2) * (t124 - t126) - g(3) * (-qJ(3) * t133 + t127) 0, 0, 0, 0, 0, 0, -g(1) * (t107 * t134 + t137) - g(2) * (-t107 * t132 + t138) - g(3) * (-t105 * t133 + t108 * t107) -g(1) * (-t105 * t134 + t90 * t107) - g(2) * (t105 * t132 + t88 * t107) - g(3) * (-t108 * t105 - t107 * t133) -t119, -g(1) * (t91 * qJ(4) + t122) - g(2) * (t89 * qJ(4) + t124) - g(3) * (t108 * pkin(3) + t127) + (-pkin(3) * t141 - g(3) * (-qJ(3) * t114 + qJ(4) * t111) - (-pkin(3) - pkin(8)) * t140) * t106, 0, 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t77 - g(3) * t80, t120, -t119, -g(1) * t117 - g(2) * t121 - g(3) * t125 + t116, 0, 0, 0, 0, 0, 0, -g(1) * (t91 * t110 + t75 * t113) - g(2) * (t89 * t110 + t77 * t113) - g(3) * (t110 * t135 + t80 * t113) -g(1) * (-t75 * t110 + t91 * t113) - g(2) * (-t77 * t110 + t89 * t113) - g(3) * (-t80 * t110 + t113 * t135) -t120, -g(1) * (t75 * pkin(5) + t74 * pkin(10) + t117) - g(2) * (t77 * pkin(5) - t76 * pkin(10) + t121) - g(3) * (t80 * pkin(5) + t79 * pkin(10) + t125) + t116;];
U_reg  = t1;
