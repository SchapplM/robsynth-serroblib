% Calculate inertial parameters regressor of potential energy for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPRR14_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:20:34
% EndTime: 2019-03-09 20:20:34
% DurationCPUTime: 0.23s
% Computational Cost: add. (277->96), mult. (600->131), div. (0->0), fcn. (727->12), ass. (0->54)
t146 = pkin(4) + pkin(9);
t111 = sin(qJ(2));
t112 = sin(qJ(1));
t114 = cos(qJ(2));
t115 = cos(qJ(1));
t139 = cos(pkin(6));
t128 = t115 * t139;
t91 = t112 * t111 - t114 * t128;
t145 = t91 * pkin(9);
t129 = t112 * t139;
t93 = t115 * t111 + t114 * t129;
t144 = t93 * pkin(9);
t143 = cos(qJ(3));
t113 = cos(qJ(5));
t142 = t113 * pkin(5) + t146;
t141 = t139 * pkin(8) + pkin(7);
t108 = sin(pkin(6));
t137 = t108 * t112;
t140 = t115 * pkin(1) + pkin(8) * t137;
t138 = t108 * t111;
t136 = t108 * t114;
t135 = t108 * t115;
t134 = pkin(2) * t138 + t141;
t133 = pkin(9) * t136;
t94 = -t111 * t129 + t115 * t114;
t132 = t94 * pkin(2) + t140;
t131 = t108 * t143;
t109 = sin(qJ(5));
t130 = pkin(5) * t109 + qJ(4);
t110 = sin(qJ(3));
t90 = t139 * t110 + t111 * t131;
t127 = t90 * pkin(3) + t134;
t85 = t110 * t137 + t94 * t143;
t126 = t85 * pkin(3) + t132;
t125 = t112 * pkin(1) - pkin(8) * t135;
t124 = g(1) * t112 - g(2) * t115;
t92 = t111 * t128 + t112 * t114;
t123 = t92 * pkin(2) + t125;
t89 = t110 * t138 - t139 * t143;
t122 = t89 * qJ(4) + t127;
t83 = -t110 * t135 + t92 * t143;
t121 = t83 * pkin(3) + t123;
t84 = t94 * t110 - t112 * t131;
t120 = t84 * qJ(4) + t126;
t82 = t92 * t110 + t115 * t131;
t119 = g(1) * t84 + g(2) * t82 + g(3) * t89;
t118 = g(1) * t85 + g(2) * t83 + g(3) * t90;
t79 = -g(1) * t93 - g(2) * t91 + g(3) * t136;
t117 = t82 * qJ(4) + t121;
t116 = -pkin(11) - pkin(10);
t107 = qJ(5) + qJ(6);
t103 = cos(t107);
t102 = sin(t107);
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t115 - g(2) * t112, t124, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t94 - g(2) * t92 - g(3) * t138, -t79, -g(3) * t139 - t124 * t108, -g(1) * t140 - g(2) * t125 - g(3) * t141, 0, 0, 0, 0, 0, 0, -t118, t119, t79, -g(1) * (t132 + t144) - g(2) * (t123 + t145) - g(3) * (-t133 + t134) 0, 0, 0, 0, 0, 0, t79, t118, -t119, -g(1) * (t120 + t144) - g(2) * (t117 + t145) - g(3) * (t122 - t133) 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t109 + t93 * t113) - g(2) * (t82 * t109 + t91 * t113) - g(3) * (t89 * t109 - t113 * t136) -g(1) * (-t93 * t109 + t84 * t113) - g(2) * (-t91 * t109 + t82 * t113) - g(3) * (t109 * t136 + t89 * t113) -t118, -g(1) * (t85 * pkin(10) + t146 * t93 + t120) - g(2) * (t83 * pkin(10) + t146 * t91 + t117) - g(3) * (t90 * pkin(10) - t146 * t136 + t122) 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t102 + t93 * t103) - g(2) * (t82 * t102 + t91 * t103) - g(3) * (t89 * t102 - t103 * t136) -g(1) * (-t93 * t102 + t84 * t103) - g(2) * (-t91 * t102 + t82 * t103) - g(3) * (t102 * t136 + t89 * t103) -t118, -g(1) * (-t85 * t116 + t130 * t84 + t142 * t93 + t126) - g(2) * (-t83 * t116 + t130 * t82 + t142 * t91 + t121) - g(3) * (-t90 * t116 + t130 * t89 - t142 * t136 + t127);];
U_reg  = t1;
