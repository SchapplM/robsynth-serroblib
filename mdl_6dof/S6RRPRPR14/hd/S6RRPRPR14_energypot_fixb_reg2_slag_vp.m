% Calculate inertial parameters regressor of potential energy for
% S6RRPRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPR14_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:36:55
% EndTime: 2019-03-09 11:36:55
% DurationCPUTime: 0.21s
% Computational Cost: add. (225->86), mult. (502->110), div. (0->0), fcn. (592->10), ass. (0->53)
t106 = cos(pkin(6));
t105 = sin(pkin(6));
t109 = sin(qJ(2));
t138 = t105 * t109;
t146 = t106 * pkin(3) + pkin(9) * t138;
t108 = sin(qJ(4));
t112 = cos(qJ(4));
t113 = cos(qJ(2));
t136 = t105 * t113;
t88 = t106 * t108 + t112 * t136;
t89 = t106 * t112 - t108 * t136;
t145 = t89 * pkin(4) + t88 * qJ(5);
t114 = cos(qJ(1));
t135 = t105 * t114;
t131 = t114 * t113;
t110 = sin(qJ(1));
t134 = t110 * t109;
t90 = -t106 * t131 + t134;
t80 = t108 * t135 + t90 * t112;
t81 = -t90 * t108 + t112 * t135;
t144 = -t81 * pkin(4) - t80 * qJ(5);
t93 = -t106 * t134 + t131;
t143 = t93 * pkin(9);
t142 = t106 * pkin(8) + pkin(7);
t137 = t105 * t110;
t139 = t114 * pkin(1) + pkin(8) * t137;
t133 = t110 * t113;
t132 = t114 * t109;
t130 = pkin(8) * t135;
t129 = pkin(2) * t138 + t142;
t128 = (-pkin(3) - pkin(8)) * t114;
t103 = t110 * pkin(1);
t91 = t106 * t132 + t133;
t127 = t91 * pkin(2) + t90 * qJ(3) + t103;
t126 = g(1) * t110 - g(2) * t114;
t92 = t106 * t133 + t132;
t125 = t93 * pkin(2) + t92 * qJ(3) + t139;
t124 = pkin(3) * t137 + t125;
t123 = -qJ(3) * t136 + t129;
t78 = t108 * t137 - t92 * t112;
t122 = g(1) * t78 - g(2) * t80 + g(3) * t88;
t79 = t92 * t108 + t112 * t137;
t121 = g(1) * t79 - g(2) * t81 + g(3) * t89;
t120 = t91 * pkin(9) + t127;
t119 = g(1) * t93 + g(2) * t91 + g(3) * t138;
t118 = -g(1) * t92 - g(2) * t90 + g(3) * t136;
t117 = t123 + t146;
t116 = t79 * pkin(4) + t78 * qJ(5) + t124;
t115 = t105 * t128 + t120;
t111 = cos(qJ(6));
t107 = sin(qJ(6));
t82 = -g(3) * t106 - t126 * t105;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t114 - g(2) * t110, t126, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t119, -t118, t82, -g(1) * t139 - g(2) * (t103 - t130) - g(3) * t142, 0, 0, 0, 0, 0, 0, t82, t119, t118, -g(1) * t125 - g(2) * (t127 - t130) - g(3) * t123, 0, 0, 0, 0, 0, 0, -t121, t122, -t119, -g(1) * (t124 + t143) - g(2) * t115 - g(3) * t117, 0, 0, 0, 0, 0, 0, -t119, t121, -t122, -g(1) * (t116 + t143) - g(2) * (t115 + t144) - g(3) * (t117 + t145) 0, 0, 0, 0, 0, 0, -g(1) * (t78 * t107 + t93 * t111) - g(2) * (-t80 * t107 + t91 * t111) - g(3) * (t88 * t107 + t111 * t138) -g(1) * (-t93 * t107 + t78 * t111) - g(2) * (-t91 * t107 - t80 * t111) - g(3) * (-t107 * t138 + t88 * t111) -t121, -g(1) * (t79 * pkin(10) + (pkin(5) + pkin(9)) * t93 + t116) - g(2) * (t91 * pkin(5) - t81 * pkin(10) + t120 + t144) - g(3) * (t89 * pkin(10) + t129 + t145 + t146) + (-g(3) * (pkin(5) * t109 - qJ(3) * t113) - g(2) * t128) * t105;];
U_reg  = t1;
