% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP9_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:49:46
% EndTime: 2019-03-09 21:49:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (311->82), mult. (727->108), div. (0->0), fcn. (903->10), ass. (0->55)
t143 = pkin(5) + pkin(10);
t108 = sin(qJ(3));
t113 = cos(qJ(1));
t106 = sin(pkin(6));
t139 = cos(qJ(3));
t131 = t106 * t139;
t109 = sin(qJ(2));
t110 = sin(qJ(1));
t112 = cos(qJ(2));
t137 = cos(pkin(6));
t129 = t113 * t137;
t94 = t109 * t129 + t110 * t112;
t82 = t94 * t108 + t113 * t131;
t142 = t82 * pkin(10);
t130 = t110 * t137;
t96 = -t109 * t130 + t113 * t112;
t84 = t96 * t108 - t110 * t131;
t141 = t84 * pkin(10);
t136 = t106 * t109;
t91 = t108 * t136 - t137 * t139;
t140 = t91 * pkin(10);
t138 = t137 * pkin(8) + pkin(7);
t135 = t106 * t110;
t134 = t106 * t112;
t133 = t106 * t113;
t132 = t113 * pkin(1) + pkin(8) * t135;
t128 = t110 * pkin(1) - pkin(8) * t133;
t127 = g(1) * t110 - g(2) * t113;
t95 = t113 * t109 + t112 * t130;
t126 = t96 * pkin(2) + t95 * pkin(9) + t132;
t85 = t108 * t135 + t96 * t139;
t125 = t85 * pkin(3) + t126;
t124 = pkin(2) * t136 - pkin(9) * t134 + t138;
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t83 = -t108 * t133 + t94 * t139;
t93 = t110 * t109 - t112 * t129;
t73 = t83 * t107 - t93 * t111;
t75 = t85 * t107 - t95 * t111;
t92 = t137 * t108 + t109 * t131;
t80 = t92 * t107 + t111 * t134;
t123 = g(1) * t75 + g(2) * t73 + g(3) * t80;
t74 = t93 * t107 + t83 * t111;
t76 = t95 * t107 + t85 * t111;
t81 = -t107 * t134 + t92 * t111;
t122 = g(1) * t76 + g(2) * t74 + g(3) * t81;
t121 = g(1) * t84 + g(2) * t82 + g(3) * t91;
t120 = t92 * pkin(3) + t124;
t119 = t94 * pkin(2) + t93 * pkin(9) + t128;
t118 = -g(1) * t95 - g(2) * t93 + g(3) * t134;
t117 = t83 * pkin(3) + t119;
t116 = t76 * pkin(4) + t75 * qJ(5) + t125;
t115 = t81 * pkin(4) + t80 * qJ(5) + t120;
t114 = t74 * pkin(4) + t73 * qJ(5) + t117;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t113 - g(2) * t110, t127, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -g(1) * t96 - g(2) * t94 - g(3) * t136, -t118, -g(3) * t137 - t127 * t106, -g(1) * t132 - g(2) * t128 - g(3) * t138, 0, 0, 0, 0, 0, 0, -g(1) * t85 - g(2) * t83 - g(3) * t92, t121, t118, -g(1) * t126 - g(2) * t119 - g(3) * t124, 0, 0, 0, 0, 0, 0, -t122, t123, -t121, -g(1) * (t125 + t141) - g(2) * (t117 + t142) - g(3) * (t120 + t140) 0, 0, 0, 0, 0, 0, -t121, t122, -t123, -g(1) * (t116 + t141) - g(2) * (t114 + t142) - g(3) * (t115 + t140) 0, 0, 0, 0, 0, 0, -t121, -t123, -t122, -g(1) * (t76 * qJ(6) + t143 * t84 + t116) - g(2) * (t74 * qJ(6) + t143 * t82 + t114) - g(3) * (t81 * qJ(6) + t143 * t91 + t115);];
U_reg  = t1;
