% Calculate inertial parameters regressor of potential energy for
% S6RRPRRP14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRP14_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP14_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:09:49
% EndTime: 2019-03-09 13:09:50
% DurationCPUTime: 0.18s
% Computational Cost: add. (251->78), mult. (570->105), div. (0->0), fcn. (686->10), ass. (0->55)
t111 = cos(pkin(6));
t144 = t111 * pkin(8) + pkin(7);
t110 = sin(pkin(6));
t114 = sin(qJ(2));
t143 = t110 * t114;
t115 = sin(qJ(1));
t142 = t110 * t115;
t118 = cos(qJ(2));
t141 = t110 * t118;
t119 = cos(qJ(1));
t140 = t110 * t119;
t139 = t115 * t114;
t138 = t115 * t118;
t137 = t119 * t114;
t136 = t119 * t118;
t135 = t119 * pkin(1) + pkin(8) * t142;
t134 = pkin(8) * t140;
t108 = t115 * pkin(1);
t95 = -t111 * t136 + t139;
t96 = t111 * t137 + t138;
t133 = t96 * pkin(2) + t95 * qJ(3) + t108;
t132 = g(1) * t115 - g(2) * t119;
t97 = t111 * t138 + t137;
t98 = -t111 * t139 + t136;
t131 = t98 * pkin(2) + t97 * qJ(3) + t135;
t130 = pkin(2) * t143 - qJ(3) * t141 + t144;
t112 = sin(qJ(5));
t116 = cos(qJ(5));
t113 = sin(qJ(4));
t117 = cos(qJ(4));
t83 = t97 * t113 + t117 * t142;
t73 = t83 * t112 - t98 * t116;
t85 = t95 * t113 - t117 * t140;
t75 = t85 * t112 - t96 * t116;
t94 = t111 * t117 - t113 * t141;
t80 = t94 * t112 - t116 * t143;
t129 = g(1) * t73 + g(2) * t75 + g(3) * t80;
t82 = t113 * t142 - t97 * t117;
t84 = t113 * t140 + t95 * t117;
t93 = t111 * t113 + t117 * t141;
t128 = g(1) * t82 - g(2) * t84 + g(3) * t93;
t127 = g(1) * t98 + g(2) * t96 + g(3) * t143;
t126 = -g(1) * t97 - g(2) * t95 + g(3) * t141;
t125 = t111 * pkin(3) + pkin(9) * t143 + t130;
t124 = pkin(3) * t142 + t98 * pkin(9) + t131;
t123 = t96 * pkin(9) + (-pkin(3) - pkin(8)) * t140 + t133;
t122 = t94 * pkin(4) + t93 * pkin(10) + t125;
t121 = t83 * pkin(4) + t82 * pkin(10) + t124;
t120 = t85 * pkin(4) - t84 * pkin(10) + t123;
t86 = -g(3) * t111 - t132 * t110;
t81 = t112 * t143 + t94 * t116;
t76 = t96 * t112 + t85 * t116;
t74 = t98 * t112 + t83 * t116;
t71 = -g(1) * t74 - g(2) * t76 - g(3) * t81;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t119 - g(2) * t115, t132, -g(3), -g(3) * pkin(7), 0, 0, 0, 0, 0, 0, -t127, -t126, t86, -g(1) * t135 - g(2) * (t108 - t134) - g(3) * t144, 0, 0, 0, 0, 0, 0, t86, t127, t126, -g(1) * t131 - g(2) * (t133 - t134) - g(3) * t130, 0, 0, 0, 0, 0, 0, -g(1) * t83 - g(2) * t85 - g(3) * t94, t128, -t127, -g(1) * t124 - g(2) * t123 - g(3) * t125, 0, 0, 0, 0, 0, 0, t71, t129, -t128, -g(1) * t121 - g(2) * t120 - g(3) * t122, 0, 0, 0, 0, 0, 0, t71, -t128, -t129, -g(1) * (t74 * pkin(5) + t73 * qJ(6) + t121) - g(2) * (t76 * pkin(5) + t75 * qJ(6) + t120) - g(3) * (t81 * pkin(5) + t80 * qJ(6) + t122);];
U_reg  = t1;
