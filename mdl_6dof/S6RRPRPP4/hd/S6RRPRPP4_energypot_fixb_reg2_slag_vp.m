% Calculate inertial parameters regressor of potential energy for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:01:39
% EndTime: 2019-03-09 10:01:39
% DurationCPUTime: 0.15s
% Computational Cost: add. (159->79), mult. (240->90), div. (0->0), fcn. (240->8), ass. (0->48)
t91 = sin(qJ(2));
t108 = qJ(3) * t91;
t92 = sin(qJ(1));
t94 = cos(qJ(2));
t112 = t92 * t94;
t121 = pkin(2) * t112 + t92 * t108;
t120 = g(3) * pkin(6);
t119 = g(3) * t94;
t118 = t91 * pkin(2) + pkin(6);
t88 = qJ(4) + pkin(9);
t81 = sin(t88);
t117 = t81 * t92;
t95 = cos(qJ(1));
t116 = t91 * t95;
t82 = cos(t88);
t115 = t92 * t82;
t90 = sin(qJ(4));
t114 = t92 * t90;
t93 = cos(qJ(4));
t113 = t92 * t93;
t111 = t93 * t95;
t110 = t94 * t95;
t109 = t95 * pkin(1) + t92 * pkin(7);
t107 = t90 * t116;
t106 = t91 * t114;
t85 = t92 * pkin(1);
t105 = t85 + t121;
t104 = -pkin(7) * t95 + t85;
t103 = -pkin(4) * t90 - qJ(3);
t102 = pkin(2) * t110 + t95 * t108 + t109;
t89 = -qJ(5) - pkin(8);
t101 = -t89 * t91 + t118;
t100 = -qJ(3) * t94 + t118;
t99 = g(1) * t95 + g(2) * t92;
t98 = t104 + t121;
t80 = pkin(4) * t93 + pkin(3);
t97 = pkin(4) * t107 + t92 * t80 + t102;
t65 = -t82 * t116 + t117;
t67 = t91 * t115 + t81 * t95;
t96 = g(1) * t65 - g(2) * t67 + t82 * t119;
t74 = g(1) * t92 - g(2) * t95;
t72 = pkin(4) * t106;
t70 = g(3) * t91 + t99 * t94;
t69 = t99 * t91 - t119;
t68 = t91 * t117 - t82 * t95;
t66 = t81 * t116 + t115;
t64 = -g(1) * t66 - g(2) * t68 + t81 * t119;
t1 = [0, 0, 0, 0, 0, 0, -t99, t74, -g(3), -t120, 0, 0, 0, 0, 0, 0, -t70, t69, -t74, -g(1) * t109 - g(2) * t104 - t120, 0, 0, 0, 0, 0, 0, -t74, t70, -t69, -g(1) * t102 - g(2) * t98 - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (t107 + t113) - g(2) * (t106 - t111) + t90 * t119, -g(1) * (t91 * t111 - t114) - g(2) * (t91 * t113 + t90 * t95) + t93 * t119, -t70, -g(1) * (t92 * pkin(3) + pkin(8) * t110 + t102) - g(2) * (pkin(8) * t112 + (-pkin(3) - pkin(7)) * t95 + t105) - g(3) * (pkin(8) * t91 + t100) 0, 0, 0, 0, 0, 0, t64, t96, -t70, -g(1) * (-t89 * t110 + t97) - g(2) * (-t89 * t112 + t72 + (-pkin(7) - t80) * t95 + t105) - g(3) * (t103 * t94 + t101) 0, 0, 0, 0, 0, 0, t64, -t70, -t96, -g(1) * (t66 * pkin(5) + t65 * qJ(6) + t97) - g(2) * (t68 * pkin(5) - t67 * qJ(6) - t80 * t95 + t72 + t98) - g(3) * t101 + (-g(3) * (-pkin(5) * t81 + qJ(6) * t82 + t103) + t99 * t89) * t94;];
U_reg  = t1;
