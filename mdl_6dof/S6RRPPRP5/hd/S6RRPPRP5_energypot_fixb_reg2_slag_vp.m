% Calculate inertial parameters regressor of potential energy for
% S6RRPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:43:56
% EndTime: 2019-03-09 08:43:56
% DurationCPUTime: 0.15s
% Computational Cost: add. (159->79), mult. (240->90), div. (0->0), fcn. (240->8), ass. (0->51)
t92 = sin(qJ(2));
t109 = qJ(3) * t92;
t93 = sin(qJ(1));
t94 = cos(qJ(2));
t116 = t93 * t94;
t124 = pkin(2) * t116 + t93 * t109;
t123 = g(3) * pkin(6);
t122 = g(3) * t94;
t121 = t92 * pkin(2) + pkin(6);
t88 = pkin(9) + qJ(5);
t81 = sin(t88);
t120 = t81 * t93;
t82 = cos(t88);
t119 = t93 * t82;
t89 = sin(pkin(9));
t118 = t93 * t89;
t90 = cos(pkin(9));
t117 = t93 * t90;
t95 = cos(qJ(1));
t115 = t94 * t95;
t114 = t95 * t81;
t113 = t95 * t82;
t112 = t95 * t89;
t111 = t95 * t90;
t110 = t95 * pkin(1) + t93 * pkin(7);
t108 = qJ(4) * t94;
t107 = t92 * t118;
t106 = t92 * t112;
t85 = t93 * pkin(1);
t105 = t85 + t124;
t104 = -t95 * pkin(7) + t85;
t103 = -pkin(4) * t89 - qJ(3);
t102 = pkin(2) * t115 + t95 * t109 + t110;
t91 = -pkin(8) - qJ(4);
t101 = -t92 * t91 + t121;
t100 = -t94 * qJ(3) + t121;
t99 = g(1) * t95 + g(2) * t93;
t98 = t104 + t124;
t80 = pkin(4) * t90 + pkin(3);
t97 = pkin(4) * t106 + t93 * t80 + t102;
t65 = -t92 * t113 + t120;
t67 = t92 * t119 + t114;
t96 = g(1) * t65 - g(2) * t67 + t82 * t122;
t75 = g(1) * t93 - g(2) * t95;
t72 = pkin(4) * t107;
t70 = g(3) * t92 + t99 * t94;
t69 = t99 * t92 - t122;
t68 = t92 * t120 - t113;
t66 = t92 * t114 + t119;
t64 = -g(1) * t66 - g(2) * t68 + t81 * t122;
t1 = [0, 0, 0, 0, 0, 0, -t99, t75, -g(3), -t123, 0, 0, 0, 0, 0, 0, -t70, t69, -t75, -g(1) * t110 - g(2) * t104 - t123, 0, 0, 0, 0, 0, 0, -t75, t70, -t69, -g(1) * t102 - g(2) * t98 - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (t106 + t117) - g(2) * (t107 - t111) + t89 * t122, -g(1) * (t92 * t111 - t118) - g(2) * (t92 * t117 + t112) + t90 * t122, -t70, -g(1) * (t93 * pkin(3) + t95 * t108 + t102) - g(2) * (t93 * t108 + (-pkin(3) - pkin(7)) * t95 + t105) - g(3) * (t92 * qJ(4) + t100) 0, 0, 0, 0, 0, 0, t64, t96, -t70, -g(1) * (-t91 * t115 + t97) - g(2) * (-t91 * t116 + t72 + (-pkin(7) - t80) * t95 + t105) - g(3) * (t103 * t94 + t101) 0, 0, 0, 0, 0, 0, t64, -t70, -t96, -g(1) * (t66 * pkin(5) + t65 * qJ(6) + t97) - g(2) * (t68 * pkin(5) - t67 * qJ(6) - t95 * t80 + t72 + t98) - g(3) * t101 + (-g(3) * (-pkin(5) * t81 + qJ(6) * t82 + t103) + t99 * t91) * t94;];
U_reg  = t1;
