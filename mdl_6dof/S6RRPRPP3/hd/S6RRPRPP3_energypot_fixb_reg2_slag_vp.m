% Calculate inertial parameters regressor of potential energy for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:57:08
% EndTime: 2019-03-09 09:57:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (200->77), mult. (268->90), div. (0->0), fcn. (278->8), ass. (0->43)
t93 = sin(qJ(1));
t94 = cos(qJ(2));
t111 = t93 * t94;
t88 = pkin(9) + qJ(4);
t82 = sin(t88);
t83 = cos(t88);
t95 = cos(qJ(1));
t66 = t82 * t111 + t83 * t95;
t109 = t95 * t82;
t67 = t83 * t111 - t109;
t122 = t67 * pkin(4) + t66 * qJ(5);
t101 = g(1) * t95 + g(2) * t93;
t121 = g(3) * pkin(6);
t92 = sin(qJ(2));
t118 = g(3) * t92;
t116 = t82 * t92;
t115 = t83 * t92;
t89 = sin(pkin(9));
t114 = t89 * t95;
t91 = -pkin(8) - qJ(3);
t113 = t91 * t92;
t112 = t93 * t89;
t110 = t94 * t95;
t108 = t95 * pkin(1) + t93 * pkin(7);
t90 = cos(pkin(9));
t80 = pkin(3) * t90 + pkin(2);
t106 = t92 * t80 + t94 * t91 + pkin(6);
t105 = t95 * t113;
t85 = t93 * pkin(1);
t104 = -pkin(7) * t95 + t85;
t103 = pkin(3) * t112 + t80 * t110 + t108;
t102 = pkin(4) * t115 + qJ(5) * t116 + t106;
t100 = pkin(2) * t94 + qJ(3) * t92;
t68 = t94 * t109 - t93 * t83;
t69 = t83 * t110 + t93 * t82;
t99 = t69 * pkin(4) + t68 * qJ(5) + t103;
t98 = g(1) * t68 + g(2) * t66 + g(3) * t116;
t97 = g(1) * t69 + g(2) * t67 + g(3) * t115;
t71 = t80 * t111;
t96 = -t93 * t113 + t71 + t85 + (-pkin(3) * t89 - pkin(7)) * t95;
t75 = g(1) * t93 - g(2) * t95;
t70 = -g(3) * t94 + t101 * t92;
t1 = [0, 0, 0, 0, 0, 0, -t101, t75, -g(3), -t121, 0, 0, 0, 0, 0, 0, -t101 * t94 - t118, t70, -t75, -g(1) * t108 - g(2) * t104 - t121, 0, 0, 0, 0, 0, 0, -g(1) * (t90 * t110 + t112) - g(2) * (t90 * t111 - t114) - t90 * t118, -g(1) * (-t89 * t110 + t93 * t90) - g(2) * (-t89 * t111 - t90 * t95) + t89 * t118, -t70, -g(1) * (t100 * t95 + t108) - g(2) * (t100 * t93 + t104) - g(3) * (pkin(2) * t92 - qJ(3) * t94 + pkin(6)) 0, 0, 0, 0, 0, 0, -t97, t98, -t70, -g(1) * (t103 - t105) - g(2) * t96 - g(3) * t106, 0, 0, 0, 0, 0, 0, -t70, t97, -t98, -g(1) * (t99 - t105) - g(2) * (t96 + t122) - g(3) * t102, 0, 0, 0, 0, 0, 0, -t70, -t98, -t97, -g(1) * (t69 * qJ(6) + t99) - g(2) * (-pkin(3) * t114 + t67 * qJ(6) + t104 + t122 + t71) - g(3) * (-pkin(5) * t94 + t102) + (-g(3) * qJ(6) * t83 - t101 * (pkin(5) - t91)) * t92;];
U_reg  = t1;
