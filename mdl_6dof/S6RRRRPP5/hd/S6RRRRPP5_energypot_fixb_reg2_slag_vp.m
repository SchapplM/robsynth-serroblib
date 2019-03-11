% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:06
% EndTime: 2019-03-09 21:09:06
% DurationCPUTime: 0.15s
% Computational Cost: add. (200->74), mult. (268->90), div. (0->0), fcn. (278->8), ass. (0->44)
t89 = qJ(3) + qJ(4);
t84 = cos(t89);
t95 = cos(qJ(1));
t112 = t95 * t84;
t92 = sin(qJ(1));
t94 = cos(qJ(2));
t114 = t92 * t94;
t83 = sin(t89);
t68 = t83 * t114 + t112;
t113 = t95 * t83;
t69 = t84 * t114 - t113;
t123 = t69 * pkin(4) + t68 * qJ(5);
t100 = g(1) * t95 + g(2) * t92;
t122 = g(3) * pkin(6);
t91 = sin(qJ(2));
t119 = g(3) * t91;
t118 = t83 * t91;
t117 = t84 * t91;
t96 = -pkin(9) - pkin(8);
t116 = t91 * t96;
t90 = sin(qJ(3));
t115 = t92 * t90;
t111 = t95 * t90;
t93 = cos(qJ(3));
t110 = t95 * t93;
t109 = t95 * pkin(1) + t92 * pkin(7);
t81 = t93 * pkin(3) + pkin(2);
t106 = t91 * t81 + t94 * t96 + pkin(6);
t105 = t95 * t116;
t86 = t92 * pkin(1);
t104 = -t95 * pkin(7) + t86;
t103 = t95 * t94 * t81 + pkin(3) * t115 + t109;
t102 = pkin(4) * t117 + qJ(5) * t118 + t106;
t101 = pkin(2) * t94 + pkin(8) * t91;
t70 = t94 * t113 - t92 * t84;
t71 = t94 * t112 + t92 * t83;
t99 = t71 * pkin(4) + t70 * qJ(5) + t103;
t98 = g(1) * t70 + g(2) * t68 + g(3) * t118;
t73 = t81 * t114;
t97 = -t92 * t116 + t73 + t86 + (-pkin(3) * t90 - pkin(7)) * t95;
t75 = g(1) * t92 - g(2) * t95;
t72 = -g(3) * t94 + t100 * t91;
t65 = -g(1) * t71 - g(2) * t69 - g(3) * t117;
t1 = [0, 0, 0, 0, 0, 0, -t100, t75, -g(3), -t122, 0, 0, 0, 0, 0, 0, -t100 * t94 - t119, t72, -t75, -g(1) * t109 - g(2) * t104 - t122, 0, 0, 0, 0, 0, 0, -g(1) * (t94 * t110 + t115) - g(2) * (t93 * t114 - t111) - t93 * t119, -g(1) * (-t94 * t111 + t92 * t93) - g(2) * (-t90 * t114 - t110) + t90 * t119, -t72, -g(1) * (t101 * t95 + t109) - g(2) * (t101 * t92 + t104) - g(3) * (t91 * pkin(2) - t94 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, t65, t98, -t72, -g(1) * (t103 - t105) - g(2) * t97 - g(3) * t106, 0, 0, 0, 0, 0, 0, t65, -t72, -t98, -g(1) * (t99 - t105) - g(2) * (t97 + t123) - g(3) * t102, 0, 0, 0, 0, 0, 0, t65, -t98, t72, -g(1) * (t71 * pkin(5) + t99) - g(2) * (-pkin(3) * t111 + t69 * pkin(5) + t104 + t123 + t73) - g(3) * (t94 * qJ(6) + t102) + (-g(3) * pkin(5) * t84 + t100 * (qJ(6) + t96)) * t91;];
U_reg  = t1;
