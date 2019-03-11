% Calculate inertial parameters regressor of potential energy for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPP4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:02:34
% EndTime: 2019-03-09 21:02:34
% DurationCPUTime: 0.17s
% Computational Cost: add. (223->85), mult. (243->107), div. (0->0), fcn. (247->10), ass. (0->44)
t121 = g(3) * pkin(6);
t103 = -pkin(9) - pkin(8);
t98 = sin(qJ(2));
t120 = g(3) * t98;
t97 = sin(qJ(3));
t119 = t97 * pkin(3);
t100 = cos(qJ(3));
t86 = t100 * pkin(3) + pkin(2);
t95 = -qJ(5) + t103;
t118 = t95 * t98;
t99 = sin(qJ(1));
t117 = t99 * t97;
t102 = cos(qJ(1));
t116 = t102 * pkin(1) + t99 * pkin(7);
t101 = cos(qJ(2));
t115 = t101 * t99;
t114 = t103 * t98;
t113 = t99 * t100;
t112 = t101 * t102;
t111 = t102 * t100;
t96 = qJ(3) + qJ(4);
t89 = cos(t96);
t79 = pkin(4) * t89 + t86;
t110 = t101 * t95 + t98 * t79 + pkin(6);
t91 = t99 * pkin(1);
t109 = -t102 * pkin(7) + t91;
t108 = pkin(2) * t101 + pkin(8) * t98;
t107 = g(1) * t102 + g(2) * t99;
t88 = sin(t96);
t80 = pkin(4) * t88 + t119;
t106 = -t102 * t118 + t79 * t112 + t99 * t80 + t116;
t87 = pkin(10) + t96;
t84 = sin(t87);
t85 = cos(t87);
t70 = t102 * t85 + t84 * t115;
t72 = t84 * t112 - t99 * t85;
t105 = g(1) * t72 + g(2) * t70 + t84 * t120;
t104 = -t99 * t118 + t79 * t115 + t91 + (-pkin(7) - t80) * t102;
t82 = g(1) * t99 - g(2) * t102;
t74 = -g(3) * t101 + t107 * t98;
t73 = t85 * t112 + t99 * t84;
t71 = -t102 * t84 + t85 * t115;
t69 = -g(1) * t73 - g(2) * t71 - t85 * t120;
t1 = [0, 0, 0, 0, 0, 0, -t107, t82, -g(3), -t121, 0, 0, 0, 0, 0, 0, -t107 * t101 - t120, t74, -t82, -g(1) * t116 - g(2) * t109 - t121, 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t111 + t117) - g(2) * (t101 * t113 - t102 * t97) - t100 * t120, -g(1) * (-t97 * t112 + t113) - g(2) * (-t97 * t115 - t111) + t97 * t120, -t74, -g(1) * (t108 * t102 + t116) - g(2) * (t108 * t99 + t109) - g(3) * (t98 * pkin(2) - t101 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t89 * t112 + t99 * t88) - g(2) * (-t102 * t88 + t89 * t115) - t89 * t120, -g(1) * (-t88 * t112 + t99 * t89) - g(2) * (-t102 * t89 - t88 * t115) + t88 * t120, -t74, -g(1) * (pkin(3) * t117 + t116) - g(2) * (-t99 * t114 + t86 * t115 + t91) - g(3) * (t101 * t103 + t98 * t86 + pkin(6)) + (-g(1) * (t101 * t86 - t114) - g(2) * (-pkin(7) - t119)) * t102, 0, 0, 0, 0, 0, 0, t69, t105, -t74, -g(1) * t106 - g(2) * t104 - g(3) * t110, 0, 0, 0, 0, 0, 0, t69, -t74, -t105, -g(1) * (t73 * pkin(5) + t72 * qJ(6) + t106) - g(2) * (t71 * pkin(5) + t70 * qJ(6) + t104) - g(3) * ((pkin(5) * t85 + qJ(6) * t84) * t98 + t110);];
U_reg  = t1;
