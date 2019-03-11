% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:31:25
% EndTime: 2019-03-10 01:31:26
% DurationCPUTime: 0.17s
% Computational Cost: add. (223->85), mult. (243->107), div. (0->0), fcn. (247->10), ass. (0->43)
t122 = g(3) * pkin(6);
t105 = -pkin(9) - pkin(8);
t99 = sin(qJ(3));
t121 = t99 * pkin(3);
t102 = cos(qJ(3));
t88 = t102 * pkin(3) + pkin(2);
t100 = sin(qJ(2));
t120 = g(3) * t100;
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t119 = t104 * pkin(1) + t101 * pkin(7);
t97 = -pkin(10) + t105;
t118 = t100 * t97;
t117 = t101 * t99;
t116 = t100 * t105;
t103 = cos(qJ(2));
t115 = t101 * t103;
t114 = t103 * t104;
t113 = t104 * t102;
t98 = qJ(3) + qJ(4);
t90 = cos(t98);
t81 = pkin(4) * t90 + t88;
t112 = t100 * t81 + t103 * t97 + pkin(6);
t93 = t101 * pkin(1);
t111 = -t104 * pkin(7) + t93;
t110 = pkin(2) * t103 + pkin(8) * t100;
t109 = g(1) * t104 + g(2) * t101;
t89 = sin(t98);
t82 = pkin(4) * t89 + t121;
t108 = t101 * t82 - t104 * t118 + t81 * t114 + t119;
t91 = qJ(5) + t98;
t86 = sin(t91);
t87 = cos(t91);
t72 = t104 * t87 + t86 * t115;
t74 = -t101 * t87 + t86 * t114;
t107 = g(1) * t74 + g(2) * t72 + t86 * t120;
t106 = -t101 * t118 + t81 * t115 + t93 + (-pkin(7) - t82) * t104;
t83 = g(1) * t101 - g(2) * t104;
t76 = -g(3) * t103 + t109 * t100;
t75 = t101 * t86 + t87 * t114;
t73 = -t104 * t86 + t87 * t115;
t71 = -g(1) * t75 - g(2) * t73 - t87 * t120;
t1 = [0, 0, 0, 0, 0, 0, -t109, t83, -g(3), -t122, 0, 0, 0, 0, 0, 0, -t109 * t103 - t120, t76, -t83, -g(1) * t119 - g(2) * t111 - t122, 0, 0, 0, 0, 0, 0, -g(1) * (t103 * t113 + t117) - g(2) * (t102 * t115 - t104 * t99) - t102 * t120, -g(1) * (t101 * t102 - t99 * t114) - g(2) * (-t99 * t115 - t113) + t99 * t120, -t76, -g(1) * (t110 * t104 + t119) - g(2) * (t110 * t101 + t111) - g(3) * (t100 * pkin(2) - t103 * pkin(8) + pkin(6)) 0, 0, 0, 0, 0, 0, -g(1) * (t101 * t89 + t90 * t114) - g(2) * (-t104 * t89 + t90 * t115) - t90 * t120, -g(1) * (t101 * t90 - t89 * t114) - g(2) * (-t104 * t90 - t89 * t115) + t89 * t120, -t76, -g(1) * (pkin(3) * t117 + t119) - g(2) * (-t101 * t116 + t88 * t115 + t93) - g(3) * (t100 * t88 + t103 * t105 + pkin(6)) + (-g(1) * (t103 * t88 - t116) - g(2) * (-pkin(7) - t121)) * t104, 0, 0, 0, 0, 0, 0, t71, t107, -t76, -g(1) * t108 - g(2) * t106 - g(3) * t112, 0, 0, 0, 0, 0, 0, t71, -t76, -t107, -g(1) * (t75 * pkin(5) + t74 * qJ(6) + t108) - g(2) * (t73 * pkin(5) + t72 * qJ(6) + t106) - g(3) * ((pkin(5) * t87 + qJ(6) * t86) * t100 + t112);];
U_reg  = t1;
