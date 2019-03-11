% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:01:04
% EndTime: 2019-03-09 16:01:04
% DurationCPUTime: 0.20s
% Computational Cost: add. (182->85), mult. (350->113), div. (0->0), fcn. (388->10), ass. (0->46)
t100 = sin(qJ(2));
t102 = cos(qJ(3));
t121 = t100 * t102;
t99 = sin(qJ(3));
t123 = t100 * t99;
t129 = pkin(3) * t121 + qJ(4) * t123;
t128 = g(3) * pkin(6);
t127 = t100 * pkin(2) + pkin(6);
t126 = g(3) * t100;
t96 = sin(pkin(10));
t125 = t96 * t99;
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t124 = t104 * pkin(1) + t101 * pkin(7);
t122 = t100 * t101;
t120 = t100 * t104;
t103 = cos(qJ(2));
t119 = t101 * t103;
t118 = t103 * t104;
t117 = t104 * t102;
t116 = pkin(5) * t96 + qJ(4);
t115 = t101 * pkin(1) - t104 * pkin(7);
t114 = t127 + t129;
t113 = pkin(2) * t118 + pkin(8) * t120 + t124;
t112 = -t103 * pkin(8) + t127;
t78 = t101 * t99 + t103 * t117;
t111 = t78 * pkin(3) + t113;
t110 = g(1) * t104 + g(2) * t101;
t109 = pkin(2) * t119 + pkin(8) * t122 + t115;
t76 = t102 * t119 - t104 * t99;
t108 = t76 * pkin(3) + t109;
t77 = -t101 * t102 + t99 * t118;
t107 = t77 * qJ(4) + t111;
t75 = t99 * t119 + t117;
t106 = g(1) * t77 + g(2) * t75 + g(3) * t123;
t105 = t75 * qJ(4) + t108;
t98 = -pkin(9) - qJ(5);
t97 = cos(pkin(10));
t95 = pkin(10) + qJ(6);
t89 = cos(t95);
t88 = sin(t95);
t86 = t97 * pkin(5) + pkin(4);
t79 = g(1) * t101 - g(2) * t104;
t72 = -g(3) * t103 + t110 * t100;
t71 = -g(1) * t78 - g(2) * t76 - g(3) * t121;
t1 = [0, 0, 0, 0, 0, 0, -t110, t79, -g(3), -t128, 0, 0, 0, 0, 0, 0, -t110 * t103 - t126, t72, -t79, -g(1) * t124 - g(2) * t115 - t128, 0, 0, 0, 0, 0, 0, t71, t106, -t72, -g(1) * t113 - g(2) * t109 - g(3) * t112, 0, 0, 0, 0, 0, 0, t71, -t72, -t106, -g(1) * t107 - g(2) * t105 - g(3) * (t112 + t129) 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t96 + t78 * t97) - g(2) * (t75 * t96 + t76 * t97) - (t102 * t97 + t125) * t126, -g(1) * (t77 * t97 - t78 * t96) - g(2) * (t75 * t97 - t76 * t96) - (-t102 * t96 + t97 * t99) * t126, t72, -g(1) * (t78 * pkin(4) - qJ(5) * t120 + t107) - g(2) * (t76 * pkin(4) - qJ(5) * t122 + t105) - g(3) * (pkin(4) * t121 + (-pkin(8) + qJ(5)) * t103 + t114) 0, 0, 0, 0, 0, 0, -g(1) * (t77 * t88 + t78 * t89) - g(2) * (t75 * t88 + t76 * t89) - (t102 * t89 + t88 * t99) * t126, -g(1) * (t77 * t89 - t78 * t88) - g(2) * (t75 * t89 - t76 * t88) - (-t102 * t88 + t89 * t99) * t126, t72, -g(1) * (t116 * t77 + t98 * t120 + t78 * t86 + t111) - g(2) * (t116 * t75 + t98 * t122 + t76 * t86 + t108) - g(3) * ((-pkin(8) - t98) * t103 + (pkin(5) * t125 + t102 * t86) * t100 + t114);];
U_reg  = t1;
