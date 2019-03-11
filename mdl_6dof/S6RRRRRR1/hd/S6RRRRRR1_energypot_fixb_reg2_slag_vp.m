% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:29:47
% EndTime: 2019-03-10 03:29:47
% DurationCPUTime: 0.10s
% Computational Cost: add. (207->61), mult. (150->77), div. (0->0), fcn. (138->12), ass. (0->37)
t115 = g(3) * pkin(6);
t102 = -pkin(8) - pkin(7);
t95 = qJ(2) + qJ(3);
t89 = qJ(4) + t95;
t86 = qJ(5) + t89;
t79 = sin(t86);
t114 = g(3) * t79;
t97 = sin(qJ(2));
t113 = t97 * pkin(2) + pkin(6);
t100 = cos(qJ(2));
t85 = t100 * pkin(2) + pkin(1);
t96 = sin(qJ(6));
t98 = sin(qJ(1));
t112 = t98 * t96;
t99 = cos(qJ(6));
t111 = t98 * t99;
t101 = cos(qJ(1));
t88 = cos(t95);
t75 = pkin(3) * t88 + t85;
t84 = cos(t89);
t74 = pkin(4) * t84 + t75;
t94 = -pkin(9) + t102;
t90 = -pkin(10) + t94;
t110 = t101 * t90 + t98 * t74;
t109 = t101 * t96;
t108 = t101 * t99;
t87 = sin(t95);
t107 = pkin(3) * t87 + t113;
t106 = t101 * t74 - t98 * t90;
t83 = sin(t89);
t105 = pkin(4) * t83 + t107;
t80 = cos(t86);
t104 = pkin(5) * t80 + pkin(11) * t79;
t103 = g(1) * t101 + g(2) * t98;
t76 = g(1) * t98 - g(2) * t101;
t71 = -g(3) * t80 + t103 * t79;
t1 = [0, 0, 0, 0, 0, 0, -t103, t76, -g(3), -t115, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t103 * t100, -g(3) * t100 + t103 * t97, -t76, -g(1) * (t101 * pkin(1) + t98 * pkin(7)) - g(2) * (t98 * pkin(1) - t101 * pkin(7)) - t115, 0, 0, 0, 0, 0, 0, -g(3) * t87 - t103 * t88, -g(3) * t88 + t103 * t87, -t76, -g(1) * (t101 * t85 - t98 * t102) - g(2) * (t101 * t102 + t98 * t85) - g(3) * t113, 0, 0, 0, 0, 0, 0, -g(3) * t83 - t103 * t84, -g(3) * t84 + t103 * t83, -t76, -g(1) * (t101 * t75 - t98 * t94) - g(2) * (t101 * t94 + t98 * t75) - g(3) * t107, 0, 0, 0, 0, 0, 0, -t103 * t80 - t114, t71, -t76, -g(1) * t106 - g(2) * t110 - g(3) * t105, 0, 0, 0, 0, 0, 0, -g(1) * (t80 * t108 + t112) - g(2) * (t80 * t111 - t109) - t99 * t114, -g(1) * (-t80 * t109 + t111) - g(2) * (-t80 * t112 - t108) + t96 * t114, -t71, -g(1) * (t104 * t101 + t106) - g(2) * (t104 * t98 + t110) - g(3) * (t79 * pkin(5) - t80 * pkin(11) + t105);];
U_reg  = t1;
