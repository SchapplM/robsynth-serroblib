% Calculate inertial parameters regressor of potential energy for
% S6RRRRRR2
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:35:46
% EndTime: 2019-03-10 03:35:47
% DurationCPUTime: 0.15s
% Computational Cost: add. (211->72), mult. (174->92), div. (0->0), fcn. (166->12), ass. (0->42)
t101 = -pkin(11) - pkin(10);
t94 = qJ(2) + qJ(3);
t88 = qJ(4) + t94;
t80 = sin(t88);
t81 = cos(t88);
t98 = cos(qJ(5));
t82 = t98 * pkin(5) + pkin(4);
t121 = -t101 * t80 + t81 * t82;
t120 = g(3) * pkin(6);
t102 = -pkin(8) - pkin(7);
t119 = g(3) * t80;
t96 = sin(qJ(2));
t118 = t96 * pkin(2) + pkin(6);
t99 = cos(qJ(2));
t83 = t99 * pkin(2) + pkin(1);
t93 = qJ(5) + qJ(6);
t84 = sin(t93);
t97 = sin(qJ(1));
t116 = t97 * t84;
t86 = cos(t93);
t115 = t97 * t86;
t95 = sin(qJ(5));
t114 = t97 * t95;
t113 = t97 * t98;
t100 = cos(qJ(1));
t87 = cos(t94);
t75 = pkin(3) * t87 + t83;
t92 = -pkin(9) + t102;
t112 = t100 * t92 + t97 * t75;
t111 = t100 * t84;
t110 = t100 * t86;
t109 = t100 * t95;
t108 = t100 * t98;
t85 = sin(t94);
t106 = pkin(3) * t85 + t118;
t74 = t100 * t75;
t105 = -t97 * t92 + t74;
t104 = pkin(4) * t81 + pkin(10) * t80;
t103 = g(1) * t100 + g(2) * t97;
t76 = g(1) * t97 - g(2) * t100;
t72 = -g(3) * t81 + t103 * t80;
t1 = [0, 0, 0, 0, 0, 0, -t103, t76, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t96 - t103 * t99, -g(3) * t99 + t103 * t96, -t76, -g(1) * (t100 * pkin(1) + t97 * pkin(7)) - g(2) * (t97 * pkin(1) - t100 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -g(3) * t85 - t103 * t87, -g(3) * t87 + t103 * t85, -t76, -g(1) * (t100 * t83 - t97 * t102) - g(2) * (t100 * t102 + t97 * t83) - g(3) * t118, 0, 0, 0, 0, 0, 0, -t103 * t81 - t119, t72, -t76, -g(1) * t105 - g(2) * t112 - g(3) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t108 + t114) - g(2) * (t81 * t113 - t109) - t98 * t119, -g(1) * (-t81 * t109 + t113) - g(2) * (-t81 * t114 - t108) + t95 * t119, -t72, -g(1) * (t104 * t100 + t105) - g(2) * (t104 * t97 + t112) - g(3) * (t80 * pkin(4) - t81 * pkin(10) + t106) 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t110 + t116) - g(2) * (t81 * t115 - t111) - t86 * t119, -g(1) * (-t81 * t111 + t115) - g(2) * (-t81 * t116 - t110) + t84 * t119, -t72, -g(1) * (t121 * t100 + t74) - g(2) * (-pkin(5) * t109 + t112) - g(3) * (t81 * t101 + t80 * t82 + t106) + (-g(1) * (pkin(5) * t95 - t92) - g(2) * t121) * t97;];
U_reg  = t1;
