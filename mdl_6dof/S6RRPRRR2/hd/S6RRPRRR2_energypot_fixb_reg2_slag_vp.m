% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:19:08
% EndTime: 2019-03-09 13:19:08
% DurationCPUTime: 0.15s
% Computational Cost: add. (211->72), mult. (174->92), div. (0->0), fcn. (166->12), ass. (0->42)
t102 = -pkin(10) - pkin(9);
t93 = qJ(2) + pkin(11);
t86 = qJ(4) + t93;
t80 = sin(t86);
t81 = cos(t86);
t99 = cos(qJ(5));
t82 = t99 * pkin(5) + pkin(4);
t121 = -t102 * t80 + t81 * t82;
t120 = g(3) * pkin(6);
t119 = g(3) * t80;
t97 = sin(qJ(2));
t118 = t97 * pkin(2) + pkin(6);
t100 = cos(qJ(2));
t83 = t100 * pkin(2) + pkin(1);
t94 = qJ(5) + qJ(6);
t87 = sin(t94);
t98 = sin(qJ(1));
t116 = t98 * t87;
t88 = cos(t94);
t115 = t98 * t88;
t96 = sin(qJ(5));
t114 = t98 * t96;
t113 = t98 * t99;
t95 = -pkin(7) - qJ(3);
t101 = cos(qJ(1));
t85 = cos(t93);
t75 = pkin(3) * t85 + t83;
t92 = -pkin(8) + t95;
t112 = t101 * t92 + t98 * t75;
t111 = t101 * t87;
t110 = t101 * t88;
t109 = t101 * t96;
t108 = t101 * t99;
t84 = sin(t93);
t106 = pkin(3) * t84 + t118;
t74 = t101 * t75;
t105 = -t98 * t92 + t74;
t104 = pkin(4) * t81 + pkin(9) * t80;
t103 = g(1) * t101 + g(2) * t98;
t76 = g(1) * t98 - g(2) * t101;
t72 = -g(3) * t81 + t103 * t80;
t1 = [0, 0, 0, 0, 0, 0, -t103, t76, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t97 - t103 * t100, -g(3) * t100 + t103 * t97, -t76, -g(1) * (t101 * pkin(1) + t98 * pkin(7)) - g(2) * (t98 * pkin(1) - t101 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -g(3) * t84 - t103 * t85, -g(3) * t85 + t103 * t84, -t76, -g(1) * (t101 * t83 - t98 * t95) - g(2) * (t101 * t95 + t98 * t83) - g(3) * t118, 0, 0, 0, 0, 0, 0, -t103 * t81 - t119, t72, -t76, -g(1) * t105 - g(2) * t112 - g(3) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t108 + t114) - g(2) * (t81 * t113 - t109) - t99 * t119, -g(1) * (-t81 * t109 + t113) - g(2) * (-t81 * t114 - t108) + t96 * t119, -t72, -g(1) * (t104 * t101 + t105) - g(2) * (t104 * t98 + t112) - g(3) * (t80 * pkin(4) - t81 * pkin(9) + t106) 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t110 + t116) - g(2) * (t81 * t115 - t111) - t88 * t119, -g(1) * (-t81 * t111 + t115) - g(2) * (-t81 * t116 - t110) + t87 * t119, -t72, -g(1) * (t121 * t101 + t74) - g(2) * (-pkin(5) * t109 + t112) - g(3) * (t81 * t102 + t80 * t82 + t106) + (-g(1) * (pkin(5) * t96 - t92) - g(2) * t121) * t98;];
U_reg  = t1;
