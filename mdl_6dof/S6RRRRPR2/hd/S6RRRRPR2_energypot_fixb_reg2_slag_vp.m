% Calculate inertial parameters regressor of potential energy for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:59:13
% EndTime: 2019-03-09 21:59:13
% DurationCPUTime: 0.17s
% Computational Cost: add. (211->72), mult. (174->92), div. (0->0), fcn. (166->12), ass. (0->42)
t96 = cos(pkin(11));
t79 = t96 * pkin(5) + pkin(4);
t94 = qJ(2) + qJ(3);
t88 = qJ(4) + t94;
t81 = sin(t88);
t82 = cos(t88);
t97 = -pkin(10) - qJ(5);
t121 = t79 * t82 - t81 * t97;
t120 = g(3) * pkin(6);
t102 = -pkin(8) - pkin(7);
t119 = g(3) * t81;
t98 = sin(qJ(2));
t118 = t98 * pkin(2) + pkin(6);
t100 = cos(qJ(2));
t83 = t100 * pkin(2) + pkin(1);
t92 = pkin(11) + qJ(6);
t84 = sin(t92);
t99 = sin(qJ(1));
t115 = t99 * t84;
t85 = cos(t92);
t114 = t99 * t85;
t95 = sin(pkin(11));
t113 = t99 * t95;
t112 = t99 * t96;
t101 = cos(qJ(1));
t87 = cos(t94);
t75 = pkin(3) * t87 + t83;
t93 = -pkin(9) + t102;
t111 = t101 * t93 + t99 * t75;
t110 = t101 * t84;
t109 = t101 * t85;
t108 = t101 * t95;
t107 = t101 * t96;
t86 = sin(t94);
t106 = pkin(3) * t86 + t118;
t74 = t101 * t75;
t105 = -t99 * t93 + t74;
t104 = g(1) * t101 + g(2) * t99;
t103 = pkin(4) * t82 + qJ(5) * t81;
t76 = g(1) * t99 - g(2) * t101;
t72 = -g(3) * t82 + t104 * t81;
t1 = [0, 0, 0, 0, 0, 0, -t104, t76, -g(3), -t120, 0, 0, 0, 0, 0, 0, -g(3) * t98 - t104 * t100, -g(3) * t100 + t104 * t98, -t76, -g(1) * (t101 * pkin(1) + t99 * pkin(7)) - g(2) * (t99 * pkin(1) - t101 * pkin(7)) - t120, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t104 * t87, -g(3) * t87 + t104 * t86, -t76, -g(1) * (t101 * t83 - t99 * t102) - g(2) * (t101 * t102 + t99 * t83) - g(3) * t118, 0, 0, 0, 0, 0, 0, -t104 * t82 - t119, t72, -t76, -g(1) * t105 - g(2) * t111 - g(3) * t106, 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t107 + t113) - g(2) * (t82 * t112 - t108) - t96 * t119, -g(1) * (-t82 * t108 + t112) - g(2) * (-t82 * t113 - t107) + t95 * t119, -t72, -g(1) * (t103 * t101 + t105) - g(2) * (t103 * t99 + t111) - g(3) * (t81 * pkin(4) - t82 * qJ(5) + t106) 0, 0, 0, 0, 0, 0, -g(1) * (t82 * t109 + t115) - g(2) * (t82 * t114 - t110) - t85 * t119, -g(1) * (-t82 * t110 + t114) - g(2) * (-t82 * t115 - t109) + t84 * t119, -t72, -g(1) * (t121 * t101 + t74) - g(2) * (-pkin(5) * t108 + t111) - g(3) * (t81 * t79 + t82 * t97 + t106) + (-g(1) * (pkin(5) * t95 - t93) - g(2) * t121) * t99;];
U_reg  = t1;
