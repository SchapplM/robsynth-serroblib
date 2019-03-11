% Calculate inertial parameters regressor of potential energy for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPRRR11_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:32:48
% EndTime: 2019-03-09 14:32:48
% DurationCPUTime: 0.22s
% Computational Cost: add. (162->90), mult. (225->107), div. (0->0), fcn. (221->10), ass. (0->51)
t93 = sin(qJ(2));
t105 = qJ(3) * t93;
t94 = sin(qJ(1));
t96 = cos(qJ(2));
t115 = t94 * t96;
t127 = pkin(2) * t115 + t94 * t105;
t126 = g(3) * pkin(6);
t98 = -pkin(9) - pkin(8);
t125 = g(3) * t96;
t92 = sin(qJ(4));
t124 = t92 * pkin(4);
t123 = t93 * pkin(2) + pkin(6);
t95 = cos(qJ(4));
t80 = t95 * pkin(4) + pkin(3);
t91 = qJ(4) + qJ(5);
t81 = sin(t91);
t72 = pkin(5) * t81 + t124;
t122 = t72 * t93;
t83 = qJ(6) + t91;
t78 = sin(t83);
t121 = t94 * t78;
t79 = cos(t83);
t120 = t94 * t79;
t119 = t94 * t81;
t82 = cos(t91);
t118 = t94 * t82;
t117 = t94 * t92;
t116 = t94 * t95;
t97 = cos(qJ(1));
t114 = t96 * t97;
t113 = t96 * t98;
t112 = t97 * t78;
t111 = t97 * t79;
t110 = t97 * t81;
t109 = t97 * t82;
t108 = t97 * t92;
t107 = t97 * t95;
t106 = t97 * pkin(1) + t94 * pkin(7);
t104 = t93 * t117;
t86 = t94 * pkin(1);
t103 = t86 + t127;
t102 = -t97 * pkin(7) + t86;
t101 = pkin(2) * t114 + t97 * t105 + t106;
t100 = -t96 * qJ(3) + t123;
t99 = g(1) * t97 + g(2) * t94;
t90 = -pkin(10) + t98;
t73 = g(1) * t94 - g(2) * t97;
t71 = pkin(5) * t82 + t80;
t70 = g(3) * t93 + t99 * t96;
t69 = t99 * t93 - t125;
t1 = [0, 0, 0, 0, 0, 0, -t99, t73, -g(3), -t126, 0, 0, 0, 0, 0, 0, -t70, t69, -t73, -g(1) * t106 - g(2) * t102 - t126, 0, 0, 0, 0, 0, 0, -t73, t70, -t69, -g(1) * t101 - g(2) * (t102 + t127) - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t108 + t116) - g(2) * (t104 - t107) + t92 * t125, -g(1) * (t93 * t107 - t117) - g(2) * (t93 * t116 + t108) + t95 * t125, -t70, -g(1) * (t94 * pkin(3) + pkin(8) * t114 + t101) - g(2) * (pkin(8) * t115 + (-pkin(3) - pkin(7)) * t97 + t103) - g(3) * (t93 * pkin(8) + t100) 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t110 + t118) - g(2) * (t93 * t119 - t109) + t81 * t125, -g(1) * (t93 * t109 - t119) - g(2) * (t93 * t118 + t110) + t82 * t125, -t70, -g(1) * (t94 * t80 + t101) - g(2) * (pkin(4) * t104 - t94 * t113 + t103) - g(3) * (-t93 * t98 + (-qJ(3) - t124) * t96 + t123) + (-g(1) * (t93 * t124 - t113) - g(2) * (-pkin(7) - t80)) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t93 * t112 + t120) - g(2) * (t93 * t121 - t111) + t78 * t125, -g(1) * (t93 * t111 - t121) - g(2) * (t93 * t120 + t112) + t79 * t125, -t70, -g(1) * (t94 * t71 + t101) - g(2) * (-t90 * t115 + t94 * t122 + t103) - g(3) * (-t93 * t90 + (-qJ(3) - t72) * t96 + t123) + (-g(1) * (-t90 * t96 + t122) - g(2) * (-pkin(7) - t71)) * t97;];
U_reg  = t1;
