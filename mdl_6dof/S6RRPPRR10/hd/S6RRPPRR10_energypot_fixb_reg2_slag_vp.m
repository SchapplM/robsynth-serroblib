% Calculate inertial parameters regressor of potential energy for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPRR10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:00
% EndTime: 2019-03-09 09:37:00
% DurationCPUTime: 0.19s
% Computational Cost: add. (162->90), mult. (225->107), div. (0->0), fcn. (221->10), ass. (0->50)
t94 = sin(qJ(2));
t105 = qJ(3) * t94;
t95 = sin(qJ(1));
t96 = cos(qJ(2));
t113 = t95 * t96;
t125 = pkin(2) * t113 + t105 * t95;
t124 = g(3) * pkin(6);
t123 = g(3) * t96;
t91 = sin(pkin(10));
t122 = t91 * pkin(4);
t121 = pkin(2) * t94 + pkin(6);
t92 = cos(pkin(10));
t79 = pkin(4) * t92 + pkin(3);
t90 = pkin(10) + qJ(5);
t80 = sin(t90);
t71 = pkin(5) * t80 + t122;
t120 = t71 * t94;
t82 = qJ(6) + t90;
t77 = sin(t82);
t119 = t95 * t77;
t78 = cos(t82);
t118 = t95 * t78;
t117 = t95 * t80;
t81 = cos(t90);
t116 = t95 * t81;
t115 = t95 * t91;
t114 = t95 * t92;
t97 = cos(qJ(1));
t112 = t97 * t77;
t111 = t97 * t78;
t110 = t97 * t80;
t109 = t97 * t81;
t108 = t97 * t91;
t107 = t97 * t92;
t93 = -pkin(8) - qJ(4);
t106 = pkin(1) * t97 + pkin(7) * t95;
t104 = qJ(4) * t96;
t103 = t94 * t115;
t86 = t95 * pkin(1);
t102 = t86 + t125;
t101 = -t97 * pkin(7) + t86;
t100 = t106 + (pkin(2) * t96 + t105) * t97;
t99 = -qJ(3) * t96 + t121;
t98 = g(1) * t97 + g(2) * t95;
t89 = -pkin(9) + t93;
t72 = g(1) * t95 - g(2) * t97;
t70 = pkin(5) * t81 + t79;
t69 = g(3) * t94 + t96 * t98;
t68 = t94 * t98 - t123;
t1 = [0, 0, 0, 0, 0, 0, -t98, t72, -g(3), -t124, 0, 0, 0, 0, 0, 0, -t69, t68, -t72, -g(1) * t106 - g(2) * t101 - t124, 0, 0, 0, 0, 0, 0, -t72, t69, -t68, -g(1) * t100 - g(2) * (t101 + t125) - g(3) * t99, 0, 0, 0, 0, 0, 0, -g(1) * (t108 * t94 + t114) - g(2) * (t103 - t107) + t91 * t123, -g(1) * (t107 * t94 - t115) - g(2) * (t114 * t94 + t108) + t92 * t123, -t69, -g(1) * (pkin(3) * t95 + t104 * t97 + t100) - g(2) * (t95 * t104 + (-pkin(3) - pkin(7)) * t97 + t102) - g(3) * (qJ(4) * t94 + t99) 0, 0, 0, 0, 0, 0, -g(1) * (t110 * t94 + t116) - g(2) * (t117 * t94 - t109) + t80 * t123, -g(1) * (t109 * t94 - t117) - g(2) * (t116 * t94 + t110) + t81 * t123, -t69, -g(1) * (t95 * t79 + t100) - g(2) * (pkin(4) * t103 - t93 * t113 + t102) - g(3) * (-t94 * t93 + (-qJ(3) - t122) * t96 + t121) + (-g(1) * (t122 * t94 - t93 * t96) - g(2) * (-pkin(7) - t79)) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t112 * t94 + t118) - g(2) * (t119 * t94 - t111) + t77 * t123, -g(1) * (t111 * t94 - t119) - g(2) * (t118 * t94 + t112) + t78 * t123, -t69, -g(1) * (t95 * t70 + t100) - g(2) * (-t89 * t113 + t95 * t120 + t102) - g(3) * (-t94 * t89 + (-qJ(3) - t71) * t96 + t121) + (-g(1) * (-t89 * t96 + t120) - g(2) * (-pkin(7) - t70)) * t97;];
U_reg  = t1;
