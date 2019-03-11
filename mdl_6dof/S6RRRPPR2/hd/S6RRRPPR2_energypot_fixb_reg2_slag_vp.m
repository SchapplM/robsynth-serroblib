% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:26:42
% EndTime: 2019-03-09 15:26:42
% DurationCPUTime: 0.11s
% Computational Cost: add. (193->63), mult. (161->75), div. (0->0), fcn. (149->10), ass. (0->38)
t92 = qJ(2) + qJ(3);
t85 = pkin(10) + t92;
t81 = sin(t85);
t105 = qJ(5) * t81;
t82 = cos(t85);
t98 = cos(qJ(1));
t111 = t82 * t98;
t116 = pkin(4) * t111 + t98 * t105;
t115 = g(3) * pkin(6);
t99 = -pkin(8) - pkin(7);
t114 = g(3) * t82;
t94 = sin(qJ(2));
t113 = t94 * pkin(2) + pkin(6);
t97 = cos(qJ(2));
t84 = t97 * pkin(2) + pkin(1);
t95 = sin(qJ(1));
t112 = t82 * t95;
t93 = sin(qJ(6));
t110 = t95 * t93;
t96 = cos(qJ(6));
t109 = t95 * t96;
t108 = t98 * t93;
t107 = t98 * t96;
t87 = cos(t92);
t73 = pkin(3) * t87 + t84;
t91 = -qJ(4) + t99;
t106 = t95 * t73 + t98 * t91;
t86 = sin(t92);
t104 = pkin(3) * t86 + t113;
t72 = t98 * t73;
t103 = -t95 * t91 + t72;
t102 = pkin(4) * t112 + t95 * t105 + t106;
t101 = g(1) * t98 + g(2) * t95;
t100 = t81 * pkin(4) - t82 * qJ(5) + t104;
t78 = g(1) * t95 - g(2) * t98;
t70 = g(3) * t81 + t101 * t82;
t69 = t101 * t81 - t114;
t1 = [0, 0, 0, 0, 0, 0, -t101, t78, -g(3), -t115, 0, 0, 0, 0, 0, 0, -g(3) * t94 - t101 * t97, -g(3) * t97 + t101 * t94, -t78, -g(1) * (t98 * pkin(1) + t95 * pkin(7)) - g(2) * (t95 * pkin(1) - t98 * pkin(7)) - t115, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t101 * t87, -g(3) * t87 + t101 * t86, -t78, -g(1) * (t98 * t84 - t95 * t99) - g(2) * (t95 * t84 + t98 * t99) - g(3) * t113, 0, 0, 0, 0, 0, 0, -t70, t69, -t78, -g(1) * t103 - g(2) * t106 - g(3) * t104, 0, 0, 0, 0, 0, 0, -t78, t70, -t69, -g(1) * (t103 + t116) - g(2) * t102 - g(3) * t100, 0, 0, 0, 0, 0, 0, -g(1) * (t81 * t108 + t109) - g(2) * (t81 * t110 - t107) + t93 * t114, -g(1) * (t81 * t107 - t110) - g(2) * (t81 * t109 + t108) + t96 * t114, -t70, -g(1) * (pkin(9) * t111 + t72 + (pkin(5) - t91) * t95 + t116) - g(2) * (-t98 * pkin(5) + pkin(9) * t112 + t102) - g(3) * (t81 * pkin(9) + t100);];
U_reg  = t1;
