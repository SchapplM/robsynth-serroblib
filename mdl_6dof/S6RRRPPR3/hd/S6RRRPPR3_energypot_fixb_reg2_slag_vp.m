% Calculate inertial parameters regressor of potential energy for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:30:06
% EndTime: 2019-03-09 15:30:06
% DurationCPUTime: 0.11s
% Computational Cost: add. (165->65), mult. (179->73), div. (0->0), fcn. (167->8), ass. (0->38)
t88 = qJ(2) + qJ(3);
t83 = sin(t88);
t104 = qJ(4) * t83;
t84 = cos(t88);
t94 = cos(qJ(1));
t110 = t84 * t94;
t116 = pkin(3) * t110 + t94 * t104;
t115 = g(3) * pkin(6);
t114 = pkin(5) * t83;
t113 = g(3) * t84;
t90 = sin(qJ(2));
t112 = t90 * pkin(2) + pkin(6);
t91 = sin(qJ(1));
t111 = t84 * t91;
t89 = sin(qJ(6));
t109 = t91 * t89;
t92 = cos(qJ(6));
t108 = t91 * t92;
t107 = t94 * t89;
t106 = t94 * t92;
t93 = cos(qJ(2));
t81 = t93 * pkin(2) + pkin(1);
t95 = -pkin(8) - pkin(7);
t105 = t91 * t81 + t94 * t95;
t103 = -qJ(5) - t95;
t102 = t83 * pkin(3) + t112;
t73 = t94 * t81;
t101 = -t91 * t95 + t73;
t100 = pkin(4) * t110 + t116 + t73;
t99 = pkin(3) * t111 + t91 * t104 + t105;
t98 = g(1) * t94 + g(2) * t91;
t97 = -t84 * qJ(4) + t102;
t96 = g(2) * (pkin(4) * t111 + t94 * qJ(5) + t99);
t78 = t83 * pkin(4);
t69 = g(1) * t91 - g(2) * t94;
t68 = g(3) * t83 + t98 * t84;
t67 = t98 * t83 - t113;
t1 = [0, 0, 0, 0, 0, 0, -t98, t69, -g(3), -t115, 0, 0, 0, 0, 0, 0, -g(3) * t90 - t98 * t93, -g(3) * t93 + t98 * t90, -t69, -g(1) * (t94 * pkin(1) + t91 * pkin(7)) - g(2) * (t91 * pkin(1) - t94 * pkin(7)) - t115, 0, 0, 0, 0, 0, 0, -t68, t67, -t69, -g(1) * t101 - g(2) * t105 - g(3) * t112, 0, 0, 0, 0, 0, 0, -t68, -t69, -t67, -g(1) * (t101 + t116) - g(2) * t99 - g(3) * t97, 0, 0, 0, 0, 0, 0, -t67, t68, t69, -g(1) * (t103 * t91 + t100) - t96 - g(3) * (t78 + t97) 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t106 - t109) - g(2) * (t83 * t108 + t107) + t92 * t113, -g(1) * (-t83 * t107 - t108) - g(2) * (-t83 * t109 + t106) - t89 * t113, -t68, -g(1) * (pkin(9) * t110 + t94 * t114 + t100) - t96 - g(3) * (t83 * pkin(9) + t78 + (-pkin(5) - qJ(4)) * t84 + t102) + (-g(1) * t103 - g(2) * (pkin(9) * t84 + t114)) * t91;];
U_reg  = t1;
