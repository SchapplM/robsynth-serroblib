% Calculate inertial parameters regressor of potential energy for
% S6RRRRRP3
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
% Datum: 2019-03-10 01:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRRRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:09:29
% EndTime: 2019-03-10 01:09:29
% DurationCPUTime: 0.14s
% Computational Cost: add. (202->77), mult. (200->96), div. (0->0), fcn. (196->10), ass. (0->44)
t117 = g(3) * pkin(6);
t97 = -pkin(10) - pkin(9);
t90 = qJ(2) + qJ(3);
t82 = sin(t90);
t116 = g(3) * t82;
t91 = sin(qJ(4));
t115 = t91 * pkin(4);
t92 = sin(qJ(2));
t114 = t92 * pkin(2) + pkin(6);
t94 = cos(qJ(4));
t78 = t94 * pkin(4) + pkin(3);
t88 = -qJ(6) + t97;
t113 = t82 * t88;
t112 = t82 * t97;
t84 = cos(t90);
t96 = cos(qJ(1));
t111 = t84 * t96;
t89 = qJ(4) + qJ(5);
t81 = sin(t89);
t93 = sin(qJ(1));
t110 = t93 * t81;
t83 = cos(t89);
t109 = t93 * t83;
t108 = t93 * t91;
t107 = t93 * t94;
t106 = t96 * t81;
t105 = t96 * t83;
t104 = t96 * t91;
t103 = t96 * t94;
t95 = cos(qJ(2));
t79 = t95 * pkin(2) + pkin(1);
t98 = -pkin(8) - pkin(7);
t102 = t93 * t79 + t96 * t98;
t76 = t96 * t79;
t101 = -t93 * t98 + t76;
t100 = pkin(3) * t84 + pkin(9) * t82;
t99 = g(1) * t96 + g(2) * t93;
t74 = g(1) * t93 - g(2) * t96;
t73 = pkin(5) * t81 + t115;
t72 = pkin(5) * t83 + t78;
t71 = -g(3) * t84 + t99 * t82;
t70 = -g(1) * (t84 * t105 + t110) - g(2) * (t84 * t109 - t106) - t83 * t116;
t69 = -g(1) * (-t84 * t106 + t109) - g(2) * (-t84 * t110 - t105) + t81 * t116;
t1 = [0, 0, 0, 0, 0, 0, -t99, t74, -g(3), -t117, 0, 0, 0, 0, 0, 0, -g(3) * t92 - t99 * t95, -g(3) * t95 + t99 * t92, -t74, -g(1) * (t96 * pkin(1) + t93 * pkin(7)) - g(2) * (t93 * pkin(1) - t96 * pkin(7)) - t117, 0, 0, 0, 0, 0, 0, -t99 * t84 - t116, t71, -t74, -g(1) * t101 - g(2) * t102 - g(3) * t114, 0, 0, 0, 0, 0, 0, -g(1) * (t84 * t103 + t108) - g(2) * (t84 * t107 - t104) - t94 * t116, -g(1) * (-t84 * t104 + t107) - g(2) * (-t84 * t108 - t103) + t91 * t116, -t71, -g(1) * (t100 * t96 + t101) - g(2) * (t100 * t93 + t102) - g(3) * (t82 * pkin(3) - t84 * pkin(9) + t114) 0, 0, 0, 0, 0, 0, t70, t69, -t71, -g(1) * (t78 * t111 - t96 * t112 + t76) - g(2) * (-pkin(4) * t104 + t102) - g(3) * (t82 * t78 + t84 * t97 + t114) + (-g(1) * (-t98 + t115) - g(2) * (t78 * t84 - t112)) * t93, 0, 0, 0, 0, 0, 0, t70, t69, -t71, -g(1) * (t72 * t111 - t96 * t113 + t76) - g(2) * (-t96 * t73 + t102) - g(3) * (t82 * t72 + t84 * t88 + t114) + (-g(1) * (t73 - t98) - g(2) * (t72 * t84 - t113)) * t93;];
U_reg  = t1;
