% Calculate inertial parameters regressor of potential energy for
% S6RPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:01:15
% EndTime: 2019-03-09 06:01:15
% DurationCPUTime: 0.14s
% Computational Cost: add. (212->74), mult. (186->95), div. (0->0), fcn. (182->10), ass. (0->40)
t89 = -pkin(9) - pkin(8);
t82 = qJ(2) + pkin(6);
t106 = g(3) * t82;
t84 = sin(qJ(3));
t105 = g(3) * t84;
t83 = sin(qJ(4));
t104 = t83 * pkin(4);
t86 = cos(qJ(4));
t70 = t86 * pkin(4) + pkin(3);
t80 = qJ(1) + pkin(10);
t71 = sin(t80);
t103 = t71 * t83;
t87 = cos(qJ(3));
t102 = t71 * t87;
t81 = qJ(4) + qJ(5);
t73 = sin(t81);
t101 = t73 * t87;
t74 = cos(t81);
t100 = t74 * t87;
t79 = -qJ(6) + t89;
t99 = t79 * t84;
t98 = t83 * t87;
t97 = t84 * t89;
t96 = t86 * t87;
t85 = sin(qJ(1));
t95 = t85 * pkin(1) + t71 * pkin(2);
t72 = cos(t80);
t88 = cos(qJ(1));
t94 = t88 * pkin(1) + t72 * pkin(2) + t71 * pkin(7);
t93 = -t72 * pkin(7) + t95;
t92 = pkin(3) * t87 + pkin(8) * t84;
t91 = g(1) * t72 + g(2) * t71;
t90 = -g(1) * t88 - g(2) * t85;
t66 = pkin(5) * t73 + t104;
t65 = pkin(5) * t74 + t70;
t64 = g(1) * t71 - g(2) * t72;
t63 = -g(3) * t87 + t91 * t84;
t62 = -g(1) * (t72 * t100 + t71 * t73) - g(2) * (t71 * t100 - t72 * t73) - t74 * t105;
t61 = -g(1) * (-t72 * t101 + t71 * t74) - g(2) * (-t71 * t101 - t72 * t74) + t73 * t105;
t1 = [0, 0, 0, 0, 0, 0, t90, g(1) * t85 - g(2) * t88, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t91, t64, -g(3), t90 * pkin(1) - t106, 0, 0, 0, 0, 0, 0, -t91 * t87 - t105, t63, -t64, -g(1) * t94 - g(2) * t93 - t106, 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t96 + t103) - g(2) * (t71 * t96 - t72 * t83) - t86 * t105, -g(1) * (t71 * t86 - t72 * t98) - g(2) * (-t71 * t98 - t72 * t86) + t83 * t105, -t63, -g(1) * (t92 * t72 + t94) - g(2) * (t92 * t71 + t93) - g(3) * (pkin(3) * t84 - pkin(8) * t87 + t82) 0, 0, 0, 0, 0, 0, t62, t61, -t63, -g(1) * (pkin(4) * t103 + t94) - g(2) * (t70 * t102 - t71 * t97 + t95) - g(3) * (t84 * t70 + t87 * t89 + t82) + (-g(1) * (t70 * t87 - t97) - g(2) * (-pkin(7) - t104)) * t72, 0, 0, 0, 0, 0, 0, t62, t61, -t63, -g(1) * (t71 * t66 + t94) - g(2) * (t65 * t102 - t71 * t99 + t95) - g(3) * (t65 * t84 + t79 * t87 + t82) + (-g(1) * (t65 * t87 - t99) - g(2) * (-pkin(7) - t66)) * t72;];
U_reg  = t1;
