% Calculate inertial parameters regressor of potential energy for
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RRPPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:15:37
% EndTime: 2019-03-09 08:15:38
% DurationCPUTime: 0.18s
% Computational Cost: add. (128->76), mult. (214->88), div. (0->0), fcn. (206->8), ass. (0->43)
t85 = sin(qJ(2));
t100 = qJ(3) * t85;
t86 = sin(qJ(1));
t87 = cos(qJ(2));
t103 = t86 * t87;
t114 = pkin(2) * t103 + t86 * t100;
t88 = cos(qJ(1));
t113 = pkin(3) * t103 + t88 * qJ(4);
t112 = g(3) * pkin(6);
t82 = sin(pkin(9));
t111 = pkin(5) * t82;
t110 = g(3) * t87;
t109 = t85 * pkin(2) + pkin(6);
t108 = t85 * t88;
t81 = pkin(9) + qJ(6);
t71 = sin(t81);
t107 = t86 * t71;
t72 = cos(t81);
t106 = t86 * t72;
t105 = t86 * t82;
t83 = cos(pkin(9));
t104 = t86 * t83;
t102 = t87 * t88;
t101 = t88 * pkin(1) + t86 * pkin(7);
t74 = t85 * pkin(3);
t99 = t74 + t109;
t77 = t86 * pkin(1);
t98 = -pkin(7) * t88 + t77;
t97 = pkin(2) * t102 + t88 * t100 + t101;
t96 = -qJ(3) * t87 + t109;
t95 = pkin(3) * t102 + t97;
t94 = g(1) * t88 + g(2) * t86;
t93 = pkin(4) * t85 + qJ(5) * t87;
t70 = pkin(5) * t83 + pkin(4);
t84 = -pkin(8) - qJ(5);
t92 = t70 * t85 - t84 * t87;
t91 = t98 + t114;
t90 = t91 + t113;
t89 = -t86 * qJ(4) + t95;
t63 = g(1) * t86 - g(2) * t88;
t62 = g(3) * t85 + t87 * t94;
t61 = t85 * t94 - t110;
t1 = [0, 0, 0, 0, 0, 0, -t94, t63, -g(3), -t112, 0, 0, 0, 0, 0, 0, -t62, t61, -t63, -g(1) * t101 - g(2) * t98 - t112, 0, 0, 0, 0, 0, 0, -t62, -t63, -t61, -g(1) * t97 - g(2) * t91 - g(3) * t96, 0, 0, 0, 0, 0, 0, -t61, t62, t63, -g(1) * t89 - g(2) * t90 - g(3) * (t74 + t96) 0, 0, 0, 0, 0, 0, -g(1) * (t83 * t108 - t105) - g(2) * (t85 * t104 + t82 * t88) + t83 * t110, -g(1) * (-t82 * t108 - t104) - g(2) * (-t85 * t105 + t83 * t88) - t82 * t110, -t62, -g(1) * (t88 * t93 + t89) - g(2) * (t86 * t93 + t90) - g(3) * (qJ(5) * t85 + (-pkin(4) - qJ(3)) * t87 + t99) 0, 0, 0, 0, 0, 0, -g(1) * (t72 * t108 - t107) - g(2) * (t85 * t106 + t71 * t88) + t72 * t110, -g(1) * (-t71 * t108 - t106) - g(2) * (-t85 * t107 + t72 * t88) - t71 * t110, -t62, -g(1) * t95 - g(2) * (t77 + t113 + t114) - g(3) * (-t85 * t84 + (-qJ(3) - t70) * t87 + t99) + (-g(1) * t92 - g(2) * (-pkin(7) + t111)) * t88 + (-g(1) * (-qJ(4) - t111) - g(2) * t92) * t86;];
U_reg  = t1;
