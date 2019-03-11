% Calculate inertial parameters regressor of potential energy for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:33
% EndTime: 2019-03-09 02:39:33
% DurationCPUTime: 0.14s
% Computational Cost: add. (211->72), mult. (160->90), div. (0->0), fcn. (152->12), ass. (0->41)
t82 = qJ(3) + pkin(10);
t73 = sin(t82);
t108 = g(3) * t73;
t87 = qJ(2) + pkin(6);
t107 = g(3) * t87;
t88 = -pkin(8) - qJ(5);
t106 = t73 * t88;
t83 = qJ(1) + pkin(9);
t74 = sin(t83);
t76 = cos(t82);
t105 = t74 * t76;
t84 = sin(pkin(11));
t104 = t74 * t84;
t85 = cos(pkin(11));
t103 = t74 * t85;
t77 = cos(t83);
t102 = t76 * t77;
t101 = t77 * t84;
t100 = t77 * t85;
t91 = cos(qJ(3));
t71 = pkin(3) * t91 + pkin(2);
t92 = cos(qJ(1));
t80 = t92 * pkin(1);
t99 = t77 * t71 + t80;
t90 = sin(qJ(1));
t79 = t90 * pkin(1);
t86 = -qJ(4) - pkin(7);
t98 = t74 * t71 + t77 * t86 + t79;
t89 = sin(qJ(3));
t97 = t89 * pkin(3) + t87;
t96 = -t74 * t86 + t99;
t95 = g(1) * t77 + g(2) * t74;
t94 = -g(1) * t92 - g(2) * t90;
t93 = pkin(4) * t76 + qJ(5) * t73;
t81 = pkin(11) + qJ(6);
t75 = cos(t81);
t72 = sin(t81);
t70 = pkin(5) * t85 + pkin(4);
t64 = g(1) * t74 - g(2) * t77;
t63 = -g(3) * t76 + t95 * t73;
t1 = [0, 0, 0, 0, 0, 0, t94, g(1) * t90 - g(2) * t92, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t95, t64, -g(3), t94 * pkin(1) - t107, 0, 0, 0, 0, 0, 0, -g(3) * t89 - t95 * t91, -g(3) * t91 + t95 * t89, -t64, -g(1) * (pkin(2) * t77 + pkin(7) * t74 + t80) - g(2) * (pkin(2) * t74 - pkin(7) * t77 + t79) - t107, 0, 0, 0, 0, 0, 0, -t95 * t76 - t108, t63, -t64, -g(1) * t96 - g(2) * t98 - g(3) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t76 * t100 + t104) - g(2) * (t76 * t103 - t101) - t85 * t108, -g(1) * (-t76 * t101 + t103) - g(2) * (-t76 * t104 - t100) + t84 * t108, -t63, -g(1) * (t93 * t77 + t96) - g(2) * (t93 * t74 + t98) - g(3) * (pkin(4) * t73 - qJ(5) * t76 + t97) 0, 0, 0, 0, 0, 0, -g(1) * (t75 * t102 + t72 * t74) - g(2) * (t75 * t105 - t72 * t77) - t75 * t108, -g(1) * (-t72 * t102 + t74 * t75) - g(2) * (-t72 * t105 - t75 * t77) + t72 * t108, -t63, -g(1) * (t70 * t102 - t77 * t106 + t99) - g(2) * (-pkin(5) * t101 + t98) - g(3) * (t70 * t73 + t76 * t88 + t97) + (-g(1) * (pkin(5) * t84 - t86) - g(2) * (t70 * t76 - t106)) * t74;];
U_reg  = t1;
