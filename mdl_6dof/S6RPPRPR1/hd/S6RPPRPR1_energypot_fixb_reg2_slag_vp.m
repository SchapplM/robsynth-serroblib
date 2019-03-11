% Calculate inertial parameters regressor of potential energy for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:40:01
% EndTime: 2019-03-09 01:40:01
% DurationCPUTime: 0.14s
% Computational Cost: add. (211->72), mult. (160->90), div. (0->0), fcn. (152->12), ass. (0->41)
t79 = pkin(10) + qJ(4);
t70 = sin(t79);
t105 = g(3) * t70;
t85 = qJ(2) + pkin(6);
t104 = g(3) * t85;
t86 = -pkin(8) - qJ(5);
t103 = t70 * t86;
t80 = qJ(1) + pkin(9);
t71 = sin(t80);
t73 = cos(t79);
t102 = t71 * t73;
t81 = sin(pkin(11));
t101 = t71 * t81;
t83 = cos(pkin(11));
t100 = t71 * t83;
t74 = cos(t80);
t99 = t73 * t74;
t98 = t74 * t81;
t97 = t74 * t83;
t84 = cos(pkin(10));
t68 = pkin(3) * t84 + pkin(2);
t89 = cos(qJ(1));
t77 = t89 * pkin(1);
t96 = t74 * t68 + t77;
t88 = sin(qJ(1));
t76 = t88 * pkin(1);
t87 = -pkin(7) - qJ(3);
t95 = t71 * t68 + t74 * t87 + t76;
t82 = sin(pkin(10));
t94 = t82 * pkin(3) + t85;
t93 = -t71 * t87 + t96;
t92 = g(1) * t74 + g(2) * t71;
t91 = -g(1) * t89 - g(2) * t88;
t90 = pkin(4) * t73 + qJ(5) * t70;
t78 = pkin(11) + qJ(6);
t72 = cos(t78);
t69 = sin(t78);
t67 = pkin(5) * t83 + pkin(4);
t61 = g(1) * t71 - g(2) * t74;
t60 = -g(3) * t73 + t92 * t70;
t1 = [0, 0, 0, 0, 0, 0, t91, g(1) * t88 - g(2) * t89, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t92, t61, -g(3), t91 * pkin(1) - t104, 0, 0, 0, 0, 0, 0, -g(3) * t82 - t92 * t84, -g(3) * t84 + t92 * t82, -t61, -g(1) * (pkin(2) * t74 + qJ(3) * t71 + t77) - g(2) * (pkin(2) * t71 - qJ(3) * t74 + t76) - t104, 0, 0, 0, 0, 0, 0, -t92 * t73 - t105, t60, -t61, -g(1) * t93 - g(2) * t95 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (t73 * t97 + t101) - g(2) * (t73 * t100 - t98) - t83 * t105, -g(1) * (-t73 * t98 + t100) - g(2) * (-t73 * t101 - t97) + t81 * t105, -t60, -g(1) * (t90 * t74 + t93) - g(2) * (t90 * t71 + t95) - g(3) * (pkin(4) * t70 - qJ(5) * t73 + t94) 0, 0, 0, 0, 0, 0, -g(1) * (t69 * t71 + t72 * t99) - g(2) * (t72 * t102 - t69 * t74) - t72 * t105, -g(1) * (-t69 * t99 + t71 * t72) - g(2) * (-t69 * t102 - t72 * t74) + t69 * t105, -t60, -g(1) * (-t74 * t103 + t67 * t99 + t96) - g(2) * (-pkin(5) * t98 + t95) - g(3) * (t67 * t70 + t73 * t86 + t94) + (-g(1) * (pkin(5) * t81 - t87) - g(2) * (t67 * t73 - t103)) * t71;];
U_reg  = t1;
