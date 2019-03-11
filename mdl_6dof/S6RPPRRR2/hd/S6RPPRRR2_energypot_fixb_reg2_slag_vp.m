% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:21:21
% EndTime: 2019-03-09 02:21:21
% DurationCPUTime: 0.16s
% Computational Cost: add. (211->71), mult. (160->87), div. (0->0), fcn. (152->12), ass. (0->43)
t87 = cos(qJ(5));
t68 = pkin(5) * t87 + pkin(4);
t78 = pkin(11) + qJ(4);
t69 = sin(t78);
t71 = cos(t78);
t89 = -pkin(9) - pkin(8);
t109 = t68 * t71 - t69 * t89;
t108 = g(3) * t69;
t83 = qJ(2) + pkin(6);
t107 = g(3) * t83;
t79 = qJ(1) + pkin(10);
t70 = sin(t79);
t80 = qJ(5) + qJ(6);
t73 = sin(t80);
t104 = t70 * t73;
t74 = cos(t80);
t103 = t70 * t74;
t85 = sin(qJ(5));
t102 = t70 * t85;
t101 = t70 * t87;
t72 = cos(t79);
t100 = t72 * t73;
t99 = t72 * t74;
t98 = t72 * t85;
t97 = t72 * t87;
t82 = cos(pkin(11));
t67 = pkin(3) * t82 + pkin(2);
t88 = cos(qJ(1));
t77 = t88 * pkin(1);
t96 = t72 * t67 + t77;
t86 = sin(qJ(1));
t76 = t86 * pkin(1);
t84 = -pkin(7) - qJ(3);
t95 = t70 * t67 + t72 * t84 + t76;
t81 = sin(pkin(11));
t94 = t81 * pkin(3) + t83;
t93 = -t70 * t84 + t96;
t92 = pkin(4) * t71 + pkin(8) * t69;
t91 = g(1) * t72 + g(2) * t70;
t90 = -g(1) * t88 - g(2) * t86;
t61 = g(1) * t70 - g(2) * t72;
t60 = -g(3) * t71 + t91 * t69;
t1 = [0, 0, 0, 0, 0, 0, t90, g(1) * t86 - g(2) * t88, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t91, t61, -g(3), t90 * pkin(1) - t107, 0, 0, 0, 0, 0, 0, -g(3) * t81 - t91 * t82, -g(3) * t82 + t91 * t81, -t61, -g(1) * (pkin(2) * t72 + qJ(3) * t70 + t77) - g(2) * (pkin(2) * t70 - qJ(3) * t72 + t76) - t107, 0, 0, 0, 0, 0, 0, -t91 * t71 - t108, t60, -t61, -g(1) * t93 - g(2) * t95 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t97 + t102) - g(2) * (t71 * t101 - t98) - t87 * t108, -g(1) * (-t71 * t98 + t101) - g(2) * (-t71 * t102 - t97) + t85 * t108, -t60, -g(1) * (t92 * t72 + t93) - g(2) * (t92 * t70 + t95) - g(3) * (pkin(4) * t69 - pkin(8) * t71 + t94) 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t99 + t104) - g(2) * (t71 * t103 - t100) - t74 * t108, -g(1) * (-t71 * t100 + t103) - g(2) * (-t71 * t104 - t99) + t73 * t108, -t60, -g(1) * (t109 * t72 + t96) - g(2) * (-pkin(5) * t98 + t95) - g(3) * (t69 * t68 + t71 * t89 + t94) + (-g(1) * (pkin(5) * t85 - t84) - g(2) * t109) * t70;];
U_reg  = t1;
