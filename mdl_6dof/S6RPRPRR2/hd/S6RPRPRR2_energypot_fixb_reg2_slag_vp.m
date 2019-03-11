% Calculate inertial parameters regressor of potential energy for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:38:56
% EndTime: 2019-03-09 03:38:57
% DurationCPUTime: 0.16s
% Computational Cost: add. (211->71), mult. (160->87), div. (0->0), fcn. (152->12), ass. (0->43)
t89 = cos(qJ(5));
t70 = pkin(5) * t89 + pkin(4);
t81 = qJ(3) + pkin(11);
t72 = sin(t81);
t74 = cos(t81);
t92 = -pkin(9) - pkin(8);
t112 = t70 * t74 - t72 * t92;
t111 = g(3) * t72;
t85 = qJ(2) + pkin(6);
t110 = g(3) * t85;
t82 = qJ(1) + pkin(10);
t73 = sin(t82);
t83 = qJ(5) + qJ(6);
t76 = sin(t83);
t107 = t73 * t76;
t77 = cos(t83);
t106 = t73 * t77;
t86 = sin(qJ(5));
t105 = t73 * t86;
t104 = t73 * t89;
t75 = cos(t82);
t103 = t75 * t76;
t102 = t75 * t77;
t101 = t75 * t86;
t100 = t75 * t89;
t90 = cos(qJ(3));
t71 = pkin(3) * t90 + pkin(2);
t91 = cos(qJ(1));
t80 = t91 * pkin(1);
t99 = t75 * t71 + t80;
t88 = sin(qJ(1));
t79 = t88 * pkin(1);
t84 = -qJ(4) - pkin(7);
t98 = t73 * t71 + t75 * t84 + t79;
t87 = sin(qJ(3));
t97 = t87 * pkin(3) + t85;
t96 = -t73 * t84 + t99;
t95 = pkin(4) * t74 + pkin(8) * t72;
t94 = g(1) * t75 + g(2) * t73;
t93 = -g(1) * t91 - g(2) * t88;
t64 = g(1) * t73 - g(2) * t75;
t63 = -g(3) * t74 + t94 * t72;
t1 = [0, 0, 0, 0, 0, 0, t93, g(1) * t88 - g(2) * t91, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t94, t64, -g(3), t93 * pkin(1) - t110, 0, 0, 0, 0, 0, 0, -g(3) * t87 - t94 * t90, -g(3) * t90 + t94 * t87, -t64, -g(1) * (pkin(2) * t75 + pkin(7) * t73 + t80) - g(2) * (pkin(2) * t73 - pkin(7) * t75 + t79) - t110, 0, 0, 0, 0, 0, 0, -t94 * t74 - t111, t63, -t64, -g(1) * t96 - g(2) * t98 - g(3) * t97, 0, 0, 0, 0, 0, 0, -g(1) * (t74 * t100 + t105) - g(2) * (t74 * t104 - t101) - t89 * t111, -g(1) * (-t74 * t101 + t104) - g(2) * (-t74 * t105 - t100) + t86 * t111, -t63, -g(1) * (t95 * t75 + t96) - g(2) * (t95 * t73 + t98) - g(3) * (pkin(4) * t72 - pkin(8) * t74 + t97) 0, 0, 0, 0, 0, 0, -g(1) * (t74 * t102 + t107) - g(2) * (t74 * t106 - t103) - t77 * t111, -g(1) * (-t74 * t103 + t106) - g(2) * (-t74 * t107 - t102) + t76 * t111, -t63, -g(1) * (t112 * t75 + t99) - g(2) * (-pkin(5) * t101 + t98) - g(3) * (t72 * t70 + t74 * t92 + t97) + (-g(1) * (pkin(5) * t86 - t84) - g(2) * t112) * t73;];
U_reg  = t1;
