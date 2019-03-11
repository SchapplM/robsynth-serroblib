% Calculate inertial parameters regressor of potential energy for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:15
% EndTime: 2019-03-09 06:55:15
% DurationCPUTime: 0.11s
% Computational Cost: add. (200->62), mult. (136->72), div. (0->0), fcn. (124->12), ass. (0->38)
t91 = -pkin(8) - pkin(7);
t83 = qJ(3) + qJ(4);
t76 = qJ(5) + t83;
t69 = sin(t76);
t104 = g(3) * t69;
t84 = qJ(2) + pkin(6);
t103 = g(3) * t84;
t89 = cos(qJ(3));
t71 = t89 * pkin(3) + pkin(2);
t81 = qJ(1) + pkin(11);
t72 = sin(t81);
t85 = sin(qJ(6));
t102 = t72 * t85;
t88 = cos(qJ(6));
t101 = t72 * t88;
t73 = cos(t81);
t100 = t73 * t85;
t99 = t73 * t88;
t75 = cos(t83);
t65 = pkin(4) * t75 + t71;
t87 = sin(qJ(1));
t78 = t87 * pkin(1);
t82 = -pkin(9) + t91;
t98 = t72 * t65 + t73 * t82 + t78;
t86 = sin(qJ(3));
t97 = t86 * pkin(3) + t84;
t74 = sin(t83);
t96 = pkin(4) * t74 + t97;
t90 = cos(qJ(1));
t80 = t90 * pkin(1);
t95 = t73 * t65 - t72 * t82 + t80;
t70 = cos(t76);
t94 = pkin(5) * t70 + pkin(10) * t69;
t93 = g(1) * t73 + g(2) * t72;
t92 = -g(1) * t90 - g(2) * t87;
t64 = g(1) * t72 - g(2) * t73;
t61 = -g(3) * t70 + t93 * t69;
t1 = [0, 0, 0, 0, 0, 0, t92, g(1) * t87 - g(2) * t90, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t93, t64, -g(3), t92 * pkin(1) - t103, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t93 * t89, -g(3) * t89 + t93 * t86, -t64, -g(1) * (pkin(2) * t73 + pkin(7) * t72 + t80) - g(2) * (pkin(2) * t72 - pkin(7) * t73 + t78) - t103, 0, 0, 0, 0, 0, 0, -g(3) * t74 - t93 * t75, -g(3) * t75 + t93 * t74, -t64, -g(1) * (t73 * t71 - t72 * t91 + t80) - g(2) * (t72 * t71 + t73 * t91 + t78) - g(3) * t97, 0, 0, 0, 0, 0, 0, -t93 * t70 - t104, t61, -t64, -g(1) * t95 - g(2) * t98 - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t70 * t99 + t102) - g(2) * (t70 * t101 - t100) - t88 * t104, -g(1) * (-t70 * t100 + t101) - g(2) * (-t70 * t102 - t99) + t85 * t104, -t61, -g(1) * (t94 * t73 + t95) - g(2) * (t94 * t72 + t98) - g(3) * (pkin(5) * t69 - t70 * pkin(10) + t96);];
U_reg  = t1;
