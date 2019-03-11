% Calculate inertial parameters regressor of potential energy for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% U_reg [1x(6*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPRRPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:44
% EndTime: 2019-03-09 04:58:44
% DurationCPUTime: 0.11s
% Computational Cost: add. (200->62), mult. (136->72), div. (0->0), fcn. (124->12), ass. (0->38)
t91 = -pkin(8) - pkin(7);
t83 = qJ(3) + qJ(4);
t74 = pkin(11) + t83;
t68 = sin(t74);
t104 = g(3) * t68;
t84 = qJ(2) + pkin(6);
t103 = g(3) * t84;
t89 = cos(qJ(3));
t71 = t89 * pkin(3) + pkin(2);
t82 = qJ(1) + pkin(10);
t72 = sin(t82);
t85 = sin(qJ(6));
t102 = t72 * t85;
t88 = cos(qJ(6));
t101 = t72 * t88;
t73 = cos(t82);
t100 = t73 * t85;
t99 = t73 * t88;
t76 = cos(t83);
t65 = pkin(4) * t76 + t71;
t87 = sin(qJ(1));
t78 = t87 * pkin(1);
t81 = -qJ(5) + t91;
t98 = t72 * t65 + t73 * t81 + t78;
t86 = sin(qJ(3));
t97 = t86 * pkin(3) + t84;
t75 = sin(t83);
t96 = pkin(4) * t75 + t97;
t90 = cos(qJ(1));
t80 = t90 * pkin(1);
t95 = t73 * t65 - t72 * t81 + t80;
t69 = cos(t74);
t94 = pkin(5) * t69 + pkin(9) * t68;
t93 = g(1) * t73 + g(2) * t72;
t92 = -g(1) * t90 - g(2) * t87;
t64 = g(1) * t72 - g(2) * t73;
t61 = -g(3) * t69 + t93 * t68;
t1 = [0, 0, 0, 0, 0, 0, t92, g(1) * t87 - g(2) * t90, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t93, t64, -g(3), t92 * pkin(1) - t103, 0, 0, 0, 0, 0, 0, -g(3) * t86 - t93 * t89, -g(3) * t89 + t93 * t86, -t64, -g(1) * (pkin(2) * t73 + pkin(7) * t72 + t80) - g(2) * (pkin(2) * t72 - pkin(7) * t73 + t78) - t103, 0, 0, 0, 0, 0, 0, -g(3) * t75 - t93 * t76, -g(3) * t76 + t93 * t75, -t64, -g(1) * (t73 * t71 - t72 * t91 + t80) - g(2) * (t72 * t71 + t73 * t91 + t78) - g(3) * t97, 0, 0, 0, 0, 0, 0, -t93 * t69 - t104, t61, -t64, -g(1) * t95 - g(2) * t98 - g(3) * t96, 0, 0, 0, 0, 0, 0, -g(1) * (t69 * t99 + t102) - g(2) * (t69 * t101 - t100) - t88 * t104, -g(1) * (-t69 * t100 + t101) - g(2) * (-t69 * t102 - t99) + t85 * t104, -t61, -g(1) * (t94 * t73 + t95) - g(2) * (t94 * t72 + t98) - g(3) * (pkin(5) * t68 - t69 * pkin(9) + t96);];
U_reg  = t1;
