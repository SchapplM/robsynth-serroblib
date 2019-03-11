% Calculate inertial parameters regressor of potential energy for
% S6RPPRRR1
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S6RPPRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:18:38
% EndTime: 2019-03-09 02:18:38
% DurationCPUTime: 0.11s
% Computational Cost: add. (200->62), mult. (136->72), div. (0->0), fcn. (124->12), ass. (0->38)
t80 = pkin(11) + qJ(4);
t74 = qJ(5) + t80;
t67 = sin(t74);
t102 = g(3) * t67;
t84 = qJ(2) + pkin(6);
t101 = g(3) * t84;
t83 = cos(pkin(11));
t69 = t83 * pkin(3) + pkin(2);
t81 = qJ(1) + pkin(10);
t71 = sin(t81);
t86 = sin(qJ(6));
t100 = t71 * t86;
t88 = cos(qJ(6));
t99 = t71 * t88;
t73 = cos(t81);
t98 = t73 * t86;
t97 = t73 * t88;
t85 = -pkin(7) - qJ(3);
t72 = cos(t80);
t63 = pkin(4) * t72 + t69;
t87 = sin(qJ(1));
t77 = t87 * pkin(1);
t79 = -pkin(8) + t85;
t96 = t71 * t63 + t73 * t79 + t77;
t82 = sin(pkin(11));
t95 = t82 * pkin(3) + t84;
t70 = sin(t80);
t94 = pkin(4) * t70 + t95;
t89 = cos(qJ(1));
t78 = t89 * pkin(1);
t93 = t73 * t63 - t71 * t79 + t78;
t68 = cos(t74);
t92 = pkin(5) * t68 + pkin(9) * t67;
t91 = g(1) * t73 + g(2) * t71;
t90 = -g(1) * t89 - g(2) * t87;
t62 = g(1) * t71 - g(2) * t73;
t59 = -g(3) * t68 + t91 * t67;
t1 = [0, 0, 0, 0, 0, 0, t90, g(1) * t87 - g(2) * t89, -g(3), -g(3) * pkin(6), 0, 0, 0, 0, 0, 0, -t91, t62, -g(3), t90 * pkin(1) - t101, 0, 0, 0, 0, 0, 0, -g(3) * t82 - t91 * t83, -g(3) * t83 + t91 * t82, -t62, -g(1) * (t73 * pkin(2) + t71 * qJ(3) + t78) - g(2) * (t71 * pkin(2) - t73 * qJ(3) + t77) - t101, 0, 0, 0, 0, 0, 0, -g(3) * t70 - t91 * t72, -g(3) * t72 + t91 * t70, -t62, -g(1) * (t73 * t69 - t71 * t85 + t78) - g(2) * (t71 * t69 + t73 * t85 + t77) - g(3) * t95, 0, 0, 0, 0, 0, 0, -t91 * t68 - t102, t59, -t62, -g(1) * t93 - g(2) * t96 - g(3) * t94, 0, 0, 0, 0, 0, 0, -g(1) * (t68 * t97 + t100) - g(2) * (t68 * t99 - t98) - t88 * t102, -g(1) * (-t68 * t98 + t99) - g(2) * (-t68 * t100 - t97) + t86 * t102, -t59, -g(1) * (t92 * t73 + t93) - g(2) * (t92 * t71 + t96) - g(3) * (t67 * pkin(5) - t68 * pkin(9) + t94);];
U_reg  = t1;
