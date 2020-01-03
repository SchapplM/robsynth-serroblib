% Calculate inertial parameters regressor of potential energy for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR15_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:07
% EndTime: 2019-12-31 20:43:07
% DurationCPUTime: 0.13s
% Computational Cost: add. (96->65), mult. (161->78), div. (0->0), fcn. (156->8), ass. (0->40)
t72 = sin(qJ(1));
t71 = sin(qJ(2));
t83 = qJ(3) * t71;
t74 = cos(qJ(2));
t91 = t72 * t74;
t100 = pkin(2) * t91 + t72 * t83;
t99 = g(3) * pkin(5);
t70 = sin(qJ(4));
t98 = pkin(4) * t70;
t97 = g(3) * t74;
t96 = t71 * pkin(2) + pkin(5);
t69 = qJ(4) + qJ(5);
t62 = sin(t69);
t95 = t72 * t62;
t63 = cos(t69);
t94 = t72 * t63;
t93 = t72 * t70;
t73 = cos(qJ(4));
t92 = t72 * t73;
t75 = cos(qJ(1));
t90 = t74 * t75;
t76 = -pkin(8) - pkin(7);
t89 = t74 * t76;
t88 = t75 * t62;
t87 = t75 * t63;
t86 = t75 * t70;
t85 = t75 * t73;
t84 = t75 * pkin(1) + t72 * pkin(6);
t82 = t71 * t93;
t66 = t72 * pkin(1);
t81 = t66 + t100;
t80 = -t75 * pkin(6) + t66;
t79 = pkin(2) * t90 + t75 * t83 + t84;
t78 = -t74 * qJ(3) + t96;
t77 = g(1) * t75 + g(2) * t72;
t61 = t73 * pkin(4) + pkin(3);
t56 = g(1) * t72 - g(2) * t75;
t55 = g(3) * t71 + t77 * t74;
t54 = t77 * t71 - t97;
t1 = [0, 0, 0, 0, 0, 0, -t77, t56, -g(3), -t99, 0, 0, 0, 0, 0, 0, -t55, t54, -t56, -g(1) * t84 - g(2) * t80 - t99, 0, 0, 0, 0, 0, 0, -t56, t55, -t54, -g(1) * t79 - g(2) * (t80 + t100) - g(3) * t78, 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t86 + t92) - g(2) * (t82 - t85) + t70 * t97, -g(1) * (t71 * t85 - t93) - g(2) * (t71 * t92 + t86) + t73 * t97, -t55, -g(1) * (t72 * pkin(3) + pkin(7) * t90 + t79) - g(2) * (pkin(7) * t91 + (-pkin(3) - pkin(6)) * t75 + t81) - g(3) * (t71 * pkin(7) + t78), 0, 0, 0, 0, 0, 0, -g(1) * (t71 * t88 + t94) - g(2) * (t71 * t95 - t87) + t62 * t97, -g(1) * (t71 * t87 - t95) - g(2) * (t71 * t94 + t88) + t63 * t97, -t55, -g(1) * (t72 * t61 + t79) - g(2) * (pkin(4) * t82 - t72 * t89 + t81) - g(3) * (-t71 * t76 + (-qJ(3) - t98) * t74 + t96) + (-g(1) * (t71 * t98 - t89) - g(2) * (-pkin(6) - t61)) * t75;];
U_reg = t1;
