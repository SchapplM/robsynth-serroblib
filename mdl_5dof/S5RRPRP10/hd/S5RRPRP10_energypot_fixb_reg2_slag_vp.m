% Calculate inertial parameters regressor of potential energy for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP10_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:04
% EndTime: 2019-12-31 20:11:05
% DurationCPUTime: 0.13s
% Computational Cost: add. (86->56), mult. (161->64), div. (0->0), fcn. (156->6), ass. (0->34)
t70 = sin(qJ(1));
t69 = sin(qJ(2));
t80 = qJ(3) * t69;
t72 = cos(qJ(2));
t85 = t70 * t72;
t92 = pkin(2) * t85 + t70 * t80;
t91 = g(3) * pkin(5);
t68 = sin(qJ(4));
t90 = pkin(4) * t68;
t89 = g(3) * t72;
t88 = t69 * pkin(2) + pkin(5);
t87 = t70 * t68;
t71 = cos(qJ(4));
t86 = t70 * t71;
t73 = cos(qJ(1));
t84 = t72 * t73;
t83 = t73 * t68;
t82 = t73 * t71;
t81 = t73 * pkin(1) + t70 * pkin(6);
t79 = t69 * t87;
t64 = t70 * pkin(1);
t78 = t64 + t92;
t77 = -t73 * pkin(6) + t64;
t76 = pkin(2) * t84 + t73 * t80 + t81;
t75 = -t72 * qJ(3) + t88;
t74 = g(1) * t73 + g(2) * t70;
t67 = -qJ(5) - pkin(7);
t61 = t71 * pkin(4) + pkin(3);
t56 = g(1) * t70 - g(2) * t73;
t55 = g(3) * t69 + t74 * t72;
t54 = t74 * t69 - t89;
t53 = -g(1) * (t69 * t83 + t86) - g(2) * (t79 - t82) + t68 * t89;
t52 = -g(1) * (t69 * t82 - t87) - g(2) * (t69 * t86 + t83) + t71 * t89;
t1 = [0, 0, 0, 0, 0, 0, -t74, t56, -g(3), -t91, 0, 0, 0, 0, 0, 0, -t55, t54, -t56, -g(1) * t81 - g(2) * t77 - t91, 0, 0, 0, 0, 0, 0, -t56, t55, -t54, -g(1) * t76 - g(2) * (t77 + t92) - g(3) * t75, 0, 0, 0, 0, 0, 0, t53, t52, -t55, -g(1) * (t70 * pkin(3) + pkin(7) * t84 + t76) - g(2) * (pkin(7) * t85 + (-pkin(3) - pkin(6)) * t73 + t78) - g(3) * (t69 * pkin(7) + t75), 0, 0, 0, 0, 0, 0, t53, t52, -t55, -g(1) * (t70 * t61 + t76) - g(2) * (pkin(4) * t79 - t67 * t85 + t78) - g(3) * (-t69 * t67 + (-qJ(3) - t90) * t72 + t88) + (-g(1) * (-t67 * t72 + t69 * t90) - g(2) * (-pkin(6) - t61)) * t73;];
U_reg = t1;
