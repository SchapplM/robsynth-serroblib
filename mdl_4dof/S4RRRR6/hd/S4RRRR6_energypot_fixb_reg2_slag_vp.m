% Calculate inertial parameters regressor of potential energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:30:50
% EndTime: 2019-12-31 17:30:50
% DurationCPUTime: 0.11s
% Computational Cost: add. (114->56), mult. (262->88), div. (0->0), fcn. (311->10), ass. (0->39)
t69 = cos(pkin(4));
t95 = t69 * pkin(6) + pkin(5);
t68 = sin(pkin(4));
t72 = sin(qJ(2));
t94 = t68 * t72;
t73 = sin(qJ(1));
t93 = t68 * t73;
t75 = cos(qJ(3));
t92 = t68 * t75;
t76 = cos(qJ(2));
t91 = t68 * t76;
t77 = cos(qJ(1));
t90 = t68 * t77;
t89 = t73 * t72;
t88 = t73 * t76;
t87 = t77 * t72;
t86 = t77 * t76;
t85 = t77 * pkin(1) + pkin(6) * t93;
t84 = t73 * pkin(1) - pkin(6) * t90;
t83 = g(1) * t73 - g(2) * t77;
t58 = t69 * t88 + t87;
t59 = -t69 * t89 + t86;
t82 = t59 * pkin(2) + t58 * pkin(7) + t85;
t81 = pkin(2) * t94 - pkin(7) * t91 + t95;
t57 = t69 * t87 + t88;
t71 = sin(qJ(3));
t48 = t57 * t71 + t75 * t90;
t50 = t59 * t71 - t73 * t92;
t54 = -t69 * t75 + t71 * t94;
t80 = g(1) * t50 + g(2) * t48 + g(3) * t54;
t56 = -t69 * t86 + t89;
t79 = t57 * pkin(2) + t56 * pkin(7) + t84;
t78 = -g(1) * t58 - g(2) * t56 + g(3) * t91;
t74 = cos(qJ(4));
t70 = sin(qJ(4));
t55 = t69 * t71 + t72 * t92;
t51 = t59 * t75 + t71 * t93;
t49 = t57 * t75 - t71 * t90;
t1 = [0, 0, 0, 0, 0, 0, -g(1) * t77 - g(2) * t73, t83, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t59 - g(2) * t57 - g(3) * t94, -t78, -g(3) * t69 - t83 * t68, -g(1) * t85 - g(2) * t84 - g(3) * t95, 0, 0, 0, 0, 0, 0, -g(1) * t51 - g(2) * t49 - g(3) * t55, t80, t78, -g(1) * t82 - g(2) * t79 - g(3) * t81, 0, 0, 0, 0, 0, 0, -g(1) * (t51 * t74 + t58 * t70) - g(2) * (t49 * t74 + t56 * t70) - g(3) * (t55 * t74 - t70 * t91), -g(1) * (-t51 * t70 + t58 * t74) - g(2) * (-t49 * t70 + t56 * t74) - g(3) * (-t55 * t70 - t74 * t91), -t80, -g(1) * (t51 * pkin(3) + t50 * pkin(8) + t82) - g(2) * (t49 * pkin(3) + t48 * pkin(8) + t79) - g(3) * (t55 * pkin(3) + t54 * pkin(8) + t81);];
U_reg = t1;
