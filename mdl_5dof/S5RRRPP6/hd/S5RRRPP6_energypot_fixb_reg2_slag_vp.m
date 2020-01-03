% Calculate inertial parameters regressor of potential energy for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:01:58
% EndTime: 2019-12-31 21:01:59
% DurationCPUTime: 0.12s
% Computational Cost: add. (129->61), mult. (184->77), div. (0->0), fcn. (187->8), ass. (0->37)
t97 = g(3) * pkin(5);
t76 = sin(qJ(2));
t96 = g(3) * t76;
t74 = -qJ(4) - pkin(7);
t95 = t74 * t76;
t75 = sin(qJ(3));
t77 = sin(qJ(1));
t94 = t77 * t75;
t79 = cos(qJ(2));
t93 = t77 * t79;
t73 = qJ(3) + pkin(8);
t67 = sin(t73);
t80 = cos(qJ(1));
t92 = t80 * t67;
t68 = cos(t73);
t91 = t80 * t68;
t90 = t80 * t75;
t78 = cos(qJ(3));
t89 = t80 * t78;
t88 = t80 * pkin(1) + t77 * pkin(6);
t66 = t78 * pkin(3) + pkin(2);
t87 = t76 * t66 + t79 * t74 + pkin(5);
t70 = t77 * pkin(1);
t86 = -t80 * pkin(6) + t70;
t85 = pkin(2) * t79 + pkin(7) * t76;
t84 = g(1) * t80 + g(2) * t77;
t83 = pkin(3) * t94 + t88 + (t66 * t79 - t95) * t80;
t54 = t67 * t93 + t91;
t56 = -t77 * t68 + t79 * t92;
t82 = g(1) * t56 + g(2) * t54 + t67 * t96;
t81 = -t77 * t95 + t66 * t93 + t70 + (-pkin(3) * t75 - pkin(6)) * t80;
t61 = g(1) * t77 - g(2) * t80;
t58 = -g(3) * t79 + t84 * t76;
t57 = t77 * t67 + t79 * t91;
t55 = t68 * t93 - t92;
t53 = -g(1) * t57 - g(2) * t55 - t68 * t96;
t1 = [0, 0, 0, 0, 0, 0, -t84, t61, -g(3), -t97, 0, 0, 0, 0, 0, 0, -t84 * t79 - t96, t58, -t61, -g(1) * t88 - g(2) * t86 - t97, 0, 0, 0, 0, 0, 0, -g(1) * (t79 * t89 + t94) - g(2) * (t78 * t93 - t90) - t78 * t96, -g(1) * (t77 * t78 - t79 * t90) - g(2) * (-t75 * t93 - t89) + t75 * t96, -t58, -g(1) * (t85 * t80 + t88) - g(2) * (t85 * t77 + t86) - g(3) * (t76 * pkin(2) - t79 * pkin(7) + pkin(5)), 0, 0, 0, 0, 0, 0, t53, t82, -t58, -g(1) * t83 - g(2) * t81 - g(3) * t87, 0, 0, 0, 0, 0, 0, t53, -t58, -t82, -g(1) * (t57 * pkin(4) + t56 * qJ(5) + t83) - g(2) * (t55 * pkin(4) + t54 * qJ(5) + t81) - g(3) * ((pkin(4) * t68 + qJ(5) * t67) * t76 + t87);];
U_reg = t1;
