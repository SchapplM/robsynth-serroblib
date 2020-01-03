% Calculate inertial parameters regressor of potential energy for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:33
% EndTime: 2019-12-31 21:54:33
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->52), mult. (143->65), div. (0->0), fcn. (138->8), ass. (0->31)
t73 = cos(qJ(4));
t61 = t73 * pkin(4) + pkin(3);
t68 = qJ(2) + qJ(3);
t64 = sin(t68);
t65 = cos(t68);
t69 = -qJ(5) - pkin(8);
t90 = t61 * t65 - t64 * t69;
t89 = g(3) * pkin(5);
t88 = g(3) * t64;
t71 = sin(qJ(2));
t87 = t71 * pkin(2) + pkin(5);
t70 = sin(qJ(4));
t72 = sin(qJ(1));
t84 = t72 * t70;
t83 = t72 * t73;
t75 = cos(qJ(1));
t82 = t75 * t70;
t81 = t75 * t73;
t74 = cos(qJ(2));
t62 = t74 * pkin(2) + pkin(1);
t76 = -pkin(7) - pkin(6);
t80 = t72 * t62 + t75 * t76;
t59 = t75 * t62;
t79 = -t72 * t76 + t59;
t78 = pkin(3) * t65 + pkin(8) * t64;
t77 = g(1) * t75 + g(2) * t72;
t57 = g(1) * t72 - g(2) * t75;
t56 = -g(3) * t65 + t77 * t64;
t55 = -g(1) * (t65 * t81 + t84) - g(2) * (t65 * t83 - t82) - t73 * t88;
t54 = -g(1) * (-t65 * t82 + t83) - g(2) * (-t65 * t84 - t81) + t70 * t88;
t1 = [0, 0, 0, 0, 0, 0, -t77, t57, -g(3), -t89, 0, 0, 0, 0, 0, 0, -g(3) * t71 - t77 * t74, -g(3) * t74 + t77 * t71, -t57, -g(1) * (t75 * pkin(1) + t72 * pkin(6)) - g(2) * (t72 * pkin(1) - t75 * pkin(6)) - t89, 0, 0, 0, 0, 0, 0, -t77 * t65 - t88, t56, -t57, -g(1) * t79 - g(2) * t80 - g(3) * t87, 0, 0, 0, 0, 0, 0, t55, t54, -t56, -g(1) * (t78 * t75 + t79) - g(2) * (t78 * t72 + t80) - g(3) * (t64 * pkin(3) - t65 * pkin(8) + t87), 0, 0, 0, 0, 0, 0, t55, t54, -t56, -g(1) * (t90 * t75 + t59) - g(2) * (-pkin(4) * t82 + t80) - g(3) * (t64 * t61 + t65 * t69 + t87) + (-g(1) * (pkin(4) * t70 - t76) - g(2) * t90) * t72;];
U_reg = t1;
