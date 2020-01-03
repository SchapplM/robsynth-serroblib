% Calculate inertial parameters regressor of potential energy for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:38
% EndTime: 2019-12-31 21:11:38
% DurationCPUTime: 0.08s
% Computational Cost: add. (119->43), mult. (120->51), div. (0->0), fcn. (113->8), ass. (0->28)
t66 = qJ(1) + qJ(2);
t61 = sin(t66);
t68 = sin(qJ(3));
t82 = qJ(4) * t68;
t71 = cos(qJ(3));
t85 = t61 * t71;
t89 = pkin(3) * t85 + t61 * t82;
t62 = cos(t66);
t78 = g(1) * t62 + g(2) * t61;
t73 = pkin(6) + pkin(5);
t86 = g(3) * t73;
t84 = t62 * t71;
t69 = sin(qJ(1));
t83 = t69 * pkin(1) + t61 * pkin(2);
t72 = cos(qJ(1));
t81 = t72 * pkin(1) + t62 * pkin(2) + t61 * pkin(7);
t80 = -t62 * pkin(7) + t83;
t79 = pkin(3) * t84 + t62 * t82 + t81;
t77 = -g(1) * t72 - g(2) * t69;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t76 = t71 * t67 - t68 * t70;
t75 = t68 * t67 + t71 * t70;
t74 = t68 * pkin(3) - t71 * qJ(4) + t73;
t52 = g(1) * t61 - g(2) * t62;
t51 = -g(3) * t68 - t78 * t71;
t50 = -g(3) * t71 + t78 * t68;
t1 = [0, 0, 0, 0, 0, 0, t77, g(1) * t69 - g(2) * t72, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t78, t52, -g(3), t77 * pkin(1) - t86, 0, 0, 0, 0, 0, 0, t51, t50, -t52, -g(1) * t81 - g(2) * t80 - t86, 0, 0, 0, 0, 0, 0, t51, -t52, -t50, -g(1) * t79 - g(2) * (t80 + t89) - g(3) * t74, 0, 0, 0, 0, 0, 0, g(3) * t76 - t78 * t75, g(3) * t75 + t78 * t76, t52, -g(1) * (pkin(4) * t84 - t61 * pkin(8) + t79) - g(2) * (pkin(4) * t85 + (-pkin(7) + pkin(8)) * t62 + t83 + t89) - g(3) * (t68 * pkin(4) + t74);];
U_reg = t1;
