% Calculate inertial parameters regressor of potential energy for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:53
% EndTime: 2020-01-03 12:05:53
% DurationCPUTime: 0.12s
% Computational Cost: add. (138->60), mult. (130->82), div. (0->0), fcn. (125->10), ass. (0->32)
t70 = pkin(6) + pkin(5);
t63 = sin(pkin(9));
t84 = pkin(7) * t63;
t83 = g(1) * t63;
t82 = g(1) * t70;
t68 = cos(qJ(1));
t81 = t68 * pkin(1);
t62 = qJ(1) + qJ(2);
t57 = sin(t62);
t64 = cos(pkin(9));
t80 = t57 * t64;
t59 = cos(t62);
t79 = t59 * t64;
t65 = sin(qJ(4));
t78 = t64 * t65;
t67 = cos(qJ(4));
t77 = t64 * t67;
t66 = sin(qJ(1));
t76 = t66 * pkin(1) + t57 * pkin(2);
t75 = pkin(4) * t65 + qJ(3);
t74 = -g(2) * t57 + g(3) * t59;
t73 = -g(2) * t66 + g(3) * t68;
t72 = -t57 * qJ(3) - t81;
t55 = t67 * pkin(4) + pkin(3);
t69 = -pkin(8) - pkin(7);
t71 = t55 * t64 - t63 * t69;
t61 = qJ(4) + qJ(5);
t58 = cos(t61);
t56 = sin(t61);
t53 = g(2) * t59 + g(3) * t57;
t52 = g(1) * t64 + t74 * t63;
t1 = [0, 0, 0, 0, 0, 0, t73, -g(2) * t68 - g(3) * t66, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t74, -t53, -g(1), t73 * pkin(1) - t82, 0, 0, 0, 0, 0, 0, t74 * t64 - t83, -t52, t53, -t82 - g(2) * (-t59 * qJ(3) + t76) - g(3) * (-t59 * pkin(2) + t72), 0, 0, 0, 0, 0, 0, -t67 * t83 - g(2) * (t57 * t77 - t59 * t65) - g(3) * (-t57 * t65 - t59 * t77), t65 * t83 - g(2) * (-t57 * t78 - t59 * t67) - g(3) * (-t57 * t67 + t59 * t78), t52, -g(1) * (t63 * pkin(3) - t64 * pkin(7) + t70) - g(2) * (pkin(3) * t80 + t57 * t84 + t76) - g(3) * t72 + (g(2) * qJ(3) - g(3) * (-pkin(3) * t64 - pkin(2) - t84)) * t59, 0, 0, 0, 0, 0, 0, -t58 * t83 - g(2) * (-t59 * t56 + t58 * t80) - g(3) * (-t57 * t56 - t58 * t79), t56 * t83 - g(2) * (-t56 * t80 - t59 * t58) - g(3) * (t56 * t79 - t57 * t58), t52, -g(1) * (t63 * t55 + t64 * t69 + t70) - g(2) * t76 + g(3) * t81 + (-g(2) * t71 + g(3) * t75) * t57 + (g(2) * t75 - g(3) * (-pkin(2) - t71)) * t59;];
U_reg = t1;
