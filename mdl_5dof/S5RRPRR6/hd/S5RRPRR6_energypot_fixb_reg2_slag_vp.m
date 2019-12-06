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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:36:09
% EndTime: 2019-12-05 18:36:10
% DurationCPUTime: 0.11s
% Computational Cost: add. (138->58), mult. (130->79), div. (0->0), fcn. (125->10), ass. (0->32)
t72 = pkin(6) + pkin(5);
t65 = sin(pkin(9));
t86 = g(1) * t65;
t85 = g(1) * t72;
t64 = qJ(1) + qJ(2);
t59 = sin(t64);
t84 = g(2) * t59;
t66 = cos(pkin(9));
t83 = t59 * t66;
t61 = cos(t64);
t82 = t61 * t66;
t67 = sin(qJ(4));
t81 = t61 * t67;
t71 = -pkin(8) - pkin(7);
t80 = t65 * t71;
t79 = t66 * t67;
t69 = cos(qJ(4));
t78 = t66 * t69;
t70 = cos(qJ(1));
t77 = t70 * pkin(1) + t61 * pkin(2) + t59 * qJ(3);
t68 = sin(qJ(1));
t76 = -t68 * pkin(1) + t61 * qJ(3);
t75 = pkin(3) * t66 + pkin(7) * t65;
t74 = -g(3) * t61 + t84;
t73 = g(2) * t68 - g(3) * t70;
t63 = qJ(4) + qJ(5);
t60 = cos(t63);
t58 = sin(t63);
t57 = t69 * pkin(4) + pkin(3);
t53 = g(2) * t61 + g(3) * t59;
t52 = g(1) * t66 + t74 * t65;
t1 = [0, 0, 0, 0, 0, 0, t73, g(2) * t70 + g(3) * t68, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t74, t53, -g(1), t73 * pkin(1) - t85, 0, 0, 0, 0, 0, 0, t74 * t66 - t86, -t52, -t53, -t85 - g(2) * (-t59 * pkin(2) + t76) - g(3) * t77, 0, 0, 0, 0, 0, 0, -t69 * t86 - g(2) * (-t59 * t78 + t81) - g(3) * (t59 * t67 + t61 * t78), t67 * t86 - g(2) * (t59 * t79 + t61 * t69) - g(3) * (t59 * t69 - t61 * t79), t52, -g(1) * (t65 * pkin(3) - t66 * pkin(7) + t72) - g(2) * t76 - g(3) * (t75 * t61 + t77) - (-pkin(2) - t75) * t84, 0, 0, 0, 0, 0, 0, -t60 * t86 - g(2) * (t61 * t58 - t60 * t83) - g(3) * (t59 * t58 + t60 * t82), t58 * t86 - g(2) * (t58 * t83 + t61 * t60) - g(3) * (-t58 * t82 + t59 * t60), t52, -g(1) * (t65 * t57 + t66 * t71 + t72) - g(2) * (pkin(4) * t81 + t76) - g(3) * (t57 * t82 - t61 * t80 + t77) + (-g(2) * (-t57 * t66 - pkin(2) + t80) - g(3) * pkin(4) * t67) * t59;];
U_reg = t1;
