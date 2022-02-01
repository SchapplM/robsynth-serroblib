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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:17:17
% EndTime: 2022-01-20 11:17:17
% DurationCPUTime: 0.12s
% Computational Cost: add. (138->59), mult. (130->79), div. (0->0), fcn. (125->10), ass. (0->32)
t74 = pkin(6) + pkin(5);
t67 = sin(pkin(9));
t88 = g(3) * t67;
t87 = g(3) * t74;
t66 = qJ(1) + qJ(2);
t59 = sin(t66);
t68 = cos(pkin(9));
t86 = t59 * t68;
t69 = sin(qJ(4));
t85 = t59 * t69;
t61 = cos(t66);
t84 = t61 * t68;
t73 = -pkin(8) - pkin(7);
t83 = t67 * t73;
t82 = t68 * t69;
t71 = cos(qJ(4));
t81 = t68 * t71;
t70 = sin(qJ(1));
t80 = t70 * pkin(1) + t59 * pkin(2);
t72 = cos(qJ(1));
t79 = t72 * pkin(1) + t61 * pkin(2) + t59 * qJ(3);
t78 = -t61 * qJ(3) + t80;
t77 = pkin(3) * t68 + pkin(7) * t67;
t76 = g(1) * t61 + g(2) * t59;
t75 = -g(1) * t72 - g(2) * t70;
t65 = qJ(4) + qJ(5);
t60 = cos(t65);
t58 = sin(t65);
t57 = t71 * pkin(4) + pkin(3);
t53 = g(1) * t59 - g(2) * t61;
t52 = -g(3) * t68 + t76 * t67;
t1 = [0, 0, 0, 0, 0, 0, t75, g(1) * t70 - g(2) * t72, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t76, t53, -g(3), t75 * pkin(1) - t87, 0, 0, 0, 0, 0, 0, -t76 * t68 - t88, t52, -t53, -g(1) * t79 - g(2) * t78 - t87, 0, 0, 0, 0, 0, 0, -g(1) * (t61 * t81 + t85) - g(2) * (t59 * t81 - t61 * t69) - t71 * t88, -g(1) * (t59 * t71 - t61 * t82) - g(2) * (-t59 * t82 - t61 * t71) + t69 * t88, -t52, -g(1) * (t77 * t61 + t79) - g(2) * (t77 * t59 + t78) - g(3) * (t67 * pkin(3) - t68 * pkin(7) + t74), 0, 0, 0, 0, 0, 0, -g(1) * (t59 * t58 + t60 * t84) - g(2) * (-t61 * t58 + t60 * t86) - t60 * t88, -g(1) * (-t58 * t84 + t59 * t60) - g(2) * (-t58 * t86 - t61 * t60) + t58 * t88, -t52, -g(1) * (pkin(4) * t85 + t79) - g(2) * (t57 * t86 - t59 * t83 + t80) - g(3) * (t67 * t57 + t68 * t73 + t74) + (-g(1) * (t57 * t68 - t83) - g(2) * (-pkin(4) * t69 - qJ(3))) * t61;];
U_reg = t1;
