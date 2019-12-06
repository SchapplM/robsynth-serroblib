% Calculate inertial parameters regressor of potential energy for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:40
% EndTime: 2019-12-05 15:54:40
% DurationCPUTime: 0.17s
% Computational Cost: add. (132->75), mult. (169->98), div. (0->0), fcn. (168->10), ass. (0->33)
t66 = sin(qJ(2));
t79 = g(3) * t66;
t61 = sin(pkin(9));
t78 = t61 * pkin(3);
t63 = cos(pkin(9));
t50 = t63 * pkin(3) + pkin(2);
t77 = g(3) * qJ(1);
t65 = -pkin(6) - qJ(3);
t59 = -pkin(7) + t65;
t76 = t59 * t66;
t62 = sin(pkin(8));
t75 = t62 * t61;
t67 = cos(qJ(2));
t74 = t62 * t67;
t64 = cos(pkin(8));
t73 = t64 * t67;
t72 = t65 * t66;
t71 = t64 * pkin(1) + t62 * pkin(5);
t60 = pkin(9) + qJ(4);
t55 = t62 * pkin(1);
t70 = -t64 * pkin(5) + t55;
t69 = g(1) * t64 + g(2) * t62;
t68 = pkin(2) * t67 + qJ(3) * t66;
t53 = qJ(5) + t60;
t52 = cos(t60);
t51 = sin(t60);
t49 = cos(t53);
t48 = sin(t53);
t47 = g(1) * t62 - g(2) * t64;
t46 = pkin(4) * t51 + t78;
t45 = pkin(4) * t52 + t50;
t44 = -g(3) * t67 + t69 * t66;
t1 = [0, 0, 0, 0, 0, 0, -t69, t47, -g(3), -t77, 0, 0, 0, 0, 0, 0, -t69 * t67 - t79, t44, -t47, -g(1) * t71 - g(2) * t70 - t77, 0, 0, 0, 0, 0, 0, -g(1) * (t63 * t73 + t75) - g(2) * (-t64 * t61 + t63 * t74) - t63 * t79, -g(1) * (-t61 * t73 + t62 * t63) - g(2) * (-t61 * t74 - t64 * t63) + t61 * t79, -t44, -g(1) * (t68 * t64 + t71) - g(2) * (t68 * t62 + t70) - g(3) * (t66 * pkin(2) - t67 * qJ(3) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t51 + t52 * t73) - g(2) * (-t64 * t51 + t52 * t74) - t52 * t79, -g(1) * (-t51 * t73 + t62 * t52) - g(2) * (-t51 * t74 - t64 * t52) + t51 * t79, -t44, -g(1) * (pkin(3) * t75 + t71) - g(2) * (t50 * t74 - t62 * t72 + t55) - g(3) * (t66 * t50 + t67 * t65 + qJ(1)) + (-g(1) * (t50 * t67 - t72) - g(2) * (-pkin(5) - t78)) * t64, 0, 0, 0, 0, 0, 0, -g(1) * (t62 * t48 + t49 * t73) - g(2) * (-t64 * t48 + t49 * t74) - t49 * t79, -g(1) * (-t48 * t73 + t62 * t49) - g(2) * (-t48 * t74 - t64 * t49) + t48 * t79, -t44, -g(1) * (t62 * t46 + t71) - g(2) * (t45 * t74 - t62 * t76 + t55) - g(3) * (t66 * t45 + t67 * t59 + qJ(1)) + (-g(1) * (t45 * t67 - t76) - g(2) * (-pkin(5) - t46)) * t64;];
U_reg = t1;
