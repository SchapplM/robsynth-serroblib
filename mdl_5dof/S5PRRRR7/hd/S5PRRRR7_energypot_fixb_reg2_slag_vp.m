% Calculate inertial parameters regressor of potential energy for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:58
% EndTime: 2019-12-05 17:12:58
% DurationCPUTime: 0.15s
% Computational Cost: add. (132->75), mult. (169->100), div. (0->0), fcn. (168->10), ass. (0->35)
t70 = -pkin(7) - pkin(6);
t67 = sin(qJ(2));
t84 = g(3) * t67;
t66 = sin(qJ(3));
t83 = t66 * pkin(3);
t68 = cos(qJ(3));
t53 = t68 * pkin(3) + pkin(2);
t82 = g(3) * qJ(1);
t62 = -pkin(8) + t70;
t81 = t62 * t67;
t64 = sin(pkin(9));
t80 = t64 * t66;
t69 = cos(qJ(2));
t79 = t64 * t69;
t65 = cos(pkin(9));
t78 = t65 * t69;
t77 = t66 * t69;
t76 = t67 * t70;
t75 = t68 * t69;
t74 = t65 * pkin(1) + t64 * pkin(5);
t63 = qJ(3) + qJ(4);
t57 = t64 * pkin(1);
t73 = -t65 * pkin(5) + t57;
t72 = pkin(2) * t69 + pkin(6) * t67;
t71 = g(1) * t65 + g(2) * t64;
t59 = qJ(5) + t63;
t55 = cos(t63);
t54 = sin(t63);
t52 = cos(t59);
t51 = sin(t59);
t50 = g(1) * t64 - g(2) * t65;
t49 = pkin(4) * t54 + t83;
t48 = pkin(4) * t55 + t53;
t47 = -g(3) * t69 + t71 * t67;
t1 = [0, 0, 0, 0, 0, 0, -t71, t50, -g(3), -t82, 0, 0, 0, 0, 0, 0, -t71 * t69 - t84, t47, -t50, -g(1) * t74 - g(2) * t73 - t82, 0, 0, 0, 0, 0, 0, -g(1) * (t65 * t75 + t80) - g(2) * (t64 * t75 - t65 * t66) - t68 * t84, -g(1) * (t64 * t68 - t65 * t77) - g(2) * (-t64 * t77 - t65 * t68) + t66 * t84, -t47, -g(1) * (t72 * t65 + t74) - g(2) * (t72 * t64 + t73) - g(3) * (t67 * pkin(2) - t69 * pkin(6) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t54 + t55 * t78) - g(2) * (-t65 * t54 + t55 * t79) - t55 * t84, -g(1) * (-t54 * t78 + t64 * t55) - g(2) * (-t54 * t79 - t65 * t55) + t54 * t84, -t47, -g(1) * (pkin(3) * t80 + t74) - g(2) * (t53 * t79 - t64 * t76 + t57) - g(3) * (t67 * t53 + t69 * t70 + qJ(1)) + (-g(1) * (t53 * t69 - t76) - g(2) * (-pkin(5) - t83)) * t65, 0, 0, 0, 0, 0, 0, -g(1) * (t64 * t51 + t52 * t78) - g(2) * (-t65 * t51 + t52 * t79) - t52 * t84, -g(1) * (-t51 * t78 + t64 * t52) - g(2) * (-t51 * t79 - t65 * t52) + t51 * t84, -t47, -g(1) * (t64 * t49 + t74) - g(2) * (t48 * t79 - t64 * t81 + t57) - g(3) * (t67 * t48 + t69 * t62 + qJ(1)) + (-g(1) * (t48 * t69 - t81) - g(2) * (-pkin(5) - t49)) * t65;];
U_reg = t1;
