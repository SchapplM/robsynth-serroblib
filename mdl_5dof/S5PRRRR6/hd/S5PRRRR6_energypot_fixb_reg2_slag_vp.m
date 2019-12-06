% Calculate inertial parameters regressor of potential energy for
% S5PRRRR6
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:11
% EndTime: 2019-12-05 17:10:12
% DurationCPUTime: 0.14s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t62 = cos(qJ(4));
t48 = t62 * pkin(4) + pkin(3);
t57 = qJ(2) + qJ(3);
t51 = sin(t57);
t53 = cos(t57);
t64 = -pkin(8) - pkin(7);
t83 = t48 * t53 - t51 * t64;
t82 = g(3) * t51;
t81 = g(3) * qJ(1);
t56 = qJ(4) + qJ(5);
t50 = sin(t56);
t58 = sin(pkin(9));
t78 = t58 * t50;
t52 = cos(t56);
t77 = t58 * t52;
t60 = sin(qJ(4));
t76 = t58 * t60;
t75 = t58 * t62;
t59 = cos(pkin(9));
t74 = t59 * t50;
t73 = t59 * t52;
t72 = t59 * t60;
t71 = t59 * t62;
t63 = cos(qJ(2));
t49 = t63 * pkin(2) + pkin(1);
t65 = -pkin(6) - pkin(5);
t70 = t58 * t49 + t59 * t65;
t61 = sin(qJ(2));
t69 = t61 * pkin(2) + qJ(1);
t45 = t59 * t49;
t68 = -t58 * t65 + t45;
t67 = pkin(3) * t53 + pkin(7) * t51;
t66 = g(1) * t59 + g(2) * t58;
t43 = g(1) * t58 - g(2) * t59;
t42 = -g(3) * t53 + t66 * t51;
t1 = [0, 0, 0, 0, 0, 0, -t66, t43, -g(3), -t81, 0, 0, 0, 0, 0, 0, -g(3) * t61 - t66 * t63, -g(3) * t63 + t66 * t61, -t43, -g(1) * (t59 * pkin(1) + t58 * pkin(5)) - g(2) * (t58 * pkin(1) - t59 * pkin(5)) - t81, 0, 0, 0, 0, 0, 0, -t66 * t53 - t82, t42, -t43, -g(1) * t68 - g(2) * t70 - g(3) * t69, 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t71 + t76) - g(2) * (t53 * t75 - t72) - t62 * t82, -g(1) * (-t53 * t72 + t75) - g(2) * (-t53 * t76 - t71) + t60 * t82, -t42, -g(1) * (t67 * t59 + t68) - g(2) * (t67 * t58 + t70) - g(3) * (t51 * pkin(3) - t53 * pkin(7) + t69), 0, 0, 0, 0, 0, 0, -g(1) * (t53 * t73 + t78) - g(2) * (t53 * t77 - t74) - t52 * t82, -g(1) * (-t53 * t74 + t77) - g(2) * (-t53 * t78 - t73) + t50 * t82, -t42, -g(1) * (t83 * t59 + t45) - g(2) * (-pkin(4) * t72 + t70) - g(3) * (t51 * t48 + t53 * t64 + t69) + (-g(1) * (pkin(4) * t60 - t65) - g(2) * t83) * t58;];
U_reg = t1;
