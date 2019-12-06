% Calculate inertial parameters regressor of potential energy for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:47
% EndTime: 2019-12-05 15:14:47
% DurationCPUTime: 0.14s
% Computational Cost: add. (131->61), mult. (143->79), div. (0->0), fcn. (138->10), ass. (0->36)
t62 = cos(qJ(4));
t47 = t62 * pkin(4) + pkin(3);
t54 = pkin(9) + qJ(3);
t48 = sin(t54);
t49 = cos(t54);
t63 = -pkin(7) - pkin(6);
t81 = t47 * t49 - t48 * t63;
t80 = g(3) * t48;
t79 = g(3) * qJ(1);
t55 = qJ(4) + qJ(5);
t50 = sin(t55);
t57 = sin(pkin(8));
t76 = t57 * t50;
t51 = cos(t55);
t75 = t57 * t51;
t61 = sin(qJ(4));
t74 = t57 * t61;
t73 = t57 * t62;
t59 = cos(pkin(8));
t72 = t59 * t50;
t71 = t59 * t51;
t70 = t59 * t61;
t69 = t59 * t62;
t58 = cos(pkin(9));
t46 = t58 * pkin(2) + pkin(1);
t60 = -pkin(5) - qJ(2);
t68 = t57 * t46 + t59 * t60;
t56 = sin(pkin(9));
t67 = t56 * pkin(2) + qJ(1);
t43 = t59 * t46;
t66 = -t57 * t60 + t43;
t65 = pkin(3) * t49 + pkin(6) * t48;
t64 = g(1) * t59 + g(2) * t57;
t41 = g(1) * t57 - g(2) * t59;
t40 = -g(3) * t49 + t64 * t48;
t1 = [0, 0, 0, 0, 0, 0, -t64, t41, -g(3), -t79, 0, 0, 0, 0, 0, 0, -g(3) * t56 - t64 * t58, -g(3) * t58 + t64 * t56, -t41, -g(1) * (t59 * pkin(1) + t57 * qJ(2)) - g(2) * (t57 * pkin(1) - t59 * qJ(2)) - t79, 0, 0, 0, 0, 0, 0, -t64 * t49 - t80, t40, -t41, -g(1) * t66 - g(2) * t68 - g(3) * t67, 0, 0, 0, 0, 0, 0, -g(1) * (t49 * t69 + t74) - g(2) * (t49 * t73 - t70) - t62 * t80, -g(1) * (-t49 * t70 + t73) - g(2) * (-t49 * t74 - t69) + t61 * t80, -t40, -g(1) * (t65 * t59 + t66) - g(2) * (t65 * t57 + t68) - g(3) * (t48 * pkin(3) - t49 * pkin(6) + t67), 0, 0, 0, 0, 0, 0, -g(1) * (t49 * t71 + t76) - g(2) * (t49 * t75 - t72) - t51 * t80, -g(1) * (-t49 * t72 + t75) - g(2) * (-t49 * t76 - t71) + t50 * t80, -t40, -g(1) * (t81 * t59 + t43) - g(2) * (-pkin(4) * t70 + t68) - g(3) * (t48 * t47 + t49 * t63 + t67) + (-g(1) * (pkin(4) * t61 - t60) - g(2) * t81) * t57;];
U_reg = t1;
