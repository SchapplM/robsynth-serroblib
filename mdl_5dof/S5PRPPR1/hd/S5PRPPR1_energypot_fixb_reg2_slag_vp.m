% Calculate inertial parameters regressor of potential energy for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:22:08
% EndTime: 2019-12-05 15:22:08
% DurationCPUTime: 0.13s
% Computational Cost: add. (138->59), mult. (130->79), div. (0->0), fcn. (125->10), ass. (0->32)
t57 = sin(pkin(8));
t77 = g(3) * t57;
t63 = pkin(5) + qJ(1);
t76 = g(3) * t63;
t55 = pkin(7) + qJ(2);
t48 = sin(t55);
t56 = sin(pkin(9));
t75 = t48 * t56;
t60 = cos(pkin(8));
t74 = t48 * t60;
t50 = cos(t55);
t73 = t50 * t60;
t72 = t56 * t60;
t62 = -pkin(6) - qJ(4);
t71 = t57 * t62;
t59 = cos(pkin(9));
t70 = t59 * t60;
t58 = sin(pkin(7));
t69 = t58 * pkin(1) + t48 * pkin(2);
t61 = cos(pkin(7));
t68 = t61 * pkin(1) + t50 * pkin(2) + t48 * qJ(3);
t67 = -t50 * qJ(3) + t69;
t66 = g(1) * t50 + g(2) * t48;
t65 = -g(1) * t61 - g(2) * t58;
t64 = pkin(3) * t60 + qJ(4) * t57;
t54 = pkin(9) + qJ(5);
t49 = cos(t54);
t47 = sin(t54);
t46 = t59 * pkin(4) + pkin(3);
t42 = g(1) * t48 - g(2) * t50;
t41 = -g(3) * t60 + t66 * t57;
t1 = [0, 0, 0, 0, 0, 0, t65, g(1) * t58 - g(2) * t61, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t66, t42, -g(3), t65 * pkin(1) - t76, 0, 0, 0, 0, 0, 0, -t66 * t60 - t77, t41, -t42, -g(1) * t68 - g(2) * t67 - t76, 0, 0, 0, 0, 0, 0, -g(1) * (t50 * t70 + t75) - g(2) * (t48 * t70 - t50 * t56) - t59 * t77, -g(1) * (t48 * t59 - t50 * t72) - g(2) * (-t48 * t72 - t50 * t59) + t56 * t77, -t41, -g(1) * (t64 * t50 + t68) - g(2) * (t64 * t48 + t67) - g(3) * (t57 * pkin(3) - t60 * qJ(4) + t63), 0, 0, 0, 0, 0, 0, -g(1) * (t48 * t47 + t49 * t73) - g(2) * (-t50 * t47 + t49 * t74) - t49 * t77, -g(1) * (-t47 * t73 + t48 * t49) - g(2) * (-t47 * t74 - t50 * t49) + t47 * t77, -t41, -g(1) * (pkin(4) * t75 + t68) - g(2) * (t46 * t74 - t48 * t71 + t69) - g(3) * (t57 * t46 + t60 * t62 + t63) + (-g(1) * (t46 * t60 - t71) - g(2) * (-pkin(4) * t56 - qJ(3))) * t50;];
U_reg = t1;
