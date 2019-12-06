% Calculate inertial parameters regressor of potential energy for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:58:06
% EndTime: 2019-12-05 14:58:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (142->70), mult. (205->97), div. (0->0), fcn. (216->10), ass. (0->40)
t48 = sin(pkin(8));
t74 = g(3) * t48;
t73 = g(3) * qJ(1);
t53 = -pkin(5) - qJ(3);
t72 = t48 * t53;
t54 = sin(qJ(5));
t71 = t48 * t54;
t55 = cos(qJ(5));
t70 = t48 * t55;
t47 = sin(pkin(9));
t49 = sin(pkin(7));
t69 = t49 * t47;
t51 = cos(pkin(8));
t68 = t49 * t51;
t46 = pkin(9) + qJ(4);
t40 = sin(t46);
t52 = cos(pkin(7));
t67 = t52 * t40;
t41 = cos(t46);
t66 = t52 * t41;
t65 = t52 * t47;
t50 = cos(pkin(9));
t64 = t52 * t50;
t63 = t52 * pkin(1) + t49 * qJ(2);
t39 = t50 * pkin(3) + pkin(2);
t62 = t48 * t39 + t51 * t53 + qJ(1);
t43 = t49 * pkin(1);
t61 = -t52 * qJ(2) + t43;
t60 = g(1) * t52 + g(2) * t49;
t59 = pkin(2) * t51 + qJ(3) * t48;
t58 = pkin(3) * t69 + t63 + (t39 * t51 - t72) * t52;
t27 = t40 * t68 + t66;
t29 = -t49 * t41 + t51 * t67;
t57 = g(1) * t29 + g(2) * t27 + t40 * t74;
t56 = -t49 * t72 + t39 * t68 + t43 + (-pkin(3) * t47 - qJ(2)) * t52;
t34 = g(1) * t49 - g(2) * t52;
t31 = -g(3) * t51 + t60 * t48;
t30 = t49 * t40 + t51 * t66;
t28 = t41 * t68 - t67;
t1 = [0, 0, 0, 0, 0, 0, -t60, t34, -g(3), -t73, 0, 0, 0, 0, 0, 0, -t60 * t51 - t74, t31, -t34, -g(1) * t63 - g(2) * t61 - t73, 0, 0, 0, 0, 0, 0, -g(1) * (t51 * t64 + t69) - g(2) * (t50 * t68 - t65) - t50 * t74, -g(1) * (t49 * t50 - t51 * t65) - g(2) * (-t47 * t68 - t64) + t47 * t74, -t31, -g(1) * (t59 * t52 + t63) - g(2) * (t59 * t49 + t61) - g(3) * (t48 * pkin(2) - t51 * qJ(3) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * t30 - g(2) * t28 - t41 * t74, t57, -t31, -g(1) * t58 - g(2) * t56 - g(3) * t62, 0, 0, 0, 0, 0, 0, -g(1) * (t30 * t55 + t52 * t71) - g(2) * (t28 * t55 + t49 * t71) - g(3) * (t41 * t70 - t51 * t54), -g(1) * (-t30 * t54 + t52 * t70) - g(2) * (-t28 * t54 + t49 * t70) - g(3) * (-t41 * t71 - t51 * t55), -t57, -g(1) * (t30 * pkin(4) + t29 * pkin(6) + t58) - g(2) * (t28 * pkin(4) + t27 * pkin(6) + t56) - g(3) * ((pkin(4) * t41 + pkin(6) * t40) * t48 + t62);];
U_reg = t1;
