% Calculate inertial parameters regressor of potential energy for
% S5PRPRR1
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:03
% EndTime: 2019-12-05 15:43:03
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->46), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t64 = pkin(5) + qJ(1);
t68 = g(3) * t64;
t61 = cos(pkin(9));
t46 = t61 * pkin(3) + pkin(2);
t63 = -pkin(6) - qJ(3);
t57 = pkin(9) + qJ(4);
t59 = sin(pkin(9));
t67 = t59 * pkin(3) + t64;
t58 = pkin(8) + qJ(2);
t48 = sin(t58);
t50 = cos(t58);
t66 = g(1) * t50 + g(2) * t48;
t60 = sin(pkin(8));
t62 = cos(pkin(8));
t65 = -g(1) * t62 - g(2) * t60;
t56 = -pkin(7) + t63;
t55 = t62 * pkin(1);
t53 = t60 * pkin(1);
t51 = qJ(5) + t57;
t49 = cos(t57);
t47 = sin(t57);
t45 = cos(t51);
t44 = sin(t51);
t42 = pkin(4) * t49 + t46;
t41 = g(1) * t48 - g(2) * t50;
t1 = [0, 0, 0, 0, 0, 0, t65, g(1) * t60 - g(2) * t62, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t66, t41, -g(3), t65 * pkin(1) - t68, 0, 0, 0, 0, 0, 0, -g(3) * t59 - t66 * t61, -g(3) * t61 + t66 * t59, -t41, -g(1) * (t50 * pkin(2) + t48 * qJ(3) + t55) - g(2) * (t48 * pkin(2) - t50 * qJ(3) + t53) - t68, 0, 0, 0, 0, 0, 0, -g(3) * t47 - t66 * t49, -g(3) * t49 + t66 * t47, -t41, -g(1) * (t50 * t46 - t48 * t63 + t55) - g(2) * (t48 * t46 + t50 * t63 + t53) - g(3) * t67, 0, 0, 0, 0, 0, 0, -g(3) * t44 - t66 * t45, -g(3) * t45 + t66 * t44, -t41, -g(1) * (t50 * t42 - t48 * t56 + t55) - g(2) * (t48 * t42 + t50 * t56 + t53) - g(3) * (pkin(4) * t47 + t67);];
U_reg = t1;
