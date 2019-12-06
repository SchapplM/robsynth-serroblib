% Calculate inertial parameters regressor of potential energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:42:04
% EndTime: 2019-12-05 17:42:04
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->45), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t63 = qJ(2) + pkin(5);
t71 = g(1) * t63;
t65 = sin(qJ(1));
t70 = t65 * pkin(1);
t62 = cos(pkin(9));
t49 = t62 * pkin(3) + pkin(2);
t64 = -pkin(6) - qJ(3);
t59 = pkin(9) + qJ(4);
t61 = sin(pkin(9));
t69 = t61 * pkin(3) + t63;
t60 = qJ(1) + pkin(8);
t51 = sin(t60);
t53 = cos(t60);
t68 = g(2) * t51 - g(3) * t53;
t66 = cos(qJ(1));
t67 = g(2) * t65 - g(3) * t66;
t58 = -pkin(7) + t64;
t57 = t66 * pkin(1);
t54 = qJ(5) + t59;
t52 = cos(t59);
t50 = sin(t59);
t48 = cos(t54);
t47 = sin(t54);
t46 = pkin(4) * t52 + t49;
t45 = g(2) * t53 + g(3) * t51;
t1 = [0, 0, 0, 0, 0, 0, t67, g(2) * t66 + g(3) * t65, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t68, t45, -g(1), t67 * pkin(1) - t71, 0, 0, 0, 0, 0, 0, -g(1) * t61 + t68 * t62, -g(1) * t62 - t68 * t61, -t45, -t71 - g(2) * (-t51 * pkin(2) + t53 * qJ(3) - t70) - g(3) * (t53 * pkin(2) + t51 * qJ(3) + t57), 0, 0, 0, 0, 0, 0, -g(1) * t50 + t68 * t52, -g(1) * t52 - t68 * t50, -t45, -g(1) * t69 - g(2) * (-t51 * t49 - t53 * t64 - t70) - g(3) * (t53 * t49 - t51 * t64 + t57), 0, 0, 0, 0, 0, 0, -g(1) * t47 + t68 * t48, -g(1) * t48 - t68 * t47, -t45, -g(1) * (pkin(4) * t50 + t69) - g(2) * (-t51 * t46 - t53 * t58 - t70) - g(3) * (t53 * t46 - t51 * t58 + t57);];
U_reg = t1;
