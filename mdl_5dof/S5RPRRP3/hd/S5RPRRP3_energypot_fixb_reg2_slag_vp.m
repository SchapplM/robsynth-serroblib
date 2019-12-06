% Calculate inertial parameters regressor of potential energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:53
% EndTime: 2019-12-05 18:03:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (110->42), mult. (86->45), div. (0->0), fcn. (73->8), ass. (0->25)
t64 = -pkin(7) - pkin(6);
t59 = qJ(2) + pkin(5);
t69 = g(1) * t59;
t61 = sin(qJ(1));
t68 = t61 * pkin(1);
t62 = cos(qJ(3));
t48 = t62 * pkin(3) + pkin(2);
t60 = sin(qJ(3));
t67 = t60 * pkin(3) + t59;
t57 = qJ(1) + pkin(8);
t49 = sin(t57);
t50 = cos(t57);
t66 = g(2) * t49 - g(3) * t50;
t63 = cos(qJ(1));
t65 = g(2) * t61 - g(3) * t63;
t58 = qJ(3) + qJ(4);
t56 = -qJ(5) + t64;
t55 = t63 * pkin(1);
t52 = cos(t58);
t51 = sin(t58);
t47 = pkin(4) * t52 + t48;
t46 = g(2) * t50 + g(3) * t49;
t45 = -g(1) * t52 - t66 * t51;
t44 = -g(1) * t51 + t66 * t52;
t1 = [0, 0, 0, 0, 0, 0, t65, g(2) * t63 + g(3) * t61, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t66, t46, -g(1), t65 * pkin(1) - t69, 0, 0, 0, 0, 0, 0, -g(1) * t60 + t66 * t62, -g(1) * t62 - t66 * t60, -t46, -t69 - g(2) * (-t49 * pkin(2) + t50 * pkin(6) - t68) - g(3) * (t50 * pkin(2) + t49 * pkin(6) + t55), 0, 0, 0, 0, 0, 0, t44, t45, -t46, -g(1) * t67 - g(2) * (-t49 * t48 - t50 * t64 - t68) - g(3) * (t50 * t48 - t49 * t64 + t55), 0, 0, 0, 0, 0, 0, t44, t45, -t46, -g(1) * (pkin(4) * t51 + t67) - g(2) * (-t49 * t47 - t50 * t56 - t68) - g(3) * (t50 * t47 - t49 * t56 + t55);];
U_reg = t1;
