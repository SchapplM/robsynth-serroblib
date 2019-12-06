% Calculate inertial parameters regressor of potential energy for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:11
% EndTime: 2019-12-05 18:48:11
% DurationCPUTime: 0.06s
% Computational Cost: add. (110->42), mult. (86->45), div. (0->0), fcn. (73->8), ass. (0->25)
t68 = pkin(6) + pkin(5);
t67 = -pkin(8) - pkin(7);
t73 = g(1) * t68;
t64 = sin(qJ(1));
t72 = t64 * pkin(1);
t65 = cos(qJ(3));
t52 = t65 * pkin(3) + pkin(2);
t63 = sin(qJ(3));
t71 = t63 * pkin(3) + t68;
t62 = qJ(1) + qJ(2);
t54 = sin(t62);
t56 = cos(t62);
t70 = g(2) * t54 - g(3) * t56;
t66 = cos(qJ(1));
t69 = g(2) * t64 - g(3) * t66;
t61 = qJ(3) + qJ(4);
t60 = -qJ(5) + t67;
t59 = t66 * pkin(1);
t55 = cos(t61);
t53 = sin(t61);
t51 = pkin(4) * t55 + t52;
t50 = g(2) * t56 + g(3) * t54;
t49 = -g(1) * t55 - t70 * t53;
t48 = -g(1) * t53 + t70 * t55;
t1 = [0, 0, 0, 0, 0, 0, t69, g(2) * t66 + g(3) * t64, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t70, t50, -g(1), t69 * pkin(1) - t73, 0, 0, 0, 0, 0, 0, -g(1) * t63 + t70 * t65, -g(1) * t65 - t70 * t63, -t50, -t73 - g(2) * (-t54 * pkin(2) + t56 * pkin(7) - t72) - g(3) * (t56 * pkin(2) + t54 * pkin(7) + t59), 0, 0, 0, 0, 0, 0, t48, t49, -t50, -g(1) * t71 - g(2) * (-t54 * t52 - t56 * t67 - t72) - g(3) * (t56 * t52 - t54 * t67 + t59), 0, 0, 0, 0, 0, 0, t48, t49, -t50, -g(1) * (pkin(4) * t53 + t71) - g(2) * (-t54 * t51 - t56 * t60 - t72) - g(3) * (t56 * t51 - t54 * t60 + t59);];
U_reg = t1;
