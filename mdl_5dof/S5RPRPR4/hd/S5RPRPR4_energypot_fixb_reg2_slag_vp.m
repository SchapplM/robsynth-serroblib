% Calculate inertial parameters regressor of potential energy for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:39:15
% EndTime: 2020-01-03 11:39:15
% DurationCPUTime: 0.10s
% Computational Cost: add. (116->44), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t65 = qJ(2) + pkin(5);
t74 = g(1) * t65;
t69 = cos(qJ(1));
t73 = t69 * pkin(1);
t68 = cos(qJ(3));
t52 = t68 * pkin(3) + pkin(2);
t64 = -qJ(4) - pkin(6);
t62 = qJ(3) + pkin(9);
t66 = sin(qJ(3));
t72 = t66 * pkin(3) + t65;
t63 = qJ(1) + pkin(8);
t54 = sin(t63);
t56 = cos(t63);
t71 = g(2) * t54 - g(3) * t56;
t67 = sin(qJ(1));
t70 = -g(2) * t67 + g(3) * t69;
t61 = -pkin(7) + t64;
t59 = t67 * pkin(1);
t57 = qJ(5) + t62;
t55 = cos(t62);
t53 = sin(t62);
t51 = cos(t57);
t50 = sin(t57);
t49 = pkin(4) * t55 + t52;
t48 = g(2) * t56 + g(3) * t54;
t1 = [0, 0, 0, 0, 0, 0, t70, -g(2) * t69 - g(3) * t67, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -t71, -t48, -g(1), t70 * pkin(1) - t74, 0, 0, 0, 0, 0, 0, -g(1) * t66 - t71 * t68, -g(1) * t68 + t71 * t66, t48, -t74 - g(2) * (t54 * pkin(2) - t56 * pkin(6) + t59) - g(3) * (-t56 * pkin(2) - t54 * pkin(6) - t73), 0, 0, 0, 0, 0, 0, -g(1) * t53 - t71 * t55, -g(1) * t55 + t71 * t53, t48, -g(1) * t72 - g(2) * (t54 * t52 + t56 * t64 + t59) - g(3) * (-t56 * t52 + t54 * t64 - t73), 0, 0, 0, 0, 0, 0, -g(1) * t50 - t71 * t51, -g(1) * t51 + t71 * t50, t48, -g(1) * (pkin(4) * t53 + t72) - g(2) * (t54 * t49 + t56 * t61 + t59) - g(3) * (-t56 * t49 + t54 * t61 - t73);];
U_reg = t1;
