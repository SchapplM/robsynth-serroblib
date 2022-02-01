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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:49:23
% EndTime: 2022-01-20 11:49:23
% DurationCPUTime: 0.06s
% Computational Cost: add. (110->43), mult. (86->45), div. (0->0), fcn. (73->8), ass. (0->25)
t70 = pkin(6) + pkin(5);
t69 = -pkin(8) - pkin(7);
t74 = g(3) * t70;
t67 = cos(qJ(3));
t53 = t67 * pkin(3) + pkin(2);
t65 = sin(qJ(3));
t73 = t65 * pkin(3) + t70;
t64 = qJ(1) + qJ(2);
t55 = sin(t64);
t57 = cos(t64);
t72 = g(1) * t57 + g(2) * t55;
t66 = sin(qJ(1));
t68 = cos(qJ(1));
t71 = -g(1) * t68 - g(2) * t66;
t63 = qJ(3) + qJ(4);
t62 = qJ(5) - t69;
t61 = t68 * pkin(1);
t59 = t66 * pkin(1);
t56 = cos(t63);
t54 = sin(t63);
t51 = pkin(4) * t56 + t53;
t50 = g(1) * t55 - g(2) * t57;
t49 = -g(3) * t54 - t72 * t56;
t48 = -g(3) * t56 + t72 * t54;
t1 = [0, 0, 0, 0, 0, 0, t71, g(1) * t66 - g(2) * t68, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t72, t50, -g(3), t71 * pkin(1) - t74, 0, 0, 0, 0, 0, 0, -g(3) * t65 - t72 * t67, -g(3) * t67 + t72 * t65, -t50, -g(1) * (t57 * pkin(2) + t55 * pkin(7) + t61) - g(2) * (t55 * pkin(2) - t57 * pkin(7) + t59) - t74, 0, 0, 0, 0, 0, 0, t49, t48, -t50, -g(1) * (t57 * t53 - t55 * t69 + t61) - g(2) * (t55 * t53 + t57 * t69 + t59) - g(3) * t73, 0, 0, 0, 0, 0, 0, t49, t48, -t50, -g(1) * (t57 * t51 + t62 * t55 + t61) - g(2) * (t55 * t51 - t57 * t62 + t59) - g(3) * (pkin(4) * t54 + t73);];
U_reg = t1;
