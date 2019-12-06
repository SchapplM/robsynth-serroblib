% Calculate inertial parameters regressor of potential energy for
% S5RPRPR3
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:31
% EndTime: 2019-12-05 17:51:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (134->45), mult. (94->57), div. (0->0), fcn. (85->10), ass. (0->29)
t71 = qJ(2) + pkin(5);
t56 = pkin(6) + t71;
t76 = g(1) * t56;
t58 = sin(pkin(9));
t75 = g(1) * t58;
t57 = qJ(1) + pkin(8);
t54 = qJ(3) + t57;
t50 = sin(t54);
t74 = g(2) * t50;
t59 = cos(pkin(9));
t60 = sin(qJ(5));
t73 = t59 * t60;
t62 = cos(qJ(5));
t72 = t59 * t62;
t53 = cos(t57);
t63 = cos(qJ(1));
t70 = t63 * pkin(1) + pkin(2) * t53;
t51 = cos(t54);
t69 = t51 * pkin(3) + t50 * qJ(4) + t70;
t52 = sin(t57);
t61 = sin(qJ(1));
t68 = -t61 * pkin(1) - pkin(2) * t52;
t67 = pkin(4) * t59 + pkin(7) * t58;
t66 = -g(3) * t51 + t74;
t65 = g(2) * t61 - g(3) * t63;
t64 = t51 * qJ(4) + t68;
t45 = g(2) * t51 + g(3) * t50;
t44 = g(1) * t59 + t66 * t58;
t1 = [0, 0, 0, 0, 0, 0, t65, g(2) * t63 + g(3) * t61, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, g(2) * t52 - g(3) * t53, g(2) * t53 + g(3) * t52, -g(1), t65 * pkin(1) - g(1) * t71, 0, 0, 0, 0, 0, 0, t66, t45, -g(1), -g(2) * t68 - g(3) * t70 - t76, 0, 0, 0, 0, 0, 0, t66 * t59 - t75, -t44, -t45, -t76 - g(2) * (-t50 * pkin(3) + t64) - g(3) * t69, 0, 0, 0, 0, 0, 0, -t62 * t75 - g(2) * (-t50 * t72 + t51 * t60) - g(3) * (t50 * t60 + t51 * t72), t60 * t75 - g(2) * (t50 * t73 + t51 * t62) - g(3) * (t50 * t62 - t51 * t73), t44, -g(1) * (t58 * pkin(4) - t59 * pkin(7) + t56) - g(2) * t64 - g(3) * (t67 * t51 + t69) - (-pkin(3) - t67) * t74;];
U_reg = t1;
