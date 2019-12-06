% Calculate inertial parameters regressor of potential energy for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:42
% EndTime: 2019-12-05 18:53:42
% DurationCPUTime: 0.06s
% Computational Cost: add. (81->28), mult. (85->39), div. (0->0), fcn. (79->10), ass. (0->26)
t64 = cos(qJ(3));
t74 = pkin(2) * t64;
t58 = qJ(3) + qJ(4);
t54 = sin(t58);
t73 = g(3) * t54;
t61 = sin(qJ(3));
t72 = g(3) * t61;
t59 = qJ(1) + qJ(2);
t55 = sin(t59);
t60 = sin(qJ(5));
t71 = t55 * t60;
t63 = cos(qJ(5));
t70 = t55 * t63;
t57 = cos(t59);
t69 = t57 * t60;
t68 = t57 * t63;
t67 = g(1) * t57 + g(2) * t55;
t62 = sin(qJ(1));
t65 = cos(qJ(1));
t66 = -g(1) * t65 - g(2) * t62;
t56 = cos(t58);
t52 = t66 * pkin(1);
t51 = g(1) * t55 - g(2) * t57;
t50 = -g(3) * t56 + t67 * t54;
t49 = -g(1) * (t65 * pkin(1) + t57 * t74) - g(2) * (t62 * pkin(1) + t55 * t74) - pkin(2) * t72;
t1 = [0, 0, 0, 0, 0, 0, t66, g(1) * t62 - g(2) * t65, -g(3), 0, 0, 0, 0, 0, 0, 0, -t67, t51, -g(3), t52, 0, 0, 0, 0, 0, 0, -t67 * t64 - t72, -g(3) * t64 + t67 * t61, -t51, t52, 0, 0, 0, 0, 0, 0, -t67 * t56 - t73, t50, -t51, t49, 0, 0, 0, 0, 0, 0, -g(1) * (t56 * t68 + t71) - g(2) * (t56 * t70 - t69) - t63 * t73, -g(1) * (-t56 * t69 + t70) - g(2) * (-t56 * t71 - t68) + t60 * t73, -t50, t49;];
U_reg = t1;
