% Calculate inertial parameters regressor of potential energy for
% S5RRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:15:18
% EndTime: 2020-01-03 12:15:18
% DurationCPUTime: 0.09s
% Computational Cost: add. (116->44), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t72 = pkin(6) + pkin(5);
t71 = -pkin(8) - pkin(7);
t77 = g(1) * t72;
t70 = cos(qJ(1));
t76 = t70 * pkin(1);
t69 = cos(qJ(3));
t55 = t69 * pkin(3) + pkin(2);
t65 = qJ(3) + qJ(4);
t67 = sin(qJ(3));
t75 = t67 * pkin(3) + t72;
t66 = qJ(1) + qJ(2);
t57 = sin(t66);
t59 = cos(t66);
t74 = g(2) * t57 - g(3) * t59;
t68 = sin(qJ(1));
t73 = -g(2) * t68 + g(3) * t70;
t64 = -pkin(9) + t71;
t62 = t68 * pkin(1);
t60 = qJ(5) + t65;
t58 = cos(t65);
t56 = sin(t65);
t54 = cos(t60);
t53 = sin(t60);
t52 = pkin(4) * t58 + t55;
t51 = g(2) * t59 + g(3) * t57;
t1 = [0, 0, 0, 0, 0, 0, t73, -g(2) * t70 - g(3) * t68, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -t74, -t51, -g(1), t73 * pkin(1) - t77, 0, 0, 0, 0, 0, 0, -g(1) * t67 - t74 * t69, -g(1) * t69 + t74 * t67, t51, -t77 - g(2) * (t57 * pkin(2) - t59 * pkin(7) + t62) - g(3) * (-t59 * pkin(2) - t57 * pkin(7) - t76), 0, 0, 0, 0, 0, 0, -g(1) * t56 - t74 * t58, -g(1) * t58 + t74 * t56, t51, -g(1) * t75 - g(2) * (t57 * t55 + t59 * t71 + t62) - g(3) * (-t59 * t55 + t57 * t71 - t76), 0, 0, 0, 0, 0, 0, -g(1) * t53 - t74 * t54, -g(1) * t54 + t74 * t53, t51, -g(1) * (pkin(4) * t56 + t75) - g(2) * (t57 * t52 + t59 * t64 + t62) - g(3) * (-t59 * t52 + t57 * t64 - t76);];
U_reg = t1;
