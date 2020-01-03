% Calculate inertial parameters regressor of potential energy for
% S5RRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:03:57
% EndTime: 2020-01-03 12:03:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->44), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t70 = pkin(6) + pkin(5);
t75 = g(1) * t70;
t69 = cos(qJ(1));
t74 = t69 * pkin(1);
t66 = cos(pkin(9));
t53 = t66 * pkin(3) + pkin(2);
t67 = -pkin(7) - qJ(3);
t63 = pkin(9) + qJ(4);
t65 = sin(pkin(9));
t73 = t65 * pkin(3) + t70;
t64 = qJ(1) + qJ(2);
t57 = sin(t64);
t58 = cos(t64);
t72 = g(2) * t57 - g(3) * t58;
t68 = sin(qJ(1));
t71 = -g(2) * t68 + g(3) * t69;
t62 = -pkin(8) + t67;
t61 = t68 * pkin(1);
t56 = qJ(5) + t63;
t55 = cos(t63);
t54 = sin(t63);
t52 = cos(t56);
t51 = sin(t56);
t50 = pkin(4) * t55 + t53;
t49 = g(2) * t58 + g(3) * t57;
t1 = [0, 0, 0, 0, 0, 0, t71, -g(2) * t69 - g(3) * t68, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -t72, -t49, -g(1), t71 * pkin(1) - t75, 0, 0, 0, 0, 0, 0, -g(1) * t65 - t72 * t66, -g(1) * t66 + t72 * t65, t49, -t75 - g(2) * (t57 * pkin(2) - t58 * qJ(3) + t61) - g(3) * (-t58 * pkin(2) - t57 * qJ(3) - t74), 0, 0, 0, 0, 0, 0, -g(1) * t54 - t72 * t55, -g(1) * t55 + t72 * t54, t49, -g(1) * t73 - g(2) * (t57 * t53 + t58 * t67 + t61) - g(3) * (-t58 * t53 + t57 * t67 - t74), 0, 0, 0, 0, 0, 0, -g(1) * t51 - t72 * t52, -g(1) * t52 + t72 * t51, t49, -g(1) * (pkin(4) * t54 + t73) - g(2) * (t57 * t50 + t58 * t62 + t61) - g(3) * (-t58 * t50 + t57 * t62 - t74);];
U_reg = t1;
