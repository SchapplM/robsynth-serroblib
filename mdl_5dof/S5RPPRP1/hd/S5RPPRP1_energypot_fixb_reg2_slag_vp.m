% Calculate inertial parameters regressor of potential energy for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:47
% EndTime: 2020-01-03 11:25:47
% DurationCPUTime: 0.12s
% Computational Cost: add. (128->50), mult. (130->65), div. (0->0), fcn. (125->8), ass. (0->29)
t56 = sin(pkin(8));
t57 = cos(pkin(8));
t77 = pkin(3) * t57 + pkin(6) * t56;
t74 = g(1) * t56;
t59 = qJ(2) + pkin(5);
t73 = g(1) * t59;
t63 = cos(qJ(1));
t72 = t63 * pkin(1);
t60 = sin(qJ(4));
t71 = t57 * t60;
t62 = cos(qJ(4));
t70 = t57 * t62;
t55 = qJ(1) + pkin(7);
t52 = sin(t55);
t61 = sin(qJ(1));
t69 = t61 * pkin(1) + t52 * pkin(2);
t68 = pkin(4) * t60 + qJ(3);
t53 = cos(t55);
t67 = -g(2) * t52 + g(3) * t53;
t66 = -g(2) * t61 + g(3) * t63;
t65 = -t52 * qJ(3) - t72;
t51 = t62 * pkin(4) + pkin(3);
t58 = -qJ(5) - pkin(6);
t64 = t51 * t57 - t56 * t58;
t49 = g(2) * t53 + g(3) * t52;
t48 = g(1) * t57 + t67 * t56;
t47 = -t62 * t74 - g(2) * (t52 * t70 - t53 * t60) - g(3) * (-t52 * t60 - t53 * t70);
t46 = t60 * t74 - g(2) * (-t52 * t71 - t53 * t62) - g(3) * (-t52 * t62 + t53 * t71);
t1 = [0, 0, 0, 0, 0, 0, t66, -g(2) * t63 - g(3) * t61, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, t67, -t49, -g(1), t66 * pkin(1) - t73, 0, 0, 0, 0, 0, 0, t67 * t57 - t74, -t48, t49, -t73 - g(2) * (-t53 * qJ(3) + t69) - g(3) * (-t53 * pkin(2) + t65), 0, 0, 0, 0, 0, 0, t47, t46, t48, -g(1) * (t56 * pkin(3) - t57 * pkin(6) + t59) - g(2) * (t77 * t52 + t69) - g(3) * t65 + (g(2) * qJ(3) - g(3) * (-pkin(2) - t77)) * t53, 0, 0, 0, 0, 0, 0, t47, t46, t48, -g(1) * (t56 * t51 + t57 * t58 + t59) - g(2) * t69 + g(3) * t72 + (-g(2) * t64 + g(3) * t68) * t52 + (g(2) * t68 - g(3) * (-pkin(2) - t64)) * t53;];
U_reg = t1;
