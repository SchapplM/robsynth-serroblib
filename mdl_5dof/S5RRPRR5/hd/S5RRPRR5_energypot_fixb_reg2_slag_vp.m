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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:02:48
% EndTime: 2022-01-20 11:02:48
% DurationCPUTime: 0.09s
% Computational Cost: add. (116->46), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t72 = pkin(6) + pkin(5);
t76 = g(3) * t72;
t68 = cos(pkin(9));
t53 = t68 * pkin(3) + pkin(2);
t69 = -pkin(7) - qJ(3);
t65 = pkin(9) + qJ(4);
t67 = sin(pkin(9));
t75 = t67 * pkin(3) + t72;
t66 = qJ(1) + qJ(2);
t58 = sin(t66);
t59 = cos(t66);
t74 = g(1) * t59 + g(2) * t58;
t70 = sin(qJ(1));
t71 = cos(qJ(1));
t73 = -g(1) * t71 - g(2) * t70;
t64 = pkin(8) - t69;
t63 = t71 * pkin(1);
t62 = t70 * pkin(1);
t57 = qJ(5) + t65;
t56 = cos(t65);
t55 = sin(t65);
t52 = cos(t57);
t51 = sin(t57);
t50 = pkin(4) * t56 + t53;
t49 = g(1) * t58 - g(2) * t59;
t1 = [0, 0, 0, 0, 0, 0, t73, g(1) * t70 - g(2) * t71, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t74, t49, -g(3), t73 * pkin(1) - t76, 0, 0, 0, 0, 0, 0, -g(3) * t67 - t74 * t68, -g(3) * t68 + t74 * t67, -t49, -g(1) * (t59 * pkin(2) + t58 * qJ(3) + t63) - g(2) * (t58 * pkin(2) - t59 * qJ(3) + t62) - t76, 0, 0, 0, 0, 0, 0, -g(3) * t55 - t74 * t56, -g(3) * t56 + t74 * t55, -t49, -g(1) * (t59 * t53 - t58 * t69 + t63) - g(2) * (t58 * t53 + t59 * t69 + t62) - g(3) * t75, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t74 * t52, -g(3) * t52 + t74 * t51, -t49, -g(1) * (t59 * t50 + t64 * t58 + t63) - g(2) * (t58 * t50 - t59 * t64 + t62) - g(3) * (pkin(4) * t55 + t75);];
U_reg = t1;
