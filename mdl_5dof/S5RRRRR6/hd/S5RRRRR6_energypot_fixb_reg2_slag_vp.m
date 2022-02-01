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
% Datum: 2022-01-20 12:09
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:08:21
% EndTime: 2022-01-20 12:08:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (116->46), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t74 = pkin(6) + pkin(5);
t73 = -pkin(8) - pkin(7);
t78 = g(3) * t74;
t71 = cos(qJ(3));
t56 = t71 * pkin(3) + pkin(2);
t67 = qJ(3) + qJ(4);
t69 = sin(qJ(3));
t77 = t69 * pkin(3) + t74;
t68 = qJ(1) + qJ(2);
t58 = sin(t68);
t60 = cos(t68);
t76 = g(1) * t60 + g(2) * t58;
t70 = sin(qJ(1));
t72 = cos(qJ(1));
t75 = -g(1) * t72 - g(2) * t70;
t66 = pkin(9) - t73;
t65 = t72 * pkin(1);
t63 = t70 * pkin(1);
t61 = qJ(5) + t67;
t59 = cos(t67);
t57 = sin(t67);
t55 = cos(t61);
t54 = sin(t61);
t52 = pkin(4) * t59 + t56;
t51 = g(1) * t58 - g(2) * t60;
t1 = [0, 0, 0, 0, 0, 0, t75, g(1) * t70 - g(2) * t72, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t76, t51, -g(3), pkin(1) * t75 - t78, 0, 0, 0, 0, 0, 0, -g(3) * t69 - t71 * t76, -g(3) * t71 + t69 * t76, -t51, -g(1) * (pkin(2) * t60 + pkin(7) * t58 + t65) - g(2) * (pkin(2) * t58 - pkin(7) * t60 + t63) - t78, 0, 0, 0, 0, 0, 0, -g(3) * t57 - t59 * t76, -g(3) * t59 + t57 * t76, -t51, -g(1) * (t56 * t60 - t58 * t73 + t65) - g(2) * (t56 * t58 + t60 * t73 + t63) - g(3) * t77, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t55 * t76, -g(3) * t55 + t54 * t76, -t51, -g(1) * (t52 * t60 + t58 * t66 + t65) - g(2) * (t52 * t58 - t60 * t66 + t63) - g(3) * (pkin(4) * t57 + t77);];
U_reg = t1;
