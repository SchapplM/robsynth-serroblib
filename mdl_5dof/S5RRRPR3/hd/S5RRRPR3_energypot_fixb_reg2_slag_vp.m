% Calculate inertial parameters regressor of potential energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:43:09
% EndTime: 2022-01-20 11:43:09
% DurationCPUTime: 0.10s
% Computational Cost: add. (116->46), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t74 = pkin(6) + pkin(5);
t78 = g(3) * t74;
t72 = cos(qJ(3));
t56 = t72 * pkin(3) + pkin(2);
t69 = -qJ(4) - pkin(7);
t67 = qJ(3) + pkin(9);
t70 = sin(qJ(3));
t77 = t70 * pkin(3) + t74;
t68 = qJ(1) + qJ(2);
t60 = sin(t68);
t61 = cos(t68);
t76 = g(1) * t61 + g(2) * t60;
t71 = sin(qJ(1));
t73 = cos(qJ(1));
t75 = -g(1) * t73 - g(2) * t71;
t66 = pkin(8) - t69;
t65 = t73 * pkin(1);
t63 = t71 * pkin(1);
t59 = qJ(5) + t67;
t58 = cos(t67);
t57 = sin(t67);
t54 = cos(t59);
t53 = sin(t59);
t52 = pkin(4) * t58 + t56;
t51 = g(1) * t60 - g(2) * t61;
t1 = [0, 0, 0, 0, 0, 0, t75, g(1) * t71 - g(2) * t73, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t76, t51, -g(3), t75 * pkin(1) - t78, 0, 0, 0, 0, 0, 0, -g(3) * t70 - t76 * t72, -g(3) * t72 + t76 * t70, -t51, -g(1) * (t61 * pkin(2) + t60 * pkin(7) + t65) - g(2) * (t60 * pkin(2) - t61 * pkin(7) + t63) - t78, 0, 0, 0, 0, 0, 0, -g(3) * t57 - t76 * t58, -g(3) * t58 + t76 * t57, -t51, -g(1) * (t61 * t56 - t60 * t69 + t65) - g(2) * (t60 * t56 + t61 * t69 + t63) - g(3) * t77, 0, 0, 0, 0, 0, 0, -g(3) * t53 - t76 * t54, -g(3) * t54 + t76 * t53, -t51, -g(1) * (t61 * t52 + t66 * t60 + t65) - g(2) * (t60 * t52 - t61 * t66 + t63) - g(3) * (pkin(4) * t57 + t77);];
U_reg = t1;
