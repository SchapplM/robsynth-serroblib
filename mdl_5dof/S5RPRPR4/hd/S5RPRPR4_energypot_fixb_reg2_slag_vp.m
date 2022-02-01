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
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:23:05
% EndTime: 2022-01-23 09:23:05
% DurationCPUTime: 0.09s
% Computational Cost: add. (116->46), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t67 = qJ(2) + pkin(5);
t75 = g(3) * t67;
t70 = cos(qJ(3));
t53 = t70 * pkin(3) + pkin(2);
t66 = -qJ(4) - pkin(6);
t64 = qJ(3) + pkin(9);
t68 = sin(qJ(3));
t74 = t68 * pkin(3) + t67;
t65 = qJ(1) + pkin(8);
t55 = sin(t65);
t57 = cos(t65);
t73 = g(1) * t57 + g(2) * t55;
t69 = sin(qJ(1));
t71 = cos(qJ(1));
t72 = -g(1) * t71 - g(2) * t69;
t63 = pkin(7) - t66;
t62 = t71 * pkin(1);
t60 = t69 * pkin(1);
t58 = qJ(5) + t64;
t56 = cos(t64);
t54 = sin(t64);
t52 = cos(t58);
t51 = sin(t58);
t49 = pkin(4) * t56 + t53;
t48 = g(1) * t55 - g(2) * t57;
t1 = [0, 0, 0, 0, 0, 0, t72, g(1) * t69 - g(2) * t71, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t73, t48, -g(3), t72 * pkin(1) - t75, 0, 0, 0, 0, 0, 0, -g(3) * t68 - t73 * t70, -g(3) * t70 + t73 * t68, -t48, -g(1) * (t57 * pkin(2) + t55 * pkin(6) + t62) - g(2) * (t55 * pkin(2) - t57 * pkin(6) + t60) - t75, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t73 * t56, -g(3) * t56 + t73 * t54, -t48, -g(1) * (t57 * t53 - t55 * t66 + t62) - g(2) * (t55 * t53 + t57 * t66 + t60) - g(3) * t74, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t73 * t52, -g(3) * t52 + t73 * t51, -t48, -g(1) * (t57 * t49 + t63 * t55 + t62) - g(2) * (t55 * t49 - t57 * t63 + t60) - g(3) * (pkin(4) * t54 + t74);];
U_reg = t1;
