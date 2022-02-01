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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:20:50
% EndTime: 2022-01-23 09:20:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (134->45), mult. (94->57), div. (0->0), fcn. (85->10), ass. (0->28)
t73 = qJ(2) + pkin(5);
t58 = pkin(6) + t73;
t77 = g(3) * t58;
t60 = sin(pkin(9));
t76 = g(3) * t60;
t61 = cos(pkin(9));
t62 = sin(qJ(5));
t75 = t61 * t62;
t64 = cos(qJ(5));
t74 = t61 * t64;
t59 = qJ(1) + pkin(8);
t53 = sin(t59);
t63 = sin(qJ(1));
t72 = t63 * pkin(1) + pkin(2) * t53;
t54 = cos(t59);
t65 = cos(qJ(1));
t71 = t65 * pkin(1) + pkin(2) * t54;
t55 = qJ(3) + t59;
t51 = sin(t55);
t52 = cos(t55);
t70 = t52 * pkin(3) + t51 * qJ(4) + t71;
t69 = pkin(4) * t61 + pkin(7) * t60;
t68 = g(1) * t52 + g(2) * t51;
t67 = -g(1) * t65 - g(2) * t63;
t66 = t51 * pkin(3) - t52 * qJ(4) + t72;
t45 = g(1) * t51 - g(2) * t52;
t44 = -g(3) * t61 + t60 * t68;
t1 = [0, 0, 0, 0, 0, 0, t67, g(1) * t63 - g(2) * t65, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t54 - g(2) * t53, g(1) * t53 - g(2) * t54, -g(3), pkin(1) * t67 - g(3) * t73, 0, 0, 0, 0, 0, 0, -t68, t45, -g(3), -g(1) * t71 - g(2) * t72 - t77, 0, 0, 0, 0, 0, 0, -t61 * t68 - t76, t44, -t45, -g(1) * t70 - g(2) * t66 - t77, 0, 0, 0, 0, 0, 0, -g(1) * (t51 * t62 + t52 * t74) - g(2) * (t51 * t74 - t52 * t62) - t64 * t76, -g(1) * (t51 * t64 - t52 * t75) - g(2) * (-t51 * t75 - t52 * t64) + t62 * t76, -t44, -g(1) * (t52 * t69 + t70) - g(2) * (t51 * t69 + t66) - g(3) * (pkin(4) * t60 - pkin(7) * t61 + t58);];
U_reg = t1;
