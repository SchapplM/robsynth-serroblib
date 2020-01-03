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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:36:20
% EndTime: 2020-01-03 11:36:20
% DurationCPUTime: 0.14s
% Computational Cost: add. (134->46), mult. (94->59), div. (0->0), fcn. (85->10), ass. (0->28)
t56 = sin(pkin(9));
t57 = cos(pkin(9));
t75 = pkin(4) * t57 + pkin(7) * t56;
t68 = qJ(2) + pkin(5);
t54 = pkin(6) + t68;
t72 = g(1) * t54;
t71 = g(1) * t56;
t58 = sin(qJ(5));
t70 = t57 * t58;
t60 = cos(qJ(5));
t69 = t57 * t60;
t55 = qJ(1) + pkin(8);
t50 = sin(t55);
t59 = sin(qJ(1));
t67 = t59 * pkin(1) + pkin(2) * t50;
t52 = qJ(3) + t55;
t48 = sin(t52);
t66 = t48 * pkin(3) + t67;
t51 = cos(t55);
t61 = cos(qJ(1));
t65 = -t61 * pkin(1) - pkin(2) * t51;
t49 = cos(t52);
t64 = -g(2) * t48 + g(3) * t49;
t63 = -g(2) * t59 + g(3) * t61;
t62 = -t48 * qJ(4) + t65;
t45 = g(2) * t49 + g(3) * t48;
t44 = g(1) * t57 + t64 * t56;
t1 = [0, 0, 0, 0, 0, 0, t63, -g(2) * t61 - g(3) * t59, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, -g(2) * t50 + g(3) * t51, -g(2) * t51 - g(3) * t50, -g(1), pkin(1) * t63 - g(1) * t68, 0, 0, 0, 0, 0, 0, t64, -t45, -g(1), -g(2) * t67 - g(3) * t65 - t72, 0, 0, 0, 0, 0, 0, t57 * t64 - t71, -t44, t45, -t72 - g(2) * (-t49 * qJ(4) + t66) - g(3) * (-pkin(3) * t49 + t62), 0, 0, 0, 0, 0, 0, -t60 * t71 - g(2) * (t48 * t69 - t49 * t58) - g(3) * (-t48 * t58 - t49 * t69), t58 * t71 - g(2) * (-t48 * t70 - t49 * t60) - g(3) * (-t48 * t60 + t49 * t70), t44, -g(1) * (pkin(4) * t56 - pkin(7) * t57 + t54) - g(2) * (t48 * t75 + t66) - g(3) * t62 + (g(2) * qJ(4) - g(3) * (-pkin(3) - t75)) * t49;];
U_reg = t1;
