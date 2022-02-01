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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:12:40
% EndTime: 2022-01-23 09:12:41
% DurationCPUTime: 0.11s
% Computational Cost: add. (128->49), mult. (130->62), div. (0->0), fcn. (125->8), ass. (0->29)
t66 = cos(qJ(4));
t53 = t66 * pkin(4) + pkin(3);
t60 = sin(pkin(8));
t61 = cos(pkin(8));
t62 = -qJ(5) - pkin(6);
t81 = t53 * t61 - t60 * t62;
t80 = g(3) * t60;
t63 = qJ(2) + pkin(5);
t79 = g(3) * t63;
t59 = qJ(1) + pkin(7);
t54 = sin(t59);
t64 = sin(qJ(4));
t77 = t54 * t64;
t75 = t61 * t64;
t74 = t61 * t66;
t65 = sin(qJ(1));
t73 = t65 * pkin(1) + t54 * pkin(2);
t55 = cos(t59);
t67 = cos(qJ(1));
t72 = t67 * pkin(1) + t55 * pkin(2) + t54 * qJ(3);
t71 = -t55 * qJ(3) + t73;
t70 = pkin(3) * t61 + pkin(6) * t60;
t69 = g(1) * t55 + g(2) * t54;
t68 = -g(1) * t67 - g(2) * t65;
t49 = g(1) * t54 - g(2) * t55;
t48 = -g(3) * t61 + t69 * t60;
t47 = -g(1) * (t55 * t74 + t77) - g(2) * (t54 * t74 - t55 * t64) - t66 * t80;
t46 = -g(1) * (t54 * t66 - t55 * t75) - g(2) * (-t54 * t75 - t55 * t66) + t64 * t80;
t1 = [0, 0, 0, 0, 0, 0, t68, g(1) * t65 - g(2) * t67, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t69, t49, -g(3), t68 * pkin(1) - t79, 0, 0, 0, 0, 0, 0, -t69 * t61 - t80, t48, -t49, -g(1) * t72 - g(2) * t71 - t79, 0, 0, 0, 0, 0, 0, t47, t46, -t48, -g(1) * (t70 * t55 + t72) - g(2) * (t70 * t54 + t71) - g(3) * (t60 * pkin(3) - t61 * pkin(6) + t63), 0, 0, 0, 0, 0, 0, t47, t46, -t48, -g(1) * (pkin(4) * t77 + t72) - g(2) * (t81 * t54 + t73) - g(3) * (t60 * t53 + t61 * t62 + t63) + (-g(1) * t81 - g(2) * (-pkin(4) * t64 - qJ(3))) * t55;];
U_reg = t1;
