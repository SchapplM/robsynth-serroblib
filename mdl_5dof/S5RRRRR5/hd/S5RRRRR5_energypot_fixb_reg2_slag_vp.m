% Calculate inertial parameters regressor of potential energy for
% S5RRRRR5
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:02:05
% EndTime: 2022-01-20 12:02:05
% DurationCPUTime: 0.12s
% Computational Cost: add. (118->42), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t72 = pkin(6) + pkin(5);
t59 = pkin(7) + t72;
t71 = g(3) * t59;
t61 = qJ(1) + qJ(2);
t53 = sin(t61);
t63 = sin(qJ(1));
t70 = t63 * pkin(1) + pkin(2) * t53;
t55 = cos(t61);
t65 = cos(qJ(1));
t69 = t65 * pkin(1) + pkin(2) * t55;
t56 = qJ(3) + t61;
t49 = sin(t56);
t50 = cos(t56);
t68 = g(1) * t50 + g(2) * t49;
t67 = -g(1) * t65 - g(2) * t63;
t66 = -pkin(9) - pkin(8);
t64 = cos(qJ(4));
t62 = sin(qJ(4));
t60 = qJ(4) + qJ(5);
t54 = cos(t60);
t52 = sin(t60);
t51 = t64 * pkin(4) + pkin(3);
t45 = g(1) * t49 - g(2) * t50;
t1 = [0, 0, 0, 0, 0, 0, t67, g(1) * t63 - g(2) * t65, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t55 - g(2) * t53, g(1) * t53 - g(2) * t55, -g(3), t67 * pkin(1) - g(3) * t72, 0, 0, 0, 0, 0, 0, -t68, t45, -g(3), -g(1) * t69 - g(2) * t70 - t71, 0, 0, 0, 0, 0, 0, -g(3) * t62 - t68 * t64, -g(3) * t64 + t68 * t62, -t45, -g(1) * (t50 * pkin(3) + t49 * pkin(8) + t69) - g(2) * (t49 * pkin(3) - t50 * pkin(8) + t70) - t71, 0, 0, 0, 0, 0, 0, -g(3) * t52 - t68 * t54, -g(3) * t54 + t68 * t52, -t45, -g(1) * (-t49 * t66 + t50 * t51 + t69) - g(2) * (t49 * t51 + t50 * t66 + t70) - g(3) * (t62 * pkin(4) + t59);];
U_reg = t1;
