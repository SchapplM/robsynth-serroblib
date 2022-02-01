% Calculate inertial parameters regressor of potential energy for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:57
% EndTime: 2022-01-20 09:48:57
% DurationCPUTime: 0.07s
% Computational Cost: add. (118->42), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t70 = qJ(2) + pkin(5);
t58 = pkin(6) + t70;
t71 = g(3) * t58;
t59 = qJ(1) + pkin(9);
t51 = sin(t59);
t62 = sin(qJ(1));
t69 = t62 * pkin(1) + pkin(2) * t51;
t52 = cos(t59);
t64 = cos(qJ(1));
t68 = t64 * pkin(1) + pkin(2) * t52;
t53 = qJ(3) + t59;
t48 = sin(t53);
t49 = cos(t53);
t67 = g(1) * t49 + g(2) * t48;
t66 = -g(1) * t64 - g(2) * t62;
t65 = -pkin(8) - pkin(7);
t63 = cos(qJ(4));
t61 = sin(qJ(4));
t60 = qJ(4) + qJ(5);
t55 = cos(t60);
t54 = sin(t60);
t50 = t63 * pkin(4) + pkin(3);
t44 = g(1) * t48 - g(2) * t49;
t1 = [0, 0, 0, 0, 0, 0, t66, g(1) * t62 - g(2) * t64, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t52 - g(2) * t51, g(1) * t51 - g(2) * t52, -g(3), t66 * pkin(1) - g(3) * t70, 0, 0, 0, 0, 0, 0, -t67, t44, -g(3), -g(1) * t68 - g(2) * t69 - t71, 0, 0, 0, 0, 0, 0, -g(3) * t61 - t67 * t63, -g(3) * t63 + t67 * t61, -t44, -g(1) * (t49 * pkin(3) + t48 * pkin(7) + t68) - g(2) * (t48 * pkin(3) - t49 * pkin(7) + t69) - t71, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t67 * t55, -g(3) * t55 + t67 * t54, -t44, -g(1) * (-t48 * t65 + t49 * t50 + t68) - g(2) * (t48 * t50 + t49 * t65 + t69) - g(3) * (t61 * pkin(4) + t58);];
U_reg = t1;
