% Calculate inertial parameters regressor of potential energy for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:31
% EndTime: 2022-01-23 09:14:31
% DurationCPUTime: 0.07s
% Computational Cost: add. (116->46), mult. (86->49), div. (0->0), fcn. (73->10), ass. (0->26)
t65 = qJ(2) + pkin(5);
t72 = g(3) * t65;
t64 = cos(pkin(9));
t50 = t64 * pkin(3) + pkin(2);
t66 = -pkin(6) - qJ(3);
t61 = pkin(9) + qJ(4);
t63 = sin(pkin(9));
t71 = t63 * pkin(3) + t65;
t62 = qJ(1) + pkin(8);
t52 = sin(t62);
t54 = cos(t62);
t70 = g(1) * t54 + g(2) * t52;
t67 = sin(qJ(1));
t68 = cos(qJ(1));
t69 = -g(1) * t68 - g(2) * t67;
t60 = pkin(7) - t66;
t59 = t68 * pkin(1);
t58 = t67 * pkin(1);
t55 = qJ(5) + t61;
t53 = cos(t61);
t51 = sin(t61);
t49 = cos(t55);
t48 = sin(t55);
t46 = pkin(4) * t53 + t50;
t45 = g(1) * t52 - g(2) * t54;
t1 = [0, 0, 0, 0, 0, 0, t69, g(1) * t67 - g(2) * t68, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t70, t45, -g(3), t69 * pkin(1) - t72, 0, 0, 0, 0, 0, 0, -g(3) * t63 - t70 * t64, -g(3) * t64 + t70 * t63, -t45, -g(1) * (t54 * pkin(2) + t52 * qJ(3) + t59) - g(2) * (t52 * pkin(2) - t54 * qJ(3) + t58) - t72, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t70 * t53, -g(3) * t53 + t70 * t51, -t45, -g(1) * (t54 * t50 - t52 * t66 + t59) - g(2) * (t52 * t50 + t54 * t66 + t58) - g(3) * t71, 0, 0, 0, 0, 0, 0, -g(3) * t48 - t70 * t49, -g(3) * t49 + t70 * t48, -t45, -g(1) * (t54 * t46 + t60 * t52 + t59) - g(2) * (t52 * t46 - t54 * t60 + t58) - g(3) * (pkin(4) * t51 + t71);];
U_reg = t1;
