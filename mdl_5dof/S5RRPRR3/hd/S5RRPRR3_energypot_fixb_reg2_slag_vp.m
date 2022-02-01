% Calculate inertial parameters regressor of potential energy for
% S5RRPRR3
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
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:15
% EndTime: 2022-01-20 10:34:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->38), mult. (63->41), div. (0->0), fcn. (50->10), ass. (0->24)
t64 = pkin(6) + pkin(5);
t60 = qJ(3) + t64;
t63 = g(3) * (pkin(7) + t60);
t51 = qJ(1) + qJ(2);
t46 = sin(t51);
t53 = sin(qJ(1));
t62 = t53 * pkin(1) + pkin(2) * t46;
t47 = cos(t51);
t55 = cos(qJ(1));
t61 = t55 * pkin(1) + pkin(2) * t47;
t45 = pkin(9) + t51;
t40 = sin(t45);
t59 = pkin(3) * t40 + t62;
t41 = cos(t45);
t58 = pkin(3) * t41 + t61;
t44 = qJ(4) + t45;
t38 = sin(t44);
t39 = cos(t44);
t57 = g(1) * t39 + g(2) * t38;
t56 = -g(1) * t55 - g(2) * t53;
t54 = cos(qJ(5));
t52 = sin(qJ(5));
t35 = g(1) * t38 - g(2) * t39;
t1 = [0, 0, 0, 0, 0, 0, t56, g(1) * t53 - g(2) * t55, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t47 - g(2) * t46, g(1) * t46 - g(2) * t47, -g(3), t56 * pkin(1) - g(3) * t64, 0, 0, 0, 0, 0, 0, -g(1) * t41 - g(2) * t40, g(1) * t40 - g(2) * t41, -g(3), -g(1) * t61 - g(2) * t62 - g(3) * t60, 0, 0, 0, 0, 0, 0, -t57, t35, -g(3), -g(1) * t58 - g(2) * t59 - t63, 0, 0, 0, 0, 0, 0, -g(3) * t52 - t57 * t54, -g(3) * t54 + t57 * t52, -t35, -g(1) * (t39 * pkin(4) + t38 * pkin(8) + t58) - g(2) * (t38 * pkin(4) - t39 * pkin(8) + t59) - t63;];
U_reg = t1;
