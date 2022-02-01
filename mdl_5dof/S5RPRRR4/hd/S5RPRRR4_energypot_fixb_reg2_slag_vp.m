% Calculate inertial parameters regressor of potential energy for
% S5RPRRR4
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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:41
% EndTime: 2022-01-23 09:34:41
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->38), mult. (63->41), div. (0->0), fcn. (50->10), ass. (0->24)
t62 = qJ(2) + pkin(5);
t59 = pkin(6) + t62;
t63 = g(3) * (pkin(7) + t59);
t50 = qJ(1) + pkin(9);
t44 = sin(t50);
t52 = sin(qJ(1));
t61 = t52 * pkin(1) + pkin(2) * t44;
t45 = cos(t50);
t54 = cos(qJ(1));
t60 = t54 * pkin(1) + pkin(2) * t45;
t46 = qJ(3) + t50;
t41 = sin(t46);
t58 = pkin(3) * t41 + t61;
t42 = cos(t46);
t57 = pkin(3) * t42 + t60;
t43 = qJ(4) + t46;
t37 = sin(t43);
t38 = cos(t43);
t56 = g(1) * t38 + g(2) * t37;
t55 = -g(1) * t54 - g(2) * t52;
t53 = cos(qJ(5));
t51 = sin(qJ(5));
t34 = g(1) * t37 - g(2) * t38;
t1 = [0, 0, 0, 0, 0, 0, t55, g(1) * t52 - g(2) * t54, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t45 - g(2) * t44, g(1) * t44 - g(2) * t45, -g(3), t55 * pkin(1) - g(3) * t62, 0, 0, 0, 0, 0, 0, -g(1) * t42 - g(2) * t41, g(1) * t41 - g(2) * t42, -g(3), -g(1) * t60 - g(2) * t61 - g(3) * t59, 0, 0, 0, 0, 0, 0, -t56, t34, -g(3), -g(1) * t57 - g(2) * t58 - t63, 0, 0, 0, 0, 0, 0, -g(3) * t51 - t56 * t53, -g(3) * t53 + t56 * t51, -t34, -g(1) * (t38 * pkin(4) + t37 * pkin(8) + t57) - g(2) * (t37 * pkin(4) - t38 * pkin(8) + t58) - t63;];
U_reg = t1;
