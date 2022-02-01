% Calculate inertial parameters regressor of potential energy for
% S5RRRPR2
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
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:48
% EndTime: 2022-01-20 11:30:48
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->38), mult. (63->41), div. (0->0), fcn. (50->10), ass. (0->24)
t62 = pkin(6) + pkin(5);
t58 = pkin(7) + t62;
t61 = g(3) * (qJ(4) + t58);
t49 = qJ(1) + qJ(2);
t43 = sin(t49);
t51 = sin(qJ(1));
t60 = t51 * pkin(1) + pkin(2) * t43;
t44 = cos(t49);
t53 = cos(qJ(1));
t59 = t53 * pkin(1) + pkin(2) * t44;
t46 = qJ(3) + t49;
t41 = sin(t46);
t57 = pkin(3) * t41 + t60;
t42 = cos(t46);
t56 = pkin(3) * t42 + t59;
t40 = pkin(9) + t46;
t34 = sin(t40);
t35 = cos(t40);
t55 = g(1) * t35 + g(2) * t34;
t54 = -g(1) * t53 - g(2) * t51;
t52 = cos(qJ(5));
t50 = sin(qJ(5));
t33 = g(1) * t34 - g(2) * t35;
t1 = [0, 0, 0, 0, 0, 0, t54, g(1) * t51 - g(2) * t53, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t44 - g(2) * t43, g(1) * t43 - g(2) * t44, -g(3), t54 * pkin(1) - g(3) * t62, 0, 0, 0, 0, 0, 0, -g(1) * t42 - g(2) * t41, g(1) * t41 - g(2) * t42, -g(3), -g(1) * t59 - g(2) * t60 - g(3) * t58, 0, 0, 0, 0, 0, 0, -t55, t33, -g(3), -g(1) * t56 - g(2) * t57 - t61, 0, 0, 0, 0, 0, 0, -g(3) * t50 - t55 * t52, -g(3) * t52 + t55 * t50, -t33, -g(1) * (t35 * pkin(4) + t34 * pkin(8) + t56) - g(2) * (t34 * pkin(4) - t35 * pkin(8) + t57) - t61;];
U_reg = t1;
