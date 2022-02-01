% Calculate inertial parameters regressor of potential energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:51:38
% EndTime: 2022-01-20 09:51:38
% DurationCPUTime: 0.07s
% Computational Cost: add. (118->42), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t69 = pkin(6) + pkin(5);
t56 = qJ(3) + t69;
t68 = g(3) * t56;
t58 = qJ(1) + qJ(2);
t52 = sin(t58);
t62 = sin(qJ(1));
t67 = t62 * pkin(1) + pkin(2) * t52;
t53 = cos(t58);
t63 = cos(qJ(1));
t66 = t63 * pkin(1) + pkin(2) * t53;
t51 = pkin(8) + t58;
t44 = sin(t51);
t45 = cos(t51);
t65 = g(1) * t45 + g(2) * t44;
t64 = -g(1) * t63 - g(2) * t62;
t61 = -pkin(7) - qJ(4);
t60 = cos(pkin(9));
t59 = sin(pkin(9));
t57 = pkin(9) + qJ(5);
t50 = cos(t57);
t49 = sin(t57);
t46 = t60 * pkin(4) + pkin(3);
t42 = g(1) * t44 - g(2) * t45;
t1 = [0, 0, 0, 0, 0, 0, t64, g(1) * t62 - g(2) * t63, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t53 - g(2) * t52, g(1) * t52 - g(2) * t53, -g(3), t64 * pkin(1) - g(3) * t69, 0, 0, 0, 0, 0, 0, -t65, t42, -g(3), -g(1) * t66 - g(2) * t67 - t68, 0, 0, 0, 0, 0, 0, -g(3) * t59 - t65 * t60, -g(3) * t60 + t65 * t59, -t42, -g(1) * (t45 * pkin(3) + t44 * qJ(4) + t66) - g(2) * (t44 * pkin(3) - t45 * qJ(4) + t67) - t68, 0, 0, 0, 0, 0, 0, -g(3) * t49 - t65 * t50, -g(3) * t50 + t65 * t49, -t42, -g(1) * (-t44 * t61 + t45 * t46 + t66) - g(2) * (t44 * t46 + t45 * t61 + t67) - g(3) * (t59 * pkin(4) + t56);];
U_reg = t1;
