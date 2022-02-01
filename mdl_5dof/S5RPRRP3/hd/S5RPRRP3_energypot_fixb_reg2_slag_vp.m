% Calculate inertial parameters regressor of potential energy for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:30:15
% EndTime: 2022-01-23 09:30:15
% DurationCPUTime: 0.06s
% Computational Cost: add. (110->43), mult. (86->45), div. (0->0), fcn. (73->8), ass. (0->25)
t66 = -pkin(7) - pkin(6);
t61 = qJ(2) + pkin(5);
t70 = g(3) * t61;
t64 = cos(qJ(3));
t49 = t64 * pkin(3) + pkin(2);
t62 = sin(qJ(3));
t69 = t62 * pkin(3) + t61;
t59 = qJ(1) + pkin(8);
t50 = sin(t59);
t51 = cos(t59);
t68 = g(1) * t51 + g(2) * t50;
t63 = sin(qJ(1));
t65 = cos(qJ(1));
t67 = -g(1) * t65 - g(2) * t63;
t60 = qJ(3) + qJ(4);
t58 = -qJ(5) + t66;
t57 = t65 * pkin(1);
t55 = t63 * pkin(1);
t53 = cos(t60);
t52 = sin(t60);
t47 = pkin(4) * t53 + t49;
t46 = g(1) * t50 - g(2) * t51;
t45 = -g(3) * t52 - t68 * t53;
t44 = -g(3) * t53 + t68 * t52;
t1 = [0, 0, 0, 0, 0, 0, t67, g(1) * t63 - g(2) * t65, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -t68, t46, -g(3), t67 * pkin(1) - t70, 0, 0, 0, 0, 0, 0, -g(3) * t62 - t68 * t64, -g(3) * t64 + t68 * t62, -t46, -g(1) * (t51 * pkin(2) + t50 * pkin(6) + t57) - g(2) * (t50 * pkin(2) - t51 * pkin(6) + t55) - t70, 0, 0, 0, 0, 0, 0, t45, t44, -t46, -g(1) * (t51 * t49 - t50 * t66 + t57) - g(2) * (t50 * t49 + t51 * t66 + t55) - g(3) * t69, 0, 0, 0, 0, 0, 0, t45, t44, -t46, -g(1) * (t51 * t47 - t50 * t58 + t57) - g(2) * (t50 * t47 + t51 * t58 + t55) - g(3) * (pkin(4) * t52 + t69);];
U_reg = t1;
