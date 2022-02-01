% Calculate inertial parameters regressor of potential energy for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:55
% EndTime: 2022-01-20 10:19:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (112->39), mult. (74->41), div. (0->0), fcn. (61->8), ass. (0->23)
t67 = pkin(6) + pkin(5);
t55 = qJ(3) + t67;
t66 = g(3) * t55;
t56 = qJ(1) + qJ(2);
t51 = sin(t56);
t59 = sin(qJ(1));
t65 = t59 * pkin(1) + pkin(2) * t51;
t52 = cos(t56);
t61 = cos(qJ(1));
t64 = t61 * pkin(1) + pkin(2) * t52;
t50 = pkin(8) + t56;
t45 = sin(t50);
t46 = cos(t50);
t63 = g(1) * t46 + g(2) * t45;
t62 = -g(1) * t61 - g(2) * t59;
t60 = cos(qJ(4));
t58 = sin(qJ(4));
t57 = -qJ(5) - pkin(7);
t49 = t60 * pkin(4) + pkin(3);
t43 = g(1) * t45 - g(2) * t46;
t42 = -g(3) * t58 - t63 * t60;
t41 = -g(3) * t60 + t63 * t58;
t1 = [0, 0, 0, 0, 0, 0, t62, g(1) * t59 - g(2) * t61, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t52 - g(2) * t51, g(1) * t51 - g(2) * t52, -g(3), t62 * pkin(1) - g(3) * t67, 0, 0, 0, 0, 0, 0, -t63, t43, -g(3), -g(1) * t64 - g(2) * t65 - t66, 0, 0, 0, 0, 0, 0, t42, t41, -t43, -g(1) * (t46 * pkin(3) + t45 * pkin(7) + t64) - g(2) * (t45 * pkin(3) - t46 * pkin(7) + t65) - t66, 0, 0, 0, 0, 0, 0, t42, t41, -t43, -g(1) * (-t45 * t57 + t46 * t49 + t64) - g(2) * (t45 * t49 + t46 * t57 + t65) - g(3) * (t58 * pkin(4) + t55);];
U_reg = t1;
