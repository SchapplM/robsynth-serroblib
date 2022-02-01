% Calculate inertial parameters regressor of potential energy for
% S5RPRRP2
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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:09
% EndTime: 2022-01-23 09:28:09
% DurationCPUTime: 0.06s
% Computational Cost: add. (112->39), mult. (74->41), div. (0->0), fcn. (61->8), ass. (0->23)
t67 = qJ(2) + pkin(5);
t56 = pkin(6) + t67;
t68 = g(3) * t56;
t57 = qJ(1) + pkin(8);
t51 = sin(t57);
t60 = sin(qJ(1));
t66 = t60 * pkin(1) + pkin(2) * t51;
t52 = cos(t57);
t62 = cos(qJ(1));
t65 = t62 * pkin(1) + pkin(2) * t52;
t53 = qJ(3) + t57;
t48 = sin(t53);
t49 = cos(t53);
t64 = g(1) * t49 + g(2) * t48;
t63 = -g(1) * t62 - g(2) * t60;
t61 = cos(qJ(4));
t59 = sin(qJ(4));
t58 = -qJ(5) - pkin(7);
t50 = t61 * pkin(4) + pkin(3);
t44 = g(1) * t48 - g(2) * t49;
t43 = -g(3) * t59 - t64 * t61;
t42 = -g(3) * t61 + t64 * t59;
t1 = [0, 0, 0, 0, 0, 0, t63, g(1) * t60 - g(2) * t62, -g(3), -g(3) * pkin(5), 0, 0, 0, 0, 0, 0, -g(1) * t52 - g(2) * t51, g(1) * t51 - g(2) * t52, -g(3), t63 * pkin(1) - g(3) * t67, 0, 0, 0, 0, 0, 0, -t64, t44, -g(3), -g(1) * t65 - g(2) * t66 - t68, 0, 0, 0, 0, 0, 0, t43, t42, -t44, -g(1) * (t49 * pkin(3) + t48 * pkin(7) + t65) - g(2) * (t48 * pkin(3) - t49 * pkin(7) + t66) - t68, 0, 0, 0, 0, 0, 0, t43, t42, -t44, -g(1) * (-t48 * t58 + t49 * t50 + t65) - g(2) * (t48 * t50 + t49 * t58 + t66) - g(3) * (t59 * pkin(4) + t56);];
U_reg = t1;
