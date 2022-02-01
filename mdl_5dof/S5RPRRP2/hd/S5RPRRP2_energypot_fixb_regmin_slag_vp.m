% Calculate minimal parameter regressor of potential energy for
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
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:12
% EndTime: 2022-01-23 09:28:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->24), mult. (49->27), div. (0->0), fcn. (43->8), ass. (0->17)
t53 = qJ(2) + pkin(5);
t45 = qJ(1) + pkin(8);
t44 = qJ(3) + t45;
t41 = sin(t44);
t42 = cos(t44);
t52 = g(1) * t42 + g(2) * t41;
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t51 = -g(1) * t50 - g(2) * t48;
t49 = cos(qJ(4));
t47 = sin(qJ(4));
t46 = -qJ(5) - pkin(7);
t43 = t49 * pkin(4) + pkin(3);
t40 = g(1) * t41 - g(2) * t42;
t39 = -g(3) * t47 - t52 * t49;
t38 = -g(3) * t49 + t52 * t47;
t1 = [0, t51, g(1) * t48 - g(2) * t50, t51 * pkin(1) - g(3) * t53, 0, -t52, t40, 0, 0, 0, 0, 0, t39, t38, t39, t38, -t40, -g(1) * (t42 * t43 - t41 * t46 + pkin(2) * cos(t45) + t50 * pkin(1)) - g(2) * (t41 * t43 + t42 * t46 + pkin(2) * sin(t45) + t48 * pkin(1)) - g(3) * (t47 * pkin(4) + pkin(6) + t53);];
U_reg = t1;
