% Calculate inertial parameters regressor of potential energy for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->33), mult. (59->36), div. (0->0), fcn. (49->8), ass. (0->20)
t52 = pkin(5) + pkin(4);
t55 = g(3) * t52;
t46 = qJ(1) + qJ(2);
t40 = sin(t46);
t42 = cos(t46);
t54 = g(1) * t42 + g(2) * t40;
t48 = sin(qJ(1));
t50 = cos(qJ(1));
t53 = -g(1) * t50 - g(2) * t48;
t51 = -pkin(7) - pkin(6);
t49 = cos(qJ(3));
t47 = sin(qJ(3));
t45 = qJ(3) + qJ(4);
t44 = t50 * pkin(1);
t43 = t48 * pkin(1);
t41 = cos(t45);
t39 = sin(t45);
t38 = t49 * pkin(3) + pkin(2);
t36 = g(1) * t40 - g(2) * t42;
t1 = [0, 0, 0, 0, 0, 0, t53, g(1) * t48 - g(2) * t50, -g(3), -g(3) * pkin(4), 0, 0, 0, 0, 0, 0, -t54, t36, -g(3), t53 * pkin(1) - t55, 0, 0, 0, 0, 0, 0, -g(3) * t47 - t54 * t49, -g(3) * t49 + t54 * t47, -t36, -g(1) * (t42 * pkin(2) + t40 * pkin(6) + t44) - g(2) * (t40 * pkin(2) - t42 * pkin(6) + t43) - t55, 0, 0, 0, 0, 0, 0, -g(3) * t39 - t54 * t41, -g(3) * t41 + t54 * t39, -t36, -g(1) * (t42 * t38 - t40 * t51 + t44) - g(2) * (t40 * t38 + t42 * t51 + t43) - g(3) * (t47 * pkin(3) + t52);];
U_reg = t1;
