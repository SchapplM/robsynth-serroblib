% Calculate inertial parameters regressor of potential energy for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:17
% EndTime: 2019-12-05 17:06:18
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->38), mult. (63->41), div. (0->0), fcn. (50->10), ass. (0->24)
t56 = pkin(5) + qJ(1);
t53 = pkin(6) + t56;
t57 = g(3) * (pkin(7) + t53);
t44 = pkin(9) + qJ(2);
t38 = sin(t44);
t45 = sin(pkin(9));
t55 = t45 * pkin(1) + pkin(2) * t38;
t39 = cos(t44);
t46 = cos(pkin(9));
t54 = t46 * pkin(1) + pkin(2) * t39;
t40 = qJ(3) + t44;
t35 = sin(t40);
t52 = pkin(3) * t35 + t55;
t36 = cos(t40);
t51 = pkin(3) * t36 + t54;
t37 = qJ(4) + t40;
t31 = sin(t37);
t32 = cos(t37);
t50 = g(1) * t32 + g(2) * t31;
t49 = -g(1) * t46 - g(2) * t45;
t48 = cos(qJ(5));
t47 = sin(qJ(5));
t28 = g(1) * t31 - g(2) * t32;
t1 = [0, 0, 0, 0, 0, 0, t49, g(1) * t45 - g(2) * t46, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t39 - g(2) * t38, g(1) * t38 - g(2) * t39, -g(3), t49 * pkin(1) - g(3) * t56, 0, 0, 0, 0, 0, 0, -g(1) * t36 - g(2) * t35, g(1) * t35 - g(2) * t36, -g(3), -g(1) * t54 - g(2) * t55 - g(3) * t53, 0, 0, 0, 0, 0, 0, -t50, t28, -g(3), -g(1) * t51 - g(2) * t52 - t57, 0, 0, 0, 0, 0, 0, -g(3) * t47 - t50 * t48, -g(3) * t48 + t50 * t47, -t28, -g(1) * (t32 * pkin(4) + t31 * pkin(8) + t51) - g(2) * (t31 * pkin(4) - t32 * pkin(8) + t52) - t57;];
U_reg = t1;
