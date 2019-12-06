% Calculate inertial parameters regressor of potential energy for
% S5PPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:12:50
% EndTime: 2019-12-05 15:12:50
% DurationCPUTime: 0.12s
% Computational Cost: add. (130->50), mult. (119->64), div. (0->0), fcn. (110->10), ass. (0->31)
t45 = pkin(9) + qJ(3);
t40 = qJ(4) + t45;
t35 = sin(t40);
t64 = g(3) * t35;
t48 = cos(pkin(9));
t37 = t48 * pkin(2) + pkin(1);
t63 = g(3) * qJ(1);
t47 = sin(pkin(8));
t51 = sin(qJ(5));
t62 = t47 * t51;
t52 = cos(qJ(5));
t61 = t47 * t52;
t49 = cos(pkin(8));
t60 = t49 * t51;
t59 = t49 * t52;
t50 = -pkin(5) - qJ(2);
t39 = cos(t45);
t31 = pkin(3) * t39 + t37;
t44 = -pkin(6) + t50;
t58 = t47 * t31 + t49 * t44;
t46 = sin(pkin(9));
t57 = t46 * pkin(2) + qJ(1);
t38 = sin(t45);
t56 = pkin(3) * t38 + t57;
t55 = t49 * t31 - t47 * t44;
t36 = cos(t40);
t54 = pkin(4) * t36 + pkin(7) * t35;
t53 = g(1) * t49 + g(2) * t47;
t32 = g(1) * t47 - g(2) * t49;
t28 = -g(3) * t36 + t53 * t35;
t1 = [0, 0, 0, 0, 0, 0, -t53, t32, -g(3), -t63, 0, 0, 0, 0, 0, 0, -g(3) * t46 - t53 * t48, -g(3) * t48 + t53 * t46, -t32, -g(1) * (t49 * pkin(1) + t47 * qJ(2)) - g(2) * (t47 * pkin(1) - t49 * qJ(2)) - t63, 0, 0, 0, 0, 0, 0, -g(3) * t38 - t53 * t39, -g(3) * t39 + t53 * t38, -t32, -g(1) * (t49 * t37 - t47 * t50) - g(2) * (t47 * t37 + t49 * t50) - g(3) * t57, 0, 0, 0, 0, 0, 0, -t53 * t36 - t64, t28, -t32, -g(1) * t55 - g(2) * t58 - g(3) * t56, 0, 0, 0, 0, 0, 0, -g(1) * (t36 * t59 + t62) - g(2) * (t36 * t61 - t60) - t52 * t64, -g(1) * (-t36 * t60 + t61) - g(2) * (-t36 * t62 - t59) + t51 * t64, -t28, -g(1) * (t54 * t49 + t55) - g(2) * (t54 * t47 + t58) - g(3) * (t35 * pkin(4) - t36 * pkin(7) + t56);];
U_reg = t1;
