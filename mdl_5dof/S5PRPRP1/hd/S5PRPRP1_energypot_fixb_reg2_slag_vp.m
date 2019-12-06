% Calculate inertial parameters regressor of potential energy for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energypot_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:38
% EndTime: 2019-12-05 15:28:38
% DurationCPUTime: 0.07s
% Computational Cost: add. (118->42), mult. (93->45), div. (0->0), fcn. (80->8), ass. (0->26)
t59 = pkin(5) + qJ(1);
t66 = g(3) * t59;
t56 = cos(pkin(8));
t44 = t56 * pkin(3) + pkin(2);
t53 = pkin(7) + qJ(2);
t46 = sin(t53);
t48 = cos(t53);
t55 = sin(pkin(7));
t50 = t55 * pkin(1);
t58 = -pkin(6) - qJ(3);
t65 = t46 * t44 + t48 * t58 + t50;
t54 = sin(pkin(8));
t64 = t54 * pkin(3) + t59;
t57 = cos(pkin(7));
t51 = t57 * pkin(1);
t63 = t48 * t44 - t46 * t58 + t51;
t62 = g(1) * t48 + g(2) * t46;
t61 = -g(1) * t57 - g(2) * t55;
t52 = pkin(8) + qJ(4);
t45 = sin(t52);
t47 = cos(t52);
t60 = pkin(4) * t47 + qJ(5) * t45;
t39 = g(1) * t46 - g(2) * t48;
t38 = -g(3) * t45 - t62 * t47;
t37 = -g(3) * t47 + t62 * t45;
t1 = [0, 0, 0, 0, 0, 0, t61, g(1) * t55 - g(2) * t57, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t62, t39, -g(3), t61 * pkin(1) - t66, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t62 * t56, -g(3) * t56 + t62 * t54, -t39, -g(1) * (t48 * pkin(2) + t46 * qJ(3) + t51) - g(2) * (t46 * pkin(2) - t48 * qJ(3) + t50) - t66, 0, 0, 0, 0, 0, 0, t38, t37, -t39, -g(1) * t63 - g(2) * t65 - g(3) * t64, 0, 0, 0, 0, 0, 0, t38, -t39, -t37, -g(1) * (t60 * t48 + t63) - g(2) * (t60 * t46 + t65) - g(3) * (t45 * pkin(4) - t47 * qJ(5) + t64);];
U_reg = t1;
