% Calculate inertial parameters regressor of potential energy for
% S5PRRRR4
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:57
% EndTime: 2019-12-05 17:07:57
% DurationCPUTime: 0.06s
% Computational Cost: add. (118->42), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t62 = pkin(5) + qJ(1);
t50 = pkin(6) + t62;
t63 = g(3) * t50;
t51 = pkin(9) + qJ(2);
t43 = sin(t51);
t53 = sin(pkin(9));
t61 = t53 * pkin(1) + pkin(2) * t43;
t44 = cos(t51);
t54 = cos(pkin(9));
t60 = t54 * pkin(1) + pkin(2) * t44;
t45 = qJ(3) + t51;
t40 = sin(t45);
t41 = cos(t45);
t59 = g(1) * t41 + g(2) * t40;
t58 = -g(1) * t54 - g(2) * t53;
t57 = -pkin(8) - pkin(7);
t56 = cos(qJ(4));
t55 = sin(qJ(4));
t52 = qJ(4) + qJ(5);
t47 = cos(t52);
t46 = sin(t52);
t42 = t56 * pkin(4) + pkin(3);
t36 = g(1) * t40 - g(2) * t41;
t1 = [0, 0, 0, 0, 0, 0, t58, g(1) * t53 - g(2) * t54, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t44 - g(2) * t43, g(1) * t43 - g(2) * t44, -g(3), t58 * pkin(1) - g(3) * t62, 0, 0, 0, 0, 0, 0, -t59, t36, -g(3), -g(1) * t60 - g(2) * t61 - t63, 0, 0, 0, 0, 0, 0, -g(3) * t55 - t59 * t56, -g(3) * t56 + t59 * t55, -t36, -g(1) * (t41 * pkin(3) + t40 * pkin(7) + t60) - g(2) * (t40 * pkin(3) - t41 * pkin(7) + t61) - t63, 0, 0, 0, 0, 0, 0, -g(3) * t46 - t59 * t47, -g(3) * t47 + t59 * t46, -t36, -g(1) * (-t40 * t57 + t41 * t42 + t60) - g(2) * (t40 * t42 + t41 * t57 + t61) - g(3) * (t55 * pkin(4) + t50);];
U_reg = t1;
