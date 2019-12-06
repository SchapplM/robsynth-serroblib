% Calculate inertial parameters regressor of potential energy for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:34
% EndTime: 2019-12-05 16:17:34
% DurationCPUTime: 0.08s
% Computational Cost: add. (134->45), mult. (94->57), div. (0->0), fcn. (85->10), ass. (0->28)
t65 = pkin(5) + qJ(1);
t50 = pkin(6) + t65;
t69 = g(3) * t50;
t52 = sin(pkin(9));
t68 = g(3) * t52;
t54 = cos(pkin(9));
t56 = sin(qJ(5));
t67 = t54 * t56;
t57 = cos(qJ(5));
t66 = t54 * t57;
t51 = pkin(8) + qJ(2);
t45 = sin(t51);
t53 = sin(pkin(8));
t64 = t53 * pkin(1) + pkin(2) * t45;
t46 = cos(t51);
t55 = cos(pkin(8));
t63 = t55 * pkin(1) + pkin(2) * t46;
t47 = qJ(3) + t51;
t43 = sin(t47);
t44 = cos(t47);
t62 = t44 * pkin(3) + t43 * qJ(4) + t63;
t61 = pkin(4) * t54 + pkin(7) * t52;
t60 = g(1) * t44 + g(2) * t43;
t59 = -g(1) * t55 - g(2) * t53;
t58 = t43 * pkin(3) - t44 * qJ(4) + t64;
t37 = g(1) * t43 - g(2) * t44;
t36 = -g(3) * t54 + t60 * t52;
t1 = [0, 0, 0, 0, 0, 0, t59, g(1) * t53 - g(2) * t55, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t46 - g(2) * t45, g(1) * t45 - g(2) * t46, -g(3), t59 * pkin(1) - g(3) * t65, 0, 0, 0, 0, 0, 0, -t60, t37, -g(3), -g(1) * t63 - g(2) * t64 - t69, 0, 0, 0, 0, 0, 0, -t60 * t54 - t68, t36, -t37, -g(1) * t62 - g(2) * t58 - t69, 0, 0, 0, 0, 0, 0, -g(1) * (t43 * t56 + t44 * t66) - g(2) * (t43 * t66 - t44 * t56) - t57 * t68, -g(1) * (t43 * t57 - t44 * t67) - g(2) * (-t43 * t67 - t44 * t57) + t56 * t68, -t36, -g(1) * (t61 * t44 + t62) - g(2) * (t61 * t43 + t58) - g(3) * (t52 * pkin(4) - t54 * pkin(7) + t50);];
U_reg = t1;
