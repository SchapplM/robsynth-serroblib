% Calculate inertial parameters regressor of potential energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:26
% EndTime: 2019-12-05 18:18:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (118->41), mult. (74->45), div. (0->0), fcn. (61->10), ass. (0->24)
t66 = pkin(6) + pkin(5);
t53 = qJ(3) + t66;
t65 = g(1) * t53;
t55 = qJ(1) + qJ(2);
t51 = cos(t55);
t60 = cos(qJ(1));
t64 = t60 * pkin(1) + pkin(2) * t51;
t50 = sin(t55);
t59 = sin(qJ(1));
t63 = -t59 * pkin(1) - pkin(2) * t50;
t49 = pkin(8) + t55;
t43 = sin(t49);
t44 = cos(t49);
t62 = g(2) * t43 - g(3) * t44;
t61 = g(2) * t59 - g(3) * t60;
t58 = -pkin(7) - qJ(4);
t57 = cos(pkin(9));
t56 = sin(pkin(9));
t54 = pkin(9) + qJ(5);
t48 = cos(t54);
t47 = sin(t54);
t45 = t57 * pkin(4) + pkin(3);
t42 = g(2) * t44 + g(3) * t43;
t1 = [0, 0, 0, 0, 0, 0, t61, g(2) * t60 + g(3) * t59, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, g(2) * t50 - g(3) * t51, g(2) * t51 + g(3) * t50, -g(1), t61 * pkin(1) - g(1) * t66, 0, 0, 0, 0, 0, 0, t62, t42, -g(1), -g(2) * t63 - g(3) * t64 - t65, 0, 0, 0, 0, 0, 0, -g(1) * t56 + t62 * t57, -g(1) * t57 - t62 * t56, -t42, -t65 - g(2) * (-t43 * pkin(3) + t44 * qJ(4) + t63) - g(3) * (t44 * pkin(3) + t43 * qJ(4) + t64), 0, 0, 0, 0, 0, 0, -g(1) * t47 + t62 * t48, -g(1) * t48 - t62 * t47, -t42, -g(1) * (t56 * pkin(4) + t53) - g(2) * (-t43 * t45 - t44 * t58 + t63) - g(3) * (-t43 * t58 + t44 * t45 + t64);];
U_reg = t1;
