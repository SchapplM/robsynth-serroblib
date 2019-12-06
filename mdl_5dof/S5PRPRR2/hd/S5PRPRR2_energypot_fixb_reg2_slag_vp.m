% Calculate inertial parameters regressor of potential energy for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:12
% EndTime: 2019-12-05 15:45:12
% DurationCPUTime: 0.12s
% Computational Cost: add. (130->50), mult. (119->64), div. (0->0), fcn. (110->10), ass. (0->31)
t49 = qJ(2) + pkin(9);
t44 = qJ(4) + t49;
t39 = sin(t44);
t68 = g(3) * t39;
t56 = cos(qJ(2));
t41 = t56 * pkin(2) + pkin(1);
t67 = g(3) * qJ(1);
t50 = sin(pkin(8));
t53 = sin(qJ(5));
t66 = t50 * t53;
t55 = cos(qJ(5));
t65 = t50 * t55;
t51 = cos(pkin(8));
t64 = t51 * t53;
t63 = t51 * t55;
t52 = -qJ(3) - pkin(5);
t43 = cos(t49);
t35 = pkin(3) * t43 + t41;
t48 = -pkin(6) + t52;
t62 = t50 * t35 + t51 * t48;
t54 = sin(qJ(2));
t61 = t54 * pkin(2) + qJ(1);
t42 = sin(t49);
t60 = pkin(3) * t42 + t61;
t59 = t51 * t35 - t50 * t48;
t40 = cos(t44);
t58 = pkin(4) * t40 + pkin(7) * t39;
t57 = g(1) * t51 + g(2) * t50;
t36 = g(1) * t50 - g(2) * t51;
t32 = -g(3) * t40 + t57 * t39;
t1 = [0, 0, 0, 0, 0, 0, -t57, t36, -g(3), -t67, 0, 0, 0, 0, 0, 0, -g(3) * t54 - t57 * t56, -g(3) * t56 + t57 * t54, -t36, -g(1) * (t51 * pkin(1) + t50 * pkin(5)) - g(2) * (t50 * pkin(1) - t51 * pkin(5)) - t67, 0, 0, 0, 0, 0, 0, -g(3) * t42 - t57 * t43, -g(3) * t43 + t57 * t42, -t36, -g(1) * (t51 * t41 - t50 * t52) - g(2) * (t50 * t41 + t51 * t52) - g(3) * t61, 0, 0, 0, 0, 0, 0, -t57 * t40 - t68, t32, -t36, -g(1) * t59 - g(2) * t62 - g(3) * t60, 0, 0, 0, 0, 0, 0, -g(1) * (t40 * t63 + t66) - g(2) * (t40 * t65 - t64) - t55 * t68, -g(1) * (-t40 * t64 + t65) - g(2) * (-t40 * t66 - t63) + t53 * t68, -t32, -g(1) * (t58 * t51 + t59) - g(2) * (t58 * t50 + t62) - g(3) * (t39 * pkin(4) - t40 * pkin(7) + t60);];
U_reg = t1;
