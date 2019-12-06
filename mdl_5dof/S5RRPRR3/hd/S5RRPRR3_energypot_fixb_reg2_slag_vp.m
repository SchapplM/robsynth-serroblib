% Calculate inertial parameters regressor of potential energy for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:30:32
% EndTime: 2019-12-05 18:30:32
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->37), mult. (63->41), div. (0->0), fcn. (50->10), ass. (0->24)
t61 = pkin(6) + pkin(5);
t58 = qJ(3) + t61;
t60 = g(1) * (pkin(7) + t58);
t48 = qJ(1) + qJ(2);
t45 = cos(t48);
t52 = cos(qJ(1));
t59 = t52 * pkin(1) + pkin(2) * t45;
t43 = pkin(9) + t48;
t40 = cos(t43);
t57 = pkin(3) * t40 + t59;
t44 = sin(t48);
t50 = sin(qJ(1));
t56 = -t50 * pkin(1) - pkin(2) * t44;
t42 = qJ(4) + t43;
t37 = sin(t42);
t38 = cos(t42);
t55 = g(2) * t37 - g(3) * t38;
t54 = g(2) * t50 - g(3) * t52;
t39 = sin(t43);
t53 = -pkin(3) * t39 + t56;
t51 = cos(qJ(5));
t49 = sin(qJ(5));
t35 = g(2) * t38 + g(3) * t37;
t1 = [0, 0, 0, 0, 0, 0, t54, g(2) * t52 + g(3) * t50, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, g(2) * t44 - g(3) * t45, g(2) * t45 + g(3) * t44, -g(1), t54 * pkin(1) - g(1) * t61, 0, 0, 0, 0, 0, 0, g(2) * t39 - g(3) * t40, g(2) * t40 + g(3) * t39, -g(1), -g(1) * t58 - g(2) * t56 - g(3) * t59, 0, 0, 0, 0, 0, 0, t55, t35, -g(1), -g(2) * t53 - g(3) * t57 - t60, 0, 0, 0, 0, 0, 0, -g(1) * t49 + t55 * t51, -g(1) * t51 - t55 * t49, -t35, -t60 - g(2) * (-t37 * pkin(4) + t38 * pkin(8) + t53) - g(3) * (t38 * pkin(4) + t37 * pkin(8) + t57);];
U_reg = t1;
