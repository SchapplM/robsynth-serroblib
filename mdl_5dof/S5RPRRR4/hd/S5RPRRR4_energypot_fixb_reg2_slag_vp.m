% Calculate inertial parameters regressor of potential energy for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% U_reg [1x(5*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_energypot_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_energypot_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:14:34
% EndTime: 2019-12-05 18:14:34
% DurationCPUTime: 0.06s
% Computational Cost: add. (117->37), mult. (63->41), div. (0->0), fcn. (50->10), ass. (0->24)
t59 = qJ(2) + pkin(5);
t57 = pkin(6) + t59;
t60 = g(1) * (pkin(7) + t57);
t47 = qJ(1) + pkin(9);
t43 = cos(t47);
t51 = cos(qJ(1));
t58 = t51 * pkin(1) + pkin(2) * t43;
t44 = qJ(3) + t47;
t40 = cos(t44);
t56 = pkin(3) * t40 + t58;
t42 = sin(t47);
t49 = sin(qJ(1));
t55 = -t49 * pkin(1) - pkin(2) * t42;
t41 = qJ(4) + t44;
t36 = sin(t41);
t37 = cos(t41);
t54 = g(2) * t36 - g(3) * t37;
t53 = g(2) * t49 - g(3) * t51;
t39 = sin(t44);
t52 = -pkin(3) * t39 + t55;
t50 = cos(qJ(5));
t48 = sin(qJ(5));
t34 = g(2) * t37 + g(3) * t36;
t1 = [0, 0, 0, 0, 0, 0, t53, g(2) * t51 + g(3) * t49, -g(1), -g(1) * pkin(5), 0, 0, 0, 0, 0, 0, g(2) * t42 - g(3) * t43, g(2) * t43 + g(3) * t42, -g(1), t53 * pkin(1) - g(1) * t59, 0, 0, 0, 0, 0, 0, g(2) * t39 - g(3) * t40, g(2) * t40 + g(3) * t39, -g(1), -g(1) * t57 - g(2) * t55 - g(3) * t58, 0, 0, 0, 0, 0, 0, t54, t34, -g(1), -g(2) * t52 - g(3) * t56 - t60, 0, 0, 0, 0, 0, 0, -g(1) * t48 + t54 * t50, -g(1) * t50 - t54 * t48, -t34, -t60 - g(2) * (-t36 * pkin(4) + t37 * pkin(8) + t52) - g(3) * (t37 * pkin(4) + t36 * pkin(8) + t56);];
U_reg = t1;
