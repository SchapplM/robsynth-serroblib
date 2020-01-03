% Calculate inertial parameters regressor of potential energy for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:04
% EndTime: 2019-12-31 16:54:04
% DurationCPUTime: 0.07s
% Computational Cost: add. (73->39), mult. (91->51), div. (0->0), fcn. (85->8), ass. (0->25)
t65 = g(3) * pkin(4);
t47 = pkin(7) + qJ(3);
t43 = sin(t47);
t64 = g(3) * t43;
t48 = sin(pkin(7));
t63 = t48 * pkin(2) + pkin(4);
t51 = sin(qJ(4));
t52 = sin(qJ(1));
t62 = t52 * t51;
t53 = cos(qJ(4));
t61 = t52 * t53;
t54 = cos(qJ(1));
t60 = t54 * t51;
t59 = t54 * t53;
t49 = cos(pkin(7));
t41 = t49 * pkin(2) + pkin(1);
t50 = -pkin(5) - qJ(2);
t58 = t52 * t41 + t54 * t50;
t57 = t54 * t41 - t52 * t50;
t44 = cos(t47);
t56 = pkin(3) * t44 + pkin(6) * t43;
t55 = g(1) * t54 + g(2) * t52;
t40 = g(1) * t52 - g(2) * t54;
t37 = -g(3) * t44 + t55 * t43;
t1 = [0, 0, 0, 0, 0, 0, -t55, t40, -g(3), -t65, 0, 0, 0, 0, 0, 0, -g(3) * t48 - t55 * t49, -g(3) * t49 + t55 * t48, -t40, -g(1) * (t54 * pkin(1) + t52 * qJ(2)) - g(2) * (t52 * pkin(1) - t54 * qJ(2)) - t65, 0, 0, 0, 0, 0, 0, -t55 * t44 - t64, t37, -t40, -g(1) * t57 - g(2) * t58 - g(3) * t63, 0, 0, 0, 0, 0, 0, -g(1) * (t44 * t59 + t62) - g(2) * (t44 * t61 - t60) - t53 * t64, -g(1) * (-t44 * t60 + t61) - g(2) * (-t44 * t62 - t59) + t51 * t64, -t37, -g(1) * (t56 * t54 + t57) - g(2) * (t56 * t52 + t58) - g(3) * (t43 * pkin(3) - t44 * pkin(6) + t63);];
U_reg = t1;
