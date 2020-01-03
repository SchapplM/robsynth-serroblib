% Calculate inertial parameters regressor of potential energy for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:45
% EndTime: 2019-12-31 16:33:45
% DurationCPUTime: 0.06s
% Computational Cost: add. (73->39), mult. (91->51), div. (0->0), fcn. (85->8), ass. (0->25)
t37 = qJ(2) + qJ(3);
t33 = sin(t37);
t55 = g(3) * t33;
t54 = g(3) * qJ(1);
t38 = sin(pkin(7));
t40 = sin(qJ(4));
t53 = t38 * t40;
t42 = cos(qJ(4));
t52 = t38 * t42;
t39 = cos(pkin(7));
t51 = t39 * t40;
t50 = t39 * t42;
t43 = cos(qJ(2));
t32 = t43 * pkin(2) + pkin(1);
t44 = -pkin(5) - pkin(4);
t49 = t38 * t32 + t39 * t44;
t41 = sin(qJ(2));
t48 = t41 * pkin(2) + qJ(1);
t47 = t39 * t32 - t38 * t44;
t34 = cos(t37);
t46 = pkin(3) * t34 + pkin(6) * t33;
t45 = g(1) * t39 + g(2) * t38;
t28 = g(1) * t38 - g(2) * t39;
t27 = -g(3) * t34 + t45 * t33;
t1 = [0, 0, 0, 0, 0, 0, -t45, t28, -g(3), -t54, 0, 0, 0, 0, 0, 0, -g(3) * t41 - t45 * t43, -g(3) * t43 + t45 * t41, -t28, -g(1) * (t39 * pkin(1) + t38 * pkin(4)) - g(2) * (t38 * pkin(1) - t39 * pkin(4)) - t54, 0, 0, 0, 0, 0, 0, -t45 * t34 - t55, t27, -t28, -g(1) * t47 - g(2) * t49 - g(3) * t48, 0, 0, 0, 0, 0, 0, -g(1) * (t34 * t50 + t53) - g(2) * (t34 * t52 - t51) - t42 * t55, -g(1) * (-t34 * t51 + t52) - g(2) * (-t34 * t53 - t50) + t40 * t55, -t27, -g(1) * (t46 * t39 + t47) - g(2) * (t46 * t38 + t49) - g(3) * (t33 * pkin(3) - t34 * pkin(6) + t48);];
U_reg = t1;
