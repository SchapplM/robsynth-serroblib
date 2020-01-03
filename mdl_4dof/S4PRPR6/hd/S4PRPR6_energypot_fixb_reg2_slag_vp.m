% Calculate inertial parameters regressor of potential energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:40
% EndTime: 2019-12-31 16:24:40
% DurationCPUTime: 0.09s
% Computational Cost: add. (71->51), mult. (115->69), div. (0->0), fcn. (113->8), ass. (0->25)
t44 = sin(qJ(2));
t55 = g(3) * t44;
t54 = g(3) * qJ(1);
t39 = sin(pkin(7));
t40 = sin(pkin(6));
t53 = t40 * t39;
t45 = cos(qJ(2));
t52 = t40 * t45;
t42 = cos(pkin(6));
t51 = t42 * t45;
t43 = -pkin(5) - qJ(3);
t50 = t43 * t44;
t49 = t42 * pkin(1) + t40 * pkin(4);
t35 = t40 * pkin(1);
t48 = -t42 * pkin(4) + t35;
t47 = g(1) * t42 + g(2) * t40;
t46 = pkin(2) * t45 + qJ(3) * t44;
t41 = cos(pkin(7));
t38 = pkin(7) + qJ(4);
t33 = cos(t38);
t32 = sin(t38);
t31 = t41 * pkin(3) + pkin(2);
t30 = g(1) * t40 - g(2) * t42;
t29 = -g(3) * t45 + t47 * t44;
t1 = [0, 0, 0, 0, 0, 0, -t47, t30, -g(3), -t54, 0, 0, 0, 0, 0, 0, -t47 * t45 - t55, t29, -t30, -g(1) * t49 - g(2) * t48 - t54, 0, 0, 0, 0, 0, 0, -g(1) * (t41 * t51 + t53) - g(2) * (-t42 * t39 + t41 * t52) - t41 * t55, -g(1) * (-t39 * t51 + t40 * t41) - g(2) * (-t39 * t52 - t42 * t41) + t39 * t55, -t29, -g(1) * (t46 * t42 + t49) - g(2) * (t46 * t40 + t48) - g(3) * (t44 * pkin(2) - t45 * qJ(3) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * (t40 * t32 + t33 * t51) - g(2) * (-t42 * t32 + t33 * t52) - t33 * t55, -g(1) * (-t32 * t51 + t40 * t33) - g(2) * (-t32 * t52 - t42 * t33) + t32 * t55, -t29, -g(1) * (pkin(3) * t53 + t49) - g(2) * (t31 * t52 - t40 * t50 + t35) - g(3) * (t44 * t31 + t45 * t43 + qJ(1)) + (-g(1) * (t31 * t45 - t50) - g(2) * (-pkin(3) * t39 - pkin(4))) * t42;];
U_reg = t1;
