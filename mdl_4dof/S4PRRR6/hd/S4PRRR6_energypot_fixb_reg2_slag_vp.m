% Calculate inertial parameters regressor of potential energy for
% S4PRRR6
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:35:06
% EndTime: 2019-12-31 16:35:06
% DurationCPUTime: 0.11s
% Computational Cost: add. (71->51), mult. (115->71), div. (0->0), fcn. (113->8), ass. (0->27)
t45 = sin(qJ(2));
t60 = g(3) * t45;
t59 = g(3) * qJ(1);
t42 = sin(pkin(7));
t44 = sin(qJ(3));
t58 = t42 * t44;
t47 = cos(qJ(2));
t57 = t42 * t47;
t43 = cos(pkin(7));
t56 = t43 * t47;
t55 = t44 * t47;
t48 = -pkin(6) - pkin(5);
t54 = t45 * t48;
t46 = cos(qJ(3));
t53 = t46 * t47;
t52 = t43 * pkin(1) + t42 * pkin(4);
t38 = t42 * pkin(1);
t51 = -t43 * pkin(4) + t38;
t50 = pkin(2) * t47 + pkin(5) * t45;
t49 = g(1) * t43 + g(2) * t42;
t41 = qJ(3) + qJ(4);
t36 = cos(t41);
t35 = sin(t41);
t34 = t46 * pkin(3) + pkin(2);
t33 = g(1) * t42 - g(2) * t43;
t32 = -g(3) * t47 + t49 * t45;
t1 = [0, 0, 0, 0, 0, 0, -t49, t33, -g(3), -t59, 0, 0, 0, 0, 0, 0, -t49 * t47 - t60, t32, -t33, -g(1) * t52 - g(2) * t51 - t59, 0, 0, 0, 0, 0, 0, -g(1) * (t43 * t53 + t58) - g(2) * (t42 * t53 - t43 * t44) - t46 * t60, -g(1) * (t42 * t46 - t43 * t55) - g(2) * (-t42 * t55 - t43 * t46) + t44 * t60, -t32, -g(1) * (t50 * t43 + t52) - g(2) * (t50 * t42 + t51) - g(3) * (t45 * pkin(2) - t47 * pkin(5) + qJ(1)), 0, 0, 0, 0, 0, 0, -g(1) * (t42 * t35 + t36 * t56) - g(2) * (-t43 * t35 + t36 * t57) - t36 * t60, -g(1) * (-t35 * t56 + t42 * t36) - g(2) * (-t35 * t57 - t43 * t36) + t35 * t60, -t32, -g(1) * (pkin(3) * t58 + t52) - g(2) * (t34 * t57 - t42 * t54 + t38) - g(3) * (t45 * t34 + t47 * t48 + qJ(1)) + (-g(1) * (t34 * t47 - t54) - g(2) * (-pkin(3) * t44 - pkin(4))) * t43;];
U_reg = t1;
