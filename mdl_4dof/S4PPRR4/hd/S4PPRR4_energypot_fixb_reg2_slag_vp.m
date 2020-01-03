% Calculate inertial parameters regressor of potential energy for
% S4PPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:39
% EndTime: 2019-12-31 16:18:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (73->39), mult. (91->51), div. (0->0), fcn. (85->8), ass. (0->25)
t32 = pkin(7) + qJ(3);
t28 = sin(t32);
t50 = g(3) * t28;
t49 = g(3) * qJ(1);
t34 = sin(pkin(6));
t38 = sin(qJ(4));
t48 = t34 * t38;
t39 = cos(qJ(4));
t47 = t34 * t39;
t36 = cos(pkin(6));
t46 = t36 * t38;
t45 = t36 * t39;
t35 = cos(pkin(7));
t27 = t35 * pkin(2) + pkin(1);
t37 = -pkin(4) - qJ(2);
t44 = t34 * t27 + t36 * t37;
t33 = sin(pkin(7));
t43 = t33 * pkin(2) + qJ(1);
t42 = t36 * t27 - t34 * t37;
t29 = cos(t32);
t41 = pkin(3) * t29 + pkin(5) * t28;
t40 = g(1) * t36 + g(2) * t34;
t23 = g(1) * t34 - g(2) * t36;
t22 = -g(3) * t29 + t40 * t28;
t1 = [0, 0, 0, 0, 0, 0, -t40, t23, -g(3), -t49, 0, 0, 0, 0, 0, 0, -g(3) * t33 - t40 * t35, -g(3) * t35 + t40 * t33, -t23, -g(1) * (t36 * pkin(1) + t34 * qJ(2)) - g(2) * (t34 * pkin(1) - t36 * qJ(2)) - t49, 0, 0, 0, 0, 0, 0, -t40 * t29 - t50, t22, -t23, -g(1) * t42 - g(2) * t44 - g(3) * t43, 0, 0, 0, 0, 0, 0, -g(1) * (t29 * t45 + t48) - g(2) * (t29 * t47 - t46) - t39 * t50, -g(1) * (-t29 * t46 + t47) - g(2) * (-t29 * t48 - t45) + t38 * t50, -t22, -g(1) * (t41 * t36 + t42) - g(2) * (t41 * t34 + t44) - g(3) * (t28 * pkin(3) - t29 * pkin(5) + t43);];
U_reg = t1;
