% Calculate inertial parameters regressor of potential energy for
% S4PRPR5
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:17
% EndTime: 2019-12-31 16:23:17
% DurationCPUTime: 0.06s
% Computational Cost: add. (73->39), mult. (91->51), div. (0->0), fcn. (85->8), ass. (0->25)
t33 = qJ(2) + pkin(7);
t29 = sin(t33);
t51 = g(3) * t29;
t50 = g(3) * qJ(1);
t34 = sin(pkin(6));
t37 = sin(qJ(4));
t49 = t34 * t37;
t39 = cos(qJ(4));
t48 = t34 * t39;
t35 = cos(pkin(6));
t47 = t35 * t37;
t46 = t35 * t39;
t40 = cos(qJ(2));
t28 = t40 * pkin(2) + pkin(1);
t36 = -qJ(3) - pkin(4);
t45 = t34 * t28 + t35 * t36;
t38 = sin(qJ(2));
t44 = t38 * pkin(2) + qJ(1);
t43 = t35 * t28 - t34 * t36;
t30 = cos(t33);
t42 = pkin(3) * t30 + pkin(5) * t29;
t41 = g(1) * t35 + g(2) * t34;
t24 = g(1) * t34 - g(2) * t35;
t23 = -g(3) * t30 + t41 * t29;
t1 = [0, 0, 0, 0, 0, 0, -t41, t24, -g(3), -t50, 0, 0, 0, 0, 0, 0, -g(3) * t38 - t41 * t40, -g(3) * t40 + t41 * t38, -t24, -g(1) * (t35 * pkin(1) + t34 * pkin(4)) - g(2) * (t34 * pkin(1) - t35 * pkin(4)) - t50, 0, 0, 0, 0, 0, 0, -t41 * t30 - t51, t23, -t24, -g(1) * t43 - g(2) * t45 - g(3) * t44, 0, 0, 0, 0, 0, 0, -g(1) * (t30 * t46 + t49) - g(2) * (t30 * t48 - t47) - t39 * t51, -g(1) * (-t30 * t47 + t48) - g(2) * (-t30 * t49 - t46) + t37 * t51, -t23, -g(1) * (t42 * t35 + t43) - g(2) * (t42 * t34 + t45) - g(3) * (t29 * pkin(3) - t30 * pkin(5) + t44);];
U_reg = t1;
