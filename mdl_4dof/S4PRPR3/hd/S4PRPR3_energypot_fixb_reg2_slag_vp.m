% Calculate inertial parameters regressor of potential energy for
% S4PRPR3
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:54
% EndTime: 2019-12-31 16:20:54
% DurationCPUTime: 0.05s
% Computational Cost: add. (69->33), mult. (59->36), div. (0->0), fcn. (49->8), ass. (0->20)
t42 = pkin(4) + qJ(1);
t45 = g(3) * t42;
t36 = pkin(6) + qJ(2);
t30 = sin(t36);
t32 = cos(t36);
t44 = g(1) * t32 + g(2) * t30;
t38 = sin(pkin(6));
t40 = cos(pkin(6));
t43 = -g(1) * t40 - g(2) * t38;
t41 = -pkin(5) - qJ(3);
t39 = cos(pkin(7));
t37 = sin(pkin(7));
t35 = pkin(7) + qJ(4);
t34 = t40 * pkin(1);
t33 = t38 * pkin(1);
t31 = cos(t35);
t29 = sin(t35);
t28 = t39 * pkin(3) + pkin(2);
t26 = g(1) * t30 - g(2) * t32;
t1 = [0, 0, 0, 0, 0, 0, t43, g(1) * t38 - g(2) * t40, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t44, t26, -g(3), t43 * pkin(1) - t45, 0, 0, 0, 0, 0, 0, -g(3) * t37 - t44 * t39, -g(3) * t39 + t44 * t37, -t26, -g(1) * (t32 * pkin(2) + t30 * qJ(3) + t34) - g(2) * (t30 * pkin(2) - t32 * qJ(3) + t33) - t45, 0, 0, 0, 0, 0, 0, -g(3) * t29 - t44 * t31, -g(3) * t31 + t44 * t29, -t26, -g(1) * (t32 * t28 - t30 * t41 + t34) - g(2) * (t30 * t28 + t32 * t41 + t33) - g(3) * (t37 * pkin(3) + t42);];
U_reg = t1;
