% Calculate inertial parameters regressor of potential energy for
% S4RPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR7_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR7_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:40
% EndTime: 2019-12-31 16:41:40
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->32), mult. (63->33), div. (0->0), fcn. (53->6), ass. (0->17)
t41 = g(3) * pkin(4);
t40 = pkin(2) + pkin(4);
t32 = sin(pkin(6));
t39 = pkin(3) * t32;
t35 = sin(qJ(1));
t36 = cos(qJ(1));
t38 = t36 * pkin(1) + t35 * qJ(2);
t29 = t35 * pkin(1);
t37 = -t36 * qJ(2) + t29;
t24 = g(1) * t35 - g(2) * t36;
t34 = -pkin(5) - qJ(3);
t33 = cos(pkin(6));
t31 = pkin(6) + qJ(4);
t27 = cos(t31);
t26 = sin(t31);
t25 = g(1) * t36 + g(2) * t35;
t1 = [0, 0, 0, 0, 0, 0, -t25, t24, -g(3), -t41, 0, 0, 0, 0, 0, 0, -g(3), t25, -t24, -g(1) * t38 - g(2) * t37 - t41, 0, 0, 0, 0, 0, 0, -g(3) * t33 - t24 * t32, g(3) * t32 - t24 * t33, -t25, -g(1) * (t36 * qJ(3) + t38) - g(2) * (t35 * qJ(3) + t37) - g(3) * t40, 0, 0, 0, 0, 0, 0, -g(3) * t27 - t24 * t26, g(3) * t26 - t24 * t27, -t25, -g(1) * (-t36 * t34 + t35 * t39 + t38) - g(2) * (-t35 * t34 + t29 + (-qJ(2) - t39) * t36) - g(3) * (t33 * pkin(3) + t40);];
U_reg = t1;
