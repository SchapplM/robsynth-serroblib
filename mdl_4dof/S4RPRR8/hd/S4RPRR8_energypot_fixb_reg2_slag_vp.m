% Calculate inertial parameters regressor of potential energy for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRR8_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:15
% EndTime: 2019-12-31 16:55:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->32), mult. (63->33), div. (0->0), fcn. (53->6), ass. (0->17)
t43 = g(3) * pkin(4);
t42 = pkin(2) + pkin(4);
t34 = sin(qJ(3));
t41 = pkin(3) * t34;
t35 = sin(qJ(1));
t37 = cos(qJ(1));
t40 = t37 * pkin(1) + t35 * qJ(2);
t31 = t35 * pkin(1);
t39 = -t37 * qJ(2) + t31;
t26 = g(1) * t35 - g(2) * t37;
t38 = -pkin(6) - pkin(5);
t36 = cos(qJ(3));
t33 = qJ(3) + qJ(4);
t29 = cos(t33);
t28 = sin(t33);
t27 = g(1) * t37 + g(2) * t35;
t1 = [0, 0, 0, 0, 0, 0, -t27, t26, -g(3), -t43, 0, 0, 0, 0, 0, 0, -g(3), t27, -t26, -g(1) * t40 - g(2) * t39 - t43, 0, 0, 0, 0, 0, 0, -g(3) * t36 - t26 * t34, g(3) * t34 - t26 * t36, -t27, -g(1) * (t37 * pkin(5) + t40) - g(2) * (t35 * pkin(5) + t39) - g(3) * t42, 0, 0, 0, 0, 0, 0, -g(3) * t29 - t26 * t28, g(3) * t28 - t26 * t29, -t27, -g(1) * (t35 * t41 - t37 * t38 + t40) - g(2) * (-t35 * t38 + t31 + (-qJ(2) - t41) * t37) - g(3) * (t36 * pkin(3) + t42);];
U_reg = t1;
