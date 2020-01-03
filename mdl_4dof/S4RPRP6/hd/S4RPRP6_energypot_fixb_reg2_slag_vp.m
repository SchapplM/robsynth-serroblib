% Calculate inertial parameters regressor of potential energy for
% S4RPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP6_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:46:13
% EndTime: 2019-12-31 16:46:13
% DurationCPUTime: 0.07s
% Computational Cost: add. (40->29), mult. (63->29), div. (0->0), fcn. (53->4), ass. (0->16)
t41 = g(3) * pkin(4);
t40 = pkin(2) + pkin(4);
t33 = sin(qJ(3));
t39 = pkin(3) * t33;
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t38 = t36 * pkin(1) + t34 * qJ(2);
t30 = t34 * pkin(1);
t37 = -t36 * qJ(2) + t30;
t27 = g(1) * t34 - g(2) * t36;
t35 = cos(qJ(3));
t32 = -qJ(4) - pkin(5);
t28 = g(1) * t36 + g(2) * t34;
t26 = g(3) * t33 - t27 * t35;
t25 = -g(3) * t35 - t27 * t33;
t1 = [0, 0, 0, 0, 0, 0, -t28, t27, -g(3), -t41, 0, 0, 0, 0, 0, 0, -g(3), t28, -t27, -g(1) * t38 - g(2) * t37 - t41, 0, 0, 0, 0, 0, 0, t25, t26, -t28, -g(1) * (t36 * pkin(5) + t38) - g(2) * (t34 * pkin(5) + t37) - g(3) * t40, 0, 0, 0, 0, 0, 0, t25, t26, -t28, -g(1) * (-t36 * t32 + t34 * t39 + t38) - g(2) * (-t34 * t32 + t30 + (-qJ(2) - t39) * t36) - g(3) * (t35 * pkin(3) + t40);];
U_reg = t1;
