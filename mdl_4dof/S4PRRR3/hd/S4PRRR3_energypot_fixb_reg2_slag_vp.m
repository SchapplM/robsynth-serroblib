% Calculate inertial parameters regressor of potential energy for
% S4PRRR3
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:38
% EndTime: 2019-12-31 16:31:38
% DurationCPUTime: 0.04s
% Computational Cost: add. (70->29), mult. (48->32), div. (0->0), fcn. (38->8), ass. (0->18)
t40 = pkin(4) + qJ(1);
t41 = g(3) * (pkin(5) + t40);
t31 = pkin(7) + qJ(2);
t25 = sin(t31);
t32 = sin(pkin(7));
t39 = t32 * pkin(1) + pkin(2) * t25;
t26 = cos(t31);
t33 = cos(pkin(7));
t38 = t33 * pkin(1) + pkin(2) * t26;
t27 = qJ(3) + t31;
t23 = sin(t27);
t24 = cos(t27);
t37 = g(1) * t24 + g(2) * t23;
t36 = -g(1) * t33 - g(2) * t32;
t35 = cos(qJ(4));
t34 = sin(qJ(4));
t20 = g(1) * t23 - g(2) * t24;
t1 = [0, 0, 0, 0, 0, 0, t36, g(1) * t32 - g(2) * t33, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t26 - g(2) * t25, g(1) * t25 - g(2) * t26, -g(3), t36 * pkin(1) - g(3) * t40, 0, 0, 0, 0, 0, 0, -t37, t20, -g(3), -g(1) * t38 - g(2) * t39 - t41, 0, 0, 0, 0, 0, 0, -g(3) * t34 - t37 * t35, -g(3) * t35 + t37 * t34, -t20, -g(1) * (t24 * pkin(3) + t23 * pkin(6) + t38) - g(2) * (t23 * pkin(3) - t24 * pkin(6) + t39) - t41;];
U_reg = t1;
