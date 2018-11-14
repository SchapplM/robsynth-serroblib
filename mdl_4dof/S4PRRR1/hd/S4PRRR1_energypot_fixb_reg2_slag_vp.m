% Calculate inertial parameters regressor of potential energy for
% S4PRRR1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energypot_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:44:29
% EndTime: 2018-11-14 13:44:29
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->28), mult. (38->31), div. (0->0), fcn. (28->8), ass. (0->17)
t35 = pkin(4) + qJ(1);
t28 = pkin(7) + qJ(2);
t23 = sin(t28);
t29 = sin(pkin(7));
t34 = t29 * pkin(1) + pkin(2) * t23;
t24 = cos(t28);
t30 = cos(pkin(7));
t33 = t30 * pkin(1) + pkin(2) * t24;
t32 = pkin(5) + t35;
t25 = qJ(3) + t28;
t31 = -g(1) * t30 - g(2) * t29;
t22 = qJ(4) + t25;
t21 = cos(t25);
t20 = sin(t25);
t17 = cos(t22);
t16 = sin(t22);
t1 = [0, 0, 0, 0, 0, 0, t31, g(1) * t29 - g(2) * t30, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t23, g(1) * t23 - g(2) * t24, -g(3), t31 * pkin(1) - g(3) * t35, 0, 0, 0, 0, 0, 0, -g(1) * t21 - g(2) * t20, g(1) * t20 - g(2) * t21, -g(3), -g(1) * t33 - g(2) * t34 - g(3) * t32, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t16, g(1) * t16 - g(2) * t17, -g(3), -g(1) * (pkin(3) * t21 + t33) - g(2) * (pkin(3) * t20 + t34) - g(3) * (pkin(6) + t32);];
U_reg  = t1;
