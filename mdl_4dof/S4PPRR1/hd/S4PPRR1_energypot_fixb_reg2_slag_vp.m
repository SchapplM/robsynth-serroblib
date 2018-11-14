% Calculate inertial parameters regressor of potential energy for
% S4PPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:40
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPRR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:40:18
% EndTime: 2018-11-14 13:40:18
% DurationCPUTime: 0.05s
% Computational Cost: add. (50->33), mult. (64->41), div. (0->0), fcn. (62->6), ass. (0->21)
t34 = g(3) * qJ(1);
t26 = sin(pkin(6));
t28 = sin(qJ(3));
t33 = t26 * t28;
t32 = -pkin(4) + qJ(1);
t27 = cos(pkin(6));
t31 = t27 * pkin(1) + t26 * qJ(2);
t23 = t26 * pkin(1);
t30 = -t27 * qJ(2) + t23;
t29 = cos(qJ(3));
t25 = qJ(3) + qJ(4);
t22 = cos(t25);
t21 = sin(t25);
t19 = t29 * pkin(3) + pkin(2);
t18 = -g(1) * t27 - g(2) * t26;
t17 = g(1) * t26 - g(2) * t27;
t16 = t26 * t29 - t27 * t28;
t15 = -t27 * t29 - t33;
t14 = -t27 * t21 + t26 * t22;
t13 = -t26 * t21 - t27 * t22;
t1 = [0, 0, 0, 0, 0, 0, t18, t17, -g(3), -t34, 0, 0, 0, 0, 0, 0, t18, -g(3), -t17, -g(1) * t31 - g(2) * t30 - t34, 0, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t16, -g(1) * t16 - g(2) * t15, g(3), -g(1) * (t27 * pkin(2) + t31) - g(2) * (t26 * pkin(2) + t30) - g(3) * t32, 0, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t14, -g(1) * t14 - g(2) * t13, g(3), -g(1) * (pkin(3) * t33 + t27 * t19 + t31) - g(2) * (t26 * t19 + t23 + (-pkin(3) * t28 - qJ(2)) * t27) - g(3) * (-pkin(5) + t32);];
U_reg  = t1;
