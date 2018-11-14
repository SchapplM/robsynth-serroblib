% Calculate inertial parameters regressor of potential energy for
% S4RPPR2
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:48
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:47:31
% EndTime: 2018-11-14 13:47:31
% DurationCPUTime: 0.05s
% Computational Cost: add. (50->33), mult. (64->41), div. (0->0), fcn. (62->6), ass. (0->21)
t39 = g(3) * pkin(4);
t31 = sin(pkin(6));
t33 = sin(qJ(1));
t38 = t33 * t31;
t37 = -qJ(3) + pkin(4);
t34 = cos(qJ(1));
t36 = t34 * pkin(1) + t33 * qJ(2);
t28 = t33 * pkin(1);
t35 = -t34 * qJ(2) + t28;
t32 = cos(pkin(6));
t30 = pkin(6) + qJ(4);
t26 = cos(t30);
t25 = sin(t30);
t24 = t32 * pkin(3) + pkin(2);
t23 = -g(1) * t34 - g(2) * t33;
t22 = g(1) * t33 - g(2) * t34;
t21 = -t34 * t31 + t33 * t32;
t20 = -t34 * t32 - t38;
t19 = -t34 * t25 + t33 * t26;
t18 = -t33 * t25 - t34 * t26;
t1 = [0, 0, 0, 0, 0, 0, t23, t22, -g(3), -t39, 0, 0, 0, 0, 0, 0, t23, -g(3), -t22, -g(1) * t36 - g(2) * t35 - t39, 0, 0, 0, 0, 0, 0, g(1) * t20 - g(2) * t21, -g(1) * t21 - g(2) * t20, g(3), -g(1) * (t34 * pkin(2) + t36) - g(2) * (t33 * pkin(2) + t35) - g(3) * t37, 0, 0, 0, 0, 0, 0, g(1) * t18 - g(2) * t19, -g(1) * t19 - g(2) * t18, g(3), -g(1) * (pkin(3) * t38 + t34 * t24 + t36) - g(2) * (t33 * t24 + t28 + (-pkin(3) * t31 - qJ(2)) * t34) - g(3) * (-pkin(5) + t37);];
U_reg  = t1;
