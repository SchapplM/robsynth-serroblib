% Calculate inertial parameters regressor of potential energy for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRRP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:43:37
% EndTime: 2018-11-14 13:43:37
% DurationCPUTime: 0.04s
% Computational Cost: add. (64->27), mult. (40->28), div. (0->0), fcn. (30->6), ass. (0->16)
t31 = pkin(4) + qJ(1);
t32 = g(3) * (pkin(5) + t31);
t25 = pkin(6) + qJ(2);
t19 = sin(t25);
t26 = sin(pkin(6));
t30 = t26 * pkin(1) + pkin(2) * t19;
t20 = cos(t25);
t27 = cos(pkin(6));
t29 = t27 * pkin(1) + pkin(2) * t20;
t28 = -g(1) * t27 - g(2) * t26;
t21 = qJ(3) + t25;
t18 = cos(t21);
t17 = sin(t21);
t14 = -g(1) * t18 - g(2) * t17;
t13 = g(1) * t17 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t28, g(1) * t26 - g(2) * t27, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -g(1) * t20 - g(2) * t19, g(1) * t19 - g(2) * t20, -g(3), t28 * pkin(1) - g(3) * t31, 0, 0, 0, 0, 0, 0, t14, t13, -g(3), -g(1) * t29 - g(2) * t30 - t32, 0, 0, 0, 0, 0, 0, t14, -g(3), -t13, -g(1) * (t18 * pkin(3) + t17 * qJ(4) + t29) - g(2) * (t17 * pkin(3) - t18 * qJ(4) + t30) - t32;];
U_reg  = t1;
