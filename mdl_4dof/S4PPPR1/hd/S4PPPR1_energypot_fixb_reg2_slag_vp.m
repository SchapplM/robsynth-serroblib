% Calculate inertial parameters regressor of potential energy for
% S4PPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PPPR1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:14
% EndTime: 2018-11-14 13:38:14
% DurationCPUTime: 0.04s
% Computational Cost: add. (37->30), mult. (52->29), div. (0->0), fcn. (46->4), ass. (0->16)
t26 = g(3) * qJ(1);
t25 = pkin(2) + qJ(1);
t18 = sin(pkin(5));
t19 = cos(pkin(5));
t24 = t19 * pkin(1) + t18 * qJ(2);
t23 = t19 * qJ(3) + t24;
t15 = t18 * pkin(1);
t22 = -t19 * qJ(2) + t15;
t21 = cos(qJ(4));
t20 = sin(qJ(4));
t12 = t18 * qJ(3);
t11 = g(1) * t19 + g(2) * t18;
t10 = g(1) * t18 - g(2) * t19;
t9 = t18 * t21 + t19 * t20;
t8 = -t18 * t20 + t19 * t21;
t1 = [0, 0, 0, 0, 0, 0, -t11, t10, -g(3), -t26, 0, 0, 0, 0, 0, 0, -g(3), t11, -t10, -g(1) * t24 - g(2) * t22 - t26, 0, 0, 0, 0, 0, 0, -t10, g(3), -t11, -g(1) * t23 - g(2) * (t12 + t22) - g(3) * t25, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t8, -g(1) * t8 - g(2) * t9, -g(3), -g(1) * (t18 * pkin(3) + t23) - g(2) * (t12 + t15 + (-pkin(3) - qJ(2)) * t19) - g(3) * (pkin(4) + t25);];
U_reg  = t1;
