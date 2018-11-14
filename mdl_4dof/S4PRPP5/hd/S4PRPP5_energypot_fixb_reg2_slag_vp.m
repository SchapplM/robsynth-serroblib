% Calculate inertial parameters regressor of potential energy for
% S4PRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:11
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP5_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP5_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP5_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP5_energypot_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:10:22
% EndTime: 2018-11-14 14:10:22
% DurationCPUTime: 0.06s
% Computational Cost: add. (28->21), mult. (32->18), div. (0->0), fcn. (22->2), ass. (0->9)
t22 = g(1) * qJ(1);
t17 = sin(qJ(2));
t18 = cos(qJ(2));
t21 = t18 * pkin(2) + t17 * qJ(3) + pkin(1);
t20 = t17 * pkin(2) - t18 * qJ(3) + qJ(1);
t19 = g(3) * pkin(4);
t13 = g(1) * t18 - g(2) * t17;
t12 = -g(1) * t17 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, -g(2), -g(3), -g(1), -t22, 0, 0, 0, 0, 0, 0, t12, -t13, g(3), -g(2) * pkin(1) + t19 - t22, 0, 0, 0, 0, 0, 0, t12, g(3), t13, -g(1) * t20 - g(2) * t21 + t19, 0, 0, 0, 0, 0, 0, t12, t13, -g(3), -g(1) * (t17 * pkin(3) + t20) - g(2) * (t18 * pkin(3) + t21) - g(3) * (qJ(4) - pkin(4));];
U_reg  = t1;
