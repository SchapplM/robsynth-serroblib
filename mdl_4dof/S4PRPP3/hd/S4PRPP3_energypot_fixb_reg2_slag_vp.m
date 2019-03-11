% Calculate inertial parameters regressor of potential energy for
% S4PRPP3
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S4PRPP3_energypot_fixb_reg2_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:19:53
% EndTime: 2019-03-08 18:19:53
% DurationCPUTime: 0.06s
% Computational Cost: add. (28->22), mult. (32->18), div. (0->0), fcn. (22->2), ass. (0->9)
t23 = g(3) * pkin(4);
t22 = g(2) * qJ(1);
t18 = sin(qJ(2));
t19 = cos(qJ(2));
t21 = t19 * pkin(2) + t18 * qJ(3) + pkin(1);
t20 = t18 * pkin(2) - t19 * qJ(3) + qJ(1);
t13 = -g(1) * t19 - g(2) * t18;
t12 = g(1) * t18 - g(2) * t19;
t1 = [0, 0, 0, 0, 0, 0, -g(1), g(3), -g(2), -t22, 0, 0, 0, 0, 0, 0, t13, t12, -g(3), -g(1) * pkin(1) - t22 - t23, 0, 0, 0, 0, 0, 0, t13, -g(3), -t12, -g(1) * t21 - g(2) * t20 - t23, 0, 0, 0, 0, 0, 0, t13, -t12, g(3), -g(1) * (t19 * pkin(3) + t21) - g(2) * (t18 * pkin(3) + t20) - g(3) * (-qJ(4) + pkin(4));];
U_reg  = t1;
