% Calculate inertial parameters regressor of potential energy for
% S4PRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPP1_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP1_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPP1_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP1_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:17:55
% EndTime: 2019-03-08 18:17:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (55->27), mult. (42->24), div. (0->0), fcn. (32->4), ass. (0->13)
t20 = pkin(4) + qJ(1);
t24 = g(3) * t20;
t17 = pkin(5) + qJ(2);
t13 = sin(t17);
t14 = cos(t17);
t19 = cos(pkin(5));
t23 = t19 * pkin(1) + t14 * pkin(2) + t13 * qJ(3);
t18 = sin(pkin(5));
t22 = t18 * pkin(1) + t13 * pkin(2) - t14 * qJ(3);
t21 = -g(1) * t19 - g(2) * t18;
t8 = g(1) * t14 + g(2) * t13;
t7 = g(1) * t13 - g(2) * t14;
t1 = [0, 0, 0, 0, 0, 0, t21, g(1) * t18 - g(2) * t19, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t8, t7, -g(3), t21 * pkin(1) - t24, 0, 0, 0, 0, 0, 0, -g(3), t8, -t7, -g(1) * t23 - g(2) * t22 - t24, 0, 0, 0, 0, 0, 0, -g(3), -t7, -t8, -g(1) * (t14 * qJ(4) + t23) - g(2) * (t13 * qJ(4) + t22) - g(3) * (pkin(3) + t20);];
U_reg  = t1;
