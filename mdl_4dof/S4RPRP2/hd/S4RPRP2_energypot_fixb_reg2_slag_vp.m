% Calculate inertial parameters regressor of potential energy for
% S4RPRP2
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPRP2_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP2_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:52
% EndTime: 2019-03-08 18:30:52
% DurationCPUTime: 0.04s
% Computational Cost: add. (42->28), mult. (64->33), div. (0->0), fcn. (62->4), ass. (0->18)
t34 = g(3) * pkin(4);
t33 = -pkin(5) + pkin(4);
t26 = sin(qJ(3));
t27 = sin(qJ(1));
t32 = t27 * t26;
t29 = cos(qJ(1));
t31 = t29 * pkin(1) + t27 * qJ(2);
t24 = t27 * pkin(1);
t30 = -t29 * qJ(2) + t24;
t28 = cos(qJ(3));
t22 = t28 * pkin(3) + pkin(2);
t21 = -g(1) * t29 - g(2) * t27;
t20 = g(1) * t27 - g(2) * t29;
t19 = -t29 * t26 + t27 * t28;
t18 = -t29 * t28 - t32;
t17 = g(1) * t18 - g(2) * t19;
t16 = -g(1) * t19 - g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t21, t20, -g(3), -t34, 0, 0, 0, 0, 0, 0, t21, -g(3), -t20, -g(1) * t31 - g(2) * t30 - t34, 0, 0, 0, 0, 0, 0, t17, t16, g(3), -g(1) * (t29 * pkin(2) + t31) - g(2) * (t27 * pkin(2) + t30) - g(3) * t33, 0, 0, 0, 0, 0, 0, t17, t16, g(3), -g(1) * (pkin(3) * t32 + t29 * t22 + t31) - g(2) * (t27 * t22 + t24 + (-pkin(3) * t26 - qJ(2)) * t29) - g(3) * (-qJ(4) + t33);];
U_reg  = t1;
