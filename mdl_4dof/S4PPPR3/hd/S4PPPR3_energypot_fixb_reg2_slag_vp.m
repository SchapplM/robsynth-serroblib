% Calculate inertial parameters regressor of potential energy for
% S4PPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta3]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PPPR3_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR3_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPPR3_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR3_energypot_fixb_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:11:16
% EndTime: 2019-03-08 18:11:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (26->21), mult. (20->18), div. (0->0), fcn. (10->4), ass. (0->10)
t19 = g(1) * qJ(2);
t18 = g(2) * qJ(1);
t17 = pkin(2) + qJ(1);
t16 = -qJ(3) - pkin(1);
t15 = cos(pkin(5));
t14 = sin(pkin(5));
t13 = pkin(5) + qJ(4);
t12 = cos(t13);
t11 = sin(t13);
t1 = [0, 0, 0, 0, 0, 0, g(3), g(1), -g(2), -t18, 0, 0, 0, 0, 0, 0, -g(2), -g(3), -g(1), g(3) * pkin(1) - t18 - t19, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t15, -g(1) * t15 + g(2) * t14, g(3), -g(2) * t17 - g(3) * t16 - t19, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t12, -g(1) * t12 + g(2) * t11, g(3), -g(1) * (t14 * pkin(3) + qJ(2)) - g(2) * (t15 * pkin(3) + t17) - g(3) * (-pkin(4) + t16);];
U_reg  = t1;
