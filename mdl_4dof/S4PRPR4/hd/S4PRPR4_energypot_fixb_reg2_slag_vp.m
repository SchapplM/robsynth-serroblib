% Calculate inertial parameters regressor of potential energy for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% U_reg [1x(4*10)]
%   inertial parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4PRPR4_energypot_fixb_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energypot_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_energypot_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energypot_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:59
% EndTime: 2019-12-31 16:21:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->27), mult. (50->28), div. (0->0), fcn. (40->6), ass. (0->15)
t25 = pkin(4) + qJ(1);
t31 = g(3) * t25;
t22 = pkin(6) + qJ(2);
t18 = sin(t22);
t19 = cos(t22);
t24 = cos(pkin(6));
t30 = t24 * pkin(1) + t19 * pkin(2) + t18 * qJ(3);
t23 = sin(pkin(6));
t29 = t23 * pkin(1) + t18 * pkin(2) - t19 * qJ(3);
t13 = g(1) * t18 - g(2) * t19;
t28 = -g(1) * t24 - g(2) * t23;
t27 = cos(qJ(4));
t26 = sin(qJ(4));
t14 = g(1) * t19 + g(2) * t18;
t1 = [0, 0, 0, 0, 0, 0, t28, g(1) * t23 - g(2) * t24, -g(3), -g(3) * qJ(1), 0, 0, 0, 0, 0, 0, -t14, t13, -g(3), t28 * pkin(1) - t31, 0, 0, 0, 0, 0, 0, -g(3), t14, -t13, -g(1) * t30 - g(2) * t29 - t31, 0, 0, 0, 0, 0, 0, -g(3) * t27 - t13 * t26, g(3) * t26 - t13 * t27, -t14, -g(1) * (t19 * pkin(5) + t30) - g(2) * (t18 * pkin(5) + t29) - g(3) * (pkin(3) + t25);];
U_reg = t1;
