% Calculate minimal parameter regressor of potential energy for
% S5PRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:55
% EndTime: 2019-12-31 17:36:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (79->28), mult. (83->35), div. (0->0), fcn. (82->8), ass. (0->18)
t112 = pkin(7) + qJ(2);
t108 = sin(t112);
t109 = cos(t112);
t121 = g(1) * t109 + g(2) * t108;
t123 = pkin(5) + qJ(1);
t122 = t109 * pkin(2) + t108 * qJ(3) + cos(pkin(7)) * pkin(1);
t120 = t108 * pkin(2) - t109 * qJ(3) + sin(pkin(7)) * pkin(1);
t113 = sin(pkin(8));
t114 = cos(pkin(8));
t119 = pkin(3) * t114 + qJ(4) * t113;
t115 = sin(qJ(5));
t116 = cos(qJ(5));
t118 = t113 * t116 - t114 * t115;
t117 = t113 * t115 + t114 * t116;
t103 = g(1) * t108 - g(2) * t109;
t102 = -g(3) * t113 - t121 * t114;
t101 = -g(3) * t114 + t121 * t113;
t1 = [-g(3) * qJ(1), 0, -t121, t103, t102, t101, -t103, -g(1) * t122 - g(2) * t120 - g(3) * t123, t102, -t103, -t101, -g(1) * (t119 * t109 + t122) - g(2) * (t119 * t108 + t120) - g(3) * (t113 * pkin(3) - t114 * qJ(4) + t123), 0, 0, 0, 0, 0, -g(3) * t118 - t121 * t117, g(3) * t117 - t121 * t118;];
U_reg = t1;
