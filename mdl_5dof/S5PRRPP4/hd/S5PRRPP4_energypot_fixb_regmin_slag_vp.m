% Calculate minimal parameter regressor of potential energy for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPP4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:41:11
% EndTime: 2019-12-31 17:41:12
% DurationCPUTime: 0.08s
% Computational Cost: add. (88->30), mult. (84->33), div. (0->0), fcn. (77->6), ass. (0->16)
t108 = sin(qJ(3));
t117 = qJ(4) * t108 + pkin(2);
t107 = pkin(7) + qJ(2);
t101 = sin(t107);
t109 = cos(qJ(3));
t115 = t101 * t109;
t102 = cos(t107);
t114 = t102 * t109;
t113 = sin(pkin(7)) * pkin(1) + pkin(3) * t115 + t117 * t101;
t112 = cos(pkin(7)) * pkin(1) + pkin(3) * t114 + t101 * pkin(6) + t117 * t102;
t111 = g(1) * t102 + g(2) * t101;
t110 = t108 * pkin(3) - t109 * qJ(4) + pkin(5) + qJ(1);
t93 = g(1) * t101 - g(2) * t102;
t92 = -g(3) * t108 - t111 * t109;
t91 = -g(3) * t109 + t111 * t108;
t1 = [-g(3) * qJ(1), 0, -t111, t93, 0, 0, 0, 0, 0, t92, t91, t92, -t93, -t91, -g(1) * t112 - g(2) * (-t102 * pkin(6) + t113) - g(3) * t110, t92, -t91, t93, -g(1) * (pkin(4) * t114 - t101 * qJ(5) + t112) - g(2) * (pkin(4) * t115 + (-pkin(6) + qJ(5)) * t102 + t113) - g(3) * (t108 * pkin(4) + t110);];
U_reg = t1;
