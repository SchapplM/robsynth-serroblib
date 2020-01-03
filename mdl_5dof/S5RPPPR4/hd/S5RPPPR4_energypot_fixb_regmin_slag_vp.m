% Calculate minimal parameter regressor of potential energy for
% S5RPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPPPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:15
% EndTime: 2019-12-31 17:45:15
% DurationCPUTime: 0.04s
% Computational Cost: add. (63->24), mult. (55->31), div. (0->0), fcn. (46->8), ass. (0->18)
t110 = qJ(2) + pkin(5);
t117 = g(3) * t110;
t107 = qJ(1) + pkin(7);
t101 = sin(t107);
t103 = cos(t107);
t112 = cos(qJ(1));
t116 = t112 * pkin(1) + t103 * pkin(2) + t101 * qJ(3);
t111 = sin(qJ(1));
t115 = t111 * pkin(1) + t101 * pkin(2) - t103 * qJ(3);
t114 = -g(1) * t101 + g(2) * t103;
t113 = -g(1) * t112 - g(2) * t111;
t109 = cos(pkin(8));
t108 = sin(pkin(8));
t106 = pkin(8) + qJ(5);
t102 = cos(t106);
t100 = sin(t106);
t96 = g(1) * t103 + g(2) * t101;
t1 = [0, t113, g(1) * t111 - g(2) * t112, t113 * pkin(1) - t117, t96, t114, -g(1) * t116 - g(2) * t115 - t117, -g(3) * t109 + t114 * t108, g(3) * t108 + t114 * t109, -t96, -g(1) * (t103 * qJ(4) + t116) - g(2) * (t101 * qJ(4) + t115) - g(3) * (pkin(3) + t110), 0, 0, 0, 0, 0, -g(3) * t102 + t114 * t100, g(3) * t100 + t114 * t102;];
U_reg = t1;
