% Calculate minimal parameter regressor of potential energy for
% S5RRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRRP3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:25
% EndTime: 2019-12-31 21:49:25
% DurationCPUTime: 0.04s
% Computational Cost: add. (81->27), mult. (55->33), div. (0->0), fcn. (52->8), ass. (0->16)
t111 = qJ(1) + qJ(2);
t110 = qJ(3) + t111;
t106 = sin(t110);
t107 = cos(t110);
t117 = g(1) * t107 + g(2) * t106;
t112 = sin(qJ(4));
t114 = cos(qJ(4));
t116 = pkin(4) * t114 + qJ(5) * t112 + pkin(3);
t115 = cos(qJ(1));
t113 = sin(qJ(1));
t109 = cos(t111);
t108 = sin(t111);
t105 = g(1) * t106 - g(2) * t107;
t104 = -g(3) * t112 - t117 * t114;
t103 = -g(3) * t114 + t117 * t112;
t1 = [0, -g(1) * t115 - g(2) * t113, g(1) * t113 - g(2) * t115, 0, -g(1) * t109 - g(2) * t108, g(1) * t108 - g(2) * t109, 0, -t117, t105, 0, 0, 0, 0, 0, t104, t103, t104, -t105, -t103, -g(1) * (t115 * pkin(1) + pkin(2) * t109) - g(2) * (t113 * pkin(1) + pkin(2) * t108) - g(3) * (t112 * pkin(4) - t114 * qJ(5) + pkin(5) + pkin(6) + pkin(7)) + (g(2) * pkin(8) - g(1) * t116) * t107 + (-g(1) * pkin(8) - g(2) * t116) * t106;];
U_reg = t1;
