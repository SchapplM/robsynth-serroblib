% Calculate minimal parameter regressor of potential energy for
% S5RPRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% U_reg [1x17]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:26:04
% EndTime: 2019-12-31 18:26:04
% DurationCPUTime: 0.04s
% Computational Cost: add. (42->27), mult. (60->38), div. (0->0), fcn. (62->8), ass. (0->19)
t110 = sin(qJ(3));
t111 = sin(qJ(1));
t117 = t111 * t110;
t114 = cos(qJ(1));
t116 = t114 * pkin(1) + t111 * qJ(2);
t108 = qJ(3) + pkin(8);
t103 = sin(t108);
t104 = cos(t108);
t115 = g(1) * (t111 * t103 + t114 * t104) + g(2) * (-t114 * t103 + t111 * t104);
t113 = cos(qJ(3));
t112 = cos(qJ(5));
t109 = sin(qJ(5));
t106 = t111 * pkin(1);
t102 = t113 * pkin(3) + pkin(2);
t101 = -g(1) * t114 - g(2) * t111;
t100 = g(1) * t111 - g(2) * t114;
t99 = -t114 * t110 + t111 * t113;
t98 = -t114 * t113 - t117;
t1 = [0, t101, t100, t101, -t100, -g(1) * t116 - g(2) * (-t114 * qJ(2) + t106) - g(3) * pkin(5), 0, g(1) * t98 - g(2) * t99, -g(1) * t99 - g(2) * t98, -g(1) * (pkin(3) * t117 + t114 * t102 + t116) - g(2) * (t111 * t102 + t106 + (-pkin(3) * t110 - qJ(2)) * t114) - g(3) * (-qJ(4) - pkin(6) + pkin(5)), 0, 0, 0, 0, 0, g(3) * t109 - t115 * t112, g(3) * t112 + t115 * t109;];
U_reg = t1;
