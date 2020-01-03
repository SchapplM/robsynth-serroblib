% Calculate minimal parameter regressor of potential energy for
% S5PRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRRPR8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR8_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:42:41
% EndTime: 2019-12-31 17:42:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->25), mult. (54->37), div. (0->0), fcn. (54->10), ass. (0->20)
t105 = qJ(2) + qJ(3);
t101 = pkin(9) + t105;
t117 = g(3) * sin(t101);
t106 = sin(pkin(8));
t108 = sin(qJ(5));
t116 = t106 * t108;
t110 = cos(qJ(5));
t115 = t106 * t110;
t107 = cos(pkin(8));
t114 = t107 * t108;
t113 = t107 * t110;
t112 = g(1) * t107 + g(2) * t106;
t111 = cos(qJ(2));
t109 = sin(qJ(2));
t104 = -qJ(4) - pkin(6) - pkin(5);
t103 = cos(t105);
t102 = sin(t105);
t100 = cos(t101);
t98 = t111 * pkin(2) + pkin(3) * t103 + pkin(1);
t1 = [-g(3) * qJ(1), 0, -g(3) * t109 - t112 * t111, -g(3) * t111 + t112 * t109, 0, -g(3) * t102 - t112 * t103, -g(3) * t103 + t112 * t102, -g(1) * (-t106 * t104 + t107 * t98) - g(2) * (t107 * t104 + t106 * t98) - g(3) * (t109 * pkin(2) + pkin(3) * t102 + qJ(1)), 0, 0, 0, 0, 0, -g(1) * (t100 * t113 + t116) - g(2) * (t100 * t115 - t114) - t110 * t117, -g(1) * (-t100 * t114 + t115) - g(2) * (-t100 * t116 - t113) + t108 * t117;];
U_reg = t1;
