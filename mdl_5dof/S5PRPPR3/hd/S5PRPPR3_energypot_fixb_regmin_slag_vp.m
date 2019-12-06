% Calculate minimal parameter regressor of potential energy for
% S5PRPPR3
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
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPPR3_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR3_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:54
% EndTime: 2019-12-05 15:26:54
% DurationCPUTime: 0.06s
% Computational Cost: add. (59->29), mult. (71->43), div. (0->0), fcn. (68->8), ass. (0->22)
t100 = qJ(2) + pkin(8);
t98 = cos(t100);
t117 = g(3) * t98;
t101 = sin(pkin(7));
t102 = cos(pkin(7));
t103 = -qJ(3) - pkin(5);
t107 = cos(qJ(2));
t96 = t107 * pkin(2) + pkin(1);
t116 = t101 * t96 + t102 * t103;
t105 = sin(qJ(2));
t115 = t105 * pkin(2) + qJ(1);
t104 = sin(qJ(5));
t114 = t101 * t104;
t106 = cos(qJ(5));
t113 = t101 * t106;
t112 = t102 * t104;
t111 = t102 * t106;
t110 = -t101 * t103 + t102 * t96;
t97 = sin(t100);
t109 = pkin(3) * t98 + qJ(4) * t97;
t108 = g(1) * t102 + g(2) * t101;
t1 = [-g(3) * qJ(1), 0, -g(3) * t105 - t108 * t107, -g(3) * t107 + t108 * t105, -g(1) * t110 - g(2) * t116 - g(3) * t115, g(3) * t97 + t108 * t98, -t108 * t97 + t117, -g(1) * (t109 * t102 + t110) - g(2) * (t109 * t101 + t116) - g(3) * (t97 * pkin(3) - t98 * qJ(4) + t115), 0, 0, 0, 0, 0, -g(1) * (t97 * t112 + t113) - g(2) * (t97 * t114 - t111) + t104 * t117, -g(1) * (t97 * t111 - t114) - g(2) * (t97 * t113 + t112) + t106 * t117;];
U_reg = t1;
