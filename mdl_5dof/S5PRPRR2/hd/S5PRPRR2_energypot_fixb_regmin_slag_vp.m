% Calculate minimal parameter regressor of potential energy for
% S5PRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% U_reg [1x15]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRR2_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR2_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR2_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR2_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:45:12
% EndTime: 2019-12-05 15:45:12
% DurationCPUTime: 0.05s
% Computational Cost: add. (50->22), mult. (51->34), div. (0->0), fcn. (51->8), ass. (0->18)
t103 = qJ(2) + pkin(9) + qJ(4);
t100 = sin(t103);
t116 = g(3) * t100;
t104 = sin(pkin(8));
t107 = sin(qJ(5));
t115 = t104 * t107;
t109 = cos(qJ(5));
t114 = t104 * t109;
t105 = cos(pkin(8));
t113 = t105 * t107;
t112 = t105 * t109;
t111 = g(1) * t105 + g(2) * t104;
t110 = cos(qJ(2));
t108 = sin(qJ(2));
t106 = -qJ(3) - pkin(5);
t102 = t110 * pkin(2) + pkin(1);
t101 = cos(t103);
t1 = [-g(3) * qJ(1), 0, -g(3) * t108 - t111 * t110, -g(3) * t110 + t111 * t108, -g(1) * (t105 * t102 - t104 * t106) - g(2) * (t104 * t102 + t105 * t106) - g(3) * (t108 * pkin(2) + qJ(1)), 0, -t111 * t101 - t116, -g(3) * t101 + t111 * t100, 0, 0, 0, 0, 0, -g(1) * (t101 * t112 + t115) - g(2) * (t101 * t114 - t113) - t109 * t116, -g(1) * (-t101 * t113 + t114) - g(2) * (-t101 * t115 - t112) + t107 * t116;];
U_reg = t1;
