% Calculate minimal parameter regressor of potential energy for
% S5PRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% U_reg [1x19]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PRPRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:39
% EndTime: 2019-12-05 15:28:39
% DurationCPUTime: 0.06s
% Computational Cost: add. (92->33), mult. (70->37), div. (0->0), fcn. (63->8), ass. (0->18)
t117 = pkin(5) + qJ(1);
t111 = pkin(7) + qJ(2);
t105 = sin(t111);
t107 = cos(t111);
t116 = g(1) * t107 + g(2) * t105;
t110 = pkin(8) + qJ(4);
t104 = sin(t110);
t106 = cos(t110);
t113 = cos(pkin(8));
t115 = t113 * pkin(3) + pkin(4) * t106 + qJ(5) * t104 + pkin(2);
t114 = -pkin(6) - qJ(3);
t112 = sin(pkin(8));
t109 = cos(pkin(7)) * pkin(1);
t108 = sin(pkin(7)) * pkin(1);
t101 = g(1) * t105 - g(2) * t107;
t100 = -g(3) * t104 - t116 * t106;
t99 = -g(3) * t106 + t116 * t104;
t1 = [-g(3) * qJ(1), 0, -t116, t101, -g(3) * t112 - t116 * t113, -g(3) * t113 + t116 * t112, -t101, -g(1) * (t107 * pkin(2) + t105 * qJ(3) + t109) - g(2) * (t105 * pkin(2) - t107 * qJ(3) + t108) - g(3) * t117, 0, 0, 0, 0, 0, t100, t99, t100, -t101, -t99, -g(1) * t109 - g(2) * t108 - g(3) * (t112 * pkin(3) + t104 * pkin(4) - t106 * qJ(5) + t117) + (-g(1) * t115 - g(2) * t114) * t107 + (g(1) * t114 - g(2) * t115) * t105;];
U_reg = t1;
