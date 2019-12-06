% Calculate minimal parameter regressor of potential energy for
% S5RPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% U_reg [1x22]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:57
% EndTime: 2019-12-05 17:59:57
% DurationCPUTime: 0.04s
% Computational Cost: add. (43->26), mult. (54->30), div. (0->0), fcn. (48->6), ass. (0->14)
t113 = sin(qJ(1));
t115 = cos(qJ(1));
t116 = t115 * pkin(1) + t113 * qJ(2);
t103 = g(1) * t113 - g(2) * t115;
t114 = cos(qJ(3));
t112 = sin(qJ(3));
t111 = qJ(3) + qJ(4);
t110 = -qJ(5) - pkin(7) - pkin(6);
t108 = t113 * pkin(1);
t106 = cos(t111);
t105 = sin(t111);
t104 = g(1) * t115 + g(2) * t113;
t102 = t112 * pkin(3) + pkin(4) * t105;
t1 = [0, -t104, t103, t104, -t103, -g(1) * t116 - g(2) * (-t115 * qJ(2) + t108) - g(3) * pkin(5), 0, 0, 0, 0, 0, -g(3) * t114 - t103 * t112, g(3) * t112 - t103 * t114, 0, 0, 0, 0, 0, -g(3) * t106 - t103 * t105, g(3) * t105 - t103 * t106, -t104, -g(1) * (t113 * t102 - t115 * t110 + t116) - g(2) * (-t113 * t110 + t108 + (-qJ(2) - t102) * t115) - g(3) * (t114 * pkin(3) + pkin(4) * t106 + pkin(2) + pkin(5));];
U_reg = t1;
