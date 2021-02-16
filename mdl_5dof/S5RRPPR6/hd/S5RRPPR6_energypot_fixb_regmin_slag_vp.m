% Calculate minimal parameter regressor of potential energy for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% U_reg [1x25]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:47:31
% EndTime: 2021-01-15 19:47:31
% DurationCPUTime: 0.08s
% Computational Cost: add. (90->45), mult. (103->64), div. (0->0), fcn. (105->12), ass. (0->31)
t98 = qJ(2) + pkin(8);
t93 = sin(t98);
t116 = g(3) * t93;
t104 = sin(qJ(2));
t115 = t104 * pkin(2) + pkin(5);
t105 = sin(qJ(1));
t95 = cos(t98);
t114 = t105 * t95;
t99 = sin(pkin(9));
t113 = t105 * t99;
t107 = cos(qJ(1));
t112 = t107 * t95;
t111 = t107 * t99;
t101 = cos(pkin(9));
t110 = t105 * t101;
t109 = t107 * t101;
t108 = g(1) * t107 + g(2) * t105;
t106 = cos(qJ(2));
t103 = pkin(6) + qJ(3);
t102 = cos(pkin(8));
t100 = sin(pkin(8));
t97 = pkin(9) + qJ(5);
t94 = cos(t97);
t92 = sin(t97);
t91 = t106 * pkin(2) + pkin(1);
t90 = t103 * t105;
t89 = t107 * t103;
t88 = g(1) * t105 - g(2) * t107;
t87 = -g(3) * t95 + t108 * t93;
t86 = (pkin(3) * t102 + qJ(4) * t100 + pkin(2)) * t106 + (-t100 * pkin(3) + qJ(4) * t102) * t104 + pkin(1);
t1 = [0, -t108, t88, 0, 0, 0, 0, 0, -g(3) * t104 - t108 * t106, -g(3) * t106 + t108 * t104, -t108 * t95 - t116, t87, -t88, -g(1) * (t107 * t91 + t90) - g(2) * (t105 * t91 - t89) - g(3) * t115, -g(1) * (t95 * t109 + t113) - g(2) * (t95 * t110 - t111) - t101 * t116, -g(1) * (-t95 * t111 + t110) - g(2) * (-t95 * t113 - t109) + t99 * t116, -t87, -g(1) * (t86 * t107 + t90) - g(2) * (t86 * t105 - t89) - g(3) * (t93 * pkin(3) - t95 * qJ(4) + t115), 0, 0, 0, 0, 0, -g(1) * (t105 * t92 + t94 * t112) - g(2) * (-t107 * t92 + t94 * t114) - t94 * t116, -g(1) * (t105 * t94 - t92 * t112) - g(2) * (-t107 * t94 - t92 * t114) + t92 * t116;];
U_reg = t1;
