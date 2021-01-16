% Calculate minimal parameter regressor of potential energy for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% U_reg [1x28]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRRPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_energypot_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:15
% EndTime: 2021-01-15 23:41:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (118->49), mult. (207->82), div. (0->0), fcn. (250->12), ass. (0->33)
t104 = qJ(3) + pkin(10);
t102 = sin(t104);
t103 = cos(t104);
t106 = cos(pkin(5));
t105 = sin(pkin(5));
t115 = cos(qJ(1));
t128 = t105 * t115;
t111 = sin(qJ(1));
t129 = t105 * t111;
t110 = sin(qJ(2));
t130 = t105 * t110;
t114 = cos(qJ(2));
t122 = t115 * t114;
t125 = t111 * t110;
t93 = t106 * t125 - t122;
t123 = t115 * t110;
t124 = t111 * t114;
t94 = -t106 * t123 - t124;
t132 = g(1) * (t102 * t129 - t93 * t103) - g(2) * (t102 * t128 + t94 * t103) + g(3) * (t106 * t102 + t103 * t130);
t109 = sin(qJ(3));
t127 = t109 * t105;
t113 = cos(qJ(3));
t126 = t110 * t113;
t121 = t109 * pkin(3) + pkin(7);
t101 = t113 * pkin(3) + pkin(2);
t107 = qJ(4) + pkin(8);
t120 = t101 * t110 - t107 * t114;
t117 = g(3) * t105 * t114 - g(1) * (t106 * t124 + t123) - g(2) * (-t106 * t122 + t125);
t112 = cos(qJ(5));
t108 = sin(qJ(5));
t92 = t101 * t114 + t107 * t110 + pkin(1);
t90 = t105 * t121 - t120 * t106;
t1 = [0, -g(1) * t115 - g(2) * t111, g(1) * t111 - g(2) * t115, 0, 0, 0, 0, 0, g(1) * t93 + g(2) * t94 - g(3) * t130, -t117, 0, 0, 0, 0, 0, -g(1) * ((-t106 * t126 + t127) * t111 + t113 * t122) - g(2) * (-t94 * t113 - t115 * t127) - g(3) * (t105 * t126 + t106 * t109), -g(1) * (t93 * t109 + t113 * t129) - g(2) * (t94 * t109 - t113 * t128) - g(3) * (t106 * t113 - t110 * t127), -t132, -g(1) * (t93 * t102 + t103 * t129) - g(2) * (t94 * t102 - t103 * t128) - g(3) * (-t102 * t130 + t106 * t103), t117, -g(1) * (t90 * t111 + t92 * t115) - g(2) * (t92 * t111 - t90 * t115) - g(3) * (t120 * t105 + t121 * t106 + pkin(6)), 0, 0, 0, 0, 0, t117 * t108 - t112 * t132, t132 * t108 + t117 * t112;];
U_reg = t1;
