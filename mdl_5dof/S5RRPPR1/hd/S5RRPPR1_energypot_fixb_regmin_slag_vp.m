% Calculate minimal parameter regressor of potential energy for
% S5RRPPR1
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
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPPR1_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energypot_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:55:59
% EndTime: 2020-01-03 11:55:59
% DurationCPUTime: 0.04s
% Computational Cost: add. (71->25), mult. (48->33), div. (0->0), fcn. (42->10), ass. (0->18)
t128 = g(1) * (qJ(3) + pkin(6) + pkin(5));
t120 = qJ(1) + qJ(2);
t115 = sin(t120);
t123 = sin(qJ(1));
t127 = t123 * pkin(1) + pkin(2) * t115;
t116 = cos(t120);
t124 = cos(qJ(1));
t126 = -t124 * pkin(1) - pkin(2) * t116;
t114 = pkin(8) + t120;
t109 = sin(t114);
t110 = cos(t114);
t125 = g(2) * t109 - g(3) * t110;
t122 = cos(pkin(9));
t121 = sin(pkin(9));
t119 = pkin(9) + qJ(5);
t113 = cos(t119);
t112 = sin(t119);
t1 = [0, -g(2) * t123 + g(3) * t124, -g(2) * t124 - g(3) * t123, 0, -g(2) * t115 + g(3) * t116, -g(2) * t116 - g(3) * t115, -g(2) * t127 - g(3) * t126 - t128, -g(1) * t121 - t125 * t122, -g(1) * t122 + t125 * t121, g(2) * t110 + g(3) * t109, -t128 - g(2) * (t109 * pkin(3) - t110 * qJ(4) + t127) - g(3) * (-t110 * pkin(3) - t109 * qJ(4) + t126), 0, 0, 0, 0, 0, -g(1) * t112 - t125 * t113, -g(1) * t113 + t125 * t112;];
U_reg = t1;
