% Calculate minimal parameter regressor of potential energy for
% S5RPRRP8
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
% U_reg [1x20]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRRP8_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_energypot_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:29
% EndTime: 2019-12-31 18:47:29
% DurationCPUTime: 0.05s
% Computational Cost: add. (58->28), mult. (108->38), div. (0->0), fcn. (124->6), ass. (0->18)
t122 = cos(qJ(3));
t121 = sin(qJ(1));
t120 = sin(qJ(3));
t114 = cos(qJ(1));
t119 = t114 * pkin(1) + t121 * qJ(2);
t118 = t121 * pkin(1) - t114 * qJ(2);
t102 = -t114 * t122 - t121 * t120;
t103 = t114 * t120 - t121 * t122;
t117 = g(1) * t103 - g(2) * t102;
t116 = g(1) * t102 + g(2) * t103;
t112 = sin(qJ(4));
t113 = cos(qJ(4));
t115 = pkin(4) * t113 + qJ(5) * t112 + pkin(3);
t105 = -g(1) * t114 - g(2) * t121;
t104 = g(1) * t121 - g(2) * t114;
t101 = g(3) * t112 + t116 * t113;
t100 = g(3) * t113 - t116 * t112;
t1 = [0, t105, t104, t105, -t104, -g(3) * pkin(5) - g(1) * t119 - g(2) * t118, 0, t116, t117, 0, 0, 0, 0, 0, t101, t100, t101, -t117, -t100, -g(1) * (t114 * pkin(2) + t119) - g(2) * (t121 * pkin(2) + t118) - g(3) * (-t112 * pkin(4) + t113 * qJ(5) + pkin(5) - pkin(6)) + (-g(1) * pkin(7) + g(2) * t115) * t103 + (g(2) * pkin(7) + g(1) * t115) * t102;];
U_reg = t1;
