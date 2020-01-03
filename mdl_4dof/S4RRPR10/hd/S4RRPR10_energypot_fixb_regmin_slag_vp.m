% Calculate minimal parameter regressor of potential energy for
% S4RRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RRPR10_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR10_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR10_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR10_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:11:55
% EndTime: 2019-12-31 17:11:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (29->24), mult. (63->33), div. (0->0), fcn. (64->6), ass. (0->17)
t110 = cos(qJ(2));
t118 = g(3) * t110;
t106 = sin(qJ(4));
t108 = sin(qJ(1));
t117 = t108 * t106;
t109 = cos(qJ(4));
t116 = t108 * t109;
t111 = cos(qJ(1));
t115 = t111 * t106;
t114 = t111 * t109;
t113 = g(1) * t111 + g(2) * t108;
t107 = sin(qJ(2));
t112 = pkin(2) * t110 + qJ(3) * t107 + pkin(1);
t105 = g(1) * t108 - g(2) * t111;
t104 = g(3) * t107 + t113 * t110;
t103 = t113 * t107 - t118;
t1 = [0, -t113, t105, 0, 0, 0, 0, 0, -t104, t103, -t105, t104, -t103, -g(3) * (t107 * pkin(2) - t110 * qJ(3) + pkin(4)) + (g(2) * pkin(5) - g(1) * t112) * t111 + (-g(1) * pkin(5) - g(2) * t112) * t108, 0, 0, 0, 0, 0, -g(1) * (t107 * t115 + t116) - g(2) * (t107 * t117 - t114) + t106 * t118, -g(1) * (t107 * t114 - t117) - g(2) * (t107 * t116 + t115) + t109 * t118;];
U_reg = t1;
