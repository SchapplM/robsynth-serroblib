% Calculate minimal parameter regressor of potential energy for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% U_reg [1x23]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RRPRR7_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR7_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:40
% EndTime: 2019-12-31 20:15:40
% DurationCPUTime: 0.03s
% Computational Cost: add. (47->19), mult. (41->25), div. (0->0), fcn. (38->8), ass. (0->13)
t106 = qJ(1) + qJ(2);
t102 = sin(t106);
t104 = cos(t106);
t99 = g(1) * t102 - g(2) * t104;
t110 = cos(qJ(1));
t109 = cos(qJ(4));
t108 = sin(qJ(1));
t107 = sin(qJ(4));
t105 = qJ(4) + qJ(5);
t103 = cos(t105);
t101 = sin(t105);
t100 = g(1) * t104 + g(2) * t102;
t1 = [0, -g(1) * t110 - g(2) * t108, g(1) * t108 - g(2) * t110, 0, -t100, t99, t100, -t99, -g(1) * (t110 * pkin(1) + t104 * pkin(2) + t102 * qJ(3)) - g(2) * (t108 * pkin(1) + t102 * pkin(2) - t104 * qJ(3)) - g(3) * (pkin(6) + pkin(5)), 0, 0, 0, 0, 0, -g(3) * t109 - t99 * t107, g(3) * t107 - t99 * t109, 0, 0, 0, 0, 0, -g(3) * t103 - t99 * t101, g(3) * t101 - t99 * t103;];
U_reg = t1;
