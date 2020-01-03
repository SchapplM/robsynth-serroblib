% Calculate minimal parameter regressor of potential energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
% 
% Output:
% U_reg [1x21]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5RPRPP5_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:41
% EndTime: 2019-12-31 18:16:41
% DurationCPUTime: 0.05s
% Computational Cost: add. (53->32), mult. (94->39), div. (0->0), fcn. (85->4), ass. (0->17)
t98 = sin(qJ(3));
t99 = sin(qJ(1));
t108 = t98 * t99;
t101 = cos(qJ(1));
t107 = t101 * pkin(1) + t99 * qJ(2);
t100 = cos(qJ(3));
t106 = qJ(4) * t100;
t94 = t99 * pkin(1);
t105 = t99 * pkin(6) + t101 * t106 + t94;
t104 = t100 * pkin(3) + t98 * qJ(4) + pkin(2) + pkin(5);
t103 = -pkin(3) * t98 - qJ(2);
t87 = g(1) * t99 - g(2) * t101;
t102 = pkin(3) * t108 + t101 * pkin(6) - t99 * t106 + t107;
t88 = g(1) * t101 + g(2) * t99;
t86 = -g(3) * t98 + t87 * t100;
t85 = -g(3) * t100 - t87 * t98;
t1 = [0, -t88, t87, t88, -t87, -g(1) * t107 - g(2) * (-t101 * qJ(2) + t94) - g(3) * pkin(5), 0, 0, 0, 0, 0, t85, -t86, t85, -t88, t86, -g(1) * t102 - g(2) * (t103 * t101 + t105) - g(3) * t104, t85, t86, t88, -g(1) * (pkin(4) * t108 + t102) - g(2) * (-t99 * qJ(5) + t105) - g(3) * (t100 * pkin(4) + t104) + (g(1) * qJ(5) - g(2) * (-pkin(4) * t98 + t103)) * t101;];
U_reg = t1;
