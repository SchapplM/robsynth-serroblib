% Calculate minimal parameter regressor of potential energy for
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% U_reg [1x18]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S4RPPR6_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_energypot_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPPR6_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_energypot_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:40:43
% EndTime: 2019-12-31 16:40:43
% DurationCPUTime: 0.08s
% Computational Cost: add. (36->23), mult. (78->32), div. (0->0), fcn. (78->6), ass. (0->16)
t104 = sin(qJ(1));
t106 = cos(qJ(1));
t110 = g(1) * t106 + g(2) * t104;
t112 = t106 * pkin(1) + t104 * qJ(2);
t111 = t104 * pkin(1) - t106 * qJ(2);
t101 = sin(pkin(6));
t102 = cos(pkin(6));
t109 = pkin(2) * t102 + qJ(3) * t101;
t103 = sin(qJ(4));
t105 = cos(qJ(4));
t108 = t101 * t105 - t102 * t103;
t107 = t101 * t103 + t102 * t105;
t96 = g(1) * t104 - g(2) * t106;
t95 = -g(3) * t101 - t110 * t102;
t94 = -g(3) * t102 + t110 * t101;
t1 = [0, -t110, t96, t95, t94, -t96, -g(3) * pkin(4) - g(1) * t112 - g(2) * t111, t95, -t96, -t94, -g(1) * (t109 * t106 + t112) - g(2) * (t109 * t104 + t111) - g(3) * (t101 * pkin(2) - t102 * qJ(3) + pkin(4)), 0, 0, 0, 0, 0, -g(3) * t108 - t110 * t107, g(3) * t107 - t110 * t108;];
U_reg = t1;
