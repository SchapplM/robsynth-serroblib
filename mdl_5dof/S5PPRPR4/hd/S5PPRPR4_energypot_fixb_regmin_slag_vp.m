% Calculate minimal parameter regressor of potential energy for
% S5PPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
% 
% Output:
% U_reg [1x16]
%   minimal parameter regressor of Potential energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function U_reg = S5PPRPR4_energypot_fixb_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_energypot_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_energypot_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_energypot_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From energy_potential_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:23
% EndTime: 2019-12-31 17:32:23
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->24), mult. (83->32), div. (0->0), fcn. (94->8), ass. (0->17)
t112 = cos(qJ(3));
t111 = sin(qJ(3));
t110 = g(3) * qJ(1);
t107 = sin(pkin(7));
t108 = cos(pkin(7));
t109 = t108 * pkin(1) + t107 * qJ(2);
t91 = -t107 * t111 - t108 * t112;
t92 = -t107 * t112 + t108 * t111;
t106 = g(1) * t92 - g(2) * t91;
t105 = g(1) * t91 + g(2) * t92;
t104 = t107 * pkin(1) - t108 * qJ(2);
t103 = cos(pkin(8));
t102 = sin(pkin(8));
t101 = pkin(8) + qJ(5);
t97 = cos(t101);
t96 = sin(t101);
t1 = [-t110, -g(1) * t109 - g(2) * t104 - t110, 0, t105, t106, g(3) * t102 + t105 * t103, g(3) * t103 - t105 * t102, -t106, -g(1) * (t108 * pkin(2) - t91 * pkin(3) + t92 * qJ(4) + t109) - g(2) * (t107 * pkin(2) - t92 * pkin(3) - t91 * qJ(4) + t104) - g(3) * (-pkin(5) + qJ(1)), 0, 0, 0, 0, 0, g(3) * t96 + t105 * t97, g(3) * t97 - t105 * t96;];
U_reg = t1;
